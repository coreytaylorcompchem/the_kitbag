import os
import logging
import requests
from pathlib import Path
from contextlib import contextmanager

import modeller
import modeller.automodel
from Bio import pairwise2, SeqIO
from Bio.PDB import PDBParser
from Bio.SeqUtils import seq1

logger = logging.getLogger(__name__)

@contextmanager
def working_directory(path):
    """Context manager hack for changing the current working directory."""
    prev_cwd = Path.cwd()
    os.chdir(path)
    try:
        yield
    finally:
        os.chdir(prev_cwd)

class ModellerBackend:
    def __init__(self, pdb_path: str = None, n_loop_models: int = 2, loop_refinement: dict = None, output_dir: str = None):
        self.env = modeller.Environ()
        self.pdb_path = pdb_path
        self.n_loop_models = n_loop_models
        self.loop_refinement = loop_refinement or {}
        self.model = None

        self.output_dir = Path(output_dir or Path.cwd())
        self.output_dir.mkdir(parents=True, exist_ok=True)

        if self.pdb_path:
            self._load_model()

    def _load_model(self):
        pdb_path = Path(self.pdb_path)
        pdb_stem = pdb_path.stem
        self.env.io.atom_files_directory.append(str(pdb_path.parent.resolve()))
        self.model = modeller.Model(self.env, file=pdb_stem)

    def _fetch_uniprot_sequence(self, uniprot_id):
        url = f"https://www.uniprot.org/uniprot/{uniprot_id}.fasta"
        response = requests.get(url)
        if not response.ok:
            raise ValueError(f"Could not fetch UniProt FASTA for {uniprot_id}")
        lines = response.text.strip().split('\n')
        return ''.join(lines[1:])  # skips header

    def _get_first_residue_and_chain(self):
        with open(self.pdb_path) as f:
            for line in f:
                if line.startswith("ATOM"):
                    chain_id = line[21].strip() or 'A'
                    res_num = int(line[22:26].strip())
                    return res_num, chain_id
        return 1, 'A'

    def _extract_pdb_sequence(self, chain_id='A'):
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure("model", self.pdb_path)
        model = structure[0]

        for chain in model:
            if chain.id == chain_id:
                seq = []
                res_ids = []
                for res in chain:
                    if res.id[0] != ' ' or 'CA' not in res:
                        continue
                    try:
                        aa = seq1(res.resname.strip(), custom_map={"MSE": "M"})
                    except Exception:
                        aa = 'X'
                    seq.append(aa)
                    res_ids.append(res.id[1])
                return ''.join(seq), res_ids

        raise ValueError(f"Chain {chain_id} not found in PDB")

    def _align_sequences(self, pdb_seq, fasta_seq):
        alignments = pairwise2.align.globalms(pdb_seq, fasta_seq, 2, -1, -10, -0.5)
        aln_pdb, aln_fasta, *_ = alignments[0]

        # Trim leading/trailing gaps in the PDB sequence
        start = next(i for i, c in enumerate(aln_pdb) if c != '-')
        end = len(aln_pdb) - next(i for i, c in enumerate(reversed(aln_pdb)) if c != '-')  # inclusive of

        trimmed_pdb = aln_pdb[start:end]
        trimmed_fasta = aln_fasta[start:end]

        return trimmed_pdb, trimmed_fasta

    def _write_pir_alignment(self, ali_path, pdb_code, pdb_aln, fasta_aln, res_ids, chain_id):
        start_res = res_ids[0]
        end_res = res_ids[-1]

        with open(ali_path, 'w') as f:
            f.write(f">P1;{pdb_code}\n")
            f.write(f"structureX:{pdb_code}:{start_res}:{chain_id}:{end_res}:{chain_id}:::0.00: 0.00\n")
            f.write(pdb_aln + "*\n")
            f.write(f">P1;{pdb_code}_seq\n")
            f.write(f"sequence:{pdb_code}_seq:::::::0.00: 0.00\n")
            f.write(fasta_aln + "*\n")

    def complete_missing_atoms(self, uniprot_id=None):
        if uniprot_id is None:
            raise ValueError("UniProt ID is required to complete missing atoms via alignment.")

        logger.info("Completing missing atoms via sequence alignment...")

        pdb_path = Path(self.pdb_path)
        pdb_stem = pdb_path.stem
        ali_file = self.output_dir / f"{pdb_stem}.ali"

        start_res, chain_id = self._get_first_residue_and_chain()
        logger.info(f"First residue: {start_res}, chain: {chain_id}")

        pdb_seq, res_ids = self._extract_pdb_sequence(chain_id)
        full_seq = self._fetch_uniprot_sequence(uniprot_id)
        pdb_aln, fasta_aln = self._align_sequences(pdb_seq, full_seq)

        self._write_pir_alignment(ali_file, pdb_stem, pdb_aln, fasta_aln, res_ids, chain_id)
        
        logger.info(f"Wrote PIR alignment to: {ali_file}")
        with working_directory(self.output_dir):
            automodel = modeller.automodel.AutoModel(self.env,
                                                    alnfile=ali_file.name,
                                                    knowns=pdb_stem,
                                                    sequence=f"{pdb_stem}_seq")
            automodel.starting_model = 1
            automodel.ending_model = 1

            automodel.make()

        generated_model = next(self.output_dir.glob(f"{pdb_stem}_seq.B9999*.pdb"))
        logger.info(f"Loading rebuilt model from: {generated_model}")
        self.model = modeller.Model(self.env, file=generated_model.name)
        self.pdb_path = str(generated_model)

    def _get_residue_range_and_chain(self):
        """Parse PDB to find the actual start/end residue numbers and chain."""
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure("model", self.pdb_path)
        model = structure[0]
        for chain in model:
            residues = [res for res in chain if res.id[0] == ' ']
            if residues:
                start = residues[0].id[1]
                end = residues[-1].id[1]
                return start, end, chain.id
        raise ValueError("No valid residues found in PDB.")

    def fix_loops(self, loop_ranges):
        logger.info(f"Starting loop modeling for {len(loop_ranges)} loop region(s)...")

        pdb_path = Path(self.pdb_path)
        pdb_stem = pdb_path.stem
        ali_file = self.output_dir / f"{pdb_stem}_loop.ali"

        start_res, end_res, chain_id = self._get_residue_range_and_chain()
        logger.info(f"Loop modelling on chain {chain_id}, residues {start_res} to {end_res}")

        pdb_seq, _ = self._extract_pdb_sequence(chain_id)

        # Write PIR alignment
        with open(ali_file, 'w') as f:
            f.write(f">P1;{pdb_stem}\n")
            f.write(f"structure:{pdb_stem}:{start_res}:{chain_id}:{end_res}:{chain_id}:::0.00: 0.00\n")
            f.write(pdb_seq + "*\n")
            f.write(f">P1;{pdb_stem}_seq\n")
            f.write(f"sequence:{pdb_stem}_seq:::::::0.00: 0.00\n")
            f.write(pdb_seq + "*\n")

        # Get residue numbers in model
        model_res_nums = [res.num for res in self.model.residues]
        model_res_range = (min(model_res_nums), max(model_res_nums))
        logger.info(f"Model residue range: {model_res_range}")

        # Filter valid loop ranges
        valid_loop_ranges = []
        for start, end in loop_ranges:
            if model_res_range[0] <= start <= model_res_range[1] and model_res_range[0] <= end <= model_res_range[1]:
                valid_loop_ranges.append((start, end))
            else:
                logger.warning(f"Skipping loop range {start}-{end}: out of bounds.")

        if not valid_loop_ranges:
            logger.warning("No valid loop regions found. Skipping loop modeling.")
            return

        class LoopModel(modeller.automodel.loopmodel):
            def select_loop_atoms(inner_self):
                selections = []
                for start, end in valid_loop_ranges:
                    try:
                        sel = inner_self.residue_range(f"{start}:{chain_id}", f"{end}:{chain_id}")
                        selections.append(sel)
                    except Exception as e:
                        logger.warning(f"Could not select residues {start}-{end} on chain {chain_id}: {e}")
                if not selections:
                    logger.warning("No loop regions could be selected. Returning empty Selection.")
                    return modeller.Selection()
                return modeller.Selection(*selections)

        with working_directory(self.output_dir):
            loop_model = LoopModel(self.env,
                                alnfile=ali_file.name,
                                knowns=pdb_stem,
                                sequence=f"{pdb_stem}_seq")
            loop_model.starting_model = 1
            loop_model.ending_model = self.n_loop_models

            try:
                loop_model.make()
            except modeller.ModellerError as e:
                if "No loop regions could be selected for refinement" in str(e):
                    logger.warning("No valid loop regions found. Skipping loop modeling.")
                    return
                else:
                    raise

            output_model = next(Path(".").glob(f"{pdb_stem}_seq.B99990001.pdb"))
            logger.info(f"Loading loop-modeled structure: {output_model}")
            self.model = modeller.Model(self.env, file=output_model.name)
            self.pdb_path = str(self.output_dir / output_model.name)


    def refine_loops(self):
        logger.info("Refining loops...")

        pdb_path = Path(self.pdb_path)
        pdb_stem = pdb_path.stem
        ali_file = self.output_dir / f"{pdb_stem}_refine.ali"

        start_res, end_res, chain_id = self._get_residue_range_and_chain()
        pdb_seq, _ = self._extract_pdb_sequence(chain_id)

        with open(ali_file, 'w') as f:
            f.write(f">P1;{pdb_stem}\n")
            f.write(f"structure:{pdb_stem}:{start_res}:{chain_id}:{end_res}:{chain_id}:::0.00: 0.00\n")
            f.write(pdb_seq + "*\n")
            f.write(f">P1;{pdb_stem}_ref\n")
            f.write(f"sequence:{pdb_stem}_ref:::::::0.00: 0.00\n")
            f.write(pdb_seq + "*\n")

        class RefinementModel(modeller.automodel.AutoModel):
            def select_atoms(inner_self):
                return modeller.Selection(inner_self.residues)

        n_models = self.loop_refinement.get("n_loop_models", 3)

        with working_directory(self.output_dir):
            refinement = RefinementModel(self.env,
                                        alnfile=ali_file.name,
                                        knowns=pdb_stem,
                                        sequence=f"{pdb_stem}_ref")
            refinement.starting_model = 1
            refinement.ending_model = n_models

            try:
                refinement.make()
            except modeller.ModellerError as e:
                logger.error(f"Loop refinement failed: {e}")
                raise

        refined_model = next(self.output_dir.glob(f"{pdb_stem}_ref.B9999*.pdb"))
        logger.info(f"Loading refined structure: {refined_model}")
        self.model = modeller.Model(self.env, file=refined_model.name)
        self.pdb_path = str(refined_model)
    
    def write_pdb(self, output_path: str):
        logger.info(f"Writing output PDB to {output_path}")
        self.model.write(file=output_path)
