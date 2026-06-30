from pathlib import Path
from collections import defaultdict
import shutil

# import MDAnalysis as mda
# from MDAnalysis.analysis import align

from Bio.PDB import MMCIFParser
from Bio.PDB import PDBParser
from Bio.PDB import PDBIO
from Bio.PDB import Superimposer

parser = MMCIFParser(QUIET=True)
pdb_parser = PDBParser(QUIET=True)

from pipeline.task_registry import register_task
from pipeline.logger import setup_logger

logger = setup_logger(__name__, debug_mode=False, simple_format=True)

def align_structure(reference_file, mobile_file, output_file):

    ref_structure = pdb_parser.get_structure(
        "ref",
        str(reference_file)
    )

    mobile_structure = pdb_parser.get_structure(
        "mobile",
        str(mobile_file)
    )

    ref_atoms = []
    mobile_atoms = []

    for ref_chain, mob_chain in zip(
        ref_structure[0],
        mobile_structure[0]
    ):
        for ref_res, mob_res in zip(ref_chain, mob_chain):

            if "CA" in ref_res and "CA" in mob_res:

                ref_atoms.append(ref_res["CA"])
                mobile_atoms.append(mob_res["CA"])

    
    logger.debug(
        f"[ALIGN] Matched {len(ref_atoms)} CA atoms"
    )

    if len(ref_atoms) < 3:
        raise ValueError(
            f"Too few atoms for alignment: {mobile_file}"
        )

    sup = Superimposer()

    sup.set_atoms(
        ref_atoms,
        mobile_atoms
    )

    sup.apply(
        mobile_structure.get_atoms()
    )

    io = PDBIO()

    io.set_structure(
        mobile_structure
    )

    io.save(
        str(output_file)
    )

    return sup.rms

def cif_to_pdb(cif_file: Path, pdb_file: Path):

    structure = parser.get_structure(
        pdb_file.stem,
        str(cif_file)
    )

    io = PDBIO()
    io.set_structure(structure)
    io.save(str(pdb_file))

    return pdb_file

def get_input_id(filename: str):

    # sstr5_model1_plddt0.88.cif
    # sstr5_model2_plddt0.83.cif

    return filename.split("_model")[0]

@register_task(
    "generate_alignment_sessions",
    category="Structure modelling",
    description="Align top structures and create PyMOL sessions."
)
def generate_alignment_sessions(backend, config, **kwargs):

    output_dir = Path(config["output_dir"])

    per_input_dir = output_dir / "top_structures_per_input"

    if not per_input_dir.exists():
        raise FileNotFoundError(
            f"{per_input_dir} does not exist"
        )

    session_root = output_dir / "structure_alignments"
    session_root.mkdir(exist_ok=True)

    grouped = defaultdict(list)

    for f in per_input_dir.iterdir():

        if f.suffix.lower() not in [".cif", ".pdb"]:
            continue

        grouped[get_input_id(f.name)].append(f)

    results = {}

    for input_id, files in grouped.items():

        logger.info(
            f"[ALIGN] {input_id}: aligning {len(files)} structures"
        )
        
        files = sorted(
            files,
            key=lambda x: int(
                x.name.split("_model")[1].split("_")[0]
            )
        )

        target_dir = session_root / input_id
        target_dir.mkdir(exist_ok=True)

        aligned_dir = target_dir / "aligned"
        aligned_dir.mkdir(exist_ok=True)

        pdb_files = []

        for f in files:

            if f.suffix.lower() == ".cif":

                pdb_file = (
                    aligned_dir /
                    f"{f.stem}.pdb"
                )

                cif_to_pdb(
                    f,
                    pdb_file
                )

                pdb_files.append(pdb_file)

            else:

                pdb_files.append(f)

        reference = pdb_files[0]

        aligned_files = []

        reference = pdb_files[0]

        reference_pdb = (
            aligned_dir /
            f"{reference.stem}_aligned.pdb"
        )

        shutil.copy(
            reference,
            reference_pdb
        )

        reference = reference_pdb

        aligned_files = [reference]

        for mobile_file in pdb_files[1:]:

            out_file = (
                aligned_dir /
                f"{mobile_file.stem}_aligned.pdb"
            )

            rms = align_structure(
                reference,
                mobile_file,
                out_file
            )

            logger.debug(
                f"[ALIGN] {mobile_file.name} RMSD={rms:.3f}"
            )

            aligned_files.append(out_file)

        # ------------------------
        # create PyMOL script
        # ------------------------

        pml_file = target_dir / f"{input_id}.pml"

        with open(pml_file, "w") as pml:

            pml.write("reinitialize\n")
            pml.write("bg_color white\n\n")

            for structure_file in aligned_files:

                obj_name = structure_file.stem

                relative_path = structure_file.relative_to(target_dir)

                pml.write(
                    f'load "{relative_path.as_posix()}", {obj_name}\n'
                )

            # visualisation defaults

            # orient on reference model
            ref_name = aligned_files[0].stem
            
            pml.write("hide everything\n")
            pml.write("show cartoon\n")

            # colour by pLDDT stored in B-factors
            pml.write("spectrum b, red_yellow_green_cyan_blue, minimum=50, maximum=100\n")

            # optional thicker cartoon in confident regions
            pml.write("set cartoon_putty, on\n")

            pml.write("orient " + ref_name + "\n")
            pml.write("zoom\n")


        logger.info(
            f"[ALIGN] Saved PyMOL script: {pml_file}"
        )

        results[input_id] = {
            "pml": str(pml_file),
            "aligned_structures": [
                str(x)
                for x in aligned_files
            ]
        }

    backend.cache["alignment_sessions"] = results

    return results