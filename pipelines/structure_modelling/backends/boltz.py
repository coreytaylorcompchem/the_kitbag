import os
import yaml
import subprocess
import logging
import json
import shutil

from collections import OrderedDict
from pathlib import Path

from backends.base import BaseStructureTool

from pipeline.logger import setup_logger

logger = setup_logger(__name__, debug_mode=False, simple_format=True)

def resolve_affinity_config(
    affinity_cfg,
    ligands,
):
    enabled = affinity_cfg.get("enabled", False)

    if not enabled:
        return None

    mode = affinity_cfg.get("mode", "auto")
    binder = affinity_cfg.get("binder", "auto")

    if mode == "never":
        return None

    if mode not in {"auto", "always"}:
        raise ValueError(
            f"Unsupported affinity.mode: {mode}. Use auto, always, or never."
        )

    has_ligand = len(ligands or []) > 0
    ligand_chain_ids = [
        f"L{i+1}"
        for i, _ in enumerate(ligands or [])
    ]

    if not has_ligand:
        if mode == "always":
            raise ValueError(
                "Affinity requested with mode=always, but no ligands are present"
            )

        logger.info("[Affinity] No ligands present; skipping affinity for this input")
        return None

    if binder == "auto":
        binder = ligand_chain_ids[0]

    if binder not in ligand_chain_ids:
        msg = (
            f"Affinity binder {binder} is not present. "
            f"Available ligand chains: {ligand_chain_ids}"
        )

        if mode == "always":
            raise ValueError(msg)

        logger.warning(f"[Affinity] {msg}; skipping affinity for this input")
        return None

    return {
        "enabled": True,
        "binder": binder,
    }

def validate_msa_policy(
    sequences: dict,
    msas: dict,
    use_msa_server: bool,
    tool_name: str,
):
    missing = [
        chain_id
        for chain_id in sequences
        if not msas or chain_id not in msas
    ]

    if missing and not use_msa_server:
        raise ValueError(
            f"{tool_name}: use_msa_server=false, but no local MSA was provided "
            f"for chains: {', '.join(missing)}"
        )
    
def normalize_boltz_token(token):
    """
    Accept either:
      [chain, residue]
      [chain, residue, atom]
      [chain, atom_name] for ligand atom-style contacts

    Return as a plain list for YAML serialization.
    """
    if not isinstance(token, (list, tuple)):
        raise ValueError(
            f"Constraint token must be a list/tuple, got: {token}"
        )

    if len(token) not in (2, 3):
        raise ValueError(
            f"Constraint token must have length 2 or 3, got: {token}"
        )

    return list(token)

def get_available_chain_ids(sequences: dict, ligands: list):
    chain_ids = set(sequences.keys())

    for i, _ in enumerate(ligands or []):
        chain_ids.add(f"L{i+1}")

    return chain_ids


def validate_protein_residue(chain_id, residue, sequences, context, mode):
    if chain_id not in sequences:
        return True

    try:
        residue = int(residue)
    except Exception:
        return True

    seq_len = len(sequences[chain_id])

    if residue < 1 or residue > seq_len:
        msg = (
            f"{context}: residue {residue} is outside chain {chain_id} "
            f"length 1-{seq_len}"
        )

        if mode == "always":
            raise ValueError(msg)

        logger.warning(f"[Constraints] Skipping constraint: {msg}")
        return False

    return True


def tokens_are_available(tokens, available_chains, sequences, mode, context):
    for token in tokens:
        if not isinstance(token, (list, tuple)) or len(token) < 2:
            continue

        chain_id = token[0]

        if chain_id not in available_chains:
            msg = (
                f"{context}: chain {chain_id} not present in this input. "
                f"Available chains: {sorted(available_chains)}"
            )

            if mode == "always":
                raise ValueError(msg)

            logger.info(f"[Constraints] Skipping constraint: {msg}")
            return False

        residue_or_atom = token[1]

        if isinstance(residue_or_atom, int):
            if not validate_protein_residue(
                chain_id=chain_id,
                residue=residue_or_atom,
                sequences=sequences,
                context=context,
                mode=mode,
            ):
                return False

    return True


def build_boltz_constraints(
    constraints_cfg: dict,
    sequences: dict,    
    ligands: list,
    row_constraints: list = None
):
    """
    Compile pipeline constraint config into native Boltz YAML constraints.

    mode:
      auto   -> skip constraints whose chains/residues are absent
      always -> raise error if a requested constraint cannot be applied
      never  -> disable constraints
    """

    if not constraints_cfg or not constraints_cfg.get("enabled", False):
        return []

    mode = constraints_cfg.get("mode", "auto")

    if mode == "never":
        return []

    if mode not in {"auto", "always"}:
        raise ValueError(
            f"Unsupported constraints.mode: {mode}. Use auto, always, or never."
        )

    available_chains = get_available_chain_ids(
        sequences=sequences,
        ligands=ligands,
    )

    compiled = []
    row_constraints = row_constraints or []

    # ------------------------------------------------------------
    # Native Boltz constraints
    # ------------------------------------------------------------
    for item in constraints_cfg.get("native", []):

        if "contact" in item:
            c = item["contact"]
            tokens = [c.get("token1"), c.get("token2")]

            if tokens_are_available(
                tokens=tokens,
                available_chains=available_chains,
                sequences=sequences,
                mode=mode,
                context="native contact constraint",
            ):
                compiled.append(item)

        elif "bond" in item:
            c = item["bond"]
            tokens = [c.get("atom1"), c.get("atom2")]

            if tokens_are_available(
                tokens=tokens,
                available_chains=available_chains,
                sequences=sequences,
                mode=mode,
                context="native bond constraint",
            ):
                compiled.append(item)

        elif "pocket" in item:
            c = item["pocket"]
            binder = c.get("binder")
            contacts = c.get("contacts", [])

            tokens = [[binder, 1]] if binder else []
            tokens.extend(contacts)

            if tokens_are_available(
                tokens=tokens,
                available_chains=available_chains,
                sequences=sequences,
                mode=mode,
                context="native pocket constraint",
            ):
                compiled.append(item)

        else:
            if mode == "always":
                raise ValueError(f"Unsupported native constraint: {item}")

            logger.warning(
                f"[Constraints] Skipping unsupported native constraint: {item}"
            )

    # ------------------------------------------------------------
    # CSV-derived constraints
    # ------------------------------------------------------------
    csv_cfg = constraints_cfg.get("csv", {})

    if csv_cfg.get("enabled", False):

        csv_type = csv_cfg.get("type", "chain_contact")
        default_max_distance = csv_cfg.get("max_distance")
        default_force = csv_cfg.get("force", False)

        for item in row_constraints:

            constraint_type = item.get("type", csv_type)

            if constraint_type != "chain_contact":
                msg = f"Unsupported CSV constraint type: {constraint_type}"

                if mode == "always":
                    raise ValueError(msg)

                logger.warning(f"[Constraints] Skipping CSV constraint: {msg}")
                continue

            chain1 = item["chain1"]
            chain2 = item["chain2"]

            if chain1 not in available_chains or chain2 not in available_chains:
                msg = (
                    f"CSV chain_contact requires chains {chain1}/{chain2}, "
                    f"but available chains are {sorted(available_chains)}"
                )

                if mode == "always":
                    raise ValueError(msg)

                logger.info(f"[Constraints] Skipping CSV constraint: {msg}")
                continue

            residues1 = item["residues1"]
            residues2 = item["residues2"]

            max_distance = item.get(
                "max_distance",
                default_max_distance,
            )

            if max_distance is None:
                raise ValueError(
                    "CSV constraint requires max_distance either in CSV constraint "
                    "or structure_prediction.constraints.csv.max_distance"
                )

            force = item.get(
                "force",
                default_force,
            )

            for res1 in residues1:
                if not validate_protein_residue(
                    chain_id=chain1,
                    residue=res1,
                    sequences=sequences,
                    context="CSV chain_contact",
                    mode=mode,
                ):
                    continue

                for res2 in residues2:
                    if not validate_protein_residue(
                        chain_id=chain2,
                        residue=res2,
                        sequences=sequences,
                        context="CSV chain_contact",
                        mode=mode,
                    ):
                        continue

                    compiled.append({
                        "contact": {
                            "token1": [chain1, int(res1)],
                            "token2": [chain2, int(res2)],
                            "max_distance": float(max_distance),
                            "force": bool(force),
                        }
                    })

    return compiled

def pdb_has_seqres(pdb_path: Path) -> bool:
    """
    Return True if the PDB contains at least one SEQRES record.
    """
    pdb_path = Path(pdb_path)

    with pdb_path.open("r", errors="replace") as handle:
        return any(
            line.startswith("SEQRES")
            for line in handle
        )


def extract_atom_residues_by_chain(pdb_path: Path) -> OrderedDict:
    """
    Extract ordered polymer residue names from ATOM records.

    Residues are identified by:
      - chain ID
      - residue number
      - insertion code

    Only ATOM records are used. HETATM records such as waters, ions,
    ligands, glycans, and cofactors are deliberately excluded.

    Alternative locations other than blank or A are ignored.
    """
    pdb_path = Path(pdb_path)

    residues_by_chain = OrderedDict()
    seen_residues = set()

    with pdb_path.open("r", errors="replace") as handle:
        for line in handle:
            if not line.startswith("ATOM  "):
                continue

            if len(line) < 27:
                logger.warning(
                    f"[Template repair] Skipping malformed ATOM line "
                    f"in {pdb_path}: {line.rstrip()}"
                )
                continue

            # PDB fixed-width fields.
            altloc = line[16].strip()
            residue_name = line[17:20].strip().upper()
            chain_id = line[21].strip()
            residue_number = line[22:26].strip()
            insertion_code = line[26].strip()

            # Ignore secondary alternative conformations.
            if altloc not in {"", "A"}:
                continue

            if not residue_name:
                continue

            # Blank chain IDs are valid in PDB files, although not ideal.
            residue_key = (
                chain_id,
                residue_number,
                insertion_code,
            )

            if residue_key in seen_residues:
                continue

            seen_residues.add(residue_key)

            residues_by_chain.setdefault(chain_id, [])
            residues_by_chain[chain_id].append(residue_name)

    return residues_by_chain


def format_seqres_records(residues_by_chain: dict) -> list:
    """
    Construct PDB SEQRES records from ordered three-letter residue names.

    PDB SEQRES records contain up to 13 residues per line.
    """
    records = []

    for chain_id, residue_names in residues_by_chain.items():
        if not residue_names:
            continue

        displayed_chain_id = chain_id if chain_id else " "
        total_residues = len(residue_names)

        chunks = [
            residue_names[i:i + 13]
            for i in range(0, total_residues, 13)
        ]

        for serial_number, chunk in enumerate(chunks, start=1):
            residue_text = " ".join(
                f"{name:>3}"
                for name in chunk
            )

            records.append(
                f"SEQRES {serial_number:>3} "
                f"{displayed_chain_id:1} "
                f"{total_residues:>4}  "
                f"{residue_text}\n"
            )

    return records


def add_seqres_to_pdb(
    input_path: Path,
    output_path: Path,
) -> Path:
    """
    Create a copy of a PDB file with SEQRES records reconstructed from
    coordinate-bearing ATOM residues.

    The original PDB is not modified.

    Important:
        The reconstructed SEQRES contains only residues present in the
        coordinate section. Unresolved residues absent from the PDB cannot
        be recovered by this procedure.
    """
    input_path = Path(input_path).expanduser().resolve()
    output_path = Path(output_path).expanduser().resolve()

    if not input_path.is_file():
        raise FileNotFoundError(
            f"Template PDB does not exist: {input_path}"
        )

    residues_by_chain = extract_atom_residues_by_chain(input_path)

    if not residues_by_chain:
        raise ValueError(
            f"Cannot reconstruct SEQRES for {input_path}: "
            "no ATOM polymer residues were found."
        )

    seqres_records = format_seqres_records(residues_by_chain)

    if not seqres_records:
        raise ValueError(
            f"Cannot reconstruct SEQRES for {input_path}: "
            "no SEQRES records could be generated."
        )

    original_lines = input_path.read_text(
        errors="replace"
    ).splitlines(keepends=True)

    # Remove any malformed or partial existing SEQRES records.
    original_lines = [
        line
        for line in original_lines
        if not line.startswith("SEQRES")
    ]

    # Place SEQRES before the first coordinate or MODEL record.
    coordinate_start = next(
        (
            index
            for index, line in enumerate(original_lines)
            if line.startswith(("ATOM  ", "HETATM", "MODEL "))
        ),
        len(original_lines),
    )

    repaired_lines = (
        original_lines[:coordinate_start]
        + seqres_records
        + original_lines[coordinate_start:]
    )

    output_path.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    output_path.write_text(
        "".join(repaired_lines)
    )

    chain_summary = {
        chain_id if chain_id else "<blank>": len(residue_names)
        for chain_id, residue_names in residues_by_chain.items()
    }

    logger.info(
        f"[Template repair] Added SEQRES to {input_path.name} "
        f"| chains={chain_summary} "
        f"| output={output_path}"
    )

    return output_path


def validate_boltz_template(template_path: Path) -> None:
    """
    Validate that Gemmi can read the template and that it contains at
    least one coordinate model, chain, and polymer residue.

    Full Boltz parsing is deliberately not performed here because the
    low-level parse_pdb() and parse_mmcif() functions require Boltz CCD
    molecule resources that are initialized by the Boltz CLI.
    """
    import gemmi

    template_path = Path(template_path).expanduser().resolve()

    if not template_path.is_file():
        raise FileNotFoundError(
            f"Boltz template does not exist: {template_path}"
        )

    if template_path.stat().st_size == 0:
        raise ValueError(
            f"Boltz template is empty: {template_path}"
        )

    suffix = template_path.suffix.lower()

    if suffix not in {".pdb", ".cif", ".mmcif"}:
        raise ValueError(
            f"Unsupported Boltz template extension "
            f"{suffix!r}: {template_path}"
        )

    try:
        structure = gemmi.read_structure(str(template_path))
    except Exception as exc:
        raise ValueError(
            f"Gemmi cannot read template: {template_path}\n"
            f"Original Gemmi error: {type(exc).__name__}: {exc}"
        ) from exc

    if len(structure) == 0:
        raise ValueError(
            f"Template contains no coordinate models recognized by Gemmi: "
            f"{template_path}"
        )

    first_model = structure[0]

    if len(first_model) == 0:
        raise ValueError(
            f"Template model 1 contains no chains: {template_path}"
        )

    chain_summary = {}
    total_polymer_residues = 0

    for chain in first_model:
        polymer_residues = [
            residue
            for residue in chain
            if residue.entity_type == gemmi.EntityType.Polymer
        ]

        chain_name = chain.name if chain.name else "<blank>"
        residue_count = len(polymer_residues)

        chain_summary[chain_name] = residue_count
        total_polymer_residues += residue_count

    if total_polymer_residues == 0:
        raise ValueError(
            f"Template contains no polymer residues recognized by Gemmi: "
            f"{template_path}\n"
            f"Chains found: {chain_summary}"
        )

    logger.info(
        f"[Template validation] Gemmi successfully parsed template "
        f"| models={len(structure)} "
        f"| chains={chain_summary} "
        f"| template={template_path}"
    )


def prepare_boltz_templates(
    templates: list,
    run_dir: Path,
    auto_add_seqres: bool = True,
) -> list:
    """
    Normalize, optionally repair, and validate Boltz template paths.

    PyMOL-exported PDB files commonly lack SEQRES records. If
    auto_add_seqres is True, a repaired copy is created under the
    current Boltz run directory.

    The original template file is never modified.
    """
    prepared_templates = []

    if not templates:
        return prepared_templates

    template_dir = Path(run_dir) / "prepared_templates"
    template_dir.mkdir(
        parents=True,
        exist_ok=True,
    )

    for index, template in enumerate(templates, start=1):
        template_path = Path(template).expanduser().resolve()

        if not template_path.is_file():
            raise FileNotFoundError(
                f"Template does not exist: {template_path}"
            )

        suffix = template_path.suffix.lower()

        if suffix not in {".pdb", ".cif", ".mmcif"}:
            raise ValueError(
                f"Unsupported template format: {template_path}"
            )

        if suffix == ".pdb" and not pdb_has_seqres(template_path):
            if not auto_add_seqres:
                raise ValueError(
                    f"Template PDB contains no SEQRES records: "
                    f"{template_path}\n"
                    "Either provide a PDB or mmCIF with complete sequence "
                    "metadata, or enable "
                    "structure_prediction.templates.auto_add_seqres."
                )

            repaired_name = (
                f"template_{index}_"
                f"{template_path.stem}_with_seqres.pdb"
            )

            prepared_path = template_dir / repaired_name

            add_seqres_to_pdb(
                input_path=template_path,
                output_path=prepared_path,
            )

        else:
            copied_name = (
                f"template_{index}_{template_path.name}"
            )

            prepared_path = template_dir / copied_name

            shutil.copy2(
                template_path,
                prepared_path,
            )

            logger.info(
                f"[Template preparation] Copied template "
                f"{template_path} to {prepared_path}"
            )

        validate_boltz_template(prepared_path)

        prepared_templates.append(
            str(prepared_path.resolve())
        )

    return prepared_templates

def generate_boltz_yaml(
    sequences: dict,
    ligands: list,
    templates: list,
    msas: dict,
    constraints: list,
    affinity: dict,
    yaml_path: Path,
):

    data = {
        "sequences": []
    }

    # templates structures, if present.
    
    if templates:
        data["templates"] = []

        for t in templates:
            template_path = Path(t).expanduser().resolve()
            ext = template_path.suffix.lower()

            if ext in {".cif", ".mmcif"}:
                entry = {
                    "cif": str(template_path)
                }
            elif ext == ".pdb":
                entry = {
                    "pdb": str(template_path)
                }
            else:
                raise ValueError(
                    f"Unsupported template format: {template_path}"
                )

            data["templates"].append(entry)

            data["templates"].append(entry)
               
        logger.info(
            f"Templates: {len(templates) if templates else 0}"
        )

    # proteins 
    for chain_id, seq in sequences.items():
        protein_entry = {
            "id": chain_id,
            "sequence": seq,
        }

        if msas and chain_id in msas:
            protein_entry["msa"] = str(msas[chain_id])

        data["sequences"].append({
            "protein": protein_entry
        })

    # ligands 
    for i, smiles in enumerate(ligands):
        lig_id = f"L{i+1}"

        data["sequences"].append({
            "ligand": {
                "id": lig_id,
                "smiles": smiles
            }
        })
    
    # constraints
    if constraints:
        data["constraints"] = constraints
        logger.info(
            f"Boltz constraints: {len(constraints)}"
        )

    # affinity property
    if affinity and affinity.get("enabled", False):
        binder = affinity.get("binder")

        if not binder:
            raise ValueError(
                "structure_prediction.affinity.enabled=true, but no binder was provided"
            )

        data["properties"] = [
            {
                "affinity": {
                    "binder": binder
                }
            }
        ]

        logger.info(
            f"Boltz affinity enabled for binder chain: {binder}"
        )

    # strict validation
    if not data["sequences"]:
        raise ValueError("No valid sequences or ligands provided for Boltz input")

    yaml_str = yaml.dump(data, sort_keys=False)

    yaml_path.write_text(yaml_str)

def log_parameters(params: dict, title="Parameters"):
    if not params:
        logger.info(f"[{title}] No parameters")
        return

    key_width = max(len(str(k)) for k in params.keys()) + 2

    logger.info(f"[{title}]")
    logger.info("─" * (key_width + 20))

    for k, v in params.items():
        logger.info(f"{k:<{key_width}} {v}")

    logger.info("─" * (key_width + 20))

class BoltzBackend(BaseStructureTool):
    def _parse_affinity(self, boltz_dir: Path):
        affinity_files = sorted(boltz_dir.rglob("*affinity*.json"))

        if not affinity_files:
            logger.info(f"[Affinity] No affinity JSON found under {boltz_dir}")
            return {}

        if len(affinity_files) > 1:
            logger.warning(
                f"[Affinity] Found multiple affinity JSON files under {boltz_dir}; "
                f"using first: {affinity_files[0]}"
            )

        affinity_file = affinity_files[0]

        logger.info(f"[Affinity] Parsing affinity JSON: {affinity_file}")

        try:
            data = json.loads(affinity_file.read_text())
        except Exception as e:
            logger.warning(
                f"Failed to parse affinity JSON {affinity_file}: {e}"
            )
            return {}

        parsed = {
            "affinity_json": str(affinity_file),
            "affinity_pred_value": data.get("affinity_pred_value"),
            "affinity_probability_binary": data.get("affinity_probability_binary"),
            "affinity_pred_value1": data.get("affinity_pred_value1"),
            "affinity_probability_binary1": data.get("affinity_probability_binary1"),
            "affinity_pred_value2": data.get("affinity_pred_value2"),
            "affinity_probability_binary2": data.get("affinity_probability_binary2"),
        }

        logger.info(
            "[Affinity] Parsed metrics: "
            f"affinity_pred_value={parsed['affinity_pred_value']} "
            f"affinity_probability_binary={parsed['affinity_probability_binary']}"
        )

        return parsed
    
    name = "boltz"
    
    def __init__(self, config=None):
            super().__init__()
            self.config = config or {}
    def run(
        self,
        run_id: int,
        device: int,
        output_dir: Path,
        sequences: dict,
        ligands: list = None,
        templates: list = None,
        msas: dict = None,
        row_constraints: list = None,
    ):

        env = os.environ.copy()
        env["CUDA_VISIBLE_DEVICES"] = str(device)

        run_dir = output_dir / f"run_{run_id}"
        run_dir.mkdir(parents=True, exist_ok=True)

        yaml_file = run_dir / "input.yaml"
        out_dir = run_dir / "results"

        # config extraction
        sp_cfg = self.config.get("structure_prediction", {})
        inf_cfg = sp_cfg.get("inference", {})
        constraints_cfg = sp_cfg.get("constraints", {})
        affinity_cfg = sp_cfg.get("affinity", {})
        templates_cfg = sp_cfg.get("templates", {})
      
        affinity_runtime = resolve_affinity_config(
            affinity_cfg,
            ligands or [],
        )

        boltz_constraints = build_boltz_constraints(
            constraints_cfg=constraints_cfg,
            sequences=sequences,
            ligands=ligands or [],
            row_constraints=row_constraints or [],
        )

        prepared_templates = prepare_boltz_templates(
            templates=templates or [],
            run_dir=run_dir,
            auto_add_seqres=templates_cfg.get(
                "auto_add_seqres",
                True,
            ),
        )

        generate_boltz_yaml(
            sequences=sequences,
            ligands=ligands or [],
            templates=prepared_templates,
            msas=msas or {},
            constraints=boltz_constraints,
            affinity=affinity_runtime,
            yaml_path=yaml_file,
        )

        cmd = [
            "boltz",
            "predict",
            str(yaml_file),
            "--out_dir", str(out_dir),
        ]

        use_msa_server = bool(inf_cfg.get("use_msa_server", False))

        validate_msa_policy(
            sequences=sequences,
            msas=msas or {},
            use_msa_server=use_msa_server,
            tool_name=self.name,
        )

        # accelerator 
        accelerator = inf_cfg.get("accelerator", "gpu")

        # Add yaml parameters to Boltz command.
        # accelerator
        cmd += ["--accelerator", accelerator]

        # optional flags
        if inf_cfg.get("use_msa_server", True):
            cmd.append("--use_msa_server")

        if inf_cfg.get("no_kernels", True):
            cmd.append("--no_kernels")            
        if inf_cfg.get("use_potentials", False):
            cmd.append("--use_potentials")

        # numeric params
        if "diffusion_samples" in inf_cfg:
            cmd += ["--diffusion_samples", str(inf_cfg["diffusion_samples"])]

        if "sampling_steps" in inf_cfg:
            cmd += ["--sampling_steps", str(inf_cfg["sampling_steps"])]

        if "step_scale" in inf_cfg:
            cmd += ["--step_scale", str(inf_cfg["step_scale"])]

        if "recycling_steps" in inf_cfg:
            cmd += ["--recycling_steps", str(inf_cfg["recycling_steps"])]

        if "sampling_steps_affinity" in inf_cfg:
            cmd += ["--sampling_steps_affinity", str(inf_cfg["sampling_steps_affinity"])]

        if "diffusion_samples_affinity" in inf_cfg:
            cmd += ["--diffusion_samples_affinity", str(inf_cfg["diffusion_samples_affinity"])]
        
        if affinity_runtime:

            if affinity_cfg.get("mw_correction", False):
                cmd.append("--affinity_mw_correction")

            if "affinity_checkpoint" in affinity_cfg:
                cmd += [
                    "--affinity_checkpoint",
                    str(affinity_cfg["affinity_checkpoint"]),
                ]

        # clean run label
        run_label = f"Boltz run={run_id} device={device}"

        # runtime parameters
        runtime_params = {
            "accelerator": accelerator,
            "use_msa_server": inf_cfg.get("use_msa_server", True),
            "no_kernels": inf_cfg.get("no_kernels", True),
            "use_potentials": inf_cfg.get("use_potentials", True),
        }

        # sampling parameters
        sampling_params = {
            "diffusion_samples": inf_cfg.get("diffusion_samples"),
            "sampling_steps": inf_cfg.get("sampling_steps"),
            "step_scale": inf_cfg.get("step_scale"),
            "recycling_steps": inf_cfg.get("recycling_steps"),
        }

        # affinity parameters 
        affinity_params = {
            "affinity_requested": affinity_cfg.get("enabled", False),
            "affinity_runtime_enabled": bool(affinity_runtime),
            "affinity_mode": affinity_cfg.get("mode", "auto"),
            "affinity_binder": affinity_runtime.get("binder") if affinity_runtime else None,
            "affinity_mw_correction": affinity_cfg.get("mw_correction", False),
            "affinity_checkpoint": affinity_cfg.get("affinity_checkpoint"),
            "diffusion_samples_affinity": inf_cfg.get("diffusion_samples_affinity"),
            "sampling_steps_affinity": inf_cfg.get("sampling_steps_affinity"),
        }

        constraints_params = {
            "constraints_enabled": constraints_cfg.get("enabled", False),
            "constraints_mode": constraints_cfg.get("mode", "auto"),
            "native_constraints": len(constraints_cfg.get("native", [])),
            "smart_constraints": len(constraints_cfg.get("smart", [])),
            "csv_constraints": len(row_constraints or []),
            "compiled_constraints": len(boltz_constraints),
        }

        # remove None values
        sampling_params = {k: v for k, v in sampling_params.items() if v is not None}
        affinity_params = {k: v for k, v in affinity_params.items() if v is not None}

        # log cleanly
        log_parameters(runtime_params, f"{run_label} | Runtime")
        log_parameters(sampling_params, f"{run_label} | Sampling")
        log_parameters(constraints_params, f"{run_label} | Constraints")

        if affinity_params:
            log_parameters(affinity_params, f"{run_label} | Affinity")

        subprocess.run(cmd, check=True, env=env)

        results = self._parse_results(out_dir)
        
        if not results:
            logger.warning(f"[WARNING] No samples returned for run {run_id}")


        return {
            "run_id": run_id,
            "device": device,
            "results": results,
            "tool": self.name,
        }
    
    def _parse_results(self, out_dir: Path):

        # find boltz_results directory
        boltz_dirs = list(out_dir.glob("boltz_results_*"))

        if not boltz_dirs:
            logger.warning(f"No boltz_results directory found in {out_dir}")
            return []

        samples = []
        affinity_data = self._parse_affinity(boltz_dirs[0])

        # recursively find JSON files
        json_files = list(boltz_dirs[0].rglob("confidence*.json"))

        if not json_files:
            logger.warning(f"No confidence JSON files found under {boltz_dirs[0]}")
            return []

        for json_file in json_files:
            logger.debug(f"Found JSON: {json_file}")

            data = json.loads(json_file.read_text())

            # find corresponding structure file

            parent_dir = json_file.parent

            # confidence_input_model_3.json
            # -> input_model_3

            model_stem = json_file.stem.replace(
                "confidence_",
                ""
            )
    
            structure_file = None

            for ext in [".cif", ".pdb"]:

                candidate = parent_dir / f"{model_stem}{ext}"

                if candidate.exists():
                    structure_file = candidate
                    break

            if structure_file is None:
                logger.warning(
                    f"No structure file found for {json_file}"
                )
                continue

            samples.append({
                "structure": str(structure_file),

                # core confidence metrics
                "confidence_score": data.get("confidence_score"),
                "plddt": (
                    data.get("complex_plddt")
                    or data.get("mean_plddt")
                    or data.get("plddt")
                ),
                "iptm": data.get("iptm"),
                "ptm": data.get("ptm"),
                "iplddt": data.get("complex_iplddt"),

                # additional Boltz metrics
                "ligand_iptm": data.get("ligand_iptm"),
                "protein_iptm": data.get("protein_iptm"),
                "complex_pde": data.get("complex_pde"),
                "complex_ipde": data.get("complex_ipde"),
                "chains_ptm": data.get("chains_ptm"),
                "pair_chains_iptm": data.get("pair_chains_iptm"),

                # affinity metrics, same per input rather than per model
                **affinity_data,
            })

        return samples

