import csv
import os
import pickle
from pipeline.task_registry import register_task
from pipeline.logger import setup_logger

logger = setup_logger(__name__, debug_mode=False, simple_format=True)

PROTON_GAS_PHASE_G = -6.28  # kcal/mol
PROTON_SOLUTION_G = -265.9  # kcal/mol

@register_task(
    'calc_pka',
    description="Calculate site pKa values from HA/A- free energies.",
    modifies_geometry=False,
    category='Property'
)
def run(backend, input_for_task, step_config, global_config=None):
    site_pka_data_path = step_config.get('site_pka_data_path')
    if site_pka_data_path is None:
        raise ValueError("Missing site_pka_data_path in step_config")

    with open(site_pka_data_path, 'rb') as f:
        site_pka_data = pickle.load(f)

    temperature = step_config.get('temperature', 298.15)

    output_dir = step_config.get("output_dir", "outputs/calc_pka")
    os.makedirs(output_dir, exist_ok=True)

    csv_path = os.path.join(output_dir, "site_pka_results.csv")
    write_header = not os.path.exists(csv_path)

    results = {}

    with open(csv_path, "a", newline="") as csvfile:
        writer = csv.writer(csvfile)
        if write_header:
            writer.writerow([
                "Molecule", "Site", "pKa", "ΔG_kcal", "Energy_HA_Hartree", "Energy_A_Hartree"
            ])

        for mol_name, sites in site_pka_data.items():
            for site_idx, data in sites.items():
                ha_xyz = data['ha_xyz']
                a_xyz = data['a_xyz']

                ha_filepath = None
                a_filepath = None

                try:
                    # Temp files for geometries
                    ha_filepath = os.path.join(output_dir, f"{mol_name}_site_{site_idx}_HA.xyz")
                    with open(ha_filepath, "w") as f:
                        f.write(ha_xyz)

                    a_filepath = os.path.join(output_dir, f"{mol_name}_site_{site_idx}_A.xyz")
                    with open(a_filepath, "w") as f:
                        f.write(a_xyz)

                    # Create log paths
                    ha_log_path = os.path.join(output_dir, f"{mol_name}_site_{site_idx}_HA.log")
                    a_log_path = os.path.join(output_dir, f"{mol_name}_site_{site_idx}_A.log")

                    # HA
                    ha_step = {**step_config}
                    ha_step.update({"charge": data.get("ha_charge", 0),
                                    "multiplicity": data.get("ha_mult", 1),
                                    "freq": True})  # request frequencies

                    if "solvent" in step_config:
                        ha_step["solvent"] = step_config["solvent"]

                    logger.info(f"Optimising HA structure for {mol_name} site {site_idx} "
                                f"from file: {ha_filepath} (charge={ha_step['charge']} mult={ha_step['multiplicity']})")

                    ha_result = backend.optimise(ha_filepath, ha_step, log_path=ha_log_path)

                    energy_ha = ha_result.get('gibbs_energy', ha_result['energy'])

                    if energy_ha is None:
                        raise ValueError(f"No energy returned for HA species in {mol_name} site {site_idx}")

                    # A–
                    a_step = {**step_config}
                    a_step.update({"charge": data.get("a_charge", -1),
                                "multiplicity": data.get("a_mult", 1),
                                "freq": True})

                    if "solvent" in step_config:
                        a_step["solvent"] = step_config["solvent"]

                    logger.info(f"Optimising A- structure for {mol_name} site {site_idx} "
                                f"from file: {a_filepath} (charge={a_step['charge']} mult={a_step['multiplicity']})")

                    a_result = backend.optimise(a_filepath, a_step, log_path=a_log_path)

                    energy_a = a_result.get('gibbs_energy', a_result['energy'])

                    if energy_a is None:
                        raise ValueError(f"No energy returned for A- species in {mol_name} site {site_idx}")

                except Exception as e:
                    logger.warning(f"Failed to calculate energies for {mol_name} site {site_idx}: {e}")
                    continue

                finally:
                    if ha_filepath and os.path.exists(ha_filepath):
                        os.remove(ha_filepath)
                    if a_filepath and os.path.exists(a_filepath):
                        os.remove(a_filepath)

                # ΔG
                hartree_to_kcal = 627.509
                deltaG = (energy_a - energy_ha) * hartree_to_kcal + PROTON_SOLUTION_G

                R = 1.9872036e-3
                ln10 = 2.302585
                pka = deltaG / (R * temperature * ln10)

                results[(mol_name, site_idx)] = {
                    'pka': pka,
                    'deltaG_kcal': deltaG,
                    'energy_ha': energy_ha,
                    'energy_a': energy_a
                }

                # Write to CSV
                writer.writerow([
                    mol_name,
                    site_idx,
                    f"{pka:.2f}",
                    f"{deltaG:.2f}",
                    f"{energy_ha:.6f}",
                    f"{energy_a:.6f}"
                ])
                csvfile.flush()

                logger.info(f"{mol_name} site {site_idx} pKa: {pka:.2f} (ΔG = {deltaG:.2f} kcal/mol)")

    return results
