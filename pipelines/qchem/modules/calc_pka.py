import tempfile
import os
import pickle
from pipeline.task_registry import register_task
from pipeline.logger import setup_logger

logger = setup_logger(__name__, debug_mode=True, simple_format=True)

PROTON_GAS_PHASE_G = -6.28  # kcal/mol, literature gas phase proton free energy (approximate)

@register_task(
    'calc_pka',
    description="Calculate site pKa values from HA/A- free energies.",
    modifies_geometry=False,
    category='Property'
)
def run(backend, input_for_task, step_config, global_config=None):
    """
    site_pka_data: dict from site_pka task containing HA and A- xyz and mol objects.
    step_config: dict with Psi4 method, basis, temperature, etc.

    Returns dict:
      {
        (molecule_name, site_idx): {
          'pka': float,
          'deltaG_kcal': float,
          'energy_ha': float,
          'energy_a': float
        },
        ...
      }
    """
    site_pka_data_path = step_config.get('site_pka_data_path')
    if site_pka_data_path is None:
        raise ValueError("Missing site_pka_data_path in step_config")

    with open(site_pka_data_path, 'rb') as f:
        site_pka_data = pickle.load(f)

    method = step_config.get('method', 'b3lyp')
    basis = step_config.get('basis', 'def2-svp')
    temperature = step_config.get('temperature', 298.15)  # Kelvin

    results = {}

    for mol_name, sites in site_pka_data.items():
        for site_idx, data in sites.items():
            ha_xyz = data['ha_xyz']
            a_xyz = data['a_xyz']

            ha_filepath = None
            a_filepath = None

            try:
                # Write HA xyz to temp file
                with tempfile.NamedTemporaryFile(mode='w', suffix='.xyz', delete=False) as ha_file:
                    ha_file.write(ha_xyz)
                    ha_filepath = ha_file.name

                # Write A- xyz to temp file
                with tempfile.NamedTemporaryFile(mode='w', suffix='.xyz', delete=False) as a_file:
                    a_file.write(a_xyz)
                    a_filepath = a_file.name

                logger.debug(f"Optimising HA structure for {mol_name} site {site_idx} from file: {ha_filepath}")
                ha_result = backend.optimise(ha_filepath, {'method': method, 'basis': basis})
                energy_ha = ha_result.get('energy')
                if energy_ha is None:
                    raise ValueError(f"No energy returned for HA species in {mol_name} site {site_idx}")

                logger.debug(f"Optimising A- structure for {mol_name} site {site_idx} from file: {a_filepath}")
                a_result = backend.optimise(a_filepath, {'method': method, 'basis': basis})
                energy_a = a_result.get('energy')
                if energy_a is None:
                    raise ValueError(f"No energy returned for A- species in {mol_name} site {site_idx}")

            except Exception as e:
                logger.warning(f"Failed to calculate energies for {mol_name} site {site_idx}: {e}")
                continue

            finally:
                # Cleanup temp files
                if ha_filepath and os.path.exists(ha_filepath):
                    os.remove(ha_filepath)
                if a_filepath and os.path.exists(a_filepath):
                    os.remove(a_filepath)

            # Convert Hartree to kcal/mol
            hartree_to_kcal = 627.509
            deltaG = (energy_a - energy_ha) * hartree_to_kcal + PROTON_GAS_PHASE_G

            # Calculate pKa: pKa = ΔG / (RT ln10)
            R = 1.9872036e-3  # kcal/(mol K)
            ln10 = 2.302585
            pka = deltaG / (R * temperature * ln10)

            results[(mol_name, site_idx)] = {
                'pka': pka,
                'deltaG_kcal': deltaG,
                'energy_ha': energy_ha,
                'energy_a': energy_a
            }

            logger.info(f"{mol_name} site {site_idx} pKa: {pka:.2f} (ΔG = {deltaG:.2f} kcal/mol)")

    return results


