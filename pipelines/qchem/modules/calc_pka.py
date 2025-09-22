import os
from pipeline.task_registry import register_task
from pipeline.logger import setup_logger

logger = setup_logger(__name__, debug_mode=False, simple_format=True)

PROTON_GAS_PHASE_G = -6.28  # kcal/mol, literature gas phase proton free energy (approximate)

@register_task(
    'calc_pka',
    description="Calculate site pKa values from HA/A- free energies.",
    modifies_geometry=False,
    category='Property'
)
def run(backend, site_pka_data, step_config, global_config=None):
    """
    site_pka_data: dict from site_pka task containing HA and A- xyz and mol objects.
    step_config: dict with Psi4 method, basis, etc.

    Returns dict:
      {
        site_idx: {
          'pka': float,
          'deltaG_kcal': float,
          'energy_ha': float,
          'energy_a': float
        },
        ...
      }
    """
    method = step_config.get('method', 'b3lyp')
    basis = step_config.get('basis', 'def2-svp')
    temperature = step_config.get('temperature', 298.15)  # Kelvin

    results = {}

    for site_idx, data in site_pka_data.items():
        ha_xyz = data['ha_xyz']
        a_xyz = data['a_xyz']

        # Run backend optimize or single_point on HA
        ha_result = backend.optimise(ha_xyz, {'method': method, 'basis': basis})
        energy_ha = ha_result.get('energy')

        # Run backend optimize or single_point on A-
        a_result = backend.optimise(a_xyz, {'method': method, 'basis': basis})
        energy_a = a_result.get('energy')

        # Convert Hartree to kcal/mol
        hartree_to_kcal = 627.509
        deltaG = (energy_a - energy_ha) * hartree_to_kcal + PROTON_GAS_PHASE_G

        # Calculate pKa: pKa = deltaG / (RT ln10)
        R = 1.9872036e-3  # kcal/(mol K)
        ln10 = 2.302585
        pka = deltaG / (R * temperature * ln10)

        results[site_idx] = {
            'pka': pka,
            'deltaG_kcal': deltaG,
            'energy_ha': energy_ha,
            'energy_a': energy_a
        }

        logger.info(f"Site {site_idx} pKa: {pka:.2f} (ΔG = {deltaG:.2f} kcal/mol)")

    return results
