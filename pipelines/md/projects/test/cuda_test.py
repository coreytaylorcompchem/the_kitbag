from openmm import Context, System, VerletIntegrator, Platform

platform = Platform.getPlatformByName('CUDA')

system = System()
system.addParticle(1.0)  # Add one particle with mass 1.0 Dalton

integrator = VerletIntegrator(1.0)

context = Context(system, integrator, platform)

cuda_version = platform.getPropertyValue(context, 'CudaToolkitVersion')
print("CUDA Toolkit version:", cuda_version)

