import pyrosetta
from pyrosetta.rosetta.core.scoring import disulfides

# List all attributes in the disulfides module
print("disulfides module contents:")
for name in dir(disulfides):
    print("  ", name)

# Try to find any function that sounds like "clear", "remove", "reset"
candidates = [nm for nm in dir(disulfides) if "clear" in nm.lower() or "remove" in nm.lower()]
print("Candidate disulfide-clear methods:", candidates)

# Also inspect pose.conformation capabilities to change cysteine disulfide state
from pyrosetta.rosetta.core.conformation import Conformation
print("Conformation methods containing 'cys' or 'disulfide':", [nm for nm in dir(Conformation) if "cys" in nm.lower()])

