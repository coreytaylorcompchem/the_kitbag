from colabfold.batch import run
from pathlib import Path
import os
os.environ["JAX_PLATFORM_NAME"] = "cpu"
# Define a test query
queries = [("test", "MSPQTETKASVGFKAGVKEYKLTYYTPEYETKDTDILAAFRVTPQPGVPPEEAGAAVAAESSTGTWTTVWTDGLTSLDRYKGRCYRIERVVGE", None)]

results = run(
    queries=queries,
    result_dir="colabfold_outputs",
    use_templates=False,
    num_models=1,  # ✅ REQUIRED
    is_complex=False,
    data_dir=Path.home() / ".colabfold" / "data",  # ✅ Safe default
)

