import subprocess
import json

def run_static_analysis(code_dir):
    results = {}

    # Cyclomatic complexity (via Radon)
    try:
        radon_out = subprocess.run(
            ["radon", "cc", "-s", "-j", code_dir],
            capture_output=True, text=True, check=True
        )
        results["complexity"] = json.loads(radon_out.stdout)
    except Exception as e:
        results["complexity_error"] = str(e)

    # Optional: run Pylint to gather potential inefficiency warnings
    try:
        pylint_out = subprocess.run(
            ["pylint", "--disable=all", "--enable=W0101,W0102,W0612", code_dir],
            capture_output=True, text=True
        )
        results["pylint_raw"] = pylint_out.stdout
    except Exception as e:
        results["pylint_error"] = str(e)

    return results
