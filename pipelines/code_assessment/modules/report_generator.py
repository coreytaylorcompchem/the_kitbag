import json
from datetime import datetime

# def generate_report(static_results, workflow_results, profile_results, output_path):
def generate_report(static_results, output_path):
    with open(output_path, "w") as f:
        f.write("<html><body>")
        f.write(f"<h1>Code Efficiency Assessment – {datetime.now()}</h1>")

        f.write("<h2>Static Analysis</h2><pre>")
        f.write(json.dumps(static_results, indent=2))
        f.write("</pre>")

        # f.write("<h2>Workflow Analysis</h2><pre>")
        # f.write(json.dumps(workflow_results["metrics"], indent=2))
        # f.write("</pre>")

        # f.write("<h2>Runtime Profiling</h2><pre>")
        # f.write(profile_results["profile_summary"])
        # f.write("</pre>")

        f.write("</body></html>")
