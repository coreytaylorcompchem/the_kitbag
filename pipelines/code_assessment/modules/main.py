import os
from datetime import datetime
import json

from modules.static_analysis import run_static_analysis
from modules.workflow_analysis import analyze_workflows
from modules.runtime_profiling import profile_pipeline
from modules.report_generator import generate_report

TOP_N = 5

def summarize_complexity(static_results, top_n=TOP_N):
    """Extract and print the top N most complex functions"""
    complexity_list = []

    for filepath, funcs in static_results.get("complexity", {}).items():
        for fn in funcs:
            complexity_list.append({
                "file": filepath,
                "function": fn["name"],
                "complexity": fn["complexity"],
                "rank": fn["rank"],
                "start_line": fn["lineno"],
                "end_line": fn["endline"]
            })

    # Sort descending by complexity
    complexity_list.sort(key=lambda x: x["complexity"], reverse=True)

    print(f"\nTop {top_n} most complex functions:\n")
    for fn in complexity_list[:top_n]:
        print(f"{fn['file']}:{fn['start_line']}-{fn['end_line']} | "
              f"{fn['function']}() | Complexity: {fn['complexity']} | Rank: {fn['rank']}")

def run_assessment(code_dir, workflows_dir, profile_entrypoint, output_dir):
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    os.makedirs(output_dir, exist_ok=True)

    print(f"[INFO] Starting assessment at {timestamp}")

    # --- Stage 1: Static Analysis ---
    print("[INFO] Running static analysis...")
    static_results = run_static_analysis(code_dir)

    # # --- Stage 2: Workflow Graph Analysis ---
    # print("[INFO] Analyzing YAML workflows...")
    # workflow_results = analyze_workflows(workflows_dir)

    # # --- Stage 3: Runtime Profiling ---
    # print("[INFO] Profiling pipeline execution...")
    # profile_results = profile_pipeline(profile_entrypoint)

    # --- Stage 4: Generate Combined Report ---
    report_path = os.path.join(output_dir, f"assessment_{timestamp}.html")
    # generate_report(static_results, workflow_results, profile_results, report_path)
    generate_report(static_results, report_path)

    # --- Save JSON ---
    json_path = os.path.join(output_dir, f"static_results_{timestamp}.json")
    with open(json_path, "w") as f:
        json.dump(static_results, f, indent=2)
    print(f"[INFO] Static results saved → {json_path}")

    summarize_complexity(static_results, top_n=TOP_N)

    print(f"[INFO] Assessment complete → {report_path}")

