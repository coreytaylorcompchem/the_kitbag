import yaml
import os
import networkx as nx

def analyze_workflows(workflows_dir):
    graph = nx.DiGraph()
    metrics = []

    for wf_file in os.listdir(workflows_dir):
        if wf_file.endswith((".yml", ".yaml")):
            path = os.path.join(workflows_dir, wf_file)
            with open(path) as f:
                wf = yaml.safe_load(f)

            tasks = wf.get("tasks", [])
            for t in tasks:
                name = t["name"]
                deps = t.get("depends_on", [])
                graph.add_node(name)
                for d in deps:
                    graph.add_edge(d, name)

    metrics.append({
        "n_tasks": graph.number_of_nodes(),
        "n_dependencies": graph.number_of_edges(),
        "max_depth": nx.dag_longest_path_length(graph) if graph else 0,
    })

    return {"graph": graph, "metrics": metrics}
