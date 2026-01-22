import importlib

import pandas as pd

from torch_geometric.loader import DataLoader
from sklearn.model_selection import train_test_split

from pipeline.task_registry import register_task
from models.cyp3a4.featurisation import mol_to_graph

from pipeline.logger import setup_logger

logger = setup_logger(__name__, debug_mode=False, simple_format=True)

@register_task("load_smiles_dataset", category="ADME")
def load_smiles_dataset(config, context):
    csv_path = config["csv_path"]

    df = pd.read_csv(csv_path)
    context["dataframe"] = df

    return {"dataframe": df}

@register_task("featurise_smiles", category="ADME")
def featurise_smiles(config, context):
    df = context["dataframe"]

    smiles_col = config.get("smiles_col", "smiles")
    label_col = config.get("label_col", "pIC50")

    # Load featuriser dynamically
    feat_cfg = config["featuriser"]
    module = importlib.import_module(feat_cfg["module"])
    featuriser = getattr(module, feat_cfg["function"])

    graphs = []
    for smi, y in zip(df[smiles_col], df[label_col]):
        g = featuriser(smi, label=y)
        if g is not None:
            graphs.append(g)

    if not graphs:
        raise ValueError("Featurisation produced zero graphs.")

    # Infer dimensions for downstream tasks
    context.update({
        "graphs": graphs,
        "input_dim": graphs[0].x.shape[1],
        "global_feat_dim": getattr(graphs[0], "global_features", None).shape[0]
    })

    return context

@register_task("split_data", category="ADME")
def split_data(config, context):
    graphs = context["graphs"]

    test_size = config.get("val_size", 0.2)
    batch_size = config.get("batch_size", 32)

    train_list, val_list = train_test_split(
        graphs,
        test_size=test_size,
        random_state=config.get("seed", 42)
    )

    train_loader = DataLoader(train_list, batch_size=batch_size, shuffle=True)
    val_loader = DataLoader(val_list, batch_size=batch_size, shuffle=False)

    context.update({
        "train_loader": train_loader,
        "val_loader": val_loader
    })

    return {
        "train_loader": train_loader,
        "val_loader": val_loader
    }

@register_task("load_model_spec", category="ADME")
def load_model_spec(config, context):
    module_path = config["module"]
    class_name = config["class"]

    module = importlib.import_module(module_path)
    ModelClass = getattr(module, class_name)

    context["model_class"] = ModelClass
    return {"model_class": ModelClass}
