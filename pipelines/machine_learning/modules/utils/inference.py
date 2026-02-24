def load_production_model(model_dir):
    import joblib
    from catboost import CatBoostRegressor
    import json
    from pathlib import Path

    model_dir = Path(model_dir)

    with open(model_dir / "metadata.json") as f:
        metadata = json.load(f)

    properties = metadata["properties"]
    n_ensemble = metadata["n_ensemble"]

    pca = joblib.load(model_dir / "pca.joblib")

    models = {}

    for prop in properties:
        prop_models = []
        for i in range(n_ensemble):
            m = CatBoostRegressor()
            m.load_model(str(model_dir / f"{prop.replace('/', '_')}_model_{i}.cbm"))
            prop_models.append(m)
        models[prop] = prop_models

    return models, pca, properties