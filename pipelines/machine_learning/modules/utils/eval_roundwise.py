from pathlib import Path
import pandas as pd

from modules.utils.eval_regression import evaluate_with_bootstrap, save_evaluation

from pathlib import Path
import pandas as pd

from modules.utils.eval_regression import evaluate_with_bootstrap, save_evaluation

def run_roundwise_evaluation(
    *,
    models,
    pca,
    X_all_pca,  # precomputed PCA embeddings for all sequences
    df_all,
    properties,
    round_col="round",
    out_dir,
    n_boot=500,
    min_n=20,
):
    """
    Runs per-round evaluation:
    - Train data: rounds <= r
    - Future data: rounds > r

    GPU-safe: uses precomputed PCA embeddings (X_all_pca)
    """
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    rounds = sorted(df_all[round_col].unique())
    summary_rows = []

    # Map sequence index for embeddings
    seq_to_idx = {seq: i for i, seq in enumerate(df_all["sequence"])}

    for r in rounds:
        df_train = df_all[df_all[round_col] <= r].reset_index(drop=True)
        df_future = df_all[df_all[round_col] > r].reset_index(drop=True)

        # Get precomputed embeddings for train/future
        train_indices = [seq_to_idx[s] for s in df_train["sequence"]]
        future_indices = [seq_to_idx[s] for s in df_future["sequence"]]

        X_train_pca = X_all_pca[train_indices]
        X_future_pca = X_all_pca[future_indices]

        # Train-set evaluation
        train_eval = evaluate_with_bootstrap(
            models=models,
            pca=None,  # embeddings are already PCA-transformed
            embed_fn=lambda seqs: X_train_pca,  # dummy, ignored
            df_eval=df_train,
            properties=properties,
            n_boot=n_boot,
            min_n=min_n,
        )
        save_evaluation(train_eval, out_dir / f"round{r}" / "train")

        # Future-only evaluation
        if len(df_future) >= min_n:
            future_eval = evaluate_with_bootstrap(
                models=models,
                pca=None,
                embed_fn=lambda seqs: X_future_pca,  # dummy
                df_eval=df_future,
                properties=properties,
                n_boot=n_boot,
                min_n=min_n,
            )
            save_evaluation(future_eval, out_dir / f"round{r}" / "future")
        else:
            future_eval = None

        # Flatten metrics for plotting
        for prop, m in train_eval.metrics.items():
            row = {
                "round": r,
                "property": prop,
                "train_spearman": m["spearman"].mean,
                "train_spearman_low": m["spearman"].low,
                "train_spearman_high": m["spearman"].high,
                "train_rmse": m["rmse"].mean,
                "train_rmse_low": m["rmse"].low,
                "train_rmse_high": m["rmse"].high,
                "n_train": m["n"],
            }

            if future_eval and prop in future_eval.metrics:
                mf = future_eval.metrics[prop]
                row.update({
                    "future_spearman": mf["spearman"].mean,
                    "future_spearman_low": mf["spearman"].low,
                    "future_spearman_high": mf["spearman"].high,
                    "n_future": mf["n"],
                })
            else:
                row.update({
                    "future_spearman": None,
                    "future_spearman_low": None,
                    "future_spearman_high": None,
                    "n_future": 0,
                })

            summary_rows.append(row)

    df_summary = pd.DataFrame(summary_rows)
    df_summary.to_csv(out_dir / "roundwise_metrics.csv", index=False)
    return df_summary
