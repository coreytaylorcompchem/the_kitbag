from pathlib import Path
from pipeline.task_registry import get_task

class BaseWorkflow:
    def __init__(self, config):
        self.config = config
        self.data = None

    def run(self):
        steps = self.config.get('workflow', [])
        for step in steps:
            print(f"\n>>> Running step: {step}")
            task_func = get_task(step)
            if not task_func:
                raise ValueError(f"Task '{step}' not found in registry.")
            self.data = task_func(self.config, self.data)
            self._log_step_output(step)

        self._save_output()
        return self.data

    def _log_step_output(self, step):
        if isinstance(self.data, dict) and "df" in self.data and hasattr(self.data["df"], "shape"):
            print(f"[{step}] → Output shape: {self.data['df'].shape}")
        elif hasattr(self.data, "shape"):
            print(f"[{step}] → Output shape: {self.data.shape}")
        else:
            print(f"[{step}] → Output data type: {type(self.data)}")

    def _save_output(self):
        output_cfg = self.config.get("output", {})
        out_dir = output_cfg.get("directory", "outputs/single_target")
        out_file = output_cfg.get("filename", "output.csv")
        overwrite = output_cfg.get("overwrite", False)

        Path(out_dir).mkdir(parents=True, exist_ok=True)
        out_path = Path(out_dir) / out_file

        if isinstance(self.data, dict) and "df" in self.data:
            df_to_save = self.data["df"]
        elif hasattr(self.data, "to_csv"):
            df_to_save = self.data
        else:
            print("No DataFrame to save.")
            return

        if out_path.exists() and not overwrite:
            print(f"File already exists at {out_path} and overwrite is False. Skipping save.")
        else:
            df_to_save.to_csv(out_path, index=False)
            print(f"Saved output to {out_path}")
