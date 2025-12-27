from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

from modules.utils.seeding import expand_seeds
from backends.jobs import InferenceJob

from pipeline.logger import setup_logger

logger = setup_logger(__name__, debug_mode=True, simple_format=True)


def _run_job(job: InferenceJob):
    """
    Runs inside a subprocess.
    """
    return job.tool_backend.run(
        run_id=job.run_id,
        seed=job.seed,
        device=job.device,
        output_dir=job.output_dir,
        refinement=job.refinement,
    )


class StructureInferenceBackend:
    def __init__(
        self,
        sequence: str,
        tools: dict,
        config: dict,
        output_dir: Path,
        available_gpus=(0, 1, 2),
    ):
        self.sequence = sequence
        self.tools = tools
        self.config = config
        self.output_dir = output_dir
        self.available_gpus = list(available_gpus)

    def _build_jobs(self):
        jobs = []

        inference_cfg = self.config["structure_prediction"]["inference"]
        refinement_cfg = self.config["structure_prediction"].get("refinement")

        global_seeds = expand_seeds(
            inference_cfg["num_runs"],
            inference_cfg.get("base_seed"),
            inference_cfg.get("seed_strategy", "increment"),
        )

        for tool_name, tool in self.tools.items():
            tool_cfg = self.config["structure_prediction"]["tools"].get(tool_name, {})
            if not tool_cfg.get("enabled", False):
                continue

            tool_inference = tool_cfg.get("inference", {})
            seeds = expand_seeds(
                tool_inference.get("num_runs", len(global_seeds)),
                tool_inference.get("base_seed", inference_cfg.get("base_seed")),
                tool_inference.get("seed_strategy", inference_cfg.get("seed_strategy")),
            )

            for run_id, seed in enumerate(seeds):
                jobs.append(
                    InferenceJob(
                        tool_name=tool_name,
                        tool_backend=tool,
                        run_id=run_id,
                        seed=seed,
                        device=None,  # assigned later
                        output_dir=self.output_dir / tool_name / f"run_{run_id}",
                        refinement=refinement_cfg if refinement_cfg and refinement_cfg.get("enabled") else None,
                    )
                )

        return jobs

    def _assign_devices(self, jobs):
        for i, job in enumerate(jobs):
            job.device = self.available_gpus[i % len(self.available_gpus)]

    def run(self):
        jobs = self._build_jobs()
        self._assign_devices(jobs)

        logger.info(f"Launching {len(jobs)} inference jobs "
                    f"across GPUs {self.available_gpus}")

        results = []

        with ProcessPoolExecutor(max_workers=len(self.available_gpus)) as pool:
            futures = {pool.submit(_run_job, job): job for job in jobs}

            for future in as_completed(futures):
                job = futures[future]
                try:
                    result = future.result()
                    results.append({
                        "tool": job.tool_name,
                        "run_id": job.run_id,
                        "seed": job.seed,
                        "device": job.device,
                        **result,
                    })
                except Exception as e:
                    logger.error(
                        f"Job failed: {job.tool_name} run {job.run_id} "
                        f"(GPU {job.device}): {e}"
                    )
                    raise

        return results
