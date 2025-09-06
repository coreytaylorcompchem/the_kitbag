from task_registry import register_task, get_task
import tempfile
import os
import json
import datetime

@register_task(
    'pipeline',
    description="Allows us to stack calculations in a workflow",
    modifies_geometry=False,
)
def run_pipeline(backend, xyz_file, config):
    steps = config.get('pipeline', [])
    if not steps:
        raise ValueError("No steps defined in 'pipeline'.")

    checkpoint_dir = config.get('checkpoint_dir', '.pipeline_checkpoints')
    os.makedirs(checkpoint_dir, exist_ok=True)

    pipeline_log_path = os.path.join(checkpoint_dir, 'pipeline.log')
    pipeline_log_lines = []
    def log_global(msg):
        print(msg)
        timestamped = f"{datetime.datetime.now().isoformat()} | {msg}"
        pipeline_log_lines.append(timestamped)

    current_xyz = xyz_file

    for i, step in enumerate(steps, 1):
        task_name = step.get('task')
        if not task_name:
            raise ValueError(f"Step {i} is missing a 'task'.")

        step_name = f"step_{i}_{task_name.lower()}"
        step_xyz = os.path.join(checkpoint_dir, f"{step_name}.xyz")
        step_log = os.path.join(checkpoint_dir, f"{step_name}.log")
        step_meta = os.path.join(checkpoint_dir, f"{step_name}.json")

        step_log_lines = []
        def log_step(msg):
            print(msg)
            timestamped = f"{datetime.datetime.now().isoformat()} | {msg}"
            step_log_lines.append(timestamped)
            pipeline_log_lines.append(timestamped)

        if os.path.exists(step_xyz) and not step.get('force', False):
            log_step(f"[SKIP] Step {i} ({task_name}) — checkpoint exists.")
            current_xyz = step_xyz
            write_log(step_log, step_log_lines)
            continue

        log_step(f"[RUN] Step {i}: {task_name}")
        task_func = get_task(task_name)
        if not task_func:
            raise ValueError(f"Task '{task_name}' is not registered.")

        step_config = config.copy()
        step_config.update(step)

        try:
            result = task_func(backend, current_xyz, step_config)
            success = True
            message = "Task completed successfully"
        except Exception as e:
            result = None
            success = False
            message = f"Task failed: {e}"
            log_step(f"[ERROR] {message}")

        if isinstance(result, str) and result.strip().startswith(str(len(result.splitlines()) - 2)):
            with open(step_xyz, 'w') as f:
                f.write(result)
            current_xyz = step_xyz
            log_step(f"[DONE] Geometry saved to {step_xyz}")
        elif success:
            log_step(f"[DONE] Task produced no new geometry.")
        else:
            log_step(f"[FAIL] Skipping remaining steps.")
            write_log(step_log, step_log_lines)
            break

        with open(step_meta, 'w') as f:
            json.dump({
                'step': i,
                'task': task_name,
                'success': success,
                'message': message,
                'config': step_config,
                'timestamp': datetime.datetime.now().isoformat()
            }, f, indent=2)

        write_log(step_log, step_log_lines)

    # Final log
    log_global(f"[Pipeline] Completed. Final geometry: {current_xyz}")
    write_log(pipeline_log_path, pipeline_log_lines)


def write_log(log_path, lines):
    with open(log_path, 'w') as f:
        f.write("\n".join(lines))

