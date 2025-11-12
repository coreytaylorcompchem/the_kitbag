from pipeline.logger import setup_logger

logger = setup_logger(
    __name__,
    debug_mode=False,
    simple_format=True
)

def update_config_input_file(config, new_input_path):
    """
    Update a workflow config so that downstream tasks automatically use a input files specified in the yaml.

    This updates:
      • Top-level 'input_file' key
      • Any nested stage dicts that contain 'input_file' (e.g. under task sections)

    Also logs which entries were updated for easier debugging.
    """
    if not new_input_path:
        logger.warning("[update_config_input_file] Called with no path — config not modified.")
        return config

    new_input_path = str(new_input_path)
    changes = []

    # Update top-level
    old_top = config.get("input_file")
    config["input_file"] = new_input_path
    if old_top != new_input_path:
        changes.append(f"top-level: '{old_top}' → '{new_input_path}'")

    # Update nested keys
    for key, value in config.items():
        if isinstance(value, dict):
            if "input_file" in value:
                old_nested = value["input_file"]
                value["input_file"] = new_input_path
                if old_nested != new_input_path:
                    changes.append(f"{key}.input_file: '{old_nested}' → '{new_input_path}'")

            # Optional: go one level deeper (if you have nested stage configs)
            for subkey, subval in value.items():
                if isinstance(subval, dict) and "input_file" in subval:
                    old_deep = subval["input_file"]
                    subval["input_file"] = new_input_path
                    if old_deep != new_input_path:
                        changes.append(f"{key}.{subkey}.input_file: '{old_deep}' → '{new_input_path}'")

    if changes:
        logger.info("[update_config_input_file] Updated input paths:\n  " + "\n  ".join(changes))
    else:
        logger.debug("[update_config_input_file] No input_file keys required updating.")

    return config

