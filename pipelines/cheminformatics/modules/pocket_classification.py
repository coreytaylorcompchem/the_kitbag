import os
import pandas as pd
import requests
from pathlib import Path
from rdkit import Chem
from pipeline.task_registry import register_task
from pipeline.logger import setup_logger
from pipeline.parallel_runner import ParallelWorkflowRunner

logger = setup_logger(__name__, debug_mode=False, simple_format=True)

def classify_pocket_location(pocket_name, tm_regions, ecd_regions):
    """
    Classify the location of the pocket as 7TM or ECD based on overlap with known regions.
    """
    if pocket_name in tm_regions:
        return "7TM"
    elif pocket_name in ecd_regions:
        return "ECD"
    else:
        return "Unknown"

def classify_functionality(pocket_name, gpcrdb_api_url):
    """
    Query GPCRdb API to classify the pocket's functionality (agonist, antagonist, etc.).
    """
    response = requests.get(f"{gpcrdb_api_url}/{pocket_name}")
    if response.status_code == 200:
        data = response.json()
        return data.get("functionality", "Unknown")
    else:
        logger.warning(f"Failed to fetch data for pocket {pocket_name}.")
        return "Unknown"

def classify_pockets(pockets_df, tm_regions, ecd_regions, gpcrdb_api_url, max_pockets=500):
    """
    Classify each pocket by location and functionality.
    """
    classified_pockets = []

    for _, row in pockets_df.iterrows():
        pocket_name = row['pocket_name']
        
        location = classify_pocket_location(pocket_name, tm_regions, ecd_regions)
        functionality = classify_functionality(pocket_name, gpcrdb_api_url)
        
        classified_pockets.append({
            "pocket_name": pocket_name,
            "location": location,
            "functionality": functionality
        })

        if len(classified_pockets) >= max_pockets:
            break

    return pd.DataFrame(classified_pockets)

class PocketClassificationRunner:
    def __init__(self, tm_regions, ecd_regions, gpcrdb_api_url, max_pockets):
        self.tm_regions = tm_regions
        self.ecd_regions = ecd_regions
        self.gpcrdb_api_url = gpcrdb_api_url
        self.max_pockets = max_pockets

    def __call__(self, task):
        pocket_chunk_file = task["pocket_chunk_file"]
        pockets_df = pd.read_csv(pocket_chunk_file)

        # Classify pockets and return the results
        return classify_pockets(
            pockets_df=pockets_df,
            tm_regions=self.tm_regions,
            ecd_regions=self.ecd_regions,
            gpcrdb_api_url=self.gpcrdb_api_url,
            max_pockets=self.max_pockets
        )

@register_task("pocket_classification", category="Pocket detection", description="Classify GPCR pockets by location and functionality.")
def gpcr_pocket_classification(config: dict, data: dict = None) -> dict:
    params = config.get("gpcr_pocket_classification", {})

    # Capture the file that fpocket generated for classification
    pocket_file = Path(params.get("input_file_pockets"))  # Should be passed correctly from fpocket
    tm_region_file = Path(params.get("input_file_7tm_regions"))
    ecd_region_file = Path(params.get("input_file_ecd_regions"))
    gpcrdb_api_url = params.get("gpcrdb_api_url")
    max_pockets = params.get("max_pockets_per_classification", 500)

    logger.info(f"Loading input files...")

    if not pocket_file.exists() or not tm_region_file.exists() or not ecd_region_file.exists():
        raise FileNotFoundError(f"Missing one or more input files: {pocket_file}, {tm_region_file}, {ecd_region_file}")

    pockets_df = pd.read_csv(pocket_file)
    tm_regions_df = pd.read_csv(tm_region_file)
    ecd_regions_df = pd.read_csv(ecd_region_file)

    logger.info(f"Loaded {len(pockets_df)} pockets, {len(tm_regions_df)} 7TM regions, {len(ecd_regions_df)} ECD regions.")

    tm_regions = tm_regions_df["region"].tolist()
    ecd_regions = ecd_regions_df["region"].tolist()

    output_cfg = config.get("output", {})
    output_dir = Path(output_cfg.get("directory", "outputs/gpcr_pocket_classification"))
    output_dir.mkdir(parents=True, exist_ok=True)
    combined_filename = output_cfg.get("filename", "classified_gpcr_pockets.csv")

    # Chunk the pockets for parallel processing (if large dataset)
    chunk_size = params.get("chunk_size", 1000)
    chunk_dir = output_dir / "pocket_chunks"
    chunk_dir.mkdir(parents=True, exist_ok=True)

    chunk_files = []
    for i in range(0, len(pockets_df), chunk_size):
        chunk = pockets_df.iloc[i:i + chunk_size].copy()
        chunk_file = chunk_dir / f"pockets_chunk_{i // chunk_size}.csv"
        chunk.to_csv(chunk_file, index=False, quoting=csv.QUOTE_ALL)
        chunk_files.append(str(chunk_file))

    task_list = []
    for chunk_file in chunk_files:
        task_list.append({"pocket_chunk_file": chunk_file})

    pocket_runner = PocketClassificationRunner(
        tm_regions=tm_regions,
        ecd_regions=ecd_regions,
        gpcrdb_api_url=gpcrdb_api_url,
        max_pockets=max_pockets
    )

    runner = ParallelWorkflowRunner(
        workflow_func=pocket_runner,
        config={"manual": task_list, **config},
        input_key="manual",
        output_key="classified_pockets",
        filename_pattern="classified_pockets_{chunk_idx}.csv",
        combined_filename=combined_filename,
        output_dir=str(output_dir),
        use_multiprocessing=True,
        reserve_cpus=4
    )

    logger.info(f"Starting parallel pocket classification over {len(task_list)} pocket chunks...")
    classified_df = runner.run()
    logger.info(f"Finished. Combined pocket classification count: {len(classified_df)}")

    return {"df": classified_df}
