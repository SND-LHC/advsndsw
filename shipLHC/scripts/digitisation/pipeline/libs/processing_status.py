from datetime import datetime
import pandas as pd
import logging

def file_modified_time(path):
    if path.exists():
        return datetime.fromtimestamp(path.stat().st_mtime)
    return None
    

def create_status_table(directories):
    rows = []

    # Skip non valid runs in tb 2026
    runs = [p for p in directories["raw"].glob("run*") if int(p.name.replace("run", "")) >= 273]
    for run_dir in sorted(runs):

        if not run_dir.is_dir():
            continue

        run_name = run_dir.name
        run_number = run_name.replace("run", "")

        # Count files
        raw_files = list(run_dir.glob("*index*.raw"))
        jsn_files = list(run_dir.glob("*index*.jsn"))

        n_raw = len(raw_files)
        n_jsn = len(jsn_files)

        # Corresponding files
        converted_file = directories["converted"] / run_name / f"{run_name}_converted.root"
        digi_file = directories["converted"] / run_name / f"{run_name}_digi_rntuple.root"
        dqm_file = directories["histos"] / f"{run_name}_dqm.root"

        # Existence checks
        converted_exists = converted_file.exists()
        digi_exists = digi_file.exists()
        dqm_exists = dqm_file.exists()

        # Timestamp comparison
        converted_mtime = file_modified_time(converted_file)
        digi_mtime = file_modified_time(digi_file)
        dqm_mtime = file_modified_time(dqm_file)

        converted_up_to_date = ((n_raw == n_jsn) and converted_exists)
        digi_up_to_date = False
        dqm_up_to_date = False

        if converted_exists and digi_exists:
            digi_up_to_date = (digi_mtime >= converted_mtime and (n_raw == n_jsn))

        if converted_exists and digi_exists and dqm_exists:
            dqm_up_to_date = (dqm_mtime >= converted_mtime and dqm_mtime >= digi_mtime and (n_raw == n_jsn))

        rows.append({
            "run": run_number,
            "n_raw": n_raw,
            "n_jsn": n_jsn,
            "converted_modified_time": converted_mtime,
            "digi_modified_time": digi_mtime,
            "dqm_modified_time": dqm_mtime,
            "converted_up_to_date": converted_up_to_date,
            "digi_up_to_date": digi_up_to_date,
            "dqm_up_to_date": dqm_up_to_date,
        })
    
    df = pd.DataFrame(rows)
    if not df.empty:
        df = df.sort_values("run")
    return df

def select_run(df):
    # There cannot be less raw than jsn files
    invalid = df[df["n_raw"] < df["n_jsn"]]

    if not invalid.empty:
        logging.error("Found invalid run with n_raw < n_jsn")
        logging.error("\n%s", invalid[["run", "n_raw", "n_jsn"]].to_string(index=False))

    candidates = df[(~df["converted_up_to_date"]) | (~df["digi_up_to_date"]) | (~df["dqm_up_to_date"])]

    if candidates.empty:
        return None
    
    # Select earliest candidate
    selected_row = candidates.loc[candidates["run"].idxmin()]
    return selected_row.to_dict()