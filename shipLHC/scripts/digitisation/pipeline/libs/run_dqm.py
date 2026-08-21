import subprocess
import logging
import time

def run_dqm(directories, run_number):
    tag = f"[run {run_number:06d}]"

    input_root_file = directories['converted'] / f"run{run_number:06d}" / f"run{run_number:06d}_digi_rntuple.root"
    output_root_file = directories['histos'] / f"run{run_number:06d}_dqm.root"

    source_advsndsw = "source /cvmfs/sndlhc.cern.ch/SNDLHC-2025/Oct7/setUp.sh"
    alienv = f"cd {directories['advsndsw']} && alienv enter advsndsw/latest"

    # Different tb have different mapping
    mapping_file = (
        "shipLHC/digitisation/rawToDigi/mapping/detector_info_tb_5_2026.csv"
        if run_number <= 386
        else "shipLHC/digitisation/rawToDigi/mapping/detector_info_tb_7_2026.csv"
    )

    command = f"""
        {source_advsndsw} &&
        {alienv} &&
        executable="$ADVSNDSW_ROOT/bin/run_real_time_monitoring" &&
        detinfo_csv="$ADVSNDSW_ROOT/{mapping_file}" &&
        "$executable" "{input_root_file}" "$detinfo_csv" "{directories['geometry']}" "{output_root_file}" 2
        """

    logging.debug("%s Running DQM command: %s", tag, command)

    start = time.perf_counter()
    
    result = subprocess.run(
        command,
        shell=True,
        executable="/bin/bash",
        capture_output=True,
        text=True
    )

    duration = time.perf_counter() - start

    if result.stdout:
        logging.info("%s DQM subprocess stdout:\n%s", tag, result.stdout)

    if result.stderr:
        logging.error("%s DQM subprocess stderr:\n%s", tag, result.stderr)

    logging.info("%s DQM finished in %.2f seconds", tag, duration)