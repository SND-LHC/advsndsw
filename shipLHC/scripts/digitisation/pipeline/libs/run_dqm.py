import subprocess
import logging
import time
import os

def run_dqm(directories, run_number):
    tag = f"[run {run_number:06d}]"

    input_root_file = directories['converted'] / f"run{run_number:06d}" / f"run{run_number:06d}_digi_rntuple.root"
    output_root_file = directories['histos'] / f"run{run_number:06d}_dqm.root"

    advsndsw_base = os.environ["ADVSNDSW_ROOT"]

    source_advsndsw = "source /cvmfs/sndlhc.cern.ch/SNDLHC-2025/Oct7/setUp.sh"
    alienv = f"cd {directories['advsndsw']} && alienv enter advsndsw/latest"

    # Different tb have different mapping
    detinfo_csv = os.path.join(advsndsw_base, "shipLHC/digitisation/rawToDigi/mapping/detector_info_tb_5_2026.csv") if run_number <= 386 else os.path.join(advsndsw_base, "shipLHC/digitisation/rawToDigi/mapping/detector_info_tb_7_2026.csv")
    executable = os.path.join(advsndsw_base, "bin/run_real_time_monitoring")

    real_time_monitoring = f"{executable} {input_root_file} {detinfo_csv} {directories['geometry']} {output_root_file} 2"

    command = (f"{source_advsndsw} && {alienv} && {real_time_monitoring}")

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