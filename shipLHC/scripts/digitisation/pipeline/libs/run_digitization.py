import subprocess
import logging
import time
import os

def run_digitization(directories, run_number):
    tag = f"[run {run_number:06d}]"

    input_root_file = (directories['converted'] / f"run{run_number:06d}" / f"run{run_number:06d}_converted.root")
    advsndsw_base = os.environ["ADVSNDSW_ROOT"]

    source_advsndsw = "source /cvmfs/sndlhc.cern.ch/SNDLHC-2025/Oct7/setUp.sh"
    alienv = f"cd {directories['advsndsw']} && alienv enter advsndsw/latest"
    # Different tb have different mapping
    detinfo_csv = os.path.join(advsndsw_base, "shipLHC/digitisation/rawToDigi/mapping/detector_info_tb_5_2026.csv") if run_number <= 386 else os.path.join(advsndsw_base, "shipLHC/digitisation/rawToDigi/mapping/detector_info_tb_7_2026.csv")
    executable = os.path.join(advsndsw_base, "bin/run_raw_to_digi")
    modes = ["rntuple", "ttree"]

    for mode in modes:
        output_root_file = (directories['converted'] / f"run{run_number:06d}" / f"run{run_number:06d}_digi_{mode}.root")

        raw_to_digi = f"{executable} {input_root_file} {detinfo_csv} {output_root_file} {mode}"

        command = (f"{source_advsndsw} && {alienv} && {raw_to_digi}")

        logging.debug("%s Running Digitization (%s): %s", tag, mode, command)

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
            logging.info("%s [%s] stdout:\n%s", tag, mode, result.stdout)

        if result.stderr:
            logging.error("%s [%s] stderr:\n%s", tag, mode, result.stderr)

        logging.info("%s [%s] finished in %.2f seconds", tag, mode, duration)