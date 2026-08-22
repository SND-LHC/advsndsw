import subprocess
import logging
import time

def run_digitization(directories, run_number):
    tag = f"[run {run_number:06d}]"

    input_root_file = (directories['converted'] / f"run{run_number:06d}" / f"run{run_number:06d}_converted.root")

    source_advsndsw = "source /opt/run4/software/setUp.sh"
    alienv = f"cd {directories['advsndsw']} && eval $(alienv load advsndsw/latest --no-refresh)"
    # Different tb have different mapping
    mapping_file = (
        "shipLHC/digitisation/rawToDigi/mapping/detector_info_tb_5_2026.csv"
        if run_number <= 386
        else "shipLHC/digitisation/rawToDigi/mapping/detector_info_tb_7_2026.csv"
    )
    modes = ["rntuple", "ttree"]

    for mode in modes:
        output_root_file = (directories['converted'] / f"run{run_number:06d}" / f"run{run_number:06d}_digi_{mode}.root")

        command = f"""
        {source_advsndsw} &&
        {alienv} &&
        executable="$ADVSNDSW_ROOT/bin/run_raw_to_digi" &&
        detinfo_csv="$ADVSNDSW_ROOT/{mapping_file}" &&
        "$executable" "{input_root_file}" "$detinfo_csv" "{output_root_file}" "{mode}"
        """

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