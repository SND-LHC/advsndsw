from pathlib import Path
import logging
import sys
import time

from libs.processing_status import select_run, create_status_table
from libs.run_conversion import run_conversion
from libs.run_digitization import run_digitization
from libs.run_dqm import run_dqm

def main():
    directories = {
        "raw" : Path("/eos/experiment/sndlhc/Run4/testbeam2026/raw_data"),
        "converted" : Path("/eos/experiment/sndlhc/Run4/testbeam2026/converted_data"),
        "histos" : Path("/eos/experiment/sndlhc/www/testbeam2026"),
        "geometry" : Path("/eos/experiment/sndlhc/Run4/testbeam2026/converted_data/geofile_testbeam_2026.root"),
        "logs" : Path("/eos/experiment/sndlhc/www/testbeam2026/logs"),
        "cmssw_src" : Path("/home/filippo/CMSSW_15_1_1/src"),
        "advsndsw" : Path("/home/filippo")
    }

    logging.basicConfig(
        level=getattr(logging, "DEBUG"),
        format="%(asctime)s [%(levelname)s] %(message)s",
        datefmt="%H:%M:%S",
        handlers=[
            logging.FileHandler(directories["logs"] / "processing_pipeline.log", mode="a"),  # append mode
            logging.StreamHandler(sys.stdout)          # still print to terminal
        ]
    )

    logging.info("Directories configuration:\n%s", "\n".join(f"[{k}] -> {v}" for k, v in directories.items()))

    POLL_INTERVAL = 60 # seconds

    while True:
        try:
            df = create_status_table(directories)
            output_csv_path = directories["logs"] / "processing_pipeline_status.csv"
            df.to_csv(output_csv_path, index=False)

            logging.debug("Current pipeline status:\n%s", df.to_string(index=False))
            logging.info("Saved pipeline status to: %s", output_csv_path)

            selected_run = select_run(df)

            if selected_run is not None:
                run_number = int(selected_run["run"])
                logging.info("Selected run: %s", run_number)
                if not selected_run["converted_up_to_date"]:
                    run_conversion(directories, run_number)
                    run_digitization(directories, run_number)
                    run_dqm(directories, run_number)
                elif not selected_run["digi_up_to_date"]:
                    run_digitization(directories, run_number)
                    run_dqm(directories, run_number)
                elif not selected_run["dqm_up_to_date"]:
                    run_dqm(directories, run_number)
                else:
                    logging.error("Run selected incorrectly")

            else:
                logging.info("Everything is up to date")
                time.sleep(POLL_INTERVAL)

        except KeyboardInterrupt:
            logging.info("Pipeline loop stopped by user")
            break

        except Exception:
            logging.exception("Error while processing pipeline loop")
            time.sleep(POLL_INTERVAL)

if __name__ == "__main__":
    main()