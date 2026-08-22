from pathlib import Path
import subprocess
import logging
import time

def run_conversion(directories, run_number):
    tag = f"[run {run_number:06d}]"
    raw_folder = (directories['raw'] / f'run{run_number:06d}')

    current_dir = Path.cwd()
    source_cms = "source /cvmfs/cms.cern.ch/cmsset_default.sh"
    cmsenv = f"cd {directories['cmssw_src']} && cmsenv && cd {current_dir}"
    # If data is on eos, this is needed to release locks on files
    clean_previous_conversion = f"rm -rf {raw_folder / '*index*.jsn'} {raw_folder / 'processing'} {raw_folder / 'open'} {raw_folder / 'mon'} {raw_folder / '*BoLS.jsn'} {raw_folder / 'fu.lock'}"
    reset_progress = f"printf '1 0' > {raw_folder / 'fu.lock'}"
    conversion = f"cmsRun cms_raw_evt_building.py rawDirectory={directories['raw']} convertedDirectory={directories['converted']} runNumber={run_number}"

    command = f"{source_cms} && {cmsenv} && {clean_previous_conversion} && {reset_progress} && {conversion}" 

    logging.debug("%s Running converion command: %s", tag, command)

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
        logging.info("%s Conversion subprocess stdout:\n%s", tag, result.stdout)

    if result.stderr:
        logging.error("%s Conversion subprocess stderr:\n%s", tag, result.stderr)

    logging.info("%s Conversion finished in %.2f seconds", tag, duration)