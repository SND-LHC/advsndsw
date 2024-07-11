advsndsw
=====

<!-- markdown-toc start - Don't edit this section. Run M-x markdown-toc-refresh-toc -->
**Table of Contents**

- [advsndsw](#advsndsw)
- [Introduction](#introduction)
    - [Contact and communication](#contact-and-communication)
    - [Branches](#branches)
- [Build instructions](#build-instructions)
    - [Introduction to `aliBuild`](#introduction-to-alibuild)
    - [On lxplus or systems with CVMFS](#on-lxplus-or-systems-with-cvmfs)
    - [On systems without access to CVMFS](#on-systems-without-access-to-cvmfs)
- [Run Instructions](#run-instructions)
    - [Use cases covered by `run_simSND.py`:](#use-cases-covered-by-run_simsndpy)
    - [Digitization of MC data:](#digitization-of-mc-data)
- [Docker Instructions](#docker-instructions)
- [Contributing Code](#contributing-code)

<!-- markdown-toc end -->

# Introduction

`advsndsw` is the software framework of the SND@LHC collaboration for the AdvSND upgrade, based on its current framework [`sndsw`](https://github.com/SND-LHC/sndsw). It is based on
the [FairShip](https://github.com/ShipSoft/FairShip) framework developed by the SHiP collaboration which in turn is
based on [FairRoot](https://github.com/FairRootGroup/FairRoot/), making use of the automatic python bindings provided by
PyROOT.

## Contact and communication

If you have questions or problems, please feel free to contact the
@SND-LHC/core-developers. For troubleshooting and development, we plan to discuss on [Mattermost](https://mattermost.web.cern.ch/sndlhc/channels/software).

The [snd-software](mailto:snd-software@cern.ch) mailing list can be used to discuss the software and report issues. Important annoucements will be made there.

The [snd-software-notifications](mailto:snd-software-notifications@cern.ch) mailing list will be used for automated notifications from Github and CI/CD etc..

Both mailing lists are self-subscribe CERN e-groups.

## Branches

<dl>
  <dt><code>master</code></dt>
  <dd>Main development branch.
      All python code is <b>required to be compatible with 3</b>
      Requires aliBuild default <code>release</code>.</dd>
</dl>

# Build instructions

The `aliBuild` family of tools developed by ALICE is used to set up `advsndsw` and
its dependencies.

## Introduction to `aliBuild`

The basic commands are the same regardless of whether CVMFS is used:

<dl>
  <dt><code>aliBuild build &lt;package-name&gt; -c snddist</code></dt>
  <dd>Build the package <code>&lt;package-name&gt;</code> (e.g. <code>advsndsw</code>) and its dependencies using the recipes and configuration provided by snddist.
  On CVMFS, it is recommended to add <code>--always-prefer-system</code> to ensure packages are used from CVMFS instead of being rebuilt.</dd>
  <dt><code>aliDoctor &lt;package-name&gt;  -c snddist</code></dt>
  <dd>Provide troubleshooting information and hints on which packages can be used from the system for <code>&lt;package-name&gt;</code> and its dependencies.</dd>
  <dt><code>alienv enter &lt;package-name&gt; /latest -c snddist</code></dt>
  <dd>Enter an environment with <code>&lt;package-name&gt;</code> and its dependencies.</dd>
</dl>

For more information on using `aliBuild`, see its
[documentation](https://alisw.github.io/alibuild/) (note: some things are ALICE
specific and will not apply to SND@LHC software).

## On lxplus or systems with CVMFS

On `lxplus` or any CC7/CC8 machine with access to CVMFS, you can do the following:

1. Make sure you can access the SNDLHC CVMFS Repository
    ```bash
    ls /cvmfs/sndlhc.cern.ch
    ```
2. Source the `setUp.sh` script
    ```bash
    source /cvmfs/sndlhc.cern.ch/SNDLHC-2024/June25/setUp.sh  # recommended latest version
    ```
3. If you don't want to modify the `advsndsw` package, skip step 3:
    ```bash
    git clone https://github.com/SND-LHC/advsndsw
    ```
   This gives you by default the `main` branch of the software. In case, you want to use a specific branch:
    ```bash
    cd advsndsw
    git checkout <branch>
    cd ..
    ```
4. Build the software using `aliBuild`
    ```bash
    aliBuild build advsndsw -c $SNDDIST --always-prefer-system --default release
    ```
If you exit your shell session and you want to go back working on it, make sure to re-execute the second step.

To load the `advsndsw` environment, after you build the software, you can simply use:

5. Load the environment
    ```bash
    alienv enter advsndsw/latest-main-release
    ```

However, this won't work if you are using HTCondor. In such case you can do:

```bash
eval $(alienv load advsndsw/latest-main-release --no-refresh)
```

If you modify `advsndsw`, simply repeat step 4 from `advsndsw`'s parent directory.

## On systems without access to CVMFS

Commands are similar to the previous case, but without access to CVMFS you need to build the required packages from source.

1. Clone the [snddist](https://github.com/SND-LHC/snddist), which containts the recipes to build `advsndsw` and it's dependencies:
    ```bash
    git clone https://github.com/SND-LHC/snddist.git
    ```
2.  If you don't want to modify the advsndsw package, skip step 2:
    ```bash
    git clone https://github.com/SND-LHC/advsndsw
    ```
    This gives you by default the master branch of the software. In case, you want to use a specific branch:
    ```bash
    cd advsndsw
    git checkout <branch>
    cd ..
    ```

3. Install [aliBuild](https://github.com/alisw/alibuild)
    ``` bash
    pipx install --user alibuild
    ```
    and make sure that it is in your $PATH.

4. Build the software using aliBuild
    ```bash
    aliBuild build advsndsw -c snddist --default release
    ```
    If you run into any problems, `aliDoctor` can help determine what the problem is.
5. Load the environment
    ```bash
    alienv enter advsndsw/latest-main-release
    ```

# Run Instructions

 Set up the bulk of the environment from CVMFS.

 ```bash
 source /cvmfs/sndlhc.cern.ch/SNDLHC-2024/June25/setUp.sh
 ```

 Load your local advsndsw environment (see above).

 ```bash
 alienv enter advsndsw/latest-main-release
 ```

  ```bash
 python $ADVSNDSW_ROOT/shipLHC/run_simSND.py  --Ntuple  -n 100 -f /eos/experiment/sndlhc/MonteCarlo/FLUKA/muons_up/version1/unit30_Nm.root  --eMin 1.0
 ```
 Output files are  sndLHC.Ntuple-TGeant4.root (events) and  geofile_full.Ntuple-TGeant4.root  (setup)

 Run the event display:

 ```bash
 python -i $ADVSNDSW_ROOT/macro/eventDisplay.py -f sndLHC.Ntuple-TGeant4.root -g geofile_full.Ntuple-TGeant4.root
 // use SHiP Event Display GUI
 Use quit() or Ctrl-D (i.e. EOF) to exit
 ```
 a) Use the GUI to display  events: SHiP actions / next event

 b) Hovering over trajectory will display additional information :

 c) At python prompt: sTree.MCTrack.Dump() will display info about all MC particles

## Use cases covered by `run_simSND.py`:

1. Transport muons, output of FLUKA simulation, to TI18 and the detector. Positive and negative muons, up and down crossing angles, exist.
Possible options are setting minimum energy for transporting particles, transport only muons, increase EM cross sections of muons.

 ```bash
 python $ADVSNDSW_ROOT/shipLHC/run_simSND.py  --Ntuple  -n nEvents  -f /eos/experiment/sndlhc/MonteCarlo/FLUKA/muons_up/version1/unit30_Nm.root  --eMin ecut
 ```

2. Muon deep inelastic scattering events, produced with pythia6, and then positioned in T18 and transported by Geant4:
 ```bash
 python $ADVSNDSW_ROOT/shipLHC/run_simSND.py  -F --MuDIS -n nEvents -f  /eos/experiment/sndlhc/MonteCarlo/Pythia6/MuonDIS/muonDis_1001.root  --eMin ecut
 ```
3. Neutrino events, produced by GENIE, and then positioned in T18 and transported by Geant4:
 ```bash
 python $ADVSNDSW_ROOT/shipLHC/run_simSND.py --Genie -n nEvents -f ...
 ```

## Digitization of MC data:

1. Convert MC points to detector hits. Input required, data from simulation together with the geometry file created when running simulation. New objects created are `Digi_ScifiHits` together with `Cluster_Scifi` and `Digi_MuFilterHit`, and in parallel objects to make the link to the original MC points, `Digi_MuFilterHits2MCPoints` and `Digi_ScifiHits2MCPoints`.

 ```bash
 python $ADVSNDSW_ROOT/shipLHC/run_digiSND.py -f sndLHC.Ntuple-TGeant4.root -g geofile_full.Ntuple-TGeant4.root
 ```

# Docker Instructions

Docker is **not** the recommended way to run `advsndsw` locally. It is ideal
for reproducing reproducible, stateless environments for debugging, HTCondor
and cluster use, or when a strict separation between `advsndsw` and the host is
desirable.

1. Build an docker image from the provided `Dockerfile`:
    ```bash
    git clone https://github.com/SND-LHC/advsndsw.git
    cd advsndsw
    docker build -t advsndsw .
    ```
2. Run the `advsndsw` docker image:
    ```bash
    docker run -i -t --rm advsndsw /bin/bash
    ```
3. Advanced docker run options:
    ```bash
    docker run -i -t --rm \
    -e DISPLAY=unix$DISPLAY -v /tmp/.X11-unix:/tmp/.X11-unix \
    -v /local\_workdir:/image\_workdir \
    advsndsw /bin/bash
    ```
    The option `-e DISPLAY=unix$DISPLAY -v /tmp/.X11-unix:/tmp/.X11-unix` forwards graphics from the docker to your local system (similar to `ssh -X`). The option `-v /local_workdir:/image_workdir` mounts `/local_workdir` on the local system as `/image_workdir` within docker.

# Contributing Code

All packages are managed in Git and GitHub. Please either use the web interface to create pull requests or issues, or [send patches via email](https://git-send-email.io/).

If your changes would also benefit [FairShip](https://github.com/ShipSoft/FairShip), please consider making a pull-request for your changes there. We can then pick up the changes from FairShip automatically.
