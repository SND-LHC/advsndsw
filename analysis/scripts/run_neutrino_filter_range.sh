#!/bin/bash

geo_file="/eos/experiment/sndlhc/convertedData/physics/2022/geofile_sndlhc_TI18_V0_2022.root"

INPUT_DIR=$1
BASE_OUT_DIR=$2
RUN_RANGE_START=$3
RUN_RANGE_END=$4
MODE=$5
CUT_SET=$6

mode_list=(STAGE1 RECO STAGE2)

if [[ ! ${mode_list[@]} =~ $MODE ]]
then
    echo Mode $MODE not available. It must be one of "${mode_list[*]}"
    echo Exitting.
    exit
fi

for i_run in `seq ${RUN_RANGE_START} ${RUN_RANGE_END}`
do

    if [ "$MODE" == "STAGE1" ]
    then
	# Check if only one digitized file exists in the input directory. Otherwise, skip.
	digi_files=(${INPUT_DIR}/${i_run}/*_digCPP.root)

	if [ "${#digi_files[@]}" -ne "1" ]
	then
	    echo "Found ${#digi_files[@]} digitized files in ${INPUT_DIR}/${i_run}. SKIPPING directory"
	    continue
	fi

	# Input file name
	input_file=${digi_files[0]}

	if [ -f "${BASE_OUT_DIR}/${i_run}/filtered_MC_00${i_run}.root" ]
	    then
	    echo "File ${BASE_OUT_DIR}/${i_run}/filtered_MC_00${i_run}.root exists. Skipping."
	    continue
	fi

    	# Set up SND environment
	if [ -z ${SNDSW_ROOT+x} ]
	then
	    echo "Setting up SNDSW"
    	    SNDBUILD_DIR=/eos/user/c/cvilela/SND_Dec16/sw/
    	    source /cvmfs/sndlhc.cern.ch/SNDLHC-2022/July14/setUp.sh
    	    eval `alienv load -w $SNDBUILD_DIR --no-refresh sndsw/latest-feature-second_stage_nu_filter-release`
	    export EOSSHIP=root://eosuser.cern.ch/
	fi

	# Run first stage filter
	neutrinoFilterGoldenSample ${input_file} filtered_MC_00${i_run}.root $CUT_SET

	# Copy output
    	mkdir -p ${BASE_OUT_DIR}/${i_run}/
    	xrdcp -f ./filtered_MC_00${i_run}.root ${BASE_OUT_DIR}/${i_run}/
    	rm -rf ./filtered_MC_00${i_run}.root
    elif [ "$MODE" == "RECO" ]
    then

	if [ -f "${BASE_OUT_DIR}/${i_run}/filtered_MC_00${i_run}_muonReco.root" ]
	    then
	    echo "File ${BASE_OUT_DIR}/${i_run}/filtered_MC_00${i_run}_muonReco.root exists. Skipping."
	    continue
	fi

	if [ ! -f "${BASE_OUT_DIR}/${i_run}/filtered_MC_00${i_run}.root" ]
	    then
	    echo "File ${BASE_OUT_DIR}/${i_run}/filtered_MC_00${i_run}.root does not exist. Skipping."
	    continue
	fi

    	# Set up SND environment
	if [ -z ${SNDSW_ROOT+x} ]
	then
	    echo "Setting up SNDSW"
#    	    SNDBUILD_DIR=/eos/user/c/cvilela/SND_Dec16/sw/
    	    SNDBUILD_DIR=/eos/user/c/cvilela/SND_153f764/sw/ # For PRL analysis use this version for RECO only
    	    source /cvmfs/sndlhc.cern.ch/SNDLHC-2022/July14/setUp.sh
    	    eval `alienv load -w $SNDBUILD_DIR --no-refresh sndsw/latest`
	    export EOSSHIP=root://eosuser.cern.ch/
	fi

	python $SNDSW_ROOT/shipLHC/run_muonRecoSND.py -f ${BASE_OUT_DIR}/${i_run}/filtered_MC_00${i_run}.root -g ${geo_file} -c passing_mu_DS -sc 1 -s ./ -hf linearSlopeIntercept -o

    	xrdcp -f ./filtered_MC_00${i_run}_muonReco.root ${BASE_OUT_DIR}/${i_run}/
	rm ./filtered_MC_00${i_run}_muonReco.root

    elif [ "$MODE" == "STAGE2" ]
    then

	if [ -f "${BASE_OUT_DIR}/${i_run}/filtered_MC_00${i_run}_stage2_noscifi2.root" ]
	    then
	    echo "File ${BASE_OUT_DIR}/${i_run}/filtered_MC_00${i_run}_stage2_noscifi2.root exists. Skipping."
	    continue
	fi

	if [ ! -f "${BASE_OUT_DIR}/${i_run}/filtered_MC_00${i_run}.root" ]
	    then
	    echo "File ${BASE_OUT_DIR}/${i_run}/filtered_MC_00${i_run}.root does not exist. Skipping."
	    continue
	fi

    	# Set up SND environment
	if [ -z ${SNDSW_ROOT+x} ]
	then
	    echo "Setting up SNDSW"
    	    SNDBUILD_DIR=/eos/user/c/cvilela/SND_Dec16/sw/
    	    source /cvmfs/sndlhc.cern.ch/SNDLHC-2022/July14/setUp.sh
    	    eval `alienv load -w $SNDBUILD_DIR --no-refresh sndsw/latest`
	    export EOSSHIP=root://eosuser.cern.ch/
	fi

	python ${SNDSW_ROOT}/analysis/neutrinoFilterGoldenSample_stage2.py -f ${BASE_OUT_DIR}/${i_run}/filtered_MC_00${i_run}.root -t ${BASE_OUT_DIR}/${i_run}/filtered_MC_00${i_run}_muonReco.root -o ./filtered_MC_00${i_run}_stage2_noscifi2.root -g ${geo_file};

    	xrdcp -f ./filtered_MC_00${i_run}_stage2_noscifi2.root ${BASE_OUT_DIR}/${i_run}/
	rm ./filtered_MC_00${i_run}_stage2_noscifi2.root
    fi
done
