#!/bin/bash -e

#-----------------------------------------------------------------------------------------
# Script information
script_name='APP - RUNNER FLOMART - MODELLED WS DISCHARGE - REALTIME'
script_version="2.0.1"
script_date='2021/11/18'

virtualenv_folder='/hydro/library/fp_libs_python3/'
virtualenv_name='virtualenv_python3'
script_folder='/hydro/library/fp_package_flomart_application/'

# Execution example:
# python3 app_flood_scenario_main.py -settings_algorithm configuration.json -time "2020-11-02 12:00"
#-----------------------------------------------------------------------------------------

#-----------------------------------------------------------------------------------------
# Get file information
script_file='/hydro/library/fp_package_flomart_application/app_flomart_main/app_flomart_main.py'
settings_file='/hydro/fp_tools_postprocessing/flomart_app_execution/Chienti/app_flomart_configuration_nwp_ecmwf0100_realtime_CHIENTI.json'

script_file_transfer='/hydro/library/fp_package_hyde/tools/tool_processing_datasets_transfer/hyde_tools_transfer_datasets.py'
settings_file_transfer='/hydro/fp_tools_postprocessing/flomart_app_execution/Chienti/app_flomart_transfer_nwp_ecmwf0100_realtime_CHIENTI.json'

# Time period execution
time_period_hour=0 # hour(s)

# Get information (-u to get gmt time)
time_now=$(date -u +"%Y-%m-%d %H:00")
#time_now="2013-11-11 04:00" # DEBUG
#-----------------------------------------------------------------------------------------

#-----------------------------------------------------------------------------------------
# Activate virtualenv
export PATH=$virtualenv_folder/bin:$PATH
source activate $virtualenv_name

# Add path to pythonpath
export PYTHONPATH="${PYTHONPATH}:$script_folder"
#-----------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------
# Info script start
echo " ==================================================================================="
echo " ==> "$script_name" (Version: "$script_version" Release_Date: "$script_date")"
echo " ==> START ..."

# Iterate over hours
time_run=$(date -d "$time_now" +'%Y-%m-%d %H:00')
for time_period_step in $(seq 0 $time_period_hour); do
    
    # Parse time information
    time_step=$(date -d "$time_run ${time_period_step} hour ago" +'%Y-%m-%d %H:00')

	# Run python script (using setting and time)
	echo " ===> COMPUTE FLOOD SCENARIOS [${time_step}] ... "
	echo " ===> COMMAND LINE: " python $script_file -settings_file $settings_file -time $time_step

	python $script_file -settings_file $settings_file -time "$time_step"
	
	echo " ===> COMPUTE FLOOD SCENARIOS [${time_step}] ... DONE"
	
	#---------------------------------------------------------------------------------
	# Run python model instance (using settings and time)
	echo " ===> RUN FLOMART HAZARD MAPS TRANSFER ... "
	# Run python script (using setting and time)
	python3 $script_file_transfer -settings_file $settings_file_transfer -time "$time_step"
	echo " ===> RUN FLOMART HAZARD MAPS TRANSFER ... DONE"
# ----------------------------------------------------------------------------------------

done
# ----------------------------------------------------------------------------------------



# ----------------------------------------------------------------------------------------
# Info script end
echo " ==> "$script_name" (Version: "$script_version" Release_Date: "$script_date")"
echo " ==> ... END"
echo " ==> Bye, Bye"
echo " ==================================================================================="
# ----------------------------------------------------------------------------------------

