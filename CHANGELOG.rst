=========
Changelog
=========

Version 2.2.1 [2024-12-15]
**************************

- APP: app_flomart_main.py
	- add method to dump the reference area mask
	- fix bugs in geographical ancillary and destination boundaries (due to update of preprocessing tool)

Version 2.2.0 [2024-11-29]
**************************

- APP: app_flomart_main.py
	- add safe memory mode to manage the ram memory usage (big domain case)
	- extend codes to manage the info domain file
	- fix bugs in applyng regional_mode (using the q/t points)
	
- PREPROCESSING: 
    - update and modify the matlab preprocessing main script (app_flomart_preprocessing.m) to manage the big domain case

Version 2.1.1 [2024-05-30]
**************************

- APP: app_flomart_main.py
	- fix the case of missing data in the observed dataset (no links mode)

- DOCS: Flomart_Manual_20221130.pdf
- DOCS: Flomart_Training_20221121.pdf

- PALETTE: palette_water_level.qml (QGIS application)

	
Version 2.1.0 [2023-11-06]
**************************
- APP: app_flomart_main.py
	- improve and refactor codes, methods and classes to adapt code to expected new cases and formats
	- add method and classes to compute tr (method_regional, method_q_t)
	- add q(t) file to extract the value of tr as function of q
	- add discharge sources file (drift, hmc) and different formats and the sections file reader (for continuum txt data)
	- fix "tr_left" boundary lower and "tr_right" boundary upper; add function to check the array of tr to try the
	new settings in all cases and configurations
	- fix weights for tr greater than tr_max
	- fix numpy errors about int and float deprecated methods 
	
	
Version 2.0.5 [2022-12-15]
**************************
- APP: app_flomart_main.py
	- add time settings for both simulated and observed datasets
	- add datasets fields "time_period_search", "time_period_left" and "time_period_right" to define
	  the time periods for simulated and observed datasets
	- add tr methods parameters in the section workspace
	- add tr methods selection in the scenario definition
	- add checking for no data case in simulated and observed discharge methods

Version 2.0.4 [2022-12-12]
**************************
- APP: app_flomart_main.py
	- Modification of the function computing the relation between Q and T for the selection of the hazard map. 
	  One unique function is used (cmp_tr_general), instead of the two functions of previous versions (one for linear and one for exponential relation).
	  added of package scipy.interpolate for function interp1d. For instance, each case study has its own relation (and hard-coded defined parameters):
	  Foglia, Chienti (Marche) and Entella (Liguria) are implemented. For Liguria domain the variables section_discharge_idx and self.correction_discharge_factor
	  are still used, while for Marche domains they are not. Im the future versions this self.correction_discharge_factor (which is hardcoded now) must be put
	  on the input configuration files. 
	- Managed the problem of maps with 0 and NaN. In particular to linearly interpolate the values of the null map and the minimum map of the abacus.
	
- PREPROCESSING: 
        - final version of matlab preprocessing main script (app_flomart_preprocessing.m)
        - added python script (merging_and_compress_rasters.py) for merging and compress the rasters of the abacus, when needed.
        - added matlab script (test_sections_indexes_continuum.m) for testing the correct sections indexes of continuum. 

- CONFIGURATION FILES:
        - added configuration and bash files for operational use in Marche domain (Chienti river) with ECMWF forecasting chain: 
          app_flomart_configuration_nwp_ecmwf0100_realtime_CHIENTI.json,
          app_flomart_manager_nwp_ecmwf0100_realtime_CHIENTI.sh,
          app_flomart_transfer_nwp_ecmwf0100_realtime_CHIENTI.json,
          
Version 2.0.3 [2022-08-01]
**************************
- APP: app_flomart_main.py
	- Bugs fixing for reading geographical file in different type and format (mat, hdf5).
	- Bugs fixing for managing settings parameters.
	- Bugs fixing for reading dynamic files in json and mat format.

Version 2.0.2 [2022-07-28]
**************************
- APP: app_flomart_preprocessing.m 
	- updated preprocessing matlab script (new save of info_{domainname}.mat) 
	- computation of drainage areas (in km2) and epsg code
	- add of possibility to generate Qindex variable from continuum outputs.
- APP: example_generation_flood_map.m	 
	- add preprocessing script matlab to test flomart (the operational generation of realtime flood map from the abacus) 
- APP: Continuum_getMap_NC.m
	- add preprocessing script matlab to read netcf continuum outputs (for the computation of Qindex). 
- APP: driver_data_io_destination.py
	- bug corrected in inputs of function self.compute_scenario_discharge() and self.compute_scenario_tr()

Version 2.0.1 [2022-07-15]
**************************
- APP: app_flomart_configuration_deterministic_example_simulated_FOGLIA.json 
	- Add of configuration file for Foglia River case study. 

- APP: driver_data_io_destination.py
	- Modification to read all epsg through domain_epsg_code (info passed directly from .mat file) in order to plot png and tiff.

- APP: app_flomart_main/driver_data_io_source.py, app_flomart_main/driver_type_io.py
	- (DriverType class) modification in order to add variables_names=self.variables_obs or self.variables_sim 

Version 2.0.0 [2022-07-05]
**************************
- APP: app_flomart_main.py
	- Operational release and refactoring of scripts and procedures
- APP: app_flomart_preprocessing.m
	- Preprocessing tools to prepare data 

Version 1.9.1 [2022-02-28]
**************************
- APP: app_flood_scenario_main.py
	- Bugs fixing for the weighted scenarios part

Version 1.9.0 [2022-02-23]
**************************
- APP: app_flood_scenario_main.py
	- Bugs fixing and code updating based on pre-operational release

Version 1.8.0 [2022-02-01]
**************************
- APP: app_flood_scenario_main.py
	- Pre-operational release
	
Version 1.7.0 [2021-10-05]
**************************
- APP: app_flood_scenario_main.py
	- Generic release for correcting bugs and managing the empty datasets for obs/mod discharges
	
Version 1.6.0 [2021-05-15]
**************************
- APP: app_flood_scenario_main.py
	- Generic release for updating the tools and the modules

Version 1.5.0 [2021-04-12]
**************************
- APP: app_flood_scenario_main.py
	- Generic release for fixing unexpected bugs

Version 1.4.0 [2021-03-19]
**************************
- APP: app_flood_scenario_main.py
	- Generic release for fixing unexpected bugs

Version 1.3.0 [2021-02-02]
**************************
- APP: app_flood_scenario_main.py
	- Generic release for fixing unexpected bugs

Version 1.2.0 [2020-12-14]
**************************
- APP: app_flood_scenario_main.py
	- Test release for testing the tools and the modules

Version 1.1.0 [2020-11-25]
**************************
- APP: app_flood_scenario_main.py
    - Generic release for including methods, apps and tools of the previous experimental library (from MatLab scripts)

Version 1.0.0 [2020-05-22]
**************************
- APP: app_flood_scenario_main.py
    - Start development and configuration of flood scenario application
    - Python 3

