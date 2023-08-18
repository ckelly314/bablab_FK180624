Analysis of FK180624 data is performed using the Python scripts in this folder.
Below, we describe the use for each script in the order that we recommend it be run. 
Raw data can be found on the BCO-DMO repository: https://www.bco-dmo.org/dataset-deployment/832393

1. clean_falkor.py 
    	Selects the measurement with 2 and 6 flags. Calculates DIC using PyCO2SYS and density using TEOS-10 GSW packages. Saves the analyzed variables and data as a .csv file.  

2. calc_oxycline_features.py
    Estimates organic matter C:N:P ratio and AOU in the oxycline using a robust least squares model between 45-95 dbar.

3. calc_reaction_stoi_and_R.py
    	Calculates and saves DNRN and denitrification reaction stoichiometric coefficients and R matrices for different organic matter compositions such as the estimate from FK180624 data, Redfield organic matter, Anderson organic matter, etc.

4. calc_falkor_outputs_func.py
	Calculates the relative reaction rates and relative contributions based on the R matrix input selected with Monte Carlo simuations. For the delta_tracer values, the script fits the tracer data from falkor_clean.csv file within each layer (defined in the script). Saves relative contributions as a text file for further plotting.

5. pH_and_reactions.py 
	Calculates the pH increase per organic carbon remineralized by denitrification for various organic matter compositions with varying carbon oxidation states for a range of DIC and TA values observed in the ETNP ODZ. Plots delta_pH versus carbon oxidation state.
	
6. pH_and_reactions_Cox.py 
	Calculates the pH increase per organic carbon remineralized by denitrification for the estimated organic matter composition based on FK180624 data. Maps delta_pH for a range of DIC and TA values.

7. calcium_carbonate_saturation.py
	Fits observed pH values from FK180624 data and creates a hypothetical pH profile. Calculates omega for using these pH profiles and in-situ temperature, alkalinity, salinity, etc. profiles. Plots both pH and omega as a function of depth. 
	

Additionally, the scripts described below are used to plot the observations from FK180624 and analysis results.
(See Cinay et al., 2022, GBC for details.)

1. O2_map_plot.py 
	Maps O2 concentration in the ETNP using Kwiecinski and Babbin 2021, GBC data products.

2. plot_station_no2_and_pH.py
	Plots observed nitrite and pH profiles for each of the 19 stations from FK180624 data as separate ridge-line plots.

3. plot_tracer_profiles.py
	Plots the vertical nitrite, pH, and oxygen profiles with respect to sigma0 for select stations. 

4. plot_st9.py
	Plots the vertical tracer profiles with respect to sigma0 for station 9. 

5. plot_multiple_p18_profiles.py
	Plots and compares vertical oxygen, nitrite, and pH profiles of station 9 to profiles from CLIVAR and GO-SHIP Cruises along P18 line at 110 W. 

6. plot_falkor_results_Cox.py
	Plots (i) the relative reaction rates (coefficients, i.e. Chi matrix in paper) as a heat map, (ii) residuals compared to the layer by layer tracer slope data, and (iii) relative contributions for a given output from calc_falkor_outputs.py. 

7. plot_Cox_outputs_compared.py
	Plots the relative contribution ratios such as Anammox % for a given range of carbon oxidation states (C_ox) as part of the sensitivity analysis. 
