"""
File: clean_f5906484.py
---------------------

Selects the measurement with 2 and 6 flags.
Calculates DIC using PyCO2SYS and density using TEOS-10 GSW packages.
Saves the analyzed variables and data as a .csv file.

INPUT:
    :bottle_data.csv: .csv file with cruise data and the following columns:
        'STNNBR'
        'CASTNO'
        'LATITUDE'
        'LONGITUDE'
        'CTDPRS'
        'CTDTMP'
        'CTDSAL'
        'CTDS_FLAG_W'
        'CTDOXY'
        'CTDOXY_FLAG_W'
        'FLOR'
        'FLOR_FLAG_W'
        'NOx_FLAG_W'
        'NITRIT_BabLab'
        'NITRIT_FLAG_W'
        'PHSPHT'
        'PHSPHT_FLAG_W'
        'NH4_FLAG_W'
        'PH_TOT'
        'PH_TOT_FLAG_W'
        'PH_TEMP'
        'TA_FLAG_W'
OUTPUT:
    :f5906484_clean.csv: .csv file with bad flagged data filtered out and calculated params:
        -absolute salinity (calculated with gsw)
        -conservative temperature (calculated with gsw)
        -DIC (calculated with pyCO2SYS)
        -pH insitu (calculated with pyCO2SYS)
"""

## Import Libraries
import pandas as pd
import numpy as np
import gsw
import PyCO2SYS as pyco2

# Import Parameters and Results
# need to add na_values so that pyco2sys doesn't try to solve carbonate system with TA=-999
f5906484 = pd.read_csv("5906484qcno2.txt",skiprows = 70, sep='\t',
    parse_dates = ['mon/day/yr'])
# add 'month' column so we can match gridded WOA data
f5906484['month'] = f5906484['mon/day/yr'].dt.strftime('%m')
# need to rename columns from BCO-DMO file to match cols in clean_f5906484.py
cols = {
    'Station': "station",
    #"CASTNO": "cast",
    'Lat [°N]': "lat",
    'Lon [°E]': "lon",
    'Pressure[dbar]': "press",
    'Temperature[°C]': "temperature",
    'Salinity[pss]': "sal",
    'QF.3': "sal_flag",
    'Oxygen[µmol/kg]': "O2",
    'QF.6': "O2_flag",
    'Nitrate[µmol/kg]':"NO3",
    'QF.8':'NO3_flag',
    'Nitrite[µmol/kg]': "NO2",
    'QF.17': "NO2_flag",
    #"PHSPHT": "phosphate",
    #"PHSPHT_FLAG_W": "phosphate_flag",
    #"NH4_FLAG_W": "NH4_flag",
    'pHinsitu[Total]': "pH_insitu",
    'QF.12': "pH_tot_flag",
    'TALK_LIAR[µmol/kg]':'TA',
    'QF.14':"TA_flag",
    'DIC_LIAR[µmol/kg]':'DIC',
    'QF.15':'DIC_flag'
}
f5906484 = f5906484.rename(columns=cols)
f5906484 = f5906484.dropna()

# Select Flagged Bad Data

flag_good = np.array([0]) # WHAT ARE THE QUALITY FLAGS FOR THESE DATA
idx_flag = np.where(
    (np.isin(f5906484["sal_flag"], flag_good))
    & (np.isin(f5906484["O2_flag"], flag_good))
    & (np.isin(f5906484["NO3_flag"], flag_good))
    #& (np.isin(f5906484["phosphate_flag"], flag_good))
    #& (np.isin(f5906484["NH4_flag"], flag_good))
    & (np.isin(f5906484["pH_tot_flag"], flag_good))
    & (np.isin(f5906484["TA_flag"], flag_good))
    & (np.isin(f5906484["NO2_flag"], flag_good))
)

# Calculate Density
SA = gsw.SA_from_SP(f5906484["sal"], f5906484["press"], f5906484["lon"], f5906484["lat"])
CT = gsw.CT_from_t(SA, f5906484["temperature"], f5906484["press"])
f5906484["sigma0"] = gsw.density.sigma0(SA, CT)
f5906484["rho"] = gsw.density.rho(SA, CT, f5906484["press"])

# optional - Calculate Carbonate System with pyCO2sys instead of using LIAR vars
'''
Z = pyco2.sys(
    par1=f5906484["TA"],
    par2=f5906484["pH_tot"],
    par1_type=1,
    par2_type=3,
    salinity=f5906484["sal"],
    temperature=25,
    temperature_out=f5906484["temperature"],
    pressure=0,
    pressure_out=f5906484["press"],
    #total_phosphate=f5906484["phosphate"],
    #total_ammonia=f5906484["NH4"],
    opt_pH_scale=1,
    opt_k_carbonic=10,
    opt_total_borate=2,
    opt_k_fluoride=2,
)
f5906484["DIC"] = Z["dic"]
f5906484["pH_insitu"] = Z["pH_out"]
'''

# Remove Bad Data and Save New Vectors
T = np.array(f5906484["temperature"])[idx_flag]  # degrees C
S = np.array(f5906484["sal"])[idx_flag]  # psu
P = np.array(f5906484["press"])[idx_flag]  # dbar
rho = np.array(f5906484["rho"])[idx_flag]  # kg/m3
sigma0 = np.array(f5906484["sigma0"])[idx_flag]  # kg/m3 #
DIC = np.array(f5906484["DIC"])[idx_flag]  # umol/kg
#DIP = np.array(f5906484["phosphate"])[idx_flag] / rho * 1000  # umol/kg
NO2 = np.array(f5906484["NO2"])[idx_flag] / rho * 1000  # umol/kg
NO3 = np.array(f5906484["NO3"])[idx_flag] / rho * 1000  # umol/kg
#NH4 = np.array(f5906484["NH4"])[idx_flag] / rho * 1000  # umol/kg
#Nstar = ((NO2 + NO3) - 16 * DIP + 2.9)  # umol/kg  ## Change the DIP coefficient (either 11.4 or 16) for sensitivity analyses.
TA = np.array(f5906484["TA"])[idx_flag]  # umol/kg
pH = np.array(f5906484["pH_insitu"])[idx_flag]
O2 = np.array(f5906484["O2"])[idx_flag]  # umol/kg
station = np.array(f5906484["station"])[idx_flag]
#cast = np.array(f5906484["cast"])[idx_flag]
lat = np.array(f5906484["lat"])[idx_flag]
lon = np.array(f5906484["lon"])[idx_flag]
#pH_tot = np.array(f5906484["pH_tot"])[idx_flag]
#pH_temp = np.array(f5906484["pH_temp"])[idx_flag]
#fluor = np.array(f5906484["fluor"])[idx_flag]
f5906484_df_data = np.array(
    [
        station,
        #cast,
        lat,
        lon,
        T,
        S,
        P,
        rho,
        sigma0,
        DIC,
        #DIP,
        NO3,
        NO2,
        #NH4,
        #Nstar,
        TA,
        pH,
        #pH_tot,
        #pH_temp,
        O2,
        #fluor,
    ]
).T
f5906484_clean = pd.DataFrame(
    data=f5906484_df_data,
    columns=[
        "Station",
        #"cast",
        "lat",
        "lon",
        "T",
        "S",
        "P",
        "rho",
        "sigma0",
        "DIC",
        #"DIP",
        "NO3",
        "NO2",
        #"NH4",
        #"Nstar",
        "TA",
        "pH",
        #"pH_tot",
        #"pH_temp",
        "O2",
        #"fluor",
    ],
)

# Save the clean version
f5906484_clean.to_csv("f5906484_clean.csv")
