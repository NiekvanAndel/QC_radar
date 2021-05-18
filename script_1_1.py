'''
##| creatieve commons attribution share alike | cc-by-sa 4.0 | Lotte de Vos |##





##   This program is distributed in the hope that it will be useful,


##   but WITHOUT ANY WARRANTY; without even the implied warranty of


##   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

PWS QC-radar | Niek van Andel | based on QC by Lotte de Vos
Script 1.1: calculating High Influx flags

Input:
-Text files with timestamp and rainfall amount columns of rainfall observations by weather stations, where the timestamp indicates the end interval of the observation and the rainfall value indicates the amount of rainfall in mm since last observation. The timeseries may not have variable time intervals.
-Text files with timestamp and rainfall amount columns  of rainfall observations of the RAP radar product, in 5 min resolution at the overlying 9 pixels of the weather stations.
-BCF_radar_raw.npy: list of the first guess bias correction factors of the RAP radar product

Output:
- HI_filter.npy: list of High Influx flags
'''
print('start file 1_1')
import numpy as np
import pandas as pd

#settings
nint=6 # more than 6 consequetive intervals of 5 min need to be dry in order to call it a faulty zero measurement (>30 min)
phi_a=0.4
phi_b=10

#Read files
#files
loc_input="/usr/people/andel/Documents/QC_radar_python/Inputfiles"
loc_input2="/usr/people/andel/Documents/QC_radar_python"
#loc_input="/Users/niek/Documents/KNMI/QC_radar/Inputfiles"
#loc_input2="/Users/niek/Documents/KNMI/QC_radar"
loc_output=loc_input2
filename_radar='avg_prec_AMS.csv'
filename_sensors='data_select.csv'
DBC_radar_raw=np.load(loc_input2+'/BCF_radar_raw.npy')
DBC_sensor=1.24

df_radar = pd.read_csv(loc_input+'/'+filename_radar)
df_radar2=np.array(df_radar)[:219168,1:].astype(float)*DBC_radar_raw


df_sensors= pd.read_csv(loc_input+'/'+filename_sensors, delimiter=',')
df_sensors2=np.array(df_sensors)[:219168,1:].astype(float)*DBC_sensor


HI_array=np.ones(np.shape(df_sensors2))*-1
HI_array[np.where(np.isnan(df_radar2)==True)]=np.nan
HI_array[np.where(df_sensors2>=0)==False]=np.nan

for i in np.arange(nint,np.shape(df_sensors2)[0]):
    print(i)
    for j in np.arange(np.shape(df_sensors2)[1]):
        if (df_radar2[i,j]>=0)==False:
            HI_array[i,j]=-1
            
        elif (df_sensors2[i,j]>=0)==False :
            HI_array[i,j]=0

        elif df_radar2[i,j] < phi_a:
            if df_sensors2[i,j]>phi_b:
                HI_array[i,j]=1
            else:
                HI_array[i,j]=0
        else:
            if df_sensors2[i,j]>phi_b/phi_a*df_radar2[i,j]:
                HI_array[i,j]=1
            else:
                HI_array[i,j]=0

np.save(loc_output+'/HI_filter.npy', HI_array)
