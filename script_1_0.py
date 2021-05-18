'''
##| creatieve commons attribution share alike | cc-by-sa 4.0 | Lotte de Vos |##





##   This program is distributed in the hope that it will be useful,


##   but WITHOUT ANY WARRANTY; without even the implied warranty of


##   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

PWS QC-radar | Niek van Andel | based on QC by Lotte de Vos
Script 1.0: calculating Faulty Zero flags

Input:
-Text files with timestamp and rainfall amount columns of rainfall observations by weather stations, where the timestamp indicates the end interval of the observation and the rainfall value indicates the amount of rainfall in mm since last observation. The timeseries may not have variable time intervals.
-Text files with timestamp and rainfall amount columns  of rainfall observations of the RAP radar product, in 5 min resolution at the overlying 9 pixels of the weather stations.
-BCF_radar_raw.npy: list of the first guess bias correction factors of the RAP radar product

Output:
- FZ_filter.npy: list of Faulty Zero flags
'''
print('start file 1_0')
import numpy as np
import pandas as pd

#settings:
nint=6 # more than 6 consequetive intervals of 5 min need to be dry in order to call it a faulty zero measurement (>30 min)

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

df_radar3=df_radar2
df_radar3[np.where(np.isnan(df_radar2)==True)]=10000
df_radar3_mask=np.ma.masked_greater(df_radar3, 9000)

df_sensors= pd.read_csv(loc_input+'/'+filename_sensors, delimiter=',')
df_sensors2=np.array(df_sensors)[:219168,1:].astype(float)

df_sensors3=df_sensors2
df_sensors3[np.where(np.isnan(df_sensors2)==True)]=10000
df_sensors3_mask=np.ma.masked_greater(df_sensors3, 9000)*DBC_sensor

##df_radar3_mask=np.ma.masked_where(np.ma.getmask(df_sensors3_mask),df_radar3_mask)
#df_sensors3_mask=np.ma.masked_where(np.ma.getmask(df_radar3_mask),df_sensors3_mask)


#faulty zero calculation
Ref_array=np.zeros(np.shape(df_sensors2))
Sensor_array=np.zeros(np.shape(df_sensors2))
FZ_array=np.ones(np.shape(df_sensors2))*-1

####fill Ref_array:
#Ref_array[np.where(np.isnan(df_radar2)==True)]=np.nan
Ref_array[np.where(df_radar2>0)]=1
Ref_array=np.ma.masked_where(np.ma.getmask(df_radar3_mask),Ref_array)

### fill Sensor_array
Sensor_array[np.where(df_sensors3_mask>0)]=1
Sensor_array[np.where(df_sensors3_mask==0)]=0
Sensor_array=np.ma.masked_where(np.ma.getmask(df_sensors3_mask),Sensor_array)



for i in np.arange(nint,np.shape(df_sensors2)[0]):
    print(i)
    for j in np.arange(np.shape(df_sensors2)[1]):
        if len(np.ma.compressed(Ref_array[i,j]))==0: #als de vergelijkende waarde niet bestaat:
            FZ_array[i,j]=-1
        
        elif len(np.ma.compressed(Sensor_array[i,j]))==0: #als de meting een NAN is, dan vorige meting overnemen
            FZ_array[i,j]=FZ_array[i-1,j]
        else: #als er genoeg goede metingen zijn:
            if Sensor_array[i,j]>0: #als er regen valt, geen FZ vlag
                FZ_array[i,j]=0
            else:
                if FZ_array[i-1,j]==1: #als er geen regen valt en vorige meting ook al FZ, dan wederom afvlaggen
                    FZ_array[i,j]=1
                elif np.sum(Sensor_array[i-nint:i+1,j]) >0: #als er geen regen valt, maar ergens in het voorgaande wel: neit afvlaggen

                    FZ_array[i,j]=0

                else: #als er een serie nullen is van nint lengte
                    if np.sum(Ref_array[i-nint:i+1,j]) <nint+1: #als de radar ook minstens af en toe nul of aangeeft: dan is alles prima, of is bepaald zonder nstat stations
                        FZ_array[i,j]=0
                    else: #dus als de sensor nul aangeeft, maar de radar niet nul (continu) in die periode
                        FZ_array[i,j]=1

np.save(loc_output+'/FZ_filter.npy', FZ_array)
