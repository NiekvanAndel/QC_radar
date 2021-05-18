'''
PWS QC-radar | Niek van Andel
Script 1.2: calculating Station Outlier flags

Input:
-Text files with timestamp and rainfall amount columns of rainfall observations by weather stations, where the timestamp indicates the end interval of the observation and the rainfall value indicates the amount of rainfall in mm since last observation. The timeseries may not have variable time intervals.
-Text files with timestamp and rainfall amount columns  of rainfall observations of the RAP radar product, in 5 min resolution at the overlying 9 pixels of the weather stations.
-BCF_radar_raw.npy: list of the first guess bias correction factors of the RAP radar product
-comparelist.npy: list of intervals at which the bias has to be corrected.
-FZ_filter & HI_filter.py: filters calculated in earlier scripts.

Output:
- SO_filter.npy: list of Station Outlier flags
'''
import script_0_2
import script_1_0
import script_1_1

import numpy as np
import pandas as pd
from scipy import stats

#settings
mint=int(12*24*7*2) # number of 5 min intervals over which the correlation and bias are calculated, here 2 weeks
mrain=100 # number of 5 min intervals that have nonzero values to calculate bias and correlation on
mmatch=200 # number of 5 min intervals that should at least overlap with neighbouring station observations for the neighbourstation to be included
gamma=0.15 # exclude station when the median of correlations with neighbouring stations is lower than corthres
nstat=5
beta=0.2


#load files
loc_input="/usr/people/andel/Documents/QC_radar_python/Inputfiles"
loc_input2="/usr/people/andel/Documents/QC_radar_python"
#loc_input="/Users/niek/Documents/KNMI/QC_radar/Inputfiles"
#loc_input2="/Users/niek/Documents/KNMI/QC_radar"
loc_output=loc_input2
filename_radar='avg_prec_AMS.csv'
filename_sensors='data_select.csv'
DBC_sensor=1.24

DBC_radar_raw=np.load(loc_input2+'/BCF_radar_raw.npy')
FZ_filter=np.load(loc_input2+'/FZ_filter.npy')
HI_filter=np.load(loc_input2+'/HI_filter.npy')
#comparelist=np.load('comparelist.npy')

default_bias_correction=np.ones(np.shape(DBC_radar_raw))*DBC_sensor

df_radar = pd.read_csv(loc_input+'/'+filename_radar)
df_radar2a=np.array(df_radar)[:219168,1:].astype(float)
df_radar2=df_radar2a*DBC_radar_raw

df_sensors= pd.read_csv(loc_input+'/'+filename_sensors, delimiter=',')
df_sensors2=np.array(df_sensors)[:219168,1:].astype(float)*DBC_sensor



#create files to calculate with
df_radar2[np.where(HI_filter==1)]=np.nan
df_radar2[np.where(FZ_filter==1)]=np.nan

df_sensors2[np.where(HI_filter==1)]=np.nan
df_sensors2[np.where(FZ_filter==1)]=np.nan

#apply masks at intervals where data is not present
df_sensors3=df_sensors2
df_sensors3[np.where(np.isnan(df_sensors2)==True)]=10000
#df_sensors3=np.nan_to_num(df_sensors2, nan=10000)
df_sensors3_mask=np.ma.masked_greater(df_sensors3, 9000)

df_radar3=df_radar2
df_radar3[np.where(np.isnan(df_radar2)==True)]=10000
#df_radar3=np.nan_to_num(df_radar2, nan=10000)
df_radar3_mask=np.ma.masked_greater(df_radar3, 9000)

df_radar3_mask=np.ma.masked_where(np.ma.getmask(df_sensors3_mask),df_radar3_mask)
df_sensors3_mask=np.ma.masked_where(np.ma.getmask(df_radar3_mask),df_sensors3_mask)

##calculate intervals of comparison: if comparelist.npy has not been calculated, calculate it here. If not: quote this part.
Nintrain=np.zeros(np.shape(df_sensors2))
Nintrain[np.where(df_sensors2>0)]=1
Nintraincum=np.cumsum(Nintrain, axis=0)

compareA=np.zeros(np.shape(df_sensors2))
for i in np.arange(np.shape(df_sensors2)[0]):
    print(i)
    for j in np.arange(np.shape(df_sensors2)[1]):
        a=np.where((Nintraincum[i,j]-mrain+1)==Nintraincum[:,j])
        if np.shape(a)[1]>0:
            compareA[i,j]=a[0][-1]
        else:
            compareA[i,j]=np.nan

compareB=np.zeros(np.shape(df_sensors2))
compareB[:mint+1,:]=np.nan
for i in np.arange(np.shape(df_sensors2)[1]):
    compareB[mint:,i]=np.arange(np.shape(df_sensors2)[0]-mint)

comparelist=np.zeros(np.shape(df_sensors2))
for i in np.arange(np.shape(df_sensors2)[0]):
    for j in np.arange(np.shape(df_sensors2)[1]):
        if compareB[i,j]==np.nan:
            comparelist[i,j]=np.nan
        elif compareA[i,j]<compareB[i,j]:
            comparelist[i,j]=compareA[i,j]
        else:
            comparelist[i,j]=compareB[i,j]

np.save(loc_output+'/comparelist2.npy',comparelist)

comparelist=np.load(loc_input2+'/comparelist2.npy')


SO_array=np.zeros(np.shape(df_sensors2))



correlation_array=np.zeros(np.shape(df_sensors2))
BCF_array2=np.ones(np.shape(df_sensors2))*DBC_sensor



#Calculate filter
for i in np.arange(200000, np.shape(df_sensors2)[0]):
    print(i)

    for j in np.arange(np.shape(df_sensors2)[1]):
            if (comparelist[i,j]>=0)==False:
                SO_array[i,j]=-1
            elif len(np.ma.compressed(df_radar3_mask[int(comparelist[i,j]):i,j]))<mmatch:
                #print('niet genoeg data')
                SO_array[i,j]=-1

            else:
                radar_array=df_radar3_mask[int(comparelist[i,j]):i,j].astype(float)
                sensor_array=df_sensors3_mask[int(comparelist[i,j]):i,j].astype(float)/DBC_sensor

                corr=stats.pearsonr(np.ma.compressed(radar_array),np.ma.compressed(sensor_array))[0]

                if np.isnan(corr)==True:
                    correlation_array[i,j]=10000
                else:
                    correlation_array[i,j]=corr


                if correlation_array[i,j]<gamma:

                    SO_array[i,j]=1

                elif correlation_array[i,j]>1:

                    SO_array[i,j]=-1
                else:
                    bias_new=np.mean(sensor_array-radar_array)/np.mean(radar_array) #bias of the original raw sensorvalues to the (new) radar
                    SO_array[i,j]=0

                    bias_corr=1/(1+bias_new)

                    if (bias_corr>=0)==False: #so if bias_corr is negative or nan
                        SO_array[i,j]=-1

                    elif bias_corr==np.inf:
                        SO_array[i,j]=-1

                    else:


                        if np.abs(np.log(bias_corr/BCF_array2[i,j]))>np.log(1+beta):

                            BCF_array2[i:,j]=bias_corr*np.ones(len(BCF_array2[i:,0]))


np.save(loc_output+'/SO_filter.npy', SO_array)
np.save(loc_output+'/BCF_array.npy', BCF_array2) #only output: this is the bias of the raw sensors(even without default bias correction)
np.save(loc_output+'/correlation_array.npy',correlation_array)
