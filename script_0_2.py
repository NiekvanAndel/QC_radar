'''
PWS QC-radar | Niek van Andel
Script 0.2: calculating first guess bias correction for RAP radar product measurements

Input:
-Text files with timestamp and rainfall amount columns of rainfall observations by weather stations, where the timestamp indicates the end interval of the observation and the rainfall value indicates the amount of rainfall in mm since last observation. The timeseries may not have variable time intervals.
-Text files with timestamp and rainfall amount columns  of rainfall observations of the RAP radar product, in 5 min resolution at the overlying 9 pixels of the weather stations.

Output:
- Comparelist.npy: list of interval lengths at each time step at which the bias correction factor is based on
- BCF_radar_raw.npy: list of bias correction factors which can be applied to the RAP radar product, so that the QC can be calculated
'''

import numpy as np
import pandas as pd

from pyproj import Proj, transform

print("start file 0.2")

#settings:
DBC_sensor=1.24#default bias correction for the sensors. For example: 1.24 for netatmo stations
DBC_radar_raw=1.42
mint=int(12*24*7*2) # number of 5 min intervals over which the correlation and bias are calculated, here 2 weeks
mrain=100 # number of 5 min intervals that have nonzero values to calculate bias and correlation on
mmatch=200 # number of 5 min intervals that should at least overlap with neighbouring station observations for the neighbourstation to be included
gamma=0.15 # exclude station when the median of correlations with neighbouring stations is lower than corthres
nstat=5

d=50000
d_use=False 

#read sensor measurements.
loc_input="/usr/people/andel/Documents/QC_radar_python/Inputfiles"
loc_input2="/usr/people/andel/Documents/QC_radar_python"
#loc_input="/Users/niek/Documents/KNMI/QC_radar/Inputfiles"
#loc_input2="/Users/niek/Documents/KNMI/QC_radar"
loc_output=loc_input2
loc_meta=loc_input

filename_sensors='data_select.csv'
lonlatlist=pd.read_csv(loc_meta+'/metadata_AMS.csv')

df_sensors= pd.read_csv(loc_input+'/'+filename_sensors, delimiter=',')
df_sensors2=np.array(df_sensors)[:219168,1:].astype(float)

#read RAP radar rainfall observations.
filename_radar_raw='avg_prec_AMS.csv'
df_radar = pd.read_csv(loc_input+'/'+filename_radar_raw)
df_radar2=np.array(df_radar)[:219168,1:].astype(float)

#apply masks at intervals where data is not present
df_sensors3=df_sensors2
df_sensors3[np.where(np.isnan(df_sensors2)==True)]=10000
#df_sensors3=np.nan_to_num(df_sensors2, nan=10000)
df_sensors3_mask=np.ma.masked_greater(df_sensors3, 9000)

df_radar3=df_radar2
df_radar3[np.where(np.isnan(df_sensors2)==True)]=10000
#df_radar3=np.nan_to_num(df_radar2, nan=10000)
df_radar3_mask=np.ma.masked_greater(df_radar3, 9000)

df_radar3_mask=np.ma.masked_where(np.ma.getmask(df_sensors3_mask),df_radar3_mask)
df_sensors3_mask=np.ma.masked_where(np.ma.getmask(df_radar3_mask),df_sensors3_mask)*DBC_sensor

# In[121]:
print('files loaded')
name=lonlatlist['id']
lonlist=lonlatlist['lon']
latlist=lonlatlist['lat']
if d_use==True:
    #calculate neighbourlist of stations

    
    #this transformes the lons and lats into rijksdriekhoekscoordinates
    m_rijks=Proj("+proj=stere +lat_0=52.1561605555555 +lon_0=5.38763888888888 +k=0.9999079 +x_0=155000 +y_0=463000 +towgs84=593.16,26.15,478.54 +ellps=bessel")
    #x_rijks,y_rijks=m_rijks(X,Y)
    
    if isinstance(lonlist, int)==True:
        lonlist=np.ones(1)*lonlist
        latlist=np.ones(1)*latlist
    
    
    xr=np.zeros(np.shape(lonlist))
    yr=np.zeros(np.shape(latlist))
    for i in np.arange(len(lonlist)):
        lon=lonlist[i]
        lat=latlist[i]
        xr[i],yr[i]=m_rijks(lon, lat)
    
    neighbourlist=[]
    q_array=np.arange(len(lonlist))
    for i in np.arange(len(lonlist)):
        dist=np.sqrt((xr-xr[i])**2+(yr-yr[i])**2)
        print(i)
        dist_mask=np.ma.masked_greater(dist,d)
        dist_mask=np.ma.masked_equal(dist_mask,0)
        q_tuple=np.ma.masked_where(np.ma.getmask(dist_mask),q_array)
    
        neighbourlist.append(q_tuple)
else:
    neighbourlist=[]
    q_array=np.arange(len(lonlist))
    for i in np.arange(len(lonlist)):
        print(i)
        q_tuple=np.ma.masked_equal(q_array,i)
        neighbourlist.append(q_tuple)
print('neighbourlist created')

# In[122]:
'''
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

np.save(loc_output+'/comparelist.npy',comparelist)
'''
print('comparelist created')

##calculate bias of the RAP radar product at sensor locations (leave one out method)
comparelist=np.load(loc_input2+'/comparelist.npy')
# In[123]:
BCF_array_radar=np.ones(np.shape(df_sensors2))*DBC_radar_raw

#q=np.arange(np.shape(df_sensors2)[1])
#comparelist=np.ma.masked_where((comparelist>=0)==False,comparelist)
for i in np.arange(np.shape(df_sensors2)[0]):
    print(i)
    biaslist=np.zeros((np.shape(df_sensors2)[1],1))*np.nan
    for j in np.arange(np.shape(df_sensors2)[1]):
        if (comparelist[i,j]>=0)==False:
            pass
        elif len(np.ma.compressed(df_radar3_mask[int(comparelist[i,j]):i,j]))<mmatch:
            pass

        else:
            radar_array=df_radar3_mask[int(comparelist[i,j]):i,j].astype(float)
            sensor_array=df_sensors3_mask[int(comparelist[i,j]):i,j].astype(float)
            biaslist[j]=np.mean(radar_array-sensor_array)/np.mean(sensor_array)

    biaslist2=np.ma.masked_where((np.abs(biaslist)>=0)==False,biaslist)
    for j in np.arange(np.shape(df_sensors2)[1]):



        #selection_tuple=tuple(q[:j])+tuple( q[j+1:])

        biaslist3=np.ma.compressed( biaslist2[neighbourlist[j].compressed(),:])

        if len(biaslist3)>nstat:
            bias_new=np.median(biaslist3)

            bias_corr=1/(1+bias_new)
            if (bias_corr>=0)==False:
                pass
            elif bias_corr>5:
                pass
            elif bias_corr<0.2:
                pass
            else:
                BCF_array_radar[i:,j]=bias_corr*np.ones(len(BCF_array_radar[i:,0]))

        else:
            pass




np.save(loc_output+'/BCF_radar_raw.npy',BCF_array_radar)
