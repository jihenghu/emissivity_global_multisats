#%%
import sys  
import os
import cdsapi
if len(sys.argv) <= 1:
    print("No sufficent args!")
    exit()
else:
    print("--------- Start download ERA5-land-"+sys.argv[1]+"_utc_"+sys.argv[2]+"00.nc ------------------------------")

# 将标准输出重定向到空设备
sys.stdout = open('/dev/null', 'w')

yyyymmdd= sys.argv[1]
UTC= sys.argv[2]
CDIR= sys.argv[3]

if not os.path.exists(CDIR+"/"+yyyymmdd):
    os.mkdir(CDIR+"/"+yyyymmdd)

c = cdsapi.Client()
c.retrieve(
    'reanalysis-era5-land',
    {
        'variable': [
            '2m_temperature', 'skin_temperature', 'snow_cover',#,'10m_u_component_of_wind', '10m_v_component_of_wind'
             'surface_pressure','volumetric_soil_water_layer_1', #'snow_depth', 'soil_temperature_level_1',
            ],
        'format': 'netcdf',
        'year': yyyymmdd[0:4],
        'month': yyyymmdd[4:6],
        'day': yyyymmdd[6:8],
        'time': UTC+':00',
        'area': [70, -180, -70, 180,],
        # 'area': [30, 120, 20, 130,],
        'grid': ['0.25','0.25']
    },
    CDIR+"/"+yyyymmdd+'/ERA5-Land-GBL-25km-'+yyyymmdd+'-'+UTC+'00.nc')
    # 'ERA5-PL-GBL-'+year+month+day+time+'.nc')

# %%
