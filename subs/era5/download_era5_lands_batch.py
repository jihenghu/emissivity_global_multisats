#%%
import h5py
import sys  
import os
import cdsapi
import netCDF4 as nc

if len(sys.argv) <= 1:
    print("No sufficent args!")
    exit()
else:
    print("--------- Start download ERA5-land-"+sys.argv[1]+"_batch.nc ------------------------------")

# 将标准输出重定向到空设备
# sys.stdout = open('/dev/null', 'w')

yyyymmdd= sys.argv[1]
CDIR= sys.argv[2]

if not os.path.exists(CDIR+"/"+yyyymmdd):
    os.mkdir(CDIR+"/"+yyyymmdd)

time= [     '00:00', '01:00', '02:00',
            '03:00', '04:00', '05:00',
            '06:00', '07:00', '08:00',
            '09:00', '10:00', '11:00',
            '12:00', '13:00', '14:00',
            '15:00', '16:00', '17:00',
            '18:00', '19:00', '20:00',
            '21:00', '22:00', '23:00',
        ]

file_path = CDIR+"/"+yyyymmdd+'/ERA5-Land-GBL-25km-'+yyyymmdd+'_batch.nc'
if os.path.exists(file_path):
    print(f"File {file_path} already exists, skipping...")
    exit()
else:
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
            'time': [
                '00:00', '01:00', '02:00',
                '03:00', '04:00', '05:00',
                '06:00', '07:00', '08:00',
                '09:00', '10:00', '11:00',
                '12:00', '13:00', '14:00',
                '15:00', '16:00', '17:00',
                '18:00', '19:00', '20:00',
                '21:00', '22:00', '23:00',
            ],
            'area': [70, -180, -70, 180,],
            # 'area': [30, 120, 20, 130,],
            'grid': ['0.25','0.25']
        },
        file_path)
        # 'ERA5-PL-GBL-'+year+month+day+time+'.nc')

# Load the original NetCDF file
data = nc.Dataset(file_path, 'r')

# Extract data dimensions
longitude = data.variables['longitude'][:]
latitude = data.variables['latitude'][:]
time = data.variables['time'][:]
t2m = data.variables['t2m'][:]
skt = data.variables['skt'][:]
snowc = data.variables['snowc'][:]
sp = data.variables['sp'][:]
swvl1 = data.variables['swvl1'][:]

# Function to create a new NetCDF file for each time step
def create_nc_file(time_index):
    filename=CDIR+"/"+yyyymmdd+'/ERA5-Land-GBL-25km-'+yyyymmdd+f'-{i:02d}00.nc'

    # Check if the file already exists
    if os.path.exists(filename):
        print(f"File {filename} already exists, skipping...")
        return
    
    with nc.Dataset(filename, 'w', format='NETCDF3_CLASSIC') as new_nc:
        # Create dimensions
        new_nc.createDimension('longitude', len(longitude))
        new_nc.createDimension('latitude', len(latitude))
        new_nc.createDimension('time', 1)
        
        # Create variables
        lon = new_nc.createVariable('longitude', 'f4', ('longitude',))
        lat = new_nc.createVariable('latitude', 'f4', ('latitude',))
        time_var = new_nc.createVariable('time', 'i4', ('time',))
        t2m_var = new_nc.createVariable('t2m', 'i2', ('time', 'latitude', 'longitude'))
        skt_var = new_nc.createVariable('skt', 'i2', ('time', 'latitude', 'longitude'))
        snowc_var = new_nc.createVariable('snowc', 'i2', ('time', 'latitude', 'longitude'))
        sp_var = new_nc.createVariable('sp', 'i2', ('time', 'latitude', 'longitude'))
        swvl1_var = new_nc.createVariable('swvl1', 'i2', ('time', 'latitude', 'longitude'))
        
        # Copy variable attributes
        lon.setncatts(data.variables['longitude'].__dict__)
        lat.setncatts(data.variables['latitude'].__dict__)
        time_var.setncatts(data.variables['time'].__dict__)
        t2m_var.setncatts(data.variables['t2m'].__dict__)
        skt_var.setncatts(data.variables['skt'].__dict__)
        snowc_var.setncatts(data.variables['snowc'].__dict__)
        sp_var.setncatts(data.variables['sp'].__dict__)
        swvl1_var.setncatts(data.variables['swvl1'].__dict__)
        
        # Write data
        lon[:] = longitude
        lat[:] = latitude
        time_var[:] = time[time_index:time_index+1]
        t2m_var[0, :, :] = t2m[time_index, :, :]
        skt_var[0, :, :] = skt[time_index, :, :]
        snowc_var[0, :, :] = snowc[time_index, :, :]
        sp_var[0, :, :] = sp[time_index, :, :]
        swvl1_var[0, :, :] = swvl1[time_index, :, :]

# Split the data into 24 files
for i in range(24):
    create_nc_file(i)

data.close()
os.remove(file_path)
os.system("touch "+file_path)

# %%
