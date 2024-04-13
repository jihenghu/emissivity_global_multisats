#%%
import sys  
import os
import cdsapi
if len(sys.argv) <= 1:
    print("ERA5-download-No sufficent args!")
    exit()
else:
    print("--------- Start download ERA5-PL-"+sys.argv[1]+"_utc_"+sys.argv[2]+"00.nc ------------------------------")

# 将标准输出重定向到空设备
sys.stdout = open('/dev/null', 'w')

yyyymmdd= sys.argv[1]
UTC= sys.argv[2]
CDIR= sys.argv[3]

if not os.path.exists(CDIR+"/"+yyyymmdd):
    os.mkdir(CDIR+"/"+yyyymmdd)

c = cdsapi.Client()
c.retrieve(
    'reanalysis-era5-pressure-levels',
    {
        'product_type': 'reanalysis',
        'format': 'netcdf',
        'variable': [
            'specific_humidity', 'temperature',
        ],
        'pressure_level': [
            '50',
            '70', '100', '125',
            '150', '175', '200',
            '225', '250', '300',
            '350', '400', '450',
            '500', '550', '600',
            '650', '700', '750',
            '775', '800', '825',
            '850', '875', '900',
            '925', '950', '975',
            '1000',
        ],
        'year': yyyymmdd[0:4],
        'month': yyyymmdd[4:6],
        'day': yyyymmdd[6:8],
        'time': UTC+':00',
        'area': [70, -180, -70, 180,],
        # 'area': [30, 120, 20, 130,],
        'grid': ['0.25','0.25']
    },
    CDIR+"/"+yyyymmdd+'/ERA5-PL-GBL-25km-'+yyyymmdd+'-'+UTC+'00.nc')
    # 'ERA5-PL-GBL-'+year+month+day+time+'.nc')

# %%
