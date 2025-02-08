# Microwave Emissivity Retrieved from Multiple Sensors Over GPM Era 

![License](https://img.shields.io/badge/license-MIT-blue.svg) <!-- Replace with your project's license -->
![Version](https://img.shields.io/badge/version-1.0-green.svg) <!-- Replace with your project's version -->

This repository contains the basic code for global clear-sky land emissivity retrieval using combined observations from multiple passive microwave sensors operating in GPM era, including GPM-GMI, Fengyun3B/C/D-MWRI, AMSR2. The project uses cloud mask data from geostationary platforms.

---

## Features
- **Instrument-Specific:** brightness temperatures from GPM-GMI, Fengyun3B/C/D-MWRI, AMSR2.
- **Multi-geostationary cloud:** Includes data from multiple geostationary satellites such as Himawari-8/9, MSG-1/2/3, and GOES-16.
- **Use ERA5 and ERA5-land analysis**.
- **Land-Only:** The retrieval is specifically designed for land surfaces.
- **Clear-Sky Conditions:** The algorithm is optimized for clear-sky conditions.
- **Stable Batch Execution:** The code is designed for stable and efficient batch processing.
- **Improved I/O:** Supports both ASCII and HDF5 formats for data input and output.
  
## Requirements
To run the code in this repository, the following dependencies are required:

- **Python3.9** and python3-pip
- **wget**：for automatically downloading the Himawari-8/9 cloud data.
- **HDF5**: For handling HDF5 file formats, you may refer to [Compiling WRF](https://www2.mmm.ucar.edu/wrf/OnLineTutorial/compilation_tutorial.php)  or my note at [RTTOV installation: dependency](https://jihenghu.github.io/research/rttov/rttov132-installlibs/)
- **NETCDF4**: For handling NetCDF file formats, refer to [Compiling WRF](https://www2.mmm.ucar.edu/wrf/OnLineTutorial/compilation_tutorial.php) or my note [RTTOV installation: dependency](https://jihenghu.github.io/research/rttov/rttov132-installlibs/).
- **RTTOV V13.2**: Radiative Transfer for TOVS (RTTOV) version 13.2. Refer to the [RTTOV documentation](https://nwp-saf.eumetsat.int/site/software/rttov/rttov-v13/) or [my installation note](https://jihenghu.github.io/research/rttov/rttov132-install/).
- **CDSAPI**: Python library for accessing the Copernicus Climate Data Store (CDS), used to retrieve ERA5 analysis, refer to the [CDS API documentation](https://cds.climate.copernicus.eu/how-to-api) to set up your condentials and the api tool.
- **EUMDAC**: Python library for accessing EUMETSAT Data Centre used to retrieve MASG cloud ata archive, refer to [eumetsat-data-access-client-eumdac-guide](https://user.eumetsat.int/resources/user-guides/eumetsat-data-access-client-eumdac-guide) to setup your account and api service.
- **Amazon AWS Service**: For accessing GOES cloud data stored on Amazon Web Services. Refer to the [AWS CLI documentation](https://docs.aws.amazon.com/cli/latest/userguide/getting-started-install.html) for installation and setup.
- **eccodes** : For handling the MSG grib2 dat, see [ecCodes installation](https://confluence.ecmwf.int/display/ECC/ecCodes+installation).

## Data archive access
- **GPM_L1C** : go to the Earth Data to creat cookies `~/.urs_cookies`, see guidance at [GPM_L1C datapool](https://gpm1.gesdisc.eosdis.nasa.gov/data/GPM_L1C/).
- **JMA HIMAWARI portal**: register at [P-Tree System](https://www.eorc.jaxa.jp/ptree/registration_top.html), you need a ftpname and passwd for the system. Replace your usename and passwd in `subs/sub_download_ahi_l2c.f90` like: 'ftp://{usename}:{passwd}@ftp.ptree.jaxa.jp/pub/himawari/L2/CLP/010/' to ensure personal valid access.
- **NOAA GOES data**: We directly extract data from [AWS S3 bucket](https://noaa-goes16.s3.amazonaws.com/index.html) for date after 20191231. If you need data before this date, you need to prepare it manually by ordering at [NOAA CLASS system](https://www.aev.class.noaa.gov/saa/products/welcome).
- **EUMSAT MSG data**: is accessed with condentials, setup it following [eumetsat-data-access-client-eumdac-guide](https://user.eumetsat.int/resources/user-guides/eumetsat-data-access-client-eumdac-guide).
- **ERA5 and ERA5-land**:  are accessed with [CDSAPI](https://cds.climate.copernicus.eu/how-to-api), the portal version we embled is currently outdated, which needs further adaption on the 'subs/read_era5*.f90' to read the lasted format.
- **Fengyun3 data**: we currently do not support batch downloading of CMA data, so you need to prepare it in advance in case of need.

## Usage
- install and set up your project locally.

```bash
# Clone the repository
git clone https://github.com/jihenghu/emissivity_global_multisats.git

# cd to the project directory
cd emissivity_global_multisats
```

Change the branch [branch-name] to work on target sensors: `gmi-dev`, `amsr2-dev`, `fy3b-dev`, `fy3c-dev`, or `fy3d-dev`.

```bash
git checkout [branch-name]
```

- Configuration

modify the `make.sh` to specify the library direction and configuration options:

```sh make.sh
gfortran -mcmodel=large -ffpe-summary=none -I/home/jihenghu/netcdf/include -I/home/jihenghu/hdf5/include -I/home/jihenghu/eccodes/include -L/home/jihenghu/netcdf/lib -lnetcdff -lnetcdf  -L/home/jihenghu/hdf5/lib -lhdf5_fortran -L/home/jihenghu/eccodes/lib64 -leccodes_f90 -leccodes  -c main_clear_retrieve_landonly.f90 -o ./main_clear_retrieve_landonly.o
gfortran -I/home/jihenghu/rttov13/mod -I/home/jihenghu/rttov13/include -fPIC -O3 -fopenmp -ffree-line-length-none  -c rttov_retrieve_gmi_emiss_clearsky.f90 -o ./rttov_retrieve_gmi_emiss_clearsky.o
gfortran -o ./main_gmi_em.exe ./main_clear_retrieve_landonly.o ./rttov_retrieve_gmi_emiss_clearsky.o \
-L/home/jihenghu/rttov13/lib -lrttov13_brdf_atlas -lrttov13_emis_atlas -lrttov13_mw_scatt -lrttov13_other -lrttov13_coef_io -lrttov13_hdf -lrttov13_parallel -lrttov13_main  \
-L/home/jihenghu/netcdf/lib -lnetcdff -L/home/jihenghu/hdf5/lib -lhdf5hl_fortran -lhdf5_hl -lhdf5_fortran -lhdf5 -lz -fopenmp -L/home/jihenghu/eccodes/lib64 -leccodes_f90 -leccodes
```
- Execution
```bash
sh make.sh
```
generates an executable file named like `./main_gmi_em.exe`, excute it with a argument to specify the date:

```bash
./main_gmi_em.exe 20150101
```

- Output
The retrieved emissivity data will be saved orbit by orbit in the output directory specified in main `main_clear_retrieve_landonly.f90`:
```f90 main_clear_retrieve_landonly.f90
EMISS_OUTDIR = '/home/jihenghu/data00/data_em/GMI_EMISSIVITY/'
```

## License
This project is licensed under the MIT License. See the LICENSE file for details.
