# emissivity_global_multisatss
Global clear-sky land emissivity retrieval using combined observations from GMI, FYs, Geostationary cloud
Basic code for global emissivity retrieval, features:
- pmw: FY3B
- cloud: H-8/9, MSG-1/2/3, GOES-16
- land only
- clear sky

# Features:
- imporved I/O, including ascii and H5;
- more GeoSats included;
- essential issures mended;
- more stable batch excuated;
- specified forFY instrument, final version;

# Requirements
- HDF5, NETCDF
- RTTOV V13.2. refer to http://home.ustc.edu.cn/~hjh18305/space/research/rttov/rttov132-column/
- cdsapi refer to https://cds.climate.copernicus.eu/api-how-to 
- eumdac refer to http://home.ustc.edu.cn/~hjh18305/space/data-n-method/MSG-SERV-download-guide/
- Amazon AWS serive refer to https://docs.aws.amazon.com/zh_cn/cli/latest/userguide/getting-started-install.html
 
