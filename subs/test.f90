     
! include 'sub_read_goes_l2c_vars.f90'
include 'sub_tools.f90'

real :: lon1, lat1
real :: x, y
character*255 ::  lon, lat
	LOGICAL :: FLAG
	
 call getarg(1,lon)
 call getarg(2,lat)
 
 read(lon,*) lon1
 read(lat,*) lat1
 
 call index_GOES_R_LonLat(lon1,lat1,x,y,FLAG)

 print*,'[x,y]',x,y, FLAG
 
 end