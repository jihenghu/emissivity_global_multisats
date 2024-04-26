
! program main_program
!     use eccodes
!   integer, dimension(3712,3712) :: cloud_mask

!   character*255 filename

  
!   filename="/data/jihenghu/data/MSG_IODC_CLM/20200101/MSG1-SEVI-MSGCLMK-0100-0100-20200101011500.000000000Z-NA.grb"
!   call read_MSG_CLM(filename,cloud_mask)

! 	PRINT*, Minval(cloud_mask),MAxval(cloud_mask)
! 	PRINT*, cloud_mask(2461,1446:1457)  
!! ===================

	! real :: lat, lon, Lambda0=41.5
	! integer::  iy, ix

	! lon= -15
	! lat=60
	! CALL geocoord2pixcoord(lat,lon, Lambda0, ix, iy)
	! print*,  lon, lat, ix, iy
	! lon= -15
	! lat=-60
	! CALL geocoord2pixcoord(lat,lon, Lambda0, ix, iy)
	! print*,  lon, lat, ix, iy

	! lon= 90
	! lat=60
	! CALL geocoord2pixcoord(lat,lon, Lambda0, ix, iy)
	! print*,  lon, lat, ix, iy

	! lon= 90
	! lat=-60
	! CALL geocoord2pixcoord(lat,lon, Lambda0, ix, iy)
	! print*,  lon, lat, ix, iy
	! lon= 41.5
	! lat= 0
	! CALL geocoord2pixcoord(lat,lon, Lambda0, ix, iy)
	! print*,  lon, lat, ix, iy  


! end program main_program


!! See examples of grib2 read using eccodes:
!! https://confluence.ecmwf.int/display/ECC/grib_print_data
!! ALSO, use grib_dump to view the meta
subroutine read_MSG_CLM(filename,cloud_mask,status)
  use eccodes
  implicit none

  character(len=*), intent(in) :: filename
  integer, dimension(3712,3712),INTENT(OUT) :: cloud_mask
  integer,INTENT(OUT) :: status
  
  real*4, dimension(:), allocatable :: clm_1d
  integer :: ifile, igrib,is,ip

  allocate(clm_1d(3712*3712)) !! radicules
  
  call codes_open_file(ifile, trim(filename),'r',status)
  IF (status/=CODES_SUCCESS) RETURN
  
  call codes_grib_new_from_file(ifile,igrib)
  call codes_get(igrib,'values',clm_1d) 
  call codes_release(igrib)	
  call codes_close_file(ifile)

  DO is=1,3712
	DO ip=1,3712
	  cloud_mask(is,ip)=int(clm_1d((ip-1)*3712+3712-is+1))
	end do
  end do

end subroutine read_MSG_CLM




! Navigating from Geodetic Latitude (φ ) and Longitude (λ ) to index, x,y
! Adapted from https://navigator.eumetsat.int/product/EO:EUM:SW:MSG:159 , the F90 program
SUBROUTINE geocoord2pixcoord( latitude,  longitude, lon0 ,column, row,FLAG)
	IMPLICIT NONE

	REAL, INTENT (IN) ::  lon0        ! Longitude of Sub-Satellite Point in degree
	REAL, INTENT (IN) :: latitude, longitude
	INTEGER, INTENT (OUT) :: column, row

	! MODULE GLOBAL
	REAL*8, PARAMETER :: PI=3.14159265359d0

	REAL*8, PARAMETER ::  SAT_HEIGHT= 42164.0d0  ! distance from Earth centre to satellite    
	REAL*8, PARAMETER ::  R_EQ = 6378.169d0      ! radius from Earth centre to equator
	REAL*8, PARAMETER ::  R_POL= 6356.5838d0     !


	! next parameters defined in REAL, for consistency with new introduced HRV values
	REAL*8, PARAMETER ::  CFAC_NONHRV = -781648343.0d0 
	REAL*8, PARAMETER ::  LFAC_NONHRV = -781648343.0d0

	! next parameters defined in REAL, since too big for INTEGER
	! REAL*8, PARAMETER ::  CFAC_HRV = -2344945030.0d0
	! REAL*8, PARAMETER ::  LFAC_HRV = -2344945030.0d0


	INTEGER, PARAMETER :: COFF_NONHRV = 1856 
	INTEGER, PARAMETER :: LOFF_NONHRV = 1856

	! INTEGER, PARAMETER ::  COFF_HRV = 5566
	! INTEGER, PARAMETER ::  LOFF_HRV = 5566
	! END MODULE GLOBAL

	INTEGER :: c=0, l=0
	INTEGER :: ccc=0, lll=0
	INTEGER :: x=0, y=0
  
	REAL*8 :: lati, longi, SUB_LON
	REAL*8 :: c_lat
	REAL*8 :: lat
	REAL*8 :: lon
	REAL*8 :: r1, r2, r3, rn, re, rl
	REAL*8 :: xx, yy, sa
	REAL*8 :: cc, ll
	REAL*8 :: dotprod
	LOGICAL :: FLAG
	FLAG=.TRUE.

	SUB_LON=lon0*PI / 180.0d0

	lati= latitude
	longi= longitude
  
	! check if the values are sane, otherwise return error value
	if (lati .LT. -90.0d0 .OR. lati .GT. 90.0d0 .OR. longi .LT. -180.0d0 .OR. longi .GT. 180.0d0 ) then
	   row = 0
	   column = 0
	   FLAG=.FALSE.
	   return
	end if
  
	! convert them to radians 
	lat = lati*PI / 180.0d0
	lon = longi *PI / 180.0d0
  
  
	! calculate the geocentric latitude from the       
	! geographic one using equations on page 24, Ref. [1] 
  
	c_lat = atan ( (0.993243d0*(sin(lat)/cos(lat)) ))
		
  
	! using c_lat calculate the length from the Earth 
	! centre to the surface of the Earth ellipsoid    
	! equations on page 23, Ref [1]                      
	
	re = R_POL / sqrt( (1.0d0 - 0.00675701d0 * cos(c_lat) * cos(c_lat) ) )
  
  
	! calculate the forward projection using equations on page 24, Ref. [1]
  
	rl = re
	r1 = SAT_HEIGHT - rl * cos(c_lat) * cos(lon - SUB_LON)
	r2 = - rl *  cos(c_lat) * sin(lon - SUB_LON)
	r3 = rl * sin(c_lat)
	rn = sqrt( r1*r1 + r2*r2 +r3*r3 )
  
	! check for visibility, whether the point on the Earth given by the
	! latitude/longitude pair is visible from the satellte or not. This 
	! is given by the dot product between the vectors of:
	! 1) the point to the spacecraft,
	! 2) the point to the centre of the Earth.
	! If the dot product is positive the point is visible otherwise it
	! is invisible.
  
	dotprod = r1*(rl * cos(c_lat) * cos(lon - SUB_LON)) - r2*r2 - r3*r3*((r_EQ/R_POL)**2)
  
	if (dotprod .LE. 0.0d0 ) then
	   column = -999
	   row = -999
	   FLAG=.FALSE.
	   return
	end if
  
	! the forward projection is x and y 
  
	xx = atan( (-r2/r1) )
	yy = asin( (-r3/rn) )
  
  
	! convert to pixel column and row using the scaling functions on 
	! page 28, Ref. [1]. And finding nearest integer value for them. 
  
	cc = DBLE(COFF_NONHRV) + xx *  2.0d0**(-16.0d0) * DBLE(CFAC_NONHRV)
	ll = DBLE(LOFF_NONHRV) + yy *  2.0d0**(-16.0d0) * DBLE(LFAC_NONHRV)
  
  
	ccc=nint(cc)
	lll=nint(ll)		
  
	column = ccc
	row = lll
  
  END SUBROUTINE geocoord2pixcoord
  
  
