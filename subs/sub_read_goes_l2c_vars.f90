! program test_read_goes
	! Character(255) :: GEOS_L2_DIR
	! Character(8):: yyyymmdd
	! Character(4) :: HHMM
	! INTEGER, DIMENSION(5424,5424) :: CLM_GOES,PHASE_GOES	
	! REAL, DIMENSION(5424,5424) :: CER_GOES
	! REAL, DIMENSION(2712,2712) ::  COT_GOES	
	! REAL, DIMENSION(1086,1086) ::  CTH_GOES	

	! LOGICAL :: A1,A2,A3,A4,A5
	
	! yyyymmdd='20200101'
	! HHMM='0100'
	! GEOS_L2_DIR='/data/jihenghu/data/GOES_R/'
	
	
    !! CALL read_GOES_CLM(yyyymmdd,HHMM,GEOS_L2_DIR,CLM_GOES,A1)
	!! IF(.NOT.A1) PRINT*,"READ GOES file error, skip ......"
	!! PRINT*,MinVal(CLM_GOES),MAXVAL(CLM_GOES)
	!! print*,CLM_GOES(:20,1)
	
	! CALL read_GOES_PHASE(yyyymmdd,HHMM,GEOS_L2_DIR,PHASE_GOES,A2)
	! IF(.NOT.A2) PRINT*,"READ GOES file error, skip ......"
	! PRINT*,MinVal(PHASE_GOES),MAXVAL(PHASE_GOES)
	 ! print*,PHASE_GOES(:20,1)
	
	! CALL read_GOES_CER(yyyymmdd,HHMM,GEOS_L2_DIR,CER_GOES,A3)  	!! 0.00152602  
	! IF(.NOT.A3) PRINT*,"READ GOES file error, skip ......"
	! PRINT*,MinVal(CER_GOES),MAXVAL(CER_GOES)
	! print*,CER_GOES(:20,1)	
	
	! CALL read_GOES_COT(yyyymmdd,HHMM,GEOS_L2_DIR,COT_GOES,A4)  
	! IF(.NOT.A4) PRINT*,"READ GOES file error, skip ......"
	! PRINT*,MinVal(COT_GOES),MAXVAL(COT_GOES)
	! print*,COT_GOES(1202,885:894)	
	

    ! CALL read_GOES_CTH(yyyymmdd,HHMM,GEOS_L2_DIR,CTH_GOES,A5) 	
	! IF(.NOT.A5) PRINT*,"READ GOES file error, skip ......"
	! PRINT*,MinVal(CTH_GOES),MAXVAL(CTH_GOES)
	! print*,CTH_GOES(330,510:521)	
! end program test_read_goes

! program test_read_goes	
    ! LOGICAL :: FLAG
	! REAL lambda,phi,x,y
	! lambda=0
	! phi=30
	! call index_GOES_R_LonLat(lambda,phi,x,y,FLAG)
	! print*,x,y
	! lambda=-130
		! phi=-30
	! call index_GOES_R_LonLat(lambda,phi,x,y,FLAG)
	! print*,x,y	
	
	! x=-0.2
	! y=0.095340 
	! call lonlat_GOES_R_xy(x,y,lambda,phi,FLAG)
	! print*,lambda,phi
	
! end program test_read_goes


!! &copy; 2023.12  Jiheng Hu @ University of Michigan, Ann Arbor
!! To read Cloud Properties from GOES-R Files

!! clear sky mask
SUBROUTINE read_GOES_CLM(yyyymmdd,HHMM,GEOS_L2_DIR,CLM_GOES,AG)
	INCLUDE 'netcdf.inc'
	Character(*), INTENT(IN) :: GEOS_L2_DIR
	Character(8), INTENT(IN) :: yyyymmdd
	Character(4), INTENT(IN) :: HHMM
	INTEGER, DIMENSION(5424,5424), INTENT(OUT) ::  CLM_GOES	
	LOGICAL, INTENT(OUT) :: AG
	Integer status,ncid,varid
	! Character*255 :: filename
	AG=.TRUE.	
	! filename=TRIM(GEOS_L2_DIR)//"/"//yyyymmdd(1:4)//"/"//yyyymmdd// &
	! "/OR_ABI-L2-ACMF-M6_G16_"//yyyymmdd//"_"//HHMM//".NC"
	
    status=nf_open(trim(adjustl(GEOS_L2_DIR)),nf_nowrite,ncid)	
	status=nf_inq_varid (ncid, 'BCM', varid) 
	status=nf_get_var_int(ncid,varid,CLM_GOES)
	status=nf_close(ncid)	
	
	IF(status.NE.0) AG=.FALSE.
	
END SUBROUTINE read_GOES_CLM

!! cloud top phase
SUBROUTINE read_GOES_PHASE(yyyymmdd,HHMM,GEOS_L2_DIR,PHASE_GOES,AG)
	INCLUDE 'netcdf.inc'
	Character(*), INTENT(IN) :: GEOS_L2_DIR
	Character(8), INTENT(IN) :: yyyymmdd
	Character(4), INTENT(IN) :: HHMM
	INTEGER, DIMENSION(5424,5424), INTENT(OUT) ::  PHASE_GOES	
	LOGICAL, INTENT(OUT) :: AG
	Integer status,ncid,varid
	! Character*255 :: filename	
	AG=.TRUE.
	
	! filename=TRIM(GEOS_L2_DIR)//"/"//yyyymmdd(1:4)//"/"//yyyymmdd// &
	! "/OR_ABI-L2-ACTPF-M6_G16_"//yyyymmdd//"_"//HHMM//".NC"
	
    status=nf_open(trim(adjustl(GEOS_L2_DIR)),nf_nowrite,ncid)	
	status=nf_inq_varid (ncid, 'Phase', varid) 
	status=nf_get_var_int(ncid,varid,PHASE_GOES)
	status=nf_close(ncid)	
	
	IF(status.NE.0)	AG=.FALSE.

END SUBROUTINE read_GOES_PHASE

!! cloud effective radius
SUBROUTINE read_GOES_CER(yyyymmdd,HHMM,GEOS_L2_DIR,CER_GOES,AG)
	INCLUDE 'netcdf.inc'
	Character(*), INTENT(IN) :: GEOS_L2_DIR
	Character(8), INTENT(IN) :: yyyymmdd
	Character(4), INTENT(IN) :: HHMM
	REAL, DIMENSION(5424,5424), INTENT(OUT)  ::  CER_GOES	
	Integer, DIMENSION(5424,5424) ::  CER	
	LOGICAL, INTENT(OUT)  :: AG
	Integer status,ncid,varid
	! Character*255 :: filename	
	AG=.TRUE.
	
	! filename=TRIM(GEOS_L2_DIR)//"/"//yyyymmdd(1:4)//"/"//yyyymmdd// &
	! "/OR_ABI-L2-CPSF-M6_G16_"//yyyymmdd//"_"//HHMM//".NC"
	
    status=nf_open(trim(adjustl(GEOS_L2_DIR)),nf_nowrite,ncid)	
	status=nf_inq_varid (ncid, 'PSD', varid) 
	status=nf_get_var_int(ncid,varid,CER)
	status=nf_close(ncid)	
		! PRINT*,MinVal(CER),MAXVAL(CER)
	where(cer<0 .and. cer.ne.-1) CER=32768*2+CER !! to unshort
	! print*,CER(2646,430:436)	
	CER_GOES=CER*0.00152602
	where(cer.eq.-1) CER_GOES=-1   !! -1 as fillvalue
	! print*,CER_GOES(2646,430:436)
	IF(status.NE.0)	AG=.FALSE.

END SUBROUTINE read_GOES_CER

!! cloud optical depth
SUBROUTINE read_GOES_COT(yyyymmdd,HHMM,GEOS_L2_DIR,COT_GOES,AG) 
	INCLUDE 'netcdf.inc'
	Character(*), INTENT(IN) :: GEOS_L2_DIR
	Character(8), INTENT(IN) :: yyyymmdd
	Character(4), INTENT(IN) :: HHMM
	REAL, DIMENSION(2712,2712), INTENT(OUT)  ::  COT_GOES	
	Integer, DIMENSION(2712,2712) ::  COD	
	LOGICAL, INTENT(OUT)  :: AG
	Integer status,ncid,varid
	! Character*255 :: filename	
	AG=.TRUE.
	
	! filename=TRIM(GEOS_L2_DIR)//"/"//yyyymmdd(1:4)//"/"//yyyymmdd// &
	! "/OR_ABI-L2-CODF-M6_G16_"//yyyymmdd//"_"//HHMM//".NC"
	
    status=nf_open(trim(adjustl(GEOS_L2_DIR)),nf_nowrite,ncid)	
	status=nf_inq_varid (ncid, 'COD', varid) 
	status=nf_get_var_int(ncid,varid,COD)
	status=nf_close(ncid)	
		! PRINT*,MinVal(COD),MAXVAL(COD)
	where(COD<0 .and. COD.ne.-1) COD=32768*2+COD !! to unshort
	! print*,COD(2646,430:436)	
	COT_GOES=COD*0.00244163
	where(COD.eq.-1) COT_GOES=-1   !! -1 as fillvalue
	! print*,COT_GOES(2646,430:436)
	IF(status.NE.0)	AG=.FALSE.

END SUBROUTINE read_GOES_COT


!! cloud top height (km)
SUBROUTINE read_GOES_CTH(yyyymmdd,HHMM,GEOS_L2_DIR,CTH_GOES,AG) 
	INCLUDE 'netcdf.inc'
	Character(*), INTENT(IN) :: GEOS_L2_DIR
	Character(8), INTENT(IN) :: yyyymmdd
	Character(4), INTENT(IN) :: HHMM
	REAL, DIMENSION(1086,1086), INTENT(OUT)  ::  CTH_GOES	
	Integer, DIMENSION(1086,1086) ::  CTH	
	LOGICAL, INTENT(OUT)  :: AG
	Integer status,ncid,varid
	! Character*255 :: filename	
	AG=.TRUE.
	
	! filename=TRIM(GEOS_L2_DIR)//"/"//yyyymmdd(1:4)//"/"//yyyymmdd// &
	! "/OR_ABI-L2-ACHAF-M6_G16_"//yyyymmdd//"_"//HHMM//".NC"
	
    status=nf_open(trim(adjustl(GEOS_L2_DIR)),nf_nowrite,ncid)	
	status=nf_inq_varid (ncid, 'HT', varid) 
	status=nf_get_var_int(ncid,varid,CTH)
	status=nf_close(ncid)	
		! PRINT*,MinVal(CTH),MAXVAL(CTH)
	where(CTH<0 .and. CTH.ne.-1) CTH=32768*2+CTH !! to unshort

	CTH_GOES=CTH* 0.3052037*0.001 !!! km
	where(CTH.eq.-1) CTH_GOES=-1   !! -1 as fillvalue
	
	IF(status.NE.0)	AG=.FALSE.

END SUBROUTINE read_GOES_CTH



!! DQF in CPSF file
SUBROUTINE read_GOES_DQF1(yyyymmdd,HHMM,GEOS_L2_DIR,DQF_GOES,AG)
	INCLUDE 'netcdf.inc'
	Character(*), INTENT(IN) :: GEOS_L2_DIR
	Character(8), INTENT(IN) :: yyyymmdd
	Character(4), INTENT(IN) :: HHMM
	INTEGER, DIMENSION(5424,5424), INTENT(OUT) ::  DQF_GOES	
	LOGICAL, INTENT(OUT) :: AG
	Integer status,ncid,varid
	! Character*255 :: filename
	AG=.TRUE.	
	! filename=TRIM(GEOS_L2_DIR)//"/"//yyyymmdd(1:4)//"/"//yyyymmdd// &
	! "/OR_ABI-L2-CPSF-M6_G16_"//yyyymmdd//"_"//HHMM//".NC"
	
    status=nf_open(trim(adjustl(GEOS_L2_DIR)),nf_nowrite,ncid)	
	status=nf_inq_varid (ncid, 'DQF', varid) 
	status=nf_get_var_int(ncid,varid,DQF_GOES)
	status=nf_close(ncid)	
	
	IF(status.NE.0) AG=.FALSE.
	
END SUBROUTINE read_GOES_DQF1

!! DQF in ACMF file
SUBROUTINE read_GOES_DQF2(yyyymmdd,HHMM,GEOS_L2_DIR,DQF_GOES,AG)
	INCLUDE 'netcdf.inc'
	Character(*), INTENT(IN) :: GEOS_L2_DIR
	Character(8), INTENT(IN) :: yyyymmdd
	Character(4), INTENT(IN) :: HHMM
	INTEGER, DIMENSION(5424,5424), INTENT(OUT) ::  DQF_GOES	
	LOGICAL, INTENT(OUT) :: AG
	Integer status,ncid,varid
	! Character*255 :: filename
	AG=.TRUE.	
	! filename=TRIM(GEOS_L2_DIR)//"/"//yyyymmdd(1:4)//"/"//yyyymmdd// &
	! "/OR_ABI-L2-ACMF-M6_G16_"//yyyymmdd//"_"//HHMM//".NC"
	
    status=nf_open(trim(adjustl(GEOS_L2_DIR)),nf_nowrite,ncid)	
	status=nf_inq_varid (ncid, 'DQF', varid) 
	status=nf_get_var_int(ncid,varid,DQF_GOES)
	status=nf_close(ncid)	
	
	IF(status.NE.0) AG=.FALSE.
	
END SUBROUTINE read_GOES_DQF2




! Navigating from Geodetic Latitude (φ ) and Longitude (λ ) to N/S Elevation Angle (y)
! and E/W Scanning Angle (x)
!
!  p4.8.2 https://www.goes-r.gov/products/docs/PUG-L2+-vol5.pdf
!
SUBROUTINE get_SZA_GOESR_LonLat(lambda16,lon,lat,x,y,FLAG)
  IMPLICIT NONE
	REAL, INTENT(IN) :: lon,lat,lambda16
	REAL, INTENT(OUT) :: x, y
	
	REAL :: lambda, phi
	REAL :: sx,sy,sz,phi_c,r_c
	REAL :: out_bound
	
	REAL, PARAMETER :: PAI = 3.14159265	
	REAL, PARAMETER :: R_eq=6378137 ! m
	REAL, PARAMETER :: R_pol=6356752.31414 ! m
	REAL, PARAMETER :: invf= 298.257222096 
	REAL            :: lambda0
	REAL, PARAMETER :: H0= 35786023 ! m
	REAL, PARAMETER :: H = 42164160 ! R_eq+H0
	REAL, PARAMETER :: e = 0.0818191910435

	LOGICAL,INTENT(OUT) :: FLAG
	
	lambda0= lambda16/180.*PAI
	phi= lat/180.*PAI
	lambda=lon/180.*PAI
	
	phi_c=ATAN(R_pol*R_pol/R_eq/R_eq*TAN(phi))
	r_c  =R_pol/sqrt(1-e*e*cos(phi_c)*cos(phi_c))
	
	sx=H-r_c*cos(phi_c)*cos(lambda-lambda0)
	sy=-1*r_c*cos(phi_c)*sin(lambda-lambda0)
	sz=r_c*sin(phi_c)
	
	out_bound=H*(H-sx)-sy*sy-sz*sz*R_eq*R_eq/R_pol/R_pol
	IF(out_bound<0) then
		FLAG= .FALSE.
		y=0.
		x=0.
		return
	ELSE 
		FLAG= .TRUE.
		y=atan(sz/sx)
		x=asin(-1*sy/sqrt(sx*sx+sy*sy+sz*sz))		
	END IF
	
END SUBROUTINE get_SZA_GOESR_LonLat


! Convert x,y to lon, lat
SUBROUTINE lonlat_GOES_R_xy(x,y,lambda,phi,FLAG)
  IMPLICIT NONE
	
	REAL ,INTENT(OUT):: lambda, phi
	REAL , INTENT(IN)::  x, y
	REAL :: sx,sy,sz,rs,a,b,c,delta

	
	REAL, PARAMETER :: PAI = 3.14159265	
	REAL, PARAMETER :: R_eq=6378137 ! m
	REAL, PARAMETER :: R_pol=6356752.31414 ! m
	REAL, PARAMETER :: invf= 298.257222096 
	REAL, PARAMETER :: lambda0= -75.0/180.*PAI
	REAL, PARAMETER :: H0= 35786023 ! m
	REAL, PARAMETER :: H = 42164160 ! R_eq+H0
	REAL, PARAMETER :: e = 0.0818191910435

	LOGICAL :: FLAG
	
	a=sin(x)*sin(x)+cos(x)*cos(x)*(cos(y)*cos(y)+R_eq*R_eq/R_pol/R_pol*sin(y)*sin(y))
	b=-2*H*cos(x)*cos(y)
	c=H*H-R_eq*R_eq
	
	! print*,a,b,c
	
	delta=b*b-4*a*c
	if(delta<0) THEN
		FLAG=.FALSE.
		lambda=0
		phi=0
	ELSE
		
		FLAG=.TRUE.
		rs=(-1*b-sqrt(delta))/2/a		
		sx=rs*cos(x)*cos(y)
		sy=-1*rs*sin(x)
		sz=rs*cos(x)*sin(y)
			! print*,rs,sx,sy,sz
		phi=atan(R_eq*R_eq/R_pol/R_pol*sz/sqrt((H-sx)*(H-sx)+sy*sy))*180/PAI 
		lambda=(lambda0-atan(sy/(H-sx)))*180/PAI 
	END IF
		
END SUBROUTINE lonlat_GOES_R_xy


