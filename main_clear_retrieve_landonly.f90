!! 2024.05.17 new branch fy3bcd-dev

!! main_clear_retrieve_landonly.f90
!! 2024.01.01 consider only clear sky, largely save download volumn and samples.
!! all subs are updated to only download and read cloud mask variable.


!! 2023.12.31 adapted from main_retrieve_gmi_emissivity_globe.f90
!! land only means we only retrieve land pixels, instead of retrieving all and do postprocesssing,
!! for the purpose of less time consuming, saving at least 2/3 time foreach

!! =============================================================================================
!!  main_retrieve_gmi_emissivity_allsky.f90
!!  Main program to retrieve land surface emissivity using GMI_Cloud_collocations：
!!	
!!  TOC:
!!  - Collocate TB and Clouds:
!!    The clouds come from three GeoSats that cover the whole globe except for high-latitudes:
!!		- Himawari-8&9/AHI over the Pacific Ocean 
!!	 	- GOES-16/ABI over the America
!!	 	- MSG-1&2&3/SERV over the EU and Africa
!!  - Collocate tb-cloud to ERA5 profiles and LST;
!!  - Call RTTOV V13.2 to retrieve multichannel allsky emissivity:
!! 	  Use the direct retrieving scheme porposed by (Baord et al., 2016), considering clearsky 
!!	  and non-precipitation cloudysky;
!!  - Post-processing to mask invalid result.  
!! 
!!  Requirements：
!!  	- NETCDF4, HDF5
!! 		- RTTOV V13.2
!!  On how to install, visit: 
!!      http://home.ustc.edu.cn/~hjh18305/space/research/rttov/rttov132-column/
!!
!!  By, Jiheng Hu, 2023.12, USTC& UMICH.
!!---------------------------------------------------------------------------------------------
!! 							DATA AVAILABLE RANGE 
!!---------------------------------------------------------------------------------------------
!! GMI L1C :  2014.03.04 - present, every orbits
!! H-8/9 AHI L2CLP : 2015.07.04 - present, every 10-min  
!! GOES-R ABI L2 Cloud : AWS (2019.339-present, 10-min)  CLASS/NCEI(2017.06.08-present) “hard”
!! MSG1/2 indian sea : 2017.2.1 - present  ONLY Cloud Mask , 15-min 
!! MSG3 0E : 2014.2.1 - present  ONLY Cloud Mask , 15-min 
!! ERA5: 1-hour, 0.25-degree, 37-layer
!! ============================================================================================

!! --------------------------------------------------------------------------------------------
!!     INCLUDE subroutines
!! --------------------------------------------------------------------------------------------
 ! Download

 INCLUDE 'subs/sub_download_ahi_l2c.f90'
 INCLUDE 'subs/sub_download_goesr_aws.f90'
 ! python 'subs/sub_download_msg_clm.py'
 ! python 'subs/sub_download_msg3_clm.py'
 ! I/O
 INCLUDE 'subs/sub_read_fy3_mwri_vars.f90'
 INCLUDE 'subs/sub_read_land_sea_mask.f90'
 INCLUDE 'subs/sub_read_ahi_l2c_vars.f90'
 INCLUDE 'subs/sub_read_goes_l2c_vars.f90'
 INCLUDE 'subs/sub_read_msg_iodc_clm.f90'

 ! Tools
 INCLUDE 'subs/sub_tools.f90'

 ! ERA5 I/O and Download  
 ! python "subs/era5/download_era5_profiles.py" 
 ! python "subs/era5/download_era5_sigles.py" 
 ! python "subs/era5/download_era5_lands.py" 
 INCLUDE 'subs/era5/sub_read_era5_vars.f90'
  
 ! RTTOV
 ! INCLUDE 'rttov_retrieve_mwri_emiss_landonly.f90'
 
 ! OUTPUT
 INCLUDE 'subs/sub_write_tb_cloud_hdf5.f90'
 INCLUDE 'subs/sub_write_fy3_emiss_hdf5.f90'
 INCLUDE 'subs/sub_write_profiles_hdf5.f90'
 
!! --------------------------------------------------------------------------------------------
!!     Main program
!! --------------------------------------------------------------------------------------------
PROGRAM main_clear_retrieve_landonly
  ! USE HDF5
  ! INCLUDE 'netcdf.inc'
  use, intrinsic :: ieee_arithmetic
  Implicit none
  CHARACTER(8) :: yyyymmdd
  CHARACTER(LEN=255) :: command
  CHARACTER(LEN=255) :: L1C_DIR,HIMA_DIR,GEOS_DIR,MSG_DIR
  INTEGER :: file_unit,text_unit
!! read FY3 Vars
  INTEGER :: nscan, npixel,nchannel, FYIO
  REAL*4, DIMENSION(:,:,:), ALLOCATABLE :: TB_swath,Emissivity_swath
  REAL*4, DIMENSION(:,:), ALLOCATABLE  :: Latitude, Longitude
  REAL*4, DIMENSION(:,:), ALLOCATABLE  :: ScanTime_swath
  REAL*4, DIMENSION(:), ALLOCATABLE  :: ScanTime
  CHARACTER(LEN=255) :: FY3_FILENAME, GOES_Sufix, FY3_Sufix
  CHARACTER(LEN=255) :: FY3_adscend,FY3_date,FY3_hhmm
  CHARACTER*4 :: FYSat

!! Land Sea Mask 10 km
  REAL*4, DIMENSION(3602,1800):: LSM_GRD

!! Himawari-8/GOES-R
  !! CLOUD FLAG:: -1:sea 0:out_of_boundary, 1:H8, 2:GOES-R, 3: MSG1,
  !!  4: Himawari missing; 5: GOES missing; 6: MSG1/2/3 missing;
  !!  7: H9; 8: MSG2; 9: MSG3
  INTEGER, DIMENSION(:,:), ALLOCATABLE  :: Cloud_Flag 
  INTEGER  :: cldflgmsg 

  REAL*4, DIMENSION(:,:), ALLOCATABLE  :: Sea_Frac  !! 0-100%
  REAL*4, DIMENSION(:,:,:), ALLOCATABLE ::CFR_swath,Clear_swath,Miss_swath 
  
  CHARACTER*255 :: AHI_FILENAME,GOES_FILENAME,AHI_FULLPATH,MSG_FILENAME,MSG_FULLPATH
  CHARACTER*255 :: record,lostfile,error_read
  CHARACTER*3 :: HIMA
  
  REAL*4, DIMENSION(2401) :: lat_h8,lon_h8

  INTEGER, DIMENSION(2401,2401) ::ctype 
  
  REAL*4, PARAMETER :: minlon_H8=80., maxlat_H8=60.
 
!! GOESR 
  REAL    :: xdeg,ydeg  !! rad

  CHARACTER*2 :: Mode
  INTEGER, DIMENSION(5424,5424)	:: CLM_GOES, DQF_CLM

  LOGICAL 						:: GOES_IN, A1, A2, A3, A4, A5, A6,A7
  LOGICAL 						:: FLAG1,FLAG2,FLAG3,FLAG4
  LOGICAL 						:: LZA, SAZ
 
  REAL :: x0
  real :: lonW,lonE,latN,latS
  REAL :: Sf,Nf,Wf,Ef
  INTEGER :: W2,E2,N2,S2,W3,E3,N3,S3,W1,E1,N1,S1
  INTEGER :: ngrid
 
  CHARACTER*255      		    :: GOES_CUR_STMAP,GOES_REC_STMAP,GOES_MISS_STMAP
  REAL*4 , PARAMETER			:: MIN_X_DEG=-0.151872, MAX_Y_DEG=0.151872 !! unit rad
  REAL*4 , PARAMETER			:: ddeg1 =MIN_X_DEG*(-2)/5424.   !! 2km nadir
  INTEGER 						:: ixg,iyg
  
!! MSG
  CHARACTER*4 :: MSGSAT
  CHARACTER*255 :: REC_MSG_STAMP,MSG_MISSING_STAMP,CURRENT_MSG_STAMP
  INTEGER,DIMENSION(3712,3712) :: CLM_MSG
  LOGICAL :: OK
  REAL :: Lambda0,lambda16
  INTEGER :: W,E,S,N, dum
 
!! general
  INTEGER:: iscan, ipixel,ilon,ilat,ix,iy
  REAL*4  :: lat, lon, stime, ls
  CHARACTER*2 :: HH, MM, MM15
  CHARACTER*3 :: DOY
  
  INTEGER , DIMENSION(11,11):: cltype,qabit,b7,b6
  INTEGER :: flag(11,11) ,ncloud,ncloud1, b13
  INTEGER :: qaflag(16)
  
  INTEGER :: STATUS
  LOGICAL :: EXISTS,EXISTSHDF,EXISTSTXT,FLAG_H8,FLAG_GOES,FLAG_MSG
  
!! Local view zenith angle, LZA
  REAL, DIMENSION(:,:),ALLOCATABLE :: LZA_swath
  REAL :: angle

!! ERA5
  CHARACTER*255 :: ERA5_HR1_STAMP,ERA5_HR2_STAMP, ERA5_DIR
  CHARACTER*255 :: req_era5_hr1, req_land_hr1, req_era5_hr2, req_land_hr2
  CHARACTER*2 :: HH1,HH2
  INTEGER :: UTC_hr1, UTC_hr2
  REAL    :: wt1,wt2
  LOGICAL :: profIO,landIO

  INTEGER,PARAMETER :: nlon=1440, nlat=561,nlevel=29, ntime=1
  REAL,PARAMETER    :: egrid=0.25
  REAL              :: elon(nlon), elat(nlat) 	


  ! grid params
  REAL, DIMENSION(nlevel) :: PLEVEL, PHALF !!hpa 	
  REAL, DIMENSION(nlon,nlat,ntime) :: lst1, psrf1, t2m1, snowc1, smc1
  REAL, DIMENSION(nlon,nlat,ntime) :: lst2, psrf2, t2m2, snowc2, smc2			
  REAL, DIMENSION(nlon,nlat,nlevel,ntime) :: qw1,ta1,zg1,rh1
  REAL, DIMENSION(nlon,nlat,nlevel,ntime) :: qw2,ta2,zg2,rh2	

  !!! swath params
  REAL, DIMENSION(:),ALLOCATABLE :: LST_land, Psrf_land, T2m_land, SnowC_land, SMC_land	
  REAL, DIMENSION(:),ALLOCATABLE :: Longitude_land,Latitude_land	
  REAL, DIMENSION(:,:),ALLOCATABLE :: LST_swath,T2m_swath,SnowC_swath,SMC_swath

  INTEGER, DIMENSION(:),ALLOCATABLE :: pixelid,scanid
  REAL, DIMENSION(:,:),ALLOCATABLE :: Vapor_land,AirTemp_land		

  REAL, DIMENSION(:,:),ALLOCATABLE :: TB_land, Emissivity_land
  LOGICAL :: alive1, alive2, alive3
  !!! tmp
  REAL, DIMENSION(:),ALLOCATABLE :: ScanTime_land,	LZA_land, Sea_land
  INTEGER, DIMENSION(:),ALLOCATABLE :: Cloudflg_land
  REAL, DIMENSION(:,:),ALLOCATABLE :: CFR_land,	Clear_land, Miss_land	
			
  !! retreve 
  CHARACTER*255 :: EMISS_OUTDIR,EMISS_HDF5,EMISS_TEXT
  CHARACTER*7 :: samples
  INTEGER :: imonth

  Integer :: num_land,iland
  
  logical :: REDO, HDF5,useERA5,ERA5land_Batchdownload

! ==================================================================================================
!	 1. SET UP VARIABLES
! ==================================================================================================	 
	FYSat="FY3B"
	
	!! Directory to save GMI_L1C HDF5 files 
	L1C_DIR = '/home/jihenghu/fy03/FY3B/ascend/'    
	HIMA_DIR = '/home/jihenghu/data04/AHI_L2/'    
	GEOS_DIR = '/home/jihenghu/data04/GOESR_CLM/'  
	MSG_DIR='/home/jihenghu/data04/MSG_CLM/'
	! ERA5_DIR = '/home/jihenghu/data00/data_em/ERA5/'
	ERA5_DIR = '/home/jihenghu/data00/ERA5/'

	!! OUTPUTs 
    EMISS_OUTDIR = '/home/jihenghu/data05/'//FYSat//'_EMISSIVITY/'    
	
	REDO=.False.  ! .True.!! redo retrieve or not?
	HDF5=.False.  ! .True.!! Output HDF orbits? Ascii format is mandatory
	
	useERA5=.False. ! True : use exist, dont download 
	ERA5land_Batchdownload=.False. !! directly download ERA5-land for 24 hours, and split
! ==================================================================================================

  ! Check if a command-line argument is provided
  IF (COMMAND_ARGUMENT_COUNT() /= 1) THEN
    WRITE(*,*) 'Error: Please provide a date string (YYYYMMDD) as a command-line argument.'
    STOP
  END IF

  ! Get the date string from the command-line argument
  CALL GETARG(1, yyyymmdd)
  IF (yyyymmdd<'20140101') THEN
    WRITE(*,*) 'Error: Please provide YYYYMMDD >= 2014.01.01.'
    STOP
  END IF  
  
  CALL system("mkdir  filelists")  
  CALL system("mkdir  logs") 
  
  CALL system("mkdir -p  "//trim(EMISS_OUTDIR))  
  CALL system("mkdir -p  "//trim(EMISS_OUTDIR)//"/"//yyyymmdd(1:4)//'/'//yyyymmdd)
  ERA5_DIR=trim(ERA5_DIR)//"/"//yyyymmdd(1:4)//'/'
  CALL system("mkdir -p  "//trim(ERA5_DIR)) 
  
! ==================================================================================================
! 	Search FY3 Files
! ==================================================================================================	 
  
  !! All FY3 L1C data Were downloaded prior to retrieval, thus no need to download...
   ! PRINT*,"Downloading FY3 files of the day ......"
   ! CALL download_GMI_L1C(yyyymmdd,L1C_DIR)
  
  ! Call the ls command to list HDF5 files in the directory
  ! FY3 sample be like: FY3B_MWRIA_GBAL_L1_20140101_0000_010KM_MS.HDF
  command ='ls '//TRIM(L1C_DIR)//"/"//yyyymmdd(1:4)//"/"//yyyymmdd//'/*.HDF'&
			// '> filelists/'//FYSat//'_filelist_'//yyyymmdd//'.txt' 

  CALL SYSTEM(command)


  ! Open a unit for reading the ls command output
  OPEN(NEWUNIT=file_unit, FILE='filelists/'//FYSat//'_filelist_'//yyyymmdd//'.txt')

!! ==================================================================================================
!!  Loop over FY3 File to retrieve
!! ==================================================================================================

176	READ(file_unit,'(A)',end=886) FY3_FILENAME
	! /data/jihenghu/data.fy3b/2014/20140701/FY3B_MWRIA_GBAL_L1_20140701_0026_010KM_MS.HDF

   	CALL GetDesiredPiece(FY3_FILENAME, '_', 2, FY3_adscend)
	CALL GetDesiredPiece(FY3_FILENAME, '_', 5, FY3_date)
	CALL GetDesiredPiece(FY3_FILENAME, '_', 6, FY3_hhmm)
	
	FY3_Sufix=FYSat//"."//trim(FY3_adscend)//'.'//trim(FY3_date)//'_'//trim(FY3_hhmm)


	EMISS_HDF5=trim(EMISS_OUTDIR)//"/"//yyyymmdd(1:4)//'/'//yyyymmdd//"/Emissivity.Clear."&
				//trim(FY3_Sufix)//".HDF5"
	EMISS_TEXT=trim(EMISS_OUTDIR)//"/"//yyyymmdd(1:4)//'/'//yyyymmdd//"/Emissivity.Clear."&
				//trim(FY3_Sufix)//".txt"
					
 	!! decide to retrieve ??	
	inquire(file=trim(EMISS_HDF5), exist=EXISTSHDF)
	inquire(file=trim(EMISS_TEXT), exist=EXISTSTXT)
	IF ((EXISTSHDF.OR.EXISTSTXT) .AND. (.NOT.REDO)) THEN
		PRINT*,"├──> Old Retrieval found, skip this orbit: "//trim(FY3_Sufix)
		GOTO 176
	END IF
	
	IF (ERA5land_Batchdownload) THEN
	 call system("python3 subs/era5/download_era5_lands_batch.py "//yyyymmdd//" "//trim(ERA5_DIR)) 
	END IF

   !! --------------------------------------------------------------------------------------------------
   !!		  read FY3 Variables                   
   !! --------------------------------------------------------------------------------------------------	
    PRINT *, '├──> Read FY3 file: '//trim(FY3_Sufix)//"_******.HDF" 
	! get HDF5 file Dimension

	CALL get_fy3bcd_swath_geometry(FY3_FILENAME,nscan,npixel,nchannel,FYIO)
	IF (FYIO.NE.0) goto 176
	
	! print*,nscan,npixel,nchannel 
	! Call the subroutine to read data	
	ALLOCATE(Latitude(npixel,nscan), Longitude(npixel,nscan), LZA_swath(npixel,nscan))	
	ALLOCATE(TB_swath(npixel,nscan,nchannel), ScanTime_swath(npixel,nscan), ScanTime(nscan))	
	CALL read_fy3bcd_Vars(FY3_FILENAME, TB_swath, Latitude, Longitude, ScanTime, nscan, npixel, nchannel,FYIO)
	
	IF (FYIO.NE.0) then 
		deallocate(TB_swath, Latitude, Longitude, ScanTime,LZA_swath,ScanTime_swath)
		goto 176
	end if
	!! read global 0.01deg land sea mask: land 0% - 25% (strict land) - 75% (coastal) - 100% strict sea
	!! longitude[3602] -0.05 ~ 360.05  ; latitude[1800]  -89.95 ~ 89.95
	CALL read_land_sea_mask("subs/IMERG_land_sea_mask.nc",LSM_GRD)
		
	! --------------------------------------------------------------------------------------------------
	!		 Allocate Himawari-8/GOES-R/MSG cloud related vars          
	! --------------------------------------------------------------------------------------------------		
	print*, "├────Cloud collocation beginning ......."
	
	ALLOCATE(CFR_swath(3,npixel,nscan))
	ALLOCATE(Miss_swath(3,npixel,nscan),Clear_swath(3,npixel,nscan))
	ALLOCATE(Sea_Frac(npixel,nscan))
	ALLOCATE(Cloud_Flag(npixel,nscan)) 
	
	
	CFR_swath=0
	Clear_swath=0
	Sea_Frac=0
	Cloud_Flag=0
	Miss_swath=100.
	LZA_swath=-999
   !! --------------------------------------------------------------------------------------------------
   !!		2. Loop over pixels and scans and match Himawari-8/GOES-R/MSG cloud             
   !! --------------------------------------------------------------------------------------------------	
	record=""
	lostfile=""
	error_read=""

	GOES_REC_STMAP="_"
	GOES_MISS_STMAP="_"

    REC_MSG_STAMP="_"
    MSG_MISSING_STAMP="_"

	
	FLAG_H8=.FALSE. 
	FLAG_GOES=.FALSE.  
	FLAG_MSG=.FALSE.  

	!!! ----------------------------do loop -----------------------------------------------------------
	do iscan=1,nscan
		
		ScanTime_swath(:,iscan)= ScanTime(iscan)
			
	  do ipixel=1,npixel
		lon=Longitude(ipixel,iscan)
		lat=Latitude(ipixel,iscan)
			
		if (lon.le.-180.or.lon.ge.180) cycle 
		if (lat.le.-90.or.lat.ge.90) cycle 
		if (ieee_is_nan(lat).or.ieee_is_nan(lon))cycle
		
		ilon=int((lon)/0.1)+1
		if (ilon.lt.0) ilon=ilon+3601	
		ilat=int((lat+90)/0.1)+1
		
		Sea_Frac(ipixel,iscan)=LSM_GRD(ilon,ilat)
	
		! if(lat.ge.70 .or. lat.le. -70) cycle    !! Out of boundary

		if(Sea_Frac(ipixel,iscan).ge.25) then
			Cloud_Flag(ipixel,iscan)=-1            !! sea/open water, will not retrieve
			cycle !! too much water   

		end if

	!! ==================================================================================================
	!!  get cloud from  H8/GOESR/MSG    
	!! ==================================================================================================
		stime=ScanTime(iscan)

		write(HH,'(i0.2)') int(stime)
		write(MM,'(i0.2)')  int(int((stime-int(stime))*60)/10)*10 
		
		! print*,stime
		!! ------------------------------------------------------------------------------------------------
		!!  Distribute the sample to branch proceduces :   
		!!  2878 :  Himawari-8 process   CLM, COT, CTH, CWP
		!!  9642 :  GOES-R process       CLM, COT, CTH, CWP
		!!  8847 :  MSG process  		 only CLM
		!!  8404 :  not found 
		!! ------------------------------------------------------------------------------------------------
		! PRINT*, lon, lat
		IF((lon.gt.80 .and. lon .le. 180).OR.((lon.ge.-180 .and. lon .lt. -160))) goto 2878   !! go to H-8		
		IF(lon.gt.-156 .and. lon .lt. -10) then
			if (yyyymmdd>'20170419') GOTO 9642     !! goto GOESR	
		end if
		IF(lon.gt.-67.5 .and. lon .lt. 113 ) GOTO  8847       !! goto MSG				
		GOTO 8404 !! goto Hell, no geostationary satellite mission region
		
	    !! ----------------------------------------------------------------------------------------------
	    !!   Himawari-8 Domain  
	    !! ----------------------------------------------------------------------------------------------
! H-8
2878    continue    		
		
		if (yyyymmdd<'20150704')  GOTO 8404
		IF (yyyymmdd>'20221212') THEN
			HIMA='H09'
		ELSE 
			HIMA='H08'
		END IF
		AHI_FILENAME="NC_"//HIMA//"_"//yyyymmdd//"_"//HH//MM//"_L2CLP010_FLDK.02401_02401.nc" 
		
		IF (trim(AHI_FILENAME) .eq. trim(error_read)) then ! file with ill
			Cloud_Flag(ipixel,iscan)=4                  
			goto 9642   !! OTHER GOSAT has it 
		END IF
		
		!! IF a totally new file needed
		IF (trim(AHI_FILENAME) .ne. trim(record)) THEN 
			record=AHI_FILENAME

			AHI_FULLPATH=trim(HIMA_DIR)//yyyymmdd//'/'//HH//'/'//trim(AHI_FILENAME)
			inquire(file=trim(AHI_FULLPATH), exist=EXISTS) !! 去指定文件夹寻找

			!! --------------------------------------------------------------------------
			!!     download Himawari-8 AHI L2CLP netcdf from online archive
			!! --------------------------------------------------------------------------
			if (.not.EXISTS) then !!! 其实不用检查，wget -c -N 会自动给检查	
				!! Not found, download
				CALL download_AHI_L2CLP(yyyymmdd,HH,AHI_FILENAME,HIMA_DIR,status)
				
				!! error download, record the lost				
				if(status.ne.0) then
					if(trim(adjustl(lostfile)).ne.trim(adjustl(AHI_FILENAME))) then           			
						lostfile=AHI_FILENAME   !!  new file download failed      
						print*, "│  │  ├── New H-8 file download failed: "//trim(AHI_FILENAME)	
					end if   
					
					Cloud_Flag(ipixel,iscan)=4
					error_read=AHI_FILENAME
					goto 9642   !!! go to GOES-R or MSG or END
				else 
					print*, "│  │  ├── New H-8 file download success :"// trim(AHI_FILENAME)
				end if
			else
				print*, "│  │  ├── Old H-8 file found, no download needed"
			endif
			
			!! --------------------------------------------------------------------------					
			!! 				Read Himawari-8/9 Cloud Proprties
			!! --------------------------------------------------------------------------
			
			CALL read_AHI_L2_cloud(AHI_FULLPATH,lat_h8,lon_h8,ctype,status) ! E80~-160; N60~-60  cer[um]  cth[km]
																		
			IF(status.ne.0) then
				print*, "│  │  └──  Error in read: "//trim(AHI_FULLPATH)  
				Cloud_Flag(ipixel,iscan)=4
				!! if read failed, must record it before cycle,
				!! to avoid being recognized as a read file and skip the read, which may result in a uncorrect use of the last readable H-8 file.
				error_read=AHI_FILENAME 
				goto 9642  !!! go to GOES-R or MSG or END
			end if
		end if 
				
		!!-------------------------------------------------------------------------------
		!!   			EFOV average  
		!!   10 GHz: 5*5;  18/23/36 GHz: 3*3;  89 GHz: 1*1;
		!!-------------------------------------------------------------------------------
		! PRINT*,"│  │  └──Calculating EFOV cloud properties..."
		
		ilon=nint((lon-minlon_H8)/0.05)+1 ! 1-2401   
		ilat=nint((maxlat_H8-lat)/0.05)+1 ! 1-2401 

		if (ilat.lt.3 .OR. ilon.lt.3 .OR. ilat.gt.2399 .OR. ilon.gt.2399) then !! H-8 boundry, abandon
			Cloud_Flag(ipixel,iscan)=0
			goto 9642  !!! go to GOES-R or MSG or END
		end if	
 
		cltype=ctype(ilon-5:ilon+5,ilat-5:ilat+5)       

		
		!! cltype:  0=Clear, 1=Ci, 2=Cs, 3=Deep convection, 
		!!          4=Ac, 5=As, 6=Ns, 7=Cu, 8=Sc, 9=St, 10=Unknown, 255=Fill
	
		!!-------------------------------------------------------------------------------
		!!   !! missing rate [primitive]: num(missing)/ngrid  
		!!-------------------------------------------------------------------------------			
		flag=1 			
		Miss_swath(1,ipixel,iscan)= sum(flag, mask=cltype.gt.20.or.cltype.lt.0)/121.*100   !%
		Miss_swath(2,ipixel,iscan)= sum(flag(3:8,3:8), &
							mask=cltype(3:8,3:8).gt.20.or.cltype(3:8,3:8).lt.0)/36.*100   !%
							
		Miss_swath(3,ipixel,iscan)= sum(flag(5:6,5:6), &
							mask=cltype(5:6,5:6).gt.20.or.cltype(5:6,5:6).lt.0)/4.*100   !%

		!!-------------------------------------------------------------------------------
		!!   cloud fraction [primitive]:   num(cloud)/ngrid , 1-cltype=0:Clear - mising
		!!   be more strict for Clear_swath=1-cloud
		!!-------------------------------------------------------------------------------	
		flag=1 
		Clear_swath(1,ipixel,iscan)= sum(flag, mask=cltype.eq.0)/121.*100   !%
		Clear_swath(2,ipixel,iscan)= sum(flag(3:8,3:8), mask=cltype(3:8,3:8).eq.0)/36.*100   !%
		Clear_swath(3,ipixel,iscan)= sum(flag(5:6,5:6), mask=cltype(5:6,5:6).eq.0)/4.*100   !%
		! if (cltype(3,3).eq.0) Clear_swath(3,ipixel,iscan)= 100.
		
		CFR_swath(:,ipixel,iscan)=100.-Clear_swath(:,ipixel,iscan)-Miss_swath(:,ipixel,iscan)
		! where(CFR_swath(:,ipixel,iscan)<0) CFR_swath(:,ipixel,iscan)=0. !! will cost lots of time !!!!
		
		IF (Miss_swath(1,ipixel,iscan)>99.99) Cloud_Flag(ipixel,iscan)=4

		!! Marked as Went Through a H8 collocation
		FLAG_H8=.TRUE.
			
		CALL calc_LZA(lon,lat,140.7,0.,angle)	
		LZA_swath(ipixel, iscan) = 	angle
		
		!! IF all missing, may be due to the night no-retrieve of Himawari-8, 
		!! wil check MSG cloud mask western-china for complementary 
		IF(Cloud_Flag(ipixel,iscan).eq. 4)  THEN			
			GOTO 8847
		END IF
		IF (HIMA.eq.'H08')  Cloud_Flag(ipixel,iscan)= 1
		IF (HIMA.eq.'H09')  Cloud_Flag(ipixel,iscan)= 7
		GOTO 8404
		

		
		!!----------------------------------------------------------------------------------------------
		!!  GOES-R Domain  
		!!        NOTE: GOESR data before 2020 can only be accessed upon order on NOAA CLASS, 
		!!              AND one should be careful of the 15-min mode3 before 20190402	
		!!        ONLY GOESR .gt. 20200101 from  AWS archive USED HERE!!!
		!! 		FullDisk mode06 ( 10-min interval: 00,10,20,30,40,50 same as Himawari-8 )
		!!        https://www.goes-r.gov/products/docs/PUG-L2+-vol5.pdf
		!!   RESOLUTIONs (NADIR):: 56 urad (2 KM), 112 urad(4KM), 280 urad(10KM)
		!!   nadir shit info: 20171214 from -89.5 to -75.2
		!!----------------------------------------------------------------------------------------------
! GOES
9642    continue 
			
		IF(yyyymmdd.gt.'20190402') Mode='M6'
		IF(yyyymmdd.lt.'20190402') Mode='M3'
		IF(yyyymmdd.eq.'20190402') then
			if (HH.lt.'16')	then
				Mode='M3'
			else
				Mode='M6'
			end if
		END IF
		
		IF (Mode=='M3') WRITE(MM,'(i0.2)') int(int((stime-int(stime))*60)/15)*15  
		IF (Mode=='M6') WRITE(MM,'(i0.2)') int(int((stime-int(stime))*60)/10)*10
		
		lambda16=-75.2
		IF (yyyymmdd.lt.'20171214') lambda16=-89.5
		
		IF(.NOT.(lon.gt.-156 .and. lon .lt. -10)) GOTO 8847		
	    IF (yyyymmdd<'20170419') THEN
			Cloud_Flag(ipixel,iscan)=5
			GOTO 8847  !! MSG has it
		ENDIF
		! print*,lon,lat
		
		! lambda17=-137.0
		CALL get_SZA_GOESR_LonLat(lambda16,lon,lat,xdeg,ydeg,GOES_IN)
		IF (.NOT.GOES_IN) THEN
			! IF (yyyymmdd>'20211201') then
				! CALL get_SZA_GOESR_LonLat(lambda17,lon,lat,xdeg,ydeg,GOES_IN)
				! IF (GOES_IN) THEN
					! USE_GOESW=.True.
				! END IF
			! END IF
			
			Cloud_Flag(ipixel,iscan)=0
			GOTO 8847  ! MSG has it
		END IF  !! outof GOES-R sight
		
		GOES_CUR_STMAP="G16_"//Mode//"_"//yyyymmdd//"_"//HH//MM
		
		IF (TRIM(GOES_CUR_STMAP).EQ.TRIM(GOES_MISS_STMAP)) THEN
			Cloud_Flag(ipixel,iscan)=5
			GOTO 8847  ! MSG has it
		END IF
		
		IF (TRIM(GOES_CUR_STMAP).NE. TRIM(GOES_REC_STMAP)) THEN  !! NEW file	
			GOES_REC_STMAP=GOES_CUR_STMAP
			
			GOES_FILENAME=TRIM(GEOS_DIR)//"/"//yyyymmdd(1:4)//"/"//yyyymmdd//&
					"/OR_ABI-L2-ACMF-"//Mode//"_G16_"//yyyymmdd//"_"//HH//MM//".NC"
			INQUIRE(FILE=TRIM(GOES_FILENAME), EXIST=EXISTS) !! 去指定文件夹寻找
			IF (.NOT.EXISTS) THEN 
				IF(yyyymmdd.gt.'20191231') THEN
					!! ---------------------------------------------------------------------
					!!    CALL Sub to download from Amazon Web Services S3 Cloud Bucket	
					!! ---------------------------------------------------------------------
					print*, '│  │  ├── GOES-R files downloading ......'				
					CALL download_GOESR_AWS_all(yyyymmdd,HH,GEOS_DIR)  !!! download for all day
				ELSE
					PRINT*,"│  │  └── GOES file notfound, "//TRIM(GOES_FILENAME)
					GOES_MISS_STMAP=GOES_CUR_STMAP
					Cloud_Flag(ipixel,iscan)=5
					GOTO 8847				
				ENDIF
			ELSE 
				print*, '│  │  ├── GOES-R file aready exists, no downloading needed'
			END IF
			
			!! -------------------------------------------------------------------------	 
			!!  Cloud Mask DIMENSION(5424,5424)
			!!	-1: NaN; 0: clear;  1:cloudy
			!! -------------------------------------------------------------------------	
			CALL read_GOES_CLM(yyyymmdd,HH//MM,TRIM(GOES_FILENAME),CLM_GOES,A1)
				! print*,MinVal(CLM_GOES),MAXVAL(CLM_GOES)  
				! print*,CLM_GOES(:20,1)

			! DQF_CLM==0 Defines the best cloud mask which is to Seattle in the northwest
			!! ![屏幕截图 2023-12-06 175659.png](https://s2.loli.net/2023/12/07/wQVL6OT58Yt3e2I.png)
			!! --------------------------------------------------------------------------
			!!    Quality flag of Cloud Mask DIMENSION(5424,5424)
			!!    valid [-1, 6] : -1 :fill
			!!    0 good_quality_qf 
			!!    1 invalid_due_to_not_geolocated_or_algorithm_non_execution_qf 
			!!	  2 degraded_due_to_LZA_threshold_exceeded_qf 
			!!    3 invalid_due_to_bad_or_missing_input_band_14_brightness_temperature_qf 
			!!    4 degraded_due_to_bad_input_band_7_pixel_qf 
			!!    5 degraded_due_to_failed_band_2_tests_qf 
			!!    6 degraded_due_to_other_bad_bands_qf
			!! --------------------------------------------------------------------------
			CALL read_GOES_DQF2(yyyymmdd,HH//MM,TRIM(GOES_FILENAME),DQF_CLM,A7)  	  
				! PRINT*,MinVal(DQF_CLM),MAXVAL(DQF_CLM)

			IF(.NOT.(A1.AND.A7)) THEN  !! SOME file MISSING
				PRINT*,"│  │  └── READ GOES file error: "//TRIM(GOES_FILENAME)
				GOES_MISS_STMAP=GOES_CUR_STMAP
				Cloud_Flag(ipixel,iscan)=5
				GOTO 8847
			END IF
							
			PRINT*,"│  │  └── READ GOES CLOUDS DONE: "//TRIM(GOES_CUR_STMAP)
			! STOP
		END IF


		ixg=int((xdeg - MIN_X_DEG)/ddeg1)+1	 
		iyg=int((MAX_Y_DEG - ydeg)/ddeg1)+1   
		IF (DQF_CLM(ixg,iyg) .ne.0) THEN
			Cloud_Flag(ipixel,iscan)=0
			goto 8847		
		END IF
		
		
!! 6.925 GHz
		lonW=lon-0.25
		lonE=lon+0.25
		latN=lat+0.25
		latS=lat-0.25
		
		CALL get_SZA_GOESR_LonLat(lambda16,lonW,lat,Wf,x0,FLAG1)
		CALL get_SZA_GOESR_LonLat(lambda16,lonE,lat,Ef,x0,FLAG2)
		CALL get_SZA_GOESR_LonLat(lambda16,lon,latN,x0,Nf,FLAG3)
		CALL get_SZA_GOESR_LonLat(lambda16,lon,latS,x0,Sf,FLAG4)
		
		IF (.not.(FLAG1.OR.FLAG2.OR.FLAG3.OR.FLAG4)) THEN
			Cloud_Flag(ipixel,iscan)=0
			GOTO 8847	
		END IF
		
		!!! 2km NADIR
		W1=int((Wf - MIN_X_DEG)/ddeg1)+1	
		E1=int((Ef - MIN_X_DEG)/ddeg1)+1	
		N1=int((MAX_Y_DEG - Nf)/ddeg1)+1         
		S1=int((MAX_Y_DEG - Sf)/ddeg1)+1  

		
		! print*,W1,E1,N1,S1
		
		IF (W1/=E1.or.N1/=S1) THEN 
			ngrid=count(CLM_GOES(W1:E1,N1:S1).gt.-999)
			
			CFR_swath(1,ipixel,iscan)= COUNT(CLM_GOES(W1:E1,N1:S1).EQ.1)*100./ngrid
									   
			Clear_swath(1,ipixel,iscan)= COUNT(CLM_GOES(W1:E1,N1:S1).EQ.0)*100./ngrid
									   
			Miss_swath(1,ipixel,iscan)= COUNT(CLM_GOES(W1:E1,N1:S1).EQ.-1)*100./ngrid
			
		ELSE  !! Wf=Ef, Sf=Nf

			IF(CLM_GOES(W1,N1) .EQ. -1) Miss_swath(:,ipixel,iscan)=100.						   
			IF(CLM_GOES(W1,N1) .EQ. 0)  Clear_swath(:,ipixel,iscan)=100.						   
			IF(CLM_GOES(W1,N1) .EQ. 1)  CFR_swath(:,ipixel,iscan)=100.	
			
			FLAG_GOES=.TRUE. 
			CALL calc_LZA(lon,lat,lambda16,0.,angle)	
			LZA_swath(ipixel, iscan) = 	angle
			Cloud_Flag(ipixel,iscan)= 2				
			GOTO 8404 
		ENDIF


!!! 18.23.36 GHZ
		lonW=lon-0.15
		lonE=lon+0.15
		latN=lat+0.15
		latS=lat-0.15
		
		CALL get_SZA_GOESR_LonLat(lambda16,lonW,lat,Wf,x0,FLAG1)
		CALL get_SZA_GOESR_LonLat(lambda16,lonE,lat,Ef,x0,FLAG2)
		CALL get_SZA_GOESR_LonLat(lambda16,lon,latN,x0,Nf,FLAG3)
		CALL get_SZA_GOESR_LonLat(lambda16,lon,latS,x0,Sf,FLAG4)
		
		IF (.not.(FLAG1.OR.FLAG2.OR.FLAG3.OR.FLAG4)) THEN
			Cloud_Flag(ipixel,iscan)=0
			GOTO 8847	
		END IF
		
		!!! 2km NADIR
		W1=int((Wf - MIN_X_DEG)/ddeg1)+1	
		E1=int((Ef - MIN_X_DEG)/ddeg1)+1	
		N1=int((MAX_Y_DEG - Nf)/ddeg1)+1         
		S1=int((MAX_Y_DEG - Sf)/ddeg1)+1  

		
		
		IF (W1/=E1.or.N1/=S1) THEN 
			ngrid=count(CLM_GOES(W1:E1,N1:S1)>-999)
			
			CFR_swath(2,ipixel,iscan)=COUNT(CLM_GOES(W1:E1,N1:S1).EQ.1)*100./ngrid
									   
			Clear_swath(2,ipixel,iscan)= COUNT(CLM_GOES(W1:E1,N1:S1).EQ.0)*100./ngrid
									   
			Miss_swath(2,ipixel,iscan)= COUNT(CLM_GOES(W1:E1,N1:S1).EQ.-1)*100./ngrid
		
		ELSE

			IF (CLM_GOES(W1,N1) .EQ. -1) Miss_swath(2:3,ipixel,iscan)=100.						   
			IF (CLM_GOES(W1,N1) .EQ. 0)  Clear_swath(2:3,ipixel,iscan)=100.						   
			IF (CLM_GOES(W1,N1) .EQ. 1)  CFR_swath(2:3,ipixel,iscan)=100.	
			
			FLAG_GOES=.TRUE.  

			Cloud_Flag(ipixel,iscan)= 2		
	
			CALL calc_LZA(lon,lat,lambda16,0.,angle)	
			LZA_swath(ipixel, iscan) = 	angle	
			GOTO 8404
		ENDIF

!!! 89 GHZ
		lonW=lon-0.08
		lonE=lon+0.08
		latN=lat+0.08
		latS=lat-0.08
		CALL get_SZA_GOESR_LonLat(lambda16,lonW,lat,Wf,x0,FLAG1)
		CALL get_SZA_GOESR_LonLat(lambda16,lonE,lat,Ef,x0,FLAG2)
		CALL get_SZA_GOESR_LonLat(lambda16,lon,latN,x0,Nf,FLAG3)
		CALL get_SZA_GOESR_LonLat(lambda16,lon,latS,x0,Sf,FLAG4)
		
		IF (.not.(FLAG1.OR.FLAG2.OR.FLAG3.OR.FLAG4)) THEN
			Cloud_Flag(ipixel,iscan)=0
			GOTO 8847	
		END IF
		
		!!! 2km NADIR
		W1=int((Wf - MIN_X_DEG)/ddeg1)+1	
		E1=int((Ef - MIN_X_DEG)/ddeg1)+1	
		N1=int((MAX_Y_DEG - Nf)/ddeg1)+1         
		S1=int((MAX_Y_DEG - Sf)/ddeg1)+1  

		IF (W1/=E1.or.N1/=S1) THEN 
			ngrid=count(CLM_GOES(W1:E1,N1:S1)>-999)
			
			! print*,ngrid,ncloud,ncloud1
			CFR_swath(3,ipixel,iscan)= COUNT(CLM_GOES(W1:E1,N1:S1).EQ.1)*100./ngrid
									   
			Clear_swath(3,ipixel,iscan)= COUNT(CLM_GOES(W1:E1,N1:S1).EQ.0)*100./ngrid
									   
			Miss_swath(3,ipixel,iscan)= COUNT(CLM_GOES(W1:E1,N1:S1).EQ.-1)*100./ngrid
			
		ELSE
			IF (CLM_GOES(W1,N1) .EQ. -1) Miss_swath(3,ipixel,iscan)=100.						   
			IF (CLM_GOES(W1,N1) .EQ. 0)  Clear_swath(3,ipixel,iscan)=100.						   
			IF (CLM_GOES(W1,N1) .EQ. 1)  CFR_swath(3,ipixel,iscan)=100.	
		
		ENDIF
	
		!! Marked as Went Through a GOES-R collocation
	
		Cloud_Flag(ipixel,iscan)= 2		
		CALL calc_LZA(lon,lat,lambda16,0.,angle)	
		LZA_swath(ipixel, iscan) = 	angle			
		FLAG_GOES=.TRUE.  
		
		GOTO 8404

		
		!!----------------------------------------------------------------------------------------------
		!!  MSG Cloud Mask 
		!!  FullDisk : 15-min interval: 00,15,30,45
		!!  Instrument replace Warning: 
		!!    2022.06.01, Meteosat-9 (45.5°E) replaced Meteosat-8(41.5° E) to provide data
		!!    20220601 8:45(MSG1) -> 20220601 09:00(MSG2)
		!!    
		!!   2018-02-20 (09:00) 	MSG4				
		!!   2018-02-20 (08:45) 	MSG3		
		!!----------------------------------------------------------------------------------------------
! MSG   WE have only one variable to deal with, ie., CLM		
8847    CONTINUE 	

		
		IF(.NOT.(lon.gt.-67.5 .and. lon .lt. 113) ) GOTO 8404  !! out of MSG domain
		IF(yyyymmdd<'20140101') THEN
			Cloud_Flag(ipixel,iscan)=0
			GOTO 8404
		END IF


		IF (yyyymmdd<'20220601') MSGSAT='MSG1'
		IF (yyyymmdd>'20220601') MSGSAT='MSG2'
		IF (yyyymmdd<'20170201') THEN
			MSGSAT='MSG3'
			IF(lon.gt.67.5) goto 8404
		END IF
		
		IF (lon.lt.-20 .and. lon.gt.-67) then
			IF (yyyymmdd>'20170201'.and. yyyymmdd<'20170608') then
				MSGSAT='MSG3'
			ENDIF
		ENDIF
			
		IF (yyyymmdd.eq.'20220601') THEN
			IF(HH.lt.'09') THEN
				MSGSAT='MSG1'
			ELSE
				MSGSAT='MSG2'
			END IF
		END IF
		
		IF (MSGSAT.EQ.'MSG1') then
			Lambda0= 41.5
			cldflgmsg=3
		ELSE IF (MSGSAT.EQ.'MSG2') THEN
			Lambda0= 45.5
			cldflgmsg=8
		ELSE IF (MSGSAT.EQ.'MSG3') THEN
			Lambda0= 0
			cldflgmsg=9
		ELSE IF (MSGSAT.EQ.'MSG4') THEN
			Lambda0= 0
		ELSE
			PRINT*, "│  │  └── ERROR: Un—Specified MSG satellite: ", MSGSAT
			Cloud_Flag(ipixel,iscan)=6
			GOTO 8404 
		END IF
		
		WRITE(MM15,'(i0.2)')  int(int((stime-int(stime))*60)/15)*15 

		CURRENT_MSG_STAMP=MSGSAT//"_"//yyyymmdd//"_"//HH//MM15
		
		IF (TRIM(CURRENT_MSG_STAMP).EQ. TRIM(MSG_MISSING_STAMP)) THEN	!! a recorded missing file
			Cloud_Flag(ipixel,iscan)=6
			GOTO 8404  
		END IF 
				
		IF (TRIM(CURRENT_MSG_STAMP).NE. TRIM(REC_MSG_STAMP)) THEN  !! we need a new file
			REC_MSG_STAMP=CURRENT_MSG_STAMP
			
			!!! search the file to open
			MSG_FULLPATH=TRIM(MSG_DIR)//yyyymmdd//'/'//MSGSAT// &
					'-SEVI-MSGCLMK-0100-0100-'//yyyymmdd//HH//MM15//'00.000000000Z-NA.grb'
			INQUIRE(FILE=TRIM(MSG_FULLPATH), EXIST=EXISTS) 
			! PRINT*,TRIM(MSG_FULLPATH)
			IF (.NOT.EXISTS) THEN 
				PRINT*,"│  │  ├── New MSG1/2/3 download needed"
				!! -----------------------------------------------------------------------------------
				!! Call python script to download from EUMETSAT DC  :: pip(>3.7) install eumdac
				!! -----------------------------------------------------------------------------------
			  IF(MSGSAT.EQ.'MSG3') THEN
				CALL SYSTEM("python3.9 subs/sub_download_msg3_clm.py "//yyyymmdd//" "//HH//MM15//" "//trim(MSG_DIR))			
			  ELSE
				CALL SYSTEM("python3.9 subs/sub_download_msg_clm.py "//yyyymmdd//" "//HH//MM15//" "//trim(MSG_DIR))			
			  END IF
			ELSE
				PRINT*,"│  │  ├── OLD MSG1/2/3 Found, no download needed"
			END IF
			!! --------------------------------------------------------------------------------------
			!!	 CLOUD MASK :  Code table 4.217
			!!  0 Clear over water; 1 Clear over land; 2 Cloud; 3 No data
			!! --------------------------------------------------------------------------------------
			!!! 手动安装eccodes:: https://confluence.ecmwf.int/display/ECC/ecCodes+installation
			!!! export LD_LIBRARY_PATH=/home/jihenghu/eccodes/lib64:$LD_LIBRARY_PATH
			!!! compile using option: -leccodes_f90 -leccodes 
			!!! or -L/home/jihenghu/eccodes/lib64/  -leccodes_f90 -leccodes  -I/home/jihenghu/eccodes/include/
			CLM_MSG=-1
			! print*,"HERE0", MSG_FULLPATH

			CALL read_MSG_CLM(MSG_FULLPATH,CLM_MSG,status)	
			
			IF (status.NE.0 .OR.(MINVAL(CLM_MSG)<0.and.MAXVAL(CLM_MSG)<0)) THEN         
				PRINT*,"│  │  └── READ MSG CLOUD ERROR**:", TRIM(CURRENT_MSG_STAMP)
				MSG_MISSING_STAMP=CURRENT_MSG_STAMP
				Cloud_Flag(ipixel,iscan)=6				    
				GOTO 8404 
			ELSE 
				PRINT*,"│  │  └── READ MSG CLOUD DONE: ", TRIM(CURRENT_MSG_STAMP)	
			END IF

			WHERE(CLM_MSG>3.OR.CLM_MSG<0) CLM_MSG=3			
		END IF
	

		!! now that we have the cloud mask,NEW OR OLD, we need to locate the sample to its disk.
		! https://www.cgms-info.org/wp-content/uploads/2021/10/cgms-lrit-hrit-global-specification-(v2-8-of-30-oct-2013).pdf
		! MSG is nominally 3km- resolutioned at nadir and increases with LZA

		! Download the MSG L1.5 geolaction calculation F90 programe to transfer lon,lat -> index ix,iy, from
		! https://navigator.eumetsat.int/product/EO:EUM:SW:MSG:159. we adapted as subrountine here.
		! This convert the ix =0, iy=0 at the south-east corner, 
		! however,  The CML data is stored with the south-western corner as the (0,0) point, with is [-5500km,-5500km].
		! so there has to be a ix=3712-ix+1, OR the disk will be west-east overturned.

		CALL geocoord2pixcoord(lat,lon, Lambda0, ix, iy, FLAG1)
		! print*, lon, lat, ix,  iy 
		ix=3712-ix+1
		IF (CLM_MSG(ix,iy).eq.3) THEN
			Cloud_Flag(ipixel,iscan)= 0	
			GOTO 8404
		END IF

		!!! ===============================================================================================================
		!! 6.925 GHz ----------------------------------------------------------------------------------
		lonW=lon-0.25
		lonE=lon+0.25
		latN=lat+0.25
		latS=lat-0.25
		
		CALL geocoord2pixcoord(lat,lonW,Lambda0, W,dum,FLAG1)   !! E<W
		CALL geocoord2pixcoord(lat,lonE,Lambda0, E,dum,FLAG2)
		CALL geocoord2pixcoord(latN,lon,Lambda0, dum,N,FLAG3)   !! S<N
		CALL geocoord2pixcoord(latS,lon,Lambda0, dum,S,FLAG4)
		
		W=3712-W+1
		E=3712-E+1  !! W<E
		
		IF (.not.(FLAG1.OR.FLAG2.OR.FLAG3.OR.FLAG4)) THEN
			Cloud_Flag(ipixel,iscan)=0
			GOTO 8404
		END IF

		IF (W/=E.or.N/=S) THEN 
			ngrid=count(CLM_MSG(W:E,S:N)<4)
			ncloud=count(CLM_MSG(W:E,S:N).EQ.2)
			
			CFR_swath(1,ipixel,iscan)= COUNT(CLM_MSG(W:E,S:N).EQ.2)*100./ngrid										   
			Clear_swath(1,ipixel,iscan)= COUNT(CLM_MSG(W:E,S:N).EQ.0.OR.CLM_MSG(W:E,S:N).EQ.1)*100./ngrid										   
			Miss_swath(1,ipixel,iscan)= COUNT(CLM_MSG(W:E,S:N).EQ.3)*100./ngrid		
			
		ELSE  !! W=E, S=N

			IF(CLM_MSG(W,N) .EQ. 3) Miss_swath(:,ipixel,iscan)=100.						   
			IF(CLM_MSG(W,N) .EQ. 1.OR. CLM_MSG(W,N) .EQ. 0)  Clear_swath(:,ipixel,iscan)=100.						   
			IF(CLM_MSG(W,N) .EQ. 2)  CFR_swath(:,ipixel,iscan)=100.	
					
			FLAG_MSG=.True. 
			! print*,FLAG_MSG
			Cloud_Flag(ipixel,iscan)= cldflgmsg
			CALL calc_LZA(lon,lat,Lambda0,0.,angle)	
			LZA_swath(ipixel, iscan) = 	angle	
			GOTO 8404 
		ENDIF

		
		!! 18.23.36 ----------------------------------------------------------------------------------
		lonW=lon-0.15
		lonE=lon+0.15
		latN=lat+0.15
		latS=lat-0.15
		
		CALL geocoord2pixcoord(lat,lonW,Lambda0, W,dum,FLAG1)   !! E<W
		CALL geocoord2pixcoord(lat,lonE,Lambda0, E,dum,FLAG2)
		CALL geocoord2pixcoord(latN,lon,Lambda0, dum,N,FLAG3)   !! S<N
		CALL geocoord2pixcoord(latS,lon,Lambda0, dum,S,FLAG4)
		W=3712-W+1
		E=3712-E+1  !! W<E
		
		IF (.not.(FLAG1.OR.FLAG2.OR.FLAG3.OR.FLAG4)) THEN
			Cloud_Flag(ipixel,iscan)=0
			GOTO 8404
		END IF

		IF (W/=E.or.N/=S) THEN 
			ngrid=count(CLM_MSG(W:E,S:N)<4)
			ncloud=count(CLM_MSG(W:E,S:N).EQ.2)
			
			CFR_swath(2,ipixel,iscan)= COUNT(CLM_MSG(W:E,S:N).EQ.2)*100./ngrid										   
			Clear_swath(2,ipixel,iscan)= COUNT(CLM_MSG(W:E,S:N).EQ.0.OR.CLM_MSG(W:E,S:N).EQ.1)*100./ngrid										   
			Miss_swath(2,ipixel,iscan)= COUNT(CLM_MSG(W:E,S:N).EQ.3)*100./ngrid		
		
		ELSE  !! W=E, S=N

			IF(CLM_MSG(W,N) .EQ. 3) Miss_swath(2:3,ipixel,iscan)=100.						   
			IF(CLM_MSG(W,N) .EQ. 1.OR. CLM_MSG(W,N) .EQ. 0)  Clear_swath(2:3,ipixel,iscan)=100.						   
			IF(CLM_MSG(W,N) .EQ. 2)  CFR_swath(2:3,ipixel,iscan)=100.	
					
			FLAG_MSG=.True. 
			! print*,FLAG_MSG
			Cloud_Flag(ipixel,iscan)= cldflgmsg
			CALL calc_LZA(lon,lat,Lambda0,0.,angle)	
			LZA_swath(ipixel, iscan) = 	angle	
			GOTO 8404 
		ENDIF

		!! 89 ----------------------------------------------------------------------------------
		lonW=lon-0.08
		lonE=lon+0.08
		latN=lat+0.08
		latS=lat-0.08
		
		CALL geocoord2pixcoord(lat,lonW,Lambda0, W,dum,FLAG1)   !! E<W
		CALL geocoord2pixcoord(lat,lonE,Lambda0, E,dum,FLAG2)
		CALL geocoord2pixcoord(latN,lon,Lambda0, dum,N,FLAG3)   !! S<N
		CALL geocoord2pixcoord(latS,lon,Lambda0, dum,S,FLAG4)
		W=3712-W+1  
		E=3712-E+1    !! W<E
		
		IF (.not.(FLAG1.OR.FLAG2.OR.FLAG3.OR.FLAG4)) THEN
			Cloud_Flag(ipixel,iscan)=0
			GOTO 8404
		END IF

		IF (W/=E.or.N/=S) THEN 
			ngrid=count(CLM_MSG(W:E,S:N)<4)
			ncloud=count(CLM_MSG(W:E,S:N).EQ.2)
			
			CFR_swath(3,ipixel,iscan)= COUNT(CLM_MSG(W:E,S:N).EQ.2)*100./ngrid										   
			Clear_swath(3,ipixel,iscan)= COUNT(CLM_MSG(W:E,S:N).EQ.0.OR.CLM_MSG(W:E,S:N).EQ.1)*100./ngrid										   
			Miss_swath(3,ipixel,iscan)= COUNT(CLM_MSG(W:E,S:N).EQ.3)*100./ngrid	
		
		ELSE  !! W=E, S=N

			IF(CLM_MSG(W,N) .EQ. 3) Miss_swath(3,ipixel,iscan)=100.						   
			IF(CLM_MSG(W,N) .EQ. 1.OR. CLM_MSG(W,N) .EQ. 0)  Clear_swath(3,ipixel,iscan)=100.						   
			IF(CLM_MSG(W,N) .EQ. 2)  CFR_swath(3,ipixel,iscan)=100.	
					
			FLAG_MSG=.True. 
			! print*,FLAG_MSG
			Cloud_Flag(ipixel,iscan)= cldflgmsg
			CALL calc_LZA(lon,lat,Lambda0,0.,angle)	
			LZA_swath(ipixel, iscan) = 	angle	
			GOTO 8404 
		ENDIF

		!!! Marked as Went Through a MSG collocation
		FLAG_MSG=.TRUE.
		! print*,FLAG_MSG
		Cloud_Flag(ipixel,iscan)= cldflgmsg		
		CALL calc_LZA(lon,lat,Lambda0,0.,angle)	
		LZA_swath(ipixel, iscan) = 	angle	
		

		
		!!----------------------------------------------------------------------------------------------
		!!    END 
		!!----------------------------------------------------------------------------------------------	 		
8404    continue
	  !! going to next [iscan, ipixel]
	  end do
	end do
	
	! print*,FLAG_MSG
	IF( .NOT. ((FLAG_H8.OR.FLAG_GOES.OR.FLAG_MSG).AND.(maxval(CFR_swath)+maxval(Clear_swath)>0)) ) THEN
		print*,"├──>> No valid sample,  goto 2333"
		GOTO 2333
	END IF
	

	!!===========================================================================================================
	!!                     3. ERA5 procedure
	!!===========================================================================================================

	
	where(CFR_swath>10.) Miss_swath=100  !!  cloudy
	
	num_land=COUNT(Miss_swath(1,:,:).lt.100)

	print*,nscan*npixel,COUNT(Sea_Frac.lt.50),num_land
	IF (num_land.lt.10) goto 2333
	!! operate in 1D londonly pixel
	print*, "├────ERA5 procedure beginning ........."	 	
	ALLOCATE(pixelid(num_land),scanid(num_land) )		
	ALLOCATE(LST_land(num_land), Psrf_land(num_land), T2m_land(num_land), SnowC_land(num_land), SMC_land(num_land))	

	ALLOCATE(Vapor_land(nlevel,num_land),AirTemp_land(nlevel,num_land))!,rh(nlevel,npixel,nscan)		

	ALLOCATE(TB_land(nchannel,num_land))

	ALLOCATE(Longitude_land(num_land),Latitude_land(num_land))
	ALLOCATE(ScanTime_land(num_land),LZA_land(num_land),Sea_land(num_land),Cloudflg_land(num_land))
	ALLOCATE(CFR_land(3,num_land),Clear_land(3,num_land),Miss_land(3,num_land))
	
	ALLOCATE(Emissivity_land(nchannel,num_land))
		Emissivity_land=0.0
	
	ERA5_HR1_STAMP="_"
	ERA5_HR2_STAMP="_"

	PLEVEL=(/50.0,70.0,100.0,125.0,150.0,175.0,200.0,225.0,250.0,300.0,350.0,400.0,450.0,500.0	&
		,550.0,600.0,650.0,700.0,750.0,775.0,800.0,825.0,850.0,875.0,900.0,925.0,950.0,975.0,1000.0/)  !! hpa
	
	PHALF=(/40.,60.,85.,112.5,137.5,162.5,187.5,212.5,237.5,275.,325.,375.,425.,475.,525.,575. &
		,625.,675.,725.,762.5,787.5,812.5,837.5,862.5,887.5,912.5,937.5,962.5,987.5/)
	
	iland=0
	do iscan=1, nscan
		do ipixel=1, npixel
			
			! IF(Sea_Frac(ipixel,iscan).ge.50) cycle
			IF(Miss_swath(1,ipixel,iscan).lt.100) then
			
			iland=iland+1
			
			pixelid(iland)=ipixel
			scanid(iland) =iscan
			
			TB_land(:,iland)=TB_swath(ipixel,iscan,:)
			Longitude_land(iland)=Longitude(ipixel,iscan)
			Latitude_land(iland) =Latitude(ipixel,iscan)
			ScanTime_land(iland)=ScanTime(iscan)
			LZA_land(iland)=LZA_swath(ipixel,iscan)
			Sea_land(iland)=Sea_Frac(ipixel,iscan)		
			Cloudflg_land(iland)=Cloud_Flag(ipixel,iscan)
			
			CFR_land(:,iland)=CFR_swath(:,ipixel,iscan)
			Clear_land(:,iland)=Clear_swath(:,ipixel,iscan)
			Miss_land(:,iland)=Miss_swath(:,ipixel,iscan)
			
			stime=ScanTime(iscan)     !! 0.33 hr    23:45  
			UTC_hr1=int(stime)        !!  0 hr      23      
			UTC_hr2=int(stime)+1      !!  1 hr      24->0
			
			wt1=(UTC_hr2-stime)*1.0/(UTC_hr2-UTC_hr1)
			wt2=(stime-UTC_hr1)*1.0/(UTC_hr2-UTC_hr1)
			
			if (UTC_hr1<0) UTC_hr1=0
			if (UTC_hr2>23) UTC_hr2=UTC_hr2-24
			
			write(HH1,'(i0.2)') UTC_hr1
			write(HH2,'(i0.2)') UTC_hr2
			
			!!--------------------------------------------------------------------------------------------------
			!!      download Global ERA5 reanaysis:  [70, -180, -70, 180,]
			!!--------------------------------------------------------------------------------------------------
			!! Hour1
			req_era5_hr1=trim(adjustl(ERA5_DIR))//"/"//yyyymmdd//"/ERA5-PL-GBL-25km-"//yyyymmdd//"-"//HH1//"00.nc"
			req_land_hr1=trim(adjustl(ERA5_DIR))//"/"//yyyymmdd//"/ERA5-Land-GBL-25km-"//yyyymmdd//"-"//HH1//"00.nc"

			IF (trim(req_era5_hr1).NE.trim(ERA5_HR1_STAMP)) THEN
				ERA5_HR1_STAMP=req_era5_hr1
				PRINT*, "│  ├──New ERA5 file need for HR1: ",HH1,":00"
					
				inquire(file=trim(adjustl(req_era5_hr1)), exist=alive1)
				inquire(file=trim(adjustl(req_land_hr1)), exist=alive2)
				
				IF (((.not.alive1).or.(.not.alive2)).and. useERA5) goto 1505

				if (.not.alive1) then 
					print*, "│  │  ├── ERA5-Plevel file Not Found, calling Python downloading script....."
					call system("python3 subs/era5/download_era5_profiles.py "//&
								yyyymmdd//" "//HH1//" "//trim(ERA5_DIR))
				ELSE
					PRINT*, "│  │  ├── Old ERA5-Plevels file found for HR1: ",HH1,":00. Extracting vars....."
				end if
				
				if (.not.alive2) then 
					print*, "│  │  ├── ERA5-Land file Not Found, calling Python downloading script....."	
					call system("python3 subs/era5/download_era5_lands.py "//&
								yyyymmdd//" "//HH1//" "//trim(ERA5_DIR))
				ELSE
					PRINT*, "│  │  ├── Old ERA5-Land file found for HR1: ",HH1,":00. Extracting vars....."					
				end if
								
				call read_ERA5_profiles(trim(req_era5_hr1),elon,elat,qw1,ta1,profIO)
				call read_ERA5_landvars(trim(req_land_hr1),lst1,psrf1,t2m1,snowc1,smc1,landIO)        !! K, hPa, K     !!  1 (top) -> 29 (bottom)
				
				IF (.NOT.(profIO.and.landIO)) goto 1505
				
				where(qw1.eq.0) qw1=1E-9			
			END IF
			!! ---------------------------------------------------------------------------------------------------------------------		
			!! Hour2
			req_era5_hr2=trim(adjustl(ERA5_DIR))//"/"//yyyymmdd//"/ERA5-PL-GBL-25km-"//yyyymmdd//"-"//HH2//"00.nc"
			req_land_hr2=trim(adjustl(ERA5_DIR))//"/"//yyyymmdd//"/ERA5-Land-GBL-25km-"//yyyymmdd//"-"//HH2//"00.nc"
			
			IF (trim(req_era5_hr2).NE.trim(ERA5_HR2_STAMP)) THEN
				ERA5_HR2_STAMP=req_era5_hr2
				PRINT*, "│  ├──New ERA5 file need for HR2: ",HH2,":00"
					
				inquire(file=trim(adjustl(req_era5_hr2)), exist=alive1)
				inquire(file=trim(adjustl(req_land_hr2)), exist=alive2)

				IF (((.not.alive1).or.(.not.alive2)).and. useERA5) goto 1505

				if (.not.alive1) then 
					print*, "│  │  ├── ERA5-Plevel file Not Found, calling Python downloading script....."
					call system("python3 subs/era5/download_era5_profiles.py "//&
								yyyymmdd//" "//HH2//" "//trim(ERA5_DIR))
				ELSE
					PRINT*, "│  │  ├── Old ERA5-Plevels file found for HR2: ",HH2,":00. Extracting vars....."
				end if
				if (.not.alive2) then 
					print*, "│  │  ├── ERA5-Land file Not Found, calling Python downloading script....."	
					call system("python3 subs/era5/download_era5_lands.py "//&
								yyyymmdd//" "//HH2//" "//trim(ERA5_DIR))
				ELSE
					PRINT*, "│  │  ├── Old ERA5-Land file found for HR2: ",HH2,":00. Extracting vars....."					
				end if
							
				call read_ERA5_profiles(trim(req_era5_hr2),elon,elat,qw2,ta2,profIO)
				call read_ERA5_landvars(trim(req_land_hr2),lst2,psrf2,t2m2,snowc1,smc2,landIO)

				IF (.NOT.(profIO.and.landIO)) goto 1505
				
				where(qw2.eq.0) qw2=1E-9	
				
			END IF
			!! -------------------------------------------------------------------------------------------------------------------------
			lon=Longitude(ipixel,iscan)
			lat=Latitude(ipixel,iscan)
			
			ilon=nint((lon-minval(elon))/egrid)+1
			ilat=nint((maxval(elat)-lat)/egrid)+1
			if (ilon>nlon) ilon=nlon
			if (ilat>nlat) ilat=nlat			
			if (ilon<1) ilon=1
			if (ilat<1) ilat=1	
		
			Vapor_land(:,iland)=qw1(ilon,ilat,:,1)*wt1+qw2(ilon,ilat,:,1)*wt2     ! kg kg-1	
			AirTemp_land(:,iland)=ta1(ilon,ilat,:,1)*wt1+ta2(ilon,ilat,:,1)*wt2     ! K
			
			
			LST_land(iland)=lst1(ilon,ilat,1)*wt1+lst2(ilon,ilat,1)*wt2     ! K
			Psrf_land(iland)=psrf1(ilon,ilat,1)*wt1+psrf2(ilon,ilat,1)*wt2  ! hpa
			
			T2m_land(iland)=t2m1(ilon,ilat,1)*wt1+t2m2(ilon,ilat,1)*wt2   ! K
			SnowC_land(iland)=snowc1(ilon,ilat,1)
			SMC_land(iland)=smc1(ilon,ilat,1)*wt1+smc2(ilon,ilat,1)*wt2   ! m3/m3
			END if
		end do
	end do
	
  ! --------------------------------------------------------------------------
  ! 4. ! Call Rttov Retrieve subruntn
  ! --------------------------------------------------------------------------

	read(yyyymmdd(5:6),'(I2)') imonth
	
	write(samples,"(i0.7)") num_land
	
	print*, "├──RTTOV retrieve begining, total samples= "//samples
	! channel 10.65 
	print*, "├──── Start RTTOV retrieve Channel 1-2 ......."

	call rttov_retrieve_mwri_emiss_clearsky(imonth,nlevel,num_land,2,53.1,(/1,2/),&
			Longitude_land,Latitude_land,PLEVEL,PHALF,&
			Vapor_land,AirTemp_land,&
			LST_land,Psrf_land,T2m_land,SnowC_land,&
			SMC_land,TB_land(1:2,:),Emissivity_land(1:2,:), FYSat) 

	! channel 18,18,23,36,36
	print*, "├──── Start RTTOV retrieve Channel 3-8 ......."
	call rttov_retrieve_mwri_emiss_clearsky(imonth,nlevel,num_land,6,53.1,(/3,4,5,6,7,8/),&
			Longitude_land,Latitude_land,PLEVEL,PHALF,&
			Vapor_land,AirTemp_land,&
			LST_land,Psrf_land,T2m_land,SnowC_land,&
			SMC_land,TB_land(3:8,:),Emissivity_land(3:8,:), FYSat) 

	! channel 89
	print*, "├──── Start RTTOV retrieve Channel 9-10 ......."
	call rttov_retrieve_mwri_emiss_clearsky(imonth,nlevel,num_land,2,53.1,(/9,10/),&
			Longitude_land,Latitude_land,PLEVEL,PHALF,&
			Vapor_land,AirTemp_land,&
			LST_land,Psrf_land,T2m_land,SnowC_land,&
			SMC_land,TB_land(9:10,:),Emissivity_land(9:10,:), FYSat) 

	!! post-process
	WHERE(Latitude_land.le.minval(elat).or.Latitude_land.ge.maxval(elat))
		Emissivity_land(1,:)=-999.9
		Emissivity_land(2,:)=-999.9
		Emissivity_land(3,:)=-999.9
		Emissivity_land(4,:)=-999.9
		Emissivity_land(5,:)=-999.9
		Emissivity_land(6,:)=-999.9
		Emissivity_land(7,:)=-999.9
		Emissivity_land(8,:)=-999.9
		Emissivity_land(9,:)=-999.9
		Emissivity_land(10,:)=-999.9
	ENDWHERE
		
	WHERE(LST_land.le.10 .or. LZA_land.lt.0)
		Emissivity_land(1,:)=-999.9
		Emissivity_land(2,:)=-999.9
		Emissivity_land(3,:)=-999.9
		Emissivity_land(4,:)=-999.9
		Emissivity_land(5,:)=-999.9
		Emissivity_land(6,:)=-999.9
		Emissivity_land(7,:)=-999.9
		Emissivity_land(8,:)=-999.9
		Emissivity_land(9,:)=-999.9
		Emissivity_land(10,:)=-999.9
		LST_land=-999.9
	ENDWHERE

	WHERE(Emissivity_land.eq.-1 ) Emissivity_land=-999.9

	!! reorgnize for HDF5 output
	IF (HDF5) THEN
		ALLOCATE(Emissivity_swath(nchannel,npixel,nscan))
		Allocate(LST_swath(npixel,nscan),T2m_swath(npixel,nscan))
		Allocate(SnowC_swath(npixel,nscan),SMC_swath(npixel,nscan))
			LST_swath  = -999.9
			T2m_swath  = -999.9
			SnowC_swath= -999.9
			SMC_swath  = -999.9
			Emissivity_swath=-999.9	
			
		DO iland=1,num_land
			iscan=scanid(iland)
			ipixel=pixelid(iland)
			Emissivity_swath(:,ipixel,iscan)=Emissivity_land(:,iland)
			LST_swath(ipixel,iscan)  = LST_land(iland)
			T2m_swath(ipixel,iscan)  = T2m_land(iland)
			SnowC_swath(ipixel,iscan)= SnowC_land(iland)
			SMC_swath(ipixel,iscan)	 = SMC_land(iland)
		END DO



		WHERE(Miss_swath(1,:,:).ge.99)
			Emissivity_swath(1,:,:)=-999.9
			Emissivity_swath(2,:,:)=-999.9
		ENDWHERE
		
		WHERE(Miss_swath(2,:,:).ge.99)	
			Emissivity_swath(3,:,:)=-999.9
			Emissivity_swath(4,:,:)=-999.9
			Emissivity_swath(5,:,:)=-999.9
			Emissivity_swath(6,:,:)=-999.9
			Emissivity_swath(7,:,:)=-999.9
			Emissivity_swath(8,:,:)=-999.9
		ENDWHERE
		
		WHERE(Miss_swath(3,:,:).ge.99)
			Emissivity_swath(9,:,:)=-999.9
			Emissivity_swath(10,:,:)=-999.9
		ENDWHERE
	END IF	
  ! --------------------------------------------------------------------------
  ! 5. output Retrieve result
  ! --------------------------------------------------------------------------
	! output to ascci text
	
	WHERE(Miss_land(1,:).ge.99)
		Emissivity_land(1,:)=-999.9
		Emissivity_land(2,:)=-999.9
	ENDWHERE
	
	WHERE(Miss_land(2,:).ge.99)	
		Emissivity_land(3,:)=-999.9
		Emissivity_land(4,:)=-999.9
		Emissivity_land(5,:)=-999.9
		Emissivity_land(6,:)=-999.9
		Emissivity_land(7,:)=-999.9
		Emissivity_land(8,:)=-999.9
	ENDWHERE
	
	WHERE(Miss_land(3,:).ge.99)
		Emissivity_land(10,:)=-999.9
		Emissivity_land(9,:)=-999.9
	ENDWHERE
	
	open(NEWUNIT=text_unit,file=trim(adjustL(EMISS_TEXT)))
	DO iland=1,num_land
		
		iscan=scanid(iland)
		ipixel=pixelid(iland)

		write(text_unit,505) nscan,npixel,iscan,ipixel,Longitude_land(iland),Latitude_land(iland),&
			ScanTime_land(iland),Emissivity_land(:,iland),TB_land(:,iland),&
			LST_land(iland),T2m_land(iland),SnowC_land(iland),SMC_land(iland), &
			LZA_land(iland),CFR_land(:,iland),Clear_land(:,iland),Miss_land(:,iland),&
			Sea_land(iland),Cloudflg_land(iland)
	
	END DO	
	PRINT*,"└── Output retrieval: "//trim(adjustL(EMISS_TEXT))
505 format(4I6, 3f10.4, 20f10.4, 5f10.4, 10f10.4, I6)
	
	close(text_unit)
	
	IF (HDF5) THEN
		CALL write_emiss_hdf5(trim(adjustL(EMISS_HDF5)),nchannel,npixel,nscan, &
				Longitude,Latitude,ScanTime_swath,TB_swath,Emissivity_swath,&
				LST_swath,T2m_swath,SnowC_swath,SMC_swath, &
				LZA_swath,CFR_swath,Clear_swath,&
				Sea_Frac,Cloud_Flag,Miss_swath,status)
		
		IF(status/=0) THEN
			PRINT*,"└── Write emissivity HDF5 error : ",trim(adjustL(EMISS_HDF5))
			Stop
		ENDIF		 
		PRINT*,"└── Output retrieval: "//trim(adjustL(EMISS_HDF5))
	END IF	
  ! --------------------------------------------------------------------------
  ! 6. free memory
  ! --------------------------------------------------------------------------	

	IF (HDF5) THEN
		DEALLOCATE(LST_swath, T2m_swath, SnowC_swath, SMC_swath)	
		DEALLOCATE(Emissivity_swath)
	END IF	

1505 continue	
	DEALLOCATE(pixelid,scanid)
	
	DEALLOCATE(LST_land, Psrf_land, T2m_land, SnowC_land, SMC_land)	
	DEALLOCATE(Vapor_land, AirTemp_land)

	DEALLOCATE(Emissivity_land)
	DEALLOCATE(TB_land)
	DEALLOCATE(Longitude_land,Latitude_land)

	DEALLOCATE(ScanTime_land,LZA_land,Sea_land,Cloudflg_land)
	DEALLOCATE(CFR_land,Clear_land,Miss_land)

2333 Continue
		! print*,"pass 2333"	
	!! Deallocate( FY3)
	DEALLOCATE(TB_swath, Latitude, Longitude, ScanTime_swath, ScanTime, LZA_swath)	
	DEALLOCATE(CFR_swath,Clear_swath,Sea_Frac,Cloud_Flag,Miss_swath)


	PRINT*,"│--> Go to next orbit ......"
	GOTO 176
	
886 continue

	PRINT*, "|| "//FYSat//"_L1C process done - "//yyyymmdd	
	! Close the filelist
	CLOSE(file_unit)
	
END PROGRAM main_clear_retrieve_landonly





