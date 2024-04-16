
!! sub to download GOES-R Cloud from AWS s3 data bucket 

!!  >  curl "https://awscli.amazonaws.com/awscli-exe-linux-x86_64.zip" -o "awscliv2.zip"
!!  >  unzip awscliv2.zip
!!  >  sudo ./aws/install

SUBROUTINE download_GOESR_AWS(yyyymmdd,HH,MM,GEOS_L2_DIR)

	CHARACTER(*), INTENT(IN) :: GEOS_L2_DIR
	CHARACTER(8), INTENT(IN) :: yyyymmdd
	CHARACTER(2), INTENT(IN) :: HH,MM
  
    ! Internal variables
    INTEGER ::  status,file_unit
    CHARACTER(3) :: DOY
    CHARACTER(256) :: command
  
    ! Parse the date string to extract year, month, and day
    CALL ParseDate(yyyymmdd,DOY)
  
  	IF(yyyymmdd.LT.'20200101') THEN
		PRINT*,'ERROR! Date before 20200101 not exist in AWS cloud bucket!'
		STOP
	END	IF
  
!! ACM Clearsky MASK 
    command = "aws s3 cp --recursive --no-sign-request --exclude '*' &
	         --include '*s"//yyyymmdd(1:4)//DOY//HH//MM//"???_e*_c*.nc' &
			 s3://noaa-goes16/ABI-L2-ACMF/"//yyyymmdd(1:4)//"/"//DOY//"/"//HH//"/ "&
			 //trim(GEOS_L2_DIR)//"/"//yyyymmdd//"/"
    CALL SYSTEM(command, status)
  
!! COD Cloud Optical Depeth
	
    command = "aws s3 cp --recursive --no-sign-request --exclude '*' &
	         --include '*s"//yyyymmdd(1:4)//DOY//HH//MM//"???_e*_c*.nc' &
			 s3://noaa-goes16/ABI-L2-CODF/"//yyyymmdd(1:4)//"/"//DOY//"/"//HH//"/ &
			 "//trim(GEOS_L2_DIR)//"/"//yyyymmdd//"/"
    CALL SYSTEM(command, status) 

!! CPS cloud partical size
	
    command = "aws s3 cp --recursive --no-sign-request --exclude '*' &
	         --include '*s"//yyyymmdd(1:4)//DOY//HH//MM//"???_e*_c*.nc' &
			 s3://noaa-goes16/ABI-L2-CPSF/"//yyyymmdd(1:4)//"/"//DOY//"/"//HH//"/ &
			 "//trim(GEOS_L2_DIR)//"/"//yyyymmdd//"/"
    CALL SYSTEM(command, status) 

!!  ACHA cloud top height

	
    command = "aws s3 cp --recursive --no-sign-request --exclude '*' &
	         --include '*s"//yyyymmdd(1:4)//DOY//HH//MM//"???_e*_c*.nc' &
			 s3://noaa-goes16/ABI-L2-ACHAF/"//yyyymmdd(1:4)//"/"//DOY//"/"//HH//"/ &
			 "//trim(GEOS_L2_DIR)//"/"//yyyymmdd//"/"
    CALL SYSTEM(command, status) 
  
!! ABI-L2-ACTPF  Cloud top phase
  
     command = "aws s3 cp --recursive --no-sign-request --exclude '*' &
	         --include '*s"//yyyymmdd(1:4)//DOY//HH//MM//"???_e*_c*.nc' &
			 s3://noaa-goes16/ABI-L2-ACTPF/"//yyyymmdd(1:4)//"/"//DOY//"/"//HH//"/ &
			 "//trim(GEOS_L2_DIR)//"/"//yyyymmdd//"/"
    CALL SYSTEM(command, status)  


    ! Check the exit status
    IF (status /= 0) THEN
      WRITE(*,*) 'Error: wget command failed with status', status
	  STOP
    ELSE
      WRITE(*,*) 'GOES-R NetCDF file Download successful'
    END IF  

	!! P.S. YOU CAN NOT IMAGINE HOW HARD IT IS TO DOWNLOAD AND USE GOES DATA UNTIL U DO IT,
	!! THE MOST HARSH THING TODO IS TO OPERATE SUCH FILES WITH TOTALLY ARBITARY NAMES......
	!! YOU CAN NOT USE WGET TO DOWNLOAD IT DIRECTLY, 
	!! NOR CAN YOU EASILY GET THESE NAMES UNLESS YOU STRUGGLE TO 
	!! CALL SYSTEM COMMAND TO GET THEIR SIGNATURES WHICH THEY ARE DIFFERENT FROM EACH OTHER
	!! AND THEN RETRIEVE THE FILENAMES, YOU GANNA GOT A LONG SHIT CODE IN THIS BLOCK.

	! CALL system("cd "//trim(GEOS_L2_DIR)//"/"//yyyymmdd//&
		! " ; ls OR_ABI-L2-ACMF-M6_G16_s"//yyyymmdd(1:4)//DOY//HH//MM//"???_e*_c*.nc &
		! >"//trim(PWD)//"/filelists/Name_ACMF_G16_s"//yyyymmdd//HH//MM//".txt")
	! OPEN(NEWUNIT=file_unit, FILE=trim(PWD)//"/filelists/Name_ACMF_G16_s"//yyyymmdd//HH//MM//".txt")
	! read(file_unit,'(A)') GOES_Sufix
	! close(file_unit)
	! print*,GOES_Sufix
	! GOES_Sufix=GOES_Sufix(23:72)  !! and length differs from variables 
	! print*,GOES_Sufix  
	
	!! INSTEAD, I Rename THEM ALL.

	CALL system("cd "//trim(GEOS_L2_DIR)//"/"//yyyymmdd//&
		" ; mv -v OR_ABI-L2-ACMF-M6_G16_s"//yyyymmdd(1:4)//DOY//HH//MM//"???_e*_c*.nc &
		OR_ABI-L2-ACMF-M6_G16_"//yyyymmdd//"_"//HH//MM//".NC")
		
	CALL system("cd "//trim(GEOS_L2_DIR)//"/"//yyyymmdd//&
		" ; mv -v OR_ABI-L2-CODF-M6_G16_s"//yyyymmdd(1:4)//DOY//HH//MM//"???_e*_c*.nc &
		OR_ABI-L2-CODF-M6_G16_"//yyyymmdd//"_"//HH//MM//".NC")	
		
	CALL system("cd "//trim(GEOS_L2_DIR)//"/"//yyyymmdd//&
		" ; mv -v OR_ABI-L2-CPSF-M6_G16_s"//yyyymmdd(1:4)//DOY//HH//MM//"???_e*_c*.nc &
		OR_ABI-L2-CPSF-M6_G16_"//yyyymmdd//"_"//HH//MM//".NC")	
		
	CALL system("cd "//trim(GEOS_L2_DIR)//"/"//yyyymmdd//&
		" ; mv -v OR_ABI-L2-ACHAF-M6_G16_s"//yyyymmdd(1:4)//DOY//HH//MM//"???_e*_c*.nc &
		OR_ABI-L2-ACHAF-M6_G16_"//yyyymmdd//"_"//HH//MM//".NC")		
		
	CALL system("cd "//trim(GEOS_L2_DIR)//"/"//yyyymmdd//&
		" ; mv -v OR_ABI-L2-ACTPF-M6_G16_s"//yyyymmdd(1:4)//DOY//HH//MM//"???_e*_c*.nc &
		OR_ABI-L2-ACTPF-M6_G16_"//yyyymmdd//"_"//HH//MM//".NC")


END SUBROUTINE download_GOESR_AWS



!! sub to download GOES-R Cloud from AWS s3 data bucket 

!!  >  curl "https://awscli.amazonaws.com/awscli-exe-linux-x86_64.zip" -o "awscliv2.zip"
!!  >  unzip awscliv2.zip
!!  >  sudo ./aws/install

SUBROUTINE download_GOESR_AWS_all(yyyymmdd,HH,GEOS_L2_DIR)

	CHARACTER(*), INTENT(IN) :: GEOS_L2_DIR
	CHARACTER(8), INTENT(IN) :: yyyymmdd
	CHARACTER(2), INTENT(IN) :: HH!,MM
  
    ! Internal variables
    INTEGER ::  status,file_unit,ih,im
    CHARACTER(3) :: DOY
    CHARACTER(256) :: command
	CHARACTER(2) :: HH0,MM0
    ! Parse the date string to extract year, month, and day
    CALL ParseDate(yyyymmdd,DOY)
  
  	IF(yyyymmdd.LT.'20200101') THEN
		PRINT*,'ERROR! Date before 20200101 not exist in AWS cloud bucket!'
		STOP
	END	IF
  
!! ACM Clearsky MASK 
    command = "aws s3 cp --quiet --recursive --no-sign-request &
			 s3://noaa-goes16/ABI-L2-ACMF/"//yyyymmdd(1:4)//"/"//DOY//"/ "&
			 //trim(GEOS_L2_DIR)//"/"//yyyymmdd(1:4)//"/"//yyyymmdd//"/"
    CALL SYSTEM(command, status)
  
    ! Check the exit status
    IF (status /= 0) THEN
      WRITE(*,*) 'Error: wget command failed with status', status
	  STOP
    ELSE
      WRITE(*,*) '│  │  ├── GOES-R NetCDF files Download successful, Renaming......'
    END IF  

	do ih=0,23
		do im=0,50,10		
		WRite(MM0,"(i0.2)") im
		WRite(HH0,"(i0.2)") ih
	CALL system("cd "//trim(GEOS_L2_DIR)//"/"//yyyymmdd(1:4)//"/"//yyyymmdd//"/"//HH0//&
		" ; mv OR_ABI-L2-ACMF-M6_G16_s"//yyyymmdd(1:4)//DOY//HH0//MM0//"???_e*_c*.nc &
		../OR_ABI-L2-ACMF-M6_G16_"//yyyymmdd//"_"//HH0//MM0//".NC")
			
		
	end do
	end do

END SUBROUTINE download_GOESR_AWS_all
