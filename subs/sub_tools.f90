

  ! Set the input string and separator
  ! inputString = 'part1.part2.part3.part4'
  ! separator = '.'

  ! Call the subroutine to extract the desired piece
  ! CALL GetDesiredPiece(inputString, separator, 3, desiredPiece)

  ! Print the result
  ! PRINT *, 'Input String:', inputString
  ! PRINT *, 'Desired Piece:', desiredPiece

SUBROUTINE GetDesiredPiece(inputString, separator, pieceIndex, desiredPiece)
  CHARACTER(LEN=*), INTENT(IN) :: inputString, separator
  INTEGER, INTENT(IN) :: pieceIndex
  CHARACTER(LEN=*), INTENT(OUT) :: desiredPiece

  INTEGER :: startPos, endPos
  INTEGER :: separatorCount, i

  ! Initialize variables
  separatorCount = 0

  startPos = 1

  ! Loop through the string to find the desired piece
  DO i = 1, LEN_TRIM(inputString)
    IF (inputString(i:i) == separator) THEN
      separatorCount = separatorCount + 1

      IF (separatorCount == pieceIndex - 1) THEN
        startPos = i + 1
      ELSE IF (separatorCount == pieceIndex) THEN
        endPos = i - 1
        EXIT
      END IF
    END IF
  END DO

  IF (separatorCount >= pieceIndex) THEN
    ! Extract the desired piece
    desiredPiece = TRIM(inputString(startPos:endPos))
  ELSE
    ! Piece not found
    desiredPiece = 'Not Found'
  END IF

END SUBROUTINE GetDesiredPiece


! Convert binary to array, and pick 7th and 6th bit
subroutine getBit_7_6(temp,b7,b6) 
	integer temp,decimals,k,quotient
	integer binarys(16),b7,b6
	decimals=temp
	if (decimals .gt.65535) then
		print*,'[Fatal Error] input exceeds maxium decimal 65535!'
		binarys=-1
		return
	end if
	k=1
	binarys=0
	do while(decimals.ne.0)
			quotient=decimals/2
			binarys(k)=decimals-quotient*2
			k=k+1
			decimals=decimals/2	
	end do
	b7=binarys(7)
	b6=binarys(6)
	! return
	! stop
end subroutine getBit_7_6


! Convert binary to array, and pick 12th bit
! (13) SOZ>80 or SAZ>70: 0=Yes, 1=No;
subroutine getBit_13(temp,b13) 
	integer temp,decimals,k,quotient
	integer binarys(16),b13
	decimals=temp
	if (decimals .gt.65535) then
		print*,'[Fatal Error] input exceeds maxium decimal 65535!'
		binarys=-1
		return
	end if
	k=1
	binarys=0
	do while(decimals.ne.0)
			quotient=decimals/2
			binarys(k)=decimals-quotient*2
			k=k+1
			decimals=decimals/2	
	end do
	b13=binarys(13)

	! return
	! stop
end subroutine getBit_13


!! convert yyyymmdd to day of year
SUBROUTINE ParseDate(date_str, doy)
  CHARACTER(8), INTENT(IN) :: date_str
  CHARACTER(3), INTENT(OUT) :: doy  
  
  INTEGER :: year, month, day, doyint, iadd
  
  ! Extract year, month, and day from the date string
  READ(date_str(1:4),*) year
  READ(date_str(5:6),*) month
  READ(date_str(7:8),*) day
  

  ! Check for leap year
  IF ((MOD(year, 4) == 0 .AND. MOD(year, 100) /= 0) .OR. MOD(year, 400) == 0) THEN
	iadd=1
	IF (month == 2 .AND. (day < 1 .OR. day > 29)) THEN
	  WRITE(*,*) 'Error: Day of month for February in a leap year should be between 1 and 29.'
	  STOP
	ENDIF
  ELSE
	iadd=0
	IF (month == 2 .AND. (day < 1 .OR. day > 28)) THEN
	  WRITE(*,*) 'Error: Day of month for February in a non-leap year should be between 1 and 28.'
	  STOP
	ENDIF
  ENDIF

  ! Validate month and day
  IF (month < 1 .OR. month > 12 .OR. day < 1 .OR. day > 31) THEN
	WRITE(*,*) 'Error: Invalid month or day.'
	STOP
  ENDIF
  

  SELECT CASE (month)
	CASE (1)  ! January
	  doyint = day
	CASE (2)  ! February
	  doyint = day + 31
	CASE (3)  ! March
	  doyint = day + 31 + 28 + iadd  ! Add 1 day for leap year
	CASE (4)  ! April
	  doyint = day + 31 + 28 + iadd + 31
	CASE (5)  ! May
	  doyint = day + 31 + 28 + iadd + 31 + 30
	CASE (6)  ! June
	  doyint = day + 31 + 28 + iadd + 31 + 30 + 31
	CASE (7)  ! July
	  doyint = day + 31 + 28 + iadd + 31 + 30 + 31 + 30
	CASE (8)  ! August
	  doyint = day + 31 + 28 + iadd + 31 + 30 + 31 + 30 + 31
	CASE (9)  ! September
	  doyint = day + 31 + 28 + iadd + 31 + 30 + 31 + 30 + 31 + 31
	CASE (10) ! October
	  doyint = day + 31 + 28 + iadd + 31 + 30 + 31 + 30 + 31 + 31 + 30
	CASE (11) ! November
	  doyint = day + 31 + 28 + iadd + 31 + 30 + 31 + 30 + 31 + 31 + 30 + 31
	CASE DEFAULT ! December
	  doyint = day + 31 + 28 + iadd + 31 + 30 + 31 + 30 + 31 + 31 + 30 + 31 + 30
  END SELECT
  
  WRITE(doy,'(i0.3)') doyint
  
END SUBROUTINE ParseDate
  
  
!! https://www.ncei.noaa.gov/sites/default/files/2021-07/GOES-R_ABI_local_zenith_angle_description_0.docx  
SUBROUTINE calc_LZA(lon,lat,lon0,lat0,LZA)

	Real, intent(IN):: lon, lat, lon0,lat0  !! degrees
	REAL, intent(OUT) :: LZA             !! degrees

    REAL,parameter :: H=42164.16 , req = 6378.137 !! km, 
	REAL :: beta, delta,lonrad,latrad , lon0rad,lat0rad
	
	latrad=lat*3.14159265/180.
	lonrad=lon*3.14159265/180.
		
	lon0rad=lon0*3.14159265/180.
	lat0rad=lat0*3.14159265/180.
	
	beta=acos(cos(latrad-lat0rad)*cos(lonrad-lon0rad))   !! degrees
	
	delta=H*sin(beta)/sqrt(H*H+req*req-2*H*req*cos(beta))
	
	IF (delta.le.1 .and. delta.ge.-1) THEN   
		LZA=asin(delta)*180/3.14159265
	ELSE 
		LZA=-999
	END IF

END SUBROUTINE calc_LZA


subroutine calculate_cwc_specific(Temperature, Pressure, VaporMixingRatio, LWC, LWC_specific)
  ! Inputs:
  real, intent(in) :: Temperature   ! Temperature in Kelvin
  real, intent(in) :: Pressure      ! Pressure in Pascals
  real, intent(in) :: VaporMixingRatio ! Vapor Mixing Ratio in kg/kg
  real, intent(in) :: LWC           ! Liquid Water Content in kg/m³
  
  ! Outputs:
  real, intent(out) :: LWC_specific ! Specific Liquid Water Content in kg/kg
  
  ! Constants
  real, parameter :: R_air = 287.05    ! Specific gas constant for dry air in J/(kg·K)
  real, parameter :: R_vapor = 461.51  ! Specific gas constant for water vapor in J/(kg·K)
  real, parameter :: M_air = 28.97     ! Molar mass of dry air in g/mol

  ! Local variables
  real:: Density_air, Density_vapor, Density_moist_air

  ! Calculate density of dry air
  Density_air = Pressure / (R_air * Temperature) * (M_air / 1000.0)

  ! Calculate density of water vapor
  Density_vapor = VaporMixingRatio * Pressure / (R_vapor * Temperature)

  ! Calculate total density of moist air
  Density_moist_air = Density_air + Density_vapor

  ! Calculate specific liquid water content
  LWC_specific = LWC / Density_moist_air

end subroutine calculate_cwc_specific

