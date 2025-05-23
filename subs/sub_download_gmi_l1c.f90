SUBROUTINE download_GMI_L1C(date_str,dir)
    CHARACTER(8), INTENT(IN) :: date_str
    CHARACTER(*), INTENT(IN) :: dir
  
    ! Internal variables
    INTEGER ::  status
    CHARACTER(3) :: doy
    CHARACTER(256) :: command
  
    ! Parse the date string to extract year, month, and day
    CALL ParseDate(date_str,doy)
  
    ! Specify the wget command
    command = 'wget -q -c -r -N -nH -nd -np'// &
	' --load-cookies ~/.urs_cookies --save-cookies ~/.urs_cookies --keep-session-cookies'// &
	' -P '//trim(adjustl(dir))//date_str(1:4)//'/'//date_str(5:6)//date_str(7:8)//&
	' -A HDF5 https://gpm1.gesdisc.eosdis.nasa.gov/data/GPM_L1C/GPM_1CGPMGMI.07/'// date_str(1:4)// '/' // doy // '/'
  
    ! Execute the wget command
    CALL SYSTEM(command, status)
  
    ! Check the exit status
    IF (status /= 0) THEN
      WRITE(*,*) 'Error: wget command failed with status', status
	  STOP
    ELSE
      WRITE(*,*) 'GMI_L1C HDF5 Download successful'
    END IF
  

  
END SUBROUTINE download_GMI_L1C

