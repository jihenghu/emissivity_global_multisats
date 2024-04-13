subroutine download_AHI_L2CLP(yyyymmdd,hh,ahi_name,ahi_L2clp_dir,status)
    CHARACTER(8), INTENT(IN) :: yyyymmdd
    INTEGER, INTENT(OUT) :: status
    CHARACTER(2), INTENT(IN) :: hh
    CHARACTER(*), INTENT(IN) :: ahi_name,ahi_L2clp_dir
    CHARACTER(255) :: command
	

    ! Specify the wget command
    command='wget -q -c -r -nH -nd -np -P '//trim(ahi_L2clp_dir)//yyyymmdd//'/'//hh// &
	 ' -A nc ftp://hjh18305_gmail.com:SP+wari8@ftp.ptree.jaxa.jp/pub/himawari/L2/CLP/010/' &
	 //yyyymmdd(1:6)//'/'//yyyymmdd(7:8)//'/'//hh//'/'//trim(ahi_name)
   
	! print*,command
	! stop
    ! Execute the wget command
    CALL SYSTEM(command, status)
  
    ! Check the exit status
    ! IF (status /= 0) THEN
      ! WRITE(*,*) 'Error: Himawari-8 NETCDF Download with status', status
    ! ELSE
      ! WRITE(*,*) 'Himawari-8 NETCDF Download successful'
    ! END IF

end subroutine download_AHI_L2CLP