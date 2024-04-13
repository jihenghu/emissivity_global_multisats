subroutine read_land_sea_mask(filename,LSM)
	include 'netcdf.inc'
	
	character*(*) :: filename
	REAL*4, DIMENSION(3602,1800) :: LSM
	
	Integer status,ncid,varid

 
    status=nf_open(trim(adjustl(filename)),nf_nowrite,ncid) !打开netcdf文件，获取文件的ID号(ncid)      
    if (status .ne. 0) then
        print*,"Open failure, sub _read_land_sea_mask"
        stop
    end if

    status=nf_inq_varid (ncid, 'landseamask', varid) 
    status=nf_get_var(ncid,varid,LSM)

	status=nf_close(ncid)
END SUBROUTINE read_land_sea_mask