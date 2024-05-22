! module fysubs   
    ! contains
! read_fy3bcd_Vars(FY3_FILENAME, TB_swath, Latitude, Longitude, ScanTime, nscan, npixel, nchannel,FYIO)
subroutine read_fy3bcd_Vars(filename,tbs,lat,lon,FYstime,nscan,npix,nchn,error)

	USE HDF5   
	CHARACTER*(*) :: filename

	INTEGER :: nscan,npix, nchn

	INTEGER(HID_T) :: file_id       ! File identifier  文件句柄
	INTEGER(HID_T) :: dset_id       ! Dataset identifier 变量句柄
	INTEGER(HID_T) :: grp_id        ! Dataset identifier  group句柄
	INTEGER(HID_T) :: dspace_id     ! Dataset identifier  只可意会，和维数和大小相关
	INTEGER     	 :: error           ! Error flag - success：0

	INTEGER:: tb(npix, nscan, nchn) !,dem(npix, nscan)
	Real*4 :: tbs(npix, nscan, nchn)
	Real*4 :: lon(npix, nscan),lat(npix, nscan)

	Real*8 :: minsec(2,nscan)
	Real*4 :: FYstime(nscan)

	! INTEGER  :: igbp(npix, nscan)   
	! INTEGER  :: lsmask(npix, nscan) 

	INTEGER(HSIZE_T), DIMENSION(3) :: dims_3D,maxdims_3D 
	INTEGER(HSIZE_T), DIMENSION(2) :: dims_2D,maxdims_2D
	INTEGER(HSIZE_T), DIMENSION(1) :: dims_1D,maxdims_1D
  
	! print*,trim(adjustl(filename))
	CALL h5open_f(error) 
	CALL h5fopen_f (trim(adjustl(filename)), H5F_ACC_RDONLY_F, file_id, error) 
			if (error.ne.0) then
			  return
			end if

	CALL h5gopen_f (file_id, "Calibration", grp_id, error)  
		!!!!! EARTH_OBSERVE_BT_10_to_89GHz
		CALL h5dopen_f(grp_id, "EARTH_OBSERVE_BT_10_to_89GHz", dset_id, error)  
		CALL h5dget_space_f(dset_id,dspace_id,error)  
		CALL h5sget_simple_extent_dims_f(dspace_id, dims_3D, maxdims_3D, error)   
		CALL h5dread_f(dset_id, H5T_NATIVE_INTEGER, tb, dims_3D, error)
		CALL h5dclose_f(dset_id, error) 
		tbs=tb*0.01+327.68
		where(tbs.gt.600 .or. tbs.lt.30) tbs=-999.9  !! Kelvin
		
		CALL h5dopen_f(grp_id, "Scan_Mscnt", dset_id, error)  
		CALL h5dget_space_f(dset_id,dspace_id,error)  
		CALL h5sget_simple_extent_dims_f(dspace_id, dims_2D, maxdims_2D, error)   
		CALL h5dread_f(dset_id, H5T_IEEE_F64LE, minsec, dims_2D, error)
		CALL h5dclose_f(dset_id, error) 
		! print*,minsec(1,:5)
		minsec=minsec*1.0  !! SLOPE: 1
		FYstime=real(minsec(1,:)/1000./60./60.)+12.  
		where(FYstime.gt.24) FYstime=FYstime-24            !! UTC	
	CALL h5gclose_f(grp_id, error) 

	CALL h5gopen_f (file_id, "Geolocation", grp_id, error)  
		CALL h5dopen_f(grp_id, "Latitude", dset_id, error)  
		CALL h5dget_space_f(dset_id,dspace_id,error)  
		CALL h5sget_simple_extent_dims_f(dspace_id, dims_2D, maxdims_2D, error)   
		CALL h5dread_f(dset_id, H5T_NATIVE_REAL, lat, dims_2D, error)
		CALL h5dclose_f(dset_id, error) 

		CALL h5dopen_f(grp_id, "Longitude", dset_id, error)  
		CALL h5dget_space_f(dset_id,dspace_id,error)  
		CALL h5sget_simple_extent_dims_f(dspace_id, dims_2D, maxdims_2D, error)   
		CALL h5dread_f(dset_id, H5T_NATIVE_REAL, lon, dims_2D, error)
		CALL h5dclose_f(dset_id, error) 		
	CALL h5gclose_f(grp_id, error) 


	CALL h5fclose_f(file_id, error)  	
    CALL h5close_f(error)				
end subroutine
! end module fysubs


subroutine get_fy3bcd_swath_geometry(filein,nscan,npix,nchanl,error)
    USE HDF5 
    implicit none 
    character*(*) :: filein
    character(len=200)::varname
    integer:: error
    integer(HID_T) :: file_id,dset_id,dspace_id
    integer(HSIZE_T) :: data_dims3(3),maxdims3(3)
    integer :: nscan ,npix, nchanl
    
    call h5open_f(error)
    call h5fopen_f(filein, H5F_ACC_RDONLY_F, file_id, error)
    if (error.ne.0) then
      return
    end if
    
    varname = "Calibration/EARTH_OBSERVE_BT_10_to_89GHz"
    CALL h5dopen_f(file_id, trim(varname), dset_id, error)
    if (error.ne.0) then
      return
    end if
    CALL h5dget_space_f(dset_id,dspace_id,error)
    CALL h5sget_simple_extent_dims_f(dspace_id, data_dims3, maxdims3, error)
    nscan = data_dims3(2) 
    npix = data_dims3(1) 
    nchanl = data_dims3(3) 
    CALL h5sclose_f(dspace_id,error)
    CALL h5dclose_f(dset_id,error)
    
    
    CALL h5fclose_f(file_id,error)
    CALL h5close_f(error)
end subroutine 