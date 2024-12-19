! program test_read_amsr
  ! character*255 filename
  ! character*2 grp
  
  ! integer status
  ! integer nscan,npixel,nchannel
  
  ! real, dimension(:,:), allocatable :: Latitude, Longitude, ScanTime_swath
  ! real, dimension(:,:,:), allocatable :: TB_swath
  ! real, dimension(:), allocatable :: ScanTime
  
  ! filename="/nfs/nuke/jihenghu/GPM_L1C/amsr2/2014/0101/1C.GCOMW1.AMSR2.XCAL2016-V.20140101-S004048-E021940.008643.V07A.HDF5"
  ! grp="S5"
  ! call get_AMSR2_primary_dims(filename,grp,nscan,npixel,nchannel,status)
  ! PRINT*, nscan,npixel,nchannel

  ! ALLOCATE(Latitude(npixel,nscan), Longitude(npixel,nscan))	
  ! ALLOCATE(TB_swath(nchannel,npixel,nscan), ScanTime_swath(npixel,nscan), ScanTime(nscan))
  ! CALL read_AMSR2_Vars(filename,grp, TB_swath, Latitude, Longitude, ScanTime, nscan, npixel, nchannel,status)
  ! print*,TB_swath(1,:10,1)
  ! print*,Latitude(:10,1)
  ! DEALLOCATE(Latitude, Longitude, ScanTime_swath)	
  ! DEALLOCATE(TB_swath, ScanTime) 
  
! end program test_read_amsr


subroutine read_AMSR2_Vars(filename,GRP, Tc, Latitude, Longitude, &
						scani_hr, nscan, npixel, nchannel,status)
  USE HDF5

  CHARACTER(LEN=*), INTENT(IN) :: filename
  CHARACTER*2 :: GRP
  REAL*4, DIMENSION(nchannel, npixel, nscan), INTENT(OUT) :: Tc
  REAL*4, DIMENSION(npixel, nscan), INTENT(OUT) :: Latitude, Longitude
  REAL*4, DIMENSION(nscan), INTENT(OUT) :: scani_hr
  INTEGER, INTENT(in) :: nscan, npixel, nchannel

  INTEGER(HID_T) ::file_id,dset_id,grp_id,dspace_id, status
  INTEGER(HSIZE_T), DIMENSION(3) :: dims_3D,maxdims_3D 
  INTEGER(HSIZE_T), DIMENSION(2) :: dims_2D,maxdims_2D
  INTEGER(HSIZE_T), DIMENSION(1) :: dims_1D,maxdims_1D
  ! tmp second
  REAL*8, DIMENSION(nscan) :: SecondOfDay

	CALL h5open_f(status) 
	CALL h5fopen_f (trim(adjustl(filename)), H5F_ACC_RDONLY_F, file_id, status) 
    if (status.ne.0) then
		WRITE(*,*) 'Error: Open hdf5 failed _sub_read_AMSR2_Vars'
		return
    end if
	!! group GRP 
	CALL h5gopen_f (file_id, GRP, grp_id, status)  
		!!!!! EARTH_OBSERVE_BT_10_to_89GHz
		CALL h5dopen_f(grp_id, "Tc", dset_id, status)  
		CALL h5dget_space_f(dset_id,dspace_id,status)  
		CALL h5sget_simple_extent_dims_f(dspace_id, dims_3D, maxdims_3D, status)   
		CALL h5dread_f(dset_id, H5T_NATIVE_REAL, Tc, dims_3D, status)
		CALL h5dclose_f(dset_id, status) 

		CALL h5dopen_f(grp_id, "Latitude", dset_id, status)  
		CALL h5dget_space_f(dset_id,dspace_id,status)  
		CALL h5sget_simple_extent_dims_f(dspace_id, dims_2D, maxdims_2D, status)   
		CALL h5dread_f(dset_id, H5T_NATIVE_REAL, Latitude, dims_2D, status)
		CALL h5dclose_f(dset_id, status) 

		CALL h5dopen_f(grp_id, "Longitude", dset_id, status)  
		CALL h5dread_f(dset_id, H5T_NATIVE_REAL, Longitude, dims_2D, status)
		CALL h5dclose_f(dset_id, status) 			
	CALL h5gclose_f(grp_id, status) 

	!! group  
	CALL h5gopen_f (file_id, GRP//"/ScanTime", grp_id, status)  
		CALL h5dopen_f(grp_id, "SecondOfDay", dset_id, status)  
		CALL h5dget_space_f(dset_id,dspace_id,status)  
		CALL h5sget_simple_extent_dims_f(dspace_id, dims_1D, maxdims_1D, status)   
		CALL h5dread_f(dset_id, H5T_NATIVE_DOUBLE, SecondOfDay, dims_1D, status)
		CALL h5dclose_f(dset_id, status) 
	CALL h5gclose_f(grp_id, status) 
	
	CALL h5fclose_f(file_id, status)  	
    CALL h5close_f(status)	

	! transfer second of day into utc hour
	scani_hr=SecondOfDay/3600.
		
end subroutine read_AMSR2_Vars

!!###############################
subroutine get_AMSR2_primary_dims(filein,grp,nscan,npixel,nchannel,status)
    USE HDF5 
    implicit none 
    character*(*), intent(in) :: filein
    character*2, intent(in) :: grp
    integer,intent(out) :: nscan,npixel,nchannel
	
    integer:: status
    integer(HID_T) :: file_id,dset_id,dspace_id
    integer(HSIZE_T) :: data_dims3(3),maxdims3(3)
    
 
    call h5open_f(status)
    call h5fopen_f(trim(filein), H5F_ACC_RDONLY_F, file_id, status)
    if (status.ne.0) then
		WRITE(*,*) 'Error: Open hdf5 failed _sub_get_AMSR2_primary_dims'
		return
    end if
    
    CALL h5dopen_f(file_id, grp//"/Tc", dset_id, status)
    if (status.ne.0) then
		WRITE(*,*) 'Error: Extract HDF Var failed _sub_get_AMSR2_primary_dims'
		return
    end if
	
    CALL h5dget_space_f(dset_id,dspace_id,status)
    CALL h5sget_simple_extent_dims_f(dspace_id, data_dims3, maxdims3, status)
    nscan = data_dims3(3) 
    npixel = data_dims3(2) 
    nchannel = data_dims3(1) 
    CALL h5sclose_f(dspace_id,status)
    CALL h5dclose_f(dset_id,status)
	
    CALL h5fclose_f(file_id,status)
    CALL h5close_f(status)
end subroutine get_AMSR2_primary_dims


