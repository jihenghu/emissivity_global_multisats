

!! subroutine sub_read_tb_cloud_prams.f90
!! to read in the hdf varibles from collocation procedure.

! subroutine sub_read_tb_cloud_prams()
	! USE HDF5
	

! end subroutine sub_read_tb_cloud_prams


subroutine get_tb_cloud_dims(filename,nscan,npixel,nchannel)
	USE HDF5 
	IMPLICIT NONE
	CHARACTER*(*),INTENT(IN) :: filename
	INTEGER,INTENT(OUT) :: nscan,npixel,nchannel

    INTEGER:: status
	INTEGER(HID_T) :: file_id,dset_id,dspace_id
	INTEGER(HSIZE_T) :: data_dims3(3),maxdims3(3)

	call h5open_f(status)
	call h5fopen_f(trim(filename), H5F_ACC_RDONLY_F, file_id, status)
	if (status.ne.0) then
		WRITE(*,*) 'Error: Open hdf5 failed _sub_get_GMI_primary_dims'
		STOP
	end if

	CALL h5dopen_f(file_id, "TB_GMI_L1C", dset_id, status)
	if (status.ne.0) then
		WRITE(*,*) 'Error: Extract HDF Var failed _sub_get_GMI_primary_dims'
		STOP
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



end subroutine get_tb_cloud_dims