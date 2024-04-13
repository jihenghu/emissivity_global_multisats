

!! ===========================================================
!! sub_write_emiss_hdf5.f90
!! subroutine to ouput emissivity retrievals in HDF5 format
!! 
!! By, J. Hu, 2023/12/23.
!! ===========================================================

subroutine write_profile_hdf5(filename,nchannel,npixel,nscan, & ! filename, geometries
							longitude,latitude, 			& ! Geolocations
							level1,level2,level3,	&					
							CFR1,LWC1,IWC1,	&					
							CFR2,LWC2,IWC2,	&					
							CFR3,LWC3,IWC3,	&					
							status)

      USE HDF5   
      IMPLICIT NONE
      CHARACTER(*), INTENT(IN) :: filename	  
	  Integer, INTENT(IN) :: nscan,npixel,nchannel
      REAL, DIMENSION(npixel,nscan), INTENT(IN) :: longitude, latitude  


	  ! REAL, DIMENSION(npixel,nscan) :: SKT,T2m,SnowC,SMC
	  
	  Integer, DIMENSION(npixel,nscan), INTENT(IN) :: level1,level2,level3
	  REAL, DIMENSION(npixel,nscan), INTENT(IN) :: CFR1,IWC1,LWC1
	  REAL, DIMENSION(npixel,nscan), INTENT(IN) :: CFR2,IWC2,LWC2
	  REAL, DIMENSION(npixel,nscan), INTENT(IN) :: CFR3,IWC3,LWC3

      INTEGER, INTENT(OUT)  :: status 	  
	  
	  
      CHARACTER*50:: dsetname != "TB_89VH_GRID25" ! Dataset name   
      ! Identifiers
      INTEGER(HID_T) :: file_id       ! File identifier
      ! INTEGER(HID_T) :: group_id      ! Group identifier
      INTEGER(HID_T) :: dset_id      ! Dataset 1 identifier
      INTEGER(HID_T) :: dspace_id    ! Dataspace 1 identifier
      INTEGER :: rank                ! Dataset rank
      INTEGER(HSIZE_T), DIMENSION(3)  :: data_dims4
      INTEGER(HSIZE_T), DIMENSION(3)  :: data_dims3
      INTEGER(HSIZE_T), DIMENSION(2)  :: data_dims2
      INTEGER(HSIZE_T), DIMENSION(1)  :: data_dims1
	
	  !! attribute
	  INTEGER(HID_T) :: attr_id     ! Attribute identifier
	  INTEGER(HID_T) :: aspace_id     ! Attribute Dataspace identifier
	  INTEGER(HID_T) :: atype_id      ! Attribute Dataspace identifier
	  INTEGER(HSIZE_T), DIMENSION(1) :: adims = (/1/) ! Attribute dimension
	  INTEGER     ::   arank = 1                      ! Attribure rank
	  INTEGER(SIZE_T) :: attrlen=80    ! Length of the attribute string
	
  ! =====================================================================
 
      ! Initialize the dset_data array 

      data_dims2 = (/npixel,nscan/)

	  
      ! Initialize Fortran interface
      CALL h5open_f(status) 
	  
	  	IF(status/=0) then
		PRINT*,"Write out HDF error : ",trim(adjustL(filename))
		Stop
		EndIF	
	  
      CALL h5fcreate_f(trim(adjustl(filename)), H5F_ACC_TRUNC_F, file_id, status)
	  
	  rank = 2 
      CALL h5screate_simple_f(rank, data_dims2, dspace_id, status)  	  

      dsetname="Cloud_Fraction1"		  
	  CALL h5dcreate_f(file_id, trim(adjustl(dsetname)), H5T_NATIVE_REAL, dspace_id,dset_id, status)
      CALL h5dwrite_f(dset_id, H5T_NATIVE_REAL, CFR1, data_dims2, status)
      CALL h5dclose_f(dset_id, status)  
	  
      dsetname="Cloud_Fraction2"		  
	  CALL h5dcreate_f(file_id, trim(adjustl(dsetname)), H5T_NATIVE_REAL, dspace_id,dset_id, status)
      CALL h5dwrite_f(dset_id, H5T_NATIVE_REAL, CFR2, data_dims2, status)
      CALL h5dclose_f(dset_id, status) 

	  dsetname="Cloud_Fraction3"		  
	  CALL h5dcreate_f(file_id, trim(adjustl(dsetname)), H5T_NATIVE_REAL, dspace_id,dset_id, status)
      CALL h5dwrite_f(dset_id, H5T_NATIVE_REAL, CFR3, data_dims2, status)
      CALL h5dclose_f(dset_id, status)  

  
      dsetname="LWC1"		  
	  CALL h5dcreate_f(file_id, trim(adjustl(dsetname)), H5T_NATIVE_REAL, dspace_id,dset_id, status)
      CALL h5dwrite_f(dset_id, H5T_NATIVE_REAL, LWC1, data_dims2, status)
      CALL h5dclose_f(dset_id, status)    
      dsetname="LWC2"		  
	  CALL h5dcreate_f(file_id, trim(adjustl(dsetname)), H5T_NATIVE_REAL, dspace_id,dset_id, status)
      CALL h5dwrite_f(dset_id, H5T_NATIVE_REAL, LWC2, data_dims2, status)
      CALL h5dclose_f(dset_id, status)    
      dsetname="LWC3"		  
	  CALL h5dcreate_f(file_id, trim(adjustl(dsetname)), H5T_NATIVE_REAL, dspace_id,dset_id, status)
      CALL h5dwrite_f(dset_id, H5T_NATIVE_REAL, LWC3, data_dims2, status)
      CALL h5dclose_f(dset_id, status)  

  
      dsetname="IWC1"		  
	  CALL h5dcreate_f(file_id, trim(adjustl(dsetname)), H5T_NATIVE_REAL, dspace_id,dset_id, status)
      CALL h5dwrite_f(dset_id, H5T_NATIVE_REAL, IWC1, data_dims2, status)
      CALL h5dclose_f(dset_id, status)    
      dsetname="IWC2"		  
	  CALL h5dcreate_f(file_id, trim(adjustl(dsetname)), H5T_NATIVE_REAL, dspace_id,dset_id, status)
      CALL h5dwrite_f(dset_id, H5T_NATIVE_REAL, IWC2, data_dims2, status)
      CALL h5dclose_f(dset_id, status)    
      dsetname="IWC3"		  
	  CALL h5dcreate_f(file_id, trim(adjustl(dsetname)), H5T_NATIVE_REAL, dspace_id,dset_id, status)
      CALL h5dwrite_f(dset_id, H5T_NATIVE_REAL, IWC3, data_dims2, status)
      CALL h5dclose_f(dset_id, status)  
	  
      dsetname="level1"		  
	  CALL h5dcreate_f(file_id, trim(adjustl(dsetname)), H5T_NATIVE_REAL, dspace_id,dset_id, status)
      CALL h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, level1, data_dims2, status)
      CALL h5dclose_f(dset_id, status)    
      dsetname="level2"		  
	  CALL h5dcreate_f(file_id, trim(adjustl(dsetname)), H5T_NATIVE_REAL, dspace_id,dset_id, status)
      CALL h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, level2, data_dims2, status)
      CALL h5dclose_f(dset_id, status)    
      dsetname="level3"		  
	  CALL h5dcreate_f(file_id, trim(adjustl(dsetname)), H5T_NATIVE_REAL, dspace_id,dset_id, status)
      CALL h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, level3, data_dims2, status)
      CALL h5dclose_f(dset_id, status)  
  
      CALL h5sclose_f(dspace_id, status)
 
 
	  rank = 2	  	  
      CALL h5screate_simple_f(rank, data_dims2, dspace_id, status)
      dsetname="Latitude"	  	  
      CALL h5dcreate_f(file_id, trim(adjustl(dsetname)), H5T_NATIVE_REAL, dspace_id,dset_id, status)
      CALL h5dwrite_f(dset_id, H5T_NATIVE_REAL, latitude, data_dims2, status)
	  	  CALL h5screate_simple_f(arank, adims, aspace_id, status)
		  CALL h5tcopy_f(H5T_NATIVE_CHARACTER, atype_id, status)
		  CALL h5tset_size_f(atype_id, attrlen, status)
		  attrlen=12
		  CALL h5tset_size_f(atype_id, attrlen, status)
		  CALL h5acreate_f(dset_id, "units", atype_id, aspace_id, attr_id, status)
		  CALL h5awrite_f(attr_id, atype_id, "degree_north", adims, status)
		  CALL h5aclose_f(attr_id, status)
		  
		  attrlen=31
		  CALL h5tset_size_f(atype_id, attrlen, status)		  
		  CALL h5acreate_f(dset_id, "description", atype_id, aspace_id, attr_id, status)
		  CALL h5awrite_f(attr_id, atype_id, "Latitude of each pixel in WGS84", adims, status)
		  CALL h5aclose_f(attr_id, status)		  
		  
		  CALL h5tclose_f(atype_id, status)
		  CALL h5sclose_f(aspace_id, status)
      CALL h5dclose_f(dset_id, status)

      dsetname="Longitude"	  	  
      CALL h5dcreate_f(file_id, trim(adjustl(dsetname)), H5T_NATIVE_REAL, dspace_id,dset_id, status)
      CALL h5dwrite_f(dset_id, H5T_NATIVE_REAL, longitude, data_dims2, status)
	  	  CALL h5screate_simple_f(arank, adims, aspace_id, status)
		  CALL h5tcopy_f(H5T_NATIVE_CHARACTER, atype_id, status)

		  attrlen=11
		  CALL h5tset_size_f(atype_id, attrlen, status)
		  CALL h5acreate_f(dset_id, "units", atype_id, aspace_id, attr_id, status)
		  CALL h5awrite_f(attr_id, atype_id, "degree_east", adims, status)
		  CALL h5aclose_f(attr_id, status)
		  
		  attrlen=32
		  CALL h5tset_size_f(atype_id, attrlen, status)		  
		  CALL h5acreate_f(dset_id, "description", atype_id, aspace_id, attr_id, status)
		  CALL h5awrite_f(attr_id, atype_id, "Longitude of each pixel in WGS84", adims, status)
		  CALL h5aclose_f(attr_id, status)		  
		  
		  CALL h5tclose_f(atype_id, status)
		  CALL h5sclose_f(aspace_id, status)
      CALL h5dclose_f(dset_id, status)	  
	   	  
      CALL h5sclose_f(dspace_id, status)
	  	  	 
     ! Close the file
     CALL h5fclose_f(file_id, status)
     ! Close FORTRAN interface
     CALL h5close_f(status)
	 


end subroutine write_profile_hdf5