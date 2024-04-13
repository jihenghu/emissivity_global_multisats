

!! ===========================================================
!! sub_write_emiss_hdf5.f90
!! subroutine to ouput emissivity retrievals in HDF5 format
!! 
!! By, J. Hu, 2023/12/23.
!! ===========================================================

subroutine write_emiss_hdf5(filename,nchannel,npixel,nscan, & ! filename, geometries
							longitude,latitude, 			& ! Geolocations
							stime, 							& ! scan time in UTC hr
							TB,								& ! TOA brightness temperatures
							emissivity,						& ! surface emissivity 
							SKT,							& ! skin emissivity , ERA5-land
							T2m,							& ! 2m air temperature, ERA5-land
							SnowC,							& ! sonow cover, ERA5-land 
							SMC,							& ! top layer soil moisture, ERA5-land
							LZA, 							& ! Local view zenith angle of geosatellite
							cfr, 							& ! cloud fraction in percentage
							clear, 							& ! clear sky percentage
							lsm,							& ! land sea mask: 0: land, 100: openwater
							flag, 							& ! cloud source flags
							missrate,						& ! missing percentage of cloud properties
							status)

      USE HDF5   
      IMPLICIT NONE
      CHARACTER(*), INTENT(IN) :: filename	  
	  Integer, INTENT(IN) :: nscan,npixel,nchannel
      REAL, DIMENSION(npixel,nscan), INTENT(IN) :: longitude, latitude  
      ! REAL, DIMENSION(nscan), INTENT(IN) 	  :: stime	  
      REAL, DIMENSION(npixel,nscan), INTENT(IN) 	  :: stime	  
	  REAL, DIMENSION(nchannel,npixel,nscan), INTENT(IN) :: TB,emissivity

	  REAL, DIMENSION(npixel,nscan) :: SKT,T2m,SnowC,SMC
	  

	  REAL, DIMENSION(3,npixel,nscan), INTENT(IN) :: cfr,clear
	  REAL, DIMENSION(npixel,nscan), INTENT(IN) :: lsm,LZA
	  REAL, DIMENSION(3,npixel,nscan), INTENT(IN) :: missrate
	  INTEGER, DIMENSION(npixel,nscan), INTENT(IN) :: flag
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
      data_dims4 = (/3,npixel,nscan/)
      data_dims3 = (/nchannel,npixel,nscan/)
      data_dims2 = (/npixel,nscan/)
      data_dims1 = (/nscan/)
	  
      ! Initialize Fortran interface
      CALL h5open_f(status) 
	  
	  	IF(status/=0) then
		PRINT*,"Write out HDF error : ",trim(adjustL(filename))
		Stop
		EndIF	
	  
      CALL h5fcreate_f(trim(adjustl(filename)), H5F_ACC_TRUNC_F, file_id, status)
	  
	  rank = 3	  
      CALL h5screate_simple_f(rank, data_dims3, dspace_id, status)  	  

	  dsetname="TB_GMI_L1C"		  
      CALL h5dcreate_f(file_id, trim(adjustl(dsetname)), H5T_NATIVE_REAL, dspace_id,dset_id, status)
      CALL h5dwrite_f(dset_id, H5T_NATIVE_REAL, TB, data_dims3, status)
	  	  CALL h5screate_simple_f(arank, adims, aspace_id, status)
		  CALL h5tcopy_f(H5T_NATIVE_CHARACTER, atype_id, status)

		  attrlen=1
		  CALL h5tset_size_f(atype_id, attrlen, status)
		  CALL h5acreate_f(dset_id, "units", atype_id, aspace_id, attr_id, status)
		  CALL h5awrite_f(attr_id, atype_id, "K", adims, status)
		  CALL h5aclose_f(attr_id, status)

		  attrlen=39
		  CALL h5tset_size_f(atype_id, attrlen, status)		  
		  CALL h5acreate_f(dset_id, "long_name", atype_id, aspace_id, attr_id, status)
		  CALL h5awrite_f(attr_id, atype_id, "GPM_L1C GMI TOA Brightness Temperatures", adims, status)
		  CALL h5aclose_f(attr_id, status)
		  attrlen=44
		  CALL h5tset_size_f(atype_id, attrlen, status)
		  CALL h5acreate_f(dset_id, "description", atype_id, aspace_id, attr_id, status)
		  CALL h5awrite_f(attr_id, atype_id, "channel: 10V,10H,18V,18H,23V,36V,36H,89V,89H", adims, status)
		  CALL h5aclose_f(attr_id, status)		  
		  
		  CALL h5tclose_f(atype_id, status)
		  CALL h5sclose_f(aspace_id, status)	 	  
      CALL h5dclose_f(dset_id, status)
	  !!    emissivity
	  dsetname="Emissivity"		  
      CALL h5dcreate_f(file_id, trim(adjustl(dsetname)), H5T_NATIVE_REAL, dspace_id,dset_id, status)
      CALL h5dwrite_f(dset_id, H5T_NATIVE_REAL, emissivity, data_dims3, status)
	  	  CALL h5screate_simple_f(arank, adims, aspace_id, status)
		  CALL h5tcopy_f(H5T_NATIVE_CHARACTER, atype_id, status)

		  ! attrlen=1
		  ! CALL h5tset_size_f(atype_id, attrlen, status)
		  ! CALL h5acreate_f(dset_id, "units", atype_id, aspace_id, attr_id, status)
		  ! CALL h5awrite_f(attr_id, atype_id, "K", adims, status)
		  ! CALL h5aclose_f(attr_id, status)

		  attrlen=22
		  CALL h5tset_size_f(atype_id, attrlen, status)		  
		  CALL h5acreate_f(dset_id, "long_name", atype_id, aspace_id, attr_id, status)
		  CALL h5awrite_f(attr_id, atype_id, "GMI surface emissivity", adims, status)
		  CALL h5aclose_f(attr_id, status)
		  attrlen=44
		  CALL h5tset_size_f(atype_id, attrlen, status)
		  CALL h5acreate_f(dset_id, "description", atype_id, aspace_id, attr_id, status)
		  CALL h5awrite_f(attr_id, atype_id, "channel: 10V,10H,18V,18H,23V,36V,36H,89V,89H", adims, status)
		  CALL h5aclose_f(attr_id, status)		  
		  
		  CALL h5tclose_f(atype_id, status)
		  CALL h5sclose_f(aspace_id, status)	 	  
      CALL h5dclose_f(dset_id, status)

      CALL h5sclose_f(dspace_id, status)

	  rank = 3	  
      CALL h5screate_simple_f(rank, data_dims4, dspace_id, status)  	  

      dsetname="Cloud_Fraction"		  
	  CALL h5dcreate_f(file_id, trim(adjustl(dsetname)), H5T_NATIVE_REAL, dspace_id,dset_id, status)
      CALL h5dwrite_f(dset_id, H5T_NATIVE_REAL, cfr, data_dims4, status)
	  	  CALL h5screate_simple_f(arank, adims, aspace_id, status)
		  CALL h5tcopy_f(H5T_NATIVE_CHARACTER, atype_id, status)

		  attrlen=1
		  CALL h5tset_size_f(atype_id, attrlen, status)
		  CALL h5acreate_f(dset_id, "units", atype_id, aspace_id, attr_id, status)
		  CALL h5awrite_f(attr_id, atype_id, "%", adims, status)
		  CALL h5aclose_f(attr_id, status)
		  attrlen=7
		  CALL h5tset_size_f(atype_id, attrlen, status)		  
		  CALL h5acreate_f(dset_id, "valid_range", atype_id, aspace_id, attr_id, status)
		  CALL h5awrite_f(attr_id, atype_id, "[0:100]", adims, status)
		  CALL h5aclose_f(attr_id, status)
		  attrlen=26
		  CALL h5tset_size_f(atype_id, attrlen, status)		  
		  CALL h5acreate_f(dset_id, "long_name", atype_id, aspace_id, attr_id, status)
		  CALL h5awrite_f(attr_id, atype_id, "Cloud Fraction in MW Beams", adims, status)
		  CALL h5aclose_f(attr_id, status)
		  
		  attrlen=50
		  CALL h5tset_size_f(atype_id, attrlen, status)
		  CALL h5acreate_f(dset_id, "description", atype_id, aspace_id, attr_id, status)
		  CALL h5awrite_f(attr_id, atype_id, "[1,:,:] 10.62 VH; [2,:,:] 18/23/36VH; [3,:,:] 89VH", adims, status)
		  CALL h5aclose_f(attr_id, status)				  
		  
		  CALL h5tclose_f(atype_id, status)
		  CALL h5sclose_f(aspace_id, status)	
      CALL h5dclose_f(dset_id, status)  

      dsetname="ClearSky_Fraction"		  
	  CALL h5dcreate_f(file_id, trim(adjustl(dsetname)), H5T_NATIVE_REAL, dspace_id,dset_id, status)
      CALL h5dwrite_f(dset_id, H5T_NATIVE_REAL, clear, data_dims4, status)
	  	  CALL h5screate_simple_f(arank, adims, aspace_id, status)
		  CALL h5tcopy_f(H5T_NATIVE_CHARACTER, atype_id, status)

		  attrlen=1
		  CALL h5tset_size_f(atype_id, attrlen, status)
		  CALL h5acreate_f(dset_id, "units", atype_id, aspace_id, attr_id, status)
		  CALL h5awrite_f(attr_id, atype_id, "%", adims, status)
		  CALL h5aclose_f(attr_id, status)
		  attrlen=7
		  CALL h5tset_size_f(atype_id, attrlen, status)		  
		  CALL h5acreate_f(dset_id, "valid_range", atype_id, aspace_id, attr_id, status)
		  CALL h5awrite_f(attr_id, atype_id, "[0:100]", adims, status)
		  CALL h5aclose_f(attr_id, status)
		  attrlen=26
		  CALL h5tset_size_f(atype_id, attrlen, status)		  
		  CALL h5acreate_f(dset_id, "long_name", atype_id, aspace_id, attr_id, status)
		  CALL h5awrite_f(attr_id, atype_id, "Clear Fraction in MW Beams", adims, status)
		  CALL h5aclose_f(attr_id, status)
		  
		  attrlen=50
		  CALL h5tset_size_f(atype_id, attrlen, status)
		  CALL h5acreate_f(dset_id, "description", atype_id, aspace_id, attr_id, status)
		  CALL h5awrite_f(attr_id, atype_id, "[1,:,:] 10.62 VH; [2,:,:] 18/23/36VH; [3,:,:] 89VH", adims, status)
		  CALL h5aclose_f(attr_id, status)		  
		  
		  CALL h5tclose_f(atype_id, status)
		  CALL h5sclose_f(aspace_id, status)	
      CALL h5dclose_f(dset_id, status)  


  

      dsetname="Cloud_Missing_Fraction"		  
	  CALL h5dcreate_f(file_id, trim(adjustl(dsetname)), H5T_NATIVE_REAL, dspace_id,dset_id, status)
      CALL h5dwrite_f(dset_id, H5T_NATIVE_REAL, missrate, data_dims4, status)
	  	  CALL h5screate_simple_f(arank, adims, aspace_id, status)
		  CALL h5tcopy_f(H5T_NATIVE_CHARACTER, atype_id, status)

		  attrlen=1
		  CALL h5tset_size_f(atype_id, attrlen, status)
		  CALL h5acreate_f(dset_id, "units", atype_id, aspace_id, attr_id, status)
		  CALL h5awrite_f(attr_id, atype_id, "%", adims, status)
		  CALL h5aclose_f(attr_id, status)
		  attrlen=7
		  CALL h5tset_size_f(atype_id, attrlen, status)		  
		  CALL h5acreate_f(dset_id, "valid_range", atype_id, aspace_id, attr_id, status)
		  CALL h5awrite_f(attr_id, atype_id, "[0:100]", adims, status)
		  CALL h5aclose_f(attr_id, status)
		  attrlen=34
		  CALL h5tset_size_f(atype_id, attrlen, status)		  
		  CALL h5acreate_f(dset_id, "long_name", atype_id, aspace_id, attr_id, status)
		  CALL h5awrite_f(attr_id, atype_id, "Cloud Missing Fraction in MW Beams", adims, status)
		  CALL h5aclose_f(attr_id, status)

		  attrlen=50
		  CALL h5tset_size_f(atype_id, attrlen, status)
		  CALL h5acreate_f(dset_id, "description", atype_id, aspace_id, attr_id, status)
		  CALL h5awrite_f(attr_id, atype_id, "[1,:,:] 10.62 VH; [2,:,:] 18/23/36VH; [3,:,:] 89VH", adims, status)
		  CALL h5aclose_f(attr_id, status)		
		  	  
		  CALL h5tclose_f(atype_id, status)
		  CALL h5sclose_f(aspace_id, status)	
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
	   	  
      dsetname="Land_Sea_Fraction"	  	  
      CALL h5dcreate_f(file_id, trim(adjustl(dsetname)), H5T_NATIVE_REAL, dspace_id,dset_id, status)
      CALL h5dwrite_f(dset_id, H5T_NATIVE_REAL, lsm, data_dims2, status)
	  	  CALL h5screate_simple_f(arank, adims, aspace_id, status)
		  CALL h5tcopy_f(H5T_NATIVE_CHARACTER, atype_id, status)
		  
		  attrlen=1
		  CALL h5tset_size_f(atype_id, attrlen, status)
		  CALL h5acreate_f(dset_id, "units", atype_id, aspace_id, attr_id, status)
		  CALL h5awrite_f(attr_id, atype_id, "%", adims, status)
		  CALL h5aclose_f(attr_id, status)		  

		  attrlen=17
		  CALL h5tset_size_f(atype_id, attrlen, status)
		  CALL h5acreate_f(dset_id, "long_name", atype_id, aspace_id, attr_id, status)
		  CALL h5awrite_f(attr_id, atype_id, "Land Sea Fraction", adims, status)
		  CALL h5aclose_f(attr_id, status)
		  
		  attrlen=77
		  CALL h5tset_size_f(atype_id, attrlen, status)
		  CALL h5acreate_f(dset_id, "description", atype_id, aspace_id, attr_id, status)
		  CALL h5awrite_f(attr_id, atype_id, &
		  "0%(strict land) - 25% (typical land) - 75% (coastal) - 100%(strict openwater)", adims, status)
		  CALL h5aclose_f(attr_id, status)		  
		  
		  CALL h5tclose_f(atype_id, status)
		  CALL h5sclose_f(aspace_id, status)
      CALL h5dclose_f(dset_id, status)	 	 
  	  
      dsetname="Local_view_Zenith_Angle"	  	  
      CALL h5dcreate_f(file_id, trim(adjustl(dsetname)), H5T_NATIVE_REAL, dspace_id,dset_id, status)
      CALL h5dwrite_f(dset_id, H5T_NATIVE_REAL, LZA, data_dims2, status)
	  	  CALL h5screate_simple_f(arank, adims, aspace_id, status)
		  CALL h5tcopy_f(H5T_NATIVE_CHARACTER, atype_id, status)
		  
		  attrlen=7
		  CALL h5tset_size_f(atype_id, attrlen, status)
		  CALL h5acreate_f(dset_id, "units", atype_id, aspace_id, attr_id, status)
		  CALL h5awrite_f(attr_id, atype_id, "degrees", adims, status)
		  CALL h5aclose_f(attr_id, status)		  

		  attrlen=50
		  CALL h5tset_size_f(atype_id, attrlen, status)
		  CALL h5acreate_f(dset_id, "long_name", atype_id, aspace_id, attr_id, status)
		  CALL h5awrite_f(attr_id, atype_id, "Local View Zenith Angle of Geostationary Satellite", adims, status)
		  CALL h5aclose_f(attr_id, status)
		  
		  attrlen=6
		  CALL h5tset_size_f(atype_id, attrlen, status)
		  CALL h5acreate_f(dset_id, "valide_range", atype_id, aspace_id, attr_id, status)
		  CALL h5awrite_f(attr_id, atype_id, "[0,90]", adims, status)
		  CALL h5aclose_f(attr_id, status)		  
		  
		  CALL h5tclose_f(atype_id, status)
		  CALL h5sclose_f(aspace_id, status)
      CALL h5dclose_f(dset_id, status)	  
	   
      dsetname="Cloud_Flag"	  	  
      CALL h5dcreate_f(file_id, trim(adjustl(dsetname)), H5T_NATIVE_INTEGER, dspace_id,dset_id, status)
      CALL h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, flag, data_dims2, status)
	  	  CALL h5screate_simple_f(arank, adims, aspace_id, status)
		  CALL h5tcopy_f(H5T_NATIVE_CHARACTER, atype_id, status)
		  
		  attrlen=25
		  CALL h5tset_size_f(atype_id, attrlen, status)
		  CALL h5acreate_f(dset_id, "long_name", atype_id, aspace_id, attr_id, status)
		  CALL h5awrite_f(attr_id, atype_id, "Cloud Product Information", adims, status)
		  CALL h5aclose_f(attr_id, status)
		  
		  attrlen=6
		  CALL h5tset_size_f(atype_id, attrlen, status)		  
		  CALL h5acreate_f(dset_id, "valid_range", atype_id, aspace_id, attr_id, status)
		  CALL h5awrite_f(attr_id, atype_id, "[-1:9]", adims, status)
		  CALL h5aclose_f(attr_id, status)
		  
		  attrlen=166
		  CALL h5tset_size_f(atype_id, attrlen, status)
		  CALL h5acreate_f(dset_id, "description", atype_id, aspace_id, attr_id, status)
		  CALL h5awrite_f(attr_id, atype_id, &
		  "CLOUD FLAG:: -1:sea, 0:NO-GEOSAT, 1:Himawari-8, 2:GOES-R, 3:MSG1/2, 4:H8 missing,"// &
		  " 5:GOES missing; 6:MSG1/2 missing, 7:H-8 LZA degrad; 8:GOESR_LZA>65; 9:MSG LZA degrad", adims, status)
		  CALL h5aclose_f(attr_id, status)		  
		  
		  CALL h5tclose_f(atype_id, status)
		  CALL h5sclose_f(aspace_id, status)
      CALL h5dclose_f(dset_id, status)	 

		
      dsetname="Scan_Time_UTC"	  	  
      CALL h5dcreate_f(file_id, trim(adjustl(dsetname)), H5T_NATIVE_REAL, dspace_id,dset_id, status)
      CALL h5dwrite_f(dset_id, H5T_NATIVE_REAL, stime, data_dims2, status)
		
		  CALL h5screate_simple_f(arank, adims, aspace_id, status)
		  CALL h5tcopy_f(H5T_NATIVE_CHARACTER, atype_id, status)
		  
		  attrlen=4
		  CALL h5tset_size_f(atype_id, attrlen, status)
		  CALL h5acreate_f(dset_id, "units", atype_id, aspace_id, attr_id, status)
		  CALL h5awrite_f(attr_id, atype_id, "a.m.", adims, status)
		  CALL h5aclose_f(attr_id, status)		  

		  
		  attrlen=20
		  CALL h5tset_size_f(atype_id, attrlen, status)
		  CALL h5acreate_f(dset_id, "long_name", atype_id, aspace_id, attr_id, status)
		  CALL h5awrite_f(attr_id, atype_id, "Sensor Scan Time UTC", adims, status)
		  CALL h5aclose_f(attr_id, status)
		  
		  attrlen=26
		  CALL h5tset_size_f(atype_id, attrlen, status)	  
		  CALL h5acreate_f(dset_id, "description", atype_id, aspace_id, attr_id, status)
		  CALL h5awrite_f(attr_id, atype_id, "hour= HH+MM/60+SS/3600+...", adims, status)
		  CALL h5aclose_f(attr_id, status)
		  
		  CALL h5tclose_f(atype_id, status)
		  CALL h5sclose_f(aspace_id, status)		  
      CALL h5dclose_f(dset_id, status)  

		
      dsetname="Skin_Temperature"	  	  
      CALL h5dcreate_f(file_id, trim(adjustl(dsetname)), H5T_NATIVE_REAL, dspace_id,dset_id, status)
      CALL h5dwrite_f(dset_id, H5T_NATIVE_REAL, SKT, data_dims2, status)
		
		  CALL h5screate_simple_f(arank, adims, aspace_id, status)
		  CALL h5tcopy_f(H5T_NATIVE_CHARACTER, atype_id, status)
		  
		  attrlen=1
		  CALL h5tset_size_f(atype_id, attrlen, status)
		  CALL h5acreate_f(dset_id, "units", atype_id, aspace_id, attr_id, status)
		  CALL h5awrite_f(attr_id, atype_id, "K", adims, status)
		  CALL h5aclose_f(attr_id, status)		  
		  
		  attrlen=36
		  CALL h5tset_size_f(atype_id, attrlen, status)
		  CALL h5acreate_f(dset_id, "long_name", atype_id, aspace_id, attr_id, status)
		  CALL h5awrite_f(attr_id, atype_id, "Surface skin temperature (ERA5-land)", adims, status)
		  CALL h5aclose_f(attr_id, status)
		  		  
		  CALL h5tclose_f(atype_id, status)
		  CALL h5sclose_f(aspace_id, status)		  
      CALL h5dclose_f(dset_id, status)  

		
      dsetname="2m_Temperature"	  	  
      CALL h5dcreate_f(file_id, trim(adjustl(dsetname)), H5T_NATIVE_REAL, dspace_id,dset_id, status)
      CALL h5dwrite_f(dset_id, H5T_NATIVE_REAL, T2m, data_dims2, status)
		
		  CALL h5screate_simple_f(arank, adims, aspace_id, status)
		  CALL h5tcopy_f(H5T_NATIVE_CHARACTER, atype_id, status)
		  
		  attrlen=1
		  CALL h5tset_size_f(atype_id, attrlen, status)
		  CALL h5acreate_f(dset_id, "units", atype_id, aspace_id, attr_id, status)
		  CALL h5awrite_f(attr_id, atype_id, "K", adims, status)
		  CALL h5aclose_f(attr_id, status)		  
		  
		  attrlen=30
		  CALL h5tset_size_f(atype_id, attrlen, status)
		  CALL h5acreate_f(dset_id, "long_name", atype_id, aspace_id, attr_id, status)
		  CALL h5awrite_f(attr_id, atype_id, "2m air temperature (ERA5-land)", adims, status)
		  CALL h5aclose_f(attr_id, status)
		  		  
		  CALL h5tclose_f(atype_id, status)
		  CALL h5sclose_f(aspace_id, status)		  
      CALL h5dclose_f(dset_id, status)  
	
      dsetname="Snow_Cover"	  	  
      CALL h5dcreate_f(file_id, trim(adjustl(dsetname)), H5T_NATIVE_REAL, dspace_id,dset_id, status)
      CALL h5dwrite_f(dset_id, H5T_NATIVE_REAL, SnowC, data_dims2, status)
		
		  CALL h5screate_simple_f(arank, adims, aspace_id, status)
		  CALL h5tcopy_f(H5T_NATIVE_CHARACTER, atype_id, status)
		  
		  attrlen=1
		  CALL h5tset_size_f(atype_id, attrlen, status)
		  CALL h5acreate_f(dset_id, "units", atype_id, aspace_id, attr_id, status)
		  CALL h5awrite_f(attr_id, atype_id, "K", adims, status)
		  CALL h5aclose_f(attr_id, status)		  
		  
		  attrlen=22
		  CALL h5tset_size_f(atype_id, attrlen, status)
		  CALL h5acreate_f(dset_id, "long_name", atype_id, aspace_id, attr_id, status)
		  CALL h5awrite_f(attr_id, atype_id, "Snow cover (ERA5-land)", adims, status)
		  CALL h5aclose_f(attr_id, status)
		  
		  attrlen=5
		  CALL h5tset_size_f(atype_id, attrlen, status)		  
		  CALL h5acreate_f(dset_id, "valid_range", atype_id, aspace_id, attr_id, status)
		  CALL h5awrite_f(attr_id, atype_id, "[0:1]", adims, status)
		  CALL h5aclose_f(attr_id, status)
		  
		  CALL h5tclose_f(atype_id, status)
		  CALL h5sclose_f(aspace_id, status)		  
      CALL h5dclose_f(dset_id, status)  

		
      dsetname="Soil_Moisture_Layer1"	  	  
      CALL h5dcreate_f(file_id, trim(adjustl(dsetname)), H5T_NATIVE_REAL, dspace_id,dset_id, status)
      CALL h5dwrite_f(dset_id, H5T_NATIVE_REAL, SMC, data_dims2, status)
		
		  CALL h5screate_simple_f(arank, adims, aspace_id, status)
		  CALL h5tcopy_f(H5T_NATIVE_CHARACTER, atype_id, status)
		  
		  attrlen=5
		  CALL h5tset_size_f(atype_id, attrlen, status)
		  CALL h5acreate_f(dset_id, "units", atype_id, aspace_id, attr_id, status)
		  CALL h5awrite_f(attr_id, atype_id, "m3/m3", adims, status)
		  CALL h5aclose_f(attr_id, status)		  
		  
		  attrlen=45
		  CALL h5tset_size_f(atype_id, attrlen, status)
		  CALL h5acreate_f(dset_id, "long_name", atype_id, aspace_id, attr_id, status)
		  CALL h5awrite_f(attr_id, atype_id, "Toper layer soil moisture content (ERA5-land)", adims, status)
		  CALL h5aclose_f(attr_id, status)
		  
		  attrlen=5
		  CALL h5tset_size_f(atype_id, attrlen, status)		  
		  CALL h5acreate_f(dset_id, "valid_range", atype_id, aspace_id, attr_id, status)
		  CALL h5awrite_f(attr_id, atype_id, "[0:1]", adims, status)
		  CALL h5aclose_f(attr_id, status)
		  
		  CALL h5tclose_f(atype_id, status)
		  CALL h5sclose_f(aspace_id, status)		  
      CALL h5dclose_f(dset_id, status)  
	  

	  
      CALL h5sclose_f(dspace_id, status)
	  	  	 
     ! Close the file
     CALL h5fclose_f(file_id, status)
     ! Close FORTRAN interface
     CALL h5close_f(status)
	 


end subroutine write_emiss_hdf5