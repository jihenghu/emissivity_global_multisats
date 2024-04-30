
use HDF5
character*6 yyyymm
character*255 filename
integer file_unit

integer nscan,npixel,iscan,ipixel,cloudflg
real flon,flat,scantime,sea,lza,lst,t2m,snow,smc
real emiss(9),tb(9)
real cfr(3),clear(3),miss(3)

integer,parameter ::  nlon=1441, nlat=721
Real mlse(nlon,nlat,9), longitude(nlon,nlat), latitude(nlon,nlat)
integer ngrid(nlon,nlat,9)

integer ilon, ilat, ifr ,is

character*255 hdfname
INTEGER(HID_T) :: file_id       ! File identifier
INTEGER(HID_T) :: dset_id     
INTEGER(HID_T) :: dspace_id    
INTEGER :: rank ,status               ! Dataset rank
INTEGER(HSIZE_T), DIMENSION(3)  :: data_dims3
INTEGER(HSIZE_T), DIMENSION(2)  :: data_dims2

INTEGER(HID_T) :: attr_id     ! Attribute identifier
INTEGER(HID_T) :: aspace_id     ! Attribute Dataspace identifier
INTEGER(HID_T) :: atype_id      ! Attribute Dataspace identifier
INTEGER(HSIZE_T), DIMENSION(1) :: adims = (/1/) ! Attribute dimension
INTEGER     ::   arank = 1                      ! Attribure rank
INTEGER(SIZE_T) :: attrlen=80    ! Length of the attribute string

mlse=0.0
ngrid=0
is=0

yyyymm="201807"
call system("ls /home/jihenghu/data00/data_em/GMI_EMISSIVITY/"//yyyymm//"??/*.txt > avg_month/filelist."//yyyymm//".txt")

open(21, file="avg_month/filelist."//yyyymm//".txt")

177 read(21,'(A)', end=404) filename
	print*, trim(filename)
	is=is+1
	! if(is>30) goto 404
	open(NEWUNIT=file_unit, file=trim(filename))	
176 read(file_unit,505, end=1001) nscan,npixel,iscan,ipixel,flon,flat,scantime,&
			emiss,tb,lst,t2m,snow,smc, lza,cfr,clear,miss,sea,cloudflg
	ilon=nint((flon+180.125)/0.25)+1
	ilat=nint((90.125-flat)/0.25)+1
	do ifr=1,9
		if(emiss(ifr).lt.0.or. emiss(ifr).gt.1.05 ) cycle
		mlse(ilon,ilat,ifr)=mlse(ilon,ilat,ifr)+emiss(ifr)
		ngrid(ilon,ilat,ifr)=ngrid(ilon,ilat,ifr)+1
	end do
		
	goto 176 !! another record line
505 format(4I6, 3f10.4, 18f10.4, 5f10.4, 10f10.4, I6)	
	
1001 continue	 
	goto 177 !! another file

404 continue 

	where(ngrid>0) mlse=mlse/ngrid
	
	do ilon=1,nlon
		longitude(ilon,:)=-180+(ilon-1)*0.25
	end do
	do ilat=1,nlat
		latitude(:,ilat)=90-(ilat-1)*0.25
	end do
	

	hdfname="GMI_Emissivity_Monthly_"//yyyymm//".HDF5"

	data_dims3 =(/nlon,nlat,9/)
	! Initialize Fortran interface
	CALL h5open_f(status) 
	CALL h5fcreate_f(trim(adjustl(hdfname)), H5F_ACC_TRUNC_F, file_id, status)	    
	
	CALL h5screate_simple_f(3, data_dims3, dspace_id, status)  	  	  
	CALL h5dcreate_f(file_id,"Emissivity", H5T_NATIVE_REAL, dspace_id,dset_id, status)
	CALL h5dwrite_f(dset_id, H5T_NATIVE_REAL, mlse, data_dims3, status)
	CALL h5dclose_f(dset_id, status)	
	CALL h5sclose_f(dspace_id, status)
	
	data_dims2=(/1441,721/)
	CALL h5screate_simple_f(2, data_dims2, dspace_id, status)
	
	CALL h5dcreate_f(file_id, "Latitude", H5T_NATIVE_REAL, dspace_id,dset_id, status)
	CALL h5dwrite_f(dset_id, H5T_NATIVE_REAL, latitude, data_dims2, status)
		  CALL h5screate_simple_f(arank, adims, aspace_id, status)
		  CALL h5tcopy_f(H5T_NATIVE_CHARACTER, atype_id, status)
		  CALL h5tset_size_f(atype_id, attrlen, status)
		  attrlen=12
		  CALL h5tset_size_f(atype_id, attrlen, status)
		  CALL h5acreate_f(dset_id, "units", atype_id, aspace_id, attr_id, status)
		  CALL h5awrite_f(attr_id, atype_id, "degree_north", adims, status)
		  CALL h5aclose_f(attr_id, status)
		    
		  CALL h5tclose_f(atype_id, status)
		  CALL h5sclose_f(aspace_id, status)
	CALL h5dclose_f(dset_id, status)

	CALL h5dcreate_f(file_id, "Longitude", H5T_NATIVE_REAL, dspace_id,dset_id, status)
	CALL h5dwrite_f(dset_id, H5T_NATIVE_REAL, longitude, data_dims2, status)
		  CALL h5screate_simple_f(arank, adims, aspace_id, status)
		  CALL h5tcopy_f(H5T_NATIVE_CHARACTER, atype_id, status)
		  CALL h5tset_size_f(atype_id, attrlen, status)
		  attrlen=11
		  CALL h5tset_size_f(atype_id, attrlen, status)
		  CALL h5acreate_f(dset_id, "units", atype_id, aspace_id, attr_id, status)
		  CALL h5awrite_f(attr_id, atype_id, "degree_east", adims, status)
		  CALL h5aclose_f(attr_id, status)
	CALL h5dclose_f(dset_id, status)	  	
	CALL h5sclose_f(dspace_id, status)
	
	CALL h5fclose_f(file_id, status)
	CALL h5close_f(status)

end