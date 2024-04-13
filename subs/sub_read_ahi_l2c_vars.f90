
subroutine read_AHI_L2_cloud(filename,latitude,longitude,ctype,status)

	include 'netcdf.inc'
	
	CHARACTER*(*) ,Intent(IN):: filename
	REAL*4, DIMENSION(2401),Intent(OUT) :: latitude,longitude

	Integer, DIMENSION(2401,2401),Intent(OUT) ::ctype

	Integer status,ncid,varid

    status=nf_open(trim(adjustl(filename)),nf_nowrite,ncid)

	status=nf_inq_varid (ncid, 'latitude', varid) 
	status=nf_get_var_real(ncid,varid,latitude)	
	
	status=nf_inq_varid (ncid, 'longitude', varid) 
	status=nf_get_var_real(ncid,varid,longitude)
	where(longitude>180) longitude= longitude-360 
		 
	! status=nf_inq_varid (ncid, 'CLER_23', varid)
	! status=nf_get_var_int(ncid,varid,icer)    
	! cer=icer*0.01 ! um
	! where(cer<0) cer=-999

	! status=nf_inq_varid (ncid, 'CLOT', varid) 
	! status=nf_get_var_int(ncid,varid,icot)
	! cot=icot*0.01 ! m-1
	! where(cot<0) cot=-999
	
	! status=nf_inq_varid (ncid, 'CLTH', varid) 
	! status=nf_get_var_int(ncid,varid,icth)
	! cth=icth*0.001 ! km
	! where(cth<0) cth=-999
	
	status=nf_inq_varid (ncid, 'CLTYPE', varid) 
	status=nf_get_var_int(ncid,varid,ctype)
	
	! status=nf_inq_varid (ncid, 'CLTT', varid)
	! status=nf_get_var_int(ncid,varid,ctt)
	! ctt=ctt*0.01+150 !K

	! status=nf_inq_varid (ncid, 'QA', varid)
	! status=nf_get_var_int(ncid,varid,qa)

	status=nf_close(ncid)
	
END SUBROUTINE read_AHI_L2_cloud

