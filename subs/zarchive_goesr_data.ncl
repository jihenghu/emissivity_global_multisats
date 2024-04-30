begin

dirs=systemfunc("ls -D -I '^m*-*'")
 print(dirs)

pwd=systemfunc("pwd")
do id=0,dimsizes(dirs)-1

    files:=systemfunc("cd "+pwd+"/"+ dirs(id)+"/001/ ; ls *.nc")    
    fulln:=pwd+"/"+ dirs(id)+"/001/"+files
    
    do iF=0,dimsizes(files)-1
        eyyyydoy=systemfunc("echo "+files(iF)+" | cut -d '_' -f4 ")
        prefix=systemfunc("echo "+files(iF)+" | cut -d '_' -f1-3 ")
        yyyydoy=str_get_cols(eyyyydoy,1,7)
        HH=str_get_cols(eyyyydoy,8,9)
        MM=str_get_cols(eyyyydoy,10,11)
        ; print(HH+""+MM)
        ; print(prefix+"")
        yyyymmdd=yyyyddd_to_yyyymmdd(toint(yyyydoy))
        ; print("move "+ fulln(iF)+" to DIR:"+ yyyymmdd)
        system("mkdir "+yyyymmdd)
        system("cp -v "+ fulln(iF)+" "+ pwd+"/"+yyyymmdd+"/"+prefix+"_"+yyyymmdd+"_"+HH+MM+".NC")
    end do
end do 
end
