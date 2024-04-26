begin

path="/home/jihenghu/data00/data_em/GMI_EMISSIVITY/"

dates=systemfunc("ls "+path)
numfiles=new(dimsizes(dates), integer)

do id=0, dimsizes(dates)-1
    filep=path+dates(id)+"/*.txt"
    files:=systemfunc("ls "+filep)
    nf = dimsizes(files)
    if (ismissing(files(0))) then nf=0 end if
    numfiles(id)=nf
end do

wks=gsn_open_wks("png","orbits")
res=True
res@gsnXYBarChartBarWidth=0.2
res@gsnXYBarChart=True

print(dates)
x=ispan(0,dimsizes(dates)-1,1)
res@tmXBLabels=dates
res@tmXBValues=x
res@tmXBMode="Explicit"
res@tmXBLabelStride=2
res@tmXBLabelFontHeightF=0.01
res@tmXBLabelAngleF=90
res@tmXBLabelJust="BottomRight"

plot=gsn_csm_xy(wks,ispan(0,dimsizes(dates)-1,1), numfiles,res)

end
