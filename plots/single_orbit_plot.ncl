begin
; file_name="Emissivity.Clear.GMI.Himawari8.GOESR.MSGIO.20150802-S000207-E013440.txt"

fillist=systemfunc("ls Emissivity.Clear.GMI.HIMA.GOESR.MSG.20150802-S061224-E074457.txt")

do iF=0,dimsizes(fillist)-1
    file_name=fillist(iF)
    nrow=numAsciiRow(file_name)
    ncol=numAsciiCol(file_name)
    io:=asciiread(file_name, (/nrow,ncol/), "float")

    flon:=io(:,4)
    flat:=io(:,5)
    lza:=io(:,ncol-12)
    mlse:=io(:,7:7+8)
    
    edgx:=new((/nrow,3/),"float")
    edgy:=new((/nrow,3/),"float")
    dx=0.075
    edgx(:,0)=flon
    edgx(:,1)=flon-dx*1.732
    edgx(:,2)=flon+dx*1.732
    edgy(:,0)=flat+dx*2
    edgy(:,1)=flat-dx
    edgy(:,2)=flat+dx

    wks:=gsn_open_wks("png", file_name+"_V")
    wks1:=gsn_open_wks("png", file_name+"_H")
    wks2:=gsn_open_wks("png", file_name+"LZA")
    gsn_define_colormap(wks, "WhiteBlueGreenYellowRed")
    gsn_define_colormap(wks1, "WhiteBlueGreenYellowRed")
    gsn_define_colormap(wks2, "WhiteBlueGreenYellowRed")
    res=True
    res@gsnFrame=False
    res@gsnDraw=False

    res@cnFillOn=True
    res@cnLinesOn = False
    res@cnLineLabelsOn = False
    res@cnInfoLabelOn = False
    res@cnFillMode = "CellFill"
    res@lbLabelBarOn   =False

    res@cnLevelSelectionMode  = "ExplicitLevels"  
    ; res@cnLevels=(/0.7,0.74,0.78,0.82,0.86,0.84,0.88,0.9,0.91,0.92,0.93,0.94,0.95,0.96,0.97,0.98,0.99/)
    res@cnLevels:=(/0.9,0.91,0.92,0.93,0.94,0.95,0.96,0.97,0.98,0.99/)

    res@mpFillOn=True
    res@gsnAddCyclic = False 
    res@mpMinLatF = -90.
    res@mpMaxLatF = 90.
    res@mpMinLonF = -180.
    res@mpMaxLonF = 180.

    res@mpMinLatF = min(flat)
    res@mpMaxLatF = max(flat)
    res@mpMinLonF = min(flon)-20
    res@mpMaxLonF = max(flon)+20

    ; res@mpFillOn = True
    res@cnFillDrawOrder      = "PostDraw"
    res@mpLandFillColor      = "lightgrey"
    res@mpOceanFillColor     = "white"  

    res@sfXArray  := flon
    res@sfYArray := flat
    res@sfXCellBounds :=edgx
    res@sfYCellBounds :=edgy

    channels=(/"10V","10H","18V","18H","23V","36V","36H","89V","89H"/)
    ps:=(/0,2,4,5,7/) ; Vpol
    ; ps=(/1,3,6,8/) ; Hpol

    plot:=new(dimsizes(ps),graphic)

    do ip=0,dimsizes(ps)-1

    res@gsnLeftString="GMI Emissivity, "+channels(ps(ip))
    plot(ip)=gsn_csm_contour_map(wks, mlse(:,ps(ip)), res)

    end do

    resP=True
    resP@gsnPanelLabelBar=True
    ; resP@
    gsn_panel(wks, plot, (/3,2/), resP)


    res@cnLevels:=(/0.7,0.74,0.78,0.82,0.86,0.84,0.88,0.9,0.91,0.92,0.93,0.94,0.95,0.96,0.97,0.98,0.99/)
    ps:=(/1,3,6,8/) ; Hpol
    plot:=new(dimsizes(ps),graphic)

    do ip=0,dimsizes(ps)-1

    res@gsnLeftString="GMI Emissivity, "+channels(ps(ip))
    plot(ip)=gsn_csm_contour_map(wks1, mlse(:,ps(ip)), res)

    end do

    resP=True
    resP@gsnPanelLabelBar=True
    ; resP@
    gsn_panel(wks1, plot, (/2,2/), resP)

    lza=where(lza.lt.0, lza@_FillValue,lza)
    lza=where(lza.gt.88, lza@_FillValue,lza)
    res@lbLabelBarOn=True
    res@gsnDraw=True
    res@gsnFrame=True
    res@cnLevels:=(/5,10,20,30,40,50,60,70,80,85/)
    res@gsnLeftString="LZA"
    plot1=gsn_csm_contour_map(wks2, lza, res)




end do
end