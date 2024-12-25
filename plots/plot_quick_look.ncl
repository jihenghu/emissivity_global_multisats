begin
; file_name="Emissivity.Clear.GMI.Himawari8.GOESR.MSGIO.20150802-S000207-E013440.txt"

; yyyymmdd=20170301
path="/home/jihenghu/data05/FY3D_EMISSIVITY/"+sprinti("%0.4i",yyyymmdd/10000)+"/"+yyyymmdd+"/"

fillist=systemfunc("cd "+path+" ; ls Emissivity.Clear.*.txt")

nrows=0
ncols=0
do iF=0,dimsizes(fillist)-1
    file_name=path+fillist(iF)
    nrow=numAsciiRow(file_name)
    nrows=nrows+nrow
end do 

fmlse=new((/nrows,10/),float,-999.9)
flon=new(nrows,float,-999.9)
flat=new(nrows,float,-999.9)
flag=new(nrows,float)

irow=0
do iF=0,dimsizes(fillist)-1
    file_name=path+fillist(iF)
    nrow=numAsciiRow(file_name)
    ncol=numAsciiCol(file_name)
    io:=asciiread(file_name, (/nrow,ncol/), "float")

    flon(irow:irow+nrow-1)=io(:,4)
    flat(irow:irow+nrow-1)=io(:,5)
    flag(irow:irow+nrow-1)=io(:,ncol-1)
    fmlse(irow:irow+nrow-1,:)=io(:,7:7+9)  

    irow=irow+nrow
end do   

fmlse=where(fmlse.lt.0,fmlse@_FillValue,fmlse)

printVarSummary(fmlse)

    edgx:=new((/nrows,3/),"float")
    edgy:=new((/nrows,3/),"float")
    dx=0.075
    edgx(:,0)=flon
    edgx(:,1)=flon-dx*1.732
    edgx(:,2)=flon+dx*1.732
    edgy(:,0)=flat+dx*2
    edgy(:,1)=flat-dx
    edgy(:,2)=flat+dx

    wks:=gsn_open_wks("png", "quick_looks/Quick_Look-"+yyyymmdd)
    gsn_define_colormap(wks, "WhiteBlueGreenYellowRed")

    res=True
    res@gsnFrame=False
    res@gsnDraw=False

    res@cnFillOn=True
    res@cnLinesOn = False
    res@cnLineLabelsOn = False
    res@cnInfoLabelOn = False
    res@cnFillMode = "CellFill"
    res@lbLabelBarOn   =True
    res@lbBoxEndCapStyle="TriangleBothEnds"
    res@lbOrientation="Vertical"
    res@pmLabelBarWidthF=0.08
    res@pmLabelBarOrthogonalPosF=0.01
    res@lbLabelAutoStride=False
    res@lbLabelStride=2

    res@gsnSpreadColorStart=20
    res@cnLevelSelectionMode  = "ExplicitLevels"  
    ; res@cnLevels=(/0.7,0.74,0.78,0.82,0.86,0.84,0.88,0.9,0.91,0.92,0.93,0.94,0.95,0.96,0.97,0.98,0.99/)
    res@cnLevels:=(/0.85,0.88,0.9,0.91,0.92,0.93,0.94,0.95,0.96,0.97,0.98,0.99/)

    res@mpFillOn=True
    res@gsnAddCyclic = False 
    res@mpMinLatF = -60.
    res@mpMaxLatF = 75.
    res@mpMinLonF = -180.
    res@mpMaxLonF = 180.
    res@tmXBLabelsOn=False
    ; res@mpMinLatF = min(flat)
    ; res@mpMaxLatF = max(flat)
    ; res@mpMinLonF = min(flon)-20
    ; res@mpMaxLonF = max(flon)+20

    ; res@mpFillOn = True
    res@cnFillDrawOrder      = "PostDraw"
    res@mpLandFillColor      = "lightgrey"
    res@mpOceanFillColor     = "white"  

    res@sfXArray  := flon
    res@sfYArray := flat
    res@sfXCellBounds :=edgx
    res@sfYCellBounds :=edgy
    
    channels=(/"10V","10H","18V","18H","23V","23H","36V","36H","89V","Flag"/)
    plot:=new(10,graphic)

    ps:=(/0,2,4,6,8/) ; Vpol
    do ip=0,dimsizes(ps)-1
        res@gsnLeftString="";+channels(ps(ip))
        plot(ip*2)=gsn_csm_contour_map(wks, fmlse(:,ps(ip)), res)
    end do

    res@cnLevels:=(/0.7,0.75,0.8,0.85,0.9,0.91,0.92,0.93,0.94,0.95,0.96,0.97,0.98,0.99/)
    ps:=(/1,3,5,7/) ; Hpol

    do ip=0,dimsizes(ps)-1
        ; res@gsnLeftString="GMI Emissivity, "+channels(ps(ip))
        ; if (ip.eq.3) then res@tmXBLabelsOn=True end if        
        plot(ip*2+1)=gsn_csm_contour_map(wks, fmlse(:,ps(ip)), res)
    end do

    ; res@gsnLeftString="Geostationary Satellite Flag "
    res@cnLevels:=(/0.5,1.5,2.5,6.5,7.5,8.5/)
    res@lbLabelAlignment="BoxCenters"
    res@lbBoxEndCapStyle="RectangleEnds"
    res@lbLabelStride=1

    res@lbLabelStrings=(/"Missing","H8","GOESR","MSG1","H9","MSG2","MSG3"/)
    cmap=read_colormap_file("grads_default")
    ; cmap=read_colormap_file("GMT_paired")

    res@cnFillPalette=cmap
    res@cnFillColors=(/15,3,4,6,7,8,9/)-2
    ; res@gsnSpreadColorEnd=6
    plot(9)=gsn_csm_contour_map(wks, flag, res)

    resP=True
    resP@gsnPanelBottom=0.01
    resP@gsnPanelTop=0.98
    resP@gsnPanelMainString="Quick Look "+yyyymmdd
    resP@gsnPanelFigureStrings=channels
    resP@gsnPanelFigureStringsJust="BottomLeft"
    gsn_panel(wks, plot, (/5,2/), resP)



; end do
end