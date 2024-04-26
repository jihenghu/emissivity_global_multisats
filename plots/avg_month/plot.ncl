begin
yyyymm=201803
    fin=addfile("GMI_Emissivity_Monthly_"+yyyymm+".HDF5","r")
    mlse=fin->Emissivity
    Latitude=fin->Latitude
    Longitude=fin->Longitude

    mlse@lat2d=Latitude
    mlse@lon2d=Longitude
    mlse@_FillValue=-999.9
    mlse=where(mlse.lt.0.5,mlse@_FillValue,mlse)

    wks:=gsn_open_wks("png", "Quick_Look-"+yyyymm)
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
    res@cnFillDrawOrder      = "PreDraw"
    res@mpLandFillColor      = "transparent"
    res@mpOceanFillColor     = "white"  

    res@cnMissingValFillColor="lightgray"
    
    channels=(/"10V","10H","18V","18H","23V","36V","36H","89V","89H"/)
    plot:=new(10,graphic)

    ps:=(/0,2,4,5,7/) ; Vpol
    do ip=0,dimsizes(ps)-1
        res@gsnLeftString="";+channels(ps(ip))
        plot(ip*2)=gsn_csm_contour_map(wks, mlse(ps(ip),:,:), res)
    end do

    res@cnLevels:=(/0.7,0.75,0.8,0.85,0.9,0.91,0.92,0.93,0.94,0.95,0.96,0.97,0.98,0.99/)
    ps:=(/1,3,6,8/) ; Hpol

    do ip=0,dimsizes(ps)-1
        ; res@gsnLeftString="GMI Emissivity, "+channels(ps(ip))
        ; if (ip.eq.3) then res@tmXBLabelsOn=True end if        
        plot(ip*2+1)=gsn_csm_contour_map(wks, mlse(ps(ip),:,:), res)
    end do

    ; res@gsnLeftString="Geostationary Satellite Flag "
    res@cnLevels:=(/0.5,1.5,2.5,6.5,7.5,8.5/)
    res@lbLabelAlignment="BoxCenters"
    res@lbBoxEndCapStyle="RectangleEnds"
    res@lbLabelStride=1

    res@lbLabelStrings=(/"Missing","H8","GOESR","MSG1","H9","MSG2","MSG3"/)
    cmap=read_colormap_file("grads_default")
    ; cmap=read_colormap_file("GMT_paired")

    ; res@cnFillPalette=cmap
    ; res@cnFillColors=(/15,3,4,6,7,8,9/)-2
    ; ; res@gsnSpreadColorEnd=6
    ; ; plot(9)=gsn_csm_contour_map(wks, flag, res)

    resP=True
    resP@gsnPanelBottom=0.01
    resP@gsnPanelTop=0.98
    resP@gsnPanelMainString="Quick Look "+yyyymm
    resP@gsnPanelFigureStrings=channels
    resP@gsnPanelFigureStringsJust="BottomLeft"
    gsn_panel(wks, plot, (/5,2/), resP)



; end do
end