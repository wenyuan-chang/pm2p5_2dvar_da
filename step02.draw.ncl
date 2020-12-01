load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"

begin
	obsval=tofloat(readAsciiTable("./RESULT/obspm25.txt",3,"float",0))
	obslon=obsval(:,0)
	obslat=obsval(:,1)
	obspmx=obsval(:,2)
	obsnum=dimsizes(obspmx)

  ;----- obspmx signle
  ;obspmx:=500.
  ;----- obspmx single

	infnm=addfile("./RESULT/da.nc","r")
	
	lat1d=infnm->lat
	lon1d=infnm->lon

	noda=infnm->pm25noda(0,:,:)
	ysda=infnm->pm25ysda(0,:,:)

	dfda=ysda-noda
	copy_VarCoords(noda,dfda)

	nodapt=linint2_points( lon1d,lat1d,noda,False,obslon,obslat,0 )
	ysdapt=linint2_points( lon1d,lat1d,ysda,False,obslon,obslat,0 )
  ;print(obspmx+" "+nodapt+" "+ysdapt)
  print(obspmx+" "+nodapt+" "+ysdapt)

	;-----------------------------------
	wks=gsn_open_wks("png","./FIGURE/da")
	cmap=read_colormap_file("WhiteBlueGreenYellowRed")
	cmap_diff=read_colormap_file("BlueRed")

	cnres=True
	cnres@gsnDraw=False
	cnres@gsnFrame=False

  cnres@mpLimitMode="LatLon"
  cnres@mpMinLatF=min(lat1d)
  cnres@mpMaxLatF=max(lat1d)
  cnres@mpMinLonF=min(lon1d)
  cnres@mpMaxLonF=max(lon1d)

  cnres@mpGridAndLimbOn=False
  cnres@mpGeophysicalLineColor="black"
  cnres@mpGeophysicalLineThicknessF=2.
  cnres@mpNationalLineThicknessF=cnres@mpGeophysicalLineThicknessF
  cnres@mpDataSetName="EArth..4"
  cnres@mpDataBaseVersion="MediumRes"
  cnres@mpOutlineOn=True
  cnres@mpOutlineSpecifiers=(/"China","Taiwan"/)
  ;cnres@mpLandFillColor="white"
  ;cnres@mpFillDrawOrder="PreDraw"
  ;cnres@mpOutlineDrawOrder="Draw"
  cnres@cnFillPalette=cmap
  cnres@gsnAddCyclic=False

	cnres@cnFillMode="RasterFill"
	cnres@cnFillOn=True
  cnres@cnLinesOn=False
  cnres@cnLineLabelsOn=False

  cnres@cnLevelSelectionMode="ManualLevels"

	plots=new( (/4/),"graphic" )
	dumpt=new( (/4,obsnum/),"graphic" )
	dumbg=new( (/4/),"graphic" )

  cnres@cnMinLevelValF=0.
  cnres@cnMaxLevelValF=200.
  cnres@cnLevelSpacingF=10.

	mkres=True
	mkres@gsMarkerIndex=16

	pid=0
	cnres@gsnLeftString="noda"
	plots(pid)=gsn_csm_contour_map(wks,noda,cnres)
	mkres@gsMarkerColor:="black"
	mkres@gsMarkerSizeF=0.008
	dumbg(pid)=gsn_add_polymarker(wks,plots(pid),obslon,obslat,mkres)
	getvalues plots@contour
	"cnLevels":cnlevels
	end getvalues
	mkres@gsMarkerSizeF=mkres@gsMarkerSizeF-0.001
	do pt=0, obsnum-1
		mkres@gsMarkerColor:=get_color_rgba(cmap,cnlevels,obspmx(pt))
		dumpt(pid,pt)=gsn_add_polymarker(wks,plots(pid),obslon(pt),obslat(pt),mkres)
	end do

  cnres@cnMinLevelValF=0.
  cnres@cnMaxLevelValF=200.
  cnres@cnLevelSpacingF=10.

	pid=1
	cnres@gsnLeftString="ysda"
	plots(pid)=gsn_csm_contour_map(wks,ysda,cnres)
	mkres@gsMarkerColor:="black"
	mkres@gsMarkerSizeF=0.008
	dumbg(pid)=gsn_add_polymarker(wks,plots(pid),obslon,obslat,mkres)
	getvalues plots@contour
	"cnLevels":cnlevels
	end getvalues
	mkres@gsMarkerSizeF=mkres@gsMarkerSizeF-0.001
	do pt=0, obsnum-1
		mkres@gsMarkerColor:=get_color_rgba(cmap,cnlevels,obspmx(pt))
		dumpt(pid,pt)=gsn_add_polymarker(wks,plots(pid),obslon(pt),obslat(pt),mkres)
	end do

  cnres@cnMinLevelValF=-10.
  cnres@cnMaxLevelValF= 10.
  cnres@cnLevelSpacingF=1.

	pid=2
  cnres@cnFillPalette:=cmap_diff
	cnres@gsnLeftString="ysda-noda"
	plots(pid)=gsn_csm_contour_map(wks,dfda,cnres)
	mkres@gsMarkerColor:="black"
	mkres@gsMarkerSizeF=0.008
	dumbg(pid)=gsn_add_polymarker(wks,plots(pid),obslon,obslat,mkres)
	;getvalues plots@contour
	;"cnLevels":cnlevels
	;end getvalues
	mkres@gsMarkerSizeF=mkres@gsMarkerSizeF-0.001
	do pt=0, obsnum-1
		;mkres@gsMarkerColor:=get_color_rgba(cmap,cnlevels,obspmx(pt))
		mkres@gsMarkerColor:="black"
		dumpt(pid,pt)=gsn_add_polymarker(wks,plots(pid),obslon(pt),obslat(pt),mkres)
	end do

	xyres=True
	xyres@gsnDraw=False
	xyres@gsnFrame=False
	xyres@xyMarkLineMode="Markers"
	xyres@xyMarker=16
	xyres@xyMarkerColors=(/"blue","red"/)

  x0=max(obspmx)
  x1=max(noda)
  x2=max(ysda)
	xmax=max((/x0,x1,x2/))
	xyres@trXMinF=0.
	xyres@trXMaxF=xmax
	xyres@trYMinF=xyres@trXMinF
	xyres@trYMaxF=xyres@trXMaxF

	xyres@vpHeightF=0.45
	xyres@vpWidthF=0.45

	df1=tostring_with_format(avg(nodapt-obspmx),"%4.2f")
	df2=tostring_with_format(avg(ysdapt-obspmx),"%4.2f")

	pid=3
	xyres@gsnLeftString="noda-obs="+df1+"  ysda-obs="+df2
	xyres@tiXAxisString="OBS"
	xyres@tiYAxisString="MOD"
	plots(pid)=gsn_csm_xy(wks,(/obspmx,obspmx/),(/nodapt,ysdapt/),xyres)
	dumln=gsn_add_polyline(wks,plots(pid),(/0,xmax/),(/0,xmax/),False)

	gsn_panel(wks,plots,(/2,2/),False)
end
