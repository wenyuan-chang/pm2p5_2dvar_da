program main
implicit none
include 'netcdf.inc'

character( len=* ), parameter           :: inpath1='./DATA/PM25.OBS/'          ! obs data
character( len=* ), parameter           :: inpath2='./DATA/GEOS_Chem/pmx.nc'   ! model background data
character( len=* ), parameter           :: oupath0='./RESULT/'

real, parameter                         :: obserr=5      ! observation error
real, parameter                         :: moderr=50.    ! model backgrond error, only useful with .true. on line 284

!----------------------------------------------------------------------
!---- Gaussian Filter
integer, parameter                      :: kn=1     ! recursive number
integer, parameter                      :: nn=3     ! nn=1 & 3 for 1st & 3rd order filter; only =1 or =3 is accepted
real, parameter                         :: R=1/5.   ! distance unit, correlation radius, sigma=R/deltX
real, parameter                         :: deltX=1. ! distance unit, grid size

!---- input nc file of model
real, parameter                         :: robsmiss=-99.

integer                                 :: ksit, kday, daynum, obsnpt, modnpt, psit
integer                                 :: ksta, ncid, varid, latid, lonid, varid2
integer                                 :: dimid_tim, dimid_lev, dimid_lat, dimid_lon
integer                                 :: ntim, nlev, nlat, nlon
integer                                 :: ko, km, km1, km2, km3, km4, klon, klat
integer                                 :: kx, ky, kk, ii
integer                                 :: ifg1, ifg2, ifg3, ifg4
real                                    :: xp, yp, yn, ys, xw, xe
real                                    :: v1, v2, v3, v4, w1, w2, w3, w4, rw
real                                    :: minlon, maxlon, minlat, maxlat
real                                    :: ss
real, dimension(98)                     :: obslon, obslat, obspmx         ! site number=98
real, dimension(98,31)                  :: obsdaypm25                     ! site number=98, daysInMonth=31
real, dimension(:), allocatable         :: okobspmx, okobslon, okobslat

real, dimension(:), allocatable         :: moderr_arr
real, dimension(:), allocatable         :: modlat, modlon
real, dimension(:,:,:,:), allocatable   :: modpmx, modpmxerr
real, dimension(:,:), allocatable       :: okmodpmx, damodpmx

real, dimension(:,:), allocatable       :: HH, DD, CC, RR, PP, PP0, PP1, PP2, HHT, HC, CH, MM
real, dimension(:), allocatable         :: DIF, bv, bv0, bv1, xv, xxv

integer                                 :: NRHS, LDA, LDB, INFO
integer                                 :: timid, vtimid, vlatid, vlonid, vid1, vid2
integer, dimension(:), allocatable      :: IPIV

!------------------------------------------ INPUT MODEL BACKGROUND DATA
! read model data
ksta=nf_open( inpath2,nf_nowrite,ncid ); call nf_check(ksta,10)

!ksta=nf_inq_dimid( ncid,'time',dimid_tim ); call nf_check(ksta,13)
!ksta=nf_inq_dimid( ncid,'lev' ,dimid_lev ); call nf_check(ksta,13)
ksta=nf_inq_dimid( ncid,'lat' ,dimid_lat ); call nf_check(ksta,13)
ksta=nf_inq_dimid( ncid,'lon' ,dimid_lon ); call nf_check(ksta,13)

!ksta=nf_inq_dimlen( ncid,dimid_tim,ntim ); call nf_check(ksta,14)
!ksta=nf_inq_dimlen( ncid,dimid_lev,nlev ); call nf_check(ksta,14)
ksta=nf_inq_dimlen( ncid,dimid_lat,nlat ); call nf_check(ksta,14)
ksta=nf_inq_dimlen( ncid,dimid_lon,nlon ); call nf_check(ksta,14)

allocate( modlat(nlat), modlon(nlon) )
!allocate( modpmx(nlon,nlat,nlev,ntim) )
allocate( modpmx(nlon,nlat,1,1) )
allocate( modpmxerr(nlon,nlat,1,1) )
allocate( okmodpmx(nlon,nlat) )

ksta=nf_inq_varid( ncid,'lat' ,latid  ); call nf_check(ksta,11)
ksta=nf_inq_varid( ncid,'lon' ,lonid  ); call nf_check(ksta,11)
ksta=nf_inq_varid( ncid,'avgpmx' ,varid  ); call nf_check(ksta,11)
ksta=nf_inq_varid( ncid,'stdpmx' ,varid2  ); call nf_check(ksta,11)

ksta=nf_get_vara_real( ncid,latid,(/1/),(/nlat/),modlat ); call nf_check(ksta,121)
ksta=nf_get_vara_real( ncid,lonid,(/1/),(/nlon/),modlon ); call nf_check(ksta,122)
ksta=nf_get_vara_real( ncid,varid,(/1,1,1,1/),(/nlon,nlat,nlev,ntim/),modpmx ); call nf_check(ksta,123)
ksta=nf_get_vara_real( ncid,varid2,(/1,1,1,1/),(/nlon,nlat,nlev,ntim/),modpmxerr ); call nf_check(ksta,123)

okmodpmx=modpmx(:,:,1,1) ! jan surface-layer
minlon=modlon(1)-0.625/2.
maxlon=modlon(nlon)+0.625/2.
minlat=modlat(1)-0.5/2.
maxlat=modlat(nlat)+0.5/2.

modnpt=nlat*nlon

!------------------------------------------ INPUT OBSERVATION DATA
! read obs longitude, latitude
open(11,file=inpath1//'zlonlat.txt',status='old',action='read')
do ksit=1, 98   ! site number=98
 read(11,*) obslon(ksit), obslat(ksit)
end do
close(11)

obspmx=0.
obsnpt=0

! read obs data: pm2.5 (ug/m3) in Jan 2016
open(12,file=inpath1//'pm25.201601.txt',status='old',action='read')
do ksit=1, 98   ! site number=98
 read(12,*) (obsdaypm25(ksit,kday),kday=1,31)
 daynum=0
 do kday=1, 31
  if( obsdaypm25(ksit,kday)/=robsmiss ) then
  obspmx(ksit)=obspmx(ksit)+obsdaypm25(ksit,kday)
  daynum=daynum+1
  end if
 end do

 ! obs data check -----------------
 ifg1=0; ifg2=0; ifg3=0; ifg4=0
 if( daynum/=0 ) ifg1=1
 ! ignore the obs data outside of the model domain
 if( obslon(ksit)>=minlon .and. obslon(ksit)<=maxlon .and. obslat(ksit)>=minlat .and. obslat(ksit)<=maxlat ) ifg2=1
 if( obslon(ksit)>modlon(1) .and. obslon(ksit)<modlon(nlon) ) ifg3=1
 if( obslat(ksit)>modlat(1) .and. obslat(ksit)<modlat(nlat) ) ifg4=1

 if( ifg1*ifg2*ifg3*ifg4==1 ) then
  obsnpt=obsnpt+1
  obspmx(ksit)=obspmx(ksit)/daynum
 else
  write(*,*) 'remove site ksit=',ksit,'  lon=', obslon(ksit), '  lat=', obslat(ksit), ' need a larger domain '
  obspmx(ksit)=robsmiss
 end if
 ! obs data check -----------------
end do
close(12)

!----------------
!obsnpt=1
!----------------
allocate( okobslon(obsnpt) )
allocate( okobslat(obsnpt) )
allocate( okobspmx(obsnpt) )
!----------------
if( .True. ) then
!----------------
 psit=0
 do ksit=1, 98  ! site number=98
  if( obspmx(ksit)/=robsmiss ) then
   psit=psit+1
   okobslon(psit)=obslon(ksit)
   okobslat(psit)=obslat(ksit)
   okobspmx(psit)=obspmx(ksit)
  end if
 end do
else
  ! signle pesudo obs for test
  okobslon(1)=115.
  okobslat(1)=35.
  okobspmx(1)=500.
end if

! output the obs data that are assimilated (for figure)
open(21,file=oupath0//'obspm25.txt',status='replace',action='write')
do ko=1, obsnpt
 write(21,*) okobslon(ko), okobslat(ko), okobspmx(ko)
end do
close(21)

!------------------------------------------ 2D-Var DA
! NOTE: The data order of 2-D okmodpmx(nlon,nlat) in H(0,:) is (klat-1)*nlon+klon 
!
write(*,*) 'HH...'
allocate( HH(obsnpt,modnpt), HHT(modnpt,obsnpt) )
allocate( DIF(obsnpt) )

! interpolate model data to obs location; set HH matrix and DIF
HH=0.
DIF=0.
do ko=1, obsnpt
 do klon=1, nlon
 if( okobslon(ko)>=modlon(klon) .and. okobslon(ko)<modlon(klon+1) ) exit
 end do
 do klat=1, nlat
 if( okobslat(ko)>=modlat(klat) .and. okobslat(ko)<modlat(klat+1) ) exit
 end do

 !--------------------------
 ! v1(xw,yn)           v12         v2(xe,yn)
 ! w1=(xe-xp)/(xe-xw)              w2=(xp-xw)/(xe-xw)
 !
 !       w=(yp-ys)/(yn-ys)
 !                  vp(xp,yp)
 !    w=(yn-yp)/(yn-ys)
 !
 ! v4(xw,ys)           v34         v3(xe,ys)
 ! w4=(xe-xp)/(xe-xw)              w3=(xp-xw)/(xe-xw)
 !--------------------------
 xp=okobslon(ko)
 yp=okobslat(ko)

 yn=modlat(klat+1)
 ys=modlat(klat)
 xw=modlon(klon)
 xe=modlon(klon+1)

 ! northwest point
 kx=klon
 ky=klat+1
 v1=okmodpmx(klon,klat+1)
 w1=(xe-xp)/(xe-xw)*(yp-ys)/(yn-ys)
 km1=(ky-1)*nlon+kx
 HH(ko,km1)=w1
 km1=(kx-1)*nlat+ky

 ! northeast point
 kx=klon+1
 ky=klat+1
 v2=okmodpmx(klon+1,klat+1)
 w2=(xp-xw)/(xe-xw)*(yp-ys)/(yn-ys)
 km2=(ky-1)*nlon+kx
 HH(ko,km2)=w2
 km2=(kx-1)*nlat+ky

 ! southeast point
 kx=klon+1
 ky=klat
 v3=okmodpmx(klon+1,klat)
 w3=(xp-xw)/(xe-xw)*(yn-yp)/(yn-ys)
 km3=(ky-1)*nlon+kx
 HH(ko,km3)=w3
 km3=(kx-1)*nlat+ky

 ! southwest point
 kx=klon
 ky=klat
 v4=okmodpmx(klon,klat)
 w4=(xe-xp)/(xe-xw)*(yn-yp)/(yn-ys)
 km4=(ky-1)*nlon+kx
 HH(ko,km4)=w4
 km4=(kx-1)*nlat+ky

 ! DIF=obs-mod, innovation
 DIF(ko)=okobspmx(ko)-(w1*v1+w2*v2+w3*v3+w4*v4)
! write(*,*) DIF(ko)

!-------------------------------------------------- DEBUG
!if( .false. ) then
! write(*,'(5f10.4,"|",6f10.4,"|",5(x,i3.3),"|",f10.4)') &
! & w1, w2, w3, w4, w1+w2+w3+w4, ys, yp, yn, xw, xp, xe, ko, km1, km2, km3, km4, DIF(ko)
!end if
!-------------------------------------------------- DEBUG
end do

!-------------------------------------------------- DEBUG
!if( .false. ) then
!write(*,*) '--------- HH obsnpt*modnpt ---------'
!do ko=1, obsnpt
!write(*,'(20f6.3)') (HH(ko,km),km=1,modnpt)
!end do
!end if
!-------------------------------------------------- DEBUG

write(*,*) 'modnpt=', modnpt, '  nlat=', nlat, '  nlon=', nlon
write(*,*) 'obsnpt=', obsnpt

!--------------------------------------------------------------
!------------ BACKGROUND ERROR --------------------------------
write(*,*) 'BE...'
allocate( RR(obsnpt,obsnpt) )
RR=0.
!----------------------------- OBS ERROR
!--- I
if( .false. ) then
 do ko=1, obsnpt
  RR(ko,ko)=obserr**2  ! obs representative error^2, ug/m3
  RR(ko,ko)=1./RR(ko,ko) ! RR inverse matrix
 end do
!--- II
else
 do ko=1, obsnpt
!!!  RR(ko,ko)=sqrt(obserr**2+(0.2*okobspmx(ko))**2)  ! obs representative error^2 + instrumental error^2, ug/m3
  RR(ko,ko)=sqrt( obserr )
  RR(ko,ko)=1./RR(ko,ko) ! RR inverse matrix
 end do
end if
!----------------------------- OBS ERROR

allocate( DD(modnpt,modnpt) )
allocate( CC(modnpt,modnpt) )
allocate( moderr_arr(modnpt) )
DD=0.
!----------------------------- MODEL ERROR
!--- I
if( .true. ) then
 write(*,*) '==> use constant moderr'
 do km=1, modnpt
  DD(km,km)=sqrt( moderr )                    ! ug/m3
  moderr_arr(km)=sqrt( moderr )
 end do
!--- II
else
 write(*,*) '==> use modpmxerr; moderr is not effective.'
 do klon=1, nlon
 do klat=1, nlat
  km=(klat-1)*nlon+klon
  !km=(klon-1)*nlat+klat
  DD(km,km)=sqrt( modpmxerr(klon,klat,1,1) ) ! stdev
  moderr_arr(km)=sqrt( modpmxerr(klon,klat,1,1) )
 end do
 end do
end if
!----------------------------- MODEL ERROR

CC=DD
do km=1, modnpt
  call GaussianFilter( CC(:,km), modnpt, kn, nn, deltX, R )
end do
do km=1, modnpt
  call GaussianFilter( CC(km,:), modnpt, kn, nn, deltX, R )
end do
if( .False. ) then
do km1=1, modnpt
  rw=DD(km1,km1)/CC(km1,km1)
  do km2=1, modnpt
    CC(km1,km2)=CC(km1,km2)*rw
  end do
end do
end if
!------------ BACKGROUND ERROR --------------------------------
!--------------------------------------------------------------

allocate( MM(nlon,nlat) )
allocate( HC(obsnpt,modnpt), CH(modnpt,obsnpt) )
do ko=1, obsnpt
  do klat=1, nlat
  do klon=1, nlon
    kk=(klat-1)*nlon+klon
    MM(klon,klat)=HH(ko,kk)
    MM(klon,klat)=MM(klon,klat)*sqrt( modpmxerr(klon,klat,1,1) )
  end do
  end do

  ! smooth 
  do klat=1, nlat
   call GaussianFilter( MM(:,klat), nlon, kn, nn, deltX, R )
  end do
  do klon=1, nlon
   call GaussianFilter( MM(klon,:), nlat, kn, nn, deltX, R )
  end do
  !MM=MM*2.  !!! amplitude adjustment factor related to normalization

  do klat=1, nlat
  do klon=1, nlon
    kk=(klat-1)*nlon+klon
    HC(ko,kk)=MM(klon,klat)
  end do
  end do
end do

write(*,*) 'CH=HC...'
do ko=1, obsnpt
do km=1, modnpt
 CH(km,ko)=HC(ko,km)
end do
end do

!-------------- bv=AA*xv
write(*,*) 'bv0=RR*DIF...'
allocate( bv(modnpt), bv0(obsnpt) )
do ko=1, obsnpt
 ss=0.
 do kk=1, obsnpt
  ss=ss+RR(ko,kk)*DIF(kk)
 end do
 bv0(ko)=ss
end do
write(*,*) 'bv=CH*bv0...'
do km=1, modnpt
 ss=0.
 do kk=1, obsnpt
  ss=ss+CH(km,kk)*bv0(kk)
 end do
 bv(km)=ss
end do
deallocate( bv0 )

write(*,*) 'PP0=RR*HC...'
allocate( PP(modnpt,modnpt), PP0(obsnpt,modnpt) )
PP0=0.
do ko=1, obsnpt
do km=1, modnpt
 ss=0.
 do kk=1, obsnpt
  ss=ss+RR(ko,kk)*HC(kk,km)
 end do
 PP0(ko,km)=ss
end do
end do
write(*,*) 'PP=CH*PP0...'
PP=0.
do km=1, modnpt
do kk=1, modnpt
 ss=0.
 do ko=1, obsnpt
  ss=ss+CH(km,ko)*PP0(ko,kk)
 end do
 PP(km,kk)=ss
 if( km==kk ) PP(km,kk)=PP(km,kk)+1
end do
end do
deallocate( PP0 )

!-------------- solver
write(*,*) 'solve...'
allocate( xv(modnpt), IPIV(modnpt) )
NRHS=1
LDA=modnpt
LDB=modnpt
call sgesv( modnpt, NRHS, PP, LDA, IPIV, bv, LDB, INFO )
xv=bv
!write(*,*) ' xv=', xv

write(*,*) 'xxv=CC*xv...'
allocate( xxv(modnpt) )
xxv=0.
do ii=1, modnpt
    ss=0.
    do kk=1, modnpt
        ss=ss+CC(ii,kk)*xv(kk)
    end do
    xxv(ii)=ss
end do
deallocate( bv )

write(*,*) '==============================='
write(*,*) '           DONE                '
allocate( damodpmx(nlon,nlat) )
damodpmx=okmodpmx
do klon=1, nlon
do klat=1, nlat
 km=(klat-1)*nlon+klon
 !km=(klon-1)*nlat+klat
 damodpmx(klon,klat)=damodpmx(klon,klat)+xxv(km)
  !write(*,*) 'klat, klon, xxv=', klat, klon, xxv(km)

 if( damodpmx(klon,klat)<=0 ) then
 write(*,*) 'ERROR: pm<=0', klon,klat,damodpmx(klon,klat)
 end if
end do
end do
write(*,*) '==============================='
!-------------- bv=AA*xv, AA=I+F*PP*F, F:GaussianFilter

! output anlaysis data
ksta=nf_create("./RESULT/da.nc",NF_CLOBBER,ncid); call nf_check(ksta,13)
ksta=nf_def_dim(ncid,"time",1,timid)
ksta=nf_def_dim(ncid,"lat",nlat,latid)
ksta=nf_def_dim(ncid,"lon",nlon,lonid)

ksta=nf_def_var(ncid,"time",nf_real,1,(/timid/),vtimid)
ksta=nf_put_att_text(ncid,vtimid,'long_name',4,'time')
ksta=nf_put_att_text(ncid,vtimid,'units',29,'months since 2000-01-01 00:00')

ksta=nf_def_var(ncid,"lat",nf_real,1,(/latid/),vlatid)
ksta=nf_put_att_text(ncid,vlatid,'long_name',8,'latitude')
ksta=nf_put_att_text(ncid,vlatid,'units',13,'degrees_north')

ksta=nf_def_var(ncid,"lon",nf_real,1,(/lonid/),vlonid)
ksta=nf_put_att_text(ncid,vlonid,'long_name',9,'longitude')
ksta=nf_put_att_text(ncid,vlonid,'units',12,'degrees_east')

ksta=nf_def_var(ncid,'pm25noda',nf_real,3,(/lonid,latid,timid/),vid1)
ksta=nf_put_att_real(ncid,vid1,'missing_value',nf_real,1,-9999.)
ksta=nf_def_var(ncid,'pm25ysda',nf_real,3,(/lonid,latid,timid/),vid2)
ksta=nf_put_att_real(ncid,vid2,'missing_value',nf_real,1,-9999.)

ksta=nf_enddef(ncid)

ntim=1
ksta=nf_put_vara_real(ncid,vtimid,1,ntim,(/1./))
ksta=nf_put_vara_real(ncid,vlatid,1,nlat,modlat)
ksta=nf_put_vara_real(ncid,vlonid,1,nlon,modlon)

ksta=nf_put_var_real(ncid,vid1,okmodpmx)
ksta=nf_put_var_real(ncid,vid2,damodpmx)
ksta=nf_close(ncid)

deallocate( xv, xxv )
deallocate( okobspmx )
deallocate( okobslon )
deallocate( okobslat )
deallocate( HH )
deallocate( HHT )
deallocate( DIF )
deallocate( DD )
deallocate( CC )
deallocate( moderr_arr )
deallocate( MM )

deallocate( modlat, modlon )
deallocate( modpmx )
deallocate( modpmxerr )
deallocate( okmodpmx )

contains
!---------------------------------- SUB 1
subroutine nf_check( ierr,iflg )
implicit none
integer, intent( in )           :: ierr, iflg
if( ierr/=nf_noerr ) then; write(*,'(i4.4,a50)') iflg, ' err: '//nf_strerror(ierr); stop; endif
end subroutine nf_check

!---------------------------------- SUB 2
subroutine GaussianFilter( yarr, nmax, kn, nn, deltX, R )
implicit none
integer, intent( in )     :: nmax, kn, nn
real, dimension(nmax), intent( inout ) :: yarr
real, intent( in )      :: deltX, R

integer         :: si, sk
real         :: sigma, E, AA, BB, ai, A1, A2, A3
real, dimension(nmax)     :: pp, ss

ss=yarr
sigma=R/deltX
!-------- Fn, n=1
if( nn==1 ) then
E=1.*kn/sigma**2
AA=1+E-sqrt(E*(E+2))    ! alpha
BB=sqrt(E*(E+2))-E      ! beta
pp=0.
do sk=1, kn
    do si=1, nmax
        if( si==1 ) then
            pp(si)=BB*ss(si)
        else
            pp(si)=BB*ss(si)+AA*pp(si-1)
        end if
    end do
    do si=nmax, 1, -1
        if( si==nmax ) then
            ss(si)=BB*pp(si)
        else
            ss(si)=BB*pp(si)+AA*ss(si+1)
        end if
    end do
end do
end if
!-------- Fn, n=3
if( nn==3 ) then
ai=3.738128+5.788982*sigma+3.382473*sigma**2+sigma**3
A1=(5.788982*sigma+6.764946*sigma**2+3*sigma**3)/ai
A2=-(3.382473*sigma**2+3*sigma**3)/ai
A3=sigma**3/ai
BB=1-(A1+A2+A3)
pp=0.
do sk=1, kn
    do si=1, nmax
        if( si==1 ) then
            pp(si)=BB*ss(si)
        else if( si==2 ) then
            pp(si)=BB*ss(si)+A1*pp(si-1)
        else if( si==3 ) then
            pp(si)=BB*ss(si)+A1*pp(si-1)+A2*pp(si-2)
        else
            pp(si)=BB*ss(si)+A1*pp(si-1)+A2*pp(si-2)+A3*pp(si-3)
        end if
    end do

    do si=nmax, 1, -1
        if( si==nmax ) then
            ss(si)=BB*pp(si)
        else if( si==nmax-1 ) then
            ss(si)=BB*pp(si)+A1*ss(si+1)
        else if( si==nmax-2 ) then
            ss(si)=BB*pp(si)+A1*ss(si+1)+A2*ss(si+2)
        else
            ss(si)=BB*pp(si)+A1*ss(si+1)+A2*ss(si+2)+A3*ss(si+3)
        end if
    end do
end do
end if
yarr=ss
end subroutine GaussianFilter

end program main
