#undef BL_LANG_CC
#define BL_LANG_FORT

#include "REAL.H"
#include "CONSTANTS.H"
#include "SPACE.H"
#include "BC_TYPES.H"
#include "ArrayLim.H"
!------------------------------------------------------------
module mmat_FVM
use Generalclass
use probcommon_module
use mof_routines_module
use MOF_pair_module

implicit none

REAL(KIND=8),PARAMETER   :: mof_tol=1.0e-4
integer,parameter        :: shapeflag = 0   ! asteroid test case 0: asteroid 1:diamand 
                                           ! 2: circle

contains


!---------------------------------------------------------
subroutine dist_fns(imat,x,y,dist,probtype_in)

implicit none

integer,intent(in)               :: imat,probtype_in
real(kind=8)                     :: x,y,dist
real(kind=8)                     :: dist1,dist2,dist3,dist4,dist5
real(kind=8)                     :: dist6,dist7,dist8
real(kind=8)                     :: dista,distb,distc,distd
real(kind=8)                     :: d1,d2,d3,d4,d5
real(kind=8)                     :: dd1,dd2,dd3,dd4,dd5
real(kind=8)                     :: xy(2)
real(kind=8)                     :: x1(2),x2(2),x3(2),x4(2),x5(2)
real(kind=8)                     :: x6(2)
real(kind=8)                     :: xxl(2),xxr(2)
integer                          :: i
real(kind=8)                     :: r1,r2,r3,r4
real(kind=8)                     :: center(2),cc(2)

real(kind=8)                     :: x0,y0
real(kind=8)                     :: c1,c2,c3,c4,tt,ttcrit,tt1,tt2
integer                          :: ttsign
real(kind=8)                     :: tcrit,tcrit1,tcrit2  !theta
real(kind=8)                     :: signtemp
real(kind=8)                     :: xtheta,ytheta
real(kind=8)                     :: xtheta1,ytheta1,xtheta2,ytheta2
integer                          :: flag
integer                          :: pp1,pp2
real(kind=8)                     :: pcurve_crit(2)
integer                          :: cenflag,cccflag

real(kind=8)                     :: vt1(2),vt2(2),vt3(2),vt4(2)
real(kind=8)                     :: vtd1,vtd2,vtd3,vtd4
real(kind=8)                     :: va(2),vb(2),vc(2),vd(2)
real(kind=8)                     :: vra(2),vrb(2),vrc(2),vrd(2)

real(kind=8)                     :: smrad
real(kind=8)                     :: rad_crit(2)
integer                          :: ppcrit
real(kind=8)                     :: pcrit1(2),pcrit2(2),pcrit(2)
real(kind=8)                     :: psig1(2),psig2(2)
real(kind=8)                     :: pclose(2),pclose1(2),pclose2(2)

real(kind=8)                     :: mx,my,c
real(kind=8)                     :: cdiff
real(kind=8)                     :: xysign
real(kind=8)                     :: xc(2),xc1(2),xc2(2),xc3(2)
real(kind=8)                     :: crossp(3)
integer                          :: choosepn
real(kind=8)                     :: r_inner  ! radius of inner circle


if(probtype_in .eq. 1 .or. probtype_in .eq. 9) then
 center(1) = 0.5d0
 center(2) = 0.5d0
 r1=radcen-radeps
 r2=radcen+radeps
 dist1= sqrt((x-center(1))**2.0d0+(y-center(2))**2.0d0)

if (imat .eq. 1) then
   dist= r1 - dist1  ! inner material
elseif(imat .eq. 2) then
   if(dist1 .le. r2 .and.  dist1 .ge. r1)then
      dist2 = r2 - dist1
      dist3 = dist1 - r1
      dist = min(dist2,dist3)
   else if (dist1.ge.r2) then
    dist=r2-dist1
   else if (dist1.le.r1) then
    dist=dist1-r1
   else
    print *,"dist1 bust"
    stop
   endif
elseif(imat .eq. 3)then  
   dist = dist1 - r2  ! outer material
else
   print *,"wrong number of mat"
   stop
endif

elseif(probtype_in .eq. 4)then

 x0 = 2.0d0*x-1.0d0
 y0 = 2.0d0*y-1.0d0

 c1 = 0.02d0*sqrt(5.0d0)
 c2 = 0.02d0*sqrt(5.0d0)
 
 tt = atan((y0-c2)/(x0-c1));  
 if((x0-c1) .ge. 0.0d0 .and. (y0-c2) .ge. 0.0d0)then
    ! do nothing
 elseif((x0-c1) .le. 0.0d0 .and. (y0-c2) .gt. 0.0d0)then
    tt = tt + pi;
 elseif((x0-c1) .lt. 0.0d0 .and. (y0-c2) .lt. 0.0d0)then
    tt = tt +pi;
 else
    tt = 2.0d0*pi + tt;
 endif

 dist1 = -(sqrt((x0-c1)**2.0d0 + (y0-c2)**2.0d0) - &
         (0.5d0 + 0.2d0*sin(5.0d0*tt)))
 if(imat .eq. 1)then
  dist = dist1
 elseif(imat .eq. 2) then
  dist = -dist1
 else
  print *,"wrong imat flag in, 130"
  stop
 endif



elseif(probtype_in .eq. 3) then

 x0 = 2.0d0*x-1.0d0
 y0 = 2.0d0*y-1.0d0

 c1 = 0.02d0*sqrt(5.0d0)
 c2 = 0.02d0*sqrt(5.0d0)
 
 tt = atan((y0-c2)/(x0-c1));  
 if((x0-c1) .ge. 0.0d0 .and. (y0-c2) .ge. 0.0d0)then
    ! do nothing
 elseif((x0-c1) .le. 0.0d0 .and. (y0-c2) .gt. 0.0d0)then
    tt = tt + pi;
 elseif((x0-c1) .lt. 0.0d0 .and. (y0-c2) .lt. 0.0d0)then
    tt = tt +pi;
 else
    tt = 2.0d0*pi + tt;
 endif

 dist1 = -(sqrt((x0-c1)**2.0d0 + (y0-c2)**2.0d0) - &
         (0.5d0 + 0.2d0*sin(5.0d0*tt)))
 dist2 = -(sqrt((x0-c1)**2.0d0 + (y0-c2)**2.0d0) - &
         (0.5d0 + pentaeps + 0.2d0*sin(5.0d0*tt)))

 if(imat .eq. 1)then
  dist = dist1  
 elseif(imat .eq. 3)then
  dist = -dist2
 elseif(imat .eq. 2)then
  if(dist1 .lt. 0.0d0 .and. dist2 .gt. 0.0d0)then
   dist = min(abs(dist1),abs(dist2))
  else
   dist = -min(abs(dist1),abs(dist2))
  endif
 else
  print *,"error 143"
 endif



elseif(probtype_in .eq. 25)then     ! asteroid 2 materials connected points
 xy(1)=x
 xy(2)=y

! center(1) = 0.02d0*sqrt(5.0d0)
! center(2) = 0.02d0*sqrt(5.0d0)
 center = 0.0d0

 cc(1)=(center(1)+1.0d0)/2.0d0
 cc(2)=(center(2)+1.0d0)/2.0d0
!  print *,"xy",xy
 call rad_cal(xy,cc,tt)
 if(tt .ge. 0.0d0 .and. tt .le. 0.5d0*pi)then
   pp1=1
   pp2=pcurve_num/4+1
 elseif(tt .le. pi)then
   pp1=pcurve_num/4+1
   pp2=pcurve_num/2+1
 elseif(tt .le. 1.5d0*pi)then
   pp1= pcurve_num/2+1
   pp2=pcurve_num/4*3+1
 elseif(tt .lt. 2.0d0*pi)then
   pp1=pcurve_num/4*3+1
   pp2=pcurve_num+1
 else
  print *,"invalid tt", tt
  stop
 endif 


 dist1=100000.0d0
 ppcrit=pp1
 do i=pp1,pp2
  call l2normd(2,xy,pcurve_ls(:,i),dist2)
  if(dist2 .lt. dist1)then
!   pcurve_crit=pcurve_ls(:,i)
   dist1=dist2
   ppcrit=i
  endif
 enddo
 !  pcrit1 >> pcrit >> pcrit2
 ! print *,"ppcrit=",ppcrit
 if(ppcrit .eq. pp1) then
  pcrit(:)=pcurve_ls(:,pp1)
  pcrit2(:)=pcurve_ls(:,pp1+1)
 elseif(ppcrit .eq. pp2)then
  pcrit(:)=pcurve_ls(:,pp2)  
  pcrit1(:)=pcurve_ls(:,ppcrit-1)
 elseif(ppcrit .gt. pp1 .and. ppcrit .lt. pp2)then
  pcrit(:)=pcurve_ls(:,ppcrit)
  pcrit1(:)=pcurve_ls(:,ppcrit-1)
  pcrit2(:)=pcurve_ls(:,ppcrit+1)   
!  print *,"pcrit", pcrit
!  print *,"pcrit1",pcrit1
!  print *,"pcrit2",pcrit2
 else
  print *,"ppcrit invalid 212"
  stop
 endif

 if(ppcrit .eq. pp1) then
!  call l2normd(2,xy,pcurve_ls(:,pp1),dist1)
  call dist_point_to_line_modify(2,pcrit,pcrit2,xy,dist1,pclose)

 elseif(ppcrit .eq. pp2)then
!  call l2normd(2,xy,pcurve_ls(:,pp2),dist1)
  call dist_point_to_line_modify(2,pcrit1,pcrit,xy,dist1,pclose)

 elseif(ppcrit .gt. pp1 .and. ppcrit .lt. pp2)then
   call dist_point_to_line_modify(2,pcrit1,pcrit,xy,dist3,pclose1)
   call dist_point_to_line_modify(2,pcrit2,pcrit,xy,dist4,pclose2)

!  print *,"pp1,pp2",pp1,pp2
!  print *,"pcrit", pcrit, dist1
!  print *,"pcrit1", pcrit1, dist3
!  print *,"pcirt2", pcrit2, dist4
!  print *,"xy",xy



  if(dist3 .lt. 0.0d0 .or. dist4 .lt. 0.0d0)then
   print *,"check sign 233"
   stop
  endif

   if(dist3 .le. dist4 )then
    dist1=dist3
    pclose=pclose1
  !  write(51,*) xy,xc1
    choosepn=-1
   elseif(dist4 .lt. dist3)then
    dist1=dist4 
    pclose=pclose2
 !   write(51,*) xy,xc2
    choosepn=+1     
   else
    print *, "err, check251",dist3,dist4,dist1
    stop
   endif
 else
   print *,"ppcrit invalid 224"
   stop
 endif

!  print *,  "dist", dist1

if(pclose(1) .eq. x .and. pclose(2) .eq. y)then
 dist=0.0d0
 crossp=0.0d0
else
 if(ppcrit .eq. pp2)then
  call cross_product(2,pcurve_ls(:,ppcrit-1),&
                     pcurve_ls(:,ppcrit),xy,crossp)
 elseif(ppcrit .eq. pp1)then
  call cross_product(2,pcurve_ls(:,ppcrit), & 
                   pcurve_ls(:,ppcrit+1),xy,crossp)
 elseif(ppcrit .gt. pp1 .and. ppcrit .lt. pp2)then
  if(choosepn .eq. -1)then
   call cross_product(2,pcurve_ls(:,ppcrit-1), & 
                   pcurve_ls(:,ppcrit),xy,crossp)  
  elseif(choosepn .eq. +1)then
    call cross_product(2,pcurve_ls(:,ppcrit), & 
                   pcurve_ls(:,ppcrit+1),xy,crossp)
  else
   print *,"check chosepn"
   stop
  endif
   
 else
  print *,"check ppcrit"
  stop
 endif

 if(crossp(3) .eq. 0.0d0)then
  dist=0.0d0
  crossp=0.0d0
!  print *,"cross_product is 0 check!"
!  print *,"xy",xy
!  print *, pcurve_ls(:,ppcrit-1)
!  print *,pcurve_ls(:,ppcrit)
!  print *,pcurve_ls(:,ppcrit+1)
!  print *,"crossp",crossp
!  stop
 endif


 if(imat .eq. 1)then
  dist=sign(dist1,crossp(3))
 elseif(imat .eq. 2)then
  dist=-1.0d0*sign(dist1,crossp(3))
 else
  print *,"wrong num of materials for test 5"
  stop
 endif

endif

elseif(probtype_in .eq. 15)then     ! diamand
 xy(1)=x
 xy(2)=y
 vt1(1)=0.75d0
 vt1(2)=0.5d0
 vt2(1)=0.5d0
 vt2(2)=0.75d0
 vt3(1)=0.25d0
 vt3(2)=0.5d0
 vt4(1)=0.5d0
 vt4(2)=0.25d0
 call dist_point_to_lined(2,vt1,vt2,xy,dist1)
 call dist_point_to_lined(2,vt2,vt3,xy,dist2) 
 call dist_point_to_lined(2,vt3,vt4,xy,dist3)
 call dist_point_to_lined(2,vt4,vt1,xy,dist4)

 dist=min(dist1,dist2,dist3,dist4)

 d1=x+y-1.25d0
 d2=x+y-0.75d0
 d3=x-y-0.25d0
 d4=x-y+0.25d0

 if(d1 .le. 0.0d0 .and. d2 .ge. 0.0d0 .and. d3 .le. 0.0d0 .and. d4 .ge. 0.0d0)then
  ! do nothing
 else
  dist = -dist
 endif

 if(imat .eq. 1)then
  ! do nothing
 elseif(imat .eq. 2)then
  dist=-dist
 endif
 
elseif(probtype_in .eq. 35)then    ! whole circle
  xy(1)=x
  xy(2)=y
  cc(1)=0.5d0
  cc(2)=0.5d0
  call l2normd(2,xy,cc,dist2)
  if(imat .eq. 1)then
   dist=0.25d0-dist2
  elseif(imat .eq. 2)then
   dist=dist2-0.25d0
  else
   print *,"Err"
   stop
  endif



elseif(probtype_in .eq. 5)then     ! asteroid 2 materials back
 xy(1)=x
 xy(2)=y

! center(1) = 0.02d0*sqrt(5.0d0)
! center(2) = 0.02d0*sqrt(5.0d0)
 center = 0.0d0

 cc(1)=(center(1)+1.0d0)/2.0d0
 cc(2)=(center(2)+1.0d0)/2.0d0
!  print *,"xy",xy
 call rad_cal(xy,cc,tt)
 if(tt .ge. 0.0d0 .and. tt .le. 0.5d0*pi)then
   pp1=1
   pp2=pcurve_num/4+1
 elseif(tt .le. pi)then
   pp1=pcurve_num/4+1
   pp2=pcurve_num/2+1
 elseif(tt .le. 1.5d0*pi)then
   pp1= pcurve_num/2+1
   pp2=pcurve_num/4*3+1
 elseif(tt .le. 2.0d0*pi)then
   pp1=pcurve_num/4*3+1
   pp2=pcurve_num+1
 else
  print *,"invalid tt", tt
  stop
 endif 


 dist1=100000.0d0
! smrad=100000.0d0
 ppcrit=pp1
 do i=pp1,pp2
  call l2normd(2,xy,pcurve_ls(:,i),dist2)
  if(dist2 .lt. dist1)then
   pcurve_crit=pcurve_ls(:,i)
   dist1=dist2
   ppcrit=i
  endif
!  if(abs(pcurve_rad(i)-tt) .lt. smrad)then
!   smrad=abs(pcurve_rad(i)-tt)
!   rad_crit= pcurve_ls(:,i)   
!  endif
 enddo
 

 !------------------------------------------------
  if(ppcrit .eq. pp2)then
   ! do nothing
  else
   call dist_point_to_lined(2,pcurve_ls(:,ppcrit), &
                           pcurve_ls(:,ppcrit+1),xy,dist5)
!   write(13,*) "dist5",dist5,dist1
  endif

 ! -------------------------------------------------
if(pcurve_ls(1,ppcrit) .eq. x .and. pcurve_ls(2,ppcrit) .eq. y)then
 dist=0.0d0
 crossp=0.0d0
else
 if(ppcrit .eq. pp2)then
  call cross_product(2,pcurve_ls(:,ppcrit-1),&
                     pcurve_ls(:,ppcrit),xy,crossp)
 else
  call cross_product(2,pcurve_ls(:,ppcrit), & 
                   pcurve_ls(:,ppcrit+1),xy,crossp) 
 endif

 if(crossp(3) .eq. 0.0d0)then
  dist=0.0d0
  crossp=0.0d0
!  print *,"cross_product is 0 check!"
!  print *,"xy",xy
!  print *, pcurve_ls(:,ppcrit-1)
!  print *,pcurve_ls(:,ppcrit)
!  print *,pcurve_ls(:,ppcrit+1)
!  print *,"crossp",crossp
!  stop
 endif

 if(imat .eq. 1)then
!  dist=sign(dist1,crossp(3))
  dist=sign(min(dist1,dist5),crossp(3))
 elseif(imat .eq. 2)then
! dist=-1.0d0*sign(dist1,crossp(3))
   dist=-1.0d0*sign(min(dist1,dist5),crossp(3))
 else
  print *,"wrong num of materials for test 5"
  stop
 endif

 endif


elseif(probtype_in .eq. 10)then
 r_inner=0.1d0
 xy(1)=x
 xy(2)=y
 cenflag=0    ! 0= axis symmetry aligned with grid
              ! 1= not aligned with grid 
 cccflag =0
 if(cenflag .eq. 1)then 
  center(1) = 0.02d0*sqrt(5.0d0)
  center(2) = 0.02d0*sqrt(5.0d0)
 elseif(cenflag .eq. 0)then
  center(1) = 0.0d0
  center(2) = 0.0d0
 else
  print *,"cenflag error"
  stop
 endif

 cc(1)=(center(1)+1.0d0)/2.0d0
 cc(2)=(center(2)+1.0d0)/2.0d0


   vra(1)=cc(1)+r_inner
   vra(2)=cc(2)
   vrb(1)=cc(1)
   vrb(2)=cc(2)+r_inner
   vrc(1)=cc(1)-r_inner
   vrc(2)=cc(2)
   vrd(1)=cc(1)
   vrd(2)=cc(2)-r_inner

   va(1)=cc(1)+0.4d0
   va(2)=cc(2)
   vb(1)=cc(1)
   vb(2)=cc(2)+0.4d0
   vc(1)=cc(1)-0.4d0
   vc(2)=cc(2)
   vd(1)=cc(1)
   vd(2)=cc(2)-0.4d0


 vt1(1)=0.6d0
 vt1(2)=0.5d0
 vt2(1)=0.5d0
 vt2(2)=0.6d0
 vt3(1)=0.4d0
 vt3(2)=0.5d0
 vt4(1)=0.5d0
 vt4(2)=0.4d0
 call dist_point_to_lined(2,vt1,vt2,xy,dd1)
 call dist_point_to_lined(2,vt2,vt3,xy,dd2) 
 call dist_point_to_lined(2,vt3,vt4,xy,dd3)
 call dist_point_to_lined(2,vt4,vt1,xy,dd4)

 dist6=min(dd1,dd2,dd3,dd4)

 d1=x+y-1.1d0
 d2=x+y-0.9d0
 d3=x-y-0.1d0
 d4=x-y+0.1d0

 if(d1 .le. 0.0d0 .and. d2 .ge. 0.0d0 .and. d3 .le. 0.0d0 .and. d4 .ge. 0.0d0)then
  ! do nothing
 else
  dist6 = -dist6
 endif






 if(cenflag .eq. 0)then              ! center aligned with grid
  if(abs(xy(1) - cc(1)) .lt. 1.0e-12 &
     .and. abs(xy(2)-cc(2)) .lt. 1.0e-12)then
   cccflag=1    
  else
   call rad_cal(xy,cc,tt)
  endif 
 elseif(cenflag .eq. 1)then    ! center not aligned with grid
  call rad_cal(xy,cc,tt)
 endif

 if(cccflag .eq. 0)then


 if(tt .ge. 0.0d0 .and. tt .le. 0.5d0*pi)then
   pp1=1
   pp2=pcurve_num/4+1
 elseif(tt .le. pi)then
   pp1=pcurve_num/4+1
   pp2=pcurve_num/2+1
 elseif(tt .le. 1.5d0*pi)then
   pp1= pcurve_num/2+1
   pp2=pcurve_num/4*3+1
 elseif(tt .le. 2.0d0*pi)then
   pp1=pcurve_num/4*3+1
   pp2=pcurve_num+1
 else
  print *,"invalid tt", tt
  stop
 endif 


 dist1=100000.0d0
! smrad=100000.0d0
 ppcrit=pp1
 do i=pp1,pp2
  call l2normd(2,xy,pcurve_ls(:,i),dist2)
  if(dist2 .lt. dist1)then
   pcurve_crit=pcurve_ls(:,i)
   dist1=dist2
   ppcrit=i
  endif
!  if(abs(pcurve_rad(i)-tt) .lt. smrad)then
!   smrad=abs(pcurve_rad(i)-tt)
!   rad_crit= pcurve_ls(:,i)   
!  endif
 enddo

if(pcurve_ls(1,ppcrit) .eq. x .and. pcurve_ls(2,ppcrit) .eq. y)then
 dist=0.0d0
 crossp=0.0d0
else
 if(ppcrit .eq. pp2)then
  call cross_product(2,pcurve_ls(:,ppcrit-1),&
                     pcurve_ls(:,ppcrit),xy,crossp)
 else
  call cross_product(2,pcurve_ls(:,ppcrit), & 
                   pcurve_ls(:,ppcrit+1),xy,crossp) 
 endif

 if(crossp(3) .eq. 0.0d0)then
  dist=0.0d0
  crossp=0.0d0
!  print *,"cross_product is 0 check!"
!  print *,"xy",xy
!  print *, pcurve_ls(:,ppcrit-1)
!  print *,pcurve_ls(:,ppcrit)
!  print *,pcurve_ls(:,ppcrit+1)
!  print *,"crossp",crossp
!  stop
 endif
endif

  if(ppcrit .eq. pp2)then
   ! do nothing
  else
   call dist_point_to_lined(2,pcurve_ls(:,ppcrit), &
                           pcurve_ls(:,ppcrit+1),xy,dist1)
!   write(13,*) "dist5",dist5,dist1
  endif



   call dist_point_to_lined(2,va,vra,xy,dist2)
   call dist_point_to_lined(2,vb,vrb,xy,dist3)
   call dist_point_to_lined(2,vc,vrc,xy,dist4)
   call dist_point_to_lined(2,vd,vrd,xy,dist5)

   call dist_point_to_lined(2,vra,vrb,xy,dista)
   call dist_point_to_lined(2,vrb,vrc,xy,distb)
   call dist_point_to_lined(2,vrc,vrd,xy,distc)
   call dist_point_to_lined(2,vrd,vra,xy,distd)

     if(dist2.lt.0.0d0 .or. dist3 .lt. 0.0d0 &
        .or. dist4.lt. 0.0d0 .or. dist5.lt. 0.0d0 )then
      print *,"check dist2345"
      stop
     endif


!   first qudrant>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

 if(tt .eq. 0.0d0)then
  if(imat .eq. 1)then
    if(crossp(3) .eq. 0.0d0)then
     dist=0.0d0
    else
     dist=sign(dist1,-crossp(3))
    endif
  elseif(imat .eq. 2)then
   if(xy(1) .lt. 0.5d0+r_inner)then
    dist=-dista
   elseif(xy(1) .gt. 0.9d0)then
    dist=-(xy(1)-0.9d0)
   else
    dist=0.0d0
   endif
  elseif(imat .eq. 3)then
   dist=-distb
  elseif(imat .eq. 4)then
   dist=-distc
  elseif(imat .eq. 5)then
   if(xy(1) .lt. 0.5d0+r_inner)then
    dist=-dista
   elseif(xy(1) .gt. 0.9d0)then
    dist=-(xy(1)-0.9d0)
   else
    dist=0.0d0
   endif
  elseif(imat .eq. 6)then
   dist=dist6
  else
    print *,"wrong material num"
    stop
  endif  

  elseif(tt .gt. 0.0d0 .and. tt .lt. 0.5d0*pi)then

   if(imat .eq. 1)then    
    if(crossp(3) .eq. 0.0d0)then
     dist=0.0d0
    else
     dist=sign(dist1,-crossp(3))
    endif
   elseif(imat .eq. 2)then
     if(crossp(3) .ge. 0.0d0 .and. dist6 .le. 0.0d0)then
     dist = min(dist1,dist2,dist3,-dist6)
     if(dist1 .lt. 0.0d0 .or.  -dist6.lt.0.0d0)then
      print *,"check dist1,6"
      stop
     endif
    else
     dist= -min(dist1,-dist6)
    endif
   elseif(imat .eq. 3)then
    dist=-min(dist3,distb)
   elseif(imat .eq. 4)then
    dist=-distc
   elseif(imat .eq. 5)then
    dist=-dist2
   elseif(imat .eq. 6)then
    dist=dist6
   else
    print *,"wrong material num"
    stop
   endif

 elseif(tt .eq. 0.5d0*pi)then 

  if(imat .eq. 1)then
    if(crossp(3) .eq. 0.0d0)then
     dist=0.0d0
    else
     dist=sign(dist1,-crossp(3))
    endif
  elseif(imat .eq. 2)then
   if(xy(2) .lt. 0.5d0+r_inner)then
    dist=-dista
   elseif(xy(2) .gt. 0.9d0)then
    dist=-(xy(2)-0.9d0)
   else
    dist=0.0d0
   endif
  elseif(imat .eq. 3)then
   if(xy(2) .lt. 0.5d0+r_inner)then
    dist=-distb
   elseif(xy(2) .gt. 0.9d0)then
    dist=-(xy(2)-0.9d0)
   else
    dist=0.0d0
   endif
  elseif(imat .eq. 4)then
   dist=-min(distc,dist4)
  elseif(imat .eq. 5)then
   dist=-min(dist2,distd)
  elseif(imat .eq. 6)then
   dist=dist6
  else
    print *,"wrong material num"
    stop
  endif  
  ! second quadrant >>>>>>>>>>>>>>>>>>>>>>>>>>>>>
 elseif(tt .lt. pi)then

   if(imat .eq. 1)then
    if(crossp(3) .eq. 0.0d0)then
     dist=0.0d0
    else
     dist=sign(dist1,-crossp(3))
    endif
   elseif(imat .eq. 2)then
    dist=-min(dist3,dista)
   elseif(imat .eq. 3)then
    if(crossp(3) .gt. 0.0d0 .and. dist6 .lt. 0.0d0)then
     dist = min(dist1,dist3,dist4,-dist6)
     if(dist1 .lt. 0.0d0 .or.  -dist6.lt.0.0d0)then
      print *,"check dist1,6"
      stop
     endif
    else
     dist= - min(dist1,-dist6)
    endif
   elseif(imat .eq. 4)then
    dist=-min(dist4,distc)
   elseif(imat .eq. 5)then 
    dist=-distd
   elseif(imat .eq. 6)then
    dist=dist6
   else
    print *,"wrong material num"
    stop
   endif
 elseif(tt .eq. pi)then
  if(imat .eq. 1)then
    if(crossp(3) .eq. 0.0d0)then
     dist=0.0d0
    else
     dist=sign(dist1,-crossp(3))
    endif
  elseif(imat .eq. 2)then
    dist=-min(dist3,dista)
  elseif(imat .eq. 3)then
   if(xy(1) .gt. 0.5d0-r_inner)then
    dist=-distb
   elseif(xy(1) .lt. 0.1d0)then
    dist=xy(1)-0.1d0
   else
    dist=0.0d0
   endif

  elseif(imat .eq. 4)then
   if(xy(1) .gt. 0.5d0-r_inner)then
    dist=-distc
   elseif(xy(1) .lt. 0.1d0)then
    dist=xy(1)-0.1d0
   else
    dist=0.0d0
   endif
  elseif(imat .eq. 5)then
   dist=-distd
  elseif(imat .eq. 6)then
   dist=dist6
  else
    print *,"wrong material num"
    stop
  endif 

  ! third quadrant
 elseif(tt .lt. 1.5d0*pi)then

   if(imat .eq. 1)then
    if(crossp(3) .eq. 0.0d0)then
     dist=0.0d0
    else
     dist=sign(dist1,-crossp(3))
    endif
   elseif(imat .eq. 2)then
    dist=-dista
   elseif(imat .eq. 3)then
    dist=-min(dist4,distb)

   elseif(imat .eq. 4)then

    if(crossp(3) .gt. 0.0d0 .and. dist6 .lt. 0.0d0)then
     dist = min(dist1,dist5,dist4,-dist6)
     if(dist1 .lt. 0.0d0 .or.  -dist6.lt.0.0d0)then
      print *,"check dist1,6"
      stop
     endif
    else
     dist= - min(dist1,-dist6)
    endif
   elseif(imat .eq. 5)then
    dist=-min(dist5,distd)
   elseif(imat .eq. 6)then
    dist=dist6
   else
    print *,"wrong material num"
    stop
   endif
 
 elseif( tt .eq. 1.5d0*pi) then

  if(imat .eq. 1)then
    if(crossp(3) .eq. 0.0d0)then
     dist=0.0d0
    else
     dist=sign(dist1,-crossp(3))
    endif
  elseif(imat .eq. 2)then
   dist = -min(dist2,dista)
  elseif(imat .eq. 3)then
   dist = -min(distb,dist4)

   elseif(imat .eq. 4)then

   if(xy(2) .gt. 0.5d0-r_inner)then
    dist=-distc
   elseif(xy(2) .lt. 0.1d0)then
    dist= xy(2)-0.1d0
   else
    dist=0.0d0
   endif

  elseif(imat .eq. 5)then
   if(xy(2) .gt. 0.5d0-r_inner)then
    dist=-distd
   elseif(xy(2) .lt. 0.1d0)then
    dist= xy(2)-0.1d0
   else
    dist=0.0d0
   endif
  elseif(imat .eq. 6)then
   dist=dist6
  else
    print *,"wrong material num"
    stop
  endif  


 elseif(tt .lt. 2.0d0*pi)then


   if(imat .eq. 1)then
    if(crossp(3) .eq. 0.0d0)then
     dist=0.0d0
    else
     dist=sign(dist1,-crossp(3))
    endif
   elseif(imat .eq. 2)then
    dist=-min(dista,dist2)

   elseif(imat .eq. 3)then
    
    dist=-distb
 
   elseif(imat .eq. 4)then
     dist=-min(distc,dist5)

   elseif(imat .eq. 5)then  

    if(crossp(3) .gt. 0.0d0 .and. dist6 .lt. 0.0d0)then
     dist = min(dist1,dist2,dist5,-dist6)
     if(dist1 .lt. 0.0d0 .or.  -dist6.lt.0.0d0)then
      print *,"check dist1,6"
      stop
     endif
    else
     dist= -min(dist1,-dist6)
    endif
   elseif(imat .eq. 6)then
    dist=dist6
   else
    print *,"wrong material num"
    stop
   endif
 else
  print *,"invalid tt", tt
  stop
 endif  
 elseif(cccflag .eq. 1)then
  if(imat .eq. 1)then
   dist=-0.1d0*sqrt(2.0)
  elseif(imat .eq. 2 .or. imat .eq. 3 .or. imat .eq. 4 &
           .or. imat .eq. 5)then
   dist=-0.1d0/sqrt(2.0)
  else
   dist=+0.1d0/sqrt(2.0)
  endif
 else
  print *,"cccflag invald 353"
  stop
 endif



elseif(probtype_in .eq. 17)then           ! backup for type 7
 xy(1)=x
 xy(2)=y
 cenflag=0    ! 0= axis symmetry aligned with grid
              ! 1= not aligned with grid
 cccflag=0    ! 0= sigular at center
              ! 1= nonsigular
 if(cenflag .eq. 1)then 
  center(1) = 0.02d0*sqrt(5.0d0)
  center(2) = 0.02d0*sqrt(5.0d0)
 elseif(cenflag .eq. 0)then
  center(1) = 0.0d0
  center(2) = 0.0d0
 else
  print *,"cenflag error"
  stop
 endif

 cc(1)=(center(1)+1.0d0)/2.0d0
 cc(2)=(center(2)+1.0d0)/2.0d0
!  print *,"center",cc
 if(cenflag .eq. 0)then              ! center aligned with grid
  if(abs(xy(1) - cc(1)) .lt. 1.0e-12 &
     .and. abs(xy(2)-cc(2)) .lt. 1.0e-12)then
   cccflag=1    
  else
   call rad_cal(xy,cc,tt)
  endif 
 elseif(cenflag .eq. 1)then    ! center not aligned with grid
  call rad_cal(xy,cc,tt)
 endif

 if(cccflag .eq. 0)then

  if(tt .ge. 0.0d0 .and. tt .le. 0.5d0*pi)then
   pp1=1
   pp2=pcurve_num/4+1
  elseif(tt .le. pi)then
   pp1=pcurve_num/4+1
   pp2=pcurve_num/2+1
  elseif(tt .le. 1.5d0*pi)then
   pp1= pcurve_num/2+1
   pp2=pcurve_num/4*3
  elseif(tt .lt. 2.0d0*pi)then
   pp1=pcurve_num/4*3+1
   pp2=pcurve_num+1
  else
   print *,"invalid tt", tt
   stop
  endif 

  dist1=100000.0d0
  smrad=100000.0d0
  do i=pp1,pp2
   call l2normd(2,xy,pcurve_ls(:,i),dist2)
   if(dist2 .lt. dist1)then
    pcurve_crit=pcurve_ls(:,i)
    dist1=dist2
   endif
   if(abs(pcurve_rad(i)-tt) .lt. smrad)then
    smrad=abs(pcurve_rad(i)-tt)
    rad_crit= pcurve_ls(:,i)   
   endif
  enddo

  call l2normd(2,cc,rad_crit,dist3) 
  call l2normd(2,cc,xy,dist2)

  x1(1)=cc(1)+0.4d0                            !          x2  
  x1(2)=cc(2)                                  ! 
  x2(1)=cc(1)                                  ! 
  x2(2)=cc(2)+0.4d0                            !    x3           x1
  x3(1)=cc(1)-0.4d0                            ! 
  x3(2)=cc(2)                                  ! 
  x4(1)=cc(1)                                  !           x4 
  x4(2)=cc(2)-0.4d0                            !  



 if(imat .eq. 1)then
  dist=-sign(dist1,(dist3-dist2))
 elseif(imat .eq. 2)then
  call dist_point_to_lined(2,cc,x1,xy,dist4)
  call dist_point_to_lined(2,cc,x2,xy,dist5)
  
  if(xy(1) .lt. cc(1) .or. xy(2) .lt. cc(2))then
   dist= -min(dist4,dist5)
  else
   dist= min(sign(dist1,(dist3-dist2)),dist4,dist5)
  endif

 elseif(imat .eq. 3)then
  call dist_point_to_lined(2,cc,x2,xy,dist4)
  call dist_point_to_lined(2,cc,x3,xy,dist5)
  
  if(xy(1) .gt. cc(1) .or. xy(2) .lt. cc(2))then
   dist= -min(dist4,dist5)
  else
   dist= min(sign(dist1,(dist3-dist2)),dist4,dist5)
  endif

 elseif(imat .eq. 4)then
  call dist_point_to_lined(2,cc,x3,xy,dist4)
  call dist_point_to_lined(2,cc,x4,xy,dist5)
  
  if(xy(1) .gt. cc(1) .or. xy(2) .gt. cc(2))then
   dist= -min(dist4,dist5)
  else
   dist= min(sign(dist1,(dist3-dist2)),dist4,dist5)
  endif

 elseif(imat .eq. 5)then
  call dist_point_to_lined(2,cc,x4,xy,dist4)
  call dist_point_to_lined(2,cc,x1,xy,dist5)
  
   if(xy(1) .lt. cc(1) .or. xy(2) .gt. cc(2))then
    dist= -min(dist4,dist5)
   else
    dist= min(sign(dist1,(dist3-dist2)),dist4,dist5)
   endif
  else
   print *, "wrong number of materials 368"
   stop
  endif

 elseif(cccflag .eq. 1)then
  if(imat .eq. 1)then
   dist=-0.1d0*sqrt(2.0)
  else
   dist=0.0d0
  endif
 else
  print *,"cccflag invald 353"
  stop
 endif




elseif(probtype_in .eq. 0) then

  if(imat .eq. 1)then
    dist = 0.3d0 - y
  elseif(imat .eq. 2)then
    dist = y - 0.3d0
  else
   print *,"imat invalid imat=",imat
   print *,"x,y= ",x,y
   stop
  endif

elseif(probtype_in .eq. 2) then

  if(imat .eq. 1)then
    dist = 0.3d0 - x
  elseif(imat .eq. 2)then
    dist = x - 0.3d0
  else
   print *,"imat invalid imat=",imat
   print *,"x,y= ",x,y
   stop
  endif

elseif(probtype_in .eq. 6)then      ! nucleate boiling set-up
 xy(1)=x
 xy(2)=y
! dist1 
 x1(1)=0.0d0
 x1(2)=0.55d0
 x2(1)=0.21d0
 x2(2)=0.55d0
 x3(1)=0.79d0
 x3(2)=0.55d0
 x4(1)=1.0d0 
 x4(2)=0.55d0
 x5(1)=0.5d0
 x5(2)=0.55d0

 xxl(2)=(0.55d0**2.0+0.5d0**2.0-(0.29d0-thermal_delta)**2.0d0)/1.1d0
 xxr(2)=(0.55d0**2.0+0.5d0**2.0-(0.29d0-thermal_delta)**2.0d0)/1.1d0
 xxl(1)= 0.5d0-sqrt(0.5d0**2.0d0-xxl(2)**2.0)
 xxr(1)= 0.5d0+sqrt(0.5d0**2.0d0-xxr(2)**2.0)

 if(y .gt. 0.55d0)then
  call dist_point_to_line(x1,x2,xy,d1)
  call dist_point_to_line(x3,x4,xy,d2)
  dist1 = -min(d1,d2)
 else
  call dist_point_to_line(x1,x2,xy,d1)
  call dist_point_to_line(x3,x4,xy,d2)
  call l2normd(2,xy,x5, r1)
  d3=r1-0.29d0
  dist1=sign(min(abs(d3),abs(d1),abs(d2)),d3)
 endif 

! dist2
  x6(1)=0.5d0
  x6(2)=0.0d0
  call l2normd(2,xy,x6,r2)
  dist2= 0.5d0-r2
  
! dist3
  call l2normd(2,xy,x5,r3)
  dist3= (0.29d0-thermal_delta)-r3

 if(imat .eq. 1)then
  dist=dist1
 elseif(imat .eq. 2)then
  call dist_point_to_arc(xy,xxr,xxl,x6,d4)
  call dist_point_to_arc(xy,xxl,xxr,x5,d5)

  if(dist2 .ge. 0.0d0 .and. dist3 .ge. 0.0d0 &
           .and. dist1 .le. 0.0d0)then
   dist=min(d4,d5)
  else
   dist=-min(d4,d5)
  endif
 elseif(imat .eq. 3)then
  call dist_point_to_arc(xy,xxr,xxl,x6,d4)
  call dist_point_to_arc(xy,xxl,xxr,x5,d5)
  dist4 = min(d4,d5,abs(dist1))
  if(dist4 .lt. 0.0)then
   print *,"error, check dist4"
   stop
  endif
  if(dist2 .ge. 0.0d0 .and. dist3 .ge. 0.0d0 &
           .and. dist1 .le. 0.0d0)then
   dist=-dist4
  elseif(dist1 .gt. 0.0d0)then
   dist=-dist4
  else
   dist=dist4
  endif
 else
  print *,"probtype6 imat invalid"
  stop
 endif
! dist2
!  x6(1)=0.5d0
!  x6(2)=0.0d0
!  call l2normd(2,xy,x6,r2)
!  dist2= 0.4d0-r2
  
! dist3
!  call l2normd(2,xy,x5,r3)
!  dist3= r3-(0.4d0-thermal_delta)

! if(imat .eq. 1)then
!  dist=dist1
! elseif(imat .eq. 2)then
!  if(dist2 .ge. 0.0d0 .and. dist3 .le. 0.0d0)then
!   dist=min(abs(dist2),abs(dist3))
!  else
!   dist=-min(abs(dist2),abs(dist3))
!  endif
! elseif(imat .eq. 3)then
!  if(dist2 .ge. 0.0d0 .and. dist3 .ge. 0.0d0 &
!           .and. dist1 .le. 0.0d0)then
!   dist=min(abs(dist3),abs(dist1))
!  else
!   dist=-min(abs(dist3),abs(dist1))
!  endif
! elseif(imat .eq. 4)then
!  dist= -max(dist1,dist2)
! endif
elseif(probtype_in .eq. 8)then



 elseif(probtype_in .eq. 20)then  ! old hypercycloid
  
 xy(1)=x
 xy(2)=y
 center(1) = 0.02d0*sqrt(5.0d0)
 center(2) = 0.02d0*sqrt(5.0d0)
 cc(1)=(center(1)+1.0d0)/2.0d0
 cc(2)=(center(2)+1.0d0)/2.0d0
!  print *,"xy",xy
  call rad_cal(xy,cc,tt)

! if(tt .ge. 0.0d0 .and. tt .le. 0.25d0*pi)then
!  print *,"p1"
!  call dist_hypocycloid(xy,0.0d0,0.25d0*pi,tcrit)
! elseif(tt .le. 0.5d0*pi)then
!  print *,"p2"
!  call dist_hypocycloid(xy,0.25d0*pi,0.5d0*pi,tcrit)
! elseif(tt .le. 0.75d0*pi)then
!  print *,"p3"
!  call dist_hypocycloid(xy,0.5d0*pi,0.75d0*pi,tcrit)
! elseif(tt .le. pi)then
!  print *,"p4"
!  call dist_hypocycloid(xy,0.75d0*pi,pi,tcrit)
! elseif(tt .le. 1.25d0*pi)then
!  print *,"p5"
!  call dist_hypocycloid(xy,pi,1.25d0*pi,tcrit)
! elseif(tt .le. 1.5d0*pi)then
!  print *,"p6"
!  call dist_hypocycloid(xy,1.25*pi,1.5d0*pi,tcrit)
! elseif(tt .le. 1.75d0*pi)then
!  print *,"p7"
!  call dist_hypocycloid(xy,1.5d0*pi,1.75d0*pi,tcrit)
! elseif(tt .le. 2.0d0*pi)then
!  print *,"p8"
!  call dist_hypocycloid(xy,1.75d0*pi,2.0d0*pi,tcrit)
! else
!  print *,"invalid tt", tt
!  stop
! endif   

 if(tt .ge. 0.0d0 .and. tt .le. 0.5d0*pi)then
  call dist_hypocycloid(xy,1,flag,tcrit)
   tcrit1=0.0d0
   tcrit2=0.5d0*pi
 elseif(tt .le. pi)then
  call dist_hypocycloid(xy,2,flag,tcrit)
  tcrit1=0.5d0*pi
  tcrit2=pi
 elseif(tt .le. 1.5d0*pi)then
  call dist_hypocycloid(xy,3,flag,tcrit)
  tcrit1=pi
  tcrit2=1.5d0*pi
 elseif(tt .le. 2.0d0*pi)then
  call dist_hypocycloid(xy,4,flag,tcrit)
  tcrit1=1.5d0*pi
  tcrit2=2.0d0*pi
 else
  print *,"invalid tt", tt
  stop
 endif 

 dist3=sqrt((x-cc(1))**2.0d0+(y-cc(2))**2.0d0)
 if(flag .eq. 0)then
  xtheta=(0.6d0*cos(tcrit)+0.2d0*cos(3.0d0*tcrit) &
           + center(1)+1.0d0)/2.0d0 
  ytheta=(0.6d0*sin(tcrit)-0.2d0*sin(3.0d0*tcrit) &
           + center(2)+1.0d0)/2.0d0
  dist1=sqrt((xtheta-x)**2.0d0+(ytheta-y)**2.0d0)
  dist2=sqrt((xtheta-cc(1))**2.0d0+(ytheta-cc(2))**2.0d0)

 elseif(flag .eq. 1)then
  xtheta1=(0.6d0*cos(tcrit1)+0.2d0*cos(3.0d0*tcrit1) &
           + center(1)+1.0d0)/2.0d0 
  ytheta1=(0.6d0*sin(tcrit1)-0.2d0*sin(3.0d0*tcrit1) &
           + center(2)+1.0d0)/2.0d0
  dist4=sqrt((xtheta1-x)**2.0d0+(ytheta1-y)**2.0d0)
  xtheta2=(0.6d0*cos(tcrit2)+0.2d0*cos(3.0d0*tcrit2) &
           + center(1)+1.0d0)/2.0d0 
  ytheta2=(0.6d0*sin(tcrit2)-0.2d0*sin(3.0d0*tcrit2) &
           + center(2)+1.0d0)/2.0d0
  dist5=sqrt((xtheta2-x)**2.0d0+(ytheta2-y)**2.0d0)  

  if(dist4 .gt. dist5)then
   xtheta=xtheta2
   ytheta=ytheta2
   dist1=dist5

  else
   xtheta=xtheta1
   ytheta=ytheta1
   dist1=dist4   
   
  endif
   dist2=sqrt((xtheta-cc(1))**2.0d0+(ytheta-cc(2))**2.0d0)
 else
  print *,"flag error 236"
  stop
 endif
 
! dist3= (xtheta-c1)*(x-xtheta) +&
!        (ytheta-c2)*(y-ytheta)

 if(imat .eq. 1)then
  if(dist3 .le. dist2)then
   dist =dist1
  else
   dist=-dist1
  endif
 elseif(imat .eq. 2) then
  if(dist3 .le. dist2)then
   dist=-dist1
  else
   dist=dist1
  endif
 else
  print *,"wrong imat flag in, 130"
  stop
 endif


!if(1 .eq. 0)then 
!elseif(probtype_in .eq. 4)then
  ! convert from (-1,1)  to (0,1)
! x0 = 2.0d0*x-1.0d0
! y0 = 2.0d0*y-1.0d0

! c1 = 0.02d0*sqrt(5.0d0)
! c2 = 0.02d0*sqrt(5.0d0)
 
! tt = atan((y0-c2)/(x0-c1));  
! if((x0-c1) .ge. 0.0d0 .and. (y0-c2) .ge. 0.0d0)then
    ! do nothing
! elseif((x0-c1) .le. 0.0d0 .and. (y0-c2) .gt. 0.0d0)then
!    tt = tt + pi;
! elseif((x0-c1) .lt. 0.0d0 .and. (y0-c2) .lt. 0.0d0)then
!    tt = tt +pi;
! else
!    tt = 2.0d0*pi + tt;
! endif

! dist1 = -(sqrt((x0-c1)**2.0d0 + (y0-c2)**2.0d0) - &
!         (0.5d0 + 0.2d0*sin(5.0d0*tt)))
! if(imat .eq. 1)then
!  dist = dist1
! elseif(imat .eq. 2) then
!  dist = -dist1
! else
!  print *,"wrong imat flag in, 130"
!  stop
! endif
!endif



else
 print *,"probtype_in invalid"
 stop
endif


end subroutine dist_fns
!----------------------------------------------------
 subroutine dist_to_boundary(xy,dist)    ! distance from center to the boudary of comutational domain
 implicit none                           ! with the same rad

 real(kind=8),intent(in)   :: xy(2)
 real(kind=8)              :: cc(2)
 real(kind=8)              :: dist
 real(kind=8)              :: dd(2)
 real(kind=8)              :: theta

 cc=0.5d0
 if(cc(1) .eq. xy(1) .and. cc(2) .eq. xy(2))then
  print *,"warning, cc=xy,934"
  stop
 endif

 call rad_cal(xy,cc,theta) 
 dd=0.0d0

 if(theta .eq. 0.0d0 .or. theta .eq. pi*0.5d0 &
     .or. theta .eq. pi .or. theta .eq. pi*1.5d0)then
   dist=0.5d0
 elseif(theta .eq. pi*0.25d0 .or. theta .eq. pi*0.75d0 &
      .or. theta .eq. pi*1.25d0  .or. theta .eq. pi*1.75d0)then
   dist=0.5d0*sqrt(2.0d0)
 elseif(theta .gt. pi*0.25d0 .and. theta .lt. pi*0.75d0)then
   dd(2)=1.0d0
   dd(1)=0.5d0/tan(theta)+0.5d0
   call l2normd(2,cc,dd,dist)
 elseif(theta .gt. pi*0.75d0 .and. theta .lt. pi*1.25d0)then
   dd(1)=0.0d0
   dd(2)=-0.5*tan(theta)+0.5d0
   call l2normd(2,cc,dd,dist)
 elseif(theta .gt. pi*1.25d0 .and. theta .lt. pi*1.75d0)then
   dd(2)=0.0d0
   dd(1)=-0.5d0/tan(theta)+0.5d0
   call l2normd(2,cc,dd,dist)
 elseif(theta .gt. 0.0d0 .and. theta .lt. pi*0.25d0)then
   dd(1)=1.0d0
   dd(2)=0.5*tan(theta)+0.5d0
   call l2normd(2,cc,dd,dist)
 elseif(theta .gt. pi*1.75d0 .and. theta .lt. pi*2.0d0)then
   dd(1)=1.0d0
   dd(2)=0.5*tan(theta)+0.5d0
   call l2normd(2,cc,dd,dist)
 else
  print *,"check theta, 960", theta
  stop
 endif

 



 end subroutine
!------------------------------------------------- 
subroutine cross_product(sdim,x1,x2,cc,crossp)
implicit none

integer, intent(in)   :: sdim
real(kind=8),intent(in)  :: x1(sdim),x2(sdim),cc(sdim)
real(kind=8)           :: crossp(3)
real(kind=8)           :: u(sdim),v(sdim)
integer                :: i

 do i=1,sdim
  u(i)=x1(i)-cc(i)
  v(i)=x2(i)-cc(i)
 enddo

 if(sdim .eq. 2)then
   crossp(1)= 0.0d0
   crossp(2)= 0.0d0
   crossp(3)= u(1)*v(2)-u(2)*v(1)
 else
  print *,"check cross_product"
  stop
 endif

end subroutine cross_product
!-------------------------------------------------
subroutine getline_from_two_points(p1,p2,mx,my,c)
implicit none

real(kind=8),intent(in)    :: p1(2),p2(2)
real(kind=8),intent(out)   :: mx,my,c
real(kind=8)                :: m

if(abs(p1(1)-p2(1)) .lt. 1.0e-10 .and. abs(p1(2)-p2(2)) .lt. 1.0e-10 ) then
  print *,"p1 and p2 are coincide, 801"
  stop
elseif(abs(p1(1)-p2(1)) .lt. 1.0e-10)then
  mx=1.0d0
  my=0.0d0
  c=-p1(1)
elseif(abs(p1(2)-p2(2)) .lt. 1.0e-10)then
  mx=0.0d0
  my=1.0d0
  c=-p1(2)
else
  mx=(p1(2)-p2(2))/(p1(1)-p2(1))
  my=-1.0d0
  c=-mx*p1(1)+p1(2)
endif

end subroutine getline_from_two_points


! --------------------------------------------------
subroutine find_cloest_2d(sflag,dist,xin,xout)
implicit none

integer     ,intent(in)    :: sflag
real(kind=8),intent(in)    :: dist
real(kind=8),intent(in)    :: xin(2)
real(kind=8)               :: xout(2)
real(kind=8)               :: ss
real(kind=8)               :: xgrad(2)
real(kind=8)               :: xnorm
integer                    :: i

if(sflag .eq. 1)then
 ss = +1.0d0
elseif(sflag .eq. -1)then
 ss = -1.0d0
else
 print *,"invalid sign flag, 717"
 stop
endif

if(dist .eq. 0.0d0)then
 print *,"dist is 0, 723"
 stop
endif

xnorm = sqrt((xin(1)-0.5d0)**2.0d0 + (xin(2)-0.5d0)**2.0d0)
do i=1,2
 xgrad(i)=(xin(i)-0.5d0)/xnorm
 xout(i)= xin(i)+ss*xgrad(i)*abs(dist) 
enddo


end subroutine find_cloest_2d
!--------------------------------------------


subroutine dist_point_to_lined(sdim,p1,p2,x,dist)
implicit none
! represent the line in parametric form,(v = f(s))
! v^x = x1 + (x2-x1)s
! v^y = y1 + (y2-y1)s
! v^z = z1 + (z2-z1)s
! 
! if the closest point on the line to point x  is outside p1 -- p2, 
! --------> return the distance from x either to p1 or p2 which is shorter.
! otherwise
! --------> return the distance from x to the cloest point

integer,intent(in)           :: sdim
real(kind=8),intent(in)      ::  p1(sdim),p2(sdim),x(sdim)
real(kind=8)                 ::  dist


real(kind=8)                 :: diff10,diff21,diffx
real(kind=8),allocatable     :: x10(:), x21(:)
integer                      :: i
real(kind=8)                 :: s


dist = 0.0d0

allocate(x10(sdim),x21(sdim))
do i = 1,sdim
 x10(i) = p1(i) - x(i)
 x21(i) = p2(i) - p1(i)
enddo

if (maxval(abs(x21)) .lt. 10d-8)then
 print *,"p1 and p2 are coincide with each other",p1,p2
 stop
endif

call l2normd(sdim, p1, x, diff10)
call l2normd(sdim, p2, p1, diff21)

!print *,"diff10",diff10
!print *,"diff21",diff21

s = -1.0d0*(dot_product(x10,x21))/(diff21**2.0d0)

!write(13,*) "s=",s

if(s .gt. 1.0d0)then
 call l2normd(sdim, p2, x,dist)
elseif(s .lt. 0.0d0)then
 call l2normd(sdim,p1,x,dist)
else
! if(abs((diff10**2.0d0 )*(diff21**2.0d0) - & 
!        (dot_product(x10,x21))**2.0d0) .lt. 1.0e-10)then
!   dist= 0.0d0
! else
  dist = sqrt(((diff10**2.0d0 )*(diff21**2.0d0) - & 
        (dot_product(x10,x21))**2.0d0)/ &
        (diff21**2.0d0))
! endif
endif

deallocate(x10,x21) 

end subroutine dist_point_to_lined
! ---------------------------------------------------------------
subroutine dist_point_to_line_modify(sdim,p1,p2,x,dist,xc)
implicit none
! represent the line in parametric form,(v = f(s))
! v^x = x1 + (x2-x1)s
! v^y = y1 + (y2-y1)s
! v^z = z1 + (z2-z1)s
! 
! if the closest point on the line to point x  is outside p1 -- p2, 
! --------> return the distance from x either to p1 or p2 which is shorter.
! otherwise
! --------> return the distance from x to the cloest point

integer,intent(in)           :: sdim
real(kind=8),intent(in)      :: p1(sdim),p2(sdim),x(sdim)
real(kind=8)                 :: dist
real(kind=8)                 :: xc(sdim)


real(kind=8)                 :: diff10,diff21,diffx
real(kind=8),allocatable     :: x10(:), x21(:)
integer                      :: i
real(kind=8)                 :: s


dist = 0.0d0

allocate(x10(sdim),x21(sdim))
do i = 1,sdim
 x10(i) = p1(i) - x(i)
 x21(i) = p2(i) - p1(i)
enddo

if (maxval(abs(x21)) .lt. 10d-8)then
 print *,"p1 and p2 are coincide with each other",p1,p2
 stop
endif

call l2normd(sdim, p1, x, diff10)
call l2normd(sdim, p2, p1, diff21)

!print *,"diff10",diff10
!print *,"diff21",diff21

s = -1.0d0*(dot_product(x10,x21))/(diff21**2.0d0)

!print *,"s=",s

if(s .gt. 1.0d0)then
 call l2normd(sdim, p2, x,dist)
 xc=p2
elseif(s .lt. 0.0d0)then
 call l2normd(sdim,p1,x,dist)
 xc=p1
else
 if(abs((diff10**2.0d0 )*(diff21**2.0d0) - & 
        (dot_product(x10,x21))**2.0d0) .lt. 1.0e-10)then
   dist= 0.0d0
   xc=x
 else
  dist = sqrt(((diff10**2.0d0 )*(diff21**2.0d0) - & 
        (dot_product(x10,x21))**2.0d0)/ &
        (diff21**2.0d0))
  do i=1,sdim
   xc(i)=p1(i)+ s*(x(i)-p1(i))
  enddo
 endif
endif

deallocate(x10,x21) 

end subroutine dist_point_to_line_modify





subroutine starshape(xt)
implicit none

real(kind=8)         :: theta(pcurve_num)
real(kind=8),intent(out):: xt(2,pcurve_num)
integer          :: num
integer :: i

num = pcurve_num

do i = 1,num
 theta(i) = (i-1)*2.0d0*pi/num
enddo

do i = 1,num
 xt(1,i) = 0.02d0*sqrt(5.0d0) + &
        (0.5d0 + 0.2d0*sin(5.0d0*theta(i)))*cos(theta(i))
 xt(2,i) = 0.02d0*sqrt(5.0d0) + &
        (0.5d0 + 0.2d0*sin(5.0d0*theta(i)))*sin(theta(i))
enddo


end subroutine starshape


subroutine starshape2(xt)
implicit none

real(kind=8)         :: theta(pcurve_num)
real(kind=8),intent(out):: xt(2,pcurve_num)
integer          :: num
integer :: i

num = pcurve_num

do i = 1,num
 theta(i) = (i-1)*2.0d0*pi/num
enddo

do i = 1,num
 xt(1,i) = 0.02d0*sqrt(5.0d0) + &
        (0.5d0 -pentaeps + 0.2d0*sin(5.0d0*theta(i)))*cos(theta(i))
 xt(2,i) = 0.02d0*sqrt(5.0d0) + &
        (0.5d0 -pentaeps + 0.2d0*sin(5.0d0*theta(i)))*sin(theta(i))
enddo


end subroutine starshape2


subroutine asteroidshape(xt)
implicit none

real(kind=8)         :: theta(pcurve_num+1)
real(kind=8),intent(out):: xt(2,pcurve_num+1)
integer          :: num
integer :: i
integer :: cenflag
real(kind=8) :: cc(2)


num = pcurve_num
 cenflag= 0

do i = 1,num+1
 theta(i) = (i-1)*2.0d0*pi/real(num,8)
enddo

if(cenflag .eq. 1)then
 cc(1) = 0.02d0*sqrt(5.0d0)
 cc(2) = 0.02d0*sqrt(5.0d0)
elseif(cenflag .eq. 0)then
 cc = 0.0d0
else
 print *,"wrong cenflag"
 stop
endif

if(shapeflag .eq. 0)then

 do i = 1,num+1

  xt(1,i) = cc(1) + &
           0.6d0*cos(theta(i)) + 0.2d0*cos(3.0d0*theta(i))
  xt(2,i) = cc(2) + &
           0.6d0*sin(theta(i)) - 0.2d0*sin(3.0d0*theta(i))
  call rad_cal(xt(:,i),cc,pcurve_rad(i))
 enddo
elseif(shapeflag .eq. 1)then
  do i = 1,pcurve_num/4+1
   xt(1,i)= 0.5d0 -(i-1)*0.5/(real(pcurve_num,8)/4.0)
   xt(2,i)= 0.0d0 +(i-1)*0.5/(real(pcurve_num,8)/4.0)

   xt(1,pcurve_num/4+i)= 0.0d0 -(i-1)*0.5/(real(pcurve_num,8)/4.0)
   xt(2,pcurve_num/4+i)= 0.5d0 -(i-1)*0.5/(real(pcurve_num,8)/4.0)

   xt(1,pcurve_num/2+i)= -0.5d0 +(i-1)*0.5/(real(pcurve_num,8)/4.0)
   xt(2,pcurve_num/2+i)= 0.0d0 -(i-1)*0.5/(real(pcurve_num,8)/4.0)
 
   xt(1,pcurve_num/4*3+i)= 0.0d0 +(i-1)*0.5/(real(pcurve_num,8)/4.0)
   xt(2,pcurve_num/4*3+i)= -0.5d0 +(i-1)*0.5/(real(pcurve_num,8)/4.0)
  enddo 
elseif(shapeflag .eq. 2)then
 
 do i=1,num+1
  xt(1,i)= cc(1) + 0.5d0*cos(theta(i))
  xt(2,i)= cc(2) + 0.5d0*sin(theta(i))
 enddo

else
 print *,"wrong shapeflag in asteroidshape,xt"
 stop
endif

end subroutine asteroidshape




subroutine dist_point_to_arc(x,arc1,arc2,arcc,dist)
! for arc of a circle
!arc1 to arc2--arc start to arcend--counterclockwise
!arcc--arc center
implicit none

real(kind=8),intent(in)   :: x(2)
real(kind=8),intent(in)   :: arc1(2),arc2(2),arcc(2)
real(kind=8)              :: dist
real(kind=8)              :: tt(3)
real(kind=8)              :: r1,r2
integer                   :: i

dist=0.0d0
call rad_cal(arc1,arcc,tt(1))
call rad_cal(arc2,arcc,tt(2))
call rad_cal(x,arcc,tt(3))

do i=1,3
 if(tt(i) .gt. 2.0d0*pi .or. tt(i) .lt. 0.0d0)then
  print *,"rad invalid 489"
  stop
 endif
enddo

if(tt(3) .lt. tt(1))then
 call l2normd(2,x,arc1,dist)
elseif(tt(3) .gt. tt(2))then
 call l2normd(2,x,arc2,dist)
else
 call l2normd(2,x,arcc,r1)
 call l2normd(2,arcc,arc1,r2)
 dist=abs(r1-r2)
endif

end subroutine dist_point_to_arc
!-------------------------------------------------------
subroutine rad_cal(x,c,theta)
! theta = [0,2pi)
implicit none

real(kind=8),intent(in)    :: x(2),c(2)
real(kind=8)               :: theta             

 theta = atan((x(2)-c(2))/(x(1)-c(1))); 
 if(x(2) .eq. c(2) .and. x(1)-c(1) .gt. 0.0d0)then
   theta=0.0d0 
 elseif((x(1)-c(1)) .gt. 0.0d0 .and. (x(2)-c(2)) .gt. 0.0d0)then
    ! do nothing
 elseif(x(1) .eq. c(1) .and. x(2) .gt. c(2))then
   theta=0.5d0*pi
 elseif((x(1)-c(1)) .lt. 0.0d0 .and. (x(2)-c(2)) .gt. 0.0d0)then
    theta = theta + pi;
 elseif(x(2) .eq. c(2) .and. (x(1)-c(1)) .lt. 0.0d0)then
   theta=pi
 elseif((x(1)-c(1)) .lt. 0.0d0 .and. (x(2)-c(2)) .lt. 0.0d0)then
    theta = theta +pi;
 elseif(x(1) .eq. c(1) .and. x(2) .lt. c(2))then
   theta=1.5d0*pi 
 elseif(x(1) .eq. c(1) .and. x(2) .eq. c(2))then
   theta=0.0d0
 else
    theta = 2.0d0*pi + theta;
 endif

end subroutine

!----------------------------------------------------------
subroutine dist_point_to_line(p1,p2,x,dist)
implicit none
! represent the line in parametric form,(v = f(s))
! v^x = x1 + (x2-x1)s
! v^y = y1 + (y2-y1)s
! v^z = z1 + (z2-z1)s
! 
! if the closest point on the line to point x  is outside p1 -- p2, 
! --------> return the distance from x either to p1 or p2 which is shorter.
! otherwise
! --------> return the distance from x to the cloest point
integer,parameter            :: SDIM = 2
real(kind=8),intent(in)      ::  p1(SDIM),p2(SDIM),x(SDIM)
real(kind=8)                 ::  dist,dist1,dist2


real(kind=8)                 :: diff10,diff21
real(kind=8)                 :: x10(SDIM), x21(SDIM)
integer                      :: i
real(kind=8)                 :: s



dist = 0.0d0

do i = 1,SDIM
 x10(i) = p1(i) - x(i)
 x21(i) = p2(i) - p1(i)
 
 if(abs(x10(i)) .lt. 10e-10)then
  x10(i) = 0.0d0
 endif
 if(abs(x21(i)) .lt. 10e-10)then
  x21(i) = 0.0d0
 endif
 
enddo



if (maxval(abs(x21)) .lt. 10d-8)then
 print *,"p1 and p2 are coincide with each other",p1,p2
 stop
endif

!if (maxval(abs(x10)) .lt. 10d-8)then
! print *,"p1 and x are coincide with each other",p1,p2
! stop
!endif

call l2normd(2,p1, x, diff10)
call l2normd(2,p2, p1, diff21)


s = -1.0d0*(dot_product(x10,x21))/(diff21**2.0d0)

if(s .le. 1.0d0 .and. s .ge. 0.0)then

 do i = 1,SDIM
  dist = dist + (x10(i) + x21(i)*s)**2.0d0
 enddo
  dist = sqrt(dist)
else
 call l2normd(2,p1,x,dist1)
 call l2normd(2,p2,x,dist2)
 dist = min(dist1,dist2)
endif

end subroutine dist_point_to_line
!---------------------------------------------------------
subroutine l2normd(sdim,x1,x2, x1x2norm)
implicit none

integer,intent(in)       :: sdim
real(kind=8),intent(in)  :: x1(sDim),x2(SDIM)
real(kind=8)             :: x1x2norm

integer                  :: i
real(kind=8),allocatable :: diff(:)

x1x2norm = 0.0d0
allocate(diff(SDIM))
do i = 1,SDIM
 diff(i) = x1(i)-x2(i)
enddo

do i = 1,SDIM
 x1x2norm = x1x2norm + diff(i)**2.0d0
enddo

x1x2norm = sqrt(x1x2norm)

deallocate(diff)

end subroutine l2normd


subroutine dist_hypocycloid_back(xin,tt1,tt2,tt)
implicit none

real(kind=8),intent(in) :: xin(2)
real(kind=8)  :: tt1,tt2

real(kind=8)            :: dist
real(kind=8)            :: c1,c2
real(kind=8)            :: x(2)
integer                 :: i
real(kind=8)            :: res
real(kind=8)            :: val1,val2,val3
real(kind=8)            :: tt
real(kind=8)            :: xtheta,ytheta


 do i=1,2
  x(i)=xin(i) 
 enddo
 print *,"xin",x

  call fund_hypocycloid(x,tt1,val1)
  call fund_hypocycloid(x,tt2,val2)
 if(abs(val1) .lt. 1.0e-8 .and. abs(val1).lt. abs(val2))then
  print *,"1",val1
  tt=tt1
 elseif(abs(val2) .lt. 1.0e-8 .and. abs(val2) .lt. abs(val1))then
  print *,"2",val2
  tt=tt2
 else 
 print *,"3"
 res=1.0e+5
 do while(res .gt. 1.0e-8)
  call fund_hypocycloid(x,tt1,val1)
  call fund_hypocycloid(x,tt2,val2)

  print *,"val1",tt1,val1
  print *,"val2",tt2,val2
  if(val1*val2 .ge. 0.0d0)then
   tt=0.5d0*(tt1+tt2)
   call fund_hypocycloid(x,tt,val3)  
   print *,"theta",tt1,tt2,tt
   print *,"val1 and val2 same sign 667",val1,val2,val3
   if(val3*val1 .lt. 0.0d0)then
    tt2=tt
   elseif(val3*val2 .lt. 0.0d0)then
    tt1=tt
   else
    print *,"check val1 and val2", val1,val2,val3
    stop
   endif
   res=val3
  else
   print *,"in bisection"
   tt=0.5d0*(tt1+tt2)
   call fund_hypocycloid(x,tt,val3)
!   print *,"val3",tt,val3
   if(val3*val1 .lt. 0.0d0)then
    tt2=tt
   elseif(val3*val2 .lt. 0.0d0)then
    tt1=tt
   else
    print *,"check val1 and val2", val1,val2,val3
    stop
   endif
   res=val3
  endif
 enddo
 endif


end subroutine dist_hypocycloid_back


subroutine dist_hypocycloid(xin,quadf,flag,tt)
implicit none

real(kind=8),intent(in) :: xin(2)
integer,intent(in)      :: quadf

real(kind=8)            :: dist
real(kind=8)            :: c1,c2
real(kind=8)            :: x(2)
integer                 :: i
real(kind=8)            :: res
real(kind=8)            :: val1,val2,val3
real(kind=8)            :: tt
real(kind=8)            :: tt1,tt2
integer                 :: flag,step

real(kind=8),parameter  :: zdelta=1.0e-8

 flag = 0
 do i=1,2
  x(i)=xin(i) 
 enddo
! print *,"xin",x
 step=0
 if(quadf .eq. 1) then
  tt1= zdelta
  tt2= pi/2.0d0-zdelta
 elseif(quadf .eq. 2)then
  tt1=pi/2.0d0+zdelta 
  tt2=pi-zdelta
 elseif(quadf .eq. 3)then
  tt1=pi+zdelta 
  tt2=1.5d0*pi-zdelta
 elseif(quadf .eq. 4)then
  tt1=1.5d0*pi+zdelta 
  tt2=2.0d0*pi-zdelta  
 else
  print *,"wrong quadrant flag"
  stop
 endif


 call fund_hypocycloid(x,tt1,val1)
 call fund_hypocycloid(x,tt2,val2)  
 if(val1 .eq. 0.0d0)then
  tt=tt1
 elseif(val2 .eq. 0.0d0)then
  tt=tt2
 elseif(val1*val2 .gt. 0.0d0)then
  flag=1
 elseif(val1*val2 .lt. 0.0d0)then
  res=1.0e+5
  step=0
  do while(res .gt. 1.0e-10 .and. step .lt. 1000)
  print *,"iter step",step 
  
  call fund_hypocycloid(x,tt1,val1)
  call fund_hypocycloid(x,tt2,val2)
   tt=0.5d0*(tt1+tt2)
   call fund_hypocycloid(x,tt,val3)
!   print *,"val3",tt,val3
   if(val3 .eq. 0.0d0)then
    !  do nothing
   elseif(val3*val1 .lt. 0.0d0)then
    tt2=tt
   elseif(val3*val2 .lt. 0.0d0)then
    tt1=tt
   else
    print *,"check val1 and val2", val1,val2,val3
    stop
   endif
   res=val3
   step=step+1 
  enddo

 else
  print *,"out of cases"
  stop
 endif  


! elseif(quadf .eq. 2) then
!  tt1=0.75d0*pi
!  tt2=0.0d0
!  res = tt1
!  do while(res .gt. 1.0e-10 .and. step .lt. 1000)
!   step = step+1
!   if(step .gt. 100)then
!   print *,"quad",quadf,"step",step,"res",res
!   endif
!   call fund_hypocycloid(x,tt1,val1)
!   call fundd_hypocycloid(x,tt1,val2)
!   tt2=tt1- 0.5d0*val1/val2 
!   if(tt2 .gt. pi .or. tt2 .le. 0.5d0*pi)then
!    flag = 1
!    exit
!   endif
!   res=abs(tt1-tt2)
!   tt1=tt2
!  enddo
! elseif(quadf .eq. 3) then
!  tt1=1.25d0*pi
!  tt2=0.0d0
!  res = tt1
!  do while(res .gt. 1.0e-10 .and. step .lt. 1000)
!   step=step+1
!   if(step .gt. 100)then
 !  print *,"quad",quadf,"step",step,"res",res
!   print *,"x",x
!   endif
!   call fund_hypocycloid(x,tt1,val1)
!   call fundd_hypocycloid(x,tt1,val2)
!   print *,"val1",val1,"val2",val2
!   tt2=tt1- 0.5d0*val1/val2 
!   print *,"tt1",tt1,"tt2",tt2

!   if(tt2 .gt. 1.5d0*pi .or. tt2 .le. pi)then
!    flag = 1
!    exit
!   endif
!   res=abs(tt1-tt2)
!   tt1=tt2
!  enddo

! elseif(quadf .eq. 4) then
!  tt1=1.75d0*pi
!  tt2=0.0d0
!  res = tt1
!  do while(res .gt. 1.0e-10 .and. step .lt. 1000)
!   step=step+1
!   if(step .gt. 100)then
!   print *,"quad",quadf,"step",step,"res",res
!   endif
!   call fund_hypocycloid(x,tt1,val1)
!   call fundd_hypocycloid(x,tt1,val2)
!   tt2=tt1- 0.5d0*val1/val2 
!   if(tt2 .gt. 2.0d0*pi .or. tt2 .le. 1.5d0*pi)then
!    flag = 1
!    exit
!   endif
!   res=abs(tt1-tt2)
!   tt1=tt2
!  enddo



end subroutine dist_hypocycloid


subroutine fund_hypocycloid(x,theta,val)
! x in domain [0,1]
implicit none

real(kind=8),intent(in)  :: x(2)
real(kind=8),intent(in)  :: theta
real(kind=8)             :: val

real(kind=8)             :: c1,c2
real(kind=8)             :: xtheta,ytheta

 c1 = 0.02d0*sqrt(5.0d0)
 c2 = 0.02d0*sqrt(5.0d0)

 xtheta=(0.6d0*cos(theta)+0.2d0*cos(3.0d0*theta) &
           + c1+1.0d0)/2.0d0 
 ytheta=(0.6d0*sin(theta)-0.2d0*sin(3.0d0*theta) &
           + c2+1.0d0)/2.0d0

 val=0.6d0*(xtheta-x(1))*(-sin(theta)-sin(3.0d0*theta))+ &
     0.6d0*(ytheta-x(2))*(cos(theta)-cos(3.0d0*theta))


end subroutine fund_hypocycloid


subroutine fundd_hypocycloid(x,theta,val)
implicit none

real(kind=8),intent(in)  :: x(2)
real(kind=8),intent(in)  :: theta
real(kind=8)             :: val

real(kind=8)             :: c1,c2
real(kind=8)             :: xtheta,ytheta

 c1 = 0.02d0*sqrt(5.0d0)
 c2 = 0.02d0*sqrt(5.0d0)

 xtheta=(0.6d0*cos(theta)+0.2d0*cos(3.0d0*theta) &
           + c1+1.0d0)/2.0d0 
 ytheta=(0.6d0*sin(theta)-0.2d0*sin(3.0d0*theta) &
           + c2+1.0d0)/2.0d0

 val=0.5d0*(-0.6d0*sin(theta)-0.6d0*sin(3.0d0*theta))**2.0d0 + &
     (xtheta-x(1))*(-0.6d0*cos(theta)-1.8d0*cos(3.0d0*theta)) + &
     0.5d0*(0.6d0*cos(theta)-0.6d0*cos(3.0d0*theta))**2.0d0 + &
     (ytheta-x(2))*(-0.6d0*sin(theta)+1.8d0*sin(3.0d0*theta)) 


end subroutine





!-----------------------------------------------------------------
subroutine update_volncen(nmat_in,dist,im,center,h,vol,cxtemp,cytemp)
implicit none

integer,intent(in)            :: nmat_in,im
real(kind=8),intent(in)            :: dist
real(kind=8),intent(in)       :: center(2)
real(kind=8),intent(in)       :: h
real(kind=8)                  :: vol(nmat_in)
real(kind=8)                  :: cxtemp(nmat_in),cytemp(nmat_in)

if(dist .gt. 0.0d0) then
  vol(im) = vol(im) + h*h
  cxtemp(im) = cxtemp(im) + h*h*center(1) 
  cytemp(im) = cytemp(im) + h*h*center(2)
endif

end subroutine update_volncen
!---------------------------------------------------------------------
subroutine renormalize_vf(nmat_in,vf)
implicit none

integer              :: nmat_in
real(kind=8)         :: vf(nmat_in)
real(kind=8)         :: vfsum
integer              :: i

vfsum = 0.0d0

do i = 1,nmat_in
  vfsum = vfsum + vf(i)
enddo
if (vfsum.le.0.0) then
 print *,"vfsum invalid"
 stop
endif
do  i =1,nmat_in
  vf(i) = vf(i)/vfsum
enddo


end subroutine
!-------------------------------------------------------
subroutine form_cell(vertex,cell)
implicit none

real(kind=8)             :: vertex(4,2)
type(polygon)            :: cell
integer                  :: i

if(cell%nodesnum .le. 2) then
  print *,"nodes num invalid in cell"
else
  
  cell%nodes(1)%val = vertex(1,:)
  cell%nodes(2)%val = vertex(2,:)
  cell%nodes(3)%val = vertex(4,:)
  cell%nodes(4)%val = vertex(3,:)  

  cell%SIDES(1)%PT(1)= cell%NODES(1)
  cell%SIDES(1)%PT(2)= cell%NODES(2)
  cell%SIDES(2)%PT(1)= cell%NODES(2)
  cell%SIDES(2)%PT(2)= cell%NODES(3)
  cell%SIDES(3)%PT(1)= cell%NODES(3)
  cell%SIDES(3)%PT(2)= cell%NODES(4)
  cell%SIDES(4)%PT(1)= cell%NODES(4)
  cell%SIDES(4)%PT(2)= cell%NODES(1)

endif


end subroutine form_cell
!-------------------------------------------------------
subroutine normalize_cutline(a,b,intercept)
implicit none

real(kind=8)              :: a, b, intercept
real(kind=8)              :: c

 c = sqrt(a*a + b*b)
 
 a = a/c
 b = b/c
 intercept = intercept/c

end subroutine normalize_cutline


!-------------------------------------------------------
subroutine AdaptQuad_sort(im,center,h,centers,dist,flag,probtype_in)
implicit none

integer                          :: flag
integer,intent(in)               :: im,probtype_in
real(kind=8),intent(in)          :: h
real(kind=8)                     :: dist
real(kind=8)                     :: center(2), centers(4,2)

  if (h.le.0.0) then
   print *,"h invalid"
   stop
  endif

  centers = 0.0d0
   flag = 0
   call dist_fns(im,center(1),center(2),dist,probtype_in)
   if(abs(dist) .gt. sqrt(2.0d0)*h) then
      ! do nothing
      flag = 0
   else
     flag = 1  
     call AdaptQuad_sub(center, h, centers)  
   endif     


end subroutine AdaptQuad_sort


!--------------------------------------------------------
subroutine AdaptQuad_sub(center, h, centers)
implicit none

real(kind=8)             :: center(2)
real(kind=8)             :: centers(4,2)
real(kind=8),intent(in)  :: h

 if (h.le.0.0) then
  print *,"h invalid"
  stop
 endif
!    
! 3  4 
! 1  2
! 
 centers(1,1)  = center(1) - h*0.25d0
 centers(1,2)  = center(2) - h*0.25d0
 centers(2,1)  = center(1) + h*0.25d0
 centers(2,2)  = center(2) - h*0.25d0
 centers(3,1)  = center(1) - h*0.25d0
 centers(3,2)  = center(2) + h*0.25d0
 centers(4,1)  = center(1) + h*0.25d0
 centers(4,2)  = center(2) + h*0.25d0


end subroutine AdaptQuad_sub
!---------------------------------------------------------
subroutine find_vertex(center, h, vertex)
implicit none

real(kind=8)           :: center(2)
real(kind=8)           :: h
real(kind=8)           :: vertex(4,2)

 if (h.le.0.0) then
  print *,"h invalid"
  stop
 endif

!    
! 3  4 
! 1  2
! 
 vertex(1,1)  = center(1) - h*0.5d0
 vertex(1,2)  = center(2) - h*0.5d0
 vertex(2,1)  = center(1) + h*0.5d0
 vertex(2,2)  = center(2) - h*0.5d0
 vertex(3,1)  = center(1) - h*0.5d0
 vertex(3,2)  = center(2) + h*0.5d0
 vertex(4,1)  = center(1) + h*0.5d0
 vertex(4,2)  = center(2) + h*0.5d0


end subroutine find_vertex

! -----------------------------------------------------

subroutine slopecal(im,center,mx,my,probtype_in)
implicit none

integer, intent(in)                  :: probtype_in
integer                              :: im
real(kind=8),intent(in)              :: center(2)
real(kind=8)                         :: mx,my

real(kind=8)                     :: x,y,dist
real(kind=8)                     :: dist1,dist2,dist3,dist4
integer                          :: i
real(kind=8)                     :: r1,r2

 
 x = center(1)
 y = center(2)



! slope is grad phi/|grad phi|

if (probtype_in.eq.0) then
 mx=0.0
 if (im.eq.1) then
  my=-1.0
 else if (im.eq.2) then
  my=1.0
 else
  print *,"im invalid 1"
  stop
 endif

else if (probtype_in.eq.2) then
 my=0.0
 if (im.eq.1) then
  mx=-1.0
 else if (im.eq.2) then
  mx=1.0
 else
  print *,"im invalid 2"
  stop
 endif

else if (probtype_in.eq.1) then

 r1=radcen-radeps
 r2=radcen+radeps

 dist1= sqrt((x-0.5d0)**2.0d0+(y-0.5d0)**2.0d0)
 
 if (im.eq. 1) then ! inner material
   dist= r1 - dist1
   mx = (x - 0.5d0)/(dist - r1)
   my = (y - 0.5d0)/(dist - r1)
 elseif(im .eq. 2) then
   if (dist1.le.radcen) then
    mx = (x - 0.5d0)/dist1
    my = (y - 0.5d0)/dist1
   else if (dist1.ge.radcen) then
    mx = -(x - 0.5d0)/dist1
    my = -(y - 0.5d0)/dist1
   else
    print *,"dist1 is bad"
    stop
   endif 
 elseif(im .eq. 3)then
   dist = dist1 - r2  ! outer material
   mx = (x - 0.5d0)/(dist + r2)
   my = (y - 0.5d0)/(dist + r2)
 else
   print *,"wrong number of mat"
 endif


else 
 print *,"probtype_in invalid"
 stop
endif

end subroutine slopecal

!--------------------------------------------------------
SUBROUTINE triangle_interface_detect(im,probtype_in,v1,v2,v3,area,cen)
implicit none

integer,intent(in)       :: probtype_in,im
real(kind=8),intent(in)  :: v1(2),v2(2),v3(2)
real(kind=8)             :: d1,d2,d3,dcheck
real(kind=8)             :: d(3)
real(kind=8)             :: v(3,2)

real(kind=8)            :: ratio,ratio1,ratio2
real(kind=8)            :: x1(2),x2(2),xx(2,2)

real(kind=8)            :: area, cen(2)
real(kind=8)            :: area1,area2
real(kind=8)            :: cen1(2),cen2(2),cencheck(2)
integer                 :: i,j
integer                 :: ct


 area = 0.0d0
 cen = 0.0d0
 x1 = 0.0d0
 x2 = 0.0d0
 xx = 0.0d0
 d = 0.0d0
 v = 0.0d0
 area1 = 0.0d0
 area2 = 0.0d0
 cen1 = 0.0d0
 cen2 = 0.0d0



 do i = 1,2
  v(1,i) = v1(i)
 enddo
 do i = 1,2
  v(2,i) = v2(i)
 enddo
 do i = 1,2
  v(3,i) = v3(i)
 enddo

!if(probtype_in .eq. 3)then
 call dist_fns(im,v1(1),v1(2),d1,probtype_in)      ! prob_type = 3
 call dist_fns(im,v2(1),v2(2),d2,probtype_in)
 call dist_fns(im,v3(1),v3(2),d3,probtype_in)

! print *,"v1",v1,v(1,:)
! print *,"v2",v2,v(2,:)
! print *,"v3",v3,v(3,:)
! print *,"d1 d2 d3",d1,d2,d3

 d(1) = d1
 d(2) = d2
 d(3) = d3

 if(abs(d1) .lt. 1.0d-10 .and. abs(d2) .lt. 1.0d-10 &                       ! 0  0  0
    .and. abs(d3) .lt. 1.0d-10)then
   print *,"invalid d1 d2 d3, check triangle_interface_detect"
!   print *, v1,d1
!   print *, v2,d2
!   print *, v3,d3
   call tri_centroid(v1,v2,v3,cencheck)
   call dist_fns(im,cencheck(1),cencheck(2),dcheck,probtype_in) 
   if(dcheck .eq. 0.0d0)then
    print *,"dcheck is still 0"
   elseif(dcheck .gt. 0.0d0)then
    call tri_area8(v1,v2,v3,area)
    call tri_centroid(v1,v2,v3,cen) 
   elseif(dcheck .lt. 0.0d0)then
    ! do nothing
   else
    print *,"dcheck in valid",dcheck
    stop
   endif
   
 elseif(d1 .gt. 0.0d0  .and. d2 .gt. 0.0d0 .and. d3 .gt. 0.0d0)then         ! +0  +0  +0
  call tri_area8(v1,v2,v3,area)
  call tri_centroid(v1,v2,v3,cen)
!  print *,"case 1", "area",area
 elseif(d1 .lt. 0.0d0  .and. d2 .lt. 0.0d0 .and. d3 .lt. 0.0d0)then         ! -0  -0  -0   
!  print *,"case 2"
    area = 0.0d0
    cen = 0.0d0
    ! do nothing

  elseif(d1 .eq. 0.0d0 .and. d2 .eq. 0.0d0 .and. d3 .ne. 0.0d0)then
   if(d3 .gt. 0.0d0)then
    call tri_area8(v1,v2,v3,area)
    call tri_centroid(v1,v2,v3,cen)
   else
    area = 0.0d0
    cen=0.0d0
   endif
  elseif(d1 .eq. 0.0d0 .and. d3 .eq. 0.0d0 .and. d2 .ne. 0.0d0)then
   if(d2 .gt. 0.0d0)then
    call tri_area8(v1,v2,v3,area)
    call tri_centroid(v1,v2,v3,cen)
   else
    area = 0.0d0
    cen=0.0d0
   endif
  elseif(d2 .eq. 0.0d0 .and. d3 .eq. 0.0d0 .and. d1 .ne. 0.0d0)then
   if(d1 .gt. 0.0d0)then
    call tri_area8(v1,v2,v3,area)
    call tri_centroid(v1,v2,v3,cen)
   else
    area = 0.0d0
    cen=0.0d0
   endif
  elseif(d1 .eq. 0.0d0 .and. d2*d3 .gt. 0.0d0)then
   if(d2 .gt. 0.0d0)then
    call tri_area8(v1,v2,v3,area)
    call tri_centroid(v1,v2,v3,cen)  
   else
    ! do nothing
   endif
  elseif(d2 .eq. 0.0d0 .and. d1*d3 .gt. 0.0d0)then
   if(d1 .gt. 0.0d0)then
    call tri_area8(v1,v2,v3,area)
    call tri_centroid(v1,v2,v3,cen)  
   else
    ! do nothing
   endif
  elseif(d3 .eq. 0.0d0 .and. d1*d2 .gt. 0.0d0)then
   if(d2 .gt. 0.0d0)then
    call tri_area8(v1,v2,v3,area)
    call tri_centroid(v1,v2,v3,cen)  
   else
    ! do nothing
   endif

 elseif(d1 .eq. 0.0d0 .and. d2*d3 .lt. 0.0d0)then
 !  print *,"case 3"
   ratio = abs(d2)/(abs(d3)+abs(d2))
   do i = 1,2
    x1(i)= v2(i) + ratio*(v3(i)-v2(i))
   enddo
   if(d2 .gt. 0.0d0) then
    call tri_area8(x1,v1,v2,area)
    call tri_centroid(x1,v1,v2,cen)
   elseif(d3 .gt. 0.0d0)then
    call tri_area8(x1,v1,v3,area)
    call tri_centroid(x1,v1,v3,cen)    
   else
    print *,"error, volate d2*d3<0,458"
    stop
   endif
 elseif(d2 .eq. 0.0d0 .and. d1*d3 .lt. 0.0d0)then
!   print *, "case4",d1,d2,d3
   ratio = abs(d1)/(abs(d1)+abs(d3))
   do i = 1,2
    x1(i)= v1(i) + ratio*(v3(i)-v1(i))
   enddo
   if(d1 .gt. 0.0d0) then
    call tri_area8(x1,v1,v2,area)
    call tri_centroid(x1,v1,v2,cen)
   elseif(d3 .gt. 0.0d0)then
    call tri_area8(x1,v2,v3,area)
    call tri_centroid(x1,v2,v3,cen)    
   else
    print *,"error, volate d2*d3<0,473"
    stop
   endif    
 elseif( d3 .eq. 0.0d0 .and. d1*d2 .lt. 0.0d0)then
 !  print *,"case5",d1,d2,d3
   ratio = abs(d1)/(abs(d1)+abs(d2))
   do i = 1,2
    x1(i)= v1(i) + ratio*(v2(i)-v1(i))
   enddo
   if(d1 .gt. 0.0d0) then
    call tri_area8(x1,v1,v3,area)
    call tri_centroid(x1,v1,v3,cen)
   elseif(d2 .gt. 0.0d0)then
    call tri_area8(x1,v2,v3,area)
    call tri_centroid(x1,v2,v3,cen)    
   else
    print *,"error, volate d2*d3<0,488"
    stop
   endif

  
 elseif(d1 .ne. 0.0d0 .and. d2 .ne. 0.0d0 .and. d3 .ne. 0.0d0)then
 ! print *,"case6"
  ct = 0
  do i = 1,2
   if(d(i)*d(i+1) .lt. 0.0d0)then
    ct = ct + 1
    ratio = abs(d(i))/(abs(d(i))+abs(d(i+1)))
    do j = 1,2
     xx(ct,j)= v(i,j) + ratio*(v(i+1,j)-v(i,j))
    enddo 
   endif
  enddo
   
  if(ct .eq. 2) then
    ! do nothing
  elseif(ct .eq. 1) then
    if(d(1)*d(3) .ge. 0.0d0)then
     print *,"ct invalid 525"
     stop
    else
     ratio = abs(d(1))/(abs(d(1))+abs(d(3)))
     do j = 1,2
      xx(2,j)= v(1,j) + ratio*(v(3,j)-v(1,j))
     enddo     
    endif
  elseif(ct .eq. 0)then
   print *,"ct invalid 526"
   stop
  else
   print *,"ct invalid 529"
   stop
  endif
 
  if(d(1)*d(2) .lt. 0.0d0 .and. d(1)*d(3) .lt. 0.0d0)then
    call tri_area8(v1,xx(1,:),xx(2,:),area1)
    call tri_centroid(v1,xx(1,:),xx(2,:),cen1)
    if(d(1) .gt. 0.0d0)then
      area = area1
      do i = 1,2
       cen(i) = cen1(i)
      enddo
    elseif(d(1) .lt. 0.0d0)then
      call tri_area8(v1,v2,v3,area2)
      call tri_centroid(v1,v2,v3,cen2)
!    print *,"area2",area2
      area = area2-area1
      do i=1,2
       cen(i)=(cen2(i)*area2 - cen1(i)*area1)/area 
      enddo
    else
     print *,"d(1) cant be zero,  565"
     stop
    endif
  elseif(d(2)*d(3) .lt. 0.0d0 .and. d(2)*d(1) .lt. 0.0d0)then
    call tri_area8(v2,xx(1,:),xx(2,:),area1)
    call tri_centroid(v2,xx(1,:),xx(2,:),cen1)
    if(d(2) .gt. 0.0d0)then
      area = area1
      do i = 1,2
       cen(i) = cen1(i)
      enddo
    elseif(d(2) .lt. 0.0d0)then
      call tri_area8(v1,v2,v3,area2)
      call tri_centroid(v1,v2,v3,cen2)
 !   print *,"area2",area2
      area = area2-area1
      do i=1,2
       cen(i)=(cen2(i)*area2 - cen1(i)*area1)/area 
      enddo
    else
     print *,"d(1) can be zero,  565"
     stop
    endif

  elseif(d(3)*d(1) .lt. 0.0d0 .and. d(3)*d(2) .lt. 0.0d0)then
    call tri_area8(v3,xx(1,:),xx(2,:),area1)
    call tri_centroid(v3,xx(1,:),xx(2,:),cen1)
    if(d(3) .gt. 0.0d0)then
      area = area1
      do i = 1,2
       cen(i) = cen1(i)
      enddo
    elseif(d(3) .lt. 0.0d0)then
      call tri_area8(v1,v2,v3,area2)
      call tri_centroid(v1,v2,v3,cen2)
!    print *,"area2",area2
      area = area2-area1
      do i=1,2
       cen(i)=(cen2(i)*area2 - cen1(i)*area1)/area
      enddo
    else
     print *,"d(1) can be zero,  565"
     stop
    endif
  else
    print *,"err, check, 611",d
    print *,"imat=",im
    print *,"vertices:"
    print *,"v1",v1
    print *,"v2",v2
    print *,"v3",v3
    stop
  endif

 else
  print *,"err out of cases", d1, d2,d3

 endif
!else
! print *,"probtype_in error in triangle_interface_detect"
! stop
!endif






end subroutine triangle_interface_detect




SUBROUTINE TRI_AREA8(V1,V2,V3,AREA)
IMPLICIT NONE

real(kind=8),INTENT(IN)                  :: V1(2),V2(2),V3(2)
REAL(KIND=8),INTENT(OUT)                 :: AREA

! AREA = 0.5* ABS(V1(1)*(V2(2)-V3(2)) &
!           & - V2(1)*(V1(2)-V3(2)) &
!           & + V3(1)*(V1(2)-V2(2)))

 AREA=0.5*ABS((v2(1)-v1(1))*(v3(2)-v1(2)) - &
                (v3(1)-v1(1))*(v2(2)-v1(2)))


END SUBROUTINE tri_area8

subroutine tri_centroid(v1,v2,v3,cen)
implicit none

real(kind=8),intent(in)      :: v1(2),v2(2),v3(2)
real(kind=8)                 :: cen(2)
integer                      :: i


do i = 1,2
 cen(i) = 1.0d0/3.0d0*(v1(i)+v2(i)+v3(i))
enddo


end subroutine tri_centroid


!subroutine form_2p_line(v1,v2,mx,my,b)
!implicit none
! 
!real(kind=8),intent(in)  :: v1(2),v2(2)
!real(kind=8)             :: mx,my,b


!end subroutine form_2p_line




!--------------------------------------------------------
!------------------------------------------------------
! Linear reconstruction normal and initial guess
!        ax+by+d=0
!  input:  G ,X, Y
!------------------------------------------------------
SUBROUTINE LS_LSF(method_flag,h,center,Gin,aout,bout,dout)
IMPLICIT NONE

! flag = 1  least square fit method   flag = 0, basic way 

integer,intent(in)                          :: method_flag
real(kind=8),intent(in)                     :: h
REAL(KIND=8), DIMENSION(2,2)                :: w
REAL(KIND=8)                                :: w1,w2
REAL(KIND=8)                                :: aa,bb,cc,dd,ee,ff,gg,hh
REAL(KIND=8)                                :: aout, bout, dout
INTEGER                                     :: i,j
REAL(KIND=8), EXTERNAL                      :: NORM_2d

REAL(KIND=8),INTENT(IN)                     :: CENTER(2)
REAL(KIND=8)                                :: VERTEX(4,2)
REAL(kind=8),intent(in)                     :: Gin(4)
real(kind=8)                                :: G(2,2)
real(kind=8)                                :: x(2),y(2)


G(1,1) = Gin(1) 
G(1,2) = Gin(2)
G(2,1) = Gin(3)
G(2,2) = Gin(4)

CALL find_vertex(center, h, vertex)

X(1) = VERTEX(1,1)
X(2) = VERTEX(4,1)
Y(1) = VERTEX(1,2)
Y(2) = VERTEX(4,2)


if (method_flag .eq. 1) then
 aa=0.0d0
 bb=0.0d0
 cc=0.0d0
 dd=0.0d0
 ee=0.0d0
 ff=0.0d0
 gg=0.0d0
 hh=0.0d0


 do i= 1, 2
  do j= 1, 2
     aa = aa + 2.0d0*((x(i)-center(1))**2.0d0)
     bb = bb + 2.0d0*(x(i)-center(1))*(y(j)-center(2))
     cc = cc + 2.0d0*(x(i)-center(1))
     dd = dd + 2.0d0*((y(j)-center(2))**2.0d0)
     ee = ee + 2.0d0*(y(j)-center(2))
     ff = ff + 2.0d0*(x(i)-center(1))*G(i,j)
     gg = gg + 2.0d0*(y(j)-center(2))*G(i,j)
     hh = hh + 2.0d0*G(i,j)
  enddo
 enddo

!!write(20,*) aa,bb,cc,dd,ee,ff,gg,hh

 if(bb .eq. 0.0d0  .and. ee .eq. 0.0d0) THEN
   bout = gg/dd
   aout = (ff-cc*hh)/(aa-cc*cc)
   dout = hh - cc*aout
 else
  dout = ((aa*hh-ff*cc)/(aa*ee-bb*cc)-(aa*gg-bb*ff)/(aa*dd-bb*bb))/ &
         &((aa-cc*cc)/(aa*ee-bb*cc)-(aa*ee-bb*cc)/(aa*dd-bb*bb))
  bout = (aa*gg-bb*ff)/(aa*dd-bb*bb)-((aa*ee-bb*cc)/(aa*dd-bb*bb))*dout
  aout = (ff-cc*dout-bb*bout)/aa
 endif


elseif(method_flag .eq. 2) then
 



endif


END SUBROUTINE LS_LSF
!----------------------------------------------------------------------
subroutine inf_lsdetect(probtype_in,im_in,h,center,flag)
implicit none

integer,intent(in)        :: probtype_in
integer,intent(in)        :: im_in
real(kind=8),intent(in)   :: h
real(kind=8),intent(in)   :: center(2)
real(kind=8)              :: vertex(4,2)

integer                   :: i,j
real(kind=8)              :: dist(4)

integer                   :: flag

 flag = 0

 call find_vertex(center, h, vertex)
 do i = 1,4
  call dist_fns(im_in,vertex(i,1),vertex(i,2),dist(i),probtype_in)
 enddo
  
 do i = 1,3
  if(dist(1)*dist(i+1) .lt. 0.0d0)then
   flag = 1
   exit
  endif
 enddo

 do i = 2,3
  if(dist(2)*dist(i+1) .lt. 0.0d0)then
   flag = 1
   exit
  endif
 enddo

  if(dist(3)*dist(4) .lt. 0.0d0)then
   flag = 1
  endif


end subroutine inf_lsdetect



!------------------------------------------------------------------------

subroutine AdaptQuad_2d(iin,jin,nmat_in,dx,center,centroid,vf,probtype_in)
implicit none

integer,intent(in)            :: nmat_in,probtype_in
real(kind=8),intent(in)       :: dx(2)
real(kind=8),intent(in)       :: center(2)
integer     ,intent(in)       :: iin,jin

integer,parameter             :: ref_num = 5
integer                       :: i,j,im,i1,i2,i3,i4,i5,dir,ilev
integer                       :: ii,jj,ic
integer                       :: flag,sign_flag
integer                       :: vcheck, lev_check   ! flags
real(kind=8)                  :: dist,h2,h3,h4,h5,h6
real(kind=8)                  :: x,y
real(kind=8)                  :: h
integer                       :: ncenter
real(kind=8)                  :: centers1(4,2),centers2(4,2)
real(kind=8)                  :: centers3(4,2),centers4(4,2),centers5(4,2)
real(kind=8)                  :: vertex(4,2)
real(kind=8)                  :: G(4)
real(kind=8)                  :: mx,my,alphadump,c

real(kind=8)                  :: vf(nmat_in)
real(kind=8)                  :: centroid(nmat_in,2)
real(kind=8)                  :: cxtemp(nmat_in),cytemp(nmat_in)
real(kind=8)                  :: vol(nmat_in)
type(polygon)                 :: cell
type(polygon)                 :: POS_PLG,NEG_PLG
type(points)                  :: cen_temp
real(kind=8)                  :: vol_temp

! check
integer                       :: t1,t2
real(kind=8)                  :: g_center,gtemp
real(kind=8)                  :: vf_check

type(points)                  :: seg(2)

real(kind=8)                  :: v1(2),v2(2),v3(2)

real(kind=8),allocatable      :: center_hold(:,:),center_hold1(:,:)
integer,allocatable           :: center_mark(:)

integer                       :: scheck


if(iin .eq. 14 .and. jin .eq. 35)then
 scheck=1
else
 scheck=0
endif


x = center(1)
y = center(2)
h = dx(1)                  !  h squared domain



if (h.le.0.0) then
 print *,"h invalid"
 stop
endif


 vf = 0.0d0
 centroid = 0.0d0
 vol = 0.0d0
 cxtemp = 0.0d0
 cytemp = 0.0d0

! inf_lsdetect(probtype_in,im_in,h,center,flag)

do im=1,nmat_in
 ilev=0
 lev_check = 0

 ncenter = 1
 allocate(center_hold(ncenter,2))
 allocate(center_mark(ncenter))
 do i=1,2
  center_hold(1,i) = center(i)
 enddo
 lev_check=1

 do while(lev_check .eq. 1 .and. ilev .lt. refine_lev+1)
  lev_check=0
  do i=1,ncenter
   call inf_lsdetect(probtype_in,im,h/(2.0d0**real(ilev,8)),center_hold(i,:),vcheck)
   if(vcheck .eq. 0)then
    call dist_fns(im,center_hold(i,1),center_hold(i,2),dist,probtype_in)
    call update_volncen(nmat_in,dist,im,center_hold(i,:), &
                        h/(2.0d0**real(ilev,8)),vol,cxtemp,cytemp)
    center_mark(i)=0
   elseif(vcheck .eq. 1)then
!     print *,"i",iin,"j",jin
    center_mark(i)=1
    lev_check=1
   else
    print *,"vcheck flag invalid 1191"
    stop
   endif
  enddo

!  if(scheck .eq. 1)then
!   write(11,*) "ncheck",ncenter
!   write(11,*) "vol",vol,vol/(h*h)
!   write(11,*) "center_mark", center_mark
!  endif

  ilev=ilev+1
  if(lev_check .eq. 1)then
   allocate(center_hold1(4**(ilev),2))
   ic=0 
   do i=1,ncenter
    if(center_mark(i) .eq. 1)then
     call AdaptQuad_sub(center_hold(i,:), h/(2.0d0**real(ilev-1,8)), centers1)
     do ii=1,4
      ic=ic+1
      do j=1,2
       center_hold1(ic,j)=centers1(ii,j)
      enddo
     enddo
    endif
   enddo
   deallocate(center_hold)
   deallocate(center_mark)
   allocate(center_hold(ic,2))
   allocate(center_mark(ic))
   do i=1,ic
   do j=1,2
    center_hold(i,j) = center_hold1(i,j)
   enddo
   enddo
   deallocate(center_hold1)
   ncenter=ic
  elseif(lev_check .eq. 0)then
   ! do nothing
  else
   print *,"lev_check flag invalid 1225"
   stop
  endif
   
 enddo  ! do while

!  if(scheck .eq. 1)then
!   write(11,*) "ncenter",ncenter
!   write(11,*) "vol",vol,vol/(h*h)
!   write(11,*) "center_mark", center_mark!
! 
!   write(11,*) "ic",ic
!   do i=1,ic
!    write(11,*) "center_hold", center_hold(i,:)
!   enddo

!  endif



 if(lev_check .ne. 0)then
  do i=1,ic
   call find_vertex(center_hold(i,:), h/(2.0d0**real(ilev,8)), vertex)

     ! triangulation
     v1 = vertex(1,:)
     v2 = vertex(2,:)
     v3 = vertex(4,:)

    if(scheck .eq. 1)then
     call dist_fns(im,vertex(1,1),vertex(1,2),dist,probtype_in)
     write(11,*) "v1",vertex(1,:),dist
     call dist_fns(im,vertex(2,1),vertex(2,2),dist,probtype_in)
     write(11,*) "v2",vertex(2,:),dist
     call dist_fns(im,vertex(3,1),vertex(3,2),dist,probtype_in)
     write(11,*) "v3",vertex(3,:),dist
     call dist_fns(im,vertex(4,1),vertex(4,2),dist,probtype_in)
     write(11,*) "v4",vertex(4,:),dist
    endif



     call triangle_interface_detect(im,probtype_in,v1,v2,v3,vol_temp,cen_temp%val)
     ! update centroid and volume
     vol(im) = vol(im) + vol_temp
     cxtemp(im) = cxtemp(im) + vol_temp*cen_temp%val(1) 
     cytemp(im) = cytemp(im) + vol_temp*cen_temp%val(2)
 

    if(scheck .eq. 1)then
     write(11,*) "vol",vol,vol/(h*h)
     write(11,*) "vol_temp",vol_temp,vol_temp/(h*h)
    endif


     !------------------------------------------------------------------------------
     vol_temp = 0.0d0
     cen_temp%val = 0.0d0
     v1 = vertex(1,:)
     v2 = vertex(4,:)
     v3 = vertex(3,:)
     call triangle_interface_detect(im,probtype_in,v1,v2,v3,vol_temp,cen_temp%val)

     ! update centroid and volume
     vol(im) = vol(im) + vol_temp
     cxtemp(im) = cxtemp(im) + vol_temp*cen_temp%val(1) 
     cytemp(im) = cytemp(im) + vol_temp*cen_temp%val(2)

    if(scheck .eq. 1)then
     write(11,*) "vol",vol,vol/(h*h)
     write(11,*) "vol_temp",vol_temp
    endif

  enddo
 endif

 deallocate(center_hold)
 deallocate(center_mark) 

enddo ! im


if (h.le.0.0) then
 print *,"h invalid"
 stop
endif

do im = 1,nmat_in
  if(vol(im) .lt. 0.0d0)then
    print *,"vol(im) is negtive, 2155"
    stop
  elseif(vol(im) .eq. 0.0)then
    centroid(im,:)=center(:)
  else
   centroid(im,1) = cxtemp(im)/vol(im)
   centroid(im,2) = cytemp(im)/vol(im)
  endif
   vf(im) = vol(im)/(h*h)     
enddo


  if (scheck .eq. 1) then
   write(11,*) iin,jin
   write(11,*) "center",center
   call find_vertex(center, h, vertex)
   do ii=1,4
    call dist_fns(1,vertex(ii,1),vertex(ii,2),dist,probtype_in)
    write(11,*) "mat1, ii,vertex,dist ",ii,vertex(ii,1),vertex(ii,2),dist
   enddo
   do ii=1,4
    call dist_fns(2,vertex(ii,1),vertex(ii,2),dist,probtype_in)
    write(11,*) "mat2, ii,vertex,dist ",ii,vertex(ii,1),vertex(ii,2),dist
   enddo

    write(11,*) "vfrac ",vf(1),vf(2)
    write(11,*) "mat1 ",vf(1),centroid(1,:)
    write(11,*) "mat2" , vf(2),centroid(2,:)
   endif



call vf_correct(iin,jin,nmat_in,vf,probtype_in)

if (probtype_in.eq.0) then                                          ! 0
 ! do nothing
else if (probtype_in.eq.2) then                                     ! 2
 ! do nothing
elseif(probtype_in .eq. 4 )then                                      ! 4
  ! do nothing
elseif(probtype_in .eq. 5)then                                       ! 5
 ! vf(2)=1.0-vf(1)
else if (probtype_in.eq.1 .or. probtype_in .eq. 3 &                  ! 1,3,9
         .or. probtype_in .eq. 9) then
 vf(2) = 1.0d0 - vf(1) -vf(3)

 if(abs(vf(2)) .lt. eps)then
  vf(2) = 0.0d0
 endif
 if(vf(2) .lt. -eps) then
  print *,"vf(2) is negative"
 endif

elseif(probtype_in .eq. 6)then                                     ! 6
 vf(3) = 1.0d0 - vf(1) -vf(2)

 if(abs(vf(3)) .lt. eps)then
  vf(3) = 0.0d0
 endif
 if(vf(3) .lt. -eps) then
  print *,"vf(3) is negative"
 endif

elseif(probtype_in .eq. 7)then                                      ! 7

 if(center(1) .gt. 0.5d0 .and. center(2) .gt. 0.5d0)then
  vf(2) = 1.0d0-vf(1)
 elseif(center(1) .lt. 0.5d0 .and. center(2) .gt. 0.5d0)then
  vf(3) = 1.0d0-vf(1)  
 elseif(center(1) .lt. 0.5d0 .and. center(2) .lt. 0.5d0)then
  vf(4) = 1.0d0-vf(1)  
 elseif(center(1) .gt. 0.5d0 .and. center(2) .lt. 0.5d0)then
  vf(5) = 1.0d0-vf(1)  
 else
  print *,"center is not aligned with grid"
  stop
 endif


! do nothing
 if(1 .eq. 0)then
 if(vf(1) .ne. 0.0d0 .and. vf(2) .ne. 0.0d0 .and. vf(3) .ne. 0.0d0  &
    .and. vf(4) .eq. 0.0d0 .and. vf(5) .eq. 0.0d0)then
   vf(3)=1.0d0-vf(1)-vf(2)
   do dir=1,2
    centroid(3,dir) =  &
     (center(dir) - vf(1)*centroid(1,dir) - vf(2)*centroid(2,dir))/vf(3)
   enddo
 elseif(vf(1) .ne. 0.0d0 .and. vf(2) .ne. 0.0d0 .and. vf(5) .ne. 0.0d0 &
    .and. vf(4) .eq. 0.0d0 .and. vf(3) .eq. 0.0d0)then
   vf(5)=1.0d0-vf(1)-vf(2)
   do dir=1,2
    centroid(5,dir) =  &
     (center(dir) - vf(1)*centroid(1,dir) - vf(2)*centroid(2,dir))/vf(5)
   enddo
 elseif(vf(1) .ne. 0.0d0 .and. vf(4) .ne. 0.0d0 .and. vf(3) .ne. 0.0d0 &
    .and. vf(2) .eq. 0.0d0 .and. vf(5) .eq. 0.0d0)then
   vf(4)=1.0d0-vf(1)-vf(3)
   do dir=1,2
    centroid(4,dir) =  &
     (center(dir) - vf(1)*centroid(1,dir) - vf(3)*centroid(3,dir))/vf(4)
   enddo
 elseif(vf(1) .ne. 0.0d0 .and. vf(4) .ne. 0.0d0 .and. vf(5) .ne. 0.0d0 &
    .and. vf(2) .eq. 0.0d0 .and. vf(3) .eq. 0.0d0)then
   vf(4)=1.0d0-vf(1)-vf(5)
   do dir=1,2
    centroid(4,dir) =  &
     (center(dir) - vf(1)*centroid(1,dir) - vf(5)*centroid(5,dir))/vf(4)
   enddo
 else
   ! do nothing
 endif
 endif

elseif(probtype_in .eq. 10)then                                     

 if(center(1) .gt. 0.5d0 .and. center(2) .gt. 0.5d0)then
  vf(2) = 1.0d0-vf(1)-vf(6)
 elseif(center(1) .lt. 0.5d0 .and. center(2) .gt. 0.5d0)then
  vf(3) = 1.0d0-vf(1)-vf(6)  
 elseif(center(1) .lt. 0.5d0 .and. center(2) .lt. 0.5d0)then
  vf(4) = 1.0d0-vf(1)-vf(6)  
 elseif(center(1) .gt. 0.5d0 .and. center(2) .lt. 0.5d0)then
  vf(5) = 1.0d0-vf(1)-vf(6)  
 else
  print *,"center is not aligned with grid"
  stop
 endif





else
 print *,"probtype_in invalid"
 stop
endif

! sanity check
do im = 1,nmat_in

   if (abs(vf(im)).le.eps) then
    vf(im)=0.0
   else if (abs(vf(im)-1.0).le.eps) then
    vf(im)=1.0
   else if ((vf(im).gt.0.0).and.(vf(im).lt.1.0)) then
    ! do nothing
   else
    print *,im,vf(im)
    print *, "Error in vf 2" 
    stop
   endif
  
enddo ! im


if (probtype_in.eq.0) then
 ! do nothing
else if (probtype_in.eq.2) then
 ! do nothing
elseif(probtype_in .eq. 4)then
  ! do nothing
elseif(probtype_in .eq. 5) then
 !  do dir=1,2
 !   centroid(2,dir) =  &
 !   (center(dir) - vf(1)*centroid(1,dir))/vf(2)
 !  enddo  

else if (probtype_in.eq.1 .or. probtype_in .eq. 3 &
         .or. probtype_in .eq. 9) then
 if(vf(2) .gt. eps)then
  do dir=1,2
   centroid(2,dir) =  &
    (center(dir) - vf(1)*centroid(1,dir) - vf(3)*centroid(3,dir))/vf(2)
  enddo
 endif
elseif(probtype_in .eq. 6)then
 if(vf(3) .gt. eps)then
  do dir=1,2
   centroid(3,dir) =  &
    (center(dir) - vf(1)*centroid(1,dir) - vf(2)*centroid(2,dir))/vf(3)
  enddo
 endif
else if (probtype_in.eq.7) then
 if(vf(2) .gt. eps .and. vf(2) .lt. 1.0d0)then
  do dir=1,2
   centroid(2,dir) =  &
    (center(dir) - vf(1)*centroid(1,dir))/vf(2)
  enddo
 elseif(vf(3) .gt. eps .and. vf(3) .lt. 1.0d0)then
  do dir=1,2
   centroid(3,dir) =  &
    (center(dir) - vf(1)*centroid(1,dir))/vf(3)
  enddo
 elseif(vf(4) .gt. eps .and. vf(4) .lt. 1.0d0)then
  do dir=1,2
   centroid(4,dir) =  &
    (center(dir) - vf(1)*centroid(1,dir))/vf(4)
  enddo
 elseif(vf(5) .gt. eps .and. vf(5) .lt. 1.0d0)then
  do dir=1,2
   centroid(5,dir) =  &
    (center(dir) - vf(1)*centroid(1,dir))/vf(5)
  enddo 

 endif

else if (probtype_in.eq.10) then
 if(vf(2) .gt. eps .and. vf(2) .lt. 1.0d0)then
   do dir=1,2
   centroid(2,dir) =  &
    (center(dir) - vf(1)*centroid(1,dir)-vf(6)*centroid(6,dir))/vf(2)
  enddo
 elseif(vf(3) .gt. eps .and. vf(3) .lt. 1.0d0)then
  do dir=1,2
   centroid(3,dir) =  &
    (center(dir) - vf(1)*centroid(1,dir)-vf(6)*centroid(6,dir))/vf(3)
  enddo
 elseif(vf(4) .gt. eps .and. vf(4) .lt. 1.0d0)then
  do dir=1,2
   centroid(4,dir) =  &
    (center(dir) - vf(1)*centroid(1,dir)-vf(6)*centroid(6,dir))/vf(4)
  enddo
 elseif(vf(5) .gt. eps .and. vf(5) .lt. 1.0d0)then
  do dir=1,2
   centroid(5,dir) =  &
    (center(dir) - vf(1)*centroid(1,dir)-vf(6)*centroid(6,dir))/vf(5)
  enddo 

 endif


else
 print *,"probtype_in invalid"
 stop
endif

!do im = 1,nmat_in
!   if(im .eq. 1)then
!    write(4,*) centroid(im,:)
!   elseif(im .eq. 2)then
!    write(5,*) centroid(im,:)
!   elseif(im .eq. 3)then
!    write(21,*) centroid(im,:)
!   elseif(im .eq. 4)then
!    write(22,*) centroid(im,:) 
!   elseif(im .eq. 5)then
!    write(23,*) centroid(im,:) 
! 
!   else
!    print *,"err 1075"
!    stop
!   endif   
!enddo


! call renormalize_vf(nmat_in,vf) 


end subroutine AdaptQuad_2d



!-------------------------------------------------------
subroutine check_sign(G,gtemp,flag)
implicit none

real(kind=8)         :: g(4)
real(kind=8)         :: gtemp
integer              :: flag
integer              :: i,j

flag = 0
do i = 1,4
  if(abs(g(i)) .lt. 1.0e-12) then
     g(i)  = 0.0d0
  endif
enddo

do i = 1,4
   do j = 1,4
      if(g(i)*g(j) .lt. 0.0d0) then
        flag = 1
      endif
   enddo      
enddo     

do i = 1,4
  if(g(i)*gtemp .lt. 0.0d0) then
     flag = 1
  endif
enddo

end subroutine check_sign





! -----------------------------------------------
!---------------------------------------------------------------------
subroutine vf_correct(iin,jin,nmat_in,vf,probtype_in)
implicit none

integer       ,intent(in)  :: nmat_in,iin,jin,probtype_in
real(kind=8)               :: vf(nmat_in)
integer                    :: i,j
real(kind=8)               :: vcheck


if ((probtype_in.eq.0).or. &
    (probtype_in.eq.2) .or. &
    (probtype_in .eq. 4) .or. &
    (probtype_in .eq. 5)) then

 if (nmat_in.ne.2) then
  print *,"nmat_in invalid"
  stop
 endif
 
 vcheck = vf(1) + vf(2)
 if(vf(1) .lt. 0.0d0 .or. vf(2) .lt. 0.0d0) then
  print *,"vf is negative"
 endif

 if(vcheck .gt. 1.0d0+eps) then
  print *,"goes into vf_correct1",vcheck
  stop
 endif

else if (probtype_in.eq.1 .or. probtype_in .eq. 3 &
         .or. probtype_in .eq. 9) then
 if (nmat_in.ne.3) then
  print *,"nmat_in invalid"
  stop
 endif
 
 vcheck = vf(1) + vf(3)
 if(vf(1) .lt. 0.0d0 .or. vf(3) .lt. 0.0d0) then
  print *,"vf is negative"
 endif

 if(vcheck .gt. 1.0d0+eps) then
  print *,"goes into vf_correct2",vcheck
  if(vf(1) .gt. 1.0d0) then
    print *,"case1"
     vf(1) = 1.0d0
     vf(3) = 0.0d0
  elseif(vf(3) .gt. 1.0d0)then
    print *,"case2"
     vf(1) = 0.0d0
     vf(3) = 1.0d0
  else
    print *,"case3",iin,jin
    print *,"overshoot", vf(1) , vf(3)
    if(vf(1) .gt. vf(3))then
       vf(1) = 1.0d0
       vf(3) = 0.0d0
    elseif(vf(1) .lt. vf(3))then
       vf(1) = 0.0d0
       vf(3) = 1.0d0
    endif
  endif
 endif
elseif(probtype_in .eq. 6)then
 vcheck = vf(1) + vf(2)
 if(vf(1) .lt. 0.0d0 .or. vf(2) .lt. 0.0d0) then
  print *,"vf is negative"
 endif

 if(vcheck .gt. 1.0d0+eps) then
  print *,"goes into vf_correct3",vcheck
  if(vf(1) .gt. 1.0d0) then
    print *,"case1"
     vf(1) = 1.0d0
     vf(2) = 0.0d0
  elseif(vf(2) .gt. 1.0d0)then
    print *,"case2"
     vf(1) = 0.0d0
     vf(2) = 1.0d0
  else
    print *,"case3",iin,jin
    print *,"overshoot", vf(1) , vf(2)
    if(vf(1) .gt. vf(2))then
       vf(1) = 1.0d0
       vf(2) = 0.0d0
    elseif(vf(1) .lt. vf(2))then
       vf(1) = 0.0d0
       vf(2) = 1.0d0
    endif
  endif
 endif
elseif(probtype_in .eq. 7)then
! vcheck = vf(1) + vf(2)+vf(3)+vf(4)+vf(5)
!  if(vcheck .gt. 1.0d0+eps) then
!  print *,"goes into vf_correct4",vcheck,vf
! endif
 do i=1,5
  if(vf(i) .lt. 0.0d0)then
   print *,"vf is negtive"
   stop
  endif
 enddo

 if(vf(1) .gt. 1.0d0)then
  vf(1)=1.0d0
  vf(2)=0.0d0
  vf(3)=0.0d0
  vf(4)=0.0d0
  vf(5)=0.0d0
 elseif(vf(2) .gt. 1.0d0)then
  vf(1)=0.0d0
  vf(2)=1.0d0
  vf(3)=0.0d0
  vf(4)=0.0d0
  vf(5)=0.0d0
 elseif(vf(3) .gt. 1.0d0)then
  vf(1)=0.0d0
  vf(2)=0.0d0
  vf(3)=1.0d0
  vf(4)=0.0d0
  vf(5)=0.0d0
 elseif(vf(4) .gt. 1.0d0)then
  vf(1)=0.0d0
  vf(2)=0.0d0
  vf(3)=0.0d0
  vf(4)=1.0d0
  vf(5)=0.0d0
 elseif(vf(5) .gt. 1.0d0)then
  vf(1)=0.0d0
  vf(2)=0.0d0
  vf(3)=0.0d0
  vf(4)=0.0d0
  vf(5)=1.0d0
 else
  ! do nothing
   
 endif
elseif(probtype_in .eq. 10)then
! vcheck = vf(1) + vf(2)+vf(3)+vf(4)+vf(5)
!  if(vcheck .gt. 1.0d0+eps) then
!  print *,"goes into vf_correct4",vcheck,vf
! endif

 if (nmat_in.ne.6) then
  print *,"nmat_in invalid"
  stop
 endif
 
 vcheck = vf(1) + vf(6)
 if(vf(1) .lt. 0.0d0 .or. vf(6) .lt. 0.0d0) then
  print *,"vf is negative"
 endif

 if(vcheck .gt. 1.0d0+eps) then
  print *,"goes into vf_correct2",vcheck
  if(vf(1) .gt. 1.0d0) then
    print *,"case1"
     vf(1) = 1.0d0
     vf(6) = 0.0d0
  elseif(vf(6) .gt. 1.0d0)then
    print *,"case2"
     vf(1) = 0.0d0
     vf(6) = 1.0d0
  else
    print *,"case3",iin,jin
    print *,"overshoot", vf(1) , vf(3)
    if(vf(1) .gt. vf(6))then
       vf(1) = 1.0d0
       vf(6) = 0.0d0
    elseif(vf(1) .lt. vf(6))then
       vf(1) = 0.0d0
       vf(6) = 1.0d0
    endif
  endif
 endif

else
 print *,"probtype_in invalid"
 stop
endif

end subroutine vf_correct



!--------------------------------------------
SUBROUTINE init_vfncen(N,cell,nmat_in,dx,centroid,vf,probtype_in)
IMPLICIT NONE

integer,intent(in)     :: N,nmat_in,probtype_in
type(polygon),intent(in)  :: cell(-1:N,-1:N)
real(kind=8),intent(in)  :: dx(2)
TYPE(POINTS),DIMENSION(-1:N,-1:N,nmat_in)  :: CENTROID
real(kind=8)                       :: center(2)
real(kind=8)                       :: cen_temp(nmat_in,2)
integer                            :: i,j,im,dir
real(kind=8)                       :: vf_temp(nmat_in)
real(kind=8)                       :: vf(-1:N,-1:N,nmat_in)

 print *,"in init_vfncen "
 print *,"N,nmat_in ",N,nmat_in
 print *,"dx ",dx(1),dx(2) 

 cen_temp = 0.0d0 ! (nmat_in,2)
 vf_temp = 0.0d0  ! (nmat_in)
 vf = 0.0d0
 do i = -1,N
  do j = -1,N
   do im= 1,nmat_in  
      centroid(i,j,im)%val = 0.0d0
   enddo
  enddo
 enddo

 do i = -1,N
   do j = -1,N
    do dir=1,2
     center(dir) = cell(i,j)%center%val(dir)
    enddo
   
!    write(3,*) "center",center

    call AdaptQuad_2d(i,j,nmat_in,dx,center,cen_temp,vf_temp,probtype_in)
    
    

    do im = 1, nmat_in
     do dir=1,2
      centroid(i,j,im)%val(dir) = cen_temp(im,dir)
     enddo
     vf(i,j,im) = vf_temp(im)

     if (1.eq.0) then
      write(3,*) "i,j,im,vf,cenx,ceny ",i,j,im,vf_temp(im), &
        cen_temp(im,1),cen_temp(im,2)
     endif  
    enddo ! im

   enddo
 enddo

END SUBROUTINE init_vfncen
!---------------------------------------------------------------------------

!--------------------------------------------------------------
subroutine init_mofdata(N,sdim,dx,nmat_in,nten,CELL,vf,CENTROID,mofdata)
implicit none

integer,intent(in)               :: N,sdim,nmat_in
TYPE(POLYGON)                    :: CELL(-1:N,-1:N)
real(kind=8),dimension(-1:N,-1:N,nmat_in) :: vf
real(kind=8),intent(in)          :: dx(sdim)
TYPE(POINTS),DIMENSION(-1:N,-1:N,nmat_in) :: CENTROID

real(kind=8)                     :: mofdata(-1:N,-1:N,(sdim*2+3)*nmat_in)

integer                          :: i,j,im,i1,j1,ii,jj,dir 
integer                          :: bfact,mof_verbose
integer,intent(in)               :: nten
integer                          :: nhalf0
real(kind=8)                     :: xsten0(-3:3,sdim)
REAL(kind=8)                     :: LS_stencil(D_DECL(-1:1,-1:1,-1:1),nmat_in+nten)
real(kind=8)                     :: xtetlist_vof(sdim+1,sdim,2000)
real(kind=8)                     :: xtetlist_cen(sdim+1,sdim,2000)
real(kind=8)                     :: multi_centroidA(nmat_in,sdim)
real(kind=8)                     :: time
integer                          :: ngeom_recon_in,vofcomp

real(kind=8)                     :: mx,my,c
type(points)                     :: seg(2)

print *,"in init_mofdata"
print *,"N,sdim,nmat_in,nten = ",N,sdim,nmat_in,nten
print *,"dx= ",dx(1),dx(2)

bfact = 1
nhalf0 = 3
do i = -3,3
  do dir = 1,sdim
    xsten0(i,dir) = i*dx(dir)*0.5d0
  enddo
enddo
mof_verbose = 0
LS_stencil = 0.0d0
xtetlist_vof = 0.0d0
xtetlist_cen = 0.0d0
ngeom_recon_in = 2*sdim+3

mofdata = 0.0d0

do  i = 0,N-1
 do j = 0,N-1
   do im = 1,nmat_in
    if (abs(vf(i,j,im)).le.eps) then
     vf(i,j,im)=0.0
    else if (abs(vf(i,j,im)-1.0).le.eps) then
     vf(i,j,im)=1.0
    else if ((vf(i,j,im).gt.0.0).and.(vf(i,j,im).lt.1.0)) then
     ! do nothing
    else
     print *,"vf invalid vf=",vf(i,j,im)
     print *,"i,j,im,N,nmat_in,sdim,ngeom_recon_in ",i,j,im,N,nmat_in,sdim, &
      ngeom_recon_in
     stop
    endif

    vofcomp=(im-1)*ngeom_recon_in+1 
    mofdata(i,j,vofcomp) = vf(i,j,im)
    
    if(vf(i,j,im) .gt. 0.0 .and. vf(i,j,im) .lt. 1.0d0)then

     do dir = 1,sdim
      mofdata(i,j,dir+vofcomp) = &
       CENTROID(i,j,im)%val(dir)-CELL(i,j)%center%val(dir)
     enddo

    else if ((vf(i,j,im).eq.0.0).or.(vf(i,j,im).eq.1.0)) then
     ! do nothing
    else
     print *,"invalid volume fraction"
     stop   
    endif
    
    mofdata(i,j,3+vofcomp) = 0.0d0  
    mofdata(i,j,4+vofcomp) = 0.0d0  
    mofdata(i,j,5+vofcomp) = 0.0d0  
    mofdata(i,j,6+vofcomp) = 0.0d0  
   enddo ! im
   
!        print *, "i=", i, "j=", j
!        print *, mofdata(i,j,:)

        ! bfact = 1   xsten0 =     nhalf0 = 3
    call multimaterial_MOF( &
        bfact,dx,xsten0,nhalf0, &
        mof_verbose, &
        0, &
        LS_stencil, &
        xtetlist_vof, &
        xtetlist_cen, &
        2000, &
        mofdata(i,j,:), &
        multi_centroidA, &
        time, &
        0, &
        0,nmat_in,nten,sdim, &
        ngeom_recon_in)

 enddo
enddo



end subroutine init_mofdata

!-------------------------------------------------------------------

subroutine get_filament_source(x_in,t_in,probtype_in,im,sdim,G_in)
IMPLICIT NONE

integer, intent(in) :: probtype_in
integer, intent(in) :: im
integer, intent(in) :: sdim
REAL*8, intent(in) :: x_in(sdim)
REAL*8, intent(in) :: t_in
REAL*8, intent(out) :: G_in
REAL*8 radius_in,theta_in,r1,r2,delx,dely
REAL*8 mypi

 mypi=4.0d0*atan(1.0d0)
 if (sdim.ne.2) then
  print *,"sdim invalid"
  stop
 endif
 if (t_in.lt.0.0) then
  print *,"t invalid"
  stop
 endif

 if (probtype_in.eq.0) then
  if (im.eq.1) then
   G_in=0.0
  else if (im.eq.2) then
   G_in=0.0
  else
   print *,"im invalid 3"
   print *,"im=",im
   print *,"sdim=",sdim
   print *,"x,y,t ",x_in(1),x_in(2),t_in
   print *,"probtype_in ",probtype_in
   stop
  endif
 else if (probtype_in.eq.2) then
  if (im.eq.1) then
   G_in=0.0
  else if (im.eq.2) then
   G_in=0.0
  else
   print *,"im invalid 4"
   stop
  endif
 else if (probtype_in.eq.1) then
  if ((im.eq.1).or.(im.eq.3)) then
   G_in=0.0
  else if (im.eq.2) then

   r1=radcen-radeps
   r2=radcen+radeps

   delx=x_in(1)-0.5d0
   dely=x_in(2)-0.5d0


   radius_in = sqrt(delx**2.0d0 +dely**2.0d0)

    ! x=r cos(theta)
    ! y=r sin(theta)
   if (radius_in.le.radeps/1000.0) then
    theta_in=0.0
   else if ((delx.ge.0.0).and.(dely.ge.0.0)) then
    theta_in=acos(delx/radius_in)
   else if ((delx.le.0.0).and.(dely.ge.0.0)) then
    theta_in=acos(abs(delx)/radius_in)
    theta_in=mypi-theta_in
   else if ((delx.le.0.0).and.(dely.le.0.0)) then
    theta_in=acos(abs(delx)/radius_in)
    theta_in=mypi+theta_in
   else if ((delx.ge.0.0).and.(dely.le.0.0)) then
    theta_in=acos(delx/radius_in)
    theta_in=2.0d0*mypi-theta_in
   else
    print *,"delx or dely invalid,1",x_in
    stop
   endif

   ! T=2+sin(theta) exp(-t/(rc^2))  alpha=1
   ! T_t - (T_rr + T_r/r + T_theta theta/r^2)=
   ! exp(-t/rc^2)(-sin(theta)/rc^2+sin(theta)/r^2)
   G_in=exp(-t_in/(radcen**2))*sin(theta_in)*(-1.0/(radcen**2)+1.0/(radius_in**2))
  else
   print *,"im invalid 5"
   stop
  endif

 else if (probtype_in.eq.3) then

!  if(1 .eq. 0)then
!  if ((im.eq.1).or.(im.eq.3)) then
!   G_in=0.0
!  else if (im.eq.2) then
!    delx=x_in(1)-0.5d0-0.02d0*sqrt(5.0)
!    dely=x_in(2)-0.5d0-0.02d0*sqrt(5.0)  
!   radius_in = sqrt(delx**2.0d0 +dely**2.0d0)

    ! x=r cos(theta)
    ! y=r sin(theta)
!   if (radius_in.le.radeps/1000.0) then
!    theta_in=0.0
!   else if ((delx.ge.0.0).and.(dely.ge.0.0)) then
!    theta_in=acos(delx/radius_in)
!   else if ((delx.le.0.0).and.(dely.ge.0.0)) then
!    theta_in=acos(abs(delx)/radius_in)
!    theta_in=mypi-theta_in
!   else if ((delx.le.0.0).and.(dely.le.0.0)) then
!    theta_in=acos(abs(delx)/radius_in)
!    theta_in=mypi+theta_in
!   else if ((delx.ge.0.0).and.(dely.le.0.0)) then
!    theta_in=acos(delx/radius_in)
!    theta_in=2.0d0*mypi-theta_in
!   else
!    print *,"delx or dely invalid"
!    stop
!   endif
   ! T=2+sin(theta) exp(-t)  alpha=1
   ! T_t - (T_rr + T_r/r + T_theta theta/r^2)=
   ! exp(-t/rc^2)(-sin(theta)/rc^2+sin(theta)/r^2)
!   G_in=exp(-t_in)*sin(theta_in)*(-1.0d0+1.0d0/(radius_in**2))
!  else
!   print *,"im invalid 5"
!   stop
!  endif
!  endif

!  if ((im.eq.1).or.(im.eq.3)) then
!   G_in=0.0
!  else if (im.eq.2) then
!   G_in= (-4.0d0-(x_in(1)**2.0d0+x_in(2)**2.0d0))*exp(-t_in)
!  endif

   G_in = 0.0d0

 elseif(probtype_in .eq.4 )then
   G_in = (-4.0d0 -(x_in(1)**2.0d0+x_in(2)**2.0d0))*exp(-t_in)

!  if(im .eq. 1)then
!   G_in = 4.0d0  
!  elseif(im .eq. 2)then
!   G_in = 16.0d0*(x_in(1)*x_in(1) + x_in(2)*x_in(2))
!  else
!   print *,"wrong im 1509"
!   stop
!  endif

! elseif(probtype_in .eq.5)then
!  if(im .eq. 1)then
!   G_in = 4.0d0  
!  elseif(im .eq. 2)then
!   G_in = 16.0d0*(x_in(1)*x_in(1) + x_in(2)*x_in(2))
!  else
!   print *,"wrong im 1509"
!   stop
!  endif
 elseif(probtype_in .eq. 5)then
   G_in = 0.0d0
 elseif(probtype_in .eq. 6)then
  ! do nothing
  !G_in = (4.0d0 -(x_in(1)**2.0d0+x_in(2)**2.0d0))*exp(-t_in)
   G_in = 0.0d0
 elseif(probtype_in .eq. 7)then
   G_in = 0.0d0
 elseif(probtype_in .eq. 9)then
   G_in = 0.0d0
 elseif(probtype_in .eq. 10)then
   G_in = 0.0d0
 else
  print *,"probtype_in invalid"
  stop
 endif

return
end subroutine get_filament_source


subroutine set_polar_2d(sdim,N,M,kappa,tau,r,z,dr,dz,u)
implicit none


integer,intent(in)          :: sdim                 
integer,intent(in)          :: N      ! discretization in r direction
integer,intent(in)          :: M       ! discretization in theta direction
!real(kind=8),parameter      :: rlo=radcen-radeps
!real(kind=8),parameter      :: rhi=radcen+radeps

real(kind=8)                :: r(0:N)
real(kind=8)                :: z(0:M)

real(kind=8)                :: dr,dz
real(kind=8)                :: tau
real(kind=8)                :: kappa

real(kind=8)                :: u(0:N,0:M)

integer                     :: i,j


if(sdim .ne. 2)then
 print *,"invalid dimension 2909"
 stop
endif
do i=0,N
do j=0,M
 u(i,j)=0.0d0
enddo
enddo

dr=(rhi-rlo)/N
dz=(2*pi)/M

do i=0,N
 r(i)=rlo+i*dr
enddo
do i=0,M
 z(i)=i*dz
enddo

tau=0.5d0*kappa*min(((dr)**2.0d0), &
       minval(r)*dr, &
      (minval(r)**2.0d0)*((dz)**2.0d0))


!do i=0,N                                                              
! do j=0,M                                                             
!  u(i,j)=2.0d0                                                      
!  u(i,j) = 1.0d0 + 10.0d0*(r(i)-rlo)                                 
! enddo                                                               
!enddo                                                                

do i=0,N
 do j=0,M
  u(i,j)=BC_T1*(rhi-r(i))/(rhi-rlo) + &
         BC_T2*(r(i)-rlo)/(rhi-rlo) + &
         100.0d0*sin(z(j))*(r(i)-rlo)*(rhi-r(i))
 enddo
enddo
                                                                
!do j=0,M                                                             
! u(0,j)=T1                                                        
! u(N,j)=3.0d0                                                         
!enddo


end subroutine set_polar_2d

!----------------------------------------------------
subroutine polar_2d_heat(sdim,N,M,kappa,tau, r,z,dr,dz,u)
! polar coordinate 2-d heat equation solver.
! forward in time, central FD discretization
implicit none

integer,intent(in)          :: sdim                 
integer,intent(in)          :: N      ! discretization in r direction
integer,intent(in)          :: M       ! discretization in theta direction

real(kind=8),intent(in)     :: r(0:N)
real(kind=8),intent(in)     :: z(0:M)
real(kind=8)                :: u(0:N,0:M),u_new(0:N,0:M)
real(kind=8)                :: dr,dz
real(kind=8),intent(in)     :: tau
real(kind=8),intent(in)     :: kappa
real(kind=8),external       :: f_src


integer                     :: i,j


 do i=1,N-1
  do j=1,M-1
   u_new(i,j)= (1.0d0+kappa*tau*(-2.0d0/(dr**2.0d0)- &
                       2.0d0/((r(i)**2.0d0)*(dz**2.0d0))))*u(i,j) &
                + kappa*tau*(  & 
                 (1.0d0/(dr**2.0d0)+1.0d0/(r(i)*2.0d0*dr))*u(i+1,j) &
                + (1.0d0/(dr**2.0d0)-1.0d0/(r(i)*2.0d0*dr))*u(i-1,j) &
                + 1.0d0/((r(i)**2.0d0)*(dz**2.0d0))*u(i,j+1) &
                + 1.0d0/((r(i)**2.0d0)*(dz**2.0d0))*u(i,j-1)) &
               + tau*f_src(r(i),z(j)) 
  enddo
 enddo
  
 do i=1,N-1
   u_new(i,0)= (1.0d0+kappa*tau*(-2.0d0/(dr**2.0d0)- &
                       2.0d0/((r(i)**2.0d0)*(dz**2.0d0))))*u(i,0) &
                + kappa*tau*( &
                + (1.0d0/(dr**2.0d0)+1.0d0/(r(i)*2.0d0*dr))*u(i+1,0) &
                + (1.0d0/(dr**2.0d0)-1.0d0/(r(i)*2.0d0*dr))*u(i-1,0) &
                + 1.0d0/((r(i)**2.0d0)*(dz**2.0d0))*u(i,1) &
                + 1.0d0/((r(i)**2.0d0)*(dz**2.0d0))*u(i,M-1)) &
               + tau*f_src(r(i),z(0)) 
 enddo

 do i=1,N-1
  u_new(i,M)=u_new(i,0)
 enddo
 do j=0,M
  u_new(0,j)=BC_T1
  u_new(N,j)=BC_T2
 enddo

! do j=0,M
!  write(91,*) u_new(0:N,j) 
! enddo


 do i=0,N
  do j=0,M
   u(i,j)=u_new(i,j)
  enddo
 enddo

!enddo



end subroutine polar_2d_heat


end module mmat_FVM



!-------------------------------------------------------

FUNCTION NORM_2d(X1,Y1,X2,Y2)
IMPLICIT NONE

REAL(KIND=8), INTENT(IN) :: X1,Y1,X2,Y2
REAL(KIND=8)             :: NORM_2d


NORM_2d = sqrt((X1-X2)**2.0d0 + (Y1-Y2)**2.0d0)

END FUNCTION norm_2d

!--------------------------------------------
!            function u,v
!---------------------------------------------
FUNCTION U(t,x,y)
USE GeneralClass
IMPLICIT NONE

REAL(KIND=8),INTENT(IN):: t,x,y
REAL(KIND=8)           :: U
REAL(KIND=8)           :: nm
!PI=4.0d0*atan(1.0d0)

!U=100.0d0*y*(1.0d0-y)*(pi/2.0d0-atan(x))
!U=2.0d0
!u=1.0d0
!if (nm .le. 0.25d0) THEN

!u= 2.0d0*pi*(2.0d0*y/100.0d0-1.0d0)*(1.0d0-((2*x/100.0d0-1.0d0)**2))
u = 2.0d0*y -100.0d0

!   u = -2.0d0*(y-0.5d0)	
!else
!  u = 0.0d0
!endif

!u=-2.0d0*pi*y*(1.0d0-x**2.0d0)

END FUNCTION 
!------------------------
FUNCTION V(t,x,y)
USE GeneralClass
IMPLICIT NONE

REAL(KIND=8),INTENT(IN):: t,x,y
REAL(KIND=8)           :: V
real(kind=8)           :: nm

!PI=4.0d0*atan(1.0d0)

!V=atan(0.1d0*x*(1.0d0-x)*y*(1.0d0-y)*(1.0d0+t))
!V=1.0d0
!v=-3.0d0

!v= 2.0d0*(x-0.5d0)*exp(-((y-0.5d0)**2+(x-0.5d0)**2)/100.0d0)


!v = 2.0d0*pi*(2*x/100.0d0-1.0d0)*(((2.0d0*y/100.0d0-1.0d0)**2)-1.0d0)

v = 100.0d0 - 2.0d0*x

!if(nm .le. 0.25d0)then
!   v= 2.0d0*(x-0.5d0)
!else
!   v=0
!endif

!v=-2.0d0*pi*x*(y**2.0d0-1.0d0)

END FUNCTION v
!------------------------------------- test_flag==2


!--------------------------------
function f_src(x,y)
! source term for subroutine polar_2d_heat
implicit none

real(kind=8)  :: x,y
real(kind=8)  :: f_src

f_src = 0.0d0


return
end function f_src

