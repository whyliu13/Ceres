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

contains

!---------------------------------------------------------
subroutine dist_concentric(imat,x,y,dist,probtype_in)
implicit none

integer,intent(in)               :: imat,probtype_in
real(kind=8)                     :: x,y,dist
real(kind=8)                     :: dist1,dist2,dist3,dist4
integer                          :: i
real(kind=8)                     :: r1,r2
real(kind=8)                     :: center(2)

real(kind=8)                     :: x0,y0
real(kind=8)                     :: c1,c2,tt   !theta

if(probtype_in .eq. 1) then
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

elseif(probtype_in .eq. 3)then
 ! convert from (-1,1)  to (0,1)
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

! if(abs(x0-c1) .lt. 1.0e-10 .and. abs(y0-c2) .lt. 1.0e-10)then 
!     beta = 0.0d0;
! end 
  dist1 = sqrt((x0-c1)**2.0d0 + (y0-c2)**2.0d0) - &
         (0.5d0 + 0.2d0*sin(5.0d0*tt))
  dist2 = sqrt((x0-c1)**2.0d0 + (y0-c2)**2.0d0) - &
         (0.5d0 -radeps + 0.2d0*sin(5.0d0*tt))
 if(imat .eq. 3) then
   dist = dist1
 elseif(imat .eq. 2) then
   if(dist1 .gt. 0.0d0 .and. dist2 .gt. 0.0d0) then
    dist = -dist1
   elseif(dist1 .lt. 0.0d0 .and. dist2 .gt. 0.0d0) then
    dist = min(abs(dist1),abs(dist2))
   else
    dist = dist2
   endif
 elseif(imat .eq. 1)then
   dist = -dist2
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

else
 print *,"probtype_in invalid"
 stop
endif


end subroutine dist_concentric
!----------------------------------------------------------


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
subroutine check_sign(G,flag)
implicit none

real(kind=8)         :: g(4)
integer              :: flag
integer              :: i,j

flag = 0
do i = 1,4
  if(abs(g(i)) .lt. eps) then
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

end subroutine check_sign
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
   call dist_concentric(im,center(1),center(2),dist,probtype_in)
   if(abs(dist) .gt. h) then
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
real(kind=8)             :: d1,d2,d3
real(kind=8)             :: d(3)
real(kind=8)             :: v(3,2)

real(kind=8)            :: ratio,ratio1,ratio2
real(kind=8)            :: x1(2),x2(2),xx(2,2)

real(kind=8)            :: area, cen(2)
real(kind=8)            :: area1,area2
real(kind=8)            :: cen1(2),cen2(2)
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
 call dist_concentric(im,v1(1),v1(2),d1,probtype_in)      ! prob_type = 3
 call dist_concentric(im,v2(1),v2(2),d2,probtype_in)
 call dist_concentric(im,v3(1),v3(2),d3,probtype_in)

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
   stop 
 elseif(d1 .ge. 0.0d0  .and. d2 .ge. 0.0d0 .and. d3 .ge. 0.0d0)then         ! +0  +0  +0
  call tri_area8(v1,v2,v3,area)
  call tri_centroid(v1,v2,v3,cen)
  print *,"case 1", "area",area
 elseif(d1 .le. 0.0d0  .and. d2 .le. 0.0d0 .and. d3 .le. 0.0d0)then         ! -0  -0  -0   
  print *,"case 2"
    ! do nothing
 elseif(abs(d1) .lt. 1.0e-10 .and. d2*d3 .lt. 0.0d0)then
  print *,"case 3"
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
 elseif(abs(d2) .lt. 1.0e-10 .and. d1*d3 .lt. 0.0d0)then
   print *, "case4"
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
 elseif(abs(d3) .lt. 1.0e-10 .and. d1*d2 .lt. 0.0d0)then
   print *,"case5"
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
 else
  print *,"case6"
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
    if(d(1)*d(3) .gt. 0.0d0)then
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
    print *,"area2",area2
      area = area2-area1
      do i=1,2
       cen(i)=(cen2(i)*area2 - cen1(i)*area1)/area 
      enddo
    else
     print *,"d(1) can be zero,  565"
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
    print *,"area2",area2
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
    print *,"area2",area2
      area = area2-area1
      do i=1,2
       cen(i)=(cen2(i)*area2 - cen1(i)*area1)/area
      enddo
    else
     print *,"d(1) can be zero,  565"
     stop
    endif
  else
    print *,"err, check, 611"
    stop
  endif
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

 AREA = 0.5* ABS(V1(1)*(V2(2)-V3(2)) &
           & - V2(1)*(V1(2)-V3(2)) &
           & + V3(1)*(V1(2)-V2(2)))


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


subroutine AdaptQuad_2d(iin,jin,nmat_in,dx,center,centroid,vf,probtype_in)
implicit none

integer,intent(in)            :: nmat_in,probtype_in
real(kind=8),intent(in)       :: dx(2)
real(kind=8),intent(in)       :: center(2)
integer     ,intent(in)       :: iin,jin

integer,parameter             :: ref_num = 5
integer                       :: i,im,i1,i2,i3,i4,i5,dir
integer                       :: ii,jj
integer                       :: flag,sign_flag
real(kind=8)                  :: dist,h2,h3,h4,h5,h6
real(kind=8)                  :: x,y
real(kind=8)                  :: h
real(kind=8)                  :: center_hold(2)
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
real(kind=8)                  :: g_center
real(kind=8)                  :: vf_check

type(points)                  :: seg(2)

real(kind=8)                  :: v1(2),v2(2),v3(2)


x = center(1)
y = center(2)
h = dx(1)

if (h.le.0.0) then
 print *,"h invalid"
 stop
endif


 cell%nodesnum = 4 
 allocate(cell%nodes(4))
 vf = 0.0d0
 centroid = 0.0d0
 vol = 0.0d0
 cxtemp = 0.0d0
 cytemp = 0.0d0
 
do im = 1,nmat_in
   call dist_concentric(im,x,y,dist,probtype_in)

   if(abs(dist) .gt. h) then
     flag = 0
     call update_volncen(nmat_in,dist,im,center,h,vol,cxtemp,cytemp)
   else
     flag = 1   
     call AdaptQuad_sub(center, h, centers1)
   endif 
 
   if (flag.eq.0) then
    ! do nothing  
   else if(flag .eq. 1) then
    do i1 = 1,4 
     do dir=1,2
      center_hold(dir)=centers1(i1,dir)
     enddo
     h2=h/2.0d0
     call AdaptQuad_sort(im,center_hold,h2,centers2,dist,flag,probtype_in)
     if(flag .eq. 1) then
      do i2 = 1,4
       do dir=1,2
        center_hold(dir)=centers2(i2,dir)
       enddo
       h3=h/4.0d0
       call AdaptQuad_sort(im,center_hold,h3,centers3,dist,flag, &
         probtype_in)
       if(flag .eq. 1) then
        do i3 = 1,4
         do dir=1,2
          center_hold(dir)=centers3(i3,dir)
         enddo
         h4=h/8.0d0
         call AdaptQuad_sort(im,center_hold,h4,centers4,dist, &
              flag,probtype_in)
         if(flag .eq. 1) then
          do i4 = 1,4
           do dir=1,2
            center_hold(dir)=centers4(i4,dir)
           enddo
           h5=h/16.0d0
           call AdaptQuad_sort(im,center_hold,h5, centers5, &
             dist,flag,probtype_in)
           if(flag .eq. 1)then
            do i5 = 1,4
             do dir=1,2
              center_hold(dir)=centers5(i5,dir)
             enddo
             h6=h/32.0d0
             call find_vertex(center_hold, h6, vertex)

             do ii = 1,4
              call dist_concentric(im,vertex(ii,1),vertex(ii,2), &
                G(ii),probtype_in)
             enddo  ! ii
             call check_sign(G,sign_flag)

             if (1.eq.0) then
              print *,"after check_sign  sign_flag=",sign_flag
              do ii=1,4
               print *,"ii,vertex,G ",ii,vertex(ii,1),vertex(ii,2),G(ii)
              enddo
             endif

             if(sign_flag .eq. 1) then
!              if(probtype_in .eq. 0 .or. &
!                 probtype_in .eq. 1 .or. &
!                 probtype_in .eq. 2) then
              if(1 .eq. 0)then
              call slopecal(im,center_hold,mx,my,probtype_in)
               !============================================================        
               ! call LS_LSF(h/32.0d0,centers5(i5,:),G,mx,my,alphadump)
               ! ax + by +alpha = 0
              call dist_concentric(im,center_hold(1),center_hold(2), &
                 G_center,probtype_in)

  ! phi=mx x + my y + c
  ! mx x' + my y' + c = phi'
  ! c=phi'-mx x' -my y'
  !              normalize 
  !              call normalize_cutline(mx,my,alpha)
  !              c = -mx*center_hold(1)-my*center_hold(2)+alpha
              c = -mx*center_hold(1)-my*center_hold(2)+ G_center

              call form_cell(vertex,cell)
              call cell_line_inter(cell,mx,my,c,seg)

              if (1.eq.0) then
               print *,"x,y,h,im,mx,my,c ",x,y,h,im,mx,my,c
              endif

              CALL VOL_FRAC_CAL(cell,mx,my,c,POS_PLG,NEG_PLG)
                       
               ! call outputplg(pos_plg)
  
              CALL POLY_CENTROID(NEG_PLG,cen_temp)
              CALL poly_area(NEG_PLG,VOL_temp) 
              CALL PLGDEL(POS_PLG)
              CALL PLGDEL(NEG_PLG)

              vol(im) = vol(im) + vol_temp
              cxtemp(im) = cxtemp(im) + vol_temp*cen_temp%val(1) 
              cytemp(im) = cytemp(im) + vol_temp*cen_temp%val(2)

              endif

            if(1 .eq. 1)then
              ! triangulation
              v1 = vertex(1,:)
              v2 = vertex(2,:)
              v3 = vertex(3,:)
              call triangle_interface_detect(im,probtype_in,v1,v2,v3,vol_temp,cen_temp%val)

              ! update centroid and volume
              vol(im) = vol(im) + vol_temp
              cxtemp(im) = cxtemp(im) + vol_temp*cen_temp%val(1) 
              cytemp(im) = cytemp(im) + vol_temp*cen_temp%val(2)
            !------------------------------------------------------------------------------
              v1 = vertex(2,:)
              v2 = vertex(3,:)
              v3 = vertex(4,:)
              call triangle_interface_detect(im,probtype_in,v1,v2,v3,vol_temp,cen_temp%val)

              ! update centroid and volume
              vol(im) = vol(im) + vol_temp
              cxtemp(im) = cxtemp(im) + vol_temp*cen_temp%val(1) 
              cytemp(im) = cytemp(im) + vol_temp*cen_temp%val(2)
            endif

              if (1.eq.0) then
               print *,"im,mx,my,c ",im,mx,my,c
               do ii=1,4
                print *,"ii,vertex,G ",ii,vertex(ii,1),vertex(ii,2),G(ii)
               enddo
               print *,"vfrac ",vol_temp/(h6*h6)
               print *,"cenx,ceny ",cen_temp%val(1),cen_temp%val(2)
              endif

             elseif(sign_flag .eq. 0) then 
              CALL dist_concentric(im,center_hold(1),center_hold(2), &
                 dist,probtype_in)
              CALL update_volncen(nmat_in,dist,im,center_hold, &
                 h/32.0d0,vol,cxtemp,cytemp)
             else
              print *,"sign_flag = ", sign_flag
              stop
             endif
            enddo  ! i5
           elseif(flag .eq. 0) then
            call update_volncen(nmat_in,dist,im,center_hold, &
              h/16.0d0,vol,cxtemp,cytemp)  
           else
            print *,"flag error in AdaptQuad_2d"
            stop
           endif
          enddo ! i4
         elseif(flag .eq. 0) then
          call update_volncen(nmat_in,dist,im,center_hold,h/8.0d0, &
            vol,cxtemp,cytemp)  
         else
          print *,"flag error in AdaptQuad_2d"
          stop
         endif
        enddo  ! i3
       elseif(flag .eq. 0) then
        call update_volncen(nmat_in,dist,im,center_hold,h/4.0d0, &
          vol,cxtemp,cytemp)  
       else
        print *,"flag error in AdaptQuad_2d" 
        stop
       endif
      enddo ! i2
     elseif(flag .eq. 0) then
      call update_volncen(nmat_in,dist,im,center_hold,h/2.0d0,vol, &
         cxtemp,cytemp)  
     else
       print *,"flag error in AdaptQuad_2d"    
     endif
    enddo ! i1

   else 
    print *,"flag invalid"
    stop
   endif

enddo ! im

if (h.le.0.0) then
 print *,"h invalid"
 stop
endif


do im = 1,nmat_in
   if(vol(im) .gt. eps)then
    centroid(im,1) = cxtemp(im)/vol(im)
    centroid(im,2) = cytemp(im)/vol(im)
   else
     vol(im) = 0.0d0
     centroid(im,:) = 0.0d0
   endif
   vf(im) = vol(im)/(h*h)

enddo


call vf_correct(iin,jin,nmat_in,vf,probtype_in)

if (probtype_in.eq.0) then
 ! do nothing
else if (probtype_in.eq.2) then
 ! do nothing
else if (probtype_in.eq.1 .or. probtype_in .eq. 3) then
 vf(2) = 1.0d0 - vf(1) -vf(3)

 if(abs(vf(2)) .lt. eps)then
  vf(2) = 0.0d0
 endif

 if(vf(2) .lt. -eps) then
  print *,"vf(2) is negative"
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
else if (probtype_in.eq.1 .or. probtype_in .eq. 3) then
 if(vf(2) .gt. eps)then
  do dir=1,2
   centroid(2,dir) =  &
    (center(dir) - vf(1)*centroid(1,dir) - vf(3)*centroid(3,dir))/vf(2)
  enddo
 endif
else
 print *,"probtype_in invalid"
 stop
endif


call renormalize_vf(nmat_in,vf) 


end subroutine AdaptQuad_2d
! -----------------------------------------------
!---------------------------------------------------------------------
subroutine vf_correct(iin,jin,nmat_in,vf,probtype_in)
implicit none

integer       ,intent(in)  :: nmat_in,iin,jin,probtype_in
real(kind=8)               :: vf(nmat_in)
integer                    :: i,j
real(kind=8)               :: vcheck

if ((probtype_in.eq.0).or. &
    (probtype_in.eq.2)) then

 if (nmat_in.ne.2) then
  print *,"nmat_in invalid"
  stop
 endif
 
 vcheck = vf(1) + vf(2)
 if(vf(1) .lt. 0.0d0 .or. vf(2) .lt. 0.0d0) then
  print *,"vf is negative"
 endif

 if(vcheck .gt. 1.0d0+eps) then
  print *,"goes into vf_correct",vcheck
  stop
 endif

else if (probtype_in.eq.1 .or. probtype_in .eq. 3) then
 if (nmat_in.ne.3) then
  print *,"nmat_in invalid"
  stop
 endif
 
 vcheck = vf(1) + vf(3)
 if(vf(1) .lt. 0.0d0 .or. vf(3) .lt. 0.0d0) then
  print *,"vf is negative"
 endif

 if(vcheck .gt. 1.0d0+eps) then
  print *,"goes into vf_correct",vcheck
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
    call AdaptQuad_2d(i,j,nmat_in,dx,center,cen_temp,vf_temp,probtype_in)

    do im = 1, nmat_in
     do dir=1,2
      centroid(i,j,im)%val(dir) = cen_temp(im,dir)
     enddo
     vf(i,j,im) = vf_temp(im)

     if (1.eq.0) then
      print *,"i,j,im,vf,cenx,ceny ",i,j,im,vf_temp(im), &
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
    print *,"delx or dely invalid"
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
    print *,"delx or dely invalid"
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



 else
  print *,"probtype_in invalid"
  stop
 endif

return
end subroutine get_filament_source

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
Function  exact_temperature(x,y,t,im,probtype_in,nmat_in,alpha,dclt_flag)
Use generalclass
implicit none

integer,intent(in)         :: im,probtype_in,nmat_in,dclt_flag
real(kind=8),intent(in)    :: alpha(nmat_in)
real(kind=8),intent(in)    :: t,x,y
real(kind=8)               :: exact_temperature
real(kind=8)              :: radius,theta,r1,r2
real(kind=8)              :: TLO,THI,yI,yHI,a1,b1,a2,b2
real(kind=8)              :: mypi,delx,dely

 mypi=4.0d0*atan(1.0d0)
 if (probtype_in.eq.1 .or. probtype_in .eq. 3) then
  if (nmat_in.ne.3) then
   print *,"nmat_in invalid"
   stop
  endif
  if (im.eq.2) then
   r1=radcen-radeps
   r2=radcen+radeps
      
   delx=x-0.5d0
   dely=y-0.5d0
   radius = sqrt(delx**2.0d0 +dely**2.0d0)
 
       ! x=r cos(theta)
       ! y=r sin(theta)
   if (radius.le.radeps/1000.0) then
    theta=0.0
   else if ((delx.ge.0.0).and.(dely.ge.0.0)) then
    theta=acos(delx/radius)
   else if ((delx.le.0.0).and.(dely.ge.0.0)) then
    theta=acos(abs(delx)/radius)
    theta=mypi-theta
   else if ((delx.le.0.0).and.(dely.le.0.0)) then
    theta=acos(abs(delx)/radius)
    theta=mypi+theta
   else if ((delx.ge.0.0).and.(dely.le.0.0)) then
    theta=acos(delx/radius)
    theta=2.0d0*mypi-theta
   else
    print *,"delx or dely invalid"
    stop
   endif

   exact_temperature=2.0d0+sin(theta)*exp(-t/(radcen**2))
  else if ((im.eq.1).or.(im.eq.3)) then
   exact_temperature=0.0
  else
   print *,"im invalid 6"
   stop
  endif

 else if ((probtype_in.eq.0).or.(probtype_in.eq.2)) then
  if (nmat_in.ne.2) then
   print *,"nmat_in invalid"
   stop
  endif
     ! T1=a1+b1 y
     ! T2=a2+b2 (y-yI)
     ! b1 * k1 = b2 * k2
     ! a1=TLO
     ! a2=a1+b1 * yI
     ! a2=TLO+b1 * yI
     ! a2+b2(yHI-yI)=THI
     ! a2+(b1 k1/k2)(yHI-yI)=THI
     ! TLO+b1 yI+(b1 k1/k2)(yHI-yI)=THI
     ! TLO k2 + b1 (yI k2 + k1 (yHI - yI))=THI k2
     ! b1(yI k2 + k1 (yHI-yI))=(THI-TLO)k2
!
!  if (dclt_flag.eq.0) then
   TLO=3.0d0
   THI=2.0d0
   yI=0.3d0
   yHI=1.0
   a1=TLO
   b1=(THI-TLO)*alpha(2)/(yI * alpha(2) +alpha(1)*(yHI-yI))
   b2=b1*alpha(1)/alpha(2)
   a2=a1+b1*yI 
!  else if (dclt_flag.eq.1) then
!   TLO=3.0d0
!   THI=2.0d0
!   a1=TLO
!   b1=(THI-TLO)/0.3d0
!   yI = 0.0d0           ! NULL the yI
!   a2=0.0d0
!   b2=0.0d0
!  else
!   print *,"dclt_flag invalid"
!   stop
!  endif

  if (probtype_in.eq.0) then
!   if(dclt_flag .eq. 0)then
!    if (im.eq.1) then
!     exact_temperature=a1+b1*y
!    else if (im.eq.2) then
!     exact_temperature=a2+b2*(y-yI)
!    else
!     print *,"im invalid 7"
!     stop
!    endif
!   elseif(dclt_flag .eq. 1)then
    if (im.eq.1) then
     exact_temperature=a1+b1*y
    else if (im.eq.2) then
     exact_temperature=a2+b2*(y-yI)
!     exact_temperature = a2
    else
     print *,"im invalid 7"
     stop
    endif
!   else
!    print *,"invalid dclt_flag in exact_temperature"
!    stop
!   endif

  else if (probtype_in.eq.2) then

   if (im.eq.1) then
    exact_temperature=a1+b1*x
   else if (im.eq.2) then
    exact_temperature=a2+b2*(x-yI)
   else
    print *,"im invalid 8"
    stop
   endif
  else
   print *,"probtype_in invalid1"
   stop
  endif

 else
  print *,"probtype_in invalid2",probtype_in
  stop
 endif

end function exact_temperature

