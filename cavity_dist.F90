
module test




contains

! From paper:  
!  Mu YT, Chen L, He YL, Kang QJ, Tao WQ. 
! Nucleate boiling performance evaluation of cavities at mesoscale level. 
! International Journal of Heat and Mass Transfer. 2017 Mar 31;106:708-19.
! 
! 
! 



! length unit ( mm )
! 
!  MATERIAL NUMBER
!  1:liquid   2: vapor   3: solid 

!  COORDIANTE TYPE
!  1 3D cartisian
!  2 R-Z   

!  SHAPE 
!  1--> rectangular
!  2--> trapzoidal
!  3--> vshape
!  4--> semicircle
!  5--> spherical reentrant



subroutine cavity_distf_13(cavity_type, coord_type, sdim, x_in, dist)
implicit none
! liquid +  solid - 

integer,intent(in)      :: coord_type
integer, intent(in)     :: cavity_type
integer,intent(in)      :: sdim 
real(kind=8),intent(in) :: x_in(sdim)
real(kind=8)            :: x(sdim)

integer                 :: dist_sign
real(kind=8)            :: dist
real(kind=8)            :: dist1,dist2,dist3
real(kind=8)            :: x1(sdim),x2(sdim),x3(sdim),x4(sdim)
real(kind=8)            :: trap,vshape,circle,radius,lcinter

dist_sign = 1.0d0
dist = 0.0d0
dist1 = 0.0d0
dist2 = 0.0d0
dist3 = 0.0d0
x1 = 0.0d0
x2 = 0.0d0
x3 = 0.0d0
x4 = 0.0d0
 
! convert from cm to mm 
do i = 1,sdim
  x(i) = x_in(i)*10.0d0
enddo


!  dimension sanity check
if(coord_type .eq. 1) then   ! 3D cartisian
 if(sdim .ne. 3)then
  print *, "coord_type is not consistent with dimension"
  stop
 endif
elseif(coord_type .eq. 2) then ! R-Z axis symmetric
 if(sdim .ne. 2) then
  print *,"coord_type is not consistent with dimension"
  stop
 endif
else
 print *,"invalid coord_type"
 stop
endif



if(cavity_type .eq. 1)then       ! 1--> rectangular
!  *                 * 
!  *       *********** 
!  *       *         * 
!  *       *         * 
!  *********         * 
!  
!*******************************************************************
if(coord_type .eq. 2)then   !  R-Z  axis symmetric
 if(x(1) .ge. 0.25d0) then   ! r >= 0.25
  dist1 = x(2) - 0.3d0        ! dist = z - 0.3   
  if(x(2) .le. 0.05d0)then
   dist2 = sqrt( (x(1)- 0.25d0)**2.0d0 + (x(2) - 0.05d0)**2.0d0)
   dist = -1.0d0*min(abs(dist1), abs(dist2))
  elseif(x(2) .lt. 0.3d0) then
   dist2 = 0.25d0 - x(1)
   dist = -1.0d0*min(abs(dist1),abs(dist2))
  else
   dist = dist1
  endif
 else  !   0 < r < 0.25
  dist1 = x(2) - 0.05
  if(x(2) .gt. 0.3d0)then
   dist2 = sqrt( (x(1)- 0.25d0)**2.0d0 + (x(2) - 0.3d0)**2.0d0)
   dist = min(dist1, dist2)
  elseif(x(2) .gt. 0.05d0)then
   dist2 = 0.25 - x(1) 
   dist = min(dist1,dist2)
  else
   dist = dist1
  endif
 endif
endif


if(coord_type .eq. 1)then   ! 3D cartisian
 if(x(1) .ge. 0.25d0 .or. x(2) .ge. 0.25d0) then
  dist1 = x(3) - 0.3d0 
  if(x(3) .gt. 0.3d0)then         ! x(3) > 0.3
   dist = dist1
  elseif(x(3) .lt. 0.05d0) then   !    x(3) < 0.05
    if(x(1) .ge. 0.25d0 .and. x(2) .le. 0.25d0)then
    x1(1) = 0.25d0
    x1(2) = 0.0d0
    x1(3) = 0.05d0
    x2(1) = 0.25d0
    x2(2) = 0.25d0
    x2(3) = 0.05d0
    call dist_point_to_line(sdim,x1,x2,x,dist2)
    dist = -1.0d0*min(dist1,dist2)
   elseif(x(1) .le. 0.25d0 .and. x(2) .ge. 0.25d0)then
    x1(1) = 0.0d0
    x1(2) = 0.25d0
    x1(3) = 0.05d0
    x2(1) = 0.25d0
    x2(2) = 0.25d0
    x2(3) = 0.05d0
    call dist_point_to_line(sdim,x1,x2,x,dist2)
    dist = -1.0d0*min(dist1,dist2)
   else
    x1(1) = 0.25d0
    x1(2) = 0.25d0
    x1(3) = 0.05d0
    call l2norm(sdim,x,x1,dist2)
    dist = -1.0d0*min(dist1,dist2)    
   endif
  else   !  0.05 <  x(3)  < 0.3
   if(x(1) .ge. 0.25d0 .and. x(2) .le. 0.25d0)then
    dist2 = 0.25d0 - x(1)
    dist = -1.0d0*min(dist1,dist2)
   elseif(x(1) .le. 0.25d0 .and. x(2) .ge. 0.25d0)then
    dist2 = 0.25d0 - x(2) 
    dist = -1.0d0*min(dist1,dist2)
   else
    x1(1) = 0.25d0
    x1(2) = 0.25d0
    x1(3) = 0.3d0
    x2(1) = 0.25d0
    x2(2) = 0.25d0
    x2(3) = 0.05d0
    call dist_point_to_line(sdim,x1,x2,x,dist2)
    dist = -1.0d0*min(dist1,dist2)    
   endif   
  endif
 else    ! x(1) .lt. 0.25d0 .and. x(2) .lt. 0.25d0
  dist1 = x(3) - 0.05d0
  if(x(3) .lt. 0.05d0 )then
   dist = dist1
  elseif(x(3) .gt. 0.3d0)then
    x1(1) = 0.25d0
    x1(2) = 0.0d0
    x1(3) = 0.3d0
    x2(1) = 0.25d0
    x2(2) = 0.25d0
    x2(3) = 0.3d0
    call dist_point_to_line(sdim,x1,x2,x,dist2) 
    x3(1) = 0.0d0
    x3(2) = 0.25d0
    x3(3) = 0.3d0
    x4(1) = 0.25d0
    x4(2) = 0.25d0
    x4(3) = 0.3d0
    call dist_point_to_line(sdim,x3,x4,x,dist3) 
    if(dist2 .lt. 0.0d0 .or. dist3 .lt. 0.0d0) then
     print *,"dist2 and dist3 should be positive"
     stop 
    endif  
    dist = min(abs(dist1),abs(dist2),abs(dist3))
  else  !   0.05 <= x(3) <= 0.3
    dist2 = 0.25d0 - x(1)   
    dist3 = 0.25d0 - x(2)
    if(dist2 .lt. 0.0d0 .or. dist3 .lt. 0.0d0) then
     print *,"dist2 and dist3 should be positive"
     stop 
    endif
    dist = min(abs(dist1),abs(dist2),abs(dist3))  
  endif
 endif
endif   
!---------------------------------------------------------------



elseif(cavity_type .eq. 2)then    ! 2--> trapzoidal
!**************************************************************
!             *                 *  
!             *         x3      x4
!             *         *********
!             *       *
!             *     *
!             *****
!            x1   x2
!**************************************************************
if(coord_type .eq. 2)then   !  R-Z  axis symmetric
 ! determine sign
 if(x(1) .ge. 0.25d0)then
  if(x(2) .ge. 0.3d0)then
   dist_sign = 1.0d0
  else
   dist_sign = -1.0d0
  endif 
 elseif(x(1) .le. 0.125d0)then
  if(x(2) .ge. 0.05d0)then
   dist_sign = 1.0d0
  else
   dist_sign = -1.0d0
  endif
 else
   trap = 2.0d0*(x(1)-0.125d0) + 0.05d0 - x(2)
   if(trap .gt. 0.0d0)then
    dist_sign = -1.0d0
   else
    dist_sign = +1.0d0
   endif
 endif
 ! determine distance
  x1(1) = 0.0d0
  x1(2) = 0.05d0 
  x2(1) = 0.125d0
  x2(2) = 0.05d0
  x3(1) = 0.25d0
  x3(2) = 0.3d0
  x4(1) = 0.6d0
  x4(2) = 0.3d0
  call dist_point_to_line(sdim,x1,x2,x,dist1)
  call dist_point_to_line(sdim,x2,x3,x,dist2) 
  call dist_point_to_line(sdim,x3,x4,x,dist3)  
  dist = dist_sign * min(dist1,dist2,dist3)

endif

if(coord_type .eq. 1)then   !  3D cartisian
 ! determine sign

 ! determin distance

endif
!------------------------------------------------------------

!  
! 
elseif(cavity_type .eq. 3) then   ! 3--> vshape
!**************************************************************
!             *            *  
!             *    x3      x4 
!             *    ********* 
!             *   *
!             * *
!             *
!            x1    
!**************************************************************
if(coord_type .eq. 2)then   !  R-Z  axis symmetric
 ! determine sign
 if(x(1) .ge. 0.25d0)then
  if(x(2) .ge. 0.3d0)then
   dist_sign = 1.0d0
  else
   dist_sign = -1.0d0
  endif 
 else
  vshape = x(1)-x(2)+0.05d0
  if(vshape .gt. 0.0d0)then
   dist_sign = -1.0d0
  else
   dist_sign = 1.0d0
  endif
 endif
 ! determine distance
  x1(1) = 0.0d0
  x1(2) = 0.05d0 
  x3(1) = 0.25d0
  x3(2) = 0.3d0
  x4(1) = 0.6d0
  x4(2) = 0.3d0
  call dist_point_to_line(sdim,x1,x3,x,dist1) 
  call dist_point_to_line(sdim,x3,x4,x,dist3)  
  dist = dist_sign * min(dist1,dist3)

endif

if(coord_type .eq. 1)then   !  3D cartisian
 ! determine sign


endif

!------------------------------------------------------------

!  
! 
elseif(cavity_type .eq. 4) then   ! 4--> semicircle
!**************************************************************
!             *                 *  
!             *         x1      x2
!             *         ********* 
!             *       *
!             *     *
!             *   *
!             *   
!**************************************************************
if(coord_type .eq. 2)then   !  R-Z  axis symmetric
 ! determine sign
 if(x(2) .ge. 0.3d0)then
   dist_sign = +1.0d0
 else
   radius = 0.25d0 - sqrt((x(1) - 0.0d0)**2.0d0 + (x(2) - 0.3d0)**2.0d0)
   if(radius .ge. 0.0d0)then
    dist_sign = +1.0d0
   else
    dist_sign = -1.0d0
   endif  
 endif

 ! determine distance
  x1(1) = 0.25d0
  x1(2) = 0.3d0 
  x2(1) = 0.6d0
  x2(2) = 0.3d0
  if(x(2) .gt. 0.3d0)then
   call dist_point_to_line(sdim,x1,x2,x,dist1)
   dist = dist_sign*dist1
  else
   dist2 = 0.25d0 - sqrt((x(1))**2.0d0 + (x(2) - 0.3d0)**2.0d0)   
   call dist_point_to_line(sdim,x1,x2,x,dist3)  
   dist = dist_sign * min(abs(dist2),abs(dist3))
  endif
endif

if(coord_type .eq. 1)then   !  3D cartisian
 ! determine sign


endif


!------------------------------------------------------------

!  
! 
elseif(cavity_type .eq. 5) then   ! 5--> spherical reentrant
!**************************************************************
!             *  x1                          x2
!             *  *****************************
!             * * 
!             * * x3
!             *   * 
!             *     *
!             *       *
!             *        *         
!             *         *
!             *         * 
!             *       *
!             *     *
!             *   *
!             * *  
!**************************************************************
if(coord_type .eq. 2)then   !  R-Z  axis symmetric
 ! determine sign
  lcinter = sqrt(0.25d0**2.0d0 - 0.15d0**2.0d0) + 0.3
 if(x(2) .ge. 0.6d0)then
   dist_sign = +1.0d0
 elseif(x(2) .lt. lcinter)then
   radius = 0.25d0 - sqrt((x(1) - 0.0d0)**2.0d0 + (x(2) - 0.3d0)**2.0d0)
   if(radius .ge. 0.0d0)then
    dist_sign = +1.0d0
   else
    dist_sign = -1.0d0
   endif  
 else
   if(x(1) .le. 0.15d0)then
    dist_sign = +1.0d0
   else
    dist_sign = -1.0d0
   endif
 endif

 ! determine distance
  x1(1) = 0.15d0
  x1(2) = 0.6d0 
  x2(1) = 0.6d0
  x2(2) = 0.6d0
  x3(1) = 0.15d0
  x3(2) = lcinter

 if(x(2) .ge. 0.6d0)then
  call dist_point_to_line(sdim,x1,x2,x,dist1)
  dist = dist_sign*dist1 
 elseif(x(2) .lt. lcinter)then
   dist2 = 0.25d0 - sqrt((x(1) - 0.0d0)**2.0d0 + (x(2) - 0.3d0)**2.0d0)
   call dist_point_to_line(sdim,x1,x2,x,dist3)  
   dist = dist_sign * min(abs(dist2),abs(dist3)) 
 else
   call dist_point_to_line(sdim,x1,x3,x,dist1)
   dist = dist_sign*dist1

   if(x(1) .ge. 0.15d0)then
    dist2 = 0.25d0 - sqrt((x(1) - 0.0d0)**2.0d0 + (x(2) - 0.3d0)**2.0d0) 
    call dist_point_to_line(sdim,x1,x2,x,dist3)  
    dist = dist_sign * min(abs(dist2),abs(dist3),abs(dist1))
   endif 
 endif

endif

if(coord_type .eq. 1)then   !  3D cartisian
 ! determine sign


endif


else
 print *, "invalid cavity type flag"
 stop
endif



! from mm to cm
dist = dist/10.0d0



end subroutine cavity_distf_13




! --------------------------------------------------
! vapor +   liquid -

subroutine cavity_distf_12(coord_type, sdim, x_in, dist)
implicit none

integer,intent(in)      :: coord_type
integer,intent(in)      :: sdim 
real(kind=8),intent(in) :: x_in(sdim)
real(kind=8)            :: x(sdim)

real(kind=8)            :: dist
real(kind=8)            :: dist1,dist2,dist3
real(kind=8)            :: x1(sdim),x2(sdim),center(sdim)


! convert from cm to mm 
do i = 1,sdim
  x(i) = x_in(i)*10.0d0
enddo


if(coord_type .eq. 1) then   ! 3D cartisian
 if(sdim .ne. 3)then
  print *, "coord_type is not consistent with dimension"
  stop
 endif
elseif(coord_type .eq. 2) then ! R-Z axis symmetric
 if(sdim .ne. 2) then
  print *,"coord_type is not consistent with dimension"
  stop
 endif
else
 print *,"invalid coord_type"
 stop
endif


if(coord_type .eq. 2)then 
  center(1) = 0.0d0
  center(2) = 0.05d0
  call l2norm(sdim, center, x , dist1)
  dist = 0.1d0 - dist1 
endif

if(coord_type .eq. 1)then 
  center(1) = 0.0d0
  center(2) = 0.0d0
  center(3) = 0.05d0
  call l2norm(sdim, center, x , dist1)
  dist = 0.1d0 - dist1  
endif




! from mm to cm
dist = dist/10.0d0

end subroutine cavity_distf_12



!------------------------------------------------------------
!       p1-------------------------p2
!               *
!               *  
!               x
!----------------------------------------------------------
subroutine dist_point_to_line(sdim,p1,p2,x,dist)
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

call l2norm(sdim, p1, x, diff10)
call l2norm(sdim, p2, p1, diff21)

s = -1.0d0*(dot_product(x10,x21))/(diff21**2.0d0)

if(s .gt. 1.0d0)then
 call l2norm(sdim, p2, x,dist)
elseif(s .lt. 0.0d0)then
 call l2norm(sdim,p1,x,dist)
else
 dist = sqrt(((diff10**2.0d0 )*(diff21**2.0d0) - & 
        (dot_product(x10,x21))**2.0d0)/ &
        (diff21**2.0d0))
endif

deallocate(x10,x21)

end subroutine dist_point_to_line
!------------------------------------------------------------
subroutine l2norm(sdim, x1,x2, x1x2norm)
implicit none

integer, intent(in)      :: sdim
real(kind=8),intent(in)  :: x1(sdim),x2(sdim)
real(kind=8)             :: x1x2norm

integer                  :: i
real(kind=8),allocatable :: diff(:)

x1x2norm = 0.0d0
allocate(diff(sdim))
do i = 1,sdim
 diff(i) = x1(i)-x2(i)
enddo

do i = 1,sdim
 x1x2norm = x1x2norm + diff(i)**2.0d0
enddo

x1x2norm = sqrt(x1x2norm)

deallocate(diff)

end subroutine l2norm



end module test






program testtest
use test
implicit none

real(kind=8),parameter            :: a = 0.6d0
real(kind=8),parameter            :: b = 0.8d0
integer,parameter                 :: N = 50

real(kind=8),allocatable          :: x(:),y(:)
integer                           :: i,j
real(kind=8)                      :: dx(2), p(2)
real(kind=8),allocatable          :: dist(:,:)

dx(1) = a/N
dx(2) = b/N 

 open(unit = 2, file = "out")

allocate(x(N+1),y(N+1),dist(N+1,N+1))

do  i = 1,N+1
 x(i) = (i-1)*dx(1)
 y(i) = (i-1)*dx(2)
enddo




 do i = 1,N+1
  do j = 1,N+1
   p(1) = x(i)
   p(2) = y(j)  
   call cavity_distf_13(5, 2, 2, p, dist(i,j))
  enddo
 enddo

do i = 1,N+1
 write(2,*) dist(:,i) 
enddo

deallocate(x,y,dist)
 close(2)

end program

