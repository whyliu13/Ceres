





! the input  x, y must within [0,1]X[0,1]



subroutine dist_concentric(imat,x,y,dist,probtype_in)
implicit none

real*8,parameter                :: mypi = 4.0d0*atan(1.0d0)
integer,parameter               :: methodtype = 1
integer,intent(in)               :: imat,probtype_in
real(kind=8)                     :: x,y,dist
real(kind=8)                     :: dist1,dist2,dist3,dist4
integer                          :: i
real(kind=8)                     :: r1,r2
real(kind=8)                     :: center(2)

real(kind=8)                     :: x0,y0
real(kind=8)                     :: c1,c2,tt   !theta
real(kind=8),allocatable         :: xt(:),yt(:)
real(kind=8)                     :: signtemp




!   method 1 
if(methodtype .eq. 1)then
 allocate(xt(1000),yt(1000))
 call starshape(1000,xt,yt)
 dist1 = 0.0d0
 do i = 1,1000
  dist2 = sqrt((x-(xt(i)+1.0d0)/2.0d0)**2.0d0 + &
               (y-(yt(i)+1.0d0)/2.0d0)**2.0d0) 
  if(dist2 .gt. dist1)then
   dist1 = dist2
  endif
 enddo

 x0 = 2.0d0*x-1.0d0
 y0 = 2.0d0*y-1.0d0

 c1 = 0.02d0*sqrt(5.0d0)
 c2 = 0.02d0*sqrt(5.0d0)
 
 tt = atan((y0-c2)/(x0-c1));  
 if((x0-c1) .ge. 0.0d0 .and. (y0-c2) .ge. 0.0d0)then
    ! do nothing
 elseif((x0-c1) .le. 0.0d0 .and. (y0-c2) .gt. 0.0d0)then
    tt = tt + mypi;
 elseif((x0-c1) .lt. 0.0d0 .and. (y0-c2) .lt. 0.0d0)then
    tt = tt +mypi;
 else
    tt = 2.0d0*mypi + tt;
 endif

 dist3 =(x0-c1)**2.0d0 + (y0-c2)**2.0d0 - &
         (0.5d0 + 0.2d0*sin(5.0d0*tt))**2.0d0
 if(imat .eq. 1)then
  dist = -sign(dist1,dist3)
 elseif(imat .eq. 2) then
  dist = sign(dist1,dist3)
 else
  print *,"wrong imat flag in, 130"
  stop
 endif
 deallocate(xt,yt)



! methodtype 2
elseif(methodtype .eq. 2)then
  x0 = 2.0d0*x-1.0d0
 y0 = 2.0d0*y-1.0d0

 c1 = 0.02d0*sqrt(5.0d0)
 c2 = 0.02d0*sqrt(5.0d0)
 
 tt = atan((y0-c2)/(x0-c1));  
 if((x0-c1) .ge. 0.0d0 .and. (y0-c2) .ge. 0.0d0)then
    ! do nothing
 elseif((x0-c1) .le. 0.0d0 .and. (y0-c2) .gt. 0.0d0)then
    tt = tt + mypi;
 elseif((x0-c1) .lt. 0.0d0 .and. (y0-c2) .lt. 0.0d0)then
    tt = tt +mypi;
 else
    tt = 2.0d0*mypi + tt;
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
endif

end subroutine

