program main
implicit none

integer,parameter           :: N=16
integer,parameter           :: M=32
real(kind=8),parameter      :: rlo=0.25d0-0.1d0
real(kind=8),parameter      :: rhi=0.25d0+0.1d0
real(kind=8),parameter      :: pi=4.0d0*atan(1.0d0)
integer,parameter           :: step=40
integer,parameter           :: testnum=3
real(kind=8),parameter      :: T1=2.0d0
real(kind=8),parameter      :: T2=2.0d0

real(kind=8)                :: r(0:N)
real(kind=8)                :: z(0:M)
real(kind=8)                :: u(0:N,0:M),u_new(0:N,0:M)
real(kind=8)                :: dr,dz
real(kind=8)                :: tau
real(kind=8)                :: kappa
real(kind=8),external       :: f_src

integer                     :: i,j,ts

dr=(rhi-rlo)/N
dz=(2*pi)/M

do i=0,N
 r(i)=rlo+i*dr
enddo
do i=0,M
 z(i)=i*dz
enddo

kappa=1.0d0
tau=0.5d0/kappa*min(((dr)**2.0d0), &
       2.0d0*minval(r)*dr, &
      (minval(r)**2.0d0)*((dz)**2.0d0))

print *,"tau=",tau

open(unit=5,file="psol.dat")



! initial cond, boundary cond

if(testnum .eq. 1)then
do i=0,N
 do j=0,M
  u(i,j)=2.0d0
!  u(i,j) = 1.0d0 + 10.0d0*(r(i)-rlo)
 enddo
enddo

do j=0,M
 u(0,j)=T1
 u(N,j)=T2
enddo

elseif(testnum .eq. 2)then

 do i=0,N
  do j= 0,M
    u(i,j)=T1*(rhi-r(i))/(rhi-rlo) + &
         T2*(r(i)-rlo)/(rhi-rlo) + &
         100.0d0*sin(z(j))*(r(i)-rlo)/(rhi-r(i))
  enddo
 enddo

elseif(testnum .eq. 3)then
 
 do i=0,N
  do j= 0,M
    u(i,j)=2.0d0*(rhi-r(i))/(rhi-rlo) + &
         2.0d0*(r(i)-rlo)/(rhi-rlo) + &
         100.0d0*sin(z(j))*(r(i)-rlo)*(rhi-r(i))
  enddo
 enddo




else
 print *,"testnum invalid"

endif


do j=0,M
 write(5,*) u(0:N,j) 
enddo


!--
do ts=1,step
 do i=1,N-1
  do j=1,M-1
   u_new(i,j)= (1.0d0+kappa*tau*(-2.0d0/dr**2.0d0- &
                       2.0d0/(r(i)**2.0d0*dz**2.0d0)))*u(i,j) &
                + kappa*tau*(  & 
                 (1.0d0/dr**2.0d0+1.0d0/(r(i)*2.0d0*dr))*u(i+1,j) &
                + (1.0d0/dr**2.0d0-1.0d0/(r(i)*2.0d0*dr))*u(i-1,j) &
                + 1.0d0/(r(i)**2.0d0*dz**2.0d0)*u(i,j+1) &
                + 1.0d0/(r(i)**2.0d0*dz**2.0d0)*u(i,j-1)) &
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
 u_new(0,j)=T1
 u_new(N,j)=T2
enddo

 do j=0,M
  write(5,*) u_new(0:N,j) 
 enddo


 do i=0,N
  do j=0,M
   u(i,j)=u_new(i,j)
  enddo
 enddo

enddo


close(5)

end program


function f_src(x,y)
implicit none

real(kind=8)  :: x,y
real(kind=8)  :: f_src

f_src = 0.0d0





return
end function
