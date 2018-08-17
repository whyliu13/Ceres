#undef BL_LANG_CC
#define BL_LANG_FORT

#include "REAL.H"
#include "CONSTANTS.H"
#include "SPACE.H"
#include "BC_TYPES.H"
#include "ArrayLim.H"




module interp_module


use Generalclass
use probcommon_module
use mof_routines_module
use MOF_pair_module
use mmat_FVM

implicit none

integer,parameter   :: sdim = 2 
 

contains
       ! centroid -> target 
       ! xprobe = point to interpolate temperature to. (dx distance from xI)
       ! xI = closest point on interface
       ! xlo = lower left hand corner of domain

      subroutine interpfabTEMP( &
       N, &
       dx, &
       xlo, &
       Temp, &
       vf, &
       centroid, &
       xprobe, &  ! normal probe position
       xI, & ! closest point on interface from a cell center.
       interpolate_temperature)
      IMPLICIT NONE

      REAL_T xlo(sdim)
      REAL_T dx(sdim)
      REAL_T xprobe(sdim)
      REAL_T xI(sdim)
      REAL_T interpolate_temperature
      INTEGER_T dir
      INTEGER_T ic,jc
      INTEGER_T i1,j1
      INTEGER_T isten,jsten
      INTEGER_T ihalf
      INTEGER_T nhalf
      INTEGER_T cell_index(sdim)
      REAL_T xsten(-3:3,sdim)
      REAL_T xsten_stencil(-3:3,sdim)
      REAL_T volcell
      REAL_T cencell(sdim)
      REAL_T T_sten(-1:1,-1:1)
      REAL_T VF_sten(-1:1,-1:1)
      REAL_T XC_sten(-1:1,-1:1,sdim)

      INTEGER       :: N
      REAL(KIND=8)  :: Temp(-1:N,-1:N)
      REAL(KIND=8)  :: vf(-1:N,-1:N)      
      REAL(KIND=8)  :: centroid(-1:N,-1:N,sdim) 


      nhalf=3

       ! find index of cell containing xprobe
       ! x=(i+1/2)dx + xlo
       ! i=(x-xlo)/dx-1/2
      do dir=1,sdim
       if (dx(dir).le.zero) then
        print *,"dx(dir) invalid"
        stop
       endif
       cell_index(dir)=NINT((xprobe(dir)-xlo(dir))/dx(dir)-half)
      enddo

      ic=cell_index(1)
      jc=cell_index(2)

      do dir=1,sdim
      do ihalf=-nhalf,nhalf
       xsten(ihalf,dir)=(cell_index(dir)+half)*dx(dir)+xlo(dir)+ &
         ihalf*dx(dir)/two
      enddo
      enddo

      do i1=-1,1
      do j1=-1,1

       isten=i1+ic
       jsten=j1+jc

       do dir=1,sdim
       do ihalf=-nhalf,nhalf
        if (dir.eq.1) then
         xsten_stencil(ihalf,dir)=(isten+0.5)*dx(dir)+xlo(dir)+ &
          ihalf*dx(dir)/two
        else if (dir.eq.2) then
         xsten_stencil(ihalf,dir)=(jsten+0.5)*dx(dir)+xlo(dir)+ &
          ihalf*dx(dir)/two
        else
         print *,"dir invalid"
         stop
        endif
       enddo
       enddo

       volcell=dx(1)*dx(2)
       do dir=1,sdim
        cencell(dir)=xsten_stencil(0,dir)
       enddo 
  !     xcenter=xsten_stencil(0,1)
  !     ycenter=xsten_stencil(0,2)

        ! temperature grid (temperature data at cell centroids)
        ! material 2 for annulus test.
!       T_sten(i1,j1)=input_temperature(isten,jsten)
        T_sten(i1,j1)=Temp(isten,jsten)

        ! volume fraction grid; (material 2 for annulus test)
!        VF_sten(i1,j1)=input_volume_fraction(isten,jsten)
        VF_sten(i1,j1)=vf(isten,jsten)


        ! centroid grid; (material 2 for annulus test)
        ! it is assumed that these centroids are relative to the cell center.
       do dir=1,sdim
!        XC_sten(i1,j1,dir)= &
!         input_centroid(isten,jsten,dir)+cencell(dir)
        XC_sten(i1,j1,dir)=centroid(isten,jsten,dir)
       enddo

      enddo
      enddo ! i1,j1

      call center_centroid_interchange( &
       dx,xlo, &
       xsten,nhalf, &
       T_sten, & 
       XC_sten, & 
       xI, &
       xprobe, &  
       VF_sten, & 
       interpolate_temperature)

      return 
      end subroutine interpfabTEMP


      subroutine center_centroid_interchange( &
       dx,xlo, & 
       xsten,nhalf, &
       T_sten, &
       XC_sten, &
       xI, &  ! closest point on interface to cell center.
       xprobe, &
       VF_sten, &
       interpolate_temperature)
      implicit none

      REAL_T dx(sdim)
      REAL_T xlo(sdim)
      REAL_T xI(sdim)
      INTEGER_T nhalf
      REAL_T xsten(-nhalf:nhalf,sdim)
      REAL_T T_sten(-1:1,-1:1)
      REAL_T XC_sten(-1:1,-1:1,sdim)
      REAL_T VF_sten(-1:1,-1:1)
      REAL_T wt_sten(-1:1,-1:1)
      REAL_T interpolate_temperature

      INTEGER_T i,j
      INTEGER_T i1,j1
      INTEGER_T dir
      REAL_T dist
      REAL_T, dimension(:,:), allocatable :: AA
      REAL_T, dimension(:,:), allocatable :: AAcopy
      REAL_T BB(sdim+1)
      REAL_T BBcopy(sdim+1)
      REAL_T xtemp(sdim)
      REAL_T xprobe(sdim)
      REAL_T xlive(sdim)
      REAL_T delx(sdim+1)
      REAL_T GRADTEMP(sdim+1)
      REAL_T TMIN,TMAX,T_test,wt_sum,VF
      INTEGER_T mat_ncomp,matstatus
      REAL_T T_hold
      INTEGER_T own_flag
      REAL_T wt_local
      REAL_T TSAT

      REAL(KIND=8),EXTERNAL :: exact_temperature                    ! use temperature on the interface

      if (nhalf.ne.3) then
       print *,"nhalf invalid"
       stop
      endif

       !   T(x)=T(xI) + grad T dot (x-xbase)

       ! this is the Dirichlet temperature condition at xI (the
       ! closest point)
  !    TSAT=temperature_function(xI)                                       ! set TSAT
  !    TSAT=exact_temperature(x,y,t,im,probtype_in,nmat_in,alpha,dclt_flag)

      TMAX=TSAT
      TMIN=TMAX

      wt_sum=zero
       
      do i=-1,1
      do j=-1,1

       xtemp(1)=xsten(2*i,1)
       xtemp(2)=xsten(2*j,2)
       do dir=1,sdim
        xlive(dir)=XC_sten(i,j,dir)
       enddo ! dir=1..sdim

       VF=VF_sten(i,j)

       own_flag=0
       if (abs(VF).le.VOFTOL) then
        own_flag=0
       else if ((VF.ge.VOFTOL).and.(VF.le.one+VOFTOL)) then
        own_flag=1
       else
        print *,"VF invalid"
        stop
       endif
 
       wt_local=zero
 
       if (own_flag.eq.0) then
        T_sten(i,j)=TSAT
        wt_local=VOFTOL
       else if (own_flag.eq.1) then
        dist=VOFTOL
        do dir=1,sdim
         dist=dist+(xlive(dir)-xprobe(dir))**2/(dx(1)**2)
        enddo ! dir
        if (dist.le.zero) then
         print *,"dist invalid"
         stop
        endif
        wt_local=one/dist
        if (VF.le.zero) then
         print *,"VF invalid"
         stop
        endif
        wt_local=wt_local*(VF**2)
       else
        print *,"own_flag invalid"
        stop
       endif

       if (wt_local.le.zero) then
        print *,"wt_local invalid"
        stop
       endif

       wt_sten(i,j)=wt_local

       T_test=T_sten(i,j)
       if (T_test.gt.TMAX) then
        TMAX=T_test
       endif
       if (T_test.lt.TMIN) then
        TMIN=T_Test
       endif
       if (abs(T_test).ge.1.0D+50) then
        print *,"T_test bust"
        stop
       endif

       wt_sum=wt_sum+wt_local

      enddo
      enddo

      if (wt_sum.le.zero) then
       print *,"wt_sum invalid"
       stop
      endif

      do i=-1,1
      do j=-1,1
       wt_sten(i,j)=wt_sten(i,j)/wt_sum
      enddo
      enddo

      mat_ncomp=sdim

      allocate(AA(mat_ncomp,mat_ncomp))
      allocate(AAcopy(mat_ncomp,mat_ncomp))

      do i=1,mat_ncomp
      do j=1,mat_ncomp
       AA(i,j)=zero
      enddo
      enddo
      do i=1,sdim+1
       BB(i)=zero
      enddo

      do i=1,mat_ncomp
      do j=1,mat_ncomp
       do i1=-1,1
       do j1=-1,1
        xtemp(1)=xsten(2*i1,1)
        xtemp(2)=xsten(2*j1,2)
        do dir=1,sdim
         xlive(dir)=XC_sten(i1,j1,dir)
        enddo ! dir

        do dir=1,sdim
         delx(dir)=(xlive(dir)-xI(dir))/dx(dir)
        enddo ! dir
        delx(sdim+1)=one

        AA(i,j)=AA(i,j)+wt_sten(i1,j1)*delx(i)*delx(j)
        if (j.eq.1) then
         BB(i)=BB(i)+wt_sten(i1,j1)*delx(i)*(T_sten(i1,j1)-TSAT)
        endif
       enddo ! j1
       enddo ! i1
       AAcopy(i,j)=AA(i,j)
       if (j.eq.1) then
        BBcopy(i)=BB(i)
       endif
      enddo   ! j
      enddo   ! i
 
      call matrix_solve(AA,GRADTEMP,BB,matstatus,mat_ncomp)

      if (matstatus.eq.1) then
       T_hold=TSAT

       do dir=1,sdim
        delx(dir)=(xprobe(dir)-xI(dir))/dx(dir)
       enddo ! dir
       delx(sdim+1)=one

       do dir=1,mat_ncomp
        T_hold=T_hold+GRADTEMP(dir)*delx(dir)
       enddo
       if (T_hold.lt.TMIN) then
        T_hold=TMIN
       endif
       if (T_hold.gt.TMAX) then
        T_hold=TMAX
       endif
       if (T_hold.le.zero) then
        print *,"temperature underflow in center_centroid_interchange"
        stop
       endif
      else if (matstatus.eq.0) then
       print *,"WARNING: 1st order failed center_centroid: doing zeroth order"
       print *,"TSAT=",TSAT
       print *,"wt_sum=",wt_sum
       print *,"mat_ncomp=",mat_ncomp
       print *,"xsten(0,dir) ",xsten(0,1),xsten(0,2)
       print *,"xprobe ",xprobe(1),xprobe(2)
       print *,"dx ",dx(1),dx(2)
       print *,"xI ",xI(1),xI(2)
       do i=1,mat_ncomp
       do j=1,mat_ncomp
        print *,"i,j,AAcopy,AA ",i,j,AAcopy(i,j),AA(i,j)
       enddo
       enddo
       do i=1,mat_ncomp
        print *,"i,BBcopy,BB ",i,BBcopy(i),BB(i)
       enddo
        
       T_hold=TSAT
      else
       print *,"matstatus has an invalid value"
       stop
      endif
 
      interpolate_temperature=T_hold

      deallocate(AA)
      deallocate(AAcopy)

      return
      end subroutine center_centroid_interchange
   
    

      subroutine polar_cart_interpolate(Np,Mp,upolar,pcenter,rlo,rhi,x_in, ux)
      implicit none
        
      integer,intent(in)          :: Np,Mp
      real(kind=8),intent(in),dimension(0:Np,0:Mp)  :: upolar
      real(kind=8),intent(in)          :: pcenter(2)
      real(kind=8),intent(in)          :: rlo,rhi
      real(kind=8),intent(in)          :: x_in(2)

      real(kind=8)                     :: ux
      integer                          :: i,j
      real(kind=8)                     :: xr,xz
      integer                          :: rl,zl      
      real(kind=8)                     :: dr,dz
      real(kind=8)                     :: r_interp1,r_interp2
      real(kind=8)                     :: alpha,beta

      dr = (rhi-rlo)/real(Np,8)     
      dz = 2.0d0*pi/real(Mp,8)

      xr = sqrt((x_in(1)-pcenter(1))**2.0d0 + (x_in(2)-pcenter(2))**2.0d0)
      rl = floor((xr-rlo)/dr) 
      alpha= (xr-(rlo+rl*dr))/dr
      
      call rad_cal(x_in,pcenter,xz)
      zl = floor(xz/dz)
      beta=(xz-zl*dz)/dz

      r_interp1=upolar(rl,zl)*(1.0d0-alpha)+upolar(rl+1,zl)*alpha
      r_interp2=upolar(rl,zl+1)*(1.0d0-alpha)+upolar(rl+1,zl+1)*alpha
      
      ux=r_interp1*(1.0d0-beta)+r_interp2*beta
   
    
      end subroutine polar_cart_interpolate

end module



