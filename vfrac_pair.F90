#undef BL_LANG_C
#define BL_LANG_FORT

#include "REAL.H"
#include "CONSTANTS.H"
#include "SPACE.H"
#include "BC_TYPES.H"
#include "ArrayLim.H"


! Generalclass is in vof_cisl.F90
! eps is defined in vof_cisl.F90
module MOF_pair_module
  USE Generalclass
  use probcommon_module
  use mof_routines_module
  use geometry_intersect_module

  implicit none




contains

  ! ==========================================================================
  ! =========================================================================
  !  Given a face of the cell, output the fraction pair along the face  
  !==========================================================================
  !==========================================================================

  subroutine vfrac_pair_along_side(nmat,sdim,frac_left,frac_right, &
       x_left,x_right,frac_pair)
    implicit none

    integer,intent(in)              :: nmat, sdim
    real(kind=8)                    :: frac_left(nmat),frac_right(nmat)
    real(kind=8)                    :: x_left(sdim,nmat),x_right(sdim,nmat)
    real(kind=8)                    :: left_sum, right_sum

    integer                         :: i,j,m
    real(kind=8)                    :: frac_pair(nmat,nmat) ! vfrac_pair(m_left, m_right)

    integer                         :: change, ml, mr, ml_crit, mr_crit
    real(kind=8)                    :: dcrit        ! >>set a dcrit value <<<<<<<
    real(kind=8)                    :: dist

    ! if F>1-eps => F=1
    ! if F< eps => F=0
    ! renormalize so that sum F_{L or R}=1

    do m = 1,nmat
       if(frac_left(m) .lt. eps) then
          frac_left(m) = 0.0d0
       elseif(frac_left(m) .gt. 1.0d0-eps)then
          frac_left(m) = 1.0d0
       endif
       if(frac_right(m) .lt. eps) then
          frac_right(m) = 0.0d0
       elseif(frac_right(m) .gt. 1.0d0-eps)then
          frac_right(m) = 1.0d0
       endif
    enddo


    left_sum = 0.0d0
    right_sum = 0.0d0
    do m = 1,nmat
       left_sum = left_sum + frac_left(m)
       right_sum = right_sum + frac_right(m)
    enddo

    do m = 1,nmat
       frac_left(m) = frac_left(m) / left_sum
       frac_right(m) = frac_right(m) / right_sum
    enddo


   !print *,"left",frac_left
    !print *,"right",frac_right


    do ml=1,nmat
    do mr=1,nmat
     frac_pair(ml,mr) = 0.0d0
    enddo
    enddo



    do m = 1,nmat
       if(frac_left(m) .eq. 0.0d0  .or. frac_right(m) .eq. 0.0d0) then
          frac_pair(m,m) = 0.0d0

       else if (abs(frac_left(m)-frac_right(m))<tol) then
          frac_pair(m,m) = frac_right(m)
          frac_right(m) = 0.0d0
          frac_left(m) = 0.0d0
       elseif(frac_left(m) .lt. frac_right(m)) then
          frac_pair(m,m) = frac_left(m)

          ! update centroid x_right(m) =  (F_R(m)*x_R(m) - F_L(m)*x_L(m))
          !                               /(F_R(m) - F_L(m))
          do i = 1,sdim
             x_right(i,m) = (frac_right(m)* x_right(i,m) - frac_left(m)* x_left(i,m)) &
                  & /(frac_right(m) - frac_left(m))
          enddo
          frac_right(m) = frac_right(m) - frac_left(m)
          frac_left(m) = 0.0d0

       elseif(frac_right(m) .lt. frac_left(m)) then
          frac_pair(m,m) = frac_right(m)
          do i = 1,sdim
             x_left(i,m) = (frac_left(m)* x_left(i,m) - frac_right(m)* x_right(i,m)) &
                  &  /(frac_left(m) - frac_right(m))
          enddo
          frac_left(m) = frac_left(m) - frac_right(m)
          frac_right(m) = 0.0d0
       else
          print *, ">>>>>>>>>>>>>>>>>>>>"
          print *, 1,"out of cases."   
          print *, frac_left(m),frac_right(m)
          print *, ">>>>>>>>>>>>>>>>>>>>"       
          stop
       endif
    enddo


    change = 1
    do while(change .eq. 1) 
       change = 0
       dcrit = 0.0d0            ! >>>>>>>>>>set a dcrit value <<<<<<<<<<<<<<<<<
       ml_crit = 0
       mr_crit = 0

       do ml = 1, nmat
          do mr = 1, nmat
             if(ml .ne. mr) then
                if(frac_left(ml) .gt. 0.0d0  .and. frac_right(mr) .gt. 0.0d0) then
                   call two_points_dist(sdim,x_left(ml,:),x_right(mr,:),dist)
                   if((change .eq. 0) .or. (dist .lt. dcrit)) then
                      dcrit = dist
                      change = 1
                      ml_crit = ml
                      mr_crit = mr
                   endif
                endif
             endif
          enddo ! mr
       enddo   !ml

       if(change .eq. 1) then
          ml = ml_crit
          mr = mr_crit
          if (abs(frac_left(ml)-frac_right(mr))<tol) then
             frac_pair(ml,mr) = frac_left(ml)
             frac_left(ml) = 0.0d0
             frac_right(mr) = 0.0d0
          elseif(frac_left(ml) .lt. frac_right(mr)) then                
             frac_pair(ml,mr) = frac_left(ml)
             do i = 1,sdim
                x_right(i,mr) = (frac_right(mr)* x_right(i,mr) - frac_left(ml)* x_left(i,ml)) &
                     /(frac_right(mr) - frac_left(ml))
             enddo
             frac_right(mr) = frac_right(mr) - frac_left(ml)
             frac_left(ml) = 0.0d0
          elseif(frac_left(ml) .gt. frac_right(mr)) then                
             frac_pair(ml,mr) = frac_right(mr)
             do i = 1,sdim
                x_left(i,ml) = (frac_left(ml)* x_left(i,ml) - frac_right(mr)* x_right(i,mr)) &
                     /(frac_left(ml) - frac_right(mr))
             enddo
             frac_left(ml) = frac_left(ml) - frac_right(mr)
             frac_right(mr) = 0.0d0
          else
             print *,2,"out of cases"
             stop
          endif
       endif

    enddo   ! end do while



  end subroutine vfrac_pair_along_side

  !----------------------------------------------------------
  subroutine two_points_dist(sdim,x1,x2,dist)
    implicit none

    integer              :: sdim
    real(kind=8)         :: x1(sdim),x2(sdim)
    integer              :: i
    real(kind=8)         :: dist

    dist = 0.0d0
    do i = 1,sdim
       dist = dist + (x1(i) - x2(i))**2.0d0
    enddo
    dist = sqrt(dist)

    return
  end subroutine two_points_dist

  !===================================================================
  !   calc all the pairs along the faces of a cell
  !====================================================================
  subroutine vfrac_pair_cell(nmat,sdim,dx, &
   ext_face_sten,thin_cen_sten,frac_pair_cell)
    implicit none

    integer,intent(in) :: sdim, nmat
    real(kind=8)       :: dx(sdim)
    real(kind=8)       :: ext_face_sten(-1:1,sdim,nmat+1,sdim,2)
    real(kind=8)       :: area(sdim)

    integer        :: i,j,k,im,im1,im2,dir,side,il,ir,ii,dir2
    real(kind=8)   :: frac_left(nmat),frac_right(nmat)
    real(kind=8)   :: x_left(sdim,nmat)
    real(kind=8)   :: x_right(sdim,nmat)
    real(kind=8)   :: frac_pair_cell(nmat,nmat,sdim,2)

    real(kind=8)   :: thin_cen_sten(-1:1,sdim,sdim,nmat,sdim,2) 
    real(kind=8)   :: frac_pair(nmat,nmat)

    do im1=1,nmat
    do im2=1,nmat
    do dir=1,sdim
    do side=1,2
     frac_pair_cell(im1,im2,dir,side) = 0.0d0 ! im2 is inside 
    enddo
    enddo
    enddo
    enddo
    do dir = 1, sdim
     do side = 1,2
      if(side .eq. 1) then
       do im = 1, nmat
        frac_left(im) = ext_face_sten(-1,dir,im,dir,2)
        frac_right(im) = ext_face_sten(0,dir,im,dir,1)  ! inside
       enddo
       do dir2=1,sdim
       do im = 1, nmat
        x_left(dir2,im)=thin_cen_sten(-1,dir,dir2,im,dir,2)
        x_right(dir2,im)=thin_cen_sten(0,dir,dir2,im,dir,1)
       enddo
       enddo
       call vfrac_pair_along_side(nmat,sdim,frac_left,frac_right, &
                x_left,x_right,frac_pair)
       do im1=1,nmat
       do im2=1,nmat
        frac_pair_cell(im1,im2,dir,side) = frac_pair(im1,im2)
       enddo
       enddo

      elseif(side .eq. 2) then
       do im = 1, nmat
        frac_left(im) = ext_face_sten(1,dir,im,dir,1)
        frac_right(im) = ext_face_sten(0,dir,im,dir,2)  ! inside
       enddo
       do dir2=1,sdim
       do im = 1, nmat
        x_left(dir2,im)=thin_cen_sten(1,dir,dir2,im,dir,1)
        x_right(dir2,im)=thin_cen_sten(0,dir,dir2,im,dir,2) ! inside
       enddo
       enddo
       call vfrac_pair_along_side(nmat,sdim,frac_left,frac_right, &
            x_left,x_right,frac_pair)
       do im1=1,nmat
       do im2=1,nmat
        frac_pair_cell(im1,im2,dir,side) = frac_pair(im1,im2)
       enddo 
       enddo 
      else
       print *,"side invalid in vfrac_pair_cell"   
       stop
      endif

     enddo ! side
    enddo ! dir



  end subroutine vfrac_pair_cell
  !/////////////////////////////////////////////////////////////////////////
  !/////////////////////////////////////////////////////////////////////////




  ! =========================================================================
  ! =========================================================================
  !   From legacy code, perturb interface and cell faces
  ! 
  ! =========================================================================
  !=========================================================================

  !------------------perturb faces of the cell-------------------
  subroutine ptb_ext(ngeom_recon, &
   nmat,sdim,dx,mofdata, &
   xsten, &
   ext_facefrac_cell, &
   multi_cen_cell)
    implicit none

    integer,         intent(in) :: ngeom_recon,nmat, sdim
    real(kind=8)                :: h       !  >__<
    real(kind=8),    intent(in) :: mofdata(ngeom_recon*nmat)
    real(kind=8),    intent(in) :: dx(sdim)
    real(kind=8)                :: mofdataproject(nmat*ngeom_recon)
    integer                     :: project_status(nmat)

    real(kind=8)                :: vcenter(nmat)
    INTEGER                     :: sorted_list(nmat)
    real(kind=8)                :: xsten(-3:3,SDIM)
    real(kind=8)                :: xsten_thin(-1:1,SDIM)
    REAL(kind=8)                :: dxthin
    integer                     :: vofcomp
    REAL(kind=8)                :: ext_facefrac_cell(nmat+1,sdim,2)

    integer                     :: im,dir,side,dir2,ivert,isten,im1
    integer                     :: im_crit,im_exclude
    integer                     :: shapeflag,nhalf, nhalf_thin, nmax
    REAL(kind=8)                :: multi_volume(nmat)
    REAL(KIND=8)                :: multi_cen(SDIM,nmat)
    REAL(KIND=8)                :: multi_area(nmat)
    REAL(KIND=8)                :: total_vol
    REAL(kind=8)                :: xtrilistuncapt(SDIM+1,SDIM,POLYGON_LIST_MAX)
    REAL(kind=8)                :: dummy_tri(SDIM+1,SDIM)  
    real(kind=8)                :: multi_cen_cell(sdim,nmat,sdim,2) 

    integer                     :: bfact

    if (ngeom_recon.ne.2*sdim+3) then
     print *,"ngeom_recon invalid"
     stop
    endif

    h = dx(1)

    bfact = 1
    nhalf=3 
    nhalf_thin=1 
    nmax=POLYGON_LIST_MAX

    do dir=1,sdim
    do im=1,nmat
    do dir2=1,sdim
    do side=1,2
     multi_cen_cell(dir,im,dir2,side) = 0.0d0
    enddo
    enddo
    enddo
    enddo
    do im=1,nmat+1
    do dir2=1,sdim
    do side=1,2
     ext_facefrac_cell(im,dir2,side) = 0.0d0
    enddo
    enddo
    enddo

    do ivert=1,SDIM+1
       do dir2=1,SDIM
          dummy_tri(ivert,dir2)=zero
       enddo
    enddo

    ! vcenter = volume fraction 
    do im=1,nmat
       vofcomp=(im-1)*ngeom_recon+1
       vcenter(im)=mofdata(vofcomp)    ! extract volume fraction
    enddo ! im
    ! sort the volume fractions.  
    im_exclude = 0
    call sort_vof(vcenter,im_exclude,sorted_list,nmat)
    im_crit=sorted_list(1)  ! vcenter(im_crit) is largest value.
    if ((im_crit.lt.1).or.(im_crit.gt.nmat)) then
       print *,"im_crit invalid"
       stop
    endif


    ! ext_facefrac_cell: nmat+1,sdim,2 components
    ! last 2*sdim components are areas of the cell face (dir,side)
    ! e.g. in 2d uniform grid: ext_facefrac_cell=h
    ! first nmat*2*sdim components are areafractions. 0<=A<=1
    ! *********cal  tot  area  of faces
    do dir=1,SDIM
     do side=1,2
      ext_facefrac_cell(nmat+1,dir,side) = h      ! 2d
     enddo
    enddo ! dir,side
    ! *************************

    if (vcenter(im_crit) .ge. one-FACETOL) then
     do dir=1,SDIM
      dxthin=FACETOL*(xsten(1,dir)-xsten(-1,dir))
      do side=1,2
       ext_facefrac_cell(im_crit,dir,side)=one
       !---------------------------------------------
       if(side .eq. 1)then
        do dir2 = 1,sdim
         if(dir2 .eq. dir)then
          multi_cen_cell(dir2,im_crit,dir,side) = xsten(-1,dir) 
         else
          multi_cen_cell(dir2,im_crit,dir,side) = xsten(0,dir2)
         endif
        enddo
       elseif(side .eq. 2)then
        do dir2 = 1,sdim
         if(dir2 .eq. dir)then
          multi_cen_cell(dir2,im_crit,dir,side) = xsten(1,dir) 
         else
          multi_cen_cell(dir2,im_crit,dir,side) = xsten(0,dir2)
         endif
        enddo
       else
        print *,"side invalid"
        stop
       endif
       !----------------------------------------------------
      enddo
     enddo  ! dir,side
       !----------------------------------------------
       !   multi_cen_cell(sdim,im_crit,sdim,2) 
       !---------------------------------------------
    else if (vcenter(im_crit).gt.zero) then
     shapeflag=0
     do dir=1,SDIM
     do side=1,2
      do isten=-1,1
       do dir2=1,SDIM
        xsten_thin(isten,dir2)=xsten(isten,dir2)
       enddo
      enddo
      dxthin=FACETOL*(xsten(1,dir)-xsten(-1,dir))
      if (side.eq.1) then
       xsten_thin(1,dir)=xsten(-1,dir)+dxthin
      else if (side.eq.2) then
       xsten_thin(-1,dir)=xsten(1,dir)-dxthin
      else
       print *,"side invalid"
       stop
      endif
      xsten_thin(0,dir)=half*(xsten_thin(-1,dir)+xsten_thin(1,dir))
      ! find volumes and areas (not scaled) of materials in
      ! xsten_thin box: xsten_thin(0,dir) = center of thin box
      ! xsten_thin(1,dir) right side in dir direction
      ! xsten_thin(-1,dir) left side 
      call project_slopes_to_face( &
       project_status, &
       bfact,dx,xsten,nhalf, &
       mofdata,mofdataproject, &
       nmat,SDIM,dir,side)

      call multi_get_volume_grid( &
       bfact,dx,xsten,nhalf, &
       mofdataproject, &  ! was mofdata
       xsten_thin,nhalf_thin, &
       dummy_tri, &
       multi_volume, &
       multi_cen, &
       multi_area, &
       levelrz, &
       xtrilistuncapt, &
       nmax,nmat,SDIM, &
       shapeflag,3)

      do im=1,nmat
       do dir2=1,SDIM
        if (dir2.eq.dir) then
         if (side.eq.1) then
          multi_cen_cell(dir2,im,dir,side) = xsten(-1,dir)
         else if (side.eq.2) then
          multi_cen_cell(dir2,im,dir,side) = xsten(1,dir)
         else
          print *,"side invalid"
          stop
         endif
        else 
         multi_cen_cell(dir2,im,dir,side) = multi_cen(dir2,im)
        endif
       enddo 
      enddo 

      total_vol=zero
      do im=1,nmat
       total_vol=total_vol+multi_volume(im)
      enddo
      if (total_vol.gt.zero) then
       do im=1,nmat
        ext_facefrac_cell(im,dir,side)=multi_volume(im)/total_vol
       enddo
      else
       print *,"total_vol invalid"
       stop
      endif
     enddo
     enddo ! dir,side

    else
     print *,"vcenter(im_crit) out of range"
     stop
    endif

    return

  end subroutine ptb_ext
  !------------------------------------------------------------------
  subroutine sort_vof(vfrac_data,im_exclude,sorted_list,nmat)
    implicit none

    integer     ,intent(in)   :: nmat,im_exclude
    real(kind=8),intent(in)   :: vfrac_data(nmat)
    integer                   :: sorted_list(nmat)

    integer                   :: im,changed,nsweeps,swap,do_swap

    do im=1,nmat
       sorted_list(im)=im
    enddo

    changed=1
    nsweeps=0
    do while ((changed.eq.1).and.(nsweeps.lt.nmat-1))
       changed=0
       do im=1,nmat-nsweeps-1
          do_swap=0
          if (sorted_list(im).eq.im_exclude) then
             do_swap=1
          else if (sorted_list(im+1).eq.im_exclude) then
             do_swap=0
          else if (vfrac_data(sorted_list(im)).lt. &
               vfrac_data(sorted_list(im+1))) then

             do_swap=1
          endif
          if (do_swap.eq.1) then
             swap=sorted_list(im)
             sorted_list(im)=sorted_list(im+1)
             sorted_list(im+1)=swap
             changed=1
          endif
       enddo
       nsweeps=nsweeps+1
    enddo

    return

  end subroutine sort_vof

  subroutine get_kappa_driver(alpha,theta1,theta2,im1,im2,kappa,nmat)
  IMPLICIT NONE

  integer im1,im2,nmat
  real(kind=8),intent(in)     :: alpha(nmat)
  real(kind=8),intent(out)    :: kappa
  real(kind=8),intent(in)     :: theta1,theta2
  real(kind=8)                :: t1,t2,theta

  t1=theta1
  t2=theta2 
  if ((t1.lt.0.0).or.(t2.lt.0.0)) then
   print *,"t1 or t2 invalid"
   stop
  endif
  theta=t1+t2
  if (theta.eq.0.0) then
   kappa=0.0
  else if (theta.gt.0.0) then
   t1=t1/theta 
   t2=t2/theta 
   if ((t1.lt.eps).or.(t2.lt.eps)) then
    kappa=0.0
   else if (alpha(im1).eq.0.0) then
    kappa=0.0
   else if (alpha(im2).eq.0.0) then
    kappa=0.0
   else if ((alpha(im1).gt.0.0).and.(alpha(im2).gt.0.0)) then
    kappa=(t1+t2)*alpha(im1)*alpha(im2)/ &
     (t1*alpha(im2)+t2*alpha(im1))
   else
    print *,"alpha invalid"
    stop
   endif
  else
   print *,"theta invalid"
   stop
  endif

  return
  end subroutine get_kappa_driver

   ! im1 is the target material
  subroutine get_kappa(dclt_test,alpha,xsten,dir,side, &
   cen1,cen2,im1,im2,kappa,nmat,sdim)
  IMPLICIT NONE

  integer im1,im2,nmat,dir,side,sdim
  integer,intent(in)     :: dclt_test 
  real(kind=8),intent(in)     :: alpha(nmat)
  real(kind=8),intent(in)     :: xsten(-3:3,sdim)
  real(kind=8),intent(in)     :: cen1(sdim)
  real(kind=8),intent(in)     :: cen2(sdim)
  real(kind=8),intent(out)    :: kappa
  real(kind=8)                :: theta1,theta2


 if(dclt_test .eq. 0)then 
  if (side.eq.1) then
   theta1=abs(xsten(-1,dir)-cen1(dir)) 
   theta2=abs(xsten(-1,dir)-cen2(dir)) 
  else if (side.eq.2) then
   theta1=abs(xsten(1,dir)-cen1(dir)) 
   theta2=abs(xsten(1,dir)-cen2(dir)) 
  else
   print *,"side invalid"
   stop
  endif
  call get_kappa_driver(alpha,theta1,theta2,im1,im2,kappa,nmat)

 elseif(dclt_test .eq. 1)then
  if(im1 .eq. im2)then
   if (side.eq.1) then
    theta1=abs(xsten(-1,dir)-cen1(dir)) 
    theta2=abs(xsten(-1,dir)-cen2(dir)) 
   else if (side.eq.2) then
    theta1=abs(xsten(1,dir)-cen1(dir)) 
    theta2=abs(xsten(1,dir)-cen2(dir)) 
   else
    print *,"side invalid"
    stop
   endif
   call get_kappa_driver(alpha,theta1,theta2,im1,im2,kappa,nmat)
  else 
   kappa = alpha(im1)
  endif
 
 else
  print *,"invalid dclt_test flag"
  stop
 endif

 ! print *,"get_kappa", kappa

  return
  end subroutine get_kappa

  subroutine get_kappa_int(alpha,dist,im1,im2,kappa,nmat)
  IMPLICIT NONE

  integer im1,im2,nmat
  real(kind=8),intent(in)     :: alpha(nmat)
  real(kind=8),intent(out)    :: kappa
  real(kind=8)                :: dist(nmat,nmat)           
  real(kind=8)                :: theta1,theta2
 
  theta1=dist(im1,im2)
  theta2=dist(im2,im1)
  call get_kappa_driver(alpha,theta1,theta2,im1,im2,kappa,nmat)

  return
  end subroutine get_kappa_int

  !---------------------------------------------------------------
  subroutine dist_point_line(sdim,nmat,dx,mofdata, &
   im_line,im_point1,im_point2,dist)
    implicit none

    integer,intent(in)          :: sdim,nmat
    integer,intent(in)          :: im_line
    integer,intent(in)          :: im_point1
    integer,intent(in)          :: im_point2
    real(kind=8),intent(in)     :: mofdata((sdim*2+3)*nmat)
    real(kind=8),    intent(in) :: dx(sdim)
    real(kind=8)                :: dist(nmat,nmat)           

    integer                     :: ngeom_recon
    integer                     :: j,dir,vofcomp
    real(kind=8)                :: m(sdim)
    real(kind=8)                :: cent(sdim)
    real(kind=8)                :: intercept
    real(kind=8)                :: modm

    ngeom_recon = sdim*2+3

    if ((im_line.lt.1).or.(im_line.gt.nmat).or. &
        (im_point1.lt.1).or.(im_point1.gt.nmat).or. &
        (im_point2.lt.1).or.(im_point2.gt.nmat)) then
     print *,"im_line or im_point invalid"
     stop
    endif
     ! vfrac,cen,order,slope,intercept
    vofcomp=ngeom_recon*(im_line-1)+1
    do dir = 1,sdim
     m(dir) = mofdata(vofcomp+sdim+1+dir)
    enddo
    intercept = mofdata(vofcomp+2*sdim+2)

    modm = 0.0d0
    do dir = 1,sdim
       modm = modm + m(dir)*m(dir)
    enddo
    modm = sqrt(modm)
    if(modm .eq. 0.0d0)then
       print *,"alert, modm = 0"
       stop
    endif

    vofcomp=ngeom_recon*(im_point1-1)+1
    cent = 0.0d0
    do dir = 1,sdim   
     cent(dir) =  mofdata(vofcomp + dir)
    enddo

    dist(im_point1,im_point2) = abs(dot_product(m,cent)+intercept)/modm    

    vofcomp=ngeom_recon*(im_point2-1)+1
    cent = 0.0d0
    do dir = 1,sdim   
     cent(dir) =  mofdata(vofcomp + dir)
    enddo

    dist(im_point2,im_point1) = abs(dot_product(m,cent)+intercept)/modm    

  end subroutine dist_point_line
  !------------------------------------------------------------
  !   perturb the interface inside the cell
  !---------------------------------------------------------
  subroutine ptb_int(ngeom_recon, &
    nmat,sdim,mofdata,dx,xsten, &
    int_facefrac_cell, &
    int_face_normal_cell, &
    dist_to_int_cell)
    implicit none

    integer,         intent(in) :: ngeom_recon,nmat, sdim
    real(kind=8),    intent(in) :: mofdata(ngeom_recon*nmat)
    real(kind=8),    intent(in) :: dx(sdim)

    real(kind=8)                :: vcenter(nmat)
    INTEGER                     :: sorted_list(nmat)
    real(kind=8)                :: xsten(-3:3,SDIM)
    real(kind=8)                :: xsten_thin(-1:1,SDIM)
    REAL(kind=8)                :: dxthin
    integer                     :: vofcomp
    REAL(kind=8)                :: int_facefrac_cell(nmat+1,nmat)
    REAL(kind=8)                :: int_face_normal_cell(nmat,nmat,sdim)

    integer                     :: im,dir,side,im_opp,dir2
    integer                     :: im1,im2
    integer                     :: im_crit,im_exclude
    integer                     :: ivert
    integer                     :: shapeflag,nhalf,nhalf_thin,nmax
    REAL(kind=8)                :: multi_volume(nmat)
    real(kind=8)                :: multi_volume_offset(nmat)
    REAL(KIND=8)                :: multi_cen(SDIM,nmat)
    REAL(KIND=8)                :: multi_cen_offset(SDIM,nmat)
    REAL(KIND=8)                :: multi_area(nmat)
    REAL(KIND=8)                :: total_vol
    REAL(kind=8)                :: xtrilistuncapt(SDIM+1,SDIM,POLYGON_LIST_MAX)
    REAL(kind=8)                :: dummy_tri(SDIM+1,SDIM) 


    real(kind=8)                :: mofdatavalid(ngeom_recon*nmat)
    REAL(kind=8)                :: multi_area_offset(nmat)
    REAL(kind=8)                :: total_face,uncaptured_volume
    INTEGER                     :: irank,testflag
    REAL(kind=8)                :: intercept
    REAL(kind=8)                :: local_facefrac(nmat)
    REAL(kind=8)                :: dist_to_int_cell(nmat,nmat)
    integer                     :: bfact
    real(kind=8)                :: dist_tol

    integer                     :: is_processed(nmat)


    if (ngeom_recon.ne.2*sdim+3) then
     print *,"ngeom_recon invalid"
     stop
    endif
    dist_tol = FACETOL*dx(1)

    bfact = 1
    nhalf=3
    nmax=POLYGON_LIST_MAX
    !print *,nmax
    !nmax = 400

    do im1=1,nmat
    do im2=1,nmat
     dist_to_int_cell(im1,im2) = 0.0d0
    enddo
    enddo

    do ivert=1,SDIM+1
     do dir2=1,SDIM
      dummy_tri(ivert,dir2)=zero
     enddo
    enddo
    ! xsten(0,dir)  dir=1,2  is center of cell
    ! xsten(1,dir)  is right coordinate in dir direction
    ! xsten(-1,dir) is left coordinate in dir direction.
    ! vcenter = volume fraction 
    do im=1,nmat
     vofcomp=(im-1)*ngeom_recon+1
     vcenter(im)=mofdata(vofcomp)
    enddo ! im

    ! sort the volume fractions.  
    im_exclude=0
    call sort_vof(vcenter,im_exclude,sorted_list,nmat)
    im_crit=sorted_list(1) ! vcenter(im_crit) is largest value.
    if ((im_crit.lt.1).or.(im_crit.gt.nmat)) then
       print *,"im_crit invalid"
       stop
    endif

    do im1=1,nmat+1
    do im2=1,nmat
     int_facefrac_cell(im1,im2)=zero
    enddo
    enddo
    do im1=1,nmat
    do im2=1,nmat
    do dir=1,sdim
     int_face_normal_cell(im1,im2,dir)=0.0
    enddo
    enddo
    enddo

    if (vcenter(im_crit).ge.one-FACETOL) then
     ! do nothing, there are no internal faces
    else if (vcenter(im_crit).gt.zero) then
     ! normalize the volume fractions so that the sum is 1.
     call make_vfrac_sum_ok(mofdata,mofdatavalid,nmat,ngeom_recon,SDIM,3)

     do im=1,nmat
      is_processed(im)=0
     enddo

     uncaptured_volume=one
     irank=1

     ! index: (im_inside-1)*(nmat+1)+im_outside
     do while ((irank.le.nmat).and.(uncaptured_volume.gt.zero))
      do im=1,nmat
       vofcomp=(im-1)*ngeom_recon+1
       testflag=NINT(mofdatavalid(vofcomp+SDIM+1))
       if (testflag.eq.irank) then

        is_processed(im)=1

        intercept=mofdatavalid(vofcomp+2*SDIM+2)
        uncaptured_volume=uncaptured_volume-mofdatavalid(vofcomp)
        if (uncaptured_volume.lt.FACETOL) then
         uncaptured_volume=zero
        endif
        if (uncaptured_volume.gt.zero) then ! we have a valid interface.
         ! dist=intercept+slopes dot (x-x0)
         ! perturb interface into the other materials
         mofdatavalid(vofcomp+2*SDIM+2)=intercept+half*FACETOL*dx(1)
         shapeflag=0

         call multi_get_volume_grid( &
          bfact,dx,xsten,nhalf, &
          mofdatavalid, &
          xsten,nhalf, &
          dummy_tri, &
          multi_volume_offset, &
          multi_cen_offset, &
          multi_area_offset, &
          levelrz, &
          xtrilistuncapt, &
          nmax,nmat,SDIM, &
          shapeflag,3)

         mofdatavalid(vofcomp+2*SDIM+2)=intercept

         call multi_get_volume_grid( &
          bfact,dx,xsten,nhalf, &
          mofdatavalid, &
          xsten,nhalf, &
          dummy_tri, &
          multi_volume, &
          multi_cen, &
          multi_area, &
          levelrz, &
          xtrilistuncapt, &
          nmax,nmat,SDIM, &
          shapeflag,3)

         if (multi_volume_offset(im).gt.multi_volume(im)) then

          if (multi_area(im).gt.zero) then

           do im_opp=1,nmat
            local_facefrac(im_opp)=zero
           enddo
           total_face=zero

           do im_opp=1,nmat
            if (im_opp.ne.im) then
             if (is_processed(im_opp).eq.0) then
              local_facefrac(im_opp)= &
               abs(multi_volume(im_opp)-multi_volume_offset(im_opp))
              total_face=total_face+local_facefrac(im_opp)
             else if (is_processed(im_opp).eq.1) then
              ! do nothing
             else
              print *,"is_processed invalid"
              stop
             endif
            endif
           enddo !im_opp

           if (total_face.gt.zero) then
            do im_opp=1,nmat
             if (im_opp.ne.im) then
              local_facefrac(im_opp)=local_facefrac(im_opp)/total_face
              int_facefrac_cell(im,im_opp)=local_facefrac(im_opp)

              if (local_facefrac(im_opp).gt.0.0) then
                ! vfrac,cen,order,slope,intercept
               vofcomp=(im-1)*ngeom_recon+1
               do dir=1,sdim
                int_face_normal_cell(im,im_opp,dir)=-mofdata(vofcomp+sdim+1+dir)
                int_face_normal_cell(im_opp,im,dir)=mofdata(vofcomp+sdim+1+dir)
               enddo
               call dist_point_line(sdim,nmat,dx,mofdata, &
                im,im,im_opp,dist_to_int_cell) 
              endif ! local_facefrac(im_opp) > 0
             endif ! im_opp<>im
            enddo ! im_opp
            int_facefrac_cell(nmat+1,im)=multi_area(im)
           else
            print *,"opposite materials disappeared"
            stop
           endif
          else
           print *,"im boundary disappeared"
           stop
          endif
         else
          print *,"im region should grow"
          stop
         endif

        endif ! uncaptured_volume>0
       endif  ! testflag=irank
      enddo ! im
      irank=irank+1
     enddo  ! while irank<=nmat and uncaptured_volume>0 

    else
     print *,"vcenter(im_crit) out of range"
     stop
    endif

    return
  end subroutine ptb_int
  !-----------------------------------------------------------
  !///////////////////////////////////////////////////////////////////////

  !==================================================================
  !===================================================================
  !  calculate the gradient and diffusion operator 
  ! 
  !==================================================================
  !===================================================================
  subroutine cell_diag_cal( &
    dclt_test, &
    ngeom_recon, &
    sdim,nmat,dx, &
    mofdata_cell, &
    alpha,  &
    mat_cen_sten, &
    im_in, &
    frac_pair_cell, &
    int_face, &
    int_face_normal, &
    dist_to_int, &
    xsten, &
    diag)
    implicit none

    integer,intent(in)       :: ngeom_recon,sdim, nmat,dclt_test
    real(kind=8)             :: mat_cen_sten(-1:1,-1:1,nmat,sdim)
    real(kind=8),intent(in)  :: frac_pair_cell(nmat,nmat,sdim,2)
    real(kind=8),intent(in)  :: int_face(nmat,nmat)
    real(kind=8),intent(in)  :: int_face_normal(nmat,nmat,sdim)
    real(kind=8),intent(in)  :: mofdata_cell(ngeom_recon*nmat)
    real(kind=8),intent(in)  :: dx(sdim)
    real(kind=8),intent(in)  :: alpha(nmat)
    real(kind=8),intent(in)  :: xsten(-3:3,sdim)
    real(kind=8),intent(in)  :: dist_to_int(nmat,nmat)
    integer     ,intent(in)  :: im_in

    integer                  :: i,j,im,im1,im2
    integer                  :: dir, side,dir3
    integer                  :: dir2,faceid,ii,jj

    real(kind=8)             :: kappa
    real(kind=8)             :: Ltemp
    real(kind=8)             :: coef,AFRAC

    real(kind=8),intent(out) :: diag

    integer                  :: nface(sdim)
    real(kind=8)             :: xface(sdim)        


    real(kind=8)             :: nf(sdim),n1(sdim)
    ! ----------
    real(kind=8)             :: dclt_ratio,dclt_ratio1
    real(kind=8)             :: dclt_diff(SDIM),dclt_diff1(SDIM)
    real(kind=8)             :: dclt_pp(SDIM),dclt_pp1(SDIM)


   
  write(2,*)  "new diag cal loop********************************************"
  write(2,*) "**************************************************************"
  write(2,*) "considering mat", im_in

    if (ngeom_recon.ne.2*sdim+3) then
     print *,"ngeom_recon invalid"
     stop
    endif
    if ((dclt_test.ne.0).and.(dclt_test.ne.1)) then
     print *,"dclt_test invalid"
     stop
    endif 
    diag = 0.0d0
    ! cell boundary(external flux) first
    im1 = im_in   ! im1 the material need to be calculated



    do im2 = 1,nmat  ! different mat im2 outside the cell
     do dir = 1,sdim
      do side = 1,2

       AFRAC=frac_pair_cell(im2,im1,dir,side)

       if (abs(AFRAC).le.eps) then
        AFRAC=0.0d0
       else if (abs(AFRAC-1.0d0).le.eps) then
        AFRAC=1.0d0
       else if ((AFRAC.gt.0.0).and.(AFRAC.lt.1.0)) then
        ! do nothing
       else
        print *,"AFRAC invalid"
        stop
       endif

        ! face: nface dot (x-xface)=0
       if ((AFRAC.gt.0.0d0).and.(AFRAC.le.1.0d0)) then
        write(2,*) "in ext_face calvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv"
       
        if((dir.eq.1).and.(side.eq.1)) then
         faceid=1
         ii=-1
         jj=0
         nface(1)=-1.0
         nface(2)=0.0
         xface(1)=xsten(-1,1)
         xface(2)=xsten(0,2)
         
         write(2,*) "side1"
        elseif((dir.eq.1).and.(side.eq.2)) then
         faceid=2
         ii=1
         jj=0
         nface(1)=1.0
         nface(2)=0.0
         xface(1)=xsten(1,1)
         xface(2)=xsten(0,2)
         write(2,*) "side2"
        elseif((dir.eq.2).and.(side.eq.1)) then
         faceid=3
         ii=0
         jj=-1
         nface(1)=0.0
         nface(2)=-1.0
         xface(1)=xsten(0,1)
         xface(2)=xsten(-1,2)
         write(2,*) "side3"
        elseif((dir.eq.2).and.(side.eq.2)) then
         faceid=4
         ii=0
         jj=1
         nface(1)=0.0
         nface(2)=1.0
         xface(1)=xsten(0,1)
         xface(2)=xsten(1,2)
         write(2,*) "side4"
        else
         print *,"dir or side invalid"
         stop
        endif

        nf(1)=ii
        nf(2)=jj

!-------------------------------------------------------------------
        call two_points_dist(sdim,mat_cen_sten(ii,jj,im2,:), &
           mat_cen_sten(0,0,im1,:),Ltemp)

        if(dclt_test .eq. 1)then
         !-------------------------------------------------
         if(im1 .eq. im2)then
          !do nothing
         else 
          dclt_ratio = abs(xsten(ii+jj,dir) - &
                           mat_cen_sten(0,0,im1,dir)) / &
                       abs(mat_cen_sten(ii,jj,im2,dir)- &
                           mat_cen_sten(0,0,im1,dir)) 
          if(dclt_ratio .lt. 1.0e-8)then
           print *,"error, dclt_ratio is invalid"
           stop
          endif
          if (dclt_ratio.gt.1.0+1.0e-8) then
           print *,"dclt_ratio invalid"
           stop
          endif
          Ltemp = Ltemp*dclt_ratio   
         endif
         !------------------------------------------------ 
        else if (dclt_test.eq.0) then
         ! do nothing
        else
         print *,"dclt_test invalid"
         stop
        endif
!---------------------------------------------------------------------

        call get_kappa(dclt_test,alpha,xsten,dir,side, &
         mat_cen_sten(0,0,im1,:), &
         mat_cen_sten(ii,jj,im2,:), &
         im1,im2,kappa,nmat,sdim)

        if (Ltemp.le.0.0) then
         print *,"Ltemp invalid"
         stop
        endif
     
       if(dclt_test .eq. 1 .and. im1 .ne. im2) then
        do dir3 = 1,sdim
            dclt_diff(dir3) = mat_cen_sten(ii,jj,im2,dir3)- &
                              mat_cen_sten(0,0,im1,dir3) 
        enddo
        do dir3 = 1,SDIM
            dclt_pp(dir3) = mat_cen_sten(0,0,im1,dir3) + &
                            dclt_ratio*dclt_diff(dir3)
        enddo
        else
          dclt_pp = 0.0d0 
       endif


       if(dclt_test .eq. 0)then
         do dir2=1,sdim
          n1(dir2) = mat_cen_sten(ii,jj,im2,dir2)-mat_cen_sten(0,0,im1,dir2)
         enddo
       elseif(dclt_test .eq. 1)then
         do dir2=1,sdim
          n1(dir2) = mat_cen_sten(ii,jj,im2,dir2)-mat_cen_sten(0,0,im1,dir2)
         enddo
         if(im1 .eq. im2)then
          ! do nothing
         else 
!          nface_connect(1)=-n1(2)
!          nface_connect(2)=n1(1)
!          do dir2=1,sdim
!           xface_connect(dir2)=mat_cen_sten(0,0,im1,dir2)  
!          enddo
           ! n1 dot (x-x1)=0
           ! n2 dot (x-x2)=0
           ! n1 dot x=n1 dot x1
           ! n2 dot x=n2 dot x2
           ! x=x1 n1 + x2 n2
           
          do dir2=1,sdim
           n1(dir2) = dclt_pp(dir2)-mat_cen_sten(0,0,im1,dir2)
          enddo
         endif 
       else 
        print *,"dclt_test invalid"
        stop
       endif

        coef=kappa*abs(dot_product(nf,n1))/(Ltemp*Ltemp)

        write(2,*) "coef in diag ext", coef
        write(2,*) "AFRAC", AFRAC
        write(2,*) "dx", dx(dir)

        write(2,*) "diag before", diag
        diag=diag+coef*AFRAC*dx(dir)
        write(2,*) "diag after", diag
       elseif(AFRAC.gt.1.0d0)then
        print *,"frac_pair_cell bust in cell_diag_cal"
        stop
       elseif(AFRAC.lt.0.0d0) then
        print *,"frac_pair_cell invalid in cell_diag_cal"
        stop

       write(2,*) "in ext_face calVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV"
       endif

      enddo ! side
     enddo ! dir
    enddo ! im2

   write(2,*) "in ext_face cal^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^"


    ! internal flux cal
    ! first, check (im1, 1), (im1, 2), (im1, 3), ...

    do im2 = 1,nmat
     AFRAC=int_face(im1,im2)

     if(AFRAC.gt.eps*dx(1)) then
      write(2,*) "in interface calvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv"
      write(2,*) "im1", im1, "im2", im2, "A*FRAC", AFRAC 
      if (im1.eq.im2) then
       print *,"im1==im2 cannot happen"
       stop
      endif
       ! points from im1 to im2.
      do dir2=1,sdim
       nf(dir2)=int_face_normal(im1,im2,dir2)
       if(abs(nf(dir2)) .lt. 1.0e-8) then
         nf(dir2) = 0.0d0
       endif
      enddo

      call two_points_dist(sdim,mat_cen_sten(0,0,im2,:), &
                           mat_cen_sten(0,0,im1,:),Ltemp)
!-------------------------------------------------------
      if(dclt_test .eq. 1)then
       !---------------------------------------
       dclt_ratio1 = dist_to_int(im1,im2)/&
                    (dist_to_int(im1,im2) + dist_to_int(im2,im1))
       if ((dclt_ratio1.lt.0.0).or.(dclt_ratio1.gt.1.0)) then
        print *,"dclt_ratio1 invalid"
        stop
       endif
       if(dist_to_int(im1,im2) .lt. 10e-8)then
        print *,"error,dist_to_int(im1,im2) is 0"
        stop
       endif
       if(dist_to_int(im2,im1) .lt. 10e-8)then
        print *,"error,dist_to_int(im2,im1) is 0"
        stop
       endif
       Ltemp = Ltemp*dclt_ratio1
       !--------------------------------------
      else if (dclt_test.eq.0) then
       ! do nothing
      else
       print *,"dclt_test invalid"
       stop
      endif


        do dir3 = 1,sdim
          dclt_diff1(dir3) = mat_cen_sten(0,0,im2,dir3)- &
                              mat_cen_sten(0,0,im1,dir3) 
        enddo
        do dir3 = 1,SDIM
          dclt_pp1(dir3) = mat_cen_sten(0,0,im1,dir3) + &
                            dclt_ratio1*dclt_diff1(dir3)
        enddo 


!-------------------------------------------------------------
   if(dclt_test .eq. 0) then
      do dir2=1,sdim
       n1(dir2) = mat_cen_sten(0,0,im2,dir2)-mat_cen_sten(0,0,im1,dir2)
       if(abs(n1(dir2)) .lt. 1.0e-8 )then
         n1(dir2) = 0.0d0
       endif
      enddo
   elseif(dclt_test .eq. 1)then
      do dir2=1,sdim
       n1(dir2) = dclt_pp1(dir2)-mat_cen_sten(0,0,im1,dir2)
       if(abs(n1(dir2)) .lt. 1.0e-8 )then
         n1(dir2) = 0.0d0
       endif
      enddo     
   else
     print *,"error"
     stop
   endif
!------------------------------------------------------------
      if(dclt_test .eq. 0)then
       call get_kappa_int(alpha,dist_to_int,im1,im2,kappa,nmat)
      elseif(dclt_test .eq. 1)then
       kappa = alpha(im1)
      else
       print *,"invalid dclt_test flag"
       stop
      endif
 !------------------------------------------------------------
     coef=kappa*abs(dot_product(nf,n1))/(Ltemp*Ltemp)
   
        write(2,*) "nf", nf
        write(2,*) "n1", n1
        write(2,*) "nfXn1", abs(dot_product(nf,n1))
        write(2,*) "Ltemp", Ltemp
        write(2,*) "coef in diag int", coef
        write(2,*) "AFRAC", AFRAC
  
        write(2,*) "diag before", diag
     diag=diag+AFRAC*coef
        write(2,*) "diag after", diag
         

   write(2,*) "in interface cal^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^"
     endif
    enddo ! im2

  write(2,*)  "new diag cal loop end********************************************"
  write(2,*) "**************************************************************"

  end subroutine cell_diag_cal

  ! -div alpha grad u
  subroutine cell_div_cal( &
   ngeom_recon, &
   sdim,nmat,dx, &
   xsten, &
   mofdata_sten, &
   rho_box, &
   alpha, &
   mat_cen_sten,&
   im_in,  &
   frac_pair_cell, &
   int_face, &
   int_face_normal, &
   dist_to_int, &
   div_tot)

    implicit none 

    integer,parameter        :: dclt_flag = 0
    integer,intent(in)       :: sdim,nmat,ngeom_recon
    real(kind=8),intent(in)  :: rho_box(-1:1,-1:1,nmat)
    real(kind=8)             :: rho(-1:1,-1:1,nmat)
    real(kind=8),intent(in)  :: xsten(-3:3,sdim)

    real(kind=8),intent(in)  :: frac_pair_cell(nmat,nmat,sdim,2)
    real(kind=8),intent(in)  :: int_face(nmat,nmat)
    real(kind=8),intent(in)  :: int_face_normal(nmat,nmat,sdim)
    real(kind=8),intent(in)  :: mofdata_sten(-1:1,-1:1,ngeom_recon*nmat)
    real(kind=8),intent(in)  :: dx(sdim)
    real(kind=8),intent(in)  :: alpha(nmat)
    real(kind=8),intent(in)  :: dist_to_int(nmat,nmat)

    integer                  :: i,j,im,im1,im2
    integer                  :: dir,side
    integer                  :: dir2
    real(kind=8)             :: mat_cen_sten(-1:1,-1:1,nmat,sdim)
    real(kind=8)             :: dist_temp

    real(kind=8)             :: kappa
    real(kind=8)             :: grad
    real(kind=8),intent(out) :: div_tot
    integer                  :: success_flag


    !integer,intent(in)       :: function_flag
    integer,intent(in)       :: im_in
    real(kind=8)             :: Ltemp
    real(kind=8)             :: coef
    real(kind=8)             :: AFRAC,vf1,vf2
    real(kind=8)             :: LOWTOL,grad_high
    real(kind=8)             :: nf(sdim),n1(sdim)

    real(kind=8),external    :: norm_2d
    type(polygon)            :: cell(-1:1,-1:1)
    integer                  :: vofcomp1,vofcomp2,debughigh
    integer                  :: faceid,ii,jj

    debughigh=0

    LOWTOL=0.01d0

    if (ngeom_recon.ne.2*sdim+3) then
     print *,"ngeom_recon invalid"
     stop
    endif

    div_tot = 0.0d0
    rho = rho_box
    success_flag=1 ! cell_div_cal

    do i=-1,1
    do j=-1,1
     cell(i,j)%center%val(1)=xsten(2*i,1) 
     cell(i,j)%center%val(2)=xsten(2*j,2) 
    enddo
    enddo

    ! cell boundary(external flux) first

    if(im_in .le. 0) then
      print *,"wrong flag"
      stop
    elseif(im_in .gt. 0) then
     im1 = im_in

     vofcomp1=(im1-1)*ngeom_recon+1
     vf1=mofdata_sten(0,0,vofcomp1)
     if ((vf1.lt.-LOWTOL).or.(vf1.gt.one+LOWTOL)) then
      print *,"vf1 invalid"
      stop
     endif
        
     do im2 = 1,nmat
      do dir = 1,sdim
       do side = 1,2

        AFRAC=frac_pair_cell(im2,im1,dir,side)

        if (abs(AFRAC).le.eps) then
         AFRAC=0.0d0
        else if (abs(AFRAC-1.0d0).le.eps) then
         AFRAC=1.0d0
        else if ((AFRAC.gt.0.0).and.(AFRAC.lt.1.0)) then
         ! do nothing
        else
         print *,"AFRAC invalid"
         stop
        endif

        if ((AFRAC.gt.0.0d0).and.(AFRAC.le.1.0d0)) then

         if((dir.eq.1).and.(side.eq.1)) then
          faceid=1
          ii=-1
          jj=0
         elseif((dir.eq.1).and.(side.eq.2)) then
          faceid=2
          ii=1
          jj=0
         elseif((dir.eq.2).and.(side.eq.1)) then
          faceid=3
          ii=0
          jj=-1
         elseif((dir.eq.2).and.(side.eq.2)) then
          faceid=4
          ii=0
          jj=1
         else
          print *,"dir or side invalid"
          stop
         endif

         nf(1)=ii
         nf(2)=jj

          ! success_flag=0 => failure
          ! success_flag=1 => success
         call two_points_dist(sdim,mat_cen_sten(ii,jj,im2,:), &
          mat_cen_sten(0,0,im1,:),Ltemp)
         do dir2=1,sdim
          n1(dir2) = mat_cen_sten(ii,jj,im2,dir2)-mat_cen_sten(0,0,im1,dir2)
         enddo
         call get_kappa(dclt_flag,alpha,xsten,dir,side, &
          mat_cen_sten(0,0,im1,:), &
          mat_cen_sten(ii,jj,im2,:), &
          im1,im2,kappa,nmat,sdim)
         grad=rho(0,0,im1)-rho(ii,jj,im2)

         if (Ltemp.le.0.0) then
          print *,"Ltemp invalid"
          stop
         endif

         coef=kappa*abs(dot_product(nf,n1))/(Ltemp*Ltemp)
         grad=-grad*coef

         vofcomp2=(im2-1)*ngeom_recon+1
         vf2=mofdata_sten(ii,jj,vofcomp2)
         if ((vf2.lt.-LOWTOL).or.(vf2.gt.one+LOWTOL)) then
          print *,"vf2 invalid"
          stop
         endif
         if  ((vf1.le.1.0d0-LOWTOL).or.(vf2.le.1.0d0-LOWTOL)) then
          call ext_grad(nmat,sdim,faceid,alpha,im2,im1,rho, &
           cell,dx,mat_cen_sten,mofdata_sten, &
           grad_high,success_flag)
          if (success_flag.eq.1) then
           if (debughigh.eq.1) then
            print *,"high order dir,side,x,y,im1,im2 ",dir,side, &
             xsten(0,1),xsten(0,2),im1,im2
           endif
           grad=grad_high
          else if (success_flag.eq.0) then
           ! do nothing
          else
           print *,"success_flag invalid"
           stop
          endif
         else if ((abs(vf1-1.0d0).le.LOWTOL).and. &
                  (abs(vf2-1.0d0).le.LOWTOL)) then
          ! do nothing
         else
          print *,"vf1 or vf2 invalid"
        stop
         endif

          ! div_tot will be an approximation to 
          ! -integral_material_boundary alpha grad rho dot n dA
          ! where n points out.
         div_tot = div_tot-grad*AFRAC*dx(dir)

        elseif(AFRAC.gt.1.0d0)then
         print *,"frac_pair_cell bust in cell_div_cal"
         stop
        elseif(AFRAC.lt.0.0d0) then
         print *,"frac_pair_cell invalid in cell_div_cal"
         stop
        endif

       enddo ! side
      enddo ! dir
     enddo ! im2

      ! internal flux cal
     do im2 = 1,nmat
      AFRAC=int_face(im1,im2)
      if(AFRAC.gt.eps*dx(1)) then
       if (im1.eq.im2) then
        print *,"im1==im2 cannot happen"
        stop
       endif

       ! points from im1 to im2.
       do dir2=1,sdim
        nf(dir2)=int_face_normal(im1,im2,dir2)
       enddo
       call two_points_dist(sdim,mat_cen_sten(0,0,im2,:), &
                            mat_cen_sten(0,0,im1,:),Ltemp)
       do dir2=1,sdim
        n1(dir2) = mat_cen_sten(0,0,im2,dir2)-mat_cen_sten(0,0,im1,dir2)
       enddo
       call get_kappa_int(alpha,dist_to_int,im1,im2,kappa,nmat)

       coef=kappa*abs(dot_product(nf,n1))/(Ltemp*Ltemp)
       grad=-coef*(rho(0,0,im1)-rho(0,0,im2))

       call int_grad(nmat,sdim,alpha,im2,im1,rho,&
        cell,dx,mat_cen_sten,mofdata_sten, &
        grad_high,success_flag)
       if (success_flag.eq.1) then
        if (debughigh.eq.1) then
         print *,"high order interior,x,y,im1,im2 ", &
           xsten(0,1),xsten(0,2),im1,im2
        endif
        grad=grad_high
       else if (success_flag.eq.0) then
        ! do nothing
       else
        print *,"success_flag invalid"
        stop
       endif
 
       div_tot = div_tot - grad*AFRAC
      endif
     enddo ! im2
       
    else
     print *,"material_in invalid" 
     stop
    endif

  end subroutine cell_div_cal

   ! div_tot = -div alpha grad u
  subroutine cell_div_cal_simple( &
   probtype, &
   dclt_test, &
   T_in, &
   ngeom_recon, &
   sdim,nmat,dx, &
   xsten, &
   mofdata_sten, &
   rho_box, &
   alpha, &
   mat_cen_sten,&
   im_in,  &
   frac_pair_cell, &
   int_face, &
   int_face_normal, &
   dist_to_int, &
   div_tot,rhs_loc)

    implicit none 
   
    integer,intent(in)       :: probtype
    integer,intent(in)       :: dclt_test
    real(kind=8)             :: T_in

    integer,intent(in)       :: sdim,nmat,ngeom_recon
    real(kind=8),intent(in)  :: rho_box(-1:1,-1:1,nmat)
    real(kind=8)             :: rho(-1:1,-1:1,nmat)
    real(kind=8),intent(in)  :: xsten(-3:3,sdim)

    real(kind=8),intent(in)  :: frac_pair_cell(nmat,nmat,sdim,2)
    real(kind=8),intent(in)  :: int_face(nmat,nmat)
    real(kind=8),intent(in)  :: int_face_normal(nmat,nmat,sdim)
    real(kind=8),intent(in)  :: mofdata_sten(-1:1,-1:1,ngeom_recon*nmat)
    real(kind=8),intent(in)  :: dx(sdim)
    real(kind=8),intent(in)  :: alpha(nmat)
    real(kind=8),intent(in)  :: dist_to_int(nmat,nmat)

    integer                  :: i,j,im,im1,im2
    integer                  :: dir,side
    integer                  :: dir2,dir3
    real(kind=8)             :: mat_cen_sten(-1:1,-1:1,nmat,sdim)
    real(kind=8)             :: dist_temp

    real(kind=8)             :: kappa,AFRAC,vf1,vf2
    real(kind=8)             :: grad
    real(kind=8),intent(out) :: div_tot

    !integer,intent(in)       :: function_flag
    integer,intent(in)       :: im_in
    real(kind=8)             :: Ltemp
    real(kind=8)             :: coef,LOWTOL
    real(kind=8)             :: nf(sdim),n1(sdim)
    integer                  :: vofcomp1,vofcomp2,debughigh
    integer                  :: faceid,ii,jj

    real(kind=8),external    :: norm_2d
    type(polygon)            :: cell(-1:1,-1:1)
    !--------------------------------------------
    real(kind=8)             :: dclt_ratio,dclt_ratio1
    real(kind=8)             :: dclt_diff(SDIM),dclt_diff1(SDIM)
    real(kind=8)             :: dclt_pp(SDIM),dclt_pp1(SDIM)
    real(kind=8),external    :: exact_temperature
   
    real(kind=8)             :: rhstemp1,rhstemp2,rhs_loc

    LOWTOL=0.01d0

  write(2,*)  "new div_cal_simple loop*************************************"
  write(2,*) "**************************************************************"
  write(2,*) "considering mat", im_in

    if (ngeom_recon.ne.2*sdim+3) then
     print *,"ngeom_recon invalid"
     stop
    endif
    if ((dclt_test.ne.0).and.(dclt_test.ne.1)) then
     print *,"dclt_test invalid"
     stop
    endif 



   rhs_loc = 0.0d0
   div_tot = 0.0d0
    write(2,*) "T profile"
    write(2,*) rho_box(0,0,1),rho_box(0,0,2)
    rho = rho_box

    do i=-1,1
    do j=-1,1
     cell(i,j)%center%val(1)=xsten(2*i,1) 
     cell(i,j)%center%val(2)=xsten(2*j,2) 
    enddo
    enddo

    ! cell boundary(external flux) first

    if(im_in .le. 0) then
      print *,"wrong flag"
      stop
    elseif(im_in .gt. 0) then
     im1 = im_in
    
     vofcomp1=(im1-1)*ngeom_recon+1
     vf1=mofdata_sten(0,0,vofcomp1)
     if ((vf1.lt.-LOWTOL).or.(vf1.gt.one+LOWTOL)) then
      print *,"vf1 invalid"
      stop
     endif
    
     do im2 = 1,nmat
      do dir = 1,sdim
       do side = 1,2
        
         rhstemp1 = 0.0d0
        AFRAC=frac_pair_cell(im2,im1,dir,side)

        if (abs(AFRAC).le.eps) then
         AFRAC=0.0d0
        else if (abs(AFRAC-1.0d0).le.eps) then
         AFRAC=1.0d0
        else if ((AFRAC.gt.0.0).and.(AFRAC.lt.1.0)) then
         ! do nothing
        else
         print *,"AFRAC invalid"
         stop
        endif

        if ((AFRAC.gt.0.0d0).and.(AFRAC.le.1.0d0)) then
        write(2,*) "in ext_face calvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv"

         if((dir.eq.1).and.(side.eq.1)) then
          faceid=1
          ii=-1
          jj=0
          write(2,*) "side1"
         elseif((dir.eq.1).and.(side.eq.2)) then
          faceid=2
          ii=1
          jj=0
          write(2,*) "side2"
         elseif((dir.eq.2).and.(side.eq.1)) then
          faceid=3
          ii=0
          jj=-1
          write(2,*) "side3"
         elseif((dir.eq.2).and.(side.eq.2)) then
          faceid=4
          ii=0
          jj=1
          write(2,*) "side4"
         else
          print *,"dir or side invalid"
          stop
         endif

         nf(1)=ii
         nf(2)=jj

         call two_points_dist(sdim,mat_cen_sten(ii,jj,im2,:), &
          mat_cen_sten(0,0,im1,:),Ltemp)
        !--------------------------------------------------
        if(dclt_test .eq. 1)then
         !-------------------------------------------------
         if(im1 .eq. im2)then
          !do nothing
         else 
          dclt_ratio = abs(xsten(ii+jj,dir) - &
                       mat_cen_sten(0,0,im1,dir)) / &
                       abs(mat_cen_sten(ii,jj,im2,dir)- &
                           mat_cen_sten(0,0,im1,dir)) 
          if(dclt_ratio .lt. 10e-8)then
           print *,"error, dclt_ratio is invalid"
           stop
          endif
          if (dclt_ratio.gt.1.0+1.0e-8) then
           print *,"dclt_ratio invalid"
           stop
          endif
          Ltemp = Ltemp*dclt_ratio   
         endif
         !------------------------------------------------ 
        else if (dclt_test.eq.0) then
         ! do nothing
        else
         print *,"dclt_test invalid"
         stop
        endif


         call get_kappa(dclt_test,alpha,xsten,dir,side, &
          mat_cen_sten(0,0,im1,:), &
          mat_cen_sten(ii,jj,im2,:), &
          im1,im2,kappa,nmat,sdim)

        !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
         grad=rho(0,0,im1)-rho(ii,jj,im2)
         write(2,*) "rho(0 0 im1)", rho(0,0,im1),"rho(ii,jj,im2)",rho(ii,jj,im2)
         if(dclt_test .eq. 1)then
     !       write(2,*) "WTF"

          if(im1 .eq. im2)then
     !       write(2,*) "WTF1"
            ! do nothing
          else
     !      write(2,*) "WTF2"
           if(dclt_ratio .lt. 10e-8)then
            print *,"error, dclt_ratio is invalid"
            stop
           endif
           do dir3 = 1,sdim
            dclt_diff(dir3) = mat_cen_sten(ii,jj,im2,dir3)- &
                              mat_cen_sten(0,0,im1,dir3) 
           enddo
           do dir3 = 1,SDIM
            dclt_pp(dir3) = mat_cen_sten(0,0,im1,dir3) + &
                            dclt_ratio*dclt_diff(dir3)
           enddo 
           grad = rho(0,0,im1) 
           rhstemp1 = exact_temperature(dclt_pp(1),dclt_pp(2), &
                       T_in,im1,probtype, nmat, alpha,dclt_test)

         write(2,*) "rho(0 0 im1)", rho(0,0,im1)
         write(2,*) "EXACT",exact_temperature(dclt_pp(1),dclt_pp(2), &
                       T_in,im1,probtype, nmat, alpha,dclt_test)
 
          endif
         endif

        write(2,*) "grad1",grad


         if (Ltemp.le.0.0) then
          print *,"Ltemp invalid"
          stop
         endif


       if(dclt_test .eq. 0)then
         do dir2=1,sdim
          n1(dir2) = mat_cen_sten(ii,jj,im2,dir2)-mat_cen_sten(0,0,im1,dir2)
         enddo
       elseif(dclt_test .eq. 1)then
         if(im1 .eq. im2)then
          do dir2=1,sdim
          n1(dir2) = mat_cen_sten(ii,jj,im2,dir2)-mat_cen_sten(0,0,im1,dir2)
         enddo
         else      
         do dir2=1,sdim
          n1(dir2) = dclt_pp(dir2)-mat_cen_sten(0,0,im1,dir2)
         enddo
         endif        
       else 
        print *,"dclt_test invalid"
        stop
       endif

       !write(2,*) "Ltemp", Ltemp
       !write(2,*) "nf", nf
       !write(2,*) "n1", n1
       !write(2,*) "kappa", kappa



         coef=kappa*abs(dot_product(nf,n1))/(Ltemp*Ltemp)
         grad=grad*coef
         if(dclt_test .eq. 1)then
          rhstemp1 = rhstemp1*coef
         endif
         vofcomp2=(im2-1)*ngeom_recon+1
         vf2=mofdata_sten(ii,jj,vofcomp2)
         if ((vf2.lt.-LOWTOL).or.(vf2.gt.one+LOWTOL)) then
          print *,"vf2 invalid"
          stop
         endif


        write(2,*) "coef in div ext", coef
        write(2,*) "AFRAC", AFRAC
        write(2,*) "dx", dx(dir)
        write(2,*) "grad2",grad
        write(2,*) "div before", div_tot
 

         div_tot = div_tot+grad*AFRAC*dx(dir)
         rhs_loc = rhs_loc + rhstemp1*AFRAC*dx(dir)
        write(2,*) "div after", div_tot

         !write(2,*) "div_tot of ext", div_tot

        elseif(AFRAC.gt.1.0d0)then
         print *,"frac_pair_cell bust in cell_div_cal"
         stop
        elseif(AFRAC.lt.0.0d0) then
         print *,"frac_pair_cell invalid in cell_div_cal"
         stop
        endif

       enddo ! side
      enddo ! dir
     enddo ! im2


! internal interface
     do im2 = 1,nmat
      AFRAC=int_face(im1,im2)
      rhstemp2 = 0.0d0
      if (1.eq.0) then
       if ((vf1.gt.0.0).and.(vf1.lt.1.0)) then
        print *,"x,y,im1,vf ",xsten(0,1),xsten(0,2),im1,vf1
       endif
      endif

      if(AFRAC.gt.eps*dx(1)) then
      write(2,*) "in interface calvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv"
      write(2,*) "im1", im1, "im2", im2, "AFRAC", AFRAC 
       if (im1.eq.im2) then
        print *,"im1==im2 cannot happen"
        stop
       endif
       ! points from im1 to im2.
       do dir2=1,sdim
        nf(dir2)=int_face_normal(im1,im2,dir2)
       enddo
       call two_points_dist(sdim,mat_cen_sten(0,0,im2,:), &
                            mat_cen_sten(0,0,im1,:),Ltemp)
!-------------------------------------------------------
      if(dclt_test .eq. 1)then
       !---------------------------------------
       dclt_ratio1 = dist_to_int(im1,im2)/&
                    (dist_to_int(im1,im2) + dist_to_int(im2,im1))
       if ((dclt_ratio1.lt.0.0).or.(dclt_ratio1.gt.1.0)) then
        print *,"dclt_ratio1 invalid"
        stop
       endif
       !write(2,*) "int_ratio", dclt_ratio1
       if(dist_to_int(im1,im2) .lt. 10e-8)then
        print *,"error,dist_to_int(im1,im2) is 0"
        stop
       endif
       if(dist_to_int(im2,im1) .lt. 10e-8)then
        print *,"error,dist_to_int(im2,im1) is 0"
        stop
       endif
       Ltemp = Ltemp*dclt_ratio1
   
       if(Ltemp .le. 0.0d0) then
        print *, "Ltemp invalid"
        stop
       endif 


       !--------------------------------------
        do dir3 = 1,sdim
          dclt_diff1(dir3) = mat_cen_sten(0,0,im2,dir3)- &
                              mat_cen_sten(0,0,im1,dir3) 
        enddo
        do dir3 = 1,SDIM
          dclt_pp1(dir3) = mat_cen_sten(0,0,im1,dir3) + &
                            dclt_ratio1*dclt_diff1(dir3)
        enddo 

       !write(2,*) "dclt_pp1", dclt_pp1 

      else if (dclt_test.eq.0) then
       ! do nothing
      else
       print *,"dclt_test invalid"
       stop

      endif

   if(dclt_test .eq. 0) then
      do dir2=1,sdim
       n1(dir2) = mat_cen_sten(0,0,im2,dir2)-mat_cen_sten(0,0,im1,dir2)
      enddo
   elseif(dclt_test .eq. 1)then
      do dir2=1,sdim
       n1(dir2) = dclt_pp1(dir2)-mat_cen_sten(0,0,im1,dir2)
      enddo     
   else
     print *,"error"
     stop
   endif

!------------------------------------------------------------
      if(dclt_test .eq. 0)then
       call get_kappa_int(alpha,dist_to_int,im1,im2,kappa,nmat)
      elseif(dclt_test .eq. 1)then
       kappa = alpha(im1)
      else
       print *,"invalid dclt_test flag"
       stop
      endif
 !------------------------------------------------------------
       coef=kappa*abs(dot_product(nf,n1))/(Ltemp*Ltemp) 

        write(2,*) "nf", nf
        write(2,*) "n1", n1
        write(2,*) "nfXn1", abs(dot_product(nf,n1))
        write(2,*) "Ltemp", Ltemp
        write(2,*) "coef in div int", coef
        write(2,*) "AFRAC", AFRAC


       !write(2,*) "rho(0,0,im1)", rho(0,0,im1)
       !write(2,*) "point", dclt_pp1, "current time", T_in, "probtype",probtype
       !write(2,*) "im1",im1, "probtype", probtype, "alpha", alpha, "flag", dclt_test
       !write(2,*) "exact", exact_temperature(dclt_pp1(1),dclt_pp1(2),T_in,&
        !          im1,probtype,nmat,alpha,dclt_test)

      if(dclt_test .eq. 0)then
       grad=coef*(rho(0,0,im1)-rho(0,0,im2))

       !write(2,*) "time", T_in, "coef",coef
       write(2,*) "Tim1", rho(0,0,im1), "Tim2", rho(0,0,im2)
       write(2,*) "grad",grad

      elseif(dclt_test .eq. 1)then

!       grad= coef*(rho(0,0,im1) - &
!                  exact_temperature(dclt_pp1(1),dclt_pp1(2),T_in,&
!                  im1,probtype,nmat,alpha,dclt_test))
       grad= coef*rho(0,0,im1)
       rhstemp2= coef*exact_temperature(dclt_pp1(1),dclt_pp1(2),T_in,&
                  im1,probtype,nmat,alpha,dclt_test)
       !write(2,*) "time", T_in, "coef",coef
       write(2,*) "Tim1", rho(0,0,im1), "Tdclt",  &
        exact_temperature(dclt_pp1(1),dclt_pp1(2),T_in,&
                 im1,probtype,nmat,alpha,dclt_test)     

       write(2,*) "grad",grad

      else
       print *,"invalid dclt_test flag"
       stop
      endif

!----------------------------------------
!      print *, "grad", grad, "Afrac", AFRAC
        write(2,*) "div before", div_tot


       div_tot = div_tot + grad*AFRAC
       rhs_loc = rhs_loc + rhstemp2*AFRAC

        write(2,*) "div after", div_tot
 
      !write(2,*) "div_tot of int", div_tot
   write(2,*) "in interface cal^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^"


      endif
     enddo ! im2
    else
     print *,"material_in invalid" 
     stop
    endif
  write(2,*)  "new div_simple cal loop end********************************************"
  write(2,*) "**************************************************************"


  end subroutine cell_div_cal_simple



  !--------------------------------------------------------------
  !////////////////////////////////////////////////////////////////////
  !////////////////////////////////////////////////////////////////////
  subroutine ext_grad(nmat,sdim,flag,kappa,mat_osd,mat_isd, T,&
         cell, dx,mat_cen_sten,mofdata, &
         grad,success_flag)
    implicit none

    integer,intent(in)            :: nmat,sdim 
    real(kind=8),intent(in)       :: mat_cen_sten(-1:1,-1:1,nmat,sdim)
    real(kind=8),intent(in)       :: kappa(nmat)

    type(polygon)                 :: cell(-1:1, -1:1)
    real(kind=8),intent(in)       :: dx(sdim)
    integer,intent(in)            :: flag
    integer,intent(in)            :: mat_osd,mat_isd
    ! material number
    real(kind=8),intent(in)       :: mofdata(-1:1,-1:1,(sdim*2+3)*nmat)
    real(kind=8),dimension(-1:1,-1:1,nmat),intent(in) :: T
    integer :: success_flag

    integer                       :: i,j,k
    integer                       :: datalen
    integer                       :: tmat

    integer,dimension(-1:1,-1:1)           :: ttag
    integer,dimension(-1:1,-1:1)           :: osd_tag,isd_tag
    real(kind=8),dimension(-1:1,-1:1,sdim) :: lsdata_osd, lsdata_isd,lsdata
    real(kind=8),dimension(-1:1,-1:1)      :: w_osd, w_isd,w
  
    real(kind=8)                  :: normal(sdim)
    type(points)                  :: fluxpoint
    type(points)                  :: x1,x2
    
    real(kind=8)                :: xdim
    real(kind=8),allocatable    :: A(:,:) , b(:), x_solved(:)
    real(kind=8), intent(out)   :: grad
    real(kind=8) ,external      :: norm_2d

    success_flag=1 ! ext_grad

    tmat = 0
    datalen = (sdim * 2 + 3)
    osd_tag = 0
    isd_tag = 0
    ttag = 0
    lsdata_isd = 0.0d0
    lsdata_osd = 0.0d0
    lsdata = 0.0d0
    normal = 0.0d0
    w_osd = 0.0d0
    w_isd = 0.0d0
    w = 0.0d0

!     print *,"mofdata", mofdata(0,0,1:7)
!     print *,"mofdata",mofdata(0,0,8:14)

    if (1.eq.0) then
     print *,"mat_cen_sten"

     do i = -1,1
     do j = -1,1
      print *,"i=",i,"j=",j
      print *, mat_cen_sten(i,j,1,:)
     enddo
     enddo
    endif

     if (flag .eq. 1) then

          normal(1) = -1
          normal(2) = 0

     elseif (flag .eq. 2) then

          normal(1) = 1
          normal(2) = 0

     elseif(flag .eq. 3)then

          normal(1) = 0
          normal(2) = -1

     elseif(flag .eq. 4)then

          normal(1) = 0
          normal(2) = 1
     else
          print *,"invalid flag in  ext_grad"
     endif


    if (mat_osd .eq. mat_isd) then

     if (flag .eq. 1) then                       
      do j = -1, 1                  
       if (mofdata(-1, j, datalen * (mat_osd - 1) + 1) .gt. eps) then
        osd_tag(-1, j) = 1
        lsdata_osd(-1, j, :) = mat_cen_sten(-1, j, mat_osd, :)
       endif
       if (mofdata(0, j, datalen * (mat_isd - 1) + 1) .gt. eps) then
        isd_tag(0, j) = 1
        lsdata_isd(0, j, :) = mat_cen_sten(0, j, mat_isd, :)
       endif
      enddo

     elseif (flag .eq. 2) then
      do j = -1, 1     
       if (mofdata(1, j, datalen * (mat_osd-1) + 1) .gt. eps) then
        osd_tag(1, j) = 1
        lsdata_osd(1, j, :) = mat_cen_sten(1, j, mat_osd, :)
       endif
       if(mofdata(0, j, datalen * (mat_isd - 1) + 1) .gt. eps) then
        isd_tag(0, j) = 1
        lsdata_isd(0, j, :) = mat_cen_sten(0, j, mat_isd, :)
       endif
      enddo

     elseif(flag .eq. 3)then
      do i = -1, 1
             if (mofdata(i, -1, datalen * (mat_osd-1) + 1) .gt. eps) then
                osd_tag(i, -1) = 1
                lsdata_osd(i, -1, :) = mat_cen_sten(i, -1, mat_osd, :)
             endif

             if(mofdata(i, 0, datalen * (mat_isd - 1) + 1) .gt. eps) then
                isd_tag(i, 0) = 1
                lsdata_isd(i, 0, :) = mat_cen_sten(i, 0, mat_isd, :)
             endif
          enddo


       elseif(flag .eq. 4)then

          do i = -1, 1     

             if (mofdata(i, 1, datalen * (mat_osd-1) + 1) .gt. eps) then
                osd_tag(i, 1) = 1
                lsdata_osd(i, 1, :) = mat_cen_sten(i, 1, mat_osd, :)
             endif

             if(mofdata(i, 0, datalen * (mat_isd - 1) + 1) .gt. eps) then
                isd_tag(i, 0) = 1
                lsdata_isd(i, 0, :) = mat_cen_sten(i, 0, mat_isd, :)
             endif

          enddo

       else
        print *,"invalid flag in  ext_grad"
        stop
       endif

  elseif(mat_osd .ne. mat_isd  .and. abs(mat_osd - mat_isd).le. 10 ) then

       do i = -1 , 1
          do j = -1, 1     
           
             if (mofdata(i, j, datalen * (mat_osd - 1) + 1) .gt. eps) then
                osd_tag(i, j) = 1
                lsdata_osd(i, j, :) = mat_cen_sten(i, j, mat_osd, :)
             endif

             if (mofdata(i, j, datalen * (mat_isd - 1) + 1) .gt. eps) then
                isd_tag(i, j) = 1
                lsdata_isd(i, j, :) = mat_cen_sten(i, j, mat_isd, :)
             endif

          enddo
        enddo

  else
    print *,"mat_osd  or mat_isd number invalid"
    stop
  endif

! find flux point,  intersection of interface and straight line through two material

  if(flag .eq. 1) then
   x1%val(1) = lsdata_osd(-1,0,1) 
   x1%val(2) = lsdata_osd(-1,0,2) 
   x2%val(1) = lsdata_isd(0,0,1) 
   x2%val(2) = lsdata_isd(0,0,2) 
   call find_fluxpoint(-1.0d0,0.0d0,cell(0,0)%center%val(1)-0.5d0*dx(1),x1,x2,fluxpoint)
    
  elseif(flag .eq. 2) then
   x1%val(1) = lsdata_osd(1,0,1) 
   x1%val(2) = lsdata_osd(1,0,2) 
   x2%val(1) = lsdata_isd(0,0,1) 
   x2%val(2) = lsdata_isd(0,0,2) 
   call find_fluxpoint(-1.0d0,0.0d0,cell(0,0)%center%val(1)+0.5d0*dx(1),x1,x2,fluxpoint)

  elseif(flag .eq. 3) then
   x1%val(1) = lsdata_osd(0,-1,1) 
   x1%val(2) = lsdata_osd(0,-1,2) 
   x2%val(1) = lsdata_isd(0,0,1) 
   x2%val(2) = lsdata_isd(0,0,2) 
   call find_fluxpoint(0.0d0,-1.0d0,cell(0,0)%center%val(2)-0.5d0*dx(2),x1,x2,fluxpoint)

  elseif(flag .eq. 4) then 
   x1%val(1) = lsdata_osd(0,1,1) 
   x1%val(2) = lsdata_osd(0,1,2) 
   x2%val(1) = lsdata_isd(0,0,1) 
   x2%val(2) = lsdata_isd(0,0,2) 
   call find_fluxpoint(0.0d0,-1.0d0,cell(0,0)%center%val(2)+0.5d0*dx(2),x1,x2,fluxpoint)
  else
    print *,"invalid flag in  ext_grad"
    stop
  endif
  
! calculate the weight
if(mat_osd .eq. mat_isd) then
  if(flag .eq. 1)then
   do j = -1, 1
    if(osd_tag(-1,j) .eq. 1) then
      CALL SDD(dx(1)/2, (NORM_2d(lsdata_osd(-1,j,1),lsdata_osd(-1,j,2), &
       lsdata_osd(-1,0,1),lsdata_osd(-1,0,2)))**2.0d0/dx(1), w(-1,j))
    elseif(osd_tag(-1,j) .eq. 0)then
      ! do nothing
    else
     print *,"osd_tag invalid when flag = 1"
     stop
    endif
    if(isd_tag(0,j) .eq. 1) then
     CALL SDD(dx(1)/2, (NORM_2d(lsdata_isd(0,j,1),lsdata_isd(0,j,2), &
       lsdata_isd(0,0,1),lsdata_isd(0,0,2)))**2.0d0/dx(1), w(0,j))
    elseif(isd_tag(0,j) .eq. 0)then
      ! do nothing
    else
     print *,"isd_tag invalid when flag = 1"
     stop
    endif
   enddo
  endif
  if(flag .eq. 2)then
   do j = -1, 1
    if(osd_tag(1,j) .eq. 1) then
      CALL SDD(dx(1)/2, (NORM_2d(lsdata_osd(1,j,1),lsdata_osd(1,j,2), &
       lsdata_osd(1,0,1),lsdata_osd(1,0,2)))**2.0d0/dx(1), w(1,j))
    elseif(osd_tag(1,j) .eq. 0)then
      ! do nothing
    else
     print *,"osd_tag invalid when flag = 2"
     stop
    endif
    if(isd_tag(0,j) .eq. 1) then
     CALL SDD(dx(1)/2, (NORM_2d(lsdata_isd(0,j,1),lsdata_isd(0,j,2), &
       lsdata_isd(0,0,1),lsdata_isd(0,0,2)))**2.0d0/dx(1), w(0,j))
    elseif(isd_tag(0,j) .eq. 0)then
      ! do nothing
    else
     print *,"isd_tag invalid when flag = 2"
     stop
    endif
   enddo
  endif

  if(flag .eq. 3)then
   do i = -1, 1
    if(osd_tag(i,-1) .eq. 1) then
      CALL SDD(dx(2)/2, (NORM_2d(lsdata_osd(i,-1,1),lsdata_osd(i,-1,2), &
       lsdata_osd(0,-1,1),lsdata_osd(0,-1,2)))**2.0d0/dx(2), w(i,-1))
    elseif(osd_tag(i,-1) .eq. 0)then
      ! do nothing
    else
     print *,"osd_tag invalid when flag = 3"
     stop
    endif
    if(isd_tag(i,0) .eq. 1) then
     CALL SDD(dx(2)/2, (NORM_2d(lsdata_isd(i,0,1),lsdata_isd(i,0,2), &
       lsdata_isd(0,0,1),lsdata_isd(0,0,2)))**2.0d0/dx(2), w(i,0))
    elseif(isd_tag(i,0) .eq. 0)then
      ! do nothing
    else
     print *,"isd_tag invalid when flag = 3"
     stop
    endif
   enddo
  endif
  if(flag .eq. 4)then
   do i = -1, 1
    if(osd_tag(i,1) .eq. 1) then
      CALL SDD(dx(2)/2, (NORM_2d(lsdata_osd(i,1,1),lsdata_osd(i,1,2), &
       lsdata_osd(0,1,1),lsdata_osd(0,1,2)))**2.0d0/dx(2), w(i,1))
    elseif(osd_tag(i,1) .eq. 0)then
      ! do nothing
    else
     print *,"osd_tag invalid when flag = 4"
     stop
    endif
    if(isd_tag(i,0) .eq. 1) then
     CALL SDD(dx(2)/2, (NORM_2d(lsdata_isd(i,0,1),lsdata_isd(i,0,2), &
       lsdata_isd(0,0,1),lsdata_isd(0,0,2)))**2.0d0/dx(2), w(i,0))
    elseif(isd_tag(i,0) .eq. 0)then
      ! do nothing
    else
     print *,"isd_tag invalid when flag = 4"
     stop
    endif
   enddo
  endif



  tmat = mat_isd
  do i = -1,1
  do j = -1,1
  do k = 1,sdim
   lsdata(i,j,k) = lsdata_isd(i,j,k) + lsdata_osd(i,j,k) 
  enddo
   ttag(i,j) = osd_tag(i,j) + isd_tag(i,j)
  enddo
  enddo

if(1 .eq. 0) then
  do i = -1 , 1
   do j = -1, 1
        print *,"weight:"
        print *,w(i,j)
        print *, "centroids:"
        print *,lsdata(i,j,:)
   enddo
  enddo
  
  print *,"fluxpoint:"
  print *, fluxpoint%val
endif

   
 if( kappa(tmat).le. eps .and. kappa(tmat) .ge. 0.0d0) then
  grad = 0.0d0
 elseif(kappa(tmat) .gt. eps) then 
  allocate(A(3,3),b(3),x_solved(3))
  call Matrix_assemble_single(sdim, nmat, dx, tmat, ttag, &
                        & lsdata, w, T, &
                        & normal, fluxpoint, kappa(tmat),&
                        & A,b)
 if(1 .eq. 0)then
  print *,"matrix_A"
  do i = 1,3
   PRINT *,A(:,i)
  enddo
  print *,"b"
  print *, b
 endif  
  
  call solve_determinant(3,A,b,x_solved,success_flag)
  grad = x_solved(2)
!  print *,"grad=", grad
  deallocate(A,b,x_solved) 
 else
  print *,"1. kappa is valid in ext_grad"
  print *,"kappa=", kappa(tmat)
  stop 
endif
  
else
  if(flag .eq. 1)then
   do j = -1, 1
    if(osd_tag(-1,j) .eq. 1) then
      CALL SDD(dx(1)/2, (NORM_2d(lsdata_osd(-1,j,1),lsdata_osd(-1,j,2), &
       lsdata_osd(-1,0,1),lsdata_osd(-1,0,2)))**2.0d0/dx(1), w_osd(-1,j))
    elseif(osd_tag(-1,j) .eq. 0)then
      ! do nothing
    else
     print *,"osd_tag invalid when flag = 1"
     stop
    endif
    if(isd_tag(0,j) .eq. 1) then
     CALL SDD(dx(1)/2, (NORM_2d(lsdata_isd(0,j,1),lsdata_isd(0,j,2), &
       lsdata_isd(0,0,1),lsdata_isd(0,0,2)))**2.0d0/dx(1), w_isd(0,j))
    elseif(isd_tag(0,j) .eq. 0)then
      ! do nothing
    else
     print *,"isd_tag invalid when flag = 1"
     stop
    endif
   enddo
  endif
  if(flag .eq. 2)then
   do j = -1, 1
    if(osd_tag(1,j) .eq. 1) then
      CALL SDD(dx(1)/2, (NORM_2d(lsdata_osd(1,j,1),lsdata_osd(1,j,2), &
       lsdata_osd(1,0,1),lsdata_osd(1,0,2)))**2.0d0/dx(1), w_osd(1,j))
    elseif(osd_tag(1,j) .eq. 0)then
      ! do nothing
    else
     print *,"osd_tag invalid when flag = 2"
     stop
    endif
    if(isd_tag(0,j) .eq. 1) then
     CALL SDD(dx(1)/2, (NORM_2d(lsdata_isd(0,j,1),lsdata_isd(0,j,2), &
       lsdata_isd(0,0,1),lsdata_isd(0,0,2)))**2.0d0/dx(1), w_isd(0,j))
    elseif(isd_tag(0,j) .eq. 0)then
      ! do nothing
    else
     print *,"isd_tag invalid when flag = 2"
     stop
    endif
   enddo
  endif

  if(flag .eq. 3)then
   do i = -1, 1
    if(osd_tag(i,-1) .eq. 1) then
      CALL SDD(dx(2)/2, (NORM_2d(lsdata_osd(i,-1,1),lsdata_osd(i,-1,2), &
       lsdata_osd(0,-1,1),lsdata_osd(0,-1,2)))**2.0d0/dx(2), w_osd(i,-1))
    elseif(osd_tag(i,-1) .eq. 0)then
      ! do nothing
    else
     print *,"osd_tag invalid when flag = 3"
     stop
    endif
    if(isd_tag(i,0) .eq. 1) then
     CALL SDD(dx(2)/2, (NORM_2d(lsdata_isd(i,0,1),lsdata_isd(i,0,2), &
       lsdata_isd(0,0,1),lsdata_isd(0,0,2)))**2.0d0/dx(2), w_isd(i,0))
    elseif(isd_tag(i,0) .eq. 0)then
      ! do nothing
    else
     print *,"isd_tag invalid when flag = 3"
     stop
    endif
   enddo
  endif
  if(flag .eq. 4)then
   do i = -1, 1
    if(osd_tag(i,1) .eq. 1) then
      CALL SDD(dx(2)/2, (NORM_2d(lsdata_osd(i,1,1),lsdata_osd(i,1,2), &
       lsdata_osd(0,1,1),lsdata_osd(0,1,2)))**2.0d0/dx(2), w_osd(i,1))
    elseif(osd_tag(i,1) .eq. 0)then
      ! do nothing
    else
     print *,"osd_tag invalid when flag = 4"
     stop
    endif
    if(isd_tag(i,0) .eq. 1) then
     CALL SDD(dx(2)/2, (NORM_2d(lsdata_isd(i,0,1),lsdata_isd(i,0,2), &
       lsdata_isd(0,0,1),lsdata_isd(0,0,2)))**2.0d0/dx(2), w_isd(i,0))
    elseif(isd_tag(i,0) .eq. 0)then
      ! do nothing
    else
     print *,"isd_tag invalid when flag = 4"
     stop
    endif
   enddo
  endif



! assemble the matrix
 if(kappa(mat_osd) .gt. eps  .and. kappa(mat_isd) .gt. eps) then
  allocate(A(4,4),b(4),x_solved(4))
  call Matrix_assemble_double(sdim, nmat, dx, mat_osd,mat_isd,osd_tag, &
    isd_tag, &
    lsdata_osd, lsdata_isd, w_osd, w_isd, T, &
    normal, fluxpoint, kappa(mat_osd), kappa(mat_isd), &
    A,b)

! use cramer's rule to solve for the flux  
  call solve_determinant(4,A,b,x_solved,success_flag)
  grad = x_solved(3)
   deallocate(A,b,x_solved)
 elseif(kappa(mat_osd) .ge. 0.0d0  .or. kappa(mat_isd) .ge. 0.0d0)then
  grad = 0.0d0
 else 
  print *,"2: kappa is invalid in ext_grad"
  print *, "kappa(mat_osd)=",kappa(mat_osd),"kappa(mat_isd)=",kappa(mat_isd)
  stop
 endif
endif

  end subroutine ext_grad

!-----------------------------------------------------------------------------
subroutine int_grad(nmat,sdim,kappa,mat_osd,mat_isd, T,&
                    & cell, dx,mat_cen_sten,mofdata, &
                    & grad,success_flag)
    implicit none

    integer                       :: success_flag
    integer,intent(in)            :: nmat,sdim
    real(kind=8),intent(in)       :: mat_cen_sten(-1:1,-1:1,nmat,sdim)
    real(kind=8),intent(in)       :: kappa(nmat)
    type(polygon)                 :: cell(-1:1,-1:1)
    real(kind=8)                  :: center(2)
    real(kind=8),intent(in)       :: dx(sdim)
    integer,intent(in)            :: mat_osd,mat_isd


    real(kind=8),intent(in)       :: mofdata(-1:1,-1:1,(sdim*2+3)*nmat)
    real(kind=8),dimension(-1:1,-1:1,nmat),intent(in) :: T

    integer                       :: i,j
    integer                       :: datalen

    integer,dimension(-1:1,-1:1)           :: osd_tag,isd_tag
    real(kind=8),dimension(-1:1,-1:1,sdim) :: lsdata_osd, lsdata_isd
    real(kind=8),dimension(-1:1,-1:1)      :: w_osd, w_isd

    real(kind=8)                  :: normal(sdim),centroid(sdim), centroid1(sdim)
    real(kind=8)                  :: c
    type(points)                  :: fluxpoint
    type(points)                  :: x1,x2

    real(kind=8),external         :: norm_2d

    real(kind=8)                :: A(4,4) , b(4), x_solved(4)
    real(kind=8), intent(out)   :: grad

    success_flag=1 ! int_grad

    datalen = (sdim * 2 + 3)
    osd_tag = 0
    isd_tag = 0
    lsdata_isd = 0.0d0
    lsdata_osd = 0.0d0
    normal = 0.0d0
    w_osd = 0.0d0
    w_isd = 0.0d0
    c = 0.0d0
    centroid = 0.0d0
    centroid1 = 0.0d0

    center(1) = cell(0,0)%center%val(1)
    center(2) = cell(0,0)%center%val(2)

  
    do i  = 1, sdim
      normal(i) = mofdata(0,0,datalen*(mat_isd - 1)+sdim+2+i)
      centroid(i) = mofdata(0,0,datalen*(mat_isd - 1) + 1+i)
      centroid(i) = centroid(i) + center(i)
      centroid1(i) = mofdata(0,0,datalen*(mat_osd - 1) + 1+i)
      centroid1(i) = centroid1(i) + center(i)
    enddo


    

    if(abs(normal(1)) .lt. eps) then
       normal(1)= 0.0d0
    endif
    if(abs(normal(2)) .lt. eps) then
       normal(2)= 0.0d0
    endif


    do  i = 1,sdim
      c =  c - normal(i)*center(i)  
    enddo
      c = c + mofdata(0,0,datalen*mat_isd)
!...........................................................................................
    do i = 1, sdim
     normal(i) = -1.0d0*normal(i)
    enddo                                                       ! ATTENSION
    C = -1.0D0*C                                           
!.......................................................................................
    x1%val = centroid
    x2%val = centroid1

    call find_fluxpoint(normal(1),normal(2),c,x1,x2,fluxpoint)

    do i  = -1 , 1
      do j = -1 , 1
         if(mofdata(i,j, datalen * (mat_osd-1) + 1) .gt. eps) then
            osd_tag(i, j) = 1
            lsdata_osd(i, j, :) = mat_cen_sten(i, j, mat_osd, :)                        
         endif
         if(mofdata(i,j, datalen * (mat_isd-1) + 1) .gt. eps) then
            isd_tag(i, j) = 1
            lsdata_isd(i, j, :) = mat_cen_sten(i, j, mat_isd, :)                        
         endif
      enddo
    enddo


! calculate the weight
  do i = -1 , 1
   do j = -1 , 1
    if(osd_tag(i,j) .ne. 0) then
      CALL SDD(dx(1)/2, (NORM_2d(lsdata_osd(i,j,1),lsdata_osd(i,j,2), &
               & fluxpoint%val(1),fluxpoint%val(2)))**2.0d0/dx(1), w_osd(i,j))
    endif

    if(isd_tag(i,j) .ne. 0) then
        CALL SDD(dx(1)/2, (NORM_2d(lsdata_isd(i,j,1),lsdata_isd(i,j,2), &
               & fluxpoint%val(1),fluxpoint%val(2)))**2.0d0/dx(1), w_isd(i,j))
    endif   

   enddo
  enddo


 if(kappa(mat_osd) .gt. eps  .and. kappa(mat_isd) .gt. eps) then

! assemble the matrix
  call Matrix_assemble_double(sdim, nmat, dx,mat_osd,mat_isd, osd_tag, isd_tag, &
                        & lsdata_osd, lsdata_isd, w_osd, w_isd, T, &
                        & normal, fluxpoint, kappa(mat_osd), kappa(mat_isd), &
                        & A,b)



! use cramer's rule to solve for the flux  
  call solve_determinant(4,A,b,x_solved,success_flag)
  grad = x_solved(3)
 elseif(kappa(mat_osd) .ge. 0.0d0  .or. kappa(mat_isd) .ge. 0.0d0)then
  grad = 0.0d0
 else 
  print *,"kappa is invalid in int_grad"
  print *, "kappa(mat_osd)=",kappa(mat_osd),"kappa(mat_isd)=",kappa(mat_isd)
  stop
 endif
 


end subroutine int_grad
!---------------------------------------------------------------------------
subroutine Matrix_assemble_single(sdim, nmat, dx, imat,tag, &
                            & lsdata, w, T, &
                            & normal, fluxpoint, k,&
                            & A,b)
implicit none

    integer,intent(in)            :: nmat,sdim 
    real(kind=8),intent(in)       :: dx(sdim)
    integer,intent(in)                                :: imat
    integer,dimension(-1:1,-1:1),intent(in)           :: tag
    real(kind=8),dimension(-1:1,-1:1,sdim),intent(in) :: lsdata
    real(kind=8),dimension(-1:1,-1:1),intent(in)      :: w
    real(kind=8),dimension(-1:1,-1:1,nmat),intent(in) :: T
   
    real(kind=8),intent(in)                  :: normal(sdim)
    real(kind=8)                  :: m(sdim)
    type(points)                  :: fluxpoint
    real(kind=8),intent(in)       :: k

    

    integer                       :: i,j
    real(kind=8)                  :: A(3,3)
    real(kind=8)                  :: b(3)
    real(kind=8)                  :: x,y,x1,y1


    m = normal
    A = 0.0d0
    b = 0.0d0
    
  !  print *,"T profile"
   if(1 .eq. 0)then
    do i = -1,1
    print *, T(:,i,imat)
    enddo
  endif
    do i = -1 , 1
     do j = -1 , 1
       if(tag(i,j) .eq. 1) then
         
         x = lsdata(i,j,1) - fluxpoint%val(1)
         y = lsdata(i,j,2) - fluxpoint%val(2) 
         A(1,1) = A(1,1) + 2.0d0*w(i,j)*((m(1)*y - m(2)*x)**2.0d0)
         A(1,2) = A(1,2) + 2.0d0*w(i,j)*((m(1)*x + m(2)*y)/k)*(m(1)*y - m(2)*x)
         A(1,3) = A(1,3) + 2.0d0*w(i,j)*(m(1)*y - m(2)*x)
         b(1) = b(1) +  2.0d0*w(i,j)*(m(1)*y - m(2)*x)*T(i,j,imat)       ! ----------T(i,j)-------
       endif
      enddo
     enddo  


    A(2,1) = A(1,2)
    do i = -1 , 1
      do j = -1 , 1
        if(tag(i,j) .eq. 1) then
         x1 = lsdata(i,j,1) - fluxpoint%val(1)
         y1 = lsdata(i,j,2) - fluxpoint%val(2) 
         A(2,2) = A(2,2) + 2.0d0*w(i,j)*((m(1)*x1 + m(2)*y1)**2.0d0)/(k*k)
         A(2,3) = A(2,3) + 2.0d0*w(i,j)*(m(1)*x1 + m(2)*y1)/k
         b(2) = b(2) + 2.0d0*w(i,j)*(m(1)*x1 + m(2)*y1)/k*T(i,j,imat)
        endif       
      enddo
    enddo


    A(3,1) = A(1,3) 
    A(3,2) = A(2,3)
    do i = -1 ,1
     do j = -1 ,1
        if(tag(i,j) .eq. 1) then          
         A(3,3) = A(3,3) + 2.0d0*w(i,j)
         b(3) = b(3) + 2.0d0*w(i,j)*T(i,j,imat)
        endif             
     enddo
    enddo

  do i = 1 , 3
   do j = 1 , 3
     if(abs(A(i,j)) .lt. eps)then
       A(i,j) = 0.0d0
     endif
   enddo
   if(abs(b(i)) .lt. eps) then
     b(i) = 0.0d0
   endif
  enddo
 

end subroutine Matrix_assemble_single


!------------------------------------------------------------------------------
 
 subroutine Matrix_assemble_double(sdim, nmat, dx, mat_osd,mat_isd,osd_tag, isd_tag, &
                            & lsdata_osd, lsdata_isd, w_osd, w_isd, T, &
                            & normal, fluxpoint, k_osd, k_isd,&
                            & A,b)
 implicit none
 
    integer,intent(in)            :: nmat,sdim 
    real(kind=8),intent(in)       :: dx(sdim)
    integer,intent(in)                                :: mat_osd,mat_isd
    integer,dimension(-1:1,-1:1),intent(in)           :: osd_tag,isd_tag
    real(kind=8),dimension(-1:1,-1:1,sdim),intent(in) :: lsdata_osd, lsdata_isd
    real(kind=8),dimension(-1:1,-1:1),intent(in)      :: w_osd, w_isd
    real(kind=8),dimension(-1:1,-1:1,nmat),intent(in) :: T
   
    real(kind=8),intent(in)                  :: normal(sdim)
    real(kind=8)                  :: m(sdim)
    type(points)                  :: fluxpoint
    real(kind=8),intent(in)       :: k_osd, k_isd
    

    integer                       :: i,j
    real(kind=8)                  :: A(4,4)
    real(kind=8)                  :: b(4)
    real(kind=8)                  :: x,y,x1,x2,y1,y2

    m = normal
    A = 0.0d0
    b = 0.0d0



    do i = -1 , 1
     do j = -1 , 1
       if(osd_tag(i,j) .eq. 1) then
         x = lsdata_osd(i,j,1) - fluxpoint%val(1)
         y = lsdata_osd(i,j,2) - fluxpoint%val(2) 
         A(1,1) = A(1,1) + 2.0d0*w_osd(i,j)*((m(1)*y - m(2)*x)**2.0d0)
         A(1,3) = A(1,3) + 2.0d0*w_osd(i,j)*((m(1)*x + m(2)*y)/k_osd)*(m(1)*y - m(2)*x)
         A(1,4) = A(1,4) + 2.0d0*w_osd(i,j)*(m(1)*y - m(2)*x)
         b(1) = b(1) +  2.0d0*w_osd(i,j)*(m(1)*y - m(2)*x)*T(i,j,mat_osd)       ! ----------T(i,j)-------
       endif
      enddo
     enddo
  
    do i = -1 , 1
     do j = -1 , 1       
       if(isd_tag(i,j) .eq. 1) then
         x = lsdata_isd(i,j,1) - fluxpoint%val(1)
         y = lsdata_isd(i,j,2) - fluxpoint%val(2)
         A(2,2) = A(2,2) + 2.0d0*w_isd(i,j)*((m(1)*y - m(2)*x)**2.0d0)
         A(2,3) = A(2,3) + 2.0d0*w_isd(i,j)*((m(1)*x + m(2)*y)/k_isd)*(m(1)*y - m(2)*x)
         A(2,4) = A(2,4) + 2.0d0*w_isd(i,j)*(m(1)*y - m(2)*x)
         b(2) = b(2) +  2.0d0*w_isd(i,j)*(m(1)*y - m(2)*x)*T(i,j,mat_isd)          
       endif
     enddo
    enddo


    A(3,1) = A(1,3)
    A(3,2) = A(2,3)
    do i = -1 , 1
      do j = -1 , 1
        if(osd_tag(i,j) .eq. 1) then
         x1 = lsdata_osd(i,j,1) - fluxpoint%val(1)
         y1 = lsdata_osd(i,j,2) - fluxpoint%val(2)           
         A(3,3) = A(3,3) + 2.0d0*w_osd(i,j)*((m(1)*x1 + m(2)*y1)**2.0d0)/((k_osd)**2.0d0)
         A(3,4) = A(3,4) + 2.0d0*w_osd(i,j)*(m(1)*x1 + m(2)*y1)/k_osd
         b(3) = b(3) + 2.0d0*w_osd(i,j)*(m(1)*x1 + m(2)*y1)/k_osd*T(i,j,mat_osd)
        endif
 
        if(isd_tag(i,j) .eq. 1) then
         x2 = lsdata_isd(i,j,1) - fluxpoint%val(1)
         y2 = lsdata_isd(i,j,2) - fluxpoint%val(2)           
         A(3,3) = A(3,3) + 2.0d0*w_isd(i,j)*((m(1)*x2 + m(2)*y2)**2.0d0)/((k_isd)**2.0d0)
         A(3,4) = A(3,4) + 2.0d0*w_isd(i,j)*(m(1)*x2 + m(2)*y2)/k_isd
         b(3) = b(3) + 2.0d0*w_isd(i,j)*(m(1)*x2 + m(2)*y2)/k_isd*T(i,j,mat_isd)
        endif        
      enddo
    enddo


    A(4,1) = A(1,4) 
    A(4,2) = A(2,4)
    A(4,3) = A(3,4)
    do i = -1 ,1
     do j = -1 ,1
        if(osd_tag(i,j) .eq. 1) then          
         A(4,4) = A(4,4) + 2.0d0*w_osd(i,j)
         b(4) = b(4) + 2.0d0*w_osd(i,j)*T(i,j,mat_osd)
        endif

        if(isd_tag(i,j) .eq. 1) then           
         A(4,4) = A(4,4) + 2.0d0*w_isd(i,j)
         b(4) = b(4) + 2.0d0*w_isd(i,j)*T(i,j,mat_isd)
        endif              
     enddo
    enddo

  do i = 1 , 4
   do j = 1 , 4
     if(abs(A(i,j)) .lt. eps)then
       A(i,j) = 0.0d0
     endif
   enddo
   if(abs(b(i)) .lt. eps) then
     b(i) = 0.0d0
   endif
  enddo




  end subroutine Matrix_assemble_double



  subroutine int_face_adjust(nmat,int_facefrac, int_face )
   implicit none
   
   integer,intent(in)            :: nmat 
   real(kind=8),intent(in)       :: int_facefrac(nmat+1,nmat)
   real(kind=8)                  :: int_face(nmat,nmat)
   real(kind=8)                  :: total_area,areafrac

   integer                       :: i , j

   do i = 1, nmat
    do j = 1, nmat
     int_face(i,j)=0.0
    enddo
   enddo

   do i = 1, nmat

     total_area=int_facefrac(nmat+1,i)

     do j = 1, nmat

      if (i.ne.j) then

       areafrac=int_facefrac(i,j)
       if (areafrac.lt.0.0) then
        print *,"areafrac invalid"
        stop
       else if (areafrac.eq.0.0) then
        ! do nothing
       else if ((areafrac.gt.0.0).and.(areafrac.le.1.0)) then
        if (total_area.le.0.0) then
         print *,"total_area invalid"
         stop
        endif
        int_face(i,j)=areafrac*total_area
        int_face(j,i)=int_face(i,j)
       else
        print *,"areafrac invalid"
        stop
       endif

      endif ! i<>j

     enddo ! j=1..nmat

   enddo ! i=1..nmat


  end subroutine int_face_adjust





  !------------------------------------------------------------
  subroutine convert_cen(nmat,sdim,N,centroid,centroid_mult)
    implicit none

    integer,intent(in)       :: nmat, sdim, N
    type(points),intent(in),dimension(-1:N,-1:N,nmat)  :: centroid
    real(kind=8),dimension(-1:N,-1:N,nmat,sdim)        :: centroid_mult

    integer                   :: i,j,im,dir

    do i = -1,N
       do j= -1,N
          do im = 1,nmat
           do dir=1,sdim
             centroid_mult(i,j,im,dir) = centroid(i,j,im)%val(dir)
           enddo
          enddo
       enddo
    enddo


  end subroutine convert_cen

  !-----------------------------------------------------------
  SUBROUTINE INIT_KAPPA(NMAT,thermal_cond,kappa)
    IMPLICIT NONE

    INTEGER,INTENT(IN)          :: NMAT
    real(kind=8),intent(in)     :: thermal_cond(nmat)
    REAL(KIND=8)                :: KAPPA(nmat,nmat)

    integer                     :: i,j

    !write(16,*) thermal_cond

    do i = 1,nmat
       do j = 1, nmat
          kappa(i,j) = 2.0d0/( 1.0d0/thermal_cond(i) + 1.0d0/thermal_cond(j) )
       enddo
    enddo



  END SUBROUTINE INIT_KAPPA







  !///////////////////////////////////////////////////////////////////
  !-------------------------------------------------------------------
  !              flux calculation
  !-------------------------------------------------------------------
  !///////////////////////////////////////////////////////////////////


  !-----------------------------------------------------------------

  subroutine find_fluxpoint(a,b,c,x1,x2,thetaI)
    ! ax + by + c = 0

    implicit none

    real(kind=8) ,intent(in)    :: a,b,c
    type(points) ,intent(in)    :: x1,x2
    integer                     :: i,j

    type(points) ,intent(out)   :: thetaI



    if(abs(a) .lt. eps  .and. abs(b) .lt. eps)then
       print *,"wrong slope in find_fluxpoint"
    elseif(abs(a) .lt. eps) then
       if(abs(x1%val(2)-x2%val(2)) .lt. eps) then
          print *,"Error, the two lines are parallel"
       else
          thetaI%val(2) = c/(-b)
          thetaI%val(1) = ((x1%val(1)-x2%val(1))*c/(-b) + x1%val(1)*(x1%val(2)-x2%val(2)) &
               & - (x1%val(1)-x2%val(1))*x1%val(2))/(x1%val(2)-x2%val(2))
       endif
    elseif(abs(b) .lt. eps) then
       if(abs(x1%val(1)-x2%val(1)) .lt. eps) then
          print *,"Error, the two lines are parallel"
       else
          thetaI%val(1) = c/(-a)
          thetaI%val(2) = ((x1%val(2)-x2%val(2))*c/(-a) + (x1%val(1)-x2%val(1))*x1%val(2) &
               & -x1%val(1)*(x1%val(2)-x2%val(2)))/(x1%val(1)-x2%val(1))
       endif
    else
       thetaI%val(1) = (-(x1%val(1)-x2%val(1))*x1%val(2) + x1%val(1)*(x1%val(2)-x2%val(2)) &
            & - c/b*(x1%val(1)-x2%val(1)))   &
            &/(x1%val(2)-x2%val(2)+(x1%val(1)-x2%val(1))*a/b)
       thetaI%val(2) = -a/b*thetaI%val(1) - c/b
    endif




  end subroutine find_fluxpoint





  !////////////////////////////////////////////////////////////
  ! use cramer's rule to solve linear system
  !////////////////////////////////////////////////////////////
  !-------------------------------------------------------
  subroutine solve_determinant(N,A_in,b_in,x,success_flag)
    implicit none

    integer, intent(in)          :: N
    real(kind=8),intent(in)      :: A_in(N,N)
    real(kind=8),intent(in)      :: b_in(N)
    real(kind=8),intent(out)     :: x(N)
    integer :: success_flag
    integer                      :: i,j,k
    real(kind=8)                :: A(N,N)
    real(kind=8)                :: b(N)


    integer                      :: mark,scl, flag

    real(kind=8)                 :: deter(N)
    real(kind=8)                 :: den

    flag = 0

    A = A_in
    b = b_in
    
    do i = 1,N
      if(b(i) .ne. 0.0d0 ) then
       flag = 1
      endif
    enddo

    if (1.eq.0) then
     print *, "A = "
     do i = 1, N
      print *, A(:,i)
     enddo
    endif

     !  print *,"b"
     !  print *,b
    call deter_cal(N,A,den)
    if(abs(den) .le. 0.0) then
      if (1.eq.0) then
       print *,"warning, denominator determinant is zero"
       print *, den
      endif
      success_flag=0
    endif

    if(flag .eq. 0) then 
     if (1.eq.0) then
      print *,"Only zero solution exists"
     endif
     x = 0.0d0
    elseif(flag .eq. 1)then
     do i = 1,N
      A(:,i) = b(:)
      call deter_cal(N,A,deter(i))
      x(i) = deter(i)/den
      A = A_in  
     enddo
    else
     print *,"flag error in solve_determinant"
    endif

  end subroutine solve_determinant
  !--------------------------------------------------
  subroutine deter_cal(N,A_in,deter)
    implicit none

    integer, intent(in)       :: N
    real(kind=8),intent(in)   :: A_in(N,N)

    real(kind=8)              :: A(N,N)
    integer                   :: i,j,k
    integer                   :: scl
    !real(kind=8)              :: 
    integer                   :: sgn
    real(kind=8)                 :: coef

    real(kind=8),intent(out)   :: deter
    integer                    :: mark

    A = A_in
    sgn = 1
    mark = 1
  !     print *, "A_in"
  !     do i = 1, N
  !      print *, A_in(:,i) 
  !     enddo

    do j = 1,N-1
       scl = j
       if(abs(A(scl,scl)) .lt. eps) then
          sgn = sgn*(-1)
          call zero_pivot(N,A,scl,mark)
    !    print *, "mark",mark
       endif
      if(mark .ne. 0) then
       do i = scl+1,N
          coef = A(i,j)/A(scl,j)
          do k = scl,N
             A(i,k) = A(i,k) - A(scl,k)*coef
          enddo
       enddo
      else
        exit
      endif
    enddo

    if(mark .ne. 0) then
     deter = sgn
     do i = 1,N
       deter = deter*A(i,i)
     enddo
    else
     deter = 0.0d0
    endif

  end subroutine deter_cal

  !--------------------------------------------------
  subroutine zero_pivot(N,A,iin,mark)
    implicit none

    integer,intent(in)           :: N
    !real(kind=8),intent(in)      :: A_in(N,N)
    real(kind=8)                 :: A(N,N)
    real(kind=8)                 :: B(N),C(N)
    integer,intent(in)           :: iin
    integer                      :: i,j,mark
!    integer                      :: flag 


    !if(abs(A(iin,iin)) .le. eps) then
    mark = 0

    do i = iin+1,N       
       if(abs(A(i,iin)) .gt. eps) then
          mark = i
          exit
       endif
    enddo
    if(mark .ne. 0) then
!       print *,"warning, the determinant is 0"
!       print *, "where:"
!       do i = 1, N
!        print *, A(:,i) 
!       enddo
!      stop
    B(:) = A(iin,:)
    C(:) = A(mark,:)
    A(iin,:) = C(:)
    A(mark,:) = B(:)
    endif


  end subroutine zero_pivot




end module MOF_pair_module






