PROGRAM test 
USE GeneralClass
USE probcommon_module 
USE MOF_pair_module
USE mmat_FVM
!USE multi_solver_module 
USE bicgstab_module
use interp_module

IMPLICIT NONE

! problem type
! 0= flat interface  
! 1= annulus  
! 2= vertical interface
! 3= star for thin filament
! 4 = star for two material canity check
! 5 = hypocycloid with 2 materials
! 6 = nucleate boiling diffusion with thin filament between vapor bubble and substrate
! 7 = hypocycloid with 5 materials
! 8 = nucleate boiling diffusion without filament


! for flat interface, interface is y=0.3.
! for dirichlet, top material has k=0 T(y=0.3)=2.0   T(y=0.0)=3.0

INTEGER,PARAMETER          :: probtype_in = 1
INTEGER,PARAMETER          :: operator_type_in = 1 !0=low,1=simple,2=least sqr
INTEGER,PARAMETER          :: dclt_test_in = 1 ! 1 = Dirichlet test  on
INTEGER,PARAMETER          :: solvtype = 1 ! 0 = CG  1 = bicgstab
INTEGER,PARAMETER          :: N =32 ,M= 1
INTEGER,PARAMETER          :: plot_int = 1
real(kind=8),parameter     :: fixed_dt = 1.25d-2  ! !!!!!!!!!!!!!!!!!!
real(kind=8),parameter     :: CFL = 0.5d0
real(kind=8),parameter     :: problo= 0.0d0, probhi= 1.0d0
integer,parameter          :: sdim_in = 2

INTEGER :: nmat_in
INTEGER :: precond_type_in
INTEGER :: dir
REAL(kind=8) :: xcen,ycen,time_init,xgrid,ygrid
REAL(kind=8) :: deltat_in
REAL(kind=8) :: tol_in
REAL(kind=8) :: current_time_in
REAL(kind=8) :: alpha_in(100)

INTEGER                    :: i,j,k,p,q,tm,imat
integer                    :: i1,i2,j1,j2
REAL(KIND=8)               :: h_in,tau
REAL(KIND=8)               :: time_n,time_np1
REAL(KIND=8),dimension(-1:N+1) :: XLINE,YLINE ! nodes
real(kind=8),dimension(-1:N) :: nox,noy       ! cell centers
REAL(KIND=8)               :: Ts(M)
real(kind=8)               :: dx_in(sdim_in)
TYPE(POLYGON),dimension(-1:N,-1:N):: CELL_FAB
real(kind=8)               :: dist
real(kind=8),external      :: u,v,exact_temperature
real(kind=8)               :: uu(0:N,0:N),vv(0:N,0:N)
character(len=70)          :: lineseg,ctd,triangles
real(kind=8)               :: vol_tot

!real(kind=8)                :: xsten(-3:3,sdim_in)
integer                     :: im,ii,im1,im2
real(kind=8)                :: sumT,sumvf

!---------------------------------------------------
REAL(KIND=8)                :: thermal_cond(100)
integer                     :: nten

real(kind=8)                :: flxtot



!----------------------------------------
INTEGER order_algorithm(1000)
INTEGER MOFITERMAX
INTEGER imaterial_override
INTEGER ngeom_recon_in

integer nx_in,ny_in,lox_in,loy_in,hix_in,hiy_in
integer hflag,vofcomp,nsteps
real(kind=8) :: sum_alpha
real(kind=8), dimension(:,:,:), allocatable :: UNEW_in
real(kind=8), dimension(:,:,:), allocatable :: UOLD_in
real(kind=8), dimension(:,:,:), allocatable :: beta_in
real(kind=8), dimension(:,:,:), allocatable :: VFRAC_MOF_in

! -1:N,-1:N,nmat
real(kind=8),dimension(:,:,:),allocatable :: vf
! -1:N,-1:N,nmat*ngeom_recon
real(kind=8),dimension(:,:,:),allocatable :: mofdata_FAB_in
! -1:N,-1:N,nmat
TYPE(POINTS),DIMENSION(:,:,:),allocatable :: CENTROID_FAB
! -1:N,-1:N,nmat,sdim
real(kind=8),dimension(:,:,:,:),allocatable :: centroid_mult   
! -1:N,-1:N,nmat
real(kind=8),dimension(:,:,:),allocatable :: T
real(kind=8),dimension(:,:,:),allocatable :: T_new


real(kind=8)  ::dtest(N+1,N+1),dtest1(N+1,N+1),dtest2(N+1,N+1)             

!----------------------------------------------------------
real(kind=8)   :: fcenter(2)
real(kind=8)   :: fprobe(2), fI(2)
real(kind=8)   :: fdist,fdist1,fdist3
integer        :: sgflag
!-----------------------------------
real(kind=8)    :: xlo(2)
real(kind=8)    :: interpolate_temperature
real(kind=8)    :: gradtemp
real(kind=8)    :: GRADERR1,GRADERR2,GRADERR3


!NAMELIST /PMTR/ probtype_in, operator_type_in, dclt_test_in, &
!                solvtype, N, M, plot_int, fixed_dt, &
!                problo, probhi, sdim_in

!OPEN (UNIT= 9, FILE='nl.dat', STATUS='UNKNOWN')
! write(9, NML=PMTR) 
! close(9)

!------------------------------------------------------------
if(dclt_test_in .eq. 1) then
 open(unit = 2, file = "out_1")
elseif(dclt_test_in .eq. 0)then
 open(unit= 2 , file= "out_0")
endif
 open(unit=3,file="vf.dat")
 open(unit=4,file="cen.dat")
 open(unit=5,file="cen1.dat")
 open(unit=21,file="cen2.dat")
 open(unit=22,file="cen3.dat")
 open(unit=23,file="cen4.dat")
 open(unit=10,file="levelset.dat")
 open(unit=31,file="levelset1.dat")
 open(unit=32,file="levelset2.dat")
 open(unit=11,file="check.dat")
 open(unit=12,file="para.dat")

 open(unit=41,file="output1.dat")
 open(unit=42,file="output2.dat")
 open(unit=43,file="output3.dat")






!if(probtype_in .eq. 4 .or. probtype_in .eq. 3)then
! pcurve_ls = 0.0d0
! call starshape(pcurve_ls)
!!   do i = 1,pcurve_num
!    write(12,*) (pcurve_ls(:,i)+1.0d0)/2.0d0
!   enddo
! call starshape2(pcurve_ls2)
if( 1 .eq. 1)then
 if(probtype_in .eq. 5 .or. probtype_in .eq. 7)then
  pcurve_ls = 0.0d0
  call asteroidshape(pcurve_ls)
  do i=1,pcurve_num+1
   pcurve_ls(1,i)=(pcurve_ls(1,i)+1.0d0)/2.0d0
   pcurve_ls(2,i)=(pcurve_ls(2,i)+1.0d0)/2.0d0  
  enddo
!   do i = 1,pcurve_num
!    write(12,*) (pcurve_ls(:,i)+1.0d0)/2.0d0
!   enddo
 endif
endif

probtype=probtype_in  ! defined in probdataf95.H (probcommon)
order_algorithm = 0
if (probtype_in.eq.0) then
 nmat_in=2
else if (probtype_in.eq.1) then
 order_algorithm(1)=1
 order_algorithm(2)=3
 order_algorithm(3)=2
 nmat_in=3
else if (probtype_in.eq.2) then
 nmat_in=2
else if(probtype_in .eq. 3)then
 order_algorithm(1)=1
 order_algorithm(2)=3
 order_algorithm(3)=2
 nmat_in=3
elseif(probtype_in .eq. 4)then
 nmat_in=2 
elseif(probtype_in .eq. 5)then
 nmat_in=2
elseif(probtype_in .eq. 6)then
 nmat_in=3
 order_algorithm(1)=1
 order_algorithm(2)=2
 order_algorithm(3)=3
elseif(probtype_in .eq. 7)then
 nmat_in=5
 order_algorithm(1)=1
 order_algorithm(2)=2
 order_algorithm(3)=3
 order_algorithm(4)=4
 order_algorithm(5)=5
else
 print *,"probtype_in invalid"
 stop
endif

MOFITERMAX=10
imaterial_override=0
ngeom_raw=1+BL_SPACEDIM
ngeom_recon=3+2*BL_SPACEDIM
ngeom_recon_in=3+2*BL_SPACEDIM

call initmof(order_algorithm,nmat_in,MOFITERMAX,imaterial_override,0,1)

! initmof(..., mof_debug_recon_in, mof_turn_off_ls_in)    in MOF.F90
!  mof_debug_recon_in = 1 , output a lot ,    
!                     = 0 , nothing
! mof_turn_off_ls_in = 1 ,  not use levelset as input,  use centroid.. ,    
!                    = 0 ,  use levelset as input
! 

print *,"in main.F90: probtype_in= ",probtype_in
print *,"N= ",N
print *,"M= ",M
print *,"fixed_dt= ",fixed_dt
print *,"radcen= ",radcen
print *,"radeps= ",radeps


h_in = (probhi-problo)/N
do dir=1,sdim_in
 dx_in(dir) = h_in
enddo

! init nodes
do i = -1 , N+1
   XLINE(i)= problo + (i)*h_in
   YLINE(i)= problo + (i)*h_in
enddo

do i= -1, N
   nox(i)= problo + (i+0.5d0)*h_in
   noy(i)= problo + (i+0.5d0)*h_in
enddo

!--intial cell
do i = -1 , N
  do j= -1 , N
     call init_cell(N,h_in,nox,noy,i,j,CELL_FAB(i,j))
  enddo
enddo

allocate(vf(-1:N,-1:N,nmat_in)) 
allocate(mofdata_FAB_in(-1:N,-1:N,ngeom_recon_in*nmat_in)) 
allocate(CENTROID_FAB(-1:N,-1:N,nmat_in))  ! type points 
allocate(centroid_mult(-1:N,-1:N,nmat_in,sdim_in)) 
allocate(T(-1:N,-1:N,nmat_in)) 
allocate(T_new(-1:N,-1:N,nmat_in)) 

! init velocity
CALL INIT_V(N,XLINE(0:N),YLINE(0:N),uu,vv)

! time step
!tau = h* CFL /max(maxval(abs(uu)),maxval(abs(vv)))

 if (fixed_dt.eq.0.0) then
  tau = h_in*0.5
 else if (fixed_dt.gt.0.0) then
  tau=fixed_dt
 else
  print *,"fixed_dt invalid"
  stop
 endif

 print *,"tau=",tau

do i = 1,M
  Ts(i) =(i-1)* tau
enddo

  ! TYPE(POINTS),DIMENSION(:,:,:),allocatable :: CENTROID_FAB 
 call init_vfncen(N,CELL_FAB,nmat_in,dx_in,CENTROID_FAB,vf,probtype_in)
 do j = 0,N
  write(3,*) vf(0:N,j,3)
 enddo


  ! real(kind=8),dimension(:,:,:,:),allocatable :: centroid_mult
 call convert_cen(nmat_in,sdim_in,N,CENTROID_FAB,centroid_mult)
 
  
! INIT_MOFdata
 nten = ( (nmat_in-1)*(nmat_in-1)+nmat_in-1 )/2
 call init_mofdata(N,sdim_in,dx_in,nmat_in,nten,CELL_FAB, &
  vf,CENTROID_FAB,mofdata_FAB_in)
 
 do i = 0,N-1
   mofdata_FAB_in(i,-1,:) = mofdata_FAB_in(i,0,:)
   mofdata_FAB_in(i,N,:) = mofdata_FAB_in(i,N-1,:)
 enddo

 do i  = -1,N
   mofdata_FAB_in(-1,i,:) = mofdata_FAB_in(0,i,:)
   mofdata_FAB_in(N,i,:) = mofdata_FAB_in(N-1,i,:)   
 enddo

! init thermal conductivity
 if (probtype_in.eq.0) then
  thermal_cond(1)=10.0d0
  thermal_cond(2)=1.0d0
  if (dclt_test_in.eq.1) then
   thermal_cond(2)=0.1d0  ! top material
  else if (dclt_test_in.eq.0) then
   ! do nothing
  else
   print *,"dclt_test_in is bad"
   stop
  endif
 else if (probtype_in.eq. 1) then
  thermal_cond(1)=0.0d0 
  thermal_cond(2)=1.0d0 
  thermal_cond(3)=0.0d0 
 else if (probtype_in.eq. 2) then
  thermal_cond(1)=1.0d0
  thermal_cond(2)=0.1d0
 elseif(probtype_in .eq. 3)then   ! penta foil with filament
  thermal_cond(1) = 0.0d0
  thermal_cond(2) = 1.0d0
  thermal_cond(3) = 0.0d0
 elseif(probtype_in .eq. 4)then
  thermal_cond(1) = 1.0d0           ! interior region
  thermal_cond(2) = 2.0d0          ! exterior region
 elseif(probtype_in .eq. 5)then
  thermal_cond(1) = 1.0d0           ! interior region
  thermal_cond(2) = 10.0d0          ! exterior region
 elseif(probtype_in .eq. 6)then
  thermal_cond(1) = 1.0d0
  thermal_cond(2) = 0.1d0 
  thermal_cond(3) = 0.01d0
 elseif(probtype_in .eq. 7)then
  thermal_cond(1) = 1.0d0
  thermal_cond(2) = 0.1d0 
  thermal_cond(3) = 1.0d0
  thermal_cond(4) = 0.1d0 
  thermal_cond(5) = 1.0d0
 else 
  print *,"probtype_in invalid"
  stop
 endif

 T = 1.0d0
 do i= -1,N
 do j= -1,N

!  xcen=centroid_mult(i,j,2,1) 
! init centroid   im1 im2 ? ? ? 
!  ycen=centroid_mult(i,j,2,2)

  xcen=centroid_mult(i,j,2,1)
  ycen=centroid_mult(i,j,2,2)

  time_init=0.0

  write(2,*) "centroid", xcen,ycen
   ! 1d problem
  if ((probtype_in.eq.0).or.(probtype_in.eq.2)) then
   T(i,j,1)=2.0
   T(i,j,2)=2.0
   if (1.eq.0) then
    do im = 1,nmat_in
     T(i,j,im)=exact_temperature(xcen,ycen,time_init,im,probtype_in, &
      nmat_in,thermal_cond,dclt_test_in)
    enddo
   endif
  else if (probtype_in.eq.1) then ! annulus problem
!   T(i,j,1)=0.0
!   T(i,j,2)=2.0
!   T(i,j,3)=0.0
  if(1 .eq. 1)then
   do im = 1,nmat_in
    T(i,j,im)=exact_temperature(xcen,ycen,time_init,im,probtype_in, &
     nmat_in,thermal_cond,dclt_test_in)
   enddo
  endif

  elseif (probtype_in .eq. 3)then
   do im = 1,nmat_in
    T(i,j,im)=exact_temperature(xcen,ycen,time_init,im,probtype_in, &
     nmat_in,thermal_cond,dclt_test_in)
   enddo
  elseif(probtype_in .eq. 4)then
   T(i,j,1)=2.0
   T(i,j,2)=2.0  
  elseif(probtype_in .eq. 5)then
   T(i,j,1)=2.0
   T(i,j,2)=2.0  
  elseif(probtype_in .eq. 6)then
   T(i,j,1)=2.0
   T(i,j,2)=2.0    
   T(i,j,3)=2.0 
  elseif(probtype_in .eq. 7)then
   T(i,j,1)=2.0
   T(i,j,2)=2.0    
   T(i,j,3)=2.0  
   T(i,j,4)=2.0    
   T(i,j,5)=2.0 

  else
   print *,"probtype_in invalid"
   stop
  endif

 enddo
 enddo

 do i= -1,N
  do j= -1,N
   do im = 1,nmat_in

    if(vf(i,j,im) .le. 1.0E-8 ) then
     sumT  = 0.0d0
     sumvf = 0.0d0
     do im1 = 1,nmat_in
      sumT  = sumT  + vf(i,j,im1)*T(i,j,im1) 
      sumvf = sumvf + vf(i,j,im1)            
     enddo
     T(i,j,im) = sumT/sumvf
    endif

   enddo ! im
  enddo
 enddo


  


do tm  = 1, M
  print *,"time_step", tm

  write(2,*) "***************************************************"
  write(2,*) "time step", tm
  write(2,*) "***************************************************"

  ngeom_recon_in=2*sdim_in+3
  nx_in=N
  ny_in=N
  tol_in=1.0E-8
  precond_type_in= 1 ! 0 M=I  1=Jacobi precond.
  hflag=0
  deltat_in=tau
  do im=1,nmat_in
   alpha_in(im)=thermal_cond(im)
  enddo
  lox_in=0
  loy_in=0
  hix_in=nx_in-1
  hiy_in=ny_in-1
  allocate(UNEW_in(lox_in-1:hix_in+1,loy_in-1:hiy_in+1,nmat_in)) 
  allocate(UOLD_in(lox_in-1:hix_in+1,loy_in-1:hiy_in+1,nmat_in)) 
  allocate(beta_in(lox_in-1:hix_in+1,loy_in-1:hiy_in+1,nmat_in)) 
  allocate(VFRAC_MOF_in(lox_in-1:hix_in+1,loy_in-1:hiy_in+1, &
     nmat_in)) 

  do i=lox_in-1,hix_in+1
  do j=loy_in-1,hiy_in+1
   do im=1,nmat_in
    xgrid=(i+0.5)*h_in
    ygrid=(j+0.5)*h_in

    sumvf=0.0
    sum_alpha=0.0
    do im1=1,nmat_in
     vofcomp=ngeom_recon_in*(im1-1)+1
     sumvf=sumvf+mofdata_FAB_in(i,j,vofcomp)
     sum_alpha=sum_alpha+mofdata_FAB_in(i,j,vofcomp)/(alpha_in(im1)+1.0E-10)
    enddo
    sum_alpha=sum_alpha/sumvf
    sum_alpha=1.0/sum_alpha 
    beta_in(i,j,im)=sum_alpha

    UNEW_in(i,j,im)=T(i,j,im)
    UOLD_in(i,j,im)=T(i,j,im)
    vofcomp=ngeom_recon_in*(im-1)+1
    VFRAC_MOF_in(i,j,im)=mofdata_FAB_in(i,j,vofcomp)

!    print *,"uold = ", uold_in(i,j,im)

   enddo ! im
  enddo
  enddo 



  write(2,*) "#####################################################################"
  write(2,*) "T initial", "  mat = 1"
  do i1 = loy_in-1,hiy_in+1
   write(2,*) Uold_in(:,i1,1) 
  enddo
   write(2,*) "T initial", "  mat = 2"
  do i1 = loy_in-1,hiy_in+1
   write(2,*) Uold_in(:,i1,2) 
  enddo
  write(2,*)"#####################################################################"


  current_time_in=Ts(tm) ! t^{n}

  


  call INIT_GLOBALS( &
   dclt_test_in, &
   operator_type_in, &
   probtype_in, &
   sdim_in,ngeom_recon_in, &
   nx_in,ny_in,lox_in,loy_in,hix_in,hiy_in, &
   UNEW_in,UOLD_in, &
   beta_in,h_in,precond_type_in,tol_in, &
   VFRAC_MOF_in,nmat_in,alpha_in,deltat_in, &
   mofdata_FAB_in,current_time_in)

  time_n=current_time_in
  time_np1=current_time_in+deltat_in

  nsteps=tm-1
  if (tm.eq.1) then
   call output_solution(UNEW_in,time_n,nsteps,plot_int)
  endif

  write(2,*)"#################################################################"
  write(2,*) "T1", "  mat = 1"
  do i1 = loy_in-1,hiy_in+1
   write(2,*) Unew_in(:,i1,1) 
   write(41,*) Unew_in(:,i1,1) 
  enddo
   write(2,*) "T1", "  mat = 2"
  do i1 = loy_in-1,hiy_in+1
   write(2,*) Unew_in(:,i1,2) 
   write(42,*) Unew_in(:,i1,2) 
  enddo
  write(2,*) "####################################################################"


 if(solvtype .eq. 1)then                        ! solver type
  call bicgstab(UNEW_in,hflag)              
 elseif(solvtype .eq. 0)then
  call cggd(UNEW_in,hflag)
 else
  print *, "solver type invalid"
  stop
 endif


  write(2,*)"#################################################################"
  write(2,*) "T2", "  mat = 1"
  do i1 = loy_in-1,hiy_in+1
   write(2,*) Unew_in(:,i1,1) 
   write(41,*) Unew_in(:,i1,1) 
  enddo
   write(2,*) "T2", "  mat = 2"
  do i1 = loy_in-1,hiy_in+1
   write(2,*) Unew_in(:,i1,2) 
   write(42,*) Unew_in(:,i1,2) 
  enddo
  write(2,*) "####################################################################"




  call set_boundary(UNEW_in,0)

  nsteps=tm
  if (((nsteps/plot_int)*plot_int.eq.nsteps).or. &
      (nsteps.eq.M)) then
   call output_solution(UNEW_in,time_np1,nsteps,plot_int)
  endif

  call DEALLOCATE_GLOBALS()
 
  deallocate(UOLD_in) 
  deallocate(beta_in) 
  deallocate(VFRAC_MOF_in) 

  do i=lox_in-1,hix_in+1
  do j=loy_in-1,hiy_in+1
   do im=1,nmat_in
    T_new(i,j,im)=UNEW_in(i,j,im)
   enddo
  enddo
  enddo

  deallocate(UNEW_in) 

 T = T_new

!  write(2,*) "#######################################################################"
!  write(2,*) "new T ", "  mat = 1"
!  do i1 = loy_in-1,hiy_in+1
!   write(2,*) T(:,i1,1) 
!  enddo
!   write(2,*) "new T ", "  mat = 2"
!  do i1 = loy_in-1,hiy_in+1
!   write(2,*) T(:,i1,2) 
!  enddo
!  write(2,*) "#######################################################################"

   
 
  
enddo ! tm=1,...,M




If(probtype_in .eq. 6)then
 flxtot=0.0d0 
 IF(N .eq. 32) then
  do i = 0,31   
   flxtot=flxtot+ (T(i,27,3)-T(i,26,3))/h
  enddo
 elseif(N .eq. 64)then
  do i = 0,63  
   flxtot=flxtot+ (T(i,54,3)-T(i,53,3))/h
  enddo 
 elseif(N .eq. 128)then
  do i = 0,127  
   flxtot=flxtot+ (T(i,108,3)-T(i,107,3))/h
  enddo 
 elseif(N .eq. 256)then
  do i = 0,255 
   flxtot=flxtot+ (T(i,216,3)-T(i,215,3))/h
  enddo   
 ENDIF

 print *,"flux total =", flxtot

endif



!do i=0,N-1
! do j = 0,N-1
!  call dist_fns(1,cell_FAB(i,j)%center%val(1),cell_FAB(i,j)%center%val(2),dtest(i+1,j+1),4)
! enddo
!enddo


if(probtype_in .eq. 6 .or. probtype_in .eq. 3)then
 do i=0,N
 do j = 0,N
  call dist_fns(1,xline(i),yline(j),dtest(i+1,j+1),probtype_in)
  call dist_fns(2,xline(i),yline(j),dtest1(i+1,j+1),probtype_in)
  call dist_fns(3,xline(i),yline(j),dtest2(i+1,j+1),probtype_in)
 enddo
 enddo
 do j=1,N+1
  write(10,*) dtest(:,j) 
 enddo
 do j=1,N+1
  write(31,*) dtest1(:,j) 
 enddo
 do j=1,N+1
  write(32,*) dtest2(:,j) 
 enddo
endif

if(1 .eq. 0)then
print *,vf(15,16,:)
print *,vf(16,15,:)
print *,vf(15,15,:)
print *,vf(16,16,:)

do i = 1,5
 print *, centroid_mult(15,16,i,:)
enddo
print *,"======================="
do i = 1,5
 print *, centroid_mult(16,15,i,:)
enddo
print *,"======================="
do i = 1,5
 print *, centroid_mult(15,15,i,:)
enddo
print *,"======================="
do i = 1,5
 print *, centroid_mult(16,16,i,:)
enddo
print *,"======================="
print *,cell_fab(15,16)%center%val
print *,cell_fab(16,15)%center%val
endif


if(probtype_in .eq. 1)then             ! flux test for probtype = 1 annulus problem
sgflag = 0
do i=1,2
 xlo(i)=problo
enddo

do i=-1,N
 do j=-1,N
  if(vf(i,j,2) .gt. 0.0001d0)then    !   material 2 cell.
   do ii=1,2
    fcenter(ii)=CELL_FAB(i,j)%center%val(ii)
   enddo
   call dist_fns(2,fcenter(1),fcenter(2), fdist , probtype_in)
   call dist_fns(1,fcenter(1),fcenter(2), fdist1 , probtype_in)
   call dist_fns(3,fcenter(1),fcenter(2), fdist3, probtype_in)
   if(abs(fdist) .lt. 2.0d0*h_in)then   
    if(abs(fdist1) .le. abs(fdist3))then 
     if(fdist1 .lt. 0.0d0)then
      sgflag=-1  
      call find_cloest_2d(-1,fdist,fcenter,xI)
     elseif(fdist1 .gt. 0.0d0)then
      sgflag=+1
      call find_cloest_2d(+1,fdist,fcenter,xI)     
     else
      XI=fcenter
     endif
    elseif(abs(fdist1) .gt. abs(fdist3))then
     if(fdist3 .lt. 0.0d0)then  
      sgflag=+1
      call find_cloest_2d(+1,fdist,fcenter,xI)
     elseif(fdist3 .gt. 0.0d0)then
      sgflag=-1
      call find_cloest_2d(-1,fdist,fcenter,xI)     
     else
      XI=fcenter
     endif    
    else
     print *,"Error 712"
     stop
    endif   ! fdist1 > < fdist3

    if(sgflag .eq. +1)then
     call find_cloest_2d(-1,h_in,xI,xprobe)
    elseif(sgflag .eq. -1)then
     call find_cloest_2d(+1,h_in,xI,xprobe)   
    else
     print *,"check"
     stop
    endif


   call interpfabTEMP( &
       N, &
       h_in, &
       problo, &
       T(:,:,2), &
       vf(:,:,2),&
       centroid_mult(:,:,2,:)
       xprobe, &  
       xI, & 
       interpolate_temperature)

    gradtemp= (interpolate_temperature - &
              exact_temperature(xI(1),xI(2),M*fixed_dt, &
              2,probtype_in,nmat_in,alpha_in,dclt_flag))/abs(fdist)
    



   endif   ! fdist > <  2*h
  
  endif   ! vf(i,j,2)
 enddo
enddo

endif




deallocate(vf)
deallocate(mofdata_FAB_in)
deallocate(CENTROID_FAB)
deallocate(centroid_mult)
deallocate(T)
deallocate(T_new)

 close(2)
 close(3)
 close(4)
 close(5)
 close(10)
 close(11)
 close(12)
 close(21)
 close(22)
 close(23)
 close(31)
 close(32)
 close(41)
 close(42)
 close(43)

END PROGRAM








