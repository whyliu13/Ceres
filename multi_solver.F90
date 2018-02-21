module multi_solver_module
implicit none

      REAL*8, PARAMETER :: AVGTOL=1.0E-8

      integer operator_type ! 0=low order 1=simple 2=high order
      integer probtypeCG ! 0=flat interface  1=annulus 2=vertical interface
      integer sdim
      integer ngeom_reconCG
      integer nx,ny,nz,lox,loy,loz,hix,hiy,hiz
      integer nmat,precond_type,nsmooth
      REAL*8 deltat
      REAL*8 bicgstab_tol
      REAL*8 h,meshvol
      REAL*8 alpha(100)
      REAL*8, dimension(:,:,:), allocatable :: beta
      REAL*8, dimension(:,:,:), allocatable :: UNEW
      REAL*8, dimension(:,:,:), allocatable :: UOLD
      REAL*8, dimension(:,:,:), allocatable :: VFRAC_MOF
      REAL*8, dimension(:,:,:), allocatable :: G
      REAL*8, dimension(:,:,:), allocatable :: DIAG_FIELD

      REAL*8, dimension(:,:,:), allocatable :: mofdata_FAB
      REAL*8                                :: current_time      

      REAL*8, dimension(:,:,:,:), allocatable :: xsten_FAB
      REAL*8, dimension(:,:,:,:), allocatable :: int_face_FAB
      REAL*8, dimension(:,:,:,:,:), allocatable :: int_face_normal_FAB
      REAL*8, dimension(:,:,:,:), allocatable :: dist_to_int_FAB
       ! nmat+1,sdim,2
      REAL*8, dimension(:,:,:,:,:), allocatable :: ext_face_FAB
      REAL*8, dimension(:,:,:,:,:,:), allocatable :: multi_cen_cell_FAB
      REAL*8, dimension(:,:,:,:,:,:), allocatable :: frac_pair_cell_FAB
      REAL*8, dimension(:,:,:,:), allocatable :: centroid_mult_FAB
      REAL*8 dx(3)

      integer dclt_test


contains
!------------------------------------------------------------
      subroutine DEALLOCATE_GLOBALS()
      IMPLICIT NONE

      DEALLOCATE(beta)
      DEALLOCATE(UNEW)
      DEALLOCATE(UOLD)
      DEALLOCATE(VFRAC_MOF)
      DEALLOCATE(G)
      DEALLOCATE(DIAG_FIELD)

      DEALLOCATE(mofdata_FAB)

      DEALLOCATE(xsten_FAB)
      DEALLOCATE(int_face_FAB)
      DEALLOCATE(int_face_normal_FAB)
      DEALLOCATE(dist_to_int_FAB)
      DEALLOCATE(ext_face_FAB)
      DEALLOCATE(multi_cen_cell_FAB)
      DEALLOCATE(frac_pair_cell_FAB)
      DEALLOCATE(centroid_mult_FAB)
      
      return
      end subroutine DEALLOCATE_GLOBALS

!-----------------------------------------------------------
     subroutine INIT_GLOBALS( &
        dclt_test_in, &
        operator_type_in, &
        probtype_in, &
        sdim_in,ngeom_recon_in, &
        nx_in,ny_in,lox_in,loy_in,hix_in,hiy_in, &
        UNEW_in,UOLD_in, &
        beta_in,h_in,precond_type_in,tol_in, &
        VFRAC_MOF_in,nmat_in,alpha_in,deltat_in,&
        mofdata_FAB_in, current_time_in)
      IMPLICIT NONE

      integer dclt_test_in


      integer operator_type_in
      integer, intent(in) :: probtype_in
      integer isten
      integer sdim_in,ngeom_recon_in
      integer nx_in,ny_in,lox_in,loy_in,hix_in,hiy_in
      integer nmat_in,precond_type_in
      integer im,i,j
      REAL*8 deltat_in,tol_in
      REAL*8 h_in
      REAL*8 alpha_in(nmat_in)
      REAL*8 beta_in(lox_in-1:hix_in+1,loy_in-1:hiy_in+1,nmat_in)
      REAL*8 UNEW_in(lox_in-1:hix_in+1,loy_in-1:hiy_in+1,nmat_in)
      REAL*8 UOLD_in(lox_in-1:hix_in+1,loy_in-1:hiy_in+1,nmat_in)
      REAL*8 VFRAC_MOF_in(lox_in-1:hix_in+1,loy_in-1:hiy_in+1,nmat_in)

      REAL*8 DIAGCOEFF,betahalf,UWALL,coeff1,coeff2
      integer dir
      integer sidesten
      integer bctype
      integer imof
      integer basecomp

      REAL*8 vf
      REAL*8 xsrc(sdim_in)
      REAL*8 xgrid,ygrid,xbc,ybc

       ! ngeom_recon_in=vfrac,centroid,order,slope,intercept=2*sdim+3
      REAL*8 mofdata_FAB_in(lox_in-1:hix_in+1,loy_in-1:hiy_in+1, &
       ngeom_recon_in*nmat_in) 
      REAL*8 current_time_in,tn,tnp1,GSRC

      operator_type=operator_type_in

      print *,"ini INIT_GLOBALS: probtype_in= ",probtype_in
      probtypeCG=probtype_in
      print *,"ini INIT_GLOBALS: probtype_in, probtypeCG= ", &
        probtype_in,probtypeCG

      sdim=sdim_in
      ngeom_reconCG=2*sdim+3
      if (ngeom_reconCG.ne.ngeom_recon_in) then
       print *,"ngeom_recon_in invalid"
       stop
      endif
      if ((sdim.ne.2).and.(sdim.ne.3)) then
       print *,"sdim invalid"
       stop
      endif
      nsmooth=8

! 

      dclt_test = dclt_test_in    
! 
      if ((dclt_test.ne.0).and.(dclt_test.ne.1)) then
       print *,"dclt_test invalid"
       stop
      endif
      nx=nx_in
      ny=ny_in
      lox=lox_in
      loy=loy_in
      hix=hix_in
      hiy=hiy_in
      nz=0
      loz=0
      hiz=0

      precond_type=precond_type_in
      deltat=deltat_in
      bicgstab_tol=tol_in
      h=h_in      
      meshvol=h*h

      nmat=nmat_in
  
      do dir=1,sdim
       dx(dir)=h
      enddo
 
      current_time = current_time_in 

      if (deltat.le.0.0) then
       print *,"deltat invalid"
       stop
      endif
      if (current_time.lt.0.0) then
       print *,"current_time invalid"
       stop
      endif
      tn=current_time
      tnp1=current_time+deltat

      do im=1,nmat
       alpha(im)=alpha_in(im)
      enddo

      allocate(beta(lox-1:hix+1,loy-1:hiy+1,nmat)) 
      allocate(UNEW(lox-1:hix+1,loy-1:hiy+1,nmat)) 
      allocate(UOLD(lox-1:hix+1,loy-1:hiy+1,nmat)) 
      allocate(VFRAC_MOF(lox-1:hix+1,loy-1:hiy+1,nmat)) 
      allocate(G(lox-1:hix+1,loy-1:hiy+1,nmat)) 
      allocate(DIAG_FIELD(lox-1:hix+1,loy-1:hiy+1,nmat)) 

      allocate(mofdata_FAB(lox-1:hix+1,loy-1:hiy+1, &
        ngeom_reconCG*nmat))

      allocate(xsten_FAB(lox-1:hix+1,loy-1:hiy+1,-3:3,sdim))
      allocate(int_face_FAB(lox-1:hix+1,loy-1:hiy+1,nmat,nmat))
      allocate(int_face_normal_FAB(lox-1:hix+1,loy-1:hiy+1, &
        nmat,nmat,sdim))
      allocate(dist_to_int_FAB(lox-1:hix+1,loy-1:hiy+1,nmat,nmat))
      allocate(ext_face_FAB(lox-1:hix+1,loy-1:hiy+1,nmat+1,sdim,2))
      allocate(multi_cen_cell_FAB(lox-1:hix+1,loy-1:hiy+1, &
        sdim,nmat,sdim,2))
      allocate(frac_pair_cell_FAB(lox-1:hix+1,loy-1:hiy+1, &
        nmat,nmat,sdim,2))
      allocate(centroid_mult_FAB(lox-1:hix+1,loy-1:hiy+1,nmat,sdim))

      do i=lox-1,hix+1
      do j=loy-1,hiy+1
      do im=1,nmat
       DIAG_FIELD(i,j,im)=1.0E+10
      enddo 
      enddo 
      enddo

    
      

      do i=lox-1,hix+1
      do j=loy-1,hiy+1

       do isten=-3,3
        xsten_FAB(i,j,isten,1)=(i+0.5)*h+isten*h*0.5
        xsten_FAB(i,j,isten,2)=(j+0.5)*h+isten*h*0.5
       enddo

       do im=1,nmat
        xgrid=xsten_FAB(i,j,0,1)
        ygrid=xsten_FAB(i,j,0,2)
      
        beta(i,j,im)=beta_in(i,j,im)
        UNEW(i,j,im)=UNEW_in(i,j,im)
        UOLD(i,j,im)=UOLD_in(i,j,im)
        VFRAC_MOF(i,j,im)=VFRAC_MOF_in(i,j,im)
        vf=VFRAC_MOF(i,j,im)

!        print *, "u_old", uold(i,j,im)

        basecomp=ngeom_reconCG*(im-1)
        do imof=1,ngeom_reconCG 
         mofdata_FAB(i,j,basecomp+imof) = &
          mofdata_FAB_in(i,j,basecomp+imof)
        enddo
        do dir=1,sdim
         centroid_mult_FAB(i,j,im,dir)= &
          mofdata_FAB(i,j,basecomp+1+dir)+xsten_FAB(i,j,0,dir)
        enddo

        do dir=1,sdim
         xsrc(dir)=centroid_mult_FAB(i,j,im,dir)
        enddo

        if (abs(vf).le.AVGTOL) then
         G(i,j,im)=(meshvol/deltat)*UOLD(i,j,im)
        else if ((vf.ge.AVGTOL).and.(vf.le.1.0+AVGTOL)) then
         if (1.eq.0) then
          print *,"probtypeCG,im,sdim before source: ",probtypeCG,im,sdim
         endif 
         call get_filament_source_local(xsrc,tnp1,im,GSRC)
         G(i,j,im)=(meshvol/deltat)*UOLD(i,j,im)+meshvol*GSRC
         if (operator_type.eq.0) then
          ! do nothing
         else if ((operator_type.eq.1).or.(operator_type.eq.2)) then
          G(i,j,im)=G(i,j,im)*vf
         else
          print *,"operator_type invalid"
          stop
         endif
        else
         print *,"vf invalid"
         stop
        endif
 
        DIAGCOEFF=meshvol/deltat


        if ((i.ge.lox).and.(i.le.hix).and. &
            (j.ge.loy).and.(j.le.hiy)) then

         xbc=xgrid
         ybc=ygrid
         betahalf=2.0*beta_in(i-1,j,im)*beta_in(i,j,im)/ &
          (beta_in(i-1,j,im)+beta_in(i,j,im))
         DIAGCOEFF=DIAGCOEFF+betahalf
         if (i.gt.lox) then
          ! do nothing
         else
          dir=1
          sidesten=1
          call getBC_TYPE(dir,sidesten,UWALL,bctype)
          xbc=0.0
          call WALLBC_DIAG(xbc,ybc,coeff2,bctype,dir,sidesten)
          DIAGCOEFF=DIAGCOEFF-betahalf*coeff2
         endif

         betahalf=2.0*beta_in(i+1,j,im)*beta_in(i,j,im)/ &
          (beta_in(i+1,j,im)+beta_in(i,j,im))
         DIAGCOEFF=DIAGCOEFF+betahalf
         if (i.lt.hix) then
          ! do nothing
         else
          dir=1
          sidesten=2
          call getBC_TYPE(dir,sidesten,UWALL,bctype)
          xbc=nx*h
          call WALLBC_DIAG(xbc,ybc,coeff2,bctype,dir,sidesten)
          DIAGCOEFF=DIAGCOEFF-betahalf*coeff2
         endif

         xbc=xgrid
         ybc=ygrid

         betahalf=2.0*beta_in(i,j-1,im)*beta_in(i,j,im)/ &
          (beta_in(i,j-1,im)+beta_in(i,j,im))
 
         DIAGCOEFF=DIAGCOEFF+betahalf
         if (j.gt.loy) then
          ! do nothing
         else
          dir=2
          sidesten=1
          call getBC_TYPE(dir,sidesten,UWALL,bctype)
          ybc=0.0
          call WALLBC_DIAG(xbc,ybc,coeff2,bctype,dir,sidesten)
          DIAGCOEFF=DIAGCOEFF-betahalf*coeff2
         endif

         betahalf=2.0*beta_in(i,j+1,im)*beta_in(i,j,im)/ &
          (beta_in(i,j+1,im)+beta_in(i,j,im))
         DIAGCOEFF=DIAGCOEFF+betahalf
         if (j.lt.hiy) then
          ! do nothing
         else
          dir=2
          sidesten=2
          call getBC_TYPE(dir,sidesten,UWALL,bctype)
          ybc=ny*h
          call WALLBC_DIAG(xbc,ybc,coeff2,bctype,dir,sidesten)
          DIAGCOEFF=DIAGCOEFF-betahalf*coeff2
         endif
        endif ! i,j in the interior

        if (DIAGCOEFF.le.0.0) then
         print *,"matrix should be positive definite"
         print *,"DIAGCOEFF=",DIAGCOEFF
         print *,"i,j,im,xgrid,ygrid,vf ",i,j,im,xgrid,ygrid,vf
         print *,"beta= ",beta(i,j,im)
         print *,"h,meshvol ",h,meshvol
         print *,"dx= ",dx(1),dx(2)
         stop
        endif
        DIAG_FIELD(i,j,im)=DIAGCOEFF
       enddo ! im
      enddo
      enddo

      call build_MM_diag()

      return
      end subroutine INIT_GLOBALS


































end module multi_solver_module
