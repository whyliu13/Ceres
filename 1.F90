      module bicgstab_module

      integer operator_type ! 0=low order 1=simple 2=high order
      integer probtype ! 0=flat interface  1=annulus
      integer sdim
      integer ngeom_reconCG
      integer nx,ny,lox,loy,hix,hiy
      integer nmat,precond_type,nsmooth
      REAL*8 deltat,diag_factor
      REAL*8 bicgstab_tol
      REAL*8 h,meshvol
      REAL*8 alpha(100)
      REAL*8, dimension(:,:,:), allocatable :: beta
      REAL*8, dimension(:,:,:), allocatable :: UNEW
      REAL*8, dimension(:,:,:), allocatable :: UOLD
      REAL*8, dimension(:,:,:), allocatable :: VFRAC_MOF
      REAL*8, dimension(:,:,:), allocatable :: G
      REAL*8, dimension(:,:,:), allocatable :: DIAG_FIELD

      REAL*8, dimension(:,:,:), allocatable :: mofdata
      REAL*8                                :: current_time      

      REAL*8, dimension(:,:,:,:), allocatable :: xsten_FAB
      REAL*8, dimension(:,:,:), allocatable :: int_face_FAB
      REAL*8, dimension(:,:,:,:,:), allocatable :: int_face_normal_FAB
      REAL*8, dimension(:,:,:,:), allocatable :: dist_to_int_FAB
       ! (nmat+1)*2*sdim
      REAL*8, dimension(:,:,:), allocatable :: ext_face_FAB
      REAL*8, dimension(:,:,:,:,:), allocatable :: multi_cen_cell_FAB
      REAL*8, dimension(:,:,:,:,:), allocatable :: frac_pair_cell_FAB
      REAL*8, dimension(:,:,:,:), allocatable :: centroid_mult_FAB
      REAL*8 dx(3)
 
      contains


      subroutine WALLBC_DIAG(x,y,coeff2,bctype,dir,sidesten)
      IMPLICIT NONE

      integer dir,sidesten
      REAL*8 x,y
      integer bctype
      REAL*8 coeff2

      if (bctype.eq.0) then
       coeff2=1.0
      else if (bctype.eq.1) then
       coeff2=-1.0
      else
       print *,"bctype invalid bctype=",bctype
       stop
      endif
      
      return
      end subroutine WALLBC_DIAG


   
! bctype=0 Neumann 
! bctype=1 Dirichlet

      subroutine DEALLOCATE_GLOBALS()
      IMPLICIT NONE

      DEALLOCATE(beta)
      DEALLOCATE(UNEW)
      DEALLOCATE(UOLD)
      DEALLOCATE(VFRAC_MOF)
      DEALLOCATE(G)
      DEALLOCATE(DIAG_FIELD)

      DEALLOCATE(mofdata)

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

      subroutine build_MM_diag()
      use MOF_pair_module

      IMPLICIT NONE

      REAL*8 mat_cen_sten(-1:1,-1:1,nmat,sdim)
      REAL*8 frac_pair_cell(nmat,nmat,sdim*2)
      REAL*8 int_face_cell(nmat*nmat)
      REAL*8 int_face_normal_cell(nmat,nmat,sdim)
      REAL*8 int_facefrac_cell((nmat+1)*nmat)
      REAL*8 ext_facefrac_cell((nmat+1)*2*sdim)
      REAL*8 dist_to_int_cell(nmat,nmat)
      REAL*8 xsten_cell(-3:3,sdim)
      REAL*8 localface((nmat+1)*2*sdim)
      REAL*8 multi_cen_cell(sdim,nmat,sdim*2)
      integer i,j,k,imof,p,q,r,im,dir,isten
      integer ii,jj,vofcomp
      REAL*8 vf,diag_local
      REAL*8 mofdata_cell(ngeom_reconCG*nmat)
      REAL*8 thin_cen_sten(-1:1,sdim,sdim,nmat,sdim*2)
      REAL*8 ext_face_sten(-1:1,sdim,(nmat+1)*2*sdim)
      integer im_in

      do i=lox-1,hix+1
      do j=loy-1,hiy+1
       do k=1,(nmat+1)*2*sdim
        ext_face_FAB(i,j,k)=0.0
       enddo
       do k=1,nmat*nmat
        int_face_FAB(i,j,k)=0.0
       enddo
       do p=1,nmat
       do q=1,nmat
        dist_to_int_FAB(i,j,p,q)=0.0
       enddo
       enddo
       do p=1,nmat
       do q=1,nmat
       do r=1,sdim
        int_face_normal_FAB(i,j,p,q,r)=0.0
       enddo
       enddo
       enddo
       do p=1,sdim
       do q=1,nmat
       do r=1,2*sdim
        multi_cen_cell_FAB(i,j,p,q,r)=0.0
       enddo
       enddo
       enddo
       do p=1,nmat
       do q=1,nmat
       do r=1,2*sdim
        frac_pair_cell_FAB(i,j,p,q,r)=0.0
       enddo
       enddo
       enddo
      enddo
      enddo ! i,j

      do i=lox-1,hix+1
      do j=loy-1,hiy+1
      
       do imof=1,ngeom_reconCG*nmat
        mofdata_cell(imof)=mofdata(i,j,imof)
       enddo
       do isten=-3,3
       do dir=1,sdim
        xsten_cell(isten,dir)=xsten_FAB(i,j,isten,dir)
       enddo
       enddo

       call ptb_int(ngeom_reconCG, &
         nmat,sdim,mofdata_cell,dx,xsten_cell, &
         int_facefrac_cell, &
         int_face_normal_cell, &
         dist_to_int_cell)

       do p=1,nmat
       do q=1,nmat
        dist_to_int_FAB(i,j,p,q)=dist_to_int_cell(p,q)
       enddo
       enddo

       do p=1,nmat
       do q=1,nmat
       do r=1,sdim
        int_face_normal_FAB(i,j,p,q,r)=int_face_normal_cell(p,q,r)
       enddo
       enddo
       enddo

       call ptb_ext(ngeom_reconCG, &
        nmat,sdim,dx,mofdata_cell, &
        xsten_cell, &
        ext_facefrac_cell, &
        multi_cen_cell)

       do k=1,(nmat+1)*2*sdim
        ext_face_FAB(i,j,k)=ext_facefrac_cell(k)
       enddo

       do p=1,sdim
       do q=1,nmat
       do r=1,2*sdim
        multi_cen_cell_FAB(i,j,p,q,r)=multi_cen_cell(p,q,r)
       enddo
       enddo
       enddo
       
       call int_face_adjust(nmat,int_facefrac_cell,int_face_cell)
       do k=1,nmat*nmat
        int_face_FAB(i,j,k)=int_face_cell(k)
       enddo

      enddo
      enddo ! i,j

      do i=lox-1,hix+1
       do k=1,(nmat+1)*2*sdim
        ext_face_FAB(i,loy-1,k)=ext_face_FAB(i,loy,k)
        ext_face_FAB(i,hiy+1,k)=ext_face_FAB(i,hiy,k)
       enddo
       do k=1,nmat*nmat
        int_face_FAB(i,loy-1,k)=int_face_FAB(i,loy,k)
        int_face_FAB(i,hiy+1,k)=int_face_FAB(i,hiy,k)
       enddo
      enddo ! i

      do j=loy-1,hiy+1
       do k=1,(nmat+1)*2*sdim
        ext_face_FAB(lox-1,j,k)=ext_face_FAB(lox,j,k)
        ext_face_FAB(hix+1,j,k)=ext_face_FAB(hix,j,k)
       enddo
       do k=1,nmat*nmat
        int_face_FAB(lox-1,j,k)=int_face_FAB(lox,j,k)
        int_face_FAB(hix+1,j,k)=int_face_FAB(hix,j,k)
       enddo
      enddo

      do i=lox,hix
      do j=loy,hiy

       do dir = 1,sdim
        do ii = -1,1
         if (dir .eq. 1) then   
          thin_cen_sten(ii,dir,:,:,:) = multi_cen_cell_FAB(i+ii,j,:,:,:)
          ext_face_sten(ii,dir,:) = ext_face_FAB(i+ii,j,:)
         elseif(dir .eq. 2)then
          thin_cen_sten(ii,dir,:,:,:) = multi_cen_cell_FAB(i,j+ii,:,:,:)
          ext_face_sten(ii,dir,:) = ext_face_FAB(i,j+ii,:)
         endif
        enddo   ! ii
       enddo   ! dir

       call vfrac_pair_cell(nmat,sdim,dx,ext_face_sten,thin_cen_sten,  &
         frac_pair_cell)

       do p=1,nmat
       do q=1,nmat
       do r=1,2*sdim
        frac_pair_cell_FAB(i,j,p,q,r)=frac_pair_cell(p,q,r)
       enddo
       enddo
       enddo
      enddo
      enddo ! i,j

      do i=lox,hix
       do p=1,nmat
       do q=1,nmat
       do r=1,2*sdim
        frac_pair_cell_FAB(i,loy-1,p,q,r) = &
         frac_pair_cell_FAB(i,loy,p,q,r)
        frac_pair_cell_FAB(i,hiy+1,p,q,r) = &
         frac_pair_cell_FAB(i,hiy,p,q,r)     
       enddo
       enddo
       enddo
      enddo

      do j=loy-1,hiy+1
       do p=1,nmat
       do q=1,nmat
       do r=1,2*sdim
        frac_pair_cell_FAB(lox-1,j,p,q,r) = &
         frac_pair_cell_FAB(lox,j,p,q,r)
        frac_pair_cell_FAB(hix+1,j,p,q,r) = &
         frac_pair_cell_FAB(hix,j,p,q,r) 
       enddo
       enddo
       enddo
      enddo

      do i=lox,hix
      do j=loy,hiy
       do imof=1,ngeom_reconCG*nmat
        mofdata_cell(imof)=mofdata(i,j,imof)
       enddo
       do ii = -1 , 1
       do jj = -1 , 1
        do im=1,nmat
        do dir=1,sdim
         mat_cen_sten(ii,jj,im,dir) =  &
           centroid_mult_FAB(i+ii,j+jj,im,dir)
        enddo
        enddo
       enddo
       enddo

       do p=1,nmat
       do q=1,nmat
       do r=1,2*sdim
        frac_pair_cell(p,q,r)=frac_pair_cell_FAB(i,j,p,q,r)
       enddo
       enddo
       enddo
       do p=1,nmat
       do q=1,nmat
        dist_to_int_cell(p,q)=dist_to_int_FAB(i,j,p,q)
       enddo
       enddo
       do k=1,nmat*nmat
        int_face_cell(k)=int_face_FAB(i,j,k)
       enddo
       do p=1,nmat
       do q=1,nmat
       do r=1,sdim
        int_face_normal_cell(p,q,r)=int_face_normal_FAB(i,j,p,q,r)
       enddo
       enddo
       enddo
       do isten=-3,3
       do dir=1,sdim
        xsten_cell(isten,dir)=xsten_FAB(i,j,isten,dir)
       enddo
       enddo

       do im_in=1,nmat
        call cell_diag_cal( &
         ngeom_reconCG, &
         sdim,nmat,dx, &
         mofdata_cell, &
         alpha,  &
         mat_cen_sten, &
         im_in, &
         frac_pair_cell, &
         int_face_cell, &
         int_face_normal_cell, &
         dist_to_int_cell, &
         xsten_cell, &
         diag_local)
        vofcomp=(im_in-1)*ngeom_reconCG+1
        vf=mofdata(i,j,vofcomp)
        if (vf.le.1.0E-8) then
         diag_local=meshvol/deltat
        else
         diag_local=(vf*meshvol/deltat)+diag_local
        endif
        if ((operator_type.eq.1).or.(operator_type.eq.2)) then
         DIAG_FIELD(i,j,im_in)=diag_local
        else if (operator_type.eq.0) then
         ! do nothing
        else
         print *,"operator_type invalid"
         stop
        endif
       enddo ! im_in

      enddo
      enddo

      return
      end subroutine build_MM_diag
 
      subroutine INIT_GLOBALS( &
        probtype_in, &
        sdim_in,ngeom_recon_in, &
        nx_in,ny_in,lox_in,loy_in,hix_in,hiy_in, &
        UNEW_in,UOLD_in, &
        beta_in,h_in,precond_type_in,tol_in, &
        VFRAC_MOF_in,nmat_in,alpha_in,deltat_in,&
        mofdata_in, current_time_in)
      USE mmat_FVM

      IMPLICIT NONE

      integer probtype_in
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
      REAL*8 DIAGCOEFF,betahalf,x,y,UWALL,coeff1,coeff2
      integer dir
      integer sidesten
      integer bctype
      integer imof
      integer basecomp

      REAL*8 AVG_TOL,vf
      REAL*8 xsrc(sdim_in)

       ! ngeom_recon=vfrac,centroid,order,slope,intercept=2*sdim+3
      REAL*8 mofdata_in(lox_in-1:hix_in+1,loy_in-1:hiy_in+1, &
       ngeom_recon_in*nmat_in) 
      REAL*8 current_time_in,tn,tnp1,GSRC

      probtype=probtype_in

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

      operator_type=0

      AVG_TOL=1.0E-8

      nx=nx_in
      ny=ny_in
      lox=lox_in
      loy=loy_in
      hix=hix_in
      hiy=hiy_in
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
      diag_factor=1.0/deltat

      allocate(beta(lox-1:hix+1,loy-1:hiy+1,nmat)) 
      allocate(UNEW(lox-1:hix+1,loy-1:hiy+1,nmat)) 
      allocate(UOLD(lox-1:hix+1,loy-1:hiy+1,nmat)) 
      allocate(VFRAC_MOF(lox-1:hix+1,loy-1:hiy+1,nmat)) 
      allocate(G(lox-1:hix+1,loy-1:hiy+1,nmat)) 
      allocate(DIAG_FIELD(lox-1:hix+1,loy-1:hiy+1,nmat)) 

      allocate(mofdata(lox-1:hix+1,loy-1:hiy+1,ngeom_reconCG*nmat))

      allocate(xsten_FAB(lox-1:hix+1,loy-1:hiy+1,-3:3,sdim))
      allocate(int_face_FAB(lox-1:hix+1,loy-1:hiy+1,nmat*nmat))
      allocate(int_face_normal_FAB(lox-1:hix+1,loy-1:hiy+1, &
        nmat,nmat,sdim))
      allocate(dist_to_int_FAB(lox-1:hix+1,loy-1:hiy+1,nmat,nmat))
      allocate(ext_face_FAB(lox-1:hix+1,loy-1:hiy+1,(nmat+1)*2*sdim))
      allocate(multi_cen_cell_FAB(lox-1:hix+1,loy-1:hiy+1, &
        sdim,nmat,sdim*2))
      allocate(frac_pair_cell_FAB(lox-1:hix+1,loy-1:hiy+1, &
        nmat,nmat,sdim*2))
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
        xsten_FAB(i,j,isten,2)=(j+0.5)*h+isten+h*0.5
       enddo

       do im=1,nmat
        x=xsten_FAB(i,j,0,1)
        y=xsten_FAB(i,j,0,2)
      
        beta(i,j,im)=beta_in(i,j,im)
        UNEW(i,j,im)=UNEW_in(i,j,im)
        UOLD(i,j,im)=UOLD_in(i,j,im)
        VFRAC_MOF(i,j,im)=VFRAC_MOF_in(i,j,im)
        vf=VFRAC_MOF(i,j,im)

        basecomp=ngeom_reconCG*(im-1)
        do imof=1,ngeom_reconCG 
         mofdata(i,j,basecomp+imof) = mofdata_in(i,j,basecomp+imof)
        enddo
        do dir=1,sdim
         centroid_mult_FAB(i,j,im,dir)=mofdata(i,j,basecomp+1+dir)+  &
          xsten_FAB(i,j,0,dir)
        enddo

        if (operator_type.eq.0) then
         G(i,j,im)=(meshvol/deltat)*UOLD(i,j,im)
        else if ((operator_type.eq.1).or.(operator_type.eq.2)) then
         if (vf.le.AVG_TOL) then
          G(i,j,im)=(meshvol/deltat)*UOLD(i,j,im)
         else
          do dir=1,sdim
           xsrc(dir)=centroid_mult_FAB(i,j,im,dir)
          enddo
          call get_filament_source(xsrc,tnp1,probtype,im,sdim,GSRC)
          G(i,j,im)=(meshvol/deltat)*vf*UOLD(i,j,im)+ &
           meshvol*vf*GSRC
         endif
        else
         print *,"operator_type invalid"
         stop
        endif
 
        DIAGCOEFF=diag_factor*meshvol

        betahalf=2.0*beta_in(i-1,j,im)*beta_in(i,j,im)/ &
          (beta_in(i-1,j,im)+beta_in(i,j,im))
        DIAGCOEFF=DIAGCOEFF+betahalf
        if (i.gt.lox) then
         ! do nothing
        else
         dir=1
         sidesten=1
         call getBC_TYPE(dir,sidesten,UWALL,bctype)
         x=0.0
         call WALLBC_DIAG(x,y,coeff2,bctype,dir,sidesten)
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
         x=nx*h
         call WALLBC_DIAG(x,y,coeff2,bctype,dir,sidesten)
         DIAGCOEFF=DIAGCOEFF-betahalf*coeff2
        endif

        betahalf=2.0*beta_in(i,j-1,im)*beta_in(i,j,im)/ &
         (beta_in(i,j-1,im)+beta_in(i,j,im))

        DIAGCOEFF=DIAGCOEFF+betahalf
        if (j.gt.loy) then
         ! do nothing
        else
         dir=2
         sidesten=1
         call getBC_TYPE(dir,sidesten,UWALL,bctype)
         y=0.0
         call WALLBC_DIAG(x,y,coeff2,bctype,dir,sidesten)
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
         y=ny*h
         call WALLBC_DIAG(x,y,coeff2,bctype,dir,sidesten)
         DIAGCOEFF=DIAGCOEFF-betahalf*coeff2
        endif
        if (DIAGCOEFF.le.0.0) then
         print *,"matrix should be positive definite"
         stop
        endif
        DIAG_FIELD(i,j,im)=DIAGCOEFF
       enddo ! im
      enddo
      enddo

      call build_MM_diag()

      return
      end subroutine INIT_GLOBALS

      subroutine WALLBC(x,y,UGHOST,UIN,UWALL,coeff1,coeff2, &
        hflag,bctype,dir,sidesten)
      IMPLICIT NONE

      integer dir,sidesten
      REAL*8 x,y
      integer hflag,bctype
      REAL*8 UGHOST,UWALL,UIN,coeff1,coeff2

      if (bctype.eq.0) then
       UGHOST=UIN
       coeff1=0.0
       coeff2=1.0
      else if (bctype.eq.1) then
       if (hflag.eq.1) then
        UGHOST=-UIN
        coeff1=0.0
        coeff2=-1.0
       else if (hflag.eq.0) then
        UGHOST=2.0*UWALL-UIN
        coeff1=2.0*UWALL
        coeff2=-1.0
       else
        print *,"hflag invalid"
        stop
       endif
      else
       print *,"bctype invalid bctype=",bctype
       stop
      endif
      
      return
      end subroutine WALLBC



      subroutine getBC_TYPE(dir,sidesten,UWALL,bctype)
      IMPLICIT NONE

      integer dir,sidesten,bctype
      REAL*8 UWALL  

      if (probtype.eq.0) then
       if ((dir.eq.1).and.(sidesten.eq.1)) then
        bctype=0 ! neumann
        UWALL=0.0
       else if ((dir.eq.1).and.(sidesten.eq.2)) then
        bctype=0 ! neumann
        UWALL=0.0
       else if ((dir.eq.2).and.(sidesten.eq.1)) then
        bctype=1 ! dirichlet
        UWALL=3.0
       else if ((dir.eq.2).and.(sidesten.eq.2)) then
        bctype=1 ! dirichlet
        UWALL=2.0
       else
        print *,"dir or sidesten invalid"
        stop
       endif
      else if (probtype.eq.1) then
       if ((dir.eq.1).and.(sidesten.eq.1)) then
        bctype=1 ! dirichlet
        UWALL=0.0
       else if ((dir.eq.1).and.(sidesten.eq.2)) then
        bctype=0 ! dirichlet
        UWALL=0.0
       else if ((dir.eq.2).and.(sidesten.eq.1)) then
        bctype=1 ! dirichlet
        UWALL=0.0
       else if ((dir.eq.2).and.(sidesten.eq.2)) then
        bctype=1 ! dirichlet
        UWALL=0.0
       else
        print *,"dir or sidesten invalid"
        stop
       endif
      else
       print *,"probtype invalid"
       stop
      endif

      return
      end subroutine getBC_TYPE

      subroutine set_boundary(U,hflag)
      IMPLICIT NONE

      integer i,j,dir,sidesten,im
      REAL*8 U(lox-1:hix+1,loy-1:hiy+1,nmat)
      REAL*8 x,y
      REAL*8 coeff1,coeff2
      integer hflag
      integer ilo,ihi,jlo,jhi,ii,jj
      integer bctype
      REAL*8 UWALL

      do im=1,nmat

      do dir=1,2
      do sidesten=1,2
       call getBC_TYPE(dir,sidesten,UWALL,bctype)
        ! xlo
       if ((dir.eq.1).and.(sidesten.eq.1)) then
        ilo=lox-1
        ihi=lox-1
        jlo=loy
        jhi=hiy
        ii=1
        jj=0
        ! xhi
       else if ((dir.eq.1).and.(sidesten.eq.2)) then
        ilo=hix+1
        ihi=hix+1
        jlo=loy
        jhi=hiy
        ii=-1
        jj=0
        ! jlo
       else if ((dir.eq.2).and.(sidesten.eq.1)) then
        jlo=loy-1
        jhi=loy-1
        ilo=lox
        ihi=hix
        ii=0
        jj=1
        ! jhi
       else if ((dir.eq.2).and.(sidesten.eq.2)) then
        jlo=hiy+1
        jhi=hiy+1
        ilo=lox
        ihi=hix
        ii=0
        jj=-1
       else
        print *,"dir or sidesten invalid"
        stop
       endif

       do i=ilo,ihi
       do j=jlo,jhi
        x=(i+0.5)*h+0.5*ii*h
        y=(j+0.5)*h+0.5*jj*h 
        call WALLBC(x,y,U(i,j,im),U(i+ii,j+jj,im),UWALL,coeff1,coeff2, &
         hflag,bctype,dir,sidesten)
       enddo 
       enddo 
      enddo 
      enddo 

      enddo ! im

      return
      end subroutine set_boundary

      subroutine DOTPROD(U,V,dsum)
      IMPLICIT NONE

      integer i,j,im
      REAL*8 U(lox-1:hix+1,loy-1:hiy+1,nmat)
      REAL*8 V(lox-1:hix+1,loy-1:hiy+1,nmat)
      REAL*8 dsum,DOT_TOL

      DOT_TOL=1.0E-8
      dsum=0.0
      do im=1,nmat
      do i=lox,hix
      do j=loy,hiy
       if (VFRAC_MOF(i,j,im).ge.DOT_TOL) then
        dsum=dsum+U(i,j,im)*V(i,j,im)
       endif
      enddo
      enddo 
      enddo 
  
      return
      end subroutine DOTPROD

      subroutine NORMPROD(U,dnorm)
      IMPLICIT NONE

      REAL*8 dnorm
      REAL*8 U(lox-1:hix+1,loy-1:hiy+1,nmat)
      REAL*8 dsum

      call DOTPROD(U,U,dsum) 
      dnorm=sqrt(dsum)

      return
      end subroutine NORMPROD

       ! W=AA*U + BB*V
      subroutine LINCOMB(U,V,W,AA,BB) 
      IMPLICIT NONE

      integer i,j,im
      REAL*8 U(lox-1:hix+1,loy-1:hiy+1,nmat)
      REAL*8 V(lox-1:hix+1,loy-1:hiy+1,nmat)
      REAL*8 W(lox-1:hix+1,loy-1:hiy+1,nmat)
      REAL*8 AA,BB,DOT_TOL

      DOT_TOL=1.0E-8
      do im=1,nmat
      do i=lox,hix
      do j=loy,hiy
       if (VFRAC_MOF(i,j,im).ge.DOT_TOL) then
        W(i,j,im)=AA*U(i,j,im)+BB*V(i,j,im)
       endif
      enddo
      enddo
      enddo

      return
      end subroutine LINCOMB

        ! V=U
      subroutine COPYVEC(U,V) 
      IMPLICIT NONE

      integer i,j,im
      REAL*8 U(lox-1:hix+1,loy-1:hiy+1,nmat)
      REAL*8 V(lox-1:hix+1,loy-1:hiy+1,nmat)
      do im=1,nmat
      do i=lox-1,hix+1
      do j=loy-1,hiy+1
       V(i,j,im)=U(i,j,im)
      enddo
      enddo
      enddo

      return
      end subroutine COPYVEC

      subroutine ZAPVEC(U) 
      IMPLICIT NONE

      integer i,j,im
      REAL*8 U(lox-1:hix+1,loy-1:hiy+1,nmat)
      do im=1,nmat
      do i=lox-1,hix+1
      do j=loy-1,hiy+1
       U(i,j,im)=0.0
      enddo
      enddo
      enddo

      return
      end subroutine ZAPVEC

       ! for k=1,2,3 ....
       ! Z^{k}=Z^{k-1}+D^{-1}(R-AZ^{k-1})
       ! e.g. A=D-L-U
       ! D Z^k = D Z^k-1 + R - (D-L-U)Z^k-1
       ! D Z^k - (L+U)Z^k-1=R 
      subroutine JACPRECOND(Z,R,hflag)
      IMPLICIT NONE

      integer i,j,im,iter
      REAL*8 R(lox-1:hix+1,loy-1:hiy+1,nmat)
      REAL*8 Z(lox-1:hix+1,loy-1:hiy+1,nmat)
      integer hflag

      REAL*8, dimension(:,:,:), allocatable :: ZN
      REAL*8, dimension(:,:,:), allocatable :: ZNP1
      REAL*8, dimension(:,:,:), allocatable :: ZSTAR

      allocate(ZN(lox-1:hix+1,loy-1:hiy+1,nmat)) 
      allocate(ZNP1(lox-1:hix+1,loy-1:hiy+1,nmat)) 
      allocate(ZSTAR(lox-1:hix+1,loy-1:hiy+1,nmat)) 

       ! ZN=Z
      call COPYVEC(Z,ZN)

      do iter=1,nsmooth
       call set_boundary(ZN,hflag)
       call RESID(ZSTAR,R,ZN,hflag)

       do im=1,nmat
       do i=lox,hix
       do j=loy,hiy
        ZNP1(i,j,im)=ZN(i,j,im)+(1.0/DIAG_FIELD(i,j,im))*ZSTAR(i,j,im)
       enddo
       enddo
       enddo
        ! ZN=ZNP1
       call COPYVEC(ZNP1,ZN)
      enddo  ! iter
        ! Z=ZN
      call COPYVEC(ZN,Z)

      deallocate(ZN) 
      deallocate(ZNP1) 
      deallocate(ZSTAR) 

      return
      end subroutine JACPRECOND

      subroutine MAKE_CONSISTENT(U,hflag)
      IMPLICIT NONE

      integer i,j,im,hflag
      REAL*8 U(lox-1:hix+1,loy-1:hiy+1,nmat)
      REAL*8 AVG_TOL,UAVG,FAVG

      AVG_TOL=1.0E-8

      call set_boundary(U,hflag)
      do i=lox-1,hix+1
      do j=loy-1,hiy+1
       UAVG=0.0
       FAVG=0.0
       do im=1,nmat
        if (VFRAC_MOF(i,j,im).gt.AVG_TOL) then
         UAVG=UAVG+U(i,j,im)*VFRAC_MOF(i,j,im)
         FAVG=FAVG+VFRAC_MOF(i,j,im)
        endif
       enddo
       if (abs(FAVG-1.0).ge.1.0E-4) then
        print *,"FAVG invalid"
        print *,"i,j ",i,j
        print *,"FAVG=",FAVG
        print *,"VFRAC(1): ",VFRAC_MOF(i,j,1)
        print *,"VFRAC(2): ",VFRAC_MOF(i,j,2)
        stop
       endif
       UAVG=UAVG/FAVG
       do im=1,nmat
        if (VFRAC_MOF(i,j,im).le.AVG_TOL) then
         U(i,j,im)=UAVG
        endif 
       enddo

      enddo 
      enddo 

      return
      end subroutine MAKE_CONSISTENT


      subroutine ATIMESU(AU,U,hflag)
      use MOF_pair_module

      IMPLICIT NONE

      integer i,j,k,im,ii,jj,imof,dir,isten
      integer p,q,r
      REAL*8 betahalfxlo
      REAL*8 betahalfxhi
      REAL*8 betahalfylo
      REAL*8 betahalfyhi
      REAL*8 DIAGCOEFF
      REAL*8 U(lox-1:hix+1,loy-1:hiy+1,nmat)
      REAL*8 AU(lox-1:hix+1,loy-1:hiy+1,nmat)
      integer hflag
      REAL*8 xsten_cell(-3:3,sdim)
      REAL*8 vf
      REAL*8 AVG_TOL
      REAL*8 rho_box(-1:1,-1:1,nmat)
      REAL*8 mofdata_sten(-1:1,-1:1,ngeom_reconCG*nmat)
      REAL*8 mat_cen_sten(-1:1,-1:1,nmat,sdim)
      REAL*8 frac_pair_cell(nmat,nmat,sdim*2)
      REAL*8 int_face_cell(nmat*nmat)
      REAL*8 int_face_normal_cell(nmat,nmat,sdim)
      REAL*8 dist_to_int_cell(nmat,nmat)
      REAL*8 div_tot
      integer success_flag

      AVG_TOL=1.0E-8

      call MAKE_CONSISTENT(U,hflag)

      call set_boundary(U,hflag)

      if (operator_type.eq.0) then

       do im=1,nmat
       do i=lox,hix
       do j=loy,hiy
 
        betahalfxlo=2.0*beta(i-1,j,im)*beta(i,j,im)/ &
         (beta(i-1,j,im)+beta(i,j,im))
        betahalfxhi=2.0*beta(i+1,j,im)*beta(i,j,im)/ &
         (beta(i+1,j,im)+beta(i,j,im))
        betahalfylo=2.0*beta(i,j-1,im)*beta(i,j,im)/ &
         (beta(i,j-1,im)+beta(i,j,im))
        betahalfyhi=2.0*beta(i,j+1,im)*beta(i,j,im)/ &
         (beta(i,j+1,im)+beta(i,j,im))
        DIAGCOEFF=diag_factor*meshvol+(betahalfxlo+betahalfxhi+ &
                 betahalfylo+betahalfyhi)
        AU(i,j,im)=DIAGCOEFF*U(i,j,im)-(betahalfxlo*U(i-1,j,im)+ &
         betahalfxhi*U(i+1,j,im)+betahalfylo*U(i,j-1,im)+ &
         betahalfyhi*U(i,j+1,im))
       enddo
       enddo
       enddo

      else if ((operator_type.eq.1).or.(operator_type.eq.2)) then

       do i=lox,hix
       do j=loy,hiy

        do isten=-3,3
        do dir=1,sdim
         xsten_cell(isten,dir)=xsten_FAB(i,j,isten,dir)
        enddo
        enddo
        do ii=-1,1
        do jj=-1,1
         do im=1,nmat
          rho_box(ii,jj,im)=U(i+ii,j+jj,im)
         enddo
         do imof=1,ngeom_reconCG*nmat
          mofdata_sten(ii,jj,imof)=mofdata(i+ii,j+jj,imof)
         enddo
        enddo
        enddo
        do ii = -1 , 1
        do jj = -1 , 1
         do im=1,nmat
         do dir=1,sdim
          mat_cen_sten(ii,jj,im,dir) = &
            centroid_mult_FAB(i+ii,j+jj,im,dir)
         enddo
         enddo
        enddo
        enddo
        do p=1,nmat
        do q=1,nmat
        do r=1,2*sdim
         frac_pair_cell(p,q,r)=frac_pair_cell_FAB(i,j,p,q,r)
        enddo
        enddo
        enddo
        do k=1,nmat*nmat
         int_face_cell(k)=int_face_FAB(i,j,k)
        enddo
        do p=1,nmat
        do q=1,nmat
        do r=1,sdim
         int_face_normal_cell(p,q,r)=int_face_normal_FAB(i,j,p,q,r)
        enddo
        enddo
        enddo
        do p=1,nmat
        do q=1,nmat
         dist_to_int_cell(p,q)=dist_to_int_FAB(i,j,p,q)
        enddo
        enddo

        do im=1,nmat
         call cell_div_cal( &
          ngeom_reconCG, &
          sdim,nmat,dx, &
          xsten_cell, &
          mofdata_sten, &
          rho_box, &
          alpha, &
          mat_cen_sten, &
          im, &
          frac_pair_cell, &
          int_face_cell, &
          int_face_normal_cell, &
          dist_to_int_cell, &
          div_tot,success_flag)

         if ((operator_type.eq.1).or.(success_flag.eq.0)) then
          call cell_div_cal_simple( &
           ngeom_reconCG, &
           sdim,nmat,dx, &
           xsten_cell, &
           rho_box, &
           alpha, &
           mat_cen_sten, &
           im, &
           frac_pair_cell, &
           int_face_cell, &
           int_face_normal_cell, &
           dist_to_int_cell, &
           div_tot)
         else if ((operator_type.eq.2).and.(success_flag.eq.1)) then
          ! do nothing
         else
          print *,"operator_type or success_flag invalid"
          stop
         endif

         vf=VFRAC_MOF(i,j,im)
         if (vf.le.AVG_TOL) then
          AU(i,j,im)=(meshvol/deltat)*U(i,j,im)
         else        
          AU(i,j,im)=(meshvol/deltat)*vf*U(i,j,im)+div_tot 
         endif

        enddo ! im
       enddo
       enddo ! i,j

      else
       print *,"operator_type invalid"
       stop
      endif
        
      return
      end subroutine ATIMESU




        ! R=RHS-A U
      subroutine RESID(R,RHS,U,hflag)
      IMPLICIT NONE

      integer i,j
      REAL*8 R(lox-1:hix+1,loy-1:hiy+1,nmat)
      REAL*8 RHS(lox-1:hix+1,loy-1:hiy+1,nmat)
      REAL*8 U(lox-1:hix+1,loy-1:hiy+1,nmat)
      REAL*8 x,y,AA,BB
      integer hflag
      REAL*8, dimension(:,:,:), allocatable :: AU

      allocate(AU(lox-1:hix+1,loy-1:hiy+1,nmat)) 

      call ATIMESU(AU,U,hflag)
      AA=1.0
      BB=-1.0
      call LINCOMB(RHS,AU,R,AA,BB)

      deallocate(AU)

      return
      end subroutine RESID


      subroutine preconditioner(Z,R,hflag)
      IMPLICIT NONE

      REAL*8 R(lox-1:hix+1,loy-1:hiy+1,nmat)
      REAL*8 Z(lox-1:hix+1,loy-1:hiy+1,nmat)
      integer hflag

      call ZAPVEC(Z)
      if (precond_type.eq.0) then
       call COPYVEC(R,Z)
      else if (precond_type.eq.1) then
       call JACPRECOND(Z,R,hflag)
      else
       print *,"precond_type invalid" 
       stop
      endif

      return
      end subroutine preconditioner

      subroutine check_fab(U,id)
      IMPLICIT NONE

      integer id
      REAL*8 U(lox-1:hix+1,loy-1:hiy+1,nmat)
      REAL*8 biggest,smallest,biginv
      integer i,j,im

      biggest=-1.0E+10
      smallest=1.0E+10
      biginv=0.0
      do im=1,nmat
      do i=lox-1,hix+1
      do j=loy-1,hiy+1
       if (U(i,j,im).gt.biggest) then
        biggest=U(i,j,im)
       endif
       if (U(i,j,im).lt.smallest) then
        smallest=U(i,j,im)
       endif
       if (U(i,j,im).ne.0.0) then
        if (1.0/abs(U(i,j,im)).gt.biginv) then
         biginv=1.0/abs(U(i,j,im))
        endif
       endif
      enddo
      enddo
      enddo

      print *,"id,biggest,smallest ",id,biggest,smallest
      print *,"id,biginv ",id,biginv

      return
      end subroutine check_fab

        ! precond_type=0 M=I, =1 Jacobi
      subroutine bicgstab(U,hflag)
      IMPLICIT NONE

      integer i,j,im,hflag,iter
      integer maxiter
      integer hflagcg,restart_flag
      REAL*8 U(lox-1:hix+1,loy-1:hiy+1,nmat)
      REAL*8 dnorm,dnorm0

      REAL*8, dimension(:,:,:), allocatable :: U0
      REAL*8, dimension(:,:,:), allocatable :: V0
      REAL*8, dimension(:,:,:), allocatable :: P0
      REAL*8, dimension(:,:,:), allocatable :: R0
      REAL*8, dimension(:,:,:), allocatable :: U1
      REAL*8, dimension(:,:,:), allocatable :: V1
      REAL*8, dimension(:,:,:), allocatable :: P1
      REAL*8, dimension(:,:,:), allocatable :: R1
      REAL*8, dimension(:,:,:), allocatable :: R0hat
      REAL*8, dimension(:,:,:), allocatable :: UINIT
      REAL*8, dimension(:,:,:), allocatable :: RHS
      REAL*8, dimension(:,:,:), allocatable :: Y
      REAL*8, dimension(:,:,:), allocatable :: Hvec
      REAL*8, dimension(:,:,:), allocatable :: S
      REAL*8, dimension(:,:,:), allocatable :: T
      REAL*8, dimension(:,:,:), allocatable :: Z

      REAL*8 rho0,w0,rho1,AA,w1,BB,a1,a2

      allocate(U0(lox-1:hix+1,loy-1:hiy+1,nmat)) 
      allocate(V0(lox-1:hix+1,loy-1:hiy+1,nmat)) 
      allocate(P0(lox-1:hix+1,loy-1:hiy+1,nmat)) 
      allocate(R0(lox-1:hix+1,loy-1:hiy+1,nmat)) 
      allocate(U1(lox-1:hix+1,loy-1:hiy+1,nmat)) 
      allocate(V1(lox-1:hix+1,loy-1:hiy+1,nmat)) 
      allocate(P1(lox-1:hix+1,loy-1:hiy+1,nmat)) 
      allocate(R1(lox-1:hix+1,loy-1:hiy+1,nmat)) 
      allocate(R0hat(lox-1:hix+1,loy-1:hiy+1,nmat)) 
      allocate(UINIT(lox-1:hix+1,loy-1:hiy+1,nmat)) 
      allocate(RHS(lox-1:hix+1,loy-1:hiy+1,nmat)) 
      allocate(Y(lox-1:hix+1,loy-1:hiy+1,nmat)) 
      allocate(Hvec(lox-1:hix+1,loy-1:hiy+1,nmat)) 
      allocate(S(lox-1:hix+1,loy-1:hiy+1,nmat)) 
      allocate(T(lox-1:hix+1,loy-1:hiy+1,nmat)) 
      allocate(Z(lox-1:hix+1,loy-1:hiy+1,nmat)) 

        ! U0=V0=P0=0 
      AA=0.0 
      do im=1,nmat
      do i=lox-1,hix+1
      do j=loy-1,hiy+1
       U0(i,j,im)=AA
       V0(i,j,im)=AA
       P0(i,j,im)=AA
      enddo
      enddo
      enddo

      if (1.eq.0) then
       call check_fab(U0,1)
       call check_fab(G,2)
       call check_fab(beta,3)
      endif

       ! R0=G-A U0
      call RESID(R0,G,U0,hflag)

      if (1.eq.0) then
       call check_fab(R0,4)
      endif
       ! R0hat=R0
      call COPYVEC(R0,R0hat)
       ! UINIT=U0
      call COPYVEC(U0,UINIT)
      call ZAPVEC(U0)
       ! RHS=R0
      call COPYVEC(R0,RHS)

       ! rho0=AA=w0=1
      rho0=1.0
      AA=1.0
      w0=1.0

      call NORMPROD(R0,dnorm0)
      print *,"initial,dnorm0 ",dnorm0
      dnorm=1.0
      maxiter=2000
      iter=0

      hflagcg=1
      do while ((dnorm.gt.bicgstab_tol).and.(iter.lt.maxiter))
       print *,"iter,dnorm ",iter,dnorm

         ! rho1= R0hat^H R0
       call DOTPROD(R0hat,R0,rho1)
       if (1.eq.0) then
        print *,"rho1=",rho1
       endif
       
       restart_flag=0
       if ((sqrt(abs(rho0)).lt.bicgstab_tol*0.01).or. &
           (sqrt(abs(w0)).lt.bicgstab_tol*0.01)) then
        restart_flag=1
       endif 

       if (restart_flag.eq.0) then
          ! (R0hat^H R0)/(R0hat^H dot R0_before)  *   (AA/w0)
        BB=(rho1/rho0)*(AA/w0)
        a1=1.0
        a2=-w0
 
         ! P1=P0-w0 V0
        call LINCOMB(P0,V0,P1,a1,a2)
         ! P1=R0+BB P1
        call LINCOMB(R0,P1,P1,a1,BB)
        ! Y=M^{-1}P1
        call preconditioner(Y,P1,hflagcg)
         
         ! V1=A Y
        call ATIMESU(V1,Y,hflagcg)

         ! AA=rho1/R0hat dot V1
        call DOTPROD(R0hat,V1,AA)

        if (sqrt(abs(AA)).lt.bicgstab_tol*0.01) then
         restart_flag=1
        endif

        if (restart_flag.eq.0) then
         AA=rho1/AA

         ! Hvec=U0+AA Y
         a1=1.0
         a2=AA
         call LINCOMB(U0,Y,Hvec,a1,a2) 
         ! U1=Hvec
         call COPYVEC(Hvec,U1)
         ! R1=RHS-A U1
         call RESID(R1,RHS,U1,hflagcg)
         call NORMPROD(R1,dnorm)
         if (1.eq.0) then
          print *,"dnorm 573 =",dnorm
          print *,"dnorm0 573 =",dnorm0
         endif
         dnorm=dnorm/dnorm0

         if (dnorm.gt.bicgstab_tol) then
          ! S=R0-AA V1
          a1=1.0
          a2=-AA
          call LINCOMB(R0,V1,S,a1,a2) 

           ! Z=M^{-1}S
          call preconditioner(Z,S,hflagcg)

           ! T=A Z
          call ATIMESU(T,Z,hflagcg)

           ! simple case is: (T,S)/(T,T)=(AZ,S)/(AZ,AZ)   (MZ=S)
          call DOTPROD(T,S,a1)
          call DOTPROD(T,T,a2)
          if (sqrt(abs(a2)).lt.bicgstab_tol*0.01) then
           restart_flag=1
          endif

          if (restart_flag.eq.0) then
           w1=a1/a2
           ! U1=Hvec+w1 Z
           a1=1.0
           a2=w1
           call LINCOMB(Hvec,Z,U1,a1,a2) 
          endif
         endif ! dnorm>bicgstab_tol
          ! R1=RHS-A U1
         call RESID(R1,RHS,U1,hflagcg)
         call NORMPROD(R1,dnorm)
         dnorm=dnorm/dnorm0
         rho0=rho1
         w0=w1
          ! R0=R1
         call COPYVEC(R1,R0) 
         call COPYVEC(P1,P0) 
         call COPYVEC(V1,V0) 
         call COPYVEC(U1,U0) 
        endif  ! restart_flag=0
       endif ! restart_flag=0

       if (restart_flag.eq.0) then
        ! do nothing
       else if (restart_flag.eq.1) then
        call RESID(R0,RHS,U0,hflagcg) ! R0=RHS-A U0
         ! R0hat=R0
        call COPYVEC(R0,R0hat)
        call COPYVEC(U0,U1) 
         ! rho0=AA=w0=1
        rho0=1.0
        AA=1.0
        w0=1.0
        call NORMPROD(R0,dnorm)
        dnorm=dnorm/dnorm0
        call ZAPVEC(V0)
        call ZAPVEC(P0)
       else
        print *,"restart_flag invalid"
        stop
       endif

       iter=iter+1
      enddo
      print *,"at the end: iter,dnorm ",iter,dnorm
       ! U=UINIT+U1
      a1=1.0
      a2=a1
      call LINCOMB(UINIT,U1,U,a1,a2)

      call MAKE_CONSISTENT(U,hflag)
 
      deallocate(U0) 
      deallocate(V0) 
      deallocate(P0) 
      deallocate(R0) 
      deallocate(U1) 
      deallocate(V1) 
      deallocate(P1) 
      deallocate(R1) 
      deallocate(R0hat) 
      deallocate(UINIT) 
      deallocate(RHS) 
      deallocate(Y) 
      deallocate(Hvec) 
      deallocate(S) 
      deallocate(T) 
      deallocate(Z) 

      return
      end subroutine bicgstab

      subroutine output_solution(UNEW)
      IMPLICIT NONE

      character*8 wavedatafile
      REAL*8 UNEW(lox-1:hix+1,loy-1:hiy+1,nmat)
      integer i,j
      REAL*8 x,y

      call set_boundary(UNEW,0)

      write(wavedatafile,'(A8)') 'wavedata'
      print *,"wavedatafile ",wavedatafile
      open(unit=11,file=wavedatafile)

      if (nmat.eq.2) then
       if (1.eq.0) then
        write(11,*) '# X,Y,U1,U2,F1,F2'
       else 
        write(11,*) 'VARIABLES="X","Y","U1","U2","F1","F2"'
        write(11,*) 'zone i=',hix-lox+3,' j=', &
         hiy-loy+3,' f=point'
       endif
      else if (nmat.eq.3) then
       if (1.eq.0) then
        write(11,*) '# X,Y,U1,U2,U3,F1,F2,F3'
       else
        write(11,*) 'VARIABLES="X","Y","U1","U2","U3","F1","F2","F3"'
        write(11,*) 'zone i=',hix-lox+3,' j=', &
         hiy-loy+3,' f=point'
       endif
      else
       print *,"nmat not supported"
       stop
      endif

      do i=lox-1,hix+1
      do j=loy-1,hiy+1
       x=(i+0.5)*h
       y=(j+0.5)*h
       if (nmat.eq.2) then
        write(11,*) x,y,UNEW(i,j,1),UNEW(i,j,2), &
         VFRAC_MOF(i,j,1), &
         VFRAC_MOF(i,j,2)
       else if (nmat.eq.3) then
        write(11,*) x,y,UNEW(i,j,1), &
         UNEW(i,j,2), &
         UNEW(i,j,3), &
         VFRAC_MOF(i,j,1), &
         VFRAC_MOF(i,j,2), &
         VFRAC_MOF(i,j,3)
       else
        print *,"nmat not supported"
        stop
       endif
      enddo
      enddo

      close(11)

      return
      end subroutine output_solution

      end module bicgstab_module

