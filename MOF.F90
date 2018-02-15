#undef  BL_LANG_CC
#ifndef BL_LANG_FORT
#define BL_LANG_FORT
#endif

#define STANDALONE 1

#include "REAL.H"
#include "CONSTANTS.H"
#include "SPACE.H"

#if (STANDALONE==0)
#include "MOF_F.H"
#endif
 
#define MLSVOFTOL (1.0E-14)
#define CENTOL (1.0E-13)
#define MAXTET (5)
#define MAXAREA (5)
#define MOFHH (1.0E-4)
#define MOF_ANGLE_MAX (0.1)
#define INTERCEPT_TOL (1.0E-10)
#define GAUSSNEWTONTOL (1.0E-8)
#define RADIUS_CUTOFF (6.0)
#define VOF_CUTOFF (1.0E-8)
#define UNCAPT_TOL (2.0E-2)
#define LSTHICK (1.0E-4)

! Author: Mark Sussman sussman@math.fsu.edu
! Department of Mathematics
! Florida State University
! Tallahassee, FL 32306
!

module geometry_intersect_module

implicit none

INTEGER_T, PARAMETER :: MAX_NUM_MATERIALS=20

INTEGER_T, PARAMETER :: maxfacelist=20
INTEGER_T, PARAMETER :: maxnodelist=20,maxtetlist=15,maxcapfacelist=5
INTEGER_T, PARAMETER :: maxmappednodes=8
INTEGER_T, PARAMETER :: n_vol_listmax=400,n_area_listmax=200

! the 3rd component of a facelist is 0 in 2D
type intersect_type
  INTEGER_T :: n_nodes,n_tet,n_faces,n_capfaces
  INTEGER_T :: n_pos_nodes
  INTEGER_T :: nodelistmap(maxnodelist,3)  
  INTEGER_T :: nodelist(maxnodelist)
  INTEGER_T :: tetlist(maxtetlist,4)
  INTEGER_T :: facelist(maxfacelist,3)
  INTEGER_T :: aligned(maxfacelist)
  INTEGER_T :: capfacelist(maxcapfacelist,3)
end type intersect_type

type mapping_type
  INTEGER_T :: mapped_nodes(maxmappednodes)
  type(intersect_type) :: intersect_geometry
end type mapping_type

! check sum and headnode determine geometry map to use
type(mapping_type), dimension(256,8)  :: hexahedron_maps
type(mapping_type), dimension(16,4)    :: tetrahedron_maps
type(mapping_type), dimension(16,4)    :: rectangle_maps
type(mapping_type), dimension(8,3)    :: triangle_maps

type(intersect_type), dimension(9) :: template_hex_plane
type(intersect_type), dimension(4) :: template_tet_plane
type(intersect_type), dimension(4) :: template_rec_plane
type(intersect_type), dimension(3) :: template_tri_plane

contains

subroutine copy_intersect_type(source,dest)
IMPLICIT NONE

type(intersect_type), intent(in) :: source
type(intersect_type), intent(out) :: dest

INTEGER_T i,dir

 dest%n_nodes=source%n_nodes
 dest%n_tet=source%n_tet
 dest%n_faces=source%n_faces
 dest%n_capfaces=source%n_capfaces
 dest%n_pos_nodes=source%n_pos_nodes
 do i=1,maxnodelist
  do dir=1,3
   dest%nodelistmap(i,dir)=0
  enddo
  dest%nodelist(i)=0
 enddo
 do i=1,maxtetlist
  do dir=1,4
   dest%tetlist(i,dir)=0
  enddo
 enddo
 do i=1,maxfacelist
  do dir=1,3
   dest%facelist(i,dir)=0
  enddo
  dest%aligned(i)=0
 enddo
 do i=1,maxcapfacelist
  do dir=1,3
   dest%capfacelist(i,dir)=0
  enddo
 enddo

 do i=1,source%n_nodes
  do dir=1,3
   dest%nodelistmap(i,dir)=source%nodelistmap(i,dir)
  enddo
  dest%nodelist(i)=source%nodelist(i)
 enddo
 do i=1,source%n_tet
  do dir=1,4
   dest%tetlist(i,dir)=source%tetlist(i,dir)
  enddo
 enddo
 do i=1,source%n_faces
  do dir=1,3
   dest%facelist(i,dir)=source%facelist(i,dir)
  enddo
  dest%aligned(i)=source%aligned(i)
 enddo
 do i=1,source%n_capfaces
  do dir=1,3
   dest%capfacelist(i,dir)=source%capfacelist(i,dir)
  enddo
 enddo

return
end subroutine copy_intersect_type


subroutine fast_copy_intersect_type(source,dest,sdim)
IMPLICIT NONE

type(intersect_type), intent(in) :: source
type(intersect_type), intent(out) :: dest
INTEGER_T :: sdim
INTEGER_T :: i,dir

 if ((sdim.ne.2).and.(sdim.ne.3)) then
  print *,"sdim invalid"
  stop
 endif

 dest%n_nodes=source%n_nodes
 dest%n_tet=source%n_tet
 dest%n_capfaces=source%n_capfaces

 do i=1,source%n_nodes
  do dir=2,3
   dest%nodelistmap(i,dir)=source%nodelistmap(i,dir)
  enddo
 enddo
 do i=1,source%n_tet
  do dir=1,sdim+1
   dest%tetlist(i,dir)=source%tetlist(i,dir)
  enddo
 enddo
 do i=1,source%n_capfaces
  do dir=1,sdim
   dest%capfacelist(i,dir)=source%capfacelist(i,dir)
  enddo
 enddo

return
end subroutine fast_copy_intersect_type



subroutine add_to_hex(gridmap,template_geom)
IMPLICIT NONE

type(intersect_type) :: template_geom
INTEGER_T :: gridmap(2,2,2)
INTEGER_T :: mapped_nodes(maxmappednodes)
INTEGER_T :: power2(maxmappednodes)
INTEGER_T :: i,j,k,nn,n_nodes,inode,checksum

 power2(1)=1
 do i=2,maxmappednodes
  power2(i)=2*power2(i-1)
 enddo

 do i=1,maxmappednodes
  mapped_nodes(i)=0
 enddo

 n_nodes=template_geom%n_nodes
 if ((n_nodes.lt.1).or.(n_nodes.gt.maxnodelist)) then
  print *,"n_nodes invalid"
  stop
 endif

 inode=1
 checksum=0
 do k=1,2
 do j=1,2
 do i=1,2
  mapped_nodes(inode)=gridmap(i,j,k)
  if ((mapped_nodes(inode).lt.1).or.(mapped_nodes(inode).gt.8)) then
   print *,"mapped nodes bust"
   stop
  endif

  do nn=1,n_nodes
   if ((template_geom%nodelistmap(nn,2).eq.inode).and. &
       (template_geom%nodelistmap(nn,3).eq.inode)) then
    checksum=checksum+power2(mapped_nodes(inode))
   endif
  enddo

  inode=inode+1
 enddo
 enddo
 enddo
 if ((checksum.lt.1).or.(checksum.gt.255)) then
  print *,"checksum invalid"
  stop
 endif
 call copy_intersect_type(template_geom, &
   hexahedron_maps(checksum,mapped_nodes(1))%intersect_geometry)
 do i=1,maxmappednodes
  hexahedron_maps(checksum,mapped_nodes(1))%mapped_nodes(i)=mapped_nodes(i)
 enddo

return
end subroutine add_to_hex


subroutine add_to_tet(linemap,template_geom)
IMPLICIT NONE

type(intersect_type) :: template_geom
INTEGER_T :: linemap(4)
INTEGER_T :: mapped_nodes(maxmappednodes)
INTEGER_T :: power2(maxmappednodes)
INTEGER_T :: i,nn,n_nodes,inode,checksum

 power2(1)=1
 do i=2,maxmappednodes
  power2(i)=2*power2(i-1)
 enddo

 do i=1,maxmappednodes
  mapped_nodes(i)=0
 enddo

 n_nodes=template_geom%n_nodes
 if ((n_nodes.lt.1).or.(n_nodes.gt.maxnodelist)) then
  print *,"n_nodes invalid"
  stop
 endif

 inode=1
 checksum=0
 do inode=1,4
  mapped_nodes(inode)=linemap(inode)
  if ((mapped_nodes(inode).lt.1).or.(mapped_nodes(inode).gt.4)) then
   print *,"mapped nodes bust"
   stop
  endif

  do nn=1,n_nodes
   if ((template_geom%nodelistmap(nn,2).eq.inode).and. &
       (template_geom%nodelistmap(nn,3).eq.inode)) then
    checksum=checksum+power2(mapped_nodes(inode))
   endif
  enddo

 enddo

 if ((checksum.lt.1).or.(checksum.gt.15)) then
  print *,"checksum invalid"
  stop
 endif
 call copy_intersect_type(template_geom, &
   tetrahedron_maps(checksum,mapped_nodes(1))%intersect_geometry)
 do i=1,maxmappednodes
  tetrahedron_maps(checksum,mapped_nodes(1))%mapped_nodes(i)=mapped_nodes(i)
 enddo

return
end subroutine add_to_tet


subroutine add_to_tri(linemap,template_geom)
IMPLICIT NONE

type(intersect_type) :: template_geom
INTEGER_T :: linemap(3)
INTEGER_T :: mapped_nodes(maxmappednodes)
INTEGER_T :: power2(maxmappednodes)
INTEGER_T :: i,nn,n_nodes,inode,checksum

 power2(1)=1
 do i=2,maxmappednodes
  power2(i)=2*power2(i-1)
 enddo

 do i=1,maxmappednodes
  mapped_nodes(i)=0
 enddo

 n_nodes=template_geom%n_nodes
 if ((n_nodes.lt.1).or.(n_nodes.gt.maxnodelist)) then
  print *,"n_nodes invalid"
  stop
 endif

 inode=1
 checksum=0
 do inode=1,3
  mapped_nodes(inode)=linemap(inode)
  if ((mapped_nodes(inode).lt.1).or.(mapped_nodes(inode).gt.3)) then
   print *,"mapped nodes bust"
   stop
  endif

  do nn=1,n_nodes
   if ((template_geom%nodelistmap(nn,2).eq.inode).and. &
       (template_geom%nodelistmap(nn,3).eq.inode)) then
    checksum=checksum+power2(mapped_nodes(inode))
   endif
  enddo

 enddo

 if ((checksum.lt.1).or.(checksum.gt.7)) then
  print *,"checksum invalid"
  stop
 endif
 call copy_intersect_type(template_geom, &
   triangle_maps(checksum,mapped_nodes(1))%intersect_geometry)
 do i=1,maxmappednodes
  triangle_maps(checksum,mapped_nodes(1))%mapped_nodes(i)=mapped_nodes(i)
 enddo

return
end subroutine add_to_tri



subroutine add_to_rec(gridmap,template_geom)
IMPLICIT NONE

type(intersect_type) :: template_geom
INTEGER_T :: gridmap(2,2)
INTEGER_T :: mapped_nodes(maxmappednodes)
INTEGER_T :: power2(maxmappednodes)
INTEGER_T :: i,j,nn,n_nodes,inode,checksum

 if (maxmappednodes.lt.4) then
  print *,"maxmappednodes invalid"
  stop
 endif

 power2(1)=1
 do i=2,maxmappednodes
  power2(i)=2*power2(i-1)
 enddo

 do i=1,maxmappednodes
  mapped_nodes(i)=0
 enddo

 n_nodes=template_geom%n_nodes
 if ((n_nodes.lt.1).or.(n_nodes.gt.maxnodelist)) then
  print *,"n_nodes invalid"
  stop
 endif

 inode=1
 checksum=0
 do j=1,2
 do i=1,2
  mapped_nodes(inode)=gridmap(i,j)
  if ((mapped_nodes(inode).lt.1).or.(mapped_nodes(inode).gt.4)) then
   print *,"mapped nodes bust"
   stop
  endif

  do nn=1,n_nodes
   if ((template_geom%nodelistmap(nn,2).eq.inode).and. &
       (template_geom%nodelistmap(nn,3).eq.inode)) then
    checksum=checksum+power2(mapped_nodes(inode))
   endif
  enddo

  inode=inode+1
 enddo
 enddo

 if ((checksum.lt.1).or.(checksum.gt.15)) then
  print *,"checksum invalid"
  stop
 endif
 if (template_geom%n_nodes.le.0) then
  print *,"n_nodes bust"
  stop
 endif
 if (1.eq.0) then
  print *,"init rectangle checksum,map(1) ",checksum,mapped_nodes(1)
  print *,"n_pos_nodes ",template_geom%n_pos_nodes
 endif

 call copy_intersect_type(template_geom, &
   rectangle_maps(checksum,mapped_nodes(1))%intersect_geometry)
 do i=1,maxmappednodes
  rectangle_maps(checksum,mapped_nodes(1))%mapped_nodes(i)=mapped_nodes(i)
 enddo

return
end subroutine add_to_rec

subroutine init_intersect_type(template_geom,n_nodes,n_faces, &
  n_capfaces,n_pos_nodes,node_array,face_array, &
  aligned_array,capface_array,sdim)
IMPLICIT NONE

INTEGER_T :: n_nodes,n_faces,n_capfaces,n_pos_nodes,sdim
INTEGER_T :: node_array(n_nodes)
INTEGER_T :: face_array(3*n_faces)
INTEGER_T :: aligned_array(n_faces)
INTEGER_T :: capface_array(3*n_capfaces)
INTEGER_T :: i,j,icomp,dir,firstnode,secondnode
INTEGER_T :: rawnode,found,itet

type(intersect_type), intent(out) :: template_geom

 if ((n_nodes.lt.1).or.(n_nodes.gt.maxnodelist)) then
  print *,"n nodes invalid"
  stop
 endif
 if ((n_faces.lt.1).or.(n_faces.gt.maxfacelist)) then
  print *,"n faces invalid"
  stop
 endif
 if ((n_capfaces.lt.0).or.(n_capfaces.gt.maxcapfacelist)) then
  print *,"n faces invalid"
  stop
 endif
 if ((n_pos_nodes.lt.1).or.(n_pos_nodes.gt.n_nodes)) then
  print *,"n nodes invalid"
  stop
 endif


 template_geom%n_nodes=n_nodes
 template_geom%n_faces=n_faces
 template_geom%n_capfaces=n_capfaces
 template_geom%n_pos_nodes=n_pos_nodes
 template_geom%n_tet=0

 if ((sdim.ne.2).and.(sdim.ne.3)) then
  print *,"sdim invalid"
  stop
 endif

 if (node_array(1).ne.11) then
  print *,"node_array(1) should be 11"
  stop
 endif

 do i=1,n_nodes
  template_geom%nodelist(i)=node_array(i)
  template_geom%nodelistmap(i,1)=node_array(i)
  firstnode=node_array(i)/10
  secondnode=node_array(i)-10*firstnode
  template_geom%nodelistmap(i,2)=firstnode
  template_geom%nodelistmap(i,3)=secondnode
 enddo

 icomp=1
 do i=1,n_faces
  do dir=1,3

   template_geom%facelist(i,dir)=0

   if (dir.le.sdim) then

    rawnode=face_array(icomp)
    found=0
    do j=1,n_nodes
     if (rawnode.eq.node_array(j)) then
      found=1
      template_geom%facelist(i,dir)=j
     endif
    enddo
    if (found.ne.1) then
     print *,"could not find node in list"
     stop
    endif

   endif ! dir<=sdim

   icomp=icomp+1
  enddo  ! dir
 enddo 

 icomp=1
 do i=1,n_capfaces
  do dir=1,3

   template_geom%capfacelist(i,dir)=0

   if (dir.le.sdim) then

    rawnode=capface_array(icomp)
    found=0
    do j=1,n_nodes
     if (rawnode.eq.node_array(j)) then
      found=1
      template_geom%capfacelist(i,dir)=j
     endif
    enddo
    if (found.ne.1) then
     print *,"could not find node in list"
     stop
    endif

   endif ! dir<=sdim

   icomp=icomp+1
  enddo  ! dir
 enddo ! i

 if (n_nodes.gt.0) then
  template_geom%n_tet=0
  itet=0

   ! first add all the tets with a capface base
  do i=1,n_capfaces
   itet=itet+1
   if ((itet.lt.1).or.(itet.gt.maxtetlist)) then
    print *,"itet invalid"
    stop
   endif
   template_geom%tetlist(itet,1)=1
   do dir=1,3
    template_geom%tetlist(itet,dir+1)=0
   enddo
   do dir=1,sdim
    template_geom%tetlist(itet,dir+1)=template_geom%capfacelist(i,dir)
   enddo
  enddo

  do i=1,n_faces
   template_geom%aligned(i)=aligned_array(i)
   if (aligned_array(i).eq.0) then
    itet=itet+1
    if ((itet.lt.1).or.(itet.gt.maxtetlist)) then
     print *,"itet invalid"
     stop
    endif
    template_geom%tetlist(itet,1)=1
    do dir=1,3
     template_geom%tetlist(itet,dir+1)=0
    enddo
    do dir=1,sdim
     template_geom%tetlist(itet,dir+1)=template_geom%facelist(i,dir)
    enddo
   else if (aligned_array(i).ne.1) then
    print *,"aligned array invalid"
    stop
   endif
  enddo  ! i=1,n_faces

  template_geom%n_tet=itet
 else 
  print *,"n_nodes invalid"
  stop
 endif

return 
end subroutine init_intersect_type
 

subroutine init_geometry_tables()
IMPLICIT NONE

INTEGER_T i,j,k,n,sdim,n_nodes,n_faces,n_capfaces,n_pos_nodes
INTEGER_T ii,jj,kk
INTEGER_T node_array3(3)
INTEGER_T node_array4(4)
INTEGER_T node_array5(5)
INTEGER_T node_array6(6)
INTEGER_T node_array8(8)
INTEGER_T node_array10(10)
INTEGER_T node_array11(11)
INTEGER_T face_array2(6)
INTEGER_T face_array3(9)
INTEGER_T face_array4(12)
INTEGER_T face_array6(18)
INTEGER_T face_array7(21)
INTEGER_T face_array10(30)
INTEGER_T face_array12(36)
INTEGER_T face_array13(39)
INTEGER_T face_array14(42)
INTEGER_T face_array15(45)
INTEGER_T capface_array1(3)
INTEGER_T capface_array2(6)
INTEGER_T capface_array3(9)
INTEGER_T capface_array4(12)
INTEGER_T aligned_array2(2)
INTEGER_T aligned_array3(3)
INTEGER_T aligned_array4(4)
INTEGER_T aligned_array6(6)
INTEGER_T aligned_array7(7)
INTEGER_T aligned_array10(10)
INTEGER_T aligned_array12(12)
INTEGER_T aligned_array13(13)
INTEGER_T aligned_array14(14)
INTEGER_T aligned_array15(15)

INTEGER_T inode
INTEGER_T basegridmap(2,2,2)
INTEGER_T basegridmap2d(2,2)
INTEGER_T baselinemap(4)
INTEGER_T rotmap(2,2,2,0:3)
INTEGER_T rotmap2d(2,2,0:3)
INTEGER_T flipx,flipy,flipz
INTEGER_T rotatex,rotatey,rotatez
INTEGER_T i1,i2,i3,i4,npos,nrot

 do i=1,256
 do j=1,8
  hexahedron_maps(i,j)%intersect_geometry%n_nodes=0
  do n=1,maxmappednodes
   hexahedron_maps(i,j)%mapped_nodes(n)=0
  enddo
 enddo
 enddo
 do i=1,16
 do j=1,4
  tetrahedron_maps(i,j)%intersect_geometry%n_nodes=0
  rectangle_maps(i,j)%intersect_geometry%n_nodes=0
  do n=1,maxmappednodes
   tetrahedron_maps(i,j)%mapped_nodes(n)=0
   rectangle_maps(i,j)%mapped_nodes(n)=0
  enddo
 enddo
 enddo
 do i=1,8
 do j=1,3
  triangle_maps(i,j)%intersect_geometry%n_nodes=0
  do n=1,maxmappednodes
   triangle_maps(i,j)%mapped_nodes(n)=0
  enddo
 enddo
 enddo
 
 sdim=3
   
 n_nodes=4
 n_faces=3
 n_capfaces=1
 n_pos_nodes=1
 node_array4=(/ 11,12,13,15 /)
 face_array3=(/ 11,13,15, 11,12,15, 11,12,13 /)
 aligned_array3=(/ 1,1,1 /)  
 capface_array1=(/ 12,13,15 /)

 call init_intersect_type(template_hex_plane(n_pos_nodes),n_nodes,n_faces, &
   n_capfaces,n_pos_nodes,node_array4,face_array3, &
   aligned_array3,capface_array1,sdim)

 n_nodes=6
 n_faces=6
 n_capfaces=2
 n_pos_nodes=2
 node_array6=(/ 11,22,13,15,24,26 /)
 face_array6=(/ 22,26,24, 11,15,13, 11,22,24, 11,13,24, &
  11,26,15, 11,26,22 /)
 aligned_array6=(/ 0,1,1,1,1,1 /)
 capface_array2=(/ 26,15,24, 15,24,13 /)

 call init_intersect_type(template_hex_plane(n_pos_nodes),n_nodes,n_faces, &
   n_capfaces,n_pos_nodes,node_array6,face_array6, &
   aligned_array6,capface_array2,sdim)

 n_nodes=8
 n_faces=6
 n_capfaces=3
 n_pos_nodes=3
 node_array8=(/ 11,22,33,15,24,34,26,37 /)
 face_array6=(/ 11,26,22, 11,26,15, 11,37,33, 11,37,15, &
   33,34,37, 24,26,22 /)
 aligned_array6=(/ 1,1,1,1,0,0 /)
 capface_array3=(/ 24,26,15, 24,34,15, 34,37,15 /)

 call init_intersect_type(template_hex_plane(n_pos_nodes),n_nodes,n_faces, &
   n_capfaces,n_pos_nodes,node_array8,face_array6, &
   aligned_array6,capface_array3,sdim)

 n_nodes=8
 n_faces=10
 n_capfaces=2
 n_pos_nodes=4
 node_array8=(/ 11,22,33,44,26,15,37,48 /)
 face_array10=(/ 11,26,22, 11,26,15, 11,37,33, 11,37,15, &
                 11,44,22, 11,44,33, 33,48,44, 33,48,37, &
                 22,48,44, 22,48,26 /)
 aligned_array10=(/ 1,1,1,1,1,1,0,0,0,0 /)
 capface_array2=(/ 15,48,26, 15,48,37 /)

 call init_intersect_type(template_hex_plane(n_pos_nodes),n_nodes,n_faces, &
   n_capfaces,n_pos_nodes,node_array8,face_array10, &
   aligned_array10,capface_array2,sdim)

 n_nodes=10
 n_faces=12
 n_capfaces=4
 n_pos_nodes=4
 node_array10=(/ 11,22,33,55,34,37,57,24,26,56 /)
 face_array12=(/ 11,34,33, 11,34,24, 11,24,22, &
     11,37,33, 11,37,57, 11,57,55, &
     11,56,55, 11,56,26, 11,26,22, &
     55,56,57, 33,34,37, 24,22,26 /)
 aligned_array12=(/ 1,1,1,1,1,1,1,1,1,0,0,0 /)
 capface_array4=(/ 34,37,57, 34,57,24, 24,26,57, 57,26,56 /)

   ! 9th element will have this extra shape
 call init_intersect_type(template_hex_plane(9),n_nodes,n_faces, &
   n_capfaces,n_pos_nodes,node_array10,face_array12, &
   aligned_array12,capface_array4,sdim)

 n_nodes=11
 n_faces=13
 n_capfaces=3
 n_pos_nodes=5
 node_array11=(/ 11,22,33,44,55,37,48,26,56,57,37 /)
 face_array13=(/ 11,44,22, 11,44,33, 11,57,55, 11,57,37, 11,37,33, &
                 11,56,55, 11,56,26, 11,26,22, 44,37,48, 44,37,33, &
                 44,26,48, 44,26,22, 55,57,56 /)

 aligned_array13=(/ 1,1,1,1,1,1,1,1,0,0,0,0,0 /)
 capface_array3=(/ 48,57,37, 48,57,56, 48,56,26 /)
 
 call init_intersect_type(template_hex_plane(n_pos_nodes),n_nodes,n_faces, &
   n_capfaces,n_pos_nodes,node_array11,face_array13, &
   aligned_array13,capface_array3,sdim)

 n_nodes=10
 n_faces=14
 n_capfaces=2
 n_pos_nodes=6
 node_array10=(/ 11,22,33,44,55,66,68,48,57,37 /)
 face_array14=(/ 11,44,33, 11,44,22, 11,66,55, 11,66,22, &
                 44,37,33, 44,37,48, 11,37,33, 11,37,57, 11,57,55, &
                 22,48,44, 22,48,68, 22,68,66, 55,68,66, 55,68,57 /)

 aligned_array14=(/ 1,1,1,1,0,0,1,1,1,0,0,0,0,0 /)
 capface_array2=(/ 48,57,37, 48,57,68 /)

 call init_intersect_type(template_hex_plane(n_pos_nodes),n_nodes,n_faces, &
   n_capfaces,n_pos_nodes,node_array10,face_array14, &
   aligned_array14,capface_array2,sdim)

 n_nodes=10
 n_faces=15
 n_capfaces=1
 n_pos_nodes=7
 node_array10=(/ 11,22,33,44,55,66,77,68,48,78 /)
 face_array15=(/ 11,44,22, 11,44,33, 55,68,66, 55,68,78, 55,78,77, &
                 11,66,55, 11,66,22, 11,77,55, 11,77,33, &
                 22,68,66, 22,68,48, 22,48,44, &
                 33,78,77, 33,78,48, 33,48,44 /)

 aligned_array15=(/ 1,1,0,0,0,1,1,1,1,0,0,0,0,0,0 /)
 capface_array1=(/ 68,48,78 /)

 call init_intersect_type(template_hex_plane(n_pos_nodes),n_nodes,n_faces, &
   n_capfaces,n_pos_nodes,node_array10,face_array15, &
   aligned_array15,capface_array1,sdim)

 n_nodes=8
 n_faces=12
 n_capfaces=0
 n_pos_nodes=8
 node_array8=(/ 11,22,33,44,55,66,77,88 /)
 face_array12=(/ 11,44,22, 11,44,33, 55,88,77, 55,88,66, &
   11,66,55, 11,66,22, 33,88,77, 33,88,44,  &
   11,77,33, 11,77,55, 22,88,44, 22,88,66 /)
 aligned_array12=(/ 1,1,0,0,1,1,0,0,1,1,0,0 /)
 capface_array1=(/ 0,0,0 /)

 call init_intersect_type(template_hex_plane(n_pos_nodes),n_nodes,n_faces, &
   n_capfaces,n_pos_nodes,node_array8,face_array12, &
   aligned_array12,capface_array1,sdim)

  ! consider all possible orderings of the 8 nodes relative to each other
  ! that does not change the "shape" of a regular hexahedron.
 do flipx=0,1
 do flipy=0,1
 do flipz=0,1

 do rotatez=0,3
 do rotatey=0,3
 do rotatex=0,3

  inode=1
  do k=1,2
  do j=1,2
  do i=1,2
   basegridmap(i,j,k)=inode
   inode=inode+1
  enddo
  enddo
  enddo

  do k=1,2
  do j=1,2
  do i=1,2
   if (flipx.eq.0) then
    ii=i
   else 
    ii=3-i
   endif
   if (flipy.eq.0) then
    jj=j
   else 
    jj=3-j
   endif
   if (flipz.eq.0) then
    kk=k
   else 
    kk=3-k
   endif
 
 
   rotmap(i,j,k,0)=basegridmap(ii,jj,kk)
  enddo
  enddo
  enddo

  ! rotate about z axis
  do k=1,2
   do nrot=1,rotatez
    rotmap(1,1,k,nrot)=rotmap(2,1,k,nrot-1)
    rotmap(2,1,k,nrot)=rotmap(2,2,k,nrot-1)
    rotmap(2,2,k,nrot)=rotmap(1,2,k,nrot-1)
    rotmap(1,2,k,nrot)=rotmap(1,1,k,nrot-1)
   enddo
  enddo

  do k=1,2
  do j=1,2
  do i=1,2
   rotmap(i,j,k,0)=rotmap(i,j,k,rotatez)
  enddo
  enddo
  enddo

  ! rotate about y axis
  do j=1,2
   do nrot=1,rotatey
    rotmap(1,j,1,nrot)=rotmap(2,j,1,nrot-1)
    rotmap(2,j,1,nrot)=rotmap(2,j,2,nrot-1)
    rotmap(2,j,2,nrot)=rotmap(1,j,2,nrot-1)
    rotmap(1,j,2,nrot)=rotmap(1,j,1,nrot-1)
   enddo
  enddo


  do k=1,2
  do j=1,2
  do i=1,2
   rotmap(i,j,k,0)=rotmap(i,j,k,rotatey)
  enddo
  enddo
  enddo

  ! rotate about x axis
  do i=1,2
   do nrot=1,rotatex
    rotmap(i,1,1,nrot)=rotmap(i,2,1,nrot-1)
    rotmap(i,2,1,nrot)=rotmap(i,2,2,nrot-1)
    rotmap(i,2,2,nrot)=rotmap(i,1,2,nrot-1)
    rotmap(i,1,2,nrot)=rotmap(i,1,1,nrot-1)
   enddo
  enddo

  do k=1,2
  do j=1,2
  do i=1,2
   basegridmap(i,j,k)=rotmap(i,j,k,rotatex)
  enddo
  enddo
  enddo

  do npos=1,9
   call add_to_hex(basegridmap,template_hex_plane(npos))
  enddo

 enddo
 enddo
 enddo  ! rotatex,y,z

 enddo
 enddo
 enddo  ! flipx,y,z

  ! now initialize the tetrahedron intersection table

 sdim=3

 n_nodes=4
 n_faces=3
 n_capfaces=1
 n_pos_nodes=1
 node_array4=(/ 11,12,13,14 /)
 face_array3=(/ 11,13,14, 11,12,14, 11,12,13 /)
 aligned_array3=(/ 1,1,1 /)  
 capface_array1=(/ 12,13,14 /)
 
 call init_intersect_type(template_tet_plane(n_pos_nodes),n_nodes,n_faces, &
   n_capfaces,n_pos_nodes,node_array4,face_array3, &
   aligned_array3,capface_array1,sdim)

 n_nodes=6
 n_faces=6
 n_capfaces=2
 n_pos_nodes=2
 node_array6=(/ 11,22,14,13,24,23 /)
 face_array6=(/ 11,24,22, 11,24,14, 11,23,22, 11,23,13, &
    11,13,14, 22,23,24 /)
 aligned_array6=(/ 1,1,1,1,1,0 /)  
 capface_array2=(/ 14,23,13, 14,23,24 /)

 call init_intersect_type(template_tet_plane(n_pos_nodes),n_nodes,n_faces, &
   n_capfaces,n_pos_nodes,node_array6,face_array6, &
   aligned_array6,capface_array2,sdim)

 n_nodes=6
 n_faces=7
 n_capfaces=1
 n_pos_nodes=3
 node_array6=(/ 11,22,33,14,34,24 /)
 face_array7=(/ 11,22,33, 22,34,33, 22,34,24, &
   11,34,33, 11,34,14, 11,24,22, 11,24,14 /)
 aligned_array7=(/ 1,0,0,1,1,1,1 /)
 capface_array1=(/ 14,34,24 /)

 call init_intersect_type(template_tet_plane(n_pos_nodes),n_nodes,n_faces, &
   n_capfaces,n_pos_nodes,node_array6,face_array7, &
   aligned_array7,capface_array1,sdim)

 n_nodes=4
 n_faces=4
 n_capfaces=0
 n_pos_nodes=4
 node_array4=(/ 11,22,33,44 /)
 face_array4=(/ 11,22,33, 11,22,44, 11,33,44, 22,33,44 /)
 aligned_array4=(/ 1,1,1,0 /)
 capface_array1=(/ 0,0,0 /)

 call init_intersect_type(template_tet_plane(n_pos_nodes),n_nodes,n_faces, &
   n_capfaces,n_pos_nodes,node_array4,face_array4, &
   aligned_array4,capface_array1,sdim)

  ! all possible orderings of the nodes with respect to each other
  ! does not change the "shape" of a tetrahedra.
 do i1=1,4
 do i2=1,4
 do i3=1,4
 do i4=1,4
  baselinemap(1)=i1
  if (i2.ne.i1) then
   baselinemap(2)=i2
   if ((i3.ne.i2).and.(i3.ne.i1)) then
    baselinemap(3)=i3
    if ((i4.ne.i3).and.(i4.ne.i2).and.(i4.ne.i1)) then
     baselinemap(4)=i4

     do npos=1,4
      call add_to_tet(baselinemap,template_tet_plane(npos))
     enddo
    endif
   endif
  endif
 enddo
 enddo
 enddo
 enddo  ! i1,i2,i3,i4
    
   ! now initialize the rectangle intersection table

 sdim=2

 n_nodes=3
 n_faces=2
 n_capfaces=1
 n_pos_nodes=1
 node_array3=(/ 11,13,12 /)
 face_array2=(/ 11,13,0, 11,12,0 /)
 aligned_array2=(/ 1,1 /)
 capface_array1=(/ 12,13,0 /)

 call init_intersect_type(template_rec_plane(n_pos_nodes),n_nodes,n_faces, &
   n_capfaces,n_pos_nodes,node_array3,face_array2, &
   aligned_array2,capface_array1,sdim)

 n_nodes=4
 n_faces=3
 n_capfaces=1
 n_pos_nodes=2
 node_array4=(/ 11,22,13,24 /)
 face_array3=(/ 11,22,0, 11,13,0, 22,24,0 /)
 aligned_array3=(/ 1,1,0 /)
 capface_array1=(/ 13,24,0 /) 

 call init_intersect_type(template_rec_plane(n_pos_nodes),n_nodes,n_faces, &
   n_capfaces,n_pos_nodes,node_array4,face_array3, &
   aligned_array3,capface_array1,sdim)

 n_nodes=5
 n_faces=4
 n_capfaces=1
 n_pos_nodes=3
 node_array5=(/ 11,22,33,34,24 /)
 face_array4=(/ 11,22,0, 11,33,0, 33,34,0, 22,24,0 /)
 aligned_array4=(/ 1,1,0,0 /)
 capface_array1=(/ 34,24,0 /)  

 call init_intersect_type(template_rec_plane(n_pos_nodes),n_nodes,n_faces, &
   n_capfaces,n_pos_nodes,node_array5,face_array4, &
   aligned_array4,capface_array1,sdim)

 n_nodes=4
 n_faces=4
 n_capfaces=0
 n_pos_nodes=4
 node_array4=(/ 11,22,33,44 /)
 face_array4=(/ 11,22,0, 22,44,0, 33,44,0, 11,33,0 /)
 aligned_array4=(/ 1,0,0,1 /)
 capface_array1=(/ 0,0,0 /)

 call init_intersect_type(template_rec_plane(n_pos_nodes),n_nodes,n_faces, &
   n_capfaces,n_pos_nodes,node_array4,face_array4, &
   aligned_array4,capface_array1,sdim)


  ! consider all possible orderings of the 4 nodes relative to each other
  ! that does not change the "shape" of a rectangle.
  !   3 4
  !   1 2
 do flipx=0,1
 do flipy=0,1

 do rotatez=0,3

  inode=1
  do j=1,2
  do i=1,2
   basegridmap2d(i,j)=inode
   inode=inode+1
  enddo
  enddo

  do j=1,2
  do i=1,2
   if (flipx.eq.0) then
    ii=i
   else 
    ii=3-i
   endif
   if (flipy.eq.0) then
    jj=j
   else 
    jj=3-j
   endif
 
   rotmap2d(i,j,0)=basegridmap2d(ii,jj)
  enddo
  enddo

  ! rotate about z axis
  do nrot=1,rotatez
    rotmap2d(1,1,nrot)=rotmap2d(2,1,nrot-1)
    rotmap2d(2,1,nrot)=rotmap2d(2,2,nrot-1)
    rotmap2d(2,2,nrot)=rotmap2d(1,2,nrot-1)
    rotmap2d(1,2,nrot)=rotmap2d(1,1,nrot-1)
  enddo

  do j=1,2
  do i=1,2
   basegridmap2d(i,j)=rotmap2d(i,j,rotatez)
  enddo
  enddo

  do npos=1,4
   call add_to_rec(basegridmap2d,template_rec_plane(npos))
  enddo

 enddo  ! rotatez

 enddo
 enddo  ! flipx,y


  ! now initialize the triangle intersection table

 sdim=2

 n_nodes=3
 n_faces=2
 n_capfaces=1
 n_pos_nodes=1
 node_array3=(/ 11,12,13 /)
 face_array2=(/ 11,13,0, 11,12,0 /)
 aligned_array2=(/ 1,1 /)  
 capface_array1=(/ 12,13,0 /)
 
 call init_intersect_type(template_tri_plane(n_pos_nodes),n_nodes,n_faces, &
   n_capfaces,n_pos_nodes,node_array3,face_array2, &
   aligned_array2,capface_array1,sdim)

 n_nodes=4
 n_faces=3
 n_capfaces=1
 n_pos_nodes=2
 node_array4=(/ 11,22,13,23 /)
 face_array3=(/ 11,22,0, 11,13,0, 22,23,0 /)
 aligned_array3=(/ 1,1,0 /)  
 capface_array1=(/ 13,23,0 /)

 call init_intersect_type(template_tri_plane(n_pos_nodes),n_nodes,n_faces, &
   n_capfaces,n_pos_nodes,node_array4,face_array3, &
   aligned_array3,capface_array1,sdim)

 n_nodes=3
 n_faces=3
 n_capfaces=0
 n_pos_nodes=3
 node_array3=(/ 11,22,33 /)
 face_array3=(/ 11,22,0, 11,33,0, 22,33,0 /)
 aligned_array3=(/ 1,1,0 /)      
 capface_array1=(/ 0,0,0 /)

 call init_intersect_type(template_tri_plane(n_pos_nodes),n_nodes,n_faces, &
   n_capfaces,n_pos_nodes,node_array3,face_array3, &
   aligned_array3,capface_array1,sdim)

  ! all possible orderings of the nodes with respect to each other
  ! does not change the "shape" of a triangle.
 do i1=1,3
 do i2=1,3
 do i3=1,3
  baselinemap(1)=i1
  if (i2.ne.i1) then
   baselinemap(2)=i2
   if ((i3.ne.i2).and.(i3.ne.i1)) then
    baselinemap(3)=i3

    do npos=1,3
     call add_to_tri(baselinemap,template_tri_plane(npos))
    enddo
   endif
  endif
 enddo
 enddo
 enddo  ! i1,i2,i3

return
end subroutine init_geometry_tables

subroutine create_xnodelist( &
  n_vol,n_area, &
  cum_volume,cum_area,cum_centroid,cum_areacentroid, &
  xnode,phinode,checksum,maxnode, &
  shapeflag,nodedomain,coord,sdim)

IMPLICIT NONE

INTEGER_T sdim

INTEGER_T n_vol,n_area
REAL_T cum_volume,cum_area
REAL_T cum_centroid(sdim)
REAL_T cum_areacentroid(sdim)

INTEGER_T checksum,maxnode,shapeflag,coord,nodedomain
REAL_T xnode(nodedomain,sdim)
REAL_T phinode(nodedomain)
INTEGER_T mapped_nodes(nodedomain)
INTEGER_T n_nodes,i,j,dir,index1,index2,n_tet,n_capfaces
REAL_T x1(sdim),x2(sdim)
REAL_T phi1,phi2
REAL_T xtet(sdim+1,sdim)
REAL_T xtri(sdim,sdim)
REAL_T local_volume,local_area
REAL_T local_centroid(sdim),local_areacentroid(sdim)
REAL_T xnodelist_array(maxnodelist,sdim)
type(intersect_type) :: template_geom

 if (sdim.eq.3) then

  if (shapeflag.eq.0) then
   call fast_copy_intersect_type( &
    hexahedron_maps(checksum,maxnode)%intersect_geometry,template_geom,sdim)
   do i=1,nodedomain
    mapped_nodes(i)=hexahedron_maps(checksum,maxnode)%mapped_nodes(i)
   enddo
  else if (shapeflag.eq.1) then
   call fast_copy_intersect_type( &
    tetrahedron_maps(checksum,maxnode)%intersect_geometry,template_geom,sdim)
   do i=1,nodedomain
    mapped_nodes(i)=tetrahedron_maps(checksum,maxnode)%mapped_nodes(i)
   enddo
  else
   print *,"shapeflag invalid"
   stop
  endif

 else if (sdim.eq.2) then

  if (shapeflag.eq.0) then
   call fast_copy_intersect_type( &
    rectangle_maps(checksum,maxnode)%intersect_geometry,template_geom,sdim)
   do i=1,nodedomain
    mapped_nodes(i)=rectangle_maps(checksum,maxnode)%mapped_nodes(i)
   enddo
  else if (shapeflag.eq.1) then
   call fast_copy_intersect_type( &
    triangle_maps(checksum,maxnode)%intersect_geometry,template_geom,sdim)
   do i=1,nodedomain
    mapped_nodes(i)=triangle_maps(checksum,maxnode)%mapped_nodes(i)
   enddo
  else
   print *,"shapeflag invalid"
   stop
  endif

 else
  print *,"sdim invalid"
  stop
 endif

 n_nodes=template_geom%n_nodes
 if (n_nodes.le.0) then
  print *,"no support for this intersection shape"
  stop
 else
  do i=1,n_nodes
   index1=template_geom%nodelistmap(i,2) 
   index2=template_geom%nodelistmap(i,3) 
   if (index1.eq.index2) then
    if ((index1.lt.1).or.(index1.gt.nodedomain)) then
     print *,"index1 invalid"
     stop
    endif
    do dir=1,sdim
     xnodelist_array(i,dir)=xnode(mapped_nodes(index1),dir)
    enddo
   else if (index1.lt.index2) then
    if ((index1.lt.1).or.(index2.gt.nodedomain)) then
     print *,"index1 or index2 invalid"
     stop
    endif
    do dir=1,sdim
     x1(dir)=xnode(mapped_nodes(index1),dir)
     x2(dir)=xnode(mapped_nodes(index2),dir)
    enddo
    phi1=phinode(mapped_nodes(index1))
    phi2=phinode(mapped_nodes(index2))
    if (((phi1.lt.zero).and.(phi2.lt.zero)).or. &
        ((phi1.ge.zero).and.(phi2.ge.zero))) then
     print *,"phi does not change sign"
     print *,"phi1,phi2 = ",phi1,phi2
     stop
    endif
    do dir=1,sdim
     xnodelist_array(i,dir)= &
      (abs(phi1)*x2(dir)+abs(phi2)*x1(dir))/ &
      (abs(phi1)+abs(phi2))
    enddo 
   else
    print *,"index1 should be < index2"
    stop
   endif  

   if (1.eq.0) then
    print *,"in create xnodelist checksum,maxnode,shapeflag,nodedomain ", &
      checksum,maxnode,shapeflag,nodedomain
    print *,"coord,sdim,n_nodes ",coord,sdim,n_nodes
    print *,"i,index1,index2 ",i,index1,index2 
    print *,"i,mapped(index1),mapped(index2) ",i, &
     mapped_nodes(index1),mapped_nodes(index2) 
    do dir=1,sdim
     print *,"i,dir,xnodelist_array ",i,dir,xnodelist_array(i,dir)
    enddo
   endif
  enddo ! initializing the nodes

  n_tet=template_geom%n_tet
  n_capfaces=template_geom%n_capfaces

  do i=1,n_tet
   n_vol=n_vol+1

   do j=1,sdim+1
    do dir=1,sdim
     xtet(j,dir)=xnodelist_array(template_geom%tetlist(i,j),dir)
    enddo
   enddo

   call tetrahedron_volume(xtet,local_volume,local_centroid,coord,sdim)

   cum_volume=cum_volume+local_volume
   do dir=1,sdim
    cum_centroid(dir)=cum_centroid(dir)+local_volume*local_centroid(dir)
   enddo

   if (1.eq.0) then
    print *,"n_tet,n_vol,cum_volume ",n_tet,n_vol,cum_volume
    do j=1,sdim+1
     print *,"i,j,tetnode ",i,j,template_geom%tetlist(i,j)
     do dir=1,sdim
      print *,"i,j,dir,xtetnode ",i,j,dir, &
       xnodelist_array(template_geom%tetlist(i,j),dir)
     enddo
    enddo
   endif

  enddo ! looping through all tets 
   
    
  do i=1,n_capfaces
   n_area=n_area+1

   do j=1,sdim
    do dir=1,sdim
     xtri(j,dir)=xnodelist_array(template_geom%capfacelist(i,j),dir)
    enddo
   enddo

   call surface_area(xtri,local_area,local_areacentroid,coord,sdim)

   cum_area=cum_area+local_area
   do dir=1,sdim
    cum_areacentroid(dir)=cum_areacentroid(dir)+ &
      local_area*local_areacentroid(dir)
   enddo
  enddo ! looping through all tris
 endif ! n_nodes>0

return
end subroutine create_xnodelist



subroutine increment_volume( &
 n_vol,n_area, &
 cum_volume,cum_area,cum_centroid,cum_areacentroid, &
 phinode,xnode,nodedomain,coord,sdim)

IMPLICIT NONE

INTEGER_T sdim

INTEGER_T n_vol,n_area
REAL_T cum_volume,cum_area
REAL_T cum_centroid(sdim)
REAL_T cum_areacentroid(sdim)

INTEGER_T nodedomain,coord
REAL_T phinode(nodedomain)
REAL_T xnode(nodedomain,sdim)
INTEGER_T maxchecksum,checksum,power2,n
REAL_T phimax
INTEGER_T maxnode,n_nodes,shapeflag

 if (nodedomain.ne.sdim+1) then
  print *,"nodedomain invalid"
  stop
 endif

 if (sdim.eq.2) then
  maxchecksum=7
 else if (sdim.eq.3) then
  maxchecksum=15
 else
  print *,"sdim invalid"
  stop
 endif

 checksum=0 
 phimax=zero
 maxnode=0
 power2=1
 do n=1,nodedomain
  if (phinode(n).ge.zero) then
   if (phinode(n).gt.phimax) then
    maxnode=n
    phimax=phinode(n)
   endif
   checksum=checksum+power2
  endif
  power2=2*power2
 enddo

 if (phimax.le.zero) then   

  ! do nothing

 else 

   if (sdim.eq.3) then
    n_nodes=tetrahedron_maps(checksum,maxnode)%intersect_geometry%n_nodes
    if (n_nodes.le.0) then
     print *,"tetra maps incomplete"
     stop
    endif
   else if (sdim.eq.2) then
    n_nodes=triangle_maps(checksum,maxnode)%intersect_geometry%n_nodes
    if (n_nodes.le.0) then
     print *,"triangle maps incomplete"
     stop
    endif
   else
    print *,"sdim invalid"
    stop
   endif

   shapeflag=1
   call create_xnodelist( &
    n_vol,n_area, &
    cum_volume,cum_area,cum_centroid,cum_areacentroid, &
    xnode,phinode,checksum,maxnode,shapeflag, &
    nodedomain,coord,sdim)

 endif ! phimax>0

return
end subroutine increment_volume

! linearcut=1 => input always intersection with plane
! linearcut=0 => input might be intersection with plane
! linearcut=-1 => always break up initial domain into tets/tris
subroutine intersection_volume( &
  cum_volume,cum_area,cum_centroid,cum_areacentroid, &
  phinode,xnode,nodedomain, &
  coord,sdim,fullelementfast,linearcut)

IMPLICIT NONE

INTEGER_T sdim

INTEGER_T n_vol,n_area
REAL_T cum_volume,cum_area
REAL_T cum_centroid(sdim)
REAL_T cum_areacentroid(sdim)

INTEGER_T nodedomain,coord,fullelementfast,linearcut,bfact
REAL_T phinode(nodedomain)
REAL_T xnode(nodedomain,sdim)
INTEGER_T maxchecksum,power2
REAL_T xsten_grid(-3:3,sdim)
REAL_T dxgrid(sdim)
INTEGER_T shapeflag,symmetry_flag,ntetbox,dir
REAL_T phimax
INTEGER_T checksum,maxnode,n,n_nodes,id,sub_nodedomain
REAL_T xx(sdim+1,sdim)
REAL_T ls(sdim+1)
INTEGER_T nhalf

 nhalf=3

 bfact=1
 symmetry_flag=0
 call get_ntetbox(ntetbox,symmetry_flag,sdim)

 if (sdim.eq.2) then
  if (nodedomain.eq.4) then ! domain is a box
   maxchecksum=15
   shapeflag=0
  else if (nodedomain.eq.3) then  ! domain is a triangle
   maxchecksum=7
   shapeflag=1
  else
   print *,"nodedomain invalid"
   stop
  endif
 else if (sdim.eq.3) then
  if (nodedomain.eq.8) then  ! domain is a regular hexahedron
   maxchecksum=255
   shapeflag=0
  else if (nodedomain.eq.4) then  ! domain is a tetrahedron
   maxchecksum=15
   shapeflag=1
  else
   print *,"nodedomain invalid"
   stop
  endif
 else
  print *,"sdim invalid"
  stop
 endif

 checksum=0 
 phimax=zero
 maxnode=0
 power2=1
 do n=1,nodedomain
  if (phinode(n).ge.zero) then
   if (phinode(n).gt.phimax) then
    maxnode=n
    phimax=phinode(n)
   endif
   checksum=checksum+power2
  endif
  power2=2*power2
 enddo

 cum_volume=zero
 cum_area=zero
 do dir=1,sdim
   cum_centroid(dir)=zero
   cum_areacentroid(dir)=zero
 enddo
 n_vol=0
 n_area=0
 
 if (phimax.le.zero) then   

  ! do nothing

 else 

  if ((checksum.eq.maxchecksum).and. &
      (fullelementfast.eq.1)) then
    if (shapeflag.eq.0) then
     do dir=1,sdim
      xsten_grid(-1,dir)=xnode(1,dir)
      xsten_grid(1,dir)=xnode(nodedomain,dir)
      dxgrid(dir)=xsten_grid(1,dir)-xsten_grid(-1,dir)
      xsten_grid(0,dir)=half*(xsten_grid(-1,dir)+xsten_grid(1,dir))
     enddo  ! dir
     call Box_volumeFAST(bfact,dxgrid,xsten_grid,nhalf, &
       cum_volume,cum_centroid,coord,sdim)
    else if (shapeflag.eq.1) then
     call tetrahedron_volume(xnode,cum_volume,cum_centroid,coord,sdim) 
    else
     print *,"shapeflag invalid"
     stop
    endif
  else
   if (sdim.eq.3) then
    if (shapeflag.eq.0) then
     n_nodes=hexahedron_maps(checksum,maxnode)%intersect_geometry%n_nodes
    else if (shapeflag.eq.1) then
     n_nodes=tetrahedron_maps(checksum,maxnode)%intersect_geometry%n_nodes
     if (n_nodes.le.0) then
      print *,"tetra maps incomplete"
      stop
     endif
    else
     print *,"shapeflag invalid"
     stop
    endif
   else if (sdim.eq.2) then
    if (shapeflag.eq.0) then
     n_nodes=rectangle_maps(checksum,maxnode)%intersect_geometry%n_nodes
    else if (shapeflag.eq.1) then
     n_nodes=triangle_maps(checksum,maxnode)%intersect_geometry%n_nodes
     if (n_nodes.le.0) then
      print *,"triangle maps incomplete"
      stop
     endif
    else
     print *,"shapeflag invalid" 
     stop
    endif
   else
    print *,"sdim invalid"
    stop
   endif

! linearcut=1 => input always intersection with plane
! linearcut=0 => input might be intersection with plane
! linearcut=-1 => always break up initial domain into tets/tris

   if ((n_nodes.eq.0).or.(linearcut.eq.-1)) then

    if ((linearcut.eq.1).and.(1.eq.0)) then
     print *,"WARNING: table entry missing"
     print *,"nodedomain=",nodedomain
     print *,"coord=",coord
     print *,"sdim=",sdim
     print *,"fullelementfast=",fullelementfast
     print *,"checksum,maxnode,phimax ",checksum,maxnode,phimax
     print *,"ntetbox,maxchecksum,shapeflag ",ntetbox,maxchecksum, &
       shapeflag
     do n=1,nodedomain 
      print *,"n,phinode ",n,phinode(n)
      do dir=1,sdim
       print *,"n,dir,xnode ",n,dir,xnode(n,dir)
      enddo
     enddo
    endif

    sub_nodedomain=sdim+1
    do id=1,ntetbox
     call extract_tet(xnode,phinode,xx,ls,id,symmetry_flag,sdim)
     call increment_volume( &
      n_vol,n_area, &
      cum_volume,cum_area,cum_centroid,cum_areacentroid, &
      ls,xx,sub_nodedomain,coord,sdim)
    enddo 

   else if (n_nodes.gt.0) then

    call create_xnodelist( &
     n_vol,n_area, &
     cum_volume,cum_area,cum_centroid,cum_areacentroid, &
     xnode,phinode,checksum,maxnode,shapeflag, &
     nodedomain,coord,sdim)

   else
    print *,"n_nodes invalid"
    stop
   endif 

   if (cum_volume.gt.zero) then
    do dir=1,sdim
     cum_centroid(dir)=cum_centroid(dir)/cum_volume
    enddo
   else
    do dir=1,sdim
     cum_centroid(dir)=zero
    enddo
   endif

   if (cum_area.gt.zero) then
    do dir=1,sdim
     cum_areacentroid(dir)=cum_areacentroid(dir)/cum_area
    enddo
   else
    do dir=1,sdim
     cum_areacentroid(dir)=zero
    enddo
   endif
  endif ! not a full element or a full element that is "tetrahedralized"

 endif ! phimax>0

return
end subroutine intersection_volume


      subroutine angle_to_slope3D(angle,nslope,sdim)
      IMPLICIT NONE

      INTEGER_T sdim
      REAL_T nslope(sdim)
      REAL_T angle(sdim-1)

      if (sdim.ne.3) then
       print *,"sdim bust angle to slope 3d"
       stop
      endif
      nslope(3)=cos(angle(2))
      nslope(1)=sin(angle(2))*cos(angle(1))
      nslope(2)=sin(angle(2))*sin(angle(1))

      return
      end subroutine angle_to_slope3D

      subroutine angle_to_slope2D(angle,nslope,sdim)
      IMPLICIT NONE

      INTEGER_T sdim
      REAL_T nslope(sdim)
      REAL_T angle(1)

      if (sdim.ne.2) then
       print *,"sdim bust angle to slope 2d"
       stop
      endif
      nslope(1)=cos(angle(1))
      nslope(2)=sin(angle(1))

      return
      end subroutine angle_to_slope2D


        ! returns centroid in absolute coordinate system

      subroutine cell_intersection_grid( &
        bfact,dxgrid,xgrid,nhalf, &
        lnode, &
        volumedark,centroiddark, &
        area,areacentroid,coord,volall,cenall,sdim)
      IMPLICIT NONE

      INTEGER_T coord,bfact

      INTEGER_T sdim,nhalf
      REAL_T lnode(4*(sdim-1))
      REAL_T xnode(4*(sdim-1),sdim)
      REAL_T xgrid(-nhalf:nhalf,sdim)
      REAL_T dxgrid(sdim)
      REAL_T volumedark
      REAL_T centroiddark(sdim)
      REAL_T volall,area
      REAL_T cenall(sdim)
      REAL_T centroididdark(sdim)
      REAL_T areacentroid(sdim)
      REAL_T areacentroiddark(sdim)
      REAL_T xx(sdim+1,sdim)
      REAL_T lsdark(sdim+1)
      REAL_T volumeiddark
      REAL_T areaiddark
      INTEGER_T i,j,k
      INTEGER_T id,dir,symmetry_flag,ntetbox
      INTEGER_T power2,inode,nplus,nminus,ndark,sumdark

      if (nhalf.lt.1) then
       print *,"nhalf invalid cell_intersection_grid"
       stop
      endif

      if ((sdim.ne.3).and.(sdim.ne.2)) then
       print *,"sdim bust cell intersection grid"
       stop
      endif
      if (bfact.lt.1) then
       print *,"bfact invalid"
       stop
      endif
 
      do j=1,sdim
       centroiddark(j)=zero
       areacentroiddark(j)=zero
       areacentroid(j)=zero
      enddo
      volumedark=zero
      area=zero

      inode=1
      power2=1
      nplus=0
      nminus=0
      sumdark=0
      ndark=0

      call Box_volumeFAST(bfact,dxgrid,xgrid,nhalf,volall,cenall,coord,sdim)

      if (volall.gt.zero) then
    
       symmetry_flag=0
       call get_ntetbox(ntetbox,symmetry_flag,sdim)
 
       if (sdim.eq.3) then
        do k=-1,1,2
        do j=-1,1,2
        do i=-1,1,2
         xnode(inode,1)=xgrid(i,1)
         xnode(inode,2)=xgrid(j,2)
         xnode(inode,sdim)=xgrid(k,sdim)

         if (lnode(inode).ge.zero) then
          nplus=nplus+1
          sumdark=sumdark+power2
          ndark=ndark+1
         endif
         if (lnode(inode).le.zero) then
          nminus=nminus+1
         endif

         inode=inode+1
         power2=power2*2
        enddo
        enddo
        enddo
       else if (sdim.eq.2) then
        do j=-1,1,2
        do i=-1,1,2
         xnode(inode,1)=xgrid(i,1)
         xnode(inode,2)=xgrid(j,2)

         if (lnode(inode).ge.zero) then
          nplus=nplus+1
          sumdark=sumdark+power2
          ndark=ndark+1
         endif
         if (lnode(inode).le.zero) then
          nminus=nminus+1
         endif

         inode=inode+1
         power2=power2*2
        enddo
        enddo
       else 
        print *,"sdim invalid"
        stop
       endif

       do id=1,ntetbox

        call extract_tet(xnode,lnode,xx,lsdark,id,symmetry_flag,sdim)

        if (sdim.eq.3) then
         call intersection_volumeXYZ(lsdark,xx,volumeiddark, &
           centroididdark,areaiddark,areacentroiddark,coord,sdim)
        else if (sdim.eq.2) then
         call int_volumeXYorRZ(lsdark,xx,volumeiddark, &
           centroididdark,areaiddark,areacentroiddark,coord,sdim)
        else
         print *,"sdim invalid"
         stop
        endif

        volumedark=volumedark+volumeiddark
        area=area+areaiddark
        do j=1,sdim
         centroiddark(j)=centroiddark(j)+centroididdark(j)*volumeiddark
         areacentroid(j)=areacentroid(j)+areacentroiddark(j)*areaiddark
        enddo

       enddo ! id
      else if (volall.eq.zero) then
       print *,"WARNING:volall cannot be zero in cell_intersection_grid"
      else if (volall.lt.zero) then
       print *,"volall invalid"
       stop
      endif

      if (area.gt.zero) then
       do j=1,sdim
        areacentroid(j)=areacentroid(j)/area
       enddo
      else
       do j=1,sdim
        areacentroid(j)=zero
       enddo
      endif

      if (volumedark.gt.zero) then
       do j=1,sdim
        centroiddark(j)=centroiddark(j)/volumedark
       enddo
      else
       volumedark=zero
       do j=1,sdim
        centroiddark(j)=zero
       enddo
      endif

      return
      end subroutine cell_intersection_grid


! finds the volume of intersection between the region phi>=0 and
! a column with triangular cross section.
! Also, obtain the centroid and surface area of intersection.

      subroutine int_volumeXYorRZ(phi,x,volumedark, &
       centroiddark,area,areacentroid,coord,sdim)
      IMPLICIT NONE

      INTEGER_T coord,sdim

      REAL_T xtrilist(sdim+1,sdim,MAXTET)
      REAL_T xarealist(sdim,sdim,MAXAREA)
      REAL_T phi(sdim+1)
      REAL_T x(sdim+1,sdim)
      REAL_T xint(sdim+1,sdim)
      REAL_T volumedark,area
      REAL_T areacentroid(sdim)
      REAL_T areacentroidlist(sdim)
      REAL_T volumelistdark,arealist
      REAL_T centroiddark(sdim)
      REAL_T centroidlistdark(sdim)

      INTEGER_T i,j,n,nlist,narea

      if (sdim.ne.2) then
       print *,"sdim invalid"
       stop
      endif

      call list_tris(phi,x,xtrilist,nlist,xarealist,narea,sdim)

      volumedark=zero
      area=zero
      do j=1,sdim
       centroiddark(j)=zero
       areacentroid(j)=zero
      enddo
      do n=1,nlist
       do i=1,sdim+1
       do j=1,sdim
        xint(i,j)=xtrilist(i,j,n)
       enddo
       enddo
       call tetrahedron_volume(xint,volumelistdark,centroidlistdark,coord,sdim)
       volumedark=volumedark+volumelistdark
       do j=1,sdim
        centroiddark(j)=centroiddark(j)+centroidlistdark(j)*volumelistdark
       enddo
      enddo
      if (volumedark.gt.zero) then
       do j=1,sdim
        centroiddark(j)=centroiddark(j)/volumedark
       enddo
      else
       do j=1,sdim
        centroiddark(j)=zero
       enddo
      endif

      do n=1,narea
       do i=1,sdim
       do j=1,sdim
        xint(i,j)=xarealist(i,j,n)
       enddo
       enddo
       do j=1,sdim
        xint(sdim+1,j)=zero
       enddo
       call areaXYorRZ(xint,1,2,arealist,areacentroidlist,coord)
       area=area+arealist
       do j=1,sdim
        areacentroid(j)=areacentroid(j)+areacentroidlist(j)*arealist
       enddo
      enddo

      if (area.gt.zero) then
       do j=1,sdim
        areacentroid(j)=areacentroid(j)/area
       enddo
      else
       do j=1,sdim
        areacentroid(j)=zero
       enddo
      endif

      return
      end subroutine int_volumeXYorRZ


! find triangles that make up the intersection of a plane (phi>0) with
! a triangle
      subroutine list_tris(phi,x,xtrilist,nlist,xarealist,narea,sdim)
      IMPLICIT NONE 
   
      INTEGER_T sdim 
      REAL_T phi(sdim+1)
      REAL_T x(sdim+1,sdim)
      REAL_T xtrilist(sdim+1,sdim,MAXTET)
      REAL_T xarealist(2,2,MAXAREA)
      INTEGER_T nlist,i,j,narea

      if (sdim.ne.2) then
       print *,"sdim invalid list_tris"
       stop
      endif

      nlist=0
      narea=0
      if ((phi(1).le.zero).and.(phi(2).le.zero).and. &
          (phi(3).le.zero)) then
       nlist=0
      else if ((phi(1).ge.zero).and.(phi(2).ge.zero).and. &
               (phi(3).ge.zero)) then
       nlist=1
       do i=1,3
       do j=1,2
        xtrilist(i,j,1)=x(i,j)
       enddo 
       enddo 
      else if ((phi(1).lt.zero).and.(phi(2).ge.zero).and. &
               (phi(3).ge.zero)) then
       call shrink_list_tri(phi,x,1,2,3,1,xtrilist,nlist,xarealist,narea)
      else if ((phi(2).lt.zero).and.(phi(1).ge.zero).and. &
               (phi(3).ge.zero)) then
       call shrink_list_tri(phi,x,2,1,3,1,xtrilist,nlist,xarealist,narea)
      else if ((phi(3).lt.zero).and.(phi(1).ge.zero).and. &
               (phi(2).ge.zero)) then
       call shrink_list_tri(phi,x,3,1,2,1,xtrilist,nlist,xarealist,narea)
      else if ((phi(1).ge.zero).and.(phi(2).lt.zero).and. &
               (phi(3).lt.zero)) then
       call shrink_list_tri(phi,x,1,2,3,0,xtrilist,nlist,xarealist,narea)
      else if ((phi(2).ge.zero).and.(phi(1).lt.zero).and. &
               (phi(3).lt.zero)) then
       call shrink_list_tri(phi,x,2,1,3,0,xtrilist,nlist,xarealist,narea)
      else if ((phi(3).ge.zero).and.(phi(1).lt.zero).and. &
               (phi(2).lt.zero)) then
       call shrink_list_tri(phi,x,3,1,2,0,xtrilist,nlist,xarealist,narea)
      else
       print *,"bust list_tris"
       print *,"phi : ",phi(1),phi(2),phi(3)
       stop
      endif

      return
      end subroutine list_tris




      subroutine get_xbounds(x,xmin,xmax,npoints,sdim)
      IMPLICIT NONE

      INTEGER_T npoints,sdim
      REAL_T x(npoints,sdim)
      REAL_T xmin(sdim)
      REAL_T xmax(sdim)
      INTEGER_T i,dir

      if ((sdim.ne.2).and.(sdim.ne.3)) then
       print *,"sdim invalid"
       stop
      endif

      do dir=1,sdim
       xmin(dir)=x(1,dir)
       xmax(dir)=x(1,dir)
      enddo
      do i=1,npoints
      do dir=1,sdim
       if (x(i,dir).lt.xmin(dir)) then
        xmin(dir)=x(i,dir)
       endif
       if (x(i,dir).gt.xmax(dir)) then
        xmax(dir)=x(i,dir)
       endif
      enddo
      enddo

      return
      end subroutine get_xbounds

! find the length of a side of a triangle; also find the centroid of the
! side.      
      subroutine areaXYorRZ(x,i1,i2,area,areacentroid,coord)
      IMPLICIT NONE

      INTEGER_T coord
      REAL_T x(3,2)
      REAL_T xx(2,2)
      REAL_T area
      REAL_T areacentroid(2)
      INTEGER_T i1,i2,dir,sdim

      sdim=2

      if ((i1.lt.1).or.(i1.gt.3).or.(i2.lt.1).or.(i2.gt.3).or. &
          (i1.eq.i2)) then
       print *,"index invalid area xy or rz"
       print *,"i1,i2 ",i1,i2
       stop
      endif

      do dir=1,sdim
       xx(1,dir)=x(i1,dir) 
       xx(2,dir)=x(i2,dir) 
      enddo
      call surface_area(xx,area,areacentroid,coord,sdim)

      return
      end subroutine areaXYorRZ

       ! in=1  (1,1,1) node
       ! in=2  (2,1,1) node
       ! in=3  (1,2,1) node
       ! in=4  (2,2,1) node
       ! in=5  (1,1,2) node
       ! in=6  (2,1,2) node
       ! in=7  (1,2,2) node
       ! in=8  (2,2,2) node
       ! in=1
       ! do k=1,2
       ! do j=1,2
       ! do i=1,2      
       !  in=in+1
       !  ....
      subroutine extract_tet(xnode,datanode,xtet,datatet, &
       id,symmetry_flag,sdim)
      IMPLICIT NONE

      INTEGER_T sdim
      REAL_T xtet(sdim+1,sdim)
      REAL_T datatet(sdim+1)
      INTEGER_T id
      INTEGER_T j,k,inode,jnode
      INTEGER_T symmetry_flag

      REAL_T datanode(4*(sdim-1))
      REAL_T xnode(4*(sdim-1),sdim)
      INTEGER_T nodelist(sdim+1)
      REAL_T xavg,davg

      if ((sdim.ne.3).and.(sdim.ne.2)) then
       print *,"sdim bust extract tet"
       stop
      endif

      if (symmetry_flag.eq.0) then

       if (sdim.eq.3) then

        if (id.eq.1) then
         nodelist(1)=3 
         nodelist(2)=1 
         nodelist(3)=2 
         nodelist(4)=5
        else if (id.eq.2) then
         nodelist(1)=5 
         nodelist(2)=6 
         nodelist(3)=2 
         nodelist(4)=8
        else if (id.eq.3) then
         nodelist(1)=3 
         nodelist(2)=5 
         nodelist(3)=7 
         nodelist(4)=8
        else if (id.eq.4) then
         nodelist(1)=3 
         nodelist(2)=4 
         nodelist(3)=2 
         nodelist(4)=8
        else if (id.eq.5) then
         nodelist(1)=5 
         nodelist(2)=8 
         nodelist(3)=3 
         nodelist(4)=2
        else
         print *,"id invalid extract_tet"
         stop
        endif

       else if (sdim.eq.2) then

        if (id.eq.1) then
         nodelist(1)=1 
         nodelist(2)=2 
         nodelist(3)=3
        else if (id.eq.2) then 
         nodelist(1)=2 
         nodelist(2)=3 
         nodelist(3)=4
        else
         print *,"id invalid"
         stop
        endif

       else
        print *,"sdim invalid"
        stop
       endif

       do j=1,sdim+1
        inode=nodelist(j)
        do k=1,sdim
         xtet(j,k)=xnode(inode,k)
        enddo
        datatet(j)=datanode(inode)
       enddo 

      else if (symmetry_flag.eq.1) then

        !  3 4
        !  1 2
       if (sdim.eq.2) then

        if (id.eq.1) then  ! dir=2 side=1
         nodelist(1)=1
         nodelist(2)=2
         nodelist(3)=0
        else if (id.eq.2) then ! dir=1 side=1
         nodelist(1)=1
         nodelist(2)=3
         nodelist(3)=0
        else if (id.eq.3) then ! dir=2 side=2
         nodelist(1)=3
         nodelist(2)=4
         nodelist(3)=0
        else if (id.eq.4) then ! dir=1 side=2
         nodelist(1)=4
         nodelist(2)=2
         nodelist(3)=0
        else
         print *,"id invalid"
         stop
        endif

        do j=1,sdim+1
         inode=nodelist(j)
         do k=1,sdim
          if (inode.eq.0) then
           xavg=zero
           do jnode=1,4
            xavg=xavg+xnode(jnode,k)
           enddo
           xtet(j,k)=xavg/four
          else
           xtet(j,k)=xnode(inode,k)
          endif
         enddo
         if (inode.eq.0) then
          davg=zero
          do jnode=1,4
           davg=davg+datanode(jnode)
          enddo
          datatet(j)=davg/four
         else
          datatet(j)=datanode(inode)
         endif
        enddo  ! j

       else if (sdim.eq.3) then
        print *,"symmetric tetrahedrazation for 3d not complete"
        stop
       else
        print *,"sdim invalid"
        stop
       endif

      else
       print *,"symmetry_flag invalid"
       stop
      endif

      return
      end subroutine extract_tet 



! find area of a face of a tetrahedron
! find the centroid of a face of a tetrahedron
      subroutine areaXYZ(x,i1,i2,i3,area,areacentroid,coord)
      IMPLICIT NONE

      INTEGER_T coord
      REAL_T x(4,3)
      REAL_T xx(3,3)
      REAL_T area
      REAL_T areacentroid(3)
      INTEGER_T i1,i2,i3 
      INTEGER_T dir
      INTEGER_T sdim

      sdim=3

      if ((i1.lt.1).or.(i1.gt.4).or.(i2.lt.1).or.(i2.gt.4).or. &
          (i3.lt.1).or.(i3.gt.4).or.(i1.eq.i2).or.(i1.eq.i3).or. &
          (i2.eq.i3)) then
       print *,"index invalid areaXYZ"
       print *,"i1,i2,i3 ",i1,i2,i3
       stop
      endif

      do dir=1,sdim
       xx(1,dir)=x(i1,dir) 
       xx(2,dir)=x(i2,dir) 
       xx(3,dir)=x(i3,dir) 
      enddo
      call surface_area(xx,area,areacentroid,coord,sdim)

      return
      end subroutine areaXYZ

! find the area of a triangular planar element in 2d or 3d
! THANK YOU JOHN BURKHARDT
      subroutine surface_area(x,area,areacen,coord,sdim)
      use LegendreNodes
      use triangle_fekete_module, only : fekete_degree,fekete_order_num, &
                                         fekete_rule
      IMPLICIT NONE

      INTEGER_T sdim,coord
      REAL_T x(sdim,sdim)
      REAL_T area
      REAL_T areacen(sdim)
      REAL_T xx(sdim-1,sdim)
      INTEGER_T i,j,dir,dircrit,itan,jtan
      REAL_T, dimension(:,:), allocatable :: xy
      REAL_T, dimension(:), allocatable :: w
      REAL_T y1_cross_y2(sdim)
      REAL_T coeff(sdim)
      REAL_T mag,DA
      INTEGER_T rule,degree,order_num  ! rule=1 degree=3
      REAL_T drho,dtheta,dz,rho,jac
      REAL_T dxpos(sdim)

      if ((sdim.ne.3).and.(sdim.ne.2)) then
       print *,"sdim bust surface area"
       stop
      endif

      do i=1,sdim-1
      do dir=1,sdim
       xx(i,dir)=x(i,dir)-x(sdim,dir)
      enddo
      enddo

      if (sdim.eq.3) then
       y1_cross_y2(1)=xx(1,2)*xx(sdim-1,sdim)-xx(sdim-1,2)*xx(1,sdim)
       y1_cross_y2(2)=-(xx(1,1)*xx(sdim-1,sdim)-xx(sdim-1,1)*xx(1,sdim))
       y1_cross_y2(sdim)=xx(1,1)*xx(sdim-1,2)-xx(sdim-1,1)*xx(1,2)
      else if (sdim.eq.2) then
       y1_cross_y2(1)=xx(1,2)
       y1_cross_y2(2)=-xx(1,1)
      else
       print *,"dimension bust"
       stop
      endif
      
      mag=zero
      do dir=1,sdim
       mag=mag+y1_cross_y2(dir)**2
      enddo
      mag=sqrt(mag)

      if (sdim.eq.3) then
       area=half*mag  ! mag is area of parallelogram
      else if (sdim.eq.2) then
       area=mag
      else
       print *,"dimension bust"
       stop
      endif

      do dir=1,sdim
       areacen(dir)=zero
       do i=1,sdim
        areacen(dir)=areacen(dir)+x(i,dir)
       enddo
       areacen(dir)=areacen(dir)/sdim
      enddo ! dir

      if (mag.gt.zero) then
       if (coord.eq.0) then
        ! do nothing
       else if ((coord.eq.1).or.(coord.eq.3)) then

        dircrit=1
        if (abs(y1_cross_y2(2)).gt.abs(y1_cross_y2(dircrit))) then
         dircrit=2
        endif
        if (sdim.eq.3) then
         if (abs(y1_cross_y2(sdim)).gt.abs(y1_cross_y2(dircrit))) then
          dircrit=sdim
         endif
        endif
        if (y1_cross_y2(dircrit).eq.zero) then
         print *,"y1_cross_y2 bust"
         stop
        endif
        do dir=1,sdim
         coeff(dir)=y1_cross_y2(dir)/y1_cross_y2(dircrit)
        enddo

        if (dircrit.eq.1) then
         itan=2
         jtan=3
        else if (dircrit.eq.2) then
         itan=1
         jtan=3
        else if ((dircrit.eq.3).and.(sdim.eq.3)) then
         itan=1
         jtan=2
        else
         print *,"dircrit invalid"
         stop
        endif

       else
        print *,"coord invalid"
        stop
       endif

       if (coord.eq.0) then
        ! do nothing
       else if (coord.eq.1) then

        area=zero
        do dir=1,sdim
         areacen(dir)=zero
        enddo
        jac=abs(y1_cross_y2(dircrit))

        if (sdim.ne.2) then
         print *,"dimension bust"
         stop
        endif
        degree=3
        order_num=2
        allocate(xy(sdim-1,order_num))
        allocate(w(order_num))

        do i=1,order_num
         xy(1,i)=(cache_gauss(order_num,i-1)+one)/two
         w(i)=cache_gauss_w(order_num,i-1)*half
        
         dxpos(itan)=xx(1,itan)*xy(1,i)
         dxpos(dircrit)=-(coeff(itan)*dxpos(itan))
         rho=x(sdim,1)+dxpos(1)
        
         ! F(rho,theta,z)=0
         ! area=integral_R |grad F|/|grad F dot p| dR
         ! p is normal to R.   
         ! dR=2 pi rho dz or
         ! dR=2 pi rho drho
         DA=two*Pi*abs(rho)*sqrt( one+coeff(itan)**2 )
 
         area=area+DA*jac*w(i)

         do dir=1,sdim
          areacen(dir)=areacen(dir)+(x(sdim,dir)+dxpos(dir))*DA*jac*w(i)
         enddo
        enddo ! i

        if (area.eq.zero) then
         print *,"area became 0 even though mag>0"
         stop
        endif
        do dir=1,sdim
         areacen(dir)=areacen(dir)/area
        enddo

        deallocate(xy)
        deallocate(w)

       else if (coord.eq.3) then ! in: surface_area

        area=zero
        do dir=1,sdim
         areacen(dir)=zero
        enddo
        jac=abs(y1_cross_y2(dircrit))

        if (sdim.eq.2) then
         degree=3
         order_num=2
         allocate(xy(sdim-1,order_num))
         allocate(w(order_num))

         do i=1,order_num
          xy(1,i)=(cache_gauss(order_num,i-1)+one)/two
          w(i)=cache_gauss_w(order_num,i-1)*half
        
          dxpos(itan)=xx(1,itan)*xy(1,i)
          dxpos(dircrit)=-(coeff(itan)*dxpos(itan))
          rho=x(sdim,1)+dxpos(1)
        
          ! F(rho,theta,z)=0
          ! area=integral_R |grad F|/|grad F dot p| dR
          ! p is normal to R.   
          ! dR=rho dtheta or
          ! dR=drho
          if (dircrit.eq.1) then
           DA=sqrt( rho**2+coeff(itan)**2 )
          else if (dircrit.eq.2) then
           DA=sqrt( one+(coeff(itan)*rho)**2 )
          else
           print *,"dircrit invalid"
           stop
          endif

          area=area+DA*jac*w(i)

          do dir=1,sdim
           areacen(dir)=areacen(dir)+(x(sdim,dir)+dxpos(dir))*DA*jac*w(i)
          enddo
         enddo ! i

         if (area.eq.zero) then
          print *,"area became 0 even though mag>0"
          stop
         endif
         do dir=1,sdim
          areacen(dir)=areacen(dir)/area
         enddo

         deallocate(xy)
         deallocate(w)

        else if (sdim.eq.3) then

         rule=1  ! degree=3
         call fekete_degree(rule,degree)
         if (degree.ne.3) then
          print *,"degree invalid"
          stop
         endif
         call fekete_order_num(rule,order_num)
         allocate(xy(sdim-1,order_num))
         allocate(w(order_num))
         call fekete_rule(rule,order_num,xy,w)

         do i=1,order_num
          dxpos(itan)=xx(1,itan)*xy(1,i)+xx(2,itan)*xy(2,i)
          dxpos(jtan)=xx(1,jtan)*xy(1,i)+xx(2,jtan)*xy(2,i)
          dxpos(dircrit)=-(coeff(itan)*dxpos(itan)+coeff(jtan)*dxpos(jtan))
          rho=x(sdim,1)+dxpos(1)
        
          ! F(rho,theta,z)=0
          ! area=integral_R |grad F|/|grad F dot p| dR
          ! p is normal to R.   
          ! dR=rho dtheta dz or
          ! dR=drho dz or
          ! dR=rho drho dtheta
          if (dircrit.eq.1) then 
           DA=sqrt( (coeff(jtan)*rho)**2+rho**2+coeff(itan)**2 )
          else if (dircrit.eq.2) then
           DA=sqrt( (coeff(jtan)*rho)**2+one+(coeff(itan)*rho)**2 )
          else if (dircrit.eq.3) then
           DA=sqrt( (coeff(jtan))**2+rho**2+(coeff(itan)*rho)**2 )
          else
           print *,"dircrit invalid"
           stop
          endif
 
          area=area+DA*jac*w(i)

          do dir=1,sdim
           areacen(dir)=areacen(dir)+(x(sdim,dir)+dxpos(dir))*DA*jac*w(i)
          enddo
         enddo ! i

         if (area.eq.zero) then
          print *,"area became 0 even though mag>0"
          stop
         endif
         do dir=1,sdim
          areacen(dir)=areacen(dir)/area
         enddo

         area=area*half

         deallocate(xy)
         deallocate(w)
        else
         print *,"dimension bust"
         stop
        endif
       else
        print *,"coord invalid"
        stop
       endif
      else if (mag.eq.zero) then
       ! do nothing
      else
       print *,"mag invalid"
       stop
      endif

      return
      end subroutine surface_area

! find the volume of a tetrahedron (in 3d) or area of triangle (2d)
      subroutine tetrahedron_volume(x,volume,centroid,coord,sdim)
      use tetrahedron_keast_module, only : keast_degree,keast_order_num, &
                                           keast_rule
      use triangle_fekete_module, only : fekete_degree,fekete_order_num, &
                                         fekete_rule

      IMPLICIT NONE

      INTEGER_T sdim,coord
      REAL_T x(sdim+1,sdim)
      REAL_T volume
      REAL_T centroid(sdim)
      REAL_T centroid_def(sdim)
      REAL_T xx(sdim,sdim)
      INTEGER_T i,j,dir
      REAL_T, dimension(:,:), allocatable :: xyz
      REAL_T, dimension(:), allocatable :: w
      INTEGER_T rule,degree,order_num
      REAL_T rho,volzero
      REAL_T xmin(sdim)
      REAL_T xmax(sdim)
      REAL_T dxpos(sdim)
      REAL_T total_weight

      if ((sdim.ne.3).and.(sdim.ne.2)) then
       print *,"sdim invalid"
       stop
      endif

      call get_xbounds(x,xmin,xmax,sdim+1,sdim)
  
      do i=1,sdim
      do dir=1,sdim
       xx(i,dir)=x(i,dir)-x(sdim+1,dir)
      enddo
      enddo

      if (sdim.eq.3) then
       volzero= &
        xx(2,1)*(xx(1,2)*xx(sdim,sdim)-xx(sdim,2)*xx(1,sdim))- &
        xx(1,1)*(xx(2,2)*xx(sdim,sdim)-xx(sdim,2)*xx(2,sdim))+ &
        xx(sdim,1)*(xx(2,2)*xx(1,sdim)-xx(1,2)*xx(2,sdim))
      else if (sdim.eq.2) then
       volzero=xx(1,1)*xx(2,2)-xx(2,1)*xx(1,2)
      else
       print *,"dimension bust"
       stop
      endif
      volzero=abs(volzero)

      do dir=1,sdim
       centroid(dir)=zero
       do i=1,sdim+1
        centroid(dir)=centroid(dir)+x(i,dir)
       enddo
       centroid(dir)=centroid(dir)/(sdim+one)
       centroid_def(dir)=centroid(dir)
      enddo ! dir

      if (coord.eq.0) then

       if (sdim.eq.3) then
        volume=volzero/six
       else if (sdim.eq.2) then
        volume=volzero/two
       else
        print *,"sdim invalid"
        stop
       endif
  
      else if (coord.eq.1) then

       if (sdim.ne.2) then
        print *,"dimension bust"
        stop
       endif

       rule=1  ! degree=3
       call fekete_degree(rule,degree)
       if (degree.ne.3) then
        print *,"degree invalid"
        stop
       endif
       call fekete_order_num(rule,order_num)
       allocate(xyz(sdim,order_num))
       allocate(w(order_num))
       call fekete_rule(rule,order_num,xyz,w)

       if (1.eq.0) then
        print *,"---------------------------"
        total_weight=zero
        do i=1,order_num
         print *,"i,x,y,w ",i,xyz(1,i),xyz(2,i),w(i)
         total_weight=total_weight+w(i)
        enddo
        print *,"total_weight: ",total_weight
        print *,"---------------------------"
       endif

       volume=zero
       do dir=1,sdim
        centroid(dir)=zero
       enddo
       do i=1,order_num
        do dir=1,sdim
         dxpos(dir)=zero
         do j=1,sdim
          dxpos(dir)=dxpos(dir)+xx(j,dir)*xyz(j,i)
         enddo
        enddo
        rho=x(sdim+1,1)+dxpos(1)
        volume=volume+two*Pi*abs(rho)*w(i)
        do dir=1,sdim
         centroid(dir)=centroid(dir)+(x(sdim+1,dir)+dxpos(dir))* &
           two*Pi*abs(rho)*w(i)
        enddo
       enddo ! i=1,order_num
       if (volume.gt.zero) then
        do dir=1,sdim
         centroid(dir)=centroid(dir)/volume
        enddo
       else
        do dir=1,sdim
         centroid(dir)=centroid_def(dir)
        enddo
       endif

       volume=abs(volume)*volzero*half

       deallocate(xyz)
       deallocate(w)

      else if (coord.eq.3) then

       if (sdim.eq.2) then

        rule=1  ! degree=3
        call fekete_degree(rule,degree)
        if (degree.ne.3) then
         print *,"degree invalid"
         stop
        endif
        call fekete_order_num(rule,order_num)
        allocate(xyz(sdim,order_num))
        allocate(w(order_num))
        call fekete_rule(rule,order_num,xyz,w)
        volume=zero
        do dir=1,sdim
         centroid(dir)=zero
        enddo
        do i=1,order_num
         do dir=1,sdim
          dxpos(dir)=zero
          do j=1,sdim
           dxpos(dir)=dxpos(dir)+xx(j,dir)*xyz(j,i)
          enddo
         enddo
         rho=x(sdim+1,1)+dxpos(1)
         volume=volume+abs(rho)*w(i)
         do dir=1,sdim
          centroid(dir)=centroid(dir)+(x(sdim+1,dir)+dxpos(dir))*abs(rho)*w(i)
         enddo
        enddo ! i=1,order_num
        if (volume.gt.zero) then
         do dir=1,sdim
          centroid(dir)=centroid(dir)/volume
         enddo
        else
         do dir=1,sdim
          centroid(dir)=centroid_def(dir)
         enddo
        endif

        volume=abs(volume)*volzero*half

        deallocate(xyz)
        deallocate(w)

       else if (sdim.eq.3) then

        rule=2  ! linear integrands exactly
        call keast_degree(rule,degree)
        if (degree.ne.1) then
         print *,"degree invalid"
         stop
        endif
        call keast_order_num(rule,order_num)
        allocate(xyz(sdim,order_num))
        allocate(w(order_num))
        call keast_rule(rule,order_num,xyz,w)
        volume=zero
        do dir=1,sdim
         centroid(dir)=zero
        enddo
        do i=1,order_num
         do dir=1,sdim
          dxpos(dir)=zero
          do j=1,sdim
           dxpos(dir)=dxpos(dir)+xx(j,dir)*xyz(j,i)
          enddo
         enddo
         rho=x(sdim+1,1)+dxpos(1)
         volume=volume+abs(rho)*w(i)
         do dir=1,sdim
          centroid(dir)=centroid(dir)+(x(sdim+1,dir)+dxpos(dir))* &
            abs(rho)*w(i)
         enddo ! dir
        enddo ! i=1,...,order_num
        if (volume.gt.zero) then
         do dir=1,sdim
          centroid(dir)=centroid(dir)/volume
         enddo
        else
         do dir=1,sdim
          centroid(dir)=centroid_def(dir)
         enddo
        endif

        volume=abs(volume)*volzero*sixth

        deallocate(xyz)
        deallocate(w)

       else
        print *,"dimension bust"
        stop
       endif

      else
       print *,"coord invalid"
       stop
      endif

      do j=1,sdim
       if ((centroid(j)+CENTOL.lt.xmin(j)).or. &
           (centroid(j)-CENTOL.gt.xmax(j))) then
        print *,"centroid invalid XYZ"
        stop
       endif
      enddo

      return
      end subroutine tetrahedron_volume



! internal routine
      subroutine shrink2D(xint,x,phi,isrc,itarg,sdim)
      IMPLICIT NONE

      INTEGER_T sdim
      REAL_T xint(sdim+1,sdim)
      REAL_T x(sdim+1,sdim)
      REAL_T phi(sdim+1)
      INTEGER_T isrc,itarg,j

      if (sdim.ne.2) then
       print *,"sdim bust shring 2d"
       stop
      endif

      if ((itarg.lt.1).or.(itarg.gt.3).or.(isrc.lt.1).or. &
          (isrc.gt.3).or.(itarg.eq.isrc)) then
       print *,"index invalid shrink2d "
       do j=1,3
        print *,"j,xint ",j,xint(j,1),xint(j,2)
        print *,"j,x ",j,x(j,1),x(j,2)
        print *,"j,phi ",j,phi(j)
       enddo
       stop
      endif
      if (((phi(itarg).ge.zero).and.(phi(isrc).ge.zero)).or. &
          ((phi(itarg).lt.zero).and.(phi(isrc).lt.zero))) then
       print *,"levelset does not change sign"
       stop
      endif

      do j=1,2
       xint(itarg,j)=(abs(phi(itarg))*x(isrc,j)+ &
                      abs(phi(isrc))*x(itarg,j))/ &
                     (abs(phi(itarg))+abs(phi(isrc)))
      enddo

      return
      end subroutine shrink2D

! internal routine
      subroutine shrink3D(xint,x,phi,isrc,itarg,sdim)
      IMPLICIT NONE

      INTEGER_T sdim
      REAL_T xint(sdim+1,sdim)
      REAL_T x(sdim+1,sdim)
      REAL_T phi(sdim+1)
      INTEGER_T isrc,itarg,j

      if (sdim.ne.3) then
       print *,"sdim bust shrink 3d"
       stop
      endif

      if ((itarg.lt.1).or.(itarg.gt.4).or.(isrc.lt.1).or. &
          (isrc.gt.4).or.(itarg.eq.isrc)) then
       print *,"index invalid shrink3d "
       do j=1,4
        print *,"j,xint ",j,xint(j,1),xint(j,2),xint(j,3)
        print *,"j,x ",j,x(j,1),x(j,2),x(j,3)
        print *,"j,phi ",j,phi(j)
       enddo
       stop
      endif
      if (((phi(itarg).ge.zero).and.(phi(isrc).ge.zero)).or. &
          ((phi(itarg).lt.zero).and.(phi(isrc).lt.zero))) then
       print *,"levelset does not change sign"
       stop
      endif

      do j=1,3
       xint(itarg,j)=(abs(phi(itarg))*x(isrc,j)+ &
                      abs(phi(isrc))*x(itarg,j))/ &
                     (abs(phi(itarg))+abs(phi(isrc)))
      enddo

      return
      end subroutine shrink3D



 
! internal routine, do not call
      subroutine shrink_list_tri(phi,x,i1,i2,i3,complement, &
        xtrilist,nlist,xarealist,narea)
      IMPLICIT NONE

      REAL_T xtrilist(3,2,MAXTET)
      REAL_T xarealist(2,2,MAXAREA)
      INTEGER_T nlist,narea
      REAL_T phi(3)
      REAL_T x(3,2)
      REAL_T xint(3,2)
      INTEGER_T i,j,i1,i2,i3,complement
      INTEGER_T sdim

      sdim=2
 
      if (complement.eq.0) then
       do i=1,sdim+1
       do j=1,sdim
        xint(i,j)=x(i,j)
       enddo
       enddo
       call shrink2D(xint,x,phi,i1,i2,sdim)
       call shrink2D(xint,x,phi,i1,i3,sdim)
       nlist=nlist+1 
       if (nlist.gt.MAXTET) then
        print *,"nlist invalid"
        stop
       endif 
       do i=1,sdim+1
       do j=1,sdim
        xtrilist(i,j,nlist)=xint(i,j)
       enddo
       enddo

       narea=narea+1
       do j=1,sdim
        xarealist(1,j,narea)=xint(i2,j)
        xarealist(2,j,narea)=xint(i3,j)
       enddo
      else if (complement.eq.1) then
       do i=1,sdim+1
       do j=1,sdim
        xint(i,j)=x(i,j)
       enddo
       enddo
       call shrink2D(xint,x,phi,i2,i1,sdim)
       nlist=nlist+1  
       if (nlist.gt.MAXTET) then
        print *,"nlist invalid"
        stop
       endif 
       do i=1,sdim+1
       do j=1,sdim
        xtrilist(i,j,nlist)=xint(i,j)
       enddo
       enddo

       do i=1,sdim+1
       do j=1,sdim
        xint(i,j)=x(i,j)
       enddo
       enddo
       call shrink2D(xint,x,phi,i1,i2,sdim)
       call shrink2D(xint,x,phi,i3,i1,sdim)
       nlist=nlist+1  
       if (nlist.gt.MAXTET) then
        print *,"nlist invalid"
        stop
       endif 
       do i=1,sdim+1
       do j=1,sdim
        xtrilist(i,j,nlist)=xint(i,j)
       enddo
       enddo

       narea=narea+1
       if (narea.gt.MAXAREA) then
        print *,"nlist invalid"
        stop
       endif 
       do j=1,sdim
        xarealist(1,j,narea)=xint(i2,j)
        xarealist(2,j,narea)=xint(i1,j)
       enddo
      else
       print *,"complement invalid"
       stop
      endif

      return
      end subroutine shrink_list_tri

! internal routine, do not call 
      subroutine shrink_list_tet(phi,x,i1,i2,i3,i4,complement, &
        xtetlist,nlist,xarealist,narea)
      IMPLICIT NONE

      REAL_T xtetlist(4,3,MAXTET)
      REAL_T xarealist(3,3,MAXAREA)
      INTEGER_T nlist,narea
      REAL_T phi(4)
      REAL_T x(4,3)
      REAL_T xint(4,3)
      INTEGER_T i,j,i1,i2,i3,i4,complement
      INTEGER_T sdim
 
      sdim=3
     
      if (complement.eq.0) then
       do i=1,sdim+1
       do j=1,sdim
        xint(i,j)=x(i,j)
       enddo
       enddo
       call shrink3D(xint,x,phi,i1,i2,sdim)  ! xint,x,phi,isrc,itarg
       call shrink3D(xint,x,phi,i1,i3,sdim)  
       call shrink3D(xint,x,phi,i1,i4,sdim)
       nlist=nlist+1  
       if (nlist.gt.MAXTET) then
        print *,"nlist invalid"
        stop
       endif
       do i=1,sdim+1
       do j=1,sdim
        xtetlist(i,j,nlist)=xint(i,j)
       enddo
       enddo
       narea=narea+1
       do j=1,sdim
        xarealist(1,j,narea)=xint(i2,j)
        xarealist(2,j,narea)=xint(i3,j)
        xarealist(3,j,narea)=xint(i4,j)
       enddo
      else if (complement.eq.1) then
       do i=1,sdim+1
       do j=1,sdim
        xint(i,j)=x(i,j)
       enddo
       enddo
       call shrink3D(xint,x,phi,i3,i1,sdim)
       nlist=nlist+1
       if (nlist.gt.MAXTET) then
        print *,"nlist invalid"
        stop
       endif
       do i=1,sdim+1
       do j=1,sdim
        xtetlist(i,j,nlist)=xint(i,j)
       enddo
       enddo

       do i=1,sdim+1
       do j=1,sdim
        xint(i,j)=x(i,j)
       enddo
       enddo
       call shrink3D(xint,x,phi,i4,i1,sdim)
       call shrink3D(xint,x,phi,i1,i3,sdim)
       nlist=nlist+1
       if (nlist.gt.MAXTET) then
        print *,"nlist invalid"
        stop
       endif
       do i=1,sdim+1
       do j=1,sdim
        xtetlist(i,j,nlist)=xint(i,j)
       enddo
       enddo

       do i=1,sdim+1
       do j=1,sdim
        xint(i,j)=x(i,j)
       enddo
       enddo
       call shrink3D(xint,x,phi,i1,i3,sdim)
       call shrink3D(xint,x,phi,i1,i4,sdim)
       call shrink3D(xint,x,phi,i2,i1,sdim)
       nlist=nlist+1
       if (nlist.gt.MAXTET) then
        print *,"nlist invalid"
        stop
       endif
       do i=1,sdim+1
       do j=1,sdim
        xtetlist(i,j,nlist)=xint(i,j)
       enddo
       enddo
       narea=narea+1
       if (narea.gt.MAXAREA) then
        print *,"narea invalid"
        stop
       endif
       do j=1,sdim
        xarealist(1,j,narea)=xint(i3,j)
        xarealist(2,j,narea)=xint(i4,j)
        xarealist(3,j,narea)=xint(i1,j)
       enddo
      else
       print *,"complement invalid"
       stop
      endif

      return
      end subroutine shrink_list_tet

! internal routine, do not call
! i1,i2 nodes are positive
! i3,i4 nodes are negative
      subroutine shrink_gableroof_list(phi,x,i1,i2,i3,i4, &
        xtetlist,nlist,xarealist,narea)
      IMPLICIT NONE

      REAL_T xtetlist(4,3,MAXTET)
      REAL_T xarealist(3,3,MAXAREA)
      INTEGER_T nlist,narea
      REAL_T phi(4)
      REAL_T x(4,3)
      REAL_T xint(4,3)
      INTEGER_T i,j,i1,i2,i3,i4

      INTEGER_T sdim

      sdim=3

      do i=1,sdim+1
      do j=1,sdim
        xint(i,j)=x(i,j)
      enddo
      enddo
      call shrink3D(xint,x,phi,i1,i3,sdim)  ! xint,x,phi,isrc,itarg
      call shrink3D(xint,x,phi,i1,i4,sdim)  
      call shrink3D(xint,x,phi,i4,i2,sdim)  
      nlist=nlist+1
      if (nlist.gt.MAXTET) then
        print *,"nlist invalid"
        stop
      endif
      do i=1,sdim+1
      do j=1,sdim
        xtetlist(i,j,nlist)=xint(i,j)
      enddo
      enddo
      narea=narea+1
      if (narea.gt.MAXAREA) then
        print *,"narea invalid"
        stop
      endif
      do j=1,sdim
       xarealist(1,j,narea)=xint(i2,j)
       xarealist(2,j,narea)=xint(i3,j)
       xarealist(3,j,narea)=xint(i4,j)
      enddo

      do i=1,sdim+1
      do j=1,sdim
        xint(i,j)=x(i,j)
      enddo
      enddo
      call shrink3D(xint,x,phi,i1,i3,sdim)
      call shrink3D(xint,x,phi,i2,i4,sdim)
      call shrink3D(xint,x,phi,i3,i2,sdim)
      nlist=nlist+1
      if (nlist.gt.MAXTET) then
        print *,"nlist invalid"
        stop
      endif
      do i=1,sdim+1
      do j=1,sdim
        xtetlist(i,j,nlist)=xint(i,j)
      enddo
      enddo
      narea=narea+1
      if (narea.gt.MAXAREA) then
        print *,"narea invalid"
        stop
      endif
      do j=1,sdim
       xarealist(1,j,narea)=xint(i2,j)
       xarealist(2,j,narea)=xint(i3,j)
       xarealist(3,j,narea)=xint(i4,j)
      enddo

      do i=1,sdim+1
      do j=1,sdim
        xint(i,j)=x(i,j)
      enddo
      enddo
      call shrink3D(xint,x,phi,i2,i4,sdim)
      call shrink3D(xint,x,phi,i2,i3,sdim)
      nlist=nlist+1
      if (nlist.gt.MAXTET) then
        print *,"nlist invalid"
        stop
      endif
      do i=1,sdim+1
      do j=1,sdim
        xtetlist(i,j,nlist)=xint(i,j)
      enddo
      enddo

      return
      end subroutine shrink_gableroof_list



   
subroutine volume_sanity_check()

use global_utility_module

IMPLICIT NONE

INTEGER_T bfact
REAL_T xsanity(3),dxgrid(3)
REAL_T xsten0(-3:3,3)
REAL_T hangle,hintercept,intercept
REAL_T angle(2)
REAL_T nslope(3)
INTEGER_T Nangle,shapeflag,sdim,Nangle2,nodedomain,coord
INTEGER_T i_int,a1,a2,i,j,k,inode,dir
REAL_T volslow,areaslow,volall
REAL_T cenall(3)
REAL_T censlow(3)
REAL_T areacenslow(3)
REAL_T xtet(4,3)
REAL_T xnode3d(8,3)
REAL_T xnode2d(4,2)
REAL_T phinode(8)
REAL_T xtarget(3)
INTEGER_T fullelementfast,linearcut,nhalf
REAL_T t1,t2

REAL_T cum_volume,cum_area
REAL_T cum_centroid(3)
REAL_T cum_areacentroid(3)

 nhalf=3
 bfact=2
 xsanity(1)=0.125
 xsanity(2)=0.25
 xsanity(3)=0.32
 dxgrid(1)=0.013
 dxgrid(2)=0.031
 dxgrid(3)=0.072
 do i=-nhalf,nhalf
  do dir=1,3
   xsten0(i,dir)=xsanity(dir)+i*dxgrid(dir)*half
  enddo
 enddo

 Nangle=512
 hangle=2*Pi/Nangle
 hintercept=two*dxgrid(3)/Nangle
 shapeflag=0

 call cpu_time(t1)

 print *,"sanity check N=",Nangle
 do sdim=2,3
  print *,"sanity check sdim= ",sdim
  if (sdim.eq.2) then
   Nangle2=0
  endif

  nodedomain=4*(sdim-1)
  do coord=0,3

   if ((coord.eq.0).or. &
       ((coord.eq.1).and.(sdim.eq.2)).or. &
       (coord.eq.3)) then
    print *,"sanity check coord=",coord
    do i_int=0,Nangle
     intercept=-dxgrid(3)+i_int*hintercept
     do a1=0,Nangle
      angle(1)=a1*hangle

      do a2=0,Nangle2
       angle(2)=a2*hangle
       if (sdim.eq.3) then
        call angle_to_slope3D(angle,nslope,sdim)
       else if (sdim.eq.2) then
        call angle_to_slope2D(angle,nslope,sdim)
       else
        print *,"sdim invalid"
        stop
       endif

       inode=1
       if (sdim.eq.3) then
        do k=-1,1,2
        do j=-1,1,2
        do i=-1,1,2
         do dir=1,sdim
          if (dir.eq.1) then
           xnode3d(inode,dir)=xsten0(i,dir)
          else if (dir.eq.2) then
           xnode3d(inode,dir)=xsten0(j,dir)
          else if (dir.eq.sdim) then
           xnode3d(inode,dir)=xsten0(k,dir)
          else
           print *,"dir invalid volume sanity check"
           stop
          endif
          xtarget(dir)=xnode3d(inode,dir)
         enddo  ! dir
         call distfunc(bfact,dxgrid,xsten0,nhalf, &
            intercept,nslope,xtarget, &
            phinode(inode),sdim)
         
         inode=inode+1
        enddo
        enddo
        enddo  ! i,j,k
       else if (sdim.eq.2) then
        do j=-1,1,2
        do i=-1,1,2
         do dir=1,sdim
          if (dir.eq.1) then
           xnode2d(inode,dir)=xsten0(i,dir)
          else if (dir.eq.2) then
           xnode2d(inode,dir)=xsten0(j,dir)
          else
           print *,"dir invalid volume sanity check 2"
           stop
          endif
          xtarget(dir)=xnode2d(inode,dir)
         enddo  ! dir
         call distfunc(bfact,dxgrid,xsten0,nhalf, &
            intercept,nslope,xtarget, &
            phinode(inode),sdim)
         
         inode=inode+1
        enddo
        enddo ! i,j
       else
        print *,"sdim invalid"
        stop
       endif

       if (inode.ne.nodedomain+1) then
        print *,"inode invalid"
        stop
       endif

       call cell_intersection_grid( &
         bfact,dxgrid,xsten0,nhalf, &
         phinode, &
         volslow,censlow,areaslow,areacenslow,coord, &
         volall,cenall,sdim)

       if (1.eq.1) then

        fullelementfast=1
          ! 1 if input is a plane, 0 if input might be plane, 
          ! -1 if force break up of domain.
        linearcut=1
        if (sdim.eq.2) then
         call intersection_volume( &
          cum_volume,cum_area,cum_centroid,cum_areacentroid, &
          phinode,xnode2d,nodedomain, &
          coord,sdim,fullelementfast,linearcut)
        else if (sdim.eq.3) then
         call intersection_volume( &
          cum_volume,cum_area,cum_centroid,cum_areacentroid, &
          phinode,xnode3d,nodedomain, &
          coord,sdim,fullelementfast,linearcut)
        else
         print *,"sdim invalid"
         stop
        endif
    
       else
        call cell_intersection_grid( &
         bfact,dxgrid,xsten0,nhalf, &
         phinode, &
         cum_volume,cum_centroid,cum_area,cum_areacentroid,coord, &
         volall,cenall,sdim)

       endif


       if (abs(cum_volume-volslow).gt.VOF_CUTOFF) then
        print *,"volume incorrect"
        print *,"cum_volume,volslow ",cum_volume,volslow
        print *,"nodedomain ",nodedomain
        do i=1,nodedomain
         print *,"i,phi ",i,phinode(i)
         do dir=1,sdim
          if (sdim.eq.2) then
           print *,"i,dir,xnode2d ",i,dir,xnode2d(i,dir)
          else
           print *,"i,dir,xnode3d ",i,dir,xnode3d(i,dir)
          endif
         enddo 
        enddo
        stop
       endif

       if (abs(cum_area-areaslow).gt.VOF_CUTOFF) then
        print *,"area incorrect"
        print *,"cum_area,areaslow ",cum_area,areaslow
        print *,"nodedomain ",nodedomain
        do i=1,nodedomain
         print *,"i,phi ",i,phinode(i)
         do dir=1,sdim
          if (sdim.eq.2) then
           print *,"i,dir,xnode2d ",i,dir,xnode2d(i,dir)
          else
           print *,"i,dir,xnode3d ",i,dir,xnode3d(i,dir)
          endif
         enddo 
        enddo
        stop
       endif

       do dir=1,sdim
        if (abs(cum_centroid(dir)-censlow(dir)).gt.VOF_CUTOFF) then
         print *,"centroid incorrect"
         stop
        endif
       enddo
       do dir=1,sdim
        if (abs(cum_areacentroid(dir)-areacenslow(dir)).gt.VOF_CUTOFF) then
         print *,"area centroid incorrect"
         stop
        endif
       enddo

      enddo ! a2
     enddo ! a1
    enddo ! i_int
   endif ! coord=0 or sdim=2
  enddo ! coord   
   
 enddo ! sdim

 call cpu_time(t2)
 print *,"elapsed time for sanity check: ",t2-t1

return
end subroutine volume_sanity_check


! find line segment that make up the intersection of the region where
! phi>0 and the given segment
      subroutine list_segment(phi,x,xsegmentlist,nseg,sdim)
      IMPLICIT NONE

      INTEGER_T nseg,sdim
      REAL_T phi(2)
      REAL_T x(2)
      REAL_T xsegmentlist(2)

      if ((sdim.ne.2).and.(sdim.ne.3)) then
       print *,"sdim invalid"
       stop
      endif

      if ((phi(1).ge.zero).and.(phi(2).ge.zero)) then
       nseg=1
       xsegmentlist(1)=x(1) 
       xsegmentlist(2)=x(2) 
      else if ((phi(1).lt.zero).and.(phi(2).lt.zero)) then
       nseg=0
      else if ((phi(1).ge.zero).and.(phi(2).lt.zero)) then
       nseg=1
       xsegmentlist(1)=x(1) 
       xsegmentlist(2)=x(1)+phi(1)*(x(2)-x(1))/(phi(1)-phi(2))
      else if ((phi(1).lt.zero).and.(phi(2).ge.zero)) then
       nseg=1
       xsegmentlist(2)=x(2) 
       xsegmentlist(1)=x(2)-phi(2)*(x(2)-x(1))/(phi(2)-phi(1))
      else
       print *,"phi invalid in list segment"
       stop
      endif

      return
      end subroutine list_segment

! OUTLINE OF YANG's routine:
! 1. initialize uncaptured_volume_fraction=1.0
! 2. for iplane=1..nmat and uncaptured_volume_fraction>0
! 3.  found=false
! 4.  for im=1..nmat, and found=false
! 5.   if (volumefraction(im)>=uncapture_volume_fraction-eps) then
! 6.    initialize uninitialized ids with "im"
! 7.    found=true 
! 8.    uncaptured_volume_fraction=0.0
! 9.   else if (order(im)=iplane) then
! 10.   cut list of tets with plane and init ids.
!         ("list_tets" cuts a tet with a plane and produces more tets)
!         (list_tets is called twice with + or - phi.)
!         (call "areaXYZ" to find areas of faces)
! 11.   uncaptured_volume_fraction=uncaptured_volume_fraction-
! 12.     volumefraction(im)
! 13.   found=true     
! 14.  else
! 15.   do nothing
! 16.  endif
! 17. enddo !im
! 18.enddo ! iplane
! FINITE VOLUME COEFFICIENTS:
! face_area(1..sdim,1..2,1..nmat+ncombine)
! face_area_internal(1..ncombine)
! e.g. if nmat=3 (11,22,33) then ncombine=12,13,23 (3)
! e.g. if nmat=4 then ncombine=12,13,14,23,24,34  (6)

! find tetrahedra that make up the intersection of a plane with
! a tetrahedron
! plane is specified by values of "phi" on the nodes of the tet.
! "x" are the coordinates of the tet.
      subroutine list_tets(phi,x,xtetlist,nlist,xarealist,narea,sdim)
      IMPLICIT NONE 
   
      INTEGER_T sdim 
      REAL_T phi(sdim+1)
      REAL_T x(sdim+1,sdim)
      REAL_T xtetlist(sdim+1,sdim,MAXTET)
      REAL_T xarealist(3,3,MAXAREA)
      INTEGER_T nlist,i,j,narea

      if (sdim.ne.3) then
       print *,"sdim invalid list_tets"
       stop
      endif

      nlist=0
      narea=0
      if ((phi(1).le.zero).and.(phi(2).le.zero).and. &
          (phi(3).le.zero).and.(phi(4).le.zero)) then
       nlist=0
      else if ((phi(1).ge.zero).and.(phi(2).ge.zero).and. &
               (phi(3).ge.zero).and.(phi(4).ge.zero)) then
       nlist=1
       do i=1,4
       do j=1,3
        xtetlist(i,j,1)=x(i,j)
       enddo 
       enddo 
      else if ((phi(1).lt.zero).and.(phi(2).ge.zero).and. &
               (phi(3).ge.zero).and.(phi(4).ge.zero)) then
       call shrink_list_tet(phi,x,1,2,3,4,1,xtetlist,nlist, &
        xarealist,narea)
      else if ((phi(2).lt.zero).and.(phi(1).ge.zero).and. &
               (phi(3).ge.zero).and.(phi(4).ge.zero)) then
       call shrink_list_tet(phi,x,2,1,3,4,1,xtetlist,nlist, &
        xarealist,narea)
      else if ((phi(3).lt.zero).and.(phi(1).ge.zero).and. &
               (phi(2).ge.zero).and.(phi(4).ge.zero)) then
       call shrink_list_tet(phi,x,3,1,2,4,1,xtetlist,nlist, &
        xarealist,narea)
      else if ((phi(4).lt.zero).and.(phi(1).ge.zero).and. &
               (phi(2).ge.zero).and.(phi(3).ge.zero)) then
       call shrink_list_tet(phi,x,4,1,2,3,1,xtetlist,nlist, &
        xarealist,narea)
      else if ((phi(1).ge.zero).and.(phi(2).lt.zero).and. &
               (phi(3).lt.zero).and.(phi(4).lt.zero)) then
       call shrink_list_tet(phi,x,1,2,3,4,0,xtetlist,nlist, &
        xarealist,narea)
      else if ((phi(2).ge.zero).and.(phi(1).lt.zero).and. &
               (phi(3).lt.zero).and.(phi(4).lt.zero)) then
       call shrink_list_tet(phi,x,2,1,3,4,0,xtetlist,nlist, &
        xarealist,narea)
      else if ((phi(3).ge.zero).and.(phi(1).lt.zero).and. &
               (phi(2).lt.zero).and.(phi(4).lt.zero)) then
       call shrink_list_tet(phi,x,3,1,2,4,0,xtetlist,nlist, &
        xarealist,narea)
      else if ((phi(4).ge.zero).and.(phi(1).lt.zero).and. &
               (phi(2).lt.zero).and.(phi(3).lt.zero)) then
       call shrink_list_tet(phi,x,4,1,2,3,0,xtetlist,nlist, &
        xarealist,narea)
      else if ((phi(1).lt.zero).and.(phi(2).lt.zero).and. &
               (phi(3).ge.zero).and.(phi(4).ge.zero)) then
       call shrink_gableroof_list(phi,x,3,4,1,2,xtetlist,nlist, &
        xarealist,narea)
      else if ((phi(1).lt.zero).and.(phi(3).lt.zero).and. &
               (phi(2).ge.zero).and.(phi(4).ge.zero)) then
       call shrink_gableroof_list(phi,x,2,4,1,3,xtetlist,nlist, &
        xarealist,narea)
      else if ((phi(1).ge.zero).and.(phi(2).ge.zero).and. &
               (phi(3).lt.zero).and.(phi(4).lt.zero)) then
       call shrink_gableroof_list(phi,x,1,2,3,4,xtetlist,nlist, &
        xarealist,narea)
      else if ((phi(1).ge.zero).and.(phi(3).ge.zero).and. &
               (phi(2).lt.zero).and.(phi(4).lt.zero)) then
       call shrink_gableroof_list(phi,x,1,3,2,4,xtetlist,nlist, &
        xarealist,narea)
      else if ((phi(3).ge.zero).and.(phi(2).ge.zero).and. &
               (phi(1).lt.zero).and.(phi(4).lt.zero)) then
       call shrink_gableroof_list(phi,x,3,2,1,4,xtetlist,nlist, &
        xarealist,narea)
      else if ((phi(1).ge.zero).and.(phi(4).ge.zero).and. &
               (phi(2).lt.zero).and.(phi(3).lt.zero)) then
       call shrink_gableroof_list(phi,x,1,4,2,3,xtetlist,nlist, &
        xarealist,narea)
      else
       print *,"bust list_tets"
       print *,"phi : ",phi(1),phi(2),phi(3),phi(4)
       stop
      endif

      return
      end subroutine list_tets

! find the triangles that make up the intersection of two triangles.
! xtrilist_old is a scratch variable
      subroutine intersect_tri(x1,x2,xtrilist_old,xtrilist,nlist,nmax,sdim)
      IMPLICIT NONE

      INTEGER_T sdim
      INTEGER_T nlist,nmax,nlist_old,i,j,iplane,n,nsub,narea,n2
      REAL_T x1(sdim+1,sdim)
      REAL_T x1old(sdim+1,sdim)
      REAL_T x2(sdim+1,sdim)
      REAL_T xtrilist_old(sdim+1,sdim,nmax)
      REAL_T xtrilist(sdim+1,sdim,nmax)
      REAL_T xsublist(sdim+1,sdim,MAXTET)
      REAL_T xarealist(2,2,MAXAREA)
      REAL_T phi1(sdim+1)
      REAL_T coeff(sdim)
      REAL_T voltest
      REAL_T centroidtest(sdim)
      REAL_T xcandidate(sdim+1,sdim)

      INTEGER_T itan(sdim)
      REAL_T sign,mag,maxside,testside1,testside2
      INTEGER_T local_coord
      
      if (sdim.ne.2) then
       print *,"sdim bust intersect tri"
       stop
      endif
      local_coord=0

      maxside=zero
      do i=2,sdim+1
       testside1=zero
       testside2=zero
       do j=1,sdim
        testside1=testside1+(x1(1,j)-x1(i,j))**2
        testside2=testside2+(x2(1,j)-x2(i,j))**2
       enddo
       testside1=sqrt(testside1)
       testside2=sqrt(testside2)
       if (testside1.gt.maxside) then
        maxside=testside1
       endif
       if (testside2.gt.maxside) then
        maxside=testside2
       endif
      enddo ! i

      if (maxside.gt.one) then
       maxside=one
      endif
      if (maxside.le.zero) then
       print *,"maxside invalid"
       stop
      endif


      nlist_old=1
      do i=1,sdim+1
      do j=1,sdim
       xtrilist_old(i,j,1)=x1(i,j)
      enddo
      enddo

! intersect the members of xtrilist_old with the planes of x2
! phi = s(n dot (x-x0)) where x0 is a point on the plane.
! s=n dot (xp-x0) where xp is not on the plane.

      do iplane=1,sdim+1

       nlist=0
       if (iplane.eq.1) then
        itan(1)=2
        itan(2)=3
       else if (iplane.eq.2) then
        itan(1)=1
        itan(2)=3
       else if (iplane.eq.3) then
        itan(1)=1
        itan(2)=2
       else
        print *,"iplane invalid"
        stop
       endif
 
       coeff(1)=-(x2(itan(2),2)-x2(itan(1),2))
       coeff(2)=(x2(itan(2),1)-x2(itan(1),1))

       mag=sqrt(coeff(1)**2+coeff(2)**2)

       if (mag.gt.MLSVOFTOL*maxside) then
        do j=1,sdim
         coeff(j)=coeff(j)/mag
        enddo
        sign=zero
        do j=1,sdim
         sign=sign+coeff(j)*(x2(iplane,j)-x2(itan(1),j))
        enddo
       
        if (abs(sign).gt.MLSVOFTOL*maxside) then

         do n=1,nlist_old

          do i=1,sdim+1
          do j=1,sdim
           x1old(i,j)=xtrilist_old(i,j,n)
          enddo
          enddo

          do i=1,sdim+1
           phi1(i)=zero
           do j=1,sdim
            phi1(i)=phi1(i)+coeff(j)*(x1old(i,j)-x2(itan(1),j))
           enddo
           phi1(i)=phi1(i)*sign
          enddo
          call list_tris(phi1,x1old,xsublist,nsub,xarealist,narea,sdim)
          do n2=1,nsub

           do i=1,sdim+1
           do j=1,sdim
            xcandidate(i,j)=xsublist(i,j,n2)
           enddo
           enddo
           call tetrahedron_volume(xcandidate,voltest, &
            centroidtest,local_coord,sdim)

           if (voltest.gt.zero) then
            nlist=nlist+1
            if (nlist.gt.nmax) then
             print *,"nlist overflow in intersect tri"
             print *,"nlist,nmax ",nlist,nmax
             stop
            endif
            do i=1,sdim+1
            do j=1,sdim
             xtrilist(i,j,nlist)=xcandidate(i,j)
            enddo
            enddo
           endif

          enddo  ! n2
         enddo  ! n
         nlist_old=nlist
         do n=1,nlist
          do i=1,sdim+1
          do j=1,sdim
           xtrilist_old(i,j,n)=xtrilist(i,j,n)
          enddo
          enddo
         enddo
        else
         print *,"sign=0 in intersect_tri"
         stop
        endif 
       else
        print *,"mag=0 in intersect_tri"
        stop
       endif
      enddo  ! iplane

      return
      end subroutine intersect_tri

! find the tetrahedra that make up the intersection of two 
! tetrahedrons.  xtetlist_old is a scratch variable.
      subroutine intersect_tet(x1,x2,xtetlist_old,xtetlist,nlist,nmax,sdim)
      IMPLICIT NONE

      INTEGER_T sdim
      INTEGER_T nlist,nmax,nlist_old,i,j,iplane,n,nsub,narea,n2
      INTEGER_T ii,jj
      REAL_T x1(sdim+1,sdim)
      REAL_T x1old(sdim+1,sdim)
      REAL_T x2(sdim+1,sdim)
      REAL_T xtetlist_old(sdim+1,sdim,nmax)
      REAL_T xtetlist(sdim+1,sdim,nmax)
      REAL_T xsublist(sdim+1,sdim,MAXTET)
      REAL_T xarealist(3,3,MAXAREA)
      REAL_T phi1(sdim+1)
      REAL_T coeff(sdim)
      REAL_T voltest
      REAL_T centroidtest(sdim)
      REAL_T xcandidate(sdim+1,sdim)

      INTEGER_T itan(sdim)
      REAL_T sign,mag
      REAL_T vec(2,sdim)
      REAL_T maxside,testside1,testside2
      INTEGER_T local_coord

      if (sdim.ne.3) then
       print *,"sdim bust intersect tet"
       stop
      endif

      maxside=zero
      do i=2,sdim+1
       testside1=zero
       testside2=zero
       do j=1,sdim
        testside1=testside1+(x1(1,j)-x1(i,j))**2
        testside2=testside2+(x2(1,j)-x2(i,j))**2
       enddo
       testside1=sqrt(testside1)
       testside2=sqrt(testside2)
       if (testside1.gt.maxside) then
        maxside=testside1
       endif
       if (testside2.gt.maxside) then
        maxside=testside2
       endif
      enddo ! i

      if (maxside.gt.one) then
       maxside=one
      endif
      if (maxside.le.zero) then
       print *,"maxside invalid"
       stop
      endif


      nlist_old=1
      do i=1,sdim+1
      do j=1,sdim
       xtetlist_old(i,j,1)=x1(i,j)
      enddo
      enddo

! intersect the members of xtetlist_old with the planes of x2
! phi = s(n dot (x-x0)) where x0 is a point on the plane.
! s=n dot (xp-x0) where xp is not on the plane.

      do iplane=1,sdim+1

       nlist=0
       if (iplane.eq.1) then
        itan(1)=2
        itan(2)=3
        itan(3)=4
       else if (iplane.eq.2) then
        itan(1)=1
        itan(2)=3
        itan(3)=4
       else if (iplane.eq.3) then
        itan(1)=1
        itan(2)=2
        itan(3)=4
       else if (iplane.eq.4) then
        itan(1)=1
        itan(2)=2
        itan(3)=3
       else
        print *,"iplane invalid"
        stop
       endif

       do i=1,2
        if ((iplane.eq.itan(i)).or.(iplane.eq.itan(3))) then
         print *,"bust intersect_tet"
         print *,"iplane=",iplane
         print *,"i=",i
         print *,"itan(i)=",itan(i)
         print *,"itan(3)=",itan(3)
         print *,"sdim=",sdim
         stop
        endif
        do j=1,sdim
         vec(i,j)=x2(itan(i),j)-x2(itan(3),j)
        enddo
       enddo
       coeff(1)=vec(1,2)*vec(2,3)-vec(1,3)*vec(2,2)
       coeff(2)=vec(1,3)*vec(2,1)-vec(1,1)*vec(2,3)
       coeff(3)=vec(1,1)*vec(2,2)-vec(2,1)*vec(1,2)
       mag=sqrt(coeff(1)**2+coeff(2)**2+coeff(3)**2)

       if (mag.gt.MLSVOFTOL*maxside) then
        do j=1,sdim
         coeff(j)=coeff(j)/mag
        enddo
        sign=zero
        do j=1,sdim
         sign=sign+coeff(j)*(x2(iplane,j)-x2(itan(1),j))
        enddo
       
        if (abs(sign).gt.MLSVOFTOL*maxside) then

         do n=1,nlist_old

          do i=1,sdim+1
          do j=1,sdim
           x1old(i,j)=xtetlist_old(i,j,n)
          enddo
          enddo

          do i=1,sdim+1
           phi1(i)=zero
           do j=1,sdim
            phi1(i)=phi1(i)+coeff(j)*(x1old(i,j)-x2(itan(1),j))
           enddo
           phi1(i)=phi1(i)*sign
          enddo
          call list_tets(phi1,x1old,xsublist,nsub,xarealist,narea,sdim)
          do n2=1,nsub

           do i=1,sdim+1
           do j=1,sdim
            xcandidate(i,j)=xsublist(i,j,n2)
           enddo
           enddo
           local_coord=0
           call tetrahedron_volume(xcandidate,voltest, &
            centroidtest,local_coord,sdim)

           if (voltest.gt.zero) then
            nlist=nlist+1
            if (nlist.gt.nmax) then
             print *,"nlist overflow in intersect_tet"
             print *,"nlist,nmax ",nlist,nmax
             stop
            endif
            do i=1,sdim+1
            do j=1,sdim
             xtetlist(i,j,nlist)=xcandidate(i,j)
            enddo
            enddo
           endif

          enddo  ! n2
         enddo  ! n
         nlist_old=nlist
         do n=1,nlist
          do i=1,sdim+1
          do j=1,sdim
           xtetlist_old(i,j,n)=xtetlist(i,j,n)
          enddo
          enddo
         enddo
        else
         print *,"sign=0 in intersect_tet"
         stop
        endif 
       else
        print *,"mag=0 in intersect_tet"
        stop
       endif
      enddo  ! iplane

      return
      end subroutine intersect_tet

! find the tets of intersection of a cube with a tet.
! xtetlist_old is a scratch variable.
      subroutine intersect_cube( &
        x1, &
        xsten,nhalf, &
        xtetlist_old,xtetlist,nlist,nmax,sdim)
      IMPLICIT NONE

      INTEGER_T sdim,nhalf
      INTEGER_T nlist,nmax,nlist_old,i,j,iplane,n,nsub,narea,n2
      REAL_T x1(sdim+1,sdim)
      REAL_T x1old(sdim+1,sdim)
      REAL_T xsten(-nhalf:nhalf,sdim)
      REAL_T x0(sdim)
      REAL_T xtetlist_old(sdim+1,sdim,nmax)
      REAL_T xtetlist(sdim+1,sdim,nmax)
      REAL_T xsublist(sdim+1,sdim,MAXTET)
      REAL_T xarealist(sdim,sdim,MAXAREA)
      REAL_T phi1(sdim+1)
      REAL_T coeff(sdim)
      REAL_T voltest
      REAL_T centroidtest(sdim)
      REAL_T xcandidate(sdim+1,sdim)
      INTEGER_T local_coord

      if (nhalf.lt.1) then
       print *,"nhalf invalid intersect cube"
       stop
      endif
      if ((sdim.ne.2).and.(sdim.ne.3)) then
       print *,"sdim bust intersect cube"
       stop
      endif

      nlist_old=1
      do i=1,sdim+1
      do j=1,sdim
       xtetlist_old(i,j,1)=x1(i,j)
      enddo
      enddo

! intersect the members of xtetlist_old with the planes of x2
! phi = (n dot (x-x0)) where x0 is a point on the plane.

      do iplane=1,2*sdim

       nlist=0

       do i=1,sdim
        coeff(i)=zero
        x0(i)=xsten(-1,i)
       enddo
 
       if (iplane.eq.1) then  ! left side
        coeff(1)=one
        x0(1)=xsten(-1,1)
       else if (iplane.eq.2) then  ! right side
        coeff(1)=-one
        x0(1)=xsten(1,1)
       else if (iplane.eq.3) then  ! front
        coeff(2)=one
        x0(2)=xsten(-1,2)
       else if (iplane.eq.4) then  ! back
        coeff(2)=-one
        x0(2)=xsten(1,2)
       else if (iplane.eq.5) then  ! bottom
        if (sdim.ne.3) then
         print *,"dimension bust"
         stop
        endif
        coeff(sdim)=one
        x0(sdim)=xsten(-1,sdim)
       else if (iplane.eq.6) then  ! top
        if (sdim.ne.3) then
         print *,"dimension bust"
         stop
        endif
        coeff(sdim)=-one
        x0(sdim)=xsten(1,sdim)
       else
        print *,"iplane invalid"
        stop
       endif

       do n=1,nlist_old

        do i=1,sdim+1
        do j=1,sdim
         x1old(i,j)=xtetlist_old(i,j,n)
        enddo
        enddo

        do i=1,sdim+1
         phi1(i)=zero
         do j=1,sdim
          phi1(i)=phi1(i)+coeff(j)*(x1old(i,j)-x0(j))
         enddo
        enddo
     
         ! in: intersect_cube 
        if (sdim.eq.3) then
         call list_tets(phi1,x1old,xsublist,nsub,xarealist,narea,sdim)
        else if (sdim.eq.2) then
         call list_tris(phi1,x1old,xsublist,nsub,xarealist,narea,sdim)
        else
         print *,"sdim invalid"
         stop
        endif


        do n2=1,nsub

         do i=1,sdim+1
         do j=1,sdim
          xcandidate(i,j)=xsublist(i,j,n2)
         enddo
         enddo
         local_coord=0
         call tetrahedron_volume(xcandidate,voltest, &
          centroidtest,local_coord,sdim)

         if (voltest.gt.zero) then
          nlist=nlist+1
          if (nlist.gt.nmax) then
           print *,"nlist overflow in intersect_cube"
           print *,"nlist,nmax ",nlist,nmax
           stop
          endif
          do i=1,sdim+1
          do j=1,sdim
           xtetlist(i,j,nlist)=xcandidate(i,j)
          enddo
          enddo
         endif

        enddo  ! n2
       enddo  ! n
       nlist_old=nlist
       do n=1,nlist
        do i=1,sdim+1
        do j=1,sdim
         xtetlist_old(i,j,n)=xtetlist(i,j,n)
        enddo
        enddo
       enddo
      enddo  ! iplane

      return
      end subroutine intersect_cube





! get volume/centroid of a triangulated region
      subroutine get_cut_geom3D(xtetlist,nlist,nmax,levelrz,volcut, &
        cencut,sdim)
      IMPLICIT NONE

      INTEGER_T nlist,nmax,levelrz,dir,sdim,i,j,n
      REAL_T volcut
      REAL_T cencut(sdim)
      REAL_T xtetlist(sdim+1,sdim,nmax)
      REAL_T volumelist
      REAL_T centroidlist(sdim)
      REAL_T xint(sdim+1,sdim)

      if ((sdim.ne.3).and.(sdim.ne.2)) then
       print *,"sdim invalid get_cut_geom3D"
       stop
      endif 

      volcut=zero
      do j=1,sdim
       cencut(j)=zero
      enddo
      do n=1,nlist
       do i=1,sdim+1
       do j=1,sdim
        xint(i,j)=xtetlist(i,j,n)
       enddo
       enddo
     
       call tetrahedron_volume(xint,volumelist,centroidlist,levelrz,sdim)

       volcut=volcut+volumelist
       do j=1,sdim
        cencut(j)=cencut(j)+centroidlist(j)*volumelist
       enddo
      enddo
      if (volcut.gt.zero) then
       do j=1,sdim
        cencut(j)=cencut(j)/volcut
       enddo
      else
       do j=1,sdim
        cencut(j)=zero
       enddo
      endif

      return
      end subroutine get_cut_geom3D


! normal points from light to dark   phi=n dot (x-x0plane) + intercept
! vof, ref centroid, order,slope,intercept  x nmat
! tetrahedra representing unfilled region
! (where LS < 0)
! routine starts with tetrahedralization of a cell,
! then intersects this region with the plane (LS<0 side)
! of already initialized materials.
      subroutine tets_box_planes( &
       bfact,dx,xsten0,nhalf0, &
       xsten_box,nhalf_box,mofdata, &
       xtetlist,nlist,nmax,nmat,sdim)
      IMPLICIT NONE

      INTEGER_T nmat,sdim,bfact,nhalf0,nhalf_box
      INTEGER_T symmetry_flag,ntetbox
      INTEGER_T nlist,nmax,nlist_old,i,j,k,iplane,n,nsub,narea,n2
      INTEGER_T local_coord,id,icrit,im,vofcomp,iorder
      INTEGER_T dir
      REAL_T dx(sdim)
      REAL_T xsten0(-nhalf0:nhalf0,sdim)
      REAL_T xsten_box(-nhalf_box:nhalf_box,sdim)
      REAL_T mofdata(nmat*(2*sdim+3))
      REAL_T xtetlist_old(sdim+1,sdim,nmax)
      REAL_T xtetlist(sdim+1,sdim,nmax)
      REAL_T xsublist(sdim+1,sdim,MAXTET)
      REAL_T xarealist(sdim,sdim,MAXAREA)
      REAL_T phi1(sdim+1),dummyphi(sdim+1)
      REAL_T coeff(sdim)
      REAL_T voltest
      REAL_T centroidtest(sdim)
      REAL_T xcandidate(sdim+1,sdim)
      REAL_T xx(sdim+1,sdim)
      REAL_T x1(sdim+1,sdim)
      REAL_T x1old(sdim+1,sdim)
      REAL_T x2(sdim+1,sdim)
      REAL_T xnode(4*(sdim-1),sdim)
      REAL_T phinode(4*(sdim-1))
      INTEGER_T inode
      REAL_T nn(sdim)
      REAL_T intercept

      if ((nhalf0.lt.1).or.(nhalf_box.lt.1)) then
       print *,"nhalf invalid tets box planes"
       stop
      endif
      if ((sdim.ne.3).and.(sdim.ne.2)) then
       print *,"sdim invalid tets_box_planes"
       stop
      endif
      if ((nmat.lt.1).or.(nmat.gt.MAX_NUM_MATERIALS)) then
       print *,"nmat invalid tets box planes"
       stop
      endif
      if (bfact.lt.1) then
       print *,"bfact invalid"
       stop
      endif
 
      local_coord=0
! YANG: get number of tetrahedra here.
      inode=1

      symmetry_flag=0
      call get_ntetbox(ntetbox,symmetry_flag,sdim)
      nlist_old=ntetbox
      nlist=ntetbox

! YANG: CONVERT xsten_box to xnode.
      if (sdim.eq.3) then
       do k=-1,1,2
       do j=-1,1,2
       do i=-1,1,2
        xnode(inode,1)=xsten_box(i,1)
        xnode(inode,2)=xsten_box(j,2)
        xnode(inode,sdim)=xsten_box(k,sdim)
        phinode(inode)=one
        inode=inode+1
       enddo
       enddo
       enddo
      else if (sdim.eq.2) then
       do j=-1,1,2
       do i=-1,1,2
        xnode(inode,1)=xsten_box(i,1)
        xnode(inode,2)=xsten_box(j,2)
        phinode(inode)=one
        inode=inode+1
       enddo
       enddo
      else
       print *,"sdim invalid"
       stop
      endif

! YANG: add these tetrahedra to your data structure.
      do id=1,ntetbox
       call extract_tet(xnode,phinode,xx,dummyphi,id,symmetry_flag,sdim)
       do i=1,sdim+1
       do j=1,sdim
        xtetlist_old(i,j,id)=xx(i,j)
        xtetlist(i,j,id)=xx(i,j)
       enddo 
       enddo 
      enddo

     
! EXAMPLE OF TRAVERSING THE PLANES IN THE PROPER ORDER.
      do iplane=1,nmat
       icrit=0
       do im=1,nmat
        vofcomp=(im-1)*(2*sdim+3)+1
        iorder=NINT(mofdata(vofcomp+sdim+1))
        if (iorder.eq.iplane) then
         icrit=im
        endif
       enddo
       if (icrit.gt.0) then
        nlist=0
        vofcomp=(icrit-1)*(2*sdim+3)+1
         ! want intersection where phi_original<0 but routine
         ! checks if phi_input>0.  So negate phi_original.
        do dir=1,sdim
         nn(dir)=-mofdata(vofcomp+sdim+1+dir)  
        enddo
        intercept=-mofdata(vofcomp+2*sdim+2)
 
        do n=1,nlist_old

         do i=1,sdim+1
         do j=1,sdim
          x1old(i,j)=xtetlist_old(i,j,n)
         enddo
         enddo

           ! LEVELSET FOR CUTTING PLANE
         do i=1,sdim+1
          phi1(i)=intercept
          do j=1,sdim
           phi1(i)=phi1(i)+nn(j)*(x1old(i,j)-xsten0(0,j))
          enddo
         enddo
          ! tetrahedras representing intersection of region where phi1>0
          ! (phi_original<0) and the tetrahedra x1old
         if (sdim.eq.3) then
          call list_tets(phi1,x1old,xsublist,nsub,xarealist,narea,sdim)
         else if (sdim.eq.2) then
          call list_tris(phi1,x1old,xsublist,nsub,xarealist,narea,sdim)
         else
          print *,"sdim invalid"
          stop
         endif

         do n2=1,nsub

          do i=1,sdim+1
          do j=1,sdim
           xcandidate(i,j)=xsublist(i,j,n2)
          enddo
          enddo
          call tetrahedron_volume(xcandidate,voltest, &
           centroidtest,local_coord,sdim)

          if (voltest.gt.zero) then
           nlist=nlist+1
           if (nlist.gt.nmax) then
            print *,"nlist overflow in tets_box_planes"
            print *,"nlist,nmax ",nlist,nmax
            stop
           endif
           do i=1,sdim+1
           do j=1,sdim
            xtetlist(i,j,nlist)=xcandidate(i,j)
           enddo
           enddo
          endif

         enddo  ! n2
        enddo  ! n
        nlist_old=nlist
        do n=1,nlist
         do i=1,sdim+1
         do j=1,sdim
          xtetlist_old(i,j,n)=xtetlist(i,j,n)
         enddo
         enddo
        enddo
       else if (icrit.eq.0) then
        nlist=nlist_old
        do n=1,nlist
         do i=1,sdim+1
         do j=1,sdim
          xtetlist(i,j,n)=xtetlist_old(i,j,n)
         enddo
         enddo
        enddo
       else
        print *,"icrit invalid"
        stop
       endif 
      enddo ! iplane

      return
      end subroutine tets_box_planes

      subroutine tets_box_planes_super( &
       bfact,dx,xsten0,nhalf0, &
       mofdata, &
       xtetlist,nlist,nmax,nmat,use_super_cell,sdim)
      IMPLICIT NONE

      INTEGER_T nmat,sdim,use_super_cell,bfact,nhalf0
      INTEGER_T isten,nhalf2
      INTEGER_T nlist,nmax,nlist_local,i1,j1,k1,itri
      INTEGER_T nn
      INTEGER_T dir
      REAL_T dx(sdim)
      REAL_T xsten0(-nhalf0:nhalf0,sdim)
      REAL_T xsten2(-1:1,sdim)
      REAL_T mofdata(nmat*(2*sdim+3))
      REAL_T xtetlist_local(sdim+1,sdim,nmax)
      REAL_T xtetlist(sdim+1,sdim,nmax)
      INTEGER_T ksten_low,ksten_high

      nhalf2=1

      if ((sdim.ne.2).and.(sdim.ne.3)) then
       print *,"sdim invalid tets_box_planes_super"
       stop
      endif
      if ((nmat.lt.1).or.(nmat.gt.MAX_NUM_MATERIALS)) then
       print *,"nmat invalid tets box planes super"
       stop
      endif
      if (bfact.lt.1) then
       print *,"bfact invalid"
       stop
      endif
      if (nhalf0.lt.3) then
       print *,"nhalf0 invalid"
       stop
      endif

       ! in: tets_box_planes_super
      if (use_super_cell.eq.0) then
       call tets_box_planes(bfact,dx,xsten0,nhalf0, &
        xsten0,nhalf0, &
        mofdata, &
        xtetlist,nlist,nmax,nmat,sdim)
      else if (use_super_cell.eq.1) then

       if (sdim.eq.3) then
        ksten_low=-1
        ksten_high=1
       else if (sdim.eq.2) then
        ksten_low=0
        ksten_high=0
       else
        print *,"sdim invalid"
        stop
       endif

       nlist=0

       do i1=-1,1
       do j1=-1,1
       do k1=ksten_low,ksten_high

        do isten=-1,1
         xsten2(isten,1)=xsten0(isten+2*i1,1)
         xsten2(isten,2)=xsten0(isten+2*j1,2)
         if (sdim.eq.3) then
          xsten2(isten,sdim)=xsten0(isten+2*k1,sdim)
         endif
        enddo ! isten

         ! in: tets_box_planes_super
        call tets_box_planes(bfact,dx,xsten0,nhalf0, &
         xsten2,nhalf2, &
         mofdata, &
         xtetlist_local,nlist_local,nmax,nmat,sdim)
        if (nlist_local+nlist.gt.nmax) then
         print *,"too many tetrahedrons in tets_box_planes_super"
         stop
        endif
        do nn=1,nlist_local
         nlist=nlist+1
         do itri=1,sdim+1
         do dir=1,sdim
          xtetlist(itri,dir,nlist)=xtetlist_local(itri,dir,nn)
         enddo
         enddo
        enddo  ! nn
          
       enddo
       enddo
       enddo ! i1,j1,k1

      else
       print *,"use_super_cell invalid"
       stop
      endif

      return
      end subroutine tets_box_planes_super


      subroutine tets_tet_planes( &
        bfact,dx,xsten0,nhalf0, &
        xtet,mofdata, &
        xtetlist,nlist,nmax,nmat,sdim)
      IMPLICIT NONE

      INTEGER_T sdim,nmat,bfact,nhalf0
      INTEGER_T nlist,nmax,nlist_old,i,j,iplane,n,nsub,narea,n2
      INTEGER_T local_coord,icrit,im,vofcomp,iorder
      INTEGER_T dir
      REAL_T xtet(sdim+1,sdim)
      REAL_T dx(sdim)
      REAL_T xsten0(-nhalf0:nhalf0,sdim)
      REAL_T mofdata(nmat*(2*sdim+3))
      REAL_T xtetlist_old(sdim+1,sdim,nmax)
      REAL_T xtetlist(sdim+1,sdim,nmax)
      REAL_T xsublist(sdim+1,sdim,MAXTET)
      REAL_T xarealist(sdim,sdim,MAXAREA)
      REAL_T phi1(sdim+1),dummyphi(sdim+1)
      REAL_T coeff(sdim)
      REAL_T voltest
      REAL_T centroidtest(sdim)
      REAL_T xcandidate(sdim+1,sdim)
      REAL_T x1(sdim+1,sdim)
      REAL_T x1old(sdim+1,sdim)
      REAL_T x2(sdim+1,sdim)
      REAL_T nn(sdim)
      REAL_T intercept

      if (nhalf0.lt.1) then
       print *,"nhalf0 invalid"
       stop
      endif
      if (bfact.lt.1) then
       print *,"bfact invalid"
       stop
      endif
      if ((sdim.ne.2).and.(sdim.ne.3)) then
       print *,"sdim invalid tets_tet_planes"
       stop
      endif
      if ((nmat.lt.1).or.(nmat.gt.MAX_NUM_MATERIALS)) then
       print *,"nmat invalid tets tet planes"
       stop
      endif
      local_coord=0

      nlist_old=1
      nlist=1
      do i=1,sdim+1
      do j=1,sdim
        xtetlist_old(i,j,1)=xtet(i,j)
        xtetlist(i,j,1)=xtet(i,j)
      enddo 
      enddo 

      do iplane=1,nmat
       icrit=0
       do im=1,nmat
        vofcomp=(im-1)*(2*sdim+3)+1
        iorder=NINT(mofdata(vofcomp+sdim+1))
        if (iorder.eq.iplane) then
         icrit=im
        endif
       enddo
       if (icrit.gt.0) then
        nlist=0
        vofcomp=(icrit-1)*(2*sdim+3)+1
         ! want intersection where phi_original<0 but routine
         ! checks if phi_input>0.  So negate phi_original.
        do dir=1,sdim
         nn(dir)=-mofdata(vofcomp+sdim+1+dir)  
        enddo
        intercept=-mofdata(vofcomp+2*sdim+2)
 
        do n=1,nlist_old

         do i=1,sdim+1
         do j=1,sdim
          x1old(i,j)=xtetlist_old(i,j,n)
         enddo
         enddo

         do i=1,sdim+1
          phi1(i)=intercept
          do j=1,sdim
           phi1(i)=phi1(i)+nn(j)*(x1old(i,j)-xsten0(0,j))
          enddo
         enddo
          ! triangles representing intersection of region where phi1>0
          ! (phi_original<0) and the triangle x1old
         if (sdim.eq.3) then
          call list_tets(phi1,x1old,xsublist,nsub,xarealist,narea,sdim)
         else if (sdim.eq.2) then
          call list_tris(phi1,x1old,xsublist,nsub,xarealist,narea,sdim)
         else
          print *,"sdim invalid"
          stop
         endif


         do n2=1,nsub

          do i=1,sdim+1
          do j=1,sdim
           xcandidate(i,j)=xsublist(i,j,n2)
          enddo
          enddo

          call tetrahedron_volume(xcandidate,voltest, &
           centroidtest,local_coord,sdim)

          if (voltest.gt.zero) then
           nlist=nlist+1
           if (nlist.gt.nmax) then
            print *,"nlist overflow in tets_tet_planes"
            print *,"nlist,nmax ",nlist,nmax
            stop
           endif
           do i=1,sdim+1
           do j=1,sdim
            xtetlist(i,j,nlist)=xcandidate(i,j)
           enddo
           enddo
          endif

         enddo  ! n2
        enddo  ! n
        nlist_old=nlist
        do n=1,nlist
         do i=1,sdim+1
         do j=1,sdim
          xtetlist_old(i,j,n)=xtetlist(i,j,n)
         enddo
         enddo
        enddo
       else if (icrit.eq.0) then
        nlist=nlist_old
        do n=1,nlist
         do i=1,sdim+1
         do j=1,sdim
          xtetlist(i,j,n)=xtetlist_old(i,j,n)
         enddo
         enddo
        enddo
       else
        print *,"icrit invalid"
        stop
       endif 
      enddo ! iplane

      return
      end subroutine tets_tet_planes
 

! find the volume/area/etc. of intersection of a plane with
! a tetrahedron
      subroutine intersection_volumeXYZ(phi,x,volumedark, &
       centroiddark,area,areacentroid,coord,sdim)
      IMPLICIT NONE

      INTEGER_T sdim,coord
      REAL_T xtetlist(sdim+1,sdim,MAXTET)
      REAL_T xarealist(sdim,sdim,MAXAREA)
      REAL_T phi(sdim+1)
      REAL_T x(sdim+1,sdim)
      REAL_T xint(sdim+1,sdim)
      REAL_T areacentroid(sdim)
      REAL_T areacentroidlist(sdim)
      REAL_T centroiddark(sdim)
      REAL_T centroidlistdark(sdim)

      REAL_T volumedark,area
      REAL_T volumelistdark,arealist

      INTEGER_T i,j,n,nlist,narea

      if (sdim.ne.3) then
       print *,"sdim invalid"
       stop
      endif

      call list_tets(phi,x,xtetlist,nlist,xarealist,narea,sdim)

      volumedark=zero
      area=zero
      do j=1,sdim
       centroiddark(j)=zero
       areacentroid(j)=zero
      enddo

      do n=1,nlist
       do i=1,sdim+1
       do j=1,sdim
        xint(i,j)=xtetlist(i,j,n)
       enddo
       enddo
       call tetrahedron_volume(xint,volumelistdark,centroidlistdark,coord,sdim)
       volumedark=volumedark+volumelistdark
       do j=1,sdim
        centroiddark(j)=centroiddark(j)+centroidlistdark(j)*volumelistdark
       enddo
      enddo
      if (volumedark.gt.zero) then
       do j=1,sdim
        centroiddark(j)=centroiddark(j)/volumedark
       enddo
      else
       do j=1,sdim
        centroiddark(j)=zero
       enddo
      endif

      do n=1,narea
       do i=1,sdim
       do j=1,sdim
        xint(i,j)=xarealist(i,j,n)
       enddo
       enddo
       do j=1,sdim
        xint(sdim+1,j)=zero
       enddo
       call areaXYZ(xint,1,2,3,arealist,areacentroidlist,coord)
       area=area+arealist
       do j=1,sdim
        areacentroid(j)=areacentroid(j)+areacentroidlist(j)*arealist
       enddo
      enddo

      if (area.gt.zero) then
       do j=1,sdim
        areacentroid(j)=areacentroid(j)/area
       enddo
      else
       do j=1,sdim
        areacentroid(j)=zero
       enddo
      endif

      return
      end subroutine intersection_volumeXYZ


      subroutine CISBOX(xsten,nhalf, &
       xlo,dx,i,j,k,lo,hi,bfact,level,vol,cen,rzflag,sdim)
      use global_utility_module

      IMPLICIT NONE

      INTEGER_T sdim
      INTEGER_T lo(sdim),hi(sdim)
      INTEGER_T bfact,rzflag,nhalf,level
      REAL_T xsten(-nhalf:nhalf,sdim)
      REAL_T dx(sdim)
      REAL_T xlo(sdim)
      INTEGER_T i,j,k
      REAL_T vol
      REAL_T cen(sdim)
      INTEGER_T dir

      if (nhalf.lt.1) then
       print *,"nhalf invalid cisbox"
       stop
      endif
      if (bfact.lt.1) then
       print *,"bfact invalid"
       stop
      endif
      if ((sdim.ne.2).and.(sdim.ne.3)) then
       print *,"sdim invalid"
       stop
      endif
      call gridsten_level(xsten,i,j,k,level,nhalf)
      call Box_volumeFAST(bfact,dx,xsten,nhalf,vol,cen,rzflag,sdim)

      return
      end subroutine CISBOX 


      subroutine CISBOXHALF(xsten,nhalf, &
       xlo,dx,i,j,k,iside,veldir, &
       lo,hi,bfact,level,vol,cen,rzflag,sdim)
      use global_utility_module

      IMPLICIT NONE

      INTEGER_T sdim
      INTEGER_T lo(sdim),hi(sdim)
      INTEGER_T bfact,rzflag,iside,veldir,nhalf,level
      REAL_T xsten(-nhalf:nhalf,sdim)
      REAL_T dx(sdim)
      REAL_T xlo(sdim)
      INTEGER_T i,j,k,ii
      REAL_T vol
      REAL_T cen(sdim)
      INTEGER_T dir

      if ((iside.ne.-1).and.(iside.ne.1)) then
       print *,"iside invalid"
       stop
      endif
      if ((veldir.lt.1).or.(veldir.gt.sdim)) then
       print *,"veldir invalid"
       stop
      endif
      if (bfact.lt.1) then
       print *,"bfact invalid"
       stop
      endif
      if (nhalf.ne.1) then
       print *,"nhalf invalid cisboxhalf"
       stop
      endif
      if ((sdim.ne.2).and.(sdim.ne.3)) then
       print *,"sdim invalid"
       stop
      endif
      call CISBOX(xsten,nhalf, &
       xlo,dx,i,j,k,lo,hi,bfact,level,vol,cen,rzflag,sdim)
      
      xsten(-iside,veldir)=xsten(0,veldir)
      xsten(0,veldir)=half*(xsten(1,veldir)+xsten(-1,veldir))

      call Box_volumeFAST(bfact,dx,xsten,nhalf,vol,cen,rzflag,sdim)

      return
      end subroutine CISBOXHALF 



       ! centroid is in absolute coordinate system 
      subroutine Box_volumeFAST(bfact,dx,xsten,nhalf, &
        volume,centroid,coord,sdim)
      IMPLICIT NONE

      INTEGER_T sdim,bfact,nhalf
      REAL_T xsten(-nhalf:nhalf,sdim)
      REAL_T dx(sdim)
      REAL_T volume
      REAL_T centroid(sdim)
      INTEGER_T coord,dir
      REAL_T rval,dr
    
      if ((sdim.ne.2).and.(sdim.ne.3)) then
       print *,"sdim invalid Box_volumeFAST"
       stop
      endif 
      if (nhalf.lt.1) then
       print *,"nhalf invalid boxvolumefast"
       stop
      endif
      if (bfact.lt.1) then
       print *,"bfact invalid"
       stop
      endif

      if (sdim.eq.2) then

       if (coord.eq.0) then
        volume=one
        do dir=1,sdim
         volume=volume*(xsten(1,dir)-xsten(-1,dir))
         centroid(dir)=half*(xsten(1,dir)+xsten(-1,dir))
        enddo
       else if (coord.eq.1) then
        volume=one
        do dir=1,sdim
         volume=volume*(xsten(1,dir)-xsten(-1,dir))
         centroid(dir)=half*(xsten(1,dir)+xsten(-1,dir))
        enddo
        rval=centroid(1)
        dr=xsten(1,1)-xsten(-1,1)
        if (rval.ne.zero) then
         volume=volume*two*Pi*abs(rval)
         centroid(1)=rval+dr*dr/(12.0*rval)
        else
         volume=zero
         centroid(1)=zero
        endif
       else if (coord.eq.3) then  ! in: Box_volumeFAST (2d)
        volume=one
        do dir=1,sdim
         volume=volume*(xsten(1,dir)-xsten(-1,dir))
         centroid(dir)=half*(xsten(1,dir)+xsten(-1,dir))
        enddo
        rval=centroid(1)
        dr=xsten(1,1)-xsten(-1,1)
        if (rval.ne.zero) then
         volume=volume*abs(rval)
         centroid(1)=rval+dr*dr/(12.0*rval)
        else
         volume=zero
         centroid(1)=zero
        endif
       else
        print *,"coord invalid"
        stop
       endif

      else if (sdim.eq.3) then
       if (coord.eq.0) then
        volume=one
        do dir=1,sdim
         volume=volume*(xsten(1,dir)-xsten(-1,dir))
         centroid(dir)=half*(xsten(1,dir)+xsten(-1,dir))
        enddo
       else if (coord.eq.3) then ! in: Box_volumeFAST (3d)
        volume=one
        do dir=1,sdim
         volume=volume*(xsten(1,dir)-xsten(-1,dir))
         centroid(dir)=half*(xsten(1,dir)+xsten(-1,dir))
        enddo
        rval=centroid(1)
        dr=xsten(1,1)-xsten(-1,1)
        if (rval.ne.zero) then
         volume=volume*abs(rval)
         centroid(1)=rval+dr*dr/(12.0*rval)
        else
         volume=zero
         centroid(1)=zero
        endif
       else
        print *,"coord invalid"
        stop
       endif
      else
       print *,"dimension bust"
       stop
      endif

      return
      end subroutine Box_volumeFAST

      subroutine Box_volume_super(bfact,dx,xsten0,nhalf0, &
       volume,centroid,sdim,levelrz)
      IMPLICIT NONE

      INTEGER_T sdim,levelrz,bfact,nhalf0,nhalf2
      REAL_T xsten0(-nhalf0:nhalf0,sdim)
      REAL_T xsten2(-1:1,sdim)
      REAL_T dx(sdim)
      REAL_T volume
      REAL_T centroid(sdim)
      INTEGER_T ksten_low,ksten_high,i1,j1,k1,dir,isten
      REAL_T volsten
      REAL_T censten(sdim)

      nhalf2=1

      if (bfact.lt.1) then
       print *,"bfact invalid"
       stop
      endif
      if (nhalf0.lt.3) then
       print *,"nhalf0 invalid"
       stop
      endif

      volume=zero
      do dir=1,sdim
       centroid(dir)=zero
      enddo

      if (sdim.eq.3) then
        ksten_low=-1
        ksten_high=1
      else if (sdim.eq.2) then
        ksten_low=0
        ksten_high=0
      else
        print *,"sdim invalid"
        stop
      endif
 
      do i1=-1,1
      do j1=-1,1
      do k1=ksten_low,ksten_high

       do isten=-1,1
        xsten2(isten,1)=xsten0(isten+2*i1,1)
        xsten2(isten,2)=xsten0(isten+2*j1,2)
        if (sdim.eq.3) then
         xsten2(isten,sdim)=xsten0(isten+2*k1,sdim)
        endif
       enddo ! isten

       call Box_volumeFAST(bfact,dx,xsten2,nhalf2,volsten,censten,levelrz,sdim)

       volume=volume+volsten
       do dir=1,sdim
        centroid(dir)=centroid(dir)+censten(dir)*volsten
       enddo
      enddo
      enddo
      enddo  ! i1,j1,k1

      if (volume.le.zero) then
        print *,"volume invalid"
        stop
      endif
      do dir=1,sdim
        centroid(dir)=centroid(dir)/volume
      enddo

      return
      end subroutine Box_volume_super


 
 
! f(x,y,z)=a+b(x-x0)+c(y-y0)+d(z-z0)
!   3  4  
!   1  2
!
!   7 8
!   5 6
!
!  the 5 tetrahedra:
!   3,1,2,5
!   5,6,2,8
!   3,5,7,8
!   3,4,2,8
!   5,8,3,2
!
! in 2d
! 1,2,3
! 2,3,4

      subroutine get_ntetbox(ntetbox,symmetry_flag,sdim)
      IMPLICIT NONE

      INTEGER_T ntetbox,symmetry_flag,sdim

      if (sdim.eq.2) then
       if (symmetry_flag.eq.0) then
        ntetbox=2
       else if (symmetry_flag.eq.1) then
        ntetbox=4
       else
        print *,"symmetry_flag invalid"
        stop
       endif
      else if (sdim.eq.3) then
       if (symmetry_flag.eq.0) then
        ntetbox=5
       else if (symmetry_flag.eq.1) then
        ntetbox=24
       else
        print *,"symmetry_flag invalid"
        stop
       endif
      else
       print *,"sdim invalid"
       stop
      endif

      return
      end subroutine get_ntetbox


        ! centroid in absolute coordinate system
      subroutine getvolume( &
        bfact,dxgrid,xsten,nhalf, &
        ldata,volume,facearea, &
        centroid,areacentroid,EBVOFTOL,sdim)
      use probcommon_module
      IMPLICIT NONE

      INTEGER_T sdim,bfact,nhalf
      REAL_T EBVOFTOL
      REAL_T ldata(D_DECL(3,3,3))
      REAL_T lnode(4*(sdim-1))
      REAL_T xsten(-nhalf:nhalf,sdim)
      REAL_T dxgrid(sdim)
      REAL_T volume,facearea
      REAL_T volcell
      REAL_T areacentroid(sdim)
      REAL_T centroid(sdim)
      REAL_T cenall(sdim)

      if (bfact.lt.1) then
       print *,"bfact invalid"
       stop
      endif
      if (nhalf.lt.1) then
       print *,"nhalf invalid getvolume"
       stop
      endif
      if (EBVOFTOL.le.zero) then
       print *,"getvolume: EBVOFTOL too small EBVOFTOL=",EBVOFTOL
       stop
      endif
      if (levelrz.eq.0) then
       ! do nothing
      else if (levelrz.eq.1) then
       if (sdim.ne.2) then
        print *,"levelrz invalid get volume"
        stop
       endif
      else if (levelrz.eq.3) then
       ! do nothing
      else
       print *,"levelrz invalid get volume 2"
       stop
      endif
      if (BL_SPACEDIM.ne.sdim) then
       print *,"dimension mismatch"
       stop
      endif
      call data_to_node(ldata,lnode,1,xsten,nhalf,sdim)

      call fast_cell_intersection_grid( &
       bfact,dxgrid,xsten,nhalf, &       
       lnode, &
       volume,centroid,facearea,areacentroid,levelrz, &
       volcell,cenall,sdim)

      if (volcell.le.zero) then
       volume=zero
      else
       volume=volume/volcell
      endif

      if (volume.le.EBVOFTOL) then
       volume=zero
      endif
      if (volume.ge.one-EBVOFTOL) then
       volume=one
      endif

      return
      end subroutine getvolume



         ! "centroid" in absolute coordinate system      
      subroutine getvolumebatch(bfact,dxgrid,xsten,nhalf, &
       ldata,volume,facearea, &
       centroid,areacentroid,nmat,EBVOFTOL,sdim)
      use probcommon_module
      use global_utility_module

      IMPLICIT NONE

      INTEGER_T sdim,bfact,nhalf
      REAL_T EBVOFTOL
      INTEGER_T im,nmat
      REAL_T ldata(D_DECL(3,3,3),nmat)
      REAL_T lnode(4*(sdim-1),nmat)
      REAL_T xsten(-nhalf:nhalf,sdim)
      REAL_T dxgrid(sdim)
      REAL_T volume(nmat)
      REAL_T facearea(nmat)
      REAL_T volcell
      REAL_T areacentroid(nmat,sdim)
      REAL_T centroid(nmat,sdim)
      REAL_T cenall(sdim)
      INTEGER_T dir,iflag

      if (bfact.lt.1) then
       print *,"bfact invalid"
       stop
      endif

      if (EBVOFTOL.le.zero) then
       print *,"getvolumebatch: EBVOFTOL too small EBVOFTOL=",EBVOFTOL
       stop
      endif
      if (nmat.ne.num_materials) then
       print *,"nmat invalid"
       stop
      endif
      if (levelrz.eq.0) then
       ! do nothing
      else if (levelrz.eq.1) then
       if (sdim.ne.2) then
        print *,"sdim invalid"
        stop
       endif
      else if (levelrz.eq.3) then
       ! do nothing
      else
       print *,"levelrz invalid get volume batch"
       stop
      endif 
      if (BL_SPACEDIM.ne.sdim) then
       print *,"dimension mismatch"
       stop
      endif
      if (nhalf.lt.3) then
       print *,"nhalf invalid getvolumebatch"
       stop
      endif

      call data_to_node(ldata,lnode,nmat,xsten,nhalf,sdim)

      call fast_cell_intersection_grid_batch( &
        bfact,dxgrid,xsten,nhalf, &
        lnode, &
        volume,centroid,facearea,areacentroid,levelrz, &
        volcell,cenall,nmat,sdim)

      do im=1,nmat
       if (volcell.le.zero) then
        volume(im)=zero
       else
        volume(im)=volume(im)/volcell
       endif

       if (volume(im).le.EBVOFTOL) then
        volume(im)=zero
       endif
       if (volume(im).ge.one-EBVOFTOL) then
        volume(im)=one
       endif
      enddo

      return
      end subroutine getvolumebatch

      subroutine data_to_node(datasten,datanode,ncomp,xsten,nhalf,sdim)
      IMPLICIT NONE

      INTEGER_T ncomp,nhalf
      INTEGER_T sdim
      REAL_T xsten(-nhalf:nhalf,sdim)
      REAL_T datasten(D_DECL(3,3,3),ncomp)
      REAL_T datanode(4*(sdim-1),ncomp)
      INTEGER_T i,j,k,i1,j1,k1,inode,im,klo,khi,dir
      REAL_T xlonode(sdim)
      REAL_T xhinode(sdim)
      REAL_T xlocell(sdim)
      REAL_T xhicell(sdim)
      REAL_T xlo,xhi,wtprod,wtsum

      if (nhalf.lt.3) then
       print *,"nhalf invalid data to node"
       stop
      endif
      if (sdim.eq.3) then
       klo=-1
       khi=1
      else if (sdim.eq.2) then
       klo=0
       khi=0
      else
       print *,"sdim bust data_to_node"
       stop
      endif

      if (ncomp.lt.1) then
       print *,"ncomp invalid data_to_node"
       stop
      endif

       ! wt=|omega_cell int omega_node|/|omega_node|

      inode=1
      do k=klo,khi,2
      do j=-1,1,2
      do i=-1,1,2
       dir=1
       xlonode(dir)=xsten(i-1,dir)
       xhinode(dir)=xsten(i+1,dir)
       dir=2
       xlonode(dir)=xsten(j-1,dir)
       xhinode(dir)=xsten(j+1,dir)
       if (sdim.eq.3) then
        dir=sdim
        xlonode(dir)=xsten(k-1,dir)
        xhinode(dir)=xsten(k+1,dir)
       endif
       do im=1,ncomp
        datanode(inode,im)=zero
       enddo
       wtsum=zero
       do i1=-1,1
       do j1=-1,1
       do k1=klo,khi
        dir=1
        xlocell(dir)=xsten(2*i1-1,dir)
        xhicell(dir)=xsten(2*i1+1,dir)
        dir=2
        xlocell(dir)=xsten(2*j1-1,dir)
        xhicell(dir)=xsten(2*j1+1,dir)
        if (sdim.eq.3) then
         dir=sdim
         xlocell(dir)=xsten(2*k1-1,dir)
         xhicell(dir)=xsten(2*k1+1,dir)
        endif
        wtprod=one
        do dir=1,sdim
         xlo=max(xlocell(dir),xlonode(dir))
         xhi=min(xhicell(dir),xhinode(dir))
         if (xhi.gt.xlo) then
          wtprod=wtprod*(xhi-xlo)
         else
          wtprod=zero
         endif
        enddo
        do im=1,ncomp
         datanode(inode,im)=datanode(inode,im)+ &
           wtprod*datasten(D_DECL(i1+2,j1+2,k1+2),im)
        enddo
        wtsum=wtsum+wtprod
       enddo
       enddo
       enddo ! i1,j1,k1
       if (wtsum.le.zero) then
        print *,"wtsum invalid in data_to_node"
        stop
       endif
       do im=1,ncomp
        datanode(inode,im)=datanode(inode,im)/wtsum
       enddo
       inode=inode+1
      enddo  
      enddo  
      enddo  ! i,j,k

      return
      end subroutine data_to_node


        ! reconstruction relative to xsten0
        ! volume domain: xsten_grid
        ! shapeflag=0 find volumes within xsten_grid
        ! shapeflag=1 find volumes within xtet

      subroutine fast_cut_cell_intersection( &
        bfact,dx,xsten0,nhalf0, &
        slope,intercept, &
        volume,centroid,area,areacentroid,coord, &
        xsten_grid,nhalf_grid,xtet,shapeflag,sdim)

      use global_utility_module

      IMPLICIT NONE

      INTEGER_T sdim

      REAL_T cum_volume,cum_area
      REAL_T cum_centroid(sdim)
      REAL_T cum_areacentroid(sdim)

      INTEGER_T coord,bfact,nhalf0,nhalf_grid

      INTEGER_T shapeflag
      REAL_T dx(sdim)
      REAL_T xsten0(-nhalf0:nhalf0,sdim)
      REAL_T xsten_grid(-nhalf_grid:nhalf_grid,sdim)
      REAL_T xtet(sdim+1,sdim)
      REAL_T slope(sdim)
      REAL_T intercept

      REAL_T volume,area
      REAL_T areacentroid(sdim)
      REAL_T centroid(sdim)
      INTEGER_T i,j,k,dir,inode
      REAL_T xtarget(sdim)
      REAL_T ls(sdim+1)

!   3  4  
!   1  2
!
!   7 8
!   5 6

      REAL_T xnode(4*(sdim-1),sdim)
      REAL_T phinode(4*(sdim-1))

      INTEGER_T linearcut,fullelementfast,nodedomain

      linearcut=1
      fullelementfast=1

      if ((nhalf0.lt.1).or.(nhalf_grid.lt.1)) then
       print *,"nhalf invalid fast cut cell intersection"
       stop
      endif
      if ((sdim.ne.3).and.(sdim.ne.2)) then
       print *,"sdim invalid fast_cut_cell_intersection"
       stop
      endif

      if (shapeflag.eq.1) then

       do i=1,sdim+1
        do j=1,sdim
         xtarget(j)=xtet(i,j)
        enddo
        call distfunc(bfact,dx,xsten0,nhalf0, &
         intercept,slope,xtarget,ls(i),sdim)
       enddo  ! i=1,sdim+1

       nodedomain=sdim+1
       call intersection_volume( &
        cum_volume,cum_area,cum_centroid,cum_areacentroid, &
        ls,xtet,nodedomain,coord, &
        sdim,fullelementfast,linearcut)

      else if (shapeflag.eq.0) then

       inode=1

       if (sdim.eq.3) then
        do k=-1,1,2
        do j=-1,1,2
        do i=-1,1,2
         do dir=1,sdim
          if (dir.eq.1) then
           xnode(inode,dir)=xsten_grid(i,dir)
          else if (dir.eq.2) then
           xnode(inode,dir)=xsten_grid(j,dir)
          else if (dir.eq.sdim) then
           xnode(inode,dir)=xsten_grid(k,dir)
          else
           print *,"dir invalid face cut cell intersection"
           stop
          endif
          xtarget(dir)=xnode(inode,dir)
         enddo  ! dir
         call distfunc(bfact,dx,xsten0,nhalf0, &
           intercept,slope,xtarget, &
           phinode(inode),sdim)
         
         inode=inode+1
        enddo
        enddo
        enddo  ! i,j,k

       else if (sdim.eq.2) then

        do j=-1,1,2
        do i=-1,1,2
         do dir=1,sdim
          if (dir.eq.1) then
           xnode(inode,dir)=xsten_grid(i,dir)
          else if (dir.eq.2) then
           xnode(inode,dir)=xsten_grid(j,dir)
          else
           print *,"dir invalid face cut cell intersection 2"
           stop
          endif
          xtarget(dir)=xnode(inode,dir)
         enddo  ! dir
         call distfunc(bfact,dx,xsten0,nhalf0, &
           intercept,slope,xtarget, &
           phinode(inode),sdim)
         
         inode=inode+1
        enddo
        enddo ! i,j
       else
        print *,"sdim invalid"
        stop
       endif

       nodedomain=4*(sdim-1)

       if (inode.ne.nodedomain+1) then
        print *,"inode invalid"
        stop
       endif
       call intersection_volume( &
        cum_volume,cum_area,cum_centroid,cum_areacentroid, &
        phinode,xnode,nodedomain,coord, &
        sdim,fullelementfast,linearcut)
      else
       print *,"shapeflag invalid"
       stop
      endif

      volume=cum_volume
      area=cum_area
      do dir=1,sdim
       centroid(dir)=cum_centroid(dir)
       areacentroid(dir)=cum_areacentroid(dir)
      enddo

      return
      end subroutine fast_cut_cell_intersection

        ! returns centroid in absolute coordinate system

      subroutine fast_cell_intersection_grid_batch( &
        bfact,dxgrid,xsten0,nhalf0, &
        lnodebatch, &
        volumedark,centroiddark, &
        area,areacentroid,coord,volall,cenall,nmat,sdim)

      IMPLICIT NONE

      INTEGER_T sdim
      INTEGER_T coord,nmat,nhalf0

      REAL_T cum_volume,cum_area
      REAL_T cum_centroid(sdim)
      REAL_T cum_areacentroid(sdim)

      INTEGER_T bfact
      REAL_T lnodebatch(4*(sdim-1),nmat)
      REAL_T lnode(4*(sdim-1))
      REAL_T xnode(4*(sdim-1),sdim)
      REAL_T xsten0(-nhalf0:nhalf0,sdim)
      REAL_T dxgrid(sdim)
      REAL_T volumedark(nmat)
      REAL_T centroiddark(nmat,sdim)
      REAL_T volall
      REAL_T area(nmat)
      REAL_T cenall(sdim)
      REAL_T areacentroid(nmat,sdim)
      INTEGER_T i,j,k
      INTEGER_T dir
      INTEGER_T inode,im
      INTEGER_T linearcut,fullelementfast,nodedomain

      if (bfact.lt.1) then
       print *,"bfact invalid"
       stop
      endif
      if ((nmat.lt.1).or.(nmat.gt.MAX_NUM_MATERIALS)) then
       print *,"nmat invalid fast cell intersection grid batch"
       stop
      endif

      linearcut=0
      fullelementfast=1
      nodedomain=4*(sdim-1)

      if ((sdim.ne.3).and.(sdim.ne.2)) then
       print *,"sdim bust fast cell intersection grid batch"
       stop
      endif
      if (nhalf0.lt.1) then
       print *,"nhalf0 invalid"
       stop
      endif

      inode=1

      call Box_volumeFAST(bfact,dxgrid,xsten0,nhalf0,volall,cenall,coord,sdim)

      if (volall.gt.zero) then
     
       if (sdim.eq.3) then
        do k=-1,1,2
        do j=-1,1,2
        do i=-1,1,2
         xnode(inode,1)=xsten0(i,1)
         xnode(inode,2)=xsten0(j,2)
         xnode(inode,sdim)=xsten0(k,sdim)
         inode=inode+1
        enddo
        enddo
        enddo
       else if (sdim.eq.2) then
        do j=-1,1,2
        do i=-1,1,2
         xnode(inode,1)=xsten0(i,1)
         xnode(inode,2)=xsten0(j,2)
         inode=inode+1
        enddo
        enddo
       else 
        print *,"sdim invalid"
        stop
       endif

       if (inode.ne.nodedomain+1) then
        print *,"inode invalid"
        stop
       endif

       do im=1,nmat
        do inode=1,nodedomain
         lnode(inode)=lnodebatch(inode,im)
        enddo        
        call intersection_volume( &
         cum_volume,cum_area,cum_centroid,cum_areacentroid, &
         lnode,xnode,nodedomain,coord, &
         sdim,fullelementfast,linearcut)

        volumedark(im)=cum_volume
        area(im)=cum_area
        do j=1,sdim
         centroiddark(im,j)=cum_centroid(j)
         areacentroid(im,j)=cum_areacentroid(j)
        enddo
       enddo  ! im

      else if (volall.eq.zero) then

       print *,"WARNING:volall cannot be zero: fast_cell_intersection_grid"
       do im=1,nmat
        volumedark(im)=zero
        area(im)=zero
        do j=1,sdim
         centroiddark(im,j)=zero
         areacentroid(im,j)=zero
        enddo
       enddo  ! im

      else if (volall.lt.zero) then
       print *,"volall invalid"
       stop
      endif

      return
      end subroutine fast_cell_intersection_grid_batch



        ! returns centroid in absolute coordinate system

      subroutine fast_cell_intersection_grid( &
        bfact,dxgrid,xsten0,nhalf0, &
        lnode, &
        volumedark,centroiddark, &
        area,areacentroid,coord,volall,cenall,sdim)

      IMPLICIT NONE

      INTEGER_T sdim
      INTEGER_T coord

      REAL_T cum_volume,cum_area
      REAL_T cum_centroid(sdim)
      REAL_T cum_areacentroid(sdim)

      INTEGER_T bfact,nhalf0
      REAL_T lnode(4*(sdim-1))
      REAL_T xnode(4*(sdim-1),sdim)
      REAL_T xsten0(-nhalf0:nhalf0,sdim)
      REAL_T dxgrid(sdim)
      REAL_T volumedark
      REAL_T centroiddark(sdim)
      REAL_T volall,area
      REAL_T cenall(sdim)
      REAL_T areacentroid(sdim)
      INTEGER_T i,j,k
      INTEGER_T dir
      INTEGER_T inode
      INTEGER_T linearcut,fullelementfast,nodedomain

      if (nhalf0.lt.1) then
       print *,"nhalf0 invalid"
       stop
      endif
      linearcut=0
      fullelementfast=1
      nodedomain=4*(sdim-1)

      if (bfact.lt.1) then
       print *,"bfact invalid"
       stop
      endif
      if ((sdim.ne.3).and.(sdim.ne.2)) then
       print *,"sdim bust fast cell intersection grid"
       stop
      endif

      inode=1

      call Box_volumeFAST(bfact,dxgrid,xsten0,nhalf0,volall,cenall,coord,sdim)

      if (volall.gt.zero) then
     
       if (sdim.eq.3) then
        do k=-1,1,2
        do j=-1,1,2
        do i=-1,1,2 
         xnode(inode,1)=xsten0(i,1)
         xnode(inode,2)=xsten0(j,2)
         xnode(inode,sdim)=xsten0(k,sdim)
         inode=inode+1
        enddo
        enddo
        enddo
       else if (sdim.eq.2) then
        do j=-1,1,2
        do i=-1,1,2
         xnode(inode,1)=xsten0(i,1)
         xnode(inode,2)=xsten0(j,2)
         inode=inode+1
        enddo
        enddo
       else 
        print *,"sdim invalid"
        stop
       endif

       if (inode.ne.nodedomain+1) then
        print *,"inode invalid"
        stop
       endif
       call intersection_volume( &
        cum_volume,cum_area,cum_centroid,cum_areacentroid, &
        lnode,xnode,nodedomain,coord, &
        sdim,fullelementfast,linearcut)

       volumedark=cum_volume
       area=cum_area
       do j=1,sdim
        centroiddark(j)=cum_centroid(j)
        areacentroid(j)=cum_areacentroid(j)
       enddo

      else if (volall.eq.zero) then
       print *,"WARNING:volall cannot be zero: fast_cell_intersection_grid"
       volumedark=zero
       area=zero
       do j=1,sdim
        centroiddark(j)=zero
        areacentroid(j)=zero
       enddo
      else if (volall.lt.zero) then
       print *,"volall invalid"
       stop
      endif

      return
      end subroutine fast_cell_intersection_grid



end module geometry_intersect_module


module MOF_routines_module

      INTEGER_T :: MOF_DEBUG_RECON_COUNT
       ! 1=>output errors when fastflag=0 or 1
       ! 2=>output errors when fastflag=0 
      INTEGER_T :: MOF_DEBUG_RECON
      INTEGER_T :: MOF_TURN_OFF_LS

contains

      subroutine get_order_algorithm(order_algorithm_out,nmat)
      use geometry_intersect_module

      IMPLICIT NONE

      INTEGER_T im,nmat
      INTEGER_T order_algorithm_out(nmat)

#include "mofdata.H"

      if ((nmat.lt.1).or.(nmat.gt.MAX_NUM_MATERIALS)) then
       print *,"nmat invalid  get order algorithm"
       print *,"nmat= ",nmat
       stop
      endif

      do im=1,nmat
       order_algorithm_out(im)=order_algorithm(im)
      enddo

      return
      end subroutine get_order_algorithm

      subroutine get_MOFITERMAX(MOFITERMAX_out)
      IMPLICIT NONE

      INTEGER_T MOFITERMAX_out

#include "mofdata.H"

      MOFITERMAX_out=MOFITERMAX

      return
      end subroutine get_MOFITERMAX

      subroutine set_MOFITERMAX(MOFITERMAX_in)
      IMPLICIT NONE

      INTEGER_T MOFITERMAX_in

#include "mofdata.H"

      MOFITERMAX=MOFITERMAX_in
      if ((MOFITERMAX.lt.0).or.(MOFITERMAX.gt.50)) then
       print *,"MOFITERMAX invalid in set mofitermax"
       stop
      endif

      return
      end subroutine set_MOFITERMAX


       ! order_algorithm=0 => try different combinations and
       ! choose combination with smallest MOF error
      subroutine set_order_algorithm(order_algorithm_in,nmat)
      use geometry_intersect_module

      IMPLICIT NONE

      INTEGER_T nmat
      INTEGER_T order_algorithm_in(nmat)
      INTEGER_T im

#include "mofdata.H"

      if ((nmat.lt.1).or.(nmat.gt.MAX_NUM_MATERIALS)) then
       print *,"nmat invalid set order algorithm"
       print *,"nmat= ",nmat
       stop
      endif

      do im=1,nmat
       order_algorithm(im)=order_algorithm_in(im)
       if (order_algorithm(im).lt.0) then
        print *,"order_alg bust"
        stop
       endif
      enddo

      return
      end subroutine set_order_algorithm



! intercept=-maxphi  => refvfrac=0
! intercept=-minphi  => refvfrac=1
      subroutine multi_phi_bounds( &
        bfact,dx,xsten,nhalf, &
        slope,xtetlist,nlist,nmax, &
        minphi,maxphi,sdim)

      use global_utility_module

      IMPLICIT NONE

      INTEGER_T nlist,nmax,sdim,bfact,nhalf
      REAL_T xtetlist(sdim+1,sdim,nmax)
      REAL_T xsten(-nhalf:nhalf,sdim)
      REAL_T slope(sdim)
      REAL_T xtarget(sdim)
      REAL_T dx(sdim)
      INTEGER_T dir,i,j,n
      REAL_T minphi,maxphi,intercept,dist

      if (nhalf.lt.1) then
       print *,"nhalf invalid multi phi bounds 3d"
       stop
      endif
      if ((sdim.ne.2).and.(sdim.ne.3)) then
       print *,"sdim invalid multi_phi_bounds"
       stop
      endif

      intercept=zero
      minphi=1.0E+10
      maxphi=-1.0E+10
      do n=1,nlist
       do i=1,sdim+1
        do dir=1,sdim
         xtarget(dir)=xtetlist(i,dir,n)
        enddo

        call distfunc(bfact,dx,xsten,nhalf, &
         intercept,slope,xtarget,dist,sdim)

        if (dist.lt.minphi) then
         minphi=dist
        endif
        if (dist.gt.maxphi) then
         maxphi=dist
        endif
       enddo  ! sweeping nodes of triangle
      enddo ! sweeping triangles

      if (((minphi.eq.zero).and.(maxphi.eq.zero)).or. &
          (minphi.ge.maxphi)) then
       print *,"cannot have zero slope"
       print *,"minphi=",minphi
       print *,"maxphi=",maxphi
       print *,"slopexyz=",slope(1),slope(2),slope(sdim)
       print *,"xcell xyz=",xsten(0,1),xsten(0,2),xsten(0,sdim)
       stop
      endif

      return
      end subroutine multi_phi_bounds 
 
! phi>0 in the material.  n points into the phi>0 region.
      subroutine closest(xint,x,nn,phi,coord,sdim)
      IMPLICIT NONE

      INTEGER_T sdim
      REAL_T xint(sdim),x(sdim)
      REAL_T nn(sdim)
      REAL_T phi
      INTEGER_T coord
      REAL_T rr
      INTEGER_T dir

      if (coord.eq.0) then
       ! do nothing
      else if (coord.eq.1) then
       if (sdim.ne.2) then
        print *,"dimension bust"
       endif
      else if (coord.eq.3) then
       ! do nothing
      else
       print *,"coord invalid"
       stop
      endif

      do dir=1,sdim
       xint(dir)=x(dir)-nn(dir)*phi
      enddo

      return
      end subroutine closest

 
      subroutine dist2or3D(x1,x2,dist,sdim)
      IMPLICIT NONE

      INTEGER_T sdim,dir
      REAL_T x1(sdim),x2(sdim)
      REAL_T dist

      if ((sdim.ne.2).and.(sdim.ne.3)) then
       print *,"dimension bust"
       stop
      endif
      dist=zero
      do dir=1,sdim
       dist=dist+(x1(dir)-x2(dir))**2
      enddo
      dist=sqrt(dist)

      return
      end subroutine dist2or3D




! find volume, centroid and area of intersection of a line with a 
! collection of triangles.
      subroutine multi_cell_intersection( &
        bfact,dx,xsten,nhalf, &
        slope,intercept, &
        volume,centroid,area,coord, &
        xtetlist,nlist,nmax,sdim)

      use global_utility_module
      use geometry_intersect_module

      IMPLICIT NONE

      INTEGER_T coord,bfact,nhalf

      INTEGER_T nlist,nmax,sdim
      REAL_T dx(sdim)
      REAL_T xsten(-nhalf:nhalf,sdim)
      REAL_T xtetlist(sdim+1,sdim,nmax)
      REAL_T slope(sdim)
      REAL_T intercept

      REAL_T volume,area
      REAL_T volumelist,arealist
      REAL_T areacentroid(sdim)
      REAL_T centroid(sdim)
      REAL_T centroidlist(sdim)
      INTEGER_T i,j,n
      REAL_T xx(sdim+1,sdim)
      REAL_T xtarget(sdim)
      REAL_T ls(sdim+1)

      if (nhalf.lt.1) then
       print *,"nhalf invalid multi cell intersection"
       stop
      endif
      if ((sdim.ne.3).and.(sdim.ne.2)) then
       print *,"sdim invalid multi_cell_intersection"
       stop
      endif

      do j=1,sdim
       centroid(j)=zero
      enddo
      area=zero
      volume=zero

      do n=1,nlist
       do i=1,sdim+1
        do j=1,sdim
         xx(i,j)=xtetlist(i,j,n)
         xtarget(j)=xx(i,j)
        enddo
        call distfunc(bfact,dx,xsten,nhalf, &
         intercept,slope,xtarget,ls(i),sdim)
       enddo

       if (sdim.eq.3) then
        call intersection_volumeXYZ(ls,xx,volumelist, &
         centroidlist,arealist,areacentroid,coord,sdim)
       else if (sdim.eq.2) then
        call int_volumeXYorRZ(ls,xx,volumelist, &
         centroidlist,arealist,areacentroid,coord,sdim)
       else
        print *,"sdim invalid"
        stop
       endif

       volume=volume+volumelist
       area=area+arealist
       do j=1,sdim
        centroid(j)=centroid(j)+centroidlist(j)*volumelist
       enddo
      enddo ! n

      if (volume.gt.zero) then
       do j=1,sdim
        centroid(j)=centroid(j)/volume
       enddo
      else
       volume=zero
       do j=1,sdim
        centroid(j)=zero
       enddo
      endif

      return
      end subroutine multi_cell_intersection




       ! centroid relative to absolute coordinate system

      subroutine multi_ff(bfact,dx,xsten0,nhalf0, &
        ff,slope,intercept, &
        continuous_mof, &
        arean,coord, &
        vtarget,xtetlist,centroid,nlist,nmax,fastflag,sdim)

      use geometry_intersect_module

      IMPLICIT NONE

      INTEGER_T continuous_mof,bfact,nhalf0
      INTEGER_T nlist,nmax,sdim,fastflag
      REAL_T xtetlist(sdim+1,sdim,nmax)
      REAL_T slope(sdim)
      REAL_T intercept
      REAL_T xsten0(-nhalf0:nhalf0,sdim)
      REAL_T xsten2(-1:1,sdim)
      INTEGER_T nhalf2
      REAL_T dx(sdim)
      INTEGER_T dir

      REAL_T ff,voln,arean,vtarget,volcell
      REAL_T cencell(sdim)
      REAL_T areacentroidn(sdim)
      INTEGER_T coord
      REAL_T centroid(sdim)
      INTEGER_T shapeflag
      REAL_T xtet(sdim+1,sdim)
      INTEGER_T ksten_low,ksten_high
      REAL_T volsten,areasten
      REAL_T censten(sdim)
      REAL_T areacentroidsten(sdim)
      INTEGER_T i1,j1,k1,isten

      nhalf2=1

      if ((sdim.ne.3).and.(sdim.ne.2)) then
       print *,"sdim invalid multi_ff"
       stop
      endif
      if ((continuous_mof.ne.0).and. &
          (continuous_mof.ne.1).and. &
          (continuous_mof.ne.2)) then
       print *,"continuous_mof invalid"
       stop
      endif
      if (nhalf0.lt.1) then
       print *,"nhalf0 invalid"
       stop
      endif

      if ((continuous_mof.eq.0).or.(continuous_mof.eq.2)) then
       call Box_volumeFAST(bfact,dx,xsten0,nhalf0,volcell, &
        cencell,coord,sdim)
      else if (continuous_mof.eq.1) then
       call Box_volume_super(bfact,dx,xsten0,nhalf0, &
        volcell,cencell,sdim,coord)
      else
       print *,"continuous_mof invalid"
       stop
      endif

      if (volcell.le.zero) then
       print *,"volcell bust"
       stop
      endif

      if (fastflag.eq.0) then
         ! xsten0 used for LS dist.
       call multi_cell_intersection(bfact,dx,xsten0,nhalf0, &
         slope,intercept,voln, &
         centroid,arean,coord,xtetlist,nlist,nmax,sdim)
      else if (fastflag.eq.1) then
       shapeflag=0

       if ((continuous_mof.eq.0).or.(continuous_mof.eq.2)) then
        call fast_cut_cell_intersection( &
         bfact,dx,xsten0,nhalf0, &
         slope,intercept, &
         voln,centroid,arean,areacentroidn,coord, &
         xsten0,nhalf0,xtet,shapeflag,sdim) 
       else if (continuous_mof.eq.1) then
        voln=zero
        arean=zero
        do dir=1,sdim
         centroid(dir)=zero
         areacentroidn(dir)=zero
        enddo 

        if (sdim.eq.3) then
         ksten_low=-1
         ksten_high=1
        else if (sdim.eq.2) then
         ksten_low=0
         ksten_high=0
        else
         print *,"sdim invalid"
         stop
        endif

        do i1=-1,1
        do j1=-1,1
        do k1=ksten_low,ksten_high

         do isten=-1,1
          xsten2(isten,1)=xsten0(isten+2*i1,1)
          xsten2(isten,2)=xsten0(isten+2*j1,2)
          if (sdim.eq.3) then
           xsten2(isten,sdim)=xsten0(isten+2*k1,sdim)
          endif
         enddo ! isten

         call fast_cut_cell_intersection( &
          bfact,dx,xsten0,nhalf0, &
          slope,intercept, &
          volsten,censten,areasten, &
          areacentroidsten,coord, &
          xsten2,nhalf2,xtet,shapeflag,sdim) 

         voln=voln+volsten
         arean=arean+areasten
         do dir=1,sdim
          centroid(dir)=centroid(dir)+volsten*censten(dir)
          areacentroidn(dir)=areacentroidn(dir)+ &
            areasten*areacentroidsten(dir) 
         enddo
        enddo
        enddo
        enddo  ! i1,j1,k1
        do dir=1,sdim
         if (voln.gt.zero) then
          centroid(dir)=centroid(dir)/voln
         else
          centroid(dir)=zero
         endif
         if (arean.gt.zero) then
          areacentroidn(dir)=areacentroidn(dir)/arean
         else
          areacentroidn(dir)=zero
         endif
        enddo ! dir 

       else
        print *,"continuous_mof invalid"
        stop
       endif
      else
       print *,"fastflag invalid multi_ff"
       stop
      endif

      ff=(voln-vtarget)/volcell

      return
      end subroutine multi_ff


      subroutine scale_MOF_variables( &
        bfact,dx,xsten,nhalf, &
        refcentroid,levelrz, &
        dx_scale,xsten_scale, &
        refcentroid_scale, &
        sdim,maxdx)
      IMPLICIT NONE

      INTEGER_T sdim,levelrz,bfact,nhalf
      REAL_T maxdx
      REAL_T refcentroid(sdim)
      REAL_T refcentroid_scale(sdim)
      REAL_T xsten(-nhalf:nhalf,sdim)
      REAL_T xsten_scale(-nhalf:nhalf,sdim)
      REAL_T dx(sdim)
      REAL_T dx_scale(sdim)
      INTEGER_T dir,i

      if (nhalf.lt.1) then
       print *,"nhalf invalid scale mof variables"
       stop
      endif

      if ((sdim.ne.3).and.(sdim.ne.2)) then
       print *,"sdim invalid"
       stop
      endif

      if ((levelrz.eq.1).or.(levelrz.eq.3)) then
       if ((levelrz.eq.1).and.(sdim.ne.2)) then
        print *,"sdim invalid"
        stop
       endif
       maxdx=one
       do dir=1,sdim
        do i=-nhalf,nhalf
         xsten_scale(i,dir)=xsten(i,dir)
        enddo
        dx_scale(dir)=dx(dir)

        refcentroid_scale(dir)=refcentroid(dir)
       enddo ! dir

      else if (levelrz.eq.0) then

       maxdx=xsten(1,1)-xsten(-1,1)
       do dir=2,sdim
        if (xsten(1,dir)-xsten(-1,dir).gt.maxdx) then
         maxdx=xsten(1,dir)-xsten(-1,dir)
        endif
       enddo
       do dir=1,sdim
        do i=-nhalf,nhalf
         xsten_scale(i,dir)=xsten(i,dir)/maxdx
        enddo
        dx_scale(dir)=dx(dir)/maxdx
        refcentroid_scale(dir)=refcentroid(dir)/maxdx
       enddo ! dir

      else
       print *,"levelrz invalid scale mof variables"
       stop
      endif
          
      return
      end subroutine scale_MOF_variables


      subroutine scale_VOF_variables( &
        bfact,dx,xsten,nhalf, &
        intercept,levelrz, &
        dx_scale,xsten_scale, &
        sdim,maxdx)
      IMPLICIT NONE

      INTEGER_T sdim,levelrz,bfact,nhalf
      REAL_T intercept,maxdx
      REAL_T dx(sdim)
      REAL_T xsten(-nhalf:nhalf,sdim)
      REAL_T dx_scale(sdim)
      REAL_T xsten_scale(-nhalf:nhalf,sdim)
      INTEGER_T dir,i

      if (nhalf.lt.1) then
       print *,"nhalf invalid scale vof variables"
       stop
      endif
      if ((sdim.ne.3).and.(sdim.ne.2)) then
       print *,"sdim invalid"
       stop
      endif

      if ((levelrz.eq.1).or.(levelrz.eq.3)) then
       if ((levelrz.eq.1).and.(sdim.ne.2)) then
        print *,"sdim bust in scale_VOF_variables"
        stop
       endif
       maxdx=one
       do dir=1,sdim
        do i=-nhalf,nhalf
         xsten_scale(i,dir)=xsten(i,dir)
        enddo
        dx_scale(dir)=dx(dir)
       enddo ! dir
      else if (levelrz.eq.0) then
       maxdx=xsten(1,1)-xsten(-1,1)
       do dir=2,sdim
        if (xsten(1,dir)-xsten(-1,dir).gt.maxdx) then
         maxdx=xsten(1,dir)-xsten(-1,dir)
        endif
       enddo
       do dir=1,sdim
        do i=-nhalf,nhalf
         xsten_scale(i,dir)=xsten(i,dir)/maxdx
        enddo
        dx_scale(dir)=dx(dir)/maxdx
       enddo
       intercept=intercept/maxdx
      else
       print *,"levelrz invalid scale vof variables"
       stop
      endif
          
      return
      end subroutine scale_VOF_variables



         ! centroid in absolute coordinate system

      subroutine multi_find_intercept( &
       bfact,dx,xsten0,nhalf0, &
       slope,intercept, &
       continuous_mof, &
       xtetlist,nlist, &
       nmax,vfrac,levelrz, &
       use_initial_guess, &
       centroid,fastflag,sdim)

      use global_utility_module
      use geometry_intersect_module

      IMPLICIT NONE

      INTEGER_T continuous_mof,bfact,nhalf0
      INTEGER_T nlist,nmax,sdim,fastflag
      REAL_T xtetlist(sdim+1,sdim,nmax)
      REAL_T slope(sdim)
      REAL_T intercept,vfrac
      REAL_T xsten0(-nhalf0:nhalf0,sdim)
      REAL_T xsten0_scale(-nhalf0:nhalf0,sdim)
      REAL_T xsten2_scale(-1:1,sdim)
      REAL_T dx(sdim)

      INTEGER_T levelrz,niter,maxiter
      INTEGER_T i,j,k,dir
      REAL_T minphi,maxphi
      REAL_T volcell
      REAL_T centroid(sdim)
      REAL_T cencell(sdim)
      REAL_T arean
      REAL_T err,moftol
      REAL_T vtarget,fc
      INTEGER_T debug_root
      INTEGER_T use_initial_guess
      REAL_T err_default,fc_default,arean_default,intercept_default
      REAL_T volcut
      REAL_T cencut(sdim)
      REAL_T xtarget(sdim)
      REAL_T vfrac_normalize
      REAL_T null_intercept,dist
      INTEGER_T klo_stencil,khi_stencil
      INTEGER_T i1,j1,k1
      REAL_T dx_scale(sdim)
      REAL_T maxdx_scale
      INTEGER_T nn,ii,isten
      REAL_T intercept_upper,intercept_lower
      REAL_T intercept_test,aa,bb,fa,fb

      if (nhalf0.lt.1) then
       print *,"nhalf0 invalid"
       stop
      endif
      if ((sdim.ne.3).and.(sdim.ne.2)) then
       print *,"sdim invalid multi_find_intercept"
       stop
      endif

      if (sdim.eq.2) then
       klo_stencil=0
       khi_stencil=0
      else if (sdim.eq.3) then
       klo_stencil=-1
       khi_stencil=1
      else
       print *,"sdim invalid"
       stop
      endif

      if ((continuous_mof.ne.0).and. &
          (continuous_mof.ne.1).and. &
          (continuous_mof.ne.2)) then
       print *,"continuous_mof invalid"
       stop
      endif
      if (bfact.lt.1) then
       print *,"bfact invalid"
       stop
      endif

       ! intercept scaled directly
      call scale_VOF_variables( &
        bfact,dx,xsten0,nhalf0, &
        intercept,levelrz, &
        dx_scale,xsten0_scale,  &
        sdim, &
        maxdx_scale)

      if (fastflag.eq.0) then
       do nn=1,nlist
       do ii=1,sdim+1
       do dir=1,sdim
        xtetlist(ii,dir,nn)=xtetlist(ii,dir,nn)/maxdx_scale
       enddo 
       enddo 
       enddo 
      else if (fastflag.eq.1) then
       ! do nothing
      else
       print *,"fastflag invalid multi_find_intercept fastflag=",fastflag
       stop
      endif

      debug_root=0

      maxiter=100
      moftol=INTERCEPT_TOL

! phi=n dot (x-x0)+int
! find max,min n dot (x-x0)
! if fastflag=0, 
!  search the vertices of all triangles that make up the "cut" domain.

      if ((continuous_mof.eq.0).or.(continuous_mof.eq.2)) then
       call Box_volumeFAST(bfact,dx_scale,xsten0_scale,nhalf0, &
        volcell,cencell,levelrz,sdim)
      else if (continuous_mof.eq.1) then
       call Box_volume_super(bfact,dx_scale,xsten0_scale,nhalf0, &
        volcell,cencell,sdim,levelrz)
      else
       print *,"continuous_mof invalid"
       stop
      endif

       ! at each node x_i one solves:
       !  n dot (x_i-x0) + b_i = 0
       !  b_i=-n dot (x_i-x0)
       !  b_Lower_bound=min b_i = - maxphi
       !  b_Upper_bound=max b_i = - minphi
      if (fastflag.eq.0) then

       call get_cut_geom3D(xtetlist,nlist,nmax,levelrz,volcut, &
         cencut,sdim)

       call multi_phi_bounds( &
         bfact,dx_scale,xsten0_scale,nhalf0, &
         slope, &
         xtetlist,nlist,nmax, &
         minphi,maxphi,sdim)

      else if (fastflag.eq.1) then
       volcut=volcell
       do dir=1,sdim
        cencut(dir)=cencell(dir)
       enddo

       minphi=1.0E+10
       maxphi=-1.0E+10
       null_intercept=zero

       if ((continuous_mof.eq.0).or. &
           (continuous_mof.eq.2)) then

        do k=klo_stencil,khi_stencil,2
        do j=-1,1,2
        do i=-1,1,2
         dir=1
         xtarget(dir)=xsten0_scale(i,dir)
         dir=2
         xtarget(dir)=xsten0_scale(j,dir)

         if (sdim.eq.3) then
          dir=sdim
          xtarget(dir)=xsten0_scale(k,dir)
         endif

         call distfunc(bfact,dx_scale,xsten0_scale,nhalf0, &
          null_intercept,slope, &
          xtarget,dist,sdim)

         if (dist.lt.minphi) then
          minphi=dist
         endif
         if (dist.gt.maxphi) then
          maxphi=dist
         endif
         
        enddo
        enddo
        enddo  ! i,j,k

        if (((minphi.eq.zero).and.(maxphi.eq.zero)).or. &
            (minphi.ge.maxphi)) then
         print *,"cannot have zero slope"
         print *,"fastflag=",fastflag
         print *,"continuous_mof=",continuous_mof
         print *,"minphi=",minphi
         print *,"maxphi=",maxphi
         print *,"slopexyz=",slope(1),slope(2),slope(sdim)
         print *,"xsten(0) xyz=",xsten0(0,1),xsten0(0,2),xsten0(0,sdim)
         print *,"xsten_scale(0) xyz=", &
          xsten0_scale(0,1),xsten0_scale(0,2),xsten0_scale(0,sdim)
         stop
        endif

       else if (continuous_mof.eq.1) then

        do i1=-1,1
        do j1=-1,1
        do k1=klo_stencil,khi_stencil

         do isten=-1,1
          xsten2_scale(isten,1)=xsten0_scale(isten+2*i1,1)
          xsten2_scale(isten,2)=xsten0_scale(isten+2*j1,2)
          if (sdim.eq.3) then
           xsten2_scale(isten,sdim)=xsten0_scale(isten+2*k1,sdim)
          endif
         enddo ! isten

         do k=klo_stencil,khi_stencil,2
         do j=-1,1,2
         do i=-1,1,2

          dir=1
          xtarget(dir)=xsten2_scale(i,dir)
          dir=2
          xtarget(dir)=xsten2_scale(j,dir)

          if (sdim.eq.3) then
           dir=sdim
           xtarget(dir)=xsten2_scale(k,dir)
          endif

          call distfunc(bfact,dx_scale,xsten0_scale,nhalf0, &
           null_intercept,slope, &
           xtarget,dist,sdim)

          if (dist.lt.minphi) then
           minphi=dist
          endif
          if (dist.gt.maxphi) then
           maxphi=dist
          endif
         enddo
         enddo
         enddo  ! i,j,k

        enddo
        enddo
        enddo  ! i1,j1,k1

       else
        print *,"continuous_mof invalid"
        stop
       endif

      else
       print *,"fastflag invalid multi_find_intercept 2"
       stop
      endif

      do dir=1,sdim
       centroid(dir)=cencut(dir)
      enddo

      intercept_lower=-maxphi
      intercept_upper=-minphi

       ! volcut is the uncaptured volume of the cell.
      if ((volcut.le.zero).and.(vfrac.gt.MLSVOFTOL)) then
       print *,"ERROR: volcut<=0 and vfrac>mslvoftol"
       stop
      endif

      if (vfrac.le.MLSVOFTOL) then
       intercept=intercept_lower
      else if (vfrac.ge.volcut/volcell-MLSVOFTOL) then
       intercept=intercept_upper
      else
       vtarget=volcell*vfrac

! solve f(xx)=0 where f(xx)=(V(n dot (x-x0)+intercept-xx)-Vtarget)/volcell
! maxphi -> (maxphi-minphi)vfrac 
! minphi -> (minphi-maxphi)(1-vfrac)

       vfrac_normalize=vfrac*volcell/volcut
       if ((vfrac_normalize.le.zero).or.(vfrac_normalize.ge.one)) then
        print *,"ERROR: vfrac_normalize out of range"
        stop
       endif

       intercept_default=intercept_lower*(one-vfrac_normalize)+ &
           intercept_upper*vfrac_normalize

         ! fc_default=(voln-vtarget)/volcell
       call multi_ff(bfact,dx_scale,xsten0_scale,nhalf0, &
        fc_default,slope,intercept_default, &
        continuous_mof, &
        arean_default,levelrz,vtarget,xtetlist,centroid,nlist,nmax, &
        fastflag,sdim)
       err_default=abs(fc_default)

       if (use_initial_guess.eq.0) then

        intercept=intercept_default
        arean=arean_default
        err=err_default
        fc=fc_default

       else if (use_initial_guess.eq.1) then

        if ((intercept.lt.intercept_lower).or. &
            (intercept.gt.intercept_upper)) then

         intercept=intercept_default
         arean=arean_default
         err=err_default
         fc=fc_default

        else

          ! fc=(voln-vtarget)/volcell
         call multi_ff(bfact,dx_scale,xsten0_scale,nhalf0, &
          fc,slope,intercept, &
          continuous_mof, &
          arean, &
          levelrz,vtarget,xtetlist,centroid,nlist,nmax, &
          fastflag,sdim)
         err=abs(fc)
         if ((err.ge.err_default).or.(arean.eq.zero)) then
          intercept=intercept_default
          arean=arean_default
          err=err_default
          fc=fc_default
         endif
 
        endif

       else
        print *,"use_initial_guess invalid multi_find_intercept"
        print *,"use_initial_guess = ",use_initial_guess
        stop
       endif

        ! err=abs(fc)=(voln-vtarget)/volcell
       if (err.gt.moftol) then

        niter=0
        do while ((niter.lt.maxiter).and.(err.gt.moftol)) 
         if (arean.gt.zero) then

          intercept_test=intercept-fc*volcell/arean
          if ((intercept_test.le.intercept_lower).or. &
              (intercept_test.ge.intercept_upper).or. &
              (niter.ge.maxiter-1)) then
           aa=intercept_lower
           bb=intercept_upper
           niter=0
           do while ((niter.lt.maxiter).and.(err.gt.moftol))

            if (niter.eq.0) then
             call multi_ff(bfact,dx_scale,xsten0_scale,nhalf0, &
              fa,slope,aa, &
              continuous_mof, &
              arean,levelrz,vtarget,xtetlist,centroid,nlist,nmax, &
              fastflag,sdim)
             call multi_ff(bfact,dx_scale,xsten0_scale,nhalf0, &
              fb,slope,bb, &
              continuous_mof, &
              arean,levelrz,vtarget,xtetlist,centroid,nlist,nmax, &
              fastflag,sdim)
            else if (niter.gt.0) then
             ! do nothing
            else
             print *,"niter invalid"
             stop
            endif

            if (abs(fa).le.moftol) then
             intercept=aa
             err=zero
            else if (abs(fb).le.moftol) then
             intercept=bb
             err=zero
            else if (fa*fb.lt.zero) then
             intercept=half*(aa+bb)
             call multi_ff(bfact,dx_scale,xsten0_scale,nhalf0, &
              fc,slope,intercept, &
              continuous_mof, &
              arean,levelrz,vtarget,xtetlist,centroid,nlist,nmax, &
              fastflag,sdim)
             err=abs(fc)
             if (fa*fc.lt.zero) then
              bb=intercept
              fb=fc
             else
              aa=intercept
              fa=fc
             endif
            else
             print *,"signs of fa and fb are inconsistent"
             stop
            endif

            niter=niter+1
            if (debug_root.eq.1) then
             print *,"bisection: niter,intercept,fc ",niter,intercept,fc
            endif  
           enddo ! bisection while
          else
           intercept=intercept_test
           call multi_ff(bfact,dx_scale,xsten0_scale,nhalf0, &
            fc,slope,intercept, &
            continuous_mof, &
            arean,levelrz,vtarget,xtetlist,centroid,nlist,nmax, &
            fastflag,sdim)
           err=abs(fc)
          endif

          niter=niter+1
          if (debug_root.eq.1) then
           print *,"newton: niter,intercept,fc ",niter,intercept,fc
          endif  
         else
          print *,"multi_find_intercept: "
          print *,"arean should not be zero"
          print *,"fastflag=",fastflag
          print *,"continuous_mof=",continuous_mof
          print *,"use_initial_guess=",use_initial_guess
          print *,"niter,vfrac,arean,volcut,fc ",niter,vfrac,arean, &
           volcut,fc
          print *,"volcell ",volcell
          print *,"vfrac_normalize ",vfrac_normalize
          print *,"minphi ",minphi
          print *,"maxphi ",maxphi
          print *,"intercept_default ",intercept_default
          stop
         endif
        enddo ! outer while statement: Newton's method

        if (niter.ge.maxiter) then
         print *,"vof recon failed in multi_find_intercept"
         print *,"niter,maxiter ",niter,maxiter
         print *,"bfact ",bfact
         print *,"vfrac ",vfrac
         do dir=1,sdim
          print *,"dir, dx, xsten0(0) ",dir,dx(dir),xsten0(0,dir)
         enddo
         do dir=1,sdim
          print *,"dir,slope ",dir,slope(dir)
         enddo
         stop
        endif
       endif ! err> moftol
      endif  ! cell has a partial vof

       ! centroid in absolute coordinate system
      do dir=1,sdim
       centroid(dir)=centroid(dir)*maxdx_scale
      enddo
      intercept=intercept*maxdx_scale 

      if (fastflag.eq.0) then
       do nn=1,nlist
       do ii=1,sdim+1
       do dir=1,sdim
        xtetlist(ii,dir,nn)=xtetlist(ii,dir,nn)*maxdx_scale
       enddo 
       enddo 
       enddo 
      else if (fastflag.eq.1) then
       ! do nothing
      else
       print *,"fastflag invalid multi find intercept 3"
       stop
      endif

      return
      end subroutine multi_find_intercept

 
        ! xcell is center of cell, not the cell centroid
        ! refcentroid_scale is passed into this routine.
      subroutine multi_rotatefunc( &
        bfact,dx,xsten0,nhalf0, &
        xtetlist_vof,nlist_vof, &
        xtetlist_cen,nlist_cen, &
        nmax, &
        refcentroid,refvfrac, &
        continuous_mof, &
        angle, &
        levelrz,ff,intercept,testcen, &
        use_initial_guess, &
        fastflag,sdim)

      use geometry_intersect_module
      use global_utility_module

      IMPLICIT NONE

      INTEGER_T continuous_mof,bfact,nhalf0
      INTEGER_T nlist_vof
      INTEGER_T nlist_cen
      INTEGER_T nmax,sdim,fastflag
      REAL_T xtetlist_vof(sdim+1,sdim,nmax)
      REAL_T xtetlist_cen(sdim+1,sdim,nmax)
      REAL_T refcentroid(sdim)
      REAL_T refcentroidT(sdim)
      REAL_T refvfrac
      INTEGER_T use_initial_guess 

      REAL_T areacentroid(sdim) 
      REAL_T testcen(sdim)
      REAL_T testcenT(sdim)
      REAL_T angle(sdim-1)
      REAL_T dx(sdim)
      REAL_T xsten0(-nhalf0:nhalf0,sdim)
      REAL_T xsten2(-1:1,sdim)
      INTEGER_T isten,nhalf2
      INTEGER_T levelrz
      REAL_T ff(sdim) 
      REAL_T volume_cut,facearea
      INTEGER_T dir
      REAL_T nslope(sdim)
      REAL_T intercept
      REAL_T volall

      INTEGER_T i1,j1,k1
      REAL_T volcell_vof
      REAL_T volcell_cen
      REAL_T cencell_vof(sdim)
      REAL_T cencell_cen(sdim)
      REAL_T xtet(sdim+1,sdim)
      INTEGER_T shapeflag
      INTEGER_T ksten_low,ksten_high
      REAL_T volsten,volstencut
      REAL_T areasten
      REAL_T areacentroidsten(sdim)
      REAL_T censten(sdim)
      REAL_T censtencut(sdim)

      nhalf2=1

      if ((sdim.ne.3).and.(sdim.ne.2)) then
       print *,"sdim invalid multi_rotatefunc"
       stop
      endif
      if (nhalf0.lt.3) then
       print *,"nhalf0 invalid"
       stop
      endif

      if ((continuous_mof.ne.0).and. &
          (continuous_mof.ne.1).and. &
          (continuous_mof.ne.2)) then
       print *,"use super_cell invalid"
       stop
      endif
      if ((use_initial_guess.ne.0).and.(use_initial_guess.ne.1)) then
       print *,"use_initial_guess  invalid multirotatefunc"
       stop
      endif

       ! if RZ, cencell can be negative

      if (continuous_mof.eq.0) then
        call Box_volumeFAST( &
         bfact,dx,xsten0,nhalf0, &
         volcell_vof, &
         cencell_vof,levelrz,sdim)
        call Box_volumeFAST( &
         bfact,dx,xsten0,nhalf0, &
         volcell_cen, &
         cencell_cen,levelrz,sdim)
      else if (continuous_mof.eq.1) then
        call Box_volume_super( &
         bfact,dx,xsten0,nhalf0, &
         volcell_vof,cencell_vof, &
         sdim,levelrz)
        call Box_volume_super( &
         bfact,dx,xsten0,nhalf0, &
         volcell_cen,cencell_cen, &
         sdim,levelrz)
      else if (continuous_mof.eq.2) then
        call Box_volumeFAST( &
         bfact,dx,xsten0,nhalf0, &
         volcell_vof, &
         cencell_vof,levelrz,sdim)
        call Box_volume_super( &
         bfact,dx,xsten0,nhalf0, &
         volcell_cen,cencell_cen, &
         sdim,levelrz)
      else
        print *,"continuous_mof invalid"
        stop
      endif

      if (sdim.eq.3) then
       call angle_to_slope3D(angle,nslope,sdim)
       ksten_low=-1
       ksten_high=1
      else if (sdim.eq.2) then
       call angle_to_slope2D(angle,nslope,sdim)
       ksten_low=0
       ksten_high=0
      else
       print *,"sdim invalid"
       stop
      endif

      do dir=1,sdim
       ff(dir)=zero
      enddo

        ! inside of multi_rotatefunc
        ! testcen in absolute coordinate system
        ! testcen is relative to the center cell.
      call multi_find_intercept( &
        bfact,dx,xsten0,nhalf0, &
        nslope,intercept, &
        continuous_mof, &
        xtetlist_vof,nlist_vof, &
        nmax,refvfrac,levelrz,use_initial_guess, &
        testcen,fastflag,sdim)

        ! testcen in absolute coordinate system
      if (fastflag.eq.0) then
         ! xcell used for LS dist.
       if (continuous_mof.eq.0) then
          ! testcen is relative to the center cell.
        call multi_cell_intersection( &
         bfact,dx,xsten0,nhalf0, &
         nslope,intercept, &
         volume_cut,testcen,facearea, &
         levelrz, &
         xtetlist_vof,nlist_vof, &
         nmax,sdim)
       else if ((continuous_mof.eq.1).or.(continuous_mof.eq.2)) then
          ! testcen is relative to the super cell.
        call multi_cell_intersection( &
         bfact,dx,xsten0,nhalf0, &
         nslope,intercept, &
         volume_cut,testcen,facearea, &
         levelrz, &
         xtetlist_cen,nlist_cen, &
         nmax,sdim)
       else
        print *,"continuous_mof invalid"
        stop
       endif
      else if (fastflag.eq.1) then

       shapeflag=0
       if (continuous_mof.eq.0) then
          ! testcen is relative to the center cell.
        call fast_cut_cell_intersection( &
         bfact,dx,xsten0,nhalf0, &
         nslope,intercept, &
         volume_cut,testcen,facearea, &
         areacentroid,levelrz, &
         xsten0,nhalf0,xtet,shapeflag,sdim) 
       else if ((continuous_mof.eq.1).or.(continuous_mof.eq.2)) then
          ! testcen is relative to the super cell.
        volume_cut=zero
        facearea=zero
        do dir=1,sdim
         testcen(dir)=zero
         areacentroid(dir)=zero
        enddo 

        if (sdim.eq.3) then
         ksten_low=-1
         ksten_high=1
        else if (sdim.eq.2) then
         ksten_low=0
         ksten_high=0
        else
         print *,"sdim invalid"
         stop
        endif

        do i1=-1,1
        do j1=-1,1
        do k1=ksten_low,ksten_high

         do isten=-1,1
          xsten2(isten,1)=xsten0(isten+2*i1,1)
          xsten2(isten,2)=xsten0(isten+2*j1,2)
          if (sdim.eq.3) then
           xsten2(isten,sdim)=xsten0(isten+2*k1,sdim)
          endif
         enddo ! isten
         call fast_cut_cell_intersection( &
          bfact,dx,xsten0,nhalf0, &
          nslope,intercept, &
          volsten,censten,areasten, &
          areacentroidsten,levelrz, &
          xsten2,nhalf2,xtet,shapeflag,sdim) 
         volume_cut=volume_cut+volsten
         facearea=facearea+areasten
         do dir=1,sdim
          testcen(dir)=testcen(dir)+volsten*censten(dir)
          areacentroid(dir)=areacentroid(dir)+ &
            areasten*areacentroidsten(dir) 
         enddo
        enddo
        enddo
        enddo  ! i1,j1,k1
        do dir=1,sdim
         if (volume_cut.gt.zero) then
          testcen(dir)=testcen(dir)/volume_cut
         else
          testcen(dir)=zero
         endif
         if (facearea.gt.zero) then
          areacentroid(dir)=areacentroid(dir)/facearea
         else
          areacentroid(dir)=zero
         endif
        enddo ! dir 

       else
        print *,"continuous_mof invalid"
        stop
       endif

      else
       print *,"fastflag invalid multi rotatefunc "
       stop
      endif

      do dir=1,sdim
       testcen(dir)=testcen(dir)-cencell_cen(dir)
      enddo

      call RT_transform_offset(refcentroid,cencell_cen,refcentroidT)
      call RT_transform_offset(testcen,cencell_cen,testcenT)

      do dir=1,sdim
       ff(dir)=(refcentroidT(dir)-testcenT(dir))
      enddo

      return
      end subroutine multi_rotatefunc



      subroutine advance_angle(angle,delangle)
      use global_utility_module
      IMPLICIT NONE
      REAL_T angle,delangle

      if (MOF_PI.eq.zero) then
       MOF_PI=four*atan(one)
      endif

      angle=angle+delangle
      if (angle.lt.-MOF_PI) then
       angle=angle+two*MOF_PI
      endif
      if (angle.gt.MOF_PI) then
       angle=angle-two*MOF_PI
      endif

      return
      end subroutine advance_angle
 
! 3D: x=cos(theta)sin(phi)  y=sin(theta)sin(phi)  z=cos(phi)
! 2D: (cos(theta),sin(theta)) 

      subroutine slope_to_angle(nn,angle,sdim)

      use global_utility_module

      IMPLICIT NONE

      INTEGER_T sdim,dir
      REAL_T nn(sdim)
      REAL_T angle(sdim-1)
      REAL_T mag
      REAL_T x,y,z

      if (MOF_PI.eq.zero) then
       MOF_PI=four*atan(one)
      endif

      mag=zero
      do dir=1,sdim
       mag=mag+nn(dir)**2
      enddo
      mag=sqrt(mag)
      if (mag.le.MLSVOFTOL) then
       print *,"slope_to_angle: invalid slope mag=",mag
       stop
      endif
      do dir=1,sdim
       nn(dir)=nn(dir)/mag
      enddo
      x=nn(1)
      y=nn(2)
      z=nn(sdim)

      if (sdim.eq.3) then

        ! tan(theta)=y/x   sin(theta)/cos(theta)=y/x
       call arctan2(y,x,angle(1))
       if ((y.eq.zero).and.(x.eq.zero)) then
        if (z.eq.zero) then
         print *,"z cannot be zero"
         stop
        else if (z.gt.zero) then
         angle(sdim-1)=zero
        else if (z.lt.zero) then
         angle(sdim-1)=MOF_PI
        else
         print *,"bust slope_to_angle"
         print *,"x,y,z ",x,y,z
         print *,"sdim=",sdim
         stop
        endif
       else if (abs(x).ge.abs(y)) then
        ! tan(phi)=x/(z cos(theta))
        ! sin(phi)/cos(phi)=x/(z cos(theta))
        ! z=cos(phi)  x=sin(phi)cos(theta)  y=x tan(theta)=sin(phi)sin(theta)
        call arctan2(x/cos(angle(1)),z,angle(sdim-1))
       else
        ! tan(phi)=y/(z sin(theta)
        ! sin(phi)/cos(phi)=y/(z sin(theta))
        ! z=cos(phi) y=sin(theta)sin(phi)  x=y/tan(theta)=cos(theta)sin(phi)
        call arctan2(y/sin(angle(1)),z,angle(sdim-1))
       endif 

      else if (sdim.eq.2) then
        ! tan(theta)=y/x   sin(theta)/cos(theta)=y/x
       call arctan2(y,x,angle(1))
 
      else
       print *,"slope_to_angle: sdim invalid"
       stop
      endif

      return
      end subroutine slope_to_angle

        
        ! refcentroid and centroidA relative to cell centroid of the
        ! super cell.
        ! xsten0(0,dir) is center of cell, not the cell centroid
        ! output: intercept,centroidA,nslope

      subroutine find_cut_geom_slope( &
        ls_mof, &
        lsnormal, &
        lsnormal_valid, &
        bfact,dx,xsten0,nhalf0, &
        refcentroid,refvfrac, &
        levelrz,npredict, &
        continuous_mof, &
        nslope,intercept, &
        xtetlist_vof,nlist_vof, &
        xtetlist_cen,nlist_cen, &
        centroidA, &
        nmax, &
        critical_material,fastflag, &
        nmat,nten,sdim)

      use geometry_intersect_module
      use global_utility_module

      IMPLICIT NONE

#include "mofdata.H"

      INTEGER_T continuous_mof,bfact,nhalf0
      INTEGER_T nlist_vof
      INTEGER_T nlist_cen
      INTEGER_T nten,nten_test
      INTEGER_T nmax,sdim,critical_material,nmat,fastflag
      REAL_T xtetlist_vof(sdim+1,sdim,nmax)
      REAL_T xtetlist_cen(sdim+1,sdim,nmax)
      REAL_T xsten0(-nhalf0:nhalf0,sdim)
      REAL_T xsten0_scale(-nhalf0:nhalf0,sdim)
      REAL_T dx(sdim)
      REAL_T refcentroid(sdim)
      REAL_T refvfrac
      INTEGER_T levelrz
      REAL_T npredict(sdim)
      REAL_T intercept
      REAL_T nslope(sdim) 

      REAL_T new_angle(sdim-1)
      REAL_T angle_base(sdim-1)
      REAL_T angle_init(sdim-1)

      REAL_T intercept_init
      REAL_T cen_derive_init(sdim)

      REAL_T intercept_array(MOFITERMAX+1)
      REAL_T cen_array(sdim,MOFITERMAX+1)
      REAL_T angle_array(sdim-1,MOFITERMAX+1)
      REAL_T f_array(sdim,MOFITERMAX+1)  
      REAL_T err_array(MOFITERMAX+1)

      INTEGER_T dir,iter,im_opp,iten
      REAL_T finit(sdim)
      REAL_T fp(sdim)
      REAL_T fm(sdim)
      REAL_T fopt(sdim)
      REAL_T fbase(sdim)

      REAL_T intp,intm,intopt,intbase
      REAL_T cenp(sdim)
      REAL_T cenm(sdim)
      REAL_T cenopt(sdim)
      REAL_T cenbase(sdim)

      REAL_T hh
      REAL_T mag,err
      REAL_T fgrad(sdim,sdim-1)  
      INTEGER_T ii,jj,iicrit
      REAL_T delangle(sdim-1)
      REAL_T RHS(sdim-1)
      REAL_T JTJ(sdim-1,sdim-1)
      REAL_T JTJINV(sdim-1,sdim-1)
      REAL_T tol,local_tol,DET
      REAL_T err_local_min
      REAL_T errinit
      REAL_T angle_initcen(sdim-1)
      REAL_T finitcen(sdim)
      REAL_T centroidA(sdim)
      INTEGER_T use_initial_guess
      REAL_T dx_scale(sdim)
      REAL_T refcentroid_scale(sdim)
      REAL_T maxdx
      INTEGER_T nn
      REAL_T dx_normalize
      INTEGER_T nguess
      REAL_T nLS(sdim)
      REAL_T LSSIGN

      REAL_T ls_mof(D_DECL(-1:1,-1:1,-1:1),nmat+nten)
      REAL_T lsnormal(nmat+nten,sdim)
      INTEGER_T lsnormal_valid(nmat+nten)

      if (nhalf0.lt.1) then
       print *,"nhalf0 invalid"
       stop
      endif
      if (bfact.lt.1) then
       print *,"bfact invalid"
       stop
      endif
      nten_test=( (nmat-1)*(nmat-1)+nmat-1 )/2
      if (nten_test.ne.nten) then
       print *,"nten invalid find cut geom sl nten nten_test ",nten,nten_test
       stop
      endif

      if ((MOFITERMAX.lt.nmat+3).or.(MOFITERMAX.gt.50)) then
       print *,"MOFITERMAX out of range find cut geom slope"
       stop
      endif

      if ((sdim.ne.3).and.(sdim.ne.2)) then
       print *,"sdim invalid find_cut_geom_slope"
       stop
      endif 

      if ((nmat.lt.1).or.(nmat.gt.MAX_NUM_MATERIALS)) then
       print *,"nmat invalid find cut geom slope"
       stop
      endif
      if ((critical_material.lt.1).or.(critical_material.gt.nmat)) then
       print *,"critical_material invalid"
       stop
      endif
      if ((continuous_mof.ne.0).and. &
          (continuous_mof.ne.1).and. &
          (continuous_mof.ne.2)) then
       print *,"continuous_mof invalid"
       stop
      endif

      call scale_MOF_variables( &
        bfact,dx,xsten0,nhalf0, &
        refcentroid,levelrz, &
        dx_scale,xsten0_scale, &
        refcentroid_scale, &
        sdim,maxdx)
 
      if (fastflag.eq.0) then

       do nn=1,nlist_vof
       do ii=1,sdim+1
       do dir=1,sdim
        xtetlist_vof(ii,dir,nn)=xtetlist_vof(ii,dir,nn)/maxdx
       enddo 
       enddo 
       enddo 

       if (continuous_mof.eq.0) then
        ! do nothing
       else if ((continuous_mof.eq.1).or.(continuous_mof.eq.2)) then
        do nn=1,nlist_cen
        do ii=1,sdim+1
        do dir=1,sdim
         xtetlist_cen(ii,dir,nn)=xtetlist_cen(ii,dir,nn)/maxdx
        enddo 
        enddo 
        enddo 
       else
        print *,"continuous_mof invalid"
        stop
       endif

      else if (fastflag.eq.1) then
       ! do nothing
      else
       print *,"fastflag invalid find cut geom slope"
       stop
      endif

      dx_normalize=dx_scale(1)
      if (dx_normalize.gt.one) then
       dx_normalize=one
      endif
 
      tol=dx_scale(1)*GAUSSNEWTONTOL
      local_tol=dx_scale(1)*tol*1.0E-2

        ! -pi < angle < pi
      call slope_to_angle(npredict,angle_init,sdim)
      nguess=1 
      do dir=1,sdim-1
       angle_array(dir,nguess)=angle_init(dir)
      enddo

      if (MOF_TURN_OFF_LS.eq.0) then

       if (lsnormal_valid(critical_material).eq.1) then
        do dir=1,sdim
         nLS(dir)=lsnormal(critical_material,dir)
        enddo
        call slope_to_angle(nLS,angle_init,sdim)
        nguess=nguess+1 
        do dir=1,sdim-1
         angle_array(dir,nguess)=angle_init(dir)
        enddo
       else if (lsnormal_valid(critical_material).eq.0) then
        ! do nothing
       else
        print *,"LSNORMAL_valid invalid1"
        print *,"critical_material,flag ",critical_material, &
         lsnormal_valid(critical_material) 
        stop
       endif
       do im_opp=1,nmat
        if (im_opp.ne.critical_material) then
         call get_iten(critical_material,im_opp,iten,nmat)
         if (lsnormal_valid(nmat+iten).eq.1) then
          LSSIGN=one
          if (im_opp.lt.critical_material) then
           LSSIGN=-one
          endif
          do dir=1,sdim
           nLS(dir)=lsnormal(nmat+iten,dir)
          enddo
          call slope_to_angle(nLS,angle_init,sdim)
          nguess=nguess+1 
          do dir=1,sdim-1
           angle_array(dir,nguess)=angle_init(dir)
          enddo
         else if (lsnormal_valid(nmat+iten).eq.0) then
          ! do nothing
         else
          print *,"LSNORMAL_valid invalid2"
          print *,"critical_material,nmat,iten,flag ", &
           critical_material,nmat,iten, &
           lsnormal_valid(nmat+iten) 
          stop
         endif
        endif
       enddo ! im_opp
 
      else if (MOF_TURN_OFF_LS.eq.1) then
       ! do nothing
      else
       print *,"MOF_TURN_OFF_LS invalid"
       stop
      endif
         
      iicrit=0
      do iter=1,nguess 
       do dir=1,sdim-1
        angle_init(dir)=angle_array(dir,iter)
       enddo

        ! find finit=xref-xact for cut domain cut by a line.
       use_initial_guess=0
       intercept_init=zero
       call multi_rotatefunc( &
        bfact,dx_scale,xsten0_scale,nhalf0, &
        xtetlist_vof,nlist_vof, &
        xtetlist_cen,nlist_cen, &
        nmax, &
        refcentroid_scale,refvfrac, &
        continuous_mof, &
        angle_init, &
        levelrz, &
        finit,intercept_init,cen_derive_init, &
        use_initial_guess,fastflag,sdim)

       errinit=zero
       do dir=1,sdim
        errinit=errinit+finit(dir)**2
       enddo
       errinit=sqrt(errinit)
       err=errinit

       do dir=1,sdim
        f_array(dir,iter)=finit(dir)
       enddo
       err_array(iter)=err
       intercept_array(iter)=intercept_init
       do dir=1,sdim
        cen_array(dir,iter)=cen_derive_init(dir)
       enddo 
       if (iicrit.eq.0) then
        iicrit=iter
       else if (err.lt.err_array(iicrit)) then
        iicrit=iter
       endif
      enddo ! iter=1..nguess

      if ((iicrit.lt.1).or.(iicrit.gt.nguess)) then
       print *,"iicrit invalid"
       stop
      endif
      do dir=1,sdim-1
       angle_array(dir,1)=angle_array(dir,iicrit)
      enddo
      do dir=1,sdim
       f_array(dir,1)=f_array(dir,iicrit)
      enddo
      err_array(1)=err_array(iicrit)
      intercept_array(1)=intercept_array(iicrit)
      do dir=1,sdim
       cen_array(dir,1)=cen_array(dir,iicrit)
      enddo 
 
      err_local_min=err_array(1)
      err=err_array(1)

      iter=0

      do while ((iter.lt.MOFITERMAX).and.(err.gt.tol).and. &
                (err_local_min.gt.local_tol))

        hh=MOFHH  ! 1 degree=pi/180
 
        do dir=1,sdim
         fbase(dir)=f_array(dir,iter+1)
        enddo

        do ii=1,sdim-1
         do dir=1,sdim-1
          angle_base(dir)=angle_array(dir,iter+1)
         enddo
 
         angle_base(ii)=angle_base(ii)+hh
         intp=intercept_array(iter+1)
         intm=intercept_array(iter+1)
         use_initial_guess=1
          ! fp=xref-cenp
         call multi_rotatefunc( &
          bfact,dx_scale,xsten0_scale,nhalf0, &
          xtetlist_vof,nlist_vof, &
          xtetlist_cen,nlist_cen, &
          nmax, &
          refcentroid_scale,refvfrac, &
          continuous_mof, &
          angle_base,levelrz, &
          fp,intp,cenp, &
          use_initial_guess, &
          fastflag,sdim)

         angle_base(ii)=angle_base(ii)-two*hh
          ! fm=xref-cenm
         call multi_rotatefunc( &
          bfact,dx_scale,xsten0_scale,nhalf0, &
          xtetlist_vof,nlist_vof, &
          xtetlist_cen,nlist_cen, &
          nmax, &
          refcentroid_scale,refvfrac, &
          continuous_mof, &
          angle_base,levelrz, &
          fm,intm,cenm, &
          use_initial_guess, &
          fastflag,sdim)
     
! jacobian matrix has:
!   f1_1  f1_2
!   f2_1  f2_2
!   f3_1  f3_2  ....

           ! fgrad ~ df/dtheta  (has dimensions of length)
         do dir=1,sdim
          fgrad(dir,ii)=(fp(dir)-fm(dir))/(two*hh)
         enddo

        enddo  ! ii

        do ii=1,sdim-1
         do jj=1,sdim-1
          JTJ(ii,jj)=zero
          do dir=1,sdim
           JTJ(ii,jj)=JTJ(ii,jj)+fgrad(dir,ii)*fgrad(dir,jj)
          enddo
         enddo
        enddo
        if (sdim.eq.3) then
         DET=JTJ(1,1)*JTJ(2,2)-JTJ(1,2)*JTJ(2,1)
        else if (sdim.eq.2) then
         DET=JTJ(1,1)
        else
         print *,"sdim invalid"
         stop
        endif 

         ! DET has dimensions of length squared
        if (abs(DET).ge.CENTOL*(dx_normalize**2)) then 
         if (sdim.eq.3) then
          JTJINV(1,1)=JTJ(2,2)
          JTJINV(2,2)=JTJ(1,1)
          JTJINV(1,2)=-JTJ(1,2)
          JTJINV(2,1)=-JTJ(2,1)
         else if (sdim.eq.2) then
          JTJINV(1,1)=one
         else
          print *,"sdim invalid"
          stop
         endif

         do ii=1,sdim-1
          do jj=1,sdim-1
           JTJINV(ii,jj)=JTJINV(ii,jj)/DET
          enddo
         enddo

         do ii=1,sdim-1
          angle_base(ii)=angle_array(ii,iter+1)
         enddo
         do ii=1,sdim-1  ! compute -JT * r
          RHS(ii)=zero
          do dir=1,sdim
           RHS(ii)=RHS(ii)-fgrad(dir,ii)*fbase(dir)
          enddo
         enddo

         err_local_min=zero
         do ii=1,sdim-1  
          err_local_min=err_local_min+RHS(ii)**2
         enddo
         err_local_min=sqrt(err_local_min)
 
         do ii=1,sdim-1  ! compute JTJ^-1 (RHS)
          delangle(ii)=zero
          do jj=1,sdim-1
           delangle(ii)=delangle(ii)+JTJINV(ii,jj)*RHS(jj)
          enddo
          if (delangle(ii).gt.MOF_ANGLE_MAX) then
           delangle(ii)=MOF_ANGLE_MAX
          else if (delangle(ii).lt.-MOF_ANGLE_MAX) then
           delangle(ii)=-MOF_ANGLE_MAX
          endif
         enddo
          ! -pi<angle<pi 
         do ii=1,sdim-1
          call advance_angle(angle_base(ii),delangle(ii))
         enddo

         intopt=intercept_array(iter+1)
         use_initial_guess=1
         call multi_rotatefunc( &
          bfact,dx_scale,xsten0_scale,nhalf0, &
          xtetlist_vof,nlist_vof, &
          xtetlist_cen,nlist_cen, &
          nmax, &
          refcentroid_scale,refvfrac, &
          continuous_mof, &
          angle_base,levelrz, &
          fopt,intopt,cenopt, &
          use_initial_guess, &
          fastflag,sdim)
        else
         err_local_min=zero

         do dir=1,sdim
          fopt(dir)=fbase(dir)
         enddo
         do dir=1,sdim-1
          angle_base(dir)=angle_array(dir,iter+1)
         enddo
         intopt=intercept_array(iter+1)
         do dir=1,sdim
          cenopt(dir)=cen_array(dir,iter+1)
         enddo
        endif  ! JTJ<=0
    
        err=zero 
        do dir=1,sdim
         f_array(dir,iter+2)=fopt(dir)
         err=err+f_array(dir,iter+2)**2
        enddo
        err=sqrt(err)

        do ii=1,sdim-1
         angle_array(ii,iter+2)=angle_base(ii) 
        enddo
        err_array(iter+2)=err

        intercept_array(iter+2)=intopt
        do dir=1,sdim
         cen_array(dir,iter+2)=cenopt(dir)
        enddo

        iter=iter+1
      enddo ! while error>tol

      iicrit=iter
      do ii=0,iter

       if ((MOF_DEBUG_RECON.eq.1).or. &
           ((MOF_DEBUG_RECON.eq.2).and.(fastflag.eq.0))) then
        print *,"iimof,eemof ",ii,err_array(ii+1)
       endif

       if (err_array(ii+1).le.err_array(iicrit+1)) then
        iicrit=ii
       endif
      enddo ! ii

      if ((MOF_DEBUG_RECON.eq.1).or. &
          ((MOF_DEBUG_RECON.eq.2).and.(fastflag.eq.0))) then
       MOF_DEBUG_RECON_COUNT=MOF_DEBUG_RECON_COUNT+1
       print *,"MOF_DEBUG_RECON_COUNT=",MOF_DEBUG_RECON_COUNT
       print *,"AFTER------------------------------- "
      endif

      do dir=1,sdim-1 
       new_angle(dir)=angle_array(dir,iicrit+1)
      enddo
      intercept=intercept_array(iicrit+1)

      if (sdim.eq.3) then
       call angle_to_slope3D(new_angle,nslope,sdim)
      else if (sdim.eq.2) then
       call angle_to_slope2D(new_angle,nslope,sdim)
      else
       print *,"sdim invalid"
       stop
      endif

      do dir=1,sdim
       centroidA(dir)=cen_array(dir,iicrit+1)
      enddo

      intercept=intercept*maxdx
      do dir=1,sdim
       centroidA(dir)=centroidA(dir)*maxdx
      enddo

      if (fastflag.eq.0) then

       do nn=1,nlist_vof
       do ii=1,sdim+1
       do dir=1,sdim
        xtetlist_vof(ii,dir,nn)=xtetlist_vof(ii,dir,nn)*maxdx
       enddo 
       enddo 
       enddo 

       if (continuous_mof.eq.0) then
        ! do nothing
       else if ((continuous_mof.eq.1).or.(continuous_mof.eq.2)) then
        do nn=1,nlist_cen
        do ii=1,sdim+1
        do dir=1,sdim
         xtetlist_cen(ii,dir,nn)=xtetlist_cen(ii,dir,nn)*maxdx
        enddo
        enddo
        enddo 
       else
        print *,"continuous_mof invalid"
        stop
       endif

      else if (fastflag.eq.1) then
       ! do nothing
      else
       print *,"fastflag invalid find_cut_geom_slope 2"
       stop
      endif


      mof_iterations(critical_material)= &
         mof_iterations(critical_material)+iter
      mof_calls(critical_material)= &
         mof_calls(critical_material)+1

      return
      end subroutine find_cut_geom_slope



#if (STANDALONE==0)

       ! cen in absolute coordinates
       ! returns a volume fraction
      subroutine find_override_vfrac_coarse( &
       bfact,dx,xsten,nhalf, &
       time,vfrac,cen, &
       volbox,cenbox,sdim)

      use global_distance_module
      use geometry_intersect_module

      IMPLICIT NONE

      INTEGER_T sdim,bfact,nhalf
      REAL_T time
      REAL_T xsten(-nhalf:nhalf,sdim)
      REAL_T dx(sdim)
      REAL_T vfrac
      REAL_T cen(sdim)
      INTEGER_T dir
      REAL_T lnode(4*(sdim-1))
      INTEGER_T inode,jnode,knode,klo,khi,isynth
      REAL_T xn,yn,zn,facearea
      REAL_T areacentroid(sdim)
      REAL_T volbox
      REAL_T cenbox(sdim)

#include "mofdata.H"

      if (bfact.lt.1) then
       print *,"bfact invalid"
       stop
      endif
      if ((sdim.ne.2).and.(sdim.ne.3)) then
       print *,"sdim invalid"
       stop
      endif
      if (nhalf.lt.1) then
       print *,"nhalf invalid find override vfrac coarse"
       stop
      endif

      isynth=1
      if (sdim.eq.2) then
        klo=0
        khi=0
      else if (sdim.eq.3) then
        klo=-1
        khi=1
      else
        print *,"sdim invalid"
        stop
      endif

      do knode=klo,khi,2
      do jnode=-1,1,2
      do inode=-1,1,2
        xn=xsten(inode,1)
        yn=xsten(jnode,2)
        if (sdim.eq.2) then
         zn=yn
        else if (sdim.eq.3) then
         zn=xsten(knode,sdim)
        else
         print *,"dimension bust calling distance_function_override"
         stop
        endif
        lnode(isynth)=distance_function_override(xn,yn,zn,time)
        isynth=isynth+1
      enddo
      enddo
      enddo

        ! returns centroid in absolute coordinate system
      call fast_cell_intersection_grid( &
        bfact,dx,xsten,nhalf, &
        lnode, &
        vfrac,cen,facearea,areacentroid, &
        levelrz,volbox,cenbox,sdim)
      if (volbox.le.zero) then
        vfrac=zero
      else
        vfrac=vfrac/volbox
      endif

      if (vfrac.le.VOF_CUTOFF) then
        vfrac=zero
        do dir=1,sdim
         cen(dir)=cenbox(dir)
        enddo
      endif
      if (vfrac.ge.one-VOF_CUTOFF) then
        vfrac=one
        do dir=1,sdim
         cen(dir)=cenbox(dir)
        enddo
      endif

      return
      end subroutine find_override_vfrac_coarse

       ! cen in absolute coordinates
      subroutine find_override_vfrac( &
       bfact,dx,xsten,nhalf, &
       nrefine,ignore_flag, &
       time,vfrac,cen,sdim)

      IMPLICIT NONE

      INTEGER_T ignore_flag,nrefine,bfact,nhalf
      INTEGER_T sdim
      REAL_T time
      REAL_T xsten(-nhalf:nhalf,sdim)
      REAL_T xsten2(-1:1,sdim)
      INTEGER_T nhalf2
      REAL_T dx(sdim)
      REAL_T vfrac
      REAL_T cen(sdim)
      INTEGER_T dir
      INTEGER_T inode,jnode,knode,khi,nside,iside
      REAL_T dxrefine(sdim)
      REAL_T vfrac_node
      REAL_T cen_node(sdim)
      REAL_T volbox,volbox_node
      REAL_T cenbox(sdim)
      REAL_T cenbox_node(sdim)
 

#include "mofdata.H"

      nhalf2=1

      if (bfact.lt.1) then
       print *,"bfact invalid"
       stop
      endif
      if (nhalf.lt.1) then
       print *,"nhalf invalid find override vfrac"
       stop
      endif
      if ((sdim.ne.2).and.(sdim.ne.3)) then
       print *,"sdim invalid"
       stop
      endif
      if ((ignore_flag.ne.0).and.(ignore_flag.ne.1)) then
       print *,"ignore_flag invalid"
       stop
      endif
      if ((nrefine.lt.0).or.(nrefine.gt.2)) then
       print *,"nrefine out of range"
       stop
      endif
      if ((imaterial_override.eq.0).and. &
          (ignore_flag.eq.0)) then

       do dir=1,sdim
        cen(dir)=xsten(0,dir)
       enddo
       vfrac=zero

      else if ((imaterial_override.gt.0).or. &
               (ignore_flag.eq.1)) then

       nside=1
       do iside=1,nrefine
        nside=nside*2
       enddo
       if (sdim.eq.2) then
        khi=1
       else if (sdim.eq.3) then
        khi=nside
       else
        print *,"sdim invalid"
        stop
       endif

       do dir=1,sdim
        cen(dir)=zero
        cenbox(dir)=zero
       enddo
       vfrac=zero
       volbox=zero

       do dir=1,sdim
        dxrefine(dir)=(xsten(1,dir)-xsten(-1,dir))/nside
       enddo

       do knode=1,khi
       do jnode=1,nside
       do inode=1,nside
        dir=1
        xsten2(-1,dir)=xsten(-1,dir)+(inode-1)*dxrefine(dir)
        dir=2
        xsten2(-1,dir)=xsten(-1,dir)+(jnode-1)*dxrefine(dir)
        if (sdim.eq.3) then
         dir=sdim
         xsten2(-1,dir)=xsten(-1,dir)+(knode-1)*dxrefine(dir)
        endif
        do dir=1,sdim
         xsten2(1,dir)=xsten2(-1,dir)+dxrefine(dir)
         xsten2(0,dir)=(xsten2(-1,dir)+xsten2(1,dir))/two
        enddo

        call find_override_vfrac_coarse( &
         bfact,dx,xsten2,nhalf2, &
         time,vfrac_node,cen_node, &
         volbox_node,cenbox_node,sdim)

        vfrac=vfrac+vfrac_node*volbox_node
        volbox=volbox+volbox_node
        do dir=1,sdim
         cenbox(dir)=cenbox(dir)+cenbox_node(dir)*volbox_node
         cen(dir)=cen(dir)+cen_node(dir)*vfrac_node*volbox_node
        enddo 
       enddo
       enddo
       enddo

       if (volbox.le.zero) then
        vfrac=zero
        do dir=1,sdim
         cenbox(dir)=xsten(0,dir)
        enddo
       else
        vfrac=vfrac/volbox
        do dir=1,sdim
         cenbox(dir)=cenbox(dir)/volbox
         if (vfrac.ge.VOF_CUTOFF) then
          cen(dir)=cen(dir)/(vfrac*volbox)
         endif
        enddo
       endif

       if (vfrac.le.VOF_CUTOFF) then
        vfrac=zero
        do dir=1,sdim
         cen(dir)=cenbox(dir)
        enddo
       endif
       if (vfrac.ge.one-VOF_CUTOFF) then
        vfrac=one
        do dir=1,sdim
         cen(dir)=cenbox(dir)
        enddo
       endif

      else
       print *,"imaterial override invalid find_override_vfrac"
       stop
      endif

      return
      end subroutine find_override_vfrac


         ! nslope points into the solid
         ! distance_function_override>0 in solid
      subroutine find_override_slope( &
       bfact,dx,xsten0,nhalf0, &
       nslope,time,sdim)

      use global_distance_module

      IMPLICIT NONE

      INTEGER_T sdim,bfact,nhalf0
      REAL_T xsten0(-nhalf0:nhalf0,sdim)
      REAL_T dx(sdim)
      REAL_T nslope(sdim)
      REAL_T LS(D_DECL(-1:1,-1:1,-1:1))
      REAL_T xn,yn,zn,time,mag
      INTEGER_T dir,inode,jnode,knode,knodelo,knodehi

#include "mofdata.H"

      if (nhalf0.lt.1) then
       print *,"nhalf0 invalid"
       stop
      endif
      if (bfact.lt.1) then
       print *,"bfact invalid"
       stop
      endif
      if ((sdim.ne.2).and.(sdim.ne.3)) then
       print *,"sdim invalid"
       stop
      endif

      if (imaterial_override.eq.0) then

       do dir=1,sdim
        nslope(dir)=zero
       enddo
       nslope(sdim)=one

      else if (imaterial_override.gt.0) then

       do dir=1,sdim
        nslope(dir)=zero
       enddo

       if (sdim.eq.2) then
        knodelo=0
        knodehi=0
       else if (sdim.eq.3) then
        knodelo=-1
        knodehi=1
       else
        print *,"dimension bust"
        stop
       endif

       do inode=-1,1,2
       do jnode=-1,1,2
       do knode=knodelo,knodehi,2
        xn=xsten0(inode,1)
        yn=xsten0(jnode,2)
        if (sdim.eq.2) then
         zn=yn
        else if (sdim.eq.3) then
         zn=xsten0(knode,sdim)
        else
         print *,"dimension bust prior to distance_function_override"
         stop
        endif
        LS(D_DECL(inode,jnode,knode))= &
          distance_function_override(xn,yn,zn,time)

        dir=1
        nslope(dir)=nslope(dir)+inode*LS(D_DECL(inode,jnode,knode))
        dir=2
        nslope(dir)=nslope(dir)+jnode*LS(D_DECL(inode,jnode,knode))

        if (sdim.eq.3) then
         dir=sdim
         nslope(dir)=nslope(dir)+knode*LS(D_DECL(inode,jnode,knode))
        else if (sdim.eq.2) then
         ! do nothing
        else
         print *,"sdim invalid"
         stop
        endif

       enddo
       enddo
       enddo

       mag=zero
       do dir=1,sdim
        nslope(dir)=nslope(dir)/ &
          (xsten0(1,dir)-xsten0(-1,dir))
        mag=mag+nslope(dir)**2
       enddo
       mag=sqrt(mag)
       if (mag.eq.zero) then
        do dir=1,sdim
         nslope(dir)=zero
        enddo
        nslope(sdim)=one
       else if (mag.gt.zero) then
        do dir=1,sdim
         nslope(dir)=nslope(dir)/mag
        enddo
       else
        print *,"mag invalid"
        stop
       endif

      else
       print *,"imaterial override invalid"
       stop
      endif

      return
      end subroutine find_override_slope

#endif


      subroutine find_predict_slope(slope,mag,cen_free,cen_ref, &
       bfact,dx,xsten,nhalf,levelrz,sdim)
      use global_utility_module

      IMPLICIT NONE

      INTEGER_T bfact,sdim,levelrz,nhalf,dir
      REAL_T xsten(-nhalf:nhalf,sdim)
      REAL_T dx(sdim)
      REAL_T slope(sdim)
      REAL_T slopeRT(sdim)
      REAL_T cen_free(sdim)
      REAL_T cen_ref(sdim)
      REAL_T cen_freeXYZ(sdim)
      REAL_T cen_refXYZ(sdim)
      REAL_T mag,RR,theta,mag_temp

      if ((sdim.ne.2).and.(sdim.ne.3)) then
       print *,"sdim invalid"
       stop
      endif
      if (bfact.lt.1) then
       print *,"bfact invalid"
       stop
      endif
      if (nhalf.lt.1) then
       print *,"nhalf invalid"
       stop
      endif
      if ((levelrz.eq.1).and.(sdim.ne.2)) then
       print *,"dimension bust"
       stop
      endif

      if ((levelrz.eq.0).or.(levelrz.eq.1)) then
       RR=one
       do dir=1,sdim
        slope(dir)=cen_ref(dir)-cen_free(dir)
       enddo
       call prepare_normal(slope,RR,mag)
      else if (levelrz.eq.3) then
       call RT_transform(cen_ref,cen_refXYZ) 
       call RT_transform(cen_free,cen_freeXYZ) 
       do dir=1,sdim
        slope(dir)=cen_refXYZ(dir)-cen_freeXYZ(dir)
       enddo
       RR=one
       call prepare_normal(slope,RR,mag)
       theta=xsten(0,2)
       slopeRT(1)=cos(theta)*slope(1)+sin(theta)*slope(2)
       slopeRT(2)=-sin(theta)*slope(1)+cos(theta)*slope(2)
       if (sdim.eq.3) then
        slopeRT(sdim)=slope(sdim)
       endif
       if (xsten(0,1).gt.zero) then
        RR=one/xsten(0,1)
        call prepare_normal(slopeRT,RR,mag_temp)
        do dir=1,sdim
         slope(dir)=slopeRT(dir)
        enddo
       else
        print *,"xsten(0,1) must be positive"
        stop
       endif
      else
       print *,"find_predict_slope: levelrz invalid"
       stop
      endif

      return
      end subroutine find_predict_slope

! this routine is called at most nmat times; after uncaptured_volume=0,
! routine is not called anymore.
! First time routine is called, the flag for all materials is 0
! reference and actual centroid relative to cell centroid.
! xcell is cell center (not cell centroid)

      subroutine individual_MOF( &
        ls_mof, &
        lsnormal, &
        lsnormal_valid, &
        bfact,dx,xsten0,nhalf0, &
        order_algorithm_in, &
        ngeom_recon, &
        xtetlist_vof, &
        xtetlist_cen, &
        nmax, &
        mofdata,levelrz, &
        imaterial_count, &
        uncaptured_volume_vof, &
        uncaptured_volume_cen, &
        multi_centroidA, &
        time, &
        continuous_mof, &
        nmat,nten,sdim)

      use geometry_intersect_module

      IMPLICIT NONE
   
      INTEGER_T sdim 
      INTEGER_T, INTENT (IN) :: nhalf0
      INTEGER_T, INTENT (IN) :: bfact
      INTEGER_T, INTENT (IN) :: nten
      REAL_T, INTENT (IN), DIMENSION(sdim) :: dx
      REAL_T, INTENT (IN), DIMENSION(-nhalf0:nhalf0,sdim) :: xsten0
      INTEGER_T ngeom_recon
      INTEGER_T continuous_mof
      INTEGER_T nmat,dir,nmax
      INTEGER_T order_algorithm_in(nmat)
      REAL_T uncaptured_volume_vof 
      REAL_T uncaptured_volume_cen
      REAL_T available_volume 
      REAL_T mofdata(nmat*ngeom_recon)
      REAL_T centroidA(sdim)
      REAL_T multi_centroidA(nmat,sdim)
      REAL_T time
      INTEGER_T levelrz,im,vofcomp
      INTEGER_T imaterial_count
      REAL_T distmax
      INTEGER_T ordermax,order_min
      INTEGER_T critical_material,override_selected
      REAL_T volcut_vof,volcell_vof
      REAL_T volcut_cen,volcell_cen
      REAL_T areacut,mag
      REAL_T cencut_vof(sdim)
      REAL_T cencell_vof(sdim)
      REAL_T cencut_cen(sdim)
      REAL_T cencell_cen(sdim)
      REAL_T xtetlist_vof(sdim+1,sdim,nmax)
      REAL_T xtetlist_cen(sdim+1,sdim,nmax)
      INTEGER_T nlist_vof
      INTEGER_T nlist_cen
      INTEGER_T single_material_takes_all
      INTEGER_T single_material_im
      INTEGER_T test_order
      REAL_T test_vfrac,max_vfrac
      INTEGER_T use_initial_guess
      REAL_T npredict(sdim)
      REAL_T refcentroid(sdim)
      REAL_T centroid_ref(sdim)
      REAL_T centroid_free(sdim)
      REAL_T refvfrac
      REAL_T single_volume
      REAL_T nslope(sdim),intercept
      INTEGER_T mat_before,vofcomp_before
      INTEGER_T fastflag,use_super_cell
      INTEGER_T nten_test
      REAL_T vofmain(nmat)

      REAL_T ls_mof(D_DECL(-1:1,-1:1,-1:1),nmat+nten)
      REAL_T lsnormal(nmat+nten,sdim)
      INTEGER_T lsnormal_valid(nmat+nten)

#include "mofdata.H"

      fastflag=1

      if (nhalf0.lt.3) then
       print *,"nhalf0 invalid"
       stop
      endif
      nten_test=( (nmat-1)*(nmat-1)+nmat-1 )/2
      if (nten_test.ne.nten) then
       print *,"nten invalid indmof nten nten_test ",nten,nten_test
       stop
      endif
      if (uncaptured_volume_vof.le.zero) then
       print *,"uncaptured_volume_vof invalid"
       stop
      endif
      if (ngeom_recon.ne.2*sdim+3) then
       print *,"ngeom_recon invalid"
       stop
      endif
      if ((sdim.ne.3).and.(sdim.ne.2)) then
       print *,"sdim invalid individual_MOF"
       stop
      endif
      if ((nmat.lt.1).or.(nmat.gt.MAX_NUM_MATERIALS)) then
       print *,"nmat invalid individual mof"
       stop
      endif
      if ((imaterial_count.lt.1).or.(imaterial_count.gt.nmat)) then
       print *,"imaterial_count invalid"
       stop
      endif
      if (nmax.lt.10) then
       print *,"nmax too small"
       stop
      endif
      if ((continuous_mof.ne.0).and. &
          (continuous_mof.ne.1).and. &
          (continuous_mof.ne.2)) then
       print *,"continuous_mof invalid"
       stop
      endif

       ! cencell_vof,cencell_cen is in absolute coordinate system

      if (continuous_mof.eq.0) then
       call Box_volumeFAST(bfact,dx,xsten0,nhalf0,volcell_vof, &
        cencell_vof,levelrz,sdim)
       call Box_volumeFAST(bfact,dx,xsten0,nhalf0,volcell_cen, &
        cencell_cen,levelrz,sdim)
      else if (continuous_mof.eq.1) then
       call Box_volume_super(bfact,dx,xsten0,nhalf0, &
         volcell_vof,cencell_vof, &
         sdim,levelrz)
       call Box_volume_super(bfact,dx,xsten0,nhalf0, &
         volcell_cen,cencell_cen, &
         sdim,levelrz)
      else if (continuous_mof.eq.2) then
       call Box_volumeFAST(bfact,dx,xsten0,nhalf0,volcell_vof, &
         cencell_vof,levelrz,sdim)
       call Box_volume_super(bfact,dx,xsten0,nhalf0, &
         volcell_cen,cencell_cen, &
         sdim,levelrz)
      else
        print *,"continuous_mof invalid"
        stop
      endif

      if ((imaterial_count.gt.1).and. &
          (imaterial_count.le.nmat)) then
       fastflag=0
      else if (imaterial_count.eq.1) then
       fastflag=1
      else 
       print *,"imaterial_count invalid"
       stop
      endif

       ! cencut is in absolute coordinate system

      if (fastflag.eq.0) then

         ! get triangulation of uncaptured space in the cell
         ! in: individual_MOF
       if (continuous_mof.eq.0) then
        use_super_cell=0
        call tets_box_planes_super(bfact,dx,xsten0,nhalf0, &
         mofdata, &
         xtetlist_vof,nlist_vof,nmax,nmat,use_super_cell,sdim)
        use_super_cell=0
        call tets_box_planes_super(bfact,dx,xsten0,nhalf0, &
         mofdata, &
         xtetlist_cen,nlist_cen,nmax,nmat,use_super_cell,sdim)
       else if (continuous_mof.eq.1) then
        use_super_cell=1
        call tets_box_planes_super(bfact,dx,xsten0,nhalf0, &
         mofdata, &
         xtetlist_vof,nlist_vof,nmax,nmat,use_super_cell,sdim)
        use_super_cell=1
        call tets_box_planes_super(bfact,dx,xsten0,nhalf0, &
         mofdata, &
         xtetlist_cen,nlist_cen,nmax,nmat,use_super_cell,sdim)
       else if (continuous_mof.eq.2) then
        use_super_cell=0
        call tets_box_planes_super(bfact,dx,xsten0,nhalf0, &
         mofdata, &
         xtetlist_vof,nlist_vof,nmax,nmat,use_super_cell,sdim)
        use_super_cell=1
        call tets_box_planes_super(bfact,dx,xsten0,nhalf0, &
         mofdata, &
         xtetlist_cen,nlist_cen,nmax,nmat,use_super_cell,sdim)
       else
        print *,"continuous_mof invalid"
        stop
       endif

       call get_cut_geom3D(xtetlist_vof,nlist_vof,nmax,levelrz, &
         volcut_vof,cencut_vof,sdim)
       call get_cut_geom3D(xtetlist_cen,nlist_cen,nmax,levelrz, &
         volcut_cen,cencut_cen,sdim)

      else if (fastflag.eq.1) then
       volcut_vof=volcell_vof
       volcut_cen=volcell_cen
       do dir=1,sdim
        cencut_vof(dir)=cencell_vof(dir)  
        cencut_cen(dir)=cencell_cen(dir)  
       enddo
      else
       print *,"fastflag invalid individual MOF"
       stop
      endif      

      available_volume=uncaptured_volume_vof

      if (abs(volcut_vof-uncaptured_volume_vof).gt.VOF_CUTOFF*volcell_vof) then
        print *,"volcut_vof invalid individual mof"
        print *,"volcut_vof ",volcut_vof
        print *,"uncaptured_volume_vof ",uncaptured_volume_vof
        stop
      endif

        ! figure out the next material to fill the unoccupied region.
      distmax=-one
      order_min=9999
      ordermax=0

      do im=1,nmat
       vofcomp=(im-1)*ngeom_recon+1
       if (NINT(mofdata(vofcomp+sdim+1)).gt.ordermax) then
        ordermax=NINT(mofdata(vofcomp+sdim+1))
       endif
       vofmain(im)=mofdata(vofcomp)
      enddo

      if (ordermax.ge.nmat) then
       print *,"all the materials already initialized"
       stop
      endif
      if (ordermax.lt.0) then
       print *,"ordermax invalid"
       stop
      endif

      critical_material=0
      single_material_takes_all=0
      single_material_im=0
      do im=1,nmat
       vofcomp=(im-1)*ngeom_recon+1
       if (mofdata(vofcomp+sdim+1).eq.zero) then
        single_volume=mofdata(vofcomp)*volcell_vof
        if (single_volume.ge.(one-VOF_CUTOFF)*available_volume) then
         if ((single_material_takes_all.ne.0).or. &
             (single_material_im.ne.0)) then
          print *,"cannot have two materials at once"
          print *,"single_material_takes_all ",single_material_takes_all
          print *,"im ",im
          print *,"single vol ",single_volume
          print *,"uncapt vol ",uncaptured_volume_vof
          print *,"uncapt volfrac ",uncaptured_volume_vof/volcell_vof
          stop
         endif
         single_material_takes_all=1
         single_material_im=im
         critical_material=im
         distmax=-one  ! tells code below not to search for a slope
        endif
       endif  ! material not already processed.
      enddo ! im

      if (single_material_takes_all.eq.1) then
       uncaptured_volume_vof=zero
       available_volume=zero
      else if (single_material_takes_all.eq.0) then
       ! do nothing
      else 
       print *,"bust individual_MOF"
       print *,"sdim,nmat ",sdim,nmat
       stop
      endif

       
         ! if uncaptured_volume_vof=0, then there is no need to find the
         ! slope since "single_material_takes_all=1". 
      if (uncaptured_volume_vof.gt.zero) then

        ! find unprocessed material whose moment is furthest from cencut. 

       override_selected=0

       do im=1,nmat
        vofcomp=(im-1)*ngeom_recon+1

        single_volume=mofdata(vofcomp)*volcell_vof

         ! only find the slope if refvfrac<>0 and refvfrac<available vol.
        if ((mofdata(vofcomp+sdim+1).eq.zero).and. &
            (single_volume.ge.VOF_CUTOFF*volcell_vof).and. &
            (single_volume.le.(one-VOF_CUTOFF)*available_volume)) then

         do dir=1,sdim
          centroid_free(dir)=cencut_cen(dir)
          centroid_ref(dir)=mofdata(vofcomp+dir)+cencell_cen(dir)
         enddo
          ! centroid_ref-centroid_free
         call find_predict_slope(npredict,mag,centroid_free,centroid_ref, &
           bfact,dx,xsten0,nhalf0,levelrz,sdim)
         
         if (mag.gt.VOF_CUTOFF*dx(1)) then

            ! order_min initialized to be 9999.
          if (order_algorithm_in(im).le.0) then
           print *,"order_algorithm_in invalid"
           stop
          else if (imaterial_override.eq.im) then ! a rigid material.
#if (STANDALONE==0)
           distmax=mag
           order_min=order_algorithm_in(im)
           critical_material=im
           override_selected=1
#else
           print *,"bust individual_MOF"
           print *,"imaterial_override ",imaterial_override
           print *,"im = ",im
           stop
#endif
          else if ((order_algorithm_in(im).lt.order_min).and. &
                   (override_selected.eq.0)) then
           distmax=mag
           order_min=order_algorithm_in(im)
           critical_material=im
          else if ((order_algorithm_in(im).eq.order_min).and. &
                   (override_selected.eq.0)) then
           if (mag.gt.distmax) then
            distmax=mag
            critical_material=im
           endif
          endif

         endif ! mag>vof_cutoff*dx(1)
        endif ! V>0 V<available vol  order=0

       enddo ! im

      else if (uncaptured_volume_vof.eq.zero) then
       if (distmax.ge.zero) then
        print *,"distmax should be negative here"
        stop
       endif
      else
       print *,"uncaptured_volume_vof invalid"
       stop
      endif 

       ! find MOF slope or find slope from prescribed distance function
      if (distmax.gt.VOF_CUTOFF*dx(1)) then

        if ((critical_material.lt.1).or. &
            (critical_material.gt.nmat)) then
         print *,"bust individual_MOF"
         print *,"critical_material=",critical_material
         print *,"nmat=",nmat
         stop
        endif

        vofcomp=(critical_material-1)*ngeom_recon+1

        do dir=1,sdim
         centroid_free(dir)=cencut_cen(dir)
         centroid_ref(dir)=mofdata(vofcomp+dir)+cencell_cen(dir)
        enddo
         ! normal points from light to dark
        call find_predict_slope(npredict,mag,centroid_free,centroid_ref, &
          bfact,dx,xsten0,nhalf0,levelrz,sdim)

        if (mag.lt.MLSVOFTOL*dx(1)) then
         print *,"mag underflow"
         stop
        endif
        do dir=1,sdim
         refcentroid(dir)=mofdata(vofcomp+dir)
        enddo
        refvfrac=mofdata(vofcomp)

        if ((imaterial_override.lt.0).or. &
            (imaterial_override.gt.nmat)) then
         print *,"imaterial_override invalid"
         stop
        endif

        if (imaterial_override.ne.critical_material) then

          ! centroidA and refcentroid relative to cell centroid of the super
          ! cell.
         call find_cut_geom_slope( &
          ls_mof, &
          lsnormal, &
          lsnormal_valid, &
          bfact,dx,xsten0,nhalf0, &
          refcentroid,refvfrac, &
          levelrz, &
          npredict, &
          continuous_mof, &
          nslope,intercept, &
          xtetlist_vof,nlist_vof, &
          xtetlist_cen,nlist_cen, &
          centroidA, &
          nmax,critical_material, &
          fastflag,nmat,nten,sdim)

        else if (imaterial_override.eq.critical_material) then
#if (STANDALONE==0)
         call find_override_slope( &
          bfact,dx,xsten0,nhalf0, &
          nslope,time,sdim) 
         use_initial_guess=0

          ! inside of individual MOF
          ! centroidA in absolute coordinate system.
         call multi_find_intercept( &
           bfact,dx,xsten0,nhalf0, &
           nslope,intercept, &
           continuous_mof, &
           xtetlist_vof,nlist_vof, &
           nmax,refvfrac,levelrz, &
           use_initial_guess,centroidA,fastflag,sdim)
          ! cencell is supercell centroid
         do dir=1,sdim
          centroidA(dir)=centroidA(dir)-cencell_cen(dir) 
         enddo
#else
         print *,"bust individual_MOF"
         print *,"sdim,nmat ",sdim,nmat
         stop
#endif
        else
         print *,"imaterial_override bust"
         stop
        endif

        mofdata(vofcomp+sdim+1)=ordermax+1
        mofdata(vofcomp+2*sdim+2)=intercept
        do dir=1,sdim
         mofdata(vofcomp+sdim+1+dir)=nslope(dir)
         multi_centroidA(critical_material,dir)=centroidA(dir)
        enddo 
        uncaptured_volume_vof=uncaptured_volume_vof-refvfrac*volcell_vof
        if (uncaptured_volume_vof.le.volcell_vof*VOF_CUTOFF) then
         uncaptured_volume_vof=zero
        endif

        ! above MOF reconstruct, below default slopes.
      else if (distmax.le.VOF_CUTOFF*dx(1)) then

         ! if single_material_takes_all=1, then
         !  distmax<0
         !  critical_material>0

        if (single_material_takes_all.eq.1) then
         critical_material=single_material_im
        else if (single_material_takes_all.eq.0) then

         max_vfrac=zero
         do im=1,nmat 
          vofcomp=(im-1)*ngeom_recon+1
          test_vfrac=mofdata(vofcomp)
          test_order=NINT(mofdata(vofcomp+sdim+1))
          if ((test_order.eq.0).and.(test_vfrac.ge.max_vfrac)) then
           max_vfrac=test_vfrac
           critical_material=im
          endif
         enddo ! im
       
        else 
         print *,"single_material_takes_all invalid"
         stop
        endif

        if ((critical_material.lt.1).or.(critical_material.gt.nmat)) then
         print *,"bust individual_MOF"
         print *,"sdim,nmat,critical_material ",sdim,nmat,critical_material
         print *,"ngeom_recon= ",ngeom_recon
         do im=1,nmat
          vofcomp=(im-1)*ngeom_recon+1
          print *,"im,vf ",im,mofdata(vofcomp)
         enddo
         stop
        endif
         
        vofcomp=(critical_material-1)*ngeom_recon+1
        test_order=NINT(mofdata(vofcomp+sdim+1))
        if (test_order.ne.0) then
         print *,"test_order invalid"
         stop
        endif

        mat_before=0
        if (ordermax.gt.0) then
         do im=1,nmat
          vofcomp_before=(im-1)*ngeom_recon+1
          if (NINT(mofdata(vofcomp_before+sdim+1)).eq.ordermax) then
           mat_before=im
          endif
         enddo
        else if (ordermax.eq.0) then
         ! do nothing
        else
         print *,"ordermax invalid"
         stop
        endif

        if ((mat_before.ge.1).and.(mat_before.le.nmat)) then
         vofcomp_before=(mat_before-1)*ngeom_recon+1
         do dir=1,sdim
          npredict(dir)=-mofdata(vofcomp_before+sdim+1+dir)
         enddo
         intercept=-mofdata(vofcomp_before+2*sdim+2)
        else if (mat_before.eq.0) then
         do dir=1,sdim
          centroid_free(dir)=cencell_cen(dir)
          centroid_ref(dir)=cencut_cen(dir)
         enddo
          ! centroid_ref-centroid_free
         call find_predict_slope(npredict,mag,centroid_free,centroid_ref, &
           bfact,dx,xsten0,nhalf0,levelrz,sdim)
         if (mag.gt.VOF_CUTOFF*dx(1)) then
          ! do nothing
         else
          do dir=1,sdim
           npredict(dir)=zero
          enddo
          npredict(1)=one
         endif
         intercept=mofdata(vofcomp)-half
        else
         print *,"mat_before invalid"
         stop
        endif

        if ((single_material_takes_all.eq.1).and. &
            (mat_before.ge.1).and.(mat_before.le.nmat)) then

         mofdata(vofcomp+sdim+1)=ordermax+1
         mofdata(vofcomp+2*sdim+2)=intercept
         do dir=1,sdim
          mofdata(vofcomp+sdim+1+dir)=npredict(dir)
           ! cencut_cen is the uncaptured centroid in absolute frame 
          centroidA(dir)=cencut_cen(dir)
          multi_centroidA(critical_material,dir)= &
           centroidA(dir)-cencell_cen(dir)
         enddo 

        else if ((single_material_takes_all.eq.0).or. &
                 (mat_before.eq.0)) then

         refvfrac=mofdata(vofcomp)
         use_initial_guess=0

          ! inside of individual MOF
          ! centroidA in absolute coordinate system.

         call multi_find_intercept( &
          bfact,dx,xsten0,nhalf0, &
          npredict,intercept, &
          continuous_mof, &
          xtetlist_vof,nlist_vof, &
          nmax,refvfrac,levelrz, &
          use_initial_guess,centroidA,fastflag,sdim)

         uncaptured_volume_vof=uncaptured_volume_vof-refvfrac*volcell_vof
         if (uncaptured_volume_vof.le.volcell_vof*VOF_CUTOFF) then
          uncaptured_volume_vof=zero
         endif
 
         if (single_material_takes_all.eq.1) then
           ! centroidA is the uncaptured centroid (absolute frame)
          do dir=1,sdim
           centroidA(dir)=cencut_cen(dir)
          enddo
         else if (single_material_takes_all.eq.0) then
          ! do nothing
         else
          print *,"single material takes all invalid"
          stop
         endif
 
         mofdata(vofcomp+sdim+1)=ordermax+1
         mofdata(vofcomp+2*sdim+2)=intercept
          ! cencell is the supercell centroid
         do dir=1,sdim
          mofdata(vofcomp+sdim+1+dir)=npredict(dir)
          multi_centroidA(critical_material,dir)= &
           centroidA(dir)-cencell_cen(dir)
         enddo 
        else
         print *,"single_material_takes_all or mat_before invalid"
         stop
        endif
      else
       print *,"distmax invalid"
       stop
      endif  

      return
      end subroutine individual_MOF


        ! called from multimaterial_MOF_intercept
        ! centroids and volume fractions are relative to center cell.
      subroutine individual_MOF_intercept( &
        ls_mof, &
        lsnormal, &
        lsnormal_valid, &
        bfact,dx,xsten0,nhalf0, &
        order_algorithm_in, &
        previous_order, &
        ngeom_recon, &
        xtetlist,nmax, &
        mofdata,levelrz, &
        imaterial_count, &
        uncaptured_volume,multi_centroidA, &
        time, &
        nmat,nten,sdim)

      use geometry_intersect_module

      IMPLICIT NONE
    
      INTEGER_T ngeom_recon,bfact,nhalf0
      INTEGER_T continuous_mof
      INTEGER_T nmat,nten,nten_test,sdim,dir,nmax
      INTEGER_T order_algorithm_in(nmat)
      INTEGER_T previous_order(nmat)
      REAL_T uncaptured_volume 
      REAL_T available_volume 
      REAL_T xsten0(-nhalf0:nhalf0,sdim)
      REAL_T dx(sdim)
      REAL_T mofdata(nmat*ngeom_recon)
      REAL_T centroidA(sdim)
      REAL_T multi_centroidA(nmat,sdim)
      REAL_T time
      INTEGER_T levelrz,im,vofcomp
      INTEGER_T imaterial_count
      REAL_T distmax
      INTEGER_T ordermax,order_min
      INTEGER_T critical_material,override_selected
      REAL_T volcut,volcell,areacut,mag
      REAL_T cencut(sdim)
      REAL_T cencell(sdim)
      REAL_T xtetlist(sdim+1,sdim,nmax)
      INTEGER_T nlist
      INTEGER_T single_material_takes_all
      INTEGER_T single_material_im
      INTEGER_T test_order
      REAL_T test_vfrac,max_vfrac
      INTEGER_T use_initial_guess
      REAL_T npredict(sdim)
      REAL_T refcentroid(sdim)
      REAL_T centroid_ref(sdim)
      REAL_T centroid_free(sdim)
      REAL_T refvfrac
      REAL_T single_volume
      REAL_T nslope(sdim),intercept
      INTEGER_T mat_before,vofcomp_before
      INTEGER_T fastflag
      REAL_T ls_mof(D_DECL(-1:1,-1:1,-1:1),nmat+nten)
      REAL_T lsnormal(nmat+nten,sdim)
      INTEGER_T lsnormal_valid(nmat+nten)

#include "mofdata.H"

      fastflag=1

      if (nhalf0.lt.1) then
       print *,"nhalf0 invalid"
       stop
      endif
      nten_test=( (nmat-1)*(nmat-1)+nmat-1 )/2
      if (nten_test.ne.nten) then
       print *,"nten invalid ind. mof int. nten nten_test ",nten,nten_test
       stop
      endif
      if (bfact.lt.1) then
       print *,"bfact invalid"
       stop
      endif
      if (uncaptured_volume.le.zero) then
       print *,"uncaptured_volume invalid"
       stop
      endif
      if (ngeom_recon.ne.2*sdim+3) then
       print *,"ngeom_recon invalid"
       stop
      endif
      if ((sdim.ne.3).and.(sdim.ne.2)) then
       print *,"sdim invalid individual_MOF_intercept"
       stop
      endif
      if ((nmat.lt.1).or.(nmat.gt.MAX_NUM_MATERIALS)) then
       print *,"nmat invalid individual mof intercept"
       stop
      endif
      if ((imaterial_count.lt.1).or.(imaterial_count.gt.nmat)) then
       print *,"imaterial_count invalid"
       stop
      endif
      if (nmax.lt.10) then
       print *,"nmax too small"
       stop
      endif

      continuous_mof=0

       ! cencell is in absolute coordinate system
      call Box_volumeFAST( &
       bfact,dx,xsten0,nhalf0, &
       volcell,cencell,levelrz,sdim)

      if ((imaterial_count.gt.1).and. &
          (imaterial_count.le.nmat)) then
       fastflag=0
      else if (imaterial_count.eq.1) then
       fastflag=1
      else 
       print *,"imaterial_count invalid"
       stop
      endif

       ! cencut is in absolute coordinate system

      if (fastflag.eq.0) then

         ! get triangulation of uncaptured space in the cell
         ! in: individual_MOF_intercept
       call tets_box_planes_super( &
         bfact,dx,xsten0,nhalf0, &
         mofdata, &
         xtetlist,nlist,nmax,nmat,continuous_mof,sdim)

       call get_cut_geom3D(xtetlist,nlist,nmax,levelrz, &
         volcut,cencut,sdim)

      else if (fastflag.eq.1) then
       volcut=volcell
       do dir=1,sdim
        cencut(dir)=cencell(dir)  ! cencell is centroid of supercell.
       enddo
      else
       print *,"fastflag invalid individual MOF_intercept"
       stop
      endif      

      available_volume=uncaptured_volume

      if (abs(volcut-uncaptured_volume).gt.VOF_CUTOFF*volcell) then
        print *,"volcut invalid individual mof intercept"
        print *,"volcut ",volcut
        print *,"uncaptured_volume ",uncaptured_volume
        stop
      endif

        ! figure out the next material to fill the unoccupied region.
      distmax=-one
      order_min=9999
      ordermax=0

      do im=1,nmat
       vofcomp=(im-1)*ngeom_recon+1
       if (NINT(mofdata(vofcomp+sdim+1)).gt.ordermax) then
        ordermax=NINT(mofdata(vofcomp+sdim+1))
       endif
      enddo

      if (ordermax.ge.nmat) then
       print *,"all the materials already initialized"
       stop
      endif
      if (ordermax.lt.0) then
       print *,"ordermax invalid"
       stop
      endif

      critical_material=0
      single_material_takes_all=0
      single_material_im=0
      do im=1,nmat
       vofcomp=(im-1)*ngeom_recon+1
       if (mofdata(vofcomp+sdim+1).eq.zero) then
        single_volume=mofdata(vofcomp)*volcell
        if (single_volume.ge.(one-VOF_CUTOFF)*available_volume) then
         if ((single_material_takes_all.ne.0).or. &
             (single_material_im.ne.0)) then
          print *,"cannot have two materials at once"
          print *,"single_material_takes_all ",single_material_takes_all
          print *,"im ",im
          print *,"single vol ",single_volume
          print *,"uncapt vol ",uncaptured_volume
          print *,"uncapt volfrac ",uncaptured_volume/volcell
          stop
         endif
         single_material_takes_all=1
         single_material_im=im
         critical_material=im
         distmax=-one  ! tells code below not to search for a slope
        endif
       endif  ! material not already processed.
      enddo ! im

      if (single_material_takes_all.eq.1) then
       uncaptured_volume=zero
       available_volume=zero
      else if (single_material_takes_all.eq.0) then
       ! do nothing
      else 
       print *,"bust individual_MOF_intercept"
       print *,"sdim,nmat ",sdim,nmat
       stop
      endif

       
         ! if uncaptured_volume=0, then there is no need to find the
         ! slope since "single_material_takes_all=1". 
      if (uncaptured_volume.gt.zero) then

        ! find unprocessed material whose moment is furthest from cencut. 

       override_selected=0

       do im=1,nmat
        vofcomp=(im-1)*ngeom_recon+1

        single_volume=mofdata(vofcomp)*volcell

         ! only find the slope if refvfrac<>0 and refvfrac<available vol.
        if ((mofdata(vofcomp+sdim+1).eq.zero).and. &
            (single_volume.ge.VOF_CUTOFF*volcell).and. &
            (single_volume.le.(one-VOF_CUTOFF)*available_volume)) then

         do dir=1,sdim
          centroid_free(dir)=cencut(dir)
          centroid_ref(dir)=mofdata(vofcomp+dir)+cencell(dir)
         enddo
          ! centroid_ref-centroid_free
         call find_predict_slope(npredict,mag,centroid_free,centroid_ref, &
           bfact,dx,xsten0,nhalf0,levelrz,sdim)

         if (mag.gt.VOF_CUTOFF*dx(1)) then

            ! order_min initialized to be 9999.
          if (order_algorithm_in(im).le.0) then
           print *,"order_algorithm_in invalid"
           stop
          else if (imaterial_override.eq.im) then ! a rigid material.
#if (STANDALONE==0)
           distmax=mag
           order_min=order_algorithm_in(im)
           critical_material=im
           override_selected=1
#else
           print *,"bust individual_MOF_intercept 2"
           print *,"sdim,nmat ",sdim,nmat
           stop
#endif
          else if ((order_algorithm_in(im).lt.order_min).and. &
                   (override_selected.eq.0)) then
           distmax=mag
           order_min=order_algorithm_in(im)
           critical_material=im
          else if ((order_algorithm_in(im).eq.order_min).and. &
                   (override_selected.eq.0)) then
           if (mag.gt.distmax) then
            distmax=mag
            critical_material=im
           endif
          endif

         endif
        endif ! V>0 V<available vol  order=0

       enddo ! im

      else if (uncaptured_volume.eq.zero) then
       if (distmax.ge.zero) then
        print *,"distmax should be negative here"
        stop
       endif
      else
       print *,"uncaptured_volume invalid"
       stop
      endif 

       ! find MOF slope or find slope from prescribed distance function
      if (distmax.gt.VOF_CUTOFF*dx(1)) then

        if ((critical_material.lt.1).or.(critical_material.gt.nmat)) then
         print *,"bust individual_MOF_intercept 3"
         print *,"sdim,nmat ",sdim,nmat
         stop
        endif

        vofcomp=(critical_material-1)*ngeom_recon+1

         ! normal points from light to dark
        if (previous_order(critical_material).eq.0) then
         do dir=1,sdim
          centroid_free(dir)=cencut(dir)
          centroid_ref(dir)=mofdata(vofcomp+dir)+cencell(dir)
         enddo
          ! centroid_ref-centroid_free
         call find_predict_slope(npredict,mag,centroid_free,centroid_ref, &
          bfact,dx,xsten0,nhalf0,levelrz,sdim)
        else if (previous_order(critical_material).gt.0) then
         mag=zero
         do dir=1,sdim
          npredict(dir)=mofdata(vofcomp+sdim+1+dir)
          mag=mag+npredict(dir)**2
         enddo
         mag=sqrt(mag)
         if (mag.gt.zero) then
          do dir=1,sdim
           npredict(dir)=npredict(dir)/mag
          enddo
         endif
        else
         print *,"previous_order invalid"
         stop
        endif 
        if (mag.lt.MLSVOFTOL*dx(1)) then
         print *,"mag underflow"
         stop
        endif
        do dir=1,sdim
         refcentroid(dir)=mofdata(vofcomp+dir)
        enddo
        refvfrac=mofdata(vofcomp)

        if ((imaterial_override.lt.0).or. &
            (imaterial_override.gt.nmat)) then
         print *,"imaterial_override invalid"
         stop
        endif

        if (imaterial_override.ne.critical_material) then
          ! centroidA and refcentroid relative to cell centroid of the super
          ! cell.

         if (previous_order(critical_material).eq.0) then
          call find_cut_geom_slope( &
           ls_mof, &
           lsnormal, &
           lsnormal_valid, &
           bfact,dx,xsten0,nhalf0, &
           refcentroid,refvfrac,levelrz, &
           npredict, &
           continuous_mof, &
           nslope,intercept, &
           xtetlist,nlist, &
           xtetlist,nlist, &
           centroidA, &
           nmax,critical_material, &
           fastflag,nmat,nten,sdim)
         else if (previous_order(critical_material).gt.0) then
          do dir=1,sdim
           nslope(dir)=npredict(dir)
          enddo

          use_initial_guess=0

           ! inside of individual MOF intercept
           ! centroidA in absolute coordinate system.
          call multi_find_intercept( &
           bfact,dx,xsten0,nhalf0, &
           nslope,intercept, &
           continuous_mof, &
           xtetlist,nlist, &
           nmax,refvfrac,levelrz, &
           use_initial_guess,centroidA,fastflag,sdim)
          do dir=1,sdim
           centroidA(dir)=centroidA(dir)-cencell(dir) 
          enddo
         else
          print *,"previous_order invalid"
          stop
         endif

        else if (imaterial_override.eq.critical_material) then
#if (STANDALONE==0)
         call find_override_slope(bfact,dx,xsten0,nhalf0, &
          nslope,time,sdim) 

         use_initial_guess=0

          ! inside of individual MOF intercept
          ! centroidA in absolute coordinate system.
         call multi_find_intercept( &
           bfact,dx,xsten0,nhalf0, &
           nslope,intercept, &
           continuous_mof, &
           xtetlist,nlist, &
           nmax,refvfrac,levelrz, &
           use_initial_guess,centroidA,fastflag,sdim)
         do dir=1,sdim
          centroidA(dir)=centroidA(dir)-cencell(dir) 
         enddo
#else
         print *,"bust individual_MOF_intercept 4"
         print *,"sdim,nmat ",sdim,nmat
         stop
#endif

        else
         print *,"imaterial_override bust"
         stop
        endif

        mofdata(vofcomp+sdim+1)=ordermax+1
        mofdata(vofcomp+2*sdim+2)=intercept
        do dir=1,sdim
         mofdata(vofcomp+sdim+1+dir)=nslope(dir)
         multi_centroidA(critical_material,dir)=centroidA(dir)
        enddo 
        uncaptured_volume=uncaptured_volume-refvfrac*volcell
        if (uncaptured_volume.le.volcell*VOF_CUTOFF) then
         uncaptured_volume=zero
        endif

        ! above MOF reconstruct, below default slopes.
      else if (distmax.le.VOF_CUTOFF*dx(1)) then

         ! if single_material_takes_all=1, then
         !  distmax<0
         !  critical_material>0

        if (single_material_takes_all.eq.1) then
         critical_material=single_material_im
        else if (single_material_takes_all.eq.0) then

         max_vfrac=zero
         do im=1,nmat
          vofcomp=(im-1)*ngeom_recon+1
          test_vfrac=mofdata(vofcomp)
          test_order=NINT(mofdata(vofcomp+sdim+1))
          if ((test_order.eq.0).and.(test_vfrac.ge.max_vfrac)) then
           max_vfrac=test_vfrac
           critical_material=im
          endif
         enddo ! im

        else 
         print *,"single_material_takes_all invalid"
         stop
        endif

        if ((critical_material.lt.1).or.(critical_material.gt.nmat)) then
         print *,"bust individual_MOF_intercept 5"
         print *,"sdim,nmat,critical_material ",sdim,nmat,critical_material
         stop
        endif
         
        vofcomp=(critical_material-1)*ngeom_recon+1
        test_order=NINT(mofdata(vofcomp+sdim+1))
        if (test_order.ne.0) then
         print *,"test_order invalid"
         stop
        endif

        mat_before=0
        if (ordermax.gt.0) then
         do im=1,nmat
          vofcomp_before=(im-1)*ngeom_recon+1
          if (NINT(mofdata(vofcomp_before+sdim+1)).eq.ordermax) then
           mat_before=im
          endif
         enddo
        else if (ordermax.eq.0) then
         ! do nothing
        else
         print *,"ordermax invalid"
         stop
        endif

        if ((mat_before.ge.1).and.(mat_before.le.nmat)) then
         vofcomp_before=(mat_before-1)*ngeom_recon+1
         do dir=1,sdim
          npredict(dir)=-mofdata(vofcomp_before+sdim+1+dir)
         enddo
         intercept=-mofdata(vofcomp_before+2*sdim+2)
        else if (mat_before.eq.0) then
         do dir=1,sdim
          centroid_free(dir)=cencell(dir)
          centroid_ref(dir)=cencut(dir)
         enddo
          ! centroid_ref-centroid_free
         call find_predict_slope(npredict,mag,centroid_free,centroid_ref, &
           bfact,dx,xsten0,nhalf0,levelrz,sdim)

         if (mag.gt.VOF_CUTOFF*dx(1)) then
          ! do nothing
         else
          do dir=1,sdim
           npredict(dir)=zero
          enddo
          npredict(1)=one
         endif

         intercept=mofdata(vofcomp)-half
        else
         print *,"mat_before invalid"
         stop
        endif

        if ((single_material_takes_all.eq.1).and. &
            (mat_before.ge.1).and.(mat_before.le.nmat)) then

         mofdata(vofcomp+sdim+1)=ordermax+1
         mofdata(vofcomp+2*sdim+2)=intercept
         do dir=1,sdim
          mofdata(vofcomp+sdim+1+dir)=npredict(dir)
           ! cencut is the uncaptured centroid in absolute frame
          centroidA(dir)=cencut(dir)
          multi_centroidA(critical_material,dir)= &
           centroidA(dir)-cencell(dir)
         enddo 

        else if ((single_material_takes_all.eq.0).or. &
                 (mat_before.eq.0)) then

         refvfrac=mofdata(vofcomp)
         use_initial_guess=0

          ! inside of individual MOF
          ! centroidA in absolute coordinate system.

         call multi_find_intercept( &
          bfact,dx,xsten0,nhalf0, &
          npredict,intercept, &
          continuous_mof, &
          xtetlist,nlist, &
          nmax,refvfrac,levelrz, &
          use_initial_guess,centroidA,fastflag,sdim)

         uncaptured_volume=uncaptured_volume-refvfrac*volcell
         if (uncaptured_volume.le.volcell*VOF_CUTOFF) then
          uncaptured_volume=zero
         endif
 
         if (single_material_takes_all.eq.1) then
           ! centroidA is the uncaptured centroid
          do dir=1,sdim
           centroidA(dir)=cencut(dir)
          enddo
         else if (single_material_takes_all.eq.0) then
          ! do nothing
         else
          print *,"single material takes all invalid"
          stop
         endif
 
         mofdata(vofcomp+sdim+1)=ordermax+1
         mofdata(vofcomp+2*sdim+2)=intercept
         do dir=1,sdim
          mofdata(vofcomp+sdim+1+dir)=npredict(dir)
          multi_centroidA(critical_material,dir)= &
           centroidA(dir)-cencell(dir)
         enddo 
        else
         print *,"single_material_takes_all or mat_before invalid"
         stop
        endif
 
      else
       print *,"distmax invalid"
       stop
      endif  

      return
      end subroutine individual_MOF_intercept


        ! n1d=1 => im material on top
        ! n1d=-1 => im material on bottom
      subroutine get_col_ht_LS( &
       bfact,dx,xsten0, &
       csten,csten_x, &
       lsdata,ht,dircrit,n1d, &
       im,im_opp,nmat,levelrz,sdim)
      use geometry_intersect_module

      IMPLICIT NONE

      INTEGER_T sdim
      INTEGER_T csten,csten_x,bfact
      REAL_T dx(sdim)
      REAL_T xsten0(-csten_x:csten_x,sdim)
      REAL_T lsdata(-csten:csten,nmat)
      REAL_T ht,n1d
      INTEGER_T dircrit
      INTEGER_T levelrz
      INTEGER_T im,im_opp,nmat

      INTEGER_T lfound,l,lmin,lmax,lcrit
      REAL_T xbottom,xtop

      INTEGER_T dir2
      REAL_T ls1,ls2,x1,x2,slope
      INTEGER_T im_crit
      REAL_T charfn_crit
      REAL_T charfn(-csten:csten)
      REAL_T LS(nmat)
      INTEGER_T imloop

      if (csten.lt.1) then
       print *,"csten invalid"
       stop
      endif
      if (csten_x.ne.2*csten+1) then
       print *,"csten_x invalid"
       stop
      endif
      if (csten.gt.4) then
       print *,"csten too big"
       stop
      endif
      if ((nmat.lt.1).or.(nmat.gt.MAX_NUM_MATERIALS)) then
       print *,"nmat invalid get col ht LS"
       print *,"nmat= ",nmat
       stop
      endif
      if (im.eq.im_opp) then
       print *,"im and im_opp invalid"
       stop
      endif
      if ((im.lt.1).or.(im.gt.nmat)) then
       print *,"im invalid"
       stop
      endif
      if ((im_opp.lt.1).or.(im_opp.gt.nmat)) then
       print *,"im_opp invalid get_col_ht_LS im_opp=",im_opp
       stop
      endif

      if (levelrz.eq.0) then
       ! do nothing
      else if (levelrz.eq.1) then
       if (sdim.ne.2) then
        print *,"sdim invalid"
        stop
       endif
       if (xsten0(0,1).le.zero) then
        print *,"cannot have r<0"
        stop
       endif
      else if (levelrz.eq.3) then
       if (xsten0(0,1).le.zero) then
        print *,"cannot have r<0"
        stop
       endif
      else 
       print *,"levelrz invalid get col ht ls"
       stop
      endif

      if ((dircrit.lt.1).or.(dircrit.gt.sdim)) then
       print *,"dircrit invalid get_col_ht_LS dircrit=",dircrit
       stop
      endif
      if ((n1d.ne.one).and.(n1d.ne.-one)) then
       print *,"n1d invalid"
       stop
      endif

      lmin=-csten
      lmax=csten
      if ((levelrz.eq.1).or.(levelrz.eq.3)) then
       if (dircrit.eq.1) then
        do while (xsten0(2*lmin,dircrit).lt.zero)
         lmin=lmin+1
        enddo
       endif
      else if (levelrz.eq.0) then
       ! do nothing
      else
       print *,"levelrz invalid get col ht ls 2"
       stop
      endif 
      xbottom=xsten0(2*lmin-1,dircrit)
      xtop=xsten0(2*lmax+1,dircrit)

      lcrit=0
      lfound=0

      do l=lmin,lmax

       do imloop=1,nmat
        LS(imloop)=lsdata(l,imloop)
       enddo
       if (LS(im).ge.LS(im_opp)) then
        charfn(l)=one
       else
        charfn(l)=-one
       endif
        
      enddo  ! l=lmin,lmax

       ! n1d>0 if im material on top and im_opp material on bottom.
       ! n1d<0 otherwise
      do l=0,csten-1
       if ((lfound.eq.0).and.(l+1.le.lmax)) then 
        ls1=charfn(l)
        ls2=charfn(l+1)
        if ((ls1*ls2.lt.zero).and.((ls2-ls1)*n1d.gt.zero)) then
         lcrit=l
         lfound=1
         ls1=lsdata(l,im)
         ls2=lsdata(l+1,im)
         x1=xsten0(2*l,dircrit)
         x2=xsten0(2*l+2,dircrit)
          ! LS=LS1+slope(x-x1)  xzero=-LS1/slope+x1
         if (ls1.eq.zero) then
          ht=x1
         else if (ls2.eq.zero) then 
          ht=x2
         else
          slope=(ls1-ls2)/(x1-x2)
          ht=x1-ls1/slope
         endif
        endif   
       endif
       if ((lfound.eq.0).and.(-(l+1).ge.lmin)) then 
        ls1=charfn(-l)
        ls2=charfn(-(l+1))
        if ((ls1*ls2.lt.zero).and.((ls1-ls2)*n1d.gt.zero)) then
         lcrit=-(l+1)
         lfound=1
         ls1=lsdata(-l,im)
         ls2=lsdata(-(l+1),im)
         x1=xsten0(-2*l,dircrit)
         x2=xsten0(-(2*l+2),dircrit)
         if (ls1.eq.zero) then
          ht=x1
         else if (ls2.eq.zero) then
          ht=x2
         else
          slope=(ls1-ls2)/(x1-x2)
          ht=x1-ls1/slope
         endif
        endif   
       endif
      enddo ! l=0,csten-1


         ! given volume for column, find the interface.
         ! for RZ, if dircrit=1, then column extends to r=0

      if (lfound.eq.1) then
       if ((lcrit.lt.lmin).or.(lcrit+1.gt.lmax)) then
        print *,"lcrit invalid"
        stop
       endif
      else if (lfound.eq.0) then
       ! do nothing
      else 
       print *,"lfound invalid"
       stop
      endif

       ! find the height of the im_crit material.

      if (n1d.eq.-one) then  ! im material on bottom
       im_crit=im
       charfn_crit=one
      else if (n1d.eq.one) then  ! im_opp material on bottom
       im_crit=im_opp
       charfn_crit=-one
      else
       print *,"n1d invalid"
       stop
      endif

      if (lfound.eq.1) then
       ! do nothing
      else if (lfound.eq.0) then
       if (lsdata(0,im_crit).ge.zero) then
        ht=xtop
       else
        ht=xbottom
       endif
      else
       print *,"lfound invalid"
       stop
      endif

      return
      end subroutine get_col_ht_LS


! normal points from (-) to (+)   phi=n dot (x-x0) + intercept
! xsten0(0,dir) is cell center (not cell centroid)
      subroutine find_cut_geom_slope_CLSVOF( &
        ls_mof, &
        lsnormal, &
        lsnormal_valid, &
        bfact,dx,xsten0,nhalf0, &
        im, &
        dxmaxLS, &
        levelrz, &
        nmat,nten,sdim)

      use geometry_intersect_module
      use global_utility_module

      IMPLICIT NONE

      INTEGER_T bfact,nhalf0,im
      INTEGER_T nmat,nten,nten_test,sdim,levelrz
      REAL_T xsten0(-nhalf0:nhalf0,sdim)
      REAL_T dx(sdim)
      REAL_T xpoint(sdim)
      REAL_T dxmaxLS,cutoff,m1,m2
      REAL_T nsimple(sdim)
      REAL_T nn(sdim)
      REAL_T distsimple,dist,LSWT
      REAL_T w(D_DECL(-1:1,-1:1,-1:1))
      REAL_T aa(sdim+1,sdim+1)
      REAL_T xx(sdim+1)
      REAL_T bb(sdim+1)
      INTEGER_T matstatus
      REAL_T wx,wy,wz
      INTEGER_T i,j,k
      INTEGER_T i1,j1,k1
      INTEGER_T klo,khi,dir

      REAL_T ls_mof(D_DECL(-1:1,-1:1,-1:1),nmat+nten)
      REAL_T lsnormal(nmat+nten,sdim)
      INTEGER_T lsnormal_valid(nmat+nten)

#include "mofdata.H"

      if (nhalf0.lt.3) then
       print *,"nhalf0 invalid"
       stop
      endif
      nten_test=( (nmat-1)*(nmat-1)+nmat-1 )/2
      if (nten_test.ne.nten) then
       print *,"nten bad cutgeomslopeclsvof nten nten_test ",nten,nten_test
       stop
      endif

      if (dxmaxLS.le.zero) then
       print *,"dxmaxLS invalid"
       stop
      endif

      if (bfact.lt.1) then
       print *,"bfact invalid"
       stop
      endif

      if ((sdim.ne.3).and.(sdim.ne.2)) then
       print *,"sdim invalid find_cut_geom_slope_CLSVOF"
       stop
      endif

      if ((im.lt.1).or.(im.gt.nmat+nten)) then
       print *,"im invalid"
       stop
      endif
      if (nmat.lt.1) then
       print *,"nmat invalid"
       stop
      endif

      if (sdim.eq.3) then
       klo=-1
       khi=1
       cutoff=sqrt(three)*dxmaxLS
      else if (sdim.eq.2) then
       klo=0
       khi=0
       cutoff=sqrt(two)*dxmaxLS
      else
       print *,"dimension bust"
       stop
      endif

      lsnormal_valid(im)=1

      do i=-1,1
      do j=-1,1
      do k=klo,khi
       if (abs(ls_mof(D_DECL(i,j,k),im)).gt.three*dxmaxLS) then
        lsnormal_valid(im)=0
       endif
      enddo
      enddo
      enddo

      if (lsnormal_valid(im).eq.1) then

       dir=1
       nsimple(dir)=(ls_mof(D_DECL(1,0,0),im)- &
                     ls_mof(D_DECL(-1,0,0),im))/ &
                    (xsten0(2,dir)-xsten0(-2,dir))
       dir=2
       nsimple(dir)=(ls_mof(D_DECL(0,1,0),im)- &
                     ls_mof(D_DECL(0,-1,0),im))/ &
                    (xsten0(2,dir)-xsten0(-2,dir))

       if (sdim.eq.3) then
        dir=sdim
        nsimple(dir)=(ls_mof(D_DECL(0,0,1),im)- &
                      ls_mof(D_DECL(0,0,-1),im))/ &
                     (xsten0(2,dir)-xsten0(-2,dir))
       endif

       distsimple=zero
       do dir=1,sdim
        distsimple=distsimple+nsimple(dir)**2
       enddo
       distsimple=sqrt(distsimple)
       if (distsimple.gt.zero) then
        do dir=1,sdim
         nsimple(dir)=nsimple(dir)/distsimple
        enddo
       else if (distsimple.eq.zero) then
        lsnormal_valid(im)=0
       else
        print *,"distsimple invalid"
        stop
       endif

       if (lsnormal_valid(im).eq.1) then
        do i=-1,1
        do j=-1,1
        do k=klo,khi
         wx=twelve
         wy=twelve
         wz=twelve
         if (i.ne.0) then
          wx=one
         endif
         if (j.ne.0) then
          wy=one
         endif
         if (k.ne.0) then
          wz=one
         endif

         wx=wx*(xsten0(2*i+1,1)-xsten0(2*i-1,1))
         wy=wy*(xsten0(2*j+1,2)-xsten0(2*j-1,2))
         if (sdim.eq.3) then
          wz=wz*(xsten0(2*k+1,sdim)-xsten0(2*k-1,sdim))
         endif

         LSWT=abs(ls_mof(D_DECL(i,j,k),im))
         w(D_DECL(i,j,k))=hsprime(LSWT,cutoff)*wx*wy*wz
        enddo
        enddo
        enddo ! i,j,k

        do i=1,sdim+1
         do j=1,sdim+1
          aa(i,j)=zero
         enddo
         bb(i)=zero
        enddo

        do i1=-1,1
        do j1=-1,1
        do k1=klo,khi
         xpoint(1)=xsten0(2*i1,1)-xsten0(0,1)
         xpoint(2)=xsten0(2*j1,2)-xsten0(0,2)
         if (sdim.eq.3) then
          xpoint(sdim)=xsten0(2*k1,sdim)-xsten0(0,sdim)
         endif
 
         do i=1,sdim+1
          if (i.eq.sdim+1) then
           m1=one
          else
           m1=xpoint(i)
          endif

          do j=1,sdim+1
           if (j.eq.sdim+1) then
            m2=one
           else
            m2=xpoint(j)
           endif
           aa(i,j)=aa(i,j)+w(D_DECL(i1,j1,k1))*m1*m2
          enddo ! j=1,sdim+1 
          bb(i)=bb(i)+ &
           w(D_DECL(i1,j1,k1))*m1*ls_mof(D_DECL(i1,j1,k1),im)
         enddo ! i=1,sdim+1
        enddo
        enddo
        enddo ! i1,j1,k1

        call matrix_solve(aa,xx,bb,matstatus,sdim+1)
        if (matstatus.eq.1) then
         dist=zero
         do dir=1,sdim
          nn(dir)=xx(dir)
          dist=dist+nn(dir)**2
         enddo
         dist=sqrt(dist)
         if (dist.gt.zero) then
          do dir=1,sdim
           nn(dir)=nn(dir)/dist
          enddo
         else if (dist.eq.zero) then
          do dir=1,sdim
           nn(dir)=nsimple(dir)
          enddo
         else
          print *,"dist invalid"
          stop
         endif
        else if (matstatus.eq.0) then
         do dir=1,sdim
          nn(dir)=nsimple(dir)
         enddo
        else
         print *,"matstatus invalid"
         stop
        endif
        do dir=1,sdim
         lsnormal(im,dir)=nn(dir)
        enddo
       else if (lsnormal_valid(im).eq.0) then
        ! do nothing
       else
        print *,"lsnormal_valid(im) invalid"
        stop
       endif
      else if (lsnormal_valid(im).eq.0) then
       ! do nothing
      else
       print *,"lsnormal_valid(im) invalid"
       stop
      endif
       
      return
      end subroutine find_cut_geom_slope_CLSVOF

      subroutine nfact(n,nf)
      IMPLICIT NONE

      INTEGER_T n,nf,i

      if (n.le.0) then
       print *,"n invalid"
       stop
      else
       nf=1
       do i=2,n
        nf=nf*i
       enddo
      endif
 
      return
      end subroutine nfact

      subroutine push_order_stack(order_stack,order_stack_count, &
       temp_order,n_orderings,n_ndef,nmat)
      IMPLICIT NONE

      INTEGER_T order_stack_count
      INTEGER_T n_orderings,n_ndef
      INTEGER_T nmat
      INTEGER_T temp_order(nmat)
      INTEGER_T order_stack(n_orderings,n_ndef)
      INTEGER_T i
     
      if (order_stack_count.ge.n_orderings) then
       print *,"order_stack_count too big - push"
       stop
      endif

      order_stack_count=order_stack_count+1
      do i=1,n_ndef
       order_stack(order_stack_count,i)=temp_order(i)
      enddo
   
      return
      end subroutine push_order_stack


      subroutine pop_order_stack(order_stack,order_stack_count, &
       temp_order,n_orderings,n_ndef,nmat)
      IMPLICIT NONE

      INTEGER_T order_stack_count
      INTEGER_T n_orderings,n_ndef
      INTEGER_T nmat
      INTEGER_T temp_order(nmat)
      INTEGER_T order_stack(n_orderings,n_ndef)
      INTEGER_T i
     
      if (order_stack_count.gt.n_orderings) then
       print *,"order_stack_count too big - pop"
       stop
      endif
      do i=1,n_ndef
       temp_order(i)=order_stack(order_stack_count,i)
      enddo
      order_stack_count=order_stack_count-1
   
      return
      end subroutine pop_order_stack


! normal points from light to dark   phi=n dot (x-x0) + intercept
! vof, ref centroid, order,slope,intercept  x nmat
! ref centroid,multi_centroidA relative to supercell centroid(not cell center).
! xcell is cell center (not cell centroid)
!
! if continuous_mof=1:
!  centroids: 3x3 super cell
!  vfrac    : 3x3 super cell
!
! if continuous_mof=2:
!  centroids: 3x3 super cell
!  vfrac    : center cell
!
! xtetlist is workspace data (use POLYGON_LIST_MAX)
!  i.e. nmax=POLYGON_LIST_MAX
! mofdata:  nmat*ngeom_recon elements; for each material:
!   vfrac (input), ref centroid (input), order (output), slope (output),
!   intercept (output)


      subroutine multimaterial_MOF( &
        bfact,dx,xsten0,nhalf0, &
        mof_verbose, &
        use_ls_data, &
        LS_stencil, &
        xtetlist_vof, &
        xtetlist_cen, &
        nmax, &
        mofdata, &
        multi_centroidA, &
        time, &
        continuous_mof, &
        levelrz,nmat,nten,sdim, &
        ngeom_recon_in)

      use global_utility_module
      use geometry_intersect_module

      IMPLICIT NONE

      INTEGER_T, INTENT (IN) :: bfact,nhalf0
      INTEGER_T, INTENT (IN) :: use_ls_data
      INTEGER_T, INTENT (IN) :: mof_verbose
      INTEGER_T, INTENT (IN) :: nmat
      INTEGER_T, INTENT (IN) :: nten
      INTEGER_T, INTENT (IN) :: sdim
      INTEGER_T, INTENT (IN) :: nmax
      INTEGER_T, INTENT (IN) :: ngeom_recon_in
      INTEGER_T, INTENT (IN) :: continuous_mof

      REAL_T, INTENT (IN) :: LS_stencil(D_DECL(-1:1,-1:1,-1:1),nmat+nten)

      REAL_T, INTENT (IN), DIMENSION(sdim) :: dx

      REAL_T, INTENT (INOUT), DIMENSION(sdim+1,sdim,nmax) :: xtetlist_vof
      REAL_T, INTENT (INOUT), DIMENSION(sdim+1,sdim,nmax) :: xtetlist_cen
      REAL_T, INTENT (IN), DIMENSION(-nhalf0:nhalf0,sdim) :: xsten0
      REAL_T, INTENT (INOUT), DIMENSION(nmat*ngeom_recon_in) :: mofdata
      REAL_T, DIMENSION(nmat*ngeom_recon_in) :: mofdata_in
      REAL_T, INTENT (OUT), DIMENSION(nmat,sdim) :: multi_centroidA
      REAL_T, INTENT (IN) :: time

      INTEGER_T imaterial2
      INTEGER_T levelrz,imaterial,vofcomp,dir
      INTEGER_T imaterial_count
      REAL_T uncaptured_volume_vof
      REAL_T uncaptured_centroid_vof(sdim)
      REAL_T uncaptured_volume_cen
      REAL_T uncaptured_centroid_cen(sdim)
      REAL_T xref_mat(sdim)
      REAL_T xact_mat(sdim)
      REAL_T xref_matT(sdim)
      REAL_T xact_matT(sdim)
      INTEGER_T single_material
      REAL_T remaining_vfrac
      REAL_T nrecon(sdim)
      INTEGER_T ngeom_recon
      INTEGER_T order_algorithm_in(nmat)
      INTEGER_T num_materials_cell
      INTEGER_T, dimension(:,:), allocatable :: order_array
      INTEGER_T, dimension(:,:), allocatable :: order_stack
      REAL_T, dimension(:,:), allocatable :: mofdata_array
      REAL_T, dimension(:,:,:), allocatable :: centroidA_array
      REAL_T, dimension(:), allocatable :: moferror_array
      INTEGER_T order_count,order_stack_count
      INTEGER_T irank,iflex,jflex,kflex,is_valid
      INTEGER_T n_ndef,place_index,n_orderings
      INTEGER_T placeholder(nmat)
      INTEGER_T placelist(nmat)
      INTEGER_T flexlist(nmat)
      INTEGER_T temp_order(nmat)
      INTEGER_T argmin_order
      REAL_T min_error,mof_err
      INTEGER_T alloc_flag
      REAL_T voftest(nmat)
      INTEGER_T mixed_cell
      INTEGER_T nten_test
      INTEGER_T i1,j1,k1,k1lo,k1hi,im_opp,iten
      INTEGER_T im_main,im_main_opp
      REAL_T dxmaxLS
      REAL_T, dimension(D_DECL(:,:,:),:), allocatable :: ls_mof
      REAL_T, dimension(:,:), allocatable :: lsnormal
      INTEGER_T, dimension(:), allocatable :: lsnormal_valid

#include "mofdata.H"

      if (bfact.lt.1) then
       print *,"bfact invalid"
       stop
      endif

      if (nhalf0.lt.3) then
       print *,"nhalf0 invalid"
       stop
      endif
      alloc_flag=0

      nten_test=( (nmat-1)*(nmat-1)+nmat-1 )/2
      if (nten_test.ne.nten) then
       print *,"nten invalid multimof nten nten_test ",nten,nten_test
       stop
      endif

      ngeom_recon=2*sdim+3
      if (ngeom_recon.ne.ngeom_recon_in) then
       print *,"ngeom_recon_in invalid"
       stop
      endif

      if ((sdim.ne.3).and.(sdim.ne.2)) then
       print *,"sdim invalid multimaterial_MOF"
       stop
      endif
      if ((continuous_mof.ne.0).and. &
          (continuous_mof.ne.1).and. &
          (continuous_mof.ne.2)) then
       print *,"continuous_mof invalid"
       stop
      endif
      if ((nmat.lt.1).or.(nmat.gt.MAX_NUM_MATERIALS)) then
       print *,"nmat invalid multimaterial mof"
       print *,"nmat= ",nmat
       stop
      endif
      if (nmax.lt.10) then
       print *,"nmax too small in multimaterial_MOF"
       stop
      endif
      if ((use_ls_data.ne.0).and.(use_ls_data.ne.1)) then
       print *,"use_ls_data invalid"
       stop
      endif

      call get_dxmaxLS(dx,bfact,dxmaxLS)

      if (sdim.eq.2) then
       k1lo=0
       k1hi=0
      else if (sdim.eq.3) then
       k1lo=-1
       k1hi=1
      else
       print *,"dimension bust"
       stop
      endif

      mixed_cell=0
      do imaterial=1,nmat
       vofcomp=(imaterial-1)*ngeom_recon+1
       voftest(imaterial)=mofdata(vofcomp)
       if ((voftest(imaterial).gt.0.1).and. &
           (voftest(imaterial).lt.0.9)) then
        mixed_cell=1
       endif
      enddo 


      allocate(ls_mof(D_DECL(-1:1,-1:1,-1:1),nmat+nten))
      allocate(lsnormal(nmat+nten,sdim))
      allocate(lsnormal_valid(nmat+nten))

      if (use_ls_data.eq.1) then
       do i1=-1,1
       do j1=-1,1
       do k1=k1lo,k1hi
       do imaterial=1,nmat+nten
        ls_mof(D_DECL(i1,j1,k1),imaterial)= &
         LS_stencil(D_DECL(i1,j1,k1),imaterial) 
       enddo
       enddo
       enddo
       enddo
       do imaterial=1,nmat+nten
        if ((imaterial.ge.1).and.(imaterial.le.nmat)) then
         im_main=imaterial
         im_main_opp=imaterial
        else if ((imaterial.gt.nmat).and.(imaterial.le.nmat+nten)) then
         iten=imaterial-nmat
         call get_inverse_iten(im_main,im_main_opp,iten,nmat)
        else
         print *,"imaterial invalid"
         stop
        endif 
        if ((voftest(im_main).le.VOF_CUTOFF).or. &
            (voftest(im_main_opp).le.VOF_CUTOFF)) then
         lsnormal_valid(imaterial)=0
        else
         call find_cut_geom_slope_CLSVOF( &
          ls_mof, &
          lsnormal, &
          lsnormal_valid, &
          bfact,dx,xsten0,nhalf0, &
          imaterial, &
          dxmaxLS, &
          levelrz, &
          nmat,nten,sdim)
        endif
       enddo ! imaterial
      else if (use_ls_data.eq.0) then
       do imaterial=1,nmat+nten
        lsnormal_valid(imaterial)=0
       enddo
      else
       print *,"use_ls_data invalid"
       stop
      endif

      if ((mof_verbose.eq.1).or. &
          ((mixed_cell.eq.1).and.(mof_verbose.eq.2))) then
       print *,"BEFORE BEFORE"
       print *,"nmax = ",nmax
       print *,"levelrz = ",levelrz
       print *,"nmat = ",nmat
       print *,"time = ",time
       print *,"sdim = ",sdim
       print *,"continuous_mof = ",continuous_mof
       print *,"ngeom_recon = ",ngeom_recon
       do imaterial=1,nmat*ngeom_recon
        print *,"i,mofdata ",imaterial,mofdata(imaterial)
       enddo
       do imaterial=1,nmat
        print *,"imaterial,order_algorithm ",imaterial, &
         order_algorithm(imaterial)
       enddo
       do dir=1,sdim
        print *,"dir,xsten0(0) ",dir,xsten0(0,dir)
        print *,"dir,xsten0(2) ",dir,xsten0(2,dir)
        print *,"dir,dx ",dir,xsten0(1,dir)-xsten0(-1,dir)
       enddo
       print *,"imaterial_override ",imaterial_override
       print *,"MOFITERMAX ",MOFITERMAX
      else if (mof_verbose.eq.0) then
       ! do nothing
      else if ((mixed_cell.eq.0).and.(mof_verbose.eq.2)) then
       ! do nothing
      else
       print *,"mof_verbose invalid in multimaterial_MOF"
       print *,"mof_verbose= ",mof_verbose
       print *,"continuous_mof=",continuous_mof
       stop
      endif

      remaining_vfrac=zero
      single_material=0
      num_materials_cell=0

        ! if F<eps for F>1-eps, then moments and vfracs are truncated.
      call make_vfrac_sum_ok(mofdata,mofdata,nmat,ngeom_recon,sdim,6)

       ! clear flag for all nmat materials.
       ! vfrac,centroid,order,slope,intercept x nmat

      do imaterial=1,nmat
       mof_iterations(imaterial)=0
       mof_calls(imaterial)=0
       vofcomp=(imaterial-1)*ngeom_recon+1
       mofdata(vofcomp+sdim+1)=zero  ! order=0
       do dir=1,sdim+1
        mofdata(vofcomp+sdim+1+dir)=zero  ! slope=0  intercept=0
       enddo
       do dir=1,sdim
        multi_centroidA(imaterial,dir)=zero
       enddo

       if (mofdata(vofcomp).gt.one-VOF_CUTOFF) then
        if (single_material.ne.0) then
         print *,"cannot have two materials at once"
         print *,"single_material ",single_material
         print *,"imaterial ",imaterial
         print *,"mofdata ",mofdata(vofcomp)
         stop
        endif
        single_material=imaterial
       else
        remaining_vfrac=remaining_vfrac+mofdata(vofcomp)
       endif

       order_algorithm_in(imaterial)=order_algorithm(imaterial)
       if (order_algorithm(imaterial).eq.0) then
        if (mofdata(vofcomp).ge.one-VOF_CUTOFF) then
         order_algorithm_in(imaterial)=nmat+1
        else if (mofdata(vofcomp).le.VOF_CUTOFF) then
         order_algorithm_in(imaterial)=nmat+1  ! was "1"
        else
         ! do nothing
        endif
       endif

       if (mofdata(vofcomp).ge.VOF_CUTOFF) then
        num_materials_cell=num_materials_cell+1
       endif

      enddo  ! imaterial

      if (num_materials_cell.le.2) then
       do imaterial=1,nmat
        if (order_algorithm(imaterial).eq.0) then
         order_algorithm_in(imaterial)=nmat+1  ! was "1"
        endif
       enddo
      endif

      n_ndef=0
      do imaterial=1,nmat
       placeholder(imaterial)=0 ! 0 if place open, 1 if place taken
       placelist(imaterial)=0  ! list of available places (>=n_ndef)
       flexlist(imaterial)=0   ! list of materials that need an ordering
      enddo
      do imaterial=1,nmat
       if (order_algorithm_in(imaterial).eq.0) then
        n_ndef=n_ndef+1
        flexlist(n_ndef)=imaterial
       else if (order_algorithm_in(imaterial).gt.nmat) then
        ! do nothing
       else if ((order_algorithm_in(imaterial).ge.1).and. &
                (order_algorithm_in(imaterial).le.nmat)) then
        placeholder(order_algorithm_in(imaterial))=1
       else
        print *,"order_algorithm_in(imaterial) invalid"
        stop
       endif
      enddo ! imaterial
        ! placelist: list of available places >= n_ndef
      place_index=0
      do imaterial=1,nmat
       if (placeholder(imaterial).eq.0) then
        place_index=place_index+1
        placelist(place_index)=imaterial
       endif
      enddo 
      if (place_index.lt.n_ndef) then
       print *,"place_index invalid"
       stop
      endif

      if (n_ndef.eq.0) then
       ! do nothing
      else if ((n_ndef.ge.1).and.(n_ndef.le.nmat)) then
       call nfact(n_ndef,n_orderings) 
        ! order_algorithm_in(flexlist(1..n_ndef))=
        !  placelist(order_array(*,1..n_ndef))
       allocate(order_array(n_orderings,n_ndef))
       allocate(order_stack(n_orderings,n_ndef))
       alloc_flag=alloc_flag+2
       order_count=0
       order_stack_count=0
       do irank=1,n_ndef
        temp_order(1)=irank
        do jflex=2,n_ndef
         temp_order(jflex)=0 
        enddo
        call push_order_stack(order_stack,order_stack_count, &
         temp_order,n_orderings,n_ndef,nmat)
       enddo  ! i
       do while (order_stack_count.gt.0)
        call pop_order_stack(order_stack,order_stack_count, &
         temp_order,n_orderings,n_ndef,nmat)
        jflex=n_ndef
        do while (temp_order(jflex).eq.0)
         jflex=jflex-1
        enddo
        if (jflex.eq.n_ndef) then
         order_count=order_count+1
         if (order_count.gt.n_orderings) then
          print *,"order_count too big"
          stop
         endif
         do iflex=1,n_ndef
          order_array(order_count,iflex)=temp_order(iflex)
         enddo
        else if ((jflex.ge.1).and.(jflex.lt.n_ndef)) then
         do irank=1,n_ndef
          is_valid=1
          do kflex=1,jflex
           if (temp_order(kflex).eq.irank) then
            is_valid=0
           endif
          enddo
          if (is_valid.eq.1) then
           temp_order(jflex+1)=irank
           call push_order_stack(order_stack,order_stack_count, &
            temp_order,n_orderings,n_ndef,nmat)
          endif
         enddo ! irank
        else
         print *,"j invalid"
         stop
        endif
       enddo ! while order_stack_count>0
    
       deallocate(order_stack) 
       alloc_flag=alloc_flag-1

       if (order_count.ne.n_orderings) then
        print *,"order_count invalid"
        stop
       endif

      else
       print *,"n_ndef invalid"
       stop
      endif

      if ((single_material.gt.0).and. &
          (remaining_vfrac.lt.VOF_CUTOFF)) then

       vofcomp=(single_material-1)*ngeom_recon+1
       mofdata(vofcomp+sdim+1)=one  ! order=1
       do dir=1,sdim
        nrecon(dir)=zero  
       enddo
       nrecon(1)=one ! null slope=(1 0 0) 
        ! phi = n dot (x-x0) + int
        ! int=-min (n dot (x-x0)) where x is a point in the cell.  
        ! x0 is cell center (xcell)
       mofdata(vofcomp+2*sdim+2)=two*bfact*dxmaxLS  ! null intercept
       do dir=1,sdim
         mofdata(vofcomp+sdim+1+dir)=nrecon(dir)
         multi_centroidA(single_material,dir)=zero  ! cell is full
       enddo 

      else if ((single_material.eq.0).or. &
               (remaining_vfrac.ge.VOF_CUTOFF)) then

          ! no need to pick an optimal ordering
       if (n_ndef.eq.0) then

        if (continuous_mof.eq.0) then
         call Box_volumeFAST(bfact,dx,xsten0,nhalf0,uncaptured_volume_vof, &
          uncaptured_centroid_vof,levelrz,sdim)
         call Box_volumeFAST(bfact,dx,xsten0,nhalf0,uncaptured_volume_cen, &
          uncaptured_centroid_cen,levelrz,sdim)
        else if (continuous_mof.eq.1) then
         call Box_volume_super( &
          bfact,dx,xsten0,nhalf0, &
          uncaptured_volume_vof,uncaptured_centroid_vof, &
          sdim,levelrz)
         call Box_volume_super( &
          bfact,dx,xsten0,nhalf0, &
          uncaptured_volume_cen,uncaptured_centroid_cen, &
          sdim,levelrz)
        else if (continuous_mof.eq.2) then
         call Box_volumeFAST(bfact,dx,xsten0,nhalf0,uncaptured_volume_vof, &
          uncaptured_centroid_vof,levelrz,sdim)
         call Box_volume_super( &
          bfact,dx,xsten0,nhalf0, &
          uncaptured_volume_cen,uncaptured_centroid_cen, &
          sdim,levelrz)
        else
         print *,"continuous_mof invalid"
         stop
        endif

        imaterial_count=1
        do while ((imaterial_count.le.nmat).and. &
                  (uncaptured_volume_vof.gt.zero))
         call individual_MOF( &
          ls_mof, &
          lsnormal, &
          lsnormal_valid, &
          bfact,dx,xsten0,nhalf0, &
          order_algorithm_in, &
          ngeom_recon, &
          xtetlist_vof, &
          xtetlist_cen, &
          nmax, &
          mofdata,levelrz, &
          imaterial_count, &
          uncaptured_volume_vof, &
          uncaptured_volume_cen, &
          multi_centroidA, &
          time, &
          continuous_mof, &
          nmat,nten,sdim)
         imaterial_count=imaterial_count+1
        enddo

         ! multiple orderings must be tested.
       else if ((n_ndef.ge.1).and.(n_ndef.le.nmat)) then

        allocate(mofdata_array(n_orderings,nmat*ngeom_recon))
        allocate(centroidA_array(n_orderings,nmat,sdim))
        allocate(moferror_array(n_orderings))
        alloc_flag=alloc_flag+3

        argmin_order=0
        min_error=0.0

        do order_count=1,n_orderings

         do dir=1,nmat*ngeom_recon
          mofdata_in(dir)=mofdata(dir)
         enddo

         do iflex=1,n_ndef
          imaterial=flexlist(iflex)
          if ((imaterial.lt.1).or.(imaterial.gt.nmat)) then
           print *,"imaterial invalid"
           stop
          endif
          irank=order_array(order_count,iflex)
          if ((irank.lt.1).or.(irank.gt.n_ndef)) then
           print *,"irank invalid"
           stop
          endif
          order_algorithm_in(imaterial)=placelist(irank)
         enddo ! iflex
        
         if (continuous_mof.eq.0) then
          call Box_volumeFAST(bfact,dx,xsten0,nhalf0, &
           uncaptured_volume_vof, &
           uncaptured_centroid_vof,levelrz,sdim)
          call Box_volumeFAST(bfact,dx,xsten0,nhalf0, &
           uncaptured_volume_cen, &
           uncaptured_centroid_cen,levelrz,sdim)
         else if (continuous_mof.eq.1) then
          call Box_volume_super(bfact,dx,xsten0,nhalf0, &
           uncaptured_volume_vof,uncaptured_centroid_vof, &
           sdim,levelrz)
          call Box_volume_super(bfact,dx,xsten0,nhalf0, &
           uncaptured_volume_cen,uncaptured_centroid_cen, &
           sdim,levelrz)
         else if (continuous_mof.eq.2) then
          call Box_volumeFAST(bfact,dx,xsten0,nhalf0, &
           uncaptured_volume_vof, &
           uncaptured_centroid_vof,levelrz,sdim)
          call Box_volume_super(bfact,dx,xsten0,nhalf0, &
           uncaptured_volume_cen,uncaptured_centroid_cen, &
           sdim,levelrz)
         else
          print *,"continuous_mof invalid"
          stop
         endif

         imaterial_count=1
         do while ((imaterial_count.le.nmat).and. &
                   (uncaptured_volume_vof.gt.zero))
          call individual_MOF( &
           ls_mof, &
           lsnormal, &
           lsnormal_valid, &
           bfact,dx,xsten0,nhalf0, &
           order_algorithm_in, &
           ngeom_recon, &
           xtetlist_vof, &
           xtetlist_cen, &
           nmax, &
           mofdata_in,levelrz, &
           imaterial_count, &
           uncaptured_volume_vof, &
           uncaptured_volume_cen, &
           multi_centroidA, &
           time, &
           continuous_mof, &
           nmat,nten,sdim)
          imaterial_count=imaterial_count+1
         enddo ! while not all of uncaptured space filled

         mof_err=zero
         do imaterial = 1,nmat
          vofcomp=(imaterial-1)*ngeom_recon+1
          do dir=1,sdim
           xref_mat(dir)=mofdata_in(vofcomp+dir)
           xact_mat(dir)=multi_centroidA(imaterial,dir)
          enddo
          call RT_transform_offset(xref_mat,uncaptured_centroid_cen,xref_matT)
          call RT_transform_offset(xact_mat,uncaptured_centroid_cen,xact_matT)
           
          do dir=1,sdim
           centroidA_array(order_count,imaterial,dir)= &
             multi_centroidA(imaterial,dir)
           mof_err = mof_err + &
            mofdata_in(vofcomp)*((xref_matT(dir)-xact_matT(dir))**2)
          enddo ! dir
         enddo ! imaterial
         do dir=1,nmat*ngeom_recon
          mofdata_array(order_count,dir)=mofdata_in(dir)
         enddo
         moferror_array(order_count)=mof_err
         if (argmin_order.eq.0) then
          argmin_order=order_count
          min_error=mof_err
         else if (mof_err.lt.min_error) then
          min_error=mof_err
          argmin_order=order_count
         endif
         
         if (1.eq.0) then
          print *,"n_ndef= ",n_ndef
          print *,"order_count=",order_count
          do imaterial=1,nmat
           print *,"imaterial,order_algorithm_in ",imaterial, &
            order_algorithm_in(imaterial)
          enddo
          print *,"mof_err=",mof_err
          print *,"argmin_order=",argmin_order
          print *,"min_error=",min_error
         endif 
        enddo ! order_count

        if ((argmin_order.lt.1).or.(argmin_order.gt.n_orderings)) then
         print *,"argmin_order invalid"
         stop
        else
         do dir=1,nmat*ngeom_recon
          mofdata(dir)=mofdata_array(argmin_order,dir)
         enddo
         do imaterial = 1,nmat
          do dir=1,sdim
           multi_centroidA(imaterial,dir)= &
            centroidA_array(argmin_order,imaterial,dir)
          enddo
         enddo
        endif

        deallocate(mofdata_array)
        deallocate(centroidA_array)
        deallocate(moferror_array)
        alloc_flag=alloc_flag-3

       else
        print *,"n_ndef invalid"
        stop
       endif

      else
       print *,"single_material or remaining_vfrac invalid"
       stop
      endif

      if (n_ndef.eq.0) then
       ! do nothing
      else if ((n_ndef.ge.1).and.(n_ndef.le.nmat)) then
       deallocate(order_array)
       alloc_flag=alloc_flag-1
      else
       print *,"n_ndef invalid"
       stop
      endif

      if (alloc_flag.ne.0) then
       print *,"alloc_flag invalid in MOF.F90"
       stop
      endif

      do imaterial=1,nmat
       vofcomp=(imaterial-1)*ngeom_recon+1
       if (abs(voftest(imaterial)-mofdata(vofcomp)).ge.1.0E-3) then
        print *,"volume fraction changed"
        print *,"imaterial,vofbefore,vofafter ",imaterial, &
          voftest(imaterial),mofdata(vofcomp)
        do imaterial2=1,nmat
         vofcomp=(imaterial2-1)*ngeom_recon+1
         print *,"imaterial2,vofbefore,vofafter ",imaterial2, &
          voftest(imaterial2),mofdata(vofcomp)
        enddo
        stop
       endif
      enddo 


      if ((mof_verbose.eq.1).or. &
          ((mixed_cell.eq.1).and.(mof_verbose.eq.2))) then
       print *,"AFTER AFTER"
       print *,"nmax = ",nmax
       print *,"levelrz = ",levelrz
       print *,"nmat = ",nmat
       print *,"time = ",time
       print *,"sdim = ",sdim
       print *,"continuous_mof = ",continuous_mof
       print *,"ngeom_recon = ",ngeom_recon
       do imaterial=1,nmat*ngeom_recon
        print *,"i,mofdata ",imaterial,mofdata(imaterial)
       enddo
       do imaterial=1,nmat
        print *,"imaterial,order_algorithm ",imaterial, &
         order_algorithm(imaterial)
       enddo
       do dir=1,sdim
        print *,"dir,xsten0(0,dir) ",dir,xsten0(0,dir)
       enddo
       print *,"imaterial_override ",imaterial_override
       print *,"MOFITERMAX ",MOFITERMAX
      else if (mof_verbose.eq.0) then
       ! do nothing
      else if ((mixed_cell.eq.0).and.(mof_verbose.eq.2)) then
       ! do nothing
      else
       print *,"mof_verbose invalid in multimaterial_MOF 2"
       print *,"mof_verbose= ",mof_verbose
       print *,"continuous_mof=",continuous_mof
       stop
      endif

      deallocate(ls_mof)
      deallocate(lsnormal)
      deallocate(lsnormal_valid)

      return
      end subroutine multimaterial_MOF


      subroutine multimaterial_MOF_intercept( &
        bfact,dx,xsten0,nhalf0, &
        use_ls_data, &
        LS_stencil, &
        xtetlist,nmax, &
        mofdata, &
        multi_centroidA, &
        time, &
        levelrz,nmat,nten,sdim)

      use geometry_intersect_module
      use global_utility_module

      IMPLICIT NONE

      INTEGER_T nmat,nten,sdim,nmax,bfact,nhalf0
      INTEGER_T, INTENT (IN) :: use_ls_data
      REAL_T, INTENT (IN) :: LS_stencil(D_DECL(-1:1,-1:1,-1:1),nmat+nten)

      REAL_T xtetlist(sdim+1,sdim,nmax)
      REAL_T xsten0(-nhalf0:nhalf0,sdim)
      REAL_T dx(sdim)
      REAL_T mofdata(nmat*(2*sdim+3))
      REAL_T mofdata_in(nmat*(2*sdim+3))
      REAL_T multi_centroidA(nmat,sdim)
      REAL_T time
      INTEGER_T levelrz,imaterial,vofcomp,dir
      INTEGER_T imaterial_count
      REAL_T uncaptured_volume
      REAL_T uncaptured_centroid(sdim)
      INTEGER_T single_material
      REAL_T remaining_vfrac
      REAL_T nrecon(sdim)
      INTEGER_T ngeom_recon
      INTEGER_T order_algorithm_in(nmat)
      INTEGER_T previous_order(nmat)
      INTEGER_T num_materials_cell
      INTEGER_T current_order
      REAL_T dxmaxLS
      INTEGER_T nten_test,i1,j1,k1,k1lo,k1hi
      INTEGER_T im_main,im_main_opp,iten
      REAL_T voftest(nmat)

      REAL_T, dimension(D_DECL(:,:,:),:), allocatable :: ls_mof
      REAL_T, dimension(:,:), allocatable :: lsnormal
      INTEGER_T, dimension(:), allocatable :: lsnormal_valid

#include "mofdata.H"

      ngeom_recon=2*sdim+3

      if (nhalf0.lt.3) then
       print *,"nhalf0 invalid"
       stop
      endif
      if (bfact.lt.1) then
       print *,"bfact invalid"
       stop
      endif
      if ((sdim.ne.3).and.(sdim.ne.2)) then
       print *,"sdim invalid multimaterial_MOF_intercept"
       stop
      endif
      if ((nmat.lt.1).or.(nmat.gt.MAX_NUM_MATERIALS)) then
       print *,"nmat invalid multimaterial mof intercept"
       print *,"nmat= ",nmat
       stop
      endif
      if (nmax.lt.10) then
       print *,"nmax too small"
       stop
      endif
      nten_test=( (nmat-1)*(nmat-1)+nmat-1 )/2
      if (nten_test.ne.nten) then
       print *,"nten invalid multimofint nten nten_test ",nten,nten_test
       stop
      endif
      if ((use_ls_data.ne.0).and.(use_ls_data.ne.1)) then
       print *,"use_ls_data invalid"
       stop
      endif

      call get_dxmaxLS(dx,bfact,dxmaxLS)

      if (sdim.eq.2) then
       k1lo=0
       k1hi=0
      else if (sdim.eq.3) then
       k1lo=-1
       k1hi=1
      else
       print *,"dimension bust"
       stop
      endif

      do imaterial=1,nmat
       vofcomp=(imaterial-1)*ngeom_recon+1
       voftest(imaterial)=mofdata(vofcomp)
      enddo 

      allocate(ls_mof(D_DECL(-1:1,-1:1,-1:1),nmat+nten))
      allocate(lsnormal(nmat+nten,sdim))
      allocate(lsnormal_valid(nmat+nten))

      if (use_ls_data.eq.1) then
       do i1=-1,1
       do j1=-1,1
       do k1=k1lo,k1hi
       do imaterial=1,nmat+nten
        ls_mof(D_DECL(i1,j1,k1),imaterial)= &
         LS_stencil(D_DECL(i1,j1,k1),imaterial) 
       enddo
       enddo
       enddo
       enddo
       do imaterial=1,nmat+nten

        if ((imaterial.ge.1).and.(imaterial.le.nmat)) then
         im_main=imaterial
         im_main_opp=imaterial
        else if ((imaterial.gt.nmat).and.(imaterial.le.nmat+nten)) then
         iten=imaterial-nmat
         call get_inverse_iten(im_main,im_main_opp,iten,nmat)
        else
         print *,"imaterial invalid"
         stop
        endif 
        if ((voftest(im_main).le.VOF_CUTOFF).or. &
            (voftest(im_main_opp).le.VOF_CUTOFF)) then
         lsnormal_valid(imaterial)=0
        else

         call find_cut_geom_slope_CLSVOF( &
          ls_mof, &
          lsnormal, &
          lsnormal_valid, &
          bfact,dx,xsten0,nhalf0, &
          imaterial, &
          dxmaxLS, &
          levelrz, &
          nmat,nten,sdim)

        endif

       enddo ! imaterial
      else if (use_ls_data.eq.0) then
       do imaterial=1,nmat+nten
        lsnormal_valid(imaterial)=0
       enddo
      else
       print *,"use_ls_data invalid"
       stop
      endif

      remaining_vfrac=zero
      single_material=0
      num_materials_cell=0

        ! if F<eps for F>1-eps, then moments and vfracs are truncated.
      call make_vfrac_sum_ok(mofdata,mofdata,nmat,ngeom_recon,sdim,6)

       ! vfrac,centroid,order,slope,intercept x nmat
       ! initialize ordering using the supercell ordering.
     
      do imaterial=1,nmat
       vofcomp=(imaterial-1)*ngeom_recon+1
       do dir=1,sdim
        multi_centroidA(imaterial,dir)=zero
       enddo

       if (mofdata(vofcomp).gt.one-VOF_CUTOFF) then
        if (single_material.ne.0) then
         print *,"cannot have two materials at once"
         print *,"single_material ",single_material
         print *,"imaterial ",imaterial
         print *,"mofdata ",mofdata(vofcomp)
         stop
        endif
        single_material=imaterial
       else
        remaining_vfrac=remaining_vfrac+mofdata(vofcomp)
       endif

       current_order=NINT(mofdata(vofcomp+sdim+1))
        ! ignore the last interface slope since it is a default slope 
        ! that is the negative of the previous.
       if (current_order.ge.nmat) then
        current_order=0
       endif 
       previous_order(imaterial)=current_order
       mofdata(vofcomp+sdim+1)=zero

       if (current_order.eq.0) then
        if (mofdata(vofcomp).ge.one-VOF_CUTOFF) then
         current_order=nmat+1
        else if (mofdata(vofcomp).le.VOF_CUTOFF) then
         current_order=1
        else
         ! do nothing
        endif
       endif
       order_algorithm_in(imaterial)=current_order

       if (mofdata(vofcomp).ge.VOF_CUTOFF) then
        num_materials_cell=num_materials_cell+1
       endif

      enddo  ! imaterial

      do imaterial=1,nmat
       if (order_algorithm_in(imaterial).eq.0) then
        order_algorithm_in(imaterial)=nmat+1
       endif
      enddo ! imaterial

      if ((single_material.gt.0).and. &
          (remaining_vfrac.lt.VOF_CUTOFF)) then
       vofcomp=(single_material-1)*ngeom_recon+1
       mofdata(vofcomp+sdim+1)=one  ! order=1
       do dir=1,sdim
        nrecon(dir)=zero  
       enddo
       nrecon(1)=one ! null slope=(1 0 0) 
        ! phi = n dot (x-x0) + int
        ! int=-min (n dot (x-x0)) where x is a point in the cell.  
        ! x0 is cell center (xcell)
       mofdata(vofcomp+2*sdim+2)=two*bfact*dxmaxLS  ! null intercept
       do dir=1,sdim
         mofdata(vofcomp+sdim+1+dir)=nrecon(dir)
         multi_centroidA(single_material,dir)=zero  ! cell is full
       enddo 
       do imaterial=1,nmat
        if (imaterial.ne.single_material) then
         vofcomp=(imaterial-1)*ngeom_recon+1
         do dir=1,sdim+1
          mofdata(vofcomp+sdim+1+dir)=zero
         enddo
        endif
       enddo ! imaterial

      else if ((single_material.eq.0).or. &
               (remaining_vfrac.ge.VOF_CUTOFF)) then

       call Box_volumeFAST(bfact,dx,xsten0,nhalf0, &
         uncaptured_volume, &
         uncaptured_centroid,levelrz,sdim)

       imaterial_count=1
       do while ((imaterial_count.le.nmat).and. &
                 (uncaptured_volume.gt.zero))
        call individual_MOF_intercept( &
          ls_mof, &
          lsnormal, &
          lsnormal_valid, &
          bfact,dx,xsten0,nhalf0, &
          order_algorithm_in, &
          previous_order, &
          ngeom_recon, &
          xtetlist,nmax, &
          mofdata,levelrz, &
          imaterial_count, &
          uncaptured_volume,multi_centroidA, &
          time, &
          nmat,nten,sdim)
        imaterial_count=imaterial_count+1
       enddo

      else
       print *,"single_material or remaining_vfrac invalid"
       stop
      endif

      deallocate(ls_mof)
      deallocate(lsnormal)
      deallocate(lsnormal_valid)

      return
      end subroutine multimaterial_MOF_intercept


      subroutine diagnostic_MOF(sdim,nmat,nmax)

      use geometry_intersect_module
      use global_utility_module


      IMPLICIT NONE

      REAL_T time
      INTEGER_T sdim,nmax,nmat
      REAL_T xtetlist(sdim+1,sdim,nmax)
      REAL_T xsten0(-3:3,sdim) 
      REAL_T dx(sdim) 
      REAL_T mofdata(nmat*(2*sdim+3))
      REAL_T multi_centroidA(nmat,sdim)
      INTEGER_T continuous_mof
      INTEGER_T levelrz
      REAL_T angle(sdim-1)
      REAL_T xpoint(sdim)
      INTEGER_T nrecon,nsamples,im,ntry,iangle,vofcomp
      INTEGER_T dir,dir2
      INTEGER_T seed
      REAL_T max_mof_error,moferror
      REAL_T nslope(sdim)
      REAL_T nslope2(sdim)
      REAL_T intercept,intercept2
      INTEGER_T shapeflag
      REAL_T volcut,volcut2,areacut,areacut2
      REAL_T cencut(sdim)
      REAL_T cencut2(sdim)
      REAL_T areacentroid(sdim)
      REAL_T areacentroid2(sdim)
      REAL_T xtet(sdim+1,sdim)
      REAL_T total_volume
      REAL_T total_centroid(sdim)
      REAL_T xref_mat(sdim)
      REAL_T xact_mat(sdim)
      REAL_T xref_matT(sdim)
      REAL_T xact_matT(sdim)
      REAL_T avgiter
      REAL_T total_calls(nmat)
      REAL_T total_iterations(nmat)
      INTEGER_T max_iterations(nmat)
      REAL_T shrink_factor
      INTEGER_T mof_verbose
      INTEGER_T use_ls_data
      INTEGER_T ngeom_recon
      INTEGER_T isten,bfact,nten
      REAL_T, dimension(D_DECL(:,:,:),:), allocatable :: LS_stencil
      INTEGER_T i1,j1,k1,k1lo,k1hi
      INTEGER_T nhalf0

#include "mofdata.H"

      print *,"in diagnostic_MOF"
      print *,"DO NOT RUN THIS TEST WITH MULTIPLE THREADS"

      nhalf0=3
      bfact=2
      shrink_factor=one/three

      ngeom_recon=2*sdim+3
      nrecon=nmat*ngeom_recon
      nten=( (nmat-1)*(nmat-1)+nmat-1 )/2
      allocate(LS_stencil(D_DECL(-1:1,-1:1,-1:1),nmat+nten))
      if (sdim.eq.2) then
       k1lo=0
       k1hi=0
      else if (sdim.eq.3) then
       k1lo=-1
       k1hi=1
      else
       print *,"dimension bust"
       stop
      endif
      do i1=-1,1
      do j1=-1,1
      do k1=k1lo,k1hi
      do im=1,nmat+nten
       LS_stencil(D_DECL(i1,j1,k1),im)=zero
      enddo
      enddo
      enddo
      enddo

      nsamples=90000
      if (nmax.ne.400) then
       print *,"nmax invalid diagnostic MOF"
       stop
      endif
      if (nmat.ne.2) then
       print *,"nmat invalid diagnostic mof"
       print *,"nmat= ",nmat
       stop
      endif
      if ((sdim.ne.2).and.(sdim.ne.3)) then
       print *,"sdim invalid"
       stop
      endif
      continuous_mof=0
      levelrz=0
      do dir=1,sdim
       dx(dir)=one
       do isten=-nhalf0,nhalf0
        xsten0(isten,dir)=isten*half*dx(dir)
       enddo
      enddo ! dir
      do im=1,nmat
       total_calls(im)=zero
       total_iterations(im)=zero
       max_iterations(im)=0
      enddo

      seed=86456
      call srand(seed)

      max_mof_error=zero

      do ntry=1,nsamples

       do iangle=1,sdim-1
#if (1==0)
        angle(iangle)=rand()
#else
        print *,"uncomment rand command"
        stop
#endif
        angle(iangle)=angle(iangle)*two*Pi 
       enddo
       if (sdim.eq.3) then
        call angle_to_slope3D(angle,nslope,sdim)
       else if (sdim.eq.2) then
        call angle_to_slope2D(angle,nslope,sdim)
       else
        print *,"sdim invalid"
        stop
       endif
        ! phi=n dot (x-xpoint)= n dot (x-xcell) +intercept
        ! intercept=n dot (xcell-xpoint)
       intercept=zero
       do dir=1,sdim
        nslope2(dir)=-nslope(dir)

#if (1==0)
        xpoint(dir)=rand()
#else
        print *,"uncomment rand command"
        stop
#endif

        xpoint(dir)=xpoint(dir)-half
        xpoint(dir)=xpoint(dir)*shrink_factor

        intercept=intercept+nslope(dir)*(xsten0(0,dir)-xpoint(dir))
       enddo  ! dir
       intercept2=-intercept
       shapeflag=0
       call fast_cut_cell_intersection( &
        bfact,dx,xsten0,nhalf0, &
        nslope,intercept, &
        volcut,cencut,areacut,areacentroid,levelrz, &
        xsten0,nhalf0,xtet,shapeflag,sdim)
       call fast_cut_cell_intersection( &
        bfact,dx,xsten0,nhalf0, &
        nslope2,intercept2, &
        volcut2,cencut2,areacut2,areacentroid2,levelrz, &
        xsten0,nhalf0,xtet,shapeflag,sdim)

       call Box_volumeFAST( &
         bfact,dx,xsten0,nhalf0, &
         total_volume, &
         total_centroid,levelrz,sdim)

       if (total_volume.lt.zero) then
        print *,"total_volume invalid"
        stop
       endif

       do dir2=1,nrecon
        mofdata(dir2)=zero
       enddo
       im=1
       vofcomp=(im-1)*(2*sdim+3)+1
       mofdata(vofcomp)=volcut/total_volume
       do dir=1,sdim
        mofdata(vofcomp+dir)=cencut(dir)-total_centroid(dir)
       enddo
       im=2
       vofcomp=(im-1)*(2*sdim+3)+1
       mofdata(vofcomp)=volcut2/total_volume
       do dir=1,sdim
        mofdata(vofcomp+dir)=cencut2(dir)-total_centroid(dir)
       enddo

        ! diagnostic_MOF
       time=zero
       mof_verbose=0
       use_ls_data=0

       call multimaterial_MOF( &
        bfact,dx,xsten0,nhalf0, &
        mof_verbose, &
        use_ls_data, &
        LS_stencil, &
        xtetlist, &
        xtetlist, &
        nmax, &
        mofdata,multi_centroidA, &
        time, &
        continuous_mof, &
        levelrz,nmat,nten,sdim, &
        ngeom_recon)

       moferror=zero
       do im=1,nmat
        vofcomp=(im-1)*(2*sdim+3)+1

        do dir=1,sdim
         xref_mat(dir)=mofdata(vofcomp+dir)
         xact_mat(dir)=multi_centroidA(im,dir)
        enddo
        call RT_transform_offset(xref_mat,total_centroid,xref_matT)
        call RT_transform_offset(xact_mat,total_centroid,xact_matT)

        do dir=1,sdim
         moferror=moferror+ &
          mofdata(vofcomp)*((xref_matT(dir)-xact_matT(dir))**2)
        enddo
       enddo  ! im
       moferror=sqrt(moferror)/dx(1)

       if (moferror.gt.max_mof_error) then
        max_mof_error=moferror
       endif

        ! diagnostic_MOF
       do im=1,nmat
        total_calls(im)=total_calls(im)+mof_calls(im)
        total_iterations(im)= &
            total_iterations(im)+mof_iterations(im)
        if (mof_iterations(im).gt.max_iterations(im)) then
         max_iterations(im)=mof_iterations(im)
        endif
       enddo

      enddo ! ntry

      print *,"sdim= ",sdim
      print *,"MOFITERMAX= ",MOFITERMAX
      print *,"nsamples= ",nsamples
      print *,"shrink_factor=",shrink_factor
      do im=1,nmat
       print *,"im= ",im
       print *,"total calls= ",total_calls(im)
       if (total_calls(im).gt.zero) then
        avgiter=total_iterations(im)/total_calls(im)
        print *,"avgiter= ",avgiter
       endif
       print *,"max iterations= ",max_iterations(im)
      enddo
      print *,"max_mof_error=",max_mof_error

      deallocate(LS_stencil)

      return
      end subroutine diagnostic_MOF

      subroutine renormalize_LS(LS_array,im_exclude,nmat)
      use geometry_intersect_module

      IMPLICIT NONE

      INTEGER_T nmat,im_exclude
      REAL_T LS_array(nmat)
      REAL_T LS_hold(nmat)
      INTEGER_T im,im_max,im2

      if ((nmat.lt.1).or.(nmat.gt.MAX_NUM_MATERIALS)) then
       print *,"nmat invalid"
       stop
      endif
      if ((im_exclude.lt.0).or.(im_exclude.gt.nmat)) then
       print *,"im_exclude invalid"
       stop
      endif

      do im=1,nmat 
       LS_hold(im)=LS_array(im)

       if (im.ne.im_exclude) then
        im_max=0
        do im2=1,nmat
         if ((im.ne.im2).and.(im2.ne.im_exclude)) then
          if (im_max.eq.0) then
           im_max=im2
          else if (LS_array(im2).gt.LS_array(im_max)) then
           im_max=im2
          endif
         endif
        enddo  ! im2
        if (im_max.eq.0) then
         ! do nothing
        else if ((im_max.ge.1).and.(im_max.le.nmat)) then
         LS_hold(im)=half*(LS_array(im)-LS_array(im_max))
        else
         print *,"im_max invalid"
         stop
        endif
       endif ! im<>im_exclude
      enddo ! im 

      do im=1,nmat 
       LS_array(im)=LS_hold(im)
      enddo

      return
      end subroutine renormalize_LS


      subroutine project_slopes_to_face( &
       project_status, &
       bfact,dx,xsten0,nhalf0, &
       mofdata,mofdataproject, &
       nmat,sdim,dir,side)

      use geometry_intersect_module

      IMPLICIT NONE

      INTEGER_T nmat,sdim,dir,side,bfact,nhalf0
      INTEGER_T project_status(nmat)
      REAL_T mofdata(nmat*(2*sdim+3))
      REAL_T mofdataproject(nmat*(2*sdim+3))
      REAL_T xsten0(-nhalf0:nhalf0,sdim)
      REAL_T dx(sdim)
      INTEGER_T vofcomp,im
      INTEGER_T ngeom_recon
      INTEGER_T dir2
      REAL_T x0_face(sdim)
      REAL_T xtilde(sdim)
      REAL_T slope(sdim)
      REAL_T slope_project(sdim)
      REAL_T intercept
      REAL_T intercept_save
      REAL_T alpha
      REAL_T mag
      REAL_T test,test_dx

      ngeom_recon=2*sdim+3

      if ((dir.lt.1).or.(dir.gt.sdim)) then
       print *,"dir invalid"
       stop
      endif
      if ((side.ne.1).and.(side.ne.2)) then
       print *,"side invalid"
       stop
      endif

      if (nhalf0.lt.1) then
       print *,"nhalf0 invalid project_slopes_to_face"
       stop
      endif
      if (bfact.lt.1) then
       print *,"bfact invalid"
       stop
      endif
      if ((sdim.ne.3).and.(sdim.ne.2)) then
       print *,"sdim invalid project_slopes_to_face"
       stop
      endif
      if ((nmat.lt.1).or.(nmat.gt.MAX_NUM_MATERIALS)) then
       print *,"nmat invalid project_slopes_to_face"
       stop
      endif

       ! F,X,order,slope,intercept
      do im=1,nmat

       project_status(im)=0

       vofcomp=(im-1)*ngeom_recon+1

       do dir2=1,ngeom_recon
        mofdataproject(vofcomp+dir2-1)=mofdata(vofcomp+dir2-1)
       enddo

       do dir2=1,sdim
        slope(dir2)=mofdata(vofcomp+sdim+1+dir2)
        slope_project(dir2)=slope(dir2)
       enddo
       intercept=mofdata(vofcomp+2*sdim+2)
       intercept_save=intercept
       slope_project(dir)=zero
       mag=zero
       do dir2=1,sdim
        mag=mag+slope_project(dir2)*slope_project(dir2)
       enddo
       mag=sqrt(mag)
       if (mag.gt.zero) then
        do dir2=1,sdim
         slope_project(dir2)=slope_project(dir2)/mag
         x0_face(dir2)=xsten0(0,dir2)
        enddo 
        if (side.eq.1) then
         x0_face(dir)=xsten0(-1,dir)
        else if (side.eq.2) then
         x0_face(dir)=xsten0(1,dir)
        else
         print *,"side invalid"
         stop
        endif
        mag=zero
        do dir2=1,sdim
         mag=mag+slope(dir2)*slope_project(dir2)
        enddo
        if (mag.gt.zero) then
         alpha=-(intercept+slope(dir)*(x0_face(dir)-xsten0(0,dir)))/mag
         intercept=zero
         do dir2=1,sdim
          xtilde(dir2)=x0_face(dir2)+alpha*slope_project(dir2)
          intercept=intercept+ &
           slope_project(dir2)*(xsten0(0,dir2)-xtilde(dir2))
         enddo
         do dir2=1,sdim
          mofdataproject(vofcomp+sdim+1+dir2)=slope_project(dir2)
         enddo
         mofdataproject(vofcomp+2*sdim+2)=intercept
         test_dx=xsten0(1,dir)-xsten0(-1,dir)
         test=intercept_save
         do dir2=1,sdim
          test=test+slope(dir2)*(xtilde(dir2)-xsten0(0,dir2))
         enddo
         if (abs(test).gt.VOF_CUTOFF*test_dx) then
          print *,"test invalid"
          stop
         endif
         project_status(im)=1
        else if (mag.eq.zero) then
         print *,"mag cannot be zero here"
         stop
        else
         print *,"mag invalid"
         stop
        endif
       else if (mag.eq.zero) then
        ! do not modify the slope or intercept
       else
        print *,"mag invalid"
        stop
       endif

      enddo ! im

      return
      end subroutine project_slopes_to_face



! vof,ref centroid,order,slope,intercept  x nmat
      subroutine make_vfrac_sum_ok(mofdata,mofdatavalid,nmat, &
        ngeom_recon,sdim,errorid)
      use geometry_intersect_module

      IMPLICIT NONE

      INTEGER_T nmat,sdim,errorid,ngeom_recon
      REAL_T mofdata(nmat*ngeom_recon)
      REAL_T mofdatavalid(nmat*ngeom_recon)

      INTEGER_T im,dir,vofcomp
      REAL_T voftotal

      if ((nmat.lt.1).or.(nmat.gt.MAX_NUM_MATERIALS)) then
       print *,"nmat bust"
       stop
      endif
      if (ngeom_recon.ne.2*sdim+3) then
       print *,"ngeom_recon invalid"
       stop
      endif
      if ((sdim.ne.2).and.(sdim.ne.3)) then
       print *,"sdim invalid"
       stop
      endif

      voftotal=zero
      do dir=1,nmat*ngeom_recon
       mofdatavalid(dir)=mofdata(dir)
      enddo

      do im=1,nmat
       vofcomp=(im-1)*ngeom_recon+1

       if (abs(mofdata(vofcomp)).le.VOF_CUTOFF) then
        mofdatavalid(vofcomp)=zero
        do dir=1,sdim
         mofdatavalid(vofcomp+dir)=zero
        enddo
       else if (abs(mofdata(vofcomp)-one).le.VOF_CUTOFF) then
        mofdatavalid(vofcomp)=one
        do dir=1,sdim
         mofdatavalid(vofcomp+dir)=zero
        enddo
       else if ((mofdata(vofcomp).gt.zero).and. &
                (mofdata(vofcomp).lt.one)) then
        ! do nothing
       else
        print *,"mofdata(vofcomp) invalid"
        print *,"im,nmat,ngeom_recon,sdim ",im,nmat,ngeom_recon,sdim
        print *,"errorid= ",errorid
        stop
       endif
       voftotal=voftotal+mofdatavalid(vofcomp)
      enddo ! im

      if (voftotal.le.zero) then
       print *,"vacuum bust in make_vfrac_sum ok"
       print *,"errorid= ",errorid
       print *,"nmat= ",nmat
       print *,"sdim= ",sdim
       print *,"voftotal= ",voftotal
       stop
      endif

      do im=1,nmat
       vofcomp=(im-1)*ngeom_recon+1
       mofdatavalid(vofcomp)=mofdatavalid(vofcomp)/voftotal
      enddo 

      return
      end subroutine make_vfrac_sum_ok



 
        ! shapeflag=0 find volumes within xsten_grid
        ! shapeflag=1 find volumes within xtet
        ! multi_cen is "absolute" (not relative to cell center)
      subroutine multi_get_volume_grid( &
       bfact,dx,xsten0,nhalf0, &
       mofdata, &
       xsten_grid,nhalf_grid, &
       xtet, &
       multi_volume,multi_cen, &
       multi_area,rz_flag,xtetlist,nmax, &
       nmat,sdim,shapeflag,caller_id)

      use geometry_intersect_module

      IMPLICIT NONE

      INTEGER_T nmat,sdim,shapeflag,caller_id,bfact,nhalf0,nhalf_grid
      REAL_T xtet(sdim+1,sdim)
      REAL_T mofdata(nmat*(2*sdim+3))
      REAL_T mofdatalocal(nmat*(2*sdim+3))
      REAL_T mofdatasave(nmat*(2*sdim+3))
      REAL_T mofdatavalid(nmat*(2*sdim+3))
      REAL_T xsten0(-nhalf0:nhalf0,sdim)
      REAL_T dx(sdim)
      REAL_T xsten_grid(-nhalf_grid:nhalf_grid,sdim)
      REAL_T multi_volume(nmat)
      REAL_T multi_cen(sdim,nmat)
      REAL_T multi_area(nmat)
      INTEGER_T rz_flag,imaterial,dir,vofcomp,im
      REAL_T uncaptured_volume
      REAL_T uncaptured_centroid(sdim)
      REAL_T volcell
      REAL_T cencell(sdim)
      REAL_T volcut,cencut(sdim)
      INTEGER_T testflag,testflag_save,nlist
      INTEGER_T nmax
      REAL_T xtetlist(sdim+1,sdim,nmax)
      INTEGER_T critical_material
      REAL_T nrecon(sdim)
      REAL_T intercept
      REAL_T voltemp,centemp(sdim),areatemp
      REAL_T areacentroidtemp(sdim)
      INTEGER_T single_material
      REAL_T remaining_vfrac
      REAL_T uncaptured_volume_fraction
      REAL_T uncaptured_volume_save
      INTEGER_T material_used(nmat)
      INTEGER_T im_test,fastflag,ngeom_recon

      ngeom_recon=2*sdim+3

      if (nmax.lt.4) then
       print *, nmax,"nmax invalid multi_get_volume_grid"
       stop
      endif
      if ((nhalf0.lt.1).or.(nhalf_grid.lt.1)) then
       print *,"nhalf invalid multi get volume grid"
       stop
      endif
      if (bfact.lt.1) then
       print *,"bfact invalid"
       stop
      endif
      if ((sdim.ne.3).and.(sdim.ne.2)) then
       print *,"sdim invalid multi_get_volume_grid"
       stop
      endif
      if ((nmat.lt.1).or.(nmat.gt.MAX_NUM_MATERIALS)) then
       print *,"nmat invalid multi get volume grid"
       stop
      endif
 
      do imaterial=1,nmat
       multi_volume(imaterial)=zero
       do dir=1,sdim
        multi_cen(dir,imaterial)=zero
       enddo
       multi_area(imaterial)=zero
      enddo

      call make_vfrac_sum_ok(mofdata,mofdatavalid,nmat,ngeom_recon,sdim,1)

      do dir=1,nmat*ngeom_recon
       mofdatalocal(dir)=mofdatavalid(dir)
       mofdatasave(dir)=mofdatavalid(dir)
      enddo

      call Box_volumeFAST(bfact,dx,xsten0,nhalf0,volcell,cencell,rz_flag,sdim)

      if (shapeflag.eq.0) then
 
       call Box_volumeFAST(bfact,dx,xsten_grid,nhalf_grid, &
         uncaptured_volume, &
         uncaptured_centroid,rz_flag,sdim)

      else if (shapeflag.eq.1) then

       call tetrahedron_volume(xtet,uncaptured_volume, &
         uncaptured_centroid,rz_flag,sdim)

      else
       print *,"shapeflag invalid"
       stop
      endif

      if (volcell.le.zero) then
       print *,"volcell invalid multigetvolume grid"
       stop
      endif
      if (uncaptured_volume.lt.zero) then
       print *,"uncaptured_volume invalid"
       stop
      endif

      if (uncaptured_volume.le.VOF_CUTOFF*volcell) then

       do imaterial=1,nmat
        vofcomp=(imaterial-1)*ngeom_recon+1
        multi_volume(imaterial)=uncaptured_volume*mofdatasave(vofcomp)
        do dir=1,sdim
         multi_cen(dir,imaterial)=uncaptured_centroid(dir)
        enddo
       enddo

      else 

       uncaptured_volume_fraction=one
       do imaterial=1,nmat
        material_used(imaterial)=0
       enddo

       imaterial=1
       do while ((imaterial.le.nmat).and.(uncaptured_volume.gt.zero).and. &
                 (uncaptured_volume_fraction.gt.zero))

        fastflag=1

        remaining_vfrac=zero
        single_material=0

        do im_test=1,nmat
         vofcomp=(im_test-1)*ngeom_recon+1

         if (material_used(im_test).eq.0) then
          if (mofdatasave(vofcomp).gt. &
              uncaptured_volume_fraction-VOF_CUTOFF) then
           if (single_material.ne.0) then
            print *,"cannot have two materials at once"
            print *,"single_material ",single_material
            print *,"im_test ",im_test
            print *,"mofdatavalid ",mofdatavalid(vofcomp)
            print *,"uncaptured_volume_fraction ", &
             uncaptured_volume_fraction
            stop
           endif
           single_material=im_test
          else
           remaining_vfrac=remaining_vfrac+mofdatasave(vofcomp)
          endif
         else if (material_used(im_test).ne.1) then
          print *,"material used bust"
          stop
         endif
        enddo  ! im_test

        if ((single_material.gt.0).and. &
            (remaining_vfrac.lt.VOF_CUTOFF)) then

         vofcomp=(single_material-1)*ngeom_recon+1
         multi_volume(single_material)=uncaptured_volume
         do dir=1,sdim
          multi_cen(dir,single_material)=uncaptured_centroid(dir)
         enddo

         uncaptured_volume=zero
         uncaptured_volume_fraction=zero
         material_used(single_material)=1

        else

         do im=1,nmat
          vofcomp=(im-1)*ngeom_recon+1
          mofdatalocal(vofcomp+sdim+1)=zero
          testflag=NINT(mofdatasave(vofcomp+sdim+1))
          if (testflag.lt.imaterial) then
           mofdatalocal(vofcomp+sdim+1)=testflag
          endif
         enddo

         if ((imaterial.gt.1).and.(imaterial.le.nmat)) then
          fastflag=0
         else if (imaterial.ne.1) then
          print *,"imaterial invalid"
          stop
         endif

         if (fastflag.eq.0) then

          if (shapeflag.eq.0) then ! volumes in a box
             ! only xsten0(0,dir) dir=1..sdim used
             ! in: multi_volume_grid
           call tets_box_planes( &
            bfact,dx,xsten0,nhalf0, &
            xsten_grid,nhalf_grid, &
            mofdatalocal, &
            xtetlist,nlist,nmax,nmat,sdim)

          else if (shapeflag.eq.1) then ! volumes in a tet.

             ! only xsten0(0,dir) dir=1..sdim used
           call tets_tet_planes( &
            bfact,dx,xsten0,nhalf0, &
            xtet,mofdatalocal, &
            xtetlist,nlist,nmax,nmat,sdim)

          else
           print *,"shapeflag invalid"
           stop
          endif

          call get_cut_geom3D(xtetlist,nlist,nmax,rz_flag, &
             volcut,cencut,sdim)

          if (abs(volcut-uncaptured_volume).gt.VOF_CUTOFF*volcell) then
           print *,"volcut invalid multi volume get volume grid "
           stop
          endif

         else if (fastflag.eq.1) then

          ! do nothing

         else 
          print *,"fastflag invalid multi get volume grid"
          stop
         endif

         critical_material=0
         do im=1,nmat
          vofcomp=(im-1)*ngeom_recon+1
          testflag=NINT(mofdatalocal(vofcomp+sdim+1))
          testflag_save=NINT(mofdatasave(vofcomp+sdim+1))
          if ((testflag_save.eq.imaterial).and.(testflag.eq.0).and. &
              (material_used(im).eq.0)) then
           critical_material=im
          endif
         enddo
         
         if (critical_material.gt.0) then        
          vofcomp=(critical_material-1)*ngeom_recon+1
          do dir=1,sdim
           nrecon(dir)=mofdatalocal(vofcomp+sdim+1+dir)
          enddo
          intercept=mofdatalocal(vofcomp+2*sdim+2)

          if (fastflag.eq.0) then
             ! only xsten0(0,dir) dir=1..sdim used
           call multi_cell_intersection( &
            bfact,dx,xsten0,nhalf0, &
            nrecon,intercept, &
            voltemp,centemp,areatemp,rz_flag, &
            xtetlist,nlist,nmax,sdim) 
          else if (fastflag.eq.1) then
             ! only xsten0(0,dir) dir=1..sdim used
           call fast_cut_cell_intersection( &
            bfact,dx,xsten0,nhalf0, &
            nrecon,intercept, &
            voltemp,centemp,areatemp, &
            areacentroidtemp,rz_flag, &
            xsten_grid,nhalf_grid,xtet,shapeflag,sdim) 
          else 
           print *,"fastflag invalid multi get volume grid 2"
           stop
          endif

          multi_volume(critical_material)=voltemp
          do dir=1,sdim
           if (voltemp.gt.zero) then
            multi_cen(dir,critical_material)=centemp(dir)
           else
            multi_cen(dir,critical_material)=zero
           endif
          enddo
          multi_area(critical_material)=areatemp

          uncaptured_volume_save=uncaptured_volume
          uncaptured_volume=uncaptured_volume-voltemp
          if (uncaptured_volume.lt.VOF_CUTOFF*volcell) then
           uncaptured_volume=zero
          endif

             ! V^{uncapt,k}=V+V^{uncapt,k+1}
             ! V^{uncapt,k}x^{uncapt,k}=V x+V^{uncapt,k+1}x^{uncapt,k+1}

          do dir=1,sdim
           if (uncaptured_volume.le.zero) then
            uncaptured_centroid(dir)=zero
           else
            uncaptured_centroid(dir)= &
              (uncaptured_volume_save*uncaptured_centroid(dir)- &
               voltemp*centemp(dir))/uncaptured_volume
           endif
          enddo
  
          uncaptured_volume_fraction=uncaptured_volume_fraction- &
            mofdatalocal(vofcomp)
          if (uncaptured_volume_fraction.lt.VOF_CUTOFF) then
           uncaptured_volume_fraction=zero
          endif
          material_used(critical_material)=1
         endif ! critical_material>0

        endif ! general case

        imaterial=imaterial+1
       enddo  ! while imaterial<=nmat and uncaptured_volume>0 

       if (uncaptured_volume.gt.UNCAPT_TOL*volcell) then
         print *,"not all volume accounted for multi get volume"
         print *,"uncaptured_volume ",uncaptured_volume
         print *,"volcell ",volcell
         print *,"fraction of uncapt volume ",uncaptured_volume/volcell
         print *,"tolerance: ",UNCAPT_TOL
         stop
       endif

      endif

      return
      end subroutine multi_get_volume_grid


        ! multi_cen is "absolute" (not relative to cell center)
      subroutine multi_get_volume_grid_default( &
       bfact,dxgrid,xsten_grid,nhalf_grid, &
       im, &
       multi_volume,multi_cen,rz_flag, &
       nmat,sdim)

      use geometry_intersect_module

      IMPLICIT NONE

      INTEGER_T nmat,sdim,bfact,nhalf_grid
      REAL_T xsten_grid(-nhalf_grid:nhalf_grid,sdim)
      REAL_T dxgrid(sdim)
      REAL_T multi_volume(nmat)
      REAL_T multi_cen(sdim,nmat)
      INTEGER_T rz_flag,im,imloop,dir
      REAL_T volcell
      REAL_T cencell(sdim)

      if (nhalf_grid.lt.1) then
       print *,"nhalf_grid invalid"
       stop
      endif
      if (bfact.lt.1) then
       print *,"bfact invalid"
       stop
      endif
      if ((sdim.ne.3).and.(sdim.ne.2)) then
       print *,"sdim invalid multi_get_volume_grid default"
       stop
      endif
      if ((nmat.lt.1).or.(nmat.gt.MAX_NUM_MATERIALS)) then
       print *,"nmat invalid multi get volume grid default"
       stop
      endif
      if ((im.lt.1).or.(im.gt.nmat)) then
       print *,"im invalid"
       stop
      endif
 
      do imloop=1,nmat
       multi_volume(imloop)=zero
       do dir=1,sdim
        multi_cen(dir,imloop)=zero
       enddo
      enddo

      call Box_volumeFAST(bfact,dxgrid,xsten_grid,nhalf_grid, &
       volcell,cencell,rz_flag,sdim)

      if (volcell.le.zero) then
       print *,"volcell invalid multigetvolume grid default"
       do dir=1,sdim
        print *,"dir,xsten_grid(0),dxgrid ", &
          dir,xsten_grid(0,dir),dxgrid(dir)
       enddo
       stop
      endif

      multi_volume(im)=volcell
      do dir=1,sdim
       multi_cen(dir,im)=cencell(dir)
      enddo

      return
      end subroutine multi_get_volume_grid_default


      subroutine compare_distance( &
       bfact,dx,xsten0,nhalf0, &
       nmat,nten,xaccept,xdonate, &
       newLS,newSLOPE,im_test,n_im, &
       slope,imslope, &
       imcell,sdim,center_stencil,im_solid, &
       donateflag, &
       truncate_volume_fractions, &
       truncate_at_junction)
      use geometry_intersect_module
      use global_utility_module
      IMPLICIT NONE

      INTEGER_T n_im,bfact,nhalf0
      INTEGER_T im_test(n_im)
      INTEGER_T imslope,imcell
      INTEGER_T nmat,nten
      INTEGER_T sdim,center_stencil,im_solid
      INTEGER_T donateflag(nmat+2)
      INTEGER_T truncate_volume_fractions(nmat)
      INTEGER_T truncate_at_junction
      REAL_T xaccept(sdim) 
      REAL_T xdonate(sdim) 
      REAL_T xsten0(-nhalf0:nhalf0,sdim) 
      REAL_T dx(sdim) 
      REAL_T slope(sdim) 
      REAL_T newLS(nmat+nten) 
      REAL_T newSLOPE(sdim*(nmat+nten)) 
      INTEGER_T nten_test,distzero,im_here,im_opp_here
      INTEGER_T im,im_opp,im3,iten,dir,nc
      REAL_T disttest,dist_compare,LSSIGN
      REAL_T slopetest(sdim)

      if (nhalf0.lt.1) then
       print *,"nhalf0 invalid"
       stop
      endif
      if (bfact.lt.1) then
       print *,"bfact invalid"
       stop
      endif
      nten_test=( (nmat-1)*(nmat-1)+nmat-1 )/2
      if (nten_test.ne.nten) then
       print *,"nten invalid multidist nten nten_test ",nten,nten_test
       stop
      endif
      if ((sdim.ne.3).and.(sdim.ne.2)) then
       print *,"sdim invalid compare_distance"
       stop
      endif
      if ((nmat.lt.1).or.(nmat.gt.MAX_NUM_MATERIALS)) then
       print *,"nmat invalid multi_get_distance"
       stop
      endif
      if ((n_im.lt.1).or.(n_im.gt.6)) then
       print *,"n_im invalid"
       stop
      endif
      do nc=1,n_im
       if ((im_test(nc).lt.1).or.(im_test(nc).gt.nmat)) then
        print *,"im_test invalid"
        stop
       endif
      enddo
      if ((imcell.lt.1).or.(imcell.gt.nmat)) then
       print *,"imcell invalid"
       stop
      endif
      if ((imslope.lt.0).or.(imslope.gt.nmat)) then
       print *,"imslope invalid"
       stop
      endif
      if ((im_solid.lt.0).or.(im_solid.gt.nmat)) then
       print *,"im_solid invalid"
       stop
      endif
      if ((center_stencil.ne.0).and.(center_stencil.ne.1)) then
       print *,"center_stencil invalid"
       stop
      endif
      if ((truncate_at_junction.ne.0).and. &
          (truncate_at_junction.ne.1)) then
       print *,"truncate_at_junction invalid"
       stop
      endif

      disttest=zero
      do dir=1,sdim
       disttest=disttest+(xaccept(dir)-xdonate(dir))**2
       slopetest(dir)=xaccept(dir)-xdonate(dir)
      enddo
      disttest=sqrt(disttest)
      distzero=0
      if (disttest.lt.VOF_CUTOFF*dx(1)) then
       distzero=1
       if (center_stencil.eq.1) then
        do dir=1,sdim
         slopetest(dir)=slope(dir)
        enddo
       else if (center_stencil.eq.0) then
        print *,"center_stencil invalid when disttest=0"
        stop
       else
        print *,"center_stencil invalid"
        stop
       endif
      else if (disttest.gt.zero) then
       ! xclosest=x -  d n
       ! n will point from xclosest to x.
       do dir=1,sdim
        slopetest(dir)=slopetest(dir)/disttest
       enddo
      else
       print *,"disttest invalid"
       stop
      endif

      do im=1,nmat
       if ((truncate_volume_fractions(im).ne.0).and. &
           (truncate_volume_fractions(im).ne.1)) then
        print *,"truncate_volume_fractions invalid"
        stop
       endif
      enddo
 
      do im=1,nmat
       if ((donateflag(im).ne.0).and. &
           (donateflag(im).ne.1)) then
        print *,"donateflag invalid"
        stop
       endif

       if (im.ne.im_solid) then

        im_here=0
        im_opp_here=0
        do nc=1,n_im
         im3=im_test(nc)
          
         if ((donateflag(im3).eq.1).or. &
             (im3.eq.im_solid).or. &
             (truncate_at_junction.eq.0)) then
          if (im3.eq.im) then
           im_here=1
          endif
          if (im3.ne.im) then
           im_opp_here=1
          endif
         endif
        enddo  ! nc

        dist_compare=abs(newLS(im))
 
         ! disttest=|xaccept-xdonate|
        if (dist_compare.gt.disttest) then
         if (center_stencil.eq.1) then
          LSSIGN=zero
          if ((imcell.eq.im).and.(im_opp_here.eq.1)) then
           LSSIGN=one
          endif
          if ((imcell.ne.im).and.(im_here.eq.1)) then
           LSSIGN=-one
          endif
          if (LSSIGN.ne.zero) then
           newLS(im)=LSSIGN*disttest  
           if (distzero.eq.1) then
            if ((imslope.ge.1).and.(imslope.le.nmat)) then
             if (imslope.eq.im) then
              LSSIGN=one
             else
              LSSIGN=-one
             endif
            else
             print *,"imslope invalid imslope= ",imslope
             print *,"nmat,nten ",nmat,nten
             print *,"xaccept ",xaccept(1),xaccept(2),xaccept(sdim) 
             print *,"xdonate ",xdonate(1),xdonate(2),xdonate(sdim) 
             print *,"xsten0 ",xsten0(0,1),xsten0(0,2),xsten0(0,sdim) 
             print *,"dx ",dx(1),dx(2),dx(sdim) 
             do nc=1,n_im
              print *,"nc,im_test ",nc,im_test(nc)
             enddo
             print *,"imcell,center_stencil ",imcell,center_stencil
             stop
            endif
           else if (distzero.eq.0) then
            ! do nothing
           else
            print *,"distzero invalid"
            stop
           endif
           do dir=1,sdim
            newSLOPE(sdim*(im-1)+dir)=LSSIGN*slopetest(dir)
           enddo
          endif
         else if (center_stencil.eq.0) then
          if (distzero.ne.0) then
           print *,"distzero invalid"
           stop
          endif
          LSSIGN=zero
          if (newLS(im).ge.zero) then
           im_here=1
          else
           im_opp_here=1
          endif
          if ((newLS(im).ge.zero).and.(im_opp_here.eq.1)) then 
           LSSIGN=one
          endif
          if ((newLS(im).lt.zero).and.(im_here.eq.1)) then
           LSSIGN=-one
          endif
          if (LSSIGN.ne.zero) then
           newLS(im)=LSSIGN*disttest  
           do dir=1,sdim
            newSLOPE(sdim*(im-1)+dir)=LSSIGN*slopetest(dir)
           enddo
          endif
         else
          print *,"center_stencil invalid"
          stop
         endif
        endif   ! dist_compare>disttest? 

        do im_opp=im+1,nmat
         if (im_opp.ne.im_solid) then
          call get_iten(im,im_opp,iten,nmat)
          dist_compare=abs(newLS(nmat+iten))
          if (dist_compare.gt.disttest) then

           im_here=0
           im_opp_here=0
           do nc=1,n_im
            im3=im_test(nc) 
            if ((donateflag(im3).eq.1).or. &
                (im3.eq.im_solid).or. &
                (truncate_at_junction.eq.0)) then
             if (im3.eq.im) then
              im_here=1
             endif
             if (im3.eq.im_opp) then
              im_opp_here=1
             endif
            endif
           enddo  ! nc

           if (center_stencil.eq.1) then
            if ((imcell.eq.im).and.(im_opp_here.eq.1)) then
             LSSIGN=one
            else if ((imcell.eq.im_opp).and.(im_here.eq.1)) then
             LSSIGN=-one
            else if ((im_here.eq.1).and.(im_opp_here.eq.1)) then
             LSSIGN=one
            else
             LSSIGN=zero
            endif
            if (LSSIGN.ne.zero) then
             newLS(nmat+iten)=LSSIGN*disttest
             if (distzero.eq.1) then
              if ((imslope.ge.1).and.(imslope.le.nmat)) then
               if (imslope.eq.im) then
                LSSIGN=one
               else
                LSSIGN=-one
               endif
              else
               print *,"imslope invalid2 imslope= ",imslope
               print *,"nmat,nten ",nmat,nten
               print *,"xaccept ",xaccept(1),xaccept(2),xaccept(sdim)
               print *,"xdonate ",xdonate(1),xdonate(2),xdonate(sdim)
               print *,"xsten0 ",xsten0(0,1),xsten0(0,2),xsten0(0,sdim)
               print *,"dx ",dx(1),dx(2),dx(sdim)
               do nc=1,n_im
                print *,"nc,im_test ",nc,im_test(nc)
               enddo
               print *,"imcell,center_stencil ",imcell,center_stencil
               stop
              endif
             else if (distzero.eq.0) then
              ! do nothing
             else
              print *,"distzero invalid"
              stop
             endif
             do dir=1,sdim
              newSLOPE(sdim*(nmat+iten-1)+dir)=LSSIGN*slopetest(dir)
             enddo
            endif  ! LSSIGN<>0
           else if (center_stencil.eq.0) then
            if (distzero.ne.0) then
             print *,"distzero invalid"
             stop
            endif
            if (newLS(im).ge.zero) then
             im_here=1
            endif
            if (newLS(im_opp).ge.zero) then
             im_opp_here=1
            endif
            if ((newLS(im).ge.zero).and.(im_opp_here.eq.1)) then
             LSSIGN=one
            else if ((newLS(im_opp).ge.zero).and.(im_here.eq.1)) then
             LSSIGN=-one
            else if ((im_here.eq.1).and.(im_opp_here.eq.1)) then
             LSSIGN=one
            else
             LSSIGN=zero
            endif
            if (LSSIGN.ne.zero) then
             newLS(nmat+iten)=LSSIGN*disttest  
             do dir=1,sdim
              newSLOPE(sdim*(nmat+iten-1)+dir)=LSSIGN*slopetest(dir)
             enddo
            endif
           else
            print *,"center_stencil invalid"
            stop
           endif
          endif   ! dist_compare>disttest? 
         endif ! im_opp.ne.im_solid
        enddo ! im_opp
       endif ! im.ne.im_solid ?
      enddo ! im
  
      return
      end subroutine compare_distance 

! find the closest point (x_cp) on the intersection of two planes (in
! 3D) or two lines (2D) to a point (x_a)
! NOTE: Assumed the vectors slope_list are "unit" normals
      subroutine closestINT( &
         bfact,dx,xsten0,nhalf0, &
         x_cp,xplus,xminus, &
         xplus2,xminus2, &
         x_a,maxdx,ilist,jlist,nlist, &
         slope_list,intercept_list,im_list,rz_flag,sdim, &
         inboxflag,nmat)
       use global_utility_module

       IMPLICIT NONE

       INTEGER_T nmat,sdim,rz_flag,nlist,bfact,nhalf0
       REAL_T maxdx
       REAL_T x_cp(sdim)
       REAL_T xplus(sdim)
       REAL_T xminus(sdim)
       REAL_T xplus2(sdim)
       REAL_T xminus2(sdim)
       REAL_T x_a(sdim)
       INTEGER_T ilist, jlist
       REAL_T slope_list(nmat,sdim)
       REAL_T intercept_list(nmat)
       INTEGER_T im_list(nmat)
       REAL_T xsten0(-nhalf0:nhalf0,sdim)
       REAL_T dx(sdim)
       INTEGER_T inboxflag

       INTEGER_T dir,islope
       REAL_T c_i,c_j,s_c
       REAL_T n_i_d_n_j
       REAL_T n_i_c_n_j(sdim)
       REAL_T dtrmn,mag
       REAL_T h(2)
       INTEGER_T indexlist(2)

       if (nhalf0.lt.1) then
        print *,"nhalf0 invalid"
        stop
       endif
       if (bfact.lt.1) then
        print *,"bfact invalid"
        stop
       endif
       if ((nlist.lt.2).or.(nlist.gt.nmat-1).or. &
           (ilist.lt.1).or.(ilist.gt.nlist).or. &
           (jlist.lt.1).or.(jlist.gt.nlist).or. &
           (ilist.eq.jlist)) then
        print *,"ilist,jlist, or nlist invalid"
        stop
       endif

       inboxflag = 0
       indexlist(1)=ilist
       indexlist(2)=jlist

        ! P_i and P_j are in the form n.(x-x_0)+d=0
        ! transform them to the form n.x=h
        ! h=n.x_0-d
       do islope=1,2
        h(islope)=zero
        do dir=1,sdim
         h(islope)=h(islope)+slope_list(indexlist(islope),dir)*xsten0(0,dir)
        enddo
        h(islope)=h(islope)-intercept_list(indexlist(islope))
       enddo ! islope

       if (sdim.eq.2) then
         ! solve the 2-by-2 linear system Ax=b
         ! A = [n_i1 n_i2]
         !     [n_j1 n_j2]
         ! x = [x_cp_1]
         !     [x_cp_2]
         ! h = [n_i1*x_01+n_i2*x_02-d_i]
         !     [n_j1*x_01+n_j2*x_02-d_j]

         dtrmn = slope_list(ilist,1)*slope_list(jlist,2)- &
                 slope_list(ilist,2)*slope_list(jlist,1)

         if (dtrmn.eq.zero) then
          inboxflag=0
         else
          x_cp(1)= &
           (slope_list(jlist,2)*h(1)- &
            slope_list(ilist,2)*h(2))/dtrmn
          x_cp(2)= &
           (-slope_list(jlist,1)*h(1)+ &
             slope_list(ilist,1)*h(2))/dtrmn
           ! in: closestINT
          call check_inbox(x_cp,xsten0,nhalf0,inboxflag)
         endif

       else if (sdim.eq.3) then

        ! The intersction line I is 
        ! I:{x|x=c_i n_i + c_j n_j + c_k (n_i cross n_j)} 
        ! c_i = (h_i-h_j(n_i.n_j))/(1-(n_i.n_j)^2)
        ! c_j = (h_j-h_i(n_i.n_j))/(1-(n_i.n_j)^2)

        n_i_d_n_j=zero
        do dir=1,sdim
         n_i_d_n_j=n_i_d_n_j+slope_list(ilist,dir)*slope_list(jlist,dir)
        enddo

        dtrmn=one-n_i_d_n_j*n_i_d_n_j
        if (dtrmn.le.zero) then
         inboxflag=0
        else
         c_i =(h(1)-h(2)*n_i_d_n_j)/dtrmn
         c_j =(h(2)-h(1)*n_i_d_n_j)/dtrmn

         n_i_c_n_j(1)=slope_list(ilist,2)*slope_list(jlist,3) &
                     -slope_list(ilist,3)*slope_list(jlist,2)
         n_i_c_n_j(2)=slope_list(ilist,3)*slope_list(jlist,1) &
                     -slope_list(ilist,1)*slope_list(jlist,3)
         n_i_c_n_j(3)=slope_list(ilist,1)*slope_list(jlist,2) &
                     -slope_list(ilist,2)*slope_list(jlist,1)


         ! The x_a in the {n_i,n_j,n_i cross n_j} basis form is
         ! x_a = s_i n_i + s_j n_j + s_c (n_i cross n_j)
         ! s_c = x_a.(n_i cross n_j)
         ! The x_cp would be
         ! x_cp = c_i n_i + c_j n_j + s_c (n_i cross n_j)

         mag=zero
         do dir=1,sdim
          mag=mag+n_i_c_n_j(dir)**2
         enddo
         mag=sqrt(mag)

         if (mag.gt.zero) then
          do dir=1,sdim
           n_i_c_n_j(dir)=n_i_c_n_j(dir)/mag
          enddo
 
          s_c=zero
          do dir=1,sdim
           s_c=s_c+x_a(dir)*n_i_c_n_j(dir)
          enddo

          do dir=1,sdim
           x_cp(dir)=c_i*slope_list(ilist,dir)+&
                     c_j*slope_list(jlist,dir)+&
                     s_c*n_i_c_n_j(dir)
          enddo
           ! in: closestINT
          call check_inbox(x_cp,xsten0,nhalf0,inboxflag)
         else
          inboxflag=0
         endif
        endif

       else
        print *,"MOF.F90::closestINT - sdim invalid!"
        stop
       endif ! sdim

       if (inboxflag.eq.0) then
        ! do nothing
       else if (inboxflag.eq.1) then
        do dir=1,sdim
         xplus(dir)=x_cp(dir)-maxdx*LSTHICK*half* &
           (slope_list(ilist,dir)+slope_list(jlist,dir)) 
         xminus(dir)=x_cp(dir)+maxdx*LSTHICK*half* &
           (slope_list(ilist,dir)+slope_list(jlist,dir)) 
         xplus2(dir)=x_cp(dir)-maxdx*LSTHICK*half* &
           (slope_list(ilist,dir)-slope_list(jlist,dir)) 
         xminus2(dir)=x_cp(dir)+maxdx*LSTHICK*half* &
           (slope_list(ilist,dir)-slope_list(jlist,dir)) 
        enddo
       else
        print *,"inboxflag invalid"
        stop
       endif

      end subroutine closestINT


      subroutine closestINT3( &
        bfact,dx,xsten0,nhalf0, &
        x_cp,xplus,xminus, &
        xplus2,xminus2, &
        x_a,maxdx,ilist,jlist,nlist, &
        slope_list,intercept_list,im_list, &
        gphi,&
        x_0side,rz_flag,sdim,inboxflag,nmat)
       use global_utility_module

       IMPLICIT NONE

       INTEGER_T nmat,sdim,rz_flag,nlist,bfact,nhalf0
       REAL_T maxdx
       REAL_T x_cp(sdim)
       REAL_T xplus(sdim)
       REAL_T xminus(sdim)
       REAL_T xplus2(sdim)
       REAL_T xminus2(sdim)
       REAL_T x_a(sdim)
       REAL_T gphi(sdim)
       REAL_T x_0side(sdim)
       INTEGER_T ilist, jlist
       REAL_T slope_list(nmat,sdim)
       REAL_T intercept_list(nmat)
       INTEGER_T im_list(nmat)
       REAL_T xsten0(-nhalf0:nhalf0,sdim)
       REAL_T dx(sdim)
       INTEGER_T inboxflag

       INTEGER_T dir,islope,dir_2,dotpr
       REAL_T c_i,c_j,s_c,c_k
       REAL_T n_i_d_n_j
       REAL_T n_i_c_n_j(sdim)
       REAL_T dtrmn,mag
       REAL_T h(2)
       INTEGER_T cp_init
       INTEGER_T indexlist(2)

       if (bfact.lt.1) then
        print *,"bfact invalid"
        stop
       endif

       if ((nlist.lt.2).or.(nlist.gt.nmat-1).or. &
           (ilist.lt.1).or.(ilist.gt.nlist).or. &
           (jlist.lt.1).or.(jlist.gt.nlist).or. &
           (ilist.eq.jlist)) then
        print *,"ilist,jlist, or nlist invalid"
        stop
       endif

       if (sdim.eq.3) then
        ! do nothing
       else
        print *,"MOF.F90::closestINT3 - sdim invalid!"
        stop
       endif

       inboxflag = 0
       indexlist(1)=ilist
       indexlist(2)=jlist

        ! P_i and P_j are in the form n.(x-x_0)+d=0
        ! transform them to the form n.x=h
        ! h=n.x_0-d
       do islope=1,2
        h(islope)=zero
        do dir=1,sdim
         h(islope)=h(islope)+slope_list(indexlist(islope),dir)*xsten0(0,dir)
        enddo
        h(islope)=h(islope)-intercept_list(indexlist(islope))
       enddo ! islope

       ! The intersction line I is 
       ! I:{x|x=c_i n_i + c_j n_j + c_k (n_i cross n_j)} 
       ! c_i = (h_i-h_j(n_i.n_j))/(1-(n_i.n_j)^2)
       ! c_j = (h_j-h_i(n_i.n_j))/(1-(n_i.n_j)^2)

       n_i_d_n_j=zero
       do dir=1,sdim
        n_i_d_n_j=n_i_d_n_j+slope_list(ilist,dir)*slope_list(jlist,dir)
       enddo

       dtrmn=one-n_i_d_n_j*n_i_d_n_j
       if (dtrmn.le.zero) then
        inboxflag=0
       else
        c_i =(h(1)-h(2)*n_i_d_n_j)/dtrmn
        c_j =(h(2)-h(1)*n_i_d_n_j)/dtrmn

        n_i_c_n_j(1)=slope_list(ilist,2)*slope_list(jlist,3) &
                    -slope_list(ilist,3)*slope_list(jlist,2)
        n_i_c_n_j(2)=slope_list(ilist,3)*slope_list(jlist,1) &
                    -slope_list(ilist,1)*slope_list(jlist,3)
        n_i_c_n_j(3)=slope_list(ilist,1)*slope_list(jlist,2) &
                    -slope_list(ilist,2)*slope_list(jlist,1)

        mag=zero
        do dir=1,sdim
         mag=mag+n_i_c_n_j(dir)**2
        enddo
        mag=sqrt(mag)

        if (mag.gt.zero) then
         do dir=1,sdim
          n_i_c_n_j(dir)=n_i_c_n_j(dir)/mag
         enddo

         ! Find the intersection of line I with a cell face 
         ! Solve for c_k
         ! I:{x|x=c_i n_i + c_j n_j + c_k (n_i cross n_j)} 
         !                 and
         ! x=x_0side(1) or y=x_0side(2) or z=x_0side(3) 
         !        {based on the direction of gphi}

         cp_init=0
         do dir=1,sdim
          if (gphi(dir).eq.zero) then
           ! do nothing
          else if (gphi(dir).eq.one) then
           if (cp_init.eq.1) then
            print *,"cp already init"
            stop
           else if (cp_init.eq.0) then
            cp_init=1
           else
            print *,"cp_init bust"
            stop
           endif
            ! Face is in this direction
           if(n_i_c_n_j(dir).ne.zero) then
            ! Line I is NOT parallel to the face
            c_k = (x_0side(dir)-c_i*slope_list(ilist,dir)&
                  -c_j*slope_list(jlist,dir))/n_i_c_n_j(dir)

            do dir_2=1,sdim
             x_cp(dir_2)=c_i*slope_list(ilist,dir_2)+&
                         c_j*slope_list(jlist,dir_2)+&
                         c_k*n_i_c_n_j(dir_2)
            enddo
             ! in: closestINT3
            call check_inbox(x_cp,xsten0,nhalf0,inboxflag)
           else
            inboxflag=0
           endif
          else
           print *,"gphi invalid"
           stop
          endif
         enddo !dir
         if (cp_init.ne.1) then
          print *,"face plane incorrectly specified"
          stop
         endif
        else
         inboxflag=0
        endif
       endif ! dtrmn<>0 ?

       if (inboxflag.eq.0) then
        ! do nothing
       else if (inboxflag.eq.1) then
        do dir=1,sdim
         xplus(dir)=x_cp(dir)-maxdx*LSTHICK*half* &
           (slope_list(ilist,dir)+slope_list(jlist,dir)) 
         xminus(dir)=x_cp(dir)+maxdx*LSTHICK*half* &
           (slope_list(ilist,dir)+slope_list(jlist,dir)) 
         xplus2(dir)=x_cp(dir)-maxdx*LSTHICK*half* &
           (slope_list(ilist,dir)-slope_list(jlist,dir)) 
         xminus2(dir)=x_cp(dir)+maxdx*LSTHICK*half* &
           (slope_list(ilist,dir)-slope_list(jlist,dir)) 
        enddo
       else
        print *,"inboxflag invalid"
        stop
       endif

       end subroutine closestINT3

         ! n dot (x-x0) + intercept = 0
         ! n dot x = d    d=n dot x0 - intercept
         ! P: ( x= a1 t1 + a2 t2 + a3 n )  a3=d/(n dot n)
         ! xaccept=b1 t1 + b2 t2 + b3 n  b3=xaccept dot n/(n dot n)
         ! xcp=b1 t1 + b2 t2 + a3 n=xaccept+(a3-b3)n
         ! if n dot n=1,
         ! xcp=xaccept+(d-xaccept dot n)n=
         !     xaccept+(n dot x0-int-xaccept dot n)n=
         !     xaccept-(n dot (xaccept-x0)+int)n
         ! xsten0 can be different from xstenbox since the LS
         ! might be projected to a lower dimension.
       subroutine closestPLANE(bfact,dx,xsten0,nhalf0, &
         xcp,xcp_plus,xcp_minus,xaccept, &
         maxdx,slope_in,intercept_in, &
         xstenbox,nhalfbox,sdim,inboxflag)
       use global_utility_module
       IMPLICIT NONE

       INTEGER_T sdim,inboxflag,dir,bfact,nhalf0,nhalfbox
       REAL_T xsten0(-nhalf0:nhalf0,sdim)
       REAL_T xstenbox(-nhalfbox:nhalfbox,sdim)
       REAL_T dx(sdim)
       REAL_T xcp(sdim)
       REAL_T xcp_plus(sdim)
       REAL_T xcp_minus(sdim)
       REAL_T xaccept(sdim)
       REAL_T maxdx
       REAL_T slope_in(sdim)
       REAL_T slope(sdim)
       REAL_T intercept_in
       REAL_T intercept,dist

       if (nhalf0.lt.1) then
        print *,"nhalf0 invalid"
        stop
       endif
       if (nhalfbox.lt.1) then
        print *,"nhalfbox invalid"
        stop
       endif
       if ((sdim.ne.2).and.(sdim.ne.3)) then
        print *,"sdim invalid"
        stop
       endif

       inboxflag=1
       dist=zero
       do dir=1,sdim
        dist=dist+slope_in(dir)**2
       enddo
       dist=sqrt(dist)
       if (dist.gt.zero) then
        do dir=1,sdim
         slope(dir)=slope_in(dir)/dist
        enddo
        intercept=intercept_in/dist
        dist=intercept
        do dir=1,sdim
         dist=dist+slope(dir)*(xaccept(dir)-xsten0(0,dir))
        enddo
        do dir=1,sdim
         xcp(dir)=xaccept(dir)-dist*slope(dir)
         xcp_plus(dir)=xcp(dir)-maxdx*LSTHICK*slope(dir)
         xcp_minus(dir)=xcp(dir)+maxdx*LSTHICK*slope(dir)
        enddo
         ! in: closestPLANE
        call check_inbox(xcp,xstenbox,nhalfbox,inboxflag)
       else
        inboxflag=0
       endif

       return
       end subroutine closestPLANE

        ! get the slope of the material whose interface is closest to
        ! the center of the cell.
       subroutine get_primary_slope( &
         bfact,dx,xsten0,nhalf0, &
         mofdata, &
         slope,imslope,nmat,sdim)
       use geometry_intersect_module
       IMPLICIT NONE

       INTEGER_T imslope,nmat,sdim,ngeom_recon,bfact,nhalf0
       REAL_T mofdata(nmat*(2*sdim+3))
       REAL_T mofdatavalid(nmat*(2*sdim+3))
       REAL_T xsten0(-nhalf0:nhalf0,sdim)
       REAL_T dx(sdim)
       REAL_T slope(sdim)
       REAL_T slope_test(sdim)
       REAL_T intercept,intercept_min
       REAL_T vfrac_data(nmat)
       INTEGER_T sorted_list(nmat)
       REAL_T uncaptured_volume
       INTEGER_T im,vofcomp,im_exclude,irank,testflag,dir

       if (bfact.lt.1) then
        print *,"bfact invalid"
        stop
       endif
       if (nhalf0.lt.1) then
        print *,"nhalf0 invalid"
        stop
       endif

       ngeom_recon=2*sdim+3
       if ((sdim.ne.3).and.(sdim.ne.2)) then
        print *,"sdim invalid get_primary_slope"
        stop
       endif
       if ((nmat.lt.1).or.(nmat.gt.MAX_NUM_MATERIALS)) then
        print *,"nmat invalid get_primary_slope"
        stop
       endif
     
       imslope=0
       do dir=1,sdim
        slope(dir)=zero
       enddo
       intercept_min=1.0e+20
 
       call make_vfrac_sum_ok(mofdata,mofdatavalid,nmat,ngeom_recon,sdim,3)

       do im=1,nmat
        vofcomp=(im-1)*ngeom_recon+1
        vfrac_data(im)=mofdatavalid(vofcomp)
       enddo
       im_exclude=0
       call sort_volume_fraction(vfrac_data,im_exclude,sorted_list,nmat)
       im=sorted_list(1)
       if (vfrac_data(im).ge.one-VOF_CUTOFF) then
        ! do nothing, there are no reconstructed interfaces in the cell.
       else if ((vfrac_data(im).ge.VOF_CUTOFF).and. &
                (vfrac_data(im).lt.one-VOF_CUTOFF)) then

        uncaptured_volume=one
        irank=1
        do while ((irank.le.nmat).and.(uncaptured_volume.gt.zero))
         do im=1,nmat
          vofcomp=(im-1)*ngeom_recon+1
          testflag=NINT(mofdatavalid(vofcomp+sdim+1))
          if (testflag.eq.irank) then
           do dir=1,sdim
            slope_test(dir)=mofdatavalid(vofcomp+sdim+1+dir)
           enddo
           intercept=mofdatavalid(vofcomp+2*sdim+2)

           uncaptured_volume=uncaptured_volume-mofdatavalid(vofcomp)
           if (uncaptured_volume.lt.VOF_CUTOFF) then
            uncaptured_volume=zero
           endif
           if (uncaptured_volume.gt.zero) then ! we have a valid interface.
            if (abs(intercept).lt.abs(intercept_min)) then
             imslope=im
             intercept_min=intercept
             do dir=1,sdim
              slope(dir)=slope_test(dir)
             enddo
            endif
           endif
          endif ! testflag=irank?
         enddo ! im
         irank=irank+1
        enddo
       else
        print *,"vfrac_data max out of range"
        stop
       endif ! vfrac(sorted_list(1))>1-eps ?

       return
       end subroutine get_primary_slope

       
      subroutine multi_get_distance( &
        bfact,dx,xsten_recon,nhalf_recon,xgrid, &
        mofdata, &
        multi_distance,multi_normal,nmat,nten,rz_flag,sdim, &
        center_stencil,im_solid,time,donateflag, &
        truncate_volume_fractions, &
        truncate_at_junction)
      use geometry_intersect_module
      use global_utility_module

      IMPLICIT NONE

      INTEGER_T nhalf_recon
      REAL_T time
      INTEGER_T im_solid,bfact
      INTEGER_T nmat,nten,sdim,dir,side,l,rz_flag,center_stencil
      INTEGER_T truncate_volume_fractions(nmat)
      INTEGER_T truncate_at_junction
      INTEGER_T donateflag(nmat+2)
      REAL_T mofdata(nmat*(2*sdim+3))
      REAL_T mofdatavalid(nmat*(2*sdim+3))
      REAL_T xsten_recon(-nhalf_recon:nhalf_recon,sdim)
      REAL_T xstenface_recon(-nhalf_recon:nhalf_recon,sdim)
      REAL_T x0side(sdim)
      REAL_T dx(sdim)
      REAL_T multi_distance(nmat+nten)
      REAL_T multi_normal(sdim*(nmat+nten))
      INTEGER_T irank,vofcomp,im,im_plus,im_minus,im0
      INTEGER_T im_plus2,im_minus2
      INTEGER_T ngeom_recon
      INTEGER_T im_exclude
      REAL_T vfrac_data(nmat)
      INTEGER_T sorted_list(nmat)
      REAL_T uncaptured_volume
      INTEGER_T testflag
      REAL_T slopes(sdim)
      REAL_T intercept,intercept_face
      INTEGER_T nten_test
      REAL_T maxdx,normgphi
      INTEGER_T inboxflag
      REAL_T xx(sdim)
      REAL_T gphi(sdim)
      REAL_T slope_list(nmat,sdim)
      REAL_T intercept_list(nmat)
      INTEGER_T im_list(nmat)
      INTEGER_T nlist
      REAL_T x0(sdim)
      INTEGER_T dir1,dir2,side1,side2,ilist,jlist
      INTEGER_T n_im
      INTEGER_T im_test(6)
      REAL_T x0face(sdim)
      INTEGER_T check_override
      REAL_T xgrid(sdim)
      REAL_T xgrid_cen(sdim)
      REAL_T xgrid_plus(sdim)
      REAL_T xgrid_minus(sdim)
      REAL_T xgrid_plus2(sdim)
      REAL_T xgrid_minus2(sdim)

      if (nhalf_recon.lt.1) then
       print *,"nhalf_recon invalid"
       stop
      endif
      nten_test=( (nmat-1)*(nmat-1)+nmat-1 )/2
      if (nten_test.ne.nten) then
       print *,"nten invalid multidist nten nten_test ",nten,nten_test
       stop
      endif
      if ((center_stencil.ne.0).and.(center_stencil.ne.1)) then
       print *,"center_stencil invalid"
       stop
      endif
      if ((im_solid.lt.0).or.(im_solid.gt.nmat)) then
       print *,"im_solid invalid"
       stop
      endif
      maxdx=dx(1)
      if (dx(2).gt.maxdx) then
       maxdx=dx(2)
      endif
      if (dx(sdim).gt.maxdx) then
       maxdx=dx(sdim)
      endif

      ngeom_recon=2*sdim+3

      if ((sdim.ne.3).and.(sdim.ne.2)) then
       print *,"sdim invalid multi_get_distance"
       stop
      endif
      if ((nmat.lt.1).or.(nmat.gt.MAX_NUM_MATERIALS)) then
       print *,"nmat invalid multi_get_distance"
       stop
      endif

      call make_vfrac_sum_ok(mofdata,mofdatavalid,nmat,ngeom_recon,sdim,3)

      do im=1,nmat
       vofcomp=(im-1)*ngeom_recon+1
       vfrac_data(im)=mofdatavalid(vofcomp)
      enddo
      im_exclude=0
      call sort_volume_fraction(vfrac_data,im_exclude,sorted_list,nmat)
      im=sorted_list(1)
       ! a full cell, so distance is either +bigdist or -bigdist,
       ! and a default normal is used.
      if (vfrac_data(im).ge.one-VOF_CUTOFF) then
       ! do nothing, there are no reconstructed interfaces in the cell.
      else if ((vfrac_data(im).ge.VOF_CUTOFF).and. &
               (vfrac_data(im).lt.one-VOF_CUTOFF)) then

        ! if check_override=1, then 
        !  if phi_override(x0)>0 then im0=im_override
       check_override=1
       do dir=1,sdim
        x0(dir)=xsten_recon(0,dir)
       enddo
       call multi_get_volumePOINT( &
        bfact,dx,xsten_recon,nhalf_recon, &
        ngeom_recon,mofdata,x0,im0,nmat,sdim, &
        time,check_override)
       check_override=0

       uncaptured_volume=one
       irank=1
       nlist=0
       do while ((irank.le.nmat).and.(uncaptured_volume.gt.zero))
        do im=1,nmat
         vofcomp=(im-1)*ngeom_recon+1
         testflag=NINT(mofdatavalid(vofcomp+sdim+1))
         if (testflag.eq.irank) then
          do dir=1,sdim
           slopes(dir)=mofdatavalid(vofcomp+sdim+1+dir)
          enddo
          intercept=mofdatavalid(vofcomp+2*sdim+2)

          uncaptured_volume=uncaptured_volume-mofdatavalid(vofcomp)
          if (uncaptured_volume.lt.VOF_CUTOFF) then
           uncaptured_volume=zero
          endif
          if (uncaptured_volume.gt.zero) then ! we have a valid interface.
           nlist=nlist+1
           if ((nlist.ge.1).and.(nlist.lt.nmat)) then
            do dir=1,sdim
             slope_list(nlist,dir)=slopes(dir)
            enddo
            intercept_list(nlist)=intercept
            im_list(nlist)=im
           else 
            print *,"nlist invalid"
            stop
           endif
           call closestPLANE(bfact,dx,xsten_recon,nhalf_recon, &
            xgrid_cen,xgrid_plus,xgrid_minus, &
            xgrid,maxdx,slopes,intercept, &
            xsten_recon,nhalf_recon,sdim,inboxflag) 

           if (inboxflag.eq.1) then
            call multi_get_volumePOINT( &
             bfact,dx,xsten_recon,nhalf_recon, &
             ngeom_recon,mofdata,xgrid_plus, &
             im_plus,nmat,sdim,time,check_override)
            call multi_get_volumePOINT( &
             bfact,dx,xsten_recon,nhalf_recon, &
             ngeom_recon,mofdata,xgrid_minus, &
             im_minus,nmat,sdim,time,check_override)
            n_im=2
            im_test(1)=im_plus
            im_test(2)=im_minus
            if (center_stencil.eq.1) then
             n_im=n_im+1
             im_test(n_im)=im0
            endif
            if (irank.eq.1) then
             n_im=n_im+1
             im_test(n_im)=im
            endif
             ! xgrid_cen is closest point on the plane.
            call compare_distance( &
             bfact,dx,xsten_recon,nhalf_recon, &
             nmat,nten,xgrid,xgrid_cen, &
             multi_distance,multi_normal,im_test,n_im,slopes, &
             im,im0,sdim,center_stencil,im_solid,donateflag, &
             truncate_volume_fractions, &
             truncate_at_junction) 
           endif

           do dir=1,sdim
           do side=-1,1,2
            do l=1,sdim
             xx(l)=xgrid(l)
             gphi(l)=slopes(l)
             x0face(l)=x0(l)
            enddo
            x0face(dir)=xsten_recon(side,dir)
            xx(dir)=x0face(dir)
            gphi(dir)=zero
            normgphi=zero
            do l=1,sdim
             normgphi=normgphi+gphi(l)**2
            enddo 
            normgphi=sqrt(normgphi)
  
            if (normgphi.ge.1.0D-12) then 
              !  n dot (x-x0)+intercept=0
              !  n dot (x-x0+x0face-x0face)+intercept=0
              !  n dot (x-x0face) + n dot (x0face-x0)+intercept=0
             intercept_face=intercept

             do dir2=1,sdim
              intercept_face=intercept_face+ &
                 slopes(dir2)*(x0face(dir2)-x0(dir2))
              xstenface_recon(0,dir2)=x0face(dir2)
             enddo
           
             call closestPLANE( &
              bfact,dx,xstenface_recon,nhalf_recon, &
              xgrid_cen,xgrid_plus,xgrid_minus, &
              xx,maxdx,gphi,intercept_face, &
              xsten_recon,nhalf_recon,sdim,inboxflag) 

             if (inboxflag.eq.1) then
              call multi_get_volumePOINT( &
                bfact,dx,xsten_recon,nhalf_recon, &
                ngeom_recon,mofdata,xgrid_plus, &
                im_plus,nmat,sdim,time,check_override)
              call multi_get_volumePOINT( &
                bfact,dx,xsten_recon,nhalf_recon, &
                ngeom_recon,mofdata,xgrid_minus, &
                im_minus,nmat,sdim,time,check_override)
              n_im=2
              im_test(1)=im_plus
              im_test(2)=im_minus
              if (center_stencil.eq.1) then
               n_im=n_im+1
               im_test(n_im)=im0
              endif
              if (irank.eq.1) then
               n_im=n_im+1
               im_test(n_im)=im
              endif

              call compare_distance( &
               bfact,dx,xsten_recon,nhalf_recon, &
               nmat,nten,xgrid,xgrid_cen, &
               multi_distance,multi_normal,im_test,n_im,slopes, &
               im,im0,sdim,center_stencil,im_solid,donateflag, &
               truncate_volume_fractions, &
               truncate_at_junction) 
             endif

             if (sdim.eq.2) then
              ! do nothing
             else if (sdim.eq.3) then
  
              do dir1=1,SDIM
              do side1=-1,1,2
               if (dir1.ne.dir) then
                do l=1,sdim
                 xx(l)=xgrid(l)
                 gphi(l)=slopes(l)
                 x0face(l)=x0(l)
                enddo
                x0face(dir)=xsten_recon(side,dir)
                x0face(dir1)=xsten_recon(side1,dir1)
                xx(dir)=x0face(dir)
                xx(dir1)=x0face(dir1)
                gphi(dir)=zero
                gphi(dir1)=zero
                normgphi=zero
                do l=1,SDIM
                 normgphi=normgphi+gphi(l)**2
                enddo 
                normgphi=sqrt(normgphi)
   
                if (normgphi.ge.1.0D-12) then
                 intercept_face=intercept
                 do dir2=1,sdim
                  intercept_face=intercept_face+ &
                   slopes(dir2)*(x0face(dir2)-x0(dir2))
                  xstenface_recon(0,dir2)=x0face(dir2)
                 enddo
                 call closestPLANE( &
                  bfact,dx,xstenface_recon,nhalf_recon, &
                  xgrid_cen,xgrid_plus,xgrid_minus, &
                  xx,maxdx,gphi,intercept_face, &
                  xsten_recon,nhalf_recon,sdim,inboxflag) 

                 if (inboxflag.eq.1) then
                  call multi_get_volumePOINT( &
                    bfact,dx,xsten_recon,nhalf_recon, &
                    ngeom_recon,mofdata,xgrid_plus, &
                    im_plus,nmat,sdim,time,check_override)
                  call multi_get_volumePOINT( &
                    bfact,dx,xsten_recon,nhalf_recon, &
                    ngeom_recon,mofdata,xgrid_minus, &
                    im_minus,nmat,sdim,time,check_override)
                  n_im=2
                  im_test(1)=im_plus
                  im_test(2)=im_minus
                  if (center_stencil.eq.1) then
                   n_im=n_im+1
                   im_test(n_im)=im0
                  endif
                  if (irank.eq.1) then
                   n_im=n_im+1
                   im_test(n_im)=im
                  endif

                  call compare_distance( &
                   bfact,dx,xsten_recon,nhalf_recon, &
                   nmat,nten,xgrid,xgrid_cen, &
                   multi_distance,multi_normal,im_test,n_im, &
                   slopes,im,im0,sdim,center_stencil,im_solid,donateflag, &
                   truncate_volume_fractions, &
                   truncate_at_junction)
                 endif
                endif ! normgphi>0 ?
               endif ! dir1<>dir
              enddo ! side1
              enddo ! dir1

             else
              print *,"sdim invalid"
              stop
             endif
            endif ! normgphi>0 ?
           enddo ! side
           enddo ! dir

          endif ! uncaptured_volume>0
         endif  ! testflag=irank
        enddo ! im
        irank=irank+1
       enddo  ! while irank<=nmat and uncaptured_volume>0 

       check_override=0
       if (nlist.ge.2) then
        do ilist=1,nlist-1
        do jlist=ilist+1,nlist
      

          ! xgrid_cen=xgrid-d n
          ! n=(xgrid-xgrid_cen)/d 
          ! if no intersection, or not in the box, then
          ! inboxflag=0.
         call closestINT( &
          bfact,dx,xsten_recon,nhalf_recon, &
          xgrid_cen,xgrid_plus,xgrid_minus, &
          xgrid_plus2,xgrid_minus2, &
          xgrid,maxdx,ilist,jlist,nlist, &
          slope_list,intercept_list,im_list,rz_flag,sdim, &
          inboxflag,nmat)
         if (inboxflag.eq.1) then
          call multi_get_volumePOINT( &
            bfact,dx,xsten_recon,nhalf_recon, &
            ngeom_recon,mofdata,xgrid_plus, &
            im_plus,nmat,sdim,time,check_override)
          call multi_get_volumePOINT( &
            bfact,dx,xsten_recon,nhalf_recon, &
            ngeom_recon,mofdata,xgrid_minus, &
            im_minus,nmat,sdim,time,check_override)
          call multi_get_volumePOINT( &
            bfact,dx,xsten_recon,nhalf_recon, &
            ngeom_recon,mofdata,xgrid_plus2, &
            im_plus2,nmat,sdim,time,check_override)
          call multi_get_volumePOINT( &
            bfact,dx,xsten_recon,nhalf_recon, &
            ngeom_recon,mofdata,xgrid_minus2, &
            im_minus2,nmat,sdim,time,check_override)
          n_im=4
          im_test(1)=im_plus
          im_test(2)=im_minus
          im_test(3)=im_plus2
          im_test(4)=im_minus2
          if (center_stencil.eq.1) then
           n_im=n_im+1
           im_test(n_im)=im0
          endif
          do dir2=1,sdim
           slopes(dir2)=slope_list(ilist,dir2)
          enddo
          im=im_list(ilist)
          call compare_distance( &
            bfact,dx,xsten_recon,nhalf_recon, &
            nmat,nten,xgrid,xgrid_cen, &
            multi_distance,multi_normal,im_test,n_im, &
            slopes,im,im0,sdim,center_stencil,im_solid,donateflag, &
            truncate_volume_fractions, &
            truncate_at_junction) 
         endif ! inboxflag=1 ?

          ! this is case where 2 planes, and a cell wall plane have
          ! a single common intersection point.
         if (sdim.eq.3) then

          do dir=1,sdim
          do side=-1,1,2
           do l=1,sdim
            xx(l)=xgrid(l)
            gphi(l)=0.0
            x0side(l)=x0(l)
           enddo
           gphi(dir)=one
           xx(dir)=xsten_recon(side,dir)
           x0side(dir)=xx(dir)

           call closestINT3( &
            bfact,dx,xsten_recon,nhalf_recon, &
            xgrid_cen,xgrid_plus,xgrid_minus, &
            xgrid_plus2,xgrid_minus2, &
            xx,maxdx,ilist,jlist,nlist, &
            slope_list,intercept_list,im_list, &
            gphi, &
            x0side,rz_flag,sdim, &
            inboxflag,nmat)

           if (inboxflag.eq.1) then
            call multi_get_volumePOINT( &
             bfact,dx,xsten_recon,nhalf_recon, &
             ngeom_recon,mofdata,xgrid_plus, &
             im_plus,nmat,sdim,time,check_override)
            call multi_get_volumePOINT( &
             bfact,dx,xsten_recon,nhalf_recon, &
             ngeom_recon,mofdata,xgrid_minus, &
             im_minus,nmat,sdim,time,check_override)
            call multi_get_volumePOINT( &
             bfact,dx,xsten_recon,nhalf_recon, &
             ngeom_recon,mofdata,xgrid_plus2, &
             im_plus2,nmat,sdim,time,check_override)
            call multi_get_volumePOINT( &
             bfact,dx,xsten_recon,nhalf_recon, &
             ngeom_recon,mofdata,xgrid_minus2, &
             im_minus2,nmat,sdim,time,check_override)

            n_im=4
            im_test(1)=im_plus
            im_test(2)=im_minus
            im_test(3)=im_plus2
            im_test(4)=im_minus2
            if (center_stencil.eq.1) then
             n_im=n_im+1
             im_test(n_im)=im0
            endif
            do dir2=1,sdim
             slopes(dir2)=slope_list(ilist,dir2)
            enddo
            im=im_list(ilist)

            call compare_distance( &
              bfact,dx,xsten_recon,nhalf_recon, &
              nmat,nten,xgrid,xgrid_cen, &
              multi_distance,multi_normal,im_test,n_im, &
              slopes,im,im0,sdim,center_stencil,im_solid,donateflag, &
              truncate_volume_fractions, &
              truncate_at_junction) 
           endif ! inboxflag=1 ?

          enddo ! side
          enddo ! dir
 
         else if (sdim.eq.2) then
           ! do nothing
         else
          print *,"sdim invalid"
          stop
         endif
        enddo 
        enddo  ! ilist,jlist
       endif ! nlist >=2

 
      else
       print *,"vfrac_data max out of range"
       stop
      endif ! vfrac(sorted_list(1))>1-eps ?

      return
      end subroutine multi_get_distance


! vof, ref centroid, order,slope,intercept  x nmat
! phi=n dot (x-x0) + intercept
! x0 is center of cell (not centroid)
! returns values -1 or 1.
      subroutine multi_get_volumePOINT( &
       bfact,dx,xsten0,nhalf0, &
       ngeom_recon_in,mofdata,xgrid, &
       im_crit,nmat,sdim,time,check_override)

#if (STANDALONE==0)
      use global_distance_module
#endif
      use geometry_intersect_module
      use global_utility_module

      IMPLICIT NONE

      INTEGER_T nmat,sdim,check_override,ngeom_recon_in,bfact,nhalf0
      REAL_T time
      REAL_T mofdata(nmat*ngeom_recon_in)
      REAL_T mofdatavalid(nmat*ngeom_recon_in)
      REAL_T xsten0(-nhalf0:nhalf0,sdim)
      REAL_T dx(sdim)
      REAL_T xgrid(sdim)
      INTEGER_T im_crit
      INTEGER_T irank,vofcomp,im
      REAL_T uncaptured_volume
      INTEGER_T testflag,dir
      REAL_T slopes(sdim)
      REAL_T intercept,ls,maxvof
      INTEGER_T im_exclude
      REAL_T vfrac_data(nmat)
      INTEGER_T vfrac_checked(nmat)
      INTEGER_T sorted_list(nmat)
      INTEGER_T im_solid_positive
      REAL_T LSSOLID

#include "mofdata.H"
 
      if (ngeom_recon_in.ne.(2*sdim+3)) then
       print *,"ngeom_recon_in invalid"
       stop
      endif

      if (bfact.lt.1) then
       print *,"bfact invalid"
       stop
      endif
      if ((sdim.ne.3).and.(sdim.ne.2)) then
       print *,"sdim invalid multi_get_volumePOINT"
       stop
      endif
      if ((nmat.lt.1).or.(nmat.gt.MAX_NUM_MATERIALS)) then
       print *,"nmat invalid multi get volume point"
       stop
      endif
      if (time.lt.zero) then
       print *,"time invalid"
       stop
      endif

      im_exclude=0 
      im_crit=0
      do im=1,nmat
       vfrac_checked(im)=0
      enddo

      call make_vfrac_sum_ok(mofdata,mofdatavalid,nmat,ngeom_recon_in,sdim,3)
      do im=1,nmat
       vofcomp=(im-1)*ngeom_recon_in+1
       vfrac_data(im)=mofdatavalid(vofcomp)
      enddo
      call sort_volume_fraction(vfrac_data,im_exclude,sorted_list,nmat)
      im_crit=sorted_list(1)
      if (vfrac_data(im_crit).ge.one-VOF_CUTOFF) then
       ! do nothing, im_crit is set.
      else
    
       im_solid_positive=0 
       if (check_override.eq.0) then
        ! do nothing
       else if (check_override.eq.1) then
        if (imaterial_override.eq.0) then
         ! do nothing
        else if ((imaterial_override.ge.1).and. &
                 (imaterial_override.le.nmat)) then
#if (STANDALONE==0)
         LSSOLID=distance_function_override(xgrid(1),xgrid(2), &
          xgrid(sdim),time)
         if (LSSOLID.ge.zero) then
          im_solid_positive=1
         endif
#else
         print *,"bust multi_get_volumePOINT"
         print *,"sdim,nmat ",sdim,nmat
         stop
#endif
        else
         print *,"imaterial_override invalid"
         stop
        endif
       else
        print *,"check_override invalid"
        stop
       endif

       if (im_solid_positive.eq.1) then

        im_crit=imaterial_override

       else if (im_solid_positive.eq.0) then

        uncaptured_volume=one
        irank=1
        do while ((irank.le.nmat).and.(uncaptured_volume.gt.zero))

         do im=1,nmat
          vofcomp=(im-1)*ngeom_recon_in+1
          testflag=NINT(mofdatavalid(vofcomp+sdim+1))
          if (testflag.eq.irank) then
           vfrac_checked(im)=1
           do dir=1,sdim
            slopes(dir)=mofdatavalid(vofcomp+sdim+1+dir)
           enddo
           intercept=mofdatavalid(vofcomp+2*sdim+2)
           call distfunc(bfact,dx,xsten0,nhalf0,intercept,slopes,xgrid,ls,sdim)

           if ((ls.ge.zero).or. &
               (mofdatavalid(vofcomp).ge.uncaptured_volume-VOF_CUTOFF)) then
            im_crit=im
            uncaptured_volume=zero
           else 
            uncaptured_volume=uncaptured_volume-mofdatavalid(vofcomp)
           endif
          endif  ! testflag=irank
         enddo ! im
         irank=irank+1
        enddo  ! while irank<=nmat and uncaptured_volume>0 

        if (uncaptured_volume.eq.zero) then
         ! do nothing
        else if (uncaptured_volume.gt.zero) then
         im_crit=0
         maxvof=zero
         do im=1,nmat
          if (vfrac_checked(im).eq.1) then
           ! do nothing
          else if (vfrac_checked(im).eq.0) then
           if (vfrac_data(im).gt.maxvof) then
            maxvof=vfrac_data(im)
            im_crit=im
           endif
          else
           print *,"vfrac_checked invalid"
           stop
          endif
         enddo ! im
         if (maxvof.le.zero) then
          print *,"failed to find material that covers point"
          do im=1,nmat
           vofcomp=(im-1)*ngeom_recon_in+1
           print *,"im,vof,flag,int ",im,mofdatavalid(vofcomp), &
            mofdatavalid(vofcomp+sdim+1), &
            mofdatavalid(vofcomp+2*sdim+2)
          enddo
          print *,"xgrid,xsten0 ",xgrid(1),xgrid(2),xsten0(0,1),xsten0(0,2)
          stop
         else
          uncaptured_volume=zero
         endif 
        else 
         print *,"uncaptured_volume invalid"
         stop
        endif

       else
        print *,"im_solid_positive invalid"
        stop
       endif

      endif ! vfrac(sorted_list(1))>1-eps ?

      return
      end subroutine multi_get_volumePOINT

 
        ! sort from largest volume fraction to smallest
      subroutine sort_volume_fraction( &
       vfrac_data,im_exclude,sorted_list,nmat)

      use geometry_intersect_module

      IMPLICIT NONE

      INTEGER_T nmat,im_exclude
      REAL_T vfrac_data(nmat)
      INTEGER_T sorted_list(nmat)
      INTEGER_T im,changed,nsweeps,swap,do_swap

      if ((nmat.lt.1).or.(nmat.gt.MAX_NUM_MATERIALS)) then
       print *,"nmat invalid sort_volume_fraction"
       print *,"nmat= ",nmat
       stop
      endif
      if ((im_exclude.lt.0).or.(im_exclude.gt.nmat)) then
       print *,"im_exclude invalid"
       stop
      endif

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
      end subroutine sort_volume_fraction
       

end module MOF_routines_module



#if (STANDALONE==0)
      subroutine FORT_INITMOF(order_algorithm_in, &
       nmat,MOFITERMAX_in, &
       imaterial_override_in, &
       MOF_DEBUG_RECON_in, &
       MOF_TURN_OFF_LS_in)
#elif (STANDALONE==1)
      subroutine initmof(order_algorithm_in, &
       nmat,MOFITERMAX_in, &
       imaterial_override_in, &
       MOF_DEBUG_RECON_in, &
       MOF_TURN_OFF_LS_in)
#else
      print *,"bust initmof"
      stop
#endif

      use geometry_intersect_module
      use MOF_routines_module

      IMPLICIT NONE

      INTEGER_T nmat,im
      INTEGER_T imaterial_override_in
      INTEGER_T order_algorithm_in(nmat)
      INTEGER_T MOFITERMAX_in
      INTEGER_T MOF_DEBUG_RECON_in
      INTEGER_T MOF_TURN_OFF_LS_in
      INTEGER_T sdim,nmat_test,nmax_test

#include "mofdata.H"

      MOF_DEBUG_RECON_COUNT=0
      MOF_DEBUG_RECON=MOF_DEBUG_RECON_in
      MOF_TURN_OFF_LS=MOF_TURN_OFF_LS_in

      call set_MOFITERMAX(MOFITERMAX_in)

      if ((nmat.lt.1).or.(nmat.gt.MAX_NUM_MATERIALS)) then
       print *,"nmat invalid init mof"
       stop
      endif
      if ((imaterial_override_in.lt.0).or. &
          (imaterial_override_in.gt.nmat)) then
       print *,"imaterial override in invalid"
       stop
      endif

      imaterial_override=imaterial_override_in

      call set_order_algorithm(order_algorithm_in,nmat)

      print *,"initializing geometry tables"

      call init_geometry_tables()
      if (1.eq.0) then
       call volume_sanity_check()
      endif
      if (1.eq.0) then
       sdim=2
       nmat_test=2
       nmax_test=400
       call diagnostic_MOF(sdim,nmat_test,nmax_test)
       stop
      endif

      return
      end


