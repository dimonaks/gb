!{\src2tex{textfont=tt}}
!!****f* ABINIT/printxsf
!! NAME
!! printxsf
!!
!! FUNCTION
!! Write a generic array in the XSF format (XCrysden format)
!!
!! COPYRIGHT
!!  Copyright (C) 2006-2010 ABINIT group (MG)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!! basis(3,3) = basis vectors of the direct real lattice or of the reciprocal lattice (fortran convention)
!!              (Bohr units if realrecip=0, Bohr^-1 if realrecip=1, see below)
!! realrecip = 0  for a plot in real space
!!             1  for a plot in reciprocal space
!! nunit   = unit number of the output file (already open by the caller, not closed here!)
!! n1=grid size along x
!! n2=grid size along y
!! n3=grid size along z
!! origin(3) = origin of the grid
!! datagrid(n1*n2*n3) = datagrid values stored using the fortran convention
!!
!! OUTPUT
!! Only write
!!
!! PARENTS
!!      mknesting
!!
!! CHILDREN
!!      leave_new,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine printxsf(n1,n2,n3,datagrid,basis,origin,natom, ntypat, typat, xcart, znucl, nunit,realrecip)

 use defs_basis
 use m_errors

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: n1,n2,n3,nunit,realrecip
 integer,intent(in) :: natom, ntypat
!arrays
 real(dp),intent(in) :: basis(3,3),datagrid(n1*n2*n3),origin(3)
 integer,intent(in) :: typat(ntypat)
 real(dp),intent(in) :: xcart(3,natom), znucl(ntypat)

!Local variables-------------------------------
!scalars
 integer :: ix,iy,iz,nslice,nsym,iatom
 real(dp) :: fact
 character(len=500) :: msg
!arrays
 real(dp) :: tau(3,natom)

! *************************************************************************

 DBG_ENTER("COLL")

 if ( ALL(realrecip/=(/0,1/)) )then
   write(msg,'(a,i6)')' The argument realrecip should be 0 or 1; however, realrecip= ',realrecip
   MSG_BUG(msg)
 end if

!conversion between ABINIT default units and XCrysden units
 fact=Bohr_Ang; if (realrecip ==1) fact=one/fact  !since we are in reciprocal space

!TODO insert crystalline structure and dummy atoms in case of reciprocal space
!need to convert basis too

 write(nunit,'(1X,A)')  'DIM-GROUP'
 write(nunit,*) '3  1'
 write(nunit,'(1X,A)') 'PRIMVEC'
 do iy = 1,3
   write(nunit,'(3(ES17.10,2X))') (Bohr_Ang*basis(ix,iy), ix=1,3)
 end do
!
!generate translated coordinates to fit origin shift
!
 do iatom = 1,natom
   tau (:,iatom) = xcart(:,iatom) - origin(:)
 end do

 write(nunit,'(1X,A)') 'PRIMCOORD'
 write(nunit,*) natom, ' 1'
 do iatom = 1,natom
   write(nunit,'(i9,3(3X,ES17.10))') NINT(znucl(typat(iatom))), &  ! WARNING alchemy not supported by XCrysden
&  Bohr_Ang*tau(1,iatom), &
&   Bohr_Ang*tau(2,iatom), &
&   Bohr_Ang*tau(3,iatom)
 end do
 write(nunit,'(1X,A)') 'ATOMS'
 do iatom = 1,natom
   write(nunit,'(i9,3(3X,ES17.10))') NINT(znucl(typat(iatom))), & ! WARNING alchemy not supported by XCrysden
&  Bohr_Ang*tau(1,iatom), &
&   Bohr_Ang*tau(2,iatom), &
&   Bohr_Ang*tau(3,iatom)
 end do

!write(nunit,'(a)')' CRYSTAL'

!if (realrecip == 1) then
!write(nunit,'(a)')' # these are primitive lattice vectors (in Angstroms)'
!write(nunit,'(a)')
!else
!write(nunit,'(a)')' # these are primitive reciprocal lattice vectors (in Angstroms -1)'
!write(nunit,'(a)')
!end if

!write(nunit,'(a)')' PRIMVEC'
!write(nunit,*)basis(:,1)*fact
!write(nunit,*)basis(:,2)*fact
!write(nunit,*)basis(:,3)*fact

!write(nunit,'(a)')'DIM-GROUP'
!write(nunit,'(a)')' 3  1'
!write(nunit,'(a)')' PRIMCOORD'
!write(nunit,'(a)')' 1  1'
!write(nunit,'(a)')'       13    0.0000000000E+00    0.0000000000E+00    0.0000000000E+00'

 write(nunit,'(a)')' BEGIN_BLOCK_DATAGRID3D'
 write(nunit,'(a)')' datagrid'
 write(nunit,'(a)')' DATAGRID_3D_DENSITY'
!NOTE: XCrysden uses aperiodical data grid
 write(nunit,*)n1+1,n2+1,n3+1
 write(nunit,*)origin
 write(nunit,*)basis(:,1)*fact
 write(nunit,*)basis(:,2)*fact
 write(nunit,*)basis(:,3)*fact

 nslice=1
 do iz=1,n3
   do iy=1,n2
     write(nunit,'(8es16.8)') datagrid(1+n1*(nslice-1):n1+n1*(nslice-1)),datagrid(1+n1*(nslice-1))
!    DEBUG
!    write(*,*)1+n1*(nslice-1),n1+n1*(nslice-1)
!    ENDDEBUG
     nslice=nslice+1
   end do
   nsym=nslice-n2
   write (nunit,'(8es16.8)') datagrid(1+n1*(nsym-1):n1+n1*(nsym-1)),datagrid(1+n1*(nsym-1))
!  DEBUG
!  write(*,*)1+n1*(nsym-1),n1+n1*(nsym-1)
!  write(*,*)' done xy plane at z = ',nslice-n2
!  ENDDEBUG
 end do

!Now write upper plane
 nslice=1
 do iy=1,n2
   write (nunit,'(8es16.8)') datagrid(1+n1*(nslice-1):n1+n1*(nslice-1)),datagrid(1+n1*(nslice-1))
   write(*,*)1+n1*(nslice-1),n1+n1*(nslice-1)
   nslice=nslice+1
 end do

 nsym=nslice-n2
 write (nunit,'(8es16.8)') datagrid(1+n1*(nsym-1):n1+n1*(nsym-1)),datagrid(1+n1*(nsym-1))
!write(*,*)1+n1*(nsym-1),n1+n1*(nsym-1)
!write(*,*)' done upper xy plane '

 write (nunit,'(a)')' END_DATAGRID_3D'
 write (nunit,'(a)')' END_BLOCK_DATAGRID3D'

 DBG_EXIT("COLL")

end subroutine printxsf

!!***
