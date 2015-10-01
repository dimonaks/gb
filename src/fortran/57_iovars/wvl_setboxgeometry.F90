!{\src2tex{textfont=tt}}
!!****f* ABINIT/wvl_setBoxGeometry
!! NAME
!! wvl_setBoxGeometry
!!
!! FUNCTION
!! When wavelets are used, the box definition needs to be changed.
!! The box size is recomputed knowing some psp informations such as
!! the radius for coarse and fine grid. Then, the atoms are translated
!! to be included in the new box. Finally the FFT grid is computed using
!! the fine wavelet mesh and a buffer characteristic of used wavelets plus
!! a buffer used to be multiple of 2, 3 or 5.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2010 ABINIT group (DC)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  me=informations about MPI parallelizationthe processor id (from mpi_enreg%me)
!!  psps <type(pseudopotential_type)>=variables related to pseudopotentials
!!  radii= the radii for each type of atoms, giving the fine and the coarse
!!         grid.
!!
!! OUTPUT
!!  rprimd(3,3)=dimensional primitive translations in real space (bohr)
!!
!! SIDE EFFECTS
!!  dtset <type(dataset_type)>=internal variables used by wavelets, describing
!!                             the box are set. The FFT grid is also changed.
!!  xred(3,natom)=reduced dimensionless atomic coordinates
!!
!! PARENTS
!!      gstate,gstateimg,wvl_memory,wvl_wfsinp_reformat
!!
!! CHILDREN
!!      leave_new,mkrdim,system_size,wrtout,xredxcart
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine wvl_setBoxGeometry(dtset, me, radii, rprimd, xred)

 use defs_basis
 use defs_datatypes
 use defs_wvltypes
 use defs_abitypes
#if defined HAVE_BIGDFT
  use BigDFT_API, only: system_size
#endif

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
 use interfaces_14_hidewrite
 use interfaces_16_hideleave
 use interfaces_42_geometry
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: me
 type(dataset_type),intent(inout) :: dtset
!arrays
 real(dp),intent(in) :: radii(dtset%ntypat,2)
 real(dp),intent(inout) :: rprimd(3,3),xred(3,dtset%natom)

!Local variables-------------------------------
!scalars
 integer :: i
 character(len=500) :: message
!arrays
 real(dp) :: rprim(3,3),acell(3)
 real(dp),allocatable :: xcart(:,:)

! *********************************************************************

#if defined HAVE_BIGDFT
 if (dtset%prtvol == 0) then
   write(message, '(a,a,a,a)' ) ch10,&
&   ' wvl_setBoxGeometry : Changing the box for wavelets computation.'
   call wrtout(std_out,message,'COLL')
 end if

!Store xcart for each atom
 allocate(xcart(3, dtset%natom))
 call xredxcart(dtset%natom, 1, rprimd, xcart, xred)

 call system_size(me, dtset%wvl%atoms, xcart, radii, dtset%wvl_crmult, &
& dtset%wvl_frmult, dtset%wvl%h(1), dtset%wvl%h(2), dtset%wvl%h(3), &
& dtset%wvl%n(1), dtset%wvl%n(2), dtset%wvl%n(3), &
& dtset%wvl%fGrid(1, 1), dtset%wvl%fGrid(1, 2), &
& dtset%wvl%fGrid(1, 3), dtset%wvl%fGrid(2, 1), &
& dtset%wvl%fGrid(2, 2), dtset%wvl%fGrid(2, 3), &
& dtset%wvl%ni(1), dtset%wvl%ni(2), dtset%wvl%ni(3))

 acell(1) = dtset%wvl%atoms%alat1
 acell(2) = dtset%wvl%atoms%alat2
 acell(3) = dtset%wvl%atoms%alat3

 if (dtset%prtvol == 0) then
   write(message, '(a,3F12.6)' ) &
&   '  | acell is now:         ', acell
   call wrtout(std_out,message,'COLL')
   write(message, '(a,2I5,a,a,2I5,a,a,2I5)' ) &
&   '  | nfl1, nfu1:           ', dtset%wvl%fGrid(:, 1), ch10, &
&   '  | nfl2, nfu2:           ', dtset%wvl%fGrid(:, 2), ch10, &
&   '  | nfl3, nfu3:           ', dtset%wvl%fGrid(:, 3)
   call wrtout(std_out,message,'COLL')
 end if

!Change the metric to orthogonal one
 rprim(:, :) = real(0., dp)
 do i = 1, 3, 1
   rprim(i, i) = real(1., dp)
 end do
 call mkrdim(acell, rprim, rprimd)

!Save shifted atom positions into xred
 call xredxcart(dtset%natom, -1, rprimd, xcart, xred)
 deallocate(xcart)

!Set the number of points for the complete system (MPI or not)
 dtset%wvl%ntot = product(dtset%wvl%ni)

 if (dtset%prtvol == 0) then
   write(message, '(a,3I12)' ) &
&   '  | box size for datas:   ', dtset%wvl%ni
   call wrtout(std_out,message,'COLL')
   write(message, '(a,3I12)' ) &
&   '  | box size for wavelets:', dtset%wvl%n
   call wrtout(std_out,message,'COLL')
 end if
 
#else
 write(message, '(a,a,a,a)' ) ch10,&
& ' wvl_setBoxGeometry : BigDFT library is not compiled.', ch10, &
& '   Action, used the flag --enable-bigdft when configuring.'
 call wrtout(std_out,message,'COLL')
 call leave_new('COLL')
#endif
end subroutine wvl_setBoxGeometry
!!***
