!{\src2tex{textfont=tt}}
!!****f* ABINIT/mati3inv
!! NAME
!! mati3inv
!!
!! FUNCTION
!! Invert and transpose orthogonal 3x3 matrix of INTEGER elements.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2010 ABINIT group (DCA, XG, GMR)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!! mm = integer matrix to be inverted
!!
!! OUTPUT
!! mit = inverse of mm input matrix
!!
!! SIDE EFFECTS
!!
!! NOTES
!! Used for symmetry operations.
!! This routine applies to ORTHOGONAL matrices only.
!! Since these form a group, inverses are also integer
!! arrays.  Returned array is TRANSPOSE of inverse, as needed.
!! Note use of integer arithmetic.
!! Also: has been designed so that mit can be same storage space as m, in
!! which case m is overwritten by resulting mit.
!!
!! TODO
!!
!! PARENTS
!!      chkgrp,classify_bands,debug_tools,get_full_kgrid,getkgrid,ingeo
!!      invars2m,m_ab6_symmetry_f90,m_crystal,m_fft_mesh,m_ptgroups,nstdy3
!!      optic,outscfcv,rdddb9,read_gkk,setsym,strainsym,symdij,symdyma,wfconv
!!
!! CHILDREN
!!      leave_new,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine mati3inv(mm,mit)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
 use interfaces_14_hidewrite
 use interfaces_16_hideleave
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!arrays
 integer,intent(in) :: mm(3,3)
 integer,intent(out) :: mit(3,3)

!Local variables-------------------------------
!scalars
 integer :: dd
 character(len=500) :: message
!arrays
 integer :: tt(3,3)

! *************************************************************************

 tt(1,1) = mm(2,2) * mm(3,3) - mm(3,2) * mm(2,3)
 tt(2,1) = mm(3,2) * mm(1,3) - mm(1,2) * mm(3,3)
 tt(3,1) = mm(1,2) * mm(2,3) - mm(2,2) * mm(1,3)
 tt(1,2) = mm(3,1) * mm(2,3) - mm(2,1) * mm(3,3)
 tt(2,2) = mm(1,1) * mm(3,3) - mm(3,1) * mm(1,3)
 tt(3,2) = mm(2,1) * mm(1,3) - mm(1,1) * mm(2,3)
 tt(1,3) = mm(2,1) * mm(3,2) - mm(3,1) * mm(2,2)
 tt(2,3) = mm(3,1) * mm(1,2) - mm(1,1) * mm(3,2)
 tt(3,3) = mm(1,1) * mm(2,2) - mm(2,1) * mm(1,2)
 dd  = mm(1,1) * tt(1,1) + mm(2,1) * tt(2,1) + mm(3,1) * tt(3,1)
!Make sure matrix is not singular
 if (dd/=0) then
   mit(:,:)=tt(:,:)/dd
 else
   write(message, '(5a,2x,9i5,a)' ) ch10,&
&   ' mati3inv : BUG -',ch10,&
&   '  Attempting to invert integer array',ch10,&
&   mm(:,:),'   ==> determinant is zero.'
   call wrtout(std_out,message,'COLL')
   call leave_new('COLL')
 end if

end subroutine mati3inv
!!***
