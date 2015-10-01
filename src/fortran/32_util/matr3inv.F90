!{\src2tex{textfont=tt}}
!!****f* ABINIT/matr3inv
!! NAME
!! matr3inv
!!
!! FUNCTION
!! Invert and transpose general 3x3 matrix of real*8 elements.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2010 ABINIT group (DCA, XG, GMR)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!! aa = 3x3 matrix to be inverted
!!
!! OUTPUT
!! ait = inverse of aa input matrix
!!
!! NOTES
!! Returned array is TRANSPOSE of inverse, as needed to get g from r.
!!
!! PARENTS
!!      berryphase,chkdilatmx,conducti_nc,electrooptic,elphon,ewald2,ewald9
!!      get_fsurf_1band,getkgrid,getspinrot,hybrid9,inwffil,metric,mkvxc3
!!      mkvxcstr3,mrgscr,newsp,optic,outwant,planeint,prtxf,rdddb9,recip,reduce
!!      relaxpol,rsiaf9,shellin,smpbz,strsym,symbrav,symdyma,symph3,symrelrot
!!      tddft,testkgrid,tetrahedron,thm9,uderiv,volumeint,xfpack,xredxcart
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine matr3inv(aa,ait)

 use defs_basis

 implicit none

!Arguments ------------------------------------
!arrays
 real(dp),intent(in) :: aa(3,3)
 real(dp),intent(out) :: ait(3,3)

!Local variables-------------------------------
!scalars
 real(dp) :: dd,t1,t2,t3

! *************************************************************************

 t1 = aa(2,2) * aa(3,3) - aa(3,2) * aa(2,3)
 t2 = aa(3,2) * aa(1,3) - aa(1,2) * aa(3,3)
 t3 = aa(1,2) * aa(2,3) - aa(2,2) * aa(1,3)
 dd  = 1.d0/ (aa(1,1) * t1 + aa(2,1) * t2 + aa(3,1) * t3)
 ait(1,1) = t1 * dd
 ait(2,1) = t2 * dd
 ait(3,1) = t3 * dd
 ait(1,2) = (aa(3,1)*aa(2,3)-aa(2,1)*aa(3,3)) * dd
 ait(2,2) = (aa(1,1)*aa(3,3)-aa(3,1)*aa(1,3)) * dd
 ait(3,2) = (aa(2,1)*aa(1,3)-aa(1,1)*aa(2,3)) * dd
 ait(1,3) = (aa(2,1)*aa(3,2)-aa(3,1)*aa(2,2)) * dd
 ait(2,3) = (aa(3,1)*aa(1,2)-aa(1,1)*aa(3,2)) * dd
 ait(3,3) = (aa(1,1)*aa(2,2)-aa(2,1)*aa(1,2)) * dd

end subroutine matr3inv
!!***
