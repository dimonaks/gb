!{\src2tex{textfont=tt}}
!!****f* ABINIT/clp
!! NAME
!! clp
!!
!! FUNCTION
!! clp(x)= x-1, if x>1/2
!!         x+1, if x<-1/2
!!
!! COPYRIGHT
!! Copyright (C) 1998-2010 ABINIT group (GZ, MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!  x= input variable
!!
!! OUTPUT
!!  clp= resulting function
!!
!! PARENTS
!!      nhatgrid
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

function clp(x)

 use defs_basis

 implicit none

!Arguments ------------------------------------
!scalars
 real(dp) :: clp
 real(dp),intent(in) :: x

!Local variables-------------------------------

! **********************************************************************

 if(x > half) then
   clp=x-one
 elseif(x < -half) then
   clp=x+one
 else
   clp=x
 end if

 end function clp
!!***
