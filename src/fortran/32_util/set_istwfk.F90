!{\src2tex{textfont=tt}}
!!****f* ABINIT/set_istwfk
!! NAME
!!  set_istwfk
!!
!! FUNCTION
!!  Returns the value of istwfk associated to the input k-point.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2010 ABINIT group
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUT 
!!  kpoint(3)=The k-point in reduced coordinates.
!!
!! OUTPUT
!!  istwfk= Integer flag internally used in the code to define the storage mode of the wavefunctions.
!!  It also define the algorithm used to apply an operator in reciprocal space as well as the FFT 
!!  algorithm used to go from G- to r-space
!!
!!   1 => time-reversal cannot be used
!!   2 => use time-reversal at the Gamma point.
!!   3 => use time-reversal symmetry for k=(1/2, 0 , 0 )
!!   4 => use time-reversal symmetry for k=( 0 , 0 ,1/2)
!!   5 => use time-reversal symmetry for k=(1/2, 0 ,1/2)
!!   6 => use time-reversal symmetry for k=( 0 ,1/2, 0 )
!!   7 => use time-reversal symmetry for k=(1/2,1/2, 0 )
!!   8 => use time-reversal symmetry for k=( 0 ,1/2,1/2)
!!   9 => use time-reversal symmetry for k=(1/2,1/2,1/2)
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

function set_istwfk(kpoint) result(istwfk)

 use defs_basis

 implicit none

!Arguments ------------------------------------
!scalars
 integer :: istwfk
!arrays
 real(dp),intent(in) :: kpoint(3)

!Local variables-------------------------------
!scalars
 integer :: bit0,ii
!arrays
 integer :: bit(3)

! *************************************************************************

 bit0=1
 
 do ii=1,3
   if (DABS(kpoint(ii))<tol10) then
     bit(ii)=0
   else if (DABS(kpoint(ii)-half)<tol10 ) then
     bit(ii)=1
   else
     bit0=0
   end if
 end do
!write(std_out,*)' bit0, bit(:)=',bit0,bit(1:3)

 if (bit0==0) then
   istwfk=1
 else
   istwfk=2+bit(1)+4*bit(2)+2*bit(3) ! Note the inversion between bit(2) and bit(3)
 end if

end function set_istwfk
!!***
