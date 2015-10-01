!{\src2tex{textfont=tt}}
!!****f* ABINIT/prmat
!!
!! NAME
!! prmat
!!
!! FUNCTION
!! This subroutine prints real*8 matrices in an attractive format.
!!
!! COPYRIGHT
!! Copyright (C) 1987-2010 ABINIT group (ZL, DCA, XG, GMR)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  mat(mi,nj)= matrix to be printed
!!  mi        = no rows of mat
!!  ni        = no rows to print
!!  nj        = no colums of mat
!!  unitm     = unit to print to, if not provided std_out is chosen
!!
!! OUTPUT
!!  (only writing)
!!
!! PARENTS
!!      newkpt
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine prmat (mat, ni, nj, mi, unitm)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in)           :: mi,ni,nj
 integer,intent(in), optional :: unitm
!arrays
 real(dp),intent(in)          :: mat(mi,nj)

!Local variables-------------------------------
!scalars
 character(len=1000)    :: message
 integer,parameter      :: nline=10
 integer                :: ii,jj,jstart,jstop,unitn

! *************************************************************************

 if (present(unitm)) then ! standard printing to std_out
   unitn=unitm
 else
   unitn=std_out
 end if


 do  jstart = 1, nj, nline
   jstop = min(nj, jstart+nline-1)
   write(message, '(3x,10(i4,8x))' ) (jj,jj=jstart,jstop)
   call wrtout(unitn,message,'COLL') 
 end do

 do ii = 1,ni
   do jstart= 1, nj, nline
     jstop = min(nj, jstart+nline-1)
     if (jstart==1) then
       write(message, '(i3,1p,10e12.4)' ) ii, (mat(ii,jj),jj=jstart,jstop)
       call wrtout(unitn,message,'COLL')
     else
       write(message, '(3x,1p,10e12.4)' )    (mat(ii,jj),jj=jstart,jstop)
       call wrtout(unitn,message,'COLL')
     end if
   end do
 end do
!

end subroutine prmat
!!***
