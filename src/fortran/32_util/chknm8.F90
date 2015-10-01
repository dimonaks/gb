!{\src2tex{textfont=tt}}
!!****f* ABINIT/chknm8
!! NAME chknm8
!! chknm8
!!
!!
!! FUNCTION
!! This small subroutine check the identity of its argument,
!! who are a6 names, and eventually send a message and stop
!! if they are found unequal
!!
!! COPYRIGHT
!! Copyright (C) 1999-2010 ABINIT group (XG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! nmfond= name which has to be checked
!! nmxpct= name expected for nmfond
!!
!! OUTPUT
!!
!!
!! SIDE EFFECTS
!!
!!
!! NOTES
!!
!! TODO
!! Describe the inputs
!!
!! PARENTS
!!      inprep8,ioddb8
!!
!! CHILDREN
!!      leave_new,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine chknm8(nmxpct,nmfond)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
 use interfaces_14_hidewrite
 use interfaces_16_hideleave
!End of the abilint section

 implicit none

!Arguments -------------------------------
!scalars
 character(len=9),intent(in) :: nmfond,nmxpct

!Local variables-------------------------------
!scalars
 character(len=500) :: message

! *********************************************************************

 if(nmxpct/=nmfond) then
   write(message, '(a,a,a,a,a,a,a,a,a,a,a,a,a)' )&
&   ' chknm8 : ERROR -',ch10,&
&   '  Reading DDB, expected name was "',nmxpct,'"',ch10,&
&   '               and name found is "',nmfond,'.',ch10,&
&   '  Likely your DDB is incorrect.',ch10,&
&   '  Action : correct your DDB, or contact the ABINIT group.'
   call wrtout(std_out,message,'COLL')
   call leave_new('COLL')
 end if

end subroutine chknm8
!!***
