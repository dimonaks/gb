!{\src2tex{textfont=tt}}
!!****f* ABINIT/appdig
!! NAME
!! appdig
!!
!! FUNCTION
!! Using input string "string" and integer "integ", make a string
!! named 'strinn' by concatenating digits of "integ" with characters
!! of "string"; return final string in "strinn".
!!
!! COPYRIGHT
!! Copyright (C) 1998-2010 ABINIT group (DCA, XG, GMR)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! integ=nonnegative integer whose digits will be appended to string
!! string=string to which digits will be appended
!!
!! OUTPUT
!! strinn=string//nn
!!
!! PARENTS
!!      berryphase_new,loop3dte,loper3,nstdy3,scfcv3,uderiv
!!
!! CHILDREN
!!      leave_new
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine appdig(integ,string,strinn)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
 use interfaces_16_hideleave
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: integ
 character(len=*),intent(in) :: string
 character(len=*),intent(out) :: strinn

!Local variables-------------------------------
!scalars
 integer :: i,length,ndig
 character(len=2) :: ncha
 character(len=8) :: form

! *************************************************************************
!
!Check that integer is nonnegative
 if (integ<0) then
   write(06, '(/,a,/,a,i10,a,/,a)' ) &
&   ' appdig: BUG -',&
&   '  Input integer=',integ,' must not be <0.',&
&   '  Argument integ was input as negative.'
   call leave_new('COLL')
 end if
!
!Fill output string initially with blanks to end of dimensioned length
 length=len(strinn)
 do i=1,length
   strinn(i:i)=' '
 end do

!Find nonwhitespace length of string
 length=len_trim(string)
!Copy input character string into first part of output string
 strinn(1:length)=string(1:length)

!Find how many digits "integ" has
 ndig=int(log10(real(integ)+0.50))+1

!Create a format for exact number of digits using internal write
 write(unit=ncha,fmt='(i2)') ndig
 form='(i'//ncha//')'
!Do internal write to get digits of integer into character string,
!placing digits into appropriate end of string.
 write(unit=strinn(1+length:1+length+ndig),fmt=form) integ
!(Note that present version writes "1" or "2" for single digit,
!not "01" or "02".  Latter may be preferable.  Can be amended.)
!
end subroutine appdig
!!***
