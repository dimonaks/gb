!{\src2tex{textfont=tt}}
!!****f* ABINIT/inupper
!! NAME
!! inupper
!!
!! FUNCTION
!! Maps all characters in string to uppercase.
!! Uses fortran90 character string manipulation but should work
!! independent of EBCDIC or ASCII assumptions--only relies on
!! 'index' intrinsic character string matching function.
!! Makes sure that the string 'lolett' remains defined as the lower
!! case 26-character alphabet string and 'uplett' remains upper case.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2010 ABINIT group (DCA, XG, GMR).
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  string= character string with arbitrary case
!!
!! OUTPUT
!!  string= same character string mapped to upper case
!!
!! SIDE EFFECTS
!!  string= (input) character string with arbitrary case
!!          (output) same character string mapped to upper case
!!
!! PARENTS
!!      anaddb,chkexi,chkvars,intagm,invars1,localorb_S,lwf,newsp,parsefile
!!      prt_cml,testlda
!!
!! CHILDREN
!!      leave_new,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine inupper(string)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 character(len=*),intent(inout) :: string

!Local variables-------------------------------
!scalars
 integer :: ii,indx,stringlen
 logical,save :: first=.true.
 character(len=1) :: cc
 character(len=500) :: message
!no_abirules
 character(len=26), parameter :: uplett='ABCDEFGHIJKLMNOPQRSTUVWXYZ',&
&                           lolett='abcdefghijklmnopqrstuvwxyz'

! *************************************************************************
!
!On first entry make sure lower case letters stayed
!lower case and upper case letters stayed upper case
 if (first) then
   do ii=1,26
!    Look for occurrence of each upper case character
!    anywhere in string of all lower case letters
     indx=index(lolett,uplett(ii:ii))
!    If found then print error message and quit
     if (indx>0) then
       write(std_out, '(a,a,a,a,a,a,a,a,a,a,a,a)' ) ch10,&
&       ' inupper: BUG -',ch10,&
&       '  Upper case string=',uplett,ch10,&
&       '  Lower case string=',lolett,ch10,&
&       '  Upper case character ',uplett(ii:ii),&
&       ' found in supposedly lower case string.'
!       call wrtout(std_out,message,'COLL')
!       call leave_new('COLL')
     end if
   end do
   first=.false.
 end if
!
 stringlen=len(string)
 do ii=1,stringlen
!  Pick off single character of string (one byte):
   cc=string(ii:ii)
!  determine whether a lowercase letter:
   indx=index(lolett,cc)
   if (indx>0) then
!    Map to uppercase:
     string(ii:ii)=uplett(indx:indx)
   end if
 end do

end subroutine inupper
!!***
