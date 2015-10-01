!{\src2tex{textfont=tt}}
!!****f* ABINIT/iofn2
!! NAME
!! iofn2
!!
!! FUNCTION
!! First, read and echo pseudopotential filenames from unit 05.
!! Store them in an array.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2010 ABINIT group (DCA, XG, GMR, FrD, AF, MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!  npsp=number of pseudopotentials
!!
!! OUTPUT
!!  pspheads(npsp)=<type pspheader_type>=all the important information from the
!!   pseudopotential file headers, as well as the psp file names
!!
!! PARENTS
!!      m_ab6_invars_f90
!!
!! CHILDREN
!!      leave_new,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

 subroutine iofn2(npsp,filnam)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
 use interfaces_14_hidewrite
 use interfaces_16_hideleave
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: npsp
!arrays
 character(len=fnlen), intent(out) :: filnam(npsp)

!Local variables-------------------------------
!scalars
 integer :: ios,ipsp
 character(len=500) :: message
 character(len=fnlen) :: filpsp
!arrays

!*************************************************************************

 do ipsp=1,npsp
!  Read the name of the psp file
   write(6, '(/,a)' ) &
&   ' iofn2 : Please give name of formatted atomic psp file'
   read (5, '(a)' , iostat=ios ) filpsp
   filnam(ipsp)=trim(filpsp)
!  It might be that a file name is missing
   if(ios/=0)then
     write(message, '(a,a,a,a,a,a,a,a,a,a)' ) ch10,&
&     ' iofn2 : ERROR -',ch10,&
&     '  There are not enough names of pseudopotentials',ch10,&
&     '  provided in the files file.',ch10,&
&     '  Action : check first the variable ntypat (and/or npsp) in the input file;',ch10,&
&     '  if they are correct, complete your files file.'
     call wrtout(std_out,message,'COLL')
     call leave_new('COLL')
   end if

   write(6, '(a,i4,a,a)' ) &
&   ' iofn2 : for atom type',ipsp,' , psp file is ',trim(filpsp)

 end do ! ipsp=1,npsp
end subroutine iofn2
!!***
