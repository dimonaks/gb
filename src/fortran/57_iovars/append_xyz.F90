!{\src2tex{textfont=tt}}
!!****f* ABINIT/append_xyz
!! NAME
!! append_xyz
!!
!! FUNCTION
!! Translate the data from a xyz file (xyz_fname),
!! and add it at the end of the usual ABINIT input data string (string),
!! taking into account the dtset (dtset_char)
!!
!! COPYRIGHT
!! Copyright (C) 2002-2010 ABINIT group (MJV).
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  dtset_char*2=possible dtset label
!!  xyz_fname = name of the xyz file
!!  strln=maximal number of characters of string, as declared in the calling routine
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!  lenstr=actual number of characters in string
!!  string*(strln)=string of characters  (upper case) to which the xyz data are appended
!!
!! PARENTS
!!      importxyz
!!
!! CHILDREN
!!      leave_new,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine append_xyz(dtset_char,lenstr,string,xyz_fname,strln)

 use defs_basis
 use m_io_tools

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
 use interfaces_14_hidewrite
 use interfaces_16_hideleave
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: strln
 integer,intent(inout) :: lenstr
 character(len=2),intent(in) :: dtset_char
 character(len=fnlen),intent(in) :: xyz_fname
 character(len=strln),intent(inout) :: string

!Local variables-------------------------------
 character :: blank=' '
!scalars
 integer :: unitxyz, iatom, natom, mu
 integer :: lenstr_new
 integer :: lenstr_old
 character(len=5) :: string5
 character(len=20) :: string20
 character(len=500) :: message
!arrays
 real(dp),allocatable :: xangst(:,:)
 character(len=2),allocatable :: elementtype(:)

!************************************************************************

 lenstr_new=lenstr

!open file with xyz data
 unitxyz = get_unit()
 open (UNIT=unitxyz, file=trim(xyz_fname), status="unknown")
 write(message, '(3a)') &
& ' importxyz : Opened file ',trim(xyz_fname),'; content stored in string_xyz'
 call wrtout(std_out,message,'COLL')

!check number of atoms is correct
 read(unitxyz,*) natom
 
 write(string5,'(i5)')natom
 lenstr_old=lenstr_new
 lenstr_new=lenstr_new+7+len_trim(dtset_char)+1+5
 string(lenstr_old+1:lenstr_new)=" _NATOM"//trim(dtset_char)//blank//string5

 allocate(xangst(3,natom), elementtype(natom))

!read dummy line
 read(unitxyz,*)

!read atomic types and positions
 do iatom = 1, natom
   read(unitxyz,*) elementtype(iatom), xangst(:,iatom)
 end do
 close (unitxyz)

!Write the element types
 lenstr_old=lenstr_new
 lenstr_new=lenstr_new+7+len_trim(dtset_char)+1
 string(lenstr_old+1:lenstr_new)=" _TYPAX"//trim(dtset_char)//blank
 do iatom=1,natom
   lenstr_old=lenstr_new
   lenstr_new=lenstr_new+3
   string(lenstr_old+1:lenstr_new)=elementtype(iatom)//blank
 end do
 lenstr_old=lenstr_new
 lenstr_new=lenstr_new+3
 string(lenstr_old+1:lenstr_new)="XX " ! end card for TYPAX

!Write the coordinates
 lenstr_old=lenstr_new
 lenstr_new=lenstr_new+8+len_trim(dtset_char)+1
 string(lenstr_old+1:lenstr_new)=" _XANGST"//trim(dtset_char)//blank

 do iatom=1,natom
   do mu=1,3
     write(string20,'(f20.12)')xangst(mu,iatom)
     lenstr_old=lenstr_new
     lenstr_new=lenstr_new+20
     string(lenstr_old+1:lenstr_new)=string20
   end do
 end do

 deallocate(elementtype,xangst)

!Check the length of the string
 if(lenstr_new>strln)then
   write(message,'(6a)')ch10,&
&   ' append_xyz : BUG -',ch10,&
&   '  The maximal size of the input variable string has been exceeded.',ch10,&
&   '  The use of a xyz file is more character-consuming than the usual input file. Sorry.'
   call wrtout(std_out,message,'COLL')
   call leave_new('COLL')
 end if

!Update the length of the string
 lenstr=lenstr_new

end subroutine append_xyz
!!***

