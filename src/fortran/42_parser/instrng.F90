!{\src2tex{textfont=tt}}
!!****f* ABINIT/instrng
!! NAME
!! instrng
!!
!! FUNCTION
!! Read the input file, and product a string of character,
!! with all data, to be analyzed in later routines. The length
!! of this string is lenstr. This number is checked to be smaller
!! than the dimension of the string of character, namely strln .
!!
!! COPYRIGHT
!! Copyright (C) 1998-2010 ABINIT group (DCA, XG, GMR).
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  filnam=name of the input file, to be read
!!  option= if 0, simple storing of the character string,
!!                 no special treatment for ABINIT (comment delimiters, checks ...)
!!          if 1, suppresses text after an ABINIT comment delimiter (! or #),
!!                 checks that a minus sign is followed by a number ...
!!  strln=maximal number of character of string, as declared in the calling routine
!!
!! OUTPUT
!!  lenstr=actual number of character in string
!!  string*(strln)=string of character
!!
!! PARENTS
!!      abinit,anaddb,importcml,localorb_S,lwf,newsp
!!
!! CHILDREN
!!      incomprs,leave_new,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine instrng (filnam,lenstr,option,strln,string)!,len_filnamin,len_string)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
! use interfaces_14_hidewrite
! use interfaces_16_hideleave
! use interfaces_42_parser, except_this_one => instrng
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: option,strln!,len_filnamin,len_string
 integer,intent(out) :: lenstr
 character(len=*),intent(in) :: filnam
 character(len=*),intent(out) :: string

!Local variables-------------------------------
 character :: blank=' '
!scalars
 integer :: ii,ii1,ii2,ij,iline,ios,lenc,mline,nline1
 character(len=1) :: string1
 character(len=3) :: string3
 character(len=500) :: message
 character(len=fnlen+20) :: line

!************************************************************************

!DEBUG
!write(6,*)' instrng : enter '
!write(6,*)' filnam =',filnam

!write(6,*)' filnam =',trim(filnam)
!stop
!ENDDEBUG

!%%%%%%%%%%%%%%%%%%%%%%%%
!read in string from file
!%%%%%%%%%%%%%%%%%%%%%%%%

!Open data file and read one line at a time, compressing data
!and concatenating into single string:
 open (unit=tmp_unit,file=trim(filnam),status='old',form='formatted')
 rewind (unit=tmp_unit)

!Initialize string to blanks
 string=blank

 lenstr=1

!Set maximum number lines to be read to some large number
 mline=50000
 do iline=1,mline

!  Keeps reading lines until end of input file
   read (unit=tmp_unit,fmt= '(a)' ,iostat=ios) line(1:fnlen+20)
!  Hello ! This is a commentary. Please, do not remove me.
!  In fact, this commentary protect tests_v4 t47 for miscopying
!  the input file into the output string. It _is_ strange.
!  The number of lines in the commentary is also resulting from
!  a long tuning..

!  DEBUG
!  write(6,*)' instrng, ios=',ios,' echo :',trim(line(1:fnlen+20))
!  ENDDEBUG

!  Exit the reading loop when arrived at the end
   if(ios/=0)then
     backspace(tmp_unit)
     read (unit=tmp_unit,fmt= '(a1)' ,iostat=ios) string1
     if(ios/=0)exit
     backspace(tmp_unit)
     read (unit=tmp_unit,fmt= '(a3)' ,iostat=ios) string3
     if(string3=='end')exit
     write(std_out, '(4a,i4,11a)' ) ch10,&
&     ' instrng : ERROR - ',ch10,&
&     '  It is observed in the input file, line number',&
&     iline,',',ch10,' that there is a non-zero IO signal.',ch10,&
&     ' This is normal when the file is completely read.',ch10,&
&     ' However, it seems that the error appears while your file has not been completely read.',ch10,&
&     ' Action : correct your file. If your file seems correct, then,',ch10,&
&     ' add the keyword ''end'' at the very beginning of the last line of your input file.'
!     call wrtout(std_out,message,'COLL')
!     call leave_new('COLL')
   end if

!  Find length of input line ignoring delimiter characters (# or !)
!  and any characters beyond it (allows for comments beyond # or !)
   ii1=index(line(1:fnlen+20),'#')
   ii2=index(line(1:fnlen+20),'!')
   if ( (ii1==0 .and. ii2==0) .or. option==0 ) then
!    delimiter character was not found on line so use full line
     ii=fnlen+20
   else if(ii1==0)then
!    ii will represent length of line up to but not including !
     ii=ii2-1
   else if(ii2==0)then
!    ii will represent length of line up to but not including #
     ii=ii1-1
   else
     ii=min(ii1,ii2)-1
   end if

!  Checks that nothing is left beyond fnlen
   if(ii>fnlen)then
     do ij=fnlen+1,ii
       if(line(ij:ij)/=' ')then
         write(std_out, '(a,a,a,a,i4,a,a,a,a,a)' ) ch10,&
&         ' instrng : ERROR - ',ch10,&
&         '  It is observed in the input file, line number',&
&         iline,',',ch10,' that more than 132 columns are used.',ch10,&
&         ' This is not allowed. Change this line of your input file.'
!         call wrtout(std_out,message,'COLL')
!         call leave_new('COLL')
       end if
     end do
   end if

   if (ii>0) then
!    Check for the occurence of a minus sign followed by a blank
     ij=index(line(1:ii),'- ')
     if (ij>0 .and. option==1) then
       write(std_out, '(a,a,a,a,i4,a,a,a,a,a,a,a,a,a,a,a)' ) ch10,&
&       ' instrng : ERROR - ',ch10,&
&       '  It is observed in the input file, line number',&
&       iline,',',ch10,'  the occurence of a minus sign followed',ch10,&
&       '  by a blank. This is forbidden.',ch10,&
&       '  If the minus sign is meaningful, do not leave a blank',ch10,&
&       '  between it and the number to which it applies.',ch10,&
&       '  Otherwise, remove it.'
!       call wrtout(std_out,message,'COLL')
!       call leave_new('COLL')
     end if
!    Check for the occurence of a tab
     ij=index(line(1:ii),char(9))
     if (ij>0 .and. option==1 ) then
       write(std_out, '(a,a,a,a,i4,a,a,a)' ) ch10,&
&       ' instrng : ERROR - ',ch10,&
&       '  The occurence of a tab, in the input file, line number',&
&       iline,',',ch10,&
&       '  is observed. This sign is confusing, and has been forbidden.'
!       call wrtout(std_out,message,'COLL')
!       call leave_new('COLL')
     end if
!    Compress: remove repeated blanks, make all ASCII characters
!    less than a blank (and '=') to become a blank.
     call incomprs(line(1:ii),lenc)

   else
!    ii=0 means line starts with #, is entirely a comment line
     lenc=0
   end if

!  Check resulting total string length
   if (lenstr+lenc>strln) then
     write(std_out, '(a,a,a,a,a,a,a,a,a,a)' ) ch10,&
&     ' instrng : ERROR --',ch10,&
&     '  The size of your input file is such that the internal',ch10,&
&     '  character string that should contain it is too small.',ch10,&
&     '  Action : decrease the size of your input file,',ch10,&
&     '  or contact the ABINIT group.'
!     call wrtout(std_out,message,'COLL')
!     call leave_new('COLL')
   end if

   if (lenc>0) then
!    Concatenate new compressed characters
!    with previous part of compressed string (unless all blank)
     string(lenstr+1:lenstr+lenc)=line(1:lenc)
   end if
!  Keep track of total string length
   lenstr=lenstr+lenc

!  If mline is reached, something is wrong
   if (iline>=mline) then
     write(std_out, '(a,a,a,a,i10,a,a,a,i10,a,a,a,a)' ) ch10,&
&     ' instrng : ERROR -',ch10,&
&     '  The number of lines already read from input file=',&
&     iline,ch10,' is equal or greater than maximum allowed',&
&     '  mline=',mline,ch10,&
&     '  Action : you could decrease the length of the input file, or',ch10,&
&     '  contact the ABINIT group.'
!     call wrtout(std_out,message,'COLL')
!     call leave_new('COLL')
   end if

!  End loop on iline. Note that there is an "exit" instruction in the loop
 end do

 nline1=iline-1
 close (unit=tmp_unit)

 write(std_out, '(a,i6,a)' ) &
& ' instrng :',nline1,' lines of input have been read'
! call wrtout(std_out,message,'COLL')

!DEBUG
!write(6,*)' instrng : string =',trim(string)
!ENDDEBUG


end subroutine instrng
!!***
