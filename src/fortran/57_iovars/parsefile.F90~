!!****f* ABINIT/parsefile
!! NAME
!! parsefile
!!
!! FUNCTION
!!  Glue function, to read the given file, put it into a string,
!!  change everything to uppercase, remove carriage returns and
!!  non significant blank characters. May also read a XML input
!!  CML file if specified. Finaly read ndtset input variable.
!!
!! INPUTS
!!  filnamin= the file to read
!!
!! OUTPUT
!!  lenstr= the length of the resulting string.
!!  ndtset= the number of declared datasets.
!!  string= contains on output the content of the file, ready for parsing.
!!
!! SIDE EFFECTS
!!
!! PARENTS
!!      abinit,m_ab6_invars_f90,ujdet
!!
!! CHILDREN
!!      importcml,importxyz,instrng,intagm,inupper,leave_new,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine parsefile(filnamin, lenstr, ndtset, string)

  use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
 use interfaces_14_hidewrite
 use interfaces_16_hideleave
 use interfaces_32_util
 use interfaces_42_parser
 use interfaces_47_xml
 use interfaces_57_iovars, except_this_one => parsefile
!End of the abilint section

  implicit none

!Arguments ------------------------------------

  character(len = *), intent(in) :: filnamin
  integer, intent(out) :: ndtset, lenstr
  character(len = strlen), intent(out) :: string

!Local variables-------------------------------
  integer :: option, marr, tread
  character(len = strlen) :: string_raw
  character(len = 30) :: token
  integer :: intarr(1)
  real(dp) :: dprarr(1)
  character(len=500) :: message

! *************************************************************************

!Read the input file, and store the information in a long string of characters
!Note : this could be done only by me=0, and then string would be BCAST
!!$ !!!!!!!!!!!!!!!!!!!!!!
!!$ DC 2007-02-06: WARNING, this preprocessing option has been removed since
!!$ the variable netcdf_input is not initialized anymore (previously done
!!$ in iofn1_ncdf() that had disapeared).
!!$
!!$#if defined HAVE_NETCDF
!!$ !If the input is in a NetCDF file, a string should be rebuilt to have the
!!$ !parser believe it was plain text.
!!$ if ( netcdf_input ) then
!!$  write(message,'(a,a,a,a)') ch10,&
!!$&  ' abinit : ERROR -',ch10,&
!!$&  ' Restarting from a NetCDF file has not been implemented yet.'
!!$  call wrtout(std_out,message,'COLL')
!!$  call leave_new('COLL')
!!$ else
!!$#endif
!!$ !!!!!!!!!!!!!!!!!!!!!!

!strlen from defs_basis module
 option=1
 call instrng (filnamin,lenstr,option,strlen,string)

!Copy original file, without change of case
 string_raw=string

!To make case-insensitive, map characters of string to upper case:
 call inupper(string(1:lenstr))

!Might import data from CML file(s) into string
!Need string_raw to deal properly with CML filenames
! call importcml(lenstr,string_raw,string,strlen)
! call importxyz(lenstr,string_raw,string,strlen)

!!$ !!!!!!!!!!!!!!!!!!!!!!
!!$#if defined HAVE_NETCDF
!!$ end if
!!$#endif
!!$ !!!!!!!!!!!!!!!!!!!!!!
 
!6) Take ndtset from the input string, then allocate
!the arrays whose dimensions depends only on ndtset and msym.
 


end subroutine parsefile
!!***
