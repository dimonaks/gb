!{\src2tex{textfont=tt}}
!!****f* ABINIT/prttagm
!!
!! NAME
!! prttagm
!!
!! FUNCTION
!! Eventually print the content of dprarr (if typevarphys='DPR','LEN', and 'ENE'),
!! or intarr (if typevarphys='INT'), arrays of effective dimensions narr and 0:ndtset_alloc
!! For the second dimension, the 0 index relates to a default.
!! Print the array only if the content for at least one value of the second
!! index is different from the default.
!! Print a generic value if the non-default values are all equal.
!! Print the detail of all values otherwise.
!! The input variable 'length' controls the print format, and, in the case
!! of the real(dp) variable, the way two numbers are determined to be
!! different or not.
!!
!! COPYRIGHT
!! Copyright (C) 1999-2010 ABINIT group (XG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  intarr(1:marr,0:ndtset_alloc), dprarr(1:marr,0:ndtset_alloc)
!!   integer or real(dp) arrays, respectively,
!!   containing the data to be printed. Use these arrays even for scalars.
!!   For the first index, only the range 1:narr is relevant.
!!  iout=unit number for echoed output
!!  jdtset_(0:ndtset_alloc)=list of dataset indices.
!!  length= if 1, short format for printing, if 2, long format for printing
!!     special formats: if 3, INT : for symrel
!!                      if 4, INT : for type
!!                      if 5, INT : for mkmem, mkqmem, mk1mem
!!                      if 3, DPR : for tnons
!!                      if 4, DPR : for wtk and znucl
!!                      if 5, DPR : for atvshift
!!     If the typevarphys is 'DPR', a negative value of 'length' will request that
!!        the equality of real(dp) numbers is determined by an ABSOLUTE
!!        difference criterion only. The absolute value of length is used
!!        to determine the format, as above.
!!
!!  ndtset_alloc=govern second dimension of intarr and dprarr
!!  marr=first dimension of the intarr and dprarr arrays, as declared in the
!!   calling subroutine.
!!  narr=actual first dimension of intarr and dprarr.
!!  token=character string for 'tag'.  Assumed no longer than 9 characters
!!  typevarphys=physical variable type (might indicate the physical meaning of
!!   for dimensionality purposes)
!!   'INT'=>integer
!!   'DPR'=>real(dp) (no special treatment)
!!   'LEN'=>real(dp) (output in bohr and angstrom)
!!   'ENE'=>real(dp) (output in hartree and eV)
!!
!! OUTPUT
!!  (only writing)
!!
!! PARENTS
!!      outvar1,outvars,pawuj_det
!!
!! CHILDREN
!!      leave_new,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine prttagm(dprarr,intarr,iout,jdtset_,length,marr,narr,ndtset_alloc,token,typevarphys)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
 use interfaces_14_hidewrite
 use interfaces_16_hideleave
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: iout,length,marr,narr,ndtset_alloc
 character(len=*),intent(in) :: token
 character(len=3),intent(in) :: typevarphys
!arrays
 integer,intent(in) :: intarr(marr,0:ndtset_alloc),jdtset_(0:ndtset_alloc)
 real(dp),intent(in) :: dprarr(marr,0:ndtset_alloc)

!Local variables-------------------------------
 character(len=*), parameter :: short_int="(1x,a16,a,1x,(t20,10i5))"
 character(len=*), parameter :: long_int ="(1x,a16,a,1x,(t20,8i8) )"
 character(len=*), parameter :: form_dpr="(1x,a16,a,1x,(t20,"
 character(len=*), parameter :: short_dpr="es16.8))"
 character(len=*), parameter :: long_dpr="es18.10))"
 character(len=*), parameter :: short_dim="es16.8),a)"
 character(len=*), parameter :: long_dim="es18.10),a)"
 character(len=*), parameter :: f_symrel="(1x,a16,a,1x,(t20,3(3i3,1x),4x,3(3i3,1x)))"
 character(len=*), parameter :: f_type  ="(1x,a16,a,1x,(t20,20i3))"
 character(len=*), parameter :: f_mem   ="(a,a16,a,1x,(t20,8i8))"
 character(len=*), parameter :: f_tnons ="(1x,a16,a,1x,(t20,3f11.7,3x,3f11.7))"
 character(len=*), parameter :: f_wtk   ="(1x,a16,a,1x,(t20,6f11.5))"
 character(len=*), parameter :: f_atvshift   ="(1x,a16,a,1x,(t20,5f11.5))"
!scalars
 integer :: iarr,idtset,jdtset,multi,ndtset_eff,print
 real(dp),parameter :: tol21=1.0d-21
 real(dp) :: diff,sum
 character(len=1) :: digit
 character(len=2) :: appen
 character(len=50) :: format_dp
 character(len=500) :: message

! *************************************************************************

!DEBUG
!write(*,*)'prttagm start'
!END DEBUG

 if(len_trim(token)>16)then
   write(message, '(a,a,a,a,i6,a,a,a,a,a)' ) ch10,&
&   ' prttagm : BUG -',ch10,&
&   '  The length of the name of the input variable ',trim(token),' is ',len_trim(token),ch10,&
&   '  This exceeds 16 characters, the present maximum in routine prttagm.'
   call wrtout(std_out,message,'COLL')
   call leave_new('COLL')
 end if

 if(ndtset_alloc<1)then
   write(message, '(a,a,a,a,i6,a,a,a,a,a)' ) ch10,&
&   ' prttagm : BUG -',ch10,&
&   '  ndtset_alloc=',ndtset_alloc,', while it should be >= 1.',ch10,&
&   '  This happened for token=',token,'.'
   call wrtout(std_out,message,'COLL')
   call leave_new('COLL')
 end if

 if(ndtset_alloc>99)then
   write(message, '(a,a,a,a,i6,a,a,a,a,a)' ) ch10,&
&   ' prttagm : BUG -',ch10,&
&   '  ndtset_alloc=',ndtset_alloc,', while it must be lower than 100.',ch10,&
&   '  This happened for token=',token,'.'
   call wrtout(std_out,message,'COLL')
   call leave_new('COLL')
 end if

 if(narr>99 .and. (typevarphys=='ENE'.or.typevarphys=='LEN'))then
   write(message, '(6a,i6,a)' ) ch10,&
&   ' prttagm : BUG -',ch10,&
&   '  typevarphys=',typevarphys,' with narr=',narr,'  is not allowed.'
   call wrtout(std_out,message,'COLL')
   call leave_new('COLL')
 end if

 if (narr>0) then

   print=1
   multi=0

!  Treat integer first
   if(typevarphys=='INT')then

!    Determine whether the different non-default occurences are all equal
     if(ndtset_alloc>1)then
       do idtset=1,ndtset_alloc
         do iarr=1,narr
           if(intarr(iarr,1)/=intarr(iarr,idtset))multi=1
         end do
       end do
     end if
!    If they are all equal, then determine whether they are equal to the default
     if(multi==0)then
       print=0
       do iarr=1,narr
         if(intarr(iarr,1)/=intarr(iarr,0))print=1
       end do
     end if
!    Print only if the values differ from the default
     if(print==1)then
       ndtset_eff=ndtset_alloc
       if(multi==0)ndtset_eff=1
       do idtset=1,ndtset_eff
         if(multi==0)then
           appen=' '
         else
           jdtset=jdtset_(idtset)
           if(jdtset<10)write(appen,'(i1)')jdtset
           if(jdtset>=10)write(appen,'(i2)')jdtset
         end if
         if(abs(length)==1)write(iout,short_int)token,appen,intarr(1:narr,idtset)
         if(abs(length)==2)write(iout,long_int)token,appen,intarr(1:narr,idtset)
         if(abs(length)==3)write(iout,f_symrel)token,appen,intarr(1:narr,idtset)
         if(abs(length)==4)write(iout,f_type)token,appen,intarr(1:narr,idtset)
         if(abs(length)==5)write(iout,f_mem)'P',token,appen,intarr(1:narr,idtset)
       end do
     end if

!    Treat now real(dp) : same structure as for integer numbers.
   else if (typevarphys=='DPR' .or. typevarphys=='LEN' .or. typevarphys=='ENE') then

     if(ndtset_alloc>1)then
       do idtset=1,ndtset_alloc
         do iarr=1,narr
!          The determination of effective equality is more difficult than in the
!          integer case :
!          - if length > 0, ask for a relative accuracy, and also include
!          the case of zero values, thanks to tol21.
!          - if length < 0, ask for absolute accuracy.
           diff=abs( dprarr(iarr,1)-dprarr(iarr,idtset) )
           if(length>0)then
             sum=abs(dprarr(iarr,1))+abs(dprarr(iarr,idtset))+10*tol21
             if(diff>sum*tol11)multi=1
           else
             if(diff>tol14)multi=1
           end if
         end do
       end do
     end if
     if(multi==0)then
       print=0
       do iarr=1,narr
         diff=abs( dprarr(iarr,1)-dprarr(iarr,0) )
         if(length>0)then
           sum=abs(dprarr(iarr,1))+abs(dprarr(iarr,0))+10*tol21
           if(diff>sum*tol11)print=1
         else
           if(diff>tol14)print=1
         end if
       end do
     end if
     if(print==1)then
!      Select the proper format
       if(abs(length)==1 .or. abs(length)==2)then
         if(typevarphys=='DPR')then
           digit='3'
           if(abs(length)==1)format_dp=form_dpr//digit//short_dpr
           if(abs(length)==2)format_dp=form_dpr//digit//long_dpr
         else if(typevarphys=='ENE' .or. typevarphys=='LEN')then
           if (narr<10) write(digit,'(i1)')narr
           if (narr> 9) write(digit,'(i2)')narr
           if(abs(length)==1)format_dp=form_dpr//digit//short_dim
           if(abs(length)==2)format_dp=form_dpr//digit//long_dim
         end if
       else
         if(abs(length)==3)format_dp=f_tnons
         if(abs(length)==4)format_dp=f_wtk
         if(abs(length)==5)format_dp=f_atvshift
       end if
       ndtset_eff=ndtset_alloc
       if(multi==0)ndtset_eff=1
       do idtset=1,ndtset_eff
         if(multi==0)then
           appen=' '
         else
           jdtset=jdtset_(idtset)
           if(jdtset<10)write(appen,'(i1)')jdtset
           if(jdtset>=10)write(appen,'(i2)')jdtset
         end if
!        DEBUG
!        if(typevarphys=='ENE')then
!        write(6,*)format_dp
!        stop
!        end if
!        ENDDEBUG
         if(typevarphys=='DPR')write(iout,format_dp)token,appen,dprarr(1:narr,idtset)
         if(typevarphys=='ENE')write(iout,format_dp)token,appen,dprarr(1:narr,idtset),' Hartree'
         if(typevarphys=='LEN')write(iout,format_dp)token,appen,dprarr(1:narr,idtset),' Bohr'
       end do
     end if

!    The type is neither 'INT' nor 'DPR','ENE','LEN'
   else

     write(message, '(a,a,a,a,a,a)' ) ch10,&
&     ' prttagm : BUG -',ch10,&
&     '  Disallowed typevarphys=',typevarphys,'.'
     call wrtout(std_out,message,'COLL')
     call leave_new('COLL')

   end if

!  End condition of narr>0
 end if

end subroutine prttagm
!!***
