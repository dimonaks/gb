!{\src2tex{textfont=tt}}
!!****f* ABINIT/symq3
!! NAME
!! symq3
!!
!! FUNCTION
!! Determines the symmetry operations by which
!! a reciprocal vector q is preserved,
!! modulo a primitive reciprocal lattice vector,
!! and the time-reversal symmetry.
!!
!! COPYRIGHT
!! Copyright (C) 1999-2010 ABINIT group (GMR)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! nsym=number of space group symmetries
!! qpt(3)= vector in reciprocal space
!! symrec(3,3,nsym)=3x3 matrices of the group symmetries (reciprocal space)
!! [prtvol]=integer flag defining the verbosity of output. =0 if no output is provided.
!!
!! OUTPUT
!! symq(4,2,nsym)= (integer) three first numbers define the G vector ;
!!     fourth number is zero if the q-vector is not preserved, is 1 otherwise
!!     second index is one without time-reversal symmetry, two with time-reversal symmetry
!! timrev=1 if the time-reversal symmetry preserves the wavevector,
!!     modulo a reciprocal lattice vector.
!!
!! NOTES
!! The condition is :  $q =  O  S(q) - G$
!! with O being either the identity or the time reversal
!! symmetry (= inversion in reciprocal space)
!! and G being a primitive vector of the reciprocal lattice
!! If the time-reversal (alone) also preserves q, modulo
!! a lattice vector, then timrev is set to 1, otherwise 0.
!!
!! TODO
!! timrev is put to 1 only for Gamma.
!! Better handling should be provided in further version.
!!
!! PARENTS
!!      bands_classification,findggp,invars2m,rdddb9,read_gkk,respfn
!!      setup_little_group,symdyma
!!
!! CHILDREN
!!      wrap2_pmhalf,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine symq3(nsym,qpt,symq,symrec,timrev,prtvol)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
 use interfaces_14_hidewrite
 use interfaces_32_util, except_this_one => symq3
!End of the abilint section

 implicit none

!Arguments -------------------------------
!scalars
 integer,intent(in) :: nsym
 integer,intent(in),optional :: prtvol
 integer,intent(out) :: timrev
!arrays
 integer,intent(in) :: symrec(3,3,nsym)
 integer,intent(out) :: symq(4,2,nsym)
 real(dp),intent(in) :: qpt(3)

!Local variables -------------------------
!scalars
 integer :: ii,isign,isym,itirev,my_prtvol
 real(dp),parameter :: tol=2.d-8
 real(dp) :: reduce
 character(len=500) :: message
!arrays
 real(dp) :: difq(3),qsym(3),shift(3)

! *********************************************************************

!By default use verbose output!
 my_prtvol=1 ; if (PRESENT(prtvol)) my_prtvol=prtvol

 do isym=1,nsym
   do itirev=1,2
!    isign is 1 without time-reversal, -1 with time-reversal
     isign=3-2*itirev

!    Initialise the array symq
     do ii=1,3
       symq(ii,itirev,isym)=0
     end do

!    Get the symmetric of the vector
     do ii=1,3
       qsym(ii)=qpt(1)*isign*symrec(ii,1,isym)&
&       +qpt(2)*isign*symrec(ii,2,isym)&
&       +qpt(3)*isign*symrec(ii,3,isym)
     end do

!    Get the difference between the symmetric and the original vector

     symq(4,itirev,isym)=1
     do ii=1,3
       difq(ii)=qsym(ii)-qpt(ii)
!      Project modulo 1 in the interval ]-1/2,1/2] such that difq = reduce + shift
       call wrap2_pmhalf(difq(ii),reduce,shift(ii))
       if(abs(reduce)>tol)symq(4,itirev,isym)=0
     end do

!    If the operation succeded, change shift from real(dp) to integer, then exit loop
     if(symq(4,itirev,isym)/=0)then
       if (my_prtvol>0) then
         if(itirev==1)write(message,'(a,i4,a)')' symq3 : found symmetry',isym,' preserves q '
         if(itirev==2)write(message,'(a,i4,a)')' symq3 : found symmetry ',isym,' + TimeReversal preserves q '
         call wrtout(std_out,message,'COLL')
       end if
!      Uses the mathematical function NINT = nearest integer
       do ii=1,3
         symq(ii,itirev,isym)=nint(shift(ii))
       end do
     end if

!    DEBUG
!    write(6,*)isym,itirev
!    write(6,*)((symrec(ii,jj,isym),ii=1,3),jj=1,3)
!    write(6,*)qsym(1),qsym(2),qsym(3)
!    write(6,*)(symq(ii,itirev,isym),ii=1,4)
!    write(6,*)
!    ENDDEBUG

   end do !itirev
 end do !isym

!Test time-reversal symmetry
 timrev=1
 do ii=1,3
!  Unfortunately, this version does not work yet ...
!  call wrap2_pmhalf(2*qpt(ii),reduce,shift(ii))
!  if(abs(reduce)>tol)timrev=0
!  So, this is left ...
   if(abs(qpt(ii))>tol)timrev=0
 end do

 if(timrev==1.and.my_prtvol>0)then
   write(message, '(a,a,a)' )&
&   ' symq3 : able to use time-reversal symmetry. ',ch10,&
&   '  (except for gamma, not yet able to use time-reversal symmetry)'
   call wrtout(std_out,message,'COLL')
 end if

end subroutine symq3
!!***
