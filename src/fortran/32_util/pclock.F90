!{\src2tex{textfont=tt}}
!!****f* ABINIT/pclock
!! NAME
!! pclock
!!
!! FUNCTION
!! Print the timing at point number itimpt
!! if itimpt=0 then start the clock, else print the cpu time elapsed since clock started
!! To have the overall time use itimpt=9999 
!!
!! COPYRIGHT
!! Copyright (C) 1999-2010 ABINIT group (GMR, VO, LR, RWG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  itimpt=See function description.
!!  unit[optional]=the unit number for output 
!!  mode_paral[optional]=either "COLL" or "PERS"
!!
!! OUTPUT
!!  (only printing)
!!
!! PARENTS
!!      cexch_haydock,m_bz_mesh,rdm
!!
!! CHILDREN
!!      timein,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine pclock(itimpt,unit,mode_paral)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
 use interfaces_14_hidewrite
 use interfaces_18_timing
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: itimpt
 integer,intent(in),optional :: unit
 character(len=4),intent(in),optional :: mode_paral

!Local variables-------------------------------
!scalars
 integer :: unt
 real(dp),save :: cstart,wstart
 real(dp) :: cpu,wall
 character(len=4) :: mode
 character(len=500) :: msg

! *************************************************************************

 unt=std_out ; if (PRESENT(unit))        unt=unit
 mode='COLL' ; if (PRESENT(mode_paral)) mode=mode_paral

 if (itimpt==0) then
!  === Start clock ===
   call timein(cstart,wstart)
   cpu=zero ; wall=zero
 else
!  === Write elapsed CPU and wall time ===
   call timein(cpu,wall)
   cpu=cpu-cstart
   wall=wall-wstart
   if (itimpt==9999) then
     write(msg,'(2a,f13.1,a,f13.1)')ch10,&
&     '+Overall time at end (sec) : cpu= ',cpu,'  wall= ',wall
     call wrtout(unt,msg,mode)
   else
     write(msg,'(2a,i4,2a,f10.2,3a,f10.2,2a)')ch10,&
&     ' timing point number ',itimpt,ch10,&
&     '-cpu time  = ',cpu, ' seconds',ch10,&
&     '-real time = ',wall,' seconds',ch10
     call wrtout(unt,msg,mode)
   end if
 end if

end subroutine pclock
!!***
