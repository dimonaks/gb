!{\src2tex{textfont=tt}}
!!****f* ABINIT/status
!! NAME
!! status
!!
!! FUNCTION
!! Routine for description of the status of the calculation
!! Eventually open the status file, write different information,
!! and close the file. The output rate and shift are governed by istat
!!
!! COPYRIGHT
!! Copyright (C) 1998-2010 ABINIT group (XG,TD)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  counter=value of the loop counter at that level
!!  file (optional argument)=name of the status file
!!  istat=gives the rate or shift of output. The status file will be opened
!!     and written only once every "istatr" calls.
!!     This variable is saved at the fourth call (just after the first
!!     call to invars0. The shift "istatshft" is saved at the fifth call.
!!     In preceeding and subsequent calls, istat has no meaning.
!!  level=number of the level of the calling subroutine (see the description later)
!!  routine=string of 14 characters indicating the status inside the level
!!
!! OUTPUT
!!  (only writing)
!!
!! NOTES
!! Warning : The string "routine" can have any size but
!! it is truncated to a size of 14.
!! because of the behaviour of some compilers, the string
!! "routine" should always have 14 characters in the calling subroutine
!!
!! PARENTS
!!      abinit,accrho3,bethe_salpeter,brdmin,calc_sigc_me,calc_sigx_me,cgwf
!!      cgwf3,delocint,driver,getgh1c,getghc,gstate,gstateimg,loop3dte,loper3
!!      m_ab6_invars_f90,moldyn,move,mv_3dte,nonlinear,resp3dte,respfn
!!      rhofermi3,scfcv,scfcv3,screening,sigma,suscep,vtorho,vtorho3,vtorhotf
!!      vtowfk,vtowfk3,wfkfermi3
!!
!! CHILDREN
!!      leave_new,timab,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine status(counter,filstat,istat,level,routine)

 use defs_basis
 use m_errors

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
 use interfaces_14_hidewrite
 use interfaces_16_hideleave
 use interfaces_18_timing
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: counter,istat,level
 character(len=*),intent(in) :: routine
 character(len=fnlen),intent(in) :: filstat

!Local variables-------------------------------
!scalars
 integer,parameter :: mcounter=2,mlevel=120
 integer,save :: first=1,output_rate=1,shift_rate=1,statnu=0
 integer :: ilevel,ios
 character(len=12) :: headwr
 character(len=500) :: message
!arrays
 integer,save :: active(mlevel),actual_counter(mlevel,mcounter)
 integer,save :: ncounter(mlevel)
 integer,save :: list_level(30)=&
&  (/1,2,100,101,102,103,104,110,111,112,113,114,10,11,12,13,14,15,16,17,18,20,30,31,40,41,50,51,52,90/)
 real(dp) :: tsec(2)
 character(len=14),save :: nm_levels(mlevel),nm_routine(mlevel)
 character(len=12),save :: nm_counter(mlevel,mcounter)

!***********************************************************************

 call timab(73,1,tsec)

!Note : all processors have their own file, so no special
!attention must be paid to the parallel case.
!Initialisation
 if(first/=0)then
   if(first==1)then
     first=4
     nm_routine(:)='              '
     active(:)=0
     actual_counter(:,:)=0

!    List of names for each level
!    Numbers from 1 to 9 are for abinit and driver
!    Numbers from 100 to 120 are for optdriver=0 routines (GS)
!    Numbers between 10 and 19 are for optdriver=1 routines (RF)
!    Numbers between 20 and 29 are for optdriver=2 routines (suscep)
!    Numbers between 30 and 39 are for optdriver=3 routines (screening)
!    Numbers between 40 and 49 are for optdriver=4 routines (sigma)
!    Numbers between 50 and 59 are for optdriver=5 routines (nonlinear)
!    When you add a level number, or modify one, do not forget to change list_level

     nm_levels(1)   ='abinit        '
     ncounter(1)=0
     nm_counter(1,1)='            '

     nm_levels(2)   ='driver        '
     ncounter(2)=1
     nm_counter(2,1)='jdtset     ='


!    Optdriver=0
     nm_levels(100)   ='gstateimg     '
     ncounter(100)=2
     nm_counter(100,1)='itimimage  ='
     nm_counter(100,2)='idynimage  ='

     nm_levels(101)   ='gstate        '
     ncounter(101)=1
     nm_counter(101,1)='itime      ='

     nm_levels(102)   ='move          '
     ncounter(102)=2
     nm_counter(102,1)='icalls     ='
     nm_counter(102,2)='itime      ='

     nm_levels(103)   ='brdmin/moldyn '
     ncounter(103)=1
     nm_counter(103,1)='itime      ='

     nm_levels(104)   ='delocint      '
     ncounter(104)=1
     nm_counter(104,1)='itime      ='

     nm_levels(110)   ='scfcv         '
     ncounter(110)=1
     nm_counter(110,1)='istep      ='

     nm_levels(111)   ='vtorho(tf)    '
     ncounter(111)=2
     nm_counter(111,1)='isppol     ='
     nm_counter(111,2)='ikpt       ='

     nm_levels(112)   ='vtowfk        '
     ncounter(112)=2
     nm_counter(112,1)='inonsc     ='
     nm_counter(112,2)='iband      ='

     nm_levels(113)   ='cgwf          '
     ncounter(113)=1
     nm_counter(113,1)='iline      ='

     nm_levels(114)   ='getghc        '
     ncounter(114)=0
     nm_counter(114,1)='            '


!    Optdriver=1
     nm_levels(10)   ='respfn        '
     ncounter(10)=0
     nm_counter(10,1)='            '

     nm_levels(11)   ='loper3        '
     ncounter(11)=1
     nm_counter(11,1)='respcase   ='

     nm_levels(12)   ='scfcv3        '
     ncounter(12)=1
     nm_counter(12,1)='istep      ='

     nm_levels(13)   ='vtorho3       '
     ncounter(13)=2
     nm_counter(13,1)='isppol     ='
     nm_counter(13,2)='ikpt       ='

     nm_levels(14)   ='vtowfk3       '
     ncounter(14)=2
     nm_counter(14,1)='inonsc     ='
     nm_counter(14,2)='iband      ='

     nm_levels(15)   ='cgwf3         '
     ncounter(15)=1
     nm_counter(15,1)='iline      ='

     nm_levels(16)   ='getgh1c       '
     ncounter(16)=0
     nm_counter(16,1)='            '

     nm_levels(17)   ='rhofermi3     '
     ncounter(17)=2
     nm_counter(17,1)='isppol     ='
     nm_counter(17,2)='ikpt       ='

     nm_levels(18)   ='wfkfermi3     '
     ncounter(18)=2
     nm_counter(18,1)='inonsc     ='
     nm_counter(18,2)='iband      ='


!    Optdriver=2
     nm_levels(20)   ='suscep        '
     ncounter(20)=0
     nm_counter(20,1)='            '


!    Optdriver=3
     nm_levels(30)   ='screening     '
     ncounter(30)=1
     nm_counter(30,1)='iqpt       ='

     nm_levels(31)   ='cchi0/cchi0q0 '
     ncounter(31)=1
!    nm_counter(31,1)='isppol  ='
     nm_counter(31,1)='ikpt       ='


!    Optdriver=4
     nm_levels(40)   ='sigma         '
     ncounter(40)=1
     nm_counter(40,1)='ikpt_gw    ='

     nm_levels(41)   ='csigme        '
     ncounter(41)=1
!    nm_counter(41,1)='isppol     ='
     nm_counter(41,1)='ikpt       ='


!    Optdriver=5
     nm_levels(50)   ='nonlinear     '
     ncounter(50)=0
     nm_counter(50,1)='            '

     nm_levels(51)   ='loop3dte      '
     ncounter(51)=2
     nm_counter(51,1)='pert1case  ='
     nm_counter(51,2)='pert3case  ='

     nm_levels(52)   ='mv_/resp3dte  '
     ncounter(52)=2
     nm_counter(52,2)='ikpt       ='

!    Optdriver=9
     nm_levels(90)   ='bethe_salpethe'
     ncounter(90)=0
     nm_counter(90,1)='            '



   else if(first==4)then

!    The true input variable "output_rate" is only available at the fourth
!    call to "status".
     if(statnu+1==4)then
       first=0
       if(istat<=0)then
         write(message, '(a,a,a,a,i7,a,a,a,a,a)' ) ch10,&
&         ' status : ERROR -',ch10,&
&         '  the value of the input variable istatr is',istat,' ,',ch10,&
&         '  while it must be a positive, non-zero number.',ch10,&
&         '  Action : change istatr in your input file.'
         call wrtout(std_out,message,'COLL')
         call leave_new('PERS')
       end if
       output_rate=istat
     end if

   end if
 end if

!The true input variable "shift_rate" is only available at the fifth call
 if(statnu+1==5)then
   if(istat<0 .or. istat>=output_rate)then
     write(message, '(4a,i7,3a,i7,2a)' ) ch10,&
&     ' status : ERROR -',ch10,&
&     '  the value of the input variable istatshft is',istat,' ,',ch10,&
&     '  while it must be a positive number smaller than istatr=',output_rate,ch10,&
&     '  Action : change istatshft in your input file.'
     call wrtout(std_out,message,'COLL')
     call leave_new('PERS')
   end if
   shift_rate=istat
 end if

!Check the value of level
 if( minval(abs(list_level(:)-level)) /= 0)then
   write(message, '(4a,i5,3a)' ) ch10,&
&   ' status : BUG -',ch10,&
&   '  The value of level in the calling routine is',level,' ,',ch10,&
&   '  which is not an allowed value.'
   call wrtout(std_out,message,'COLL')
   call leave_new('PERS')
 end if

!Assign the info about the actual routine
 write(unit=nm_routine(level),fmt='(a14)') routine
 if(trim(nm_routine(level))=='exit')then
!  The value 2 will be changed to 0 at the end of the routine.
   active(level)=2
 else if(trim(nm_routine(level))=='')then
   active(level)=0
 else
   active(level)=1
 end if

!Assign the info about the actual counter
 if(counter>=0)then
   if(ncounter(level)==1)then
     actual_counter(level,1)=counter
   else if(ncounter(level)==2)then
     actual_counter(level,2)=counter/100
!    The counter number 1 does not allow more than 99 passes
     actual_counter(level,1)=counter-(counter/100)*100
   end if
 end if

!============================================================

!After treatment of present information, output of the status
 statnu=statnu+1

!DEBUG
!write(6,*)' status : statnu=',statnu
!ENDDEBUG

 if( mod(statnu,output_rate)==shift_rate .or. statnu==10 )then

   open (tmp_unit,file=filstat,form='formatted',status='unknown',iostat=ios)
   if (ios/=0) then
     message = " Opening file: "//TRIM(filstat)
     MSG_PERS_ERROR(message)
   end if

   rewind tmp_unit
   write(tmp_unit,*)

   headwr='(a,i4,a,i6 )'
   if(statnu>=100000)   headwr='(a,i4,a,i9 )'
   if(output_rate>=1000)headwr='(a,i6,a,i6 )'
   if(statnu>=100000 .and. output_rate>=1000)   headwr='(a,i6,a,i9 )'
   if(statnu>=100000000)headwr='(a,i6,a,i12)'
   write(tmp_unit,headwr)&
&   ' Status file, with repetition rate',output_rate,&
&   ', status number',statnu
   write(tmp_unit,*)

!  Treat every level one after the other
   do ilevel=1,mlevel
!    This level must be activated in order to have a corresponding output
     if(active(ilevel)==1 .or. active(ilevel)==2)then

       write(tmp_unit,'(4a)')&
&       '  Level ',nm_levels(ilevel),' : ',nm_routine(ilevel)

!      Is there a counter for this level ?
       if(ncounter(ilevel)>=1)then

         if(actual_counter(ilevel,1)>0)then
           write(tmp_unit,'(a,a,i5)')'  ',nm_counter(ilevel,1),&
&           actual_counter(ilevel,1)
         end if
         if(ncounter(ilevel)==2)then
           if(actual_counter(ilevel,2)>0)then
             write(tmp_unit,'(a,a,i5)')'  ',nm_counter(ilevel,2),&
&             actual_counter(ilevel,2)
           end if
         end if

       end if

!      End of the check on activation of the level
     end if

!    End of the loop on the levels
   end do

   close (tmp_unit)

!  End of the repetition rate check
 end if

 if (active(level)==2)then
   active(level)=0
   nm_routine(level)='              '
 end if

 call timab(73,2,tsec)

end subroutine status
!!***
