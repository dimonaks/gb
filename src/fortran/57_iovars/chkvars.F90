!{\src2tex{textfont=tt}}
!!****f* ABINIT/chkvars
!! NAME
!! chkvars
!!
!! FUNCTION
!! Examine the input string, to check whether all names are allowed
!!
!! COPYRIGHT
!! Copyright (C) 2007-2010 ABINIT group (XG).
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  string*(*)=string of character
!!   the string (with upper case) from the input file, to which the CML data are appended
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!
!! PARENTS
!!      abinit
!!
!! CHILDREN
!!      inupper,leave_new,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine chkvars (string)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
 use interfaces_14_hidewrite
 use interfaces_16_hideleave
 use interfaces_32_util
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 character(len=*),intent(in) :: string

!Local variables-------------------------------
 character :: blank=' '
!scalars
 integer :: index_blank,index_current,index_endword,index_endwordnow,index_list_vars
 character(len=100) :: list_logicals,list_strings
 character(len=10000) :: list_vars
 character(len=500) :: message

!************************************************************************

!DEBUG
!write(6,*)' chkvars : enter '
!write(6,*)trim(string)
!stop
!ENDDEBUG

!Here, list all admitted variable names (max 10 per line, to fix the ideas)
!A
 list_vars=                 ' accesswff acell algalch alpha amu angdeg atvshift awtr'
!B
 list_vars=trim(list_vars)//' bandpp bdberry bdeigrf bdgw berryopt bfield bmass boxcenter boxcutmin brvltt bxctmindg'
!BEGIN VARIABLES FOR @Bethe-Salpeter
 list_vars=trim(list_vars)//' bs_algorithm bs_haydock_niter bs_haydock_tol bs_exchange_term bs_coulomb_term '
 list_vars=trim(list_vars)//' bs_calctype bs_coupling bs_haydock_tol bs_eh_basis_set bs_eh_cutoff bs_freq_mesh '
!END VARIABLES FOR @Bethe-Salpeter.
!C
 list_vars=trim(list_vars)//' charge chkexit chkprim chksymbreak cpus cpum cpuh'
!D
 list_vars=trim(list_vars)//' delayperm densty diecut diegap dielam dielng diemac'
 list_vars=trim(list_vars)//' diemix diemixmag dilatmx dmatpuopt dmatudiag dmatpawu dmft_iter '
 list_vars=trim(list_vars)//' dmft_mxsf dmft_nwli dmft_nwlo dmft_rslf dmft_solv dmftbandi dmftbandf'
 list_vars=trim(list_vars)//' dmftcheck dosdeltae dtion dynimage'
!E
 list_vars=trim(list_vars)//' ecut ecuteps ecutsigx ecutsm ecutwfn effmass efield'
 list_vars=trim(list_vars)//' enunit eshift esmear etsfgroups etsfmain exchmix exchn2n3d'
!F
 list_vars=trim(list_vars)//' pawfatbnd fband fftalg fftcache fftgw fft_opt_lob fixmom'
 list_vars=trim(list_vars)//' freqremax freqspmax freqsusin freqsuslo friction frzfermi fxcartfactor '
!G
 list_vars=trim(list_vars)//' genafm getbseig getcell getsuscep getddk getden getkss getocc getqps getscr'
 list_vars=trim(list_vars)//' getvel getwfk getwfq getxcart getxred get1den get1wf'
 list_vars=trim(list_vars)//' gwcalctyp gwcomp gwencomp gwgamma gwmem gwpara gwrpacorr'
 list_vars=trim(list_vars)//' gw_nqlwl gw_nstep gw_qlwl gw_sctype gw_sctype_name gw_sigxcore gw_toldfeig '
 list_vars=trim(list_vars)//' gw_EET gw_EET_nband'

!I
 list_vars=trim(list_vars)//' iatcon iatfix iatfixx iatfixy iatfixz iatsph'
 list_vars=trim(list_vars)//' iboxcut icoulomb icutcoul idyson ieig2rf ikhxc' ! XG091017 : iextrapwf has been removed from this list,
!because no doc or test were provided ... 
 list_vars=trim(list_vars)//' imgmov inclvkb intexact intxc ionmov iprcch'
 list_vars=trim(list_vars)//' iprcel iprctfvw iprcfc irdbseig irdddk irdden irdqps'
 list_vars=trim(list_vars)//' irdkss irdscr irdsuscep irdwfk irdwfq ird1wf iscf'
 list_vars=trim(list_vars)//' isecur istatr istatshft istwfk ixc ixcpositron'
!J
 list_vars=trim(list_vars)//' jdtset jellslab jpawu'
!K
 list_vars=trim(list_vars)//' kberry kpt kptbounds kptgw'
 list_vars=trim(list_vars)//' kptnrm kptopt kptrlatt kptrlen kssform'
!L
 list_vars=trim(list_vars)//' ldgapp lexexch localrdwf lpawu'
!M
 list_vars=trim(list_vars)//' macro_uj maxnsym mdftemp mditemp mdwall mffmem mixalch'
 list_vars=trim(list_vars)//' mkmem mkqmem mk1mem mqgrid mqgriddg'
!N
 list_vars=trim(list_vars)//' natcon natfix natfixx natfixy natfixz'
 list_vars=trim(list_vars)//' natom natrd natsph natvshift nband nbandkss'
 list_vars=trim(list_vars)//' nbandsus nbdblock nbdbuf nberry nconeq'
 list_vars=trim(list_vars)//' nctime ndivk ndivsm ndtset ndyson'
 list_vars=trim(list_vars)//' nfreqim nfreqre nfreqsp nfreqsus ngfft ngfftdg'
 list_vars=trim(list_vars)//' ngkpt ngroup_rf nkpt nkptgw nimage nline nloalg nomegasi nnos nnsclo'
 list_vars=trim(list_vars)//' nobj nomegasf nomegasrd normpawu noseinert npband'
 list_vars=trim(list_vars)//' npfft npimage npkpt npsp npulayit npweps npwkss npwsigx'
 list_vars=trim(list_vars)//' npwwfn nqpt nqptdm nscforder nsheps nshiftk'
 list_vars=trim(list_vars)//' nshsigx nshwfn nspden nspinor nsppol nstep nsym'
 list_vars=trim(list_vars)//' ntime ntimimage ntypalch ntypat nwfshist nshift'
!O
 list_vars=trim(list_vars)//' objaat objbat objaax objbax objan objbn objarf'
 list_vars=trim(list_vars)//' objbrf objaro objbro objatr objbtr occ'
 list_vars=trim(list_vars)//' occopt omegasimax omegasrdmax optcell optdriver optforces'
 list_vars=trim(list_vars)//' optfreqsus optnlxccc optstress ortalg'
!P
 list_vars=trim(list_vars)//' paral_kgb paral_rf pawcpxocc pawecutdg pawlcutd pawlmix pawmixdg'
 list_vars=trim(list_vars)//' pawnhatxc pawnphi pawntheta pawnzlm pawovlp pawoptmix'
 list_vars=trim(list_vars)//' pawprtden pawprtdos pawprtvol  pawprtwf pawspnorb'
 list_vars=trim(list_vars)//' pawstgylm pawujat pawujrad pawujv pawusecp pawxcdev pitransform'
 list_vars=trim(list_vars)//' positron posnstep posocc postoldfe postoldff ppmfrq ppmodel prepanl prepgkk'
 list_vars=trim(list_vars)//' prtbbb prtbltztrp prtcif prtcml prtden papiopt'
 list_vars=trim(list_vars)//' prtdensph prtdipole prtdos prtdosm prtefg prteig prtelf'
 list_vars=trim(list_vars)//' prtfc prtfsurf prtgden prtgeo prtgkk prtkden prtkpt prtlden' 
 list_vars=trim(list_vars)//' prtnabla prtnest prtposcar prtpot prtspcur prtstm prtvha prtvhxc'
 list_vars=trim(list_vars)//' prtvol prtvxc prtwant prtwf prtxangst prtxcart prtxml prtxred prt1dm pspso ptcharge'
!Q
 list_vars=trim(list_vars)//' qmass qprtrb qpt qptdm qptnrm quadmom'
!R
 list_vars=trim(list_vars)//' ratsph rcut'  !  XG091017 : rdmnb has been removed from this list, because no doc or test were provided ... 
 list_vars=trim(list_vars)//' recefermi recgratio recnpath recnrec recptrott recrcut rectesteg rectolden'
 list_vars=trim(list_vars)//' restartxf rfasr rfatpol rfddk rfdir rfelfd rfmgfd rfmeth rfphon'
 list_vars=trim(list_vars)//' rfstrs rfuser rf1atpol rf1dir rf1elfd'
 list_vars=trim(list_vars)//' rf1phon rf2atpol rf2dir rf2elfd rf2phon'
 list_vars=trim(list_vars)//' rf3atpol rf3dir rf3elfd rf3phon rhoqpmix rprim '
!S
 list_vars=trim(list_vars)//' scalecart sciss scphon_temp shiftk signperm '
 list_vars=trim(list_vars)//' slabwsrad slabzbeg slabzend smdelta soenergy so_psp '
 list_vars=trim(list_vars)//' spbroad spgaxor spgorig spgroup spgroupma spmeth spinat '
 list_vars=trim(list_vars)//' spnorbscl stmbias strfact strprecon strtarget scphon_supercell supercell'
 list_vars=trim(list_vars)//' suskxcrs symafm symchi symrel symmorphi symsigma'
!T
 list_vars=trim(list_vars)//' td_maxene td_mexcit tfkinfunc tl_nprccg tl_radius timopt '
 list_vars=trim(list_vars)//' tnons toldfe toldff tolimg tolmxf tolrff tolsym tolvrs'
 list_vars=trim(list_vars)//' tolwfr tphysel tsmear typat'
!U
 list_vars=trim(list_vars)//' udtset upawu usedmatpu usedmft useexexch usekden usepawu usewvl'
 list_vars=trim(list_vars)//' useria userib useric userid userie'
 list_vars=trim(list_vars)//' userra userrb userrc userrd userre usexcnhat useylm'
!V
 list_vars=trim(list_vars)//' vaclst vacnum vacuum vacwidth vcutgeo vdw_xc vdw_nwan vdw_supercell vel vis vprtrb'
!W
 list_vars=trim(list_vars)//' wfoptalg wtatcon wtk'
 list_vars=trim(list_vars)//' wvl_cpmult wvl_crmult wvl_fpmult wvl_frmult wvl_hgrid wvl_nprccg'
 list_vars=trim(list_vars)//' w90iniprj w90prtunk'
!X
 list_vars=trim(list_vars)//' xangst xcart xred xyzfile '
!Y
!Z
 list_vars=trim(list_vars)//' zcut zeemanfield znucl'

!Extra token, also admitted :
 list_vars=trim(list_vars)//' au Angstr Angstrom Angstroms Bohr Bohrs eV Ha Hartree Hartrees K Ry Rydberg Rydbergs T Tesla sqrt end'
!SIESTA strings
 list_vars=trim(list_vars)//' LatticeConstant '

!Logical input variables
 list_logicals=' SpinPolarized '

!String input variables
 list_strings=' cmlfile XCname '

!Transform to upper case
 call inupper(list_vars)
 call inupper(list_logicals)
 call inupper(list_strings)

!DEBUG
!write(6,*)' chkvars : len_trim(string)=',len_trim(string)
!write(6,*)' chkvars : trim(string)=',trim(string)
!ENDDEBUG

 index_current=1
 do ! Infinite do-loop, to identify the presence of each potential variable names

   if(len_trim(string)<=index_current)exit
   index_blank=index(string(index_current:),blank)+index_current-1

   if(index('ABCDEFGHIJKLMNOPQRSTUVWXYZ',string(index_current:index_current))/=0)then

!    Skip characters like : + or the digits at the end of the word
!    Start from the blank that follows the end of the word
     do index_endword=index_blank-1,index_current,-1
       if(index('ABCDEFGHIJKLMNOPQRSTUVWXYZ',string(index_endword:index_endword))/=0)exit
     end do

!    Find the index of the potential variable name in the list of variables
     index_list_vars=index(list_vars,blank//string(index_current:index_endword)//blank)

!    Treat the complications due to the possibility of images
     if(index_list_vars==0)then

!      Treat possible LASTIMG appendix
       if(index_endword-6>=1)then
         if(string(index_endword-6:index_endword)=='LASTIMG')index_endword=index_endword-7
       end if

!      Treat possible IMG appendix
       if(index_endword-2>=1)then
         if(string(index_endword-2:index_endword)=='IMG')index_endword=index_endword-3
       end if

       index_endwordnow=index_endword

!      Again skip characters like : + or the digits before IMG
!      Start from the blank that follows the end of the word
       do index_endword=index_endwordnow,index_current,-1
         if(index('ABCDEFGHIJKLMNOPQRSTUVWXYZ',string(index_endword:index_endword))/=0)exit
       end do

!      Find the index of the potential variable name in the list of variables
       index_list_vars=index(list_vars,blank//string(index_current:index_endword)//blank)

     end if

!    DEBUG
!    write(6,*)' chkvars :'
!    write(6,*)' index_current, index_blank, string',index_current,index_blank,string(index_current:index_endword)
!    write(6,*)' index_list_vars=',index_list_vars
!    if(index_list_vars/=0)then
!    write(6,*)' list_vars(index_list_vars:index_list_vars+12)=',&
!    &    list_vars(index_list_vars:min(index_list_vars+12,len_trim(list_vars)))
!    else 
!    write(6,*)' not in list_vars, might be in list_logicals or list_strings '
!    endif
!    write(6,*)
!    ENDDEBUG

     if(index_list_vars==0)then

!      DEBUG
!      write(6,*)' list_logicals,index=:',list_logicals,':',index(list_strings,string(index_current:index_endword)//blank)
!      write(6,*)' list_strings,index=:',list_strings,':',index(list_strings,string(index_current:index_endword)//blank)
!      ENDDEBUG

!      Treat possible logical input variables
       if(index(list_logicals,blank//string(index_current:index_endword)//blank)/=0)then
         index_blank=index(string(index_current:),blank)+index_current-1
         if(index(' F T ',string(index_blank:index_blank+2))==0)then
           write(message, '(11a)' ) ch10,&
&           ' chkvars : ERROR - ',ch10,&
&           '  Found the token ',string(index_current:index_endword),' in the input file.',ch10,&
&           '  This variable should be given a logical value (T or F), but the following string was found :',&
&           string(index_blank:index_blank+2),ch10,&
&           '  Action : check your input file. You likely misused the input variable.'
           call wrtout(std_out,message,'COLL')
           call leave_new('COLL')
         else
           index_blank=index_blank+2
         end if
!        Treat possible string input variables
       else if(index(list_strings,blank//string(index_current:index_endword)//blank)/=0)then
!        Every following string is accepted
         index_current=index(string(index_current:),blank)+index_current
         index_blank=index(string(index_current:),blank)+index_current-1

!        If still not admitted, then there is a problem
       else
         write(message, '(10a)' ) ch10,&
&         ' chkvars : ERROR - ',ch10,&
&         '  Found the token ',string(index_current:index_endword),' in the input file.',ch10,&
&         '  This name is not one of the registered input variable names (see the Web list of input variables).',ch10,&
&         '  Action : check your input file. You likely mistyped the input variable.'
         call wrtout(std_out,message,'COLL')
         call leave_new('COLL')
       end if
     end if
   end if

   index_current=index_blank+1

 end do

!DEBUG
!write(6,*)' chkvars : exit '
!stop
!ENDDEBUG

end subroutine chkvars
!!***
