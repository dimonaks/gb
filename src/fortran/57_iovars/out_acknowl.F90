!{\src2tex{textfont=tt}}
!!****f* ABINIT/out_acknowl
!! NAME
!! out_acknowl
!!
!! FUNCTION
!! Echo acknowledgments for the ABINIT code.
!!
!! COPYRIGHT
!! Copyright (C) 2009-2010 ABINIT group (XG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  iout=unit number for echoed output
!!  dtsets(0:ndtset_alloc)=<type datafiles_type>contains all input variables
!!  ndtset_alloc=number of datasets, corrected for allocation of at least
!!   one data set. Use for most dimensioned arrays.
!!  npsp=number of pseudopotentials
!!  pspheads(npsp)=<type pspheader_type>=all the important information from the
!!   pseudopotential file headers, as well as the psp file names
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!
!! PARENTS
!!      abinit
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine out_acknowl(dtsets,iout,ndtset_alloc,npsp,pspheads) 

 use defs_basis
 use defs_datatypes
 use defs_abitypes

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: iout,npsp,ndtset_alloc
 type(pspheader_type),intent(in) :: pspheads(npsp)
!arrays
 type(dataset_type),intent(in) :: dtsets(0:ndtset_alloc)

!Local variables-------------------------------
 integer :: idtset,iprior,iref,ncited,nrefs,ipsp,print_optional
 integer, allocatable :: cite(:),priority(:)
 character(len=600), allocatable :: ref(:)
 character(len=500), allocatable :: comment(:)

! *************************************************************************

!DEBUG
!write(6,*)' out_acknowl : enter '
!ENDDEBUG

!Allocate and initialize, for each possible reference, the flag for citation,
!the priority of the citation, the reference, and the comment. 
 nrefs=30
 allocate(cite(nrefs),ref(nrefs),comment(nrefs))

 allocate(priority(nrefs))
!The highest, the best, except that one from -1 and -2 should be cited. 
!0 means, cite if there are less than five papers total, otherwise forget, and any case, mention that it is optional.
!1-19 means specific papers, that must be cited. However, they might not appear in the top list of papers.
!20 means papers that should appear in the top list (usually, the most specific papers).

 cite(:)=0

 ref(1)=' ABINIT : First-principles approach of materials and nanosystem properties.'//ch10//&
& ' X. Gonze, B. Amadon, P.-M. Anglade, J.-M. Beuken, F. Bottin, P. Boulanger, F. Bruneval,'//ch10//&
& ' D. Caliste, R. Caracas, M. Cote, T. Deutsch, L. Genovese, Ph. Ghosez, M. Giantomassi'//ch10//&
& ' S. Goedecker, D.R. Hamann, P. Hermet, F. Jollet, G. Jomard, S. Leroux, M. Mancini, S. Mazevet,'//ch10//&
& ' M.J.T. Oliveira, G. Onida, Y. Pouillon, T. Rangel, G.-M. Rignanese, D. Sangalli, R. Shaltaf,'//ch10//&
& ' M. Torrent, M.J. Verstraete, G. Zerah, J.W. Zwanziger'//ch10//&
& ' Computer Phys. Comm. 180, 2582-2615 (2009).'
 comment(1)=' Comment : the third generic paper describing the ABINIT project.'//ch10//&
& ' Note that a version of this paper, that is not formatted for Computer Phys. Comm. '//ch10//&
& ' is available at http://www.abinit.org/about/ABINIT_CPC_v10.pdf .'//ch10//&
& ' The licence allows the authors to put it on the Web.'
 priority(1)=2

 ref(2)=' A brief introduction to the ABINIT software package.'//ch10//&
& ' X. Gonze, G.-M. Rignanese, M. Verstraete, J.-M. Beuken, Y. Pouillon, R. Caracas, F. Jollet,'//ch10//&
& ' M. Torrent, G. Zerah, M. Mikami, Ph. Ghosez, M. Veithen, J.-Y. Raty, V. Olevano, F. Bruneval,'//ch10//&
& ' L. Reining, R. Godby, G. Onida, D.R. Hamann, and D.C. Allan.'//ch10//&
& ' Z. Kristallogr. 220, 558-562 (2005).'
 comment(2)=' Comment : the second generic paper describing the ABINIT project. Note that this paper'//ch10//&
& ' should be cited especially if you are using the GW part of ABINIT, as several authors'//ch10//& 
& ' of this part are not in the list of authors of the first or third paper.'//ch10//&
& ' The .pdf of the latter paper is available at http://www.abinit.org/about/zfk_0505-06_558-562.pdf.'//ch10//&
& ' Note that it should not redistributed (Copyright by Oldenburg Wissenshaftverlag,'//ch10//&
& ' the licence allows the authors to put it on the Web).'
 priority(2)=1

 ref(3)=' First-principles computation of material properties : the ABINIT software project. '//ch10//&
& ' X. Gonze, J.-M. Beuken, R. Caracas, F. Detraux, M. Fuchs, G.-M. Rignanese, L. Sindic,'//ch10//&
& ' M. Verstraete, G. Zerah, F. Jollet, M. Torrent, A. Roy, M. Mikami, Ph. Ghosez, J.-Y. Raty, D.C. Allan.'//ch10//&
& ' Computational Materials Science 25, 478-492 (2002). http://dx.doi.org/10.1016/S0927-0256(02)00325-7'
 comment(3)=' Comment : the original paper describing the ABINIT project.'
 priority(3)=0

 ref(4)=' Fast radix 2, 3, 4 and 5 kernels for Fast Fourier Transformations'//ch10//&
& ' on computers with overlapping multiply-add instructions.'//ch10//&
& ' S. Goedecker, SIAM J. on Scientific Computing 18, 1605 (1997).'
 comment(4)=' '
 priority(4)=0 

 ref(5)=' Towards a potential-based conjugate gradient algorithm for order-N self-consistent'//ch10//&
& ' total energy calculations.'//ch10//&
& ' X. Gonze, Phys. Rev. B 54, 4383 (1996).'
 comment(5)=' Comment : The potential-based conjugate-gradient algorithm, used when iscf=5, is not published.'//ch10//&
& ' However, many elements of this algorithm have been explained in the paper above.'
 priority(5)=0

 ref(6)=' First-principles responses of solids to atomic displacements and homogeneous electric fields:,'//ch10//&
& ' implementation of a conjugate-gradient algorithm. X. Gonze, Phys. Rev. B55, 10337 (1997).'
 comment(6)=' Comment : Non-vanishing rfphon and/or rfelfd, in the norm-conserving case.'
 priority(6)=3

 ref(7)=' Dynamical matrices, Born effective charges, dielectric permittivity tensors, and ,'//ch10//&
& ' interatomic force constants from density-functional perturbation theory,'//ch10//&
& ' X. Gonze and C. Lee, Phys. Rev. B55, 10355 (1997).'
 comment(7)=' Comment : Non-vanishing rfphon and/or rfelfd, in the norm-conserving case.'
 priority(7)=3

 ref(8)=' Metric tensor formulation of strain in density-functional perturbation theory, '//ch10//&
& ' D. R. Hamann, X. Wu, K. M. Rabe, and D. Vanderbilt, Phys. Rev. B71, 035117 (2005).'
 comment(8)=' Comment : Non-vanishing rfstrs. Strong suggestion to cite this paper in your publications.'
 priority(8)=20

 ref(9)=' Ab initio pseudopotentials for electronic structure calculations of poly-atomic systems, '//ch10//&
& ' using density-functional theory.'//ch10//&
& ' M. Fuchs, M. Scheffler, Comput. Phys. Commun. 119, 67 (1999).'
 comment(9)=' Comment : Some pseudopotential generated using the FHI code were used.'
 priority(9)=3

 ref(10)=' Nonlinear optical susceptibilities, Raman efficiencies, and electrooptic tensors'//ch10//&
& ' from first principles density functional theory.'//ch10//&
& ' M. Veithen, X. Gonze, and Ph. Ghosez, Phys. Rev. B 71, 125107 (2005).'
 comment(10)=' Comment : to be cited for non-linear response calculations, with optdriver=5.'
 priority(10)=20

 ref(11)=' Effect of self-consistency on quasiparticles in solids'//ch10//&
& ' F. Bruneval, N. Vast, L. Reining, Phys. Rev. B 74, 045102 (2006).'
 comment(11)=' Comment : in case gwcalctyp >= 10.'
 priority(11)=18

 ref(12)=' Accurate GW self-energies in a plane-wave basis using only a few empty states:'//ch10//&
& ' towards large systems. F. Bruneval, X. Gonze, Phys. Rev. B 78, 085125 (2008).'
 comment(12)=' Comment : to be cited for non-vanishing gwcomp. Strong suggestion to cite this paper in your publications.'
 priority(12)=20

 ref(13)=' Large scale ab initio calculations based on three levels of parallelization'//ch10//&
& ' F. Bottin, S. Leroux, A. Knyazev, G. Zerah, Comput. Mat. Science 42, 329, (2008).'
 comment(13)=' Comment : in case paral_kgb is non-zero. Strong suggestion to cite this paper in your publications.'//ch10//&
& ' This paper is also available at http://www.arxiv.org/abs/0707.3405'
 priority(13)=20

 ref(14)=' Implementation of the Projector Augmented-Wave Method in the ABINIT code.'//ch10//&
& ' M. Torrent, F. Jollet, F. Bottin, G. Zerah, and X. Gonze Comput. Mat. Science 42, 337, (2008).'
 comment(14)=' Comment : PAW calculations. Strong suggestiong to cite this paper.'
 priority(14)=15

 ref(15)=' Gamma and beta cerium: LDA+U calculations of ground-state parameters.'//ch10//&
& ' B. Amadon, F. Jollet and M. Torrent, Phys. Rev. B 77, 155104 (2008).'
 comment(15)=' Comment : LDA+U calculations, usepawu/=0. Strong suggestion to cite this paper.'
 priority(15)=20

 ref(16)=' Preconditioning of self-consistent-field cycles in density functional theory : the extrapolar method'//ch10//& 
& ' P.-M. Anglade, X. Gonze, Phys. Rev. B 78, 045126 (2008).'
 comment(16)=' Comment : to be cited in case the extrapolar conditioner is used, i.e. non-vanishing iprcel.'
 priority(16)=10

 ref(17)=' Sharing electronic structure and crystallographic data with ETSF_IO'//ch10//&
& ' D. Caliste, Y. Pouillon, M.J. Verstraete, V. Olevano, X. Gonze,'//ch10//&
& ' Comput. Physics Communications 179, 748 (2008).'
 comment(17)=' Comment : to be cited in case the ETSF_IO file format is used, i.e. accesswff=3.'
 priority(17)=20

 ref(18)=' Daubechies wavelets as a basis set for density functional pseudopotential calculations.'//ch10//&
& ' L. Genovese, A. Neelov, S. Goedecker, T. Deutsch, S.A. Ghasemi, A. Willand,'// &
& ' D. Caliste, O. Zilberberg, M. Rayson, A. Bergman et R. Schneider,'//ch10//&
& ' J. Chem. Phys. 129, 014109 (2008).'
 comment(18)=' Comment : to be cited in case BigDFT project is used, i.e. usewvl=1.'
 priority(18)=5

 ref(19)=' Calculations of the transport properties within the PAW formalism.'//ch10//&
& ' S. Mazevet, M. Torrent, V. Recoules, F. Jollet,'// &
& ' High Energy Density Physics, 6, 84-88 (2010).'
 comment(19)=' Comment : to be cited in case output for transport properties calculation within PAW is used,'//ch10//&
& '           i.e. prtnabla>0 and usepaw=1.'
 priority(19)=20

 ref(20)=' Plane-wave based electronic strcture calculations for correlated materials.'//ch10//&
& ' using dynamical mean-field theory and projected local orbitals,'//ch10// &
& ' B. Amadon, F. Lechermann, A. Georges, F. Jollet, T.O. Wehling, A.I. Lichenstein,'//ch10// &
& ' Phys. Rev. B 77, 205112 (2008).'
 comment(20)=' Comment : to be cited in case the computation of overlap operator'// &
& ' for Wannier90 interface within PAW is used,'//ch10//&
& '           i.e. prtwant=2 and usepaw=1.'
 priority(20)=15

 ref(21)=' First-principles calculation of electric field gradients in metals, semiconductors, and insulators.'//ch10//&
& ' J.W. Zwanziger, M. Torrent,'// &
& ' Applied Magnetic Resonance 33, 447-456 (2008).'
 comment(21)=' Comment : to be cited in case the computation of electric field gradient is used, i.e. prtefg=1 and usepaw=1.'
 priority(21)=20

!---------------------------------------------------------------------------------------------
!Determine the papers to be cited

!Generic papers, not subject to conditions for citations
 cite(1:4)=1

!Go through the datasets
 do idtset=1,ndtset_alloc

!  If iscf=5 or iscf=15 used, cite Gonze96
   if(dtsets(idtset)%iscf==5)cite(5)=1
   if(dtsets(idtset)%iscf==15)cite(5)=1

!  If rfphon/=0 or rfelfd/=0, cite Gonze97a
   if(dtsets(idtset)%rfphon/=0)cite(6)=1
   if(dtsets(idtset)%rfelfd/=0)cite(6)=1

!  If rfphon/=0 or rfelfd/=0, cite Gonze97b
   if(dtsets(idtset)%rfphon/=0)cite(7)=1
   if(dtsets(idtset)%rfelfd/=0)cite(7)=1

!  If rfstrs/=0, cite Hamann05
   if(dtsets(idtset)%rfstrs/=0)cite(8)=1

!  If optdriver==5, cite Veithen2005
   if(dtsets(idtset)%optdriver==5)cite(10)=1

!  If gwcalctyp>=10, cite Bruneval2006
   if(dtsets(idtset)%gwcalctyp>=10)cite(11)=1

!  If gwcomp/=0, cite Bruneval2008
   if(dtsets(idtset)%gwcomp/=0)cite(12)=1

!  If paral_kgb/=0, cite Bottin2008
   if(dtsets(idtset)%paral_kgb/=0)cite(13)=1

!  If usepaw/=0, cite Torrent2008
   if(dtsets(idtset)%usepaw/=0)cite(14)=1

!  If usepawu/=0, cite Amadon2008
   if(dtsets(idtset)%usepawu/=0)cite(15)=1

!  If iprcel/=0, cite Anglade2008
   if(dtsets(idtset)%iprcel/=0)cite(16)=1

!  If accesswff==IO_MODE_ETSF, cite Caliste2008
   if(dtsets(idtset)%accesswff==IO_MODE_ETSF)cite(17)=1

!  If usewvl/=0, cite Genovese2008
   if(dtsets(idtset)%usewvl/=0)cite(18)=1

!  If prtnabla/=0, cite Mazevet2010
   if(dtsets(idtset)%usepaw==1.and.dtsets(idtset)%prtnabla>0)cite(19)=1

!  If prtnabla/=0, cite Amadon2008
   if(dtsets(idtset)%usepaw==1.and.dtsets(idtset)%prtwant==2)cite(20)=1

!  If prtefg/=0, cite Zwanziger2008
   if(dtsets(idtset)%usepaw==1.and.dtsets(idtset)%prtefg>0)cite(21)=1

 end do

!Go through the pseudopotentials
 do ipsp=1,npsp

!  If FHI pseudopotential, cite Fuchs 1999
   if(pspheads(ipsp)%pspcod==6)cite(9)=1
 end do

!-------------------------------------------------------------------------------------------
!Assemble the acknowledgment notice

 write(iout, '(30a)' )ch10,&
& '================================================================================',ch10,ch10,&
& ' Suggested references for the acknowledgment of ABINIT usage.',ch10,ch10,&
& ' The users of ABINIT have little formal obligations with respect to the ABINIT group',ch10,&
& ' (those specified in the GNU General Public License, http://www.gnu.org/copyleft/gpl.txt).',ch10,&
& ' However, it is common practice in the scientific literature,',ch10,&
& ' to acknowledge the efforts of people that have made the research possible.',ch10,&
& ' In this spirit, please find below suggested citations of work written by ABINIT developers,',ch10,&
& ' corresponding to implementations inside of ABINIT that you have used in the present run.',ch10,&
& ' Note also that it will be of great value to readers of publications presenting these results,',ch10,&
& ' to read papers enabling them to understand the theoretical formalism and details',ch10,&
& ' of the ABINIT implementation.',ch10,&
& ' For information on why they are suggested, see also http://www.abinit.org/about/?text=acknowledgments.'

 ncited=0 
 print_optional=1

 do iprior=20,0,-1
   do iref=1,nrefs
     if(cite(iref)==1)then
       if(priority(iref)==iprior)then
         if(priority(iref)>0 .or. &
&         (priority(iref)==1 .and. ncited<5) .or. (priority(iref)==0 .and. ncited<5)) then
           ncited=ncited+1
           cite(iref)=0
           if(priority(iref)==0 .and. print_optional==1)then
             print_optional=0
             write(iout,'(a)')ch10,' And optionally :'
           end if
           if(len_trim(comment(iref))/=0)then
             if(ncited<10)write(iout, '(2a,i1,4a)')ch10,' [',ncited,']',trim(ref(iref)),ch10,trim(comment(iref))
             if(ncited>=10)write(iout, '(2a,i2,4a)')ch10,' [',ncited,']',trim(ref(iref)),ch10,trim(comment(iref))
           else
             if(ncited<10)write(iout, '(2a,i1,4a)')ch10,' [',ncited,']',trim(ref(iref))
             if(ncited>=10)write(iout, '(2a,i2,4a)')ch10,' [',ncited,']',trim(ref(iref))
           end if
         end if
       end if
       if(priority(iref)==0 .and. ncited>=5)cite(iref)=0
     end if
   end do
 end do

!-------------------------------------------------------------------------------------------
!Cleaning

 deallocate(cite,ref,comment,priority)

!DEBUG
!write(6,*)' out_acknowl : end of subroutine '
!ENDDEBUG

end subroutine out_acknowl
!!***
