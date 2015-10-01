! This file has been automatically generated, do not modify.

  subroutine ab6_invars_get_integer(dtsetsId, value, id, idtset, errno)
    integer, intent(in) :: dtsetsId
    integer, intent(in) :: id, idtset
    integer, intent(out) :: value
    integer, intent(out) :: errno

    type(dtsets_list), pointer :: token

    call get_token(token, dtsetsId)
    if (.not. associated(token)) then
      errno = AB6_ERROR_OBJ
      return
    end if
    if (idtset < 0 .or. idtset > size(token%dtsets)) then
      errno = AB6_ERROR_INVARS_ID
      return
    end if
    
    errno = AB6_NO_ERROR
    select case (id)
    case (ab6_invars_symchi)
      value = token%dtsets(idtset)%symchi
    case (ab6_invars_get1wf)
      value = token%dtsets(idtset)%get1wf
    case (ab6_invars_rfstrs)
      value = token%dtsets(idtset)%rfstrs
    case (ab6_invars_exchn2n3d)
      value = token%dtsets(idtset)%exchn2n3d
    case (ab6_invars_nstep)
      value = token%dtsets(idtset)%nstep
    case (ab6_invars_wvl_nprccg)
      value = token%dtsets(idtset)%wvl_nprccg
    case (ab6_invars_gwgamma)
      value = token%dtsets(idtset)%gwgamma
    case (ab6_invars_nsheps)
      value = token%dtsets(idtset)%nsheps
    case (ab6_invars_prt1dm)
      value = token%dtsets(idtset)%prt1dm
    case (ab6_invars_prepanl)
      value = token%dtsets(idtset)%prepanl
    case (ab6_invars_nspinor)
      value = token%dtsets(idtset)%nspinor
    case (ab6_invars_prtxred)
      value = token%dtsets(idtset)%prtxred
    case (ab6_invars_ndtset)
      value = token%dtsets(idtset)%ndtset
    case (ab6_invars_usepawu)
      value = token%dtsets(idtset)%usepawu
    case (ab6_invars_mpw)
      value = token%dtsets(idtset)%mpw
    case (ab6_invars_occopt)
      value = token%dtsets(idtset)%occopt
    case (ab6_invars_prtxcart)
      value = token%dtsets(idtset)%prtxcart
    case (ab6_invars_npspalch)
      value = token%dtsets(idtset)%npspalch
    case (ab6_invars_td_mexcit)
      value = token%dtsets(idtset)%td_mexcit
    case (ab6_invars_frzfermi)
      value = token%dtsets(idtset)%frzfermi
    case (ab6_invars_rfphon)
      value = token%dtsets(idtset)%rfphon
    case (ab6_invars_prtden)
      value = token%dtsets(idtset)%prtden
    case (ab6_invars_gwpara)
      value = token%dtsets(idtset)%gwpara
    case (ab6_invars_npfft)
      value = token%dtsets(idtset)%npfft
    case (ab6_invars_recptrott)
      value = token%dtsets(idtset)%recptrott
    case (ab6_invars_jdtset)
      value = token%dtsets(idtset)%jdtset
    case (ab6_invars_bs_haydock_niter)
      value = token%dtsets(idtset)%bs_haydock_niter
    case (ab6_invars_pawprtden)
      value = token%dtsets(idtset)%pawprtden
    case (ab6_invars_bs_algorithm)
      value = token%dtsets(idtset)%bs_algorithm
    case (ab6_invars_gwmem)
      value = token%dtsets(idtset)%gwmem
    case (ab6_invars_npkpt)
      value = token%dtsets(idtset)%npkpt
    case (ab6_invars_iextrapwf)
      value = token%dtsets(idtset)%iextrapwf
    case (ab6_invars_usekden)
      value = token%dtsets(idtset)%usekden
    case (ab6_invars_nberry)
      value = token%dtsets(idtset)%nberry
    case (ab6_invars_chksymbreak)
      value = token%dtsets(idtset)%chksymbreak
    case (ab6_invars_usewvl)
      value = token%dtsets(idtset)%usewvl
    case (ab6_invars_useexexch)
      value = token%dtsets(idtset)%useexexch
    case (ab6_invars_getden)
      value = token%dtsets(idtset)%getden
    case (ab6_invars_tfkinfunc)
      value = token%dtsets(idtset)%tfkinfunc
    case (ab6_invars_irdbseig)
      value = token%dtsets(idtset)%irdbseig
    case (ab6_invars_ionmov)
      value = token%dtsets(idtset)%ionmov
    case (ab6_invars_npwwfn)
      value = token%dtsets(idtset)%npwwfn
    case (ab6_invars_mffmem)
      value = token%dtsets(idtset)%mffmem
    case (ab6_invars_bs_exchange_term)
      value = token%dtsets(idtset)%bs_exchange_term
    case (ab6_invars_mk1mem)
      value = token%dtsets(idtset)%mk1mem
    case (ab6_invars_nfreqsus)
      value = token%dtsets(idtset)%nfreqsus
    case (ab6_invars_prtkpt)
      value = token%dtsets(idtset)%prtkpt
    case (ab6_invars_prtcs)
      value = token%dtsets(idtset)%prtcs
    case (ab6_invars_imgmov)
      value = token%dtsets(idtset)%imgmov
    case (ab6_invars_ndyson)
      value = token%dtsets(idtset)%ndyson
    case (ab6_invars_userec)
      value = token%dtsets(idtset)%userec
    case (ab6_invars_spmeth)
      value = token%dtsets(idtset)%spmeth
    case (ab6_invars_mkmem)
      value = token%dtsets(idtset)%mkmem
    case (ab6_invars_natpawu)
      value = token%dtsets(idtset)%natpawu
    case (ab6_invars_prtfc)
      value = token%dtsets(idtset)%prtfc
    case (ab6_invars_usedmatpu)
      value = token%dtsets(idtset)%usedmatpu
    case (ab6_invars_paral_kgb)
      value = token%dtsets(idtset)%paral_kgb
    case (ab6_invars_irdden)
      value = token%dtsets(idtset)%irdden
    case (ab6_invars_getbseig)
      value = token%dtsets(idtset)%getbseig
    case (ab6_invars_getxred)
      value = token%dtsets(idtset)%getxred
    case (ab6_invars_w90prtunk)
      value = token%dtsets(idtset)%w90prtunk
    case (ab6_invars_intexact)
      value = token%dtsets(idtset)%intexact
    case (ab6_invars_rfmgfd)
      value = token%dtsets(idtset)%rfmgfd
    case (ab6_invars_useylm)
      value = token%dtsets(idtset)%useylm
    case (ab6_invars_macro_uj)
      value = token%dtsets(idtset)%macro_uj
    case (ab6_invars_dmft_rslf)
      value = token%dtsets(idtset)%dmft_rslf
    case (ab6_invars_npimage)
      value = token%dtsets(idtset)%npimage
    case (ab6_invars_getcell)
      value = token%dtsets(idtset)%getcell
    case (ab6_invars_maxnsym)
      value = token%dtsets(idtset)%maxnsym
    case (ab6_invars_kssform)
      value = token%dtsets(idtset)%kssform
    case (ab6_invars_wfoptalg)
      value = token%dtsets(idtset)%wfoptalg
    case (ab6_invars_rf1elfd)
      value = token%dtsets(idtset)%rf1elfd
    case (ab6_invars_irdwfk)
      value = token%dtsets(idtset)%irdwfk
    case (ab6_invars_ntimimage)
      value = token%dtsets(idtset)%ntimimage
    case (ab6_invars_pawlcutd)
      value = token%dtsets(idtset)%pawlcutd
    case (ab6_invars_prteig)
      value = token%dtsets(idtset)%prteig
    case (ab6_invars_enunit)
      value = token%dtsets(idtset)%enunit
    case (ab6_invars_prtnest)
      value = token%dtsets(idtset)%prtnest
    case (ab6_invars_irdwfq)
      value = token%dtsets(idtset)%irdwfq
    case (ab6_invars_pawstgylm)
      value = token%dtsets(idtset)%pawstgylm
    case (ab6_invars_iprcel)
      value = token%dtsets(idtset)%iprcel
    case (ab6_invars_ldgapp)
      value = token%dtsets(idtset)%ldgapp
    case (ab6_invars_restartxf)
      value = token%dtsets(idtset)%restartxf
    case (ab6_invars_irdddk)
      value = token%dtsets(idtset)%irdddk
    case (ab6_invars_spgroup)
      value = token%dtsets(idtset)%spgroup
    case (ab6_invars_npwsigx)
      value = token%dtsets(idtset)%npwsigx
    case (ab6_invars_nomegasi)
      value = token%dtsets(idtset)%nomegasi
    case (ab6_invars_nomegasf)
      value = token%dtsets(idtset)%nomegasf
    case (ab6_invars_nline)
      value = token%dtsets(idtset)%nline
    case (ab6_invars_gw_eet_nband)
      value = token%dtsets(idtset)%gw_eet_nband
    case (ab6_invars_rfmeth)
      value = token%dtsets(idtset)%rfmeth
    case (ab6_invars_useria)
      value = token%dtsets(idtset)%useria
    case (ab6_invars_nfreqsp)
      value = token%dtsets(idtset)%nfreqsp
    case (ab6_invars_useric)
      value = token%dtsets(idtset)%useric
    case (ab6_invars_userid)
      value = token%dtsets(idtset)%userid
    case (ab6_invars_userie)
      value = token%dtsets(idtset)%userie
    case (ab6_invars_ntypalch)
      value = token%dtsets(idtset)%ntypalch
    case (ab6_invars_localrdwf)
      value = token%dtsets(idtset)%localrdwf
    case (ab6_invars_prtdosm)
      value = token%dtsets(idtset)%prtdosm
    case (ab6_invars_nimage)
      value = token%dtsets(idtset)%nimage
    case (ab6_invars_irdkss)
      value = token%dtsets(idtset)%irdkss
    case (ab6_invars_rectesteg)
      value = token%dtsets(idtset)%rectesteg
    case (ab6_invars_ntypat)
      value = token%dtsets(idtset)%ntypat
    case (ab6_invars_icoulomb)
      value = token%dtsets(idtset)%icoulomb
    case (ab6_invars_ikhxc)
      value = token%dtsets(idtset)%ikhxc
    case (ab6_invars_prtxml)
      value = token%dtsets(idtset)%prtxml
    case (ab6_invars_ndynimage)
      value = token%dtsets(idtset)%ndynimage
    case (ab6_invars_nbdbuf)
      value = token%dtsets(idtset)%nbdbuf
    case (ab6_invars_prtgeo)
      value = token%dtsets(idtset)%prtgeo
    case (ab6_invars_rf1phon)
      value = token%dtsets(idtset)%rf1phon
    case (ab6_invars_iscf)
      value = token%dtsets(idtset)%iscf
    case (ab6_invars_bdeigrf)
      value = token%dtsets(idtset)%bdeigrf
    case (ab6_invars_pawujat)
      value = token%dtsets(idtset)%pawujat
    case (ab6_invars_ixc)
      value = token%dtsets(idtset)%ixc
    case (ab6_invars_delayperm)
      value = token%dtsets(idtset)%delayperm
    case (ab6_invars_nscforder)
      value = token%dtsets(idtset)%nscforder
    case (ab6_invars_prtpmp)
      value = token%dtsets(idtset)%prtpmp
    case (ab6_invars_rf3elfd)
      value = token%dtsets(idtset)%rf3elfd
    case (ab6_invars_suskxcrs)
      value = token%dtsets(idtset)%suskxcrs
    case (ab6_invars_iprcfc)
      value = token%dtsets(idtset)%iprcfc
    case (ab6_invars_usexcnhat)
      value = token%dtsets(idtset)%usexcnhat
    case (ab6_invars_dmft_solv)
      value = token%dtsets(idtset)%dmft_solv
    case (ab6_invars_prepgkk)
      value = token%dtsets(idtset)%prepgkk
    case (ab6_invars_prtvxc)
      value = token%dtsets(idtset)%prtvxc
    case (ab6_invars_nfftdg)
      value = token%dtsets(idtset)%nfftdg
    case (ab6_invars_paral_rf)
      value = token%dtsets(idtset)%paral_rf
    case (ab6_invars_nbdblock)
      value = token%dtsets(idtset)%nbdblock
    case (ab6_invars_irdscr)
      value = token%dtsets(idtset)%irdscr
    case (ab6_invars_usepaw)
      value = token%dtsets(idtset)%usepaw
    case (ab6_invars_nsym)
      value = token%dtsets(idtset)%nsym
    case (ab6_invars_rfuser)
      value = token%dtsets(idtset)%rfuser
    case (ab6_invars_prtfsurf)
      value = token%dtsets(idtset)%prtfsurf
    case (ab6_invars_npband)
      value = token%dtsets(idtset)%npband
    case (ab6_invars_nfft)
      value = token%dtsets(idtset)%nfft
    case (ab6_invars_rfasr)
      value = token%dtsets(idtset)%rfasr
    case (ab6_invars_accesswff)
      value = token%dtsets(idtset)%accesswff
    case (ab6_invars_fftgw)
      value = token%dtsets(idtset)%fftgw
    case (ab6_invars_prtkden)
      value = token%dtsets(idtset)%prtkden
    case (ab6_invars_rfelfd)
      value = token%dtsets(idtset)%rfelfd
    case (ab6_invars_prtvhxc)
      value = token%dtsets(idtset)%prtvhxc
    case (ab6_invars_mgfft)
      value = token%dtsets(idtset)%mgfft
    case (ab6_invars_dmatudiag)
      value = token%dtsets(idtset)%dmatudiag
    case (ab6_invars_nbandsus)
      value = token%dtsets(idtset)%nbandsus
    case (ab6_invars_recgratio)
      value = token%dtsets(idtset)%recgratio
    case (ab6_invars_iprctfvw)
      value = token%dtsets(idtset)%iprctfvw
    case (ab6_invars_nshwfn)
      value = token%dtsets(idtset)%nshwfn
    case (ab6_invars_pawspnorb)
      value = token%dtsets(idtset)%pawspnorb
    case (ab6_invars_npwkss)
      value = token%dtsets(idtset)%npwkss
    case (ab6_invars_gw_nqlwl)
      value = token%dtsets(idtset)%gw_nqlwl
    case (ab6_invars_prtcml)
      value = token%dtsets(idtset)%prtcml
    case (ab6_invars_gwcalctyp)
      value = token%dtsets(idtset)%gwcalctyp
    case (ab6_invars_irdsuscep)
      value = token%dtsets(idtset)%irdsuscep
    case (ab6_invars_chkprim)
      value = token%dtsets(idtset)%chkprim
    case (ab6_invars_prtefg)
      value = token%dtsets(idtset)%prtefg
    case (ab6_invars_bs_coupling)
      value = token%dtsets(idtset)%bs_coupling
    case (ab6_invars_getvel)
      value = token%dtsets(idtset)%getvel
    case (ab6_invars_dmftbandf)
      value = token%dtsets(idtset)%dmftbandf
    case (ab6_invars_pawxcdev)
      value = token%dtsets(idtset)%pawxcdev
    case (ab6_invars_dmftcheck)
      value = token%dtsets(idtset)%dmftcheck
    case (ab6_invars_prtwf)
      value = token%dtsets(idtset)%prtwf
    case (ab6_invars_gwcomp)
      value = token%dtsets(idtset)%gwcomp
    case (ab6_invars_mqgriddg)
      value = token%dtsets(idtset)%mqgriddg
    case (ab6_invars_inclvkb)
      value = token%dtsets(idtset)%inclvkb
    case (ab6_invars_gw_nstep)
      value = token%dtsets(idtset)%gw_nstep
    case (ab6_invars_tl_nprccg)
      value = token%dtsets(idtset)%tl_nprccg
    case (ab6_invars_pawprtwf)
      value = token%dtsets(idtset)%pawprtwf
    case (ab6_invars_ird1wf)
      value = token%dtsets(idtset)%ird1wf
    case (ab6_invars_rf2phon)
      value = token%dtsets(idtset)%rf2phon
    case (ab6_invars_rfddk)
      value = token%dtsets(idtset)%rfddk
    case (ab6_invars_nshsigx)
      value = token%dtsets(idtset)%nshsigx
    case (ab6_invars_nqpt)
      value = token%dtsets(idtset)%nqpt
    case (ab6_invars_pawlmix)
      value = token%dtsets(idtset)%pawlmix
    case (ab6_invars_prtlden)
      value = token%dtsets(idtset)%prtlden
    case (ab6_invars_optstress)
      value = token%dtsets(idtset)%optstress
    case (ab6_invars_positron)
      value = token%dtsets(idtset)%positron
    case (ab6_invars_smdelta)
      value = token%dtsets(idtset)%smdelta
    case (ab6_invars_getwfk)
      value = token%dtsets(idtset)%getwfk
    case (ab6_invars_diismemory)
      value = token%dtsets(idtset)%diismemory
    case (ab6_invars_getwfq)
      value = token%dtsets(idtset)%getwfq
    case (ab6_invars_nbandkss)
      value = token%dtsets(idtset)%nbandkss
    case (ab6_invars_gw_sigxcore)
      value = token%dtsets(idtset)%gw_sigxcore
    case (ab6_invars_pitransform)
      value = token%dtsets(idtset)%pitransform
    case (ab6_invars_pawmixdg)
      value = token%dtsets(idtset)%pawmixdg
    case (ab6_invars_prtdipole)
      value = token%dtsets(idtset)%prtdipole
    case (ab6_invars_istatshft)
      value = token%dtsets(idtset)%istatshft
    case (ab6_invars_npweps)
      value = token%dtsets(idtset)%npweps
    case (ab6_invars_bandpp)
      value = token%dtsets(idtset)%bandpp
    case (ab6_invars_awtr)
      value = token%dtsets(idtset)%awtr
    case (ab6_invars_nkptgw)
      value = token%dtsets(idtset)%nkptgw
    case (ab6_invars_npulayit)
      value = token%dtsets(idtset)%npulayit
    case (ab6_invars_dmft_nwli)
      value = token%dtsets(idtset)%dmft_nwli
    case (ab6_invars_rdmnb)
      value = token%dtsets(idtset)%rdmnb
    case (ab6_invars_rf3phon)
      value = token%dtsets(idtset)%rf3phon
    case (ab6_invars_nnos)
      value = token%dtsets(idtset)%nnos
    case (ab6_invars_gwrpacorr)
      value = token%dtsets(idtset)%gwrpacorr
    case (ab6_invars_dmatpuopt)
      value = token%dtsets(idtset)%dmatpuopt
    case (ab6_invars_natvshift)
      value = token%dtsets(idtset)%natvshift
    case (ab6_invars_signperm)
      value = token%dtsets(idtset)%signperm
    case (ab6_invars_prtgden)
      value = token%dtsets(idtset)%prtgden
    case (ab6_invars_prtstm)
      value = token%dtsets(idtset)%prtstm
    case (ab6_invars_berryopt)
      value = token%dtsets(idtset)%berryopt
    case (ab6_invars_ixcpositron)
      value = token%dtsets(idtset)%ixcpositron
    case (ab6_invars_getqps)
      value = token%dtsets(idtset)%getqps
    case (ab6_invars_timopt)
      value = token%dtsets(idtset)%timopt
    case (ab6_invars_mgfftdg)
      value = token%dtsets(idtset)%mgfftdg
    case (ab6_invars_getocc)
      value = token%dtsets(idtset)%getocc
    case (ab6_invars_prtposcar)
      value = token%dtsets(idtset)%prtposcar
    case (ab6_invars_nomegasrd)
      value = token%dtsets(idtset)%nomegasrd
    case (ab6_invars_iboxcut)
      value = token%dtsets(idtset)%iboxcut
    case (ab6_invars_pawntheta)
      value = token%dtsets(idtset)%pawntheta
    case (ab6_invars_usedmft)
      value = token%dtsets(idtset)%usedmft
    case (ab6_invars_nshiftk)
      value = token%dtsets(idtset)%nshiftk
    case (ab6_invars_prtelf)
      value = token%dtsets(idtset)%prtelf
    case (ab6_invars_prtcif)
      value = token%dtsets(idtset)%prtcif
    case (ab6_invars_intxc)
      value = token%dtsets(idtset)%intxc
    case (ab6_invars_bs_coulomb_term)
      value = token%dtsets(idtset)%bs_coulomb_term
    case (ab6_invars_ngroup_rf)
      value = token%dtsets(idtset)%ngroup_rf
    case (ab6_invars_spgaxor)
      value = token%dtsets(idtset)%spgaxor
    case (ab6_invars_get1den)
      value = token%dtsets(idtset)%get1den
    case (ab6_invars_gw_eet)
      value = token%dtsets(idtset)%gw_eet
    case (ab6_invars_dmftbandi)
      value = token%dtsets(idtset)%dmftbandi
    case (ab6_invars_irdqps)
      value = token%dtsets(idtset)%irdqps
    case (ab6_invars_optnlxccc)
      value = token%dtsets(idtset)%optnlxccc
    case (ab6_invars_jellslab)
      value = token%dtsets(idtset)%jellslab
    case (ab6_invars_ppmodel)
      value = token%dtsets(idtset)%ppmodel
    case (ab6_invars_pawoptmix)
      value = token%dtsets(idtset)%pawoptmix
    case (ab6_invars_optcell)
      value = token%dtsets(idtset)%optcell
    case (ab6_invars_pawprtdos)
      value = token%dtsets(idtset)%pawprtdos
    case (ab6_invars_pawfatbnd)
      value = token%dtsets(idtset)%pawfatbnd
    case (ab6_invars_gw_sctype)
      value = token%dtsets(idtset)%gw_sctype
    case (ab6_invars_kptopt)
      value = token%dtsets(idtset)%kptopt
    case (ab6_invars_w90iniprj)
      value = token%dtsets(idtset)%w90iniprj
    case (ab6_invars_nsppol)
      value = token%dtsets(idtset)%nsppol
    case (ab6_invars_nfreqre)
      value = token%dtsets(idtset)%nfreqre
    case (ab6_invars_ptgroupma)
      value = token%dtsets(idtset)%ptgroupma
    case (ab6_invars_pawnhatxc)
      value = token%dtsets(idtset)%pawnhatxc
    case (ab6_invars_pawcpxocc)
      value = token%dtsets(idtset)%pawcpxocc
    case (ab6_invars_optfreqsus)
      value = token%dtsets(idtset)%optfreqsus
    case (ab6_invars_pawusecp)
      value = token%dtsets(idtset)%pawusecp
    case (ab6_invars_prtxangst)
      value = token%dtsets(idtset)%prtxangst
    case (ab6_invars_userib)
      value = token%dtsets(idtset)%userib
    case (ab6_invars_idyson)
      value = token%dtsets(idtset)%idyson
    case (ab6_invars_prtpot)
      value = token%dtsets(idtset)%prtpot
    case (ab6_invars_optforces)
      value = token%dtsets(idtset)%optforces
    case (ab6_invars_brvltt)
      value = token%dtsets(idtset)%brvltt
    case (ab6_invars_ntime)
      value = token%dtsets(idtset)%ntime
    case (ab6_invars_prtvha)
      value = token%dtsets(idtset)%prtvha
    case (ab6_invars_npsp)
      value = token%dtsets(idtset)%npsp
    case (ab6_invars_symsigma)
      value = token%dtsets(idtset)%symsigma
    case (ab6_invars_prtvol)
      value = token%dtsets(idtset)%prtvol
    case (ab6_invars_nqptdm)
      value = token%dtsets(idtset)%nqptdm
    case (ab6_invars_iprcch)
      value = token%dtsets(idtset)%iprcch
    case (ab6_invars_prtbltztrp)
      value = token%dtsets(idtset)%prtbltztrp
    case (ab6_invars_isecur)
      value = token%dtsets(idtset)%isecur
    case (ab6_invars_natsph)
      value = token%dtsets(idtset)%natsph
    case (ab6_invars_getddk)
      value = token%dtsets(idtset)%getddk
    case (ab6_invars_nwfshist)
      value = token%dtsets(idtset)%nwfshist
    case (ab6_invars_nctime)
      value = token%dtsets(idtset)%nctime
    case (ab6_invars_nkpt)
      value = token%dtsets(idtset)%nkpt
    case (ab6_invars_nconeq)
      value = token%dtsets(idtset)%nconeq
    case (ab6_invars_vdw_xc)
      value = token%dtsets(idtset)%vdw_xc
    case (ab6_invars_ieig2rf)
      value = token%dtsets(idtset)%ieig2rf
    case (ab6_invars_getsuscep)
      value = token%dtsets(idtset)%getsuscep
    case (ab6_invars_icutcoul)
      value = token%dtsets(idtset)%icutcoul
    case (ab6_invars_prtspcur)
      value = token%dtsets(idtset)%prtspcur
    case (ab6_invars_natom)
      value = token%dtsets(idtset)%natom
    case (ab6_invars_prtdos)
      value = token%dtsets(idtset)%prtdos
    case (ab6_invars_spgorig)
      value = token%dtsets(idtset)%spgorig
    case (ab6_invars_nspden)
      value = token%dtsets(idtset)%nspden
    case (ab6_invars_symmorphi)
      value = token%dtsets(idtset)%symmorphi
    case (ab6_invars_prtnabla)
      value = token%dtsets(idtset)%prtnabla
    case (ab6_invars_pawnzlm)
      value = token%dtsets(idtset)%pawnzlm
    case (ab6_invars_nfreqim)
      value = token%dtsets(idtset)%nfreqim
    case (ab6_invars_dmft_nwlo)
      value = token%dtsets(idtset)%dmft_nwlo
    case (ab6_invars_prtwant)
      value = token%dtsets(idtset)%prtwant
    case (ab6_invars_mband)
      value = token%dtsets(idtset)%mband
    case (ab6_invars_getscr)
      value = token%dtsets(idtset)%getscr
    case (ab6_invars_optdriver)
      value = token%dtsets(idtset)%optdriver
    case (ab6_invars_prtbbb)
      value = token%dtsets(idtset)%prtbbb
    case (ab6_invars_prtgkk)
      value = token%dtsets(idtset)%prtgkk
    case (ab6_invars_pawprtvol)
      value = token%dtsets(idtset)%pawprtvol
    case (ab6_invars_pawnphi)
      value = token%dtsets(idtset)%pawnphi
    case (ab6_invars_recnrec)
      value = token%dtsets(idtset)%recnrec
    case (ab6_invars_bs_calctype)
      value = token%dtsets(idtset)%bs_calctype
    case (ab6_invars_posnstep)
      value = token%dtsets(idtset)%posnstep
    case (ab6_invars_recnpath)
      value = token%dtsets(idtset)%recnpath
    case (ab6_invars_istatr)
      value = token%dtsets(idtset)%istatr
    case (ab6_invars_natrd)
      value = token%dtsets(idtset)%natrd
    case (ab6_invars_getkss)
      value = token%dtsets(idtset)%getkss
    case (ab6_invars_nnsclo)
      value = token%dtsets(idtset)%nnsclo
    case (ab6_invars_xclevel)
      value = token%dtsets(idtset)%xclevel
    case (ab6_invars_ntyppure)
      value = token%dtsets(idtset)%ntyppure
    case (ab6_invars_mqgrid)
      value = token%dtsets(idtset)%mqgrid
    case (ab6_invars_prtdensph)
      value = token%dtsets(idtset)%prtdensph
    case (ab6_invars_dmft_iter)
      value = token%dtsets(idtset)%dmft_iter
    case (ab6_invars_getxcart)
      value = token%dtsets(idtset)%getxcart
    case (ab6_invars_chkexit)
      value = token%dtsets(idtset)%chkexit
    case (ab6_invars_vacnum)
      value = token%dtsets(idtset)%vacnum
    case (ab6_invars_ortalg)
      value = token%dtsets(idtset)%ortalg
    case (ab6_invars_fft_opt_lob)
      value = token%dtsets(idtset)%fft_opt_lob
    case (ab6_invars_mkqmem)
      value = token%dtsets(idtset)%mkqmem
    case (ab6_invars_rf2elfd)
      value = token%dtsets(idtset)%rf2elfd

    case default
      errno = AB6_ERROR_INVARS_ATT
    end select
  end subroutine ab6_invars_get_integer

  subroutine ab6_invars_get_real(dtsetsId, value, id, idtset, errno)
    integer, intent(in) :: dtsetsId
    integer, intent(in) :: id, idtset
    real(dp), intent(out) :: value
    integer, intent(out) :: errno

    type(dtsets_list), pointer :: token

    call get_token(token, dtsetsId)
    if (.not. associated(token)) then
      errno = AB6_ERROR_OBJ
      return
    end if
    if (idtset < 0 .or. idtset > size(token%dtsets)) then
      errno = AB6_ERROR_INVARS_ID
      return
    end if
    
    errno = AB6_NO_ERROR
    select case (id)
    case (ab6_invars_gwencomp)
      value = token%dtsets(idtset)%gwencomp
    case (ab6_invars_slabwsrad)
      value = token%dtsets(idtset)%slabwsrad
    case (ab6_invars_gw_toldfeig)
      value = token%dtsets(idtset)%gw_toldfeig
    case (ab6_invars_tolmxf)
      value = token%dtsets(idtset)%tolmxf
    case (ab6_invars_charge)
      value = token%dtsets(idtset)%charge
    case (ab6_invars_postoldff)
      value = token%dtsets(idtset)%postoldff
    case (ab6_invars_postoldfe)
      value = token%dtsets(idtset)%postoldfe
    case (ab6_invars_bxctmindg)
      value = token%dtsets(idtset)%bxctmindg
    case (ab6_invars_kptnrm)
      value = token%dtsets(idtset)%kptnrm
    case (ab6_invars_slabzbeg)
      value = token%dtsets(idtset)%slabzbeg
    case (ab6_invars_vis)
      value = token%dtsets(idtset)%vis
    case (ab6_invars_userre)
      value = token%dtsets(idtset)%userre
    case (ab6_invars_mdftemp)
      value = token%dtsets(idtset)%mdftemp
    case (ab6_invars_ecutsigx)
      value = token%dtsets(idtset)%ecutsigx
    case (ab6_invars_wvl_hgrid)
      value = token%dtsets(idtset)%wvl_hgrid
    case (ab6_invars_exchmix)
      value = token%dtsets(idtset)%exchmix
    case (ab6_invars_recefermi)
      value = token%dtsets(idtset)%recefermi
    case (ab6_invars_vacwidth)
      value = token%dtsets(idtset)%vacwidth
    case (ab6_invars_td_maxene)
      value = token%dtsets(idtset)%td_maxene
    case (ab6_invars_qptnrm)
      value = token%dtsets(idtset)%qptnrm
    case (ab6_invars_zcut)
      value = token%dtsets(idtset)%zcut
    case (ab6_invars_bmass)
      value = token%dtsets(idtset)%bmass
    case (ab6_invars_freqsuslo)
      value = token%dtsets(idtset)%freqsuslo
    case (ab6_invars_tolsym)
      value = token%dtsets(idtset)%tolsym
    case (ab6_invars_rectolden)
      value = token%dtsets(idtset)%rectolden
    case (ab6_invars_userrd)
      value = token%dtsets(idtset)%userrd
    case (ab6_invars_userra)
      value = token%dtsets(idtset)%userra
    case (ab6_invars_userrc)
      value = token%dtsets(idtset)%userrc
    case (ab6_invars_userrb)
      value = token%dtsets(idtset)%userrb
    case (ab6_invars_dmft_mxsf)
      value = token%dtsets(idtset)%dmft_mxsf
    case (ab6_invars_pawujv)
      value = token%dtsets(idtset)%pawujv
    case (ab6_invars_diemix)
      value = token%dtsets(idtset)%diemix
    case (ab6_invars_posocc)
      value = token%dtsets(idtset)%posocc
    case (ab6_invars_ecutsm)
      value = token%dtsets(idtset)%ecutsm
    case (ab6_invars_freqremax)
      value = token%dtsets(idtset)%freqremax
    case (ab6_invars_dielam)
      value = token%dtsets(idtset)%dielam
    case (ab6_invars_pawovlp)
      value = token%dtsets(idtset)%pawovlp
    case (ab6_invars_cpus)
      value = token%dtsets(idtset)%cpus
    case (ab6_invars_ecuteps)
      value = token%dtsets(idtset)%ecuteps
    case (ab6_invars_scphon_temp)
      value = token%dtsets(idtset)%scphon_temp
    case (ab6_invars_rcut)
      value = token%dtsets(idtset)%rcut
    case (ab6_invars_diemixmag)
      value = token%dtsets(idtset)%diemixmag
    case (ab6_invars_dielng)
      value = token%dtsets(idtset)%dielng
    case (ab6_invars_mdwall)
      value = token%dtsets(idtset)%mdwall
    case (ab6_invars_sciss)
      value = token%dtsets(idtset)%sciss
    case (ab6_invars_friction)
      value = token%dtsets(idtset)%friction
    case (ab6_invars_dosdeltae)
      value = token%dtsets(idtset)%dosdeltae
    case (ab6_invars_kptrlen)
      value = token%dtsets(idtset)%kptrlen
    case (ab6_invars_esmear)
      value = token%dtsets(idtset)%esmear
    case (ab6_invars_wvl_crmult)
      value = token%dtsets(idtset)%wvl_crmult
    case (ab6_invars_strfact)
      value = token%dtsets(idtset)%strfact
    case (ab6_invars_nelect)
      value = token%dtsets(idtset)%nelect
    case (ab6_invars_effmass)
      value = token%dtsets(idtset)%effmass
    case (ab6_invars_pawujrad)
      value = token%dtsets(idtset)%pawujrad
    case (ab6_invars_fixmom)
      value = token%dtsets(idtset)%fixmom
    case (ab6_invars_fband)
      value = token%dtsets(idtset)%fband
    case (ab6_invars_freqspmax)
      value = token%dtsets(idtset)%freqspmax
    case (ab6_invars_diegap)
      value = token%dtsets(idtset)%diegap
    case (ab6_invars_strprecon)
      value = token%dtsets(idtset)%strprecon
    case (ab6_invars_spnorbscl)
      value = token%dtsets(idtset)%spnorbscl
    case (ab6_invars_omegasimax)
      value = token%dtsets(idtset)%omegasimax
    case (ab6_invars_ecutwfn)
      value = token%dtsets(idtset)%ecutwfn
    case (ab6_invars_wvl_cpmult)
      value = token%dtsets(idtset)%wvl_cpmult
    case (ab6_invars_dilatmx)
      value = token%dtsets(idtset)%dilatmx
    case (ab6_invars_pawecutdg)
      value = token%dtsets(idtset)%pawecutdg
    case (ab6_invars_tolimg)
      value = token%dtsets(idtset)%tolimg
    case (ab6_invars_diemac)
      value = token%dtsets(idtset)%diemac
    case (ab6_invars_wvl_fpmult)
      value = token%dtsets(idtset)%wvl_fpmult
    case (ab6_invars_soenergy)
      value = token%dtsets(idtset)%soenergy
    case (ab6_invars_omegasrdmax)
      value = token%dtsets(idtset)%omegasrdmax
    case (ab6_invars_noseinert)
      value = token%dtsets(idtset)%noseinert
    case (ab6_invars_freqsusin)
      value = token%dtsets(idtset)%freqsusin
    case (ab6_invars_diecut)
      value = token%dtsets(idtset)%diecut
    case (ab6_invars_tolwfr)
      value = token%dtsets(idtset)%tolwfr
    case (ab6_invars_tphysel)
      value = token%dtsets(idtset)%tphysel
    case (ab6_invars_wvl_frmult)
      value = token%dtsets(idtset)%wvl_frmult
    case (ab6_invars_boxcutmin)
      value = token%dtsets(idtset)%boxcutmin
    case (ab6_invars_spbroad)
      value = token%dtsets(idtset)%spbroad
    case (ab6_invars_rhoqpmix)
      value = token%dtsets(idtset)%rhoqpmix
    case (ab6_invars_ecut)
      value = token%dtsets(idtset)%ecut
    case (ab6_invars_recrcut)
      value = token%dtsets(idtset)%recrcut
    case (ab6_invars_bs_haydock_tol)
      value = token%dtsets(idtset)%bs_haydock_tol
    case (ab6_invars_alpha)
      value = token%dtsets(idtset)%alpha
    case (ab6_invars_dtion)
      value = token%dtsets(idtset)%dtion
    case (ab6_invars_mditemp)
      value = token%dtsets(idtset)%mditemp
    case (ab6_invars_tolvrs)
      value = token%dtsets(idtset)%tolvrs
    case (ab6_invars_tsmear)
      value = token%dtsets(idtset)%tsmear
    case (ab6_invars_tl_radius)
      value = token%dtsets(idtset)%tl_radius
    case (ab6_invars_fxcartfactor)
      value = token%dtsets(idtset)%fxcartfactor
    case (ab6_invars_eshift)
      value = token%dtsets(idtset)%eshift
    case (ab6_invars_tolrff)
      value = token%dtsets(idtset)%tolrff
    case (ab6_invars_slabzend)
      value = token%dtsets(idtset)%slabzend
    case (ab6_invars_ppmfrq)
      value = token%dtsets(idtset)%ppmfrq
    case (ab6_invars_toldff)
      value = token%dtsets(idtset)%toldff
    case (ab6_invars_toldfe)
      value = token%dtsets(idtset)%toldfe
    case (ab6_invars_stmbias)
      value = token%dtsets(idtset)%stmbias

    case default
      errno = AB6_ERROR_INVARS_ATT
    end select
  end subroutine ab6_invars_get_real

  subroutine ab6_invars_get_integer_array(dtsetsId, values, n, id, idtset, errno)
    integer, intent(in) :: dtsetsId
    integer, intent(in) :: id, idtset, n
    integer, intent(out) :: values(n)
    integer, intent(out) :: errno

    type(dtsets_list), pointer :: token
    integer :: n_dt

    call get_token(token, dtsetsId)
    if (.not. associated(token)) then
      errno = AB6_ERROR_OBJ
      return
    end if
    if (idtset < 0 .or. idtset > size(token%dtsets)) then
      errno = AB6_ERROR_INVARS_ID
      return
    end if
    
    errno = AB6_NO_ERROR
    select case (id)
    case (ab6_invars_algalch)
      n_dt = product(shape(token%dtsets(idtset)%algalch))
      if (n_dt /= n) then
        errno = AB6_ERROR_INVARS_SIZE
      else
        values = reshape(token%dtsets(idtset)%algalch, (/ n_dt /))
      end if
    case (ab6_invars_kptrlatt)
      n_dt = product(shape(token%dtsets(idtset)%kptrlatt))
      if (n_dt /= n) then
        errno = AB6_ERROR_INVARS_SIZE
      else
        values = reshape(token%dtsets(idtset)%kptrlatt, (/ n_dt /))
      end if
    case (ab6_invars_nloalg)
      n_dt = product(shape(token%dtsets(idtset)%nloalg))
      if (n_dt /= n) then
        errno = AB6_ERROR_INVARS_SIZE
      else
        values = reshape(token%dtsets(idtset)%nloalg, (/ n_dt /))
      end if
    case (ab6_invars_qprtrb)
      n_dt = product(shape(token%dtsets(idtset)%qprtrb))
      if (n_dt /= n) then
        errno = AB6_ERROR_INVARS_SIZE
      else
        values = reshape(token%dtsets(idtset)%qprtrb, (/ n_dt /))
      end if
    case (ab6_invars_nband)
      n_dt = product(shape(token%dtsets(idtset)%nband))
      if (n_dt /= n) then
        errno = AB6_ERROR_INVARS_SIZE
      else
        values = reshape(token%dtsets(idtset)%nband, (/ n_dt /))
      end if
    case (ab6_invars_iatfix)
      n_dt = product(shape(token%dtsets(idtset)%iatfix))
      if (n_dt /= n) then
        errno = AB6_ERROR_INVARS_SIZE
      else
        values = reshape(token%dtsets(idtset)%iatfix, (/ n_dt /))
      end if
    case (ab6_invars_dynimage)
      n_dt = product(shape(token%dtsets(idtset)%dynimage))
      if (n_dt /= n) then
        errno = AB6_ERROR_INVARS_SIZE
      else
        values = reshape(token%dtsets(idtset)%dynimage, (/ n_dt /))
      end if
    case (ab6_invars_iatsph)
      n_dt = product(shape(token%dtsets(idtset)%iatsph))
      if (n_dt /= n) then
        errno = AB6_ERROR_INVARS_SIZE
      else
        values = reshape(token%dtsets(idtset)%iatsph, (/ n_dt /))
      end if
    case (ab6_invars_bdberry)
      n_dt = product(shape(token%dtsets(idtset)%bdberry))
      if (n_dt /= n) then
        errno = AB6_ERROR_INVARS_SIZE
      else
        values = reshape(token%dtsets(idtset)%bdberry, (/ n_dt /))
      end if
    case (ab6_invars_symafm)
      n_dt = product(shape(token%dtsets(idtset)%symafm))
      if (n_dt /= n) then
        errno = AB6_ERROR_INVARS_SIZE
      else
        values = reshape(token%dtsets(idtset)%symafm, (/ n_dt /))
      end if
    case (ab6_invars_typat)
      n_dt = product(shape(token%dtsets(idtset)%typat))
      if (n_dt /= n) then
        errno = AB6_ERROR_INVARS_SIZE
      else
        values = reshape(token%dtsets(idtset)%typat, (/ n_dt /))
      end if
    case (ab6_invars_ngfftdg)
      n_dt = product(shape(token%dtsets(idtset)%ngfftdg))
      if (n_dt /= n) then
        errno = AB6_ERROR_INVARS_SIZE
      else
        values = reshape(token%dtsets(idtset)%ngfftdg, (/ n_dt /))
      end if
    case (ab6_invars_rf3atpol)
      n_dt = product(shape(token%dtsets(idtset)%rf3atpol))
      if (n_dt /= n) then
        errno = AB6_ERROR_INVARS_SIZE
      else
        values = reshape(token%dtsets(idtset)%rf3atpol, (/ n_dt /))
      end if
    case (ab6_invars_bs_eh_basis_set)
      n_dt = product(shape(token%dtsets(idtset)%bs_eh_basis_set))
      if (n_dt /= n) then
        errno = AB6_ERROR_INVARS_SIZE
      else
        values = reshape(token%dtsets(idtset)%bs_eh_basis_set, (/ n_dt /))
      end if
    case (ab6_invars_rf1dir)
      n_dt = product(shape(token%dtsets(idtset)%rf1dir))
      if (n_dt /= n) then
        errno = AB6_ERROR_INVARS_SIZE
      else
        values = reshape(token%dtsets(idtset)%rf1dir, (/ n_dt /))
      end if
    case (ab6_invars_rfatpol)
      n_dt = product(shape(token%dtsets(idtset)%rfatpol))
      if (n_dt /= n) then
        errno = AB6_ERROR_INVARS_SIZE
      else
        values = reshape(token%dtsets(idtset)%rfatpol, (/ n_dt /))
      end if
    case (ab6_invars_rf1atpol)
      n_dt = product(shape(token%dtsets(idtset)%rf1atpol))
      if (n_dt /= n) then
        errno = AB6_ERROR_INVARS_SIZE
      else
        values = reshape(token%dtsets(idtset)%rf1atpol, (/ n_dt /))
      end if
    case (ab6_invars_bdgw)
      n_dt = product(shape(token%dtsets(idtset)%bdgw))
      if (n_dt /= n) then
        errno = AB6_ERROR_INVARS_SIZE
      else
        values = reshape(token%dtsets(idtset)%bdgw, (/ n_dt /))
      end if
    case (ab6_invars_scphon_supercell)
      n_dt = product(shape(token%dtsets(idtset)%scphon_supercell))
      if (n_dt /= n) then
        errno = AB6_ERROR_INVARS_SIZE
      else
        values = reshape(token%dtsets(idtset)%scphon_supercell, (/ n_dt /))
      end if
    case (ab6_invars_vdw_nwan)
      n_dt = product(shape(token%dtsets(idtset)%vdw_nwan))
      if (n_dt /= n) then
        errno = AB6_ERROR_INVARS_SIZE
      else
        values = reshape(token%dtsets(idtset)%vdw_nwan, (/ n_dt /))
      end if
    case (ab6_invars_rf3dir)
      n_dt = product(shape(token%dtsets(idtset)%rf3dir))
      if (n_dt /= n) then
        errno = AB6_ERROR_INVARS_SIZE
      else
        values = reshape(token%dtsets(idtset)%rf3dir, (/ n_dt /))
      end if
    case (ab6_invars_vdw_supercell)
      n_dt = product(shape(token%dtsets(idtset)%vdw_supercell))
      if (n_dt /= n) then
        errno = AB6_ERROR_INVARS_SIZE
      else
        values = reshape(token%dtsets(idtset)%vdw_supercell, (/ n_dt /))
      end if
    case (ab6_invars_rf2dir)
      n_dt = product(shape(token%dtsets(idtset)%rf2dir))
      if (n_dt /= n) then
        errno = AB6_ERROR_INVARS_SIZE
      else
        values = reshape(token%dtsets(idtset)%rf2dir, (/ n_dt /))
      end if
    case (ab6_invars_rfdir)
      n_dt = product(shape(token%dtsets(idtset)%rfdir))
      if (n_dt /= n) then
        errno = AB6_ERROR_INVARS_SIZE
      else
        values = reshape(token%dtsets(idtset)%rfdir, (/ n_dt /))
      end if
    case (ab6_invars_kberry)
      n_dt = product(shape(token%dtsets(idtset)%kberry))
      if (n_dt /= n) then
        errno = AB6_ERROR_INVARS_SIZE
      else
        values = reshape(token%dtsets(idtset)%kberry, (/ n_dt /))
      end if
    case (ab6_invars_supercell)
      n_dt = product(shape(token%dtsets(idtset)%supercell))
      if (n_dt /= n) then
        errno = AB6_ERROR_INVARS_SIZE
      else
        values = reshape(token%dtsets(idtset)%supercell, (/ n_dt /))
      end if
    case (ab6_invars_normpawu)
      n_dt = product(shape(token%dtsets(idtset)%normpawu))
      if (n_dt /= n) then
        errno = AB6_ERROR_INVARS_SIZE
      else
        values = reshape(token%dtsets(idtset)%normpawu, (/ n_dt /))
      end if
    case (ab6_invars_symrel)
      n_dt = product(shape(token%dtsets(idtset)%symrel))
      if (n_dt /= n) then
        errno = AB6_ERROR_INVARS_SIZE
      else
        values = reshape(token%dtsets(idtset)%symrel, (/ n_dt /))
      end if
    case (ab6_invars_so_psp)
      n_dt = product(shape(token%dtsets(idtset)%so_psp))
      if (n_dt /= n) then
        errno = AB6_ERROR_INVARS_SIZE
      else
        values = reshape(token%dtsets(idtset)%so_psp, (/ n_dt /))
      end if
    case (ab6_invars_istwfk)
      n_dt = product(shape(token%dtsets(idtset)%istwfk))
      if (n_dt /= n) then
        errno = AB6_ERROR_INVARS_SIZE
      else
        values = reshape(token%dtsets(idtset)%istwfk, (/ n_dt /))
      end if
    case (ab6_invars_rf2atpol)
      n_dt = product(shape(token%dtsets(idtset)%rf2atpol))
      if (n_dt /= n) then
        errno = AB6_ERROR_INVARS_SIZE
      else
        values = reshape(token%dtsets(idtset)%rf2atpol, (/ n_dt /))
      end if
    case (ab6_invars_lexexch)
      n_dt = product(shape(token%dtsets(idtset)%lexexch))
      if (n_dt /= n) then
        errno = AB6_ERROR_INVARS_SIZE
      else
        values = reshape(token%dtsets(idtset)%lexexch, (/ n_dt /))
      end if
    case (ab6_invars_lpawu)
      n_dt = product(shape(token%dtsets(idtset)%lpawu))
      if (n_dt /= n) then
        errno = AB6_ERROR_INVARS_SIZE
      else
        values = reshape(token%dtsets(idtset)%lpawu, (/ n_dt /))
      end if
    case (ab6_invars_ngfft)
      n_dt = product(shape(token%dtsets(idtset)%ngfft))
      if (n_dt /= n) then
        errno = AB6_ERROR_INVARS_SIZE
      else
        values = reshape(token%dtsets(idtset)%ngfft, (/ n_dt /))
      end if

    case default
      errno = AB6_ERROR_INVARS_ATT
    end select
  end subroutine ab6_invars_get_integer_array

  subroutine ab6_invars_get_real_array(dtsetsId, values, n, id, idtset, errno)
    integer, intent(in) :: dtsetsId
    integer, intent(in) :: id, idtset, n
    real(dp), intent(out) :: values(n)
    integer, intent(out) :: errno

    type(dtsets_list), pointer :: token
    integer :: n_dt

    call get_token(token, dtsetsId)
    if (.not. associated(token)) then
      errno = AB6_ERROR_OBJ
      return
    end if
    if (idtset < 0 .or. idtset > size(token%dtsets)) then
      errno = AB6_ERROR_INVARS_ID
      return
    end if
    
    errno = AB6_NO_ERROR
    select case (id)
    case (ab6_invars_boxcenter)
      n_dt = product(shape(token%dtsets(idtset)%boxcenter))
      if (n_dt /= n) then
        errno = AB6_ERROR_INVARS_SIZE
      else
        values = reshape(token%dtsets(idtset)%boxcenter, (/ n_dt /))
      end if
    case (ab6_invars_qptdm)
      n_dt = product(shape(token%dtsets(idtset)%qptdm))
      if (n_dt /= n) then
        errno = AB6_ERROR_INVARS_SIZE
      else
        values = reshape(token%dtsets(idtset)%qptdm, (/ n_dt /))
      end if
    case (ab6_invars_xred_orig)
      n_dt = product(shape(token%dtsets(idtset)%xred_orig))
      if (n_dt /= n) then
        errno = AB6_ERROR_INVARS_SIZE
      else
        values = reshape(token%dtsets(idtset)%xred_orig, (/ n_dt /))
      end if
    case (ab6_invars_kptns)
      n_dt = product(shape(token%dtsets(idtset)%kptns))
      if (n_dt /= n) then
        errno = AB6_ERROR_INVARS_SIZE
      else
        values = reshape(token%dtsets(idtset)%kptns, (/ n_dt /))
      end if
    case (ab6_invars_bs_freq_mesh)
      n_dt = product(shape(token%dtsets(idtset)%bs_freq_mesh))
      if (n_dt /= n) then
        errno = AB6_ERROR_INVARS_SIZE
      else
        values = reshape(token%dtsets(idtset)%bs_freq_mesh, (/ n_dt /))
      end if
    case (ab6_invars_ptcharge)
      n_dt = product(shape(token%dtsets(idtset)%ptcharge))
      if (n_dt /= n) then
        errno = AB6_ERROR_INVARS_SIZE
      else
        values = reshape(token%dtsets(idtset)%ptcharge, (/ n_dt /))
      end if
    case (ab6_invars_qmass)
      n_dt = product(shape(token%dtsets(idtset)%qmass))
      if (n_dt /= n) then
        errno = AB6_ERROR_INVARS_SIZE
      else
        values = reshape(token%dtsets(idtset)%qmass, (/ n_dt /))
      end if
    case (ab6_invars_gw_qlwl)
      n_dt = product(shape(token%dtsets(idtset)%gw_qlwl))
      if (n_dt /= n) then
        errno = AB6_ERROR_INVARS_SIZE
      else
        values = reshape(token%dtsets(idtset)%gw_qlwl, (/ n_dt /))
      end if
    case (ab6_invars_mixalch)
      n_dt = product(shape(token%dtsets(idtset)%mixalch))
      if (n_dt /= n) then
        errno = AB6_ERROR_INVARS_SIZE
      else
        values = reshape(token%dtsets(idtset)%mixalch, (/ n_dt /))
      end if
    case (ab6_invars_ratsph)
      n_dt = product(shape(token%dtsets(idtset)%ratsph))
      if (n_dt /= n) then
        errno = AB6_ERROR_INVARS_SIZE
      else
        values = reshape(token%dtsets(idtset)%ratsph, (/ n_dt /))
      end if
    case (ab6_invars_vcutgeo)
      n_dt = product(shape(token%dtsets(idtset)%vcutgeo))
      if (n_dt /= n) then
        errno = AB6_ERROR_INVARS_SIZE
      else
        values = reshape(token%dtsets(idtset)%vcutgeo, (/ n_dt /))
      end if
    case (ab6_invars_tnons)
      n_dt = product(shape(token%dtsets(idtset)%tnons))
      if (n_dt /= n) then
        errno = AB6_ERROR_INVARS_SIZE
      else
        values = reshape(token%dtsets(idtset)%tnons, (/ n_dt /))
      end if
    case (ab6_invars_kpt)
      n_dt = product(shape(token%dtsets(idtset)%kpt))
      if (n_dt /= n) then
        errno = AB6_ERROR_INVARS_SIZE
      else
        values = reshape(token%dtsets(idtset)%kpt, (/ n_dt /))
      end if
    case (ab6_invars_densty)
      n_dt = product(shape(token%dtsets(idtset)%densty))
      if (n_dt /= n) then
        errno = AB6_ERROR_INVARS_SIZE
      else
        values = reshape(token%dtsets(idtset)%densty, (/ n_dt /))
      end if
    case (ab6_invars_corecs)
      n_dt = product(shape(token%dtsets(idtset)%corecs))
      if (n_dt /= n) then
        errno = AB6_ERROR_INVARS_SIZE
      else
        values = reshape(token%dtsets(idtset)%corecs, (/ n_dt /))
      end if
    case (ab6_invars_atvshift)
      n_dt = product(shape(token%dtsets(idtset)%atvshift))
      if (n_dt /= n) then
        errno = AB6_ERROR_INVARS_SIZE
      else
        values = reshape(token%dtsets(idtset)%atvshift, (/ n_dt /))
      end if
    case (ab6_invars_genafm)
      n_dt = product(shape(token%dtsets(idtset)%genafm))
      if (n_dt /= n) then
        errno = AB6_ERROR_INVARS_SIZE
      else
        values = reshape(token%dtsets(idtset)%genafm, (/ n_dt /))
      end if
    case (ab6_invars_rprimd_orig)
      n_dt = product(shape(token%dtsets(idtset)%rprimd_orig))
      if (n_dt /= n) then
        errno = AB6_ERROR_INVARS_SIZE
      else
        values = reshape(token%dtsets(idtset)%rprimd_orig, (/ n_dt /))
      end if
    case (ab6_invars_qpt)
      n_dt = product(shape(token%dtsets(idtset)%qpt))
      if (n_dt /= n) then
        errno = AB6_ERROR_INVARS_SIZE
      else
        values = reshape(token%dtsets(idtset)%qpt, (/ n_dt /))
      end if
    case (ab6_invars_dmatpawu)
      n_dt = product(shape(token%dtsets(idtset)%dmatpawu))
      if (n_dt /= n) then
        errno = AB6_ERROR_INVARS_SIZE
      else
        values = reshape(token%dtsets(idtset)%dmatpawu, (/ n_dt /))
      end if
    case (ab6_invars_zeemanfield)
      n_dt = product(shape(token%dtsets(idtset)%zeemanfield))
      if (n_dt /= n) then
        errno = AB6_ERROR_INVARS_SIZE
      else
        values = reshape(token%dtsets(idtset)%zeemanfield, (/ n_dt /))
      end if
    case (ab6_invars_kptgw)
      n_dt = product(shape(token%dtsets(idtset)%kptgw))
      if (n_dt /= n) then
        errno = AB6_ERROR_INVARS_SIZE
      else
        values = reshape(token%dtsets(idtset)%kptgw, (/ n_dt /))
      end if
    case (ab6_invars_wtatcon)
      n_dt = product(shape(token%dtsets(idtset)%wtatcon))
      if (n_dt /= n) then
        errno = AB6_ERROR_INVARS_SIZE
      else
        values = reshape(token%dtsets(idtset)%wtatcon, (/ n_dt /))
      end if
    case (ab6_invars_occ_orig)
      n_dt = product(shape(token%dtsets(idtset)%occ_orig))
      if (n_dt /= n) then
        errno = AB6_ERROR_INVARS_SIZE
      else
        values = reshape(token%dtsets(idtset)%occ_orig, (/ n_dt /))
      end if
    case (ab6_invars_efield)
      n_dt = product(shape(token%dtsets(idtset)%efield))
      if (n_dt /= n) then
        errno = AB6_ERROR_INVARS_SIZE
      else
        values = reshape(token%dtsets(idtset)%efield, (/ n_dt /))
      end if
    case (ab6_invars_bfield)
      n_dt = product(shape(token%dtsets(idtset)%bfield))
      if (n_dt /= n) then
        errno = AB6_ERROR_INVARS_SIZE
      else
        values = reshape(token%dtsets(idtset)%bfield, (/ n_dt /))
      end if
    case (ab6_invars_jpawu)
      n_dt = product(shape(token%dtsets(idtset)%jpawu))
      if (n_dt /= n) then
        errno = AB6_ERROR_INVARS_SIZE
      else
        values = reshape(token%dtsets(idtset)%jpawu, (/ n_dt /))
      end if
    case (ab6_invars_wtk)
      n_dt = product(shape(token%dtsets(idtset)%wtk))
      if (n_dt /= n) then
        errno = AB6_ERROR_INVARS_SIZE
      else
        values = reshape(token%dtsets(idtset)%wtk, (/ n_dt /))
      end if
    case (ab6_invars_ziontypat)
      n_dt = product(shape(token%dtsets(idtset)%ziontypat))
      if (n_dt /= n) then
        errno = AB6_ERROR_INVARS_SIZE
      else
        values = reshape(token%dtsets(idtset)%ziontypat, (/ n_dt /))
      end if
    case (ab6_invars_acell_orig)
      n_dt = product(shape(token%dtsets(idtset)%acell_orig))
      if (n_dt /= n) then
        errno = AB6_ERROR_INVARS_SIZE
      else
        values = reshape(token%dtsets(idtset)%acell_orig, (/ n_dt /))
      end if
    case (ab6_invars_shiftk)
      n_dt = product(shape(token%dtsets(idtset)%shiftk))
      if (n_dt /= n) then
        errno = AB6_ERROR_INVARS_SIZE
      else
        values = reshape(token%dtsets(idtset)%shiftk, (/ n_dt /))
      end if
    case (ab6_invars_vel_orig)
      n_dt = product(shape(token%dtsets(idtset)%vel_orig))
      if (n_dt /= n) then
        errno = AB6_ERROR_INVARS_SIZE
      else
        values = reshape(token%dtsets(idtset)%vel_orig, (/ n_dt /))
      end if
    case (ab6_invars_strtarget)
      n_dt = product(shape(token%dtsets(idtset)%strtarget))
      if (n_dt /= n) then
        errno = AB6_ERROR_INVARS_SIZE
      else
        values = reshape(token%dtsets(idtset)%strtarget, (/ n_dt /))
      end if
    case (ab6_invars_amu)
      n_dt = product(shape(token%dtsets(idtset)%amu))
      if (n_dt /= n) then
        errno = AB6_ERROR_INVARS_SIZE
      else
        values = reshape(token%dtsets(idtset)%amu, (/ n_dt /))
      end if
    case (ab6_invars_rprim_orig)
      n_dt = product(shape(token%dtsets(idtset)%rprim_orig))
      if (n_dt /= n) then
        errno = AB6_ERROR_INVARS_SIZE
      else
        values = reshape(token%dtsets(idtset)%rprim_orig, (/ n_dt /))
      end if
    case (ab6_invars_bs_eh_cutoff)
      n_dt = product(shape(token%dtsets(idtset)%bs_eh_cutoff))
      if (n_dt /= n) then
        errno = AB6_ERROR_INVARS_SIZE
      else
        values = reshape(token%dtsets(idtset)%bs_eh_cutoff, (/ n_dt /))
      end if
    case (ab6_invars_spinat)
      n_dt = product(shape(token%dtsets(idtset)%spinat))
      if (n_dt /= n) then
        errno = AB6_ERROR_INVARS_SIZE
      else
        values = reshape(token%dtsets(idtset)%spinat, (/ n_dt /))
      end if
    case (ab6_invars_znucl)
      n_dt = product(shape(token%dtsets(idtset)%znucl))
      if (n_dt /= n) then
        errno = AB6_ERROR_INVARS_SIZE
      else
        values = reshape(token%dtsets(idtset)%znucl, (/ n_dt /))
      end if
    case (ab6_invars_upawu)
      n_dt = product(shape(token%dtsets(idtset)%upawu))
      if (n_dt /= n) then
        errno = AB6_ERROR_INVARS_SIZE
      else
        values = reshape(token%dtsets(idtset)%upawu, (/ n_dt /))
      end if
    case (ab6_invars_quadmom)
      n_dt = product(shape(token%dtsets(idtset)%quadmom))
      if (n_dt /= n) then
        errno = AB6_ERROR_INVARS_SIZE
      else
        values = reshape(token%dtsets(idtset)%quadmom, (/ n_dt /))
      end if
    case (ab6_invars_qptn)
      n_dt = product(shape(token%dtsets(idtset)%qptn))
      if (n_dt /= n) then
        errno = AB6_ERROR_INVARS_SIZE
      else
        values = reshape(token%dtsets(idtset)%qptn, (/ n_dt /))
      end if
    case (ab6_invars_vprtrb)
      n_dt = product(shape(token%dtsets(idtset)%vprtrb))
      if (n_dt /= n) then
        errno = AB6_ERROR_INVARS_SIZE
      else
        values = reshape(token%dtsets(idtset)%vprtrb, (/ n_dt /))
      end if

    case default
      errno = AB6_ERROR_INVARS_ATT
    end select
  end subroutine ab6_invars_get_real_array

  subroutine ab6_invars_get_shape(dtsetsId, dims, ndim, id, idtset, errno)
    integer, intent(in) :: dtsetsId
    integer, intent(in) :: id, idtset
    integer, intent(out) :: dims(7), ndim
    integer, intent(out) :: errno

    type(dtsets_list), pointer :: token

    call get_token(token, dtsetsId)
    if (.not. associated(token)) then
      errno = AB6_ERROR_OBJ
      return
    end if
    if (idtset < 0 .or. idtset > size(token%dtsets)) then
      errno = AB6_ERROR_INVARS_ID
      return
    end if
    
    errno = AB6_NO_ERROR
    select case (id)
    case (ab6_invars_boxcenter)
      ndim = size(shape(token%dtsets(idtset)%boxcenter))
      dims(1:ndim) = shape(token%dtsets(idtset)%boxcenter)
    case (ab6_invars_qptdm)
      ndim = size(shape(token%dtsets(idtset)%qptdm))
      dims(1:ndim) = shape(token%dtsets(idtset)%qptdm)
    case (ab6_invars_algalch)
      ndim = size(shape(token%dtsets(idtset)%algalch))
      dims(1:ndim) = shape(token%dtsets(idtset)%algalch)
    case (ab6_invars_kptrlatt)
      ndim = size(shape(token%dtsets(idtset)%kptrlatt))
      dims(1:ndim) = shape(token%dtsets(idtset)%kptrlatt)
    case (ab6_invars_nloalg)
      ndim = size(shape(token%dtsets(idtset)%nloalg))
      dims(1:ndim) = shape(token%dtsets(idtset)%nloalg)
    case (ab6_invars_xred_orig)
      ndim = size(shape(token%dtsets(idtset)%xred_orig))
      dims(1:ndim) = shape(token%dtsets(idtset)%xred_orig)
    case (ab6_invars_kptns)
      ndim = size(shape(token%dtsets(idtset)%kptns))
      dims(1:ndim) = shape(token%dtsets(idtset)%kptns)
    case (ab6_invars_bs_freq_mesh)
      ndim = size(shape(token%dtsets(idtset)%bs_freq_mesh))
      dims(1:ndim) = shape(token%dtsets(idtset)%bs_freq_mesh)
    case (ab6_invars_ptcharge)
      ndim = size(shape(token%dtsets(idtset)%ptcharge))
      dims(1:ndim) = shape(token%dtsets(idtset)%ptcharge)
    case (ab6_invars_qmass)
      ndim = size(shape(token%dtsets(idtset)%qmass))
      dims(1:ndim) = shape(token%dtsets(idtset)%qmass)
    case (ab6_invars_gw_qlwl)
      ndim = size(shape(token%dtsets(idtset)%gw_qlwl))
      dims(1:ndim) = shape(token%dtsets(idtset)%gw_qlwl)
    case (ab6_invars_qprtrb)
      ndim = size(shape(token%dtsets(idtset)%qprtrb))
      dims(1:ndim) = shape(token%dtsets(idtset)%qprtrb)
    case (ab6_invars_mixalch)
      ndim = size(shape(token%dtsets(idtset)%mixalch))
      dims(1:ndim) = shape(token%dtsets(idtset)%mixalch)
    case (ab6_invars_ratsph)
      ndim = size(shape(token%dtsets(idtset)%ratsph))
      dims(1:ndim) = shape(token%dtsets(idtset)%ratsph)
    case (ab6_invars_vcutgeo)
      ndim = size(shape(token%dtsets(idtset)%vcutgeo))
      dims(1:ndim) = shape(token%dtsets(idtset)%vcutgeo)
    case (ab6_invars_tnons)
      ndim = size(shape(token%dtsets(idtset)%tnons))
      dims(1:ndim) = shape(token%dtsets(idtset)%tnons)
    case (ab6_invars_nband)
      ndim = size(shape(token%dtsets(idtset)%nband))
      dims(1:ndim) = shape(token%dtsets(idtset)%nband)
    case (ab6_invars_iatfix)
      ndim = size(shape(token%dtsets(idtset)%iatfix))
      dims(1:ndim) = shape(token%dtsets(idtset)%iatfix)
    case (ab6_invars_dynimage)
      ndim = size(shape(token%dtsets(idtset)%dynimage))
      dims(1:ndim) = shape(token%dtsets(idtset)%dynimage)
    case (ab6_invars_iatsph)
      ndim = size(shape(token%dtsets(idtset)%iatsph))
      dims(1:ndim) = shape(token%dtsets(idtset)%iatsph)
    case (ab6_invars_kpt)
      ndim = size(shape(token%dtsets(idtset)%kpt))
      dims(1:ndim) = shape(token%dtsets(idtset)%kpt)
    case (ab6_invars_bdberry)
      ndim = size(shape(token%dtsets(idtset)%bdberry))
      dims(1:ndim) = shape(token%dtsets(idtset)%bdberry)
    case (ab6_invars_symafm)
      ndim = size(shape(token%dtsets(idtset)%symafm))
      dims(1:ndim) = shape(token%dtsets(idtset)%symafm)
    case (ab6_invars_densty)
      ndim = size(shape(token%dtsets(idtset)%densty))
      dims(1:ndim) = shape(token%dtsets(idtset)%densty)
    case (ab6_invars_corecs)
      ndim = size(shape(token%dtsets(idtset)%corecs))
      dims(1:ndim) = shape(token%dtsets(idtset)%corecs)
    case (ab6_invars_typat)
      ndim = size(shape(token%dtsets(idtset)%typat))
      dims(1:ndim) = shape(token%dtsets(idtset)%typat)
    case (ab6_invars_atvshift)
      ndim = size(shape(token%dtsets(idtset)%atvshift))
      dims(1:ndim) = shape(token%dtsets(idtset)%atvshift)
    case (ab6_invars_genafm)
      ndim = size(shape(token%dtsets(idtset)%genafm))
      dims(1:ndim) = shape(token%dtsets(idtset)%genafm)
    case (ab6_invars_ngfftdg)
      ndim = size(shape(token%dtsets(idtset)%ngfftdg))
      dims(1:ndim) = shape(token%dtsets(idtset)%ngfftdg)
    case (ab6_invars_rf3atpol)
      ndim = size(shape(token%dtsets(idtset)%rf3atpol))
      dims(1:ndim) = shape(token%dtsets(idtset)%rf3atpol)
    case (ab6_invars_rprimd_orig)
      ndim = size(shape(token%dtsets(idtset)%rprimd_orig))
      dims(1:ndim) = shape(token%dtsets(idtset)%rprimd_orig)
    case (ab6_invars_bs_eh_basis_set)
      ndim = size(shape(token%dtsets(idtset)%bs_eh_basis_set))
      dims(1:ndim) = shape(token%dtsets(idtset)%bs_eh_basis_set)
    case (ab6_invars_qpt)
      ndim = size(shape(token%dtsets(idtset)%qpt))
      dims(1:ndim) = shape(token%dtsets(idtset)%qpt)
    case (ab6_invars_dmatpawu)
      ndim = size(shape(token%dtsets(idtset)%dmatpawu))
      dims(1:ndim) = shape(token%dtsets(idtset)%dmatpawu)
    case (ab6_invars_rf1dir)
      ndim = size(shape(token%dtsets(idtset)%rf1dir))
      dims(1:ndim) = shape(token%dtsets(idtset)%rf1dir)
    case (ab6_invars_rfatpol)
      ndim = size(shape(token%dtsets(idtset)%rfatpol))
      dims(1:ndim) = shape(token%dtsets(idtset)%rfatpol)
    case (ab6_invars_zeemanfield)
      ndim = size(shape(token%dtsets(idtset)%zeemanfield))
      dims(1:ndim) = shape(token%dtsets(idtset)%zeemanfield)
    case (ab6_invars_rf1atpol)
      ndim = size(shape(token%dtsets(idtset)%rf1atpol))
      dims(1:ndim) = shape(token%dtsets(idtset)%rf1atpol)
    case (ab6_invars_kptgw)
      ndim = size(shape(token%dtsets(idtset)%kptgw))
      dims(1:ndim) = shape(token%dtsets(idtset)%kptgw)
    case (ab6_invars_bdgw)
      ndim = size(shape(token%dtsets(idtset)%bdgw))
      dims(1:ndim) = shape(token%dtsets(idtset)%bdgw)
    case (ab6_invars_wtatcon)
      ndim = size(shape(token%dtsets(idtset)%wtatcon))
      dims(1:ndim) = shape(token%dtsets(idtset)%wtatcon)
    case (ab6_invars_occ_orig)
      ndim = size(shape(token%dtsets(idtset)%occ_orig))
      dims(1:ndim) = shape(token%dtsets(idtset)%occ_orig)
    case (ab6_invars_efield)
      ndim = size(shape(token%dtsets(idtset)%efield))
      dims(1:ndim) = shape(token%dtsets(idtset)%efield)
    case (ab6_invars_bfield)
      ndim = size(shape(token%dtsets(idtset)%bfield))
      dims(1:ndim) = shape(token%dtsets(idtset)%bfield)
    case (ab6_invars_scphon_supercell)
      ndim = size(shape(token%dtsets(idtset)%scphon_supercell))
      dims(1:ndim) = shape(token%dtsets(idtset)%scphon_supercell)
    case (ab6_invars_jpawu)
      ndim = size(shape(token%dtsets(idtset)%jpawu))
      dims(1:ndim) = shape(token%dtsets(idtset)%jpawu)
    case (ab6_invars_wtk)
      ndim = size(shape(token%dtsets(idtset)%wtk))
      dims(1:ndim) = shape(token%dtsets(idtset)%wtk)
    case (ab6_invars_vdw_nwan)
      ndim = size(shape(token%dtsets(idtset)%vdw_nwan))
      dims(1:ndim) = shape(token%dtsets(idtset)%vdw_nwan)
    case (ab6_invars_ziontypat)
      ndim = size(shape(token%dtsets(idtset)%ziontypat))
      dims(1:ndim) = shape(token%dtsets(idtset)%ziontypat)
    case (ab6_invars_rf3dir)
      ndim = size(shape(token%dtsets(idtset)%rf3dir))
      dims(1:ndim) = shape(token%dtsets(idtset)%rf3dir)
    case (ab6_invars_vdw_supercell)
      ndim = size(shape(token%dtsets(idtset)%vdw_supercell))
      dims(1:ndim) = shape(token%dtsets(idtset)%vdw_supercell)
    case (ab6_invars_acell_orig)
      ndim = size(shape(token%dtsets(idtset)%acell_orig))
      dims(1:ndim) = shape(token%dtsets(idtset)%acell_orig)
    case (ab6_invars_rf2dir)
      ndim = size(shape(token%dtsets(idtset)%rf2dir))
      dims(1:ndim) = shape(token%dtsets(idtset)%rf2dir)
    case (ab6_invars_rfdir)
      ndim = size(shape(token%dtsets(idtset)%rfdir))
      dims(1:ndim) = shape(token%dtsets(idtset)%rfdir)
    case (ab6_invars_shiftk)
      ndim = size(shape(token%dtsets(idtset)%shiftk))
      dims(1:ndim) = shape(token%dtsets(idtset)%shiftk)
    case (ab6_invars_kberry)
      ndim = size(shape(token%dtsets(idtset)%kberry))
      dims(1:ndim) = shape(token%dtsets(idtset)%kberry)
    case (ab6_invars_vel_orig)
      ndim = size(shape(token%dtsets(idtset)%vel_orig))
      dims(1:ndim) = shape(token%dtsets(idtset)%vel_orig)
    case (ab6_invars_strtarget)
      ndim = size(shape(token%dtsets(idtset)%strtarget))
      dims(1:ndim) = shape(token%dtsets(idtset)%strtarget)
    case (ab6_invars_supercell)
      ndim = size(shape(token%dtsets(idtset)%supercell))
      dims(1:ndim) = shape(token%dtsets(idtset)%supercell)
    case (ab6_invars_normpawu)
      ndim = size(shape(token%dtsets(idtset)%normpawu))
      dims(1:ndim) = shape(token%dtsets(idtset)%normpawu)
    case (ab6_invars_amu)
      ndim = size(shape(token%dtsets(idtset)%amu))
      dims(1:ndim) = shape(token%dtsets(idtset)%amu)
    case (ab6_invars_symrel)
      ndim = size(shape(token%dtsets(idtset)%symrel))
      dims(1:ndim) = shape(token%dtsets(idtset)%symrel)
    case (ab6_invars_so_psp)
      ndim = size(shape(token%dtsets(idtset)%so_psp))
      dims(1:ndim) = shape(token%dtsets(idtset)%so_psp)
    case (ab6_invars_istwfk)
      ndim = size(shape(token%dtsets(idtset)%istwfk))
      dims(1:ndim) = shape(token%dtsets(idtset)%istwfk)
    case (ab6_invars_rprim_orig)
      ndim = size(shape(token%dtsets(idtset)%rprim_orig))
      dims(1:ndim) = shape(token%dtsets(idtset)%rprim_orig)
    case (ab6_invars_rf2atpol)
      ndim = size(shape(token%dtsets(idtset)%rf2atpol))
      dims(1:ndim) = shape(token%dtsets(idtset)%rf2atpol)
    case (ab6_invars_bs_eh_cutoff)
      ndim = size(shape(token%dtsets(idtset)%bs_eh_cutoff))
      dims(1:ndim) = shape(token%dtsets(idtset)%bs_eh_cutoff)
    case (ab6_invars_spinat)
      ndim = size(shape(token%dtsets(idtset)%spinat))
      dims(1:ndim) = shape(token%dtsets(idtset)%spinat)
    case (ab6_invars_znucl)
      ndim = size(shape(token%dtsets(idtset)%znucl))
      dims(1:ndim) = shape(token%dtsets(idtset)%znucl)
    case (ab6_invars_upawu)
      ndim = size(shape(token%dtsets(idtset)%upawu))
      dims(1:ndim) = shape(token%dtsets(idtset)%upawu)
    case (ab6_invars_lexexch)
      ndim = size(shape(token%dtsets(idtset)%lexexch))
      dims(1:ndim) = shape(token%dtsets(idtset)%lexexch)
    case (ab6_invars_quadmom)
      ndim = size(shape(token%dtsets(idtset)%quadmom))
      dims(1:ndim) = shape(token%dtsets(idtset)%quadmom)
    case (ab6_invars_qptn)
      ndim = size(shape(token%dtsets(idtset)%qptn))
      dims(1:ndim) = shape(token%dtsets(idtset)%qptn)
    case (ab6_invars_lpawu)
      ndim = size(shape(token%dtsets(idtset)%lpawu))
      dims(1:ndim) = shape(token%dtsets(idtset)%lpawu)
    case (ab6_invars_vprtrb)
      ndim = size(shape(token%dtsets(idtset)%vprtrb))
      dims(1:ndim) = shape(token%dtsets(idtset)%vprtrb)
    case (ab6_invars_ngfft)
      ndim = size(shape(token%dtsets(idtset)%ngfft))
      dims(1:ndim) = shape(token%dtsets(idtset)%ngfft)

    case default
      errno = AB6_ERROR_INVARS_ATT
    end select
  end subroutine ab6_invars_get_shape
