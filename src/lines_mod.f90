module lines_mod
!---------------------------------------------------------------------------
! MoCHII: collisional line luminosities on the converged state (Tier 2).
!
! For every registry (element, stage) with a Tier-2 file
! data/atomic/nlevel_<el>_<stage>.txt, solve the n-level atom per leaf at
! the converged (T_e, n_e), multiply by the cascade stage density, and
! integrate over the grid:  L_line = sum_leaf j_line n_ion V  [erg/s].
! H beta from the case-B effective recombination coefficient
! alpha_eff(Hbeta) = 3.03e-14 (T/1e4)^-0.9 cm^3/s (accurate to a few % at
! 5-20 kK; Storey & Hummer tables arrive with the emission-map stage).
! Writes '<base>_lines.txt': element, stage, wavelength, L, L/L(Hbeta).
!
! With par%emis_output, also writes '<base>_emis' (HDF5/FITS): the same
! line set leaf by leaf as emissivity blocks emis_h (SH95 H I lines) and
! emis_<el>_<stage> [erg/s/cm^3; sum(emis*V) = L_line; row order matches
! the text table], wl_* wavelength arrays, and the state needed to build
! maps without the rates file: LeafXYZ, n_H, n_e, T_e, x_HI/HeI/HeII,
! and the x_<el>_stages fraction blocks.
!---------------------------------------------------------------------------
  use define
  use octree_mod,    only : amr_grid, leaf_half, leaf_cx, leaf_cy, leaf_cz
  use gas_state_mod, only : gas_nH, gas_xHI, gas_xHeI, gas_xHeII, &
                            gas_ne, gas_Te, gas_nleaf
  use species_mod,   only : n_elements, elem_name, elem_nstage, elem_abund, &
                            species_fractions, species_write
  use nlevel_mod
  use sh95_mod
  implicit none
  private

  public :: lines_write

contains

  !=========================================================================
  subroutine lines_write()
    use utility, only : get_base_name
    use iofile_mod
    implicit none
    type(nlevel_atom_type) :: atom
    type(io_file_type) :: efile
    character(len=256) :: fname
    character(len=192) :: outname
    character(len=64)  :: extname
    character(len=8)   :: nlsuf
    real(kind=wp) :: Lline(600), emis(600), frac(6), wl(600)
    real(kind=wp) :: T, ne, nH, vol, nion, LHb, t4, nHI, nHII
    real(kind=wp), allocatable :: em(:,:)
    integer :: ie, ist, il, k, unit, nl, status, ic
    logical :: ok, do_emis

    if (mpar%p_rank /= 0) return
    do_emis = par%emis_output
    status  = 0

    !--- the SH95 H I / He II tables now follow par%case_ab (case-A tables
    !--- are loaded for case A, case-B for case B), so they are consistent
    !--- with the ionization balance.  The He I Porter table is case B only
    !--- (no case-A Porter set); note the residual mismatch under case A.
    if (trim(par%case_ab) == 'A') write(*,'(a)') &
       ' LINE: note - He I Porter lines are case B (no case-A Porter'// &
       ' table); H I/He II lines use case A.'

    !--- emissivity file: open and write the state blocks first.
    if (do_emis) then
       call io_open_new(efile, trim(get_base_name(par%out_file))//'_emis'// &
                        trim(io_file_extension(par%file_format)), status)
       call emis_write_state(efile)
    end if

    !--- H I recombination lines (Storey & Hummer 1995 case B).
    call sh95_setup()
    block
      real(kind=wp) :: LH(12), wlh(12)
      integer :: kk, iHb, nlh
      nlh = sh95_nlines()
      LH = 0.0_wp;  iHb = 0
      do kk = 1, nlh
         if (trim(sh95_label(kk)) == 'Hbeta') iHb = kk
         wlh(kk) = sh95_wl(kk)
      end do
      if (do_emis) then
         allocate(em(nlh, gas_nleaf))
         em = 0.0_wp
      end if
      do il = 1, gas_nleaf
         nH = gas_nH(il)
         if (nH <= 0.0_wp) cycle
         T  = gas_Te(il);  ne = gas_ne(il)
         vol = (2.0_wp*leaf_half(il)*par%distance2cm)**3
         do kk = 1, nlh
            emis(kk) = ne*nH*(1.0_wp - gas_xHI(il))*sh95_emis(kk, T, ne)
            LH(kk)   = LH(kk) + emis(kk)*vol
         end do
         if (do_emis) em(1:nlh, il) = emis(1:nlh)
      end do
      LHb = LH(iHb)
      if (do_emis) then
         call io_append_image(efile, em, status, bitpix=-64)
         call io_put_keyword(efile,'EXTNAME','emis_h', &
              'SH95 H I line emissivities [erg/s/cm^3]',status)
         call io_append_image(efile, wlh(1:nlh), status, bitpix=-64)
         call io_put_keyword(efile,'EXTNAME','wl_h', &
              'H I line wavelengths [A]',status)
         deallocate(em)
      end if

      outname = trim(get_base_name(par%out_file))//'_lines.txt'
      open(newunit=unit, file=trim(outname), status='replace')
      write(unit,'(a)') '# MoCHII line luminosities (converged state)'
      write(unit,'(a)') '# H I recombination lines: Storey & Hummer (1995) '// &
         'case '//trim(par%case_ab)//', bilinear (log T, log ne) interpolation'
      write(unit,'(a)') '# (SH95 H I/He II lines use case '// &
         trim(par%case_ab)//'; He I Porter lines are case B)'
      write(unit,'(a)') '# metals: Tier-2 n-level solve per leaf'
      write(unit,'(a,es14.6,a)') '# L(Hbeta) = ', LHb, ' erg/s'
      write(unit,'(a)') '# elem stage  lambda[A]      L[erg/s]     L/L(Hbeta)'
      do kk = 1, nlh
         write(unit,'(a4,i4,f12.2,2es14.5,2x,a)') 'h', 1, wlh(kk), &
            LH(kk), LH(kk)/LHb, trim(sh95_label(kk))
      end do
    end block

    !--- He II recombination lines (SH95 Z = 2 case B; needs the He III
    !--- zone — the luminosities weight n_e n_HeIII).
    if (sh95_nlines(ion=2) > 0) then
       block
         real(kind=wp) :: LHe(12), wlhe(12), xHeIII
         integer :: kk, nlhe
         nlhe = sh95_nlines(ion=2)
         LHe = 0.0_wp
         do kk = 1, nlhe
            wlhe(kk) = sh95_wl(kk, ion=2)
         end do
         if (do_emis) then
            allocate(em(nlhe, gas_nleaf))
            em = 0.0_wp
         end if
         do il = 1, gas_nleaf
            nH = gas_nH(il)
            if (nH <= 0.0_wp) cycle
            xHeIII = max(0.0_wp, 1.0_wp - gas_xHeI(il) - gas_xHeII(il))
            if (xHeIII <= 1.0e-30_wp) cycle
            T  = gas_Te(il);  ne = gas_ne(il)
            vol = (2.0_wp*leaf_half(il)*par%distance2cm)**3
            do kk = 1, nlhe
               emis(kk) = ne*nH*par%He_abund*xHeIII &
                          *sh95_emis(kk, T, ne, ion=2)
               LHe(kk)  = LHe(kk) + emis(kk)*vol
            end do
            if (do_emis) em(1:nlhe, il) = emis(1:nlhe)
         end do
         do kk = 1, nlhe
            write(unit,'(a4,i4,f12.2,2es14.5,2x,a)') 'he', 2, wlhe(kk), &
               LHe(kk), LHe(kk)/LHb, trim(sh95_label(kk, ion=2))
         end do
         if (do_emis) then
            call io_append_image(efile, em, status, bitpix=-64)
            call io_put_keyword(efile,'EXTNAME','emis_heii', &
                 'SH95 He II line emissivities [erg/s/cm^3]',status)
            call io_append_image(efile, wlhe(1:nlhe), status, bitpix=-64)
            call io_put_keyword(efile,'EXTNAME','wl_heii', &
                 'He II line wavelengths [A]',status)
            deallocate(em)
         end if
       end block
    end if

    !--- He I recombination lines (Porter case B; emitted by He+ recombining,
    !--- so the luminosities weight n_e n_HeII).
    if (sh95_nlines(ion=3) > 0) then
       block
         real(kind=wp) :: LHe1(12), wlhe1(12)
         integer :: kk, nlhe1
         nlhe1 = sh95_nlines(ion=3)
         LHe1 = 0.0_wp
         do kk = 1, nlhe1
            wlhe1(kk) = sh95_wl(kk, ion=3)
         end do
         if (do_emis) then
            allocate(em(nlhe1, gas_nleaf))
            em = 0.0_wp
         end if
         do il = 1, gas_nleaf
            nH = gas_nH(il)
            if (nH <= 0.0_wp) cycle
            if (gas_xHeII(il) <= 1.0e-30_wp) cycle
            T  = gas_Te(il);  ne = gas_ne(il)
            vol = (2.0_wp*leaf_half(il)*par%distance2cm)**3
            do kk = 1, nlhe1
               emis(kk) = ne*nH*par%He_abund*gas_xHeII(il) &
                          *sh95_emis(kk, T, ne, ion=3)
               LHe1(kk) = LHe1(kk) + emis(kk)*vol
            end do
            if (do_emis) em(1:nlhe1, il) = emis(1:nlhe1)
         end do
         do kk = 1, nlhe1
            write(unit,'(a4,i4,f12.2,2es14.5,2x,a)') 'he', 1, wlhe1(kk), &
               LHe1(kk), LHe1(kk)/LHb, trim(sh95_label(kk, ion=3))
         end do
         if (do_emis) then
            call io_append_image(efile, em, status, bitpix=-64)
            call io_put_keyword(efile,'EXTNAME','emis_hei', &
                 'Porter He I line emissivities [erg/s/cm^3]',status)
            call io_append_image(efile, wlhe1(1:nlhe1), status, bitpix=-64)
            call io_put_keyword(efile,'EXTNAME','wl_hei', &
                 'He I line wavelengths [A]',status)
            deallocate(em)
         end if
       end block
    end if

    !--- optional collisional H I Balmer component (par%h_coll_effects):
    !--- solve the 25-level H I atom with collisional excitation from the
    !--- ground level and integrate the Balmer emission grouped into
    !--- Halpha/Hbeta/Hgamma windows.  Written as separate 'hc' rows so it
    !--- never contaminates the SH95 recombination H lines above.  Physics
    !--- caveat: nlevel_emissivities uses the atom's raw optically-thin
    !--- (case-A) radiative cascade; under case B the trapped Lyman photons
    !--- would enhance the Balmer collisional component somewhat.  This is a
    !--- first-cut collisional estimate, output separately by design.
    if (par%h_coll_effects) then
       call nlevel_load(trim(par%atomic_dir)//'/nlevel_h_1.txt', atom, ok)
       if (.not. ok) then
          write(*,'(a)') ' LINE: note - collisional H I file '// &
             'nlevel_h_1.txt not found; skipping the hc rows.'
       else
          block
            real(kind=wp) :: Lc(3), gemis(3), wlk, wlhc(3)
            integer :: kk, ig
            real(kind=wp), allocatable :: emc(:,:)
            nl = nlevel_nlines(atom)
            Lc = 0.0_wp
            wlhc = [6564.60_wp, 4862.65_wp, 4341.66_wp]  ! Halpha/Hbeta/Hgamma
            if (do_emis) then
               allocate(emc(3, gas_nleaf))
               emc = 0.0_wp
            end if
            do il = 1, gas_nleaf
               nH = gas_nH(il)
               if (nH <= 0.0_wp) cycle
               T   = gas_Te(il);  ne = gas_ne(il)
               nHI = nH*gas_xHI(il)
               if (nHI <= 1.0e-30_wp) cycle
               vol = (2.0_wp*leaf_half(il)*par%distance2cm)**3
               call nlevel_emissivities(atom, T, ne, emis)
               gemis = 0.0_wp
               do kk = 1, nl
                  wlk = nlevel_line_ident(atom, kk)
                  if      (wlk >= 6520.0_wp .and. wlk <= 6600.0_wp) then
                     ig = 1
                  else if (wlk >= 4840.0_wp .and. wlk <= 4880.0_wp) then
                     ig = 2
                  else if (wlk >= 4320.0_wp .and. wlk <= 4360.0_wp) then
                     ig = 3
                  else
                     cycle
                  end if
                  gemis(ig) = gemis(ig) + emis(kk)
               end do
               Lc(1:3) = Lc(1:3) + nHI*gemis(1:3)*vol
               if (do_emis) emc(1:3, il) = nHI*gemis(1:3)
            end do
            write(unit,'(a)') '# hc: collisional H I Balmer component'// &
               ' (excitation of neutral H from the ground level;'// &
               ' separate from the SH95 recombination H lines)'
            write(unit,'(a4,i4,f12.2,2es14.5,2x,a)') 'hc', 1, wlhc(1), &
               Lc(1), Lc(1)/LHb, 'Halpha_c'
            write(unit,'(a4,i4,f12.2,2es14.5,2x,a)') 'hc', 1, wlhc(2), &
               Lc(2), Lc(2)/LHb, 'Hbeta_c'
            write(unit,'(a4,i4,f12.2,2es14.5,2x,a)') 'hc', 1, wlhc(3), &
               Lc(3), Lc(3)/LHb, 'Hgamma_c'
            if (do_emis) then
               call io_append_image(efile, emc, status, bitpix=-64)
               call io_put_keyword(efile,'EXTNAME','emis_hc', &
                    'collisional H I Balmer emissivities [erg/s/cm^3]',status)
               call io_append_image(efile, wlhc, status, bitpix=-64)
               call io_put_keyword(efile,'EXTNAME','wl_hc', &
                    'collisional H I Balmer wavelengths [A]',status)
               deallocate(emc)
            end if
          end block
       end if
    end if

    do ie = 1, n_elements
       do ist = 1, elem_nstage(ie)
          !--- optional expanded Fe II/III model (par%fe_levels_full).
          nlsuf = ''
          if (par%fe_levels_full .and. trim(elem_name(ie)) == 'fe') &
             nlsuf = '_full'
          write(fname,'(a,a,a,i0,a,a)') trim(par%atomic_dir)//'/nlevel_', &
             trim(elem_name(ie)), '_', ist, trim(nlsuf), '.txt'
          call nlevel_load(trim(fname), atom, ok)
          if (.not. ok) cycle
          nl = nlevel_nlines(atom)
          Lline(1:nl) = 0.0_wp
          do k = 1, nl
             wl(k) = nlevel_line_ident(atom, k)
          end do
          if (do_emis) then
             allocate(em(nl, gas_nleaf))
             em = 0.0_wp
          end if
          do il = 1, gas_nleaf
             nH = gas_nH(il)
             if (nH <= 0.0_wp) cycle
             T = gas_Te(il);  ne = gas_ne(il)
             nHI = nH*gas_xHI(il);  nHII = nH - nHI
             call species_fractions(ie, il, T, ne, nHI, nHII, frac)
             nion = elem_abund(ie)*nH*frac(ist)
             if (nion <= 1.0e-30_wp) cycle
             vol = (2.0_wp*leaf_half(il)*par%distance2cm)**3
             call nlevel_emissivities(atom, T, ne, emis)
             Lline(1:nl) = Lline(1:nl) + emis(1:nl)*nion*vol
             if (do_emis) em(1:nl, il) = emis(1:nl)*nion
          end do
          do k = 1, nl
             write(unit,'(a4,i4,f12.2,2es14.5)') trim(elem_name(ie)), ist, &
                wl(k), Lline(k), Lline(k)/LHb
          end do
          if (do_emis) then
             write(extname,'(a,a,a,i0)') 'emis_', trim(elem_name(ie)), '_', ist
             call io_append_image(efile, em, status, bitpix=-64)
             call io_put_keyword(efile,'EXTNAME',trim(extname), &
                  'line emissivities [erg/s/cm^3]; sum(emis*V)=L',status)
             write(extname,'(a,a,a,i0)') 'wl_', trim(elem_name(ie)), '_', ist
             call io_append_image(efile, wl(1:nl), status, bitpix=-64)
             call io_put_keyword(efile,'EXTNAME',trim(extname), &
                  'line wavelengths [A]',status)
             deallocate(em)
          end if
       end do
    end do
    close(unit)
    write(*,'(2a)') ' LINE: Tier-2 line luminosities written to: ', trim(outname)
    if (do_emis) then
       call io_close(efile, status)
       write(*,'(2a)') ' LINE: leaf emissivities + state written to: ', &
          trim(get_base_name(par%out_file))//'_emis'// &
          trim(io_file_extension(par%file_format))
    end if
  end subroutine lines_write

  !=========================================================================
  ! State blocks of the emissivity file: everything a map needs without
  ! opening the rates file (positions, densities, temperature, fractions).
  !=========================================================================
  subroutine emis_write_state(efile)
    use iofile_mod
    implicit none
    type(io_file_type), intent(inout) :: efile
    real(kind=wp), allocatable :: tmp(:), lxyz(:,:)
    integer :: il, ic, status

    status = 0
    allocate(lxyz(gas_nleaf,3))
    do il = 1, gas_nleaf
       lxyz(il,1) = leaf_cx(il)
       lxyz(il,2) = leaf_cy(il)
       lxyz(il,3) = leaf_cz(il)
    end do
    call io_append_image(efile, lxyz, status, bitpix=-64)
    call io_put_keyword(efile,'EXTNAME','LeafXYZ', &
         'leaf center x,y,z (code units)',status)
    call io_put_keyword(efile,'DIST_CM', par%distance2cm, &
         'distance unit (cm)', status)
    deallocate(lxyz)

    allocate(tmp(gas_nleaf))
    do il = 1, gas_nleaf
       tmp(il) = 2.0_wp*leaf_half(il)        ! leaf edge length, code units
    end do
    call io_append_image(efile, tmp, status, bitpix=-64)
    call io_put_keyword(efile,'EXTNAME','LeafSize', &
         'leaf edge length (code units); V = (size*DIST_CM)^3',status)
    tmp = gas_nH(1:gas_nleaf)
    call io_append_image(efile, tmp, status, bitpix=-64)
    call io_put_keyword(efile,'EXTNAME','n_H','H density [cm^-3]',status)
    tmp = gas_ne(1:gas_nleaf)
    call io_append_image(efile, tmp, status, bitpix=-64)
    call io_put_keyword(efile,'EXTNAME','n_e','electron density [cm^-3]',status)
    tmp = gas_Te(1:gas_nleaf)
    call io_append_image(efile, tmp, status, bitpix=-64)
    call io_put_keyword(efile,'EXTNAME','T_e','electron temperature [K]',status)
    tmp = gas_xHI(1:gas_nleaf)
    call io_append_image(efile, tmp, status, bitpix=-64)
    call io_put_keyword(efile,'EXTNAME','x_HI','n_HI/n_H',status)
    tmp = gas_xHeI(1:gas_nleaf)
    call io_append_image(efile, tmp, status, bitpix=-64)
    call io_put_keyword(efile,'EXTNAME','x_HeI','n_HeI/n_He',status)
    tmp = gas_xHeII(1:gas_nleaf)
    call io_append_image(efile, tmp, status, bitpix=-64)
    call io_put_keyword(efile,'EXTNAME','x_HeII','n_HeII/n_He',status)
    deallocate(tmp)

    !--- metal stage fractions (same blocks as the rates file).
    if (par%use_metals) call species_write(efile)
  end subroutine emis_write_state

end module lines_mod
