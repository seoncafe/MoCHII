module lines_mod
!---------------------------------------------------------------------------
! MoCHII: collisional line luminosities on the converged state.
!
! For every registry (element, stage) with an n-level file
! data/atomic/excitation/<par%metal_excitation_model>/nlevel_<el>_<stage>.txt
! (H I: data/atomic/nlevel_h_1.txt), solve the n-level atom per leaf at
! the converged (T_e, n_e), multiply by the cascade stage density, and
! integrate over the grid:  L_line = sum_leaf j_line n_ion V  [erg/s].
! H beta (the normalization of every L/L(Hbeta) column) and the other H I
! and He II recombination lines from the Storey & Hummer (1995) tables
! (sh95_mod), He I from Porter et al. (2012, 2013).  The C I, C II and O I
! ground-term fine structure also carries atomic-hydrogen impact
! (nlevel_mod, hcoll_<el>_<stage>.txt).
! Writes '<base>_lines.txt': element, stage, wavelength, L, L/L(Hbeta).
!
! With par%emis_output, also writes '<base>_emis' (HDF5/FITS): the same
! line set leaf by leaf as emissivity blocks emis_h (SH95 H I lines) and
! emis_<el>_<stage> [erg/s/cm^3; sum(emis*V) = L_line; for the metals the
! lines with L >= par%emis_line_floor L(Hbeta), in text-table order], wl_*
! wavelength arrays, and the state needed to build
! maps without the rates file: LeafXYZ, n_H, n_e, T_e, x_HI/HII/HeI/HeII/HeIII,
! and the x_<el>_stages fraction blocks.
!---------------------------------------------------------------------------
  use define
  use octree_mod,    only : amr_grid, leaf_half, leaf_cx, leaf_cy, leaf_cz
  use gas_state_mod, only : gas_nH, gas_xHI, gas_xHII, gas_xHeI, gas_xHeII, &
                            gas_xHeIII, &
                            gas_ne, gas_Te, gas_nleaf, &
                            gas_leaf_owner_first, gas_leaf_owner_stride
  use mpi
  use species_mod,   only : n_elements, elem_name, elem_nstage, elem_abund, &
                            species_fractions, species_write
  use nlevel_mod
  use sh95_mod
  use hydrogenic_recombination_mod, only : lyman_line_case
  implicit none
  private

  public :: lines_write

contains

  !=========================================================================
  subroutine lines_write()
    use utility, only : get_base_name, fatal_error
    use iofile_mod
    implicit none
    type(nlevel_atom_type) :: atom
    type(io_file_type) :: efile
    character(len=PATH_LEN) :: fname
    character(len=PATH_LEN) :: outname
    character(len=64)  :: extname
    real(kind=wp) :: frac(MAX_ION_STAGES)
    real(kind=wp), allocatable :: Lline(:), emis(:), wl(:)
    real(kind=wp) :: T, ne, nH, vol, nion, LHb, t4, nHI, nHII
    real(kind=wp), allocatable :: em(:,:), emk(:,:)
    integer :: ie, ist, il, k, unit, nl, status, ic, nsel, ierr
    integer :: i0, nk, nchunk
    !--- leaf emissivities of the metal lines are summed onto rank 0 in
    !--- chunks of at most this many values, so only rank 0 holds the whole
    !--- (line, leaf) table (every rank held it until 2026-10-12: ~5 GB per
    !--- rank for ~300 lines on 2.1M leaves)
    integer, parameter :: EMIS_CHUNK_VALUES = 2**24
    integer, allocatable :: sel(:)
    logical :: ok, do_emis

    do_emis = par%emis_output
    status  = 0
    LHb = 0.0_wp
    !--- the recombination lines and every file operation on rank 0; the
    !--- metal n-level lines below on all ranks, leaves divided
    if (mpar%p_rank == 0) then

    !--- the SH95 H I / He II tables follow par%lyman_lines, the fate of the
    !--- Lyman-LINE photons, which is deliberately independent of
    !--- par%recombination_continuum: a nebula whose diffuse continuum is
    !--- transported still traps its Lyman lines.  The He I Porter table is
    !--- case B only (no case-A Porter set).
    if (trim(par%lyman_lines) == 'escaping') &
       write(*,'(a)') &
       ' LINE: note - H I/He II lines use the case-A SH95 tables'// &
       ' (He I Porter lines are case B; no case-A Porter table).'

    !--- emissivity file: open and write the state blocks first.
    if (do_emis) then
       call io_open_new(efile, trim(get_base_name(par%out_file))//'_emis'// &
                        trim(io_file_extension(par%file_format)), status)
       call emis_write_state(efile, status)
    end if

    !--- H I recombination lines (Storey & Hummer 1995 case B).
    call sh95_setup()
    allocate(emis(12))
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
            emis(kk) = ne*nH*gas_xHII(il)*sh95_emis(kk, T, ne)
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
         'case '//trim(lyman_line_case())//', bilinear (log T, log ne) '// &
         'interpolation'
      write(unit,'(a)') '# (SH95 H I/He II lines use case '// &
         trim(lyman_line_case())//'; He I Porter lines are case B)'
      write(unit,'(a)') '# metals: n-level solve leaf by leaf'
      write(unit,'(a,es14.6,a)') '# L(Hbeta) = ', LHb, ' erg/s'
      write(unit,'(a)') '# elem stage  lambda[A]      L[erg/s]     L/L(Hbeta)'
      do kk = 1, nlh
         write(unit,'(a4,i4,f12.2,2es15.5e3,2x,a)') 'h', 1, wlh(kk), &
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
            xHeIII = gas_xHeIII(il)
            if (xHeIII <= 1.0e-30_wp) cycle
            T  = gas_Te(il);  ne = gas_ne(il)
            vol = (2.0_wp*leaf_half(il)*par%distance2cm)**3
            do kk = 1, nlhe
               emis(kk) = ne*nH*par%abund_He*xHeIII &
                          *sh95_emis(kk, T, ne, ion=2)
               LHe(kk)  = LHe(kk) + emis(kk)*vol
            end do
            if (do_emis) em(1:nlhe, il) = emis(1:nlhe)
         end do
         do kk = 1, nlhe
            write(unit,'(a4,i4,f12.2,2es15.5e3,2x,a)') 'he', 2, wlhe(kk), &
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

    !--- He I lines (Porter case-B collisional-radiative emissivities: the
    !--- table already includes collisional excitation out of the 2^3S
    !--- metastable, so it is density-dependent, not a pure recombination
    !--- cascade; luminosities weight n_e n_HeII, He+ setting the triplet
    !--- reservoir through recombination).
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
               emis(kk) = ne*nH*par%abund_He*gas_xHeII(il) &
                          *sh95_emis(kk, T, ne, ion=3)
               LHe1(kk) = LHe1(kk) + emis(kk)*vol
            end do
            if (do_emis) em(1:nlhe1, il) = emis(1:nlhe1)
         end do
         do kk = 1, nlhe1
            write(unit,'(a4,i4,f12.2,2es15.5e3,2x,a)') 'he', 1, wlhe1(kk), &
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
            deallocate(emis)
            allocate(emis(atom%nrad))
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
            write(unit,'(a4,i4,f12.2,2es15.5e3,2x,a)') 'hc', 1, wlhc(1), &
               Lc(1), Lc(1)/LHb, 'Halpha_c'
            write(unit,'(a4,i4,f12.2,2es15.5e3,2x,a)') 'hc', 1, wlhc(2), &
               Lc(2), Lc(2)/LHb, 'Hbeta_c'
            write(unit,'(a4,i4,f12.2,2es15.5e3,2x,a)') 'hc', 1, wlhc(3), &
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

    end if
    call MPI_BCAST(LHb, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

    !--- metal lines: every rank solves the atom in its own leaves
    !--- (gas_leaf_owner_first/stride, as the thermal solve) and the line
    !--- luminosities are summed on rank 0.  With par%emis_output, the leaf
    !--- emissivities of the lines with L >= par%emis_line_floor L(Hbeta)
    !--- are recomputed in a second pass and summed on rank 0 (a full atom
    !--- carries up to ~1e4 lines; all of them leaf by leaf would be 22 GB
    !--- for Fe III on a 262144-leaf grid).
    do ie = 1, n_elements
       do ist = 1, elem_nstage(ie)
          fname=metal_excitation_data_file('nlevel',trim(elem_name(ie)),ist)
          call nlevel_load(trim(fname), atom, ok)
          if (.not. ok) cycle
          nl = nlevel_nlines(atom)
          if (allocated(Lline)) deallocate(Lline, emis, wl)
          if (allocated(emis)) deallocate(emis)
          allocate(Lline(atom%nrad), emis(atom%nrad), wl(atom%nrad))
          Lline(1:nl) = 0.0_wp
          do k = 1, nl
             wl(k) = nlevel_line_ident(atom, k)
          end do
          do il = gas_leaf_owner_first(), gas_nleaf, gas_leaf_owner_stride()
             if (.not. metal_ion_density(il, nion, T, ne, nHI, vol)) cycle
             call atom_emissivities(T, ne, nHI)
             Lline(1:nl) = Lline(1:nl) + emis(1:nl)*nion*vol
          end do
          call MPI_ALLREDUCE(MPI_IN_PLACE, Lline, nl, MPI_DOUBLE_PRECISION, &
                             MPI_SUM, MPI_COMM_WORLD, ierr)
          if (mpar%p_rank == 0) then
             do k = 1, nl
                write(unit,'(a4,i4,f12.2,2es15.5e3)') trim(elem_name(ie)), ist, &
                   wl(k), Lline(k), Lline(k)/LHb
             end do
          end if
          if (do_emis) then
             nsel = count(Lline(1:nl) >= par%emis_line_floor*LHb .and. Lline(1:nl) > 0.0_wp)
             if (nsel > 0) then
                allocate(sel(nsel))
                sel = pack([(k, k = 1, nl)], Lline(1:nl) >= par%emis_line_floor*LHb &
                           .and. Lline(1:nl) > 0.0_wp)
                if (mpar%p_rank == 0) then
                   allocate(em(nsel, gas_nleaf))
                else
                   allocate(em(1, 1))
                end if
                nchunk = max(1, min(gas_nleaf, EMIS_CHUNK_VALUES/nsel))
                allocate(emk(nsel, nchunk))
                do i0 = 1, gas_nleaf, nchunk
                   nk = min(nchunk, gas_nleaf - i0 + 1)
                   emk(:, 1:nk) = 0.0_wp
                   !--- this rank's leaves in [i0, i0+nk-1]
                   do il = i0 + modulo(gas_leaf_owner_first() - i0, gas_leaf_owner_stride()), &
                           i0 + nk - 1, gas_leaf_owner_stride()
                      if (.not. metal_ion_density(il, nion, T, ne, nHI, vol)) cycle
                      call atom_emissivities(T, ne, nHI)
                      emk(:, il - i0 + 1) = emis(sel)*nion
                   end do
                   !--- one nonzero contribution per leaf: the sum is exact
                   if (mpar%p_rank == 0) then
                      call MPI_REDUCE(emk, em(1, i0), nsel*nk, MPI_DOUBLE_PRECISION, &
                                      MPI_SUM, 0, MPI_COMM_WORLD, ierr)
                   else
                      call MPI_REDUCE(emk, emk, nsel*nk, MPI_DOUBLE_PRECISION, &
                                      MPI_SUM, 0, MPI_COMM_WORLD, ierr)
                   end if
                end do
                deallocate(emk)
                if (mpar%p_rank == 0) then
                   write(extname,'(a,a,a,i0)') 'emis_', trim(elem_name(ie)), '_', ist
                   call io_append_image(efile, em, status, bitpix=-64)
                   call io_put_keyword(efile,'EXTNAME',trim(extname), &
                        'line emissivities [erg/s/cm^3]; sum(emis*V)=L',status)
                   write(extname,'(a,a,a,i0)') 'wl_', trim(elem_name(ie)), '_', ist
                   call io_append_image(efile, wl(sel), status, bitpix=-64)
                   call io_put_keyword(efile,'EXTNAME',trim(extname), &
                        'line wavelengths [A] (lines with L >= emis_line_floor L(Hbeta))',status)
                end if
                deallocate(sel, em)
             end if
          end if
       end do
    end do
    if (mpar%p_rank /= 0) return
    close(unit)
    write(*,'(2a)') ' LINE: line luminosities written to: ', trim(outname)
    if (do_emis) then
       call io_close(efile, status)
       if (status /= 0) call fatal_error('failed to write emissivity file')
       write(*,'(2a)') ' LINE: leaf emissivities + state written to: ', &
          trim(get_base_name(par%out_file))//'_emis'// &
          trim(io_file_extension(par%file_format))
    end if
  contains

    !--- n(ion) and the leaf state for (element ie, stage ist) in leaf il;
    !--- .false. where the leaf is empty or the ion absent
    logical function metal_ion_density(il, nion, T, ne, nHI, vol) result(present_ion)
      integer,       intent(in)  :: il
      real(kind=wp), intent(out) :: nion, T, ne, nHI, vol
      real(kind=wp) :: nH, nHII
      present_ion = .false.
      nion = 0.0_wp;  T = 0.0_wp;  ne = 0.0_wp;  nHI = 0.0_wp;  vol = 0.0_wp
      nH = gas_nH(il)
      if (nH <= 0.0_wp) return
      T = gas_Te(il);  ne = gas_ne(il)
      nHI = nH*gas_xHI(il);  nHII = nH*gas_xHII(il)
      call species_fractions(ie, il, T, ne, nHI, nHII, frac)
      nion = elem_abund(ie)*nH*frac(ist)
      if (nion <= 1.0e-30_wp) return
      vol = (2.0_wp*leaf_half(il)*par%distance2cm)**3
      present_ion = .true.
    end function metal_ion_density

    !--- emis(1:nl) of the loaded atom; nHI adds atomic-hydrogen impact on
    !--- the ground-term fine structure of C I, C II, O I
    !--- (hcoll_<el>_<stage>.txt)
    subroutine atom_emissivities(T, ne, nHI)
      real(kind=wp), intent(in) :: T, ne, nHI
      if (par%neutral_hydrogen_excitation) then
         call nlevel_emissivities(atom, T, ne, emis, nHI)
      else
         call nlevel_emissivities(atom, T, ne, emis)
      end if
    end subroutine atom_emissivities

  end subroutine lines_write

  !=========================================================================
  ! State blocks of the emissivity file: everything a map needs without
  ! opening the rates file (positions, densities, temperature, fractions).
  !=========================================================================
  subroutine emis_write_state(efile, status)
    use iofile_mod
    implicit none
    type(io_file_type), intent(inout) :: efile
    integer, intent(inout) :: status
    real(kind=wp), allocatable :: tmp(:), lxyz(:,:)
    integer :: il, ic

    if (status /= 0) return
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
    tmp = gas_xHII(1:gas_nleaf)
    call io_append_image(efile, tmp, status, bitpix=-64)
    call io_put_keyword(efile,'EXTNAME','x_HII','n_HII/n_H',status)
    tmp = gas_xHeI(1:gas_nleaf)
    call io_append_image(efile, tmp, status, bitpix=-64)
    call io_put_keyword(efile,'EXTNAME','x_HeI','n_HeI/n_He',status)
    tmp = gas_xHeII(1:gas_nleaf)
    call io_append_image(efile, tmp, status, bitpix=-64)
    call io_put_keyword(efile,'EXTNAME','x_HeII','n_HeII/n_He',status)
    tmp = gas_xHeIII(1:gas_nleaf)
    call io_append_image(efile, tmp, status, bitpix=-64)
    call io_put_keyword(efile,'EXTNAME','x_HeIII','n_HeIII/n_He',status)
    deallocate(tmp)

    !--- metal stage fractions (same blocks as the rates file).
    if (par%use_metals) call species_write(efile, status)
  end subroutine emis_write_state

end module lines_mod
