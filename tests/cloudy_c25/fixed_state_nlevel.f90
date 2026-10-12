program fixed_state_nlevel
!---------------------------------------------------------------------------
! Line cooling coefficients of MoCHII's ions on a list of fixed states, for
! the fixed-state cooling comparison against Cloudy c25.00
! (tests/cloudy_c25/fixed_state_cooling.py writes the input and reads the
! output; build and run with tests/cloudy_c25/run_fixed_state.sh).
!
! Input  (argument 1): first non-comment line nzone, then nzone rows
!                      T [K]  n_e [cm^-3]  n_HI [cm^-3]
! Output (argument 2): one row for each zone; for each ion (ION_EL, ION_ST) four
!                      columns, all from the compiled MoCHII modules:
! Model  (argument 3, optional): par%metal_excitation_model, the default when
!                      omitted.  Every metal file -- the zero-density fit, the
!                      n-level atom, the (T, n_e) table built from it, and the
!                      hydrogen-impact rates beside it -- is the one
!                      nlevel_mod::metal_excitation_data_file resolves for
!                      that model, i.e. the files a run with that selection
!                      reads.  H I is not selected by the model: its fit and
!                      atom are data/atomic/cooling_h_1.txt and nlevel_h_1.txt,
!                      as in cooling_mod.  The model is written in the output
!                      header.
!   lam_fit     Tier-1 zero-density fit, cooling_mod::cooling_eval of
!               cooling_<ion>.txt  [erg cm^3 s^-1; volume cooling
!               n_e n_ion lam_fit]  (par%cooling_model = 'low_density').
!               An ion the model gives no fit (mocassin_v2: C I, S III) is
!               evaluated as species_mod::metal_cooling evaluates it: the
!               n-level table, at its lowest density column here and at the
!               local n_e in lam_localne.
!   lam_localne the thermal-loop coefficient under par%cooling_model =
!               'local_ne': nlevel_cooling_mod::nlevel_cooling_apply on the
!               (T, n_e) table built from the same n-level atom
!               [erg cm^3 s^-1]; equals lam_fit for an ion with no n-level atom
!   P_e         sum of the n-level line emissivities, electron impact only,
!               nlevel_mod::nlevel_emissivities(atom, T, n_e, emis)
!               [erg s^-1 per ion]
!   P_eH        the same with the atomic-hydrogen impact rates of
!               hcoll_<ion>.txt in the same statistical
!               equilibrium, nlevel_emissivities(atom, T, n_e, emis, n_HI)
!               [erg s^-1 per ion]; equals P_e for an ion without them
! Run from the repository root (the atomic data are read from data/atomic).
!---------------------------------------------------------------------------
  use define
  use cooling_mod,        only : cooling_fit_type, cooling_load, cooling_eval
  use nlevel_mod,         only : nlevel_atom_type, nlevel_load, &
                                 nlevel_emissivities, metal_excitation_data_file
  use nlevel_cooling_mod, only : nlevel_cooling_add, nlevel_cooling_apply, &
                                 nlevel_cooling_lam_se
  implicit none

  integer, parameter :: NION = 15
  character(len=2), parameter :: ION_EL(NION) = [character(len=2) :: &
       'h', 'c', 'c', 'c', 'n', 'n', 'n', 'o', 'o', 'o', 'ne', 'ne', 's', 's', 's']
  integer, parameter :: ION_ST(NION) = [1, 1, 2, 3, 1, 2, 3, 1, 2, 3, 2, 3, 2, 3, 4]

  type(cooling_fit_type) :: fit(NION)
  type(nlevel_atom_type), allocatable :: atom(:)
  integer  :: idx(NION)
  logical  :: have_nl(NION), have_fit(NION)
  real(kind=wp) :: T, ne, nHI, lf
  real(kind=wp), allocatable :: out(:,:), emis(:)
  character(len=512) :: fin, fout, line, fname, model
  integer :: nz, iz, k, u, ios

  call get_command_argument(1, fin)
  call get_command_argument(2, fout)
  call get_command_argument(3, model)
  if (len_trim(fin) == 0 .or. len_trim(fout) == 0) then
     write(*,'(a)') 'usage: fixed_state_nlevel.x <state.txt> <coeff.txt> [chianti|mocassin_v2]'
     stop 1
  end if
  if (len_trim(model) > 0) par%metal_excitation_model = trim(model)

  par%atomic_dir    = 'data/atomic'
  par%cooling_model = 'local_ne'
  mpar%p_rank       = 1          ! quiet: no setup banners

  allocate(atom(NION))
  do k = 1, NION
     if (trim(ION_EL(k)) == 'h') then
        write(fname,'(a,i0,a)') trim(par%atomic_dir)//'/cooling_h_', ION_ST(k), '.txt'
     else
        fname = metal_excitation_data_file('cooling', trim(ION_EL(k)), ION_ST(k))
     end if
     inquire(file=trim(fname), exist=have_fit(k))
     if (have_fit(k)) call cooling_load(trim(fname), fit(k))
     if (trim(ION_EL(k)) == 'h') then
        write(fname,'(a,i0,a)') trim(par%atomic_dir)//'/nlevel_h_', ION_ST(k), '.txt'
     else
        fname = metal_excitation_data_file('nlevel', trim(ION_EL(k)), ION_ST(k))
     end if
     call nlevel_load(trim(fname), atom(k), have_nl(k))
     if (.not. have_nl(k)) then
        write(*,'(a)') 'ERROR: no n-level atom '//trim(fname)
        stop 1
     end if
     idx(k) = nlevel_cooling_add(trim(ION_EL(k)), ION_ST(k))
  end do

  open(newunit=u, file=trim(fin), status='old', action='read')
  do
     read(u,'(a)') line
     line = adjustl(line)
     if (len_trim(line) > 0 .and. line(1:1) /= '#') exit
  end do
  read(line,*) nz
  allocate(out(4*NION, nz))
  do iz = 1, nz
     read(u,*,iostat=ios) T, ne, nHI
     if (ios /= 0) then
        write(*,'(a,i0)') 'ERROR: short state file at zone ', iz
        stop 1
     end if
     do k = 1, NION
        if (allocated(emis)) deallocate(emis)
        allocate(emis(atom(k)%nrad))
        if (have_fit(k)) then
           lf = cooling_eval(fit(k), T)
           out(4*k-3, iz) = lf
           out(4*k-2, iz) = nlevel_cooling_apply(idx(k), lf, T, ne)
        else
           out(4*k-3, iz) = nlevel_cooling_lam_se(idx(k), T, tiny(1.0_wp))
           out(4*k-2, iz) = nlevel_cooling_lam_se(idx(k), T, ne)
        end if
        call nlevel_emissivities(atom(k), T, ne, emis)
        out(4*k-1, iz) = sum(emis(1:atom(k)%nrad))
        call nlevel_emissivities(atom(k), T, ne, emis, nHI)
        out(4*k,   iz) = sum(emis(1:atom(k)%nrad))
     end do
  end do
  close(u)

  open(newunit=u, file=trim(fout), status='replace', action='write')
  write(u,'(2a)') '# metal_excitation_model ', trim(par%metal_excitation_model)
  write(u,'(a)', advance='no') '# zone'
  do k = 1, NION
     write(line,'(a,a,i0)') trim(ION_EL(k)), '_', ION_ST(k)
     write(u,'(4(1x,a))', advance='no') trim(line)//':lam_fit', &
          trim(line)//':lam_localne', trim(line)//':P_e', trim(line)//':P_eH'
  end do
  write(u,*)
  do iz = 1, nz
     write(u,'(i5,*(1x,es24.16))') iz, out(:, iz)
  end do
  close(u)
end program fixed_state_nlevel
