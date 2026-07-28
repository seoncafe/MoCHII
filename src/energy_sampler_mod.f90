module energy_sampler_mod
!---------------------------------------------------------------------------
! Continuous photon-energy samplers.
!
! A sampler stores a piecewise-linear, non-negative energy-luminosity density
! and its analytically integrated CDF.  Inversion is therefore exact for a
! tabulated spectrum's stated linear interpolation model.  Planck spectra use
! the same representation after adaptive refinement, independently of the
! diagnostic ionizing-band grid.
!---------------------------------------------------------------------------
  use, intrinsic :: iso_fortran_env, only : real64
  implicit none
  private

  integer, parameter :: wp = real64
  integer, parameter :: MAX_PLANCK_KNOTS = 200000
  real(wp), parameter :: KB_EV = 8.617333262145e-5_wp
  real(wp), parameter :: ZETA4 = 1.0823232337111381915_wp

  type, public :: energy_sampler_type
     real(wp), allocatable :: energy(:)
     real(wp), allocatable :: density(:)
     real(wp), allocatable :: cdf(:)
     real(wp) :: total_integral = 0.0_wp
  end type energy_sampler_type

  public :: build_tabulated_sampler, build_tabulated_band_sampler
  public :: build_planck_sampler
  public :: sample_energy_cdf, energy_cdf_value
  public :: planck_energy_shape, sample_full_planck_bc

contains

  subroutine build_tabulated_sampler(sampler, energy, density)
    type(energy_sampler_type), intent(out) :: sampler
    real(wp), intent(in) :: energy(:), density(:)
    integer :: i, n
    real(wp) :: area, total

    n = size(energy)
    if (n < 2 .or. size(density) /= n) &
       error stop 'energy sampler: table needs at least two matching columns'
    if (any(density < 0.0_wp)) &
       error stop 'energy sampler: negative spectral density'
    do i = 1, n - 1
       if (energy(i+1) <= energy(i)) &
          error stop 'energy sampler: energies must be strictly ascending'
    end do

    allocate(sampler%energy(n), sampler%density(n), sampler%cdf(n))
    sampler%energy = energy
    sampler%density = density
    sampler%cdf(1) = 0.0_wp
    do i = 1, n - 1
       area = 0.5_wp * (density(i) + density(i+1)) * &
              (energy(i+1) - energy(i))
       sampler%cdf(i+1) = sampler%cdf(i) + area
    end do
    total = sampler%cdf(n)
    if (total <= 0.0_wp) error stop 'energy sampler: zero integrated luminosity'
    sampler%total_integral = total
    sampler%cdf = sampler%cdf / total
    sampler%cdf(n) = 1.0_wp
  end subroutine build_tabulated_sampler

  subroutine build_tabulated_band_sampler(sampler, energy, density, emin, emax)
    type(energy_sampler_type), intent(out) :: sampler
    real(wp), intent(in) :: energy(:), density(:), emin, emax
    real(wp), allocatable :: clipped_e(:), clipped_f(:)
    real(wp) :: lo, hi
    integer :: i, n, nc

    n = size(energy)
    if (n < 2 .or. size(density) /= n) &
       error stop 'energy sampler: table needs at least two matching columns'
    if (any(density < 0.0_wp)) &
       error stop 'energy sampler: negative spectral density'
    do i = 1, n - 1
       if (energy(i+1) <= energy(i)) &
          error stop 'energy sampler: energies must be strictly ascending'
    end do
    if (emax <= emin) error stop 'energy sampler: invalid requested band'
    lo = max(emin, energy(1))
    hi = min(emax, energy(n))
    if (hi <= lo) error stop 'energy sampler: table does not overlap requested band'

    nc = 2 + count(energy > lo .and. energy < hi)
    allocate(clipped_e(nc), clipped_f(nc))
    clipped_e(1) = lo
    clipped_f(1) = linear_table_value(energy, density, lo)
    nc = 1
    do i = 1, n
       if (energy(i) > lo .and. energy(i) < hi) then
          nc = nc + 1
          clipped_e(nc) = energy(i)
          clipped_f(nc) = density(i)
       end if
    end do
    nc = nc + 1
    clipped_e(nc) = hi
    clipped_f(nc) = linear_table_value(energy, density, hi)
    call build_tabulated_sampler(sampler, clipped_e(1:nc), clipped_f(1:nc))
    deallocate(clipped_e, clipped_f)
  end subroutine build_tabulated_band_sampler

  subroutine build_planck_sampler(sampler, temperature, emin, emax, rtol, thresholds)
    type(energy_sampler_type), intent(out) :: sampler
    real(wp), intent(in) :: temperature, emin, emax
    real(wp), intent(in), optional :: rtol
    real(wp), intent(in), optional :: thresholds(:)
    real(wp), allocatable :: mandatory(:), work_e(:), work_f(:)
    real(wp) :: tol, scale, span
    integer :: i, j, nmandatory, nwork

    if (temperature <= 0.0_wp .or. emin < 0.0_wp .or. emax <= emin) &
       error stop 'energy sampler: invalid Planck temperature or energy band'
    tol = 1.0e-8_wp
    if (present(rtol)) tol = rtol
    if (tol <= 0.0_wp) error stop 'energy sampler: Planck rtol must be positive'

    nmandatory = 2
    if (present(thresholds)) then
       do i = 1, size(thresholds)
          if (thresholds(i) > emin .and. thresholds(i) < emax) &
             nmandatory = nmandatory + 1
       end do
    end if
    allocate(mandatory(nmandatory))
    mandatory(1) = emin
    j = 1
    if (present(thresholds)) then
       do i = 1, size(thresholds)
          if (thresholds(i) > emin .and. thresholds(i) < emax) then
             j = j + 1
             mandatory(j) = thresholds(i)
          end if
       end do
    end if
    mandatory(nmandatory) = emax
    call sort_unique(mandatory, nmandatory)

    allocate(work_e(MAX_PLANCK_KNOTS), work_f(MAX_PLANCK_KNOTS))
    scale = rough_planck_integral(temperature, emin, emax)
    span = emax - emin
    nwork = 1
    work_e(1) = mandatory(1)
    work_f(1) = planck_energy_shape(mandatory(1), temperature)
    do i = 1, nmandatory - 1
       call refine_planck_interval(temperature, mandatory(i), mandatory(i+1), &
            planck_energy_shape(mandatory(i), temperature), &
            planck_energy_shape(mandatory(i+1), temperature), tol, scale, span, &
            0, work_e, work_f, nwork)
    end do

    call build_tabulated_sampler(sampler, work_e(1:nwork), work_f(1:nwork))
    deallocate(mandatory, work_e, work_f)
  end subroutine build_planck_sampler

  recursive subroutine refine_planck_interval(temperature, e0, e1, f0, f1, &
       tol, scale, span, depth, work_e, work_f, nwork)
    real(wp), intent(in) :: temperature, e0, e1, f0, f1, tol, scale, span
    integer, intent(in) :: depth
    real(wp), intent(inout) :: work_e(:), work_f(:)
    integer, intent(inout) :: nwork
    real(wp) :: em, fm, coarse, fine, err, allowance

    em = 0.5_wp * (e0 + e1)
    fm = planck_energy_shape(em, temperature)
    coarse = 0.5_wp * (f0 + f1) * (e1 - e0)
    fine = 0.25_wp * (f0 + 2.0_wp*fm + f1) * (e1 - e0)
    err = abs(fine - coarse) / 3.0_wp
    allowance = tol * (abs(fine) + scale * (e1 - e0) / span)

    if (err <= allowance .or. depth >= 50) then
       nwork = nwork + 1
       if (nwork > size(work_e)) error stop 'energy sampler: Planck knot limit'
       work_e(nwork) = e1
       work_f(nwork) = f1
    else
       call refine_planck_interval(temperature, e0, em, f0, fm, tol, scale, &
            span, depth+1, work_e, work_f, nwork)
       call refine_planck_interval(temperature, em, e1, fm, f1, tol, scale, &
            span, depth+1, work_e, work_f, nwork)
    end if
  end subroutine refine_planck_interval

  real(wp) function sample_energy_cdf(sampler, u) result(energy)
    type(energy_sampler_type), intent(in) :: sampler
    real(wp), intent(in) :: u
    integer :: lo, hi, mid
    real(wp) :: target, interval_area, local_area, dx, y0, y1, slope
    real(wp) :: disc, offset, slope_scale

    if (.not. allocated(sampler%cdf)) error stop 'energy sampler: uninitialized CDF'
    target = min(max(u, 0.0_wp), 1.0_wp)
    if (target <= 0.0_wp) then
       energy = sampler%energy(1)
       return
    else if (target >= 1.0_wp) then
       energy = sampler%energy(size(sampler%energy))
       return
    end if

    lo = 1
    hi = size(sampler%cdf)
    do while (hi - lo > 1)
       mid = (lo + hi) / 2
       if (sampler%cdf(mid) <= target) then
          lo = mid
       else
          hi = mid
       end if
    end do

    interval_area = sampler%cdf(hi) - sampler%cdf(lo)
    if (interval_area <= 0.0_wp) then
       energy = sampler%energy(hi)
       return
    end if
    dx = sampler%energy(hi) - sampler%energy(lo)
    y0 = sampler%density(lo)
    y1 = sampler%density(hi)
    local_area = (target - sampler%cdf(lo)) / interval_area * &
                 (0.5_wp * (y0 + y1) * dx)
    slope = (y1 - y0) / dx
    slope_scale = max(abs(y0), abs(y1), tiny(1.0_wp))
    if (abs(slope*dx) <= 64.0_wp*epsilon(1.0_wp)*slope_scale) then
       if (y0 > 0.0_wp) then
          offset = local_area / y0
       else
          offset = 0.0_wp
       end if
    else
       disc = sqrt(max(0.0_wp, y0*y0 + 2.0_wp*slope*local_area))
       offset = 2.0_wp * local_area / (y0 + disc)
    end if
    energy = sampler%energy(lo) + min(max(offset, 0.0_wp), dx)
  end function sample_energy_cdf

  real(wp) function energy_cdf_value(sampler, energy) result(value)
    type(energy_sampler_type), intent(in) :: sampler
    real(wp), intent(in) :: energy
    integer :: lo, hi, mid
    real(wp) :: dx, offset, y0, slope, interval_total, partial

    if (energy <= sampler%energy(1)) then
       value = 0.0_wp
       return
    else if (energy >= sampler%energy(size(sampler%energy))) then
       value = 1.0_wp
       return
    end if
    lo = 1
    hi = size(sampler%energy)
    do while (hi - lo > 1)
       mid = (lo + hi) / 2
       if (sampler%energy(mid) <= energy) then
          lo = mid
       else
          hi = mid
       end if
    end do
    dx = sampler%energy(hi) - sampler%energy(lo)
    offset = energy - sampler%energy(lo)
    y0 = sampler%density(lo)
    slope = (sampler%density(hi) - y0) / dx
    interval_total = 0.5_wp * (y0 + sampler%density(hi)) * dx
    partial = y0*offset + 0.5_wp*slope*offset*offset
    if (interval_total > 0.0_wp) then
       value = sampler%cdf(lo) + (sampler%cdf(hi) - sampler%cdf(lo)) * &
               partial / interval_total
    else
       value = sampler%cdf(lo)
    end if
  end function energy_cdf_value

  real(wp) function planck_energy_shape(energy, temperature) result(shape)
    real(wp), intent(in) :: energy, temperature
    real(wp) :: x
    x = energy / (KB_EV * temperature)
    if (x <= 0.0_wp) then
       shape = 0.0_wp
    else if (x < 1.0e-3_wp) then
       shape = energy**3 / (x * (1.0_wp + 0.5_wp*x + x*x/6.0_wp))
    else if (x < 700.0_wp) then
       shape = energy**3 / (exp(x) - 1.0_wp)
    else
       shape = energy**3 * exp(-x)
    end if
  end function planck_energy_shape

  real(wp) function sample_full_planck_bc(temperature, u) result(energy)
    real(wp), intent(in) :: temperature
    real(wp), intent(in) :: u(5)
    real(wp) :: target, cumulative, product
    integer :: ell

    if (temperature <= 0.0_wp) error stop 'energy sampler: invalid Planck temperature'
    target = min(max(u(1), 0.0_wp), 1.0_wp-epsilon(1.0_wp)) * ZETA4
    cumulative = 0.0_wp
    ell = 0
    do
       ell = ell + 1
       cumulative = cumulative + 1.0_wp / real(ell,wp)**4
       if (cumulative >= target) exit
       if (ell >= 1000000) error stop 'energy sampler: Planck mixture did not converge'
    end do
    product = product_clamped(u(2:5))
    energy = KB_EV * temperature * (-log(product)) / real(ell,wp)
  end function sample_full_planck_bc

  real(wp) function product_clamped(u) result(product)
    real(wp), intent(in) :: u(:)
    integer :: i
    product = 1.0_wp
    do i = 1, size(u)
       product = product * min(max(u(i), tiny(1.0_wp)), 1.0_wp)
    end do
  end function product_clamped

  real(wp) function rough_planck_integral(temperature, emin, emax) result(total)
    real(wp), intent(in) :: temperature, emin, emax
    integer, parameter :: n = 256
    integer :: i
    real(wp) :: e0, e1
    total = 0.0_wp
    do i = 1, n
       e0 = emin + (emax-emin)*real(i-1,wp)/real(n,wp)
       e1 = emin + (emax-emin)*real(i,wp)/real(n,wp)
       total = total + 0.5_wp * (planck_energy_shape(e0,temperature) + &
                                 planck_energy_shape(e1,temperature)) * (e1-e0)
    end do
    if (total <= 0.0_wp) error stop 'energy sampler: zero Planck band luminosity'
  end function rough_planck_integral

  real(wp) function linear_table_value(energy, density, value) result(result)
    real(wp), intent(in) :: energy(:), density(:), value
    integer :: lo, hi, mid
    real(wp) :: fraction
    if (value <= energy(1)) then
       result = density(1)
       return
    else if (value >= energy(size(energy))) then
       result = density(size(density))
       return
    end if
    lo = 1
    hi = size(energy)
    do while (hi-lo > 1)
       mid = (lo+hi)/2
       if (energy(mid) <= value) then
          lo = mid
       else
          hi = mid
       end if
    end do
    fraction = (value-energy(lo))/(energy(hi)-energy(lo))
    result = density(lo) + fraction*(density(hi)-density(lo))
  end function linear_table_value

  subroutine sort_unique(values, n)
    real(wp), allocatable, intent(inout) :: values(:)
    integer, intent(inout) :: n
    real(wp), allocatable :: tmp(:)
    real(wp) :: hold
    integer :: i, j, nkeep

    do i = 2, n
       hold = values(i)
       j = i - 1
       do while (j >= 1)
          if (values(j) <= hold) exit
          values(j+1) = values(j)
          j = j - 1
       end do
       values(j+1) = hold
    end do
    nkeep = 1
    do i = 2, n
       if (values(i) > values(nkeep)) then
          nkeep = nkeep + 1
          values(nkeep) = values(i)
       end if
    end do
    if (nkeep < size(values)) then
       allocate(tmp(nkeep))
       tmp = values(1:nkeep)
       call move_alloc(tmp, values)
    end if
    n = nkeep
  end subroutine sort_unique

end module energy_sampler_mod
