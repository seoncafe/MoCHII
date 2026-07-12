! MoCHII: copied from MoCafe_v2.00/src/define.f90 (2026-07-11)
module define
  use, intrinsic :: iso_fortran_env, only: real32, real64, real128, int32, int64
public

!! kind parameter for single and double precision (6 and 15 significant decimal digits)
!! precision of 32-, 64-, and 128-bit reals
!  integer, parameter :: sp = selected_real_kind(6, 37)
!  integer, parameter :: dp = selected_real_kind(15, 307)
!  integer, parameter :: qp = selected_real_kind(33, 4931)
!! precision of 32-, 64- bit integers
!  integer, parameter :: i4b = selected_int_kind(r=9)
!  integer, parameter :: i8b = selected_int_kind(r=18)
  integer, parameter :: sp  = real32
  integer, parameter :: dp  = real64
  integer, parameter :: qp  = real128
  integer, parameter :: i4b = int32
  integer, parameter :: i8b = int64

! working precision
! Note that (1) iwp (integer working precision) applies only to seed of KISS random number generator.
!           (2) nphotons is always 64-bit integer.
!           (3) other integer variables are always 32-bit integers.
!  integer, parameter :: wp  = sp
!  integer, parameter :: iwp = i4b
  integer, parameter :: wp  = dp
  integer, parameter :: iwp = i8b

! numerical values
  real(kind=wp), parameter :: pi      = 3.141592653589793238462643383279502884197_wp
  real(kind=wp), parameter :: twopi   = 6.283185307179586476925286766559005768394_wp
  real(kind=wp), parameter :: fourpi  = 12.56637061435917295385057353311801153679_wp
  real(kind=wp), parameter :: halfpi  = pi/2.0_wp
  real(kind=wp), parameter :: sqrttwo = sqrt(2.0_wp)

! conversion factor (radian to degree, degree to radian)
  real(kind=wp), parameter :: rad2deg = 180.0_wp/pi
  real(kind=wp), parameter :: deg2rad = pi/180.0_wp

! speed of light in km/s
  real(kind=wp), parameter :: speedc  = 2.99792458e5_wp

! h_planck = Planck's constant in m^2 kg/s
! massH    = Hydrogen mass in kg
  real(kind=wp), parameter :: h_planck = 6.62607004e-34_wp
  real(kind=wp), parameter :: massH    = 1.6737236e-27_wp
  real(kind=wp), parameter :: g0       = h_planck/massH

! distances
  real(kind=wp), parameter :: pc2cm  = 3.0856776e18_wp
  real(kind=wp), parameter :: kpc2cm = pc2cm * 1e3_wp
  real(kind=wp), parameter :: au2cm  = 1.4960e13_wp
  real(kind=wp), parameter :: ang2m  = 1.0e-10_wp
  real(kind=wp), parameter :: ang2km = 1.0e-13_wp
  real(kind=wp), parameter :: um2m   = 1.0e-6_wp
  real(kind=wp), parameter :: um2km  = 1.0e-9_wp

! MoCHII: cgs constants for the ionizing band and gas microphysics
  real(kind=wp), parameter :: h_planck_cgs = 6.62607015e-27_wp   ! erg s
  real(kind=wp), parameter :: clight_cgs   = 2.99792458e10_wp    ! cm/s
  real(kind=wp), parameter :: ev2erg       = 1.602176634e-12_wp  ! erg/eV
  real(kind=wp), parameter :: kboltz_cgs   = 1.380649e-16_wp     ! erg/K
  real(kind=wp), parameter :: eth_HI       = 13.598_wp           ! eV, H I ionization threshold
  real(kind=wp), parameter :: eth_HeI      = 24.587_wp           ! eV, He I
  real(kind=wp), parameter :: eth_HeII     = 54.416_wp           ! eV, He II

! dust mass per hydrogen nuclei (g/H)
  real(kind=wp), parameter :: h2dust = 1.87e-26_wp

! kappa_v : dust extinction cross-section (cm^2/g) at V-band
  real(kind=wp), parameter :: kappa_v = 25992.6_wp

! maximum number of observers
  integer, parameter :: MAX_OBSERVERS = 181

! maximum number of albedo/asymmetry-factor values in a single-run (a,g) scan
  integer, parameter :: MAX_SCAN = 32

! maximum number of stellar source components (multi-population SED, Stage 6)
  integer, parameter :: MAX_SRC = 16

! tinest = the smallest positive number
! eps    = the least positive number
!          that added to 1 returns a number that is greater than 1
! hugest = the hugest positive number
  real(kind=wp), parameter :: tinest = tiny(0.0_wp)
  real(kind=wp), parameter :: eps    = epsilon(0.0_wp)
  real(kind=sp), parameter :: eps_sp = epsilon(0.0_sp)
  real(kind=dp), parameter :: eps_dp = epsilon(0.0_dp)
  real(kind=wp), parameter :: hugest = huge(1.0_wp)

! define NaN
  real(real64),  parameter :: nan64 = transfer(-2251799813685248_int64, 1._real64)
  real(real32),  parameter :: nan32 = transfer(-4194304_int32, 1._real32)

! photon type
  type photon_type
     real(kind=wp) :: x,y,z
     real(kind=wp) :: kx,ky,kz
     real(kind=wp) :: nx,ny,nz
     real(kind=wp) :: mx,my,mz
     integer       :: icell,jcell,kcell
     integer       :: icell_clump = 0    ! current/owner clump index (0 = vacuum; clump grid only)
     integer       :: icell_amr   = 0    ! current leaf index (0 = unknown/outside; amr grid only)
     integer       :: nscatt
     real(kind=wp) :: wgt
     logical       :: inside
     real(kind=wp) :: albedo, hgg
     !--- SED (multi-wavelength) mode: wavelength bin index, wavelength [um],
     !--- and the grey-rescaling factor s_ext = C_ext(lambda)/C_ext(lambda_ref).
     !--- s_ext = 1 in monochromatic mode (keeps the mono path bit-identical).
     integer       :: il     = 0
     real(kind=wp) :: lambda = 0.0_wp
     real(kind=wp) :: s_ext  = 1.0_wp
     !--- MoCHII: ionizing-band frequency bin index (0 = not an ionizing
     !--- photon).  The ionizing band carries its own opacity array
     !--- kap_ion(inu,leaf); the grey s_ext rescale does not apply there.
     integer       :: inu    = 0
     !--- physical luminosity carried by the packet [erg/s], used by the
     !--- energy (J_lambda) tally so stellar and dust-emission photons with
     !--- different packet energies combine correctly (1.0 in mono mode).
     real(kind=wp) :: Lpacket = 1.0_wp
     ! Stokes parameters
     real(kind=wp) :: I,Q,U,V
  end type photon_type

! cylindrical grid type
! nr, np, nz          - number of cells
! rface, pface, zface - locations of cell faces
! opacity             - opacity per unit length
  type grid_type
     integer :: nx,ny,nz
     integer :: i0 = 0, j0 = 0, k0 = 0
     real(kind=wp) :: rmax
     real(kind=wp) :: xmin,xmax,xrange
     real(kind=wp) :: ymin,ymax,yrange
     real(kind=wp) :: zmin,zmax,zrange
     real(kind=wp) :: dx,dy,dz
     real(kind=wp), pointer :: xface(:)      => null()
     real(kind=wp), pointer :: yface(:)      => null()
     real(kind=wp), pointer :: zface(:)      => null()
     real(kind=wp), pointer :: rhokap(:,:,:) => null()
  end type grid_type

! input parameters
  type params_type
     integer(kind=int64) :: nphotons   = 1e5
     integer             :: nprint     = 1e7
     real(kind=wp)       :: no_photons = 1e5
     real(kind=wp)       :: no_print   = 1e7
     integer       :: iseed        = 0
     real(kind=wp) :: luminosity   = 1.0
     real(kind=wp) :: tauhomo      = -999.0
     real(kind=wp) :: taumax       = -999.0
     real(kind=wp) :: lambda0      = 6563.0_wp
     !--- grid geometry parameters
     character(len=128) :: geometry = ''
     logical       :: xyz_symmetry    = .false.
     logical       :: xy_periodic     = .false.
     logical       :: z_symmetry      = .false.
     logical       :: sightline_tau   = .false.
     integer       :: nx   = 1
     integer       :: ny   = 1
     integer       :: nz   = 11
     real(kind=wp) :: xmax = 1.0_wp
     real(kind=wp) :: ymax = 1.0_wp
     real(kind=wp) :: zmax = 1.0_wp
     !--- to simulate spherical geometry
     integer       :: nr   = -999
     real(kind=wp) :: rmax = -999.0
     real(kind=wp) :: rmin = 0.0
     !real(kind=wp) :: rsource = 0.0_wp
     !--- continuum parameters
     real(kind=wp) :: source_zscale          = nan64
     !--- radial scale length of an 'exponential' stellar disk source; when set,
     !--- x,y are drawn radially exponential (r ~ r*exp(-r/source_rscale))
     !--- instead of uniform over the plane.
     real(kind=wp) :: source_rscale          = nan64
     real(kind=wp) :: density_zscale         = nan64
     real(kind=wp) :: density_rscale         = nan64
     real(kind=wp) :: density_powerlaw_index = nan64
     !--- disk vertical profile ('exp' = exp(-|z|/zscale), 'sech' = sech^2(z/zscale))
     !--- and log-spiral arm modulation (1 + spiral_amp*sin(m*(ln r/tan(pitch)-phi)))
     !--- applied to the cylinder dust density and to 'exp_spiral' disk sources.
     !--- spiral_m = 0 disables the arms.
     character(len=8)  :: density_zprofile = 'exp'
     integer           :: spiral_m         = 0
     real(kind=wp)     :: spiral_pitch     = 15.0_wp    ! pitch angle [deg]
     real(kind=wp)     :: spiral_amp       = 0.5_wp     ! arm amplitude
     real(kind=wp) :: distance2cm          = kpc2cm
     character(len=128) :: distance_unit   = 'kpc'
     character(len=128) :: source_geometry = 'point'
     character(len=128) :: base_name       = ''
     character(len=128) :: out_file        = ''
     integer            :: out_bitpix      = -32
     logical            :: out_bitpix_force = .false.   ! .true. = keep out_bitpix, skip the float32->float64 auto-promote
     character(len=8)   :: file_format     = 'hdf5'
     !--- location of a point source.
     real(kind=wp) :: xs_point = 0.0_wp
     real(kind=wp) :: ys_point = 0.0_wp
     real(kind=wp) :: zs_point = 0.0_wp
     !--- shared memory & master-slave algorithm
     integer       :: num_send_at_once   = 10000
     logical       :: use_shared_memory  = .false.
     logical       :: use_master_slave   = .true.
     logical       :: use_reduced_wgt    = .true.
     !--- Dust-related parameters
     real(kind=wp) :: hgg           = 0.6761
     real(kind=wp) :: albedo        = 0.3253
     real(kind=wp) :: cext_dust     = 1.6059e-21
     logical       :: use_stokes    = .true.
     character(len=128) :: scatt_mat_file = ''
     !--- (albedo, asymmetry-factor) scan: one run, many (a,g).  See Seon 2010, PKAS, 25, 177.
     !--- When use_ag_list = .true., par%hgg is the simulated g0; the scattered image
     !--- becomes a 4-D array scatt(x,y,albedo_list,hgg_list).  Empty (NaN) lists are
     !--- auto-filled with the canonical paper grids (a=0.1..1.0, g=0.0..0.9).
     !--- SED (multi-wavelength) mode (MoCafe v2.00, Stage 1).  When
     !--- use_sed = .true., a single run transports photons over a log-spaced
     !--- wavelength grid [lambda_min, lambda_max] (um) with nlambda bins.
     !--- Dust properties C_ext/albedo/g vs lambda are read from kext_file
     !--- (e.g. SEDust calc_kext_astrodust.x output: lambda, albedo, <cos>,
     !--- C_ext/H columns; '#' comments).  The source spectrum is either a
     !--- 2-column file source_spectrum (lambda[um], L_lambda[arb]) or a
     !--- Planck function with temperature tstar (K).  The grid opacity and
     !--- par%taumax/tauhomo refer to the reference wavelength lambda_ref.
     logical            :: use_sed         = .false.
     integer            :: nlambda         = 128
     real(kind=wp)      :: lambda_min      = 0.0912_wp
     real(kind=wp)      :: lambda_max      = 2000.0_wp
     real(kind=wp)      :: lambda_ref      = 0.55_wp
     character(len=128) :: kext_file       = ''
     character(len=128) :: source_spectrum = ''
     real(kind=wp)      :: tstar           = -999.0_wp
     !--- Stage 2: mean-intensity tally J_lambda(x,y,z) in each cell (Lucy 1999
     !--- pathlength estimator); writes '<base>_jlam.<ext>'.  Requires
     !--- use_sed and the plain Cartesian grid.
     logical            :: save_jlam       = .false.
     !--- Stage 3: dust thermal emission via the SEDust library (Mode 1, Lucy
     !--- 1999; non-iterative first version).  Requires use_sed + save_jlam.
     !--- par%luminosity must be the physical stellar luminosity in erg/s for
     !--- absolute dust temperatures.  Writes '<base>_dustsed.<ext>'.
     logical            :: use_dustemis     = .false.
     !--- dust emission method: 'lucy' = Lucy 1999 + SEDust (equilibrium +
     !--- stochastic + PAH; iterable); 'bw01' = Bjorkman & Wood 2001 immediate
     !--- reemission (approximate, equilibrium mixture-mean opacity, no PAH).
     character(len=8)   :: dust_emission_method = 'lucy'
     character(len=16)  :: dust_model_sed   = 'astrodust'
     !--- SEDust optics/size-distribution paths, resolved relative to
     !--- par%sed_workdir (build_astrodust is called from there so that
     !--- SEDust's own '../data/dielectric/...' relative reads resolve too).
     character(len=256) :: sed_qtable       = '../tmatrix/output/q_astrodust_P0.20_Fe0.00_1.400.dat'
     character(len=256) :: sed_sizedist     = '../data/release/size_distribution.dat'
     !--- SEDust sed/ directory (self-contained under the MoCafe tree).
     !--- Blank (default) = auto-resolve to <executable dir>/SEDust/sed at
     !--- run time (see read_input), so a fresh checkout works from any path.
     !--- Set an explicit absolute path here only to override that.
     character(len=256) :: sed_workdir      = ''
     integer            :: sed_NT           = 200
     real(kind=wp)      :: sed_Tlo          = 2.7_wp
     real(kind=wp)      :: sed_Thi          = 5.0e3_wp
     !--- DL07 (Draine & Li 2007; dust_model_sed='dl07') reuses sed_qtable and
     !--- sed_sizedist (same files as astrodust).  sed_dl07_sdindex = WD01
     !--- size-distribution index (7 = MW R_V=3.1, b_C=6e-5); sed_dl07_uisrf =
     !--- reference radiation-field scaling U (1 = MMP83 diffuse ISM).
     integer            :: sed_dl07_sdindex = 7
     real(kind=wp)      :: sed_dl07_uisrf   = 1.0_wp
     !--- Zubko (ZDA 2004 BARE-GR-S; dust_model_sed='zubko'): the ZDA config
     !--- file and the DustEM optics/calorimetry directory, resolved from
     !--- sed_workdir (defaults point at the copied SEDust/data/zubko).
     character(len=256) :: sed_zubko_config = '../data/zubko/ZDA_BARE_GR_S_Config.dat'
     character(len=256) :: sed_zubko_dir    = '../data/zubko/'
     !--- Lucy iteration for dust self-absorption (Stage 3 follow-up #1).
     !--- dust_niter = max iterations (1 = non-iterative); dust_nphotons =
     !--- dust-emission photons per iteration; dust_tol = relative-change
     !--- convergence on the total emitted (=absorbed) luminosity.
     integer            :: dust_niter       = 1
     real(kind=wp)      :: dust_no_photons  = 1.0e6_wp
     real(kind=wp)      :: dust_tol         = 1.0e-3_wp
     !--- fast emission-table path (Stage 3 follow-up #3): interpolate the
     !--- emission spectrum in the field-intensity scaling U over a fixed
     !--- reference shape, instead of an exact SEDust solve per cell.  Much
     !--- faster on large grids; approximate where the local spectral shape
     !--- departs from the reference (hardened fields near sources).
     logical            :: dust_fast_table  = .false.
     integer            :: dust_nU          = 40
     !--- single-temperature (equilibrium-only) SEDust solve: one mixture
     !--- equilibrium T per cell (dust_emission_single_teq).  Faster than the
     !--- full stochastic solve and gives the true Teq, but drops PAH /
     !--- stochastic-heating features.  More accurate than the B&W mixture mean.
     logical            :: dust_single_teq  = .false.
     !--- multiple stellar source components (Stage 6).  When par%nsource > 1,
     !--- each component i has its own spectrum (Planck src_tstar(i) or file
     !--- src_spectrum(i)), luminosity src_lum(i) [erg/s], geometry
     !--- src_geometry(i) ('point'|'uniform'|'gaussian'|'exponential') and
     !--- geometry parameters (src_x/y/z for 'point', src_zscale/src_rscale
     !--- for disks).  Sources are sampled in proportion to src_lum; the total
     !--- luminosity par%luminosity is set to sum(src_lum).
     integer            :: nsource                 = 1
     character(len=16)  :: src_geometry(MAX_SRC)   = 'point'
     real(kind=wp)      :: src_tstar(MAX_SRC)      = -999.0_wp
     character(len=128) :: src_spectrum(MAX_SRC)   = ''
     real(kind=wp)      :: src_lum(MAX_SRC)        = -999.0_wp
     real(kind=wp)      :: src_x(MAX_SRC)          = 0.0_wp
     real(kind=wp)      :: src_y(MAX_SRC)          = 0.0_wp
     real(kind=wp)      :: src_z(MAX_SRC)          = 0.0_wp
     real(kind=wp)      :: src_zscale(MAX_SRC)     = -999.0_wp
     real(kind=wp)      :: src_rscale(MAX_SRC)     = -999.0_wp
     !--- Sersic bulge (src_geometry = 'sersic'): index n, effective radius Re,
     !--- and vertical axial ratio (oblate flattening z = r*cos*axratio).
     real(kind=wp)      :: src_sersic_index(MAX_SRC) = 4.0_wp
     real(kind=wp)      :: src_reff(MAX_SRC)         = -999.0_wp
     real(kind=wp)      :: src_axial_ratio(MAX_SRC)  = 1.0_wp
     !--- boxiness exponent for 'boxy'/'bar'/'xbar' bulges (2 = ellipsoidal,
     !--- >2 = boxy, <2 = disky); src_reff sets the scale, src_sersic_index is
     !--- the 'boxy' power-law index, src_axial_ratio flattens z.
     real(kind=wp)      :: src_boxiness(MAX_SRC)     = 2.0_wp
     !--- HEALPix all-sky interior observer (Milky-Way case).  When
     !--- par%allsky = .true. the observer sits at par%allsky_x/y/z INSIDE the
     !--- grid and each event peels toward it, binned into the HEALPix pixel
     !--- (RING scheme, npix = 12*allsky_nside^2) of the incoming sky direction.
     logical            :: allsky        = .false.
     integer            :: allsky_nside  = 64       ! HEALPix nside (power of 2)
     real(kind=wp)      :: allsky_x      = 0.0_wp   ! interior observer position (code units)
     real(kind=wp)      :: allsky_y      = 0.0_wp
     real(kind=wp)      :: allsky_z      = 0.0_wp
     !--- Modified Random Walk (Stage 8) for very optically thick cells; used
     !--- in the Lucy energy passes.  mrw_gamma is the trigger threshold on
     !--- R0*rho*kappa (nearest-wall optical depth); ~2 is typical.
     logical            :: use_mrw       = .false.
     real(kind=wp)      :: mrw_gamma     = 2.0_wp
     logical       :: use_ag_list   = .false.
     real(kind=wp) :: albedo_list(MAX_SCAN) = nan64
     real(kind=wp) :: hgg_list(MAX_SCAN)    = nan64
     !--- polychromatic optical-depth (tau) scan: one run, many tau.  See Jonsson 2006,
     !--- MNRAS, 372, 2 (eqs. 29-32).  When use_tau_list = .true., par%taumax is the
     !--- simulated reference tau0 and the scattered image gains a tau axis
     !--- scatt(x,y,albedo,hgg,tau) while the direct image becomes direc(x,y,tau).
     !--- Composes with use_ag_list (the a/g axes collapse to length 1 when it is off).
     !--- tau_list holds target taumax values; the reference par%taumax is auto-inserted
     !--- if absent.  An empty (NaN) list is auto-filled with sqrt(2)-spaced
     !--- tau/taumax in [0.5, 2].
     logical       :: use_tau_list  = .false.
     real(kind=wp) :: tau_list(MAX_SCAN)     = nan64
     !--- grid-type selector (modeled on LaRT v2.10).  'car' = Cartesian (default);
     !--- 'clump' = clumpy spherical dust medium (Part B); 'amr' reserved (Part A).
     !--- Legacy booleans use_clump_medium / use_amr_grid are mapped to grid_type
     !--- in read_input for backward compatibility.
     character(len=8) :: grid_type        = 'car'
     logical          :: use_clump_medium = .false.
     logical          :: use_amr_grid     = .false.
     !--- AMR octree grid (dust only; velocity/temperature deferred).  See
     !--- AMR_CLUMPS_PLAN.md Part A.  Leaf data are read from a generic AMR
     !--- file (FITS/HDF5/text) produced by the Python builders/converters;
     !--- amr_type='ramses' is rejected (convert to 'generic' first).  The
     !--- dust opacity of each leaf follows dust_model: 'global_dgr'
     !--- (nH*cext_dust*DGR), 'from_file' (ndust column), or 'laursen09'
     !--- ((Z/Z_ref)*(nHI+f_ion*nHII), requires an xHI column).
     character(len=32)  :: amr_type   = 'generic'
     character(len=128) :: amr_file   = ''
     character(len=32)  :: dust_model = 'global_dgr'
     real(kind=wp)      :: DGR        = 1.0e-2_wp
     real(kind=wp)      :: Z_global   = 0.0134_wp
     real(kind=wp)      :: Z_ref      = 0.0134_wp
     real(kind=wp)      :: f_ion_dust = 0.01_wp
     !--- clumpy-medium parameters (dust only; velocity/temperature deferred).
     !--- See AMR_CLUMPS_PLAN.md Part B and docs/MoCafe_clump.tex.  A sphere of
     !--- radius par%rmax holds N non-overlapping dust clumps of radius
     !--- clump_radius; the inter-clump medium is vacuum (par%rmin carves an
     !--- inner cavity).  Clump count: one of clump_f_cov / clump_f_vol /
     !--- clump_N_clumps.  Dust opacity of each clump kappa_c: one of clump_tau0
     !--- (dimensionless, center->surface) / clump_ndust / clump_nH; otherwise
     !--- back-solved from the system target par%taumax or par%tauhomo.
     real(kind=wp) :: clump_radius   = -1.0_wp
     real(kind=wp) :: clump_f_cov    = -1.0_wp
     real(kind=wp) :: clump_f_vol    = -1.0_wp
     integer       :: clump_N_clumps = -1
     real(kind=wp) :: clump_tau0     = -1.0_wp
     real(kind=wp) :: clump_ndust    = -1.0_wp
     real(kind=wp) :: clump_nH       = -1.0_wp
     logical       :: clump_fully_inside = .true.
     !--- optional spatial profiles (geometry only): constant|powerlaw|gaussian|exponential|file
     character(len=16) :: clump_radius_profile  = 'constant'
     character(len=16) :: clump_density_profile = 'constant'
     character(len=16) :: clump_number_profile  = 'constant'
     real(kind=wp) :: clump_radius_alpha  = 0.0_wp
     real(kind=wp) :: clump_radius_r0     = -1.0_wp
     real(kind=wp) :: clump_density_alpha = 0.0_wp
     real(kind=wp) :: clump_density_r0    = -1.0_wp
     real(kind=wp) :: clump_number_alpha  = 0.0_wp
     real(kind=wp) :: clump_number_r0     = -1.0_wp
     character(len=128) :: clump_profile_file = ''
     real(kind=wp) :: cone_opening    = 0.0_wp
     character(len=128) :: clump_input_file = ''
     logical       :: save_clump_info = .false.
     !--- output-related parameters
     logical       :: save_direc0   = .false.
     character(len=10)  :: output_normalization = 'luminosity'
     !--- distribution function for external illumination
     character(len=128) :: radiation_angular_PDF_file = ''
     !--- density file
     character(len=128) :: density_file  = ''
     integer            :: reduce_factor = 1
     integer            :: centering     = 0
     !--- parameters for PEELING-OFF
     integer :: nobs = 1
     integer :: nxim = 129
     integer :: nyim = 129
     real(kind=wp) :: distance          = nan64
     real(kind=wp) :: inclination_angle(MAX_OBSERVERS) = nan64
     real(kind=wp) :: position_angle(MAX_OBSERVERS)    = nan64
     real(kind=wp) :: phase_angle(MAX_OBSERVERS)       = nan64
     real(kind=wp) :: alpha(MAX_OBSERVERS) = nan64
     real(kind=wp) :: beta(MAX_OBSERVERS)  = nan64
     real(kind=wp) :: gamma(MAX_OBSERVERS) = nan64
     real(kind=wp) :: obsx(MAX_OBSERVERS)  = nan64
     real(kind=wp) :: obsy(MAX_OBSERVERS)  = nan64
     real(kind=wp) :: obsz(MAX_OBSERVERS)  = nan64
     real(kind=wp) :: dxim  = nan64
     real(kind=wp) :: dyim  = nan64
     !--- MoCHII: ionizing frequency band (nu >= 13.6 eV) and gas leaf state.
     !--- Stage G0: fixed neutral fractions, no equilibrium feedback; the band
     !--- transports source packets and tallies J_nu per leaf, from which the
     !--- photoionization rates Gamma_i and photoheating H_i follow.
     logical            :: use_ion_band = .false.
     integer            :: nnu_ion      = 16          ! frequency bins
     real(kind=wp)      :: eion_min     = 13.598_wp   ! band lower edge [eV]
     real(kind=wp)      :: eion_max     = 100.0_wp    ! band upper edge [eV]
     !--- source spectrum in the band: 2-column file (E[eV], L_E[arb]) or,
     !--- when empty, a Planck function with par%tstar [K].
     character(len=128) :: ion_spectrum = ''
     !--- FUV extension of the band: nnu_fuv extra log bins on
     !--- [efuv_min, eion_min] BELOW the nnu_ion ionizing bins, so the
     !--- ionizing bin resolution is unchanged.  par%luminosity keeps its
     !--- meaning as the [eion_min, eion_max] luminosity; the FUV part is
     !--- added on top from the same source spectrum.  H/He cross sections
     !--- vanish below threshold, so the FUV bins are absorbed by dust and
     !--- (with ion_metal_abs) low-threshold metal stages only.
     logical            :: add_fuv      = .false.
     real(kind=wp)      :: efuv_min     = 6.0_wp      ! FUV lower edge [eV]
     integer            :: nnu_fuv      = 8           ! FUV bins
     !--- further registry elements (default 0 = off; pass GAS-PHASE
     !--- values — Si and especially Ca are strongly depleted onto
     !--- grains, Cl is nearly undepleted; typical Orion gas-phase:
     !--- Si 3e-6, Cl 3e-7, Ca 2e-8).  Five stages each; Cl uses the
     !--- VY95 outer-shell photoionization form (absent from VFKY96).
     real(kind=wp)      :: abund_Si     = 0.0_wp
     real(kind=wp)      :: abund_Cl     = 0.0_wp
     real(kind=wp)      :: abund_Ca     = 0.0_wp
     !--- metal photoionization absorption in kap_ion (active with
     !--- use_metals): sum over stages of n_el frac_i sigma_VFKY96,i(E).
     !--- Negligible next to H/He above 13.6 eV; the only GAS opacity in
     !--- the FUV bins (Mg I 7.65, C I 11.26, S I 10.36, Fe I 7.90 eV).
     logical            :: ion_metal_abs = .true.
     !--- PDR-zone physics (each defaults off so the recorded gates
     !--- reproduce; turn all three on together with add_fuv for
     !--- meaningful PDR temperatures):
     !--- metal_ne:    electrons from the metal cascade join the n_e
     !---              closure (the only electrons beyond the I-front);
     !--- metal_heat:  metal photoheating joins the thermal balance
     !---              (heating integrals accumulated with the metal
     !---              Gamma's);
     !--- grain_pe:    grain photoelectric heating + grain recombination
     !---              cooling (Bakes & Tielens 1994 fits) driven by the
     !---              local FUV field in Habing units (needs add_fuv);
     !---              scaled by pe_scale x DGR relative to MW dust.
     logical            :: metal_ne     = .false.
     logical            :: metal_heat   = .false.
     logical            :: grain_pe     = .false.
     real(kind=wp)      :: pe_scale     = 1.0_wp
     !--- He I EXCITED-level recombination radiation as a fourth diffuse
     !--- channel (needs diffuse_field): ~60% of He II recombinations end
     !--- in excited levels whose decay photons ionize H — AGN3 low-density
     !--- branching, 3/4 -> 2^3S (19.82 eV), 1/4 x 2/3 -> 2^1P (584 A,
     !--- 21.22 eV), 1/4 x 1/3 -> 2^1S (two-photon; probability 0.56 of an
     !--- H-ionizing photon, energy sampled flat in [13.6, 20.62] eV).
     logical            :: hei_diffuse  = .false.
     !--- free-free Gaunt factors from the van Hoof et al. (2014) table
     !--- (data/gauntff_vh14.dat) instead of the ported Hummer (1988)
     !--- getGauntFF — the PLAN item-3 swap; off by default so the
     !--- recorded gates reproduce.
     logical            :: gaunt_vh14   = .false.
     !--- peel-off imaging of the band: after convergence, one extra
     !--- transport pass peels direct (stellar + diffuse) and
     !--- dust-scattered contributions toward the observers
     !--- (par%obsx/y/z or par%alpha/beta/gamma, par%nxim/nyim,
     !--- par%dxim/dyim auto-sized) -> '<base>_image'.
     logical            :: ion_peel     = .false.
     !--- bin-resolved peel images: accumulate one image for each band
     !--- bin -> 'direct_cube'/'scatt_cube' (nxim, nyim, nnu_ion) blocks
     !--- (hardness maps, synthetic EUV/FUV photometry); the
     !--- band-integrated channel images are still written (derived from
     !--- the cube).
     logical            :: peel_bins    = .false.
     !--- convergence criterion: 'cell' = cell-max tests (max|dx_HII| <
     !--- gas_tol, max|dTe|/Te < gas_tol_te — stall at the front-cell MC
     !--- noise floor and at no-heating cells); 'vol' = volume-integrated
     !--- tests against the SAME tolerances:
     !---   sum V|dx_HII| / sum V x_HII         < gas_tol
     !---   sum V ne|dTe| / sum V ne Te         < gas_tol_te
     !--- (n_e weight removes vacuum / no-heating cells).  Both measures
     !--- are printed each iteration.
     character(len=8)   :: conv_crit    = 'cell'
     !--- leaf-by-leaf line emissivities -> '<base>_emis' (HDF5/FITS):
     !--- emis_h (SH95 H I lines) and emis_<el>_<stage> blocks
     !--- [erg/s/cm^3; sum(emis*V) = L_line, rows match <base>_lines.txt]
     !--- with wl_* wavelength arrays, plus the state needed for maps
     !--- (LeafXYZ, n_H, n_e, T_e, x_HI/HeI/HeII, x_<el>_stages).
     logical            :: emis_output  = .false.
     !--- G0 gas-state initialization (overridden by AMR file columns where
     !--- present): ion fractions and the He abundance n_He/n_H.
     real(kind=wp)      :: xHI_init     = 1.0_wp
     real(kind=wp)      :: xHeI_init    = 1.0_wp
     real(kind=wp)      :: xHeII_init   = 0.0_wp
     real(kind=wp)      :: He_abund     = 0.1_wp
     !--- add dust absorption into the ionizing band (dusty Stromgren
     !--- competition, PLAN section 7 item 2).  The leaf grey rhokap is
     !--- scaled per bin by (1-albedo) C_ext(E)/C_ext(lambda_ref) read
     !--- from the kext table ion_dust_kext (must cover the EUV, e.g.
     !--- kext_albedo_WD_MW_3.1_60_D03).  EUV dust scattering deferred.
     logical            :: ion_add_dust = .false.
     character(len=256) :: ion_dust_kext = ''
     !--- sample dust scattering in the band (forced first interaction
     !--- after the analytic direct walk; HG phase with g(E) from the
     !--- kext table; albedo weight + Russian roulette).  When off, the
     !--- scattered fraction simply passes (absorption-only treatment).
     logical            :: ion_dust_scatter = .false.
     !--- SEDust stochastic dust emission at output time (needs
     !--- ion_add_dust + sed_workdir/sed_qtable/sed_sizedist paths).
     logical            :: dust_sed     = .false.
     !--- dust-band leaf emissivities at these wavelengths [um] (finite
     !--- entries are active) -> '<base>_dustemis' with emis_dust /
     !--- wl_dust blocks [erg/s/cm^3/um]; the Python map maker projects
     !--- them like line emissivities (IR optically thin, so a column
     !--- map IS the image).
     real(kind=wp)      :: dust_emis_bands(8) = nan64
     !--- x_HII-dependent PAH survival in the EMISSION model: the leaf's
     !--- SEDust spectrum is assembled from the model's output channels
     !--- with the PAH channel weighted by (xHI + f_ion_pah xHII)
     !--- [divided by the dust survival (xHI + f_ion_dust xHII) under
     !--- laursen09_live, whose rhokap already carries it].  Requires a
     !--- model with an explicit 'PAH' channel (astrodust).
     logical            :: sed_pah_live = .false.
     !--- PAH share of the reference extinction and the PAH survival in
     !--- ionized gas (laursen09_live: PLAN section 7 item 3 option).
     real(kind=wp)      :: f_pah        = 0.0_wp
     real(kind=wp)      :: f_ion_pah    = 0.0_wp
     !--- Stage G1: ionization equilibrium + opacity feedback iteration.
     !--- gas_niter = 0 keeps the G0 single-pass behavior (no solve).
     !--- Convergence: max |delta x_HII| < gas_tol.  te_fixed [K] is the
     !--- fixed electron temperature until the G2 thermal balance.
     !--- ion_relax < 1 under-relaxes x between iterations (sharp I-fronts).
     integer            :: gas_niter    = 0
     real(kind=wp)      :: gas_tol      = 1.0e-3_wp
     real(kind=wp)      :: te_fixed     = 1.0e4_wp
     real(kind=wp)      :: ion_relax    = 1.0_wp
     character(len=1)   :: case_ab      = 'B'
     !--- Stage G2: thermal balance.  solve_te replaces the fixed-Te solve
     !--- with heating = cooling bisection in [te_min, te_max]; convergence
     !--- additionally requires max |delta Te|/Te < gas_tol_te.  Tier-1
     !--- cooling coefficient files are read from atomic_dir.
     logical            :: solve_te     = .false.
     real(kind=wp)      :: te_min       = 1.0e3_wp
     real(kind=wp)      :: te_max       = 5.0e4_wp
     real(kind=wp)      :: gas_tol_te   = 1.0e-3_wp
     character(len=256) :: atomic_dir   = 'data/atomic'
     !--- Stage G2b: trace-metal registry.  An element is active when its
     !--- abundance n(X)/n(H) > 0; defaults are the MOCASSIN HII20/HII40
     !--- benchmark abundances.  use_metals enables the cascade + Tier-1
     !--- metal cooling inside the thermal balance.
     !--- Stage G4: solution-driven I-front re-refinement.  At iteration
     !--- refine_iter the octree is rebuilt: level refine_lmax where
     !--- refine_eps < x_HI < 1-refine_eps or a face-neighbor x_HI jump
     !--- exceeds refine_dx; refine_lbase elsewhere.  State maps by
     !--- position; the iteration then continues on the new grid.
     logical            :: refine_front = .false.
     integer            :: refine_iter  = 20
     integer            :: refine_lbase = 5
     integer            :: refine_lmax  = 7
     real(kind=wp)      :: refine_eps   = 1.0e-2_wp
     real(kind=wp)      :: refine_dx    = 0.2_wp
     !--- Stage G3: explicit diffuse ionizing field.  Ground-recombination
     !--- packets are emitted from every leaf and transported; forces
     !--- case A rates (case B double-counts the on-the-spot absorption).
     logical            :: diffuse_field = .false.
     logical            :: use_metals   = .false.
     real(kind=wp)      :: abund_C      = 2.2e-4_wp
     real(kind=wp)      :: abund_N      = 4.0e-5_wp
     real(kind=wp)      :: abund_O      = 3.3e-4_wp
     real(kind=wp)      :: abund_Ne     = 5.0e-5_wp
     real(kind=wp)      :: abund_S      = 9.0e-6_wp
     !--- G5: argon, magnesium, iron (not in the MOCASSIN benchmark set;
     !--- default off; Mg/Fe are strongly depleted onto grains — use
     !--- gas-phase values, e.g. Mg ~ 3e-6, Fe ~ 3e-7 in HII regions).
     real(kind=wp)      :: abund_Ar     = 0.0_wp
     real(kind=wp)      :: abund_Mg     = 0.0_wp
     real(kind=wp)      :: abund_Fe     = 0.0_wp
     !--- number of scatterings, these are not input parameters.
     real(kind=wp) :: nscatt_tot  = 0.0_wp
     real(kind=wp) :: exetime     = 0.0_wp
  end type params_type

  !--- observer type
  type observer_type
     real(kind=wp) :: x,y,z
     real(kind=wp) :: inclination_angle = 0.0_wp
     real(kind=wp) :: position_angle    = 0.0_wp
     real(kind=wp) :: phase_angle       = 0.0_wp
     real(kind=wp) :: alpha             = nan64
     real(kind=wp) :: beta              = nan64
     real(kind=wp) :: gamma             = nan64
     real(kind=wp) :: distance          = -999.9_wp
     real(kind=wp) :: rmatrix(3,3)
     integer       :: nxim,nyim
     real(kind=wp) :: dxim,dyim
     real(kind=wp) :: steradian_pix
     real(kind=wp), pointer :: tau(:,:)    => null()
     real(kind=wp), pointer :: scatt(:,:)  => null()
     !--- 4-D scattered image (nxim,nyim,n_albedo,n_hgg) used only when par%use_ag_list = .true.
     real(kind=wp), pointer :: scatt_ag(:,:,:,:) => null()
     !--- 5-D scattered image (nxim,nyim,n_albedo,n_hgg,n_tau) used when par%use_tau_list = .true.
     real(kind=wp), pointer :: scatt_agt(:,:,:,:,:) => null()
     !--- 3-D direct image (nxim,nyim,n_tau) used when par%use_tau_list = .true.
     real(kind=wp), pointer :: direc_t(:,:,:) => null()
     !--- 3-D wavelength-resolved images (nxim,nyim,nlambda) used when par%use_sed = .true.
     real(kind=wp), pointer :: scatt_sed(:,:,:)  => null()
     real(kind=wp), pointer :: direc_sed(:,:,:)  => null()
     real(kind=wp), pointer :: direc0_sed(:,:,:) => null()
     real(kind=wp), pointer :: direc(:,:)  => null()
     real(kind=wp), pointer :: direc0(:,:) => null()
     real(kind=wp), pointer :: I(:,:)      => null()
     real(kind=wp), pointer :: Q(:,:)      => null()
     real(kind=wp), pointer :: U(:,:)      => null()
     real(kind=wp), pointer :: V(:,:)      => null()
  end type observer_type

! scattering related global variables
  type scattering_matrix_type
     integer :: nPDF
     real(kind=wp), pointer :: coss(:)      => null()
     real(kind=wp), pointer :: S11(:)       => null()
     real(kind=wp), pointer :: S12(:)       => null()
     real(kind=wp), pointer :: S33(:)       => null()
     real(kind=wp), pointer :: S34(:)       => null()
     real(kind=wp), pointer :: phase_PDF(:) => null()
     integer,       pointer :: alias(:)     => null()
  end type scattering_matrix_type

  type params_mpi
     !--- p_rank    = process (thread) rank
     !--- h_rank    = host rank (rank defined for each process of a node)
     !--- nproc     = number of total processes
     !--- hostcomm  = communicator for host
     !--- SAME_HRANK_COMM  = communicator for processes with the same h_rank
     !--- SAME_HRANK_NPROC = number of processes with the same h_rank
     integer :: p_rank
     integer :: h_rank
     integer :: nproc
     integer :: hostcomm
     integer :: SAME_HRANK_COMM
     integer :: SAME_HRANK_NPROC
  end type params_mpi

  ! globar parameters
  type(params_type) :: par
  type(params_mpi)  :: mpar
  type(scattering_matrix_type) :: scatt_mat
  type(observer_type), allocatable :: observer(:)

  !--- MoCHII: scattered-peel hook (set to ion_peel_mod::ion_peel_scatter
  !--- during the imaging pass; a procedure pointer avoids the circular
  !--- dependence raytrace_amr <-> ion_peel).
  abstract interface
     subroutine peel_hook_intf(photon)
        import :: photon_type
        type(photon_type), intent(in) :: photon
     end subroutine peel_hook_intf
  end interface
  procedure(peel_hook_intf), pointer :: ion_peel_scatter_hook => null()

  !------------------------------------------------------------
  ! define procedure pointers
  procedure(raytrace_tau), pointer :: raytrace_to_tau => null()
  abstract interface
     subroutine raytrace_tau(photon,grid,tau_in)
     import
     type(photon_type), intent(inout) :: photon
     type(grid_type),   intent(inout) :: grid
     real(kind=wp),     intent(in)    :: tau_in
     end subroutine raytrace_tau
  end interface

  procedure(raytrace_edge), pointer :: raytrace_to_edge => null()
  abstract interface
     subroutine raytrace_edge(photon,grid,tau)
     import
     type(photon_type), intent(in)  :: photon
     type(grid_type),   intent(in)  :: grid
     real(kind=wp),     intent(out) :: tau
   end subroutine raytrace_edge
  end interface

  procedure(peeling_direct), pointer :: peeling_direct_photon => null()
  abstract interface
     subroutine peeling_direct(photon,grid)
     import
     type(photon_type), intent(in) :: photon
     type(grid_type),   intent(in) :: grid
   end subroutine peeling_direct
  end interface

  procedure(peeling_scattered), pointer :: peeling_scattered_photon => null()
  abstract interface
     subroutine peeling_scattered(photon,grid)
     import
     type(photon_type), intent(in) :: photon
     type(grid_type),   intent(in) :: grid
   end subroutine peeling_scattered
  end interface

  procedure(scatter), pointer :: scattering => null()
  abstract interface
     subroutine scatter(photon,grid)
     import
     type(photon_type), intent(inout) :: photon
     type(grid_type),   intent(inout) :: grid
     end subroutine scatter
  end interface

  procedure(run_sim), pointer :: run_simulation => null()
  abstract interface
     subroutine run_sim(grid)
     import
     type(grid_type),   intent(inout) :: grid
     end subroutine run_sim
  end interface

  procedure(write_out), pointer :: write_output => null()
  abstract interface
     subroutine write_out(filename,grid)
     import
     character(len=*),  intent(in) :: filename
     type(grid_type),   intent(in) :: grid
     end subroutine write_out
  end interface
  !------------------------------------------------------------
end module define
