! nameofvariablesuperscriptsubscript_functionof
!cr - core, st - sheath, k - constituent

module globalvars

implicit none

!Define pi
double precision, parameter :: pi = 2.0D0*dacos(0.0D0)

!Time step and space discretization variables
integer :: days, nts, s!time, steps/day, total steps
double precision :: dt, steps

!Number of constituents
integer :: k, void_state

!Mass density prescribed vars
double precision :: mu, rhob, rhopga, rhoseal !actual rho of each constituent
double precision :: vfracpga, poresizepga, fiberdiameterpga, vfracseal, poresizeseal, fiberdiameterseal
double precision, allocatable, dimension(:) :: rho_s,  capM_s, scaffoldporesize_s, capMfiller_s
double precision, allocatable, dimension(:,:) :: rhop_s, vfracp_s, poresize_s, fiberdiameter_s
double precision, allocatable, dimension(:,:) ::  capQ_s, capQrho_s, phik_s, capMk_s, mk_tau, qk_tau !calculated rho over time
integer :: infl_switch

!Deformation prescribed variables
double precision, allocatable, dimension(:)  :: Gkh ! homeostatic prestretch
double precision, allocatable, dimension(:,:) :: lk_s, lk2_s, alphak_tau

!Vessel geometry
double precision :: Ri, H, Rmid
double precision, allocatable, dimension(:) :: ril_s, hl_s, Ri_s, H_s, R0_s, rmidl_s, Rmid_s !Evolving radial subdomains over time
double precision, allocatable, dimension(:) :: lamda_s, lamdaul_s !Evolving axial stretch over time

!Constituent mechanical parameters
double precision, allocatable, dimension(:) :: c2, c3 !tensile elastic parameters, polymers have c3 = 0
double precision, allocatable, dimension(:) :: cc2, cc3 !compressive elastic parameters, polymers have cc3 = 0
double precision, allocatable, dimension(:) :: alphakh
double precision :: gmmod, gamma_mod2, c2_col, c3_col, c2_smc, c3_smc, c2_pga, c2_seal !Stiffening for inflamed constituents
double precision :: c2_ground, phiground, gamma_mod3, normfd, maxfd

!Kinetic parameters
double precision, allocatable, dimension(:) :: kq, kqrho, kinflstatus !degradation rates
double precision, allocatable, dimension(:) :: Kkalpha, Kkwound, macalpha, macbeta !infl production parameters
double precision, allocatable, dimension(:) :: Kkbstrs, Kkinfmax, Kktau, Kksigma
double precision, allocatable, dimension(:) :: mkb !production rates
double precision, allocatable, dimension(:) :: tnot, tnotrho !sigmoid midpoint
double precision :: normpore, maxpore, ECMmax, mat_rate !Pore size for normalization

!Loading conditions
double precision :: Ph, Tauwh, Qh, Sigmatth

!Experimental Comparisons
integer :: expcount, no_tests, strs_switch
double precision, allocatable, dimension(:) :: exptimes, Pexp, Area, RComp
double precision, allocatable, dimension(:,:) :: PDdata, Stiffness

!Construct stress and stretch history
double precision, allocatable, dimension(:,:) :: cauchy, tension, stretch
double precision, allocatable, dimension(:) :: Pcalc, Fcalc, Tauw_s, J_s, Pstore


end module globalvars
