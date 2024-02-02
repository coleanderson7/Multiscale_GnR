res
program main

use globalvars
use stress
use kinetics
use outputs
use experiments

integer :: ts, strs_flag, i, j, ldng_flag, t1, t2, clock_rate, clock_max, etime
double precision :: temp, t, dF, degmod, invar1, invar2, invar3, invar4
double precision :: invar5, invar6, invar7, invar8, invar9, invar10

call system_clock(t1, clock_rate, clock_max)

!Units - standard SI
open(unit=10, file = 'general.dat', status='replace', action='write')
open(unit=11, file = 'arealmasshist.dat', status='replace', action='write')
open(unit=12, file = 'massfractions.dat', status='replace', action='write')
open(unit=13, file = 'mkratehistory.dat', status='replace', action='write')
open(unit=14, file = 'qkratehistory.dat', status='replace', action='write')
open(unit=15, file = 'polymerlog.dat', status='replace', action='write')
open(unit=16, file = 'PDdata.dat', status='replace', action='write')

open(unit=19, file = 'GnRinput.dat', status='old', action='read')
read(19,*) invar1, invar2, invar3, invar4, invar5

!Introduce simulation length so we can allocate storage
days = 734 !time 0.25
steps = 1 !steps/day 732
dt = 1/real(steps)
nts = days*steps !total steps
s = 1 !Define start time index
t0 = 0 !Start time

!Establish number of constituents
k = 12
!k1 = PGA
!k2 = PCL-PLA
!k3-6 = strs-mediated col (ax, circ, diag, diag+90)
!k7-10 = infl-mediated col (ax, circ, diag, diag+90)
!k11 = strs SMCs
!k12 = infl SMCs

allocate(Gkh(k))
!Establish prestretches
Gkh(1) = 0 !Polymer felt
Gkh(2) = 0 !Polymer sealant
Gkh(3:6) = 1.08 !Stress mediated collagen families
Gkh(7:10) = 1.04 !Inflammation mediated collagen families
Gkh(11) = 1.14 !Stress/Infl SMCs
Gkh(12) = 1.14 !Stress/Infl SMCs

allocate(c2(k),c3(k),cc2(k),cc3(k))
!Establish consituent mechanical properities
!c2 is inifitesimal shear modulus for polymers
!young's modulus for tissue
!c3 is nonlinearity parameter
!General material parameters
c2_col = 2.00D4
c3_col = 0.50
c2_smc = 3.15D3
c3_smc = 3.15
c2_pga = 7E9 !PGA felt
c2_seal = invar1 !PCL-PLA sealant
gamma_mod2 = 3.27 !Enhancement for inflammatory collagen
gamma_mod3 = 0.50
c2_ground = 2E3 !Stiffness of ground matrix

!For this graft choice
!In tension
c2(1) = c2_pga
c2(2) = c2_seal
c2(3:6) = c2_col
c2(7:10) = gamma_mod2*c2_col
c2(11) = c2_smc
c2(12) = gamma_mod2*c2_smc

!Non-linear terms
c3(1) = 0.0
c3(2) = 0.0
c3(3:6) = c3_col
c3(7:10) = gamma_mod3*c3_col
c3(11) = c3_smc
c3(12) = gamma_mod3*c3_smc

!TC switch
cc2(3:12) = 0.3*c2(3:12)
cc3(3:12) = 0.3*c3(3:12)

!Initialize additional directional properties for ECM components
!Assume for now they do not vary directions with time
allocate(alphakh(k), alphak_tau(k,nts))
!Homeostatic orientations
alphakh(1) = 0*(pi/180.0)
alphakh(2) = 0*(pi/180.0)
alphakh(3:6) = [0, 90, 45, -45]*(pi/180.0) !Col
alphakh(7:10) = [0, 90, 45, -45]*(pi/180.0) !Infl Col
alphakh(11:12) = [90, 90]*(pi/180.0) !SMC, Infl SMC

!Inflamatory constituents only produced ~isotropcially
!Polymers never aligned
alphak_tau(1,:) = alphakh(1)
alphak_tau(2,:) = alphakh(2)
!mech collagen families deposited according to stress state
alphak_tau(3,:) = alphakh(3)
alphak_tau(4,:) = alphakh(4)
alphak_tau(5,:) = alphakh(5)
alphak_tau(6,:) = alphakh(6)
!infl collagen families always depsotied in homeostatic directions
alphak_tau(7,:) = alphakh(7)
alphak_tau(8,:) = alphakh(8)
alphak_tau(9,:) = alphakh(9)
alphak_tau(10,:) = alphakh(10)
!SMCs always circumferential
alphak_tau(11,:) = alphakh(11)
alphak_tau(12,:) = alphakh(12)

!Initialize Geometry
allocate(ril_s(nts), hl_s(nts), Ri_s(nts), H_s(nts), lamda_s(nts), lamdaul_s(nts), &
         R0_s(nts), rmidl_s(nts), Rmid_s(nts))

Ri = 450D-6
H = 290D-6
Rmid = Ri + H/2

Ri_s(:) = Ri
Rmid_s(:) = Rmid
R0_s(:) = Ri

lamda_s(:) = 1.0 !Assuming no stretch for now
lamdaul_s(:) = 1.0 !Unloading stretch initialization

!Initialize mass variables
rhob = 1.05D3 !mass density of blood/tissue

!Initialize polymer characteristics
rhopga = 1.53D3 !mass density of polymer (PGA in this case)
rhoseal = 1.23D3 !mass density of polymer (PCL-PLA in this case)
vfracpga = 0.115 !0.115 !1 - Porosity PGA
vfracseal = invar2 !0.10 !1 - Porosity Sealant
fiberdiameterpga = 16*1D-6 !16.0D-6 !Fiber diameter PGA
fiberdiameterseal = invar3*1D-6 !1.0D-6 !Fiber diameter PCL-PLA

!void state is trigger for filling up void space
void_state = 0

!Pore Size inflammation parameters
normpore = 10D-6 !Optimal pore size for cell infiltration
normfd = 5D-6
maxpore = 100D-6 !Maximum pore size detected as a pore rather than void
maxfd = 50D-6

!Calculate pore size
poresizepga = poresize_calc(fiberdiameterpga,vfracpga)
poresizeseal = poresize_calc(fiberdiameterseal,vfracseal)

allocate(rho_s(nts), rhop_s(2,nts), capM_s(nts), capMfiller_s(nts)) !Storage for  mixture density
allocate(vfracp_s(2,nts), poresize_s(2,nts), fiberdiameter_s(2,nts)) !Storage for evolving polymer properties
allocate(phik_s(k,nts), capMk_s(k,nts), mk_tau(k,nts), qk_tau(k,nts)) !Storage for evolving production/degradation rates
allocate(capQ_s(2,nts), capQrho_s(2,nts), scaffoldporesize_s(nts))
vfracp_s(1,:) = vfracpga
vfracp_s(2,:) = vfracseal

rhop_s(1,1) = rhopga
rhop_s(2,1) = rhoseal
rho_s(1) = vfracp_s(1,1)*rhop_s(1,1) + vfracp_s(2,1)*rhop_s(2,1) +&
           (1 - vfracp_s(1,1) - vfracp_s(2,1))*rhob !mixture density

fiberdiameter_s(1,1) = fiberdiameterpga !fiber diameter
fiberdiameter_s(2,1) = fiberdiameterseal
poresize_s(1,1) = poresizepga !pore size
poresize_s(2,1) = poresizeseal

!Smallest pore size is scaffold pore size
if (fiberdiameter_s(1,1) <= fiberdiameter_s(2,1)) then
    scaffoldporesize_s(s) = poresize_calc(fiberdiameterpga,vfracpga + vfracseal)
else
    scaffoldporesize_s(s) = poresize_calc(fiberdiameterseal,vfracpga + vfracseal)
end if

phik_s(1,:) = vfracp_s(1,1)*rhop_s(1,1)/rho_s(1)
phik_s(2,:) = vfracp_s(2,1)*rhop_s(2,1)/rho_s(1)
phik_s(3:k,:) = 0

capMk_s(1,:) = vfracp_s(1,1)*rhop_s(1,1)*H
capMk_s(2,:) = vfracp_s(2,1)*rhop_s(2,1)*H
capMk_s(3:k,:) = 0 !current area mass density of ECM constituent

capMfiller_s(:) = (1 - vfracp_s(1,1) - vfracp_s(2,1))*rhob*H

capM_s(:) = 0
capM_s(1) = sum(capMk_s(:,1)) +  capMfiller_s(1) !current total area mass density of constituent

mk_tau(:,:) = 0
qk_tau(:,:) = 1 !Degradation of produced constituents
capQ_s(1:2,:) = 1 !Surface degradation of initial constituents
capQrho_s(1:2,:) = 1 !Bulk degradation of initial constituents

rhop_s(1,:) = rhopga*capQrho_s(1,:)
rhop_s(2,:) = rhoseal*capQrho_s(2,:)

!Initialize phenomenologic kinetic parameters
allocate(kq(k), kqrho(2), kinflstatus(k))
allocate(Kkalpha(k), Kkwound(k), macalpha(k), macbeta(k))
allocate(Kkbstrs(k), Kkinfmax(k), mkb(k), tnot(2), tnotrho(2))
allocate(Kksigma(k), Kktau(k))

!List whether constituents are stress or inflammation mediated
!0 - Stress, 1 - Inflammation
kinflstatus(1:2) = 0
kinflstatus(3:6) = 0
kinflstatus(7:10) = 1
kinflstatus(11) = 0
kinflstatus(12) = 1

!Degradation rates
!Sigmoidal rate constants for Pol
kq(1) = 0.128*(4.01)!*(scaffoldporesize_s(1)/normpore + 1)) !For PGA
kq(2) = invar4*(4.01)!*(scaffoldporesize_s(1)/normpore + 1)) !For sealant
!Exponential rate constants for Tissue
kq(3:12) = 1/real(80) !For ECM
kq(7:10) = 1/real(80)
kq(12) = 1/real(80)
!Maximum life-span of tissue
ECMmax = 2000

!Density degradation rate for polymers
!Assuming surface erosion
kqrho(1:2) = 0

!Sigmoidal offsets for polymers
tnot(1) = 28.5/(4.01)!*(scaffoldporesize_s(1)/normpore + 1))
tnot(2) = invar5/(4.01)!*(scaffoldporesize_s(1)/normpore + 1))
!No bulk degradation
tnotrho(1:2) = 0

!Mass Production Parameters
!Basal mass production rates
mkb(1:2) = 0
mkb(3:4) = 1.225D-4 !1.25D-4*1.0!1.25D-4*5
mkb(5:6) = 0.6125D-4 !0.75D-4*1.0!0.75D-4*5
mkb(7:8) = 1.225D-4 !1.25D-4*1.0
mkb(9:10) = 0.6125D-4 !0.75D-4*1.0
mkb(11) = 1.3125D-4 !1.75D-4*1.0!1.75D-4*5
mkb(12) = 1.3125D-4 !1.75D-4

!Maturation rate
mat_rate = 1.01 !0.025

!Other inflammation mediated parameters
Kkalpha(1:k) = 3.31
Kkwound(1:k) = 1.0
macalpha(1:k) = 0.05
macbeta(1:k) = 3.08
Kkinfmax(1:k) = 3.73

!Other strs parameters
Kkbstrs(1:k) = 1.0 !1.0
Kksigma = 1
Kktau = 10 !1.0!20.0

!Declare loading variables
Ph = 2.0*133.33 !Pressure in Pa
Tauwh = 550D-6 !0.001 !WSS in Pa
mu = 0.0037 !Viscosity in kg/(m*s)
Qh = Tauwh*pi*Ri**3/(4*mu) !Flow in m^3/s
Sigmatth = 7000

!Determine if we're using phenomenologic model 0, or mechanistic model 1
infl_switch = 0

!Initialize load outputs storage
allocate(cauchy(2,nts),tension(2,nts),stretch(2,nts), Stiffness(2,nts), J_s(nts))
allocate(Pcalc(nts),Fcalc(nts),Tauw_s(nts), Pstore(nts), Area(nts), RComp(nts))
Pstore = 0
Stiffness = 0
J_s = 1
RComp = 0

!Set up for numerical experiments
no_tests = 2
allocate(Pexp(no_tests))
do i = 1, no_tests
    Pexp(i) = 1.0*133.33/1*(i - 1) + 2*133.33
end do

!Initialize some deformation vars
allocate(lk_s(k,nts),lk2_s(k,nts))
lk_s(:,:) = 1
lk2_s = 1

!Find initial loaded configuration with Secant Method
strs_switch = 0
call equil_gr()

!Find loaded thickness
hl_s(:) = findh(rmidl_s(s))
ril_s(s) = rmidl_s(s) - hl_s(s)/2

!Find loading quantities
strs_flag = 1
ldng_flag = 0
temp = 0
call loaded(temp,rmidl_s(s),dF,strs_flag,ldng_flag)
expcount = 1

!Ready to start G&R
do s = 1, nts - 1

    t = dt*(s - 1)
    call update_kinetics

    call run_experiments
    call update_outputs

end do

call system_clock(t2, clock_rate, clock_max)
etime = t2 - t1
!write(*,*) etime

end program main
