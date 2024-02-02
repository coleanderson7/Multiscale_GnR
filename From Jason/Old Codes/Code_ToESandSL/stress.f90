
module stress

implicit none
contains

subroutine equil_gr()

    use globalvars

    implicit none

    !Solving for current configuration of TEVG
    !Need to do a Newton loop on change in radius

    double precision :: r, rg, Frg, dFrg !returned radius
    double precision :: max_tol, beta, dF
    integer :: iter, max_iter, strs_flag, ldng_flag

    open(unit=7, file = 'NewtonLog.dat', status='old', action='write')

    !convergence criteria
    max_tol = 1D-8
    max_iter = 100

    if (strs_switch == 1) then
        strs_flag = 1
    else
        strs_flag = 0
    end if

    ldng_flag = 0
    beta = 1.00

    iter = 1
    if (s > 1) then
        rg = rmidl_s(s - 1)*1.00
        if (strs_switch == 1) then
            rg = rmidl_s(s)
        end if
    else
        rg = Rmid
    end if

    dFrg = 1
    Frg = 10000

    do while (iter < max_iter .and. dabs(Frg) > max_tol)
        call loaded(Frg,rg,dF,strs_flag,ldng_flag)
        dFrg = Frg/dF !Frg1*(rg1 - rg)/(Frg1 - Frg) !Evaluate derivative

        write(7,*) dt*(s - 1), iter, dFrg, Frg, rg

        rg = rg - beta*dFrg !Find next guess

        iter = iter + 1
    end do

    write(7,*) '-------------------------'

    if (iter == max_iter) then
        !write(7,*) 'Failed radial convergence'
    end if

    rmidl_s(s) = rg !return converged guess

end subroutine equil_gr

subroutine loaded(obj, rg, dF, strs_flag, ldng_flag)

    use globalvars

    implicit none

    integer :: dir, strs_flag, ks, ldng_flag
    double precision :: rg, obj, lt, lz, P, dF
    double precision :: sum_dWkdl_nh, sum_dWkdl_exp, sum_ddWkddl_nh, sum_ddWkddl_exp
    double precision :: ddWkddC, hcurr

    !open(unit=8, file = 'NewtonCalcLog.dat', status='old', action='write')

    dir = 2

    if (ldng_flag == 0) then
        !For G&R equilibrium
        lt = rg/Rmid
        lz = lamda_s(s)
        hcurr = hl_s(s)
    elseif (ldng_flag == 1) then
        !For solving the unloaded configuration
        lt = rg/Rmid
        lz = lamdaul_s(s)
        hcurr = H_s(s)
    end if

    do ks = 3,k
        lk_s(ks,s) = sqrt(lz**2*cos(alphak_tau(ks,s))**2 + &
                     lt**2*sin(alphak_tau(ks,s))**2)
    end do

    sum_dWkdl_nh = 0
    sum_dWkdl_exp = 0
    sum_ddWkddl_nh = 0
    sum_ddWkddl_exp = 0

    !For polymers
    sum_dWkdl_nh = dWkdl_nH(1,lt,lz,dir) + dWkdl_nH(2,lt,lz,dir)
    sum_ddWkddl_nh = ddWkddl_nH(1,lt,lz,dir) + ddWkddl_nH(2,lt,lz,dir)

    !For ECM/cells
    if (s > 1) then
        do ks = 3,k
            sum_dWkdl_exp = sum_dWkdl_exp + dWkdl_exp(ks,lt,lz,dir)
            sum_ddWkddl_exp = sum_ddWkddl_exp + ddWkddl_exp(ks,lt,lz,dir)
        end do
    end if

    P = 1/lz*(sum_dWkdl_nh + sum_dWkdl_exp)/(rg - hcurr)
    dF = 1/lz*(sum_ddWkddl_nh + sum_ddWkddl_exp)/(rg - hcurr) - Ph

    !write(8,*) dt*(s - 1), P, lk_s(4,s)/lk_s(4,s - 1), &
               !sum_dWkdl_nh, sum_dWkdl_exp

    obj = 1/lz*(sum_dWkdl_nh + sum_dWkdl_exp) - Ph*(rg - hcurr)

    if (strs_flag == 1) then

        !Calculate and store loading metrics
        dir = 2
        stretch(1,s) = lt
        tension(1,s) = 1/lz*(sum_dWkdl_nh + sum_dWkdl_exp)
        cauchy(1,s) = tension(1,s)/hcurr
        Pcalc(s) = P

        !Need to solve for strain energy in the other direction
        dir = 3
        sum_dWkdl_nh = 0
        sum_dWkdl_exp = 0
        sum_ddWkddl_nh = 0
        sum_ddWkddl_exp = 0

        !For polymers
        sum_dWkdl_nh = dWkdl_nH(1,lt,lz,dir) + dWkdl_nH(2,lt,lz,dir)
        sum_ddWkddl_nh = ddWkddl_nH(1,lt,lz,dir) + ddWkddl_nH(2,lt,lz,dir)

        !For ECM/cells
        if (s > 1) then
            do ks = 3,k
                sum_dWkdl_exp = sum_dWkdl_exp + dWkdl_exp(ks,lt,lz,dir)
                sum_ddWkddl_exp = sum_ddWkddl_exp + ddWkddl_exp(ks,lt,lz,dir)
            end do
        end if

        J_s(s) = lt*lz
        stretch(2,s) = lz
        tension(2,s) = 1/lt*(sum_dWkdl_nh + sum_dWkdl_exp)
        cauchy(2,s) = tension(2,s)/hcurr
        Fcalc(s) = pi*(2*rg + hcurr)*tension(2,s)

        Tauw_s(s) = (4*mu*Qh)/(pi*(rg - hcurr)**3)

    end if

end subroutine loaded

function dWkdl_nH(ks,lt,lz,dir)

    use globalvars

    implicit none

    integer, intent(in) :: dir, ks
    double precision, intent(in) :: lt, lz
    double precision :: dWkdl_nH, dWnativedl_nH
    !Returns derivative of neo-Hookean strain energy function

    if (dir == 2) then !circumferential

        dWkdl_nH = c2(ks)*poldeg(ks)*(lt - 1/(lz**2*lt**3))*&
                   capMk_s(ks,s)/rho_s(s)

        !Polymer becomes ground matrix
        if (c2(ks)*poldeg(ks) < c2_ground) then
            dWkdl_nH = c2_ground*(lt - 1/(lz**2*lt**3))*&
                       capMk_s(ks,s)/rho_s(s)
        end if


    elseif (dir == 3) then

        dWkdl_nH = c2(ks)*poldeg(ks)*(lz - 1/(lz**3*lt**2))*&
                   capMk_s(ks,s)/rho_s(s)

        if (c2(ks)*poldeg(ks) < c2_ground) then
            dWkdl_nH = c2_ground*(lz - 1/(lz**3*lt**2))*&
                       capMk_s(ks,s)/rho_s(s)
        end if

    end if

end function dWkdl_nH

function ddWkddl_nH(ks,lt,lz,dir)

    use globalvars

    implicit none

    integer, intent(in) :: dir, ks
    double precision, intent(in) :: lt, lz
    double precision :: ddWkddl_nH
    !Returns derivative of neo-Hookean strain energy function

    if (dir == 2) then !circumferential

        ddWkddl_nH = c2(ks)*poldeg(ks)*(1 + 3/(lz**2*lt**4))*&
                   capMk_s(ks,s)/rho_s(s)

        if (c2(ks)*poldeg(ks) < c2_ground) then
            ddWkddl_nH = c2_ground*(1 + 3/(lz**2*lt**4))*&
                       capMk_s(ks,s)/rho_s(s)
        end if


    elseif (dir == 3) then

        ddWkddl_nH = c2(ks)*poldeg(ks)*(1 + 3/(lz**4*lt**2))*&
                   capMk_s(ks,s)/rho_s(s)

        if (c2(ks)*poldeg(ks) < c2_ground) then
            ddWkddl_nH = c2_ground*(1 + 3/(lz**4*lt**2))*&
                       capMk_s(ks,s)/rho_s(s)
        end if

    end if

end function

function dWkdl_exp(ks,lt,lz,dir)

    use globalvars

    implicit none

    integer, intent(in) :: dir, ks
    integer :: ntau
    double precision, intent(in) :: lt, lz
    double precision :: dWkdl_exp, str_energ_cohort

    !Calculate the constituent specific stress of a cohort produced at time ntau
    dWkdl_exp = 0
    str_energ_cohort = 0

    if (dir == 2) then

        do ntau = 1,s

            str_energ_cohort = Gkh(ks)/lk_s(ks,ntau)*dWkdlkntau_exp(ks,ntau)* &
                               lt*sin(alphak_tau(ks,ntau))**2/lk_s(ks,s)

            dWkdl_exp = dWkdl_exp + mk_tau(ks,ntau)*qk_tau(ks,ntau)*str_energ_cohort*dt

        end do

        dWkdl_exp = dWkdl_exp/rho_s(s)

    elseif (dir == 3) then

        do ntau = 1,s

            str_energ_cohort = Gkh(ks)/lk_s(ks,ntau)*dWkdlkntau_exp(ks,ntau)* &
                               lz*cos(alphak_tau(ks,ntau))**2/lk_s(ks,s)

            dWkdl_exp = dWkdl_exp + mk_tau(ks,ntau)*qk_tau(ks,ntau)*str_energ_cohort*dt

        end do

        dWkdl_exp = dWkdl_exp/rho_s(s)

    end if

end function dWkdl_exp

function ddWkddl_exp(ks,lt,lz,dir)

use globalvars

    implicit none

    integer, intent(in) :: dir, ks
    integer :: ntau
    double precision, intent(in) :: lt, lz
    double precision :: ddWkddl_exp, str_energ_cohort

    !Calculate the constituent specific stress of a cohort produced at time ntau
    ddWkddl_exp = 0
    str_energ_cohort = 0

    if (dir == 2) then

        do ntau = 1,s

            str_energ_cohort = dWkdlkntau_exp(ks,ntau)*Gkh(ks)/lk_s(ks,ntau)* &
                               (lz*sin(alphak_tau(ks,ntau))*cos(alphak_tau(ks,ntau)))**2 &
                               /lk_s(ks,s)**3 + ddWkddlkntau_exp(ks,ntau)* &
                               (Gkh(ks)/lk_s(ks,ntau)*lt*sin(alphak_tau(ks,ntau))**2/lk_s(ks,s))**2

            ddWkddl_exp = ddWkddl_exp + mk_tau(ks,ntau)*qk_tau(ks,ntau)*str_energ_cohort*dt

        end do

        ddWkddl_exp = ddWkddl_exp/rho_s(s)

    elseif (dir == 3) then

        do ntau = 1,s

            str_energ_cohort = dWkdlkntau_exp(ks,ntau)*Gkh(ks)/lk_s(ks,ntau)* &
                               (lt*cos(alphak_tau(ks,ntau))*sin(alphak_tau(ks,ntau)))**2 &
                               /lk_s(ks,s)**3 + ddWkddlkntau_exp(ks,ntau)* &
                               (Gkh(ks)/lk_s(ks,ntau)*lz*cos(alphak_tau(ks,ntau))**2/lk_s(ks,s))**2

            ddWkddl_exp = ddWkddl_exp + mk_tau(ks,ntau)*qk_tau(ks,ntau)*str_energ_cohort*dt

        end do

        ddWkddl_exp = ddWkddl_exp/rho_s(s)

    end if

end function ddWkddl_exp

function dWkdlkntau_exp(ks,ntau)

    use globalvars

    implicit none

    integer, intent(in) :: ks, ntau
    double precision :: dWkdlkntau_exp, lkntau_s, t, tau
    double precision :: modc2, modc3, modcc2, modcc3
    t = dt*(s - 1)
    tau = dt*(ntau - 1)
    modc2 = xlink(t,tau)
    modc3 = xlink(t,tau)
    modcc2 = modc2
    modcc3 = modc3

    !Possibility for TC switch
    lkntau_s = Gkh(ks)*(lk_s(ks,s)/lk_s(ks,ntau))

    if (lkntau_s >= 1.00) then
        dWkdlkntau_exp = modc2*c2(ks)*lkntau_s*(lkntau_s**2 - 1)* &
                         exp(modc3*c3(ks)*(lkntau_s**2 - 1)**2)
    else
        dWkdlkntau_exp = modcc2*cc2(ks)*lkntau_s*(lkntau_s**2 - 1)* &
                         exp(modcc3*cc3(ks)*(lkntau_s**2 - 1)**2)
    end if

end function dWkdlkntau_exp

function ddWkddlkntau_exp(ks,ntau)

    use globalvars

    implicit none

    integer, intent(in) :: ks, ntau
    double precision :: ddWkddlkntau_exp, lkntau_s, t, tau
    double precision :: modc2, modc3, modcc2, modcc3
    t = dt*(s - 1)
    tau = dt*(ntau - 1)
    modc2 = xlink(t,tau)
    modc3 = xlink(t,tau)
    modcc2 = modc2
    modcc3 = modc3

    !Possibility for TC switch
    lkntau_s = Gkh(ks)*(lk_s(ks,s)/lk_s(ks,ntau))

    if (lkntau_s >= 1.00) then
        ddWkddlkntau_exp = modc2*c2(ks)*(3*lkntau_s**2 - 1 + &
                           4*modc3*c3(ks)*(lkntau_s*(lkntau_s**2 - 1))**2)* &
                          exp(modc3*c3(ks)*(lkntau_s**2 - 1)**2)
    else
        ddWkddlkntau_exp = modcc2*cc2(ks)*(3*lkntau_s**2 - 1 + &
                           4*modcc3*cc3(ks)*(lkntau_s*(lkntau_s**2 - 1))**2)* &
                          exp(modcc3*cc3(ks)*(lkntau_s**2 - 1)**2)
    end if

end function ddWkddlkntau_exp

function poldeg(ks)

    use globalvars

    implicit none

    integer, intent(in) :: ks
    double precision :: swellingfactor
    double precision :: poldeg

    !Returns changes to polymer due to degradation/porosity
    !Surface degradation
    swellingfactor = 0.03
    poldeg = 3*swellingfactor*vfracp_s(ks,s)**2

    !if (poldeg < 0.02) then
        !poldeg = 0.02
    !end if

end function poldeg

function findh(rg)

    use globalvars

    implicit none

    double precision, intent(in) :: rg
    double precision :: lt, lz, J
    double precision :: findh

    lt = rg/Rmid
    lz = lamda_s(s)

    J = lt*lz

    findh = (capM_s(s))/(J*rho_s(s))

end function findh

subroutine update_ECMalignments()

    !This function finds the new angle at which mechanically mediated collagen
    !is placed from one time step to the next
    !Governed by principal stresses according to Baek 2006

    use globalvars

    implicit none

    integer :: i, ks
    double precision :: f1, f2, lz, lt, sigma1, sigma2

    lt = stretch(1,s - 1)
    lz = stretch(2,s - 1)

    if (dabs(cauchy(1,s - 1)) > dabs(cauchy(2,s - 1))) then
        sigma1 = cauchy(1,s - 1)
        sigma2 = cauchy(2,s - 1)
    else
        sigma1 = cauchy(2,s - 1)
        sigma2 = cauchy(1,s - 1)
    end if

    f1 = sigma2/sqrt(sigma1**2 + sigma2**2)
    f2 = sigma1/sqrt(sigma1**2 + sigma2**2)

    alphak_tau(5,s) = atan(f1/f2)
    alphak_tau(6,s) = -alphak_tau(5,s)

end subroutine update_ECMalignments

function xlink(t,tau)

    !Controls rate of maturation of collagen

    use globalvars

    implicit none

    double precision :: t, tau, t2, xlink, init

    t2 = t - tau
    init = 0.0

    xlink = init + (1 - init)*(1 - exp(-mat_rate*t2))

end function xlink



end module stress
