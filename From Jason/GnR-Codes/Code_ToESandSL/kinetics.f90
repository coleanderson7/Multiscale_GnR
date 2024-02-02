
module kinetics

implicit none
contains

subroutine update_kinetics()

    use globalvars
    use stress
    use experiments

    implicit none

    integer :: ks, ts, strs_flag, ldng_flag, i
    double precision :: capMkintegral, temp, dF

    if (s == 1) then

        call mk_tau_fun()

        call qk_tau_fun()

    else

        !Predict kinetic relations
        call predict_kinetics()

        !Correct mass estimates
        do i = 1, 1
            call correct_kinetics()
        end do

        !Pstore(s) = Pcalc(s)

        !Find the unloaded vessel configuration
        call find_unloaded_config()

        !Find the linearized stiffness quantities
        call find_linearized_stiffness()

        !Find in vivo loading quantities
        call equil_gr()
        strs_flag = 1
        ldng_flag = 0
        temp = 0
        call loaded(temp,rmidl_s(s),dF,strs_flag,ldng_flag)

        if (ril_s(s) .NE. ril_s(s)) then
            stop
        end if

        if (Ri_s(s) .NE. Ri_s(s)) then
            stop
        end if

        if (hl_s(s) .NE. hl_s(s) .or. hl_s(s) > 1) then
            stop
        end if

        if (Stiffness(1,s) .NE. Stiffness(1,s) .or. Stiffness(1,s) > 10D15) then
            stop
        end if

        call update_ECMalignments()

    end if

end subroutine update_kinetics

subroutine predict_kinetics()

    use globalvars
    use outputs
    use stress

    implicit none

    integer :: ks, ts, strs_flag, ldng_flag
    double precision :: capMkintegral, temp, dF, t
    double precision :: matrix_frac, KKinf, kkq

    !Produce mass based on kinetic quantitites based on rates from the previous
    !time steps

    t = dt*(s - 1)

    !Polymers present initially
    call update_polymer()
    capM_s(s) = capMk_s(1,s) + capMk_s(2,s)

    !Update constituent alignments
    call update_ECMalignments()

    !Enhancement to initial matrix degradation
    kkq = kq(7)

    !ECM produced after time 0
    capMfiller_s(s) = capMfiller_s(s - 1)
    do ks = 3, k
        capMkintegral = 0.0
        do ts = 1, s - 1
            capMkintegral = capMkintegral + mk_tau(ks,ts)*qk_tau(ks,ts)*dt!/2
        end do
        capMk_s(ks,s) = capMkintegral
        capM_s(s) = capM_s(s) + capMk_s(ks,s)

        !Update the amount of interstitial water
        if (t > 14) then
            matrix_frac = 0.10 !Fraction of initial ground matrix that persists in the graft
            capMfiller_s(s) = capMfiller_s(1)*((1 - matrix_frac)*exp(-kkq*(t - 14)) + matrix_frac)
                              !- *mk_tau(ks,s - 1)*dt
        end if

    end do

    if (capMfiller_s(s) > 0 .and. void_state .ne. 1) then
        capMfiller_s(s) = capMfiller_s(s)
    else
        capMfiller_s(s) = 0
        void_state = 1
    end if

    capM_s(s) = sum(capMk_s(:,s)) + capMfiller_s(s)

    if (dabs(capM_s(s) - capM_s(s - 1)) < 0.001*capM_s(s)) then
        capM_s(s) = capM_s(s - 1)
    end if

    !Update mass fractions for strain energy calculation
    do ks = 1,k
        phik_s(ks,s) = capMk_s(ks,s)/capM_s(s)
    end do

    !Save some polymer as ground matrix
    if (phik_s(1,s)<=0.01) then
        capMk_s(1,s) = 0.01*capM_s(s)
    end if

    if (phik_s(2,s)<=0.01) then
        capMk_s(2,s) = 0.01*capM_s(s)
    end if

    !Recalculate mass fractions in case polymer became ground substance
    do ks = 1,k
        phik_s(ks,s) = capMk_s(ks,s)/capM_s(s)
    end do

    rho_s(s) = vfracp_s(1,s)*rhop_s(1,s) &
                + vfracp_s(2,s)*rhop_s(2,s) + &
                (1 - vfracp_s(1,s) - vfracp_s(2,s))*rhob

    if (rho_s(s)/rhob - 1 < 0.005) then
        rho_s(s) = rhob
    end if

    !Solve for mechanics based on newly produced mass
    !Find new equilibrium configuration
    call equil_gr()

    if (ril_s(s) .NE. ril_s(s)) then
        write(*,*) s, 5
        call update_outputs
        stop
    end if

    !Update loaded thickness
    hl_s(s) = findh(rmidl_s(s))
    ril_s(s) = rmidl_s(s) - hl_s(s)/2

    !Find loading quantities
    strs_flag = 1
    ldng_flag = 0
    temp = 0
    call loaded(temp,rmidl_s(s),dF,strs_flag,ldng_flag)

    call mk_tau_fun()

    call qk_tau_fun()

end subroutine predict_kinetics

subroutine correct_kinetics()

    use globalvars
    use outputs
    use stress

    implicit none

    integer :: ks, ts, strs_flag, ldng_flag
    double precision :: capMkintegral, temp, dF, t
    double precision :: matrix_frac, Kkinf, kkq

    t = dt*(s - 1)

    call update_polymer()
    capM_s(s) = capMk_s(1,s) + capMk_s(2,s)

    !Update constituent alignments
    call update_ECMalignments

    Kkinf = Kkalpha(7)*2*scaffoldporesize_s(1)/normpore + Kkwound(7)
    kkq = kq(7)!*(1.0 + Kkinf/Kkinfmax(7))

    !ECM produced after time 0
    capMfiller_s(s) = capMfiller_s(s - 1)
    do ks = 3, k
        capMkintegral = 0.0
        do ts = 2, s
            capMkintegral = capMkintegral + (mk_tau(ks,ts)*qk_tau(ks,ts) + &
                                             mk_tau(ks,ts - 1)*qk_tau(ks,ts - 1))*dt/2
        end do
        capMk_s(ks,s) = capMkintegral
        capM_s(s) = capM_s(s) + capMk_s(ks,s)

        !Update the amount of interstitial water
        if (t > 14) then
            matrix_frac = 0.10 !Fraction of initial ground matrix that persists in the graft
            capMfiller_s(s) = capMfiller_s(1)*((1 - matrix_frac)*exp(-kkq*(t - 14)) + matrix_frac)
        end if

    end do

    if (capMfiller_s(s) > 0 .and. void_state .ne. 1) then
        capMfiller_s(s) = capMfiller_s(s)
    else
        capMfiller_s(s) = 0
        void_state = 1
    end if

    capM_s(s) = sum(capMk_s(:,s)) + capMfiller_s(s)

    if (dabs(capM_s(s) - capM_s(s - 1)) < 0.001*capM_s(s)) then
        capM_s(s) = capM_s(s - 1)
    end if

    !Update mass fractions for strain energy calculation
    do ks = 1,k
        phik_s(ks,s) = capMk_s(ks,s)/capM_s(s)
    end do

    !Save some polymer as ground matrix
    if (phik_s(1,s)<=0.01) then
        capMk_s(1,s) = 0.01*capM_s(s)
    end if

    if (phik_s(2,s)<=0.01) then
        capMk_s(2,s) = 0.01*capM_s(s)
    end if

    !Recalculate mass fractions in case polymer became ground substance
    do ks = 1,k
        phik_s(ks,s) = capMk_s(ks,s)/capM_s(s)
    end do

    rho_s(s) = vfracp_s(1,s)*rhop_s(1,s) &
               + vfracp_s(2,s)*rhop_s(2,s) + &
              (1 - vfracp_s(1,s) - vfracp_s(2,s))*rhob

    if (rho_s(s)/rhob - 1 < 0.005) then
        rho_s(s) = rhob
    end if

    !Solve for mechanics based on newly produced mass
    !Find new equilibrium configuration
    call equil_gr()

    if (ril_s(s) .NE. ril_s(s)) then
        write(*,*) s, 6
        call update_outputs
        stop
    end if

    !Update loaded thickness
    hl_s(s) = findh(rmidl_s(s))
    ril_s(s) = rmidl_s(s) - hl_s(s)/2

    !Find loading quantities
    strs_flag = 1
    ldng_flag = 0
    temp = 0
    call loaded(temp,rmidl_s(s),dF,strs_flag,ldng_flag)

    call mk_tau_fun()

    call qk_tau_fun()

end subroutine correct_kinetics

!set mkr for each radial point based on the current time
subroutine mk_tau_fun()

    use globalvars

    implicit none

    double precision :: t, t0, deltaSigma, deltaTauw, strsshielding, taumod, peffect
    double precision :: steady_inf, pore_fun, fd_fun
    double precision, dimension(k) :: Kkstrs, mkstrs, Kkinf, mkinf
    integer :: i, j, ks

    !open(unit=14, file = 'infllog.dat', status='old', action='write')

    !Actual time
    t = (s - 1)*dt
    t0 = 0

    !Calculate mechanically mediated matrix production
    Kkstrs(3:k) = 1!0.1*(2*scaffoldporesize_s(1)/normpore + 1)

    !Higher than normal circumferential stress causes production
    deltaSigma = cauchy(1,s)/Sigmatth - 1

    !Higher than normal wall stress will decrease production
    deltaTauw = Tauwh**3/(ril_s(s))**3 - 1 !Tauw_s(s)/Tauwh - 1

    !Add choice of polmyers for stress shielding term

    mkstrs(3:k) = (Kksigma(3:k)*deltaSigma - Kktau(3:k)*deltaTauw)
                  !(1 - capQ_s(1,s)**4)*Kkstrs(3:k)*&Kkstrs(3:k)*Kkstrs(3:k)*

    fd_fun = 1.0
    pore_fun = 1.0

    !Calculate inflammatory mediated matrix production
    if (s > 1) then
        pore_fun = scaffoldporesize_s(2)/normpore
        fd_fun = fiberdiameter_s(2,2)/normfd
        Kkinf(3:k) = Kkalpha(3:k)*(pore_fun*fd_fun) &
                     + Kkwound(3:k) !
    else
        Kkinf(3:k) = 0
    end if

    mkinf(3:k) = Kkinf(3:k)*macalpha**macbeta*t**(macbeta - 1) &
                 *exp(-macalpha*t)/&
                 (macalpha*(macbeta - 1)**(macbeta - 1)*exp(1 - macbeta))

    if ((fd_fun + 1/pore_fun - 2) > 0) then
        steady_inf = (1 - exp(-macalpha(7)*t))*&
                     Kkalpha(7)*(fd_fun + 1/pore_fun - 2)
    else
        steady_inf = 0
    end if

    mkinf(3:k) = mkinf(3:k) + steady_inf

    !mkinf(3:k) = (Kkinf(3:k)*macalpha*t*exp(1 - macalpha*t))**1

    !Seperate constituents into stress or inflammatory mediated production
    !(1 - exp(-s)) for delay during cell infiltration
    do ks = 3,k

        if (kinflstatus(ks) == 1) then
            mk_tau(ks,s) = mkb(ks)*(mkinf(ks)) ! + capQ_s(2,s)) !!*(1 - exp(-t))
            if (mk_tau(ks,s) < 1D-10) then
                mk_tau(ks,s) = 0
            end if
        elseif (mkstrs(ks) + 1 > 0.1) then
            mk_tau(ks,s) = mkb(ks)*(mkstrs(ks) + 1) !*(1 - exp(-t))
        else
            mk_tau(ks,s) = mkb(ks)/10 ! mkb(ks)*(1 - exp(-t))*(1 - exp(-t))
        end if

    end do

end subroutine mk_tau_fun

!set qkr for each radial point based on the current time
! and the production time
subroutine qk_tau_fun()

    use globalvars

    implicit none

    integer :: nt, ntau, i, j, js, ks, rs, ts, tss
    double precision :: t, tau, kkq, pore_fun, fd_fun
    double precision, dimension(k) :: Kkinf

    t = (s - 1)*dt
    if (s > 1) then
        pore_fun = scaffoldporesize_s(2)/normpore
        fd_fun = fiberdiameter_s(2,2)/normfd
        Kkinf(3:k) = Kkalpha(3:k)*(pore_fun*fd_fun) &
                     + Kkwound(3:k) !
        !if (Kkinf(7) > Kkinfmax(ks)) then
    else
        Kkinf(3:k) = 0
    end if

    !Need to integrate from deposition time to current time
    do ntau = 1, s - 1
        tau = (ntau - 1)*dt
        do ks = 3,k

            !qk_tau(ks,ntau) = 1

            if (t - tau > ECMmax) then
                qk_tau(ks,ntau) = 0.0
            elseif (kinflstatus(ks) == 1) then

                qk_tau(ks,ntau) = 0.0
                do ts = ntau, s - 1
                        kkq = kq(ks)*(1.0 + Kkinf(ks)/Kkinfmax(ks))
                        qk_tau(ks,ntau) = qk_tau(ks,ntau) + kkq*dt

                end do
                qk_tau(ks,ntau) = exp(-qk_tau(ks,ntau))

            else

                qk_tau(ks,ntau) = 0.0
                do ts = ntau, s - 1
                        kkq = kq(ks)*(1.0 + (lk_s(ks,s)/lk_s(ks,ts) - 1)**2)
                        qk_tau(ks,ntau) = qk_tau(ks,ntau) + kkq*dt
                end do
                qk_tau(ks,ntau) = exp(-qk_tau(ks,ntau)) !exp(-kq(ks)*(t - tau)) !

            end if
        end do

    end do



end subroutine qk_tau_fun

subroutine update_polymer()
    use globalvars

    implicit none

    integer :: nt, i, rs, ks, ts
    double precision :: t, t0, Q1, Q2, Kkinf, kkq, mkinf

    !We need to update the degradation of the polymer constituents
    !Final "output" is the updated areal mass densities
    !Also need to update relevant microstructural characteristics for
    !mechanical properties

    t = (s - 1)*dt !We need an actual time
    t0 = 0

    !We prescribe the fiber diameter degradation profile as Q
    !Inflammation dependent polymer degradation
    !Q1 = 0
    !Q2 = 0
    !do ts = 1, s
        !Kkinf = Kkalpha(1)*2*scaffoldporesize_s(ts)/normpore + Kkwound(1)
        !mkinf = (Kkinf*macalpha(1)*t*exp(1 - macalpha(1)*t))**macbeta(1)
        !kkq = 1.0 !+ mkinf)
        !Q1 = Q1 + kq(1)*kkq*dt
        !Q2 = Q2 + kq(2)*kkq*dt
    !end do
    !capQ_s(1,s) = exp(-Q1)
    !capQ_s(2,s) = exp(-Q2)

    capQ_s(1,s) = sigmoidal(t,t0,kq(1),tnot(1))
    capQ_s(2,s) = sigmoidal(t,t0,kq(2),tnot(2))

    capQrho_s(1,s) = sigmoidal(t,t0,kqrho(1),tnotrho(1))
    capQrho_s(2,s) = sigmoidal(t,t0,kqrho(2),tnotrho(2))

    !Assume some residual polymer
    if (capQ_s(1,s) < 0.01) then
        capQ_s(1,s) = 0.01
    end if

    if (capQ_s(2,s) < 0.01) then
        capQ_s(2,s) = 0.01
    end if

    !Need to update porosity, fiber diameter and calculate updated pore size
    vfracp_s(1,s) = vfracp_s(1,1)*capQ_s(1,s)**2
    vfracp_s(2,s) = vfracp_s(2,1)*capQ_s(2,s)**2

    fiberdiameter_s(1,s) = fiberdiameter_s(1,1)*capQ_s(1,s)
    fiberdiameter_s(2,s) = fiberdiameter_s(2,1)*capQ_s(2,s)

    poresize_s(1,s) = poresize_calc(fiberdiameter_s(1,s),vfracp_s(1,s))
    poresize_s(2,s) = poresize_calc(fiberdiameter_s(2,s),vfracp_s(2,s))

    if (fiberdiameter_s(1,1) <= fiberdiameter_s(2,1)) then
        scaffoldporesize_s(s) = poresize_calc(fiberdiameter_s(1,s),vfracp_s(1,s) + vfracp_s(2,s))
    else
        scaffoldporesize_s(s) = poresize_calc(fiberdiameter_s(2,s),vfracp_s(1,s) + vfracp_s(2,s))
    end if

    if (scaffoldporesize_s(s) < 5D-7) then
        scaffoldporesize_s(s) = 5D-7
    end if

    rhop_s(1,s) = rhop_s(1,1)*capQrho_s(1,s)
    rhop_s(2,s) = rhop_s(2,1)*capQrho_s(2,s)

    !Need to update polymer areal masses
    capMk_s(1,s) = vfracp_s(1,s)*rhop_s(1,s)*hl_s(s - 1)
    capMk_s(2,s) = vfracp_s(2,s)*rhop_s(2,s)*hl_s(s - 1)

end subroutine update_polymer

pure function sigmoidal(t,t0,rate,off_set)

    use globalvars

    implicit none

    double precision, intent(in) :: t, t0, rate, off_set
    double precision :: sigmoidal, init
    init = (1.0/(1.0 + exp(rate*(t0 - off_set))))
    sigmoidal = (1.0/(1.0 + exp(rate*(t - off_set))))/init

end function sigmoidal

pure function poresize_calc(fiberdiameter, vfrac)

    use globalvars

    implicit none

    double precision, intent(in) :: fiberdiameter, vfrac
    double precision :: poresize_calc

    !Function to calculate pore size of construct

    poresize_calc = -sqrt(pi)/4*(1 + pi/(2*log(1 - vfrac)))*fiberdiameter

    if (poresize_calc > maxpore .or. poresize_calc /= poresize_calc) then
        poresize_calc = maxpore
    end if

end function poresize_calc

end module kinetics
