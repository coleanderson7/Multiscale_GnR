module experiments

implicit none
contains

subroutine run_experiments()

    use globalvars
    use stress
    implicit none

    double precision :: t, dF, temp, R0, hnew
    integer :: i, j, strs_flag, ldng_flag

    t = dt*(s - 1)
    strs_switch = 1

    !write(16,*) s, R0, Ph

    do i = 1, no_tests
        Ph = Pexp(i)
        call equil_gr()
        hl_s(s) = findh(rmidl_s(s))
        ril_s(s) = rmidl_s(s) - hl_s(s)/2
        if (i .eq. 1) then
            R0 = ril_s(s)
        end if
        !hnew = sqrt((R0 + hl_s(s))**2 - R0**2 + ril_s(s)**2) - ril_s(s)
    end do

    write(16,*) s, ril_s(s), R0

    RComp(s) = (ril_s(s)) - (R0)

    Ph = 2.0*133.33
    call equil_gr()
    hl_s(s) = findh(rmidl_s(s))
    ril_s(s) = rmidl_s(s) - hl_s(s)/2
    strs_switch = 0

    !Find in vivo loading quantities
    strs_flag = 1
    ldng_flag = 0
    temp = 0
    call loaded(temp,rmidl_s(s),dF,strs_flag,ldng_flag)

end subroutine run_experiments

subroutine find_unloaded_config()

    !Routine to find unloaded configuration of the construct
    !Adapted from Biaxial data analysis code from Matt Bersi

    use globalvars
    use stress
    implicit none

    double precision :: rulg, tol, obj, temp, m21, norm
    double precision :: dummy1, dummy2, dummy3, dX_1, dX_2
    double precision, dimension(2) :: X, F, X1, F1, X_new
    double precision, dimension(2,2) :: J, J0
    integer :: nmax, n, strs_flag, ldng_flag

    dummy1 = 0
    dummy2 = 0
    dummy3 = 0

    rulg = Rmid
    if (s > 1) then
        rulg = Rmid_s(s - 1)
        lamdaul_s(s) = lamdaul_s(s - 1)
    end if

    n = 1
    nmax = 100
    tol = 1D-10
    obj = 1

    strs_flag = 1
    ldng_flag = 1

    J0 = reshape([1, 0, 0, 1], shape(J0))
    J = J0

    X = [rulg, lamdaul_s(s)]

    do while (n < nmax .and. obj > tol)

        !Find Pressure and Force with the given guess
        lamdaul_s(s) = X(2)
        !Update thickness
        H_s(s) = sqrt(((ril_s(s) + hl_s(s))**2 - ril_s(s)**2)*1/X(2) + X(1)**2) - X(1)
        call loaded(dummy1, X(1) + H_s(s)/2, dummy2, strs_flag, ldng_flag)

        !Get global store local
        F(1) = Pcalc(s)
        F(2) = Fcalc(s)

        !Find the Jacobian matrix
        J = J0 !Reinitialize

        !Specify increments
        dX_1 = 1D-8*dabs(X(1)) + 1D-8
        dX_2 = 1D-8*dabs(X(2)) + 1D-8

        !Find Numerical "Partial Derivative" wrt R
        X1(1) = X(1) + dX_1
        X1(2) = X(2)

        !Update guess
        lamdaul_s(s) = X1(2)
        H_s(s) = sqrt(((ril_s(s) + hl_s(s))**2 - ril_s(s)**2)*1/X1(2) + X1(1)**2) - X1(1)
        call loaded(dummy1, X1(1) + H_s(s)/2, dummy2, strs_flag, ldng_flag)
        F1(1) = Pcalc(s)
        F1(2) = Fcalc(s)

        !Update Jacobian
        J(1,:) = (F1 - F)/dX_1

        !Find Numerical "Partial Derivative wrt Z
        X1(1) = X(1)
        X1(2) = X(2) + dX_2

        !Update guess
        lamdaul_s(s) = X1(2)
        H_s(s) = sqrt(((ril_s(s) + hl_s(s))**2 - ril_s(s)**2)*1/X1(2) + X1(1)**2) - X1(1)
        call loaded(dummy1, X1(1) + H_s(s)/2, dummy2, strs_flag, ldng_flag)
        F1(1) = Pcalc(s)
        F1(2) = Fcalc(s)

        !Update Jacobian
        J(2,:) = (F1 - F)/dX_2

        !Gaussian Elimination /w partial pivoting
        if (dabs(J(1,2)) > dabs(J(1,1))) then

            !Swap the rows of J
            temp = J(1,1)
            J(1,1) = J(1,2)
            J(1,2) = temp
            temp = J(2,1)
            J(2,1) = J(2,2)
            J(2,2) = temp

            !Swap elements of F
            temp = F(1)
            F(1) = F(2)
            F(2) = temp

        end if

        !write(*,*) J(1,1), J(1,2), J(2,1), J(2,2)

        !Mulitplier and row multiplications
        m21 = J(1,2)/J(1,1)
        J(1,2) = 0
        J(2,2) = J(2,2) - m21*J(2,1)
        F(2) = F(2) - m21*F(1)

        !Backward Substituiotn
        X_new(2) = -F(2)/J(2,2)
        X_new(1) = (-F(1) - J(2,1)*X_new(2))/J(1,1)

        !New solution
        X_new(1) = X(1) + X_new(1)
        X_new(2) = X(2) + X_new(2)

        !Evaluate convergence
        obj = dabs(X_new(1) - X(1))
        if (dabs(X_new(2) - X(2)) > obj) then
            obj = dabs(X_new(2) - X(2))
        end if
        norm = dabs(X(1))
        if (dabs(X(2)) > norm) then
            norm = dabs(X(2))
        end if

        obj = obj/norm
        if (obj < tol) then
            exit
        else
            X = X_new
        end if

        n = n + 1

    end do

    !write(*,*) s, n, F(1), F(2)

    Ri_s(s) = X_new(1)
    lamdaul_s(s) = X_new(2)
    H_s(s) = sqrt(((ril_s(s) + hl_s(s))**2 - ril_s(s)**2)*1/lamdaul_s(s) + Ri_s(s)**2) - Ri_s(s)
    Rmid_s(s) = Ri_s(s) + H_s(s)/2

end subroutine find_unloaded_config

subroutine find_linearized_stiffness()

    use globalvars
    use stress
    implicit none

    !Calculate linearized stiffness in the axial
    !and circumferential directions
    !Adapted from Baek's small on large paper (2006)

    integer :: strs_flag, ldng_flag, ks, ts, dir
    double precision :: lt, lz, sum_dWkdl2, sum_ddWkddl2, lkntau_s, t, tau
    double precision :: tex, Ctttt, Czzzz, tcalc, rl_temp, Rl, modc2, modc3, Q1, Q2

    !Deformation from unloaded to loaded configuration
    lt = rmidl_s(s)/Rmid_s(s)
    lz = 1/lamdaul_s(s)

    !Update current constituent deformations
    do ks = 3,k
        lk2_s(ks,s) = sqrt(lz**2*cos(alphak_tau(ks,s))**2 + &
                     lt**2*sin(alphak_tau(ks,s))**2)
    end do

    !Circumferential direction
    dir = 2
    tex = 0
    !Find extra stress in polymers
    !Find extra stress in ECM constituents
    sum_dWkdl2 = 0
    sum_ddWkddl2 = 0

    do ks = 1,2
        sum_dWkdl2 = sum_dWkdl2 + c2(ks)/2*poldeg(ks)*capMk_s(ks,s)/rho_s(s)
        sum_ddWkddl2 = sum_ddWkddl2 + 0
    end do

    t = (s - 1)*dt
    if (s > 1) then
            do ks = 3,k
                do ts = 1,s
                    tau = (ts - 1)*dt
                    modc2 = xlink(t,tau)
                    modc3 = xlink(t,tau)
                    lkntau_s = Gkh(ks)*(lk2_s(ks,s)/lk2_s(ks,ts))

                    if (lkntau_s > 1.00) then
                        Q1 = (lkntau_s**2 - 1)
                        Q2 = modc3*c3(ks)*Q1**2
                        sum_dWkdl2 = sum_dWkdl2 + &
                                     modc2*c2(ks)/2*Q1*exp(Q2)*&
                                     sin(alphak_tau(ks,ts))**2*&
                                     mk_tau(ks,ts)*qk_tau(ks,ts)/rho_s(s)*dt
                        sum_ddWkddl2 = sum_ddWkddl2 + &
                                       modc2*c2(ks)/2*(1+2*Q2)*exp(Q2)*&
                                       sin(alphak_tau(ks,ts))**4*&
                                       mk_tau(ks,ts)*qk_tau(ks,ts)/rho_s(s)*dt
                    else
                        Q1 = (lkntau_s**2 - 1)
                        Q2 = modc3*cc3(ks)*Q1**2
                        sum_dWkdl2 = sum_dWkdl2 + &
                                     modc2*cc2(ks)/2*Q1*exp(Q2)*&
                                     sin(alphak_tau(ks,ts))**2*&
                                     mk_tau(ks,ts)*qk_tau(ks,ts)/rho_s(s)*dt
                        sum_ddWkddl2 = sum_ddWkddl2 + &
                                       modc2*cc2(ks)/2*(1+2*Q2)*exp(Q2)*&
                                       sin(alphak_tau(ks,ts))**4*&
                                       mk_tau(ks,ts)*qk_tau(ks,ts)/rho_s(s)*dt
                    end if

                end do
            end do
    end if

    tex = 2*lt**2*sum_dWkdl2/hl_s(s)
    Ctttt = 2*tex + 4*lt**4*sum_ddWkddl2/hl_s(s)

    !Axial direction
    !Circumferential direction
    dir = 3
    tex = 0
    !Find extra stress in polymers
    !Find extra stress in ECM constituents
    sum_dWkdl2 = 0
    sum_ddWkddl2 = 0

    do ks = 1,2
        sum_dWkdl2 = sum_dWkdl2 + c2(ks)/2*poldeg(ks)*capMk_s(ks,s)/rho_s(s)
        sum_ddWkddl2 = sum_ddWkddl2 + 0
    end do

    t = (s - 1)*dt
    if (s > 1) then
            do ks = 3,k
                do ts = 1,s
                    tau = (ts - 1)*dt
                    modc2 = xlink(t,tau)
                    modc3 = xlink(t,tau)
                    lkntau_s = Gkh(ks)*(lk2_s(ks,s)/lk2_s(ks,ts))

                    if (lkntau_s > 1.00) then
                        Q1 = (lkntau_s**2 - 1)
                        Q2 = modc3*c3(ks)*Q1**2
                        sum_dWkdl2 = sum_dWkdl2 + &
                                     modc2*c2(ks)/2*Q1*exp(Q2)*&
                                     cos(alphak_tau(ks,ts))**2*&
                                     mk_tau(ks,ts)*qk_tau(ks,ts)/rho_s(s)*dt
                        sum_ddWkddl2 = sum_ddWkddl2 + &
                                       modc2*c2(ks)/2*(1+2*Q2)*exp(Q2)*&
                                       cos(alphak_tau(ks,ts))**4*&
                                       mk_tau(ks,ts)*qk_tau(ks,ts)/rho_s(s)*dt
                    else
                        Q1 = (lkntau_s**2 - 1)
                        Q2 = modc3*cc3(ks)*Q1**2
                        sum_dWkdl2 = sum_dWkdl2 + &
                                     modc2*cc2(ks)/2*Q1*exp(Q2)*&
                                     cos(alphak_tau(ks,ts))**2*&
                                     mk_tau(ks,ts)*qk_tau(ks,ts)/rho_s(s)*dt
                        sum_ddWkddl2 = sum_ddWkddl2 + &
                                       modc2*cc2(ks)/2*(1+2*Q2)*exp(Q2)*&
                                       cos(alphak_tau(ks,ts))**4*&
                                       mk_tau(ks,ts)*qk_tau(ks,ts)/rho_s(s)*dt
                    end if

                end do
            end do
    end if

    tex = 2*lz**2*sum_dWkdl2/hl_s(s)
    Czzzz = 2*tex + 4*lz**4*sum_ddWkddl2/hl_s(s)

    Stiffness(1,s) = Ctttt
    Stiffness(2,s) = Czzzz

end subroutine find_linearized_stiffness

end module experiments

