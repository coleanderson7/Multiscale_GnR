module outputs

implicit none
contains

subroutine update_outputs()

    use globalvars
    implicit none

    double precision :: t

    t = dt*(s - 1)

    !write(10,*) t, ril_s(s), hl_s(s), rho_s(s), cauchy(1,s), Ri_s(s), lamdaul_s(s)
    write(10,*) t, ril_s(s), rmidl_s(s), hl_s(s), cauchy(1,s), Ri_s(s), H_s(s), lamdaul_s(s), &
                Stiffness(1,s), Stiffness(2,s), rho_s(s), Pcalc(s), RComp(s)

    write(11,*) capMk_s(1,s), capMk_s(2,s), capMk_s(3,s), capMk_s(4,s), capMk_s(5,s), &
               capMk_s(6,s), capMk_s(7,s), capMk_s(8,s), capMk_s(9,s), capMk_s(10,s), &
               capMk_s(11,s), capMk_s(12,s), capMfiller_s(s), sum(capMk_s(:,s)), &
               (sum(capMk_s(:,s))+ capMfiller_s(s))/hl_s(s)!+ capMfiller_s(s)
    write(12,*) phik_s(1,s), phik_s(2,s), phik_s(3,s), phik_s(4,s), phik_s(5,s), &
               phik_s(6,s), phik_s(7,s), phik_s(8,s), phik_s(9,s), phik_s(10,s), &
               phik_s(11,s), phik_s(12,s), sum(phik_s(:,s)), sum(phik_s(3:6,s)), sum(phik_s(7:10,s))
    write(13,*) mk_tau(1,s), mk_tau(2,s), mk_tau(3,s), mk_tau(4,s), mk_tau(5,s), &
               mk_tau(6,s), mk_tau(7,s), mk_tau(8,s), mk_tau(9,s), mk_tau(10,s), &
               mk_tau(11,s), mk_tau(12,s)
    write(15,*) scaffoldporesize_s(s), poresize_s(1,s), poresize_s(2,s), &
               fiberdiameter_s(1,s), fiberdiameter_s(2,s), vfracp_s(1,s), vfracp_s(2,s)

end subroutine update_outputs

end module outputs
