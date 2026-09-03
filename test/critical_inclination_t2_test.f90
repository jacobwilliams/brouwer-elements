
!------------------------------------------------------------------------------
!>
!  Test of [[critical_inclination_t2]].

program critical_inclination_t2_test

    use brouwer_module, only: wp => brouwer_module_wp, critical_inclination_t2
    use pyplot_module, only: pyplot

    implicit none

    integer :: i
    real(wp) :: inc, t2, eqn, cosinc
    real(wp),dimension(:), allocatable :: theta_vec, t2_vec, eqn_vec
    type(pyplot) :: plt

    real(wp),dimension(3),parameter :: blue = [0, 0, 1]   ! blue
    real(wp),dimension(3),parameter :: red = [1, 0, 0]    ! red
    real(wp),dimension(3),parameter :: black = [0, 0, 0]  ! black
    integer,dimension(2),parameter :: figsize = [6, 4]
    real(wp),parameter :: deg2rad = acos(-1.0_wp) / 180.0_wp
    real(wp),parameter :: rad2deg = 1.0_wp / deg2rad

    allocate(theta_vec(0), t2_vec(0), eqn_vec(0))
    do i = 0, 180, 1
        inc = real(i, wp) * deg2rad
        cosinc = cos(inc)
        t2 = critical_inclination_t2(cosinc)  ! the approximation of the function
        eqn =  1.0_wp / (1.0_wp - 5.0_wp * cosinc**2)  ! the function being approximated
        theta_vec = [theta_vec, inc * rad2deg]
        t2_vec    = [t2_vec, t2]
        eqn_vec   = [eqn_vec, eqn]
        write(*,'(a, i5, a, f18.10, a, f18.10, a, f18.10)') "inc = ", int(inc * rad2deg), " deg, t2 = ", t2, ", eqn = ", eqn, ', diff = ', t2 - eqn
    end do

    ! the approximation should be very close to the actual function,
    ! except at the singularities at 63.4349 deg and 116.5651 deg
    call plt%initialize(title="critical_inclination_t2 Comparison", &
                        figsize=figsize, &
                        xlabel="Inclination (deg)", &
                        ylabel="", &
                        grid=.true., &
                        legend = .true.)
    call plt%add_plot(theta_vec, t2_vec,  label="t2", linestyle="-", color=blue, linewidth=2)
    call plt%add_plot(theta_vec, eqn_vec, label="$(1 - 5 cos^2(i))^{-1}$", linestyle=":", color=red, linewidth=2)
    call plt%add_plot([63.4349_wp, 63.4349_wp],   [-50.0_wp, 50.0_wp],  label="63.4349 deg", linestyle="--", color=black, linewidth=1)
    call plt%add_plot([116.5651_wp, 116.5651_wp], [-50.0_wp, 50.0_wp],  label="116.5651 deg", linestyle="--", color=black, linewidth=1)
    call plt%savefig("critical_inclination_t2_comparison.png")

end program critical_inclination_t2_test