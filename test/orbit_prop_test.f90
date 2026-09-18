!------------------------------------------------------------------------------
!>
!  Numerical orbit propagation test and visualization for Brouwer elements.
!  Integrates an Earth orbit under zonal harmonics (J2, J3, J4, J5) using DDEABM,
!  computes osculating and Brouwer mean elements along the trajectory,
!  and generates comparative plots with pyplot-fortran.

program orbit_prop_test

    use iso_fortran_env, only: error_unit, int64
    use brouwer_module, wp => brouwer_module_wp
    use ddeabm_module, only: ddeabm_class
    use pyplot_module, only: pyplot

    implicit none

    ! Gravitational and planetary parameters (Earth)
    real(wp), parameter :: mu_earth  = 398600.4415_wp           ! km^3/s^2
    real(wp), parameter :: req_earth = 6378.1363_wp             ! km
    real(wp), parameter :: j2_earth  = 1.082626925638815e-3_wp
    real(wp), parameter :: j3_earth  = -0.2532307818191774e-5_wp
    real(wp), parameter :: j4_earth  = -0.1620429990000000e-5_wp
    real(wp), parameter :: j5_earth  = -0.2270711043920343e-6_wp

    ! Orbit simulation parameters: LEO orbit, ~2 days (approx 30 orbits)
    integer, parameter :: n_steps = 1000
    real(wp), parameter :: t_final = 0.1_wp * 86400.0_wp        ! seconds
    ! real(wp), parameter :: t_final = 1.0_wp * 86400.0_wp        ! seconds
    ! real(wp), parameter :: t_final = 80.0_wp * 86400.0_wp        ! seconds
    real(wp), parameter :: dt = t_final / real(n_steps, wp)

    ! Initial orbital elements: [sma (km), ecc, inc (deg), raan (deg), aop (deg), ta (deg)]
    real(wp), dimension(6) :: kep_0, cart_0, cart_state, kep_osc, bl_short, bl_long
    real(wp), dimension(6) :: cart_prop, kep_prop
    real(wp), dimension(n_steps + 1) :: t_hrs
    real(wp), dimension(n_steps + 1) :: sma_osc, sma_short, sma_long, sma_prop
    real(wp), dimension(n_steps + 1) :: ecc_osc, ecc_short, ecc_long, ecc_prop
    real(wp), dimension(n_steps + 1) :: inc_osc, inc_short, inc_long, inc_prop
    real(wp), dimension(n_steps + 1) :: aop_osc, aop_short, aop_long, aop_prop
    real(wp), dimension(n_steps + 1) :: raan_osc, raan_short, raan_long, raan_prop
    real(wp), dimension(n_steps + 1) :: ma_osc, ma_short, ma_long, ma_prop

    type(ddeabm_class) :: solver
    type(pyplot) :: plt
    real(wp) :: t, t_out
    integer :: i, stat, idid
    character(len=100) :: file_suffix
    real(wp), dimension(3) :: acc1, acc2

    ! colors:
    real(wp),dimension(3),parameter :: c0 = [0.0_wp, 0.4470_wp, 0.7410_wp]
    real(wp),dimension(3),parameter :: c1 = [0.8500_wp, 0.3250_wp, 0.0980_wp]
    real(wp),dimension(3),parameter :: c2 = [0.9290_wp, 0.6940_wp, 0.1250_wp]
    real(wp),dimension(3),parameter :: c3 = [0.4940_wp, 0.1840_wp, 0.5560_wp]

    integer,dimension(2),parameter :: figsize = [10,5]

    print *, "=========================================================="
    print *, " Orbit Propagation & Brouwer Mean Element Plotting Test"
    print *, "=========================================================="

    write(file_suffix, '(i10)') int(t_final / 3600.0_wp) ! hrs
    file_suffix = '_TF='//trim(adjustl(file_suffix))//'h'

    ! 1. Initialize Orbit
    kep_0 = [6800.0_wp, 0.02_wp, 51.6_wp, 30.0_wp, 40.0_wp, 0.0_wp] ! ISS-like LEO orbit
    ! kep_0 = [7500.0_wp, 0.02_wp, 63.4349488_wp, 50.0_wp, 40.0_wp, 80.0_wp] ! critical inclination
    ! kep_0 = [6800.0_wp, 0.02_wp, 33.3_wp, 30.0_wp, 40.0_wp, 0.0_wp] ! lower inc case
    call keplerian_to_cartesian(mu_earth, kep_0, anomaly_type="TA", stat=stat, cart=cart_0)
    if (stat /= BROUWER_SUCCESS) error stop "Error converting initial Keplerian to Cartesian state."

    cart_state = cart_0
    t = 0.0_wp

    ! Record initial step (t = 0)
    t_hrs(1) = 0.0_wp
    call cartesian_to_keplerian(mu_earth, cart_state, anomaly_type="MA", kepl=kep_osc, stat=stat)
    call cartesian_to_brouwer_mean_short(mu_earth, req_earth, j2_earth, cart_state, stat=stat, blms=bl_short)
    call cartesian_to_brouwer_mean_long(mu_earth, req_earth, j2_earth, j3_earth, j4_earth, j5_earth, cart_state, stat=stat, blml=bl_long)
    call brouwer_lyddane_propagate(mu_earth, req_earth, j2_earth, j3_earth, j4_earth, j5_earth, cart_0, 0.0_wp, stat, cart_prop)
    call cartesian_to_keplerian(mu_earth, cart_prop, anomaly_type="MA", kepl=kep_prop, stat=stat)

    sma_osc(1) = kep_osc(1);  sma_short(1) = bl_short(1);  sma_long(1) = bl_long(1);  sma_prop(1) = kep_prop(1)
    ecc_osc(1) = kep_osc(2);  ecc_short(1) = bl_short(2);  ecc_long(1) = bl_long(2);  ecc_prop(1) = kep_prop(2)
    inc_osc(1) = kep_osc(3);  inc_short(1) = bl_short(3);  inc_long(1) = bl_long(3);  inc_prop(1) = kep_prop(3)
    raan_osc(1) = kep_osc(4); raan_short(1) = bl_short(4); raan_long(1) = bl_long(4); raan_prop(1) = kep_prop(4)
    aop_osc(1) = kep_osc(5);  aop_short(1) = bl_short(5);  aop_long(1) = bl_long(5);  aop_prop(1) = kep_prop(5)
    ma_osc(1) = kep_osc(6);  ma_short(1) = bl_short(6);  ma_long(1) = bl_long(6);  ma_prop(1) = kep_prop(6)

    ! 2. Initialize Integrator
    call solver%initialize(6,maxnum=10000000,df=grav_derivs,rtol=[1.0e-12_wp],atol=[1.0e-12_wp])

    print *, "Propagating orbit with J2, J3, J4, J5 gravity perturbation..."
    do i = 1, n_steps
        t_out = real(i, wp) * dt
        call solver%integrate(t, cart_state, t_out, idid=idid)
        if (idid < 1) then
            print*, 'step = ', i, 'idid = ', idid
            error stop "Integrator error"
        end if

        t_hrs(i + 1) = t / 3600.0_wp
        call cartesian_to_keplerian(mu_earth, cart_state, anomaly_type="MA", kepl=kep_osc, stat=stat)
        call cartesian_to_brouwer_mean_short(mu_earth, req_earth, j2_earth, cart_state, stat=stat, blms=bl_short)
        call cartesian_to_brouwer_mean_long(mu_earth, req_earth, j2_earth, j3_earth, j4_earth, j5_earth, cart_state, stat=stat, blml=bl_long)
        call brouwer_lyddane_propagate(mu_earth, req_earth, j2_earth, j3_earth, j4_earth, j5_earth, cart_0, t_out, stat, cart_prop)
        call cartesian_to_keplerian(mu_earth, cart_prop, anomaly_type="MA", kepl=kep_prop, stat=stat)

        sma_osc(i + 1) = kep_osc(1);  sma_short(i + 1) = bl_short(1);  sma_long(i + 1) = bl_long(1);  sma_prop(i + 1) = kep_prop(1)
        ecc_osc(i + 1) = kep_osc(2);  ecc_short(i + 1) = bl_short(2);  ecc_long(i + 1) = bl_long(2);  ecc_prop(i + 1) = kep_prop(2)
        inc_osc(i + 1) = kep_osc(3);  inc_short(i + 1) = bl_short(3);  inc_long(i + 1) = bl_long(3);  inc_prop(i + 1) = kep_prop(3)
        raan_osc(i + 1) = kep_osc(4); raan_short(i + 1) = bl_short(4); raan_long(i + 1) = bl_long(4); raan_prop(i + 1) = kep_prop(4)
        aop_osc(i + 1) = kep_osc(5);  aop_short(i + 1) = bl_short(5);  aop_long(i + 1) = bl_long(5);  aop_prop(i + 1) = kep_prop(5)
        ma_osc(i + 1) = kep_osc(6);  ma_short(i + 1) = bl_short(6);  ma_long(i + 1) = bl_long(6);  ma_prop(i + 1) = kep_prop(6)
    end do

    print *, "Propagation completed successfully."

    ! 3. Speed Comparison: DDEABM Numerical Integration (grav_derivs force
    !    model) vs Brouwer-Lyddane Analytic Propagation.
    !
    !    Each method is used to answer the same set of "state at time t"
    !    queries from a fixed epoch. The analytic side converts the epoch
    !    Cartesian state to Brouwer mean elements ONCE (that osculating ->
    !    mean conversion is an iterative fixed-point solve, done once per
    !    epoch in any real use case) and then reuses the O(1), non-iterative
    !    propagate_brouwer_mean_long / brouwer_mean_long_to_osculating /
    !    keplerian_to_cartesian chain for every query - the actual cost of
    !    "propagate to a new time" once you're set up at an epoch.
    block
        type(ddeabm_class) :: speed_solver
        real(wp), dimension(6) :: cart_speed_state, blml_speed0, blml_speedt, kepl_speed, cart_speed_prop
        real(wp) :: t_speed, elapsed_numerical, elapsed_analytic, speedup
        integer(int64) :: clock_rate, clock_start, clock_end
        integer :: j, n_queries, speed_idid, speed_stat

        print *, "=========================================================="
        print *, " Speed Comparison: DDEABM Integration vs Brouwer-Lyddane Propagation"
        print *, "=========================================================="

        n_queries = n_steps * 1

        call system_clock(count_rate=clock_rate)

        call speed_solver%initialize(6, maxnum=10000000, df=grav_derivs, rtol=[1.0e-12_wp], atol=[1.0e-12_wp])

        elapsed_numerical = 0.0_wp
        t_out = real(n_queries, wp) * dt
        call system_clock(clock_start)
        do j = 1, n_queries
            cart_speed_state = cart_0
            t_speed = 0.0_wp
            call speed_solver%first_call()
            call speed_solver%integrate(t_speed, cart_speed_state, t_out, idid=speed_idid)
            if (speed_idid < 1) error stop "Integrator error in speed test"
        end do
        call system_clock(clock_end)
        elapsed_numerical = elapsed_numerical + real(clock_end - clock_start, wp) / real(clock_rate, wp)

        elapsed_analytic = 0.0_wp
        t_out = real(n_queries, wp) * dt
        call system_clock(clock_start)
        do j = 1, n_queries
            call cartesian_to_brouwer_mean_long(mu_earth, req_earth, j2_earth, j3_earth, j4_earth, j5_earth, &
                                                cart_0, stat=speed_stat, blml=blml_speed0)
            call propagate_brouwer_mean_long(mu_earth, req_earth, j2_earth, j4_earth, blml_speed0, t_out, speed_stat, blml_speedt)
            call brouwer_mean_long_to_osculating(mu_earth, req_earth, j2_earth, j3_earth, j4_earth, j5_earth, &
                                                  blml_speedt, speed_stat, kepl_speed)
            call keplerian_to_cartesian(mu_earth, kepl_speed, anomaly_type="MA", stat=speed_stat, cart=cart_speed_prop)
        end do
        call system_clock(clock_end)
        elapsed_analytic = elapsed_analytic + real(clock_end - clock_start, wp) / real(clock_rate, wp)

        speedup = elapsed_numerical / elapsed_analytic

        print '(A,I0)',      " Number of propagation queries       : ", n_queries
        print '(A,ES12.5,A)', " DDEABM (grav_derivs) total time     : ", elapsed_numerical, " s"
        print '(A,ES12.5,A)', " Brouwer-Lyddane propagate total time: ", elapsed_analytic, " s"
        print '(A,F12.1,A)',  " Speedup (analytic vs numerical)     : ", speedup, "x"
        if (speedup > 1.0_wp) then
            print *, "Brouwer-Lyddane propagation was faster for this workload."
        else
            print *, "DDEABM propagation was faster for this workload."
        end if

        print *, "=========================================================="
    end block

    ! 4. Generate Comparative Plots using pyplot-fortran
    print *, "Generating element plots..."

    ! Plot 1: Semi-major Axis Comparison
    call plt%initialize(title="Semi-Major Axis: Osculating vs Brouwer Mean", &
                        figsize=figsize, &
                        xlabel="Time (hours)", &
                        ylabel="Semi-major axis $a$ (km)", &
                        legend = .true.)
    call plt%add_plot(t_hrs, sma_osc,   label="Osculating", linestyle="-", color=c0, linewidth=1)
    call plt%add_plot(t_hrs, sma_short, label="Brouwer Short-Period Mean", linestyle="--", color=c1)
    call plt%add_plot(t_hrs, sma_long,  label="Brouwer Long-Period Mean", linestyle=":", color=c2)
    call plt%add_plot(t_hrs, sma_prop,  label="Brouwer-Lyddane Propagation", linestyle="-.", color=c3)
    call plt%savefig("brouwer_sma_comparison"//trim(file_suffix)//".png")

    ! Plot 2: Eccentricity Comparison
    call plt%initialize(title="Eccentricity: Osculating vs Brouwer Mean", &
                        figsize=figsize, &
                        xlabel="Time (hours)", &
                        ylabel="Eccentricity $e$", &
                        legend = .true.)
    call plt%add_plot(t_hrs, ecc_osc,   label="Osculating", linestyle="-", color=c0, linewidth=1)
    call plt%add_plot(t_hrs, ecc_short, label="Brouwer Short-Period Mean", linestyle="--", color=c1)
    call plt%add_plot(t_hrs, ecc_long,  label="Brouwer Long-Period Mean", linestyle=":", color=c2)
    call plt%add_plot(t_hrs, ecc_prop,  label="Brouwer-Lyddane Propagation", linestyle="-.", color=c3)
    call plt%savefig("brouwer_ecc_comparison"//trim(file_suffix)//".png" )

    ! Plot 3: Inclination Comparison
    call plt%initialize(title="Inclination: Osculating vs Brouwer Mean", &
                        figsize=figsize, &
                        xlabel="Time (hours)", &
                        ylabel="Inclination $i$ (deg)", &
                        legend = .true.)
    call plt%add_plot(t_hrs, inc_osc,   label="Osculating", linestyle="-", color=c0, linewidth=1)
    call plt%add_plot(t_hrs, inc_short, label="Brouwer Short-Period Mean", linestyle="--", color=c1)
    call plt%add_plot(t_hrs, inc_long,  label="Brouwer Long-Period Mean", linestyle=":", color=c2)
    call plt%add_plot(t_hrs, inc_prop,  label="Brouwer-Lyddane Propagation", linestyle="-.", color=c3)
    call plt%savefig("brouwer_inc_comparison"//trim(file_suffix)//".png")

    ! Plot 4: Argument of Periapsis Comparison
    call plt%initialize(title="Argument of Periapsis: Osculating vs Brouwer Mean", &
                        figsize=figsize, &
                        xlabel="Time (hours)", &
                        ylabel="Argument of Periapsis $\\omega$ (deg)", &
                        legend = .true.)
    call plt%add_plot(t_hrs, aop_osc,   label="Osculating", linestyle="-", color=c0, linewidth=1)
    call plt%add_plot(t_hrs, aop_short, label="Brouwer Short-Period Mean", linestyle="--", color=c1)
    call plt%add_plot(t_hrs, aop_long,  label="Brouwer Long-Period Mean", linestyle=":", color=c2)
    call plt%add_plot(t_hrs, aop_prop,  label="Brouwer-Lyddane Propagation", linestyle="-.", color=c3)
    call plt%savefig("brouwer_aop_comparison"//trim(file_suffix)//".png")

    ! Plot 5: RAAN Comparison
    call plt%initialize(title="RAAN: Osculating vs Brouwer Mean", &
                        figsize=figsize, &
                        xlabel="Time (hours)", &
                        ylabel="RAAN $\\Omega$ (deg)", &
                        legend = .true.)
    call plt%add_plot(t_hrs, raan_osc,   label="Osculating", linestyle="-", color=c0, linewidth=1)
    call plt%add_plot(t_hrs, raan_short, label="Brouwer Short-Period Mean", linestyle="--", color=c1)
    call plt%add_plot(t_hrs, raan_long,  label="Brouwer Long-Period Mean", linestyle=":", color=c2)
    call plt%add_plot(t_hrs, raan_prop,  label="Brouwer-Lyddane Propagation", linestyle="-.", color=c3)
    call plt%savefig("brouwer_raan_comparison"//trim(file_suffix)//".png")

    ! Plot 6: Mean Anomaly Comparison
    call plt%initialize(title="Mean Anomaly: Osculating vs Brouwer Mean", &
                        figsize=figsize, &
                        xlabel="Time (hours)", &
                        ylabel="Mean Anomaly $M$ (deg)", &
                        legend = .true.)
    call plt%add_plot(t_hrs, ma_osc,   label="Osculating", linestyle="-", color=c0, linewidth=1)
    call plt%add_plot(t_hrs, ma_short, label="Brouwer Short-Period Mean", linestyle="--", color=c1)
    call plt%add_plot(t_hrs, ma_long,  label="Brouwer Long-Period Mean", linestyle=":", color=c2)
    call plt%add_plot(t_hrs, ma_prop,  label="Brouwer-Lyddane Propagation", linestyle="-.", color=c3)
    call plt%savefig("brouwer_ma_comparison"//trim(file_suffix)//".png")

    print *, "All plots generated successfully: brouwer_*.png"
    print *, "=========================================================="

    print*, ''
    print *, "=========================================================="
    print *, "Gravity unit test: comparing gravity_j2_j3_j4_j5 vs gravity_j2_j3_j4"
    print *, "=========================================================="

    ! test of the two grav routines:
    print*, 'general test of gravity_j2_j3_j4_j5 vs gravity_j2_j3_j4'
    call gravity_j2_j3_j4_j5(cart_0(1:3),mu_earth,req_earth,j2_earth,j3_earth,j4_earth,j5_earth,acc1)
    call gravity_j2_j3_j4(   cart_0(1:3),mu_earth,req_earth,j2_earth,j3_earth,j4_earth,         acc2)
    print*, 'gravity_j2_j3_j4_j5 acc = ', acc1
    print*, 'gravity_j2_j3_j4 acc    = ', acc2
    if (any(abs(acc1 - acc2) > 1.0e-5_wp)) then ! should be close
        error stop 'Error: gravity_j2_j3_j4_j5 and gravity_j2_j3_j4 not close!'
    else
        print*, 'Success: gravity_j2_j3_j4_j5 and gravity_j2_j3_j4 are close.'
    end if

    ! if j5=0, then these should return the same acceleration vector. Let's test that:
    print*, 'test when j5=0'
    call gravity_j2_j3_j4_j5(cart_0(1:3),mu_earth,req_earth,j2_earth,j3_earth,j4_earth,0.0_wp,acc1)
    call gravity_j2_j3_j4(   cart_0(1:3),mu_earth,req_earth,j2_earth,j3_earth,j4_earth,       acc2)
    print*, 'gravity_j2_j3_j4_j5 acc = ', acc1
    print*, 'gravity_j2_j3_j4 acc    = ', acc2
    if (any(abs(acc1 - acc2) > 1.0e-12_wp)) then
        error stop 'Error: gravity_j2_j3_j4_j5 and gravity_j2_j3_j4 do not match when j5=0!'
    else
        print*, 'Success: gravity_j2_j3_j4_j5 and gravity_j2_j3_j4 match when j5=0.'
    end if
    print *, "=========================================================="

contains

    subroutine grav_derivs(me, t, y, dydt)
        !! Right-hand-side equations of motion with zonal harmonics J2, J3, J4, J5
        class(ddeabm_class),intent(inout) :: me
        real(wp), intent(in) :: t !! time [s] - not used here
        real(wp), dimension(:), intent(in) :: y      !! [r,v]
        real(wp), dimension(:), intent(out) :: dydt  !! [v,a]

        real(wp),dimension(3) :: acc

        ! Kinematics: dr/dt = v
        dydt(1:3) = y(4:6)

        ! Dynamics: dv/dt = a
        call gravity_j2_j3_j4_j5(y(1:3),mu_earth,req_earth,j2_earth,j3_earth,j4_earth,j5_earth,acc)
        ! call gravity_j2_j3_j4(y(1:3),mu_earth,req_earth,j2_earth,j3_earth,j4_earth,acc)
        dydt(4:6) = acc(1:3)

    end subroutine grav_derivs

!*****************************************************************************************
!>
!  Gravitational acceleration due to simplified spherical harmonic
!  expansion (only the J2-J5 terms are used).
!
!@note This is an AI-generated routine.

    subroutine gravity_j2_j3_j4_j5(r,mu,req,j2,j3,j4,j5,acc)

    real(wp),dimension(3),intent(in)  :: r   !! satellite position vector [km]
    real(wp),intent(in)               :: mu  !! central body gravitational parameter [km^3/s^2]
    real(wp),intent(in)               :: req !! body equatorial radius [km]
    real(wp),intent(in)               :: j2  !! j2 coefficient
    real(wp),intent(in)               :: j3  !! j3 coefficient
    real(wp),intent(in)               :: j4  !! j4 coefficient
    real(wp),intent(in)               :: j5  !! j5 coefficient
    real(wp),dimension(3),intent(out) :: acc !! gravity acceleration vector [km/s^2]

    real(wp) :: rmag, r2, r3, z, z_r, z2_r2, z3_r3, z4_r4
    real(wp) :: re_r, re_r2, re_r3, re_r4, re_r5
    real(wp) :: f_r, f_z, mu_r3

    r2 = dot_product(r, r)
    rmag = sqrt(r2)
    if (rmag==0.0_wp) error stop 'Error in gravity_j2_j3_j4_j5: spacecraft at center of body.'

    r3 = rmag * r2
    z = r(3)
    z_r = z / rmag
    z2_r2 = z_r * z_r
    z3_r3 = z2_r2 * z_r
    z4_r4 = z2_r2 * z2_r2

    re_r = req / rmag
    re_r2 = re_r * re_r
    re_r3 = re_r2 * re_r
    re_r4 = re_r3 * re_r
    re_r5 = re_r4 * re_r

    mu_r3 = mu / r3

    ! Potential partial derivatives: a = - (mu/rmag^3)*r_vec + a_pert
    ! Common factors for zonal perturbation acceleration components
    f_r = 1.0_wp + 1.5_wp * j2 * re_r2 * (1.0_wp - 5.0_wp * z2_r2) &
            + 2.5_wp * j3 * re_r3 * (3.0_wp * z_r - 7.0_wp * z3_r3) &
            - 0.625_wp * j4 * re_r4 * (3.0_wp - 42.0_wp * z2_r2 + 63.0_wp * z4_r4) &
            - 2.625_wp * j5 * re_r5 * (5.0_wp * z_r - 30.0_wp * z3_r3 + 33.0_wp * z4_r4 * z_r)

    f_z = -3.0_wp * j2 * re_r2 * z_r &
        + 0.5_wp * j3 * re_r3 * (3.0_wp - 15.0_wp * z2_r2) &
        + 2.5_wp * j4 * re_r4 * (3.0_wp * z_r - 7.0_wp * z3_r3) &
        - 1.875_wp * j5 * re_r5 * (1.0_wp - 14.0_wp * z2_r2 + 21.0_wp * z4_r4)

    acc(1) = -mu_r3 * r(1) * f_r
    acc(2) = -mu_r3 * r(2) * f_r
    acc(3) = -mu_r3 * (r(3) * f_r - rmag * f_z)

    end subroutine gravity_j2_j3_j4_j5
!*****************************************************************************************

!*****************************************************************************************
!>
!  Gravitational acceleration due to simplified spherical harmonic
!  expansion (only the J2-J4 terms are used).
!
!### Reference
!  * http://www.ni.com/pdf/manuals/370762a.pdf
!
!@note This is from the Fortran Astrodynamics Toolkit.

    subroutine gravity_j2_j3_j4(r,mu,req,j2,j3,j4,acc)

    real(wp),dimension(3),intent(in)  :: r   !! satellite position vector [km]
    real(wp),intent(in)               :: mu  !! central body gravitational parameter [km^3/s^2]
    real(wp),intent(in)               :: req !! body equatorial radius [km]
    real(wp),intent(in)               :: j2  !! j2 coefficient
    real(wp),intent(in)               :: j3  !! j3 coefficient
    real(wp),intent(in)               :: j4  !! j4 coefficient
    real(wp),dimension(3),intent(out) :: acc !! gravity acceleration vector [km/s^2]

    real(wp) :: mmor3,reqor,reqor2,reqor3,reqor4,&
                rmag,rmag2,rmag3,rzor,rzor2,rzor3,rzor4,c,d

    rmag2 = dot_product(r,r)
    rmag  = sqrt(rmag2)

    if (rmag==0.0_wp) error stop 'Error in gravity_j2_j3_j4: spacecraft at center of body.'

    rmag3  = rmag*rmag2
    mmor3  = -mu/rmag3
    reqor  = req/rmag
    reqor2 = reqor*reqor
    reqor3 = reqor2*reqor
    reqor4 = reqor3*reqor
    rzor   = r(3)/rmag
    rzor2  = rzor*rzor
    rzor3  = rzor2*rzor
    rzor4  = rzor3*rzor

    c = mmor3 * (1.0_wp - 1.5_wp*J2*reqor2*(5.0_wp*rzor2-1.0_wp) + &
                    2.5_wp*J3*reqor3*(3.0_wp*rzor-7.0_wp*rzor3) - &
                    0.625_wp*J4*reqor4*(3.0_wp-42.0_wp*rzor2+63.0_wp*rzor4))
    d = mmor3 * (r(3) + 1.5_wp*J2*reqor2*(3.0_wp-5.0_wp*rzor2)*r(3) + &
                0.5_wp*J3*reqor3*(30.0_wp*rzor*r(3)-35.0_wp*rzor3*r(3)-3.0_wp*rmag) - &
                0.625_wp*J4*reqor4*(15.0_wp-70.0_wp*rzor2+63.0_wp*rzor4)*r(3))

    acc(1) = c * r(1)
    acc(2) = c * r(2)
    acc(3) = d

    end subroutine gravity_J2_J3_J4
!*****************************************************************************************

end program orbit_prop_test
