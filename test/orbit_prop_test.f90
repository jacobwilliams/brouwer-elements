!------------------------------------------------------------------------------
!>
!  Numerical orbit propagation test and visualization for Brouwer elements.
!  Integrates an Earth orbit under zonal harmonics (J2, J3, J4, J5) using DDEABM,
!  computes osculating and Brouwer mean elements along the trajectory,
!  and generates comparative plots with pyplot-fortran.

program orbit_prop_test

    use iso_fortran_env, only: wp => real64, error_unit
    use brouwer_module
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
    real(wp), parameter :: t_final = 1.0_wp * 86400.0_wp        ! seconds
    real(wp), parameter :: dt = t_final / real(n_steps, wp)

    ! Initial orbital elements: [sma (km), ecc, inc (deg), raan (deg), aop (deg), ta (deg)]
    real(wp), dimension(6) :: kep_0, cart_0, cart_state, kep_osc, bl_short, bl_long
    real(wp), dimension(n_steps + 1) :: t_hrs
    real(wp), dimension(n_steps + 1) :: sma_osc, sma_short, sma_long
    real(wp), dimension(n_steps + 1) :: ecc_osc, ecc_short, ecc_long
    real(wp), dimension(n_steps + 1) :: inc_osc, inc_short, inc_long
    real(wp), dimension(n_steps + 1) :: aop_osc, aop_short, aop_long
    real(wp), dimension(n_steps + 1) :: raan_osc, raan_short, raan_long

    type(ddeabm_class) :: solver
    type(pyplot) :: plt
    real(wp) :: t, t_out
    integer :: i, stat, idid

    ! colors:
    real(wp),dimension(3) :: c0 = [0.0_wp, 0.4470_wp, 0.7410_wp]
    real(wp),dimension(3) :: c1 = [0.8500_wp, 0.3250_wp, 0.0980_wp]
    real(wp),dimension(3) :: c2 = [0.9290_wp, 0.6940_wp, 0.1250_wp]

    print *, "=========================================================="
    print *, " Orbit Propagation & Brouwer Mean Element Plotting Test"
    print *, "=========================================================="

    ! 1. Initialize Orbit (e.g. ISS-like LEO orbit)
    kep_0 = [6800.0_wp, 0.02_wp, 51.6_wp, 30.0_wp, 40.0_wp, 0.0_wp]
    cart_0 = keplerian_to_cartesian(mu_earth, kep_0, anomaly_type="TA", stat=stat)
    if (stat /= 0) error stop "Error converting initial Keplerian to Cartesian state."

    cart_state = cart_0
    t = 0.0_wp

    ! Record initial step (t = 0)
    t_hrs(1) = 0.0_wp
    kep_osc = cartesian_to_keplerian(mu_earth, cart_state, anomaly_type="MA")
    bl_short = cartesian_to_brouwer_mean_short(mu_earth, req_earth, j2_earth, cart_state)
    bl_long  = cartesian_to_brouwer_mean_long(mu_earth, req_earth, j2_earth, j3_earth, j4_earth, j5_earth, cart_state)

    sma_osc(1) = kep_osc(1);  sma_short(1) = bl_short(1);  sma_long(1) = bl_long(1)
    ecc_osc(1) = kep_osc(2);  ecc_short(1) = bl_short(2);  ecc_long(1) = bl_long(2)
    inc_osc(1) = kep_osc(3);  inc_short(1) = bl_short(3);  inc_long(1) = bl_long(3)
    raan_osc(1) = kep_osc(4); raan_short(1) = bl_short(4); raan_long(1) = bl_long(4)
    aop_osc(1) = kep_osc(5);  aop_short(1) = bl_short(5);  aop_long(1) = bl_long(5)

    ! 2. Initialize Integrator
    call solver%initialize(6,maxnum=10000,df=grav_derivs,rtol=[1.0e-12_wp],atol=[1.0e-12_wp])

    print *, "Propagating orbit with J2, J3, J4, J5 gravity perturbation..."
    do i = 1, n_steps
        t_out = real(i, wp) * dt
        call solver%integrate(t, cart_state, t_out, idid=idid)
        if (idid < 1) error stop "Integrator error"

        t_hrs(i + 1) = t / 3600.0_wp
        kep_osc  = cartesian_to_keplerian(mu_earth, cart_state, anomaly_type="MA")
        bl_short = cartesian_to_brouwer_mean_short(mu_earth, req_earth, j2_earth, cart_state)
        bl_long  = cartesian_to_brouwer_mean_long(mu_earth, req_earth, j2_earth, j3_earth, j4_earth, j5_earth, cart_state)

        sma_osc(i + 1) = kep_osc(1);  sma_short(i + 1) = bl_short(1);  sma_long(i + 1) = bl_long(1)
        ecc_osc(i + 1) = kep_osc(2);  ecc_short(i + 1) = bl_short(2);  ecc_long(i + 1) = bl_long(2)
        inc_osc(i + 1) = kep_osc(3);  inc_short(i + 1) = bl_short(3);  inc_long(i + 1) = bl_long(3)
        raan_osc(i + 1) = kep_osc(4); raan_short(i + 1) = bl_short(4); raan_long(i + 1) = bl_long(4)
        aop_osc(i + 1) = kep_osc(5);  aop_short(i + 1) = bl_short(5);  aop_long(i + 1) = bl_long(5)
    end do

    print *, "Propagation completed successfully."

    ! 3. Generate Comparative Plots using pyplot-fortran
    print *, "Generating element plots..."

    ! Plot 1: Semi-major Axis Comparison
    call plt%initialize(title="Semi-Major Axis: Osculating vs Brouwer Mean", &
                        xlabel="Time (hours)", &
                        ylabel="Semi-major axis $a$ (km)", &
                        legend = .true.)
    call plt%add_plot(t_hrs, sma_osc,   label="Osculating", linestyle="-", color=c0)
    call plt%add_plot(t_hrs, sma_short, label="Brouwer Short-Period Mean", linestyle="--", color=c1)
    call plt%add_plot(t_hrs, sma_long,  label="Brouwer Long-Period Mean", linestyle=":", color=c2)
    call plt%savefig("brouwer_sma_comparison.png")

    ! Plot 2: Eccentricity Comparison
    call plt%initialize(title="Eccentricity: Osculating vs Brouwer Mean", &
                     xlabel="Time (hours)", &
                     ylabel="Eccentricity $e$", &
                     legend = .true.)
    call plt%add_plot(t_hrs, ecc_osc,   label="Osculating", linestyle="-", color=c0)
    call plt%add_plot(t_hrs, ecc_short, label="Brouwer Short-Period Mean", linestyle="--", color=c1)
    call plt%add_plot(t_hrs, ecc_long,  label="Brouwer Long-Period Mean", linestyle=":", color=c2)
    call plt%savefig("brouwer_ecc_comparison.png" )

    ! Plot 3: Inclination Comparison
    call plt%initialize(title="Inclination: Osculating vs Brouwer Mean", &
                        xlabel="Time (hours)", &
                        ylabel="Inclination $i$ (deg)", &
                        legend = .true.)
    call plt%add_plot(t_hrs, inc_osc,   label="Osculating", linestyle="-", color=c0)
    call plt%add_plot(t_hrs, inc_short, label="Brouwer Short-Period Mean", linestyle="--", color=c1)
    call plt%add_plot(t_hrs, inc_long,  label="Brouwer Long-Period Mean", linestyle=":", color=c2)
    call plt%savefig("brouwer_inc_comparison.png")

    ! Plot 4: Argument of Periapsis Comparison
    call plt%initialize(title="Argument of Periapsis: Osculating vs Brouwer Mean", &
                        xlabel="Time (hours)", &
                        ylabel="Argument of Periapsis $\\omega$ (deg)", &
                        legend = .true.)
    call plt%add_plot(t_hrs, aop_osc,   label="Osculating", linestyle="-", color=c0)
    call plt%add_plot(t_hrs, aop_short, label="Brouwer Short-Period Mean", linestyle="--", color=c1)
    call plt%add_plot(t_hrs, aop_long,  label="Brouwer Long-Period Mean", linestyle=":", color=c2)
    call plt%savefig("brouwer_aop_comparison.png")

    print *, "All plots generated successfully: brouwer_*.png"
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
        ! call gravity_j2_j3_j4_j5(y(1:3),mu_earth,req_earth,j2_earth,j3_earth,j4_earth,j5_earth,acc)
        call gravity_j2_j3_j4(y(1:3),mu_earth,req_earth,j2_earth,j3_earth,j4_earth,acc)
        dydt(4:6) = acc(1:3)

    end subroutine grav_derivs

!     subroutine gravity_j2_j3_j4_j5(r,mu,req,j2,j3,j4,j5,acc)

!     !!@warning This is an AI-generated routine that seems to be wrong?

!     real(wp),dimension(3),intent(in)  :: r   !! satellite position vector [km]
!     real(wp),intent(in)               :: mu  !! central body gravitational parameter [km^3/s^2]
!     real(wp),intent(in)               :: req !! body equatorial radius [km]
!     real(wp),intent(in)               :: j2  !! j2 coefficient
!     real(wp),intent(in)               :: j3  !! j3 coefficient
!     real(wp),intent(in)               :: j4  !! j4 coefficient
!     real(wp),intent(in)               :: j5  !! j5 coefficient
!     real(wp),dimension(3),intent(out) :: acc !! gravity acceleration vector [km/s^2]

!     real(wp) :: rmag, r2, r3, r_xy2, z, z_r, z2_r2, z3_r3, z4_r4
!     real(wp) :: re_r, re_r2, re_r3, re_r4, re_r5
!     real(wp) :: f_r, f_z, mu_r3

!     r2 = r(1)**2 + r(2)**2 + r(3)**2
!     rmag = sqrt(r2)
!     r3 = rmag * r2
!     z = r(3)
!     z_r = z / rmag
!     z2_r2 = z_r * z_r
!     z3_r3 = z2_r2 * z_r
!     z4_r4 = z2_r2 * z2_r2

!     re_r = req / rmag
!     re_r2 = re_r * re_r
!     re_r3 = re_r2 * re_r
!     re_r4 = re_r3 * re_r
!     re_r5 = re_r4 * re_r

!     mu_r3 = mu / r3

!     ! Potential partial derivatives: a = - (mu/rmag^3)*r_vec + a_pert
!     ! Common factors for zonal perturbation acceleration components
!     f_r = 1.0_wp + 1.5_wp * j2 * re_r2 * (1.0_wp - 5.0_wp * z2_r2) &
!             + 2.5_wp * j3 * re_r3 * (3.0_wp * z_r - 7.0_wp * z3_r3) &
!             - 0.625_wp * j4 * re_r4 * (3.0_wp - 42.0_wp * z2_r2 + 63.0_wp * z4_r4) &
!             - 2.625_wp * j5 * re_r5 * (5.0_wp * z_r - 30.0_wp * z3_r3 + 33.0_wp * z4_r4 * z_r)

!     f_z = 3.0_wp * j2 * re_r2 * z_r &
!             - 0.5_wp * j3 * re_r3 * (3.0_wp - 7.0_wp * z2_r2) &
!             - 2.5_wp * j4 * re_r4 * (3.0_wp * z_r - 7.0_wp * z3_r3) &
!             + 0.625_wp * j5 * re_r5 * (1.0_wp - 14.0_wp * z2_r2 + 21.0_wp * z4_r4)

!     acc(1) = -mu_r3 * r(1) * f_r
!     acc(2) = -mu_r3 * r(2) * f_r
!     acc(3) = -mu_r3 * (r(3) * f_r - rmag * f_z)

!     end subroutine gravity_j2_j3_j4_j5
! !*****************************************************************************************

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
