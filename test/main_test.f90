program main_test

    use iso_fortran_env, only: wp => real64
    use brouwer_module

    implicit none

    real(wp), dimension(6) :: cart_orig, blms, cart_roundtrip, blms_roundtrip
    real(wp), dimension(6) :: blml, cart_long, blml_roundtrip
    real(wp), dimension(6) :: kep_orig, cart_test, kep_res
    real(wp) :: ta, ma, ea, ha, diff
    integer :: stat
    logical :: all_passed

    all_passed = .true.

    print *, "=========================================================="
    print *, " Testing Modern Fortran Brouwer Element Conversions"
    print *, "=========================================================="

    ! -------------------------------------------------------------
    ! Test 1: Standard LEO Orbit Roundtrips
    ! -------------------------------------------------------------
    kep_orig = [7000.0_wp, 0.01_wp, 28.5_wp, 45.0_wp, 30.0_wp, 15.0_wp]
    cart_orig = keplerian_to_cartesian(mu_earth, kep_orig, anomaly_type="TA", stat=stat)
    call check(stat == 0, "keplerian_to_cartesian LEO")

    blms = cartesian_to_brouwer_mean_short(mu_earth, cart_orig, stat=stat)
    call check(stat == 0, "cartesian_to_brouwer_mean_short LEO")
    cart_roundtrip = brouwer_mean_short_to_cartesian(mu_earth, blms, stat=stat)
    call check(stat == 0, "brouwer_mean_short_to_cartesian LEO")
    diff = norm2(cart_orig - cart_roundtrip) / norm2(cart_orig)
    call check(diff < 1.0e-7_wp, "Brouwer Short LEO roundtrip")

    blml = cartesian_to_brouwer_mean_long(mu_earth, cart_orig, stat=stat)
    call check(stat == 0, "cartesian_to_brouwer_mean_long LEO")
    cart_long = brouwer_mean_long_to_cartesian(mu_earth, blml, stat=stat)
    call check(stat == 0, "brouwer_mean_long_to_cartesian LEO")
    diff = norm2(cart_orig - cart_long) / norm2(cart_orig)
    call check(diff < 1.0e-7_wp, "Brouwer Long LEO roundtrip")

    blms_roundtrip = cartesian_to_brouwer_mean_short(mu_earth, cart_roundtrip, stat=stat)
    diff = norm2(blms - blms_roundtrip) / norm2(blms)
    call check(diff < 1.0e-7_wp, "Brouwer Short elements roundtrip")

    blml_roundtrip = cartesian_to_brouwer_mean_long(mu_earth, cart_long, stat=stat)
    diff = norm2(blml - blml_roundtrip) / norm2(blml)
    call check(diff < 1.0e-7_wp, "Brouwer Long elements roundtrip")

    ! -------------------------------------------------------------
    ! Test 2: High Inclination / Polar Orbit
    ! -------------------------------------------------------------
    kep_orig = [7200.0_wp, 0.005_wp, 87.0_wp, 120.0_wp, 45.0_wp, 200.0_wp]
    cart_orig = keplerian_to_cartesian(mu_earth, kep_orig, anomaly_type="TA", stat=stat)
    blms = cartesian_to_brouwer_mean_short(mu_earth, cart_orig, stat=stat)
    cart_roundtrip = brouwer_mean_short_to_cartesian(mu_earth, blms, stat=stat)
    diff = norm2(cart_orig - cart_roundtrip) / norm2(cart_orig)
    call check(diff < 1.0e-7_wp, "High inclination Short roundtrip")

    blml = cartesian_to_brouwer_mean_long(mu_earth, cart_orig, stat=stat)
    cart_long = brouwer_mean_long_to_cartesian(mu_earth, blml, stat=stat)
    diff = norm2(cart_orig - cart_long) / norm2(cart_orig)
    call check(diff < 1.0e-7_wp, "High inclination Long roundtrip")

    ! -------------------------------------------------------------
    ! Test 3: Retrograde Orbit (inc > 90 deg)
    ! -------------------------------------------------------------
    kep_orig = [8000.0_wp, 0.05_wp, 105.0_wp, 210.0_wp, 75.0_wp, 315.0_wp]
    cart_orig = keplerian_to_cartesian(mu_earth, kep_orig, anomaly_type="TA", stat=stat)
    blms = cartesian_to_brouwer_mean_short(mu_earth, cart_orig, stat=stat)
    cart_roundtrip = brouwer_mean_short_to_cartesian(mu_earth, blms, stat=stat)
    diff = norm2(cart_orig - cart_roundtrip) / norm2(cart_orig)
    call check(diff < 1.0e-7_wp, "Retrograde Short roundtrip")

    blml = cartesian_to_brouwer_mean_long(mu_earth, cart_orig, stat=stat)
    cart_long = brouwer_mean_long_to_cartesian(mu_earth, blml, stat=stat)
    diff = norm2(cart_orig - cart_long) / norm2(cart_orig)
    call check(diff < 1.0e-7_wp, "Retrograde Long roundtrip")

    ! -------------------------------------------------------------
    ! Test 4: Critical Inclination Orbit (i ~ 63.435 deg and 116.565 deg)
    ! -------------------------------------------------------------
    kep_orig = [7500.0_wp, 0.02_wp, 63.4349488_wp, 50.0_wp, 40.0_wp, 80.0_wp]
    cart_orig = keplerian_to_cartesian(mu_earth, kep_orig, anomaly_type="TA", stat=stat)
    blml = cartesian_to_brouwer_mean_long(mu_earth, cart_orig, stat=stat)
    cart_long = brouwer_mean_long_to_cartesian(mu_earth, blml, stat=stat)
    call check(stat == 0, "Critical inclination prograde Long conversion")

    kep_orig = [7500.0_wp, 0.02_wp, 116.5650512_wp, 50.0_wp, 40.0_wp, 80.0_wp]
    cart_orig = keplerian_to_cartesian(mu_earth, kep_orig, anomaly_type="TA", stat=stat)
    blml = cartesian_to_brouwer_mean_long(mu_earth, cart_orig, stat=stat)
    cart_long = brouwer_mean_long_to_cartesian(mu_earth, blml, stat=stat)
    call check(stat == 0, "Critical inclination retrograde Long conversion")

    ! -------------------------------------------------------------
    ! Test 5: Pseudostate Orbit (i > 175 deg)
    ! -------------------------------------------------------------
    kep_orig = [7100.0_wp, 0.01_wp, 177.0_wp, 60.0_wp, 45.0_wp, 100.0_wp]
    cart_orig = keplerian_to_cartesian(mu_earth, kep_orig, anomaly_type="TA", stat=stat)
    blms = cartesian_to_brouwer_mean_short(mu_earth, cart_orig, stat=stat)
    cart_roundtrip = brouwer_mean_short_to_cartesian(mu_earth, blms, stat=stat)
    diff = norm2(cart_orig - cart_roundtrip) / norm2(cart_orig)
    call check(diff < 1.0e-7_wp, "Pseudostate (i > 175) Short roundtrip")

    blml = cartesian_to_brouwer_mean_long(mu_earth, cart_orig, stat=stat)
    cart_long = brouwer_mean_long_to_cartesian(mu_earth, blml, stat=stat)
    diff = norm2(cart_orig - cart_long) / norm2(cart_orig)
    call check(diff < 1.0e-7_wp, "Pseudostate (i > 175) Long roundtrip")

    ! -------------------------------------------------------------
    ! Test 6: Osculating direct calls with special cases
    ! -------------------------------------------------------------
    ! Circular mean elements (ecc <= 1e-11)
    blms = [7000.0_wp, 0.0_wp, 28.5_wp, 45.0_wp, 30.0_wp, 50.0_wp]
    kep_res = brouwer_mean_short_to_osculating(mu_earth, blms, stat=stat)
    call check(stat == 0, "Short osculating circular inclined")

    ! Circular equatorial (ecc <= 1e-11, inc <= 1e-7)
    blms = [7000.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 50.0_wp]
    kep_res = brouwer_mean_short_to_osculating(mu_earth, blms, stat=stat)
    call check(stat == 0, "Short osculating circular equatorial")

    ! Negative eccentricity input
    blms = [7000.0_wp, -0.01_wp, 28.5_wp, 45.0_wp, 30.0_wp, 50.0_wp]
    kep_res = brouwer_mean_short_to_osculating(mu_earth, blms, stat=stat)
    call check(stat == 0, "Short osculating negative eccentricity")

    ! Long osculating circular inclined & circular equatorial
    blml = [7000.0_wp, 0.0_wp, 28.5_wp, 45.0_wp, 30.0_wp, 50.0_wp]
    kep_res = brouwer_mean_long_to_osculating(mu_earth, blml, stat=stat)
    call check(stat == 0, "Long osculating circular inclined")

    blml = [7000.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 50.0_wp]
    kep_res = brouwer_mean_long_to_osculating(mu_earth, blml, stat=stat)
    call check(stat == 0, "Long osculating circular equatorial")

    ! Long osculating inclined with near-equatorial inc <= 1e-7 and ecc > 0
    blml = [7000.0_wp, 0.01_wp, 0.0_wp, 0.0_wp, 30.0_wp, 50.0_wp]
    kep_res = brouwer_mean_long_to_osculating(mu_earth, blml, stat=stat)
    call check(stat == 0, "Long osculating elliptic equatorial")

    ! -------------------------------------------------------------
    ! Test 7: Cartesian to Keplerian Orbital Regimes
    ! -------------------------------------------------------------
    ! Case 1: Non-circular, Inclined with quadrant flips
    kep_orig = [8000.0_wp, 0.1_wp, 45.0_wp, 200.0_wp, 210.0_wp, 220.0_wp]
    cart_test = keplerian_to_cartesian(mu_earth, kep_orig, anomaly_type="TA", stat=stat)
    kep_res = cartesian_to_keplerian(mu_earth, cart_test, anomaly_type="TA", stat=stat)
    diff = norm2(kep_orig - kep_res) / norm2(kep_orig)
    call check(diff < 1.0e-7_wp, "Cart2Kep Case 1 quadrant flips")

    ! Case 2: Non-circular, Equatorial Prograde (inc = 0)
    kep_orig = [8000.0_wp, 0.1_wp, 0.0_wp, 0.0_wp, 40.0_wp, 50.0_wp]
    cart_test = keplerian_to_cartesian(mu_earth, kep_orig, anomaly_type="TA", stat=stat)
    kep_res = cartesian_to_keplerian(mu_earth, cart_test, anomaly_type="TA", stat=stat)
    call check(abs(kep_res(1) - 8000.0_wp) < 1.0e-4_wp .and. abs(kep_res(2) - 0.1_wp) < 1.0e-6_wp, "Cart2Kep Case 2 prograde")

    ! Case 2b: Non-circular, Equatorial Prograde with pos.vel < 0
    kep_orig = [8000.0_wp, 0.1_wp, 0.0_wp, 0.0_wp, 200.0_wp, 220.0_wp]
    cart_test = keplerian_to_cartesian(mu_earth, kep_orig, anomaly_type="TA", stat=stat)
    kep_res = cartesian_to_keplerian(mu_earth, cart_test, anomaly_type="TA", stat=stat)
    call check(abs(kep_res(1) - 8000.0_wp) < 1.0e-4_wp, "Cart2Kep Case 2 prograde quadrant 2")

    ! Case 2c: Non-circular, Equatorial Retrograde (inc = 180)
    kep_orig = [8000.0_wp, 0.1_wp, 180.0_wp, 0.0_wp, 40.0_wp, 50.0_wp]
    cart_test = keplerian_to_cartesian(mu_earth, kep_orig, anomaly_type="TA", stat=stat)
    kep_res = cartesian_to_keplerian(mu_earth, cart_test, anomaly_type="TA", stat=stat)
    call check(abs(kep_res(1) - 8000.0_wp) < 1.0e-4_wp, "Cart2Kep Case 2 retrograde")

    ! Case 3: Circular, Inclined (e = 0, inc = 45)
    kep_orig = [8000.0_wp, 0.0_wp, 45.0_wp, 30.0_wp, 0.0_wp, 60.0_wp]
    cart_test = keplerian_to_cartesian(mu_earth, kep_orig, anomaly_type="TA", stat=stat)
    kep_res = cartesian_to_keplerian(mu_earth, cart_test, anomaly_type="TA", stat=stat)
    call check(abs(kep_res(1) - 8000.0_wp) < 1.0e-4_wp .and. abs(kep_res(3) - 45.0_wp) < 1.0e-5_wp, "Cart2Kep Case 3")

    ! Case 3b: Circular, Inclined with pos(3) < 0 and nodeVec(2) < 0
    kep_orig = [8000.0_wp, 0.0_wp, 45.0_wp, 200.0_wp, 0.0_wp, 220.0_wp]
    cart_test = keplerian_to_cartesian(mu_earth, kep_orig, anomaly_type="TA", stat=stat)
    kep_res = cartesian_to_keplerian(mu_earth, cart_test, anomaly_type="TA", stat=stat)
    call check(abs(kep_res(1) - 8000.0_wp) < 1.0e-4_wp, "Cart2Kep Case 3b")

    ! Case 4: Circular, Equatorial Prograde (e = 0, inc = 0)
    kep_orig = [8000.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 45.0_wp]
    cart_test = keplerian_to_cartesian(mu_earth, kep_orig, anomaly_type="TA", stat=stat)
    kep_res = cartesian_to_keplerian(mu_earth, cart_test, anomaly_type="TA", stat=stat)
    call check(abs(kep_res(1) - 8000.0_wp) < 1.0e-4_wp, "Cart2Kep Case 4 prograde")

    ! Case 4b: Circular, Equatorial with pos(2) < 0
    kep_orig = [8000.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 220.0_wp]
    cart_test = keplerian_to_cartesian(mu_earth, kep_orig, anomaly_type="TA", stat=stat)
    kep_res = cartesian_to_keplerian(mu_earth, cart_test, anomaly_type="TA", stat=stat)
    call check(abs(kep_res(1) - 8000.0_wp) < 1.0e-4_wp, "Cart2Kep Case 4b")

    ! Case 4c: Circular, Equatorial Retrograde (e = 0, inc = 180)
    kep_orig = [8000.0_wp, 0.0_wp, 180.0_wp, 0.0_wp, 0.0_wp, 45.0_wp]
    cart_test = keplerian_to_cartesian(mu_earth, kep_orig, anomaly_type="TA", stat=stat)
    kep_res = cartesian_to_keplerian(mu_earth, cart_test, anomaly_type="TA", stat=stat)
    call check(abs(kep_res(1) - 8000.0_wp) < 1.0e-4_wp, "Cart2Kep Case 4c retrograde")

    ! Convert to Mean Anomaly mode in CartesianToKeplerian & KeplerianToCartesian
    kep_orig = [7500.0_wp, 0.05_wp, 35.0_wp, 45.0_wp, 60.0_wp, 75.0_wp]
    cart_test = keplerian_to_cartesian(mu_earth, kep_orig, anomaly_type="MA", stat=stat)
    kep_res = cartesian_to_keplerian(mu_earth, cart_test, anomaly_type="MA", stat=stat)
    diff = norm2(kep_orig - kep_res) / norm2(kep_orig)
    call check(diff < 1.0e-7_wp, "Cart2Kep and Kep2Cart with MA anomaly type")

    ! -------------------------------------------------------------
    ! Test 8: Anomaly Conversions (Elliptic, Parabolic, Hyperbolic)
    ! -------------------------------------------------------------
    ! Elliptic: True <-> Eccentric <-> Mean
    ta = 45.0_wp * deg2rad
    ea = true_to_eccentric_anomaly(ta, 0.1_wp)
    ma = true_to_mean_anomaly(ta, 0.1_wp)
    ta = mean_to_true_anomaly(ma, 0.1_wp, tol=1.0e-12_wp, stat=stat)
    call check(abs(ta - 45.0_wp * deg2rad) < 1.0e-10_wp .and. stat == 0, "Elliptic anomaly roundtrip")

    ! Special Elliptic angle E = pi, M = pi
    ta = pi
    ea = true_to_eccentric_anomaly(ta, 0.2_wp)
    ma = true_to_mean_anomaly(ta, 0.2_wp)
    ta = mean_to_true_anomaly(ma, 0.2_wp, stat=stat)
    call check(abs(ta - pi) < 1.0e-10_wp, "Anomaly at pi")

    ! Parabolic case
    ma = true_to_mean_anomaly(45.0_wp * deg2rad, 1.0_wp)
    call check(ma == 0.0_wp, "Parabolic true_to_mean_anomaly returns 0")

    ! Hyperbolic case
    ha = true_to_hyperbolic_anomaly(0.5_wp, 1.5_wp)
    call check(ha > 0.0_wp, "Hyperbolic true_to_hyperbolic_anomaly")
    ma = true_to_mean_anomaly(0.5_wp, 1.5_wp)
    call check(ma > 0.0_wp, "Hyperbolic true_to_mean_anomaly")
    ta = mean_to_true_anomaly(ma, 1.5_wp, stat=stat)
    call check(abs(ta - 0.5_wp) < 1.0e-8_wp .and. stat == 0, "Hyperbolic mean_to_true_anomaly roundtrip")

    ! -------------------------------------------------------------
    ! Test 9: Error Handling / Status Flags
    ! -------------------------------------------------------------
    ! Invalid mu for Brouwer Short
    blms = cartesian_to_brouwer_mean_short(100.0_wp, cart_orig, stat=stat)
    call check(stat == 1, "Short error on invalid mu")

    ! Invalid mu for Brouwer Long
    blml = cartesian_to_brouwer_mean_long(100.0_wp, cart_orig, stat=stat)
    call check(stat == 1, "Long error on invalid mu")

    ! Invalid mu for Osculating Short & Long
    kep_res = brouwer_mean_short_to_osculating(100.0_wp, blms, stat=stat)
    call check(stat == 1, "Short osculating error on invalid mu")
    cart_test = brouwer_mean_short_to_cartesian(100.0_wp, blms, stat=stat)
    call check(stat == 1, "Short cartesian error on invalid mu")

    kep_res = brouwer_mean_long_to_osculating(100.0_wp, blml, stat=stat)
    call check(stat == 1, "Long osculating error on invalid mu")
    cart_test = brouwer_mean_long_to_cartesian(100.0_wp, blml, stat=stat)
    call check(stat == 1, "Long cartesian error on invalid mu")

    ! Invalid inclination > 180
    cart_test = [10000.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 6.0_wp, 0.0_wp]
    kep_res = [7000.0_wp, 0.01_wp, 185.0_wp, 0.0_wp, 0.0_wp, 0.0_wp]
    kep_res = brouwer_mean_short_to_osculating(mu_earth, kep_res, stat=stat)
    call check(stat == 2, "Short osculating error on inc > 180")

    kep_res = [7000.0_wp, 0.01_wp, 185.0_wp, 0.0_wp, 0.0_wp, 0.0_wp]
    kep_res = brouwer_mean_long_to_osculating(mu_earth, kep_res, stat=stat)
    call check(stat == 2, "Long osculating error on inc > 180")

    ! Invalid eccentricity >= 0.99
    kep_res = [1000000.0_wp, 0.995_wp, 28.5_wp, 0.0_wp, 0.0_wp, 0.0_wp]
    kep_res = brouwer_mean_short_to_osculating(mu_earth, kep_res, stat=stat)
    call check(stat == 3, "Short osculating error on ecc >= 0.99")

    kep_res = [1000000.0_wp, 0.995_wp, 28.5_wp, 0.0_wp, 0.0_wp, 0.0_wp]
    kep_res = brouwer_mean_long_to_osculating(mu_earth, kep_res, stat=stat)
    call check(stat == 3, "Long osculating error on ecc >= 0.99")

    ! Invalid periapsis radius < 3000 km
    kep_res = [3500.0_wp, 0.2_wp, 28.5_wp, 0.0_wp, 0.0_wp, 0.0_wp]
    kep_res = brouwer_mean_short_to_osculating(mu_earth, kep_res, stat=stat)
    call check(stat == 4, "Short osculating error on radper < 3000 km")

    kep_res = [3500.0_wp, 0.2_wp, 28.5_wp, 0.0_wp, 0.0_wp, 0.0_wp]
    kep_res = brouwer_mean_long_to_osculating(mu_earth, kep_res, stat=stat)
    call check(stat == 4, "Long osculating error on radper < 3000 km")

    ! Cartesian to Brouwer errors on radper < 3000 km and ecc >= 0.99
    cart_test = keplerian_to_cartesian(mu_earth, [3500.0_wp, 0.2_wp, 28.5_wp, 0.0_wp, 0.0_wp, 0.0_wp])
    blms = cartesian_to_brouwer_mean_short(mu_earth, cart_test, stat=stat)
    call check(stat == 4, "Cart2BrouwerShort error on radper < 3000 km")
    blml = cartesian_to_brouwer_mean_long(mu_earth, cart_test, stat=stat)
    call check(stat == 4, "Cart2BrouwerLong error on radper < 3000 km")

    cart_test = keplerian_to_cartesian(mu_earth, [1000000.0_wp, 0.995_wp, 28.5_wp, 0.0_wp, 0.0_wp, 0.0_wp])
    blms = cartesian_to_brouwer_mean_short(mu_earth, cart_test, stat=stat)
    call check(stat == 3, "Cart2BrouwerShort error on ecc >= 0.99")
    blml = cartesian_to_brouwer_mean_long(mu_earth, cart_test, stat=stat)
    call check(stat == 3, "Cart2BrouwerLong error on ecc >= 0.99")

    ! Zero / invalid Cartesian inputs
    cart_test = [0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp]
    kep_res = cartesian_to_keplerian(mu_earth, cart_test, stat=stat)
    call check(stat == 6, "Cart2Kep error on zero pos/vel")

    kep_res = cartesian_to_keplerian(0.0_wp, cart_orig, stat=stat)
    call check(stat == 6, "Cart2Kep error on zero mu")

    ! Singular radius in Keplerian to Cartesian
    kep_res = [0.0_wp, 1.0_wp, 28.5_wp, 0.0_wp, 0.0_wp, 0.0_wp]
    cart_test = keplerian_to_cartesian(mu_earth, kep_res, stat=stat)
    call check(stat == 6, "Kep2Cart error on singular p")

    print *, "=========================================================="
    if (all_passed) then
        print *, " ALL TESTS PASSED SUCCESSFULLY!"
    else
        stop 1
    end if
    print *, "=========================================================="

contains

    subroutine check(condition, test_name)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: test_name
        if (condition) then
            print *, "PASS: ", test_name
        else
            print *, "FAIL: ", test_name
            all_passed = .false.
        end if
    end subroutine check

end program main_test
