program main_test

    use brouwer_module, wp => brouwer_module_wp

    implicit none

    real(wp), parameter :: pi = acos(-1.0_wp)
    real(wp), parameter :: deg2rad = pi / 180.0_wp

    ! Earth parameters
    real(wp), parameter :: mu_earth  = 398600.4415_wp           ! km^3/s^2
    real(wp), parameter :: req_earth = 6378.1363_wp             ! km
    real(wp), parameter :: j2_earth  = 1.082626925638815e-3_wp
    real(wp), parameter :: j3_earth  = -0.2532307818191774e-5_wp
    real(wp), parameter :: j4_earth  = -0.1620429990000000e-5_wp
    real(wp), parameter :: j5_earth  = -0.2270711043920343e-6_wp

    ! Mars parameters
    real(wp), parameter :: mu_mars  = 42828.375214_wp           ! km^3/s^2
    real(wp), parameter :: req_mars = 3396.19_wp                ! km
    real(wp), parameter :: j2_mars  = 1.96045e-3_wp
    real(wp), parameter :: j3_mars  = 3.15e-5_wp
    real(wp), parameter :: j4_mars  = -1.54e-5_wp
    real(wp), parameter :: j5_mars  = 0.0_wp

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
    ! Test 1: Standard LEO Orbit Roundtrips (Earth)
    ! -------------------------------------------------------------
    kep_orig = [7000.0_wp, 0.01_wp, 28.5_wp, 45.0_wp, 30.0_wp, 15.0_wp]
    call keplerian_to_cartesian(mu_earth, kep_orig, anomaly_type="TA", stat=stat, cart=cart_orig)
    call check(stat == 0, "keplerian_to_cartesian LEO")

    call cartesian_to_brouwer_mean_short(mu_earth, req_earth, j2_earth, cart_orig, stat=stat, blms=blms)
    call check(stat == 0, "cartesian_to_brouwer_mean_short LEO")
    call brouwer_mean_short_to_cartesian(mu_earth, req_earth, j2_earth, blms, stat=stat, cart=cart_roundtrip)
    call check(stat == 0, "brouwer_mean_short_to_cartesian LEO")
    diff = norm2(cart_orig - cart_roundtrip) / norm2(cart_orig)
    call check(diff < 1.0e-7_wp, "Brouwer Short LEO roundtrip")

    call cartesian_to_brouwer_mean_long(mu_earth, req_earth, j2_earth, j3_earth, j4_earth, j5_earth, cart_orig, stat=stat, blml=blml)
    call check(stat == 0, "cartesian_to_brouwer_mean_long LEO")
    call brouwer_mean_long_to_cartesian(mu_earth, req_earth, j2_earth, j3_earth, j4_earth, j5_earth, blml, stat=stat, cart=cart_long)
    call check(stat == 0, "brouwer_mean_long_to_cartesian LEO")
    diff = norm2(cart_orig - cart_long) / norm2(cart_orig)
    call check(diff < 1.0e-7_wp, "Brouwer Long LEO roundtrip")

    call cartesian_to_brouwer_mean_short(mu_earth, req_earth, j2_earth, cart_roundtrip, stat=stat, blms=blms_roundtrip)
    diff = norm2(blms - blms_roundtrip) / norm2(blms)
    call check(diff < 1.0e-7_wp, "Brouwer Short elements roundtrip")

    call cartesian_to_brouwer_mean_long(mu_earth, req_earth, j2_earth, j3_earth, j4_earth, j5_earth, cart_long, stat=stat, blml=blml_roundtrip)
    diff = norm2(blml - blml_roundtrip) / norm2(blml)
    call check(diff < 1.0e-7_wp, "Brouwer Long elements roundtrip")

    ! -------------------------------------------------------------
    ! Test 2: High Inclination / Polar Orbit
    ! -------------------------------------------------------------
    kep_orig = [7200.0_wp, 0.005_wp, 87.0_wp, 120.0_wp, 45.0_wp, 200.0_wp]
    call keplerian_to_cartesian(mu_earth, kep_orig, anomaly_type="TA", stat=stat, cart=cart_orig)
    call cartesian_to_brouwer_mean_short(mu_earth, req_earth, j2_earth, cart_orig, stat=stat, blms=blms)
    call brouwer_mean_short_to_cartesian(mu_earth, req_earth, j2_earth, blms, stat=stat, cart=cart_roundtrip)
    diff = norm2(cart_orig - cart_roundtrip) / norm2(cart_orig)
    call check(diff < 1.0e-7_wp, "High inclination Short roundtrip")

    call cartesian_to_brouwer_mean_long(mu_earth, req_earth, j2_earth, j3_earth, j4_earth, j5_earth, cart_orig, stat=stat, blml=blml)
    call brouwer_mean_long_to_cartesian(mu_earth, req_earth, j2_earth, j3_earth, j4_earth, j5_earth, blml, stat=stat, cart=cart_long)
    diff = norm2(cart_orig - cart_long) / norm2(cart_orig)
    call check(diff < 1.0e-7_wp, "High inclination Long roundtrip")

    ! -------------------------------------------------------------
    ! Test 3: Retrograde Orbit (inc > 90 deg)
    ! -------------------------------------------------------------
    kep_orig = [8000.0_wp, 0.05_wp, 105.0_wp, 210.0_wp, 75.0_wp, 315.0_wp]
    call keplerian_to_cartesian(mu_earth, kep_orig, anomaly_type="TA", stat=stat, cart=cart_orig)
    call cartesian_to_brouwer_mean_short(mu_earth, req_earth, j2_earth, cart_orig, stat=stat, blms=blms)
    call brouwer_mean_short_to_cartesian(mu_earth, req_earth, j2_earth, blms, stat=stat, cart=cart_roundtrip)
    diff = norm2(cart_orig - cart_roundtrip) / norm2(cart_orig)
    call check(diff < 1.0e-7_wp, "Retrograde Short roundtrip")

    call cartesian_to_brouwer_mean_long(mu_earth, req_earth, j2_earth, j3_earth, j4_earth, j5_earth, cart_orig, stat=stat, blml=blml)
    call brouwer_mean_long_to_cartesian(mu_earth, req_earth, j2_earth, j3_earth, j4_earth, j5_earth, blml, stat=stat, cart=cart_long)
    diff = norm2(cart_orig - cart_long) / norm2(cart_orig)
    call check(diff < 1.0e-7_wp, "Retrograde Long roundtrip")

    ! -------------------------------------------------------------
    ! Test 4: Critical Inclination Orbit (i ~ 63.435 deg and 116.565 deg)
    ! -------------------------------------------------------------
    ! Brouwer Short-period is non-singular at critical inclination:
    kep_orig = [7500.0_wp, 0.02_wp, 63.4349488_wp, 50.0_wp, 40.0_wp, 80.0_wp]
    call keplerian_to_cartesian(mu_earth, kep_orig, anomaly_type="TA", stat=stat, cart=cart_orig)
    call cartesian_to_brouwer_mean_short(mu_earth, req_earth, j2_earth, cart_orig, stat=stat, blms=blms)
    call brouwer_mean_short_to_cartesian(mu_earth, req_earth, j2_earth, blms, stat=stat, cart=cart_roundtrip)
    diff = norm2(cart_orig - cart_roundtrip) / norm2(cart_orig)
    ! Brouwer Long-period has a singularity at critical inclination (1 - 5*cos^2(i) = 0),
    ! triggering non-convergence flag (stat = 5):
    call cartesian_to_brouwer_mean_long(mu_earth, req_earth, j2_earth, j3_earth, j4_earth, j5_earth, cart_orig, stat=stat, blml=blml)
    call check(stat == 6, "Critical inclination prograde Long returns stat=6 (singularity)")

    kep_orig = [7500.0_wp, 0.02_wp, 116.5650512_wp, 50.0_wp, 40.0_wp, 80.0_wp]
    call keplerian_to_cartesian(mu_earth, kep_orig, anomaly_type="TA", stat=stat, cart=cart_orig)
    call cartesian_to_brouwer_mean_long(mu_earth, req_earth, j2_earth, j3_earth, j4_earth, j5_earth, cart_orig, stat=stat, blml=blml)
    call check(stat == 6, "Critical inclination retrograde Long returns stat=6 (singularity)")

    ! -------------------------------------------------------------
    ! Test 5: Pseudostate Orbit (i > 175 deg)
    ! -------------------------------------------------------------
    kep_orig = [7100.0_wp, 0.01_wp, 177.0_wp, 60.0_wp, 45.0_wp, 100.0_wp]
    call keplerian_to_cartesian(mu_earth, kep_orig, anomaly_type="TA", stat=stat, cart=cart_orig); if (stat /= 0) error stop "Error converting pseudostate Keplerian to Cartesian state."
    call cartesian_to_brouwer_mean_short(mu_earth, req_earth, j2_earth, cart_orig, stat=stat, blms=blms); if (stat /= 0) error stop "Error converting pseudostate Cartesian to Brouwer Short state."
    call brouwer_mean_short_to_cartesian(mu_earth, req_earth, j2_earth, blms, stat=stat, cart=cart_roundtrip); if (stat /= 0) error stop "Error converting pseudostate Brouwer Short to Cartesian state."
    diff = norm2(cart_orig - cart_roundtrip) / norm2(cart_orig)
    call check(diff < 1.0e-7_wp, "Pseudostate (i > 175) Short roundtrip")

    call cartesian_to_brouwer_mean_long(mu_earth, req_earth, j2_earth, j3_earth, j4_earth, j5_earth, cart_orig, stat=stat, blml=blml); if (stat /= 0) error stop "Error converting pseudostate Cartesian to Brouwer Long state."
    call brouwer_mean_long_to_cartesian(mu_earth, req_earth, j2_earth, j3_earth, j4_earth, j5_earth, blml, stat=stat, cart=cart_long); if (stat /= 0) error stop "Error converting pseudostate Brouwer Long to Cartesian state."
    diff = norm2(cart_orig - cart_long) / norm2(cart_orig)
    call check(diff < 1.0e-7_wp, "Pseudostate (i > 175) Long roundtrip")

    ! -------------------------------------------------------------
    ! Test 6: Mars Satellite Orbit Roundtrip
    ! -------------------------------------------------------------
    kep_orig = [4000.0_wp, 0.02_wp, 45.0_wp, 70.0_wp, 25.0_wp, 60.0_wp]
    call keplerian_to_cartesian(mu_mars, kep_orig, anomaly_type="TA", stat=stat, cart=cart_orig)
    call check(stat == 0, "keplerian_to_cartesian Mars")

    call cartesian_to_brouwer_mean_short(mu_mars, req_mars, j2_mars, cart_orig, stat=stat, blms=blms)
    call check(stat == 0, "cartesian_to_brouwer_mean_short Mars")
    call brouwer_mean_short_to_cartesian(mu_mars, req_mars, j2_mars, blms, stat=stat, cart=cart_roundtrip)
    call check(stat == 0, "brouwer_mean_short_to_cartesian Mars")
    diff = norm2(cart_orig - cart_roundtrip) / norm2(cart_orig)
    call check(diff < 1.0e-7_wp, "Mars Short roundtrip")

    call cartesian_to_brouwer_mean_long(mu_mars, req_mars, j2_mars, j3_mars, j4_mars, j5_mars, cart_orig, stat=stat, blml=blml)
    call check(stat == 0, "cartesian_to_brouwer_mean_long Mars")
    call brouwer_mean_long_to_cartesian(mu_mars, req_mars, j2_mars, j3_mars, j4_mars, j5_mars, blml, stat=stat, cart=cart_long)
    call check(stat == 0, "brouwer_mean_long_to_cartesian Mars")
    diff = norm2(cart_orig - cart_long) / norm2(cart_orig)
    call check(diff < 1.0e-7_wp, "Mars Long roundtrip")

    ! -------------------------------------------------------------
    ! Test 7: Osculating direct calls with special cases
    ! -------------------------------------------------------------
    ! Circular mean elements (ecc <= 1e-11)
    blms = [7000.0_wp, 0.0_wp, 28.5_wp, 45.0_wp, 30.0_wp, 50.0_wp]
    call brouwer_mean_short_to_osculating(mu_earth, req_earth, j2_earth, blms, stat=stat, kepl=kep_res)
    call check(stat == 0, "Short osculating circular inclined")

    ! Circular equatorial (ecc <= 1e-11, inc <= 1e-7)
    blms = [7000.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 50.0_wp]
    call brouwer_mean_short_to_osculating(mu_earth, req_earth, j2_earth, blms, stat=stat, kepl=kep_res)
    call check(stat == 0, "Short osculating circular equatorial")

    ! Negative eccentricity input
    blms = [7000.0_wp, -0.01_wp, 28.5_wp, 45.0_wp, 30.0_wp, 50.0_wp]
    call brouwer_mean_short_to_osculating(mu_earth, req_earth, j2_earth, blms, stat=stat, kepl=kep_res)
    call check(stat == 0, "Short osculating negative eccentricity")

    ! Long osculating circular inclined & circular equatorial
    blml = [7000.0_wp, 0.0_wp, 28.5_wp, 45.0_wp, 30.0_wp, 50.0_wp]
    call brouwer_mean_long_to_osculating(mu_earth, req_earth, j2_earth, j3_earth, j4_earth, j5_earth, blml, stat=stat, kepl=kep_res)
    call check(stat == 0, "Long osculating circular inclined")

    blml = [7000.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 50.0_wp]
    call brouwer_mean_long_to_osculating(mu_earth, req_earth, j2_earth, j3_earth, j4_earth, j5_earth, blml, stat=stat, kepl=kep_res)
    call check(stat == 0, "Long osculating circular equatorial")

    ! Long osculating inclined with near-equatorial inc <= 1e-7 and ecc > 0
    blml = [7000.0_wp, 0.01_wp, 0.0_wp, 0.0_wp, 30.0_wp, 50.0_wp]
    call brouwer_mean_long_to_osculating(mu_earth, req_earth, j2_earth, j3_earth, j4_earth, j5_earth, blml, stat=stat, kepl=kep_res)
    call check(stat == 0, "Long osculating elliptic equatorial")

    ! -------------------------------------------------------------
    ! Test 8: Cartesian to Keplerian Orbital Regimes
    ! -------------------------------------------------------------
    ! Case 1: Non-circular, Inclined with quadrant flips
    kep_orig = [8000.0_wp, 0.1_wp, 45.0_wp, 200.0_wp, 210.0_wp, 220.0_wp]
    call keplerian_to_cartesian(mu_earth, kep_orig, anomaly_type="TA", stat=stat, cart=cart_test)
    call cartesian_to_keplerian(mu_earth, cart_test, anomaly_type="TA", stat=stat, kepl=kep_res)
    diff = norm2(kep_orig - kep_res) / norm2(kep_orig)
    call check(diff < 1.0e-7_wp, "Cart2Kep Case 1 quadrant flips")

    ! Case 2: Non-circular, Equatorial Prograde (inc = 0)
    kep_orig = [8000.0_wp, 0.1_wp, 0.0_wp, 0.0_wp, 40.0_wp, 50.0_wp]
    call keplerian_to_cartesian(mu_earth, kep_orig, anomaly_type="TA", stat=stat, cart=cart_test)
    call cartesian_to_keplerian(mu_earth, cart_test, anomaly_type="TA", stat=stat, kepl=kep_res)
    call check(abs(kep_res(1) - 8000.0_wp) < 1.0e-4_wp .and. abs(kep_res(2) - 0.1_wp) < 1.0e-6_wp, "Cart2Kep Case 2 prograde")

    ! Case 2b: Non-circular, Equatorial Prograde with pos.vel < 0
    kep_orig = [8000.0_wp, 0.1_wp, 0.0_wp, 0.0_wp, 200.0_wp, 220.0_wp]
    call keplerian_to_cartesian(mu_earth, kep_orig, anomaly_type="TA", stat=stat, cart=cart_test)
    call cartesian_to_keplerian(mu_earth, cart_test, anomaly_type="TA", stat=stat, kepl=kep_res)
    call check(abs(kep_res(1) - 8000.0_wp) < 1.0e-4_wp, "Cart2Kep Case 2 prograde quadrant 2")

    ! Case 2c: Non-circular, Equatorial Retrograde (inc = 180)
    kep_orig = [8000.0_wp, 0.1_wp, 180.0_wp, 0.0_wp, 40.0_wp, 50.0_wp]
    call keplerian_to_cartesian(mu_earth, kep_orig, anomaly_type="TA", stat=stat, cart=cart_test)
    call cartesian_to_keplerian(mu_earth, cart_test, anomaly_type="TA", stat=stat, kepl=kep_res)
    call check(abs(kep_res(1) - 8000.0_wp) < 1.0e-4_wp, "Cart2Kep Case 2 retrograde")

    ! Case 3: Circular, Inclined (e = 0, inc = 45)
    kep_orig = [8000.0_wp, 0.0_wp, 45.0_wp, 30.0_wp, 0.0_wp, 60.0_wp]
    call keplerian_to_cartesian(mu_earth, kep_orig, anomaly_type="TA", stat=stat, cart=cart_test)
    call cartesian_to_keplerian(mu_earth, cart_test, anomaly_type="TA", stat=stat, kepl=kep_res)
    call check(abs(kep_res(1) - 8000.0_wp) < 1.0e-4_wp .and. abs(kep_res(3) - 45.0_wp) < 1.0e-5_wp, "Cart2Kep Case 3")

    ! Case 3b: Circular, Inclined with pos(3) < 0 and nodeVec(2) < 0
    kep_orig = [8000.0_wp, 0.0_wp, 45.0_wp, 200.0_wp, 0.0_wp, 220.0_wp]
    call keplerian_to_cartesian(mu_earth, kep_orig, anomaly_type="TA", stat=stat, cart=cart_test)
    call cartesian_to_keplerian(mu_earth, cart_test, anomaly_type="TA", stat=stat, kepl=kep_res)
    call check(abs(kep_res(1) - 8000.0_wp) < 1.0e-4_wp, "Cart2Kep Case 3b")

    ! Case 4: Circular, Equatorial Prograde (e = 0, inc = 0)
    kep_orig = [8000.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 45.0_wp]
    call keplerian_to_cartesian(mu_earth, kep_orig, anomaly_type="TA", stat=stat, cart=cart_test)
    call cartesian_to_keplerian(mu_earth, cart_test, anomaly_type="TA", stat=stat, kepl=kep_res)
    call check(abs(kep_res(1) - 8000.0_wp) < 1.0e-4_wp, "Cart2Kep Case 4 prograde")

    ! Case 4b: Circular, Equatorial with pos(2) < 0
    kep_orig = [8000.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 220.0_wp]
    call keplerian_to_cartesian(mu_earth, kep_orig, anomaly_type="TA", stat=stat, cart=cart_test)
    call cartesian_to_keplerian(mu_earth, cart_test, anomaly_type="TA", stat=stat, kepl=kep_res)
    call check(abs(kep_res(1) - 8000.0_wp) < 1.0e-4_wp, "Cart2Kep Case 4b")

    ! Case 4c: Circular, Equatorial Retrograde (e = 0, inc = 180)
    kep_orig = [8000.0_wp, 0.0_wp, 180.0_wp, 0.0_wp, 0.0_wp, 45.0_wp]
    call keplerian_to_cartesian(mu_earth, kep_orig, anomaly_type="TA", stat=stat, cart=cart_test)
    call cartesian_to_keplerian(mu_earth, cart_test, anomaly_type="TA", stat=stat, kepl=kep_res)
    call check(abs(kep_res(1) - 8000.0_wp) < 1.0e-4_wp, "Cart2Kep Case 4c retrograde")

    ! Convert to Mean Anomaly mode in CartesianToKeplerian & KeplerianToCartesian
    kep_orig = [7500.0_wp, 0.05_wp, 35.0_wp, 45.0_wp, 60.0_wp, 75.0_wp]
    call keplerian_to_cartesian(mu_earth, kep_orig, anomaly_type="MA", stat=stat, cart=cart_test)
    call cartesian_to_keplerian(mu_earth, cart_test, anomaly_type="MA", stat=stat, kepl=kep_res)
    diff = norm2(kep_orig - kep_res) / norm2(kep_orig)
    call check(diff < 1.0e-7_wp, "Cart2Kep and Kep2Cart with MA anomaly type")

    ! -------------------------------------------------------------
    ! Test 9: Anomaly Conversions (Elliptic, Parabolic, Hyperbolic)
    ! -------------------------------------------------------------
    ! Elliptic: True <-> Eccentric <-> Mean
    ta = 45.0_wp * deg2rad
    ea = true_to_eccentric_anomaly(ta, 0.1_wp)
    ma = true_to_mean_anomaly(ta, 0.1_wp)
    call mean_to_true_anomaly(ma, 0.1_wp, tol=1.0e-12_wp, stat=stat, ta=ta)
    call check(abs(ta - 45.0_wp * deg2rad) < 1.0e-10_wp .and. stat == 0, "Elliptic anomaly roundtrip")

    ! Special Elliptic angle E = pi, M = pi
    ta = pi
    ea = true_to_eccentric_anomaly(ta, 0.2_wp)
    ma = true_to_mean_anomaly(ta, 0.2_wp)
    call mean_to_true_anomaly(ma, 0.2_wp, stat=stat, ta=ta)
    call check(abs(ta - pi) < 1.0e-10_wp, "Anomaly at pi")

    ! Parabolic case
    ma = true_to_mean_anomaly(45.0_wp * deg2rad, 1.0_wp)
    call check(ma == 0.0_wp, "Parabolic true_to_mean_anomaly returns 0")

    ! Hyperbolic case
    ha = true_to_hyperbolic_anomaly(0.5_wp, 1.5_wp)
    call check(ha > 0.0_wp, "Hyperbolic true_to_hyperbolic_anomaly")
    ma = true_to_mean_anomaly(0.5_wp, 1.5_wp)
    call check(ma > 0.0_wp, "Hyperbolic true_to_mean_anomaly")
    call mean_to_true_anomaly(ma, 1.5_wp, stat=stat, ta=ta)
    call check(abs(ta - 0.5_wp) < 1.0e-8_wp .and. stat == 0, "Hyperbolic mean_to_true_anomaly roundtrip")

    ! -------------------------------------------------------------
    ! Test 10: Error Handling / Status Flags
    ! -------------------------------------------------------------
    ! Invalid mu or req <= 0
    call cartesian_to_brouwer_mean_short(-1.0_wp, req_earth, j2_earth, cart_orig, stat=stat, blms=blms)
    call check(stat == 1, "Short error on negative mu")

    call cartesian_to_brouwer_mean_long(mu_earth, -100.0_wp, j2_earth, j3_earth, j4_earth, j5_earth, cart_orig, stat=stat, blml=blml)
    call check(stat == 1, "Long error on negative req")

    ! Invalid mu for Osculating Short & Long
    call brouwer_mean_short_to_osculating(-1.0_wp, req_earth, j2_earth, blms, stat=stat, kepl=kep_res)
    call check(stat == 1, "Short osculating error on invalid mu")
    call brouwer_mean_short_to_cartesian(-1.0_wp, req_earth, j2_earth, blms, stat=stat, cart=cart_test)
    call check(stat == 1, "Short cartesian error on invalid mu")

    call brouwer_mean_long_to_osculating(mu_earth, 0.0_wp, j2_earth, j3_earth, j4_earth, j5_earth, blml, stat=stat, kepl=kep_res)
    call check(stat == 1, "Long osculating error on zero req")
    call brouwer_mean_long_to_cartesian(mu_earth, 0.0_wp, j2_earth, j3_earth, j4_earth, j5_earth, blml, stat=stat, cart=cart_test)
    call check(stat == 1, "Long cartesian error on zero req")

    ! Invalid inclination > 180
    cart_test = [10000.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 6.0_wp, 0.0_wp]
    kep_orig = [7000.0_wp, 0.01_wp, 185.0_wp, 0.0_wp, 0.0_wp, 0.0_wp]
    call brouwer_mean_short_to_osculating(mu_earth, req_earth, j2_earth, kep_orig, stat=stat, kepl=kep_res)
    call check(stat == 2, "Short osculating error on inc > 180")

    call brouwer_mean_long_to_osculating(mu_earth, req_earth, j2_earth, j3_earth, j4_earth, j5_earth, kep_orig, stat=stat, kepl=kep_res)
    call check(stat == 2, "Long osculating error on inc > 180")

    ! Invalid eccentricity >= 0.99
    kep_orig = [1000000.0_wp, 0.995_wp, 28.5_wp, 0.0_wp, 0.0_wp, 0.0_wp]
    call brouwer_mean_short_to_osculating(mu_earth, req_earth, j2_earth, kep_orig, stat=stat, kepl=kep_res)
    call check(stat == 3, "Short osculating error on ecc >= 0.99")

    call brouwer_mean_long_to_osculating(mu_earth, req_earth, j2_earth, j3_earth, j4_earth, j5_earth, kep_orig, stat=stat, kepl=kep_res)
    call check(stat == 3, "Long osculating error on ecc >= 0.99")

    ! Invalid periapsis radius <= 0 km
    kep_orig = [-100.0_wp, 0.2_wp, 28.5_wp, 0.0_wp, 0.0_wp, 0.0_wp]
    call brouwer_mean_short_to_osculating(mu_earth, req_earth, j2_earth, kep_orig, stat=stat, kepl=kep_res)
    call check(stat == 4, "Short osculating error on radper <= 0 km")

    call brouwer_mean_long_to_osculating(mu_earth, req_earth, j2_earth, j3_earth, j4_earth, j5_earth, kep_orig, stat=stat, kepl=kep_res)
    call check(stat == 4, "Long osculating error on radper <= 0 km")

    ! Cartesian to Brouwer errors on ecc >= 0.99
    call keplerian_to_cartesian(mu_earth, [1000000.0_wp, 0.995_wp, 28.5_wp, 0.0_wp, 0.0_wp, 0.0_wp], stat=stat, cart=cart_test)
    call cartesian_to_brouwer_mean_short(mu_earth, req_earth, j2_earth, cart_test, stat=stat, blms=blms)
    call check(stat == 3, "Cart2BrouwerShort error on ecc >= 0.99")
    call cartesian_to_brouwer_mean_long(mu_earth, req_earth, j2_earth, j3_earth, j4_earth, j5_earth, cart_test, stat=stat, blml=blml)
    call check(stat == 3, "Cart2BrouwerLong error on ecc >= 0.99")

    ! Zero / invalid Cartesian inputs
    cart_test = [0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp]
    call cartesian_to_keplerian(mu_earth, cart_test, stat=stat, kepl=kep_res)
    call check(stat == 6, "Cart2Kep error on zero pos/vel")

    call cartesian_to_keplerian(0.0_wp, cart_orig, stat=stat, kepl=kep_res)
    call check(stat == 6, "Cart2Kep error on zero mu")

    ! Singular radius in Keplerian to Cartesian
    kep_res = [0.0_wp, 1.0_wp, 28.5_wp, 0.0_wp, 0.0_wp, 0.0_wp]
    call keplerian_to_cartesian(mu_earth, kep_res, stat=stat, cart=cart_test)
    call check(stat == 6, "Kep2Cart error on singular p")

    print *, "=========================================================="
    if (.not. all_passed) error stop 'TEST FAILURE: One or more tests failed.'
    print *, " ALL TESTS PASSED SUCCESSFULLY!"
    print *, "=========================================================="

contains

    subroutine check(condition, test_name)
        !! Check a condition and print the result of a test.
        logical, intent(in) :: condition
        character(len=*), intent(in) :: test_name
        if (condition) then
            print *, "PASS: ", test_name
        else
            print *, "FAIL: ", test_name
            all_passed = .false.
            error stop
        end if
    end subroutine check

end program main_test
