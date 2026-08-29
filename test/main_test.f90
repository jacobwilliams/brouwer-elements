program main_test

    use iso_fortran_env, only: wp => real64
    use brouwer_module

    implicit none

    real(wp), dimension(6) :: cart_orig, blms, cart_roundtrip, blms_roundtrip
    real(wp), dimension(6) :: blml, cart_long, blml_roundtrip
    real(wp), dimension(6) :: kep_orig, cart_from_kep
    integer :: stat
    real(wp) :: diff
    logical :: all_passed

    all_passed = .true.

    print *, "=========================================================="
    print *, " Testing Modern Fortran Brouwer Element Conversions"
    print *, "=========================================================="

    ! Define an Earth orbit test case (e.g. LEO/MEO satellite)
    ! sma = 7000 km, ecc = 0.01, inc = 28.5 deg, raan = 45 deg, aop = 30 deg, ta = 15 deg
    kep_orig = [7000.0_wp, 0.01_wp, 28.5_wp, 45.0_wp, 30.0_wp, 15.0_wp]
    cart_orig = keplerian_to_cartesian(mu_earth, kep_orig, anomaly_type="TA", stat=stat)
    print *, "cart_orig stat:", stat
    if (stat /= 0) then
        print *, "FAIL: keplerian_to_cartesian returned error:", stat
        all_passed = .false.
    end if

    print *, "Initial Cartesian state:"
    print '(6(1X,F15.6))', cart_orig

    ! -------------------------------------------------------------
    ! Test 1: Cartesian -> Brouwer Mean Short -> Cartesian roundtrip
    ! -------------------------------------------------------------
    blms = cartesian_to_brouwer_mean_short(mu_earth, cart_orig, stat=stat)
    if (stat /= 0) then
        print *, "FAIL: cartesian_to_brouwer_mean_short returned error:", stat
        all_passed = .false.
    end if

    print *, "Brouwer Mean Short elements:"
    print '(6(1X,F15.6))', blms

    cart_roundtrip = brouwer_mean_short_to_cartesian(mu_earth, blms, stat=stat)
    if (stat /= 0) then
        print *, "FAIL: brouwer_mean_short_to_cartesian returned error:", stat
        all_passed = .false.
    end if

    diff = norm2(cart_orig - cart_roundtrip)
    print *, "Brouwer Short roundtrip position/velocity norm difference:", diff
    if (diff < 1.0e-5_wp) then
        print *, "PASS: Brouwer Mean Short roundtrip successful!"
    else
        print *, "FAIL: Brouwer Mean Short roundtrip difference too high!"
        all_passed = .false.
    end if

    ! -------------------------------------------------------------
    ! Test 2: Cartesian -> Brouwer Mean Long -> Cartesian roundtrip
    ! -------------------------------------------------------------
    blml = cartesian_to_brouwer_mean_long(mu_earth, cart_orig, stat=stat)
    print *, "cartesian_to_brouwer_mean_long returned stat =", stat
    if (stat /= 0) then
        print *, "FAIL: cartesian_to_brouwer_mean_long returned error:", stat
        all_passed = .false.
    end if

    print *, "Brouwer Mean Long elements:"
    print '(6(1X,F15.6))', blml

    cart_long = brouwer_mean_long_to_cartesian(mu_earth, blml, stat=stat)
    if (stat /= 0) then
        print *, "FAIL: brouwer_mean_long_to_cartesian returned error:", stat
        all_passed = .false.
    end if

    diff = norm2(cart_orig - cart_long)
    print *, "Brouwer Long roundtrip absolute diff (km):", diff
    print *, "Brouwer Long roundtrip relative diff:", diff / norm2(cart_orig)
    if (diff / norm2(cart_orig) < 1.0e-7_wp) then
        print *, "PASS: Brouwer Mean Long roundtrip successful!"
    else
        print *, "FAIL: Brouwer Mean Long roundtrip difference too high!"
        all_passed = .false.
    end if

    ! -------------------------------------------------------------
    ! Test 3: Brouwer Elements -> Cartesian -> Brouwer Elements
    ! -------------------------------------------------------------
    blms_roundtrip = cartesian_to_brouwer_mean_short(mu_earth, cart_roundtrip, stat=stat)
    diff = norm2(blms - blms_roundtrip)
    print *, "Brouwer Mean Short elements roundtrip diff:", diff
    print *, "Brouwer Mean Short elements relative diff:", diff / norm2(blms)
    if (diff / norm2(blms) < 1.0e-7_wp) then
        print *, "PASS: Element roundtrip for short-period successful!"
    else
        print *, "FAIL: Element roundtrip diff too high!"
        all_passed = .false.
    end if

    blml_roundtrip = cartesian_to_brouwer_mean_long(mu_earth, cart_long, stat=stat)
    diff = norm2(blml - blml_roundtrip)
    print *, "Brouwer Mean Long elements roundtrip diff:", diff
    print *, "Brouwer Mean Long elements relative diff:", diff / norm2(blml)
    if (diff / norm2(blml) < 1.0e-7_wp) then
        print *, "PASS: Element roundtrip for long-period successful!"
    else
        print *, "FAIL: Element roundtrip diff too high!"
        all_passed = .false.
    end if

    ! -------------------------------------------------------------
    ! Test 4: Polar / High inclination orbit test
    ! -------------------------------------------------------------
    kep_orig = [7200.0_wp, 0.005_wp, 87.0_wp, 120.0_wp, 45.0_wp, 200.0_wp]
    cart_orig = keplerian_to_cartesian(mu_earth, kep_orig, anomaly_type="TA", stat=stat)
    blms = cartesian_to_brouwer_mean_short(mu_earth, cart_orig, stat=stat)
    cart_roundtrip = brouwer_mean_short_to_cartesian(mu_earth, blms, stat=stat)
    diff = norm2(cart_orig - cart_roundtrip) / norm2(cart_orig)
    if (diff < 1.0e-7_wp) then
        print *, "PASS: High inclination Short roundtrip successful! (rel diff =", diff, ")"
    else
        print *, "FAIL: High inclination Short roundtrip failed! (rel diff =", diff, ")"
        all_passed = .false.
    end if

    blml = cartesian_to_brouwer_mean_long(mu_earth, cart_orig, stat=stat)
    cart_long = brouwer_mean_long_to_cartesian(mu_earth, blml, stat=stat)
    diff = norm2(cart_orig - cart_long) / norm2(cart_orig)
    if (diff < 1.0e-7_wp) then
        print *, "PASS: High inclination Long roundtrip successful! (rel diff =", diff, ")"
    else
        print *, "FAIL: High inclination Long roundtrip failed! (rel diff =", diff, ")"
        all_passed = .false.
    end if

    ! -------------------------------------------------------------
    ! Test 5: Retrograde orbit test (inc > 90 deg)
    ! -------------------------------------------------------------
    kep_orig = [8000.0_wp, 0.05_wp, 105.0_wp, 210.0_wp, 75.0_wp, 315.0_wp]
    cart_orig = keplerian_to_cartesian(mu_earth, kep_orig, anomaly_type="TA", stat=stat)
    blms = cartesian_to_brouwer_mean_short(mu_earth, cart_orig, stat=stat)
    cart_roundtrip = brouwer_mean_short_to_cartesian(mu_earth, blms, stat=stat)
    diff = norm2(cart_orig - cart_roundtrip) / norm2(cart_orig)
    if (diff < 1.0e-7_wp) then
        print *, "PASS: Retrograde Short roundtrip successful! (rel diff =", diff, ")"
    else
        print *, "FAIL: Retrograde Short roundtrip failed! (rel diff =", diff, ")"
        all_passed = .false.
    end if

    blml = cartesian_to_brouwer_mean_long(mu_earth, cart_orig, stat=stat)
    cart_long = brouwer_mean_long_to_cartesian(mu_earth, blml, stat=stat)
    diff = norm2(cart_orig - cart_long) / norm2(cart_orig)
    if (diff < 1.0e-7_wp) then
        print *, "PASS: Retrograde Long roundtrip successful! (rel diff =", diff, ")"
    else
        print *, "FAIL: Retrograde Long roundtrip failed! (rel diff =", diff, ")"
        all_passed = .false.
    end if

    print *, "=========================================================="
    if (all_passed) then
        print *, " ALL TESTS PASSED SUCCESSFULLY!"
    else
        stop 1
    end if
    print *, "=========================================================="

end program main_test
