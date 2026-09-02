!------------------------------------------------------------------------------
!>
!  Modern Fortran implementation of Brouwer-Lyddane Mean Elements.
!
!  Based on the GMAT code `StateConversionUtil.cpp`.
!
!  Long-period terms use Phipps' T2 regularization to avoid the
!  critical-inclination singularity at approximately 63.4 degrees.
!
!## References
!
!  * Brouwer, D., "Solution of the Problem of Artificial Satellite Theory without Drag,"
!    *Astronomical Journal*, Vol. 64, Nov. 1959, pp. 378-397.
!  * Lyddane, R. H., "Small Eccentricities or Inclinations in the Brouwer
!    Theory of the Artificial Satellite," *Astronomical Journal*,
!    Vol. 68, Oct. 1963, pp. 555-558.
!  * NASA General Mission Analysis Tool (GMAT), `StateConversionUtil.cpp`.
!  * Phipps Jr., W., "Parallelization of the Navy Space Surveillance Center
!    (NAVSPASUR) Satellite Model," Naval Postgraduate School, 1992.

module brouwer_module

    use, intrinsic :: iso_fortran_env

    implicit none

    private

#ifdef REAL32
    integer,parameter,public :: brouwer_module_wp = real32   !! Real working precision [4 bytes]
#elif REAL64
    integer,parameter,public :: brouwer_module_wp = real64   !! Real working precision [8 bytes]
#elif REAL128
    integer,parameter,public :: brouwer_module_wp = real128  !! Real working precision [16 bytes]
#else
    integer,parameter,public :: brouwer_module_wp = real64   !! Real working precision if not specified [8 bytes]
#endif
    integer,parameter :: wp = brouwer_module_wp !! real kind to use in this module

    ! Error codes:
    integer,parameter,public :: BROUWER_SUCCESS                    = 0 !! Successful execution
    integer,parameter,public :: BROUWER_INVALID_MU                 = 1 !! Invalid gravitational parameter or equatorial radius
    integer,parameter,public :: BROUWER_INVALID_INCLINATION        = 2 !! Invalid inclination
    integer,parameter,public :: BROUWER_INVALID_ECCENTRICITY       = 3 !! Eccentricity outside [0, 0.99)
    integer,parameter,public :: BROUWER_PERIAPSIS_TOO_LOW          = 4 !! Periapsis radius < 1 km
    integer,parameter,public :: BROUWER_ITERATION_DID_NOT_CONVERGE = 5 !! Iteration did not converge
    integer,parameter,public :: BROUWER_DEGENERATE_STATE           = 6 !! Degenerate Cartesian state (e.g., zero position and velocity)
    integer,parameter,public :: BROUWER_PARABOLIC_ORBIT            = 7 !! Parabolic orbit (eccentricity close to 1)
    integer,parameter,public :: BROUWER_SINGULARITY                = 8 !! some other error (in mean_to_true_anomaly)

    ! Mathematical constants
    real(wp), parameter :: pi = acos(-1.0_wp)
    real(wp), parameter :: two_pi = 2.0_wp * pi
    real(wp), parameter :: deg2rad = pi / 180.0_wp
    real(wp), parameter :: rad2deg = 180.0_wp / pi
    real(wp), parameter :: min_brouwer_radper = 1.0_wp !! Minimum periapsis radius for Brouwer-Lyddane Mean Elements (km)

    ! Numerical tolerances
    real(wp), parameter :: parabolic_tol = 1.0e-10_wp !! Tolerance for Keplerian elements

    ! Public API
    public :: cartesian_to_brouwer_mean_short
    public :: brouwer_mean_short_to_cartesian
    public :: cartesian_to_brouwer_mean_long
    public :: brouwer_mean_long_to_cartesian

    ! Supporting routines
    public :: brouwer_mean_short_to_osculating
    public :: brouwer_mean_long_to_osculating
    public :: cartesian_to_keplerian
    public :: keplerian_to_cartesian
    public :: true_to_mean_anomaly
    public :: mean_to_true_anomaly
    public :: true_to_eccentric_anomaly
    public :: true_to_hyperbolic_anomaly

contains

    !--------------------------------------------------------------------------
    !>
    !  Converts Cartesian state to Brouwer-Lyddane Mean Elements (short-period terms only).

    pure subroutine cartesian_to_brouwer_mean_short(mu, req, j2, cartesian, stat, blms)
        real(wp), intent(in) :: mu !! Central body gravitational parameter (km^3/s^2)
        real(wp), intent(in) :: req !! Central body equatorial radius (km)
        real(wp), intent(in) :: j2 !! Central body J2 zonal harmonic coefficient
        real(wp), dimension(6), intent(in) :: cartesian !! Cartesian state vector [x, y, z, vx, vy, vz] (km, km/s)
        integer, intent(out) :: stat !! Status: 0 success; /=0 failure.
        real(wp), dimension(6), intent(out) :: blms !! Brouwer mean elements [sma(km), ecc, inc(deg), raan(deg), aop(deg), ma(deg)]

        real(wp), parameter :: tol = 1.0e-8_wp
        integer, parameter :: maxiter = 75

        real(wp), dimension(6) :: cart, kep, kep2, blmean, blmean2, &
                                  aeq, aeq2, aeqmean, aeqmean2, tmp, cart2
        real(wp) :: radper, emag, emag_old, sum_sq_diff, sum_sq_cart, inc_arg
        integer :: pseudostate, ii

        stat = BROUWER_SUCCESS
        blms = 0.0_wp

        if (mu <= 0.0_wp .or. req <= 0.0_wp) then
            stat = BROUWER_INVALID_MU
            return
        end if

        cart = cartesian
        call cartesian_to_keplerian(mu, cart, anomaly_type="TA", stat=stat, kepl=kep)
        if (stat /= BROUWER_SUCCESS) return

        if (kep(3) > 180.0_wp) then
            stat = BROUWER_INVALID_INCLINATION
            return
        end if

        if (kep(2) >= 0.99_wp .or. kep(2) < 0.0_wp) then
            stat = BROUWER_INVALID_ECCENTRICITY
            return
        end if

        radper = kep(1) * (1.0_wp - kep(2))
        if (radper < min_brouwer_radper) then
            stat = BROUWER_PERIAPSIS_TOO_LOW
            return
        end if

        ! Convert true anomaly (deg) to mean anomaly (deg)
        kep(6) = true_to_mean_anomaly(kep(6) * deg2rad, kep(2)) * rad2deg

        pseudostate = 0
        if (kep(3) > 175.0_wp) then
            kep(3) = 180.0_wp - kep(3)
            kep(4) = -kep(4)
            call keplerian_to_cartesian(mu, kep, anomaly_type="MA", stat=stat, cart=cart)
            if (stat /= BROUWER_SUCCESS) return
            pseudostate = 1
        end if

        blmean = kep
        call brouwer_mean_short_to_osculating(mu, req, j2, kep, stat=stat, kepl=kep2)
        if (stat /= BROUWER_SUCCESS) return

        ! Convert to alternate equinoctial elements
        aeq(1) = kep(1)
        aeq(2) = kep(2) * sin((kep(5) + kep(4)) * deg2rad)
        aeq(3) = kep(2) * cos((kep(5) + kep(4)) * deg2rad)
        aeq(4) = sin(kep(3) * 0.5_wp * deg2rad) * sin(kep(4) * deg2rad)
        aeq(5) = sin(kep(3) * 0.5_wp * deg2rad) * cos(kep(4) * deg2rad)
        aeq(6) = kep(4) + kep(5) + kep(6)

        aeq2(1) = kep2(1)
        aeq2(2) = kep2(2) * sin((kep2(5) + kep2(4)) * deg2rad)
        aeq2(3) = kep2(2) * cos((kep2(5) + kep2(4)) * deg2rad)
        aeq2(4) = sin(kep2(3) * 0.5_wp * deg2rad) * sin(kep2(4) * deg2rad)
        aeq2(5) = sin(kep2(3) * 0.5_wp * deg2rad) * cos(kep2(4) * deg2rad)
        aeq2(6) = kep2(4) + kep2(5) + kep2(6)

        aeqmean(1) = blmean(1)
        aeqmean(2) = blmean(2) * sin((blmean(5) + blmean(4)) * deg2rad)
        aeqmean(3) = blmean(2) * cos((blmean(5) + blmean(4)) * deg2rad)
        aeqmean(4) = sin(blmean(3) * 0.5_wp * deg2rad) * sin(blmean(4) * deg2rad)
        aeqmean(5) = sin(blmean(3) * 0.5_wp * deg2rad) * cos(blmean(4) * deg2rad)
        aeqmean(6) = blmean(4) + blmean(5) + blmean(6)

        aeqmean2 = aeqmean + (aeq - aeq2)

        emag = 0.9_wp
        emag_old = 1.0_wp
        ii = 0

        do while (emag > tol)
            blmean2(1) = aeqmean2(1)
            blmean2(2) = sqrt(aeqmean2(2)**2 + aeqmean2(3)**2)

            inc_arg = aeqmean2(4)**2 + aeqmean2(5)**2
            blmean2(3) = acos(1.0_wp - 2.0_wp * min(1.0_wp, inc_arg)) * rad2deg

            blmean2(4) = atan2(aeqmean2(4), aeqmean2(5)) * rad2deg
            if (blmean2(4) < 0.0_wp) blmean2(4) = blmean2(4) + 360.0_wp

            blmean2(5) = atan2(aeqmean2(2), aeqmean2(3)) * rad2deg - blmean2(4)
            if (blmean2(5) < 0.0_wp) blmean2(5) = blmean2(5) + 360.0_wp

            blmean2(6) = aeqmean2(6) - atan2(aeqmean2(2), aeqmean2(3)) * rad2deg

            call brouwer_mean_short_to_osculating(mu, req, j2, blmean2, stat=stat, kepl=kep2)
            if (stat /= BROUWER_SUCCESS) return
            call keplerian_to_cartesian(mu, kep2, anomaly_type="MA", stat=stat, cart=cart2)
            if (stat /= BROUWER_SUCCESS) return

            tmp = cart - cart2
            sum_sq_diff = sum(tmp**2)
            sum_sq_cart = sum(cart**2)
            emag = sqrt(sum_sq_diff) / sqrt(sum_sq_cart)

            aeq2(1) = kep2(1)
            aeq2(2) = kep2(2) * sin((kep2(5) + kep2(4)) * deg2rad)
            aeq2(3) = kep2(2) * cos((kep2(5) + kep2(4)) * deg2rad)
            aeq2(4) = sin(kep2(3) * 0.5_wp * deg2rad) * sin(kep2(4) * deg2rad)
            aeq2(5) = sin(kep2(3) * 0.5_wp * deg2rad) * cos(kep2(4) * deg2rad)
            aeq2(6) = kep2(4) + kep2(5) + kep2(6)

            if (emag_old > emag) then
                emag_old = emag
                aeqmean = aeqmean2
                aeqmean2 = aeqmean + (aeq - aeq2)
            else
                ! Not converging
                stat = BROUWER_ITERATION_DID_NOT_CONVERGE
                exit
            end if

            ii = ii + 1
            if (ii > maxiter) then
                stat = BROUWER_ITERATION_DID_NOT_CONVERGE
                exit
            end if
        end do

        blmean(1) = aeqmean2(1)
        blmean(2) = sqrt(aeqmean2(2)**2 + aeqmean2(3)**2)

        inc_arg = aeqmean2(4)**2 + aeqmean2(5)**2
        blmean(3) = acos(1.0_wp - 2.0_wp * min(1.0_wp, inc_arg)) * rad2deg

        blmean(4) = atan2(aeqmean2(4), aeqmean2(5)) * rad2deg
        blmean(5) = atan2(aeqmean2(2), aeqmean2(3)) * rad2deg - blmean(4)
        blmean(6) = aeqmean2(6) - atan2(aeqmean2(2), aeqmean2(3)) * rad2deg

        if (pseudostate /= 0) then
            blmean(3) = 180.0_wp - blmean(3)
            blmean(4) = -blmean(4)
        end if

        call wrap_0_360(blmean(4))
        call wrap_0_360(blmean(5))
        call wrap_0_360(blmean(6))

        blms = blmean

    end subroutine cartesian_to_brouwer_mean_short

    !--------------------------------------------------------------------------
    !>
    !  Converts Brouwer-Lyddane Mean Elements (short-period terms only) to
    !  Osculating Keplerian Elements.

    pure subroutine brouwer_mean_short_to_osculating(mu, req, j2, blms, stat, kepl)
        real(wp), intent(in) :: mu !! Central body gravitational parameter (km^3/s^2)
        real(wp), intent(in) :: req !! Central body equatorial radius (km)
        real(wp), intent(in) :: j2 !! Central body J2 zonal harmonic coefficient
        real(wp), dimension(6), intent(in) :: blms !! Brouwer mean elements [sma(km), ecc, inc(deg), raan(deg), aop(deg), ma(deg)]
        integer, intent(out) :: stat !! Status: 0 success; /=0 failure.
        real(wp), dimension(6), intent(out) :: kepl !! Osculating Keplerian elements [sma(km), ecc, inc(deg), raan(deg), aop(deg), ma(deg)]

        real(wp) :: smap, eccp, incp, raanp, aopp, meanAnom, radper
        real(wp) :: eta, theta, p, k2, gm2, gm2p, tap, rp, adr
        real(wp) :: sma1, decc, dinc, draan, aop1, ma1, lgh, eccpdl
        real(wp) :: ecosl, esinl, ecc1, sinhalfisinh, sinhalficosh, inc1, raan1, sqr_inc
        integer :: pseudostate

        stat = BROUWER_SUCCESS
        kepl = 0.0_wp

        if (mu <= 0.0_wp .or. req <= 0.0_wp) then
            stat = BROUWER_INVALID_MU
            return
        end if

        smap = blms(1) / req
        eccp = blms(2)
        incp = blms(3) * deg2rad
        raanp = blms(4) * deg2rad
        aopp = blms(5) * deg2rad
        meanAnom = blms(6) * deg2rad

        if (incp < 0.0_wp .or. incp > pi) then
            stat = BROUWER_INVALID_INCLINATION
            return
        end if

        radper = blms(1) * (1.0_wp - blms(2))
        if (radper < min_brouwer_radper) then
            stat = BROUWER_PERIAPSIS_TOO_LOW
            return
        end if

        if (eccp < 0.0_wp) then
            eccp = -eccp
            meanAnom = meanAnom - pi
            aopp = aopp + pi
        end if

        if (eccp > 0.99_wp) then
            stat = BROUWER_INVALID_ECCENTRICITY
            return
        end if

        pseudostate = 0
        if (incp > 175.0_wp * deg2rad) then
            incp = pi - incp
            raanp = -raanp
            pseudostate = 1
        end if

        call wrap_0_2pi(raanp)
        call wrap_0_2pi(aopp)
        call wrap_0_2pi(meanAnom)

        eta = sqrt(max(0.0_wp, 1.0_wp - eccp**2))
        theta = cos(incp)
        p = smap * eta**2
        k2 = 0.5_wp * j2
        gm2 = k2 / (smap**2)
        gm2p = gm2 / (eta**4)

        call mean_to_true_anomaly(meanAnom, eccp, 1.0e-8_wp, stat=stat, ta=tap)
        if (stat /= BROUWER_SUCCESS) return
        if (tap < 0.0_wp) tap = tap + two_pi

        rp = p / (1.0_wp + eccp * cos(tap))
        adr = smap / rp

        sma1 = smap + smap * gm2 * ((adr**3 - 1.0_wp / (eta**3)) * (-1.0_wp + 3.0_wp * theta**2) &
               + 3.0_wp * (1.0_wp - theta**2) * adr**3 * cos(2.0_wp * aopp + 2.0_wp * tap))

        decc = (eta**2 / 2.0_wp) * ((3.0_wp / (eta**6) * gm2 * (1.0_wp - theta**2) * cos(2.0_wp * aopp + 2.0_wp * tap) &
               * (3.0_wp * eccp * cos(tap)**2 + 3.0_wp * cos(tap) + (eccp**2) * (cos(tap)**3) + eccp)) &
               - gm2p * (1.0_wp - theta**2) * (3.0_wp * cos(2.0_wp * aopp + tap) + cos(3.0_wp * tap + 2.0_wp * aopp)) &
               + (3.0_wp * theta**2 - 1.0_wp) * gm2 / (eta**6) * (eccp * eta + eccp / (1.0_wp + eta) &
               + 3.0_wp * eccp * cos(tap)**2 + 3.0_wp * cos(tap) + (eccp**2) * (cos(tap)**3)))

        dinc = (gm2p / 2.0_wp) * theta * sin(incp) * (3.0_wp * cos(2.0_wp * aopp + 2.0_wp * tap) &
               + 3.0_wp * eccp * cos(2.0_wp * aopp + tap) + eccp * cos(2.0_wp * aopp + 3.0_wp * tap))

        draan = -(gm2p / 2.0_wp) * theta * (6.0_wp * (tap - meanAnom + eccp * sin(tap)) &
                - 3.0_wp * sin(2.0_wp * aopp + 2.0_wp * tap) - 3.0_wp * eccp * sin(2.0_wp * aopp + tap) &
                - eccp * sin(2.0_wp * aopp + 3.0_wp * tap))

        if (eccp > 1.0e-11_wp) then ! Avoid singularity at zero eccentricity
            aop1 = aopp + 3.0_wp * j2 / (2.0_wp * p**2) * ((2.0_wp - 2.5_wp * sin(incp)**2) * (tap - meanAnom + eccp * sin(tap)) &
               + (1.0_wp - 1.5_wp * sin(incp)**2) * ((1.0_wp / eccp) * (1.0_wp - 0.25_wp * eccp**2) * sin(tap) &
               + 0.5_wp * sin(2.0_wp * tap) + (eccp / 12.0_wp) * sin(3.0_wp * tap)) &
               - (1.0_wp / eccp) * (0.25_wp * sin(incp)**2 + (0.5_wp - (15.0_wp / 16.0_wp) * sin(incp)**2) * eccp**2) &
               * sin(tap + 2.0_wp * aopp) + (eccp / 16.0_wp) * sin(incp)**2 * sin(tap - 2.0_wp * aopp) &
               - 0.5_wp * (1.0_wp - 2.5_wp * sin(incp)**2) * sin(2.0_wp * tap + 2.0_wp * aopp) &
               + (1.0_wp / eccp) * ((7.0_wp / 12.0_wp) * sin(incp)**2 - (1.0_wp / 6.0_wp) &
               * (1.0_wp - (19.0_wp / 8.0_wp) * sin(incp)**2) * eccp**2) * sin(3.0_wp * tap + 2.0_wp * aopp) &
                + (3.0_wp / 8.0_wp) * sin(incp)**2 * sin(4.0_wp * tap + 2.0_wp * aopp) &
                + (eccp / 16.0_wp) * sin(incp)**2 * sin(5.0_wp * tap + 2.0_wp * aopp))

            ma1 = meanAnom + 3.0_wp * j2 * eta / (2.0_wp * eccp * p**2) * (-(1.0_wp - 1.5_wp * sin(incp)**2) &
              * ((1.0_wp - 0.25_wp * eccp**2) * sin(tap) + (eccp / 2.0_wp) * sin(2.0_wp * tap) + (eccp**2 / 12.0_wp) * sin(3.0_wp * tap)) &
              + sin(incp)**2 * (0.25_wp * (1.0_wp + 1.25_wp * eccp**2) * sin(tap + 2.0_wp * aopp) &
              - (eccp**2 / 16.0_wp) * sin(tap - 2.0_wp * aopp) - (7.0_wp / 12.0_wp) * (1.0_wp - (eccp**2 / 28.0_wp)) &
              * sin(3.0_wp * tap + 2.0_wp * aopp) - (3.0_wp * eccp / 8.0_wp) * sin(4.0_wp * tap + 2.0_wp * aopp) &
              - (eccp**2 / 16.0_wp) * sin(5.0_wp * tap + 2.0_wp * aopp)))
        else
            aop1 = 0.0_wp
            ma1 = 0.0_wp
        end if

        lgh = raanp + aopp + meanAnom + (gm2p / 4.0_wp) * (6.0_wp * (-1.0_wp - 2.0_wp * theta + 5.0_wp * theta**2) &
              * (tap - meanAnom + eccp * sin(tap)) + (3.0_wp + 2.0_wp * theta - 5.0_wp * theta**2) &
              * (3.0_wp * sin(2.0_wp * aopp + 2.0_wp * tap) + 3.0_wp * eccp * sin(2.0_wp * aopp + tap) &
              + eccp * sin(2.0_wp * aopp + 3.0_wp * tap))) + (gm2p / 4.0_wp) * (eta**2 / (eta + 1.0_wp)) * eccp &
              * (3.0_wp * (1.0_wp - theta**2) * (sin(3.0_wp * tap + 2.0_wp * aopp) * ((1.0_wp / 3.0_wp) + adr**2 * eta**2 + adr) &
              + sin(2.0_wp * aopp + tap) * (1.0_wp - adr**2 * eta**2 - adr)) &
              + 2.0_wp * sin(tap) * (3.0_wp * theta**2 - 1.0_wp) * (1.0_wp + adr**2 * eta**2 + adr))

        eccpdl = -(eta**3 / 4.0_wp) * gm2p * (2.0_wp * (-1.0_wp + 3.0_wp * theta**2) * (adr**2 * eta**2 + adr + 1.0_wp) * sin(tap) &
                 + 3.0_wp * (1.0_wp - theta**2) * ((-adr**2 * eta**2 - adr + 1.0_wp) * sin(2.0_wp * aopp + tap) &
                 + (adr**2 * eta**2 + adr + (1.0_wp / 3.0_wp)) * sin(2.0_wp * aopp + 3.0_wp * tap)))

        ecosl = (eccp + decc) * cos(meanAnom) - eccpdl * sin(meanAnom)
        esinl = (eccp + decc) * sin(meanAnom) + eccpdl * cos(meanAnom)
        ecc1 = sqrt(ecosl**2 + esinl**2)

        if (ecc1 < 1.0e-11_wp) then
            ma1 = 0.0_wp
        else
            ma1 = atan2(esinl, ecosl)
            if (ma1 < 0.0_wp) ma1 = ma1 + two_pi
        end if

        sinhalfisinh = (sin(0.5_wp * incp) + cos(0.5_wp * incp) * 0.5_wp * dinc) * sin(raanp) &
                       + 0.5_wp * (sin(incp) / cos(incp * 0.5_wp)) * draan * cos(raanp)
        sinhalficosh = (sin(0.5_wp * incp) + cos(0.5_wp * incp) * 0.5_wp * dinc) * cos(raanp) &
                       - 0.5_wp * (sin(incp) / cos(incp * 0.5_wp)) * draan * sin(raanp)

        sqr_inc = sqrt(sinhalfisinh**2 + sinhalficosh**2)
        inc1 = 2.0_wp * asin(min(1.0_wp, sqr_inc))

        if (inc1 == 0.0_wp .or. abs(inc1 - pi) < 1.0e-14_wp) then
            raan1 = 0.0_wp
            aop1 = lgh - ma1 - raan1
        else
            raan1 = atan2(sinhalfisinh, sinhalficosh)
            if (raan1 < 0.0_wp) raan1 = raan1 + two_pi
            aop1 = lgh - ma1 - raan1
        end if

        call wrap_0_2pi(aop1)
        call wrap_0_2pi(raan1)

        kepl(1) = sma1 * req
        kepl(2) = ecc1
        kepl(3) = inc1 * rad2deg
        kepl(4) = raan1 * rad2deg
        kepl(5) = aop1 * rad2deg
        kepl(6) = ma1 * rad2deg

        if (pseudostate /= 0) then
            kepl(3) = 180.0_wp - kepl(3)
            kepl(4) = 360.0_wp - kepl(4)
        end if

    end subroutine brouwer_mean_short_to_osculating

    !--------------------------------------------------------------------------
    !>
    !  Converts Brouwer-Lyddane Mean Elements (short-period terms only) to Cartesian.

    pure subroutine brouwer_mean_short_to_cartesian(mu, req, j2, blms, stat, cart)
        real(wp), intent(in) :: mu !! Central body gravitational parameter (km^3/s^2)
        real(wp), intent(in) :: req !! Central body equatorial radius (km)
        real(wp), intent(in) :: j2 !! Central body J2 zonal harmonic coefficient
        real(wp), dimension(6), intent(in) :: blms !! Brouwer mean elements [sma(km), ecc, inc(deg), raan(deg), aop(deg), ma(deg)]
        integer, intent(out) :: stat !! Status propagated from [[brouwer_mean_short_to_osculating]] or [[keplerian_to_cartesian]] conversion
        real(wp), dimension(6), intent(out) :: cart !! Cartesian state vector [x, y, z, vx, vy, vz] (km, km/s)

        real(wp), dimension(6) :: kepl

        call brouwer_mean_short_to_osculating(mu, req, j2, blms, stat, kepl)
        if (stat /= BROUWER_SUCCESS) then
            cart = 0.0_wp
            return
        end if

        call keplerian_to_cartesian(mu, kepl, anomaly_type="MA", stat=stat, cart=cart)

    end subroutine brouwer_mean_short_to_cartesian

    !--------------------------------------------------------------------------
    !>
    !  Converts Cartesian state to Brouwer-Lyddane Mean Elements (short and long period terms).

    pure subroutine cartesian_to_brouwer_mean_long(mu, req, j2, j3, j4, j5, cartesian, stat, blml)
        real(wp), intent(in) :: mu !! Central body gravitational parameter (km^3/s^2)
        real(wp), intent(in) :: req !! Central body equatorial radius (km)
        real(wp), intent(in) :: j2 !! Central body J2 zonal harmonic coefficient
        real(wp), intent(in) :: j3 !! Central body J3 zonal harmonic coefficient
        real(wp), intent(in) :: j4 !! Central body J4 zonal harmonic coefficient
        real(wp), intent(in) :: j5 !! Central body J5 zonal harmonic coefficient
        real(wp), dimension(6), intent(in) :: cartesian !! Cartesian state vector [x, y, z, vx, vy, vz] (km, km/s)
        integer, intent(out) :: stat !! Status: 0 success; /=0 failure.
        real(wp), dimension(6), intent(out) :: blml !! Brouwer mean elements [sma(km), ecc, inc(deg), raan(deg), aop(deg), ma(deg)]

        real(wp), parameter :: tol = 1.0e-8_wp
        integer, parameter :: maxiter = 75

        real(wp), dimension(6) :: cart, kep, kep2, blmean, blmean2
        real(wp), dimension(6) :: aeq, aeq2, aeqmean, aeqmean2, tmp, cart2
        real(wp) :: radper, emag, emag_old, sum_sq_diff, sum_sq_cart, inc_arg
        integer :: pseudostate, ii

        stat = BROUWER_SUCCESS
        blml = 0.0_wp

        if (mu <= 0.0_wp .or. req <= 0.0_wp) then
            stat = BROUWER_INVALID_MU
            return
        end if

        cart = cartesian
        call cartesian_to_keplerian(mu, cart, anomaly_type="TA", stat=stat, kepl=kep)
        if (stat /= BROUWER_SUCCESS) return

        if (kep(2) >= 0.99_wp .or. kep(2) < 0.0_wp) then
            stat = BROUWER_INVALID_ECCENTRICITY
            return
        end if

        radper = kep(1) * (1.0_wp - kep(2))
        if (radper < min_brouwer_radper) then
            stat = BROUWER_PERIAPSIS_TOO_LOW
            return
        end if

        if (kep(3) > 180.0_wp) then
            stat = BROUWER_INVALID_INCLINATION
            return
        end if

        ! Convert TA to MA
        kep(6) = true_to_mean_anomaly(kep(6) * deg2rad, kep(2)) * rad2deg

        pseudostate = 0
        if (kep(3) > 175.0_wp) then
            kep(3) = 180.0_wp - kep(3)
            kep(4) = -kep(4)
            call keplerian_to_cartesian(mu, kep, anomaly_type="MA", stat=stat, cart=cart)
            if (stat /= BROUWER_SUCCESS) return
            pseudostate = 1
        end if

        blmean = kep
        call brouwer_mean_long_to_osculating(mu, req, j2, j3, j4, j5, kep, stat=stat, kepl=kep2)
        if (stat /= BROUWER_SUCCESS) return

        ! Alternate equinoctial elements
        aeq(1) = kep(1)
        aeq(2) = kep(2) * sin((kep(5) + kep(4)) * deg2rad)
        aeq(3) = kep(2) * cos((kep(5) + kep(4)) * deg2rad)
        aeq(4) = sin(kep(3) * 0.5_wp * deg2rad) * sin(kep(4) * deg2rad)
        aeq(5) = sin(kep(3) * 0.5_wp * deg2rad) * cos(kep(4) * deg2rad)
        aeq(6) = kep(4) + kep(5) + kep(6)

        aeq2(1) = kep2(1)
        aeq2(2) = kep2(2) * sin((kep2(5) + kep2(4)) * deg2rad)
        aeq2(3) = kep2(2) * cos((kep2(5) + kep2(4)) * deg2rad)
        aeq2(4) = sin(kep2(3) * 0.5_wp * deg2rad) * sin(kep2(4) * deg2rad)
        aeq2(5) = sin(kep2(3) * 0.5_wp * deg2rad) * cos(kep2(4) * deg2rad)
        aeq2(6) = kep2(4) + kep2(5) + kep2(6)

        aeqmean(1) = blmean(1)
        aeqmean(2) = blmean(2) * sin((blmean(5) + blmean(4)) * deg2rad)
        aeqmean(3) = blmean(2) * cos((blmean(5) + blmean(4)) * deg2rad)
        aeqmean(4) = sin(blmean(3) * 0.5_wp * deg2rad) * sin(blmean(4) * deg2rad)
        aeqmean(5) = sin(blmean(3) * 0.5_wp * deg2rad) * cos(blmean(4) * deg2rad)
        aeqmean(6) = blmean(4) + blmean(5) + blmean(6)

        aeqmean2 = aeqmean + (aeq - aeq2)

        emag = 0.9_wp
        emag_old = 1.0_wp
        ii = 0

        do while (emag > tol)
            blmean2(1) = aeqmean2(1)
            blmean2(2) = sqrt(aeqmean2(2)**2 + aeqmean2(3)**2)

            inc_arg = aeqmean2(4)**2 + aeqmean2(5)**2
            blmean2(3) = acos(1.0_wp - 2.0_wp * min(1.0_wp, inc_arg)) * rad2deg

            blmean2(4) = atan2(aeqmean2(4), aeqmean2(5)) * rad2deg
            if (blmean2(4) < 0.0_wp) blmean2(4) = blmean2(4) + 360.0_wp

            blmean2(5) = atan2(aeqmean2(2), aeqmean2(3)) * rad2deg - blmean2(4)
            if (blmean2(5) < 0.0_wp) blmean2(5) = blmean2(5) + 360.0_wp

            blmean2(6) = aeqmean2(6) - atan2(aeqmean2(2), aeqmean2(3)) * rad2deg

            call brouwer_mean_long_to_osculating(mu, req, j2, j3, j4, j5, blmean2, stat=stat, kepl=kep2)
            if (stat /= BROUWER_SUCCESS) return
            call keplerian_to_cartesian(mu, kep2, anomaly_type="MA", stat=stat, cart=cart2)
            if (stat /= BROUWER_SUCCESS) return

            tmp = cart - cart2
            sum_sq_diff = sum(tmp**2)
            sum_sq_cart = sum(cart**2)
            emag = sqrt(sum_sq_diff) / sqrt(sum_sq_cart)

            if (emag_old > emag) then
                emag_old = emag

                aeq2(1) = kep2(1)
                aeq2(2) = kep2(2) * sin((kep2(5) + kep2(4)) * deg2rad)
                aeq2(3) = kep2(2) * cos((kep2(5) + kep2(4)) * deg2rad)
                aeq2(4) = sin(kep2(3) * 0.5_wp * deg2rad) * sin(kep2(4) * deg2rad)
                aeq2(5) = sin(kep2(3) * 0.5_wp * deg2rad) * cos(kep2(4) * deg2rad)
                aeq2(6) = kep2(4) + kep2(5) + kep2(6)

                aeqmean = aeqmean2
                aeqmean2 = aeqmean + (aeq - aeq2)
            else
                stat = BROUWER_ITERATION_DID_NOT_CONVERGE
                exit
            end if

            ii = ii + 1
            if (ii > maxiter) then
                stat = BROUWER_ITERATION_DID_NOT_CONVERGE
                exit
            end if
        end do

        blmean(1) = aeqmean(1)
        blmean(2) = sqrt(aeqmean(2)**2 + aeqmean(3)**2)

        inc_arg = aeqmean(4)**2 + aeqmean(5)**2
        blmean(3) = acos(1.0_wp - 2.0_wp * min(1.0_wp, inc_arg)) * rad2deg

        blmean(4) = atan2(aeqmean(4), aeqmean(5)) * rad2deg
        if (blmean(4) < 0.0_wp) blmean(4) = blmean(4) + 360.0_wp

        blmean(5) = atan2(aeqmean(2), aeqmean(3)) * rad2deg - blmean(4)
        blmean(6) = aeqmean(6) - atan2(aeqmean(2), aeqmean(3)) * rad2deg

        if (pseudostate /= 0) then
            blmean(3) = 180.0_wp - blmean(3)
            blmean(4) = -blmean(4)
        end if

        call wrap_0_360(blmean(4))
        call wrap_0_360(blmean(5))
        call wrap_0_360(blmean(6))

        blml = blmean

    end subroutine cartesian_to_brouwer_mean_long

    !--------------------------------------------------------------------------
    !>
    !  Converts Brouwer-Lyddane Mean Elements (short and long period terms) to
    !  Osculating Keplerian Elements.

    pure subroutine brouwer_mean_long_to_osculating(mu, req, j2, j3, j4, j5, blml, stat, kepl)
        real(wp), intent(in) :: mu !! Central body gravitational parameter (km^3/s^2)
        real(wp), intent(in) :: req !! Central body equatorial radius (km)
        real(wp), intent(in) :: j2 !! Central body J2 zonal harmonic coefficient
        real(wp), intent(in) :: j3 !! Central body J3 zonal harmonic coefficient
        real(wp), intent(in) :: j4 !! Central body J4 zonal harmonic coefficient
        real(wp), intent(in) :: j5 !! Central body J5 zonal harmonic coefficient
        real(wp), dimension(6), intent(in) :: blml !! Brouwer mean elements [sma(km), ecc, inc(deg), raan(deg), aop(deg), ma(deg)]
        integer, intent(out) :: stat !! Status: 0 success; /=0 failure.
        real(wp), dimension(6), intent(out) :: kepl !! Osculating Keplerian elements [sma(km), ecc, inc(deg), raan(deg), aop(deg), ma(deg)]

        real(wp) :: smadp, eccdp, incdp, raandp, aopdp, meanAnom, radper, t2, &
                    bk2, bk3, bk4, bk5, eccdp2, cn2, cn, gm2, gmp2, gm4, gmp4, &
                    theta, theta2, theta4, gm3, gmp3, gm5, gmp5, &
                    g3dg2, g4dg2, g5dg2, sinMADP, cosMADP, sinraandp, cosraandp, &
                    tadp, rp, adr, sinta, costa, cs2gta, adr2, adr3, costa2, &
                    a1, a2, a3, a4, a5, a6, sinI, a10, a7, a8p, a8, &
                    b13, b14, b15, a11, a12, a13, a14, a17, a15, a16, a18, a19, a21, a22, &
                    sinI2, cosI2, tanI2, a26, a27, &
                    b1, b2, b3, b4, b5, b6, b7, b8, b9, b10, b11, b12, &
                    sma, sn2gta, snf2gd, csf2gd, sn2gd, cs2gd, sin3gd, cs3gd, sn3fgd, cs3fgd, &
                    sinGD, cosGD, bisubc, blghp, eccdpdl, dltI, sinDH, dlt1e, &
                    blgh, dlte, eccdpdl2, eccdpde2, ecc, sinDH2, squar, sqrI, inc, &
                    ma, raan, aop, arg1, arg2
        integer :: pseudostate

        stat = BROUWER_SUCCESS
        kepl = 0.0_wp

        if (mu <= 0.0_wp .or. req <= 0.0_wp) then
            stat = BROUWER_INVALID_MU
            return
        end if

        pseudostate = 0
        smadp = blml(1) / req
        eccdp = blml(2)
        incdp = blml(3) * deg2rad
        raandp = blml(4) * deg2rad
        aopdp = blml(5) * deg2rad
        meanAnom = blml(6) * deg2rad

        if (incdp > 175.0_wp * deg2rad) then
            incdp = pi - incdp
            raandp = -raandp
            pseudostate = 1
        end if

        if (eccdp > 0.99_wp) then
            stat = BROUWER_INVALID_ECCENTRICITY
            return
        end if

        radper = blml(1) * (1.0_wp - blml(2))
        if (radper < min_brouwer_radper) then
            stat = BROUWER_PERIAPSIS_TOO_LOW
            return
        end if

        if (blml(3) > 180.0_wp) then
            stat = BROUWER_INVALID_INCLINATION
            return
        end if

        call wrap_0_2pi(raandp)
        call wrap_0_2pi(aopdp)
        call wrap_0_2pi(meanAnom)

        bk2 = 0.5_wp * j2
        bk3 = -j3
        bk4 = -(3.0_wp / 8.0_wp) * j4
        bk5 = -j5

        eccdp2 = eccdp * eccdp
        cn2 = 1.0_wp - eccdp2
        cn = sqrt(max(0.0_wp, cn2))
        gm2 = bk2 / (smadp**2)
        gmp2 = gm2 / (cn2**2)
        gm4 = bk4 / (smadp**4)
        gmp4 = gm4 / (cn**8)
        theta = cos(incdp)
        theta2 = theta * theta
        theta4 = theta2 * theta2
        t2 = critical_inclination_t2(theta)

        gm3 = bk3 / (smadp**3)
        gmp3 = gm3 / (cn2**3)
        gm5 = bk5 / (smadp**5)
        gmp5 = gm5 / (cn**10)

        g3dg2 = gmp3 / gmp2
        g4dg2 = gmp4 / gmp2
        g5dg2 = gmp5 / gmp2

        sinMADP = sin(meanAnom)
        cosMADP = cos(meanAnom)
        sinraandp = sin(raandp)
        cosraandp = cos(raandp)

        call mean_to_true_anomaly(meanAnom, eccdp, 1.0e-12_wp, stat=stat, ta=tadp)
        if (stat /= BROUWER_SUCCESS) return

        rp = smadp * (1.0_wp - eccdp2) / (1.0_wp + eccdp * cos(tadp))
        adr = smadp / rp
        sinta = sin(tadp)
        costa = cos(tadp)
        cs2gta = cos(2.0_wp * aopdp + 2.0_wp * tadp)
        adr2 = adr * adr
        adr3 = adr2 * adr
        costa2 = costa * costa

        a1 = ((1.0_wp / 8.0_wp) * gmp2 * cn2) * (1.0_wp - theta2 * (11.0_wp + 40.0_wp * theta2 * t2))
        a2 = ((5.0_wp / 12.0_wp) * g4dg2 * cn2) * (1.0_wp - theta2 * (3.0_wp + 8.0_wp * theta2 * t2))
        a3 = g5dg2 * (3.0_wp * eccdp2 + 4.0_wp)
        a4 = g5dg2 * (1.0_wp - theta2 * (9.0_wp + 24.0_wp * theta2 * t2))
        a5 = (g5dg2 * (3.0_wp * eccdp2 + 4.0_wp)) * (1.0_wp - theta2 * (9.0_wp + 24.0_wp * theta2 * t2))
        a6 = g3dg2 * 0.25_wp
        sinI = sin(incdp)
        a10 = cn2 * sinI
        a7 = a6 * a10
        a8p = g5dg2 * eccdp * (1.0_wp - theta2 * (5.0_wp + 16.0_wp * theta2 * t2))
        a8 = a8p * eccdp

        b13 = eccdp * (a1 - a2)
        b14 = a7 + (5.0_wp / 64.0_wp) * a5 * a10
        b15 = a8 * a10 * (35.0_wp / 384.0_wp)

        a11 = 2.0_wp + eccdp2
        a12 = 3.0_wp * eccdp2 + 2.0_wp
        a13 = theta2 * a12
        a14 = (5.0_wp * eccdp2 + 2.0_wp) * theta4 * t2
        a17 = theta4 * t2**2
        a15 = eccdp2 * theta4 * theta2 * t2**2
        a16 = theta2 * t2
        a18 = eccdp * sinI
        a19 = a18 / (1.0_wp + cn)
        a21 = eccdp * theta
        a22 = eccdp2 * theta
        sinI2 = sin(incdp * 0.5_wp)
        cosI2 = cos(incdp * 0.5_wp)
        tanI2 = tan(incdp * 0.5_wp)
        a26 = 16.0_wp * a16 + 40.0_wp * a17 + 3.0_wp
        a27 = a22 * 0.125_wp * (11.0_wp + 200.0_wp * a17 + 80.0_wp * a16)

        b1 = cn * (a1 - a2) - ((a11 - 400.0_wp * a15 - 40.0_wp * a14 - 11.0_wp * a13) * (1.0_wp / 16.0_wp) &
             + (11.0_wp + 200.0_wp * a17 + 80.0_wp * a16) * a22 * (1.0_wp / 8.0_wp)) * gmp2 &
             + ((-80.0_wp * a15 - 8.0_wp * a14 - 3.0_wp * a13 + a11) * (5.0_wp / 24.0_wp) + (5.0_wp / 12.0_wp) * a26 * a22) * g4dg2
        b2 = a6 * a19 * (2.0_wp + cn - eccdp2) + (5.0_wp / 64.0_wp) * a5 * a19 * cn2 - (15.0_wp / 32.0_wp) * a4 * a18 * cn * cn2 &
             + ((5.0_wp / 64.0_wp) * a5 + a6) * a21 * tanI2 + (9.0_wp * eccdp2 + 26.0_wp) * (5.0_wp / 64.0_wp) * a4 * a18 &
             + (15.0_wp / 32.0_wp) * a3 * a21 * a26 * sinI * (1.0_wp - theta)
        b3 = ((80.0_wp * a17 + 5.0_wp + 32.0_wp * a16) * a22 * sinI * (theta - 1.0_wp) * (35.0_wp / 576.0_wp) * g5dg2 * eccdp) &
             - ((a22 * tanI2 + (2.0_wp * eccdp2 + 3.0_wp * (1.0_wp - cn2 * cn)) * sinI) * (35.0_wp / 1152.0_wp) * a8p)
        b4 = cn * eccdp * (a1 - a2)
        b5 = ((9.0_wp * eccdp2 + 4.0_wp) * a10 * a4 * (5.0_wp / 64.0_wp) + a7) * cn
        b6 = (35.0_wp / 384.0_wp) * a8 * cn2 * cn * sinI
        b7 = cn2 * a18 * t2 * ((1.0_wp / 8.0_wp) * gmp2 * (1.0_wp - 15.0_wp * theta2) &
             + (1.0_wp - 7.0_wp * theta2) * g4dg2 * (-(5.0_wp / 12.0_wp)))
        b8 = (5.0_wp / 64.0_wp) * a3 * cn2 * (1.0_wp - theta2 * (9.0_wp + 24.0_wp * theta2 * t2)) + a6 * cn2
        b9 = a8 * (35.0_wp / 384.0_wp) * cn2
        b10 = sinI * (a22 * a26 * g4dg2 * (5.0_wp / 12.0_wp) - a27 * gmp2)
        b11 = a21 * (a5 * (5.0_wp / 64.0_wp) + a6 + a3 * a26 * (15.0_wp / 32.0_wp) * sinI * sinI)
        b12 = -((80.0_wp * a17 + 32.0_wp * a16 + 5.0_wp) * (a22 * eccdp * sinI * sinI * (35.0_wp / 576.0_wp) * g5dg2) &
              + (a8 * a21 * (35.0_wp / 1152.0_wp)))

        sma = smadp * (1.0_wp + gm2 * ((3.0_wp * theta2 - 1.0_wp) * (eccdp2 / (cn2**3)) * (cn + (1.0_wp / (1.0_wp + cn))) &
              + ((3.0_wp * theta2 - 1.0_wp) / (cn2**3)) * (eccdp * costa) * (3.0_wp + 3.0_wp * eccdp * costa + eccdp2 * costa2) &
              + 3.0_wp * (1.0_wp - theta2) * adr3 * cs2gta))

        sn2gta = sin(2.0_wp * aopdp + 2.0_wp * tadp)
        snf2gd = sin(2.0_wp * aopdp + tadp)
        csf2gd = cos(2.0_wp * aopdp + tadp)
        sn2gd = sin(2.0_wp * aopdp)
        cs2gd = cos(2.0_wp * aopdp)
        sin3gd = sin(3.0_wp * aopdp)
        cs3gd = cos(3.0_wp * aopdp)
        sn3fgd = sin(3.0_wp * tadp + 2.0_wp * aopdp)
        cs3fgd = cos(3.0_wp * tadp + 2.0_wp * aopdp)
        sinGD = sin(aopdp)
        cosGD = cos(aopdp)

        bisubc = t2**2 * ((25.0_wp * theta4 * theta) * (gmp2 * eccdp2))

        if (bisubc >= 0.001_wp) then
            dlt1e = 0.0_wp
            blghp = 0.0_wp
            eccdpdl = 0.0_wp
            dltI = 0.0_wp
            sinDH = 0.0_wp
        else
            blghp = raandp + aopdp + meanAnom + b3 * cs3gd + b1 * sn2gd + b2 * cosGD
            call wrap_0_2pi(blghp)

            dlt1e = b14 * sinGD + b13 * cs2gd - b15 * sin3gd
            eccdpdl = b4 * sn2gd - b5 * cosGD + b6 * cs3gd - 0.25_wp * cn2 * cn * gmp2 &
                      * (2.0_wp * (3.0_wp * theta2 - 1.0_wp) * (adr2 * cn2 + adr + 1.0_wp) * sinta &
                      + 3.0_wp * (1.0_wp - theta2) * ((-adr2 * cn2 - adr + 1.0_wp) * snf2gd &
                      + (adr2 * cn2 + adr + (1.0_wp / 3.0_wp)) * sn3fgd))
            dltI = 0.5_wp * theta * gmp2 * sinI * (eccdp * cs3fgd + 3.0_wp * (eccdp * csf2gd + cs2gta)) &
                   - (a21 / cn2) * (b8 * sinGD + b7 * cs2gd - b9 * sin3gd)
            sinDH = (1.0_wp / cosI2) * (0.5_wp * (b12 * cs3gd + b11 * cosGD + b10 * sn2gd &
                    - 0.5_wp * gmp2 * theta * sinI * (6.0_wp * (eccdp * sinta - meanAnom + tadp) &
                    - (3.0_wp * (sn2gta + eccdp * snf2gd) + eccdp * sn3fgd))))
        end if

        blgh = blghp + ((1.0_wp / (cn + 1.0_wp)) * 0.25_wp * eccdp * gmp2 * cn2 &
               * (3.0_wp * (1.0_wp - theta2) * (sn3fgd * ((1.0_wp / 3.0_wp) + adr2 * cn2 + adr) &
               + snf2gd * (1.0_wp - (adr2 * cn2 + adr))) + 2.0_wp * sinta * (3.0_wp * theta2 - 1.0_wp) &
               * (adr2 * cn2 + adr + 1.0_wp))) + gmp2 * 1.5_wp * ((-2.0_wp * theta - 1.0_wp + 5.0_wp * theta2) &
               * (eccdp * sinta + tadp - meanAnom)) + (3.0_wp + 2.0_wp * theta - 5.0_wp * theta2) &
               * (gmp2 * 0.25_wp * (eccdp * sn3fgd + 3.0_wp * (sn2gta + eccdp * snf2gd)))
        call wrap_0_2pi(blgh)

        dlte = dlt1e + (0.5_wp * cn2 * ((3.0_wp / (cn2**3) * gm2 * (1.0_wp - theta2) * cs2gta &
               * (3.0_wp * eccdp * costa2 + 3.0_wp * costa + eccdp2 * costa * costa2 + eccdp)) &
               - (gmp2 * (1.0_wp - theta2) * (3.0_wp * csf2gd + cs3fgd)) &
               + (3.0_wp * theta2 - 1.0_wp) * gm2 * (1.0_wp / (cn2**3)) * (eccdp * cn + (eccdp / (1.0_wp + cn)) &
               + 3.0_wp * eccdp * costa2 + 3.0_wp * costa + eccdp2 * costa * costa2)))

        eccdpdl2 = eccdpdl * eccdpdl
        eccdpde2 = (eccdp + dlte)**2
        ecc = sqrt(eccdpdl2 + eccdpde2)
        sinDH2 = sinDH * sinDH
        squar = (dltI * cosI2 * 0.5_wp + sinI2)**2
        sqrI = sqrt(sinDH2 + squar)
        inc = 2.0_wp * asin(min(1.0_wp, max(-1.0_wp, sqrI)))
        call wrap_0_2pi(inc)

        if (ecc <= 1.0e-11_wp) then
            aop = 0.0_wp
            if (inc <= 1.0e-7_wp) then
                raan = 0.0_wp
                ma = blgh
            else
                arg1 = sinDH * cosraandp + sinraandp * (0.5_wp * dltI * cosI2 + sinI2)
                arg2 = cosraandp * (0.5_wp * dltI * cosI2 + sinI2) - sinDH * sinraandp
                raan = atan2(arg1, arg2)
                ma = blgh - aop - raan
            end if
        else
            arg1 = eccdpdl * cosMADP + (eccdp + dlte) * sinMADP
            arg2 = (eccdp + dlte) * cosMADP - eccdpdl * sinMADP
            ma = atan2(arg1, arg2)
            call wrap_0_2pi(ma)

            if (inc <= 1.0e-7_wp) then
                raan = 0.0_wp
                aop = blgh - raan - ma
            else
                arg1 = sinDH * cosraandp + sinraandp * (0.5_wp * dltI * cosI2 + sinI2)
                arg2 = cosraandp * (0.5_wp * dltI * cosI2 + sinI2) - sinDH * sinraandp
                raan = atan2(arg1, arg2)
                aop = blgh - ma - raan
            end if
        end if

        call wrap_0_2pi(raan)
        call wrap_0_2pi(aop)

        kepl(1) = sma * req
        kepl(2) = ecc
        kepl(3) = inc * rad2deg
        kepl(4) = raan * rad2deg
        kepl(5) = aop * rad2deg
        kepl(6) = ma * rad2deg

        if (pseudostate /= 0) then
            kepl(3) = 180.0_wp - kepl(3)
            kepl(4) = 360.0_wp - kepl(4)
        end if

    end subroutine brouwer_mean_long_to_osculating

    !--------------------------------------------------------------------------
    !>
    !  Converts Brouwer-Lyddane Mean Elements (short and long period terms) to Cartesian.

    pure subroutine brouwer_mean_long_to_cartesian(mu, req, j2, j3, j4, j5, blml, stat, cart)
        real(wp), intent(in) :: mu !! Central body gravitational parameter (km^3/s^2)
        real(wp), intent(in) :: req !! Central body equatorial radius (km)
        real(wp), intent(in) :: j2 !! Central body J2 zonal harmonic coefficient
        real(wp), intent(in) :: j3 !! Central body J3 zonal harmonic coefficient
        real(wp), intent(in) :: j4 !! Central body J4 zonal harmonic coefficient
        real(wp), intent(in) :: j5 !! Central body J5 zonal harmonic coefficient
        real(wp), dimension(6), intent(in) :: blml !! Brouwer mean elements [sma(km), ecc, inc(deg), raan(deg), aop(deg), ma(deg)]
        integer, intent(out) :: stat !! Status propagated from [[brouwer_mean_long_to_osculating]] or [[keplerian_to_cartesian]] conversion
        real(wp), dimension(6), intent(out) :: cart !! Cartesian state vector [x, y, z, vx, vy, vz] (km, km/s)

        real(wp), dimension(6) :: kepl

        call brouwer_mean_long_to_osculating(mu, req, j2, j3, j4, j5, blml, stat, kepl)
        if (stat /= BROUWER_SUCCESS) then
            cart = 0.0_wp
            return
        end if

        call keplerian_to_cartesian(mu, kepl, anomaly_type="MA", stat=stat, cart=cart)

    end subroutine brouwer_mean_long_to_cartesian

    !--------------------------------------------------------------------------
    !>
    !  Phipps' bounded approximation to \( (1 - 5 cos^2(i))^{-1} \).
    !
    ! This method is used to avoid singularity at the critical inclination (i = 63.4deg).
    !
    ! This method, based on Warren Phipps's 1992 thesis (Eq. 2.47 and 2.48),
    ! approximates the factor that causes the singularity by a function,
    ! named T2 in the thesis.
    !
    !### See also
    !  * The Orekit code `BrouwerLyddanePropagator.java`.

    pure function critical_inclination_t2(cos_inc) result(t2)
        real(wp), intent(in) :: cos_inc !! cosine of the mean inclination
        real(wp) :: t2 !! approximation to \( (1 - 5 cos^2(i))^{-1} \)

        real(wp) :: x, x2, series, product, beta_power, factorial
        integer :: ii

        real(wp), parameter :: beta = 100.0_wp / 2.0_wp**11 !! Phipps T2 regularization parameter

        x = 1.0_wp - 5.0_wp * cos_inc**2
        x2 = x**2
        series = 0.0_wp
        beta_power = 1.0_wp
        factorial = 1.0_wp
        do ii = 0, 12
            factorial = factorial * real(ii + 1, wp)
            if (mod(ii, 2) == 0) then
                series = series + beta_power * x2**ii / factorial
            else
                series = series - beta_power * x2**ii / factorial
            end if
            beta_power = beta_power * beta
        end do

        product = 1.0_wp
        do ii = 0, 10
            product = product * (1.0_wp + exp(-beta * x2 * 2.0_wp**ii))
        end do
        t2 = beta * x * series * product

    end function critical_inclination_t2

    !--------------------------------------------------------------------------
    !>
    !  Converts Cartesian state to Keplerian elements.

    pure subroutine cartesian_to_keplerian(mu, cartesian, anomaly_type, stat, kepl)
        real(wp), intent(in) :: mu !! Central body gravitational parameter (km^3/s^2)
        real(wp), dimension(6), intent(in) :: cartesian !! Cartesian state vector [x, y, z, vx, vy, vz] (km, km/s)
        character(len=*), intent(in), optional :: anomaly_type !! Anomaly type: "TA" (True Anomaly, default) or "MA" (Mean Anomaly)
        integer, intent(out) :: stat !! Status: 0 success; /=0 failure.
        real(wp), dimension(6), intent(out) :: kepl !! Keplerian state [sma(km), ecc, inc(deg), raan(deg), aop(deg), anom(deg)]

        real(wp), dimension(3) :: pos, vel, angMomentum, nodeVec, eccVec
        real(wp) :: h, n, posMag, velMag, e, zeta, sma, inc, raan, argPeriapsis, trueAnom, anom
        character(len=2) :: anom_type

        real(wp), parameter :: singular_tol = 0.001_wp !! minimum allowed rp (km)

        stat = BROUWER_SUCCESS
        kepl = 0.0_wp

        if (present(anomaly_type)) then
            anom_type = anomaly_type
        else
            anom_type = "TA"
        end if

        if (abs(mu) < 1.0e-30_wp) then
            stat = BROUWER_INVALID_MU
            return
        end if

        pos = cartesian(1:3)
        vel = cartesian(4:6)

        posMag = norm2(pos)
        velMag = norm2(vel)

        if (posMag == 0.0_wp .or. velMag == 0.0_wp) then
            stat = BROUWER_DEGENERATE_STATE
            return
        end if

        ! Angular momentum vector h = pos x vel
        angMomentum = cross(pos, vel)
        h = norm2(angMomentum)

        if (h == 0.0_wp) then
            stat = BROUWER_DEGENERATE_STATE
            return
        end if

        ! Line of nodes vector n = [0,0,1] x angMomentum
        nodeVec(1) = -angMomentum(2)
        nodeVec(2) = angMomentum(1)
        nodeVec(3) = 0.0_wp
        n = norm2(nodeVec)

        ! Eccentricity vector
        eccVec = (1.0_wp / mu) * ((velMag**2 - mu / posMag) * pos - dot_product(pos, vel) * vel)
        e = norm2(eccVec)

        ! Specific energy zeta
        zeta = 0.5_wp * velMag**2 - mu / posMag
        if (zeta == 0.0_wp .or. abs(1.0_wp - e) <= parabolic_tol) then
            stat = BROUWER_PARABOLIC_ORBIT
            return
        end if

        sma = -mu / (2.0_wp * zeta)

        if (abs(sma * (1.0_wp - e)) < singular_tol) then
            stat = BROUWER_DEGENERATE_STATE
            return
        end if

        inc = acos(min(1.0_wp, max(-1.0_wp, angMomentum(3) / h)))

        raan = 0.0_wp
        argPeriapsis = 0.0_wp
        trueAnom = 0.0_wp

        if (e >= 1.0e-11_wp .and. (inc >= 1.0e-11_wp .and. inc <= (pi - 1.0e-11_wp))) then
            ! Case 1: Non-circular, Inclined Orbit
            if (n == 0.0_wp) then
                stat = BROUWER_DEGENERATE_STATE
                return
            else
                raan = acos(min(1.0_wp, max(-1.0_wp, nodeVec(1) / n)))
                if (nodeVec(2) < 0.0_wp) raan = two_pi - raan

                argPeriapsis = acos(min(1.0_wp, max(-1.0_wp, dot_product(nodeVec, eccVec) / (n * e))))
                if (eccVec(3) < 0.0_wp) argPeriapsis = two_pi - argPeriapsis
            end if

            trueAnom = acos(min(1.0_wp, max(-1.0_wp, dot_product(eccVec, pos) / (e * posMag))))
            if (dot_product(pos, vel) < 0.0_wp) trueAnom = two_pi - trueAnom

        else if (e >= 1.0e-11_wp .and. (inc < 1.0e-11_wp .or. inc > (pi - 1.0e-11_wp))) then
            ! Case 2: Non-circular, Equatorial Orbit
            raan = 0.0_wp
            argPeriapsis = acos(min(1.0_wp, max(-1.0_wp, eccVec(1) / e)))
            if (eccVec(2) < 0.0_wp) argPeriapsis = two_pi - argPeriapsis
            if (inc > (pi - 1.0e-11_wp)) argPeriapsis = -argPeriapsis
            if (argPeriapsis < 0.0_wp) argPeriapsis = argPeriapsis + two_pi

            trueAnom = acos(min(1.0_wp, max(-1.0_wp, dot_product(eccVec, pos) / (e * posMag))))
            if (dot_product(pos, vel) < 0.0_wp) trueAnom = two_pi - trueAnom

        else if (e < 1.0e-11_wp .and. (inc >= 1.0e-11_wp .and. inc <= (pi - 1.0e-11_wp))) then
            ! Case 3: Circular, Inclined Orbit
            if (n == 0.0_wp) then
                stat = BROUWER_DEGENERATE_STATE
                return
            else
                raan = acos(min(1.0_wp, max(-1.0_wp, nodeVec(1) / n)))
                if (nodeVec(2) < 0.0_wp) raan = two_pi - raan
                argPeriapsis = 0.0_wp
                trueAnom = acos(min(1.0_wp, max(-1.0_wp, dot_product(nodeVec, pos) / (n * posMag))))
                if (pos(3) < 0.0_wp) trueAnom = two_pi - trueAnom
            end if

        else
            ! Case 4: Circular, Equatorial Orbit
            raan = 0.0_wp
            argPeriapsis = 0.0_wp
            trueAnom = acos(min(1.0_wp, max(-1.0_wp, pos(1) / posMag)))
            if (pos(2) < 0.0_wp) trueAnom = two_pi - trueAnom
            if (inc > (pi - 1.0e-11_wp)) trueAnom = -trueAnom
            if (trueAnom < 0.0_wp) trueAnom = trueAnom + two_pi
        end if

        anom = trueAnom * rad2deg
        if (anom_type == "MA") then
            anom = true_to_mean_anomaly(trueAnom, e) * rad2deg
        end if

        kepl(1) = sma
        kepl(2) = e
        kepl(3) = inc * rad2deg
        kepl(4) = raan * rad2deg
        kepl(5) = argPeriapsis * rad2deg
        kepl(6) = anom

    end subroutine cartesian_to_keplerian

    !--------------------------------------------------------------------------
    !>
    !  Converts Keplerian elements to Cartesian state vector.

    pure subroutine keplerian_to_cartesian(mu, keplerian, anomaly_type, stat, cart)
        real(wp), intent(in) :: mu !! Central body gravitational parameter (km^3/s^2)
        real(wp), dimension(6), intent(in) :: keplerian !! Keplerian state [sma(km), ecc, inc(deg), raan(deg), aop(deg), anom(deg)]
        character(len=*), intent(in), optional :: anomaly_type !! Anomaly type: "TA" (True Anomaly, default) or "MA" (Mean Anomaly)
        integer, intent(out) :: stat !! Status: 0 success; /=0 failure.
        real(wp), dimension(6), intent(out) :: cart !! Cartesian state [x, y, z, vx, vy, vz] (km, km/s)

        real(wp) :: sma, ecc, inc, raan, per, anom, true_anom, p, onePlusECos, rad, &
                    cosPerAnom, sinPerAnom, cosInc, sinInc, cosRaan, sinRaan, &
                    sqrtGravP, cosAnomPlusE, sinAnom, cosPer, sinPer
        character(len=2) :: anom_type

        stat = BROUWER_SUCCESS
        cart = 0.0_wp

        if (present(anomaly_type)) then
            anom_type = anomaly_type
        else
            anom_type = "TA"
        end if

        sma = keplerian(1)
        ecc = keplerian(2)
        inc = keplerian(3) * deg2rad
        raan = keplerian(4) * deg2rad
        per = keplerian(5) * deg2rad
        anom = keplerian(6) * deg2rad

        if (anom_type == "MA") then
            call mean_to_true_anomaly(anom, ecc, 1.0e-8_wp, stat=stat, ta=true_anom)
            if (stat /= BROUWER_SUCCESS) return
            anom = true_anom
        end if

        p = sma * (1.0_wp - ecc**2)
        if (abs(p) < 1.0e-30_wp) then
            stat = BROUWER_DEGENERATE_STATE
            return
        end if

        onePlusECos = 1.0_wp + ecc * cos(anom)
        if (onePlusECos < 1.0e-10_wp) then
            stat = BROUWER_DEGENERATE_STATE
            return
        end if

        rad = p / onePlusECos

        cosPerAnom = cos(per + anom)
        sinPerAnom = sin(per + anom)
        cosInc = cos(inc)
        sinInc = sin(inc)
        cosRaan = cos(raan)
        sinRaan = sin(raan)
        sqrtGravP = sqrt(mu / p)
        cosAnomPlusE = cos(anom) + ecc
        sinAnom = sin(anom)
        cosPer = cos(per)
        sinPer = sin(per)

        cart(1) = rad * (cosPerAnom * cosRaan - cosInc * sinPerAnom * sinRaan)
        cart(2) = rad * (cosPerAnom * sinRaan + cosInc * sinPerAnom * cosRaan)
        cart(3) = rad * sinPerAnom * sinInc

        cart(4) = sqrtGravP * cosAnomPlusE * (-sinPer * cosRaan - cosInc * sinRaan * cosPer) &
                  - sqrtGravP * sinAnom * (cosPer * cosRaan - cosInc * sinRaan * sinPer)
        cart(5) = sqrtGravP * cosAnomPlusE * (-sinPer * sinRaan + cosInc * cosRaan * cosPer) &
                  - sqrtGravP * sinAnom * (cosPer * sinRaan + cosInc * cosRaan * sinPer)
        cart(6) = sqrtGravP * (cosAnomPlusE * sinInc * cosPer - sinAnom * sinInc * sinPer)

    end subroutine keplerian_to_cartesian

    !--------------------------------------------------------------------------
    !>
    !  Computes Mean Anomaly from True Anomaly.

    pure function true_to_mean_anomaly(ta_radians, ecc) result(ma)
        real(wp), intent(in) :: ta_radians !! True anomaly in radians
        real(wp), intent(in) :: ecc !! Eccentricity
        real(wp) :: ma !! Mean anomaly in radians [0, 2*pi)

        real(wp) :: ea, ha

        if (ecc < (1.0_wp - parabolic_tol)) then
            ea = true_to_eccentric_anomaly(ta_radians, ecc)
            ma = ea - ecc * sin(ea)
        else if (ecc > (1.0_wp + parabolic_tol)) then
            ha = true_to_hyperbolic_anomaly(ta_radians, ecc)
            ma = ecc * sinh(ha) - ha
        else
            ma = 0.0_wp
        end if
        call wrap_0_2pi(ma)
    end function true_to_mean_anomaly

    !--------------------------------------------------------------------------
    !>
    !  Computes Eccentric Anomaly from True Anomaly.

    pure function true_to_eccentric_anomaly(ta_radians, ecc) result(ea)
        real(wp), intent(in) :: ta_radians !! True anomaly in radians
        real(wp), intent(in) :: ecc !! Eccentricity
        real(wp) :: ea !! Eccentric anomaly in radians [0, 2*pi)

        real(wp) :: cosTa, eccCosTa, sinEa, cosEa

        ea = 0.0_wp
        if (ecc <= (1.0_wp - 1.0e-11_wp)) then
            cosTa = cos(ta_radians)
            eccCosTa = ecc * cosTa
            sinEa = (sqrt(max(0.0_wp, 1.0_wp - ecc**2)) * sin(ta_radians)) / (1.0_wp + eccCosTa)
            cosEa = (ecc + cosTa) / (1.0_wp + eccCosTa)
            ea = atan2(sinEa, cosEa)
            call wrap_0_2pi(ea)
        end if
    end function true_to_eccentric_anomaly

    !--------------------------------------------------------------------------
    !>
    !  Computes Hyperbolic Anomaly from True Anomaly.

    pure function true_to_hyperbolic_anomaly(ta_radians, ecc) result(ha)
        real(wp), intent(in) :: ta_radians !! True anomaly in radians
        real(wp), intent(in) :: ecc !! Eccentricity (> 1)
        real(wp) :: ha !! Hyperbolic anomaly in radians

        real(wp) :: tanhHa2

        ha = 0.0_wp
        if (ecc >= (1.0_wp + parabolic_tol)) then
            tanhHa2 = tan(ta_radians * 0.5_wp) * sqrt((ecc - 1.0_wp) / (ecc + 1.0_wp))
            ha = 2.0_wp * atanh(tanhHa2)
        end if
    end function true_to_hyperbolic_anomaly

    !--------------------------------------------------------------------------
    !>
    !  Computes True Anomaly from Mean Anomaly via Newton-Raphson iteration.

    pure subroutine mean_to_true_anomaly(ma_radians, ecc, tol, stat, ta)
        real(wp), intent(in) :: ma_radians !! Mean anomaly in radians
        real(wp), intent(in) :: ecc !! Eccentricity
        real(wp), intent(in), optional :: tol !! Optional convergence tolerance (default = 1.0e-8)
        integer, intent(out) :: stat !! Status: 0 success; /=0 failure.
        real(wp), intent(out) :: ta !! True anomaly in radians [0, 2*pi)

        real(wp), parameter :: ztol = 1.0e-30_wp !! zero tolerance to avoid division by zero
        integer, parameter :: max_iter = 1000 !! maximum number of iterations for Newton-Raphson

        real(wp) :: tol_val, rm, e, e1, e2, temp, temp2, c, f, g, f1, f2
        integer :: iter
        logical :: done

        stat = BROUWER_SUCCESS
        ta = 0.0_wp

        if (present(tol)) then
            tol_val = tol
        else
            tol_val = 1.0e-8_wp
        end if

        rm = ma_radians
        iter = 0

        if (ecc <= 1.0_wp) then
            ! Elliptical orbit
            e2 = rm + ecc * sin(rm)
            done = .false.

            do while (.not. done)
                iter = iter + 1
                if (iter > max_iter) then
                    stat = BROUWER_ITERATION_DID_NOT_CONVERGE
                    return
                end if

                temp = 1.0_wp - ecc * cos(e2)
                if (abs(temp) < ztol) then
                    stat = BROUWER_SINGULARITY
                    return
                end if

                e1 = e2 - (e2 - ecc * sin(e2) - rm) / temp

                if (abs(e2 - e1) < tol_val) then
                    done = .true.
                end if
                e2 = e1
            end do

            e = e2
            if (e < 0.0_wp) e = e + two_pi

            c = abs(e - pi)
            if (c >= 1.0e-8_wp) then
                temp = 1.0_wp - ecc
                if (abs(temp) < ztol) then
                    stat = BROUWER_PARABOLIC_ORBIT
                    return
                else
                    temp2 = (1.0_wp + ecc) / temp
                    if (temp2 < 0.0_wp) then
                        stat = BROUWER_SINGULARITY
                        return
                    else
                        f = sqrt(temp2)
                        g = tan(e * 0.5_wp)
                        ta = 2.0_wp * atan(f * g)
                    end if
                end if
            else
                ta = e
            end if

            call wrap_0_2pi(ta)

        else
            ! Hyperbolic orbit
            f2 = 0.0_wp
            done = .false.

            do while (.not. done)
                iter = iter + 1
                if (iter > max_iter) then
                    stat = BROUWER_ITERATION_DID_NOT_CONVERGE
                    return
                end if

                temp = ecc * cosh(f2) - 1.0_wp
                if (abs(temp) < ztol) then
                    stat = BROUWER_SINGULARITY
                    return
                end if

                f1 = f2 - (ecc * sinh(f2) - f2 - rm) / temp
                if (abs(f2 - f1) < tol_val) then
                    done = .true.
                end if
                f2 = f1
            end do

            f = f2
            temp = ecc - 1.0_wp
            if (abs(temp) < ztol) then
                stat = BROUWER_PARABOLIC_ORBIT
                return
            else
                temp2 = (ecc + 1.0_wp) / temp
                if (temp2 < 0.0_wp) then
                    stat = BROUWER_SINGULARITY
                    return
                else
                    e = sqrt(temp2)
                    g = tanh(f * 0.5_wp)
                    ta = 2.0_wp * atan(e * g)
                end if
            end if

            call wrap_0_2pi(ta)
        end if

    end subroutine mean_to_true_anomaly

!*****************************************************************************************
    pure subroutine wrap_0_2pi(angle)
        !! Wraps an angle in radians to the range [0, 2*pi).
        real(wp), intent(inout) :: angle
        angle = modulo(angle, two_pi)
        if (angle < 0.0_wp) angle = angle + two_pi
    end subroutine wrap_0_2pi
!*****************************************************************************************

!*****************************************************************************************
    pure subroutine wrap_0_360(angle)
        !! Wraps an angle in degrees to the range [0, 360).
        real(wp), intent(inout) :: angle
        angle = modulo(angle, 360.0_wp)
        if (angle < 0.0_wp) angle = angle + 360.0_wp
    end subroutine wrap_0_360
!*****************************************************************************************

!*****************************************************************************************
!> author: Jacob Williams
!
!  Cross product of two 3x1 vectors

    pure function cross(r,v) result(rxv)

    real(wp),dimension(3),intent(in) :: r
    real(wp),dimension(3),intent(in) :: v
    real(wp),dimension(3)            :: rxv

    rxv = [r(2)*v(3) - v(2)*r(3), &
           r(3)*v(1) - v(3)*r(1), &
           r(1)*v(2) - v(1)*r(2) ]

    end function cross
!*****************************************************************************************

end module brouwer_module