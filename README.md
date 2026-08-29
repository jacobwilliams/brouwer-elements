Modern Fortran version of Brouwer-Lyddane mean orbital element conversion routines.

## Overview

The `brouwer_module` provides modern Fortran functions to convert between Cartesian state vectors $[x, y, z, v_x, v_y, v_z]$ and Brouwer-Lyddane mean orbital elements $[a, e, i, \Omega, \omega, M]$ (semi-major axis, eccentricity, inclination, right ascension of ascending node, argument of periapsis, and mean anomaly).

Both short-period perturbation models and full short- plus long-period perturbation models ($J_2, J_3, J_4, J_5$) are supported.

## Features

- **Four Primary Conversion Routines**:
  - `cartesian_to_brouwer_mean_short`: Cartesian $\rightarrow$ Brouwer Mean Elements (short-period terms only)
  - `brouwer_mean_short_to_cartesian`: Brouwer Mean Elements (short-period) $\rightarrow$ Cartesian
  - `cartesian_to_brouwer_mean_long`: Cartesian $\rightarrow$ Brouwer Mean Elements (short and long period terms)
  - `brouwer_mean_long_to_cartesian`: Brouwer Mean Elements (short and long period) $\rightarrow$ Cartesian
- **Supporting Routines**:
  - `brouwer_mean_short_to_osculating` / `brouwer_mean_long_to_osculating`
  - `cartesian_to_keplerian` / `keplerian_to_cartesian` (supports True Anomaly `"TA"` and Mean Anomaly `"MA"`)
  - `true_to_mean_anomaly` / `mean_to_true_anomaly`
  - `true_to_eccentric_anomaly` / `true_to_hyperbolic_anomaly`
- **Modern Fortran**:
  - Standard Fortran 2008+ using `iso_fortran_env: real64`.
  - Pure functions where applicable, with optional `stat` error status arguments.
  - Compatible with Fortran Package Manager (`fpm`).

## Building and Testing

### Using `fpm` directly

```bash
fpm build
fpm test
```

## Usage Example

```fortran
program example
    use iso_fortran_env, only: wp => real64
    use brouwer_module

    implicit none

    ! Central body parameters (e.g. Earth)
    real(wp), parameter :: mu_earth  = 398600.4415_wp           ! km^3/s^2
    real(wp), parameter :: req_earth = 6378.1363_wp             ! km
    real(wp), parameter :: j2_earth  = 1.082626925638815e-3_wp
    real(wp), parameter :: j3_earth  = -0.2532307818191774e-5_wp
    real(wp), parameter :: j4_earth  = -0.1620429990000000e-5_wp
    real(wp), parameter :: j5_earth  = -0.2270711043920343e-6_wp

    real(wp), dimension(6) :: cartesian, blms, blml
    integer :: stat

    ! Cartesian state [x, y, z, vx, vy, vz] in km and km/s
    cartesian = [420.040413_wp, 6512.298337_wp, 2338.986311_wp, &
                 -7.156474_wp,  -0.443318_wp,   2.577366_wp]

    ! Convert to Brouwer Mean Short elements [a (km), e, i (deg), raan (deg), aop (deg), ma (deg)]
    blms = cartesian_to_brouwer_mean_short(mu_earth, req_earth, j2_earth, cartesian, stat=stat)
    if (stat == 0) then
        print *, "Brouwer Mean Short:", blms
    end if

    ! Convert to Brouwer Mean Long elements
    blml = cartesian_to_brouwer_mean_long(mu_earth, req_earth, j2_earth, j3_earth, j4_earth, j5_earth, cartesian, stat=stat)
    if (stat == 0) then
        print *, "Brouwer Mean Long:", blml
    end if

    ! Convert back to Cartesian
    cartesian = brouwer_mean_long_to_cartesian(mu_earth, req_earth, j2_earth, j3_earth, j4_earth, j5_earth, blml, stat=stat)

end program example
```

## Adding to Your Project

Add `brouwer-elements` to your `fpm.toml`:

```toml
[dependencies]
brouwer_elements = { git = "https://github.com/jacobwilliams/brouwer-elements.git" }
```

## References

 * Brouwer, D., "[Solution of the Problem of Artificial Satellite Theory without Drag](https://articles.adsabs.harvard.edu/pdf/1959AJ.....64..378B)," *Astronomical Journal*, Vol. 64, Nov. 1959, pp. 378-397.
 * Lyddane, R. H., "[Small Eccentricities or Inclinations in the Brouwer Theory of the Artificial Satellite](https://articles.adsabs.harvard.edu/pdf/1963AJ.....68..555L)," *Astronomical Journal*, Vol. 68, Oct. 1963, pp. 555-558.
 * [NASA General Mission Analysis Tool](https://github.com/nasa/GMAT) (GMAT), `StateConversionUtil.cpp`.
