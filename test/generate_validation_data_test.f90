!------------------------------------------------------------------------------
!>
!  Generates a CSV file of test case states across a grid of orbit sizes,
!  eccentricities, and inclinations, for validation against other tools.
!
!  For each case, the following are written to `validation_data.csv`:
!
!  * The input Keplerian elements (true anomaly)
!  * The corresponding Cartesian state vector
!  * The osculating Keplerian elements (recovered from the Cartesian state)
!  * The Brouwer mean elements (short-period theory)
!  * The Brouwer mean elements (long-period theory)

program generate_validation_data_test

    use iso_fortran_env, only: error_unit
    use brouwer_module, wp => brouwer_module_wp

    implicit none

    ! Earth parameters
    real(wp), parameter :: mu_earth  = 398600.4415_wp           ! km^3/s^2
    real(wp), parameter :: req_earth = 6378.1363_wp             ! km
    real(wp), parameter :: j2_earth  = 1.082626925638815e-3_wp
    real(wp), parameter :: j3_earth  = -0.2532307818191774e-5_wp
    real(wp), parameter :: j4_earth  = -0.1620429990000000e-5_wp
    real(wp), parameter :: j5_earth  = -0.2270711043920343e-6_wp

    ! Grid of test case parameters
    real(wp), dimension(*), parameter :: sma_grid = [7000.0_wp, 8000.0_wp, &
                                                       12000.0_wp, 26560.0_wp, &
                                                       42164.0_wp]                ! km
    ! Note: 1.0e-6_wp (not exactly 0) avoids a singularity in the
    ! osculating-to-mean iteration, where the argument of periapsis is
    ! recovered via atan2 of two near-zero equinoctial components.
    real(wp), dimension(*), parameter :: ecc_grid = [1.0e-6_wp, 0.01_wp, 0.1_wp, &
                                                       0.3_wp, 0.6_wp]
    ! Note: 1.0e-6_wp (not exactly 0) avoids the equatorial analog of the
    ! same singularity, where RAAN is recovered via atan2 of two
    ! near-zero equinoctial inclination-vector components.
    real(wp), dimension(*), parameter :: inc_grid = [1.0e-6_wp, 28.5_wp, 51.6_wp, &
                                                       63.4349488_wp, 90.0_wp, &
                                                       98.7_wp, 116.5650512_wp, &
                                                       150.0_wp]                  ! deg
    real(wp), dimension(*), parameter :: raan_grid = [0.0_wp, 120.0_wp]          ! deg
    real(wp), dimension(*), parameter :: aop_grid = [0.0_wp, 45.0_wp]            ! deg
    real(wp), dimension(*), parameter :: ta_grid = [0.0_wp, 90.0_wp, 200.0_wp]   ! deg

    character(len=*), parameter :: csv_file = 'validation_data.csv'
    integer :: iunit, i_sma, i_ecc, i_inc, i_raan, i_aop, i_ta
    integer :: case_number, n_written
    real(wp), dimension(6) :: kep_in, cart, kep_osc, blms, blml
    real(wp), dimension(6) :: cart_short_rt, cart_long_rt
    real(wp) :: short_rt_relerr, long_rt_relerr
    integer :: stat_cart, stat_osc, stat_short, stat_long, stat_short_rt, stat_long_rt

    print *, "=========================================================="
    print *, " Generating Brouwer Elements Validation Data"
    print *, "=========================================================="

    open(newunit=iunit, file=csv_file, status='replace', action='write')

    write(iunit, '(A)') &
        'case,'//&
        'in_sma_km,in_ecc,in_inc_deg,in_raan_deg,in_aop_deg,in_ta_deg,'//&
        'x_km,y_km,z_km,vx_kmps,vy_kmps,vz_kmps,'//&
        'osc_sma_km,osc_ecc,osc_inc_deg,osc_raan_deg,osc_aop_deg,osc_ma_deg,'//&
        'short_sma_km,short_ecc,short_inc_deg,short_raan_deg,short_aop_deg,short_ma_deg,'//&
        'long_sma_km,long_ecc,long_inc_deg,long_raan_deg,long_aop_deg,long_ma_deg,'//&
        'short_rt_relerr,long_rt_relerr,'//&
        'stat_cart,stat_osc,stat_short,stat_long,stat_short_rt,stat_long_rt'

    case_number = 0
    n_written = 0

    do i_sma = 1, size(sma_grid)
        do i_ecc = 1, size(ecc_grid)
            do i_inc = 1, size(inc_grid)
                do i_raan = 1, size(raan_grid)
                    do i_aop = 1, size(aop_grid)
                        do i_ta = 1, size(ta_grid)

                            case_number = case_number + 1

                            ! skip physically invalid combinations: periapsis
                            ! radius must clear the atmosphere/surface
                            if (sma_grid(i_sma) * (1.0_wp - ecc_grid(i_ecc)) < req_earth + 150.0_wp) cycle

                            kep_in = [sma_grid(i_sma), ecc_grid(i_ecc), inc_grid(i_inc), &
                                      raan_grid(i_raan), aop_grid(i_aop), ta_grid(i_ta)]

                            call keplerian_to_cartesian(mu_earth, kep_in, anomaly_type="TA", &
                                                         stat=stat_cart, cart=cart)

                            if (stat_cart == BROUWER_SUCCESS) then
                                call cartesian_to_keplerian(mu_earth, cart, anomaly_type="MA", &
                                                             stat=stat_osc, kepl=kep_osc)
                            else
                                stat_osc = -1
                                kep_osc = 0.0_wp
                            end if

                            if (stat_cart == BROUWER_SUCCESS) then
                                call cartesian_to_brouwer_mean_short(mu_earth, req_earth, j2_earth, &
                                                                      cart, stat=stat_short, blms=blms)
                            else
                                stat_short = -1
                                blms = 0.0_wp
                            end if

                            if (stat_cart == BROUWER_SUCCESS) then
                                call cartesian_to_brouwer_mean_long(mu_earth, req_earth, j2_earth, &
                                                                     j3_earth, j4_earth, j5_earth, &
                                                                     cart, stat=stat_long, blml=blml)
                            else
                                stat_long = -1
                                blml = 0.0_wp
                            end if

                            ! Round-trip: map the (possibly not-fully-converged) mean
                            ! elements back to Cartesian and compare against the
                            ! original state, to see whether they are still "good
                            ! enough" even when the internal iteration didn't
                            ! strictly converge.
                            if (stat_short == BROUWER_SUCCESS) then
                                call brouwer_mean_short_to_cartesian(mu_earth, req_earth, j2_earth, &
                                                                      blms, stat=stat_short_rt, cart=cart_short_rt)
                                short_rt_relerr = norm2(cart - cart_short_rt) / norm2(cart)
                            else
                                stat_short_rt = -1
                                short_rt_relerr = -1.0_wp
                            end if

                            if (stat_long == BROUWER_SUCCESS) then
                                call brouwer_mean_long_to_cartesian(mu_earth, req_earth, j2_earth, &
                                                                     j3_earth, j4_earth, j5_earth, &
                                                                     blml, stat=stat_long_rt, cart=cart_long_rt)
                                long_rt_relerr = norm2(cart - cart_long_rt) / norm2(cart)
                            else
                                stat_long_rt = -1
                                long_rt_relerr = -1.0_wp
                            end if

                            write(iunit, '(I0,32(",",ES24.16E3),6(",",I0))') &
                                case_number, kep_in, cart, kep_osc, blms, blml, &
                                short_rt_relerr, long_rt_relerr, &
                                stat_cart, stat_osc, stat_short, stat_long, stat_short_rt, stat_long_rt

                            n_written = n_written + 1

                        end do
                    end do
                end do
            end do
        end do
    end do

    close(iunit)

    print '(A,I0,A,I0,A)', ' Wrote ', n_written, ' of ', case_number, ' generated cases.'
    print *, ' Output file: '//csv_file
    print *, "=========================================================="

end program generate_validation_data_test
