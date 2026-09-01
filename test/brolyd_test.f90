program brolyd_test

    use iso_fortran_env, only: wp => real64
    use brolyd_module, only: brolyd
    use brouwer_module, only: keplerian_to_cartesian, &
                              cartesian_to_brouwer_mean_long, &
                              cartesian_to_brouwer_mean_short

    implicit none

    real(wp), parameter :: pi = acos(-1.0_wp)
    real(wp), parameter :: deg2rad = pi / 180.0_wp
    real(wp), parameter :: rad2deg = 1.0_wp / deg2rad
    real(wp), parameter :: mu_earth  = 398600.8_wp
    real(wp), parameter :: req_earth = 6378.135_wp
    real(wp), parameter :: j2_earth  = 0.10826158e-2_wp
    real(wp), parameter :: j3_earth  = -0.25388100e-5_wp
    real(wp), parameter :: j4_earth  = -0.16559700e-5_wp
    real(wp), parameter :: j5_earth  = -0.21848266e-6_wp

    real(wp), dimension(6) :: Oscele !! OUTPUT OSCULATING ELEMENTS AT TIME TTO
    real(wp), dimension(6) :: Dpele  !! INPUT IS OSCULATING ELEMENTS AT EPOCH IF IDMEAN = 0
    real(wp), dimension(6) :: osculating_degrees, cartesian
    real(wp), dimension(6) :: mean_short, mean_long
    integer :: Ipert, Ipass, Idmean
    integer :: stat
    real(wp), dimension(5) :: Orbel

    print*, "BROLYD test program"

    ! set inputs
    Idmean = 0 ! INPUT IS OSCULATING ELEMENTS AT EPOCH IF IDMEAN = 0
    Dpele = [7000.0_wp, 0.001_wp, 33.3_wp*deg2rad, 10.0_wp*deg2rad, 10.0_wp*deg2rad, 10.0_wp*deg2rad] !! example osculating elements

    ! IPERT =0, NO PERTURBATIONS DUE TO OBLATENESS COMPUTED
    ! =1, SECULAR TERMS COMPUTED
    ! =2, SECULAR + LONG PERIODIC + SHORT PERIODIC TERMS
    Ipert = 2 !! full secular, long-periodic, and short-periodic model

    write(*,'(a30,1x,*(f18.10,1x))') "Osculating elements: ", Dpele(1), Dpele(2), Dpele(3)*rad2deg, Dpele(4)*rad2deg, Dpele(5)*rad2deg, Dpele(6)*rad2deg

    osculating_degrees = [Dpele(1), Dpele(2), Dpele(3:6) * rad2deg]
    call keplerian_to_cartesian(mu_earth, osculating_degrees, anomaly_type='MA', stat=stat, cart=cartesian)
    if (stat /= 0) error stop 'Unable to form Cartesian state from osculating elements.'

    ! call main routine
    ! IPASS =1, COMPUTE CONSTANTS NEEDED IN COMPUTATION OF OSCULATING ELEMENTS
    !       =2, UPDATE OSCULATING ELEMENT TO OBSERVATION TIME WITHOUT UPDATING
    Ipass = 1
    call brolyd(Oscele,Dpele,Ipert,Ipass,Idmean,Orbel)
    write(*,'(a30,1x,*(f18.10,1x))') "brolyd mean elements: ", Oscele(1), Oscele(2), Oscele(3)*rad2deg, Oscele(4)*rad2deg, Oscele(5)*rad2deg, Oscele(6)*rad2deg

    ! compare to cartesian_to_brouwer_mean_short and cartesian_to_brouwer_mean_long routines
    call cartesian_to_brouwer_mean_short(mu_earth, req_earth, j2_earth, cartesian, stat, mean_short)
    if (stat /= 0) error stop 'Short-period Brouwer mean conversion failed.'

    call cartesian_to_brouwer_mean_long(mu_earth, req_earth, j2_earth, j3_earth, j4_earth, j5_earth, &
                                        cartesian, stat, mean_long)
    if (stat /= 0) error stop 'Long-period Brouwer mean conversion failed.'

    write(*,'(a30,1x,*(f18.10,1x))') "Brouwer mean short [km, -, deg]: ", mean_short
    write(*,'(a30,1x,*(f18.10,1x))') "Brouwer mean long  [km, -, deg]: ", mean_long
    write(*,'(a30,1x,*(f18.10,1x))') "Long minus short   [km, -, deg]: ", mean_long - mean_short
    write(*,'(a30,1x,*(f18.10,1x))') "brolyd - Brouwer mean short [km, -, deg]: ", [Oscele(1), Oscele(2), Oscele(3:6) * rad2deg] - mean_short
    write(*,'(a30,1x,*(f18.10,1x))') "brolyd - Brouwer mean long  [km, -, deg]: ", [Oscele(1), Oscele(2), Oscele(3:6) * rad2deg] - mean_long



end program brolyd_test