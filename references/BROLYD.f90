
module brolyd_module

  use iso_fortran_env, only: wp => real64

  implicit none

  real(wp),parameter :: pi2 = 2.0_wp * acos(-1.0_wp) !6.283185307179586

  private

  public :: brolyd

contains

SUBROUTINE brolyd(Oscele,Dpele,Ipert,Ipass,Idmean,Orbel)
!*******************************************************************************************
!*    REF. "BROUWER-LYDDANE ORBIT GENERATOR ROUTINE"                *
!*            (X-553-70-223)                                        *
!*        BY E.A. GALBREATH 1970                                    *
!*------------------------------------------------------------------
!*    MODIFIED 7/31/74 VIONA BROWN AND R.A. GORDON TO INTERFACE     *
!*    WITH GTDS                                                     *
!********************************************************************************************
   IMPLICIT NONE

   REAL*8 a , a0 , a1 , a10 , a11 , a12 , a13 , a14 , a15 , a16 , a17 , a18 , a19 , a2 , a21 , a22 , a26 , a27 , a3 , a4
   REAL*8 a5 , a6 , a7 , a8 , a8p , adp , Ae , anu , arg1 , arg2 , b1 , b10 , b11 , b12 , b13 , b14 , b15 , b2 , b3 , b4
   REAL*8 b5 , b6 , b7 , b8 , b9 , bi , bi0 , bidp , bisubc , Bj2 , Bj3 , Bj4 , Bj5 , bk2 , bk3 , bk4 , bk5 , bksubc , bl , bl0
   REAL*8 bldot , bldp , blgh , blghp , bmu , cn , cn2 , cosde , cosfd , cosfd2 , cosgd , coshdp , cosi2 , cosldp , cs2gd ,        &
        & cs2gfd , cs3fgd , cs3gd , csf2gd , dadr
   REAL*8 dadr2 , dadr3 , delt , dlt1e , dlte , dlti , Dpele , e , e0 , ea , eadp , edp , edp2 ,        &
        & edpde2 , edpdl , edpdl2 , ek
   REAL*8 Esq , f15d16 , f15d32 , f1d16 , f1d2 , f1d3 , f1d4 , f1d8 , f35384 , f35576 , f35d52 , f3d2 , f3d32 , f3d8 , f5d12 ,     &
        & f5d16 , f5d24 , f5d4 , f5d64 , fdp
   REAL*8 Fltinv , g , g0 , g3dg2 , g4dg2 , g5dg2 , gdot , gdp , Gm , gm2 , gm3 , gm4 , gm5 , gmp2 , gmp3 , gmp4 , gmp5 , h , h0 , &
        & hdot
   REAL*8 hdp , Orbel , Oscele , R , re , sinde , sindh , sindh2 , sinfd , singd , sinhdp , sini , sini2 , sinldp , sn2gd ,        &
        & sn2gfd , sn3fgd , sn3gd , snf2gd , sqri
   REAL*8 squar , tan12 , tani2 , theta , theta2 , theta4 , Tto , Xke
   INTEGER id8 , Idmean , if , iflg , Ipass , Ipert , nn

   DIMENSION Oscele(6) , Dpele(6) , Orbel(5)
   COMMON /blcnst/ Tto , R , Ae , Gm , Bj2 , Bj3 , Bj4 , Bj5 , Fltinv , Xke , Esq
   INTEGER :: todo
   DATA bmu , re/1.0D0 , 1.0D0/ , bksubc/0.01D0/

   todo = 1
   main: DO
      SELECT CASE (todo)
      CASE (1)
         ek = dsqrt(Gm/Ae**3)
         delt = ek*Tto
         IF ( Ipass==2 ) THEN
            IF ( Ipert==0 ) THEN
               todo = 4
               CYCLE main
            ENDIF
            IF ( Idmean/=0 ) THEN
               todo = 2
               CYCLE main
            ENDIF
            adp = Dpele(1)/Ae
            edp = Dpele(2)
            bidp = Dpele(3)
            hdp = Dpele(4)
            gdp = Dpele(5)
            bldp = Dpele(6)
         ELSE
!
! EPOCH ELEMENTS AT EPOCH TIME
!
            adp = Dpele(1)/Ae
            edp = Dpele(2)
            bidp = Dpele(3)
            hdp = Dpele(4)
            gdp = Dpele(5)
            bldp = Dpele(6)
            a0 = adp
            e0 = edp
            bi0 = bidp
            h0 = hdp
            g0 = gdp
            bl0 = bldp
            iflg = 0
!
! COMPUTE MEAN MOTION
!
            anu = dsqrt(bmu/a0**3)
!
! COMPUTE FRACTIONS
!
            f3d8 = 3.0D0/8.0D0
            f1d2 = 1.0D0/2.0D0
            f3d2 = 3.0D0/2.0D0
            f1d4 = 1.0D0/4.0D0
            f5d4 = 5.0D0/4.0D0
            f1d8 = 1.0D0/8.0D0
            f5d12 = 5.0D0/12.0D0
            f1d16 = 1.0D0/16.0D0
            f15d16 = 15.0D0/16.0D0
            f5d24 = 5.0D0/24.0D0
            f3d32 = 3.0D0/32.0D0
            f15d32 = 15.0D0/32.0D0
            f5d64 = 5.0D0/64.0D0
            f35384 = 35.0D0/384.0D0
            f35576 = 35.0D0/576.0D0
            f35d52 = 35.0D0/1152.0D0
            f1d3 = 1.0D0/3.0D0
            f5d16 = 5.0D0/16.0D0
            bk2 = -f1d2*(Bj2*re*re)
            bk3 = Bj3*re**3
            bk4 = f3d8*(Bj4*re**4)
            bk5 = Bj5*re**5
         ENDIF
         edp2 = edp*edp
         cn2 = 1.0 - edp2
         cn = dsqrt(cn2)
         gm2 = bk2/adp**2
         gmp2 = gm2/(cn2*cn2)
         gm4 = bk4/adp**4
         gmp4 = gm4/cn**8
         theta = dcos(bidp)
         theta2 = theta*theta
         theta4 = theta2*theta2
         todo = 2
      CASE (2)
         IF ( Idmean/=0 ) THEN
            IF ( Ipass==2 ) THEN
               todo = 3
               CYCLE main
            ENDIF
!
! COMPUTE LDOT,GDOT,HDOT
!
            bldot = cn*anu*(gmp2*(f3d2*(3.0*theta2-1)+gmp2*f3d32*(theta2*(-96.0*cn+30.0-90.0*cn2)+(16.0*cn+25.0*cn2-15.0)          &
                  & +theta4*(144.0*cn+25.0*cn2+105.0)))+edp2*gmp4*f15d16*(3.0+35.0*theta4-30.0*theta2))
            gdot = anu*(f5d16*gmp4*((theta2*(126.0*cn2-270.0)+theta4*(385.0-189.0*cn2))-9.0*cn2+21.0)                              &
                 & +gmp2*(f3d32*gmp2*(theta4*(45.0*cn2+360.0*cn+385.0)+theta2*(90.0-192.0*cn-126.0*cn2)+(24.0*cn+25.0*cn2-35))     &
                 & +f3d2*(5*theta2-1)))
            hdot = anu*(gmp4*f5d4*theta*(3.0-7.0*theta2)*(5.0-3.0*cn2)                                                             &
                 & +gmp2*(gmp2*f3d8*(theta*(12.0*cn+9.0*cn2-5.0)-theta*theta2*(5.0*cn2+36.0*cn+35.0))-3*theta))
         ENDIF
         IF ( iflg/=1 ) THEN
!
!  COMPUTE ISUBC TO TEST CRITICAL INCLINATION
!
            bisubc = ((1.0-5.0*theta2)**(-2))*((25.0*theta4*theta)*(gmp2*edp2))
            iflg = 1
!
!  FIRST CHECK FOR CRITICAL INCLINATION
!
            IF ( bisubc>bksubc ) THEN
!
! MODIFICATIONS FOR CRITICAL INCLINATION
!
               dlt1e = 0.0
               blghp = 0.0
               edpdl = 0.0
               dlti = 0.0
               sindh = 0.0
               ASSIGN 40 TO id8
               todo = 3
               CYCLE main
            ELSE
               ASSIGN 20 TO id8
            ENDIF
!
!  IS THERE CRITICAL INCLINATION?
!
         ELSEIF ( bisubc>bksubc ) THEN
            todo = 3
            CYCLE main
         ENDIF
         IF ( Ipert/=1 ) THEN
            gm3 = bk3/adp**3
            gmp3 = gm3/(cn2*cn2*cn2)
            gm5 = bk5/adp**5
            gmp5 = gm5/cn**10
            g3dg2 = gmp3/gmp2
            g4dg2 = gmp4/gmp2
            g5dg2 = gmp5/gmp2
!
! COMPUTE A1-A8
!
            a1 = (f1d8*gmp2*cn2)*(1.0-11.0*theta2-((40.0*theta4)/(1.0-5.0*theta2)))
            a2 = (f5d12*g4dg2*cn2)*(1.0-((8.0*theta4)/(1.0-5.0*theta2))-3.0*theta2)
            a3 = g5dg2*((3.0*edp2)+4.0)
            a4 = g5dg2*(1.0-(24.0*theta4)/(1.0-5.0*theta2)-9.0*theta2)
            a5 = (g5dg2*(3.0*edp2+4.0))*(1.0-(24.0*theta4)/(1.0-5.0*theta2)-9.0*theta2)
            a6 = g3dg2*f1d4
            sini = dsin(bidp)
            a10 = cn2*sini
            a7 = a6*a10
            a8p = g5dg2*edp*(1.0-(16.0*theta4)/(1.0-5.0*theta2)-5.0*theta2)
            a8 = a8p*edp
!
!   COMPUTE B13-B15
!
            b13 = edp*(a1-a2)
            b14 = a7 + f5d64*a5*a10
            b15 = a8*a10*f35384
!
!   COMPUTE A11-A27
!
            a11 = 2.0 + edp2
            a12 = 3.0*edp2 + 2.0
            a13 = theta2*a12
            a14 = (5.0*edp2+2.0)*(theta4/(1.0-5.0*theta2))
            a17 = theta4/((1.0-5.0*theta2)*(1.0-5.0*theta2))
            a15 = (edp2*theta4*theta2)/((1.0-5.0*theta2)*(1.0-5.0*theta2))
            a16 = theta2/(1.0-5.0*theta2)
            a18 = edp*sini
            a19 = a18/(1.0+cn)
            a21 = edp*theta
            a22 = edp2*theta
            sini2 = dsin(bidp/2.0)
            cosi2 = dcos(bidp/2.0)
            tani2 = dtan(bidp/2.0)
            a26 = 16.0*a16 + 40.0*a17 + 3.0
            a27 = a22*f1d8*(11.0+200.0*a17+80.0*a16)
!
!  COMPUTE B1-B12
!
            b1 = cn*(a1-a2) - ((a11-400.0*a15-40.0*a14-11.0*a13)*f1d16+(11.0+200.*a17+80.0*a16)*a22*f1d8)                          &
               & *gmp2 + ((-80.0*a15-8.0*a14-3.0*a13+a11)*f5d24+f5d12*a26*a22)*g4dg2
            b2 = a6*a19*(2.0+cn-edp2) + f5d64*a5*a19*cn2 - f15d32*a4*a18*cn*cn2 + (f5d64*a5+a6)*a21*tan12 + (9.0*edp2+26.0)        &
               & *f5d64*a4*a18 + f15d32*a3*a21*a26*sini*(1.0-theta)
            b3 = ((80.0*a17+5.0+32.0*a16)*a22*sini*(theta-1.0)*f35576*g5dg2*edp) - ((a22*tani2+(2.0*edp2+3.0*(1.0-cn2*cn))*sini)   &
               & *f35d52*a8p)
            b4 = cn*edp*(a1-a2)
            b5 = ((9.0*edp2+4.0)*a10*a4*f5d64+a7)*cn
            b6 = f35384*a8*cn2*cn*sini
            b7 = ((cn2*a18)/(1.0-5.0*theta2))*(f1d8*gmp2*(1.0-15.0*theta2)+(1.0-7.0*theta2)*g4dg2*(-f5d12))
            b8 = f5d64*(a3*cn2*(1.0-9.0*theta2-(24.0*theta4/(1.0-5.0*theta2)))) + a6*cn2
            b9 = a8*f35384*cn2
            b10 = sini*(a22*a26*g4dg2*f5d12-a27*gmp2)
            b11 = a21*(a5*f5d64+a6+a3*a26*f15d32*sini*sini)
            b12 = -((80.0*a17+32.0*a16+5.0)*(a22*edp*sini*sini*f35576*g5dg2)+(a*a21*f35d52))
         ENDIF
         todo = 3
      CASE (3)
         IF ( Ipert==0 ) THEN
            todo = 4
            CYCLE main
         ENDIF
         IF ( Idmean/=0 ) THEN
!
!  COMPUTE SECULAR TERMS
!  "MEAN" MEAN ANOMALY
!
            bldp = anu*delt + bldot*delt + bl0
            bldp = dmod(bldp,pi2)
            IF ( bldp<0.0D0 ) bldp = bldp + pi2
!
! MEAN ARGUMENT OF PERIGEE
!
            gdp = gdot*delt + g0
            gdp = dmod(gdp,pi2)
            IF ( gdp<0.0D0 ) gdp = gdp + pi2
!
!  MEAN LONGITUDE OF ASCENDING NODE
!
            hdp = hdot*delt + h0
            hdp = dmod(hdp,pi2)
            IF ( hdp<0.0D0 ) hdp = hdp + pi2
         ENDIF
         DO nn = 1 , 6
            Oscele(nn) = Dpele(nn)
         ENDDO
         a = adp
         e = edp
         bi = bidp
         h = hdp
         g = gdp
         bl = bldp
!
! COMPUTE TRUE ANOMALY (DOUBLE PRIMED)
!
         eadp = dkeplr(bldp,edp)
         sinde = dsin(eadp)
         cosde = dcos(eadp)
         sinfd = cn*sinde
         cosfd = cosde - edp
         fdp = datan0(sinfd,cosfd)
         IF ( Ipert==1 ) THEN
            todo = 4
            CYCLE main
         ENDIF
         dadr = (1.0-edp*cosde)**(-1)
         sinfd = sinfd*dadr
         cosfd = cosfd*dadr
         cs2gfd = dcos(2.0*gdp+2.0*fdp)
         dadr2 = dadr*dadr
         dadr3 = dadr2*dadr
         cosfd2 = cosfd*cosfd
!
! COMPUTE A (SEMI-MAJOR AXIS)
!
         a = adp*(1.0+gm2*((3.0*theta2-1.0)*(edp2/(cn2*cn2*cn2))*(cn+(1.0/(1+cn)))+((3.0*theta2-1.0)/(cn2*cn2*cn2))*(edp*cosfd)    &
           & *(3.0+3.0*edp*cosfd+edp2*cosfd2)+3.0*(1.0-theta2)*dadr3*cs2gfd))
         sn2gfd = dsin(2.0*gdp+2.0*fdp)
         snf2gd = dsin(2.0*gdp+fdp)
         csf2gd = dcos(2.0*gdp+fdp)
         sn2gd = dsin(2.0*gdp)
         cs2gd = dcos(2.0*gdp)
         sn3gd = dsin(3.0*gdp)
         cs3gd = dcos(3.0*gdp)
         sn3fgd = dsin(3.0*fdp+2.0*gdp)
         cs3fgd = dcos(3.0*fdp+2.0*gdp)
         singd = dsin(gdp)
         cosgd = dcos(gdp)
         GOTO id8
 20      dlt1e = b14*singd + b13*cs2gd - b15*sn3gd
!
! COMPUTE (L+G+H) PRIMED
!
         blghp = hdp + gdp + bldp + b3*cs3gd + b1*sn2gd + b2*cosgd
         blghp = dmod(blghp,pi2)
         IF ( blghp<0.0D0 ) blghp = blghp + pi2
         edpdl = b4*sn2gd - b5*cosgd + b6*cs3gd - f1d4*cn2*cn*gmp2*(2.0*(3.0*theta-1.0)*(dadr2*cn2+dadr+1.0)*sinfd+3.0*(1.0-theta2)&
               & *((-dadr2*cn2-dadr+1.0)*snf2gd+(dadr2*cn2+dadr+f1d3)*sn3fgd))
         dlti = f1d2*theta*gmp2*sini*(edp*cs3fgd+3.0*(edp*csf2gd+cs2gfd)) - (a21/cn2)*(b8*singd+b7*cs2gd-b9*sn3gd)
         sindh = (1.0/cosi2)                                                                                                       &
               & *(f1d2*(b12*cs3gd+b11*cosgd+b10*sn2gd-(f1d2*gmp2*theta*sini*(6.0*(edp*sinfd-bldp+fdp)-(3.0*(sn2gfd+edp*snf2gd)    &
               & +edp*sn3fgd)))))
!
! COMPUTE (L+G+H)
!
 40      blgh = blghp +                                                                                                            &
              & ((1.0/(cn+1.0))*f1d4*edp*gmp2*cn2*(3.0*(1.0-theta2)*(sn3fgd*(f1d3+dadr2*cn2+dadr)+snf2gd*(1.0-(dadr2*cn2+dadr)))   &
              & +2.0*sinfd*(3.0*theta2-1.0)*(dadr2*cn2+dadr+1.0))) + gmp2*f3d2*((-2.0*theta-1.0+5.0*theta2)*(edp*sinfd+fdp-bldp))  &
              & + (3.0+2.0*theta-5.0*theta2)*(gmp2*f1d4*(edp*sn3fgd+3.0*(sn2gfd+edp*snf2gd)))
         blgh = dmod(blgh,pi2)
         IF ( blgh<0.0D0 ) blgh = blgh + pi2
         dlte = dlt1e +                                                                                                            &
              & (f1d2*cn2*((3.0*(1.0/(cn2*cn2*cn2))*gm2*(1.0-theta2)*cs2gfd*(3.0*edp*cosfd2+3.0*cosfd+edp2*cosfd*cosfd2+edp))      &
              & -(gmp2*(1.0-theta2)*(3.0*csf2gd+cs3fgd))+(3.0*theta2-1.0)*gm2*(1.0/(cn2*cn2*cn2))                                  &
              & *(edp*cn+(edp/(1.0+cn))+3.0*edp*cosfd2+3.0*cosfd+edp2*cosfd*cosfd2)))
         edpdl2 = edpdl*edpdl
         edpde2 = (edp+dlte)*(edp+dlte)
!
! COMPUTE E (ECCENTRICITY)
!
         e = dsqrt(edpdl2+edpde2)
         sindh2 = sindh*sindh
         squar = (dlti*cosi2*f1d2+sini2)*(dlti*cosi2*f1d2+sini2)
         sqri = dsqrt(sindh2+squar)
!
! COMPUTE BI (INCLINATION)
!
         bi = asin(sqri)  ! this was darsin in the original code
         bi = 2.0*bi
         bi = dmod(bi,pi2)
         IF ( bi<0.0D0 ) bi = bi + pi2
!
! CHECK FOR E (ECCENTRICITY)=0
!
         IF ( e/=0.0 ) THEN
            sinldp = dsin(bldp)
            cosldp = dcos(bldp)
            sinhdp = dsin(hdp)
            coshdp = dcos(hdp)
!
! COMPUTE L (MEAN ANOMALY)
!
            arg1 = edpdl*cosldp + (edp+dlte)*sinldp
            arg2 = (edp+dlte)*cosldp - (edpdl*sinldp)
            bl = datan2(arg1,arg2)
            bl = dmod(bl,pi2)
            IF ( bl<0.0D0 ) bl = bl + pi2
         ELSE
            bl = 0.0
         ENDIF
!
! CHECK FOR BI (INCLINATION)=0
!
         IF ( bi/=0.0 ) THEN
!
! COMPUTE H (LONGITUDE OF ASCENDING NODE)
!
            arg1 = sindh*coshdp + sinhdp*(f1d2*dlti*cosi2+sini2)
            arg2 = coshdp*(f1d2*dlti*cosi2+sini2) - (sindh*sinhdp)
            h = datan2(arg1,arg2)
            h = dmod(h,pi2)
            IF ( h<0.0D0 ) h = h + pi2
         ELSE
            h = 0.0
         ENDIF
!
! COMPUTE G (ARGUMENT OF PERIGEE)
!
         g = blgh - bl - h
         g = dmod(g,pi2)
         IF ( g<0.0D0 ) g = g + pi2
!
! COMPUTE TRUE ANOMALY
!
         ea = dkeplr(bl,e)
         arg1 = dsin(ea)*dsqrt(1.0-e**2)
         arg2 = dcos(ea) - e
         if = datan0(arg1,arg2)
         Oscele(1) = a*Ae
         Oscele(2) = e
         Oscele(3) = bi
         Oscele(4) = h
         Oscele(5) = g
         Oscele(6) = bl
         todo = 4
      CASE (4)
         Dpele(1) = adp*Ae
         Dpele(2) = edp
         Dpele(3) = bidp
         Dpele(4) = hdp
         Dpele(5) = gdp
         Dpele(6) = bldp
         IF ( Ipert==0 ) bl = dmod(anu*delt,pi2)
         Orbel(1) = eadp
         Orbel(2) = gdp + fdp
         Orbel(3) = gdp
         Orbel(4) = ek*(anu+bldot)
         Orbel(5) = fdp
         R = a*Ae*(1.0D0-e*dcos(ea))
         EXIT main
      END SELECT
   ENDDO main
END SUBROUTINE brolyd

FUNCTION dkeplr(M,E)
   !! SUBROUTINE TO SOLVE KEPLER'S EQUATION
   !! KEPLER'S EQUATION RELATES GEOMETRY OR POSITION IN ORBIT PLANE TO TIME.
   !! M - MEAN ANOMALY (0<M<2PI)
   !! E - ECCENTRICITY
   !! EA - ECCENTRIC ANOMALY
   IMPLICIT NONE

   REAL*8 delea , dkeplr , E , ea , fe , oldea
   INTEGER i

   REAL*8 M , tol/0.5D-15/

   ea = 0
   IF ( M/=0 ) THEN
      ea = M + E*dsin(M)
      DO i = 1 , 12
         oldea = ea
         fe = ea - E*dsin(ea) - M
         ea = ea - fe/(1-E*dcos(ea-0.5D0*fe))
         ! TEST FOR CONVERGENCE
         delea = dabs(ea-oldea)
         IF ( delea<=tol ) EXIT
      END DO
   END IF
   ea = dmod(ea,pi2)
   dkeplr = ea
END FUNCTION dkeplr

REAL*8 FUNCTION DATAN0(ARG1,ARG2)
!  VERSION OF 03/10/71
!
!  FORTRAN IV FUNCTION SUBROUTINE FOR THE IBM-360
!  PURPOSE
!  COMPUTE A VALUE FOR THE ARCTAN BETWEEN 0 AND 2 PI WHERE THE
!  TANGENT IS DEFINED BY THE TWO INPUT ARGUMENTS AS ARG1/ARG2
!
!  CALLING SEQUENCE
!  NONE
!  INPUT
!  ARG1 - FIST ARGUMENT OF THE ARC TANGENT
!  ARG2 - SECOND ARGUMENT OF THE ARC TANGENT
!
!  OUTPUT
!  A DOUBLE PRECISION ARC TANGENT (+ VALUE BETWEEN 0 AND 2PI)
!
!  METHOD
!  USES FORTRAN MATH SUBROUTINE DATAN2 WHICH RETURNS A VALUE
!  BETWEEN -PI AND PI, GIVEN TWO ARGUMENTS
!
!  REQUIRED SUBROUTINES
!  1- FUNCTION SUBROUTINE DATAN2
!  PROGRAMMER
!  R. E. GILLIAN - COMPUTING AND SOFTWARE
!
!  COMPUTE ARCTAN BETWEEN -PI AND PI
   REAL*8,intent(in) :: ARG1,ARG2

   DATAN0=DATAN2(ARG1,ARG2)
   IF (DATAN0 >= 0.0d0) return
   ! IF ARCTAN IS NEGATIVE, ADD 2PI TO THE RESULT
   DATAN0 = DATAN0 + pi2

   END FUNCTION DATAN0

end module brolyd_module