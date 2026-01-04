       module mod_specialfun
       use mod_kinds, only:dp 

       implicit none
       private
       public:: dlgamma,dbeta, beta,psi
       
     contains

!CS    REAL FUNCTION ALGAMA(X)
!C    DOUBLE PRECISION FUNCTION DLGAMA(X)
!C      real(dp)  FUNCTION QLGAMMA(X)
        real(dp)  FUNCTION DLGAMMA(X)
!C----------------------------------------------------------------------
!C
!C This routine calculates the LOG(GAMMA) function for a positive real
!C   argument X.  Computation is based on an algorithm outlined in
!C   references 1 and 2.  The program uses rational functions that
!C   theoretically approximate LOG(GAMMA) to at least 18 significant
!C   decimal digits.  The approximation for X > 12 is from reference
!C   3, while approximations for X < 12.0 are similar to those in
!C   reference 1, but are unpublished.  The accuracy achieved depends
!C   on the arithmetic system, the compiler, the intrinsic functions,
!C   and proper selection of the machine-dependent constants.
!C
!C
!C*********************************************************************
!C*********************************************************************
!C
!C Explanation of machine-dependent constants
!C
!C beta   - radix for the floating-point representation
!C maxexp - the smallest positive power of beta that overflows
!C XBIG   - largest argument for which LN(GAMMA(X)) is representable
!C          in the machine, i.e., the solution to the equation
!C                  LN(GAMMA(XBIG)) = beta**maxexp
!C XINF   - largest machine representable floating-point number;
!C          approximately beta**maxexp.
!C EPS    - The smallest positive floating-point number such that
!C          1.0+EPS .GT. 1.0
!C FRTBIG - Rough estimate of the fourth root of XBIG
!C
!C
!C     Approximate values for some important machines are:
!C
!C                           beta      maxexp         XBIG
!C
!C CRAY-1        (S.P.)        2        8191       9.62E+2461
!C Cyber 180/855
!C   under NOS   (S.P.)        2        1070       1.72E+319
!C IEEE (IBM/XT,
!C   SUN, etc.)  (S.P.)        2         128       4.08E+36
!C IEEE (IBM/XT,
!C   SUN, etc.)  (D.P.)        2        1024       2.55D+305
!C IBM 3033      (D.P.)       16          63       4.29D+73
!C VAX D-Format  (D.P.)        2         127       2.05D+36
!C VAX G-Format  (D.P.)        2        1023       1.28D+305
!C
!C
!C                           XINF        EPS        FRTBIG
!C
!C CRAY-1        (S.P.)   5.45E+2465   7.11E-15    3.13E+615
!C Cyber 180/855
!C   under NOS   (S.P.)   1.26E+322    3.55E-15    6.44E+79
!C IEEE (IBM/XT,
!C   SUN, etc.)  (S.P.)   3.40E+38     1.19E-7     1.42E+9
!C IEEE (IBM/XT,
!C   SUN, etc.)  (D.P.)   1.79D+308    2.22D-16    2.25D+76
!C IBM 3033      (D.P.)   7.23D+75     2.22D-16    2.56D+18
!C VAX D-Format  (D.P.)   1.70D+38     1.39D-17    1.20D+9
!C VAX G-Format  (D.P.)   8.98D+307    1.11D-16    1.89D+76
!C
!C**************************************************************
!C**************************************************************
!C
!C Error returns
!C
!C  The program returns the value XINF for X .LE. 0.0 or when
!C     overflow would occur.  The computation is believed to 
!C     be free of underflow and overflow.
!C
!C
!C Intrinsic functions required are:
!C
!C      LOG
!C
!C
!C References:
!C
!C  1) W. J. Cody and K. E. Hillstrom, 'Chebyshev Approximations for
!C     the Natural Logarithm of the Gamma Function,' Math. Comp. 21,
!C     1967, pp. 198-203.
!C
!C  2) K. E. Hillstrom, ANL/AMD Program ANLC366S, DGAMMA/DLGAMA, May,
!C     1969.
!C 
!C  3) Hart, Et. Al., Computer Approximations, Wiley and sons, New
!C     York, 1968.
!C
!C
!C  Authors: W. J. Cody and L. Stoltz
!C           Argonne National Laboratory
!C
!C  Latest modification: June 16, 1988
!C
!C----------------------------------------------------------------------
!c     implicit real*8 (a-h, o-z)
      use mod_kinds, only:dp 
      implicit none
      
      INTEGER I
!CS    REAL      
!CD    DOUBLE PRECISION
       real(dp)::    &
         C,CORR,D1,D2,D4,EPS,FRTBIG,FOUR,HALF,ONE,PNT68,P1,P2,P4,&
         Q1,Q2,Q4,RES,SQRTPI,THRHAL,TWELVE,TWO,X,XBIG,XDEN,XINF, &
         XM1,XM2,XM4,XNUM,Y,YSQ,ZERO
      DIMENSION C(7),P1(8),P2(8),P4(8),Q1(8),Q2(8),Q4(8)
!C----------------------------------------------------------------------
!C  Mathematical constants
!C----------------------------------------------------------------------
!CS    DATA ONE,HALF,TWELVE,ZERO/1.0E0,0.5E0,12.0E0,0.0E0/,
!CS   1     FOUR,THRHAL,TWO,PNT68/4.0E0,1.5E0,2.0E0,0.6796875E0/,
!CS   2     SQRTPI/0.9189385332046727417803297E0/
      DATA ONE,HALF,TWELVE,ZERO/1.0D0,0.5D0,12.0D0,0.0D0/,   &
          FOUR,THRHAL,TWO,PNT68/4.0D0,1.5D0,2.0D0,0.6796875D0/, &
          SQRTPI/0.9189385332046727417803297D0/
!C----------------------------------------------------------------------
!C  Machine dependent parameters
!C----------------------------------------------------------------------
!CS    DATA XBIG,XINF,EPS,FRTBIG/4.08E36,3.401E38,1.19E-7,1.42E9/
      DATA XBIG,XINF,EPS,FRTBIG/2.55D305,1.79D308,2.22D-16,2.25D76/
!C----------------------------------------------------------------------
!C  Numerator and denominator coefficients for rational minimax
!C     approximation over (0.5,1.5).
!C----------------------------------------------------------------------
!CS    DATA D1/-5.772156649015328605195174E-1/
!CS    DATA P1/4.945235359296727046734888E0,2.018112620856775083915565E2,
!CS   1        2.290838373831346393026739E3,1.131967205903380828685045E4,
!CS   2        2.855724635671635335736389E4,3.848496228443793359990269E4,
!CS   3        2.637748787624195437963534E4,7.225813979700288197698961E3/
!CS    DATA Q1/6.748212550303777196073036E1,1.113332393857199323513008E3,
!CS   1        7.738757056935398733233834E3,2.763987074403340708898585E4,
!CS   2        5.499310206226157329794414E4,6.161122180066002127833352E4,
!CS   3        3.635127591501940507276287E4,8.785536302431013170870835E3/
      DATA D1/-5.772156649015328605195174D-1/
      DATA P1/4.945235359296727046734888D0,2.018112620856775083915565D2,&
             2.290838373831346393026739D3,1.131967205903380828685045D4,&
             2.855724635671635335736389D4,3.848496228443793359990269D4,&
             2.637748787624195437963534D4,7.225813979700288197698961D3/
      DATA Q1/6.748212550303777196073036D1,1.113332393857199323513008D3,&
             7.738757056935398733233834D3,2.763987074403340708898585D4,&
             5.499310206226157329794414D4,6.161122180066002127833352D4,&
             3.635127591501940507276287D4,8.785536302431013170870835D3/
!C----------------------------------------------------------------------
!C  Numerator and denominator coefficients for rational minimax
!C     Approximation over (1.5,4.0).
!C----------------------------------------------------------------------
!CS    DATA D2/4.227843350984671393993777E-1/
!CS    DATA P2/4.974607845568932035012064E0,5.424138599891070494101986E2,
!CS   1        1.550693864978364947665077E4,1.847932904445632425417223E5,
!CS   2        1.088204769468828767498470E6,3.338152967987029735917223E6,
!CS   3        5.106661678927352456275255E6,3.074109054850539556250927E6/
!CS    DATA Q2/1.830328399370592604055942E2,7.765049321445005871323047E3,
!CS   1        1.331903827966074194402448E5,1.136705821321969608938755E6,
!CS   2        5.267964117437946917577538E6,1.346701454311101692290052E7,
!CS   3        1.782736530353274213975932E7,9.533095591844353613395747E6/
      DATA D2/4.227843350984671393993777D-1/
      DATA P2/4.974607845568932035012064D0,5.424138599891070494101986D2,&
             1.550693864978364947665077D4,1.847932904445632425417223D5,&
             1.088204769468828767498470D6,3.338152967987029735917223D6,&
             5.106661678927352456275255D6,3.074109054850539556250927D6/
      DATA Q2/1.830328399370592604055942D2,7.765049321445005871323047D3,&
             1.331903827966074194402448D5,1.136705821321969608938755D6,&
             5.267964117437946917577538D6,1.346701454311101692290052D7,&
             1.782736530353274213975932D7,9.533095591844353613395747D6/
!C----------------------------------------------------------------------
!C  Numerator and denominator coefficients for rational minimax
!C     Approximation over (4.0,12.0).
!C----------------------------------------------------------------------
!CS    DATA D4/1.791759469228055000094023E0/
!CS    DATA P4/1.474502166059939948905062E4,2.426813369486704502836312E6,
!CS   1        1.214755574045093227939592E8,2.663432449630976949898078E9,
!CS   2      2.940378956634553899906876E10,1.702665737765398868392998E11,
!CS   3      4.926125793377430887588120E11,5.606251856223951465078242E11/
!CS    DATA Q4/2.690530175870899333379843E3,6.393885654300092398984238E5,
!CS   2        4.135599930241388052042842E7,1.120872109616147941376570E9,
!CS   3      1.488613728678813811542398E10,1.016803586272438228077304E11,
!CS   4      3.417476345507377132798597E11,4.463158187419713286462081E11/
      DATA D4/1.791759469228055000094023D0/
      DATA P4/1.474502166059939948905062D4,2.426813369486704502836312D6,&
             1.214755574045093227939592D8,2.663432449630976949898078D9,&
           2.940378956634553899906876D10,1.702665737765398868392998D11,&
           4.926125793377430887588120D11,5.606251856223951465078242D11/
      DATA Q4/2.690530175870899333379843D3,6.393885654300092398984238D5,&
             4.135599930241388052042842D7,1.120872109616147941376570D9,&
           1.488613728678813811542398D10,1.016803586272438228077304D11,&
           3.417476345507377132798597D11,4.463158187419713286462081D11/
!C----------------------------------------------------------------------
!C  Coefficients for minimax approximation over (12, INF).
!C----------------------------------------------------------------------
!CS    DATA C/-1.910444077728E-03,8.4171387781295E-04,
!CS   1     -5.952379913043012E-04,7.93650793500350248E-04,
!CS   2     -2.777777777777681622553E-03,8.333333333333333331554247E-02,
!CS   3      5.7083835261E-03/
      DATA C/-1.910444077728D-03,8.4171387781295D-04,  &
          -5.952379913043012D-04,7.93650793500350248D-04,  &
          -2.777777777777681622553D-03,8.333333333333333331554247D-02, &
           5.7083835261D-03/
!C----------------------------------------------------------------------
      Y = X
      IF ((Y .GT. ZERO) .AND. (Y .LE. XBIG)) THEN
            IF (Y .LE. EPS) THEN
                  RES = -LOG(Y)
               ELSE IF (Y .LE. THRHAL) THEN
!C----------------------------------------------------------------------
!C  EPS .LT. X .LE. 1.5
!C----------------------------------------------------------------------
                  IF (Y .LT. PNT68) THEN
                        CORR = -LOG(Y)
                        XM1 = Y
                     ELSE
                        CORR = ZERO
                        XM1 = (Y - HALF) - HALF
                  END IF
                  IF ((Y .LE. HALF) .OR. (Y .GE. PNT68)) THEN
                        XDEN = ONE
                        XNUM = ZERO
                        DO 140 I = 1, 8
                           XNUM = XNUM*XM1 + P1(I)
                           XDEN = XDEN*XM1 + Q1(I)
  140                   CONTINUE
                        RES = CORR + (XM1 * (D1 + XM1*(XNUM/XDEN)))
                     ELSE
                        XM2 = (Y - HALF) - HALF
                        XDEN = ONE
                        XNUM = ZERO
                        DO 220 I = 1, 8
                           XNUM = XNUM*XM2 + P2(I)
                           XDEN = XDEN*XM2 + Q2(I)
  220                   CONTINUE
                        RES = CORR + XM2 * (D2 + XM2*(XNUM/XDEN))
                  END IF
               ELSE IF (Y .LE. FOUR) THEN
!C----------------------------------------------------------------------
!C  1.5 .LT. X .LE. 4.0
!C----------------------------------------------------------------------
                  XM2 = Y - TWO
                  XDEN = ONE
                  XNUM = ZERO
                  DO 240 I = 1, 8
                     XNUM = XNUM*XM2 + P2(I)
                     XDEN = XDEN*XM2 + Q2(I)
  240             CONTINUE
                  RES = XM2 * (D2 + XM2*(XNUM/XDEN))
               ELSE IF (Y .LE. TWELVE) THEN
!C----------------------------------------------------------------------
!C  4.0 .LT. X .LE. 12.0
!C----------------------------------------------------------------------
                  XM4 = Y - FOUR
                  XDEN = -ONE
                  XNUM = ZERO
                  DO 340 I = 1, 8
                     XNUM = XNUM*XM4 + P4(I)
                     XDEN = XDEN*XM4 + Q4(I)
  340             CONTINUE
                  RES = D4 + XM4*(XNUM/XDEN)
               ELSE 
!C----------------------------------------------------------------------
!C  Evaluate for argument .GE. 12.0,
!C----------------------------------------------------------------------
                  RES = ZERO
                  IF (Y .LE. FRTBIG) THEN
                        RES = C(7)
                        YSQ = Y * Y
                        DO 450 I = 1, 6
                           RES = RES / YSQ + C(I)
  450                   CONTINUE
                  END IF
                  RES = RES/Y
                  CORR = LOG(Y)
                  RES = RES + SQRTPI - HALF*CORR
                  RES = RES + Y*(CORR-ONE)
            END IF
         ELSE
!C----------------------------------------------------------------------
!C  Return for bad arguments
!C----------------------------------------------------------------------
            RES = XINF
      END IF
!C----------------------------------------------------------------------
!C  Final adjustments and return
!C----------------------------------------------------------------------
!CS    ALGAMA = RES
      DLGAMMA = RES
      RETURN
!C ---------- Last line of DLGAMA ----------
      END FUNCTION DLGAMMA




      
!C
!C       ==========================================
!C       Purpose: Compute the density of beta distribution
!C       Input :  p  --- Parameter  ( p > 0 )
!C                q  --- Parameter  ( q > 0 )
!C       Output:   --- x^{p-1}*(1-x)^{q-1}/B(p,q)
!C       Routine called: GAMMA for computing â(x)
!C       ==========================================
!C

        real(dp) FUNCTION dbeta(x, p, q)
!c        implicit real*8 (a-h, o-z)
        use mod_kinds, only:dp 
        implicit none
        real(dp)::x, p, q
!c        real(dp),external:: beta       
        
        dbeta = x**(p-1d0)*(1d0-x)**(q-1d0)/beta(p,q)

        return 
        end  FUNCTION dbeta

        real(dp) FUNCTION BETA(P,Q)
!C
!C       ==========================================
!C       Purpose: Compute the beta function B(p,q)
!C       Input :  p  --- Parameter  ( p > 0 )
!C                q  --- Parameter  ( q > 0 )
!C       Output:  BT --- B(p,q)
!C       Routine called: GAMMA for computing â(x)
!C       ==========================================
!C
!c        IMPLICIT real*8 (A-H,O-Z)
        use mod_kinds, only:dp 
        implicit none
        real(dp):: p, q
!C        real(dp),external::QLGAMMA
        
        beta=exp(DLGAMMA(P)+DLGAMMA(Q)-DLGAMMA(P+Q))

        RETURN
        END FUNCTION BETA


!C=======================================================================
!C
! 


      real(dp) FUNCTION PSI(XX)
!C----------------------------------------------------------------------
!C
!C This function program evaluates the logarithmic derivative of the
!C   gamma function, 
!C
!C      psi(x) = d/dx (gamma(x)) / gamma(x) = d/dx (ln gamma(x))
!C
!C   for real x, where either
!C
!C          -xmax1 < x < -xmin (x not a negative integer), or
!C            xmin < x.
!C
!C   The calling sequence for this function is 
!C
!C                  Y = PSI(X)
!C
!C   The main computation uses rational Chebyshev approximations
!C   published in Math. Comp. 27, 123-127 (1973) by Cody, Strecok and
!C   Thacher.  This transportable program is patterned after the
!C   machine-dependent FUNPACK program PSI(X), but cannot match that
!C   version for efficiency or accuracy.  This version uses rational
!C   approximations that are theoretically accurate to 20 significant
!C   decimal digits.  The accuracy achieved depends on the arithmetic
!C   system, the compiler, the intrinsic functions, and proper selection
!C   of the machine-dependent constants.
!C
!C*******************************************************************
!C*******************************************************************
!C
!C Explanation of machine-dependent constants
!C
!C   XINF   = largest positive machine number
!C   XMAX1  = beta ** (p-1), where beta is the radix for the
!C            floating-point system, and p is the number of base-beta
!C            digits in the floating-point significand.  This is an
!C            upper bound on non-integral floating-point numbers, and
!C            the negative of the lower bound on acceptable negative
!C            arguments for PSI.  If rounding is necessary, round this
!C            value down.
!C   XMIN1  = the smallest in magnitude acceptable argument.  We
!C            recommend XMIN1 = MAX(1/XINF,xmin) rounded up, where
!C            xmin is the smallest positive floating-point number.
!C   XSMALL = absolute argument below which  PI*COTAN(PI*X)  may be
!C            represented by 1/X.  We recommend XSMALL < sqrt(3 eps)/pi,
!C            where eps is the smallest positive number such that
!C            1+eps > 1. 
!C   XLARGE = argument beyond which PSI(X) may be represented by
!C            LOG(X).  The solution to the equation
!C               x*ln(x) = beta ** p
!C            is a safe value.
!C
!C     Approximate values for some important machines are
!C
!C                        beta  p     eps     xmin       XINF  
!C
!C  CDC 7600      (S.P.)    2  48  7.11E-15  3.13E-294  1.26E+322
!C  CRAY-1        (S.P.)    2  48  7.11E-15  4.58E-2467 5.45E+2465
!C  IEEE (IBM/XT,
!C    SUN, etc.)  (S.P.)    2  24  1.19E-07  1.18E-38   3.40E+38
!C  IEEE (IBM/XT,
!C    SUN, etc.)  (D.P.)    2  53  1.11D-16  2.23E-308  1.79D+308
!C  IBM 3033      (D.P.)   16  14  1.11D-16  5.40D-79   7.23D+75
!C  SUN 3/160     (D.P.)    2  53  1.11D-16  2.23D-308  1.79D+308
!C  VAX 11/780    (S.P.)    2  24  5.96E-08  2.94E-39   1.70E+38
!C                (D.P.)    2  56  1.39D-17  2.94D-39   1.70D+38
!C   (G Format)   (D.P.)    2  53  1.11D-16  5.57D-309  8.98D+307
!C
!C                         XMIN1      XMAX1     XSMALL    XLARGE
!C
!C  CDC 7600      (S.P.)  3.13E-294  1.40E+14  4.64E-08  9.42E+12
!C  CRAY-1        (S.P.)  1.84E-2466 1.40E+14  4.64E-08  9.42E+12
!C  IEEE (IBM/XT,
!C    SUN, etc.)  (S.P.)  1.18E-38   8.38E+06  1.90E-04  1.20E+06
!C  IEEE (IBM/XT,
!C    SUN, etc.)  (D.P.)  2.23D-308  4.50D+15  5.80D-09  2.71D+14
!C  IBM 3033      (D.P.)  1.39D-76   4.50D+15  5.80D-09  2.05D+15
!C  SUN 3/160     (D.P.)  2.23D-308  4.50D+15  5.80D-09  2.71D+14
!C  VAX 11/780    (S.P.)  5.89E-39   8.38E+06  1.35E-04  1.20E+06
!C                (D.P.)  5.89D-39   3.60D+16  2.05D-09  2.05D+15
!C   (G Format)   (D.P.)  1.12D-308  4.50D+15  5.80D-09  2.71D+14
!C
!C*******************************************************************
!C*******************************************************************
!C
!C Error Returns
!C
!C  The program returns XINF for  X < -XMAX1, for X zero or a negative
!C    integer, or when X lies in (-XMIN1, 0), and returns -XINF
!C    when X lies in (0, XMIN1).
!C
!C Intrinsic functions required are:
!C
!C     ABS, AINT, DBLE, INT, LOG, REAL, TAN
!C
!C
!C  Author: W. J. Cody
!C          Mathematics and Computer Science Division 
!C          Argonne National Laboratory
!C          Argonne, IL 60439
!C
!C  Latest modification: June 8, 1988
!C
!C----------------------------------------------------------------------
!c     implicit real*8 (a-h,o-z)

      use mod_kinds,only: dp
      implicit none
      INTEGER I,N,NQ
!CS    REAL
!CD    DOUBLE PRECISION
      real(dp)::  &
        AUG,CONV,DEN,FOUR,FOURTH,HALF,ONE,P1,P2,PIOV4,Q1,Q2, &
        SGN,THREE,XLARGE,UPPER,W,X,XINF,XMAX1,XMIN1,XSMALL,X01,&
        X01D,X02,XX,Z,ZERO
      DIMENSION P1(9),P2(7),Q1(8),Q2(6)
!C----------------------------------------------------------------------
!C  Mathematical constants.  PIOV4 = pi / 4
!C----------------------------------------------------------------------
!CS    DATA ZERO,FOURTH,HALF,ONE/0.0E0,0.25E0,0.5E0,1.0E0/
!CS    DATA THREE,FOUR/3.0E0,4.0E0/,PIOV4/7.8539816339744830962E-01/
      DATA ZERO,FOURTH,HALF,ONE/0.0D0,0.25D0,0.5D0,1.0D0/
      DATA THREE,FOUR/3.0D0,4.0D0/,PIOV4/7.8539816339744830962D-01/
!C----------------------------------------------------------------------
!C  Machine-dependent constants
!C----------------------------------------------------------------------
!CS    DATA XINF/1.70E+38/, XMIN1/5.89E-39/, XMAX1/8.38E+06/,
!CS   1     XSMALL/1.35E-04/, XLARGE/1.20E+06/
      DATA XINF/1.70D+38/, XMIN1/5.89D-39/, XMAX1/3.60D+16/, &
          XSMALL/2.05D-09/, XLARGE/2.04D+15/
!C----------------------------------------------------------------------
!C  Zero of psi(x)
!C----------------------------------------------------------------------
!CS    DATA X01/187.0E0/,X01D/128.0E0/,X02/6.9464496836234126266E-04/
      DATA X01/187.0D0/,X01D/128.0D0/,X02/6.9464496836234126266D-04/
!C----------------------------------------------------------------------
!C  Coefficients for approximation to  psi(x)/(x-x0)  over [0.5, 3.0]
!C----------------------------------------------------------------------
!CS    DATA P1/4.5104681245762934160E-03,5.4932855833000385356E+00,
!CS   1        3.7646693175929276856E+02,7.9525490849151998065E+03,
!CS   2        7.1451595818951933210E+04,3.0655976301987365674E+05,
!CS   3        6.3606997788964458797E+05,5.8041312783537569993E+05,
!CS   4        1.6585695029761022321E+05/
!CS    DATA Q1/9.6141654774222358525E+01,2.6287715790581193330E+03,
!CS   1        2.9862497022250277920E+04,1.6206566091533671639E+05,
!CS   2        4.3487880712768329037E+05,5.4256384537269993733E+05,
!CS   3        2.4242185002017985252E+05,6.4155223783576225996E-08/
      DATA P1/4.5104681245762934160D-03,5.4932855833000385356D+00,&
             3.7646693175929276856D+02,7.9525490849151998065D+03,&
             7.1451595818951933210D+04,3.0655976301987365674D+05,&
             6.3606997788964458797D+05,5.8041312783537569993D+05,&
             1.6585695029761022321D+05/
      DATA Q1/9.6141654774222358525D+01,2.6287715790581193330D+03,&
             2.9862497022250277920D+04,1.6206566091533671639D+05,&
             4.3487880712768329037D+05,5.4256384537269993733D+05,&
             2.4242185002017985252D+05,6.4155223783576225996D-08/
!C----------------------------------------------------------------------
!C  Coefficients for approximation to  psi(x) - ln(x) + 1/(2x) 
!C     for  x > 3.0
!C----------------------------------------------------------------------
!CS    DATA P2/-2.7103228277757834192E+00,-1.5166271776896121383E+01,
!CS   1        -1.9784554148719218667E+01,-8.8100958828312219821E+00,
!CS   2        -1.4479614616899842986E+00,-7.3689600332394549911E-02,
!CS   3        -6.5135387732718171306E-21/
!CS    DATA Q2/ 4.4992760373789365846E+01, 2.0240955312679931159E+02,
!CS   1         2.4736979003315290057E+02, 1.0742543875702278326E+02,
!CS   2         1.7463965060678569906E+01, 8.8427520398873480342E-01/
      DATA P2/-2.7103228277757834192D+00,-1.5166271776896121383D+01,&
             -1.9784554148719218667D+01,-8.8100958828312219821D+00,&
             -1.4479614616899842986D+00,-7.3689600332394549911D-02,&
             -6.5135387732718171306D-21/
      DATA Q2/ 4.4992760373789365846D+01, 2.0240955312679931159D+02,&
              2.4736979003315290057D+02, 1.0742543875702278326D+02,&
              1.7463965060678569906D+01, 8.8427520398873480342D-01/
!C----------------------------------------------------------------------
!CS    CONV(I) = REAL(I)
      CONV(I) = DBLE(I)
      X = XX
      W = ABS(X)
      AUG = ZERO
!C----------------------------------------------------------------------
!C  Check for valid arguments, then branch to appropriate algorithm
!C----------------------------------------------------------------------
      IF ((-X .GE. XMAX1) .OR. (W .LT. XMIN1)) THEN
            GO TO 410
         ELSE IF (X .GE. HALF) THEN
            GO TO 200
!C----------------------------------------------------------------------
!C  X < 0.5, use reflection formula: psi(1-x) = psi(x) + pi * cot(pi*x)
!C     Use 1/X for PI*COTAN(PI*X)  when  XMIN1 < |X| <= XSMALL.  
!C----------------------------------------------------------------------
         ELSE IF (W .LE. XSMALL) THEN
            AUG = -ONE / X
            GO TO 150
      END IF
!C----------------------------------------------------------------------
!C  Argument reduction for cot
!C----------------------------------------------------------------------
      IF (X .LT. ZERO) THEN
            SGN = PIOV4
         ELSE
            SGN = -PIOV4
      END IF
      W = W - AINT(W)
      NQ = INT(W * FOUR)
      W = FOUR * (W - CONV(NQ) * FOURTH)
!C----------------------------------------------------------------------
!C  W is now related to the fractional part of  4.0 * X.
!C     Adjust argument to correspond to values in the first
!C     quadrant and determine the sign.
!C----------------------------------------------------------------------
      N = NQ / 2
      IF ((N+N) .NE. NQ) W = ONE - W
      Z = PIOV4 * W
      IF (MOD(N,2) .NE. 0) SGN = - SGN
!C----------------------------------------------------------------------
!C  determine the final value for  -pi * cotan(pi*x)
!C----------------------------------------------------------------------
      N = (NQ + 1) / 2
      IF (MOD(N,2) .EQ. 0) THEN
!C----------------------------------------------------------------------
!C  Check for singularity
!C----------------------------------------------------------------------
            IF (Z .EQ. ZERO) GO TO 410
            AUG = SGN * (FOUR / TAN(Z))
         ELSE
            AUG = SGN * (FOUR * TAN(Z))
      END IF
  150 X = ONE - X
  200 IF (X .GT. THREE) GO TO 300
!C----------------------------------------------------------------------
!C  0.5 <= X <= 3.0
!C----------------------------------------------------------------------
      DEN = X
      UPPER = P1(1) * X
      DO 210 I = 1, 7
         DEN = (DEN + Q1(I)) * X
         UPPER = (UPPER + P1(I+1)) * X
  210 CONTINUE
      DEN = (UPPER + P1(9)) / (DEN + Q1(8))
      X = (X-X01/X01D) - X02
      PSI = DEN * X + AUG
      GO TO 500
!C----------------------------------------------------------------------
!C  3.0 < X 
!C----------------------------------------------------------------------
  300 IF (X .LT. XLARGE) THEN
         W = ONE / (X * X)
         DEN = W
         UPPER = P2(1) * W
         DO 310 I = 1, 5
            DEN = (DEN + Q2(I)) * W
            UPPER = (UPPER + P2(I+1)) * W
  310    CONTINUE
         AUG = (UPPER + P2(7)) / (DEN + Q2(6)) - HALF / X + AUG
      END IF
      PSI = AUG + LOG(X)
      GO TO 500
!C----------------------------------------------------------------------
!C  Error return
!C----------------------------------------------------------------------
  410 PSI = XINF
      IF (X .GT. ZERO) PSI = -XINF
  500 RETURN
!C---------- Last card of PSI ----------
    END FUNCTION PSI

    end module mod_specialfun
