MODULE X_MOD
!-------------------------------------------------------------------------------
!X_MOD is an F90 module of cross-section and reaction rate routines.
!
!References:
!  W.A.Houlberg, F90 free format 8/2004
!
!Contains PUBLIC routines:
!  X_SIG_FUS           -fusion cross-sections
!  X_SIG_PE            -e impact ionization cross-section for H0,D0,T0
!  X_SIG_PI            -i impact ionization cross-section for H0,D0,T0
!  X_SIG_PX            -cx cross-section between H+,D+,T+ and H0,D0,T0
!  X_SIG_ZI            -cx + i impact ionization cross-section for Z>1 on
!                       H+,D+,T+
!  X_SIGV_FUS          -fusion reaction rates for Maxwellian distributions
!  X_SIGV_PE           -e impact ionization reaction rate for Maxwellian e,
!                       and monoenergetic H0,D0,T0
!  X_SIGV_PI           -i impact ionization reaction rate between Maxwellian
!                       H+,D+,T+ and monoenergetic H0,D0,T0
!  X_SIGV_PX           -cx reaction rate between Maxwellian H+,D+,T+ and
!                       monoenergetic H0,D0,T0
!  X_SIGV_RC           -recombination reaction rate for H+,D+,T+
!
!Comments:
!  Units are contained in [] in the description of variables
!-------------------------------------------------------------------------------
USE SPEC_KIND_MOD
IMPLICIT NONE 

!-------------------------------------------------------------------------------
! Private data
!-------------------------------------------------------------------------------
!Physical constants and conversion factors
REAL(KIND=rspec), PRIVATE, PARAMETER :: &  
  z_emass=9.1095e-31,  & !electron mass [kg]
  z_j7kv=1.6022e-16,   & !unit conversion factor [joules/keV]
  z_pi=3.141592654,    & !pi [-]
  z_pmass=1.6726e-27,  & !proton mass [kg]
  z_eion=0.0136          !H ionization energy [keV]

!Abscissas and weight factors for 10th order Gaussian integration
REAL(KIND=rspec), PRIVATE :: &
  x(1:10)=(/1.3779347e-1, 7.2945455e-1, 1.8083429e+0, 3.4014337e+0, &
            5.5524961e+0, 8.3301527e+0, 1.1843786e+1, 1.6279258e+1, &
            2.1996586e+1, 2.9920697e+1/),                           &
  w(1:10)=(/3.0844112e-1, 4.0111993e-1, 2.1806829e-1, 6.2087456e-2, &
            9.5015170e-3, 7.5300839e-4, 2.8259234e-5, 4.2493140e-7, &
            1.8395648e-9, 9.9118272e-13/)


!-------------------------------------------------------------------------------
! Procedures
!-------------------------------------------------------------------------------
CONTAINS
      
FUNCTION X_SIG_FUS(vrel,k_reaction) 
!-------------------------------------------------------------------------------
!X_SIG_FUS gives the fusion cross-section for several fusion reactions, given
!  the relative velocity of the incident ions
!
!References:
!  Pacific Northwest Laboratory Annual Report on Controlled Thermonuclear
!    Technology-1972 BNWL-1685
!  W.A.Houlberg, F90 free format 8/2004
!
!Comments:
!  The energies for the products are inversely proportional to their masses
!  The branching ratio for reaction 5 is energy sensitive
!-------------------------------------------------------------------------------

!Declaration of input variables
INTEGER, INTENT(IN) :: &
  k_reaction             !index denoting reaction [-]
                         !=1   D+D =>3He(  817keV)+  n( 2450keV)
                         !=2   D+D =>  T( 1008keV)+  p( 3024keV)
                         !=3   D+T =>4He( 3517keV)+  n(14069keV)
                         !=4   T+T =>4He( 1259keV)+  n( 5034keV)+n( 5034keV)
                         !=5 3He+T =>  D( 9546keV)+4He( 4773keV)              43%
                         !         =>  p( 5374keV)+4He( 1344keV)+n( 5374keV)  51%
                         !         =>  p(10077keV)+4He(  403keV)+n( 1612keV)   6%
                         !=6 3He+D =>  p(14681keV)+4He( 3670keV)

REAL(KIND=rspec), INTENT(IN) :: &
  vrel                   !relative velocity of reacting ions [m/s]

!Declaration of function output
REAL(KIND=rspec) :: &
  X_SIG_FUS              !fusion cross-section [m**2]

!Declaration of local variables
REAL(KIND=rspec), SAVE :: &
  a1(1:6)=(/47.877,  46.096,  45.948,  38.390,  123.11,  89.271 /), &                
  a2(1:6)=(/4.82e2,  3.72e2,  5.02e4,  4.48e2,  1.125e4, 2.59e4 /), &
  a3(1:6)=(/3.08e-4, 4.36e-4, 1.368e-2,1.02e-3, 0.0,     3.98e-3/), &                   
  a4(1:6)=(/1.177,   1.220,   1.076,   2.09,    0.0,     1.297  /), &
  a5(1:6)=(/0.0,     0.0,     4.09e2,  0.0,     0.0,     6.47e2 /), &
  am(1:6)=(/2.0,     2.0,     2.0,     3.0,     3.0,     2.0    /)

REAL(KIND=rspec) :: &
  ekev,ex0,sigb

!-------------------------------------------------------------------------------
!Evaluate cross-section
!-------------------------------------------------------------------------------
!Initialization
X_SIG_FUS=0

!Conversion of relative velocity to keV
ekev=am(k_reaction)*z_pmass*vrel**2/2/z_j7kv

IF(ekev >= 0.03) THEN

  ex0=a1(k_reaction)/SQRT(ekev)

  IF(ex0 > LOG(HUGE(1.0_rspec))) THEN
 
    !Set cross-section to 0 if expansion term is too large
    X_SIG_FUS=0
 
  ELSE
 
    sigb=(a5(k_reaction)+a2(k_reaction)/(1.0+(a3(k_reaction)*ekev              &
         -a4(k_reaction))**2))/(ekev*(EXP(ex0)-1.0))

    IF(sigb < 1.0e-10) THEN

      !Set cross-section to 0 if less than machine precision
      X_SIG_FUS=0

    ELSE

      X_SIG_FUS=1.0e-28*sigb

    ENDIF

  ENDIF

ENDIF

END FUNCTION X_SIG_FUS

FUNCTION X_SIG_PE(e)
!-------------------------------------------------------------------------------
!X_SIG_PE gives the cross-section for electron impact ionization of neutral
!  hydrogenic atoms
!
!References:
!  Freeman, Jones, Culham Laboratory Report, CLM-R137 (1974)
!  Gryzinski, Phys Rev 138, A305, A322, A336 (1965)
!  W.A.Houlberg, F90 free format 8/2004
!-------------------------------------------------------------------------------

!Declaration of input variables
REAL(KIND=rspec), INTENT(IN) :: &
  e                      !electron energy [keV]

!Declaration of function output
REAL(KIND=rspec):: &
  X_SIG_PE               !electron impact ionization cross-section [m**2]

!Declaration of local variables
REAL(KIND=rspec) :: &
  er

!-------------------------------------------------------------------------------
!Evaluate cross-section
!-------------------------------------------------------------------------------
!Initialization
X_SIG_PE=0

!Normalize energy to ionization energy
er=e/z_eion

IF(er > 1.0) THEN

  X_SIG_PE=3.513e-20/er*((er-1.0)/(er+1.0))**1.5                               &
           *(1.0+2.0/3.0*(1.0-0.5/er)*LOG(2.7+SQRT(er-1.0)))

ENDIF

END FUNCTION X_SIG_PE

FUNCTION X_SIG_PI(e,ai)
!-------------------------------------------------------------------------------
!X_SIG_PI gives the cross-section for ion impact ionization of neutral
!  hydrogenic atoms
!
!References:
!  Freeman, Jones, Culham Laboratory Report, CLM-R137 (1974)
!  W.A.Houlberg, F90 free format 8/2004
!-------------------------------------------------------------------------------

!Declaration of input variables
REAL(KIND=rspec), INTENT(IN) :: &
  ai,                  & !atomic mass number of ion [1-3]
  e                      !ion energy (stationary neutral) [keV]

!Declaration of function output
REAL(KIND=rspec) :: &
  X_SIG_PI               !ion impact ionization cross-section [m**2]

!Declaration of local variables
REAL(KIND=rspec), SAVE :: &
  a(1:7)=(/-5.124343e+1,  3.557321e+0, -1.045134e+0,  3.139238e-1, &
           -7.454475e-2,  8.459113e-3, -3.495444e-4/)

REAL(KIND=rspec) :: &
  ep,epl

!-------------------------------------------------------------------------------
!Evaluate cross-section
!-------------------------------------------------------------------------------
!Initialization
X_SIG_PI=0

!Define effective proton energy
ep=e/ai

IF(e > z_eion .AND. ep <= 1.0) THEN

  !Fit valid for z_eion < ep <= 1 keV
  X_SIG_PI=5.56e-23*ep**2.55

ELSEIF(ep > 1.0 .AND. ep <= 5.0e2) THEN

  !Fit from freeman and jones valid for 0.1 keV < ep <= 500 keV
  epl=LOG(ep)
  X_SIG_PI=EXP(a(1)+epl*(a(2)+epl*(a(3)+epl*(a(4)+epl*(a(5) &
               +epl*(a(6)+epl*a(7)))))))

ELSEIF(ep > 5.0e2) THEN

  !Fit valid for ep > 500 keV
  X_SIG_PI=7.37e-19/ep**0.86

ENDIF

END FUNCTION X_SIG_PI

FUNCTION X_SIG_PX(e,ai)
!-------------------------------------------------------------------------------
!X_SIG_PX gives the cross-section for charge exchange between hydrogenic species
!
!References:
!  Freeman, Jones, Culham Laboratory Report, CLM-R137 (1974)
!  W.A.Houlberg, F90 free format 8/2004
!-------------------------------------------------------------------------------

!Declaration of input variables
REAL(KIND=rspec), INTENT(IN):: &
  ai,                  & !atomic mass number of ion [1-3]
  e                      !ion energy (stationary neutral) [keV]

!Declaration of function output
REAL(KIND=rspec) :: &
  X_SIG_PX               !charge exchange cross-section [m**2]

!Declaration of local variables
REAL(KIND=rspec) :: &
  ep,epe

!-------------------------------------------------------------------------------
!Evaluate cross-section
!-------------------------------------------------------------------------------
!Initialization
X_SIG_PX=0

!Define effective proton energy
ep=e/ai

IF(ep <= 1.0e-2) THEN

  !Fit valid for ep <= .01 keV
  X_SIG_PX=2.55e-19/ep**0.146

ELSEIF(ep <= 1.0e2) THEN

  !Fit from Freeman and Jones valid for .001 keV < ep <= 100 keV
  epe=1.0e3*ep
  X_SIG_PX=6.937e-19*(1.0-6.73e-2*LOG(epe))**2/(1.0+1.112e-15*epe**3.3)

ELSE

  !Fit valid for ep > 100 keV
  X_SIG_PX=2.28e-14/ep**3.67

ENDIF

END FUNCTION X_SIG_PX

FUNCTION X_SIG_ZI(e,an,zi)
!-------------------------------------------------------------------------------
!X_SIG_ZI calculates the cross section for charge exchange plus ion impact
!  ionization of species Z>1 with atomic hydrogen
!
!References:
!  Olson, et al, PRL 41, 163
!  W.A.Houlberg, F90 free format 8/2004
!-------------------------------------------------------------------------------

!Declaration of input variables
REAL(KIND=rspec), INTENT(IN) :: &
  an,                  & !atomic mass number of neutral [1-3]
  e,                   & !relative energy [keV]
  zi                     !charge state of impurity [-]

!Declaration of function output
REAL(KIND=rspec) :: &
  X_SIG_ZI               !charge exchange plus ion impact cross-section [m**2]

!Declaration of local variables
REAL(KIND=rspec) :: &
  con1,con2

!-------------------------------------------------------------------------------
!Evaluate cross-section
!-------------------------------------------------------------------------------
!Initialization
X_SIG_ZI=0

con1=e/(32*zi*an)

IF(con1 > 87.0) THEN

  con2=1.0

ELSE

  con2=1.0-EXP(-con1)

ENDIF

X_SIG_ZI=4.6e-20*zi*con2/con1

END FUNCTION X_SIG_ZI

FUNCTION X_SIGV_FUS(ti,k_reaction)
!-------------------------------------------------------------------------------
!X_SIGV_FUS gives <sigma*v> integrated over two Maxwellian distributions each
!  with temperature ti, for several fusion reactions
!
!References:
!  McNally, Rothe, Sharp, ORNL/TM-6914 (1979)
!  Pacific Northwest Laboratory Annual Report on Controlled
!   Thermonuclear Technology, BNWL-1685 (1972)
!  UCRL-70552 (1967)
!  Sov Phys JETP 12 (1961) 163
!  W.A.Houlberg, F90 free format 8/2004
!
!Comments:
!  The energies for the products are inversely proportional to the
!  mass of each product
!  The branching ratio for reaction 5 is energy sensitive
!  Cross-sections tabulated in steps of .2 keV for   1 <t< 10  keV
!                                       2. keV for  10 <t<100  keV
!                                      20. keV for 100 <t<200  keV
!  D-T reaction below 1 keV valid to 10**-56 m**3/s
!-------------------------------------------------------------------------------

!Declaration of input variables
INTEGER, INTENT(IN) :: &
  k_reaction             !index denoting reaction [-]
                         !=1   D+D =>3He(  817keV)+  n( 2450keV)
                         !=2   D+D =>  T( 1008keV)+  p( 3024keV)
                         !=3   D+T =>4He( 3517keV)+  n(14069keV)
                         !=4   T+T =>4He( 1259keV)+  n( 5034keV)+n( 5034keV)
                         !=5 3He+T =>  D( 9546keV)+4He( 4773keV)             41%
                         !         =>  p( 5374keV)+4He( 1344keV)+n( 5374keV) 55%
                         !         =>  p(10077keV)+4He(  403keV)+n( 1612keV)  4%
                         !=6 3He+D =>  p(14681keV)+4He( 3670keV)

REAL(KIND=rspec), INTENT(IN) :: &
  ti                     !Maxwellian temperature [keV]

!Declaration of function output
REAL(KIND=rspec) :: &
  X_SIGV_FUS             !fusion rate <sigma*v> [m**3/s]

!Declaration of local variables
INTEGER :: &
  ipt,ix

REAL(KIND=rspec) :: &
  dtx,fsv,fsvp1,svb,tkev,tl,tx

REAL(KIND=rspec), SAVE :: &
  fsv1(1:96)=(/9.65e-29, 2.61e-28, 5.75e-28, 1.10e-27, 1.90e-27, &
               3.04e-27, 4.58e-27, 6.57e-27, 9.06e-27, 1.21e-26, &
               1.57e-26, 2.00e-26, 2.49e-26, 3.04e-26, 3.67e-26, &
               4.37e-26, 5.14e-26, 5.98e-26, 6.90e-26, 7.90e-26, &
               8.97e-26, 1.01e-25, 1.13e-25, 1.26e-25, 1.40e-25, &
               1.55e-25, 1.70e-25, 1.86e-25, 2.03e-25, 2.21e-25, & 
               2.39e-25, 2.58e-25, 2.78e-25, 2.99e-25, 3.20e-25, &
               3.42e-25, 3.64e-25, 3.88e-25, 4.12e-25, 4.37e-25, &
               4.62e-25, 4.88e-25, 5.15e-25, 5.42e-25, 5.70e-25, &
               5.99e-25, 9.18e-25, 1.29e-24, 1.71e-24, 2.16e-24, &
               2.65e-24, 3.17e-24, 3.71e-24, 4.27e-24, 4.85e-24, &
               5.44e-24, 6.04e-24, 6.66e-24, 7.28e-24, 7.90e-24, &
               8.54e-24, 9.18e-24, 9.82e-24, 1.05e-23, 1.11e-23, &
               1.18e-23, 1.24e-23, 1.30e-23, 1.37e-23, 1.43e-23, &
               1.50e-23, 1.56e-23, 1.63e-23, 1.69e-23, 1.76e-23, &
               1.82e-23, 1.88e-23, 1.95e-23, 2.01e-23, 2.07e-23, &
               2.13e-23, 2.20e-23, 2.26e-23, 2.32e-23, 2.38e-23, &
               2.44e-23, 2.50e-23, 2.56e-23, 2.62e-23, 2.68e-23, &
               2.74e-23, 3.32e-23, 3.86e-23, 4.37e-23, 4.86e-23, &
               5.32e-23/)

REAL(KIND=rspec), SAVE :: &
  fsv2(1:96)=(/9.66e-29, 2.61e-28, 5.76e-28, 1.10e-27, 1.91e-27, &
               3.05e-27, 4.59e-27, 6.58e-27, 9.07e-27, 1.21e-26, &
               1.57e-26, 1.99e-26, 2.48e-26, 3.04e-26, 3.66e-26, &
               4.35e-26, 5.11e-26, 5.95e-26, 6.86e-26, 7.84e-26, &
               8.90e-26, 1.00e-25, 1.12e-25, 1.25e-25, 1.39e-25, &
               1.53e-25, 1.68e-25, 1.83e-25, 2.00e-25, 2.17e-25, &
               2.35e-25, 2.53e-25, 2.72e-25, 2.92e-25, 3.12e-25, &
               3.33e-25, 3.55e-25, 3.77e-25, 4.00e-25, 4.23e-25, &
               4.48e-25, 4.72e-25, 4.97e-25, 5.23e-25, 5.49e-25, &
               5.76e-25, 8.72e-25, 1.21e-24, 1.58e-24, 1.99e-24, &
               2.41e-24, 2.85e-24, 3.31e-24, 3.79e-24, 4.27e-24, &
               4.76e-24, 5.25e-24, 5.75e-24, 6.26e-24, 6.77e-24, &
               7.28e-24, 7.79e-24, 8.30e-24, 8.81e-24, 9.32e-24, &
               9.84e-24, 1.03e-23, 1.09e-23, 1.14e-23, 1.19e-23, &
               1.24e-23, 1.29e-23, 1.34e-23, 1.39e-23, 1.44e-23, &
               1.49e-23, 1.54e-23, 1.59e-23, 1.64e-23, 1.69e-23, &
               1.73e-23, 1.78e-23, 1.83e-23, 1.88e-23, 1.93e-23, &
               1.97e-23, 2.02e-23, 2.07e-23, 2.12e-23, 2.16e-23, &
               2.21e-23, 2.66e-23, 3.09e-23, 3.51e-23, 3.90e-23, &
               4.29e-23/)

REAL(KIND=rspec), SAVE :: &
  fsv3(1:96)=(/6.27e-27, 1.86e-26, 4.44e-26, 9.11e-26, 1.67e-25, &
               2.83e-25, 4.47e-25, 6.72e-25, 9.67e-25, 1.34e-24, &
               1.81e-24, 2.38e-24, 3.06e-24, 3.86e-24, 4.79e-24, &
               5.86e-24, 7.07e-24, 8.43e-24, 9.95e-24, 1.16e-23, &
               1.35e-23, 1.55e-23, 1.77e-23, 2.00e-23, 2.26e-23, &
               2.53e-23, 2.81e-23, 3.12e-23, 3.44e-23, 3.78e-23, &
               4.14e-23, 4.52e-23, 4.91e-23, 5.31e-23, 5.74e-23, &
               6.17e-23, 6.63e-23, 7.09e-23, 7.57e-23, 8.07e-23, &
               8.57e-23, 9.09e-23, 9.62e-23, 1.02e-22, 1.07e-22, &
               1.13e-22, 1.74e-22, 2.39e-22, 3.06e-22, 3.70e-22, &
               4.31e-22, 4.88e-22, 5.39e-22, 5.86e-22, 6.28e-22, &
               6.65e-22, 6.98e-22, 7.27e-22, 7.52e-22, 7.74e-22, &
               7.93e-22, 8.09e-22, 8.23e-22, 8.35e-22, 8.45e-22, &
               8.54e-22, 8.61e-22, 8.66e-22, 8.70e-22, 8.74e-22, &
               8.76e-22, 8.77e-22, 8.78e-22, 8.78e-22, 8.77e-22, &
               8.76e-22, 8.74e-22, 8.72e-22, 8.70e-22, 8.67e-22, &
               8.64e-22, 8.61e-22, 8.57e-22, 8.54e-22, 8.50e-22, &
               8.46e-22, 8.42e-22, 8.37e-22, 8.33e-22, 8.29e-22, &
               8.24e-22, 7.77e-22, 7.31e-22, 6.89e-22, 6.50e-22, &
               6.16e-22/)

REAL(KIND=rspec), SAVE :: &
  fsv4(1:96)=(/3.28e-29, 1.02e-28, 2.52e-28, 5.29e-28, 9.84e-28, &
               1.68e-27, 2.66e-27, 4.00e-27, 5.75e-27, 7.96e-27, &
               1.07e-26, 1.40e-26, 1.79e-26, 2.24e-26, 2.76e-26, &
               3.35e-26, 4.01e-26, 4.75e-26, 5.56e-26, 6.45e-26, &
               7.42e-26, 8.47e-26, 9.60e-26, 1.08e-25, 1.21e-25, &
               1.35e-25, 1.49e-25, 1.64e-25, 1.80e-25, 1.97e-25, &
               2.14e-25, 2.32e-25, 2.51e-25, 2.71e-25, 2.91e-25, &
               3.12e-25, 3.34e-25, 3.56e-25, 3.79e-25, 4.03e-25, &
               4.27e-25, 4.52e-25, 4.77e-25, 5.03e-25, 5.30e-25, &
               5.57e-25, 8.54e-25, 1.19e-24, 1.56e-24, 1.96e-24, &
               2.37e-24, 2.80e-24, 3.24e-24, 3.68e-24, 4.13e-24, &
               4.59e-24, 5.04e-24, 5.50e-24, 5.96e-24, 6.41e-24, &
               6.87e-24, 7.33e-24, 7.78e-24, 8.23e-24, 8.68e-24, &
               9.12e-24, 9.57e-24, 1.00e-23, 1.04e-23, 1.09e-23, &
               1.13e-23, 1.17e-23, 1.22e-23, 1.26e-23, 1.30e-23, &
               1.34e-23, 1.39e-23, 1.43e-23, 1.47e-23, 1.51e-23, &
               1.55e-23, 1.59e-23, 1.63e-23, 1.67e-23, 1.71e-23, &
               1.75e-23, 1.79e-23, 1.83e-23, 1.87e-23, 1.91e-23, &
               1.95e-23, 2.35e-23, 2.79e-23, 3.28e-23, 3.78e-23, &
               4.27e-23/)

REAL(KIND=rspec), SAVE :: &
  fsv5(1:96)=(/0.00e-00, 0.00e-00, 0.00e-00, 0.00e-00, 0.00e-00, &
               0.00e-00, 0.00e-00, 0.00e-00, 0.00e-00, 0.00e-00, &
               0.00e-00, 0.00e-00, 0.00e-00, 0.00e-00, 0.00e-00, &
               0.00e-00, 0.00e-00, 0.00e-00, 0.00e-00, 0.00e-00, &
               0.00e-00, 0.00e-00, 0.00e-00, 0.00e-00, 0.00e-00, &
               0.00e-00, 0.00e-00, 0.00e-00, 0.00e-00, 0.00e-00, &
               0.00e-00, 0.00e-00, 0.00e-00, 0.00e-00, 0.00e-00, &
               0.00e-00, 0.00e-00, 0.00e-00, 0.00e-00, 0.00e-00, &
               0.00e-00, 0.00e-00, 0.00e-00, 0.00e-00, 0.00e-00, &
               1.16e-26, 2.50e-26, 5.12e-26, 9.74e-26, 1.69e-25, &
               2.62e-25, 3.77e-25, 5.24e-25, 7.04e-25, 9.12e-25, &
               1.13e-24, 1.38e-24, 1.68e-24, 2.03e-24, 2.41e-24, &
               2.81e-24, 3.22e-24, 3.69e-24, 4.20e-24, 4.74e-24, &
               5.29e-24, 5.86e-24, 6.48e-24, 7.14e-24, 7.82e-24, &
               8.52e-24, 9.23e-24, 9.99e-24, 1.08e-23, 1.16e-23, &
               1.24e-23, 1.32e-23, 1.41e-23, 1.50e-23, 1.60e-23, &
               1.69e-23, 1.78e-23, 1.87e-23, 1.97e-23, 2.07e-23, &
               2.18e-23, 2.29e-23, 2.39e-23, 2.50e-23, 2.61e-23, &
               2.71e-23, 3.81e-23, 5.09e-23, 6.46e-23, 7.86e-23, &
               9.18e-23/)

REAL(KIND=rspec), SAVE :: &
  fsv6(1:96)=(/3.10e-32, 1.80e-31, 7.28e-31, 2.31e-30, 6.11e-30, &
               1.41e-29, 2.94e-29, 5.63e-29, 1.01e-28, 1.70e-28, &
               2.73e-28, 4.22e-28, 6.29e-28, 9.09e-28, 1.28e-27, &
               1.75e-27, 2.36e-27, 3.11e-27, 4.03e-27, 5.14e-27, &
               6.46e-27, 8.03e-27, 9.87e-27, 1.20e-26, 1.45e-26, &
               1.73e-26, 2.05e-26, 2.41e-26, 2.81e-26, 3.26e-26, &
               3.76e-26, 4.32e-26, 4.93e-26, 5.61e-26, 6.35e-26, &
               7.15e-26, 8.03e-26, 8.98e-26, 1.00e-25, 1.11e-25, &
               1.23e-25, 1.36e-25, 1.50e-25, 1.65e-25, 1.81e-25, &
               1.97e-25, 4.33e-25, 8.17e-25, 1.39e-24, 2.19e-24, &
               3.26e-24, 4.61e-24, 6.28e-24, 8.27e-24, 1.06e-23, &
               1.32e-23, 1.62e-23, 1.95e-23, 2.30e-23, 2.68e-23, &
               3.09e-23, 3.52e-23, 3.96e-23, 4.42e-23, 4.89e-23, &
               5.37e-23, 5.87e-23, 6.36e-23, 6.87e-23, 7.37e-23, &
               7.88e-23, 8.38e-23, 8.88e-23, 9.38e-23, 9.87e-23, &
               1.04e-22, 1.08e-22, 1.13e-22, 1.18e-22, 1.22e-22, &
               1.27e-22, 1.31e-22, 1.36e-22, 1.40e-22, 1.44e-22, &
               1.48e-22, 1.52e-22, 1.56e-22, 1.60e-22, 1.63e-22, &
               1.67e-22, 1.97e-22, 2.19e-22, 2.34e-22, 2.45e-22, &
               2.52e-22/) 

REAL(KIND=rspec), SAVE :: &
  a(1:4)=(/-4.137748e+0,-5.416161e+0, 9.860669e-1,-2.094218e-1/)

!-------------------------------------------------------------------------------
!Evaluate reaction rate
!-------------------------------------------------------------------------------
!Initialization
X_SIGV_FUS=0
tkev=ti

IF(ti < 1.0) THEN
  
  !Lower cutoff on rates at 1 keV except D-T 
  IF(k_reaction == 3) THEN

    !D-T reaction below 1 keV valid to 10**-56 m**3/s
    tl=LOG(ti)
    svb=EXP(-(a(1)+tl*(a(2)+tl*(a(3)+tl*a(4)))))
    X_SIGV_FUS=1.0e-28*svb

  ENDIF

ELSE
  
  !Interpolate tabulated values
  IF(ti < 10.0) THEN

    ! 1<ti<10 tabulated every 0.2 keV
    dtx=0.2
    ipt=-4

  ELSEIF(ti < 100.0) THEN

    ! 10<ti<100 tabulated every 2 keV
    dtx=2.0
    ipt=41

  ELSE

    ! 100<ti<200 tabulated every 20 keV
    dtx=20.0
    ipt=86

    IF(ti >= 199.999) THEN

      ! 200<ti cutoff at 200 keV
      tkev=199.999

    ENDIF
     
  ENDIF

  ix=INT(tkev/dtx)+ipt
  tx=dtx*(ix-ipt)

  IF(k_reaction == 1) THEN

    fsv=fsv1(ix)
    fsvp1=fsv1(ix+1)

  ELSEIF(k_reaction == 2) THEN

    fsv=fsv2(ix)
    fsvp1=fsv2(ix+1)

  ELSEIF(k_reaction == 3) THEN

    fsv=fsv3(ix)
    fsvp1=fsv3(ix+1)

  ELSEIF(k_reaction == 4) THEN

    fsv=fsv4(ix)
    fsvp1=fsv4(ix+1)

  ELSEIF(k_reaction == 5) THEN

    fsv=fsv5(ix)
    fsvp1=fsv5(ix+1)

  ELSEIF(k_reaction == 6) THEN

    fsv=fsv6(ix)
    fsvp1=fsv6(ix+1)

  ENDIF

  X_SIGV_FUS=fsv+(fsvp1-fsv)*(tkev-tx)/dtx

ENDIF

END FUNCTION X_SIGV_FUS

FUNCTION X_SIGV_PE(te,en,an)
!-------------------------------------------------------------------------------
!X_SIGV_PE gives <sigma*v> for electron impact ionization of neutral hydrogenic
!  atoms by integrating over a Maxwellian electron distribution and a
!  monoenergetic neutral distribution
!
!References:
!  Freeman, Jones, Culham Laboratory Report, CLM-R137 (1974)
!  Gryzinski, Phys Rev 138 a305,a322,a336 (1965)
!  W.A.Houlberg, F90 free format 8/2004
!-------------------------------------------------------------------------------

!Declaration of input variables
REAL(KIND=rspec), INTENT(IN) :: &
  an,                  & !atomic mass number of neutral [1-3]
  en,                  & !neutral hydrogen atom energy [keV]
  te                     !electron temperature [keV]

!Declaration of function output
REAL(KIND=rspec) :: &
  X_SIGV_PE              !<sigma*v> [m**3/s]

!-------------------------------------------------------------------------------
!Declaration of local variables
INTEGER :: &
  i

REAL(KIND=rspec) :: &
  er,svint,tl,u,ve,vn,x12       

REAL(KIND=rspec), SAVE :: & 
  a(1:7)=(/-3.140212e+1,-3.024379e-1,-5.616546e-2, 7.902886e-3, &
            1.246713e-3, 2.217222e-4,-9.486967e-5/)                

!-------------------------------------------------------------------------------
!Evaluate reaction rate
!-------------------------------------------------------------------------------
!Initialize reaction rate
X_SIGV_PE=0

!Physical and conversion constants
ve=SQRT(2*te*z_j7kv/z_emass)
vn=SQRT(2*en*z_j7kv/(z_pmass*an))
u=vn/ve

!Check dimensionless velocities
IF(u < 0.3) THEN

  !vn < 0.3*ve, ignore vn
  IF(te <= 1.0e2) THEN

    !Fit from Freeman and Jones valid for .001 keV < te < 100 keV
    tl=LOG(te)
    X_SIGV_PE=EXP(a(1)+tl*(a(2)+tl*(a(3)+tl*(a(4)+tl*(a(5)+tl*(a(6) &
                  +tl*a(7)))))))

  ELSEIF(te > 1.0e2) THEN

    !Integrate for te > 100 keV
    svint=0

    DO i=1,10

      er=te*x(i)
      svint=svint*w(i)*x(i)*X_SIG_PE(er)

    ENDDO
        
    X_SIGV_PE=svint*2*ve/SQRT(z_pi)

  ENDIF

ELSEIF(u < 4.0) THEN

  !Comparable velocities, 0.3*ve < vn < 4.0*ve
  svint=0

  DO i=1,10

    er=te*x(i)
    x12=SQRT(x(i))
    svint=svint+w(i)*x12*SINH(2*u*x12)*X_SIG_PE(er)

  ENDDO   

  X_SIGV_PE=svint*ve*EXP(-u**2)/SQRT(z_pi)/u

ELSE

  !Ignore ve, vn > 4.0*ve
  er=en*z_emass/(an*z_pmass)
  X_SIGV_PE=X_SIG_PE(er)*vn

ENDIF

END FUNCTION X_SIGV_PE

FUNCTION X_SIGV_PI(ti,en,ai,an)
!-------------------------------------------------------------------------------
!X_SIGV_PI gives <sigma*v> for ion impact ionization of neutral hydrogenic atoms
!  by integrating over a Maxwellian distribution for the ions and a
!  monoenergetic distribution for the neutrals
!
!References:
!  Freeman, Jones, Culham Laboratory Report, CLM-R137 (1974)
!  W.A.Houlberg, F90 free format 8/2004
!-------------------------------------------------------------------------------

!Declaration of input variables
REAL(KIND=rspec), INTENT(IN) :: &
  ai,                  & !atomic mass number of maxwellian component [1-3]
  an,                  & !atomic mass number of monoenergetic component [1-3]
  en,                  & !neutral hydrogen atom energy [keV]
  ti                     !ion temperature [keV]

!Declaration of function output
REAL(KIND=rspec) :: &
  X_SIGV_PI              !<sigma*v> [m**3/s]

!Declaration of local variables
INTEGER :: &
  i

REAL(KIND=rspec) :: &
  er,svint,u,vi,vn,x12,y

REAL(KIND=rspec), SAVE :: &
  a(1:4)=(/-5.900469e+1, 8.728930e+0,-5.545365e-1, 9.956948e-3/)
    
!-------------------------------------------------------------------------------
!Evaluate reaction rate
!-------------------------------------------------------------------------------
!Initialize function
X_SIGV_PI=0

!Physical and conversion constants
vi=SQRT(2*ti*z_j7kv/(z_pmass*ai))
vn=SQRT(2*en*z_j7kv/(z_pmass*an))
u=vn/vi

!Check dimensionless velocities
IF(u < 0.2) THEN

  !Ignore vn, vn < 0.2*vi
  IF((ti/ai) >= 5.0e-1 .AND. (ti/ai) <= 1.0e2) THEN

    !Fit valid for 0.5 keV < ti/ai < 100 keV
    y=LOG(1.0e3*ti/ai)
    X_SIGV_PI=1.0e-6*EXP(a(1)+y*(a(2)+y*(a(3)+y*a(4))))

  ELSE

    !Integrate for ti/ai < 0.5 keV and ti/ai > 100 keV
    svint=0

    DO i=1,10

      er=ti*x(i)
      svint=svint+w(i)*x(i)*X_SIG_PI(er,ai)

    ENDDO
   
    X_SIGV_PI=svint*2*vi/SQRT(z_pi)

  ENDIF

ELSEIF(u < 5.0) THEN

  !Comparable velocities, 0.2*vi < vn < 5.0*vi
  svint=0

  DO i=1,10

    er=ti*x(i)
    x12=SQRT(x(i))
    svint=svint+w(i)*x12*SINH(2*u*x12)*X_SIG_PI(er,ai)

  ENDDO
   
  X_SIGV_PI=svint*vi*EXP(-u**2)/SQRT(z_pi)/u

ELSE

  !Ignore vi, 5.0*vi < vn
  X_SIGV_PI=X_SIG_PI(en,an)*vn

ENDIF

END FUNCTION X_SIGV_PI

FUNCTION X_SIGV_PX(ti,en,ai,an)
!-------------------------------------------------------------------------------
!X_SIGV_PX gives <sigma*v> for charge exchange between hydrogenic species by
!  integrating over a Maxwellian distribution for the ions and a monoenergetic
!  distribution for the neutrals
!
!References:
!  Freeman, Jones, Culham Laboratory Report, CLM-R137 (1974)
!  W.A.Houlberg, F90 free format 8/2004
!-------------------------------------------------------------------------------

!Declaration of input variables
REAL(KIND=rspec), INTENT(IN) :: &
  ai,                  & !atomic mass number of Maxwellian component [1-3]
  an,                  & !atomic mass number of monoenergetic component [1-3]
  en,                  & !neutral hydrogen atom energy [keV]
  ti                     !Maxwellian ion temperature [keV]

!Declaration of function output
REAL(KIND=rspec) :: &
  X_SIGV_PX              !<sigma*v> [m**3/s]

!Declaration of local variables
INTEGER :: &
  i

REAL(KIND=rspec) :: &
  er,svint,u,vi,vn,x12,y

REAL(KIND=rspec), SAVE :: & 
  a(1:9)=(/-3.220815e+1, 3.479700e-1,-1.196432e-2, 7.548590e-3, &
            5.301262e-5,-9.399207e-4, 2.215189e-4,-1.952323e-5, &
            5.947583e-7/)            
     
!-------------------------------------------------------------------------------
!Evaluate reaction rate
!-------------------------------------------------------------------------------
!Initialize reaction rate
X_SIGV_PX=0

!Physical and conversion constants
vi=SQRT(2*ti*z_j7kv/(z_pmass*ai))
vn=SQRT(2*en*z_j7kv/(z_pmass*an))
u=vn/vi
      
!Check dimensionless velocities
IF(u < 0.2) THEN

  !Ignore vn, vn < vi/5
  IF((ti/ai) >= 1.0e-4 .AND. (ti/ai) <= 1.0e2) THEN

    !Fit valid for .0001 keV < ti/ai < 100 keV
    y=LOG(1.0e3*ti/ai)
    X_SIGV_PX=EXP(a(1)+y*(a(2)+y*(a(3)+y*(a(4)+y*(a(5)+y*(a(6)+y*(a(7) &
                  +y*(a(8)+y*a(9)))))))))

  ELSE

    !Integrate for ti/ai < 0.0001 keV and ti/ai > 100 keV
    svint=0

    DO i=1,10

      er=ti*x(i)
      svint=svint+w(i)*x(i)*X_SIG_PX(er,ai)

    ENDDO
   
    X_SIGV_PX=svint*2*vi/SQRT(z_pi)

  ENDIF

ELSEIF(u < 5.0) THEN 

  !Comparable velocities, vi/5 < vn < 5*vi
  svint=0

  DO i=1,10

    er=ti*x(i)
    x12=SQRT(x(i))
    svint=svint+w(i)*x12*SINH(2*u*x12)*X_SIG_PX(er,ai)

  ENDDO
   
  X_SIGV_PX=svint*vi*EXP(-u**2)/SQRT(z_pi)/u

ELSE

  !Ignore vi for vi < vn/5
  X_SIGV_PX=X_SIG_PX(en,an)*vn

ENDIF

END FUNCTION X_SIGV_PX

FUNCTION X_SIGV_RC(te)
!-------------------------------------------------------------------------------
!X_SIGV_RC gives <sigma*v> for recombination of hydrogen with electrons
!
!References:
!  Gordeev, et al. JETP Lett, 25, 204
!  W.A.Houlberg, F90 free format 8/2004
!-------------------------------------------------------------------------------

!Declaration of input variables
REAL(KIND=rspec), INTENT(IN) :: &
  te                     !Maxwellian electron temperature [keV]

!Declaration of function output
REAL(KIND=rspec) :: &
  X_SIGV_RC              !<sigma*v> [m**3/s]

!-------------------------------------------------------------------------------
!Evaluate reaction rate
!-------------------------------------------------------------------------------
IF(te > 0.4) THEN

  X_SIGV_RC=3.58e-22/te**1.388

ELSE

  X_SIGV_RC=1.48e-20/(SQRT(te)*(1.0+43.4*te))

ENDIF

END FUNCTION X_SIGV_RC

END MODULE X_MOD 
