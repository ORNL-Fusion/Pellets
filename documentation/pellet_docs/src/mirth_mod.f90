MODULE MIRTH_MOD
!-------------------------------------------------------------------------------
!MIRTH_MOD is an F90 module that contains procedures for fast ions
!
!References:
!
!  W.A.Houlberg, F90 free format 8/2004
!
!Contains PUBLIC routines:
!
!  MIRTH_SS            -solves for steady-state fast ion distribution  
!-------------------------------------------------------------------------------
USE SPEC_KIND_MOD
USE X_MOD
IMPLICIT NONE

!-------------------------------------------------------------------------------
! Procedures
!-------------------------------------------------------------------------------
CONTAINS

SUBROUTINE MIRTH_SS(af,izf,e0,p,nr,te,den,zbr,dn0,dvol,h,ng,eg, &
                    dg, &
                    vc)
!-------------------------------------------------------------------------------
!MIRTH_SS populates the density profiles in energy groups for fast ions from a
!  monoenergetic source (neutral beam injection or fusion products) using a
!  steady-state solution of the Fokker-Planck equation
!
!References:
!  S.E. Attenberger, W.A. Houlberg, Nucl Technol/Fusion 4 (1983) 129
!  W.A.Houlberg, F90 free format 8/2004
!-------------------------------------------------------------------------------

!Declaration of input variables
INTEGER, INTENT(IN) :: &
  izf,                 & !charge number of fast ions [-]
  nr,                  & !number of radial cells in plasma [-]
  ng                     !number of fast ion energy groups [-]

REAL(KIND=rspec), INTENT(IN) :: &
  den(:),              & !electron density [/m**3]
  dn0(:),              & !neutral density [/m**3]
  dvol(:),             & !volume of plasma cells [m**3]
  eg(:),               & !fast ion energy group boundaries, eg(j+1)<eg(j) [keV]
  h(:),                & !source profile shape for fast ions [-]
  te(:),               & !electron temperature [keV]
  zbr(:)                 !sum(Z_j**2*n_j/A_j)/ne in v_crit calculation [-]

REAL(KIND=rspec), INTENT(IN) :: &
  af,                  & !atomic mass number of fast ions [-]
  e0,                  & !initial energy of fast ions [keV]
  p                      !power in fast ion source [w]

!Declaration of input/output variables
REAL(KIND=rspec), INTENT(INOUT) :: &
  dg(:,:)                !fast ion density [/m**3]

!Declaration of output variables
REAL(KIND=rspec), INTENT(OUT) :: &
  vc(:)                  !critical velocity of fast ions [m/s]

!-------------------------------------------------------------------------------
!Declaration of local variables
INTEGER :: &
  i,j

REAL(KIND=rspec) :: &
  con1,con2,deldg,e0j,egj,sav,sgj,tcxm,ts,tslow,vc3,ve,vlo,vup,xln

REAL(KIND=rspec), PARAMETER :: &
  z_pi=3.141592654, &
  z_coulomb=1.6022e-19, &
  z_electronmass=9.1095e-31, &
  z_epsilon0=8.8542e-12, &
  z_protonmass=1.6726e-27, &
  z_j7kv=1.6022e-16

!-------------------------------------------------------------------------------
!Initialization
!-------------------------------------------------------------------------------
!Constants
e0j=e0*z_j7kv
tcxm=0

!Constant for Ecrit calculation
con1=(3*SQRT(z_pi)/4*(z_electronmass/z_protonmass))**(1.0/3.0)

!Constant for characteristic momentum relaxation time on electrons
con2=3*(af*z_protonmass)/4/SQRT(2*z_pi*z_electronmass) &
     *(4*z_pi*z_epsilon0/z_coulomb/izf)**2

!-------------------------------------------------------------------------------
!Fast ion density
!-------------------------------------------------------------------------------
!Average source amplitude
sav=p/e0j/SUM(dvol(1:nr))

!Set group densities using steady state Fokker-Planck solution
DO i=1,nr !Over radial nodes

  xln=37.8-LOG(SQRT(den(i))/te(i))
  ve=SQRT(2*te(i)*z_j7kv/z_electronmass)
  vc(i)=con1*zbr(i)**(1.0/3.0)*ve
  ts=con2*(te(i)*z_j7kv)**1.5/z_coulomb/(z_coulomb*den(i))/xln
  vc3=vc(i)**3
  egj=eg(1)*z_j7kv
  IF(egj > e0j) egj=e0j
  vlo=SQRT(2*egj/(af*z_protonmass))
  sgj=sav*h(i)

  DO j=1,ng !Over energy groups

    egj=eg(j+1)*z_j7kv

    IF(egj < 0.99999*e0j) THEN

      vup=vlo
      vlo=SQRT(2.0_rspec*egj/(af*z_protonmass))
      tslow=ts/3*LOG((vup**3+vc3)/(vlo**3+vc3))

      IF(izf == 1) THEN

        tcxm=dn0(i)*(vup+vlo)/2*X_SIG_PX((eg(j)+eg(j+1))/2,af)

      ENDIF

      deldg=sgj/(tcxm+1.0_rspec/tslow)
      dg(i,j)=dg(i,j)+deldg
      sgj=deldg/tslow

    ENDIF

  ENDDO !Over energy groups

ENDDO !Over radial nodes

END SUBROUTINE MIRTH_SS

END MODULE MIRTH_MOD
