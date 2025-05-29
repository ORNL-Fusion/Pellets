SUBROUTINE SETUP_AJAX(k_equil,n_eq,r0,a0,s0,e0,e1,d1,bt0,q0,q1,n_rho,rho,      &
                      nrho_ajax,ntheta_ajax,nzeta_ajax, &
                      iflag,message)
!-------------------------------------------------------------------------------
!SETUP_AJAX sets up the equilibrium information for AJAX
!
!References:
!  W.A.Houlberg, F90 free format 8/2004
!-------------------------------------------------------------------------------
USE SPEC_KIND_MOD
USE AJAX_MOD
IMPLICIT NONE

!Declaration of input variables
INTEGER, INTENT(IN) :: &
  k_equil,             & !MHD equilibrium/geometry option [-]
                         !=1 read inverse coordinate data file
                         !=else use simple approximation to 2-D plasma
  n_eq,                & !input unit number [-]
  n_rho,               & !no. radial nodes for output [-]
  nrho_ajax,           & !radial nodes to use in internal AJAX data [-]
  ntheta_ajax,         & !poloidal nodes to use in internal AJAX data [-]
  nzeta_ajax             !toroidal nodes to use in 3-D internal AJAX data [-]

REAL(KIND=rspec), INTENT(IN) :: &
  r0,                  & !major radius of geometric center [m]
  a0,                  & !minor radius in midplane [m]
  s0,                  & !axis shift normalized to a0 [-]
  e0,                  & !axis elongation normalized to a0 [-]
  e1,                  & !edge elongation normalized to a0 [-]
  d1,                  & !edge triangularity normalized to a0 [-]
  bt0,                 & !toroidal field at r0 [T]
  q0,                  & !safety factor on axis [-]
  q1,                  & !safety factor at edge [-]
  rho(1:n_rho)           !radial grid [rho]

!Declaration of output variables
CHARACTER(len=*), INTENT(OUT) :: &
  message                !warning or error message [character]

INTEGER, INTENT(OUT) :: &
  iflag                  !error and warning flag [-]
                         !=-1 warning
                         !=0 none
                         !=1 error

!-------------------------------------------------------------------------------
!Declaration of local variables
INTEGER :: &
  k_pflx

REAL(kind=rspec) :: &
  phitot,q(1:n_rho),v1(1:n_rho),v2(1:n_rho)

!Physical and conversion constants
REAL(KIND=rspec), PARAMETER :: &
  z_pi=3.141592654

!-------------------------------------------------------------------------------
!Get data
!-------------------------------------------------------------------------------
!Read inverse coordinate expansion
IF(k_equil == 1) THEN

  iflag=0
  message=''
  CALL READEQ(n_eq,nrho_ajax,ntheta_ajax,nzeta_ajax, &
              iflag,message)

  !Check messages
  IF(iflag > 0) THEN

    message='SETUP_AJAX(1)/'//message
    GOTO 9999

  ENDIF

ELSE

  !Load boundary values of R,Z and fill in with approximations
  iflag=0
  message=''
  CALL AJAX_LOAD_RZBDY(r0,a0,s0,e0,e1,d1, &
                       iflag,message, &
                       NRHO_AJAX=nrho_ajax, &
                       NTHETA_AJAX=ntheta_ajax)

  !Check messages
  IF(iflag > 0) THEN

    message='SETUP_AJAX(2)/'//message
    GOTO 9999

  ENDIF

  !Get metrics to evaluate total toroidal magnetic flux
  iflag=0
  message=''
  CALL AJAX_FLUXAV_G(n_rho,rho, &
                     iflag,message, &
                     RM2_R=v1, &
                     VP_R=v2)
  phitot=r0*bt0/4/z_pi*rho(n_rho)*v1(n_rho)*v2(n_rho)

  !Check messages
  IF(iflag > 0) THEN

    message='SETUP_AJAX(3)/'//message
    GOTO 9999

  ENDIF

  k_pflx=0
  q(1:n_rho)=q0+(q1-q0)*rho(1:n_rho)**2
  iflag=0
  message=''
  CALL AJAX_LOAD_MAGFLUX(phitot,k_pflx,n_rho,rho,q, &
                         iflag,message)

  !Check messages
  IF(iflag > 0) THEN

    message='SETUP_AJAX(4)/'//message
    GOTO 9999

  ENDIF

ENDIF

!-------------------------------------------------------------------------------
!Cleanup and exit
!-------------------------------------------------------------------------------
9999 CONTINUE

END SUBROUTINE SETUP_AJAX

SUBROUTINE READEQ(n_eq,nrho_ajax,ntheta_ajax,nzeta_ajax, &
                  iflag,message)
!-------------------------------------------------------------------------------
!READEQ reads a simplified form of a 3-D VMEC equilbrium solution
!
!References:
!  W.A.Houlberg, F90 free format 8/2004
!-------------------------------------------------------------------------------
USE SPEC_KIND_MOD
USE AJAX_MOD
IMPLICIT NONE

!Declaration of input variables
INTEGER, INTENT(IN) :: &
  n_eq,                & !input unit number [-]
  nrho_ajax,           & !radial nodes to use in internal AJAX data [-]
  ntheta_ajax,         & !poloidal nodes to use in internal AJAX data [-]
  nzeta_ajax             !toroidal nodes to use in 3-D internal AJAX data [-]

!Declaration of output variables
CHARACTER(len=*), INTENT(OUT) :: &
  message                !warning or error message [character]

INTEGER, INTENT(OUT) :: &
  iflag                  !error and warning flag [-]
                         !=-1 warning
                         !=0 none
                         !=1 error

!-------------------------------------------------------------------------------
!Declaration of local variables
LOGICAL :: &
  l_axisymmetric

INTEGER :: &
  i,k, &
  k_grid,k_pflx, &
  nr_iota,nr_rzl,nk_rzl

INTEGER, ALLOCATABLE :: &
  m(:),n(:)

REAL(KIND=rspec) :: &
  phitot

REAL(KIND=rspec), ALLOCATABLE :: &
  rho_iota(:),iotabar(:),rho_rz(:),rho_lam(:),rmn(:,:),zmn(:,:),lmn(:,:)

!-------------------------------------------------------------------------------
!Read data
!-------------------------------------------------------------------------------
!phitot [Wb]
READ(n_eq,*) phitot

!iotabar [-], rho~sqrt(toroidal flux) and normalized to unity at boundary
READ(n_eq,*) nr_iota
ALLOCATE(rho_iota(1:nr_iota), &
         iotabar(1:nr_iota))

  READ(n_eq,*) (rho_iota(i),i=1,nr_iota)
  READ(n_eq,*) (iotabar(i),i=1,nr_iota)

!R,Z and lambda grids, rho~sqrt(toroidal flux) and arbitrary normalization
k_grid=0
READ(n_eq,*) nr_rzl
ALLOCATE(rho_rz(1:nr_rzl), &
         rho_lam(1:nr_rzl))

  READ(n_eq,*) (rho_rz(i),i=1,nr_rzl)
  READ(n_eq,*) (rho_lam(i),i=1,nr_rzl)

!R,Z and lambda modes
READ(n_eq,*) nk_rzl
ALLOCATE(m(nk_rzl), &
         n(nk_rzl), &
         rmn(nr_rzl,nk_rzl), &
         zmn(nr_rzl,nk_rzl), &
         lmn(nr_rzl,nk_rzl))

  !Check to see whether this is a tokamak or stellarator
  l_axisymmetric=.TRUE.

  DO k=1,nk_rzl

    READ(n_eq,*) m(k),n(k)
    IF(n(k) /= 0 .AND. l_axisymmetric) l_axisymmetric=.FALSE.
    READ(n_eq,*) (rmn(i,k),i=1,nr_rzl)
    READ(n_eq,*) (zmn(i,k),i=1,nr_rzl)
    READ(n_eq,*) (lmn(i,k),i=1,nr_rzl)

  ENDDO

!-------------------------------------------------------------------------------
!Load data into AJAX
!-------------------------------------------------------------------------------
SELECT CASE (l_axisymmetric)

CASE (.TRUE.)

  !Axisymmetric plasma, calculate lambdas in AJAX on the internal grid
  iflag=0
  message=''
  CALL AJAX_LOAD_RZLAM(nr_rzl,nk_rzl,rho_rz,m,n,rmn,zmn, &
                       iflag,message, &
                       K_GRID=k_grid, &
                       L_MFILTER_AJAX=.TRUE., &
                       NRHO_AJAX=nrho_ajax, &
                       NTHETA_AJAX=ntheta_ajax)

  !Check messages
  IF(iflag > 0) THEN

    message='READEQ(1)/'//message
    GOTO 9999

  ENDIF

CASE (.FALSE.)

  !Non-axisymmetric plasma, use the lambdas from the data file
  iflag=0
  message=''
  CALL AJAX_LOAD_RZLAM(nr_rzl,nk_rzl,rho_rz,m,n,rmn,zmn, &
                       iflag,message, &
                       K_GRID=k_grid, &
                       NRHO_AJAX=nrho_ajax, &
                       NTHETA_AJAX=ntheta_ajax, &
                       NZETA_AJAX=nzeta_ajax, &
                       NR_LAM=nr_rzl, &
                       NK_LAM=nk_rzl, &
                       RHO_LAM=rho_lam, &
                       LAM=lmn)

  !Check messages
  IF(iflag > 0) THEN

    message='READEQ(2)/'//message
    GOTO 9999

  ENDIF

END SELECT

iflag=0
message=''
k_pflx=1
CALL AJAX_LOAD_MAGFLUX(phitot,k_pflx,nr_iota,rho_iota,iotabar, &
                       iflag,message)

!Check messages
IF(iflag > 0) message='READEQ(3)/'//message

!-------------------------------------------------------------------------------
!Cleanup and exit
!-------------------------------------------------------------------------------
9999 CONTINUE

DEALLOCATE(rho_iota, &
           iotabar)

DEALLOCATE(rho_rz, &
           rho_lam, &
           m, &
           n, &
           rmn, &
           zmn, &
           lmn)

END SUBROUTINE READEQ

