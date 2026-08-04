MODULE EQDSK_MOD
!-------------------------------------------------------------------------------
! EQDSK_MOD is a module of routines that read geometric and equilibrium
!   details from an eqdsk/geqdsk file.
!
! Contains PUBLIC routines:
!
!   READ_EFIT_EQDSK     -reads eqdsk file and returns constants, fields defined
!                        on the psi(x,y) grid.
!                        *moved from pellet_dr.f90*
!   PELLET_EFIT         -generates plasma geometry information for PELLET by 
!                        reading an eqdsk/geqdsk file
!                        *moved from pellet_dr.f90*
!
!
!-------------------------------------------------------------------------------

USE SPEC_KIND_MOD
USE FLUXAV_MOD
IMPLICIT NONE

CONTAINS

SUBROUTINE READ_EFIT_EQDSK(nin,cnin,mxnx_xy,mxny_xy,mxn_lim, &
                           bt0,cur,psimag,psilim,r0,rmag,zmag, &
                           nx_xy,ny_xy,x_xy,y_xy,psi_xy, &
                           f_x,ffp_x,psi_x,q_x,rhop_x, &
                           n_lim,x_lim,y_lim, &
                           iflag,message)
!-------------------------------------------------------------------------------
!READ_EFIT_EQDSK reads an EQDSK file from EFIT
!
!References:
!  W.A.Houlberg, F90 free format 8/2004
!
!Comments:
!  This routine reads and discards some of the stored EIFT data (e.g., the
!    boundary points that are not presently used with this application).
!  It also constructs the relevant 1-D and 2-D grids that are implicit in the
!    stored data.
!-------------------------------------------------------------------------------
USE SPEC_KIND_MOD
IMPLICIT NONE

!Declaration of input variables
CHARACTER(len=*), INTENT(IN) :: &
  cnin                   !input file name [character]

INTEGER,INTENT(IN) :: &
  mxnx_xy,             & !maximum number of x points on psi(x,y) grid [-]
  mxny_xy,             & !maximum number of y points on psi(x,y) grid [-]
  mxn_lim,             & !maximum number of points on limiter surface [-]
  nin                    !input unit number [-]

!Declaration of output variables
CHARACTER(len=*), INTENT(OUT) :: &
  message                !warning or error message [character]

INTEGER, INTENT(OUT) :: &
  iflag,               & !error and warning flag [-]
                         !=-1 warning
                         !=0 none
                         !=1 error
  nx_xy,               & !number of x points on psi(x,y) grid [-]
  ny_xy,               & !number of y points on psi(x,y) grid [-]
  n_lim                  !number of points on limiter [-]

REAL(KIND=rspec), INTENT(OUT) :: &
  bt0,                 & !toroidal field at r0 [T]
  cur,                 & !toroidal plasma current [A]
  psimag,              & !poloidal flux/(2*pi) at axis [Wb/rad]
  psilim,              & !poloidal flux/(2*pi) at limiter/separatrix [Wb/rad]
  r0,                  & !reference major radius, center of limiter surface [m]
  rmag,                & !horizontal position of magnetic axis [m]
  zmag                   !vertical position of magnetic axis [m]

REAL(KIND=rspec), INTENT(OUT) :: &
  x_xy(mxnx_xy),             & !vertical grid for 2-D poloidal flux [m]
  y_xy(mxny_xy),             & !horizontal grid for 2-D poloidal flux [m]
  psi_xy(mxnx_xy,mxny_xy),   & !poloidal flux/(2*pi) on 2-D grid [Wb/rad]
  f_x(mxnx_xy),              & !F=R*B_t on equilibrium psi grid [m*T]
  ffp_x(mxnx_xy),            & !F*dF/dpsi on equilibrium psi grid [rad*T]
  psi_x(mxnx_xy),            & !poloidal flux/(2*pi) = equilibrium psi grid [Wb/rad]
  q_x(mxnx_xy),              & !safety factor on equilibrium psi grid [-]
  rhop_x(mxnx_xy),           & !normalized poloidal flux grid proportional to psi [-]
  x_lim(mxn_lim),            & !horizontal positions of limiter points [m]
  y_lim(mxn_lim)               !vertical positions of limiter points [m]              

!-------------------------------------------------------------------------------
!Declaration of local variables
!Input from EQDSK file that is not retained
INTEGER :: &
  n_bdry                 !number of points on plasma boundary [-]


REAL(KIND=rspec) :: &
  rmin,                & !horizontal inside of computational domain [m]
  zmid,                & !vertical center of comoputational domain [m]
  rdim,                & !width of computational domain [m]
  zdim,                & !height of computational domain [m]
  x_bdry,              & !horizontal positions of boundary points [m]
  y_bdry,              & !vertical positions of boundary points [m]
  p_x(mxnx_xy),        & !plasma kinetic pressure [N/m**2]
  pp_x(mxnx_xy)          !dp/dpsi on equilibrium psi grid [rad*N/m**2/Wb]

!Other
INTEGER :: &
  i,j

REAL(KIND=rspec) :: &
  dum

CHARACTER(len=256) :: line

!-------------------------------------------------------------------------------
!Initialization
!-------------------------------------------------------------------------------
!Null output
iflag=0
message=''

! PRINT *, "Reading EQDSK."

!Open the EQDSK file
OPEN(UNIT=nin, &
     STATUS='old', & 
     FILE=cnin, &
     FORM='formatted')

!-------------------------------------------------------------------------------
!Read the EQDSK file
!-------------------------------------------------------------------------------
!Point data - dum values are duplicate information or not used
READ(nin,'(52x,2i4)') nx_xy,ny_xy

!Check if x dimension is exceeded
IF(nx_xy > mxnx_xy) THEN

  !Horizontal grid points exceed dimensions set by parameters
  iflag=1
  message='READ_EFIT_EQDSK(1)/ERROR:x grid dimension exceeded'
  GOTO 9999

ENDIF

!Check if y dimension is exceeded
IF(ny_xy > mxny_xy) THEN

  !Vertical grid points exceed dimensions set by parameters
  iflag=1
  message='READ_EFIT_EQDSK(2)/ERROR:y grid dimension exceeded'
  GOTO 9999

ENDIF

READ(nin,'(5e16.9)') rdim,zdim,r0,rmin,zmid
READ(nin,'(5e16.9)') rmag,zmag,psimag,psilim,bt0
READ(nin,'(5e16.9)') cur
READ(nin,'(5e16.9)') dum


!Read 1-D and 2-D data, radial grid is equally spaced in poloidal flux (1:nx_xy)
READ(nin,'(5e16.9)') (f_x(i),i=1,nx_xy)
READ(nin,'(5e16.9)') (p_x(i),i=1,nx_xy)
READ(nin,'(5e16.9)') (ffp_x(i),i=1,nx_xy)
READ(nin,'(5e16.9)') (pp_x(i),i=1,nx_xy)
READ(nin,'(5e16.9)') ((psi_xy(i,j),i=1,nx_xy),j=1,ny_xy)
READ(nin,'(5e16.9)') (q_x(i),i=1,nx_xy)

!Boundary and limiter data
READ(nin,'(2i5)') n_bdry,n_lim

!Check if boundary dimension is exceeded
!IF(n_bdry > mxn_bdry) THEN
!
!  !Boundary points exceed dimensions set by parameters
!  iflag=1
!  message='READ_EFIT_EQDSK(3)/ERROR:bdry grid dim exceeded'
!  GOTO 9999
!
!ENDIF

!Check if limiter dimension is exceeded
IF(n_lim > mxn_lim) THEN

  !Limiter points exceed dimensions set by parameters
  iflag=1
  message='READ_EFIT_EQDSK(4)/ERROR:lim grid dim exceeded'
  GOTO 9999

ENDIF

READ(nin,'(5e16.9)') (x_bdry,y_bdry,i=1,n_bdry)
READ(nin,'(5e16.9)') (x_lim(i),y_lim(i),i=1,n_lim)

!Construct implied grids
!2D grid
x_xy(1:nx_xy)=rmin+rdim*(/ (i-1,i=1,nx_xy) /)/(nx_xy-1)
y_xy(1:ny_xy)=zmid-zdim/2+zdim*(/ (i-1,i=1,ny_xy) /)/(ny_xy-1)

!1D radial grid and poloidal flux
psi_x(1:nx_xy)=psimag+(psilim-psimag)*(/ (i-1,i=1,nx_xy) /)/(nx_xy-1)
rhop_x(1:nx_xy)=(psi_x(1:nx_xy)-psi_x(1))/(psi_x(nx_xy)-psi_x(1))

!-------------------------------------------------------------------------------
!Cleanup and exit
!-------------------------------------------------------------------------------
9999 CONTINUE

!Close the EQDSK file
CLOSE(unit=nin)

END SUBROUTINE READ_EFIT_EQDSK

SUBROUTINE PELLET_EFIT(nin,cnin,nr_r,rhot_r, &
                       r0,a0,bt0,s0,e0,e1,d1,q0,q1,iflag,message)
!-------------------------------------------------------------------------------
!PELLET_EFIT generates plasma geometry information for PELLET by reading an
!  EQDSK file, calling FLUXAV_LOAD to load the data in the FLUXAV module,
!  calling FLUXAV to generate needed flux surface quantities.
!
!References:
!  W.A.Houlberg, F90 free format 8/2004
!-------------------------------------------------------------------------------

CHARACTER(len=*), INTENT(IN) :: &
  cnin                   !name of unit for data [character]

INTEGER, INTENT(IN) :: &
  nin,                 & !unit number for data [-]
  nr_r                   !number of radial points [-]

REAL(KIND=rspec), INTENT(IN) :: &
  rhot_r(nr_r)           !normalized tor flux grid proportional to (Phi)**0.5 [-]

!Declaration of output variables
INTEGER, INTENT(OUT) :: &
  iflag                  !error and warning flag [-]
                         !=-1 warning
                         !=0 none
                         !=1 error

CHARACTER(len=*), INTENT(OUT) :: &
  message                !warning or error message [character]

REAL(KIND=rspec), INTENT(OUT) :: &
  a0,                  & !minor radius, half diameter of boundary flux surface [m]
  bt0,                 & !toroidal field at r0 [T]
  d1,                  & !edge triangularity normalized to a0 [-]
  e0,                  & !axis elongation normalized to a0 [-]
  e1,                  & !edge elongation normalized to a0 [-]
  q0,                  & !axial safety factor [-]
  q1,                  & !edge sagety factor [-]
  r0,                  & !major radius, center of boundary flux suface [m]
  s0                     !axis shift normalized to a0 [-]
  

!-------------------------------------------------------------------------------
!Declaration of local variables
INTEGER, PARAMETER :: &
  mxnx_xy=300, &
  mxny_xy=300,  &
  mxnr_r=300, &
  mxn_lim=500, &
  mxn_bdry=300

REAL(KIND=rspec), PARAMETER :: &
  z_mu0=1.2566e-06, &
  z_pi=3.141592654

INTEGER :: &
  k_grid

INTEGER :: &
  n_lim,nx_xy,ny_xy

REAL(KIND=rspec) :: &
  cur,psimag,psilim,rmag,zmag

REAL(KIND=rspec) :: &
  b2_r(nr_r),          & !<B**2> [T**2]
  bm2_r(nr_r),         & !<1/B**2> [/T**2]
  bpout_r(nr_r),       & !poloidal field at rout_r(i) [T]
  btout_r(nr_r),       & !toroidal field at rout_r(i) [T]
  elong_r(nr_r),       & !elongation [-]
  triang_r(nr_r),      & !triangularity [-]
  f_r(nr_r),           & !2*pi*R*B_t/mu0 [A]
  fhat_r(nr_r),        & !
  fm_r(3,nr_r),        & !geometric factor [-]
  ftrap_r(nr_r),       & !trapped fraction [-]
  gph_r(nr_r),         & !poloidal flux metric  [-]
  grho1_r(nr_r),       & !a0*<|grad(rhot_r)|> [-]
  grho2_r(nr_r),       & !a0**2*<|grad(rhot_r)|**2> [-]
  gr2bm2_r(nr_r),      & !a0**2*<|grad(rhot_r)|**2/B**2> [/T**2]
  grth_r(nr_r),        & !<n.grad(theta)> [/m]
  gth_r(nr_r),         & !toroidal flux metric  [-]
  phit_r(nr_r),        & !toroidal flux [Wb]
  psi_r(nr_r),         & !poloidal flux [Wb/rad]
  q_r(nr_r),           & !safety factor [-]
  r2_r(1:mxnr_r),      & !<R**2> [m**2]
  rm2_r(1:mxnr_r),     & !<1/R**2> [/m**2]
  rhop_r(nr_r),        & !normalized poloidal flux grid proportional to psi [-]
  rin_r(nr_r),         & !major radius grid on inside of torus in axis plane [m]
  rout_r(nr_r),        & !major radius grid on outside of torus in axis plane [m]
  vol_r(nr_r),         & !volume enclosed [m**3]
  vp_r(nr_r)             !d vol_r/d rhot_r/a0 [m**2]

REAL(KIND=rspec) :: &
  f_x(mxnx_xy),ffp_x(mxnx_xy),psi_x(mxnx_xy),q_x(mxnx_xy),rhop_x(mxnx_xy), &
  x_xy(1:mxnx_xy),y_xy(1:mxny_xy),psi_xy(1:mxnx_xy,1:mxny_xy), &
  x_lim(1:mxn_lim),y_lim(1:mxn_lim)

! PRINT *, "Entering PELLET_EFIT."

!-------------------------------------------------------------------------------
!Initialization
!-------------------------------------------------------------------------------
!Local arrays
x_xy(:)=0
y_xy(:)=0
psi_xy(:,:)=0
f_x(:)=0
ffp_x(:)=0
rhop_x(:)=0
psi_x(:)=0
q_x(:)=0
rhop_x(:)=0
x_lim(:)=0
y_lim(:)=0

r2_r(:)=0
rm2_r(:)=0

!-------------------------------------------------------------------------------
!Get EFIT data
!-------------------------------------------------------------------------------
CALL READ_EFIT_EQDSK(nin,cnin,mxnx_xy,mxny_xy,mxn_lim, &
                     bt0,cur,psimag,psilim,r0,rmag,zmag, &
                     nx_xy,ny_xy,x_xy,y_xy,psi_xy, &
                     f_x,ffp_x,psi_x,q_x,rhop_x, &
                     n_lim,x_lim,y_lim, &
                     iflag,message)


! PRINT *, "Finished loading EQDSK."
!Check messages
IF(iflag /= 0) THEN

  message='FORCEBAL_EFIT(1)/'//message
  IF(iflag > 0) GOTO 9999

ENDIF

! r0 = a0*2.31907
! bt0 = bt0*(1.0/(r0 - 0.4)**2)


!-------------------------------------------------------------------------------
!Call FLUXAV to generate metrics from EFIT MHD equilibrium
!-------------------------------------------------------------------------------
CALL FLUXAV_LOAD(cur,r0,rmag,zmag,psimag,psilim, &
                 nx_xy,ny_xy,x_xy,y_xy,psi_xy,f_x,ffp_x,rhop_x,q_x, &
                 n_lim,x_lim,y_lim, &
                 iflag,message)

!Check messages
IF(iflag /= 0) THEN

  message='FORCEBAL_EFIT(2)/'//message
  IF(iflag > 0) GOTO 9999

ENDIF

k_grid=0

CALL FLUXAV(k_grid,nr_r,rhot_r, &
            a0,b2_r,bm2_r,bpout_r,btout_r,elong_r,triang_r,f_r,fhat_r,fm_r, &
            ftrap_r,gph_r,gr2bm2_r,grho1_r,grho2_r,grth_r,gth_r,phit_r, &
            psi_r,q_r,r2_r,rin_r,rm2_r,rout_r,vol_r,vp_r, &
            iflag,message)

!Check messages
IF(iflag /= 0) THEN

  message='FORCEBAL_EFIT(3)/'//message
  IF(iflag > 0) GOTO 9999

ENDIF

! PRINT *, "q_r: ", q_r

!Set 0-D quantities
r0=(rout_r(nr_r)+rin_r(nr_r))/2
s0=(rout_r(1)-r0)/a0
e0=elong_r(1)
e1=elong_r(nr_r)
d1=triang_r(nr_r)
q0=q_r(1)
q1=q_r(nr_r)

!Change grid normalization to a0
fhat_r(:)=fhat_r(:)*a0
gph_r(:)=gph_r(:)/a0
gr2bm2_r(:)=gr2bm2_r(:)*a0**2
grho1_r(:)=grho1_r(:)*a0
grho2_r(:)=grho2_r(:)*a0**2
vp_r(:)=vp_r(:)/a0

!Define rhop_r
rhop_r(1:nr_r)=(psi_r(1:nr_r)-psi_r(1))/(psi_r(nr_r)-psi_r(1))

!At this point q_r is always positive, correct sign for coordinate consistency
q_r(1:nr_r)=q_r(1:nr_r)*SIGN(1.0_rspec,bpout_r(nr_r)*btout_r(nr_r))
!-------------------------------------------------------------------------------
!Cleanup and exit
!-------------------------------------------------------------------------------
9999 CONTINUE

END SUBROUTINE PELLET_EFIT

END MODULE EQDSK_MOD