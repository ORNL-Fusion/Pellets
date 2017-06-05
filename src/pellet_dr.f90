PROGRAM PELLET_DR
!-------------------------------------------------------------------------------
!PELLET_DR is a program to run the PELLET module
!
!References:
!  W.A.Houlberg, L.R.Baylor, F90 free format 8/2004
!
!Comments:
!  Units are contained in [] in the description of variables
!-------------------------------------------------------------------------------
USE SPEC_KIND_MOD
USE AJAX_MOD
USE MIRTH_MOD
USE PELLET_MOD
USE TRACK_MOD
USE WRITE_MOD
IMPLICIT NONE

!-------------------------------------------------------------------------------
!Declaration of namelist input variables
CHARACTER(len=25) :: &
  cn_eq,               & !name of input file for MHD equilibrium data [character]
                         !required only for k_equil=1
  cn_prof,             & !name of input file for plasma profiles [character]
                         !required only for k_readd=1
  cn_runid,            & !Run id to use in output filename
                         !
  cn_device              !name of device for pellet experiment [character]
                         !i.e. JET, D3D, etc.


INTEGER :: &
  k_equil,             & !option for MHD equilibrium/geometry [-]
                         !=0    simple approximation to 2D plasma using AJAX
                         !      (specify r0,a0,s0,e0,e1,d1,bt0,cur,q0)
                         !=1    read 2D or 3D inverse coordinate data file using AJAX
                         !      (specify file name in cn_eq)
                         !=else not allowed
  k_pel,               & !option for ablation model [hydrogenic 0-5, other 10-11]
                         !=0    NGS:  reference model
                         !            electron distribution
                         !            elliptical neutral shield
                         !=1    NGS:  original Milora model
                         !            single electron energy
                         !            spherical neutral shield
                         !=2    NGPS: 1mm neutral layer thickness
                         !            electron distribution
                         !=3    Macaulay neutral gas shielding
                         !=4    Kuteev neutral gas shielding
                         !=5    Parks neutral gas shielding
                         !=6    Parks neutral gas shielding arbitrary Q
                         !=10   Parks impurity model
                         !=11   Kuteev impurity model
                         !=12   Parks new impurity model 2012
                         !=else failure
  k_readd,             & !option to read profile data set [-]
                         !=1    read an input text file for profiles
                         !=else do not read data file, use input profile factors.
  k_seg_p,             & !option for coordinates of the pellet path segment [-]
                         !=1    flux coordinates (rho,theta,zeta)
                         !=2    Cartesian coordinates (x,y,z)
                         !=else default cylindrical coordinates (R,phi,Z)
  ncplas,              & !number of radial cells in plasma [-]
  ncsol,               & !number of radial cells in scrape-off layer [-]
  nhorz,               & !number of horizontal chords for line densities [0-5]
  nvert,               & !number of vertical chords for line densities [0-5]
  nprlcld,             & !number of cloudlets to be used in PRL modeling
  iprlcld,             & !id number of cloudlet to output PRL diagnostic data
  k_prl,               & !switch to use PRL deposition model, 1=on
  k_ped                  !switch for pedestal Te model, 1= on


REAL(KIND=rspec) :: &
  amu_pel,             & !atomic mass number of pellet atoms [-]
  r_pel,               & !equivalent spherical radius of pellet [m]
  v_pel,               & !pellet velocity [m/s]
  rseg_p(3,2),         & !coordinates defining the segment [see k_seg_p]
                         !(:,1) is launch point
                         !(:,2) is target point
                         !required only for k_equil = 0
  r0,                  & !major radius of geometric center [m]
  a0,                  & !minor radius in midplane [m]
  dsol,                & !SOL thickness in midplane normalized to a0 [-]
  s0,                  & !axis shift normalized to a0 [-]
  e0,                  & !axis elongation normalized to a0 [-]
  e1,                  & !edge elongation normalized to a0 [-]
  d1,                  & !edge triangularity normalized to a0 [-]
  bt0,                 & !toroidal field at R0 [T]
  cur,                 & !toroidal plasma current [MA]
  q0,                  & !axial safety factor [-]
                         !required only for k_readd /= 1
                         !te(x)=te1+(te0-te1)*(1-x**px_te)**qx_te
  te0,                 & !electron temperature on axis [kev]
  te1,                 & !electron temperature at separatrix/limiter [kev]
  px_te,               & !electron temperature profile exponent [-]
  qx_te,               & !electron temperature profile exponent [-]
  rl_te,               & !electron temperature exponential decay length in SOL [m]
                         !required only for k_readd /= 1
                         !den(x)=den1+(den0-den1)*(1-x**px_den)**qx_den
  den0,                & !electron density on axis [/m**3]
  den1,                & !electron density at separatrix/limiter [/m**3]
  px_den,              & !electron density profile exponent [-]
  qx_den,              & !electron density profile exponent [-]
  rl_den,              & !electron density exponential decay length in SOL [m]
  rvert(5),            & !major radii of vertical chords [m]
  zhorz(5),            & !vertical positions of normal horizontal chords [m]
  pb(3),               & !power in full, half, and third beam components [w]
                         !hb(x)=(1-x**px_hb(j))**qx_hb(j))
  px_hb(3),            & !fast beam ion hj(r) profile exponent
  qx_hb(3),            & !fast beam ion hj(r) profile exponent
  amu_b,               & !atomic mass number of beam ions [1-3]
  eb0,                 & !neutral beam injection energy [keV]
  amu_i,               & !atomic mass number of plasma ions [1-3]
  dn01,                & !neutral density at separatrix/limiter [/m**3]
  rl_dn0,              & !neutral density exp falloff in plasma [m]
  pa,                  & !power in fast alphas [w]
                         !ha_r(x)=(1-x**px_ha)**qx_ha
  px_ha,               & !fast alpha h(r) profile exponent
  qx_ha,               & !fast alpha h(r) profile exponent
! PRL specific namelist variables
  pedte,               & !pedestal Te value
  pedne,               & !pedestal ne value
  pedwid,              & !pedestal Te width
  fpelprl,             & !fraction of pellet mass unaffected by PRL drift model 
  prlinjang,           & !PRL injection angle
  prlq0,               & !PRL q0 value
  prlqa,               & !PRL qa value 
  prlqf                  !PRL q profile exponent q(r) = (qa-q0) + q0*(1-(r/a)^qf)
!
!
!


!-------------------------------------------------------------------------------
!Declaration of other variables
LOGICAL :: &
  l_fa,l_fb(3)

CHARACTER(len=45) :: &
  cn_msg,cn_sum,cn_tmp

CHARACTER(len=120) :: &
  label,message

CHARACTER(len=1) :: &
  cb

INTEGER :: &
  n_msg,n_sum,n_tmp

INTEGER :: &
  i,i0,idum(5),ii,j,n,n_c,n_p,nc,nefmax,nf,iflag, &
  nrho_ajax,ntheta_ajax,nzeta_ajax

INTEGER, ALLOCATABLE :: &
  irho_c(:),izone_c(:), &
  irho_p(:),izone_p(:), &
  iz_f(:),ne_f(:)

INTEGER, PARAMETER :: &
  nbg=17, &
  nag=20

REAL(KIND=rspec) :: &
  drp,drso,frpen,fvpen,hnorm,pel_ions,ptemp,q1,raxis,rpen,rsoff,rtot, &
  rdum(5),r_flx(3),r_cyl(3),rseg_c(3,2),volp,volpen,dndlv(4,5),dndlh(4,5), &
  timndl(4)

REAL(KIND=rspec), ALLOCATABLE :: &
  cur_rm(:),den_r(:),dvol_r(:),pden_r(:),rho_r(:),rho_rm(:),te_r(:), &
  dn0_r(:),zbr_r(:),amu_f(:),e0_f(:),p_f(:),e_ef(:,:),h_rf(:,:), &
  vc_rf(:,:),den_ref(:,:,:), &
  t_p(:),s_p(:),te0_p(:),te1_p(:),den0_p(:),den1_p(:),rpel1_p(:), &
  src_p(:),rcyl_p(:,:),rflx_p(:,:), &
  s_c(:), pedte1(:), pedte2(:), pedne1(:), pedne2(:), prlq_r(:), prldep(:)

!Physical constants, mathematical constants, conversion factors
REAL(KIND=rspec), PARAMETER :: &
  z_eion=32.6e-3, &
  z_pi=3.141592654, &
  z_mu0=4.0e-7*z_pi  

!PRL added variables
REAL(KIND=rspec) :: &
   pedlambda, xped, xmid, thill

!Output arrays
INTEGER :: &
  npro

INTEGER, PARAMETER :: &
  mx_npro=100

CHARACTER(len=15) :: &
  namepro(mx_npro),unitpro(mx_npro)

CHARACTER(len=95) :: &
  descpro(mx_npro)

REAL(KIND=rspec), ALLOCATABLE :: &
  valpro(:,:)

!-------------------------------------------------------------------------------
NAMELIST/indata/cn_eq,cn_prof, cn_runid, cn_device, &
                k_equil,k_pel,k_readd,k_seg_p, &
                ncplas,ncsol,nhorz,nvert, &
                amu_pel,r_pel,v_pel, &
                rseg_p, &
                r0,a0,dsol,s0,e0,e1,d1,bt0,cur,q0, &
                te0,te1,px_te,qx_te,rl_te, &
                den0,den1,px_den,qx_den,rl_den, &
                rvert,zhorz, &
                pb,px_hb,qx_hb,amu_b,eb0,amu_i,dn01,rl_dn0, &
                pa,px_ha,qx_ha, k_prl, nprlcld, iprlcld, &
				pedte, pedne, pedwid, k_ped, fpelprl, prlinjang, prlq0, prlqa, prlqf
               

!-------------------------------------------------------------------------------
!Initialization
!-------------------------------------------------------------------------------
!Warning/error flag and message
iflag=0
message=''

!AJAX grids
nrho_ajax=31
ntheta_ajax=33
nzeta_ajax=21

!Namelist input - character
cn_eq=''
cn_prof=''
cn_device=''

!Namelist input - integer
k_equil=0
k_pel=0
k_readd=0
k_seg_p=0
k_prl=0     ! PRL model flag  LRB
nprlcld=0
ncplas=0
ncsol=0
nvert=0
nhorz=0
iprlcld=0   ! PRL diag flag - cloudlet number to dump data  LRB

!Namelist input - real
amu_pel=0
r_pel=0
v_pel=0
rseg_p(:,:)=0
r0=0
a0=0
dsol=0
s0=0
e0=0
e1=0
d1=0
bt0=0
cur=0
q0=0
te0=0
te1=0
px_te=0
qx_te=0
rl_te=0
den0=0
den1=0
px_den=0
qx_den=0
rl_den=0
rvert(:)=0
zhorz(:)=0
pb(:)=0
px_hb(:)=0
qx_hb(:)=0
amu_b=0
eb0=0
amu_i=0
dn01=0
rl_dn0=0
pa=0
px_ha=0
qx_ha=0

pedte=0
pedne=0
pedwid=0
fpelprl=0
prlinjang=0
prlq0=0
prlqa=0
prlqf=0


!-------------------------------------------------------------------------------
!Set the input namelist unit, open, read and close file
!-------------------------------------------------------------------------------
n_tmp=20
cn_tmp='nml_pellet.dat'
OPEN(UNIT=n_tmp, &
     FILE=cn_tmp, &
     STATUS='old', &
     ACCESS='sequential')

READ(n_tmp,indata)

CLOSE(UNIT=n_tmp)

!-------------------------------------------------------------------------------
!Open output files
!-------------------------------------------------------------------------------
!Summary
n_sum=10
cn_sum='sum_pellet'//cn_runid//'.dat'
OPEN(UNIT=n_sum, &
     FILE=cn_sum, &
     STATUS='unknown', &
     FORM='formatted')

!Messages
n_msg=11
cn_msg='msg_pellet'//cn_runid//'.dat'
OPEN(UNIT=n_msg, &
     FILE=cn_msg, &
     STATUS='unknown', &
     FORM='formatted')

!-------------------------------------------------------------------------------
!Set radial grid
!-------------------------------------------------------------------------------
!This section sets up two grids - main grid (rho_r)
! - interface/mid grid (rho_rm)
!The interface grid lies half-way between the main grid points
!Densities and temperatures are on the main grid and represent the values
!  in the cells that are bounded by the interface/mid grid points
!There are both main and interface/mid nodes at the axis
!The last interface/mid node is the plasma limiter/separatrix or wall
!The last main grid node is for ghost values of density, temperature, etc
!Metric quantities are usually calculated on the interface/mid grid
!Convention for number of nodes and cells:
!  ncplas-number of radial cells in the plasma
!  ncsol-number of radial cells on the SOL
!  nc-total number of radial cells (=ncplas+ncsol)
!  n-total number of grid points (=nc+1)
!-------------------
!Schematic w/o SOL:
!
!  Main grid        |     |     | ... |     |
!    Label          1     2     3    n-1    n
!
!  Interface grid   |  |     | ... |     |
!    Label          1  2     3    n-1    n
!
!                   ^                    ^
!                 Axis                Sep/lim
!                rho_r=0              rho_r=1
!-------------------
!Schematic w/ SOL:
!
!  Main grid        |     |     |   ... | ...  |     |
!    Label          1     2     3   ncplas+1   n-1   n
!
!  Interface grid   |  |     | ... | ...   ||
!    Label          1  2     3   ncplas+1  n-1     n
!
!                   ^              ^               ^ 
!                 Axis          Sep/lim         Wall/edge
!               rho_r=0         rho_r=1          rho_r>1
!-------------------------------------------------------------------------------
!Number of radial cells
IF(ncsol < 0) ncsol=0
IF(dsol <= 0.0) ncsol=0
IF(ncsol == 0) dsol=0
nc=ncplas+ncsol
n=nc+1

!Allocate arrays with radial dependence
ALLOCATE(cur_rm(n), &
         den_r(n), &
         dvol_r(n), &
         te_r(n), &
         pden_r(n), &
         rho_r(n), &
         rho_rm(n), &
		 pedte1(n),  &
		 pedte2(n), &
		 pedne1(n),  &
		 pedne2(n), &
		 prlq_r(n), &
		 prldep(n))

  cur_rm(:)=0
  den_r(:)=0
  dvol_r(:)=0
  te_r(:)=0
  pden_r(:)=0
  rho_r(:)=0
  rho_rm(:)=0
  pedte1(:)=0
  pedte2(:)=0
  pedne1(:)=0
  pedne2(:)=0
  prlq_r(:)=0
  prldep(:)=0

!Set cell grid sizes and first node outside core (ghost or SOL)
IF(ncsol == 0) THEN

  !Core only
  drp=1/(REAL(ncplas,rspec)-0.5)
  rho_r(ncplas+1)=1.0+drp/2

ELSE

  !Core + SOL
  drso=dsol/ncsol
  drp=(1.0-drso/2)/(ncplas-1)
  rho_r(ncplas+1)=1.0+drso/2

ENDIF

!Main grid
!Core
rho_r(1:ncplas)=(/ (i-1,i=1,ncplas) /)*drp

!SOL
IF(ncsol /= 0) THEN

  rho_r(ncplas+2:n)=rho_r(ncplas+1)+(/ (i,i=1,ncsol) /)*drso

ENDIF

!Interface/mid grid
rho_rm(1)=0

DO i=2,n !Over radial nodes

  rho_rm(i)=(rho_r(i)+rho_r(i-1))/2

ENDDO !Over radial nodes

!-------------------------------------------------------------------------------
!Set up MHD equilibrium information
!-------------------------------------------------------------------------------
IF(k_equil == 0 .OR. &
   k_equil == 1) THEN

  !AJAX option for MHD equilibrium interface
  IF(k_equil == 1) THEN

    !Open data file for reading
    OPEN(UNIT=n_tmp, &
         FILE=cn_eq, &
         STATUS='old', &
         ACCESS='sequential')

  ENDIF

  q1=4*z_pi*a0**2*SQRT(e0)*bt0/(z_mu0*cur*r0)
  ptemp=2

  DO WHILE(ABS(ptemp-1.0) > 1.0e-3)

    CALL SETUP_AJAX(k_equil,n_tmp,r0,a0,s0,e0,e1,d1,bt0,q0,q1,ncplas+1,rho_rm, &
                    nrho_ajax,ntheta_ajax,nzeta_ajax, &
                    iflag,message)

    !Check messages
    IF(iflag /= 0) THEN

      CALL WRITE_LINE(n_msg,message,1,1)
      IF(iflag > 0) GOTO 9999
      iflag=0
      message=''

    ENDIF

    IF(k_equil == 1) THEN

      !Current determined by q profile in equilibrium file
      ptemp=1

    ELSE

      !Iterate on q at edge to match desired current
      CALL AJAX_I(ncplas+1,rho_rm, &
                  iflag,message, &
                  CUR_I_R=cur_rm)

      !Check messages
      IF(iflag /= 0) THEN

        CALL WRITE_LINE(n_msg,message,1,1)
        IF(iflag > 0) GOTO 9999
        iflag=0
        message=''

      ENDIF

      IF(ncsol > 0) cur_rm(ncplas+2:n)=cur_rm(ncplas+1)

      ptemp=cur_rm(n)/cur
      q1=q1*ptemp

    ENDIF

  ENDDO

  IF(k_equil == 1) CLOSE(UNIT=n_tmp)

ELSEIF(k_equil == 2) THEN

  !Open data file for reading
  OPEN(UNIT=n_tmp, &
       FILE=cn_eq, &
       STATUS='old', &
       ACCESS='sequential')

  CALL PELLET_EFIT(n_tmp,cn_eq,ncplas+1,rho_rm, &
                   r0,a0,bt0,s0,e0,e1,d1,q0,q1,iflag,message)

  CLOSE(UNIT=n_tmp)

  !Check messages
  IF(iflag /= 0) THEN

    CALL WRITE_LINE(n_msg,message,1,1)
    IF(iflag > 0) GOTO 9999
    iflag=0
    message=''

  ENDIF

  CALL SETUP_AJAX(k_equil,n_tmp,r0,a0,s0,e0,e1,d1,bt0,q0,q1,ncplas+1,rho_rm, &
                  nrho_ajax,ntheta_ajax,nzeta_ajax, &
                  iflag,message)

  !Check messages
  IF(iflag /= 0) THEN

    CALL WRITE_LINE(n_msg,message,1,1)
    IF(iflag > 0) GOTO 9999
    iflag=0
    message=''

  ENDIF

ELSE

  !Illegal choice of k_equil
  iflag=1
  message='PELLET_DR/ERROR:illegal k_equil (use 0,1)'
  CALL WRITE_LINE(n_msg,message,1,1)
  GOTO 9999

ENDIF

!Get geometric quantities
iflag=0
message=''
CALL AJAX_FLUXAV_G(n,rho_rm, &
                   iflag,message, &
                   DVOL_R=dvol_r)

!Check messages
IF(iflag /= 0) THEN

  CALL WRITE_LINE(n_msg,message,1,1)
  IF(iflag > 0) GOTO 9999

ENDIF

!-------------------------------------------------------------------------------
!Plasma profiles
!-------------------------------------------------------------------------------
IF(k_readd == 1) THEN

  !Get electron profiles from file
  CALL READ_PELLET_PRO(n_tmp,cn_prof,n,rho_r, &
                       den_r,te_r,iflag,message)

ELSE

  !Set electron profiles from profile parameters
  !Core
  te_r(1:ncplas)=te1+(te0-te1)*(1.0-rho_r(1:ncplas)**px_te)**qx_te
  den_r(1:ncplas)=den1+(den0-den1)*(1.0-rho_r(1:ncplas)**px_den)**qx_den


!
! Pedestal profile for Te and ne
!	
  if (k_ped .gt. 0) then  ! Implement pedestal profile using tanh model
	   xped = 1-pedwid
	   xmid = 1-pedwid/2
	   pedlambda = (pedte-te1)/(1+tanh(1.0))
	   thill = te0-pedte
       do i=1,ncplas
	      pedte1(i) = te1+pedlambda*  &
                      (tanh(1.0)-tanh(2*(rho_r(i)-xmid)/pedwid))
	      if (rho_r(i) .lt. (1.0-pedwid)) then
               pedte2(i) = thill*((1-(rho_r(i)/xped)**px_te)**qx_te)
	      else
	         pedte2(i) = 0.0
	      endif
          te_r(i)=(pedte1(i)+pedte2(i))
       enddo
	   pedlambda = (pedne-den1)/(1+tanh(1.0))
	   thill = den0-pedne
       do i=1,ncplas
	      pedne1(i) = den1+pedlambda*  &
                      (tanh(1.0)-tanh(2*(rho_r(i)-xmid)/pedwid))
	      if (rho_r(i) .lt. (1.0-pedwid)) then
               pedne2(i) = thill*((1-(rho_r(i)/xped)**px_den)**qx_den)
	      else
	         pedne2(i) = 0.0
	      endif
          den_r(i)=(pedne1(i)+pedne2(i))
       enddo
  endif

!
! q profile for PRL model
!
  if (k_prl .gt. 0) then  ! Implement q profile using polynomial model
       do i=1,ncplas
	      prlq_r(i) = prlq0 + (prlqa-prlq0)*rho_r(i)**prlqf
       enddo
  endif


  !SOL
  DO i=ncplas+1,n

    te_r(i)=te1*EXP(-(rho_r(i)-1.0)*a0/rl_te)
    den_r(i)=den1*EXP(-(rho_r(i)-1.0)*a0/rl_den)

  ENDDO

ENDIF

!------------------------------------------------------------------------------
!Fast ion species
!------------------------------------------------------------------------------
!Determine number of species and energy groups
nf=0
nefmax=0
l_fb(:)=.FALSE.
l_fa=.FALSE.

DO j=1,3

  IF(pb(j) > 1.0e3) THEN

    l_fb(j)=.TRUE.
    nf=nf+1
    nefmax=nbg

  ENDIF

ENDDO

IF(pa > 1.0e3) THEN

  l_fa=.TRUE.
  nf=nf+1
  IF(nag > nefmax) nefmax=nag

ENDIF

!Allocate arrays
IF(nf > 0) THEN

  ALLOCATE(dn0_r(n), &
           zbr_r(n), &
           p_f(nf), &
           ne_f(nf), &
           iz_f(nf), &
           amu_f(nf), &
           e0_f(nf), &
           e_ef(nefmax+1,nf), &
           h_rf(n,nf), &
           vc_rf(n,nf), &
           den_ref(n,nefmax,nf))

    dn0_r(:)=0
    zbr_r(:)=0
    p_f(:)=0
    ne_f(:)=0
    iz_f(:)=0
    amu_f(:)=0
    e0_f(:)=0
    e_ef(:,:)=0
    h_rf(:,:)=0
    vc_rf(:,:)=0
    den_ref(:,:,:)=0

  !Neutral density
  !Core
  dn0_r(1:ncplas)=dn01*EXP((rho_r(1:ncplas)-1.0)*a0/rl_dn0)

  !SOL
  dn0_r(ncplas+1:n)=dn01

  ![Z] approximated for critical velocity/energy
  zbr_r(1:n)=1/amu_i

  !Fast ion density
  i=0

  DO j=1,3

    IF(l_fb(j)) THEN

      !Neutral beam ions
      i=i+1
      p_f(i)=pb(j)
      ne_f(i)=nbg
      iz_f(i)=1
      amu_f(i)=amu_b
      e0_f(i)=eb0/j
      e_ef(1:nbg+1,i)=eb0*(1.0-REAL((/ (i-1,i=1,nbg+1) /),rspec) &
                      /REAL(nbg+1,rspec))

      !H(r) normalized so that integral H(r)dV=V over plasma volume
      h_rf(1:ncplas,i)=(1.0-rho_r(1:ncplas)**px_hb(j))**qx_hb(j)
      hnorm=SUM(dvol_r(1:ncplas))/SUM(dvol_r(1:ncplas)*h_rf(1:ncplas,i))
      h_rf(:,i)=hnorm*h_rf(:,i)

      !Steady-state slowing down distribution
      CALL MIRTH_SS(amu_f(i),iz_f(i),e0_f(i),p_f(i),n,te_r,den_r,zbr_r,dn0_r, &
                    dvol_r,h_rf(:,i),ne_f(i),e_ef(:,i), &
                    den_ref(:,:,i), &
                    vc_rf(:,i))


    ENDIF

  ENDDO

  IF(l_fa) THEN

    !Fusion alphas
    i=i+1
    p_f(i)=pa
    ne_f(i)=nag
    iz_f(i)=2
    amu_f(i)=4
    e0_f(i)=3.52e3
    e_ef(1:nag+1,i)=e0_f(i)*(1.0-REAL((/ (i-1,i=1,nag+1) /),rspec)/(nag+1))

    !H(r) normalized so that integral H(r)dV=V over plasma volume
    h_rf(1:ncplas,i)=(1.0-rho_r(1:ncplas)**px_ha)**qx_ha
    h_rf(1+ncplas:n,i)=0
    hnorm=SUM(dvol_r(1:ncplas))/SUM(dvol_r(1:ncplas)*h_rf(1:ncplas,i))
    h_rf(:,i)=hnorm*h_rf(:,i)

    !Steady-state slowing down distribution
    CALL MIRTH_SS(amu_f(i),iz_f(i),e0_f(i),p_f(i),n,te_r,den_r,zbr_r,dn0_r, &
                  dvol_r,h_rf(:,i),ne_f(i),e_ef(:,i), &
                  den_ref(:,:,i), &
                  vc_rf(:,i))

  ENDIF

ENDIF

!-------------------------------------------------------------------------------
!Allocate arrays along pellet path and call TRACK for pellet path
!-------------------------------------------------------------------------------
!Dimensioned to allow up to 6 times as many intersections as surfaces
ALLOCATE(irho_p(6*n), &
         izone_p(6*n), &
         s_p(6*n), &
         t_p(6*n), &
         te0_p(6*n), &
         te1_p(6*n), &
         den0_p(6*n), &
         den1_p(6*n), &
         src_p(6*n), &
         rpel1_p(6*n), &
         rcyl_p(3,6*n), &
         rflx_p(3,6*n))

  irho_p(:)=0
  izone_p(:)=0
  s_p(:)=0
  t_p(:)=0
  te0_p(:)=0
  te1_p(:)=0
  den0_p(:)=0
  den1_p(:)=0
  src_p(:)=0
  rpel1_p(:)=0
  rcyl_p(:,:)=0
  rflx_p(:,:)=0

!Get pellet path
CALL TRACK(n,rho_rm,2,rseg_p, &
           n_p,irho_p,s_p,iflag,message, &
           K_SEG=k_seg_p, &
           IZONE_INT=izone_p, &
           RCYL_INT=rcyl_p, &
           RFLX_INT=rflx_p)

!Check messages
IF(iflag /= 0) THEN

  CALL WRITE_LINE(n_msg,message,1,1)
  IF(iflag > 0) GOTO 9999
  iflag=0
  message=''

ENDIF

!-------------------------------------------------------------------------------
!Call PELLET
!-------------------------------------------------------------------------------
CALL PELLET(k_pel,amu_pel,r_pel,v_pel,nc,dvol_r,den_r,te_r,n_p-1,izone_p,s_p, &
            pden_r,iflag,message, &
            NF=nf, &
            NE_F=ne_f, &
            IZ_F=iz_f, &
            AMU_F=amu_f, &
            E_EF=e_ef, &
            VC_RF=vc_rf, &
            DEN_REF=den_ref, &
            PEL_IONS=pel_ions, &
            T_P=t_p, &
            RPEL1_P=rpel1_p, &
            SRC_P=src_p, &
            DEN0_P=den0_p, &
            DEN1_P=den1_p, &
            TE0_P=te0_p, &
            TE1_P=te1_p, &
			R0=r0, &
			A0=a0, &
			BT0=bt0,  &
			NCSOL=ncsol, &
			K_PRL=k_prl, &
			NPRLCLD=nprlcld, &
			IPRLCLD=iprlcld, &
			FPELPRL=fpelprl, &
		    PRLINJANG=prlinjang, &
			PRLQ_R=prlq_r,   &
			PRLDEP=prldep)

!Check messages
IF(iflag /= 0) THEN

  CALL WRITE_LINE(n_msg,message,1,1)
  IF(iflag > 0) GOTO 9999
  iflag=0
  message=''

ENDIF

!-------------------------------------------------------------------------------
!Fractional radius and volume penetrated
!-------------------------------------------------------------------------------
!Cylindrical coordinates of axis
r_flx(:)=0
r_cyl(:)=0
CALL AJAX_FLX2CYL(r_flx, &
                  r_cyl,iflag,message)

!Check messages
IF(iflag /= 0) THEN

  CALL WRITE_LINE(n_msg,message,1,1)
  IF(iflag > 0) GOTO 9999
  iflag=0
  message=''

ENDIF

raxis=r_cyl(1)

!Cylindrical coordinates of edge
r_flx(1)=1
CALL AJAX_FLX2CYL(r_flx, &
                  r_cyl,iflag,message)

!Check messages
IF(iflag /= 0) THEN

  CALL WRITE_LINE(n_msg,message,1,1)
  IF(iflag > 0) GOTO 9999
  iflag=0
  message=''

ENDIF

rsoff=r_cyl(1)
rtot=rsoff-raxis
volp=SUM(dvol_r(1:ncplas))
i0=0

DO i=1,ncplas !Over core zones

  IF(pden_r(i) < 1.0e10) i0=i

ENDDO !Over core zones

IF(i0 > 0) THEN

  !Partial penetration
  r_flx(1)=rho_rm(i0+1)
  CALL AJAX_FLX2CYL(r_flx, &
                    r_cyl,iflag,message)
  rpen=rsoff-r_cyl(1)
  volpen=volp-SUM(dvol_r(1:i0))

ELSE

  !Penetration to axis
  rpen=rtot
  volpen=volp

ENDIF

frpen=rpen/rtot
fvpen=volpen/volp

!------------------------------------------------------------------------------
!Allocate arrays for chordal densities and set arbitrary times
!------------------------------------------------------------------------------
!Dimensioned to allow up to 6 times as many intersections as surfaces
IF(nvert > 0 .OR. &
   nhorz > 0) THEN

  ALLOCATE(irho_c(6*n), &
           izone_c(6*n), &
           s_c(6*n))

    irho_c(:)=0
    izone_c(:)=0
    s_c(:)=0

  timndl(1)=-1.0
  timndl(2)=-0.05
  timndl(3)=0.05
  timndl(4)=1.0

ENDIF

!------------------------------------------------------------------------------
!Line densities for vertical chords
!------------------------------------------------------------------------------
dndlv(:,:)=0

IF(nvert > 0) THEN

  DO j=1,nvert !Over vertical chords

    !Set endpoints below and above plasma
    rseg_c(1,1:2)=rvert(j)
    rseg_c(2,1:2)=0
    rseg_c(3,1)=-1.1*e1*a0*rho_rm(n)
    rseg_c(3,2)=1.1*e1*a0*rho_rm(n)

    !Get intersections with plasma
    n_c=0
    CALL TRACK(n,rho_rm,2,rseg_c, &
               n_c,irho_c,s_c,iflag,message, &
               IZONE_INT=izone_c)

    !Check messages
    IF(iflag /= 0) THEN

      CALL WRITE_LINE(n_msg,message,1,1)
      IF(iflag > 0) GOTO 9999
      iflag=0
      message=''

    ENDIF

    IF(n_c > 1) THEN

      !Integrate vertical line densities
      DO ii=1,n_c-1 !Over chord segments

        i=izone_c(ii)

        IF(i > 0 .AND. &
           i <= nc) THEN

          !Inside plasma
          dndlv(2,j)=dndlv(2,j)+den_r(i)*(s_c(ii+1)-s_c(ii))
          dndlv(3,j)=dndlv(3,j)+(den_r(i)+pden_r(i))*(s_c(ii+1)-s_c(ii))

        ENDIF

      ENDDO !Over chord segments

    ENDIF

    dndlv(1,j)=dndlv(2,j)
    dndlv(4,j)=dndlv(3,j)

  ENDDO !Over vertical chords

ENDIF

!------------------------------------------------------------------------------
!Line densities for horizontal chords
!------------------------------------------------------------------------------
dndlh(:,:)=0

IF(nhorz > 0) THEN

  DO j=1,nhorz !Over horizontal chords

    !Set endpoints outside and inside plasma
    rseg_c(1,1)=r0+1.1*a0*rho_rm(n)
    rseg_c(1,2)=r0-1.1*a0*rho_rm(n)
    rseg_c(2,1:2)=0
    rseg_c(3,1:2)=zhorz(j)

    !Get intersections with plasma
    n_c=0
    CALL TRACK(n,rho_rm,2,rseg_c, &
               n_c,irho_c,s_c,iflag,message, &
               IZONE_INT=izone_c)

    !Check messages
    IF(iflag /= 0) THEN

      CALL WRITE_LINE(n_msg,message,1,1)
      IF(iflag > 0) GOTO 9999
      iflag=0
      message=''

    ENDIF

    IF(n_c > 1) THEN

      !Integrate horizontal line densities
      DO ii=1,n_c-1 !Over chord segments

        i=izone_c(ii)

        IF(i > 0 .AND. &
           i <= nc) THEN

          !Inside Plasma
          dndlh(2,j)=dndlh(2,j)+den_r(i)*(s_c(ii+1)-s_c(ii))
          dndlh(3,j)=dndlh(3,j)+(den_r(i)+pden_r(i))*(s_c(ii+1)-s_c(ii))

        ENDIF

      ENDDO !Over chord segments
   
    ENDIF

    dndlh(1,j)=dndlh(2,j)
    dndlh(4,j)=dndlh(3,j)

  ENDDO !Over horizontal chords

ENDIF

!-------------------------------------------------------------------------------
!Print to summary file
!-------------------------------------------------------------------------------
!
! MAke listing of namelist variables
!
write(n_sum,indata)

!Global information
label='*** Global Data ***'
CALL WRITE_LINE(n_sum,label,1,1)

label='Atomic mass number of pellet ions [-] ='
rdum(1)=amu_pel
CALL WRITE_LINE_IR(n_sum,label,0,idum,1,rdum,15,1)

label='Pellet effective spherical radius [m] ='
rdum(1)=r_pel
CALL WRITE_LINE_IR(n_sum,label,0,idum,1,rdum,15,1)

label='Pellet velocity [m/s] ='
rdum(1)=v_pel
CALL WRITE_LINE_IR(n_sum,label,0,idum,1,rdum,15,1)

label='Distance penetrated [m] ='
rdum(1)=rpen
CALL WRITE_LINE_IR(n_sum,label,0,idum,1,rdum,15,1)

label='Distance to axis [m] ='
rdum(1)=rtot
CALL WRITE_LINE_IR(n_sum,label,0,idum,1,rdum,15,1)

label='Fractional radius penetrated [-] ='
rdum(1)=frpen
CALL WRITE_LINE_IR(n_sum,label,0,idum,1,rdum,15,1)

label='Volume penetrated [m**3] ='
rdum(1)=volpen
CALL WRITE_LINE_IR(n_sum,label,0,idum,1,rdum,15,1)

label='Plasma volume [m**3] ='
rdum(1)=volp
CALL WRITE_LINE_IR(n_sum,label,0,idum,1,rdum,15,1)

label='Fractional volume penetrated [-] ='
rdum(1)=fvpen
CALL WRITE_LINE_IR(n_sum,label,0,idum,1,rdum,15,1)

!Profiles
!Allocate and initialize output radial arrays
ALLOCATE(valpro(n,mx_npro))

  namepro(:)=''
  unitpro(:)=''
  descpro(:)=''
  valpro(:,:)=0
  npro=0

!Grids
npro=npro+1
namepro(npro)='rho_t'
unitpro(npro)='-'
descpro(npro)='Normalized toroidal flux grid - ' &
              //'proportional to square root toroidal flux'
valpro(:,npro)=rho_r(:)

!Geometry
npro=npro+1
namepro(npro)='dVol'
unitpro(npro)='m**3'
descpro(npro)='Cell volume'
valpro(:,npro)=dvol_r(:)

!Temperatures
npro=npro+1
namepro(npro)='Te(tpel-)'
unitpro(npro)='keV'
descpro(npro)='Initial electron temperature'
valpro(:,npro)=te_r(:)

!Densities
npro=npro+1
namepro(npro)='ne(tpel-)'
unitpro(npro)='/m**3'
descpro(npro)='Initial electron density'
valpro(:,npro)=den_r(:)


!Fast ions and neutrals
IF(nf > 0) THEN

  npro=npro+1
  namepro(npro)='n0'
  unitpro(npro)='/m**3'
  descpro(npro)='Neutral density'
  valpro(:,npro)=dn0_r(:)

  DO j=1,nf !Over fast ion components

    IF(iz_f(j) == 1) THEN

      !Fast beam ions
      WRITE(cb,'(i1)') j

      npro=npro+1
      namepro(npro)='Hb('//cb//')'
      unitpro(npro)='-'
      descpro(npro)='H(r) for beam component '//cb
      valpro(:,npro)=h_rf(:,j)

      npro=npro+1
      namepro(npro)='nb('//cb//')'
      unitpro(npro)='/m**3'
      descpro(npro)='Density of fast beam component '//cb

      DO i=1,n

        valpro(i,npro)=SUM(den_ref(i,1:ne_f(j),j))

      ENDDO

    ELSE

      !Fast alphas
      npro=npro+1
      namepro(npro)='Halpha'
      unitpro(npro)='-'
      descpro(npro)='H(r) for fast alphas'
      valpro(:,npro)=h_rf(:,j)

      npro=npro+1
      namepro(npro)='nalpha'
      unitpro(npro)='/m**3'
      descpro(npro)='Density of fast alphas'

      DO i=1,n

        valpro(i,npro)=SUM(den_ref(i,1:ne_f(j),j))

      ENDDO

    ENDIF

  ENDDO !Over beam components

ENDIF


!Perturbations to plasma profles from pellet
npro=npro+1
namepro(npro)='delta_ne'
unitpro(npro)='/m**3'
descpro(npro)='Electron density perturbation'
valpro(:,npro)=pden_r(:)

npro=npro+1
namepro(npro)='Te(tpel+)'
unitpro(npro)='keV'
descpro(npro)='Final electron temperature'
valpro(:,npro)=(den_r(:)*te_r(:)-pden_r(:)*(2.0/3.0) &
                 *z_eion)/(den_r(:)+pden_r(:))

npro=npro+1
namepro(npro)='ne(tpel+)'
unitpro(npro)='/m**3'
descpro(npro)='Final electron density'
valpro(:,npro)=den_r(:)+pden_r(:)

!PRL Deposition
if (k_prl > 0) then
    npro=npro+1
    namepro(npro)='ne(prl)'
    unitpro(npro)='/m**3'
    descpro(npro)='PRL electron density perturbation'
    valpro(:,npro)=prldep(:)
endif

npro=npro+1
namepro(npro)='delta_ne/ne'
unitpro(npro)='-'
descpro(npro)='Fractional electron density perturbation'
valpro(:,npro)=pden_r(:)/den_r(:)

label='*** Profiles ***'
CALL WRITE_OUT1(n_sum,'sum',n,npro,valpro,namepro,unitpro, &
                descpro,1,LABEL=label)

!Print radial profiles to 1D data file
n_tmp=20
cn_tmp='1d_pellet_r'//cn_runid//'.dat'
OPEN(UNIT=n_tmp, &
     FILE=cn_tmp, &
     STATUS='unknown', &
     FORM='formatted')

CALL WRITE_OUT1(n_tmp,'1d',n,npro,valpro,namepro,unitpro,descpro,-1)

CLOSE(UNIT=n_tmp)

DEALLOCATE(valpro)

!Pellet path
!Allocate and initialize output radial arrays
ALLOCATE(valpro(n_p,mx_npro))

  namepro(:)=''
  unitpro(:)=''
  descpro(:)=''
  valpro(:,:)=0
  npro=0

npro=npro+1
namepro(npro)='t'
unitpro(npro)='s'
descpro(npro)='Time at entry to cell'
valpro(1:n_p,npro)=t_p(1:n_p)

npro=npro+1
namepro(npro)='rho'
unitpro(npro)='-'
descpro(npro)='Normalized minor radius at entry to cell'
valpro(1:n_p,npro)=rflx_p(1,1:n_p)

npro=npro+1
namepro(npro)='R'
unitpro(npro)='m'
descpro(npro)='Major radius at entry to cell'
valpro(1:n_p,npro)=rcyl_p(1,1:n_p)

npro=npro+1
namepro(npro)='Z'
unitpro(npro)='m'
descpro(npro)='Elevation at entry to cell'
valpro(1:n_p,npro)=rcyl_p(3,1:n_p)

npro=npro+1
namepro(npro)='s'
unitpro(npro)='m'
descpro(npro)='Distance at entry to cell'
valpro(1:n_p,npro)=s_p(1:n_p)

npro=npro+1
namepro(npro)='src'
unitpro(npro)='/s'
descpro(npro)='Pellet source rate in cell'
valpro(1:n_p,npro)=src_p(1:n_p)

npro=npro+1
namepro(npro)='rpel1'
unitpro(npro)='m'
descpro(npro)='Pellet radius at exit from cell'
valpro(1:n_p,npro)=rpel1_p(1:n_p)

npro=npro+1
namepro(npro)='Te0'
unitpro(npro)='keV'
descpro(npro)='Electron temperature at entry to cell'
valpro(1:n_p,npro)=te0_p(1:n_p)

npro=npro+1
namepro(npro)='Te1'
unitpro(npro)='keV'
descpro(npro)='Electron temperature at exit from cell'
valpro(1:n_p,npro)=te1_p(1:n_p)

npro=npro+1
namepro(npro)='ne0'
unitpro(npro)='/m**3'
descpro(npro)='Electron density at entry to cell'
valpro(1:n_p,npro)=den0_p(1:n_p)

npro=npro+1
namepro(npro)='ne1'
unitpro(npro)='/m**3'
descpro(npro)='Electron density at exit from cell'
valpro(1:n_p,npro)=den1_p(1:n_p)

label='*** Pellet path ***'
CALL WRITE_OUT1(n_sum,'sum',n_p,npro,valpro,namepro,unitpro,descpro,1, &
                LABEL=label)

!Print path profiles to 1D data file
n_tmp=20
cn_tmp='1d_pellet_p.dat'
OPEN(UNIT=n_tmp, &
     FILE=cn_tmp, &
     STATUS='unknown', &
     FORM='formatted')

CALL WRITE_OUT1(n_tmp,'1d',n_p,npro,valpro,namepro,unitpro,descpro,-1)

CLOSE(UNIT=n_tmp)

!-------------------------------------------------------------------------------
!Exit
!-------------------------------------------------------------------------------
9999 CONTINUE

END PROGRAM PELLET_DR

SUBROUTINE READ_PELLET_PRO(nin,cin,n_rho,rho_r, &
                           den_r,te_r,iflag,message)
!-------------------------------------------------------------------------------
!READ_PELLET_PRO reads a data file of profile information for the PELLET code
!
!References:
!  W.A.Houlberg, L.R.Baylor, F90 free format 8/2004
!
!Comments:
!  Data file contains rho, ne, Te, Ti only
!  Make sure the radial grids (external rho and internal xr) are defined in the
!    same way (presently proportional to square root of the toroidal flux and
!    unity at the separatrix/limiter)
!  If the input nodes are outside the data range, edge values are used in the
!    extrapolation (LINEAR1_INTERP)
!-------------------------------------------------------------------------------
USE SPEC_KIND_MOD
USE LINEAR1_MOD
IMPLICIT NONE

!Declaration of input variables
CHARACTER(len=*), INTENT(IN) :: &
  cin                    !name of unit for data [character]

INTEGER, INTENT(IN) :: &
  nin,                 & !unit number for data [-]
  n_rho                  !no. of radial grid points [-]

 REAL(KIND=rspec), INTENT(IN) :: &
  rho_r(n_rho)           !radial nodes [rho]

!Declaration of output variables
CHARACTER(len=*), INTENT(OUT) :: &
  message                !warning or error message [character]

INTEGER, INTENT(OUT) :: &
  iflag                  !error and warning flag [-]
                         !=-1 warning
                         !=0 none
                         !=1 error

REAL(KIND=rspec), INTENT(OUT) :: &
  den_r(n_rho),        & !electron density [/m**3]
  te_r(1:n_rho)          !electron temperature [keV]

!-------------------------------------------------------------------------------
!Declaration of local variables
INTEGER :: &
  i,nxr

CHARACTER(len=15) :: &
  cdummy

REAL(KIND=rspec) :: &
  tixr

REAL(KIND=rspec), ALLOCATABLE :: &
  xr(:),denxr(:),texr(:)

!-------------------------------------------------------------------------------
!Open file, get number of radial points, allocate arrays
!-------------------------------------------------------------------------------
OPEN(UNIT=nin, &
     FILE=cin, &
     STATUS='old')

!Read text header line
READ(nin,*) cdummy

!Read number of data points
READ(nin,*) nxr

!Allocate data arrays
ALLOCATE(xr(nxr), &
         denxr(nxr), &
         texr(nxr))

  xr(:)=0
  denxr(:)=0
  texr(:)=0

!-------------------------------------------------------------------------------
!Read profile data and interpolate to external grid
!-------------------------------------------------------------------------------
READ(nin,*) (xr(i),denxr(i),texr(i),tixr, i=1,nxr)

!Change density to 10^19 /m^3 units
denxr(:)=denxr(:)*(1e19)

CALL LINEAR1_INTERP(nxr,xr,denxr,n_rho,rho_r, &
                    den_r,iflag,message)

IF(iflag /= 0) THEN

  message='READ_PELLET_PRO den interpolation/'//message
  IF(iflag > 0) GOTO 9999

ENDIF

CALL LINEAR1_INTERP(nxr,xr,texr,n_rho,rho_r, &
                    te_r,iflag,message)

IF(iflag /= 0) THEN

  message='READ_PELLET_PRO te interpolation/'//message
  IF(iflag > 0) GOTO 9999

ENDIF

!-------------------------------------------------------------------------------
!Cleanup and exit
!-------------------------------------------------------------------------------
9999 CONTINUE

CLOSE(UNIT=nin)
DEALLOCATE(xr,denxr,texr)

END SUBROUTINE READ_PELLET_PRO

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
USE SPEC_KIND_MOD
USE FLUXAV_MOD
IMPLICIT NONE

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
  mxnx_xy=130, &
  mxny_xy=130,  &
  mxnr_r=300, &
  mxn_lim=200, &
  mxn_bdry=1500

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

!Check messages
IF(iflag /= 0) THEN

  message='FORCEBAL_EFIT(1)/'//message
  IF(iflag > 0) GOTO 9999

ENDIF

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

!-------------------------------------------------------------------------------
!Initialization
!-------------------------------------------------------------------------------
!Null output
iflag=0
message=''

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

!1D radial grid and and poloidal flux
psi_x(1:nx_xy)=psimag+(psilim-psimag)*(/ (i-1,i=1,nx_xy) /)/(nx_xy-1)
rhop_x(1:nx_xy)=(psi_x(1:nx_xy)-psi_x(1))/(psi_x(nx_xy)-psi_x(1))

!-------------------------------------------------------------------------------
!Cleanup and exit
!-------------------------------------------------------------------------------
9999 CONTINUE

!Close the EQDSK file
CLOSE(unit=nin)

END SUBROUTINE READ_EFIT_EQDSK
