MODULE PRL_MOD
!-------------------------------------------------------------------------------
!PRL_MOD is an F90 module of routines that support the calculation of ExB pellet
!  cloud drift in the PRL code. 
!
!References:
!
!  L.R.Baylor 9/2004
!  W.A.Houlberg, L.R.Baylor, Standardization for NTCC Module Library 1/2005
!
!Contains PUBLIC routines:
!
!  PRL                 -calculates drift of a pellet generated cloudlet
!-------------------------------------------------------------------------------
USE SPEC_KIND_MOD
IMPLICIT NONE

!-------------------------------------------------------------------------------
!Private procedures
!-------------------------------------------------------------------------------
PRIVATE :: &
  PRL_BESSI0, &
  PRL_BESSI1, &
  PRL_BESSK, &
  PRL_BESSK0, &
  PRL_BESSK1, &
  PRL_DTIME, &
  PRL_GETNORMDEP, &
  PRL_INIT, &
  PRL_INTEGRATEF, &
  PRL_LAGSETUP, &
  PRL_LFUNC, &
  PRL_OUT, &
  PRL_TRAPZD

!-------------------------------------------------------------------------------
!Private data
!-------------------------------------------------------------------------------
!Higher (double) precision required for selected calculations
INTEGER, PRIVATE, PARAMETER :: &
  dpspec=SELECTED_REAL_KIND(p=12,r=100)

!Integer parameters
INTEGER, PRIVATE, PARAMETER :: &
  na_pr=200,           & !Axial cells [-]
  nr_pr=100,           & !Radial zones [-]
  ns_pr=500,           & !Saved time points [-]
  nout = 12              !Diagnostic output lun 

!Physical and conversion constants
REAL(KIND=rspec), PRIVATE, PARAMETER :: &   
  z_pi=3.141592654       !pi [-]

!Parameters
REAL(KIND=rspec), PARAMETER :: &
  gamma=5.0/3.0 , &      !gamma = 5/3
  mu0 = 4.0*z_pi*1e-7, & !magnetic permeability in H/m
  gammasqrt=1.291        !sqrt(gamma)

!Logical switches
LOGICAL, PRIVATE :: &
  l_chi_pr  =.TRUE.,   & ! [default=.TRUE.]
  l_cltem_pr=.TRUE.,   & !Limit cloud temp from reaching ambient [default=.TRUE.]
  l_mcor_pr =.TRUE.,   & !Mach number correction [default=.TRUE.]
!  l_mcor_pr =.FALSE.,   & !Mach number correction [default=.TRUE.]
!  l_mom_pr  =.TRUE.,   & !Momentum flux correction [default=.TRUE.]
  l_mom_pr  =.FALSE.,   & !Momentum flux correction [default=.TRUE.]
  l_msnew_pr=.FALSE.,   & !New mass shedding model [default=.TRUE.]
  l_msold_pr=.FALSE.,  & !Old mass shedding model [default=.FALSE.]
  l_out_pr  =.FALSE.,  & !Output profiles to file [default=.FALSE.]
  l_debug_pr=.TRUE.,  &  !Output debug info to screen [default=.FALSE.]
  l_rad_pr  =.TRUE.,   & !Include radiation loss term [default=.TRUE.]
  l_tesig_pr=.TRUE.,   & !Electron temperature dependent sigma [default=.TRUE.]
  l_awdamp_pr=.FALSE., & !Alfven wave damping ala Pegourie [default=.FALSE.]
  l_i2switch_pr=.FALSE.  !

!Real input stored for general use in module
REAL(KIND=rspec), PRIVATE :: &
  bt_pr,rm_pr,am_pr,rpel_pr,rhop_pr

!Time dependent data
REAL(kind=dpspec), PRIVATE, DIMENSION(0:ns_pr) :: &
  psip_pr,ptime_pr,vcld_pr

!Radial data
REAL(KIND=rspec), PRIVATE, DIMENSION(na_pr) :: &
  etaplus_pr,etaminus_pr,gplus_pr,gminus_pr,uplus_pr,uminus_pr,uplusint_pr, &
  zplus_pr,zminus_pr, energy_pr, vol_pr

REAL(KIND=dpspec), PRIVATE, DIMENSION(0:na_pr+1) :: &
  delq_pr,dencon_pr,deni_pr,denn_pr,deno_pr,etahat_pr,gdiff_pr,ghat_pr, &
  khat_pr,presn_pr,preso_pr,tempi_pr,tempn_pr,tempo_pr,xcenn_pr,xceno_pr, &
  theta_pr, psis, delbcs

REAL(KIND=dpspec), PRIVATE, DIMENSION(na_pr+1) :: &
  vhatn_pr,vhato_pr,xinti_pr,xintn_pr,xinto_pr

!Scalar data
CHARACTER(LEN=16), PRIVATE :: &
  cidrun_pr

INTEGER, PRIVATE :: &
  icycle_pr, mcycle_pr, ocycle_pr, isave

REAL(KIND=dpspec), PRIVATE :: &
afact_pr,cs0_pr,delt_pr,gacc_pr,inert_pr,lc_pr,n0_pr,presbc_pr,pratio_pr, &
qhat_pr,rin_pr,sigma_pr,sigma0_pr,tau_pr,taubar_pr,tcloud_pr,wpel_pr, loglambda, sigma_spitz

REAL(KIND=rspec), PRIVATE :: &
chi0_pr,kappac_pr,neinf_pr,teinf_pr


!Storage arrays
REAL(KIND=rspec), PRIVATE, DIMENSION(ns_pr) :: &
  presbc_arr,presNC,tau_arr,u_arr,rhotw_arr, u1_arr,u1i2_arr, u2_arr,u3_arr,psi_arr, cos_arr, &
  Gradn_arr, Ac_arr, ut_arr, uchi_arr, chin_arr, Fluxn_arr, Ndiff_arr, i2_arr, i1_arr, nu1_arr, nu2_arr, &
  n1_arr, n2_arr, n3_arr, T1_arr, T2_arr, T3_arr, LC3_arr

REAL(KIND=rspec), PRIVATE, DIMENSION(na_pr) :: &
   Lc1, Lc2, Lc3, Pr1, Pr2,Pr3

REAL(KIND=rspec), PRIVATE, DIMENSION(ns_pr,na_pr) :: &
  beta_arr,dtw_arr,Lc_arr,Mnum_arr,phi_arr,rho_arr,Rmtw_arr,rtw_arr,thtw_arr, &
  vhat_arr, Q_arr

!-------------------------------------------------------------------------------
!Procedures
!-------------------------------------------------------------------------------
CONTAINS

SUBROUTINE PRL(irunid,idiag,Bt,Rm,am,chi0,rpel,rhop,pnum,vpel,n_r,te_r,ne_r,q_r, &
               tea_r,nea_r,rho_r,nump)
!-------------------------------------------------------------------------------
!PRL is a 1D Lagrangian finite difference calculation of the pressure relaxation
!  model of Parks for the drift of a pellet generated cloudlet
!
!References:
!  W.D.Sessions, Tennessee Tech Univ, original coding
!	 L.R. Baylor, PRL_LAGSETUP included as part of code  
!    Mach # effect
!    Radiation
!    Momentum flux
!    Te dependent (time dependent) incident heat flux
!    Mass shedding
!    Enhanced Heating test 5/2003
!    New Mass shedding 11/2003
!    New new mass shedding from APS poster 3/2004
!    Updated variable names and using new notation 4/2004
!    Completed analytical test for 1st shed cell 4/2004
!    Modified for F90 inclusion in latest version of PELLET 9/2004
!  W.A.Houlberg, L.R.Baylor, Standardization for NTCC Module Library 1/2005
!    Added energy in cloudlet that is lost when exiting plasma. LRB 5/2008
!
!Comments:
!  Uses a duct profile with time-varying boundary conditions (pres_ambient)
!  Key variables used:
!    denn,deno,deni - density of cloud cells (n=new, o=old, i=intial)
!    tempn,tempo,tempi - temperature of cloud cells (n=new, o=old, i=intial)
!    presbc - pressure at boundary
!    lc - initial cloud half length
!    rhotw,rhotwo - r/a (o=old)
!    chin,chiold - cloud poloidal coordinate
!    vhatn,vhato - parallel velocity of cloud cells normalized to c0bar.
!    unew,uold - velocity of pellet cloud (normalized) in minor radius direction
!    uchinew, uchiold - velocity of pellet cloud (normalized) in poloidal direction
!    uchiint - integral of u_chi/rho dt used in unew eqn of motion.
!    Mtw - normalized mass of cloud relative to original mass
!    Erpl - radial electric field in plasma (V/m)
!-------------------------------------------------------------------------------

!Declaration of input variables
INTEGER, INTENT(IN) :: &
  irunid,              & !run id number (0 if stand alone code call) [-]
  idiag,               & !level of diagnostic output written to file [-]
  n_r                    !number of points in profile arrays [-]

REAL(KIND=rspec), INTENT(IN) :: &
  Bt,                  & !magnetic field on axis [T]
  Rm,                  & !major radius [cm]
  am,                  & !minor radius for circular approximation [cm]
  chi0,                & !angle of injection trajectory (HFS=pi, LFS=0) [rad]
  rpel,                & !pellet radius [cm]
  rhop,                & !rho location of pellet cloudlet [-]
  pnum,                & !particles in pellet cloudlet [-]
  vpel                   !pellet velocity [cm/s]

REAL(KIND=rspec), INTENT(IN) :: &
  rho_r(:),            & !radial grid normalized to am (axis=0, edge = 1) [-]
  te_r(:),             & !initial electron temperature profile [eV]
  ne_r(:),             & !initial electron density profile [/cm**3]
  q_r(:)                 !safety factor profile [-]

!Declaration of output variables
REAL(KIND=rspec), INTENT(OUT) :: &
  nump                   !number of particles calculated in cloudlet 

REAL(KIND=rspec), INTENT(OUT) :: &
  tea_r(:),            & !electron temperature profile after cloudlet motion [eV]
  nea_r(:)               !electron density profile after cloudlet motion [/cm**3]

!-------------------------------------------------------------------------------
!Declaration of local variables     
!
CHARACTER(len=2) :: &
  runid2

CHARACTER(len=16) :: &
  tmprid

INTEGER :: &
  i,ixloc,j,jjmin,jsave,jsp,jsp2,k,na,Numshed, dummy,ii

REAL(KIND=dpspec) :: &
  delbc,delxn,delxo,psi,taumax

REAL(KIND=rspec) :: &
  alphar,chin,chiold,cltemfact,delpress,delr,delxc,dendelm,eps,Erpl,fr, &
  fracpsi,fracion,h,hfact,Ke,lambda_ef,lambda_ef0,Lsloc,Mcorfact,momterm,Mtw, &
  Ncl,nudrag,presbc0,q0,qa,qfact,qhatt,qloc,radterm,rcloud,rho_0,rhotw,rhotwo, &
  rpresr,shedlevel,taubarn,te_inf0,uchiint,uchinew,uchiold,uexb,unew,uold, &
  vexb,xloc, unew1, unew1i2, unew2, unew3, unewtot, uchi1, uchi2, uchi3, duchi, cAinf, &
  etot, i2, i2new, ntw, Utot, ne_inf0, etap, Rp, Lp, Za, nu1, nu2, tstar, &
  I1, i_1, Istar, Vsource

!Radial arrays
REAL(KIND=rspec), DIMENSION(nr_pr) :: &
  Ls,frad,massprof,ne_prl,q_prl,r_prl,shat,te_prl

!Axial arrays
REAL(KIND=dpspec), DIMENSION(0:na_pr+1) :: &
  Mnum

REAL(KIND=rspec), DIMENSION(na_pr) :: &
  chiang,delxint,deltatw,massloc,q,Rmtw,rperp,rtw

!-------------------------------------------------------------------------------
!Initialization
!-------------------------------------------------------------------------------
na=na_pr                              !Start with all axial cells being used
Numshed=0                             !Initialize the mass shedding index
IF(ocycle_pr == 0) ocycle_pr=2000     !Default to 1000 steps for output cycle
IF(mcycle_pr == 0) mcycle_pr=10000000 !Default to 10M steps for output cycle
IF(taumax == 0) taumax=25             !Default to max of 100 tau steps
IF(shedlevel == 0) shedlevel=1.0      !Default shedding factor
tau_pr=0
jsave=0
icycle_pr=0
psi=0                                 !Initialize psi integral
hfact=1                               !Heating enhancement factor (for testing)
isave = 100                           !Save increment for diagnostic output


IF(idiag>0) then
    l_out_pr = .TRUE.         !Turn on diagnostic output
    l_debug_pr = .TRUE.       !Turn on debugging output
endif

!Parameters for perpendicular transport
rpresr=1/4.9
Erpl = 0         ! Eradial in units of V/m
Ke = 0.73

!Load inputs into PRL variables
bt_pr=Bt
chi0_pr=chi0

!chi0_pr = 3.1415-rhop*0.975    !Test for FTU/D3D_V+1

rm_pr=Rm
am_pr=am
rpel_pr=rpel
rhop_pr=rhop

!Linear interpolation to get Te, ne, and q profile on PRL radial grid
jjmin=1

DO i=1,nr_pr

  fr=(REAL(i)/nr_pr)
  j=n_r

  DO WHILE (j >= jjmin)
    IF(rho_r(j) >= fr) THEN
      te_prl(i)=te_r(j)
      ne_prl(i)=ne_r(j)
      q_prl(i)=q_r(j)
      IF(rho_r(j) == fr) jjmin=j
      j=j-1
    ELSE

      te_prl(i)=te_r(j)-((te_r(j)-te_r(j+1))/(rho_r(j)-rho_r(j+1))) &
                        *(rho_r(j)-fr)
      ne_prl(i)=ne_r(j)-((ne_r(j)-ne_r(j+1))/(rho_r(j)-rho_r(j+1))) &
                        *(rho_r(j)-fr)
      q_prl(i)=q_r(j)-((q_r(j)-q_r(j+1))/(rho_r(j)-rho_r(j+1))) &
                      *(rho_r(j)-fr)
      jjmin=j
      j=j-1

    ENDIF
  ENDDO

ENDDO

!Plasma Te and ne values at starting location in plasma
DO i=1,nr_pr

  IF(REAL(i)/nr_pr < rhop) THEN
    teinf_pr=te_prl(i)
    neinf_pr=ne_prl(i)
    qloc=q_prl(i)
  ENDIF

ENDDO

!Calculate s and Ls on PRL radial grid for mass shedding model
delr=am_pr/nr_pr

DO i=2,nr_pr-1
  r_prl(i)=i*delr
  shat(i)=(r_prl(i)/q_prl(i))*((q_prl(i+1)-q_prl(i-1))/(2*delr))
  Ls(i)=q_prl(i)*(rm_pr+COS(chi0)*r_prl(i))/shat(i)
  IF (Ls(i) > 10000) Ls(i)=10000.0
ENDDO

r_prl(1)=0
shat(1)=shat(2)
Ls(1)=Ls(2)
r_prl(nr_pr)=r_prl(nr_pr-1)
shat(nr_pr)=shat(nr_pr-1)
Ls(nr_pr)=Ls(nr_pr-1)

!***
!neinf_pr = 1.0e14
!teinf_pr = 5000.0  !5 kev


!Call PRL_LAGSETUP to define cell bounday positions and heat flux  LRB  Apr-2001
CALL PRL_LAGSETUP(neinf_pr,teinf_pr)

uold=vpel/cs0_pr/gammasqrt       !Initialize velocity to pellet velocity (neg inward velocity)
unew=uold
Utot=unew
uchiold=0
uchinew=uchiold
chiold=chi0
eps=am_pr/rin_pr
rhotw=rhop                       !Initial minor radius of cloudlet
rhotwo=rhotw
h=1                              !Parks FAX 5-Apr-2004
rcloud=10*rpel_pr                !Estimated cloud radius = 10*rpel  WDS
rcloud=kappac_pr*rpel_pr         !From Parks PoP 2000
Ncl=2*lc_pr*n0_pr*z_pi*rcloud**2 !Particles in initial cloud
nump=Ncl

i2 = 0                           !Initial Alfven short circuit - Parks Oct2008

CALL PRL_INIT
rpresr=1/pratio_pr 
te_inf0=teinf_pr
ne_inf0=neinf_pr
presbc0=preso_pr(na_pr+1)
k=1

!Calculate heat source term from heat deposition  (Parks, PoP Eqn 13)
khat_pr(1:na_pr)=sigma_pr*qhat_pr*ghat_pr(1:na_pr)

!? write stuff
!IF(l_out_pr) THEN      !If profile output desired, open output file
!  CALL PRL_OUT(0)      !Open output file
!  CALL PRL_OUT(1)      !Write data to output file
!ENDIF

!Determine run identification and display starting info
IF(irunid > 0) THEN    !Build new runid string

  WRITE(runid2,'(i2)') irunid
  jsp=INDEX(cidrun_pr,CHAR(0))
  IF((jsp >0) .OR. (irunid>0)) cidrun_pr = 'Run      '
  jsp=INDEX(cidrun_pr,CHAR(32))
  tmprid=cidrun_pr(1:14)
  IF(jsp>0) tmprid=cidrun_pr(1:jsp-1)
  jsp2=INDEX(runid2,CHAR(32))
  cidrun_pr=tmprid(1:jsp-1)//runid2(jsp2+1:2)

ENDIF

!? write stuff
if(l_out_pr) then
   write(*,*) '  '
   write(*,*) '  PRLagrangian output for runid: '//cidrun_pr
   write(*,*) '  '
endif
open(unit=15, file = 'PRLout_'//cidrun_pr//'.dat')
write(15,*) '  PRLagrangian output for runid: '//cidrun_pr
write(15,*) ' '
write(15,*) ' '
write(15,*) ' Scaling factors: '
write(15,150) 'Lc(cm): ', lc_pr
write(15,150) 'c0(cm/s): ', cs0_pr/gammasqrt
write(15,150) 'Lc/c0(s): ', lc_pr/(cs0_pr/gammasqrt)
write(15,150)  'Ncl: ', Ncl
write(15,150)  'Pnum: ', pnum
write(*,150)  'Ncl(PRL): ', Ncl
write(*,150)  'Pnum(PELLET): ', pnum
write(15,*) ' '
write(15,*) ' '
write(15,*)'tau            psi            L           u',  &
           '         presbc         r/a'
write(*,*)'tau            psi            L           u',  &
          '         presbc         r/a'
150 format(1x,a10,e12.5)

!Set new values to old ones for first time step:
xintn_pr(1:na+1)=xinto_pr(1:na+1)
vhatn_pr(1:na+1)=vhato_pr(1:na+1)        
denn_pr(1:na)=deno_pr(1:na)
tempn_pr(1:na)=tempo_pr(1:na)
presn_pr(1:na)=preso_pr(1:na)
xcenn_pr(1:na)=xceno_pr(1:na)

!Initial rperp for mass shedding
rperp(1:na)=0

!Initial debug variables to zero
rho_arr(1:ns_pr,1:na_pr)=0       ! Density
phi_arr(1:ns_pr,1:na_pr)=0       ! Temp
beta_arr(1:ns_pr,1:na_pr)=0      ! Press
vhat_arr(1:ns_pr,1:na_pr)=0      ! Vel_axial
Lc_arr(1:ns_pr,1:na_pr)=0        ! Length Cloud
dtw_arr(1:ns_pr,1:na_pr)=0
thtw_arr(1:ns_pr,1:na_pr)=0
rtw_arr(1:ns_pr,1:na_pr)=0
Rmtw_arr(1:ns_pr,1:na_pr)=0
Mnum_arr(1:ns_pr,1:na_pr)=0      ! Mach num
u_arr(1:ns_pr)=0                 ! V rad
u1_arr(1:ns_pr)=0                ! V1 rad
u2_arr(1:ns_pr)=0                ! V2 rad
u3_arr(1:ns_pr)=0                ! V3 rad
ut_arr(1:ns_pr)=0                ! Vt rad
uchi_arr(1:ns_pr)=0              ! 
chin_arr(1:ns_pr)=0
cos_arr(1:ns_pr)=0               ! V rad
psi_arr(1:ns_pr)=0               ! psi 
presNC(1:ns_pr)=0                ! Press cloud end
presbc_arr(1:ns_pr)=0            ! Press plasma boundary condition
tau_arr(1:ns_pr)=0               ! tau
rhotw_arr(1:ns_pr) = 0			 ! r/a cloud location
i2_arr(1:ns_pr) = 0				 ! i2 sc current
nu1_arr(1:ns_pr) = 0			 ! nu1 resistive drag term
nu2_arr(1:ns_pr) = 0			 ! nu2 radiation resistance drag term


!-------------------------------------------------------------------------------
!Outer loop
!-------------------------------------------------------------------------------
LOOP_J: DO j=1,mcycle_pr
         
  IF(tau_pr > taumax .OR. &
     !na < 5 .OR. xintn_pr(na) > 60) THEN !Finished
     na < 5 ) THEN !Finished

    !Store final location of cells in rho space      
    massloc(1:na)=rhotw

    !Return normalized mass depos profile to caller in nea_r array
    CALL PRL_GETNORMDEP(n_r,pnum,massloc,massprof,frad,rho_r,nea_r)

    !? write stuff
    write(15,160) tau_pr,psi,xintn_pr(na+1),unew,presbc_pr,rhotw
    160 format(e12.5,1x,e12.5,1x,e12.5,1x,e12.5,1x,e12.5,1x,e12.5)
    161 format(e12.5,1x,e12.5,1x,e12.5,1x,e12.5,1x,e12.5,1x,e12.5,1x,e12.5)
 
    IF(l_out_pr) THEN
      CALL PRL_OUT(idiag,nout)  !Write output file
      !CALL PRL_OUT(2)  !Close file
    ENDIF

    EXIT LOOP_J
                   
  ENDIF  

  !Calculate max time step via Courant condition
  CALL PRL_DTIME(na)

  !Find q at current cloudlet location
  xloc=rhotw ! r/a
  ixloc=INT(xloc*100.0)
  IF(ixloc < 1) ixloc=1

  !Determine starting point for each cloudlet for mass shedding
  IF(l_msold_pr .AND. &
     j == 1) THEN

    rho_0=xloc !rho_0 = r(0)/a

    DO i=1,na

      Rmtw(i)=(rm_pr+xloc*am_pr*COS(theta_pr(i)))/rin_pr
      theta_pr(i)=chi0+i*COS(chi0)*lc_pr/na_pr/(Rmtw(i)*rin_pr*qloc)
      rtw(i)=1
      deltatw(i)=0
      q(i)=qloc

    ENDDO

  ENDIF

  !Determine theta angle and perp distance of each cell 
  IF(l_msold_pr .AND. &
     j > 1) THEN

    !Calculate new r, delta, theta, R 
    DO i=1,na

      rtw(i)=rtw(i)+((afact_pr/(eps*rho_0))*unew*COS(theta_pr(i)))*delt_pr

      IF(chi0 > z_pi/2) THEN

        theta_pr(i)=theta_pr(i)-(afact_pr*vhatn_pr(i)/(q(i)*Rmtw(i)))*delt_pr &
                    -(afact_pr/(eps*rho_0)*unew*SIN(theta_pr(i))/rtw(i))*delt_pr
        deltatw(i)=deltatw(i)+2/(taubar_pr*gacc_pr)*unew*COS(theta_pr(i)/2) &
                   *delt_pr

      ELSE

        theta_pr(i)=theta_pr(i)+(afact_pr*vhatn_pr(i)/(q(i)*Rmtw(i)))*delt_pr &
                    -(afact_pr/(eps*rho_0)*unew*SIN(theta_pr(i))/rtw(i))*delt_pr
        deltatw(i)=deltatw(i)+2/(taubar_pr*gacc_pr)*unew*SIN(theta_pr(i)/2) &
                   *delt_pr

      ENDIF

      !remove 0.0
      IF(chi0 > z_pi/2) THEN

        Rmtw(i)=Rmtw(i)+afact_pr*unew*delt_pr+eps*afact_pr*rho_0 &
                *rtw(i)/(q(i)*Rmtw(i))*vhatn_pr(i)*SIN(theta_pr(i))*delt_pr

      ELSE

        Rmtw(i)=Rmtw(i)+afact_pr*unew*delt_pr-eps*afact_pr*rho_0 &
                *rtw(i)/(q(i)*Rmtw(i))*vhatn_pr(i)*SIN(theta_pr(i))*delt_pr

      ENDIF

      !**q0, qa, and qfact have not been given values, obsolete coding
      !q(i)=q0+qa*(rtw(i)*xloc)**qfact  !q profile is known
      !q(i)=260.0/(Rmtw(i)*Rin)  !Keep Rq constant

    ENDDO

  ENDIF

  !Check if new mass shedding scheme is to be used
  IF(l_msnew_pr .AND. &
     j == 1) THEN

    rho_0=xloc !rho_0 = r(0)/a

    DO i=1,na

      Rmtw(i)=(rm_pr+COS(chi0)*xloc*am_pr)/rin_pr  !Assume at pellet location
      rtw(i)=1
      deltatw(i)=0
      chiang(i)=0
      q(i)=qloc

    ENDDO

  ENDIF

  !New shedding criterion - Mar 2004 from APS 2003  
  IF(l_msnew_pr .AND. &
     j > 1) THEN

    !Calculate new rperp 
    DO i=1,na

      !Update R location of cloudlet, Assume at pellet location
      Rmtw(i)=(rm_pr+COS(chi0)*xloc*am_pr*COS(chiang(i)))/rin_pr

      !Determine new delta shedding parameter, was ABS(unew) Sep 2004 LRB
      deltatw(i)=deltatw(i)+(lc_pr**2)/(2*rcloud*rm_pr)  &
                 *ABS(shat(ixloc)/(q_prl(ixloc))*xcenn_pr(i))*delt_pr*(-unew)  

    ENDDO

  ENDIF

  !Do mass shedding if deltatw >= shedlevel (nominally 1)
  IF(l_msnew_pr .OR. &
     l_msold_pr) THEN

    DO i=1,na

      IF(deltatw(i) >= shedlevel) THEN

        massloc(i)=rhotw ! r/a value of shed mass
        Numshed=Numshed + 1

        IF(na > 4) THEN

          !Shed a cell
          na=na-1
          Write(*,*) ' Cell shed: ', Numshed,j,tau_pr
          Write(*,*) '     shed location: ', rhotw
!          Write(15,*) ' Cell shed: ',Numshed,j,tau_pr
!          Write(15,*) '     shed location: ', rhotw

        ENDIF

      ENDIF

    ENDDO

  ENDIF

  !Momentum equation       
  DO i=2,na

    delpress=preso_pr(i)- preso_pr(i-1)             
    dendelm=(dencon_pr(i)+dencon_pr(i-1))/2
    vhatn_pr(i)=vhato_pr(i)-(delt_pr/dendelm)*(delpress) 

    IF(l_mom_pr) THEN

      !Correct momentum flux
      vhatn_pr(i)=vhatn_pr(i) &
                  +presbc_pr*(1.0+Ke)/2*EXP(-1.8)/4*sigma_pr*gdiff_pr(i)*delt_pr

    ENDIF 

    !? debug stuff
    IF((vhatn_pr(i) > 1.0e6) .and. l_debug_pr) THEN
      WRITE(*,*) 'vhatn(i) > 1e6', i, vhatn_pr(i)
    ENDIF

  ENDDO
        
  vhatn_pr(1)=0
  vhatn_pr(na+1)=vhato_pr(na+1)-(delt_pr/dencon_pr(na+1)) &
                                *(presbc_pr-preso_pr(na))

  IF(l_mom_pr) THEN

    !Momentum flux correction 
    vhatn_pr(na+1)=vhatn_pr(na+1)+presbc_pr*(1.0+Ke)/2*EXP(-1.8)/4 &
                                  *sigma_pr*gdiff_pr(na+1)*delt_pr

  ENDIF 

  !Update interface positions via new velocities:
  xintn_pr(1:na+1)=xinto_pr(1:na+1)+delt_pr*vhatn_pr(1:na+1)          

  !Now check to see which cells are under compression, if any!
  DO i=1,na        

    delxo=xinto_pr(i+1)-xinto_pr(i)
    delxn=xintn_pr(i+1)-xintn_pr(i)
	vol_pr(i) = 2*z_pi*rcloud**2 * Lc_pr * delxn   ! 2 is for both sides of cloudlet (cm^3)  LRB 5/2008
    vol_pr(na) = vol_pr(na-1)   ! Fix last cell volume
          
    IF(delxn >= delxo) THEN          

      delq_pr(i)=0            

    ELSE          

      delq_pr(i)=1            

    ENDIF        
 
  ENDDO

  !Now redo the momentum equation to add numerical viscosity for compressed cells
  DO i=2,na              

    IF(delq_pr(i-1) == 1.0) THEN              

      delq_pr(i-1)=2*((vhato_pr(i)-vhato_pr(i-1))**2)*deno_pr(i-1)                 

    ELSE               

      delq_pr(i)=0                 

    ENDIF              

    delpress=preso_pr(i)-preso_pr(i-1)             
    dendelm=(dencon_pr(i)+dencon_pr(i-1))/2
    vhatn_pr(i)=vhato_pr(i)-(delt_pr/dendelm)*(delpress+delq_pr(i-1))                

    IF(l_mom_pr) THEN

      !Momentum flux correction 
      momterm=presbc_pr*(1.0+Ke)/2*EXP(-1.8)/4*sigma_pr*gdiff_pr(i)*delt_pr
      vhatn_pr(i)=vhatn_pr(i)+momterm             

    ENDIF 

    IF(vhatn_pr(i) < 0.0) vhatn_pr(i)=0

  ENDDO
        
  vhatn_pr(1)=0
              
  IF(delq_pr(na) == 1.0) THEN                
    delq_pr(na)=2*((vhato_pr(na+1)-vhato_pr(na))**2)*deno_pr(na)                
    vhatn_pr(na+1)=vhato_pr(na+1)-(delt_pr/dencon_pr(na+1)) &          
                                  *((presbc_pr-preso_pr(na))+delq_pr(na))
  ELSE              
    vhatn_pr(na+1)=vhato_pr(na+1)-(delt_pr/dencon_pr(na+1)) &
                                  *(presbc_pr-preso_pr(na))                         
  ENDIF

  IF(l_mom_pr) THEN
    !Momentum flux correction 
    momterm=presbc_pr*(1.0+Ke)/2*EXP(-1.8)/4*sigma_pr*gdiff_pr(na+1)*delt_pr
    vhatn_pr(na+1)=vhatn_pr(na+1)+momterm
  ENDIF 

  !Calculate Mach number for corrections to j_perp  (Parks Note 8Jun2000)  LRB
  !Divide by T^1/2 to get the local cs, not cs0.  Jun2003 LRB      
  Mnum(1:na+1)=vhatn_pr(1:na+1)/SQRT(tempo_pr(1:na+1))  
              
  !Update interface positions via new velocities:
  xintn_pr(1:na+1)=xinto_pr(1:na+1)+delt_pr*vhatn_pr(1:na+1)          
              
  !Cell-centered positions
  DO i=1,na      
    xcenn_pr(i)=(xintn_pr(i+1)+xintn_pr(i))/2     
  ENDDO
   
  !Calculate new densities via change in delxint:
  DO i=1,na        
    delxint(i)=xintn_pr(i+1)-xintn_pr(i)
    denn_pr(i)=dencon_pr(i)/delxint(i) 

    IF(denn_pr(i) < 0.0) THEN  ! Make sure density not negative!

      denn_pr(i)=deno_pr(i)
      !? debug stuff
      if (l_debug_pr) then
	     write(*,*) 'denn < 0', j, i
      endif
    ENDIF
  ENDDO

  !Energy equation (cell-centered equation):
  !Radiation term
  alphar=3.8e-20*n0_pr*lc_pr*SQRT(wpel_pr)/tcloud_pr

  DO i=1,na          

    !Radiation term - fractional ionization - Parks FAX 26-Jun-2002
    fracpsi=(3.e21/(denn_pr(i)*n0_pr))*EXP(-13.25/(tempo_pr(i)*tempi_pr(i)))

    IF(fracpsi < 0.0) THEN   ! Make sure not less than zero

      !? debug stuff
      if (l_debug_pr) then
         write (*,*) 'Fracpsi,i = ', fracpsi,i
      endif
      fracpsi=0

    ENDIF

    IF((fracpsi**2 + 4*fracpsi) > 0.0) THEN
      fracion=(SQRT(fracpsi**2+4*fracpsi)-fracpsi)/2

      IF(fracion < 0.0) fracion=0
    ELSE
      fracion=0
    ENDIF

    IF(tempo_pr(i) < 0.0) tempo_pr(i)=0

    radterm=2*alphar/3*(denn_pr(i)**2)*SQRT(tempo_pr(i))*fracion &
            *(1.0+ 33.0/(tempo_pr(i)*tempi_pr(i)))*delt_pr

    !Momentum flux term - Parks FAX 18-Jun-2002
    momterm=EXP(-1.8)/4*presbc_pr*(1.0+Ke)/2*etahat_pr(i) &
            *(denn_pr(i)-deno_pr(i))/denn_pr(i)

    !Check whether cloud temp should be limited from reaching ambient
    cltemfact=1
    IF(l_cltem_pr) cltemfact=(1.0-(tempo_pr(i)*tempi_pr(i))/teinf_pr)
    presn_pr(i)=preso_pr(i)+((5*preso_pr(i)/3+2*delq_pr(i)/3)/denn_pr(i)) &
                *(denn_pr(i)-deno_pr(i))+2*delt_pr/3*denn_pr(i)*khat_pr(i) &
                *cltemfact
                
    !Check whether to include radiation loss term
    IF(l_rad_pr) presn_pr(i)=presn_pr(i)-radterm

    !Check whether to include momentum flux heating term
    IF(l_mom_pr) presn_pr(i)=presn_pr(i)+momterm

  ENDDO         
        
  !Update cycle counter 
  icycle_pr=icycle_pr+1

  !Update time
  tau_pr = tau_pr+delt_pr

  ! Calculate new temps from denn and presn values:
  DO i=1,na      
    tempn_pr(i) = presn_pr(i)/denn_pr(i)
    IF(tempn_pr(i) < 0.0) tempn_pr(i) = tempn_pr(i-1)
  ENDDO      

  !Calculate psi function (actually psi twiddle - psi/(p0*Lc)    
  !(adjusted for Mach num correction  LRB)    
  psi=0
  uchiint=uchiint+(lc_pr/am_pr)*uchinew/rhotw   ! lc_pr/am_pr = h ???

  DO i=1,na            

    delxc=xinti_pr(i+1)-xinti_pr(i)		!distance between cells

    IF(l_mcor_pr) THEN

      Mcorfact = 1.0+(vhatn_pr(i)**2)/2/tempo_pr(i) ! Same as 1+(Mnum(i)**2)/2
      delbc = presn_pr(i)*Mcorfact-presbc_pr     

    ELSE            

      delbc = presn_pr(i)-presbc_pr    

    ENDIF             

    psi = psi + (deni_pr(i)/denn_pr(i))*delbc*delxc 
	psis(i) = psi
	delbcs(i) = delbc
	
	!***   &
        !*COS((lc_pr/(q_prl(ixloc)*rm_pr))*xcenn_pr(i))  
    !***    *COS((lc_pr*xintn_pr(i)/(q_prl(ixloc)*rm_pr))*xcenn_pr(i))  
		! Eqn 25 of PoP 2000
		! Note this is psi_twiddle - see PoP2000 Eqn 25.
                      
  ENDDO 
  
  if (psi < 0.0) then
     psi = psi
  endif       
                   
  !Update old values to new ones:
  xinto_pr(1:na+1)=xintn_pr(1:na+1)
  vhato_pr(1:na+1)=vhatn_pr(1:na+1)        
  deno_pr(1:na)=denn_pr(1:na)
  tempo_pr(1:na)=tempn_pr(1:na)
  preso_pr(1:na)=presn_pr(1:na)
  xceno_pr(1:na)=xcenn_pr(1:na)        
      
  !Set BC for increasing pressure, i.e. cloudlet moves thru ambient plasma
  ! and progress perp transport variables.
  rhotw=rhotwo + uold*delt_pr*lc_pr/am_pr

  h=1                         !April 5, 2004  Parks FAX
  nudrag = 1.0/inert_pr       !Same as Lc/Lskid definition from Jun2004

  cAinf = 2.18e11*(2.**(-0.5))*(neinf_pr**(-0.5))*bt_pr*10**4  ! Alfven velocity from NRL PF in cm/s
  nudrag = (2.0/taubar_pr)*(neinf_pr/n0_pr)*(cAinf/cs0_pr/1.265)  !from PRL 2005 formula

  ! Alfven drag terms
  loglambda = 17                           ! Approx for ln Lambda
  sigma_spitz = (teinf_pr**1.5)/(loglambda*1.65e-9)   ! (ohm-m)^-1  from Wesson 2.10
  etap = 1/sigma_spitz    ! Parallel resistivity from Spitzer (ohm-m)
  
  Rp = etap*2*z_pi*am_pr*rm_pr/rcloud**3   ! Pegourie resistive circuit element formulation
  Lp = 1*z_pi*qloc*mu0*rm_pr			   ! Parks self inductance formula
  Za = (z_pi/16)*mu0*cAinf
  Vsource = (z_pi/2)*Bt*rcloud*Utot*cs0_pr
  Istar = 8*cs0_pr*rcloud*Bt/(cAinf*mu0)
  I1 = Vsource/Za
  i_1 = I1/Istar
  
  tstar = lc_pr/(cs0_pr/gammasqrt)		   ! time normalization factor
  nu1 = tstar*Rp/Lp
  nu2 = tstar*Za/Lp
  ntw = neinf_pr/ne_inf0				   ! Density normalized to starting plasma density


!
! Turn on i2 when tau reaches 5 microseconds
!
  if (tau_pr*tstar > 5e-5) then
      l_i2switch_pr = .TRUE.
  else
      i2 = 0.0
  endif

  if (l_i2switch_pr) then
	  i2new = i2 + delt_pr*(nu2*Utot/ntw**0.5 - nu1*i2)
	  i2 = i2new 	  
	  
  endif
 

  chin=chi0
  IF(l_chi_pr) chin=chiold + uchiold*delt_pr*(lc_pr/am_pr)/rhotw
  Mtw=REAL(na)/na_pr
  vexb=(Erpl/bt_pr)*100       !cm/s from m/s
  uexb=vexb/cs0_pr/gammasqrt  !Normalize to c0

  !Include new centripetal term - Jun2004
  unew=(1.0-delt_pr*nudrag*(1+i2/(Utot*ntw**0.5))/Mtw)*uold + (delt_pr/Mtw)*gacc_pr*psi*COS(chin) &
       + delt_pr*(lc_pr/am_pr)*(uchinew**2)/rhotw

  unew1= (-delt_pr*nudrag*(1+i2/(Utot*ntw**0.5))/Mtw)*uold     !Save individual terms for debugging
  unew1i2 = i2/(Utot*ntw**0.5)
  unew2= (delt_pr/Mtw)*gacc_pr*psi*COS(chin)
  unew3=  delt_pr*(lc_pr/am_pr)*(uchinew**2)/rhotw

!  unew4=  Pcon*(pi*rpel**2*sigma_spitzer)/Lcon

  unewtot = unew1+unew2+unew3


  !Include new coriolis term  - Jun2004 - Fixed - sign on psi term  
  uchinew=uchiold - delt_pr*nudrag*(1+i2/(Utot*ntw**0.5))/Mtw*(uchiold-uexb) - (delt_pr/Mtw)*gacc_pr &
          *psi*SIN(chin) - delt_pr*(lc_pr/am_pr)*(uchinew*unew)/rhotw
  uchi1=-delt_pr*nudrag*(1+i2/(Utot*ntw**0.5))/Mtw*(uchiold-uexb)
  uchi2=(delt_pr/Mtw)*gacc_pr*psi*SIN(chin)
  uchi3=-delt_pr*(lc_pr/am_pr)*(uchinew*unew)/rhotw

  duchi = unew*Lc_pr*xintn_pr(na)/Lsloc

  Utot = (unew**2 + (uchinew-uexb)**2)**0.5


  IF(rhotw <= 0.0 ) THEN !HFS

    !? debug stuff
    if (l_debug_pr) then
        write(*,*) ' Blob has reached plasma center ', na
        !CALL PRL_OUT(1,nout)
        !CALL PRL_OUT(2)
!        write(15,160) tau_pr,psi,xintn_pr(na+1),unew,presbc_pr,rhotw
!        write(*,160) tau_pr,psi,xintn_pr(na+1),unew,presbc_pr,rhotw
!        write(15,*) ' Blob has reached plasma center ', na
!        close(15)     
    endif

    !Store final location of cells in rho space
    massloc(1:na)=rhotw  

    !Return normalized mass depos profile to caller in nea_r array
    CALL PRL_GETNORMDEP(n_r,pnum,massloc,massprof,frad,rho_r,nea_r)
    EXIT LOOP_J

  ENDIF

  IF(rhotw >= 1.01) THEN   ! LFS

    !? write stuff
    if (l_debug_pr) then
        write(*,*) ' Blob has exited plasma ', na,  etot
        !CALL PRL_OUT(1,nout)
        !CALL PRL_OUT(2)
!        write(15,160) tau_pr,psi,xintn_pr(na+1),unew,presbc_pr,rhotw  
!        write(*,160) tau_pr,psi,xintn_pr(na+1),unew,presbc_pr,rhotw
!        write(15,*) ' Blob has exited plasma ',na
!        close(15)
    endif

    !Store final location of cells in rho space
    massloc(1:na)=rhotw

    !Return normalized mass depos profile to caller in nea_r array
    CALL PRL_GETNORMDEP(n_r,pnum,massloc,massprof,frad,rho_r,nea_r)
    EXIT LOOP_J

  ENDIF

  !Save new values to old variables          
  rhotwo=rhotw
  uold=unew
  uchiold=uchinew
  chiold=chin

  !Get new local values from stored profiles
  DO i=1,nr_pr

    IF(REAL(i)/nr_pr < rhotw) THEN

      teinf_pr=te_prl(i)
      neinf_pr=ne_prl(i)
      qloc=q_prl(i)
      Lsloc=Ls(i)

    ENDIF

  ENDDO


!***
!teinf_pr = 5000.0
!neinf_pr = 1.0e14

  etot = 0.0
  presbc_pr=presbc0*teinf_pr/te_inf0     ! Assumes constant plasma density profile ---  needs improving.

  IF(l_tesig_pr) THEN

    !Calculate new sigma value - temp dependent (Parks FAX 2-Jul-2002)
    lambda_ef=1.35e10*teinf_pr/(n0_pr**0.5)
    lambda_ef0=1.35e10*te_inf0/(n0_pr**0.5)
    sigma_pr=sigma0_pr*((te_inf0/teinf_pr)**2)*LOG(lambda_ef)/LOG(lambda_ef0)
    IF(sigma_pr < 0.0) sigma_pr=sigma0_pr
    taubarn=na*taubar_pr/na_pr


 
    DO i=1,na

      uplus_pr(i)=sigma_pr*(taubarn+uplusint_pr(i))   
      uminus_pr(i)=sigma_pr*(taubarn-uplusint_pr(i))

      IF(uminus_pr(i) < 0.0) THEN
        !? debug stuff
        if (l_debug_pr) then
!           write(*,*) 'Tedepend uminus <0', j, i, uminus_pr(i)
		endif
        uminus_pr(i)=0.0001
      ENDIF

      IF(uplus_pr(i) < 0.0) THEN
        !? debug stuff
        if (l_debug_pr) then
!           write(*,*) 'Tedepend uplus <0', j, i, uminus_pr(i)
		endif
        uplus_pr(i)=0.0001
      ENDIF

      zplus_pr(i)=SQRT(uplus_pr(i))
      zminus_pr(i)=SQRT(uminus_pr(i))
      gplus_pr(i)=SQRT(zplus_pr(i))*PRL_BESSK1(SQRT(zplus_pr(i)))/4
      gminus_pr(i)=SQRT(zminus_pr(i))*PRL_BESSK1(SQRT(zminus_pr(i)))/4
      ghat_pr(i)=gplus_pr(i)+gminus_pr(i)
      gdiff_pr(i)=gplus_pr(i)-gminus_pr(i)
      etaplus_pr(i)=uplus_pr(i)*PRL_BESSK(2,zplus_pr(i))/2
      etaminus_pr(i)=uminus_pr(i)*PRL_BESSK(2,zminus_pr(i))/2
      etahat_pr(i)=etaplus_pr(i)+etaminus_pr(i)

      !Calculate new value of heating quotient (PoP Eqn 14 scaled by (Te/Te0)^3/2)
      qhatt=qhat_pr*(teinf_pr/te_inf0)**1.5
!***          qhatt=qhat_pr
      khat_pr(i)=sigma_pr*qhatt*ghat_pr(i)*hfact  !hfact is enhanced heating factor


! Energy in cloudlet - LRB   5/2008
	  energy_pr(i) = 1.5*n0_pr*denn_pr(i)*tcloud_pr*tempn_pr(i)*vol_pr(i)*1.602e-19      !  3/2nkT * Vol  (Joules) 
      etot = etot + energy_pr(i)

    ENDDO

  ENDIF

  !Store output into arrays for viewing
  IF((j == (jsave+1)*isave) .AND. &
     (jsave < 500)) THEN

    !Save into arrays
    jsave=jsave+1
    rho_arr(jsave,1:na)=denn_pr(1:na)        ! Density
    phi_arr(jsave,1:na)=tempn_pr(1:na)       ! Temp
    beta_arr(jsave,1:na)=presn_pr(1:na)      ! Press
    Q_arr(jsave,1:na)=khat_pr(1:na)          ! Heat Q
    vhat_arr(jsave,1:na)=vhato_pr(1:na)      ! Vel_axial
    Lc_arr(jsave,1:na)=xintn_pr(1:na)        ! Length Cloud
    dtw_arr(jsave,1:na)=deltatw(1:na)
    thtw_arr(jsave,1:na)=theta_pr(1:na)
    rtw_arr(jsave,1:na)=rtw(1:na)
    Rmtw_arr(jsave,1:na)=Rmtw(1:na)
    Mnum_arr(jsave,1:na)=Mnum(1:na)          ! Mach num
    u_arr(jsave)=unew                        ! V rad total
    u1_arr(jsave)=unew1                      ! V rad 1
    u1i2_arr(jsave)=unew1i2                  ! V rad 1 i2 factor
    u2_arr(jsave)=unew2                      ! V rad 2
    u3_arr(jsave)=unew3                      ! V rad 3
    ut_arr(jsave)=unewtot                    ! V rad 3
    uchi_arr(jsave)=uchinew
	chin_arr(jsave)=chin                    
    cos_arr(jsave)=cos(chin)                 ! cos(chin)
	psi_arr(jsave)=psi
    presNC(jsave)=presn_pr(na)               ! Press cloud end
    presbc_arr(jsave)=presbc_pr              ! Press plasma boundary condition
    tau_arr(jsave)=tau_pr                    ! tau
	rhotw_arr(jsave) = rhotw			     ! r/a cloud location
	Ac_arr(jsave) = 2*Lc_pr*xintn_pr(na)*2*z_pi*rcloud    ! Cloud area cm^2 - first 2 is cloud halves
	Gradn_arr(jsave) = (n0_pr*denn_pr(50) - neinf_pr)/rcloud  ! # cm^-4
	Fluxn_arr(jsave) = Gradn_arr(jsave)*10000  ! #/ (cm^2 sec)  - assume 1 m^2/sec diffusion coefficient
	Ndiff_arr(jsave) = Fluxn_arr(jsave)*Ac_arr(jsave)   ! #/sec   - diffusive loss

    n1_arr(jsave)=denn_pr(1)         ! Density 1
    n2_arr(jsave)=denn_pr(100)       ! Density 100
    n3_arr(jsave)=denn_pr(na-1)      ! Density 199
    T1_arr(jsave)=tempn_pr(1)        ! Temp 1
    T2_arr(jsave)=tempn_pr(100)      ! Temp 100
    T3_arr(jsave)=tempn_pr(na-1)     ! Temp 199
    LC3_arr(jsave)=xintn_pr(na-1)    ! Lc 199


    i2_arr(jsave) = i2               ! i2 sc current  Jul2009
    i1_arr(jsave) = i_1              ! i1 current  Jul2009
    nu1_arr(jsave) = nu1
    nu2_arr(jsave) = nu2

    if (jsave==20) then
	    Pr1=tempn_pr(1:na)
	    Lc1=xintn_pr(1:na)
    endif
    if (jsave==50) then
	    Pr2=tempn_pr(1:na)
	    Lc2=xintn_pr(1:na)
    endif
    if (jsave>=200) then
	    Pr3=tempn_pr(1:na)
	    Lc3=xintn_pr(1:na)

        open(unit=15, file = 'pressout.dat')

        DO ii=1,na 
		    write(15,160) Lc1(ii),Pr1(ii),Lc2(ii),Pr2(ii),Lc3(ii),Pr3(ii)
!    160     format(e12.5,1x,e12.5,1x,e12.5,1x,e12.5,1x,e12.5,1x,e12.5)
		
		enddo           

        close(15)

    endif


    





  ENDIF        

  !Write output and store at appropriate times
  IF((icycle_pr/ocycle_pr)*ocycle_pr == icycle_pr) THEN

    !? write stuff
!    IF(l_out_pr) THEN    
!      CALL PRL_OUT(1,nout)  !Write axial profile output data
!    ENDIF
         
!    WRITE(15,161) tau_pr,psi,xintn_pr(na+1),unew,presbc_pr,rhotw,tempo_pr(1)
!    WRITE(*,161) tau_pr,psi,xintn_pr(na+1),unew,presbc_pr,rhotw,tempo_pr(1)

    !Save data to global common variables  (LRB)
    ptime_pr(k)=tau_pr
    vcld_pr(k)=unew
    psip_pr(k)=psi
    k=k+1

  ENDIF

ENDDO LOOP_J

!CLOSE(15)dummy = 0
!IF(l_out_pr) THEN
!    CALL PRL_OUT(2) ! Close output file
!ENDIF
      
    k = k+1    !dummy test location

END SUBROUTINE PRL




SUBROUTINE PRL_LAGSETUP(nex,tex)
!-------------------------------------------------------------------------------
!PRL_LAGSETUP sets up cells for 1D Lagrangian 
!
!References:
!  L.R.Baylor, Fortran version 14-Jun-2000
!              Ported to subroutine for PRLag 4-Apr-2001
!              Modified for momentum flux terms 24-Jun-2002
!              Modified for use inside PRL_sub  10-Jun-2003
!              Modified for F90 version of PRL  10-Sep-2004
!  W.A.Houlberg, L.R.Baylor, Standardization for NTCC Module Library 1/2005
!
!Comments:
!  Finite Difference Code
!  Pressure relaxation model (Parks and Sessions)
!  Duct Profile w/ time-varying BC(pres_ambient)
!-------------------------------------------------------------------------------

!Declaration of input variables
REAL, INTENT(IN) :: &
  tex,                 & !electron temperature [keV]
  nex                    !electron density [/cm**3]

!-------------------------------------------------------------------------------
!Declaration of local variables
INTEGER :: &
  i

REAL(KIND=rspec) :: &
  Bpel,cs1,Itotal,dI,dxi,gam,lambda_eb,lambda_ef,loglambda,M0,mion,muE,R,Rout, &
  theta0,xlo,xhi,xint(na_pr+1),xcen(na_pr)

!-------------------------------------------------------------------------------
!Initialization
!-------------------------------------------------------------------------------
muE = 0.64
!*** M0 = 0.8
M0 = 0.8
wpel_pr = 2
tcloud_pr = 2

!-------------------------------------------------------------------------------
!Calculate Sigma0 and n0 from Parks, et al., 2001
!-------------------------------------------------------------------------------
R = rm_pr+COS(chi0_pr)*am_pr*rhop_pr  !major radius of pellet cloud formation
rin_pr = rm_pr-am_pr
Rout = rm_pr+am_pr
Bpel = bt_pr*(rm_pr/R)                !magnetic field at pellet cloud
gam = 5.0/3.0                         !ratio of specific heats Cp/Cv
mion = 3.3e-27                        !mass of ion in kg
lambda_eb = 2.52*tex/7.5              !7.5eV mean electronic excitation energy of H

cs1=SQRT(2*(5.0/3.0)*tcloud_pr*1.6e-19/mion)*100  !convert to cm/s
cs0_pr=981761.0*SQRT(2*(5.0/3.0)*tcloud_pr/wpel_pr)

kappac_pr=(2.17226e4*Tex**(1.0/6.0)*SQRT(4*gam*tcloud_pr+15.8)) &
          /((nex*rpel_pr*LOG(lambda_eb))**(1.0/3.0)*wpel_pr**(1.0/6.0)) 
!kappac_pr=2.17226e4*teinf_pr**(1.0/6.0)/((wpel_pr**(1.0/6.0)) &
!          *(neinf_pr*rpel_pr*LOG(lambda_eb))**(1.0/3.0))
!          *SQRT(6.67*tcloud_pr+15.8)


lc_pr=SQRT(kappac_pr*rpel_pr*R)
n0_pr=(1.304e15*Tex**(11.0/6.0)/(M0*kappac_pr**2*cs0_pr &
      *LOG(lambda_eb)**(2.0/3.0)))*(nex/(wpel_pr*rpel_pr**2))**(1.0/3.0)
lambda_ef=1.35e10*tex/(n0_pr**0.5)
loglambda=LOG(lambda_ef)

!Changed from 652 on 4-Mar-2003 from Parks memo
sigma0_pr=679*R**0.5/(M0*kappac_pr**1.5*wpel_pr**(1.0/3.0)*cs0_pr) &
          *(nex**2/(Tex*rpel_pr))**(1.0/6.0) &
          *LOG(lambda_ef)/(LOG(lambda_eb))**(2.0/3.0)

!Store value
sigma_pr=sigma0_pr

pratio_pr=(406.0*cs0_pr*(Tex**(5.0/6.0))/(M0*kappac_pr**2)) &
          *(wpel_pr/(nex*rpel_pr*LOG(lambda_eb)))**(2.0/3.0)


!****  Test reduced pratio  LRB Jan2010
pratio_pr = pratio_pr/2.

!****
!sigma_pr = 5.0
!sigma0_pr = 5.0
!pratio_pr = 5.0


!Calculate value of heating quotient (PoP Eqn 15)
theta0=((1.0+gam)**2)*(M0**2)/(1.0+gam*M0**2)**2
qhat_pr=2*gam**(3.0/2.0)*(((1.0+gam*M0**2)**2)/(M0*(1.0+gam)**2)) &
        *(1.0+(15.8*theta0)/(4*gam*tcloud_pr))

!***
!qhat_pr = 5.0

!-------------------------------------------------------------------------------
!Display  sigma0 and n0 values
!-------------------------------------------------------------------------------
!? debug stuff
if (l_debug_pr) then
    write(*,*) ' PRLagrangian Code  -  setup parameters '
    write(*,*) '   '
    write(*,*) 'Sigma0 and Qhat and Kappac: ', sigma0_pr, qhat_pr ,kappac_pr
    write(*,*) 'Beta0/BetaPl: ', pratio_pr
    write(*,*) 'loglambda_eb: ', loglambda
    write(*,*) ' '
    write(*,150) ' Lc (cm): ', lc_pr
    write(*,150) ' c0 (cm/s): ', cs0_pr/1.265        ! 1.265 = (5./3)^0.5
    write(*,150) ' Lc/c0 (s): ', lc_pr/(cs0_pr/1.265)
    write(*,*) ' '

  150 format(1x,a10,e12.4)
endif

!-------------------------------------------------------------------------------
!Calculate Itotal - integrate function PRL_LFUNC from 0.0 to 1.0    
!-------------------------------------------------------------------------------
!***Itotal=PRL_INTEGRATEF(3,0.0,1.0)
Itotal=PRL_INTEGRATEF(3,0.0,1.0)
taubar_pr=Itotal
dI=Itotal/na_pr
dxi=dI
xlo=0
xhi=0

!Calculate dimensionless cloudlet inertia and acceleration parameter
!Modified from 0.229 11/14/2003 LRB
inert_pr=0.231*taubar_pr*(tex**(11.0/6.0)) &
         *((wpel_pr/nex)**(1.0/6.0))/((rpel_pr**(2.0/3.0))*Bpel &
         *(M0*kappac_pr**2)*(LOG(lambda_eb))**(2.0/3.0))

gacc_pr=(2/taubar_pr)*SQRT(kappac_pr*rpel_pr/R)
afact_pr=lc_pr/rin_pr

!Loop to determine xhi
DO i=2,na_pr+1

  xhi=xhi+dxi
  xint(i)=xhi
  xlo=xhi

ENDDO

xint(na_pr+1)=1.0

!Loop to get xcentered values
DO i=1,na_pr

  xcen(i)=(xint(i+1)+xint(i))/2

ENDDO

!Save output for PRLag LRB
xinto_pr(1:na_pr+1) = xint(1:na_pr+1)
xceno_pr(1:na_pr) = xcen(1:na_pr)


!-------------------------------------------------------------------------------
!Main loop -  calculation of incident heat flux
!-------------------------------------------------------------------------------
DO i=1,na_pr


!***  uplusint_pr(i) = PRL_INTEGRATEF(3,0.0,xcen(i))
  uplusint_pr(i) = PRL_INTEGRATEF(3,0.0,xcen(i))
  uplus_pr(i) = sigma0_pr*(taubar_pr+uplusint_pr(i))   
  uminus_pr(i) = sigma0_pr*(taubar_pr-uplusint_pr(i))


  write(*,*) i, uplus_pr(i), uminus_pr(i)
  if (uminus_pr(i) < 0.0) uminus_pr(i) = 0.0

  zplus_pr(i) = SQRT(uplus_pr(i))
  zminus_pr(i) = SQRT(uminus_pr(i))

  write(*,*)  i, zplus_pr(i), zminus_pr(i)

  gplus_pr(i) = SQRT(zplus_pr(i))*PRL_BESSK1(SQRT(zplus_pr(i)))/4
    if (zminus_pr(i) <= 0.0) then
	   gminus_pr(i) = 0.0
	else
       gminus_pr(i) = SQRT(zminus_pr(i))*PRL_BESSK1(SQRT(zminus_pr(i)))/4
	endif
  ghat_pr(i) = gplus_pr(i)+gminus_pr(i)
  gdiff_pr(i) = gplus_pr(i)-gminus_pr(i)

  etaplus_pr(i) = uplus_pr(i)*PRL_BESSK(2,zplus_pr(i))/2
    if (zminus_pr(i) <= 0.0) then
	   etaminus_pr(i) = 0.0
	else
       etaminus_pr(i) = uminus_pr(i)*PRL_BESSK(2,zminus_pr(i))/2
	endif
  etahat_pr(i) = etaplus_pr(i)+etaminus_pr(i)

ENDDO
       
END SUBROUTINE PRL_LAGSETUP


      
SUBROUTINE PRL_INIT
!-------------------------------------------------------------------------------
!PRL_INIT initializes longitudinal profiles along the cloud axis
!
!References:
!  W.A.Houlberg, L.R.Baylor, Standardization for NTCC Module Library 1/2005
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
!Declaration of local variables
INTEGER :: i
      
REAL(KIND=dpspec) :: &
  delxint,denbc

REAL(KIND=rspec) :: &
  rhol,rhor,x

!-------------------------------------------------------------------------------
!Initialize cell-centered profiles
xinti_pr(1:na_pr+1)=xinto_pr(1:na_pr+1)

!Old velocity
vhato_pr(1:na_pr+1)=0

!Old temperature
tempo_pr(1:na_pr+1)=1

!Initial temperature
tempi_pr(1:na_pr+1)=tempo_pr(1:na_pr+1)*tcloud_pr

rhol=1
rhor=1/pratio_pr
      
!Initial density
DO i=1, na_pr+1

  x=xinto_pr(i)
  deni_pr(i)=(rhol+rhor)/2+(rhol-rhor)*TANH(10*(1.0-x))/2
  deni_pr(i)=deni_pr(i)-((1.0-deni_pr(i))/(1.0-(rhor+rhol)/2)) &
                        *((1.0-rhor)-(1.0-(rhor+rhol)/2))

!***  Dummy profile for test
!***   deni_pr(i) = rhol - x * (rhol-rhor)
ENDDO

!Old density                                     
deno_pr(1:na_pr+1)=deni_pr(1:na_pr+1)

!Old pressure
preso_pr(1:na_pr+1)=deno_pr(1:na_pr+1)*tempo_pr(1:na_pr+1)

!Boundary conditions
denbc=deno_pr(na_pr+1)
presbc_pr=preso_pr(na_pr+1)

!? debug stuff
if (l_debug_pr) then
    write(*,*) '  '
    write(*,*) 'press_ratio= ', preso_pr(1)/preso_pr(na_pr+1)
    write(*,*) 'I, g = ', inert_pr, gacc_pr
    write(*,*) '  '
endif
        
!Calculate conserved mass variable, i.e. rho * delxint:
DO i=1,na_pr

  delxint=xinto_pr(i+1)-xinto_pr(i)
  dencon_pr(i)=deno_pr(i)*delxint

ENDDO    
      
dencon_pr(na_pr+1)=denbc*(xinto_pr(na_pr+1)-xinto_pr(na_pr))
pratio_pr=preso_pr(1)/preso_pr(na_pr+1)
   
END SUBROUTINE PRL_INIT

 

SUBROUTINE PRL_GETNORMDEP(n_r,pnum,massloc,massprof,frad,rho_r,nea_r)
!-------------------------------------------------------------------------------
!PRL_GETNORMDEP calculates the normalized deposition result
!
!References:
!  W.A.Houlberg, L.R.Baylor, Standardization for NTCC Module Library 1/2005
!-------------------------------------------------------------------------------

!Declaration of input variables
INTEGER, INTENT(IN) :: &
  n_r                   !

REAL(KIND=rspec), INTENT(IN) :: &
  pnum,                & !
  rho_r(:),            & !
  massloc(:)

!Declaration of input/output variables
REAL(KIND=rspec), INTENT(INOUT) :: &
  massprof(:)            !

!Declaration of output variables
REAL(KIND=rspec), INTENT(OUT) :: &
  frad(:),             & !
  nea_r(:)                 !

!-------------------------------------------------------------------------------
!Declaration of local variables
CHARACTER(len=24) :: &
  fname

INTEGER :: &
  i,jj,jjmin,jsp,l

REAL(KIND=rspec) :: &
  fmtot,masstot

!-------------------------------------------------------------------------------
! Initialization
!-------------------------------------------------------------------------------
fmtot=0

!-------------------------------------------------------------------------------
!Determine normalized mass deposition profile
!-------------------------------------------------------------------------------
DO i=1,na_pr

  l=massloc(i)*nr_pr      !Radial cell where mass is deposited

  IF(l > 0 .AND. &
     l <= nr_pr) THEN

    massprof(l)=massprof(l)+1.0
    fmtot=fmtot+1.0

  ENDIF

ENDDO

!-------------------------------------------------------------------------------
!Write output to nprof file
!-------------------------------------------------------------------------------
!? debug stuff
if (l_debug_pr) then

    jsp=INDEX(cidrun_pr,' ')
    fname='PRL_'//cidrun_pr(1:jsp-1)//'_nprof'//'.dat'
    OPEN(UNIT=12, FILE=fname)

    DO i=1,nr_pr
        WRITE(12,*) i,massprof(i)
    ENDDO

    CLOSE(12)

endif

!-------------------------------------------------------------------------------
!Normalized mass deposition profile
!-------------------------------------------------------------------------------
jjmin=1

DO i=1,n_r !Go through resulting massprof array

  jj=nr_pr

  DO WHILE(jj > jjmin)

    frad(jj)=REAL(jj,rspec)/nr_pr

    IF(frad(jj) >= rho_r(i)) THEN

      nea_r(i)=massprof(jj)
      jj=jj-1
      IF(frad(jj) == rho_r(i)) jjmin=jj+1

    ELSE

      nea_r(i)=massprof(jj)+((massprof(jj)-massprof(jj+1))/  &
               (frad(jj)-frad(jj+1)))*(frad(jj)-rho_r(i))
      jjmin=jj+1
      jj=jj-1

    ENDIF

  ENDDO

ENDDO

!Scale the nea_r output array
masstot=SUM(nea_r(1:n_r))
IF(masstot > 0.0) nea_r(1:n_r)=nea_r(1:n_r)*(pnum/masstot)*fmtot/na_pr

END SUBROUTINE PRL_GETNORMDEP  

SUBROUTINE PRL_DTIME(na)
!-------------------------------------------------------------------------------
!PRL_DTIME calculates the maximum time step via Courant condition
!
!References:
!  W.A.Houlberg, L.R.Baylor, Standardization for NTCC Module Library 1/2005
!-------------------------------------------------------------------------------

!Declaration of input variables
INTEGER, INTENT(IN) :: &
  na

!-------------------------------------------------------------------------------
!Declaration of local variables
INTEGER :: &
  i
      
REAL(KIND=dpspec) :: &
  cspd,cspd2,rdtime,rdx,rdtx,vcavg
      
!-------------------------------------------------------------------------------
!Initialization
rdtime=0

!Time step calculation
DO i=2,na-1

  vcavg=(vhato_pr(i)+vhato_pr(i-1))/2          
  cspd2=2*ABS(tempo_pr(i))
  cspd=SQRT(cspd2)
  rdx=1/(xinto_pr(i)-xinto_pr(i-1))        
  rdtx=rdx*(ABS(vcavg)+cspd)        
  rdtime=MAX(rdtime,rdtx)

ENDDO
      
delt_pr=0.1/rdtime
   
END SUBROUTINE PRL_DTIME 




                  
SUBROUTINE PRL_OUT(oflag, nout)
!-------------------------------------------------------------------------------
!PRL_OUT writes data to a output file for debugging
!
!References:
!  W.A.Houlberg, L.R.Baylor, Standardization for NTCC Module Library 1/2005
!
! 19-Aug-2005   Changed to netcdf output files.   LRB
!
!-------------------------------------------------------------------------------
!-----------------------------------------------------------------------
! NOTE:      
! This subroutine requires netcdf v. 3.5 or above w. Robert Pincus f90
! interface
!-----------------------------------------------------------------------

USE netcdf
      
!Declaration of input variables
INTEGER, INTENT(IN) :: &
     nout,   &          ! msg output unit 
     oflag              ! option for file operation

!-------------------------------------------------------------------------------
!Declaration of local variables
INTEGER :: &
  i, istat, ierr
                     
REAL(KIND=dpspec) :: &
  x

CHARACTER(len=26) :: &
        cnetcdf                         ! output filename

INTEGER :: &
       NCID,      &                      ! NETCDF ID number
       VAR_ID(25), &                     ! List of variable id's
       VAR_ID_RHO,  &                    ! rho grid id  
       VAR_ID_TIME,  &                   ! time grid id'
       ID_RHO, &
       ID_TIME, &
       ID_DEN, &
       ID_TEMP, &
       ID_PRES, &
       ID_VHAT, &
       ID_LCARR, &
       ID_MACH, &
       ID_VRAD, &
       ID_PRESN, &
       ID_PRESBC, &
       ID_TAU, &
       ID_RMIN, &
	   ID_LC, &
	   ID_CS0

!-------------------------------------------------------------------------------
!Perform requested file interaction
!-------------------------------------------------------------------------------

IF(oflag==1) THEN

!-----------------------------------------------------------------------
! Open the NETCDF file, overwriting any existing file with the same name
!-----------------------------------------------------------------------
      cnetcdf = 'PRL_'//trim(cidrun_pr)//'_profiles.ncd'
      ierr=nf90_create(trim(cnetcdf),nf90_clobber,ncid)

      if (ierr/=nf90_noerr) then
        write(nout,*)
        write(nout,*) ' Warning: NETCDF file error:'
        write(nout,*)   trim(nf90_strerror(ierr))
        return
      else
        write(nout,*)
        write(nout,*) 'netCDF v.'//trim(NF90_INQ_LIBVERS())//   &
                     ' output file :',trim(cnetcdf)
      end if


!-----------------------------------------------------------------------
! Define the grid to store things on
!-----------------------------------------------------------------------

! Define the base grid - rho (nr)        
      ! define dimension  
      istat=nf90_def_dim(ncid,'rho',na_pr,id_rho)
      call check_err(nout,istat)

      istat=nf90_def_dim(ncid,'time',ns_pr,id_time)
      call check_err(nout,istat)      

      ! Define variable id 
      istat=nf90_def_var(ncid,'rho',NF90_FLOAT,(/id_rho/),var_id_rho)
      call check_err(nout,istat)

      istat=nf90_def_var(ncid,'time',NF90_FLOAT,(/id_time/),var_id_time)
      call check_err(nout,istat)

      ! Set attributes - units
      istat=nf90_put_att(ncid,var_id_rho,'units','')
      call check_err(nout,istat)
      ! Set attributes - define long names
      istat=nf90_put_att(ncid,var_id_rho,'long_name','sqrt toroidal flux')
      call check_err(nout,istat)

      ! Set attributes - units
      istat=nf90_put_att(ncid,var_id_time,'units','s')
      call check_err(nout,istat)
      ! Set attributes - define long names
      istat=nf90_put_att(ncid,var_id_time,'long_name','time')              
      call check_err(nout,istat)

      ! Define variable ids 
      istat=nf90_def_var(ncid,'den_pr',NF90_FLOAT,(/id_time, id_rho/),id_den)
      call check_err(nout,istat)
      istat=nf90_def_var(ncid,'temp_pr',NF90_FLOAT,(/id_time, id_rho/),id_temp)
      call check_err(nout,istat)
      istat=nf90_def_var(ncid,'pres_pr',NF90_FLOAT,(/id_time, id_rho/),id_pres)
      call check_err(nout,istat)
      istat=nf90_def_var(ncid,'vhat_pr',NF90_FLOAT,(/id_time, id_rho/),id_vhat)
      call check_err(nout,istat)
      istat=nf90_def_var(ncid,'lc_pr',NF90_FLOAT,(/id_time, id_rho/),id_lcarr)
      call check_err(nout,istat)
      istat=nf90_def_var(ncid,'mach_pr',NF90_FLOAT,(/id_time, id_rho/),id_mach)
      call check_err(nout,istat)
      istat=nf90_def_var(ncid,'vrad',NF90_FLOAT,(/id_time/),id_vrad)
      call check_err(nout,istat)
      istat=nf90_def_var(ncid,'presn',NF90_FLOAT,(/id_time/),id_presn)
      call check_err(nout,istat)
      istat=nf90_def_var(ncid,'presbc',NF90_FLOAT,(/id_time/),id_presbc)
      call check_err(nout,istat)
      istat=nf90_def_var(ncid,'tau',NF90_FLOAT,(/id_time/),id_tau)
      call check_err(nout,istat)
      istat=nf90_def_var(ncid,'rmin',NF90_FLOAT,(/id_time/),id_rmin)
      call check_err(nout,istat)

      ! Define scalar variables
      istat=nf90_def_var(ncid,'Lc',NF90_FLOAT,id_lc)
      call check_err(nout,istat)
      istat=nf90_def_var(ncid,'cs0',NF90_FLOAT,id_cs0)
      call check_err(nout,istat)

      ! Set attributes - units
      istat=nf90_put_att(ncid,id_lc,'units','cm')
      call check_err(nout,istat)
      ! Set attributes - define long names
      istat=nf90_put_att(ncid,id_lc,'long_name','Cloud Length')
      call check_err(nout,istat)

      ! Set attributes - units
      istat=nf90_put_att(ncid,id_cs0,'units','cm/s')
      call check_err(nout,istat)
      ! Set attributes - define long names
      istat=nf90_put_att(ncid,id_cs0,'long_name','Sound Speed')              
      call check_err(nout,istat)

!-----------------------------------------------------------------------
!Close the definition section of the file 
!-----------------------------------------------------------------------

      istat=nf90_ENDDEF(ncid)

      call check_err(nout,istat)

!-----------------------------------------------------------------------
! Define the independent rho grid 
!-----------------------------------------------------------------------
      istat=nf90_PUT_VAR(ncid,var_id_rho,xintn_pr(1:na_pr))
      call check_err(nout,istat)

      ! The above statements only needs to be done once for every file
!-----------------------------------------------------------------------
! Define the independent time variables
!-----------------------------------------------------------------------
      istat=nf90_PUT_VAR(ncid,var_id_time,tau_arr)
      call check_err(nout,istat)
!-----------------------------------------------------------------------
! Write data to the file
!-----------------------------------------------------------------------

      istat=nf90_PUT_VAR(ncid,id_den,rho_arr)
      call check_err(nout,istat)
      istat=nf90_PUT_VAR(ncid,id_temp,phi_arr)
      call check_err(nout,istat)
      istat=nf90_PUT_VAR(ncid,id_pres,beta_arr)
      call check_err(nout,istat)
      istat=nf90_PUT_VAR(ncid,id_vhat,vhat_arr)
      call check_err(nout,istat)
      istat=nf90_PUT_VAR(ncid,id_lcarr,lc_arr)
      call check_err(nout,istat)
      istat=nf90_PUT_VAR(ncid,id_mach,Mnum_arr)
      call check_err(nout,istat)
      istat=nf90_PUT_VAR(ncid,id_vrad,u_arr)
      call check_err(nout,istat)
      istat=nf90_PUT_VAR(ncid,id_presn,presNC)
      call check_err(nout,istat)
      istat=nf90_PUT_VAR(ncid,id_presbc,presbc_arr)
      call check_err(nout,istat)
      istat=nf90_PUT_VAR(ncid,id_tau,tau_arr)
      call check_err(nout,istat)
      istat=nf90_PUT_VAR(ncid,id_rmin,rhotw_arr)
      call check_err(nout,istat)

      istat=nf90_PUT_VAR(ncid,id_lc,lc_pr)
      call check_err(nout,istat)
      istat=nf90_PUT_VAR(ncid,id_cs0,cs0_pr/gammasqrt)
      call check_err(nout,istat)

!-----------------------------------------------------------------------
! Close the file and return
!-----------------------------------------------------------------------      
      
      call check_err(nout,nf90_close(ncid))


endif

      CONTAINS

      subroutine check_err(ilun,ierr)
      integer  ::  ilun,ierr
      if (ierr/=nf90_noerr) then
        write(ilun,*)
        write(ilun,*)'   Status: '//trim(nf90_strerror(ierr))
      end if  
      end SUBROUTINE check_err
      















!  OPEN(UNIT=11, &
!       FILE='PRL_'//cidrun_pr//'_profiles'//'.dat')

!ELSEIF(oflag==1) THEN

!  write(11,*) ' '  ! Blank line for delimiting                           
!  write(11,*) tau_pr
!  write(11,*) 'x              n              T           p  ', &
!              '         M'

!  DO i=1,na_pr,5   ! Every 5th pt
!    x=xceno_pr(i)  ! Display cell centered positions
!    write(11,170) x,deno_pr(i),tempo_pr(i),preso_pr(i),vhato_pr(i)
!    170 format(e12.4,1x,e12.4,1x,e12.4,1x,e12.4,1x,e12.4)
!  ENDDO

!ELSEIF(oflag==2) THEN

!  CLOSE(11)

!ENDIF

END SUBROUTINE PRL_OUT 


              
 
FUNCTION PRL_INTEGRATEF(k_func,a,b)
!-------------------------------------------------------------------------------
!PRL_INTEGRATEF integrate function PRL_LFUNC from a to b using trapezoidal rule
!
!References:
!  W.A.Houlberg, L.R.Baylor, Standardization for NTCC Module Library 1/2005
!
!Comments:
!  See Numerical Recipes for algorithm used
!-------------------------------------------------------------------------------

!Declaration of input varaibles
INTEGER, INTENT(IN) :: &
  k_func                 !option for cloud density profile function [-]

REAL(KIND=rspec), INTENT(IN) :: &
  a,                   & !lower bound of integral [arb]
  b                      !upper bound of integral [arb]

!Declaration of output variables
REAL(KIND=rspec) :: &
  PRL_INTEGRATEF

!-------------------------------------------------------------------------------
!Declaration of local variables
INTEGER, PARAMETER :: &
  jmax=20

REAL(KIND=rspec), PARAMETER :: &
  eps=1.0e-6

INTEGER :: &
  j

REAL(KIND=rspec) :: &
  s,olds

!-------------------------------------------------------------------------------
!Initialization
olds=-1.0e30

!Return integral of function func from a to b
!Integration is performed to accuracy eps by at trapeziodal rule using at most
!  2^(jmax-1) steps
DO j=1,jmax

  CALL PRL_TRAPZD(k_func,a,b,j,s)
  PRL_INTEGRATEF=s
  IF(ABS(s-olds) < eps*ABS(olds)) EXIT
  olds=s

ENDDO

END FUNCTION PRL_INTEGRATEF




SUBROUTINE PRL_TRAPZD(k_func,a,b,n,s)
!-------------------------------------------------------------------------------
!PRL_BESSI0 gives the value of the I0 Bessel function
!
!References:
!  W.A.Houlberg, L.R.Baylor, Standardization for NTCC Module Library 1/2005
!  See Abramowitz and Stegun, Chapter 9, Dover 1964.
!-------------------------------------------------------------------------------

!Declaration of input variables
INTEGER, INTENT(IN) :: &
  k_func,              & !option for cloud density profile function [-]
  n                      !number of intervals

REAL(KIND=rspec), INTENT(IN) :: &
  a,                   & !lower bound of integral [arb]
  b                      !upper bound of integral [arb]

!Declaration of input/output variables
REAL, INTENT(INOUT) :: &
  s                      !integral veraged with previous value [arb]

!-------------------------------------------------------------------------------
!Declaration of local variables
INTEGER, SAVE :: &
  it

INTEGER :: &
  j

REAL(KIND=rspec) :: &
  x,sum,del, tnm

!-------------------------------------------------------------------------------
IF(n == 1) THEN

  s=(b-a)*(PRL_LFUNC(k_func,a)+PRL_LFUNC(k_func,b))/2
  it=1

ELSE

  tnm = it

  del = (b-a)/it
  x = a + del/2
  sum = 0.0

  DO j=1,it

    sum=sum+PRL_LFUNC(k_func,x)
    x=x+del

  ENDDO

  s=(s+(b-a)*sum/tnm)/2
  it=2*it

ENDIF

END SUBROUTINE PRL_TRAPZD
 
FUNCTION PRL_LFUNC(k_func,x)
!-------------------------------------------------------------------------------
!PRL_LFUNC provides the cloud axial density profile
!
!References:
!  W.A.Houlberg, L.R.Baylor, Standardization for NTCC Module Library 1/2005
!-------------------------------------------------------------------------------

!Declaration of input variables
INTEGER, INTENT(IN) :: &
  k_func                 !option for cloud density profile function [-]

REAL(KIND=rspec), INTENT(IN) :: &
  x                      !abscissa where density function is evaluated [arb]

!Declaration of output variables
REAL(KIND=rspec) :: &
  PRL_LFUNC

!-------------------------------------------------------------------------------
!Declaration of local variables
REAL(KIND=rspec) :: &
  a1,a2,rhor,rhol

!-------------------------------------------------------------------------------
!Set function value 
IF(k_func == 1) THEN

  !Constant cloud axial density profile case
  PRL_LFUNC=1

ELSEIF(k_func == 2) THEN

  !Variable cloud axial density profile case
  !**The profile factors are given temporary arbitrary values here
  a1=1
  a2=1
  PRL_LFUNC=(2/z_pi)*ATAN(EXP((-1/a2)*(x-a1)))

ELSEIF(k_func == 3) THEN

  !Variable cloud axial density profile case
  rhol=1
  rhor=1/pratio_pr
  PRL_LFUNC=(rhor+rhol)/2+(rhol-rhor)/2*tanh(10*(1.0-x))
  PRL_LFUNC=PRL_LFUNC-((1-PRL_LFUNC)/(1-(rhor+rhol)/2)) &
                      *((1-rhor)-(1-(rhor+rhol)/2))

ENDIF

END FUNCTION PRL_LFUNC

FUNCTION PRL_BESSI0(x)
!-------------------------------------------------------------------------------
!PRL_BESSI0 gives the value of the I0 Bessel function
!
!References:
!  W.A.Houlberg, L.R.Baylor, Standardization for NTCC Module Library 1/2005
!  See Abramowitz and Stegun, Chapter 9, Dover 1964.
!-------------------------------------------------------------------------------

!Declaration of input variables
REAL(KIND=rspec), INTENT(IN) :: &
  x

!Declaration of output variables
REAL(KIND=rspec) :: &
  PRL_BESSI0

!-------------------------------------------------------------------------------
!Declaration of local variables
REAL(KIND=dpspec), PARAMETER :: &
  p1= 1.0,          p2= 3.5156229,    p3= 3.0899424,    p4= 1.2067492,    &
  p5= 0.2659732,    p6= 0.360768e-1,  p7= 0.45813e-2,                     &
  q1= 0.39894228e0, q2= 0.1328592e-1, q3= 0.225319e-2,  q4=-0.157565e-2,  &
  q5= 0.916281e-2,  q6=-0.2057706e-1, q7= 0.2635537e-1, q8=-0.1647633e-1, &
  q9= 0.392377e-2

REAL(KIND=dpspec) :: &
  y

REAL(KIND=rspec) :: &
  ax

!-------------------------------------------------------------------------------
!Set function value
IF(ABS(x) < 3.75) THEN

  y=(x/3.75)**2
  PRL_BESSI0=p1+y*(p2+y*(p3+y*(p4+y*(p5+y*(p6+y*p7)))))

ELSE

  ax=ABS(x)
  y=3.75/ax
  PRL_BESSI0=(EXP(ax)/SQRT(ax)) &
             *(q1+y*(q2+y*(q3+y*(q4+y*(q5+y*(q6+y*(q7+y*(q8+y*q9))))))))

ENDIF

END FUNCTION PRL_BESSI0
  
FUNCTION PRL_BESSI1(x)
!-------------------------------------------------------------------------------
!PRL_BESSI1 gives the value of the I1 Bessel function
!
!References:
!  W.A.Houlberg, L.R.Baylor, Standardization for NTCC Module Library 1/2005
!  See Abramowitz and Stegun, Chapter 9, Dover 1964.
!-------------------------------------------------------------------------------

!Declaration of input variables
REAL(KIND=rspec), INTENT(IN) :: &
  x

!Declaration of output variables
REAL(KIND=rspec) :: &
  PRL_BESSI1

!-------------------------------------------------------------------------------
!Declaration of local variables
REAL(KIND=dpspec), PARAMETER :: &
  p1= 0.5,          p2= 0.87890594,   p3= 0.51498869,   p4= 0.15084934,   &
  p5= 0.2658733e-1, p6= 0.301532e-2,  p7= 0.32411e-3,                     &
  q1= 0.39894228,   q2=-0.3988024e-1, q3=-0.362018e-2,  q4= 0.163801e-2,  &
  q5=-0.1031555e-1, q6= 0.2282967e-1, q7=-0.2895312e-1, q8=0.1787654e-1,  &
  q9=-0.420059e-2

REAL(KIND=dpspec) :: &
  y

REAL(KIND=rspec) :: &
  ax

!-------------------------------------------------------------------------------
!Set function value
IF(ABS(x) < 3.75) THEN

  y=(x/3.75)**2
  PRL_BESSI1=x*(p1+y*(p2+y*(p3+y*(p4+y*(p5+y*(p6+y*p7))))))

ELSE

  ax=ABS(x)
  y=3.75/ax
  PRL_BESSI1=(EXP(ax)/SQRT(ax))   &
             *(q1+y*(q2+y*(q3+y*(q4+y*(q5+y*(q6+y*(q7+y*(q8+y*q9))))))))
  IF(x < 0.0) PRL_BESSI1=-PRL_BESSI1

ENDIF

END FUNCTION PRL_BESSI1

FUNCTION PRL_BESSK(n,x)
!-------------------------------------------------------------------------------
!PRL_BESSK gives the value of the K Bessel function
!
!References:
!  W.A.Houlberg, L.R.Baylor, Standardization for NTCC Module Library 1/2005
!  See Abramowitz and Stegun, Chapter 9, Dover 1964.
!
!Comments:
!  USES PRL_BESSK0 and PRL_BESSK1
!-------------------------------------------------------------------------------

!Declaration of input variables
INTEGER, INTENT(IN) :: &
  n

REAL(KIND=rspec), INTENT(IN) :: &
  x

!Declaration of output variables
REAL(KIND=rspec) :: &
  PRL_BESSK

!-------------------------------------------------------------------------------
!Declaration of local variables
INTEGER :: &
  j

REAL(KIND=rspec) :: &
  bk,bkm,bkp,tox

!-------------------------------------------------------------------------------
!Set function value
tox=2/x
bkm=PRL_BESSK0(x)
bk=PRL_BESSK1(x)

DO j=1,n-1

  bkp=bkm+j*tox*bk
  bkm=bk
  bk=bkp

ENDDO

PRL_BESSK=bk

END FUNCTION PRL_BESSK

FUNCTION PRL_BESSK0(x)
!-------------------------------------------------------------------------------
!PRL_BESSK0 gives the value of the K0 Bessel function
!
!References:
!  W.A.Houlberg, L.R.Baylor, Standardization for NTCC Module Library 1/2005
!  See Abramowitz and Stegun, Chapter 9, Dover 1964.
!
!Comments:
!  USES PRL_BESSI0
!-------------------------------------------------------------------------------

!Declaration of input variables
REAL(KIND=rspec) :: & 
  x

!Declaration of output variables
REAL (KIND=rspec) :: &
  PRL_BESSK0

!-------------------------------------------------------------------------------
!Declaration of local variables
REAL(KIND=dpspec), PARAMETER :: &
  p1=-0.57721566,   p2= 0.42278420,   p3= 0.23069756,   p4= 0.3488590e-1, &
  p5= 0.262698e-2,  p6= 0.10750e-3,   p7= 0.74e-5,                        &
  q1= 1.25331414,   q2=-0.7832358e-1, q3= 0.2189568e-1, q4=-0.1062446e-1, &
  q5= 0.587872e-2,  q6=-0.251540e-2,  q7= 0.53208e-3

REAL(KIND=dpspec) :: &
  y

!-------------------------------------------------------------------------------
!Set function value
IF(x <= 2.0) THEN

  y=x*x/4
  PRL_BESSK0=(-LOG(x/2)*PRL_BESSI0(x)) &
             +(p1+y*(p2+y*(p3+y*(p4+y*(p5+y*(p6+y*p7))))))

ELSE

  y=(2/x)
  PRL_BESSK0=(EXP(-x)/SQRT(x)) &
             *(q1+y*(q2+y*(q3+y*(q4+y*(q5+y*(q6+y*q7))))))

ENDIF

END FUNCTION PRL_BESSK0

FUNCTION PRL_BESSK1(x)
!-------------------------------------------------------------------------------
!PRL_BESSK1 gives the value of the K1 Bessel function
!
!References:
!  W.A.Houlberg, L.R.Baylor, Standardization for NTCC Module Library 1/2005
!  See Abramowitz and Stegun, Chapter 9, Dover 1964.
!
!Comments:
!  USES PRL_BESSI1
!-------------------------------------------------------------------------------

!Declaration of input variables
REAL(KIND=rspec) :: & 
  x

!Declaration of output variables
REAL (KIND=rspec) :: &
  PRL_BESSK1

!-------------------------------------------------------------------------------
!Declaration of local variables
REAL(KIND=dpspec), PARAMETER :: &
  p1= 1.0,          p2= 0.15443144,   p3=-0.67278579,   p4=-0.18156897,   &
  p5=-0.1919402e-1, p6=-0.110404e-2,  p7=-0.4686e-4,                      &
  q1= 1.25331414,   q2= 0.23498619,   q3=-0.3655620e-1, q4= 0.1504268e-1, &
  q5=-0.780353e-2,  q6= 0.325614e-2,  q7=-0.68245e-3

REAL(KIND=dpspec) :: &
  y

!-------------------------------------------------------------------------------
!Set function value
IF(x <= 2.0) THEN

  y=x*x/4
  PRL_BESSK1=(LOG(x/2)*PRL_BESSI1(x))  &
             +(1/x)*(p1+y*(p2+y*(p3+y*(p4+y*(p5+y*(p6+y*p7))))))

ELSE

  y=2/x
  PRL_BESSK1=(EXP(-x)/SQRT(x))*(q1+y*(q2+y*(q3+y*(q4+y*(q5+y*(q6+y*q7))))))

ENDIF

END FUNCTION PRL_BESSK1

END MODULE PRL_MOD
