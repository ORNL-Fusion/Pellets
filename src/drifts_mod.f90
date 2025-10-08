MODULE DRIFTS_MOD
!!-------------------------------------------------------------------------------
!! DRIFTS_MOD is an F90 module of routines that calculates the change in 
!! pellet penetraion depth based on three drift scaling laws.	
!!
!! References:
!!
!!  Parks 2000
!!  Baylor 2007
!!  Koechl
!!  
!!-------------------------------------------------------------------------------

USE SPEC_KIND_MOD
USE LINEAR1_MOD
USE AJAX_MOD
IMPLICIT NONE
                     
CONTAINS

SUBROUTINE PARKS_DRIFT(W,r_p,ne_inf,Te_inf,R,B,M0, &
                       T0,cA_inf,beta_inf,kap_c, &
                       beta_ratio, &
                       Sigma_0,c_0_bar,Psi_int,DelR)

  !!------------------------------------------------------------------
  !! PARKS_DRIFT calculates the change in pellet position Delta R
  !! based on an ad-hoc scaling law.
  !! See Parks et al. (2000)
  !!
  !! We use the same convention as Parks but also include the
  !! corresponding PELLET parameters along with their units.
  !!------------------------------------------------------------------

REAL(KIND=RSPEC), PARAMETER :: &
!> free space permeability [kg/s**2/A**2]
  mu0 = 1.25663706e-06,                 &
!> Coloumb charge (eV --> J)
  e0 = 1.60217663e-19,                  &
!> ion mass [kg]
!> PELLET: NA
  mi = 1.6735575e-27,                   &
!> gas constant [-]
  gam = 5./3.,                          &
!> latent energy of ionization [eV/ion]
  eps_ion = 13.6,                       &
!> dissociation energy              
  eps_diss = 2.2,                       &
!> heat flux attenutation
  muE = 0.5,                            &
!> scale /cm**3 --> /m**3
  ne_scale = 1.0e6

REAL(KIND=RSPEC) :: &
!> Coulomb log for electron-neutral interactions (?)
  lnLam_en,                        &
!> Coloumb log for electron-electron interactions
  lnLam_ee,                        &
!> beta at the sonic radius
  beta_star_prime,                 &
!> temperature at the sonic radius [eV]
  T_star_prime,                    &
!> temperature at the location where M=1
  T_star,                          &
!> sound speed at the channel entrance
  cs_0,                            &
!> beta at the channel entrance [-]
  beta_0,                          &
!> background plasma density (ion) [/m**3]
!> PELLET: NA
  ni_inf,                          &
!>------------------------------------------------------------------
!> These "termX" variable names could
!> stand to be a little more descriptive
!>------------------------------------------------------------------
!> separating terms for Sigma_0 calc 
!> to keep the code a bit cleaner
  Sigma_0_term1, Sigma_0_term2,    &
  Sigma_0_term3,                   &
!> separating terms for beta_ratio calc 
!> to keep the code a bit cleaner
  beta_ratio_term1,                & 
  beta_ratio_term2,                &
!> separating terms for kap_c calc 
!> to keep the code a bit cleaner
  kap_c_term1, kap_c_term2,        &
  kap_c_term3,                     & 
  !> quant is a scalar multiple of many variables that
!> pops up a lot, so useful to have defined
!> (yes, it desperately needs a different name but I don't
!> have that computing capacity right now.)
  quant 

REAL(KIND=RSPEC), INTENT(IN) :: &
  W,                            &
    !! pellet mass in amu [-]
    !! PELLET: amu_pel [-] 
  r_p,                          &
    !! pellet radius [cm]
    !! PELLET: r0 [m]
    !! scale input by 1.0E2
  ne_inf,                       &
    !! background plasma density [/cm**3]
    !! PELLET: den0 [/m**3]
    !! scale input by 1.0E6
  Te_inf,                       &
    !! background plasma temperature [keV]
    !! PELLET: te0 [keV]
    !! no scaling
  R,                            &
    !! major radius [m]
    !! PELLET: r0 [m]
    !! no scaling
  B,                            &
    !! toroidal magnetic field [T]
    !! PELLET: bt0 [T]
    !! no scaling
  M0,                           &
    !! mach number at the channel entrance [-]
    !! PELLET: NA
  T0                             
    !! temperature at the channel entrance [eV]
    !! PELLET: NA

!>------------------------------------------------------------------#
!> Currently these out quantities are just to compare with
!> the calculations in Parks et al. (2000)
!>------------------------------------------------------------------#

REAL(KIND=RSPEC), INTENT(OUT) :: &
!> background plasma beta
  beta_inf,                      &
!> 
  c_0_bar,                       &
!> Alfven velocity of background plasma
  cA_inf,                        &
!> normalized cloud opacity (?) (dimensionless)
  Sigma_0,                       &
!> ratio of channel entrance and  background plasma betas
  beta_ratio,                    &
!> kap_c == r_perp/r_p
  kap_c,                         &
!> toroidal drive integral
  Psi_int,                       &
!> final penetration depth
  DelR    

!>------------------------------------------------------------------#
!> Calculating scalar quantities
!>------------------------------------------------------------------#

lnLam_en = log(2.0*Te_inf/7.5)
lnLam_ee = 23.5 - log(((ne_inf)**(1./2.))*((Te_inf)**(-5./6.)))
quant = ne_inf*r_p*lnLam_en
beta_inf = 4.0*mu0*(ne_inf*ne_scale)*(Te_inf*e0)/B**2.0

!>------------------------------------------------------------------
!> Calculate sonic radius quantities
!>
!> Beta and temperatures at the sonic radius r_star_prime, which is the
!> transition from sub- to super-sonic velocity, r_star_prime ~ 1.6r_p.
!>------------------------------------------------------------------

beta_star_prime = 4.3e3*((W/quant)**(1./3.))*(Te_inf**(2./3.))*beta_inf
PRINT *, "beta_star_prime: ", beta_star_prime
T_star_prime = 1.88e-9*((W/Te_inf)**(1./3.))*(quant**(2./3.))

!>------------------------------------------------------------------
!> temperature where M = 1
!>------------------------------------------------------------------
T_star = T0*(((1.0 + gam*(M0**2.0))**2.0)/(((1.0 + gam)**2.0)*(M0**2)))

!>------------------------------------------------------------------
!> background plasma Alfven velocity
!>------------------------------------------------------------------

!> assuming quasi-neutrality
ni_inf = ne_inf

cA_inf = B/((mu0*(W*mi)*(ni_inf*ne_scale)))**(1./2.)
cA_inf = cA_inf*1.0e2

!>------------------------------------------------------------------
!> Quantities at channel entrance
!>------------------------------------------------------------------

cs_0 = (2.0*gam*(T0*e0)/(W*mi))**(1./2.)
c_0_bar = ((cs_0**2.0)/gam)**(1./2.)
cs_0 = cs_0*1.0e2
c_0_bar = c_0_bar*1.0e2

kap_c_term1 = 1.54e4*(Te_inf**(1./6.))
kap_c_term2 = ((1.0 - muE)**(1./2.))*(W**(1./6.))*(quant**(1./3.))
kap_c_term3 = (4.0*gam*T_star + eps_diss + eps_ion)**(1./2.)
kap_c = (kap_c_term1/kap_c_term2)*kap_c_term3

!>------------------------------------------------------------------
!> Dimensionless quantities at channel entrance
!>------------------------------------------------------------------

!> For the life of me I cannot get Sigma_0 to match the paper value...
!> which is leading to a final DelR 50% larger than it should be.
Sigma_0_term1 = ((582.0*(R**(1./2.)))/(M0*((kap_c)**(3./2.))*(W**(1./3.))*cs_0))
Sigma_0_term2 = ((((ne_inf)**2.0)/(Te_inf*r_p))**(1./6.))
Sigma_0_term3 = (lnLam_ee/(lnLam_en**(2./3.)))
Sigma_0 = Sigma_0_term1*Sigma_0_term2*Sigma_0_term3

beta_ratio_term1 = 406.0*cs_0*(Te_inf**(5./6.))/(M0*(kap_c**2.0))
beta_ratio_term2 = (W/quant)**(2./3.)
beta_ratio = beta_ratio_term1*beta_ratio_term2

!>------------------------------------------------------------------
!> Toroidal drive integral
!>------------------------------------------------------------------

!> This expression is valid only for beta_ratio < 10
!> Can add extended definition later on

Psi_int = 0.036*(Sigma_0**1.1)*((beta_ratio - 1.0)**2.64)

!>------------------------------------------------------------------
!> Penetration depth (Eq. 28)
!>------------------------------------------------------------------

beta_0 = beta_ratio*beta_inf
DelR = 0.5*beta_0*kap_c*r_p*(cA_inf/c_0_bar)*Psi_int

OPEN (UNIT=10,FILE="parks_2000_test.txt",ACTION="write")
WRITE(10,*) kap_c, &
            beta_ratio,beta_0, &
            Sigma_0,Psi_int,DelR
CLOSE (10)

END SUBROUTINE PARKS_DRIFT

SUBROUTINE BAYLOR_DRIFT(B,Te0,Teped,r_pel,qa, &
                        DelR)

REAL(KIND=RSPEC), intent(in) :: &
!> magnetic field [T] 
  B,                            &
!> central electron temp [keV]
  Te0,                          &
!> pedestal electron temp [keV]
  Teped,                        &
!> pellet radius [m]
  r_pel,                        &
!> safety factor at the edge    
  qa

real(kind=rspec), intent(out) :: &
!> drift radius
  DelR 

DelR = B**(-0.15)*Te0**(-0.13)*Teped**(0.5)*r_pel**(0.76)*qa**(-0.15)

! open (unit=10,file="baylor_2007_test.txt",action="write")
! write(10,*) DelR
! close (10)

END SUBROUTINE

SUBROUTINE HPI2_DRIFT(v_p,r_p,ne0,Te0,alpha,Lambda,a0,r0, &
                      B0,kappa,Del_drift)

real(kind=rspec), parameter :: &
!> constants (see Table 4 in Koechl for Delta 1)  
  C1  = 0.116,    &
  C2  = 0.120,    &
  C3  = 0.368,    &
  C4  = 0.041,    &
  C5  = 0.015,    &
  C6  = 1.665,    &
  C7  = 0.439,    &
  C8  = 0.217,    &
  C9  = -0.038,   &
  C10 = 0.493,    &
  C11 = 0.193,    &
  C12 = -0.346,   &
  C13 = -0.204

real(kind=rspec), intent(in) :: &
  v_p,      &
!> pellet velocity [m/s]
!> input velocity from PELLET in [m/s]
!> no scaling needed
  ne0,      &
! !> axial electron density [1E19 1/m**3]
! !> input ne0 from PELLET in [1/m**3]
! !> so scale input by 1.0E-19.
  r_p,      &      
! !> pellet radius [mm]
! !> input pellet radius from PELLET in [m]
! !> so scale input by 1.0E3.
  Te0,      &
!> axial electron temperature [keV]
!> input Te0 from PELLET [keV]
!> no scaling needed
  alpha,    &
!> pellet injection angle w.r.t the horizontal
!> outward direction in range [-pi, pi]
!> input from PELLET NA.
  Lambda,   &
!> "impact parameter of the pellet trajectory" [-]
!> input from PELLET NA.
  a0,       &
!> minor radius [m]
!> input a0 from PELLET in [m]
!> no scaling needed
  r0,       &
!> major radius [m]
!> input r0 from PELLET in [m]
!> no scaling needed
  B0,       &
!> toroidal field strength [T]
!> input B0 from PELLET in [T]
!> no scaling needed
  kappa
!> plasma elongation close to the separatrix [-]

! real(kind=rspec), intent(in) :: & 
!   ne0,      &
! !> axial electron density [1E19 1/m**3]
! !> input ne0 from PELLET in [1/m**3]
! !> so scale input by 1.0E-19.
!   r_p      
! !> pellet radius [mm]
! !> input pellet radius from PELLET in [m]
! !> so scale input by 1.0E3.

real(kind=rspec), intent(out) :: &
  Del_drift

Del_drift = C1*((v_p/100)**C2)*(r_p**C3)*(ne0**C4) &
            *(Te0**C5)*((ABS(ABS(alpha) - C6) + C8)**C7) &
            *((1.0 - Lambda)**C9)*(a0**C10)*(R0**C11)*(B0**C12) &
            *(kappa**C13)

! open (unit=50,file="HPI2_2012_test.txt",action="write")
! write(50,*) "v_p, r_p, te0, ne0, alpha, lambda, kappa, a0, r0, B0\n"
! write(50,*) v_p, r_p, Te0, ne0, alpha, lambda, kappa, a0, r0, B0
! write(50,*) "Del_drift\n"
! write(50,*) Del_drift
! close (50)

END SUBROUTINE

SUBROUTINE APPLY_DRIFTS(ncplas,n,rho_r,rcyl_p,rflx_p,pden_r,drift,pden_r_shifted,message,iflag)

INTEGER, INTENT(IN) :: &
  ncplas,              &
  !> number of points in plasma
  n
  !> number of points (plasma+sol)

REAL(KIND=RSPEC), INTENT(IN) :: &
  rho_r(:),                     & 
  !> normalised rho array
  rcyl_p(:,:),                  & 
  !> cylindrical coordinates of pellet trajectory
  rflx_p(:,:),                  &
  !> flux coordinates of pellet trajectory
  pden_r(:),                    &
  !> density change due to injected pellet
  drift         
  !> total gradB drift in R

!!Declaration of output variables
CHARACTER(len=*), INTENT(OUT) :: &
  message                !!warning or error message [character]

INTEGER, INTENT(OUT) :: &
  iflag                  !!error and warning flag [-]
                         !!=-1 warning
                         !!=0 no warnings or errors
                         !!=1 error

REAL(KIND=RSPEC), INTENT(OUT) :: &
  pden_r_shifted(:)
  !> density change shifted due to gradB drift

REAL(KIND=RSPEC) :: &
  new_R                          
  !> shfted R value

INTEGER :: &
  ii, jj, kk
  !> counters

REAL(KIND=RSPEC), ALLOCATABLE :: &
  interp_cyl(:,:),               &
  !> cylindrical coords at rho of deposition depth
  new_cyl(:),                    &
  !> cylindrical coordinates of shifted deposition
  new_flx(:),                    &
  !> flux coordinates of shifted deposition
  rho_r_map(:)
  !> rho mapped to new deposition array

ALLOCATE(interp_cyl(3,n), &
          new_cyl(3),    &
          new_flx(3),    &
          rho_r_map(n))
          !> Allocate array sizes

interp_cyl(:,:)=0
new_cyl(:)=0
new_flx(:)=0
rho_r_map(:)=0
pden_r_shifted(:)=0
!> Initialise arrays

DO ii=1,n
  IF (pden_r(ii).NE.0) EXIT
ENDDO
!> Find the first value of dne (density change) =/= 0...
!> This will be the deposition depth.

PRINT *, "RHO_R(ii): ", rho_r(ii)

CALL LINEAR1_INTERP(n,rflx_p(1,n:1:-1),rcyl_p(1,n:1:-1),n,rho_r(n:1:-1),interp_cyl(1,n:1:-1),iflag,message)
!> Interpolate to find the R coord of rho

PRINT *, "R(ii): ", interp_cyl(1,ii)

new_R = interp_cyl(1,ii) + drift
!> Add the drift distance to R

PRINT *, "Rnew: ", new_R

CALL LINEAR1_INTERP(n,rcyl_p(1,:),rcyl_p(3,:),n,interp_cyl(1,n:1:-1),interp_cyl(3,n:1:-1),iflag,message)
!> Interpolate to find the Z coord for the (R,Z) pair

new_cyl(1) = new_R
new_cyl(3) = interp_cyl(3,ii)
!> Assign new cylindrical coords using shifted R but SAME Z

CALL AJAX_CYL2FLX(new_cyl,new_flx,iflag,message)
!> Calculate the rho coord for the shifted R location

PRINT *, "CYL_NEW: ", new_cyl
PRINT *, "FLX_NEW: ", new_flx

DO jj=1,n
  IF (rho_r(jj)>new_flx(1)) EXIT
ENDDO
!> Find closest rho in original array that matches the new shifted rho
!> Could interp again to get exact...

rho_r_map(ncplas:n) = rho_r(ncplas:n)
!> Mapped SOL rho matches original rho array

DO kk=jj,ncplas
  rho_r_map(kk) = rho_r(ncplas) + ((rho_r(kk) - rho_r(ncplas))/(rho_r(jj) - rho_r(ncplas)))*(rho_r(ii) - rho_r(ncplas))
ENDDO
!> Mapping the shifted rho array

CALL LINEAR1_INTERP(ncplas,rho_r,pden_r,ncplas-jj+1,rho_r_map(jj:ncplas),pden_r_shifted(jj:ncplas),iflag,message)
!> Interpolating dne onto the mapped rho array

pden_r_shifted(:) = (rho_r(ii) - rho_r(ncplas)) / (rho_r(jj) - rho_r(ncplas))*pden_r_shifted(:)
!> Normalised shifted dne

PRINT *, "Old density change: ", pden_r
PRINT *, "New density change: ", pden_r_shifted

END SUBROUTINE

END MODULE DRIFTS_MOD