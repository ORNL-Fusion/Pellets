module DRIFTS_MOD
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

use SPEC_KIND_MOD
IMPLICIT NONE

!>------------------------------------------------------------------
!> Private data
!>------------------------------------------------------------------

real(kind=rspec), private, save :: &
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

real(kind=rspec), private, parameter :: &
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

contains

SUBROUTINE PARKS_DRIFT(W,r_p,ne_inf,Te_inf,R,B,M0, &
                       T0,cA_inf,beta_inf,kap_c, &
                       beta_ratio, &
                       Sigma_0,c_0_bar,Psi_int,DelR)

!>------------------------------------------------------------------
!> PARKS_DRIFT calculates the change in pellet position Delta R
!> based on an ad-hoc scaling law.
!> See Parks et al. (2000)
!>
!> We use the same convention as Parks but also include the
!> corresponding PELLET parameters along with their units.
!>------------------------------------------------------------------

real(kind=rspec), intent(in) :: &
!> pellet mass in amu <br />
!> PELLET: amu_pel [-] 
  W,                            &
!> pellet radius [cm]
!> PELLET: r0 [m]
  r_p,                          &
!> background plasma density [/cm**3]
!> PELLET: d0 [/m**3]
  ne_inf,                       &
!> background plasma temperature [keV]
!> PELLET: te0 [keV]
  Te_inf,                       &
!> major radius [m]
!> PELLET: r0 [m]
  R,                            &
!> toroidal magnetic field [T]
!> PELLET: bt0 [T]
  B,                            &
!> mach number at the channel entrance [-]
!> PELLET: NA ? 
  M0,                           &
!> temperature at the channel entrance [eV]
!> PELLET: NA ?
  T0                             

!>------------------------------------------------------------------#
!> Currently these out quantities are just to compare with
!> the calculations in Parks et al. (2000)
!>------------------------------------------------------------------#

real(kind=rspec), intent(out) :: &
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

open (unit=10,file="test.txt",action="write")
write(10,*) kap_c, &
            beta_ratio,beta_0, &
            Sigma_0,Psi_int,DelR
close (10)

END SUBROUTINE PARKS_DRIFT

SUBROUTINE BAYLOR_DRIFT(B,Te0,Teped,r_pel,qa, &
                        DelR)

real(kind=rspec), intent(in) :: &
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


END SUBROUTINE

END MODULE DRIFTS_MOD