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
IMPLICIT NONE

!>------------------------------------------------------------------
!> Private data
!>------------------------------------------------------------------

REAL(KIND=rspec), PRIVATE, SAVE :: &
!> Coulomb log for electron-neutral interactions (?)
  lnLam_en,                        &
!> beta at the sonic radius
  beta_star_prime,                 &
!> temperature at the sonic radius [eV]
  T_star_prime,                    &
!> normalized cloud opacity (?) (dimensionless)
  Sigma_0,                         &
!> ratio of channel entrance and  background plasma betas
  beta_ratio,                      &
!> kap_c == r_perp/r_p
  kap_c

CONTAINS

SUBROUTINE PARKS_DRIFT(W,r_p,ne_inf,te_inf,R,B,M0, &
                       T0)

!>------------------------------------------------------------------
!> PARKS_DRIFT calculates the change in pellet position Delta R
!> based on an ad-hoc scaling law.
!> See Parks et al. (2000)
!>
!> We use the same convention as Parks but also include the
!> corresponding PELLET parameter along with their units.
!>------------------------------------------------------------------

REAL(KIND=rspec), INTENT(IN) :: &
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
  te_inf,                       &
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
!> Defining some useful quantities that will come up a lot
!>------------------------------------------------------------------#

!! Coulomb log for electron-neutral interactions (?)
lnLam_en = 2.0*T0/7.5


!> Beta and temperatures at the sonic radius r_star_prime, which is the
!> transition from sub- to super-sonic velocity, r_star_prime ~ 1.6r_p.
  
!beta_star_prime = 

END SUBROUTINE PARKS_DRIFT

END MODULE DRIFTS_MOD