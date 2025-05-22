MODULE DRIFTS_MOD
!-------------------------------------------------------------------------------
! DRIFTS_MOD is an F90 module of routines that calculates the change in 
! pellet penetraion depth based on three drift scaling laws.	
!
!References:
!
!  Parks 2000
!  Baylor 2007
!  Koechl
!  
!-------------------------------------------------------------------------------

USE SPEC_KIND_MOD
IMPLICIT NONE

CONTAINS

SUBROUTINE PARKS_DRIFT(r_pel)

! What needs to be included as input?
! DEFINED IN PAPER : DEFINED IN PELLET : 
! >>> DEFINITION
! r_p              : r_pel             : 
! >>> pellet radius
! t_e_inf          :                   : 
! >>> electron temperature of background plasma
! n_e_inf          :                   :
! >>> electron density of background plasma
! 

REAL(KIND=rspec), INTENT(IN) :: &
  r_pel

END SUBROUTINE PARKS_DRIFT

END MODULE DRIFTS_MOD