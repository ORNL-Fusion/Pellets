MODULE PRLDATA_MOD
!-------------------------------------------------------------------------------
!PRL_MOD is an F90 module of data variables and routines that support the retrieval of
! internal data from the calculation of ExB pellet
!  cloud drift in the PRL code. 
!
!References:
!
!  L.R.Baylor 9/2004
!  W.A.Houlberg, L.R.Baylor, Standardization for NTCC Module Library 1/2005
!
!Contains PUBLIC routines:
!
!  PRLDATA  - retrieves internal data from PRL code
!  PRL_PRDATA - retreives profile data and temp data from PRL code for a particular cloudlet
!-------------------------------------------------------------------------------
USE SPEC_KIND_MOD
IMPLICIT NONE


!-------------------------------------------------------------------------------
!Public data
!-------------------------------------------------------------------------------


! Proile data variables 
!
! IF pr output turned on then allocate variables and have prl populate them
!x,deno_pr(i),tempo_pr(i),preso_pr(i),vhato_pr(i), crunid


! IF temporal output turned on then allocate variables and have prl populate them
! temporal data of cloud 
!tau , psi integral,cloud Length ,  uradial


!-------------------------------------------------------------------------------
!Procedures
!-------------------------------------------------------------------------------
CONTAINS

SUBROUTINE PRLDATA()
!-------------------------------------------------------------------------------
!PRL is a 1D Lagrangian finite difference calculation of the pressure relaxation
!  model of Parks for the drift of a pellet generated cloudlet
!
! optional data to return to user depending on what he wants
!
!
!Comments:
!  Original, 2-Jul-2005, L.R. Baylor
!-------------------------------------------------------------------------------

!Declaration of input variables
INTEGER, INTENT(IN) :: &
  irunid,              & !run id number (0 if stand alone code call) [-]
  n_r                    !number of points in profile arrays [-]


REAL(KIND=rspec), INTENT(IN) :: &
  rho_r(:),            & !radial grid normalized to am (axis=0, edge = 1) [-]
  te_r(:),             & !initial electron temperature profile [eV]
  ne_r(:),             & !initial electron density profile [/cm**3]
  q_r(:)                 !safety factor profile [-]

!Declaration of output variables
REAL(KIND=rspec), INTENT(OUT) :: &
  nump                   !number of particles calculated in cloudlet 

!-------------------------------------------------------------------------------
!Declaration of local variables     

INTEGER :: &
  i,ixloc,j,jjmin,jsave,jsp,jsp2,k,na,Numshed

REAL(KIND=dpspec) :: &
  delbc,delxn,delxo,psi,taumax


!-------------------------------------------------------------------------------
!Initialization
!-------------------------------------------------------------------------------
pressure = ne_r*te_r
pq = q_r*pressure
nump = rho_r * pq * irunid * n_r
       
END SUBROUTINE PRLDATA
      

END MODULE PRLDATA_MOD
