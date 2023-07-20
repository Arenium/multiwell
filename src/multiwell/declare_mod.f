!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
!                                                                      !   
!  MultiWell: a code for master equation simulations.                  !   
!  Copyright (C) 2021 John R. Barker                            !   
!                                                                      !   
!  John R. Barker                                                      !   
!  jrbarker@umich.edu                                                  !   
!  University of Michigan                                              !   
!  Ann Arbor, MI 48109-2143                                            !   
!  (734) 763 6239                                                      !   
!                                                                      !   
!  This program is free software; you can redistribute it and/or       !   
!  modify it under the terms of the GNU General Public License         !   
!  (version 2) as published by the Free Software Foundation.           !   
!                                                                      !   
!  This program is distributed in the hope that it will be useful,     !   
!  but WITHOUT ANY WARRANTY; without even the implied warranty of      !   
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the        !   
!  GNU General Public License for more details.                        !   
!                                                                      !   
!  See the 'ReadMe' file for a copy of the GNU General Public License, !   
!  or contact:                                                         !   
!                                                                      !   
!  Free Software Foundation, Inc.                                      !   
!  59 Temple Place - Suite 330                                         !   
!  Boston, MA 02111-1307, USA.                                         !   
!                                                                      !   
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
 
      MODULE declare_mod
!
!     For declaring and sharing data

      IMPLICIT NONE
      SAVE
!
!  *** SETTINGS CAN BE CHANGED BY CHANGING PARAMETER DECLARATIONS ***
!
      INTEGER, PARAMETER :: Imax = 100000       ! upper limit size of double array
      INTEGER, PARAMETER :: MaxWell = 200       !  maximum number of wells
      INTEGER, PARAMETER :: MaxProd = 200       !  maximum of product sets
      INTEGER, PARAMETER :: MaxChan = 200       !  maximum number of reaction channels from any given well
      INTEGER, PARAMETER :: MaxSup = 10         ! maximum allowed supplementary rxns for each Well (index Mol)
      REAL(8), PARAMETER :: tminlogt = 1.d-16   ! minimum time point for logt output 
      REAL(8), PARAMETER :: tunthresh = 1.0e-12 ! minimum tunneling probability; if the tunneling probability is less
                                                ! than tunthresh, calculation of the tunneling corrected sum of states
                                                ! for the transition state is stopped

      INTEGER, PARAMETER :: Ntime = 101         ! number of time-steps for reporting results in standard output
      INTEGER, PARAMETER :: Ndist = 1001        ! max number of energy bins for reporting vibrational distributions
      INTEGER, PARAMETER :: Mtime = 101         ! number of time-steps for reporting vibrational distributions

      REAL(8), PARAMETER :: Qthresh = 1.0e-12   ! when the relative contribution to the partition function
                                                ! becomes less than Qthresh we stop evaluating the partition function

      REAL(8), PARAMETER :: FRACT = 0.5d+00     ! Fraction of step size used in COLSTEP and COLNORM numerical integrations
      REAL(8), PARAMETER :: ERRo = 1.d-05       ! for convergence in COLSTEP and COLNORM numerical integration

!  ---------------------------------------------------------------------------------------------------------------------
!     ALLOCATABLE ARRAYS 

      INTEGER :: istat
      
      REAL(8), ALLOCATABLE, DIMENSION(:,:)           :: Dens     ! density of states for molecule Mol (MaxWell,Imax)
      REAL(8), ALLOCATABLE, DIMENSION(:,:)           :: COLLUP   ! probability of an up-step (MaxWell,Imax)
      REAL(8), ALLOCATABLE, DIMENSION(:,:)           :: CNORM    ! normalization factor for energy transfer collision model (MaxWell,Imax)
      REAL(8), ALLOCATABLE, DIMENSION(:,:)           :: DensTS   ! density of states for a specific transition state; the first part of the
      REAL(8), ALLOCATABLE, DIMENSION(:,:,:)         :: Rate     ! k(E) for a given molecule and channel (no tunneling) (MaxWell,MaxChan,Imax)
      REAL(8), ALLOCATABLE, DIMENSION(:,:,:)         :: TSsum    ! sum of states for a given molecule and channel (no tunneling) (MaxWell,MaxChan,Imax)
      REAL(8), ALLOCATABLE, DIMENSION(:,:)           :: TSdensBK ! density of TS states for one channel, including tunneling (MaxWell,2*Imax)
      REAL(8), ALLOCATABLE, DIMENSION(:,:,:)         :: TRate    ! k(E) with tunneling correction (MaxWell,MaxChan,Imax*2)
      REAL(8), ALLOCATABLE, DIMENSION(:,:,:)         :: KRate    ! k(E) bookkeeping array containing Rate and TRate; used for output purposes only (MaxWell,MaxChan,Imax*2)
      REAL(8), ALLOCATABLE,DIMENSION( : )            :: Pstart   ! Normalized cumulative distribution of starting energies
      INTEGER, ALLOCATABLE,DIMENSION( : )            :: Nstart   ! actual cumulative probability distribution selected by Monte Carlo
      INTEGER                                        :: Jsize    ! size of the cumulative initial energy distribution
      INTEGER, ALLOCATABLE, DIMENSION(:,:)           :: WARN1    ! Warning that k(E) from reactant an empty energy grain in product well
      INTEGER, ALLOCATABLE, DIMENSION(:)             :: IMol     ! index number of any given molecule
      REAL(8), ALLOCATABLE, DIMENSION(:)             :: HMol     ! enthalpy (actually internal energy) of a molecule or product set
      INTEGER, ALLOCATABLE, DIMENSION(:)             :: Molsym   ! symmetry number of a given molecule
      INTEGER, ALLOCATABLE, DIMENSION(:)             :: Molopt   ! number of optical isomers of a given molecule
      REAL(8), ALLOCATABLE, DIMENSION(:)             :: MolMom   ! rotation constant for 2-D external rotation of a given molceule
      REAL(8), ALLOCATABLE, DIMENSION(:)             :: Molele   ! electronic partition function of a given molecule
      REAL(8), ALLOCATABLE, DIMENSION(:)             :: Viblo    ! lowest vibrational frequency in Mol
      REAL(8), ALLOCATABLE, DIMENSION(:)             :: Sig      ! Lennard-Jones sigma for Mol
      REAL(8), ALLOCATABLE, DIMENSION(:)             :: Eps      ! Lennard-Jones epsilon for Mol
      REAL(8), ALLOCATABLE, DIMENSION(:)             :: Sigma    ! Lennard-Jones, combined
      REAL(8), ALLOCATABLE, DIMENSION(:)             :: Epsil    ! Lennard-Jones, combined
      REAL(8), ALLOCATABLE, DIMENSION(:)             :: Mass     ! reduced mass
      REAL(8), ALLOCATABLE, DIMENSION(:)             :: klj      ! Lennard-Jones collision rate constant
      REAL(8), ALLOCATABLE, DIMENSION(:)             :: kqm      ! quantum total bimolecular collision frequency rate constant
      REAL(8), ALLOCATABLE, DIMENSION(:)             :: CollK    ! bimolecular collision frequency rate constant
      INTEGER, ALLOCATABLE, DIMENSION(:)             :: Nchan    ! number of rxn channels from a given molecule
      INTEGER, ALLOCATABLE, DIMENSION(:,:)           :: Jrev     ! flag for reversible reaction
      INTEGER, ALLOCATABLE, DIMENSION(:,:)           :: Jto      ! index number of product molecule
      INTEGER, ALLOCATABLE, DIMENSION(:,:)           :: TSsym    ! symmetry number of a transiton state
      REAL(8), ALLOCATABLE, DIMENSION(:,:)           :: TSele    ! electronic partition function of a transition state
      INTEGER, ALLOCATABLE, DIMENSION(:,:)           :: TSopt    ! number of optical isomers of a transition state
      INTEGER, ALLOCATABLE, DIMENSION(:,:)           :: Jread    ! flag for how k(E)s are determined: ILT, SUM, CRP, or from an external file
      REAL(8), ALLOCATABLE, DIMENSION(:,:)           :: Path     ! rxn path degeneracy
      REAL(8), ALLOCATABLE, DIMENSION(:,:)           :: Afac     ! A-factor used with ILT method (Afac<0 for Reactant conc (molecule/cc) for pseudo-1st-order bimol rxn 
      REAL(8), ALLOCATABLE, DIMENSION(:,:)           :: Eo       ! critical energy for a reaction channel
      REAL(8), ALLOCATABLE, DIMENSION(:,:)           :: Eor      ! critical energy for a reaction channel, corrected for rotational energy
      REAL(8), ALLOCATABLE, DIMENSION(:,:)           :: Hts      ! enthalpy (internal energy) of a transition state
      INTEGER, ALLOCATABLE, DIMENSION(:,:)           :: NCENT    ! flag for applying centrifugal correction
      REAL(8), ALLOCATABLE, DIMENSION(:,:)           :: TSmom    ! moment of inertia for 2-D external rotation of a transition state
      REAL(8), ALLOCATABLE, DIMENSION(:,:)           :: Press    ! for microcanonical pseudo-first order rxn
      REAL(8), ALLOCATABLE, DIMENSION(:)             :: Var
      REAL(8), ALLOCATABLE, DIMENSION(:)             :: FracMol
      REAL(8), ALLOCATABLE, DIMENSION(:)             :: AveE
      INTEGER, ALLOCATABLE, DIMENSION(:)             :: Nreaindex
      INTEGER, ALLOCATABLE, DIMENSION(:)             :: Nproindex
      INTEGER, ALLOCATABLE, DIMENSION(:,:)           :: iivr     ! flag to calculate IVR (iivr=1), or not (iivr=0)
      REAL(8), ALLOCATABLE, DIMENSION(:,:)           :: vivr     ! characteristic reaction frequency (as in RRK theory) (cm-1)
      REAL(8), ALLOCATABLE, DIMENSION(:,:)           :: vave     ! average molecular frequency (cm-1)
      REAL(8), ALLOCATABLE, DIMENSION(:,:)           :: pivr     ! bimolecular collision rate constant for collision-induced IVR
      REAL(8), ALLOCATABLE, DIMENSION(:,:)           :: tivr     ! IVR threshold energy (cm-1), measured from the  reaction threshold (Eor)
      REAL(8), ALLOCATABLE, DIMENSION(:,:,:)         :: civr     ! coefficients for 2nd order fit of 
      INTEGER, ALLOCATABLE, DIMENSION(:)             :: nsmax    ! number of supplementary rxns for a Well (index Mol) 
      INTEGER, ALLOCATABLE, DIMENSION(:,:)           :: nsup     ! index number for supplentary rxn for a Well (index Mol)
      INTEGER, ALLOCATABLE, DIMENSION(:,:)           :: sto      ! supplentary rxn product index
      INTEGER, ALLOCATABLE, DIMENSION(:,:)           :: norder   ! supplentary rxn reaction order
      REAL(8), ALLOCATABLE, DIMENSION(:,:)           :: Afc      ! supplentary rxn Arrhenius A-factor (s-1 or cm3 s-1)
      REAL(8), ALLOCATABLE, DIMENSION(:,:)           :: Bsr      ! supplentary rxn E(activation)/R
      REAL(8), ALLOCATABLE, DIMENSION(:,:)           :: dde      ! Quench energy (energy lost to environment)
      REAL(8), ALLOCATABLE, DIMENSION(:,:)           :: ksup     ! Supplemntary rate constant = Afc*exp(-Bsr/T)
      REAL(8), ALLOCATABLE, DIMENSION(:)             :: CNORMref ! reference normalization for collision frequency at rxn Eor ----> NEW UP-STEP COLLISION METHOD
      REAL(8), ALLOCATABLE, DIMENSION(:)             :: Eref     ! reference energy for CNORMref                              ----> NEW UP-STEP COLLISION METHOD
      INTEGER, ALLOCATABLE, DIMENSION(:)             :: iset
      INTEGER, ALLOCATABLE, DIMENSION(:)             :: iset2
      REAL(8), ALLOCATABLE, DIMENSION(:,:)           :: V1       ! energy barrier in Eckart?
      REAL(8), ALLOCATABLE, DIMENSION(:,:)           :: vimag    ! imaginary frequency
      INTEGER, ALLOCATABLE, DIMENSION(:,:)           :: itun     ! flag designating tunneling (=1), or not (0)
      INTEGER, ALLOCATABLE, DIMENSION(:)             :: ITYPE    ! flag designating TYPE of collisional energy transfer model
      REAL(8), ALLOCATABLE, DIMENSION(:,:)           :: DC       ! coefficients defining a energy transfer collision model
      REAL(8), ALLOCATABLE, DIMENSION(:)             :: Rab      ! Rab(i) are the reaction rates and collision rate calculated in STEPPER and used in BOOK1
      REAL(8), ALLOCATABLE, DIMENSION(:,:,:)         :: Flux
      CHARACTER(len=2),  ALLOCATABLE, DIMENSION(:)   :: LJQM     ! dummy label for type of collision rate constant
      CHARACTER(len=100),ALLOCATABLE, DIMENSION(:,:) :: LI       ! comment lines associated with each energy transfer type
      CHARACTER(len=10), ALLOCATABLE, DIMENSION(:)   :: MolName  ! names of Wells and product-sets
      CHARACTER(len=10), ALLOCATABLE, DIMENSION(:,:) :: TSname   ! name of a transition state
      CHARACTER(len=10), ALLOCATABLE, DIMENSION(:,:) :: SupName  ! name of supplementary reaction

c      REAL(8), DIMENSION(MaxWell,Imax)          :: Dens     ! density of states for molecule Mol (MaxWell,Imax)
c      REAL(8), DIMENSION(MaxWell,Imax)          :: COLLUP   ! probability of an up-step (MaxWell,Imax)
c      REAL(8), DIMENSION(MaxWell,Imax)          :: CNORM    ! normalization factor for energy transfer collision model (MaxWell,Imax)
c      REAL(8), DIMENSION(MaxChan*MaxWell,Imax)  :: DensTS   ! density of states for a TS; the first part of the
c      REAL(8), DIMENSION(MaxWell,MaxChan,Imax)  :: Rate     ! k(E) for a given molecule and channel (no tunneling) (MaxWell,MaxChan,Imax)
c      REAL(8), DIMENSION(MaxWell,MaxChan,Imax)  :: TSsum    ! sum of states for a given molecule and channel (no tunneling) (MaxWell,MaxChan,Imax)
c      REAL(8), DIMENSION(MaxChan*MaxWell,2*Imax):: TSdensBK ! sum of states for a given molecule and channel, including tunneling (MaxWell,2*Imax)
c      REAL(8), DIMENSION(MaxWell,MaxChan,2*Imax):: TRate    ! k(E) with tunneling correction (MaxWell,MaxChan,Imax*2)
c      REAL(8), DIMENSION(MaxWell,MaxChan,2*Imax):: KRate    ! k(E) bookkeeping array containing Rate and TRate; used for output purposes only (MaxWell,MaxChan,Imax*2)

c      INTEGER, DIMENSION(MaxWell,MaxChan)       :: WARN1    ! Warning that k(E) from reactant an empty energy grain in product well
c      INTEGER, DIMENSION(MaxWell+MaxProd)       :: IMol     ! index number of any given molecule
c      REAL(8), DIMENSION(MaxWell+MaxProd)       :: HMol     ! enthalpy (actually internal energy) of a molecule or product set
c      INTEGER, DIMENSION(MaxWell)               :: Molsym   ! symmetry number of a given molecule
c      INTEGER, DIMENSION(MaxWell)               :: Molopt   ! number of optical isomers of a given molecule
c      REAL(8), DIMENSION(MaxWell)               :: MolMom   ! rotation constant for 2-D external rotation of a given molceule
c      REAL(8), DIMENSION(MaxWell)               :: Molele   ! electronic partition function of a given molecule
c      REAL(8), DIMENSION(MaxWell)               :: Viblo    ! lowest vibrational frequency in Mol
c      REAL(8), DIMENSION(MaxWell)               :: Sig      ! Lennard-Jones sigma for Mol
c      REAL(8), DIMENSION(MaxWell)               :: Eps      ! Lennard-Jones epsilon for Mol
c      REAL(8), DIMENSION(MaxWell)               :: Sigma    ! Lennard-Jones, combined
c      REAL(8), DIMENSION(MaxWell)               :: Epsil    ! Lennard-Jones, combined
c      REAL(8), DIMENSION(MaxWell)               :: Mass     ! reduced mass
c      REAL(8), DIMENSION(MaxWell)               :: klj      ! Lennard-Jones collision rate constant
c      REAL(8), DIMENSION(MaxWell)               :: kqm      ! quantum total bimolecular collision frequency rate constant
c      REAL(8), DIMENSION(MaxWell)               :: CollK    ! bimolecular collision frequency rate constant
c      CHARACTER(len=2), DIMENSION(MaxWell)      :: LJQM     ! dummy label for type of collision rate constant
c      INTEGER, DIMENSION(MaxWell)               :: Nchan    ! number of rxn channels from a given molecule
c      INTEGER, DIMENSION(MaxWell,MaxChan)       :: Jrev     ! flag for reversible reaction
c      INTEGER, DIMENSION(MaxWell,MaxChan)       :: Jto      ! index number of product molecule
c      INTEGER, DIMENSION(MaxWell,MaxChan)       :: TSsym    ! symmetry number of a transiton state
c      REAL(8), DIMENSION(MaxWell,MaxChan)       :: TSele    ! electronic partition function of a transition state
c      INTEGER, DIMENSION(MaxWell,MaxChan)       :: TSopt    ! number of optical isomers of a transition state
c      INTEGER, DIMENSION(MaxWell,MaxChan)       :: Jread    ! flag for how k(E)s are determined: ILT, SUM, CRP, or from an external file
c      REAL(8), DIMENSION(MaxWell,MaxChan)       :: Path     ! rxn path degeneracy
c      REAL(8), DIMENSION(MaxWell,MaxChan)       :: Afac     ! A-factor used with ILT method (Afac<0 for Reactant conc (molecule/cc) for pseudo-1st-order bimol rxn 
c      REAL(8), DIMENSION(MaxWell,MaxChan)       :: Eo       ! critical energy for a reaction channel
c      REAL(8), DIMENSION(MaxWell,MaxChan)       :: Eor      ! critical energy for a reaction channel, corrected for rotational energy
c      REAL(8), DIMENSION(MaxWell,MaxChan)       :: Hts      ! enthalpy (internal energy) of a transition state
c      INTEGER, DIMENSION(MaxWell,MaxChan)       :: NCENT    ! flag for applying centrifugal correction
c      REAL(8), DIMENSION(MaxWell,MaxChan)       :: TSmom    ! moment of inertia for 2-D external rotation of a transition state
c      REAL(8), DIMENSION(MaxWell,MaxChan)       :: Press    ! for microcanonical pseudo-first order rxn
c      REAL(8), DIMENSION( MaxWell + MaxProd )   :: Var
c      REAL(8), DIMENSION( MaxWell + MaxProd )   :: FracMol
c      REAL(8), DIMENSION( MaxWell + MaxProd )   :: AveE
c      INTEGER, DIMENSION( MaxWell*MaxChan )     :: Nreaindex
c      INTEGER, DIMENSION( MaxWell*MaxChan )     :: Nproindex
c      INTEGER, DIMENSION( MaxWell , MaxChan )   :: iivr     ! flag to calculate IVR (iivr=1), or not (iivr=0)
c      REAL(8), DIMENSION( MaxWell , MaxChan )   :: vivr     ! characteristic reaction frequency (as in RRK theory) (cm-1)
c      REAL(8), DIMENSION( MaxWell , MaxChan )   :: vave     ! average molecular frequency (cm-1)
c      REAL(8), DIMENSION( MaxWell , MaxChan )   :: pivr     ! bimolecular collision rate constant for collision-induced IVR
c      REAL(8), DIMENSION( MaxWell , MaxChan )   :: tivr     ! IVR threshold energy (cm-1), measured from the  reaction threshold (Eor)
c      REAL(8), DIMENSION( MaxWell , MaxChan, 3 ):: civr     ! coefficients for 2nd order fit of 
c      INTEGER, DIMENSION(MaxWell)               :: nsmax    ! number of supplementary rxns for a Well (index Mol) 
c      INTEGER, DIMENSION(MaxWell,MaxSup)        :: nsup     ! index number for supplentary rxn for a Well (index Mol)
c      INTEGER, DIMENSION(MaxWell,MaxSup)        :: sto      ! supplentary rxn product index
c      INTEGER, DIMENSION(MaxWell,MaxSup)        :: norder   ! supplentary rxn reaction order
c      REAL(8), DIMENSION(MaxWell,MaxSup)        :: Afc      ! supplentary rxn Arrhenius A-factor (s-1 or cm3 s-1)
c      REAL(8), DIMENSION(MaxWell,MaxSup)        :: Bsr      ! supplentary rxn E(activation)/R
c      REAL(8), DIMENSION(MaxWell,MaxSup)        :: dde      ! Quench energy (energy lost to environment)
c      REAL(8), DIMENSION(MaxWell,MaxSup)        :: ksup     ! Supplemntary rate constant = Afc*exp(-Bsr/T)
c      REAL(8), DIMENSION(MaxWell)               :: CNORMref ! reference normalization for collision frequency at rxn Eor ----> NEW UP-STEP COLLISION METHOD
c      REAL(8), DIMENSION(MaxWell)               :: Eref     ! reference energy for CNORMref                              ----> NEW UP-STEP COLLISION METHOD
c      INTEGER, DIMENSION(MaxWell)               :: iset, iset2
c      REAL(8), DIMENSION(MaxWell,MaxChan)       :: V1       ! energy barrier in Eckart?
c      REAL(8), DIMENSION(MaxWell,MaxChan)       :: vimag    ! imaginary frequency
c      INTEGER, DIMENSION(MaxWell,MaxChan)       :: itun     ! flag designating tunneling (=1), or not (0)
c      INTEGER, DIMENSION(MaxWell)               :: ITYPE    ! flag designating TYPE of collisional energy transfer model
c      REAL(8), DIMENSION(MaxWell,8)             :: DC       ! coefficients defining a energy transfer collision model
c      CHARACTER(len=100), DIMENSION(MaxWell,5)  :: LI       ! comment lines associated with each energy transfer type
c      REAL(8), DIMENSION(MaxChan+1)             :: Rab      ! Rab(i) are the reaction rates and collision rate calculated in STEPPER and used in BOOK1
c      REAL(8), DIMENSION( MaxWell , MaxChan , Ntime)  :: Flux
c      CHARACTER(len=10), DIMENSION(MaxWell+MaxProd)   :: MolName  ! names of Wells and product-sets
c      CHARACTER(len=10), DIMENSION(MaxWell,MaxChan)   :: TSname   ! name of a transition state
c      CHARACTER(len=10), DIMENSION(MaxWell,MaxSup)    :: SupName  ! name of supplementary reaction


!  ---------------------------------------------------------------------------------------------------------------------
!
!      NIST Physical Constants (1998, 2000)

      REAL(8), PARAMETER :: clight = 2.99792458d+10       ! cm/s
      REAL(8), PARAMETER :: hor = 1.4387752d+00           ! h / R (cm-1)
      REAL(8), PARAMETER :: jtocm = 83.59347d+00          ! kJ/mol to cm-1
      REAL(8), PARAMETER :: caltoj = 4.184d+00            ! kcal/mol to kJ/mol  (thermochemical calorie)
      REAL(8), PARAMETER :: caltocm = 349.7551d+00        ! caltocm = kcal/mol to cm-1    (thermochemical calorie)
      REAL(8), PARAMETER :: plank = 6.62606896d-34        ! h; Planck's constant
      REAL(8), PARAMETER :: avogadro = 6.02214179d+23     ! Avogadro's number

      REAL(8), PARAMETER :: dpi = 3.1415926535897932384626433  ! pi

!
!      Input/output unit numbers, etc.

      INTEGER, PARAMETER :: KIN   = 2     ! Data file
      INTEGER, PARAMETER :: KSTD  = 6     ! standard output (screen)
      INTEGER, PARAMETER :: KSMM  = 7     ! Short summary of output file 
      INTEGER, PARAMETER :: KARY  = 8     ! Summary of array variables
      INTEGER, PARAMETER :: KOUT  = 9     ! General time-dependent output
      INTEGER, PARAMETER :: KRAT  = 10    ! average rate 'constants'
      INTEGER, PARAMETER :: KFLX  = 11    ! relative reactive flux via a given path
      INTEGER, PARAMETER :: KDIS  = 12    ! vib distribution

      CHARACTER(len=46)             :: dumcut
      CHARACTER(len=46), PARAMETER :: 
     &          cut = '**************INPUT DATA SUMMARY**************'

!
!     Flags

      INTEGER :: ntunflag                 ! flag for tunneling
      INTEGER :: nivrflag                 ! flag for slow IVR
      INTEGER :: Pressflag                ! flag for reactant [B]  (for microcanonical pseudo-first order rxn)
!
!      Energy grains and maximum energies for double arrays (cm-1)

      REAL(8)     :: Egrain1
      REAL(8)     :: Egrain2
      REAL(8)     :: Emax1
      REAL(8)     :: Emax2
      INTEGER     :: imax1
!
!     Wells and products

      INTEGER                             :: Isize     ! size of double array in a given calculation
      INTEGER                             :: NWells      ! number of wells in a given calculation
      INTEGER                             :: NProds      ! number of product sets in a given calculation
!
!     Collision rate constants

      REAL(8)                                   :: SigM     ! LJ sigma of collider gas
      REAL(8)                                   :: EpsM     ! LJ epsilon of collider gas
      REAL(8)                                   :: AmuM     ! Molecular weight of collider gas
      REAL(8)                                   :: Amu            ! Molecular weight of Mol
      CHARACTER (len=2), PARAMETER              :: LJ = 'LJ'      ! label for Lennard-Jones collison rate constant
      CHARACTER (len=2), PARAMETER              :: QM = 'QM'      ! label for total quantum mechanical rate constant
      
!     Reactions and Transition States
!
      INTEGER                                   :: Nreac    ! total number of forward and reverse reactions
      INTEGER                                   :: largest  ! largest value of Nchan among all of the Wells
      INTEGER                             :: nsux     ! total number of supplementary reactions
!
! in MAIN, STEPPER and BOOK1 subroutines

      CHARACTER(len=4)  :: Tspec          ! simulated timescale: 'time', 'coll', or 'logt'
!
!     Other variables


      REAL(8) :: Statest      ! test to determine if a state is present in an energy grain
      REAL(8) :: CE
      INTEGER :: IDUM   ! random number seed

      CHARACTER(len=5)              :: ETKEY    ! dummy variable for key-word
      CHARACTER(len=5), PARAMETER         :: oldet = 'OLDET'      ! traditional treatment of energy transfer at low E
      CHARACTER(len=5), PARAMETER         :: newet = 'NEWET'      ! Barker's 'New Approach" to energy transfer [J.R. Barker, Int. J. Chem. Kinetics, 41, 748-763 (2009)]
      CHARACTER(len=5), PARAMETER         :: xxxet = 'XXXET'      ! if not specified by user, then use Barker's 'New Approach"


      REAL(8) :: Hstart       ! Energy step size for initial energy distribution Pstart
      REAL(8) :: Edels


!  ---------------------------------------------------------------------
!  ---------------------------------------------------------------------


      END MODULE declare_mod
