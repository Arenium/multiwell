!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
!                                                                      !   
!  MultiWell: a code for master equation simulations.                  !   
!  Copyright (C) 2001 - 2021 John R. Barker                            !   
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

      MODULE bookstep_mod
      
!     The following procedures are contained in this module:
!
!     SUBROUTINE STEPPER
!     SUBROUTINE BOOK1
!     SUBROUTINE INITIAL
!     FUNCTION rivr

      USE declare_mod
      IMPLICIT NONE
      SAVE
      
      REAL(8)               :: Tlim                  ! maximum simulated time duration
      REAL(8)               :: Tstep                 ! time step
      REAL(8)               :: DistStep              ! energy step size for binning the energy distribution

      REAL(8), ALLOCATABLE, DIMENSION(:,:)      :: EMol       ! reported average vibrational energy
      REAL(8), ALLOCATABLE, DIMENSION(:,:,:)    :: EDIST      ! vibrational energy distribution (MaxWell,Ndist,Mtime)
      REAL(8), ALLOCATABLE, DIMENSION(:,:,:)    :: kuni       ! reported average unimolecular rate constant
      INTEGER, ALLOCATABLE, DIMENSION(:,:)      :: Num        ! number of trials which access a given time bin
      INTEGER, ALLOCATABLE, DIMENSION(:,:)      :: Mum        ! number of trials which access a given time bin for EDIST
      INTEGER, ALLOCATABLE, DIMENSION(:,:)      :: Nuni       ! number of kuni samples for a given time increment

c      REAL(8), DIMENSION(MaxWell,Ntime)               :: EMol       ! reported average vibrational energy
c      REAL(8), DIMENSION(MaxWell,Ndist,Mtime)         :: EDIST      ! vibrational energy distribution (MaxWell,Ndist,Mtime)
c      REAL(8), DIMENSION(MaxWell,MaxChan,Ntime)       :: kuni       ! reported average unimolecular rate constant
c      INTEGER, DIMENSION(MaxWell+MaxProd,Ntime)       :: Num        ! number of trials which access a given time bin
c      INTEGER, DIMENSION(MaxWell,Mtime)               :: Mum        ! number of trials which access a given time bin for EDIST
c      INTEGER, DIMENSION(MaxWell,Ntime)               :: Nuni       ! number of kuni samples for a given time increment

      CONTAINS

! ---------------------------------------------------------------------------------------------------------------------
      SUBROUTINE STEPPER(Mol,E,NumbDens,IP,Delt)
 
c       Carries out Gillespie'st Stochastic Simulation Algorithm for one step

      USE declare_mod
      USE utili_mod
      USE subs_mod
      IMPLICIT NONE
      SAVE 

      INTEGER, INTENT(IN)    :: Mol               ! index of Well
      REAL(8), INTENT(INOUT) :: E                 ! energy in Well
      REAL(8), INTENT(IN)    :: numbdens          ! number density of collider gas
      INTEGER, INTENT(OUT)   :: IP                ! index of selected process
      REAL(8), INTENT(OUT)   :: Delt              ! time step (s)
      
      INTEGER :: i , IE  , IPtot
      REAL(8) :: Atot , RX , arg , X , Atemp , 
     &           Y1 , Y2 , Y3 , Aref , CC , arg2 , rsup

      IF ( ETKEY.EQ.newet ) THEN                      ! Barker's NEW energy transfer treatment <<<<<<<<<<
        CC = fxnofe(E , Mol , CNORM)
        CC = CC/CNORMref(Mol)       
      ELSEIF ( ETKEY.EQ.oldet ) THEN                  ! Traditional energy transfer treatment
        CC = 1.0d+00                                  
      ENDIF

      IPtot = Nchan(Mol) + nsmax(Mol) + 1                 ! rxn channels, supplemenary rxns, and collisions for a given Mol
      Rab(IPtot) = CollK(Mol)*numbdens*CC                 ! Inelastic collision frequency
      Atot = Rab(IPtot)                                   ! Collision frequency
      IF ( Nchan(Mol).GT.0 ) THEN                         ! If reaction channels
        DO i = 1 , Nchan(Mol)
          Rab(i) = 0.0D+00
          arg = E - Eor(Mol,i)                             ! Active energy in TS to determine k(E); centrifugal correction
          arg2 =-imax1*Egrain1                             ! minimum active energy in TS for the tunneling case
          IF( (itun(Mol,i) .eq. 1 ).and.( arg .ge. arg2) ) THEN     ! TUNNELING
            call NTERP2(arg,IE,X)
            Y1 = TRate(Mol,i,IE-1)              ! tunneling k(E)
            Y2 = TRate(Mol,i,IE)                ! tunneling k(E)
            Y3 = TRate(Mol,i,IE+1)              ! tunneling k(E)
            Rab(i) = EXP( XLINT(Y1,Y2,Y3,X) )     ! interpolated k(E)
          ELSEIF ( arg .GE. 0.0D+00 ) THEN            ! when NO TUNNELING and above rxn threshold energy
            call NTERP(arg,IE,X)
            Y1 = rivr( Mol, i, IE-1, numbdens )     ! IVR transmission coefficient
            Y2 = rivr( Mol, i, IE, numbdens )       ! IVR transmission coefficient
            Y3 = rivr( Mol, i, IE+1, numbdens )     ! IVR transmission coefficient
            Rab(i) = XLINT(Y1,Y2,Y3,X) * rxnofe( arg, Mol, i, Rate )  ! interpolated rivr * interpolated k(E)
          ENDIF
          Atot = Atot + Rab(i)
        END DO  ! Nchan
      ENDIF  ! if Nchan
      
      IF ( nsmax(Mol) .GT. 0 ) THEN                            ! If canonical supplementary reactions
         DO i = 1, nsmax(Mol)
           IF ( norder(Mol,i) .EQ. 1) THEN                     ! order of reaction
              rsup = ksup( Mol, i )                            ! first-order 
            ELSEIF ( norder(Mol,i) .EQ. 2) THEN
              rsup = numbdens*ksup( Mol, i )                     ! PSEUDO-first-order (2nd order rate constant)
           ENDIF
           Rab( Nchan(Mol) + i ) = rsup
           Atot = Atot + rsup                                      ! includes supplementary rxn channels
        END DO
      ENDIF   ! if supplementary rxns nsmax(Mol)

      IF ( Atot .GT. 1.e-15 ) THEN      ! Probably small enough for most applications.....
        RX = RAN1(IDUM)
        Aref = RX*Atot
        Atemp = 0.0D+00        
        IP = IPtot + 1                                   ! includes all processes
        DO WHILE ( (Atemp.LE.Aref) .AND. (IP .GT. 1 ) )       ! select process; IP = process number (index)
          IP = IP - 1
          Atemp = Atemp + Rab(IP)
        ENDDO         
        RX = RAN1(IDUM)                                  ! fresh random number
        IF ( RX .LE. 1.00D-10 ) RX = RAN1(IDUM)          ! if random number too small, call another one
        Delt = -LOG(RX)/Atot                             ! Determine time step
      ELSE
        IP =  IPtot                     ! i.e., select a collision
        Delt = Tlim
      ENDIF

      RETURN
      END SUBROUTINE stepper

! ---------------------------------------------------------------------------------------------------------------------
      SUBROUTINE BOOK1(Mol,E,Temp,IP,Time,Delt,Igo)
      
c       Bookkeeping routine, works in concert with Subroutine STEPPER
c
      USE declare_mod
      USE utili_mod
      USE subs_mod

      IMPLICIT NONE
      SAVE 
      INTEGER, INTENT(INOUT)        :: Mol              ! index for Well
      REAL(8), INTENT(INOUT)        :: E              ! Energy
      REAL(8), INTENT(IN)               :: Temp       ! Temperature
      INTEGER, INTENT(IN)               :: IP              ! = rxn channel number, or = (Nchan+1) for collision
      REAL(8), INTENT(INOUT)        :: Time       ! Time (s)
      REAL(8), INTENT(IN)               :: Delt       ! Time step (s)
      INTEGER, INTENT(OUT)              :: Igo        ! Process: contineu stepping (=1), stop (=0)

      INTEGER :: IPtot , IPP , ifin , i , Last , N , ifind , Dlast 
      INTEGER :: k , kstep , klast
      REAL(8) :: Step , rat
      REAL(8) :: Enew, ddel

      kstep = Ntime                                          ! number of time steps for collision sampling

      N = 1 + INT(E/DistStep)                            ! Bin energy index for binned distributions
      IF ( N .GT. Ndist ) N = Ndist
c
c       First, record the current step
c
      IF ( Time.EQ.0.0D+00 ) THEN                      ! First step of a trial
         Last = 1                                      ! index of final bin for last step
         Dlast = 1                                     ! index of final bin for vib distribution
         klast = 1
         Num(Mol,1) = Num(Mol,1) + 1                   ! Number of times 1st time bin accessed (for averages)
         Mum(Mol,1) = Mum(Mol,1) + 1                   ! Number of times 1st time bin accessed (for EDIST)
         EMol(Mol,1) = EMol(Mol,1) + E                 ! For calculating average vibrational energy
         EDIST(Mol,N,1) = EDIST(Mol,N,1) + 1.d+00      ! For calculating vibrational distribution
         DO k = 1 , Nchan(Mol) + nsmax(Mol)            ! Rab(k)= k(E) for kth rxn channel calculated by STEPPER [includes IVR effects]; includes supplementary rxns
            kuni(Mol,k,1) = Rab(k) + kuni(Mol,k,1)     ! for computing average rate constant [includes IVR effects]
         END DO  ! k
      ENDIF

      IF ( Tspec.EQ.'TIME' .OR. Tspec.EQ.'COLL') THEN
          ifin = 1 + INT(((Time+Delt)/Tlim)*(Ntime-1.0) )         ! index of  time-bin at the end of this step
        ELSEIF ( Tspec.EQ.'LOGT' ) THEN
          rat = log(( Time + Delt + tminlogt )/tminlogt ) 
          rat = rat / log( Tlim/tminlogt)
          ifin = Ntime*rat
      ENDIF
      
      IF ( Time+Delt .GE. Tlim ) THEN
         ifin = Ntime
         Igo = 0                                       ! Stop stepping
      ENDIF

      IF ( ifin.GT.Last ) THEN                         ! snapshot at end of each time interval
         DO i = Last+1 , ifin                          ! i = time bin; ifin is .LE. Ntime (see earlier lines)
            Num(Mol,i) = Num(Mol,i) + 1                ! Number of times ith time bin accessed
            EMol(Mol,i) = EMol(Mol,i) + E              ! for average vibrational energy
            DO k = 1 , Nchan(Mol)                      ! Rab(k)= k(E) for kth rxn channel calculated by STEPPER [includes IVR effects]
               kuni(Mol,k,i) = Rab(k) + kuni(Mol,k,i)  ! for computing average rate constant [includes IVR effects]
            END DO ! k
         END DO ! i
         Last = ifin                                   ! re-set Last
      ENDIF  ! ifin

      ifind = 1 + (Time+Delt)*(Mtime-1)/Tlim           ! Vib distribution final time-bin index
      IF ( ifind.GT.Mtime ) ifind = Mtime
  
      IF ( ifind.GT.Dlast ) THEN                       ! snapshot at end of each time interval
         DO i = Dlast + 1 , Mtime                      ! ifind is LE Mtime (see earlier lines)
           Mum(Mol,i) = Mum(Mol,i) + 1                 ! Number of times ith time bin accessed
           EDIST(Mol,N,i) = EDIST(Mol,N,i) + 1.d+00    ! Energy-dependent vibrational distribution (N is index of energy)
         END DO ! i
         Dlast = ifind                                 ! re-set Dlast
      ENDIF  ! ifind

c--------This block for sampling molecular energies    <------------------------- FOR ENERGIES
c      *** Normally, commented-out ***
c         kfin = 1 + (Time+Delt)*(kstep-1)/Tlim          ! time-step index for collisional sampling
c         IF ( kfin.GT.klast ) THEN                      ! snapshot at end of each time interval
c            DO i = klast + 1 , MIN(kfin,kstep)
c               WRITE (20,*) E
c            END DO
c            klast = kfin
c         ENDIF
c------End energy sampling block-------------

c
c       Now apply the stochastic changes
c 
      Time = Time + Delt
      IPtot = 1 + Nchan(Mol) + nsmax(Mol)              ! = Collisions + reaction channels + supplementary reactions of Mol

      IF ( Igo.EQ.0 ) THEN                              ! Stop stepping
         RETURN
      ENDIF 
      
      IF ( IP .EQ. IPtot ) THEN                          ! *** COLLISION ***
         Step = COLSTEP(Mol,E,Temp)                     ! Collision step size
         E = E + Step                                   ! Energy after collision
c         write(21,*) Step   !         <----------------un-comment this line FOR WRITING COLLISION STEP-SIZES (in Unit 21)

      ELSEIF ( IP .LE. Nchan(Mol) ) THEN                      ! *** UNIMOLECULAR REACTION ***
         IF ( Jto(Mol,IP) .LE. Nwells ) THEN
           Enew = E + HMol(Mol) - HMol(Jto(Mol,IP))
           Mol = Jto(Mol,IP)                           ! New WELL index
           E = Enew                                    ! Energy of new Well
         ELSEIF ( Jto(Mol,IP) .GT. Nwells ) THEN
           Mol = Jto(Mol,IP)                           ! New PRODUCT index
         ENDIF  ! Jto
         
      ELSEIF ( (IP .GT. Nchan(Mol) ) .AND. (IP .LT. IPtot)  ) THEN       ! *** SUPPLEMENTARY REACTION ***
         IPP = IP - Nchan(Mol)                            ! channel number for supplementary reaction
         IF ( sto(Mol,IPP) .LE. Nwells ) THEN
           ddel = dde( Mol,IPP )
           Enew = E + HMol( Mol ) - HMol( sto(Mol,IPP) )
           Mol = sto(Mol,IPP)                          ! replace Mol with NEW WELL index
           E = Enew - ddel                             ! Energy of new Well, after subtracting quenching energy
         ELSEIF ( sto(Mol,IPP) .GT. Nwells ) THEN
           Mol = sto(Mol,IPP)                          ! New irreversible PRODUCT index
         ENDIF
      ENDIF  ! IP
     
      IF ( Mol.GT.NWells .AND. Last.LT.Ntime ) THEN    ! If irreversible product
         Igo = 0                                       ! Stop stepping
         DO  i = Last + 1 , Ntime
            Num(Mol,i) = Num(Mol,i) + 1                ! Number of times ith bin accessed
         END DO
      ENDIF

      RETURN
      END SUBROUTINE book1

!  ---------------------------------------------------------------------
      SUBROUTINE INITIAL
c
c       Initialize bookkeeping arrays between trials
c
      USE declare_mod

      IMPLICIT NONE
      INTEGER:: i , j , k , n
      SAVE 
 
      DO 200 i = 1 , NWells
 
         DO 50 j = 1 , Ntime
            EMol(i,j) = 0.0D+00               ! average vibrational energy
            Num(i,j) = 0                      ! number of samples for molecular properties
            Nuni(i,j) = 0                     ! number of samples for reaction rates
            DO 20 n = 1 , Nchan(i)
               kuni(i,n,j) = 0.0D+00       ! average unimolecular rates
 20         CONTINUE
 50      CONTINUE
 
         DO 100 j = 1 , Ndist
            DO 60 k = 1 , Mtime
               EDIST(i,j,k) = 0 ! vibrational distribution
 60         CONTINUE
 100     CONTINUE
 
         DO 110 j = 1 , Mtime
            Mum(i,j) = 0
 110     CONTINUE

 200  CONTINUE
       DO 300 i = 1 , Jsize
         Nstart(i) = 0          ! distribution of selected starting energies
 300  CONTINUE
 
 
      IF ( NProds.GT.0 ) THEN
         DO 350 i = NWells + 1 , NWells + NProds
            DO 320 j = 1 , Ntime
               Num(i,j) = 0     ! number of samples for molecular properties
 320        CONTINUE
 350     CONTINUE
      ENDIF
 
      RETURN
      END SUBROUTINE INITIAL


! ---------------------------------------------------------------------------------------------------------------------
      REAL(8) FUNCTION rivr( Mol , i , IE , numbdens )
      USE declare_mod

      IMPLICIT NONE
      INTEGER, INTENT(IN) :: Mol          ! Well
      INTEGER, INTENT(IN) :: i            ! rxn channel
      INTEGER, INTENT(IN) :: IE           ! energy grain number based on imax1 and imax2
      REAL(8), INTENT(IN) :: numbdens     ! number density for collision-induced ivr

      REAL(8) ::  EE , colls , rativr , vfrq , rmax
      SAVE
c
c     vivr(Mol,i)   : characteristic reaction frequency (as in RRK theory) (cm-1)
c     vave(Mol,i)   : average molecular frequency (cm-1)
c     pivr(Mol,i)   : bimolecular collision rate constant for collision-induced IVR
c     tivr(Mol,i)   : IVR threshold energy (cm-1), measured from the  reaction threshold (Eor)
c     civr(Mol,i,j) : coefficients for 2nd order fit of 
c                     k(ivr) = civr(-,-,1) + civr(-,-,2)*EE + cirv(-,-,3)*EE*EE
c                     where EE = (E - Eor) is the Energy measured from the rxn threshold
c     iivr(Mol,i)   : flag to calculate IVR (iivr=1), or not (iivr=0)
c
      rivr = 1.0d+00
      IF ( iivr(Mol , i) .EQ. 0 ) THEN              ! Neglect IVR
         
         RETURN

      ELSE          ! Calculate IVR transmission coefficient

        rmax = 2.d+00*2.9979d+10*vave( Mol , i )    ! Upper limit to IVR rate

         IF ( IE.LE.imax1 ) THEN                    ! EE = (E - Eor) is the Energy measured from the rxn threshold
            EE = (IE-1)*Egrain1
           ELSE
            EE = (IE-imax1-1)*Egrain2
         ENDIF
c
c        Calculate IVR rate constant = f(EE)
c
         IF ( EE .GT. tivr( Mol , i ) ) THEN                          ! Energy must be .GE. IVR threshold energy
           rativr = civr( Mol , i , 1 ) +                             ! IVR rate constant = f(EE) [second order fit]
     &        (civr( Mol , i , 2 ) + civr( Mol , i , 3 )*EE )*EE
           IF ( rativr .LT. 0.0d+00 ) rativr = 0.0d+00
           IF ( rativr .GT. rmax )    rativr = rmax                   ! upper limit corresponding to vave(Mol,i) cm-1
         ELSE
           rativr = 0.0d+00
         ENDIF

         vfrq = vivr(Mol , i)*2.9979d+10                              ! Convert from cm-1 to frequency (s-1)
         colls = pivr( Mol , i )*numbdens                             ! Collision-induced IVR pseudo-first-order rate

         rivr = (rativr + colls) / ( rativr + colls + vfrq )          ! Wolynes-Leitner IVR transmission coefficient
      ENDIF

      IF (rivr.LT.1.0d-10) rivr = 1.0d-10           ! ensures that log(rivr) will not be NAN when rivr is very small
      
      RETURN
      END FUNCTION rivr
      
! ---------------------------------------------------------------------------------------------------------------------
           
      END MODULE bookstep_mod
      
