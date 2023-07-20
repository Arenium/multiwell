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

      MODULE subs_mod

!     The following procedures are contained in this module:
!
!     SUBROUTINE DensArray      Copyright 2001 - 2021 John R. Barker
!     SUBROUTINE DensArrayTS    Copyright 2001 - 2021 John R. Barker
!     SUBROUTINE RateArray(T)   Copyright 2001 - 2021 John R. Barker, Philip Stimac
!     SUBROUTINE COLNORM        Copyright 2001 - 2021 John R. Barker
!     FUNCTION COLSTEP          Copyright 2001 - 2021 John R. Barker
!     SUBROUTINE PSTEP          Copyright 2001 - 2021 John R. Barker
!     FUNCTION Estart           Copyright 2001 - 2021 John R. Barker
!     FUNCTION Etherm           Copyright 2001 - 2021 John R. Barker
!     SUBROUTINE QKinf          Copyright 2001 - 2021 John R. Barker
!     SUBROUTINE QKinftun       Copyright 2001 - 2021 John R. Barker, Philip Stimac
!     SUBROUTINE ECKART         Copyright 2001 - 2021 John R. Barker, Philip Stimac
!     SUBROUTINE rotunits       Copyright 2001 - 2021 John R. Barker
!     SUBROUTINE UnitTest       Copyright 2001 - 2021 John R. Barker

      CONTAINS
      
!  ---------------------------------------------------------------------
      SUBROUTINE DensArray
c      Generates densities of states in 'double' array
c      Options: Whitten-Rabinovitch, or read densities from external file
c      SMol = 0 : A product molecule and do not calculate densities
c      Note that File name and path in OPEN statements may need revision
c            for non-Macintosh

      USE declare_mod
      USE utili_mod

      IMPLICIT NONE
      INTEGER :: i , JJ , j , k
      REAL(8) :: A , B , Test1 , Test2
      CHARACTER(len=100) :: TITLE
      SAVE  
 
      WRITE(*,*) 'Entering Subroutine DensArray'
 
      DO 100 i = 1 , NWells
         Viblo(i) = 100.                        ! lowest vibrational frequency
            OPEN (12,FILE='DensData/'//
     &      MolName(i)(1:lenstr(MolName(i)))//'.dens',STATUS='OLD')

            READ (12,99030) dumcut                 ! Check to see if file format includes input data summary
            IF ( dumcut .EQ. cut ) THEN            ! Format starting with v.2008.1
               dumcut = 'x'
               DO WHILE (dumcut .NE. cut)
                 READ (12,99030) dumcut
               END DO
            ELSE
               BACKSPACE 12                        ! If format from prior to v.2008.1
            ENDIF

            READ (12,*) TITLE
            READ (12,*) TITLE
c
c                   Egrain1 , imax1 , Emax2 , Isize , Viblo(i)
            READ (12,*) A , JJ , B , k , Viblo(i)
            READ (12,*) TITLE
            Test1 = ABS(A-Egrain1)/1.D-06
            Test2 = ABS(B-Emax2)/1.D-06
            IF ( Test1.LT.1.0 .AND. JJ.EQ.imax1 .AND. Test2.LT.1.0 .AND. 
     &           k.EQ.Isize ) THEN              ! Test for proper array bounds, etc.
               DO 10 k = 1 , Isize
                  READ (12,*) j , A , Dens(i,k) , B
                  IF ( Dens(i,k).GT.0.0 ) THEN
                     Dens(i,k) = LOG(Dens(i,k)) ! Store natural log.
                  ELSE
                     Dens(i,k) = -70.
                  ENDIF
 10            CONTINUE
            ELSE                        ! If grain or energy limit does not match
               WRITE (KSTD,99002) MolName(i)
               STOP
            ENDIF
            CLOSE (12)
 100  CONTINUE
 
      WRITE(*,*) 'Leaving Subroutine DensArray'

      RETURN
 
c99001 FORMAT (A100)
99002 FORMAT (/'Energy grain or energy limit does not match: ',
     &        A10/'****** FATAL ERROR: Run Terminated *******')
99030 FORMAT (A46)

      END SUBROUTINE DensArray


!  ---------------------------------------------------------------------
      SUBROUTINE DensArrayTS(Mol,i)
c      same code as densarray.f except this version generates the double
c      array for the transition state to be used later in subroutine 
c      ECKART to calculate the tunneling-corrected sum of states for 
c      the transition state

      USE declare_mod
      USE utili_mod

      IMPLICIT NONE
      INTEGER, INTENT(IN) :: Mol        ! index of Well
      INTEGER, INTENT(IN) :: i          ! index of rxn channel
      INTEGER :: j , k
      REAL(8) :: A , B
      
      SAVE  

      statest = 0.1/Egrain1

      OPEN (12,FILE='DensData/'//
     &    TSname(Mol,i)(1:lenstr(TSname(Mol,i)))//'.dens',STATUS='OLD')
       DO 10 k = 1 , Isize
          READ (12,*) j , A , DensTS(i,k) , B
          IF ( DensTS(i,k) .GT. statest ) THEN
             DensTS(i,k) = LOG(DensTS(i,k)) ! Store natural log.
          ELSE
             DensTS(i,k) = -70.
          ENDIF
 10    CONTINUE
      CLOSE(12)
      RETURN
 
      END SUBROUTINE DensArrayTS

!  ---------------------------------------------------------------------
      SUBROUTINE RateArray(T)

c       KEYWORDS:
c
c         'NOCENT'   no centrifugal correction             [NCENT =0]
c         'CENT1'    1-D adiabatic rotor centrifugal       [NCENT= 1]
c         'CENT2'    2-D adiabatic rotor centrifugal       [NCENT= 2]
c         'CENT3'    3-D adiabatic rotor centrifugal       [NCENT= 3]
c         'CENTX'    LEGACY centrifugal correction         [NCENT=-2]
c                    (not recommended)
c
c         'NOREV'    neglecting the reverse reaction        [Jrev = Jr]
c         'REV'      calculating reverse rxn rate; 
c                    afterwards, Jrev=1 is changed to Jrev=2
c
c         'FAST'     neglecting limitations due to IVR
c         'SLOW'     including IVR limitations; a following line contains 
c                        'SLOW', pivr(Mol,i), vivr(Mol,i), vave(Mol,i), 
c                             tivr(Mol,i), civr(Mol,i,1), civr(Mol,i,2), 
c                             civr(Mol,i,3)
c
c         'NOTUN'    for neglecting tunneling
c         'TUN'      for including tunneling via unsymmetrical Eckart barrier;
c                       a following line contains: 'TUN', tunnel(Mol,i)
c
c          'ILT'     Inverse laplace transform method for k(E)
c          'SUM'     File containing sums of states (usually generated by Densum)
c          'RKE'     File containing k(E), generated externally
c          'REVRKE'  k(E) for reverse reaction calculated from k(E) file for forward rxn
c
c         WARNINGS:
c          WARN1    = warning that k(E) = 0 was imposed because Product Well has empty grain corresponding to E

      USE declare_mod
      USE utili_mod

      IMPLICIT NONE

      REAL(8), INTENT(IN) :: T

      INTEGER i , k , Mol , JJ , jdum
      REAL(8) E , EE , DD , A , B , DUM , ER , Test1 , Test2 ,
     &        cfac , f , rat , Eto , Dto , psum , absbar
      CHARACTER TITLE*100
      SAVE 

      WRITE(*,*) 'Entering Subroutine RateArray'

      statest = 0.1/Egrain1        ! statest: Density criterion to test if a state is present in a given energy bin 

      DO Mol = 1 , NWells
            DO i = 1 , Nchan(Mol)
 
               Path(Mol,i) = (Molsym(Mol)*TSele(Mol,i)*TSopt(Mol,i))
     &                       /(TSsym(Mol,i)*Molele(Mol)*Molopt(Mol))
                          ! Rxn path degeneracy X electronic partition fxn ratio used for all k(E)

               IF ( (MolMom(Mol) .GT. 0.0 )
     &            .AND. (TSmom(Mol,i) .GT. 0.0)  )  THEN
                  rat = TSmom(Mol,i)/MolMom(Mol)
                 ELSE
                  rat = 1.0D+00
               ENDIF

               IF ( NCENT(Mol,i) .EQ. 0 ) THEN
                  cfac = 1.0d+00
                  Eor(Mol,i) = Eo(Mol,i)
               ELSEIF ( NCENT(Mol,i).GT.0 .AND. NCENT(Mol,i).LE.3) THEN 
                  f = 0.5*NCENT(Mol,i)*( 1.0d+00 - rat )      ! takes into account the number of adiabatic rotors
                  Eor(Mol,i) = Eo(Mol,i) + f*T/hor
                  cfac = rat*EXP(f)
                  IF ( NCENT(Mol,i) .EQ. -2 ) cfac = 1.d+00
               ENDIF

               IF ( ( Eor(Mol,i) .LT. 0.0 ) .AND. 
     &              ( Jread(Mol,i) .EQ. 1 )  ) THEN
                 write(*,*) ' '
                 write(*,*) '********** WARNING **********'
                 write(*,*) ' Eocent should not be <0 for SUM input '
                 write(*,*) ' Centrifugal correction is too large.'
                 write(*,*) '   Reaction from ....... to '
                 write(*,*) '   ', Mol, Jto(Mol,i)
                 write(*,*) '********** WARNING **********'
                 write(KOUT,*) ' '
                 write(KOUT,*) '********** WARNING **********'
                 write(KOUT,*) ' Eocent should not be <0 for SUM input '
                 write(KOUT,*) ' Centrifugal correction is too large'
                 write(KOUT,*) '   Reaction from ....... to '  
                 write(KOUT,*) '   ', Mol, Jto(Mol,i)
                 write(KOUT,*) '********** WARNING **********'
                 write(KSMM,*) ' '
                 write(KSMM,*) '********** WARNING **********'
                 write(KSMM,*) ' Eocent should not be <0 for SUM input '
                 write(KSMM,*) ' Centrifugal correction is too large'
                 write(KSMM,*) '   Reaction from ....... to '   
                 write(KSMM,*) '   ', Mol, Jto(Mol,i)
                 write(KSMM,*) '********** WARNING **********'
                 write(KSMM,*) ' '
               ENDIF                                    ! ++++++++++++++END STANDARD CENTRIFUGAL CORRECTION++++++++++++

               IF ( Jread(Mol,i).EQ.0 ) THEN            ! ++++++++++++++Inverse-Laplace-Transform++++++++++++++
                   
                  DO 5 k = 1 , Isize
                     IF ( k.LE.imax1 ) THEN             ! E = energy from top of barrier
                        E = (k-1)*Egrain1
                     ELSE
                        E = (k-imax1-1)*Egrain2
                     ENDIF
                     EE = E + Eo(Mol,i)                 ! Does NOT include centrifugal correction
                     DD = fxnexp(EE, Mol, Dens)
                     IF ( EE .LE. Emax2 .AND. DD.GT.statest ) THEN   ! calculate k(E) only up to Emax2, since dens is undefined for higher E
                       Rate(Mol,i,k) = Afac(Mol,i)*
     &                                     fxnexp(E, Mol, Dens)/DD    ! k(E)
                       Rate(Mol,i,k) = LOG(Rate(Mol,i,k))             ! k(E) stored as log
                     ELSE
                       Rate(Mol,i,k) = LOG( 1.e-20 )                  ! k(E) stored as log
                     ENDIF
                       
 5                CONTINUE
 
               ELSEIF ( Jread(Mol,i).EQ.1 ) THEN        ! ++++++++++++++Sums of states from external file++++++++++++++
                  OPEN (12,FILE='DensData/'//
     &            TSname(Mol,i)(1:lenstr(TSname(Mol,i)))//'.dens',
     &            STATUS='OLD')                         ! TST sum of states

                  READ (12,99030) dumcut                ! Check to see if file format includes input data summary
                  IF ( dumcut .EQ. cut ) THEN           ! Format starting with v.2008.1
                    dumcut = 'x'
                    DO WHILE (dumcut .NE. cut)
                      READ (12,99030) dumcut
                    END DO
                 ELSE
                    BACKSPACE 12                        ! If format from prior to v.2008.1
                 ENDIF
                 READ (12,*) TITLE
                 READ (12,*) TITLE
                 READ (12 , * )  A , JJ , B , k , DUM
                 READ (12,*) TITLE
                 Test1 = ABS(A-Egrain1)/1.D-06
                 Test2 = ABS(B-Emax2)/1.D-06
                 IF ( Test1.LT.1.0 .AND. JJ.EQ.imax1 .AND.
     &              Test2.LT.1.0 .AND. k.EQ.Isize ) THEN        ! Test for proper array bounds, etc.

                      absbar=HMol(Mol)+Eor(Mol,i)   ! absolute energy of barrier needed to determine
                                                    ! back reaction barrier

                    V1(Mol,i)=absbar-HMol(Jto(Mol,i))   ! HMol(Jto(Mol,i)) is the enthalpy of the product


                    if(itun(Mol,i).eq.1)then     ! check to make sure the back reaction barrier is
                                                 ! greater than zero, which must be the case if tunneling
                                                 ! is to be considered
                      if(V1(Mol,i).lt.0)then
                         write(6,99005)Jto(Mol,i),Mol
                         itun(Mol,i)=0
                      end if
                    end if

      if(itun(Mol,i) .eq. 0) then    ! ++++++++++++++ NO tunneling section ++++++++++++++

         DO k = 1 , Isize
           READ (12,*) jdum , A , B , TSsum(Mol,i,k)    ! Read sum of states  (jdum, A, and B are dummy variables)
         END DO  ! k
         
      else                        ! ++++++++++++++ TUNNELING section *************************************

        write(*,99003) Mol        ! Entering tunneling routine write statement

        CALL DensArrayTS(Mol,i)   ! Reads input file TSname to get densTS array so that the DENSITY (not TSsum) of TS states
                                  ! of the transition state can be interpolated later
        do k=1,Isize+imax1        ! loop over energies for which k(E) will be calculated

          IF ( k .LE. 2*imax1 ) THEN                
            E = (k-1)*Egrain1                      ! energy relative to (imax1-1)*Egrain1 below the barrier top
            EE = Eor(Mol,i) - (imax1-1)*Egrain1 + E  ! relative to the bottom of the Well 
          ELSE 
            E = (k-imax1*2-1)*Egrain2
            EE = Eor(Mol,i) + E  ! relative to bottom of the well ; ensures 2nd part
          ENDIF                  ! of double array starts at 0 cm-1, the top of the barrier

          if((Jrev(Mol,Jto(Mol,i))).eq.1)then  ! if the reverse rxn from well Mol to product Jto(Mol,i)
                                               ! is considered than store the density of states of the 
                                               ! product calculated for the forward rxn in TSdensBK  and 
                                               ! use this quantity to calculate the reverse rxn rate 

             iset(Mol)=Jto(Mol,i)              ! the reverse rxn will have iset(Mol) as the well and 
             iset2(Mol)=Mol                    ! iset2(Mol) as the product

          else if ( (Mol .eq. iset(Mol)) .and.
     &          ( Jto(Mol,i) .eq. iset2(Mol) ) ) then

             if ( TSdensBK(iset2(Mol),k) .lt. 1e-30 ) then
                TRate(Mol,i,k)=-70.0d0       
                goto 19                      
             end if

             if( fxnexp(EE, Mol, Dens) .lt. statest) then
                TRate(Mol,i,k)=-70.0d0
                goto 19
             end if

             IF ( EE .LE. Emax2 ) THEN                       ! If molecular energy is not too high, claculate k(E)
                TRate(Mol,i,k)=TSdensBK(iset2(Mol),k)*
     &            clight*Path(Mol,i)/fxnexp(EE, Mol, Dens)   !    k(E) with tunneling correction
                if ( TRate(Mol,i,k) .GT. 1.e-30 ) then
                    TRate(Mol,i,k) = LOG(TRate(Mol,i,k))     !    log( k(E) )
                else
                     TRate(Mol,i,k) = -70.
                endif
                goto 19
             ELSE
                TRate(Mol,i,k) = -70.0                  ! Do not calculate k(E) if molecular energy >Emax2
             ENDIF
          end if 

          CALL ECKART(Mol,i,EE,k,psum)       ! calculate modified sum of states including tunneling correction (unsymmetrical Eckart potential)
              
          if((Jrev(Mol,Jto(Mol,i))).eq.1)then   ! if the reverse rxn  Jto(Mol,i)-->Mol is considered 
             TSdensBK(Mol,k)=psum               ! save the denisty of states of the TS for the forward
          end if                                ! rxn and use them to calculate the reverse rxn rate, i.e.
                                                ! don't perform the calculation twice

          if(psum.lt.1e-30)then         ! this if statement is necessary because if the sum of 
             TRate(Mol,i,k)=-70.0d0     ! is zero taking LOG(TRate(Mol,i,k) below results in NAN 
             goto 19                    ! this circumvents that problem
          end if

          if( fxnexp(EE, Mol, Dens) .lt. statest) then
            TRate(Mol,i,k)=-70.0d0
            goto 19
          end if

          IF ( EE .LE. Emax2 ) THEN                    ! If molecular energy is not too high, claculate k(E)
             TRate(Mol,i,k) = psum*clight*Path(Mol,i) /
     &                         fxnexp(EE, Mol, Dens)           ! k(E) with tunneling correction
             TRate(Mol,i,k) = LOG(TRate(Mol,i,k))              ! k(E) stored as log
          ELSE
             TRate(Mol,i,k) = -70.0                    ! Do not calculate k(E) if molecular energy >Emax2
          ENDIF

19           CONTINUE
        end do  ! k

                  write(*,99004)
                  
        end if   ! end of no-tunneling/tunneling section ------------------------------------------------------------
           
         ELSE                               ! FATAL ERRORS if grain or energy limit does not match
           WRITE (KSTD,99002) TSname(Mol,i)
           WRITE (KOUT,99002) TSname(Mol,i)
           STOP
         ENDIF
         CLOSE (12)
 
               ELSEIF ( Jread(Mol,i).EQ.3 ) THEN        ! ++++++++++++++Read CRP file++++++++++++++
                 OPEN (12,FILE='DensData/'//
     &           TSname(Mol,i)(1:lenstr(TSname(Mol,i)))//'.crp',
     &           STATUS='OLD')                          ! CRP computed using SCTST of states
                  READ (12,99030) dumcut                ! Check to see if file format includes input data summary
                  IF ( dumcut .EQ. cut ) THEN           ! Format starting with v.2008.1
                    dumcut = 'x'
                    DO WHILE (dumcut .NE. cut)
                      READ (12,99030) dumcut            ! Check to see if file format includes input data summary
                    END DO
                 ELSE
                    BACKSPACE 12                        ! If format from prior to v.2008.1
                 ENDIF
                 READ (12,*) TITLE
                 READ (12,*) TITLE
                 READ (12 , * )  A , JJ , B , k , DUM
                 READ (12,*) TITLE
                 Test1 = ABS(A-Egrain1)/1.D-06
                 Test2 = ABS(B-Emax2)/1.D-06
                 IF ( Test1.LT.1.0 .AND. JJ.EQ.imax1 .AND. ! Test for proper array bounds, etc.
     &              Test2.LT.1.0 .AND. k.EQ.Isize ) THEN

                       DO k = 1 , Isize
                         READ (12,*) jdum , A , B , TSsum(Mol,i,k)   ! Read sum of states  (here, A and B are dummy variables)
                       END DO   ! k
                       
                 ENDIF
                 CLOSE (12)
                 

               ELSEIF ( Jread(Mol,i).EQ.2 ) THEN        ! ++++++++++++++Read k(E)'s from external file ++++++++++++++ (not including path degeneracy) 
                  OPEN (12,FILE='DensData/'//
     &            TSname(Mol,i)(1:lenstr(TSname(Mol,i)))//'.rke',
     &            STATUS='OLD')                   ! k(E)'s
                  READ (12,*) TITLE
                  READ (12,*) TITLE
c
c                           Egrain1 , imax1 , Emax2 , Isize
                  READ (12 , * )  A , JJ , B , k , DUM
                  READ (12,*) TITLE
                  Test1 = ABS(A-Egrain1)/1.D-06
                  Test2 = ABS(B-Emax2)/1.D-06
                  IF ( Test1.LT.1.0 .AND. JJ.EQ.imax1 .AND. 
     &                 Test2.LT.1.0 .AND. k.EQ.Isize ) THEN     ! Test for proper array bounds, etc.
                     DO k = 1 , Isize
                        READ (12,*) jdum , A , B , Rate(Mol,i,k)    ! Read k(E)'s
                        IF ( Afac(Mol,i).lt.0 ) Rate(Mol,i,k) =     ! Pseudo-first-order bimolecular rate constant
     &                    Rate(Mol,i,k)*ABS(Afac(Mol,i))            !  concentration [Reactant] = -Afac; thus, : kI(E) = k(E)*[Reactant] = -Afac*k(E)
                        IF ( Rate(Mol,i,k).GT.1.D-20 ) THEN
                           Rate(Mol,i,k) = LOG(Rate(Mol,i,k))       ! k(E) stored as log
                        ELSE
                           Rate(Mol,i,k) = LOG(1.D-20)              ! k(E) stored as log
                        ENDIF
                     END DO   ! k
                  ELSE                  ! If grain or energy limit does not match
                     WRITE (KSTD,99002) TSname(Mol,i)
                     WRITE (KOUT,99002) TSname(Mol,i)
                     STOP
                  ENDIF
                  CLOSE (12)
 
               ELSEIF ( Jread(Mol,i).EQ.-2 ) THEN       ! ++++++++++++++Read k(E)'s for REVERSE rxn (not including path degeneracy) ++++++++++++++
                  OPEN (12,FILE='DensData/'//
     &               TSname(Mol,i)(1:lenstr(TSname(Mol,i)))//'.rke',
     &               STATUS='OLD')                                    ! open k(E)'s previously computed for what is now the REVERSE rxn
                  READ (12,*) TITLE
                  READ (12,*) TITLE
c
c                           Egrain1 , imax1 , Emax2 , Isize
                  READ (12 , * )  A , JJ , B , k , DUM
                  READ (12,*) TITLE
                  Test1 = ABS(A-Egrain1)/1.D-06
                  Test2 = ABS(B-Emax2)/1.D-06
                  IF ( Test1 .LT. 1.0 .AND. JJ .EQ. imax1 .AND. 
     &                 Test2 .LT .1.0 .AND. k .EQ. Isize ) THEN
                                                ! Test for proper array bounds, etc.
                     DO k = 1 , Isize
                        READ (12,*) jdum , A , B , Rate(Mol,i,k)  ! Read k(E)'s for REVERSE reaction (not including path degeneracy) and use them to fill this array for the forward rxn
 
                        IF ( k.LE.imax1 ) THEN                  ! E = energy from top of barrier
                           E = (k-1)*Egrain1                    ! Energy in forward direction
                        ELSE
                           E = (k-imax1-1)*Egrain2              ! Energy in forward direction
                        ENDIF
                        E = E + Eor(Mol,i)                      ! With centrifugal correction
                        ER = E - (HMol(Jto(Mol,i))-HMol(Mol))   ! Energy in reverse direction
                        DD = fxnexp(ER, Jto(Mol,i), Dens)/
     &                                    fxnexp(E, Mol, Dens)  ! Density of states ratio
                        Rate(Mol,i,k) = Rate(Mol,i,k)*DD*Path(Mol,i)
                                                                ! now use DD and Path degeneracy to obtain k(E) in forward direction
                        IF ( Rate(Mol,i,k).GT.1.D-20 ) THEN
                           Rate(Mol,i,k) = LOG(Rate(Mol,i,k))
                                                        ! k(E) stored as log
                        ELSE
                           Rate(Mol,i,k) = LOG(1.D-20)  ! k(E) approximately =0 stored as log
                        ENDIF

                     END DO  ! k
                  ELSE                  ! If grain or energy limit does not match
                     WRITE (KSTD,99002) TSname(Mol,i)
                     WRITE (KOUT,99002) TSname(Mol,i)
                     STOP
                  ENDIF
                  CLOSE (12)
 
               ENDIF
               
               IF ( Jread(Mol,i).EQ.1 .OR. Jread(Mol,i).EQ.3 ) THEN      ! ++++++++ COMPUTE RATE CONSTANTS ++++++++
c              Jread     = 0  use inverse Laplace transform method for channel
c                        = 1  read sums of states for channel from external file
c                        = 2  read k(E) for this channel from external file
c                        = 3  read cumulative rxn prob from crp file
c                                                                        ! +++++++++IF STANDARD CENTRIFUGAL CORRECTION+++++
c                'NOCENT'   no centrifugal correction             [NCENT=  0]
c                'CENT1'    1-D adiabatic rotor centrifugal       [NCENT=  1]
c                'CENT2'    2-D adiabatic rotor centrifugal       [NCENT=  2]
c                'CENT3'    3-D adiabatic rotor centrifugal       [NCENT=  3]
c                'CENTX'    LEGACY centrifugal correction         [NCENT= -2]  (not recommended)

                   DO k = 1 , Isize
                     IF ( k.LE.imax1 ) THEN                    ! E = active energy in TS from top of barrier+ZPE
                       E = (k-1)*Egrain1                       ! Active energy in Transition State
                     ELSE
                       E = (k-imax1-1)*Egrain2                 ! Active energy in Transition State
                     ENDIF
                     EE = E + Eor(Mol,i)                       ! Active energy in REACTANT Well, includes centrifugal correction

                     IF ( EE.LE.Emax2 ) THEN                   ! calculate k(E) only up to Emax2, since Wells only have states up to Emax2
                       DD = fxnexp(EE, Mol, Dens)              ! density of states in reactant Well at energy EE
                       IF ( DD.GE.statest) THEN         
                         DD = TSsum(Mol,i,k)/DD 
                         Rate(Mol,i,k) = clight*DD*Path(Mol,i)     ! k(E) w/o tunneling
                         Rate(Mol,i,k) = cfac*Rate(Mol,i,k)        ! ++++++++ STANDARD CENTRIFUGAL CORRECTION +++++++
                         Rate(Mol,i,k) = LOG(Rate(Mol,i,k)+1.e-20) ! k(E) w/o tunneling correction stored as log
                       ELSE
                         Rate(Mol,i,k)= LOG( 1.d-20 )   ! If no state in this grain (empty grain),  rate = 0.
                       ENDIF
                     ELSE
                       Rate(Mol,i,k)= LOG( 1.d-20 )     ! if EE .GE. Emax2
                     ENDIF
                   END DO   ! k
                 ENDIF  ! Jread

c    FOR WELL-TO-WELL RXNS, CHECK TO MAKE SURE A STATE IS PRESENT IN THE PRODUCT WELL, regardless of how k(E) was calculated.

               WARN1(Mol,i) = 0                                   ! For counting empty grains (that contain no states)
               IF ( Jto(Mol,i) .LE. NWells ) THEN
                   DO k = 1 , Isize
                     IF ( k.LE.imax1 ) THEN                       ! E = energy in TS from top of barrier+ZPE
                       E = (k-1)*Egrain1                          ! energy in Transition State
                     ELSE
                       E = (k-imax1-1)*Egrain2                    ! energy in Transition State
                     ENDIF
                     EE = E + Eor(Mol,i)                          ! energy in REACTANT Well, includes possible centrifugal correction
                     Eto = EE + HMol(Mol) - HMol(Jto(Mol,i))      ! energy of PRODUCT WELL

                     IF ( Eto.GE.0.d+00 .AND. EE.GE.0.d+00 ) THEN
                        IF ( EE.LE.Emax2 .AND. Eto.LE.Emax2 ) THEN      ! calculate k(E) only up to Emax2, since Wells only have states up to Emax2
                          DD = fxnexp(EE, Mol, Dens)                    ! density of states in reactant Well at energy EE
                          IF ( DD.GE.statest) THEN         
                             Dto = fxnexp( Eto, Jto(Mol,i), Dens )      ! Test to see if a state exists in PRODUCT WELL energy grain
                             IF ( Dto .LT. statest ) THEN
                               WARN1(Mol,i) = WARN1(Mol,i) + 1          ! count number of empty grains
                               Rate(Mol,i,k)= LOG( 1.d-21 )             ! IF NO STATE IS PRESENT in this PRODUCT grain, then set k(E) = 0
                             ENDIF
                          ENDIF
                        ELSE
                          Rate(Mol,i,k)= LOG( 1.d-22 )
                        ENDIF
                     ELSE
                        WARN1(Mol,i) = WARN1(Mol,i) + 1            ! count number of empty grains
                        Rate(Mol,i,k)= LOG( 1.d-23 )               ! IF Eto < 0 or EE < 0, then set k(E) = 0
                     ENDIF
                   END DO   ! k
               ENDIF        ! end of test for empty grains

           END DO  ! index i (rxn channels)
        END DO  ! index Mol (loop over Wells)
         
      WRITE(*,*) 'Leaving Subroutine RateArray'

        RETURN
 
c99001   FORMAT (A100)
99002   FORMAT (/'Energy grain or energy limit does not match: ',
     &        A10/'****** Run Terminated *******')
99003   FORMAT ('Entering Tunneling Subroutine: Molecule #',I2)
99030   FORMAT (A46)
99004   FORMAT ('Leaving Tunneling Subroutine')
99005   FORMAT (/'The reaction barrier from product',i3,' to well', i3,
     &          /' is less than the enthalpy of this product.',
     &   /'Calculations will proceed as if tunneling were not invoked.')
 
        END SUBROUTINE RateArray

!  ---------------------------------------------------------------------
      SUBROUTINE COLNORM(Mol,Temp)
c      Normalization for a down-step model given in FUNCTION PSTEP(E,EE,Mol,Temp,Hs)
c
c            Mol      = index of molecule
c            Temp      = translational temperature
c             E      = initial energy
c            EE      = final energy
c            Hs      = characteristic stepsize at energy E (calculated by PSTEP)
c
c      Note that COLLUP and CNORM are 'double arrays' with the first
c      imax1+1 elements starting at E = 0 with Egrain1 spacing; the remaining
c      elements are at Egrain2 spacing, starting at E = 0.
c
c      Coefficients stored in array DC
c      Probability of up-steps stored in COLLUP(MaxWell,Isize)
c      Normalization factor (up + down) stored in CNORM(MaxWell,Isize)
c
c       ERRo    = Convergence criterion declared in Declare3.inc include file
c      EDIFF       = Minimum integration interval (cm-1)
c      FRACT        = Step-size relative to 'Hs'
c      ITERS        = Iterations on up-steps

      USE declare_mod
      USE utili_mod
      
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: Mol            ! index number of Well
      REAL(8), INTENT(IN) :: Temp            ! Temperature
 
      REAL(8) :: RATIO , DELMIN, B , DENI , E , EE , Hs , PROB , H ,  
     &         DENJ , TERM , SPAN , TLAST , TEST , a ,  AA
      INTEGER :: I , Ilow , itemp
      REAL(8), DIMENSION(Isize) :: CN

      WRITE(*,*) 'Entering Subroutine COLNORM: Molecule #',Mol
 
      Statest = 0.5d+00/Egrain1      ! Criterion for whether a state is present within an energy grain

c      Integrate over UP steps (temporarily store up-step integrals in COLLUP)
c
       DO 100 I = 1 , Isize
         IF ( I.LE.imax1 ) THEN
            E = (I-1)*Egrain1
           ELSE
            E = (I-imax1-1)*Egrain2
         ENDIF

         DENI = fxnexp(E, Mol, Dens)                                  ! Density at E
                
         IF ( I. EQ. 1 ) DENI = 1.0                          ! Special case for Up-steps at E=0

         IF ( E.GE.0.0d+00 .AND. DENI.GE.Statest ) THEN      !  ONLY IF STATE, OTHERWISE SKIP
            TEST = 1.0
            COLLUP(Mol,I) = 0.0D+00     ! Temporary storage of UP-normalization
            EE = E
            CALL PSTEP(EE,E,Mol,Temp, PROB, Hs)
            TLAST = 1.0d+00                                     ! Initial Hs and 'old' TERM
            DELMIN = Viblo(Mol) + 0.5*Hs*(Temp/hor)/(Hs + (Temp/hor))

            DO WHILE ( TEST.GT.ERRo )                            ! ********** UP STEPS **********
               IF ( EE.GT.Emax1 ) THEN
                  H = FRACT*Hs*(Temp/hor)/(Hs + (Temp/hor))    
                ELSE
                  H = Egrain1
               ENDIF

               EE = EE + H                                     ! UP-steps

               DENJ = fxnexp(EE, Mol, Dens)
               
               IF ( DENJ .GE. Statest ) THEN
                  RATIO = DENJ/DENI
                  B = EXP(-(EE-E)*hor/Temp)                     ! EE > E
                  CALL PSTEP(EE,E,Mol,Temp, PROB, Hs)
                  TERM = PROB*B*RATIO
               ELSE
                  TERM = 0.0D+00
               ENDIF
               
               IF ( EE .GT. Emax1 ) THEN
                  a = -LOG(TERM/TLAST)/H                        ! Exponential form
                  AA = TLAST                                    ! Exponential form
                  SPAN = (AA/a)*(1.0D+00-EXP(-a*H))             ! Exponential form: integral
               ELSE
                  SPAN = 0.5D+00*H*( TERM + TLAST )             ! Trapezoidal Rule: integral
               ENDIF
 
               COLLUP(Mol,I) = COLLUP(Mol,I) + SPAN
c99     FORMAT(3(1pe10.3,3x))
c          IF ( I .EQ. 1 ) write(50,99) E,EE-E, COLLUP(Mol,I)  !  +++++++++++++++++++++++++++++++++++++++++++++
 
               IF ( SPAN.GT.0.0 .AND. ABS(E-EE).GT.DELMIN )
     &              TEST = ABS(SPAN/COLLUP(Mol,I))
               TLAST = TERM                                     ! 'old' TERM
            ENDDO  ! while
         ENDIF   ! if state present
 100  CONTINUE

      WRITE (*,*) '    Finished up steps, starting down-steps'
 
c            ! ********** DOWN STEPS ***********
c*************************************************************
c 
      Ilow = INT(Emax1/Egrain2) + 2 + imax1                 ! index of lowest element in upper half not overlapped with lower half

      CNORM(Mol,1) = COLLUP(Mol,1)
      CN(1) = 0.0d+00
 
        DO 250 I = 2 , Isize

        IF ( I.LE.imax1 ) THEN
           E = (I-1)*Egrain1
          ELSE
           E = (I-imax1-1)*Egrain2
        ENDIF

        DENI = EXP( Dens(Mol,I) )                   ! Density at E(I)

      IF ( DENI.LT.Statest ) THEN
          CNORM(Mol,I) = 0.0d+00
      ELSE                                        ! PROCEED ONLY IF STATE PRESENT

            EE = E
            CALL PSTEP(E,EE,Mol,Temp, PROB, Hs)                    ! to get Hs
            DELMIN = Viblo(Mol) + 0.5*Hs
            TLAST = 1.0d+00                                  ! Initial Hs and 'old' TERM
            TEST = 1.0
            CN(I) = 0.0d+00

            DO WHILE ( TEST.GT.ERRo .AND. EE .GE. Egrain1 )
              IF ( EE.GT.Emax1 ) THEN
                H = FRACT*Hs
                H = MIN( H , EE-Emax1 )
               ELSE
                H = Egrain1
              ENDIF
              IF ( H.GT.EE ) H = EE

              EE = EE - H                                  ! DOWN-steps
              DENJ = fxnexp(EE, Mol, Dens)
              IF ( DENJ.GE.Statest ) THEN
                CALL PSTEP(E,EE,Mol,Temp, PROB, Hs)                  
                TERM = PROB
               ELSE
                TERM = 0.0D+00
              ENDIF
 
              IF ( EE-H .GE. Emax1 ) THEN
                a = -LOG(TERM/TLAST)/H                     ! Exponential form
                AA = TLAST                                 ! Exponential form
                SPAN = (AA/a)*(1.0D+00-EXP(-a*H))          ! Exponential form
               ELSE
                SPAN = 0.5D+00*H*( TERM + TLAST )          ! Trapezoidal Rule
              ENDIF
              TLAST = TERM                          ! 'old' value of TERM
              CN(I) = CN(I) + SPAN
              
c               IF ( I .EQ. 1 ) write(50,99) E, EE-E, CN(I)  !  +++++++++++++++++++++++++++++++++++++++++++++

              IF ( SPAN.GT.0.0 .AND. ABS(E-EE).GT.DELMIN ) 
     &             TEST = ABS(SPAN/CN(I))

            ENDDO  ! while

            IF ( ( COLLUP(Mol,I)+CN(I) ) .GT. 0.0 ) THEN
              CNORM(Mol,I) = ( COLLUP(Mol,I) + CN(I) )         
             ELSE
              CNORM(Mol,I) = 0.0d+00
            ENDIF

        ENDIF           ! End test for state present
 
 250  CONTINUE          ! end energy steps

      CNORM(Mol,1) = COLLUP(Mol,1) 
      CNORM(Mol,imax1+1) = CNORM(Mol,1)
 
c      Now calculate final COLLUP array

      DO 500 I = 1 , Isize
        IF ( CNORM(Mol,I) .GT. 0.d+00 ) THEN
          COLLUP(Mol,I) = COLLUP(Mol,I)/CNORM(Mol,I)             ! Evaluate up-probability array
         ELSE
          COLLUP(Mol,I) = 0.0d+00
        ENDIF
 500  CONTINUE

c      Now use values for lower half of double array to replace overlapped values in upper half

      DO I = imax1+2 , Ilow-1
         EE = (I-imax1-1)*Egrain2
         itemp = 1 + NINT( EE/Egrain1 )
         CNORM(Mol,I)  = CNORM(Mol,itemp)
         COLLUP(Mol,I) = COLLUP(Mol,itemp)
      END DO

      COLLUP(Mol,1) = 1.0d+00
      COLLUP(Mol,imax1+1) = COLLUP(Mol,1)
      CNORM(Mol,imax1+1) = CNORM(Mol,1)

      WRITE(*,*) 'Leaving Subroutine COLNORM'
 
      RETURN
 
      END subroutine colnorm

!  ---------------------------------------------------------------------
      REAL(8) FUNCTION COLSTEP(Mol,E,Temp)
      
C     FROM RANDOM NUMBER RON AND ENERGY E COMPUTE A COLLISIONAL
C     STEP SIZE FROM THE NORMALIZED CUMULATIVE DISTRIBUTION 
C     FUNCTION CORRESPONDING TO THE COLLISIONAL DOWN-STEP MODEL:
C
C      IFLAG2 = -1 FOR DOWN-STEPS
C      IFLAG2 = +1 FOR UP-STEPS
C
c      Step-size distribution given in FUNCTION PSTEP(E,EE,Mol,Temp,Hs)
c
c            Mol      = index of molecule
c            Temp      = translational temperature
c             E      = initial energy
c            EE      = final energy
c            DC      = coefficients for stepsize distribution
c            Hs      = characteristic stepsize at energy E (calculated by PSTEP
c
c      Note that COLLUP and CNORM are 'double arrays' with the first
c      imax1+1 elements starting at E = 0 with Egrain1 spacing; the remaining
c      elements are at Egrain2 spacing, starting at E = 0.
c
c      Coefficients stored in array DC
c      Probability of up-steps stored in COLLUP(MaxWell,Isize)
c      Log of Normalization factor (up + down) stored in CNORM(MaxWell,Isize)
c
c      ERRo      = Convergence criterion (declared in "declare3.inc")
c      EDIFF       = Minimum integration interval (cm-1)
c      FRACT        = Step-size relative to 'Hs'
c      ITERS        = Iterations on up-steps
 
      USE declare_mod
      USE utili_mod

      IMPLICIT NONE
      INTEGER, INTENT(IN) :: Mol                  ! index number of species
      REAL(8), INTENT(IN) :: E, Temp                  ! Energy, Temperature
      REAL(8) :: CNORMI , 
     &                 RATIO , CD 
      REAL(8) :: UP , RON , TOT , ADOWN , AUP , B , a , 
     &                 AA , Ht , arg , DELMIN , TEST 

      INTEGER :: IFLAG2 , Igg , Ntrap

 
      REAL(8) :: DENI , EE , Hs , PROB , H , DENJ , TERM , 
     &          SPAN , TLAST , Tiggp , Tiggm
c      SAVE 

      IF ( E .GT. Emax2 ) THEN
        WRITE(*,*) '**********************************************'
        WRITE(*,*) ' FATAL ERROR at point #1 in colstep.f!'
        WRITE(*,*) ' E is greater than Emax2'
        WRITE(*,*) ' Suggest you increase Emax2'
        WRITE(*,*) '**********************************************'
        STOP 
      ENDIF

      Statest = 0.1/Egrain1

      DENI = fxnexp(E, Mol, Dens)                              ! DENI = Density at E
      CNORMI = fxnofe(E , Mol , CNORM)
      UP = fxnofe(E , Mol , COLLUP)                         ! Probability of an up-step 
      IF ( UP .GT. 1.0d+00 ) UP = 1.0d+00

      CALL PSTEP(E,E,Mol,Temp, PROB, Hs)     ! Evaluate to get initial Hs
      
      RON = RAN1(IDUM)
      IF ( RON.LE.UP ) THEN             ! Select up-steps, or down-steps
         IFLAG2 = 1                     ! Up-step
        ELSE
         IFLAG2 = -1                    ! Down-step
      ENDIF
      IF ( E .LE. 0.0 ) IFLAG2 = 1
 

c ******************** Down steps ***************
      IF ( IFLAG2.EQ.-1 ) THEN                  !  Integrate over down steps
 
         DELMIN = Viblo(Mol) + 0.5*Hs
         CD = CNORMI*(1.d+00-UP)
 50      RON = RAN1(IDUM)
         IF ( RON .GT. (1.0d+00-ERRo) ) RON = 1.00D+00 - ERRo
 
         ADOWN = RON*CD !  Random-number-selected integral of down-steps

         TLAST = 1.0d+00                        ! 'old' TERM
         TOT = 0.0
         EE = E
         TEST = 1.0 
 
         DO WHILE (TOT.LT.ADOWN .AND. EE.GE.Egrain1 .AND. TEST.GT.ERRo)
              IF ( EE.GT.Emax1 ) THEN
                H = FRACT*Hs
                H = MIN( H , EE-Emax1 )
              ELSE
                H = Egrain1
              ENDIF
              IF ( H.GT.EE ) H = EE                        ! for down-steps

              EE = EE - H                                  ! DOWN-steps

              DENJ = fxnexp(EE, Mol, Dens)
              IF ( DENJ.GE.Statest ) THEN
                CALL PSTEP(E,EE,Mol,Temp, PROB, Hs)                      
                TERM = PROB
               ELSE
                TERM = 0.0D+00
              ENDIF
 
              IF ( EE-H .GE. Emax1 ) THEN
                Ntrap = 0
                a = -LOG(TERM/TLAST)/H                     ! Exponential form
                AA = TLAST                                 ! Exponential form
                SPAN = (AA/a)*(1.0D+00-EXP(-a*H))          ! Exponential form
               ELSE
                Ntrap = 1
                SPAN = 0.5D+00*H*( TERM + TLAST )          ! Trapezoidal Rule
              ENDIF
              TLAST = TERM                          ! 'old' value of TERM
              TOT = TOT + SPAN                      ! Down-step normalization integral
 
            IF ( DENJ .GE. Statest ) THEN
               IF ( (E-EE) .GE. DELMIN ) TEST = SPAN/TOT
               IF ( TOT .GE. ADOWN ) TEST = 0.0
            ENDIF
            
         ENDDO

         IF ( TOT.LT.ADOWN ) THEN       ! Invoked when estimated CNORMI is too large
            CD = TOT*(1.0d+00-ERRo)     ! Use new cumulative sum
            GOTO 50                     ! Repeat search with a fresh random number
         ENDIF
c                                       Begin correction for overshoot ( TOT > ADOWN )

         IF ( Ntrap .EQ. 1) THEN               ! Trapezoidal rule
             EE = EE + (TOT - ADOWN)*H/SPAN
           ELSE                                ! Exponential form
             arg = 1.0D+00 - a*( ADOWN - TOT + SPAN )/AA
             Ht = -LOG( arg )/a
             EE = EE + H - Ht           ! Correction for overshoot             
          ENDIF
          
c ******************** Up steps ***************
      ELSE                                              ! Up-collisions
 
         CD = CNORMI*UP
 100     RON = RAN1(IDUM)
         IF ( RON .GT. (1.0d+00-ERRo) ) RON = 1.00D+00 - ERRo
         AUP = RON*CD                   ! Random number-selected
         CALL PSTEP(E,E,Mol,Temp, PROB, Hs)
         DELMIN = Viblo(Mol) + 0.5*Hs*(Temp/hor)/(Hs + (Temp/hor))
         TLAST = 1.0d+00                                ! 'old' TERM
         TOT = 0.0d+00
         EE = E
         TEST = 1.0
 
         DO WHILE ( TOT.LT.AUP .AND. TEST.GT.ERRo )
               IF ( EE.GE.Emax1 ) THEN
                  H = FRACT*Hs*(Temp/hor)/(Hs + (Temp/hor))
                ELSE
                  H = Egrain1
               ENDIF

               EE = EE + H                                     ! UP-steps

               DENJ = fxnexp(EE, Mol, Dens)
               IF ( DENJ .GE. Statest ) THEN
                  RATIO = DENJ/DENI
                  B = EXP(-(EE-E)*hor/Temp)                     ! EE > E
                  CALL PSTEP(EE,E,Mol,Temp, PROB, Hs)               
                  TERM = PROB*B*RATIO
                ELSE
                  TERM = 0.0D+00
               ENDIF
 
               IF ( EE .GT. Emax1 ) THEN
                  Ntrap = 0
                  a = -LOG(TERM/TLAST)/H                        ! Exponential form
                  AA = TLAST                                    ! Exponential form
                  SPAN = (AA/a)*(1.0D+00-EXP(-a*H))             ! Exponential form: integral
               ELSE
                  Ntrap = 1
                  SPAN = 0.5D+00*H*( TERM + TLAST )             ! Trapezoidal Rule: integral
               ENDIF
 
            TLAST = TERM                        ! Old value for TERM
            TOT = TOT + SPAN                    ! Up-step normalization integral
 
            IF ( DENJ .GE. Statest ) THEN
               IF ( (EE-E) .GE. DELMIN ) TEST = SPAN/TOT
               IF ( TOT .GE. AUP ) TEST = 0.0
            ENDIF

      IF ( EE .GT. Emax2 ) THEN
        WRITE(*,*) '**********************************************'
        WRITE(*,*) ' FATAL ERROR at point #2 in colstep.f!'
        WRITE(*,*) ' EE is greater than Emax2'
        WRITE(*,*) ' Suggest you increase Emax2'
        WRITE(*,*) '**********************************************'
        STOP 
      ENDIF

         ENDDO
 
         IF ( TOT.LT.AUP ) THEN                 ! Invoked when estimated CNORMI was too large
            CD = TOT*(1.0-ERRo)                 ! Use new cumulative sum
            GOTO 100                            ! Repeat search
         ENDIF
c                                       Begin correction for overshoot ( TOT > AUP )

         IF ( Ntrap .EQ. 1) THEN               ! Trapezpidal rule
             EE = EE - (TOT - AUP)*H/SPAN
           ELSE                                ! Exponential form
             arg = 1.0D+00 - a*( AUP - TOT + SPAN )/AA
             Ht = -LOG( arg )/a
             EE = EE - H +  Ht           ! Correction for overshoot
          ENDIF

      ENDIF  ! up or down steps
 
      IF ( EE.LT.0.0 ) EE = 0.0d+00

      IF ( EE .LE. Emax1) THEN               ! align with energy grains
        Igg = NINT( EE/Egrain1 )
        EE = Egrain1*Igg
        IF ( fxnexp(EE, Mol, Dens) .LT. Statest ) THEN  ! if a state is not present, move up or down one grain
           Tiggp = fxnexp(EE+Egrain1, Mol, Dens)
           Tiggm = fxnexp(EE-Egrain1, Mol, Dens)
           IF ( Tiggp .GT. Tiggm ) THEN
               EE = EE + Egrain1
             ELSE
               EE = EE - Egrain1
           ENDIF
        ENDIF

        IF ( EE.LT.0.0 ) EE = 0.0d+00

      ENDIF
        
      COLSTEP = EE - E
      
      RETURN

      END FUNCTION COLSTEP

!  ---------------------------------------------------------------------
      SUBROUTINE PSTEP(E,EE,Mol,TEMPI,PROB,Hs)
c
c      PROB      = step-size probability density for DOWN-steps:  E > EE, not normalized,
c      E      = initial energy
c      EE      = final energy
c      Mol      = index number for a molecule
c      TEMPI      = translational temperature
c      Hs      = characteristic 'average' step size corresponding to PROB;
c            used in SUBROUTINEs COLNORM and COLSTEP to determine integration
c            step-sizes.
c
c      ITYPE (stored in a COMMON) is a flag used to designate collision model:
c
c      ITYPE = 1  for Biexponential Model
c      ITYPE = 2  for Density-weighted Biexponential Model
c      ITYPE = 3  for Off-set Gaussian with constant offset and E-dependent width
c      ITYPE = 4  for Biexponential Model with energy-dependent fraction
c      ITYPE = 5  for Generalized Gaussian with energy-dependent exponent
c      ITYPE = 6  for Generalized Gaussian plus Exponential term
c      ITYPE = 7  for Weibull Model
c      ITYPE = 8  for Lorentzian Step-Ladder Model
c      ITYPE = 9  for Exponential+Elastic Model
c      ITYPE = 10 for Klaus Luther's empirical function
c      ITYPE = 11 for Triplet state empirical function
c      ITYPE = 12 for Exponential Model with Alpha(E)=linear + exponential
c      ITYPE = 13 for Exponential Model with alpha(E) switching function
c      ITYPE = 14 for Density-weighted Exponential Model, motivated by Barker & Weston [J.Phys.Chem. A, 114, 10619-10633 (2010)]

      USE declare_mod
      USE utili_mod
 
      IMPLICIT NONE
      REAL(8), INTENT(IN) :: E
      REAL(8), INTENT(IN) :: EE
      REAL(8), INTENT(IN) :: TEMPI
      INTEGER, INTENT(IN) :: Mol
      REAL(8), INTENT(OUT) :: Hs
      REAL(8), INTENT(OUT) :: PROB
      REAL(8) :: tfac
      REAL(8) :: B , DENE , DENEE
      REAL(8) :: TERM1 , TERM11 , TERM2 , TERM22 , ALPHA , ALP ,
     &           ALPHA1 , ALPHA2 , EN , Fractt , A , BET , BETA , 
     &           X , ARG , ARG1 , ARG2 , WIDTH , F , TT , Rob_3
 
      IF ( E.LE.-0.01 ) THEN     ! Indicates that we just want the title lines
 
         LI(Mol,1) = '  '
         LI(Mol,2) = '  '
         LI(Mol,3) = '  '
         LI(Mol,4) = '  '
         LI(Mol,5) = '  ' 

         IF ( ITYPE(Mol).EQ.1 ) THEN
         LI(Mol,1) = ' Biexponential Model'
         LI(Mol,2) = ' Alpha1= [C(1) + E*C(2) + E*E*C(3)]*(T/300)**C(8)'
         LI(Mol,3) = ' Alpha2= [C(5) + E*C(6) + E*E*C(7)]*(T/300)**C(8)'
         LI(Mol,4) = 
     &' Pdown = (1-C(4))*EXP(-(E-EE)/Alpha1) + C(4)*EXP(-(E-EE)/Alpha2)'
         
       
       ELSEIF ( ITYPE(Mol).EQ.2 ) THEN
               LI(Mol,1) = ' Density-weighted Biexponential Model'
         LI(Mol,2) = ' Alpha1= [C(1) + E*C(2) + E*E*C(3)]*(T/300)**C(8)'
         LI(Mol,3) = ' Alpha2= [C(5) + E*C(6) + E*E*C(7)]*(T/300)**C(8)'
               LI(Mol,4) = 
     &' Pdown = rho*((1-C(4))*EXP(-(E-EE)/Alpha1) + C(4)*EXP(-(E-EE)/Alp
     &ha2))'
       
       ELSEIF ( ITYPE(Mol).EQ.3 ) THEN
               LI(Mol,1) = 
     &                ' Constant off-set, variable width Gaussian Model'
               LI(Mol,2) = 
     &              '  Alpha1= [C(1) + E*C(2) + E*E*C(3)]*(T/300)**C(8)'
               LI(Mol,3) = 
     &              '  C(4) = constant off-set, Alpha1 is the std. dev.'
               LI(Mol,4) = ' Pdown = EXP(-(0.5*(E-EE-C(4))/Alpha1)**2)'
       
       ELSEIF ( ITYPE(Mol).EQ.4 ) THEN
               LI(Mol,1) = 
     &             ' Biexponential Model with Energy-Dependent Fraction'
               LI(Mol,2) = '  Alpha1 = C(1) + [E*C(2) + E*E*C(3)]'
               LI(Mol,3) = '  Alpha2 = C(5) + [E*C(6) + E*E*C(7)]'
               LI(Mol,4) = 
     &' Pdown = (1-[C(4)+E*C(8)])*EXP(-(E-EE)/Alpha1) + [C(4)+E*C(8)]*EX
     &P(-(E-EE)/Alpha2)'
       
       ELSEIF ( ITYPE(Mol).EQ.5 ) THEN
               LI(Mol,1) = 
     &            'Generalized Gaussian Model with E-Dependent Exponent'
               LI(Mol,2) = '  Alpha = C(1) + [E*C(2) + E*E*C(3)]'
               LI(Mol,3) = '  exponent = C(5) + [E*C(6) + E*E*C(7)]'
               LI(Mol,4) = ' Pdown = EXP(-[(E-EE)/Alpha]**Exponent)'
      
      ELSEIF ( ITYPE(Mol).EQ.6 ) THEN
               LI(Mol,1) = ' Generalized Gaussian + Exponential'
               LI(Mol,2) = 
     &'  Alpha1   = C(1) + [E*C(2) + E*E*C(3)]    Alpha2 = C(7)+ E*C(8)'
               LI(Mol,3) = '  Exponent = C(4) + E*C(5)'
               LI(Mol,4) = '  F        = C(6)'
               LI(Mol,5) = 
     &' PDOWN     = (1-C(6))*EXP(-[(E-EE)/Alpha1]**Exponent) + C(6)*EXP(
     &-(E-EE)/Alpha2)'
       
       ELSEIF ( ITYPE(Mol).EQ.7 ) THEN
               LI(Mol,1) = ' Weibull Model'
               LI(Mol,2) = 
     &                    '  Alpha = {C(1) + [E*C(2) + E*E*C(3)]}**Beta'
               LI(Mol,3) = '  Beta  = C(4) + E*C(5)'
               LI(Mol,4) = 
     &' PDOWN  = (Beta/Alpha)*{(E-EE)**(beta-1)}*exp(-((E-EE)/Alpha)**Be
     &ta)'
      
      ELSEIF ( ITYPE(Mol).EQ.8 ) THEN
               LI(Mol,1) = ' Lorentzian Step-Ladder Model'
               LI(Mol,2) = '  Alpha = C(1) + [E*C(2) + E*E*C(3)]'
               LI(Mol,3) = '  Width = C(5) + [E*C(6) + E*E*C(7)]'
               LI(Mol,4) = ' Pdown = 1 / [(E-EE-Alpha)^2 + Width^2]'
      
      ELSEIF ( ITYPE(Mol).EQ.9 ) THEN
               LI(Mol,1) = ' Exponential+Elastic Model'
               LI(Mol,2) = 
     &              '  Alpha = [C(1) + E*C(2) + E*E*C(3)]*(T/300)**C(8)'
               LI(Mol,3) = 
     &      '      F = [C(5) + (E/C(6))**C(7)]*(T/300)**C(8); when E=EE'
               LI(Mol,4) = 
     &              ' Pdown = [F/(F+C(4))]*EXP(-(E-EE)/Alpha) + elastic'
       
       ELSEIF ( ITYPE(Mol).EQ.10 ) THEN 
               LI(Mol,1) = ' Klaus Luther s P(EE,E) Empirical Function'
               LI(Mol,2) = '  Alpha = C(1) + E*C(2)'
               LI(Mol,3) = '  Beta  = C(3)'
               LI(Mol,4) = ' Pdown = EXP {-[(E-EE)/Alpha]^Beta}'
       
       ELSEIF ( ITYPE(Mol).EQ.11 ) THEN 
               LI(Mol,1) = ' Triplet State P(EE,E) Empirical Function'
               LI(Mol,2) = 
     &                   '  Alpha = C(1)*{1-exp[-(E/C(2))^C(3)] + C(4)}'
               LI(Mol,3) = '  PROB = EXP [-(E-EE)/Alpha]'
               LI(Mol,4) = ' '

       ELSEIF ( ITYPE(Mol).EQ.12 ) THEN 
               LI(Mol,1) = 
     &             'Exponential Model with alpha(E)= Linear+exponential'
               LI(Mol,2) = ' Alpha1 = C(1) + E*C(2) + C(3)*exp(-E/C(4))'
               LI(Mol,3) =  ' Pdown = EXP(-(E-EE)/Alpha1)'

       ELSEIF ( ITYPE(Mol).EQ.13 ) THEN 
               LI(Mol,1) = ' Exponential Model with Switching Function'
               LI(Mol,2) = 
     &           '  Alpha1 = C(1) + E*C(2);      Alpha2 = C(3) + E*C(4)'
               LI(Mol,3) = 
     &' Alp = Alpha1 + 0.5*(Alpha2 - Alpha1)*(1. - TANH((C(5)-E)/C(6)))'
               LI(Mol,4) =  ' Pdown = EXP(-(E-EE)/Alp)'

      ELSEIF ( ITYPE(Mol).EQ.14 ) THEN
        LI(Mol,1) = ' Boltzmann-weighted, after Barker & Weston [2010]'
               LI(Mol,2) = ' Pdown = B(T;EE,E)*EXP(-(E-EE)/Alp)'
               LI(Mol,3) = ' B(T;EE,E) = SQRT(rho(EE)*exp(-(EE-E)/RT))'
               LI(Mol,4) = ' Alp = C(1) ; rho(EE) = density of states '
               
          ENDIF  ! ITYPE

       ELSE    ! we want to calculate probability
 
       IF ( ITYPE(Mol).EQ.1 ) THEN
         tfac = (TEMPI/300.0)**DC(Mol,8)
         ALPHA1 = ( ABS(DC(Mol,1) + E*( DC(Mol,2)+E*DC(Mol,3) ) )*tfac)
         ALPHA2 = ( ABS(DC(Mol,5) + E*( DC(Mol,6)+E*DC(Mol,7) ) )*tfac)
  
         IF ( ALPHA1.GT.0.0 ) THEN
            TERM1 = EXP(-ABS(E-EE)/ALPHA1)
            TERM11 = TERM1/ALPHA1
         ELSE
            TERM1 = 0.0D+00
            TERM11 = 0.0D+00
         ENDIF
  
         IF ( ALPHA2.GT.0.0 ) THEN
            TERM2 = EXP(-ABS(E-EE)/ALPHA2)
            TERM22 = TERM2/ALPHA2
         ELSE
            TERM2 = 0.0D+00
            TERM22 = 0.0D+00
         ENDIF
  
         PROB = (1.0D00-DC(Mol,4))*TERM1 + DC(Mol,4)*TERM2
         Hs = PROB/((1.0D00-DC(Mol,4))*TERM11+DC(Mol,4)*TERM22)

         
       ELSEIF ( ITYPE(Mol).EQ.2 ) THEN
         tfac = (TEMPI/300.0)**DC(Mol,8)
         ALPHA1 = ( ABS(DC(Mol,1) + E*( DC(Mol,2)+E*DC(Mol,3) ) )*tfac)
         ALPHA2 = ( ABS(DC(Mol,5) + E*( DC(Mol,6)+E*DC(Mol,7) ) )*tfac)
 
            IF ( ALPHA1.GT.0.0 ) THEN
               TERM1 = EXP(-ABS(E-EE)/ALPHA1)
               TERM11 = TERM1/ALPHA1
            ELSE
               TERM1 = 0.0D+00
               TERM11 = 0.0D+00
            ENDIF
 
            IF ( ALPHA2.GT.0.0 ) THEN
               TERM2 = EXP(-ABS(E-EE)/ALPHA2)
               TERM22 = TERM2/ALPHA2
            ELSE
               TERM2 = 0.0D+00
               TERM22 = 0.0D+00
            ENDIF
 
            DENE = fxnexp(E, Mol, Dens)
            PROB = DENE*((1.0D00-DC(Mol,4))*TERM1+DC(Mol,4)*TERM2)
            Hs = DENE*((1.0D00-DC(Mol,4))*TERM11+DC(Mol,4)*TERM22)
            Hs = PROB/Hs
            IF ( Hs.GT.MAX(ALPHA1,ALPHA2) ) THEN
               Hs = MAX(ALPHA1,ALPHA2)
c            ELSEIF ( Hs.LT.MIN(ALPHA1,ALPHA2) ) THEN
c               Hs = MIN(ALPHA1,ALPHA2)
c               Hs = MAX(Hs,Egrain1)
            ENDIF
  
                        
         ELSEIF ( ITYPE(Mol).EQ.3 ) THEN
           tfac = (TEMPI/300.0)**DC(Mol,8)
           ALPHA1= ( ABS(DC(Mol,1) + E*( DC(Mol,2)+E*DC(Mol,3) ) )*tfac)
 
            IF ( ALPHA1.GT.0.0 ) THEN
               TERM1 = EXP(-(0.5*(E-EE-DC(Mol,4))/ALPHA1)**2)
               TERM11 = -(E-EE-DC(Mol,4))*TERM1/(ALPHA1**2)
            ELSE
               TERM1 = 0.0D+00
               TERM11 = 0.0D+00
            ENDIF
 
            PROB = TERM1
            Hs = ABS(TERM1/TERM11)
            IF ( Hs.GT.0.5*ALPHA1 ) Hs = 0.5*ALPHA1
 
            
         ELSEIF ( ITYPE(Mol).EQ.4 ) THEN
            ALPHA1 = ABS(DC(Mol,1)+(E*(DC(Mol,2)+E*DC(Mol,3))))
            ALPHA2 = ABS(DC(Mol,5)+(E*(DC(Mol,6)+E*DC(Mol,7))))
 
            IF ( ALPHA1.GT.0.0 ) THEN
               TERM1 = EXP(-ABS(E-EE)/ALPHA1)
               TERM11 = TERM1/ALPHA1
            ELSE
               TERM1 = 0.0D+00
               TERM11 = 0.0D+00
            ENDIF
 
            IF ( ALPHA2.GT.0.0 ) THEN
               TERM2 = EXP(-ABS(E-EE)/ALPHA2)
               TERM22 = TERM2/ALPHA2
            ELSE
               TERM2 = 0.0D+00
               TERM22 = 0.0D+00
            ENDIF
 
            Fractt = DC(Mol,4) + DC(Mol,8)*E
 
            PROB = (1.0D00-Fractt)*TERM1 + Fractt*TERM2
            Hs = PROB/((1.0D00-Fractt)*TERM11+Fractt*TERM22)
            IF ( Hs.GT.MAX(ALPHA1,ALPHA2) ) THEN
               Hs = MAX(ALPHA1,ALPHA2)
c            ELSEIF ( Hs.LT.MIN(ALPHA1,ALPHA2) ) THEN
c               Hs = MIN(ALPHA1,ALPHA2)
c               Hs = MAX(Hs,Egrain1)
            ENDIF
            
         ELSEIF ( ITYPE(Mol).EQ.5 ) THEN
            ALPHA = ABS(DC(Mol,1)+(E*(DC(Mol,2)+E*DC(Mol,3))))
            EN = ABS(DC(Mol,5)+(E*(DC(Mol,6)+E*DC(Mol,7))))
            ARG = (ABS(E-EE)/ALPHA)**EN
 
            IF ( ARG.GT.100 ) THEN
               PROB = 0.
               Hs = 10000.      ! Take giant steps, because PROB=0, anyway
            ELSE
               PROB = EXP(-ARG)
               IF ( ALPHA.GT.0.0 ) THEN
                  Hs = ALPHA/EN
               ELSE
                  Hs = Egrain1
               ENDIF
            ENDIF
            
         ELSEIF ( ITYPE(Mol).EQ.6 ) THEN
            ALPHA1 = ABS(DC(Mol,1)+(E*(DC(Mol,2)+E*DC(Mol,3))))
            EN = ABS(DC(Mol,4)+E*(DC(Mol,5)))
            Fractt = DC(Mol,6)
            ALPHA2 = ABS(DC(Mol,7)+E*(DC(Mol,8)))
 
            IF ( ALPHA1.GT.0.0 ) THEN
               ARG = (ABS(E-EE)/ALPHA1)**EN
               IF ( ARG.GT.100. ) ARG = 100.
               TERM1 = EXP(-ARG)
               TERM11 = EN/ALPHA1
            ELSE
               TERM11 = 1.00D+00/Egrain1
            ENDIF
 
            IF ( ALPHA2.GT.0.0 ) THEN
               TERM2 = EXP(-ABS(E-EE)/ALPHA2)
               TERM22 = TERM2/ALPHA2
            ELSE
               TERM2 = 0.0D+00
               TERM22 = 0.0D+00
            ENDIF
 
            A = 1.0/TERM11
            B = 1.0/TERM22
 
            PROB = (1.0D00-Fractt)*TERM1 + Fractt*TERM2
            Hs = PROB/((1.0D00-Fractt)*TERM11+Fractt*TERM22)
            IF ( Hs.GT.MAX(A,B) ) THEN
               Hs = MAX(A,B)
c            ELSEIF ( Hs.LT.MIN(A,B) ) THEN
c               Hs = MIN(A,B)
c               Hs = MAX(Hs,Egrain1)
            ENDIF
 
            
         ELSEIF ( ITYPE(Mol).EQ.7 ) THEN
            ALPHA = ABS(DC(Mol,1)+(E*(DC(Mol,2)+E*DC(Mol,3))))  ! Energy units
            BET = DC(Mol,4) + E*DC(Mol,5)
            TT = ALPHA*BET
            IF ( TT.LE.0.0 ) THEN
               WRITE(*,*) 
     &          '*** Weibull Model requires ALPHA & BETA > 0 ***'
               WRITE(KOUT,*)  
     &          '*** Weibull Model requires ALPHA & BETA > 0 ***'
               STOP
            ENDIF
 
            Rob_3 = 1.0
            Hs = MAX(ALPHA/BET,Rob_3)   ! Energy units
            ALPHA = ALPHA**BET          ! units for Weibull Distribution
            X = ABS(E-EE)
 
            IF ( X.GT.0.0 ) THEN        ! if X >0
               ARG2 = (X**(BET-2.))/ALPHA
               ARG1 = X*ARG2
               ARG = X*ARG1
               IF ( ARG.GT.100. ) THEN  ! then ARG1 and ARG2 must be even bigger
                  PROB = 0.0
                  Hs = 10000.   ! Take giant steps, because PROB = 0 anyway
                ELSE
                  PROB = (BET/ALPHA)*ARG1*EXP(-ARG)      ! This is General Weibull Model
                ENDIF
            ELSEIF ( ABS(BET-1.0).LT.1.0E-06 ) THEN     ! if X=0 and BET=1, exponential model
               PROB = 1.0
             ELSE                        ! if X=0 and BET≠1
               PROB = 0.0
            ENDIF

         ELSEIF ( ITYPE(Mol).EQ.8 ) THEN
            ALPHA = ABS(DC(Mol,1)+(E*(DC(Mol,2)+E*DC(Mol,3))))
            WIDTH = ABS(DC(Mol,5)+(E*(DC(Mol,6)+E*DC(Mol,7))))
            X = (E-EE-ALPHA)
            PROB = WIDTH*WIDTH/(X*X+WIDTH*WIDTH)
            Hs = MIN(WIDTH/PROB,0.1*ALPHA)
            Hs = MAX(Hs,Egrain1)
 
            
         ELSEIF ( ITYPE(Mol).EQ.9 ) THEN
            tfac = (TEMPI/300.0)**DC(Mol,8)
            ALPHA= ( ABS(DC(Mol,1) + E*( DC(Mol,2)+E*DC(Mol,3) ) )*tfac)
            F = (DC(Mol,5) + (E/DC(Mol,6))**DC(Mol,7))*tfac
 
            IF ( ALPHA.GT.0.0 ) THEN
               TERM1 = EXP(-ABS(E-EE)/ALPHA)
            ELSE
               TERM1 = 0.0D+00
            ENDIF
 
            IF ( ABS(EE-E).LT.1.0E-06 ) THEN                    ! Elastic collisions
               PROB = F*TERM1/(DC(Mol,4)+F) + DC(Mol,4)/(DC(Mol,4)+F)
            ELSE
               PROB = F*TERM1/(DC(Mol,4)+F)
            ENDIF
 
            Hs = MAX(ALPHA,Egrain1)
 
            
         ELSEIF ( ITYPE(Mol).EQ.10 ) THEN
            ALPHA = ABS( DC(Mol,1) + E*DC(Mol,2) )
            BETA = ABS( DC(Mol,3) )
 
            IF ( ALPHA.GT.0.0 ) THEN
               PROB = EXP( -( ABS(E-EE)/ALPHA )**BETA )
               Hs = ALPHA
            ELSE
               PROB = 0.0D+00
               Hs = 1.0
            ENDIF
 
            
         ELSEIF ( ITYPE(Mol).EQ.11 ) THEN
             ALPHA = ABS(DC(Mol,1)*(1-EXP(-((E/DC(Mol,2))**DC(Mol,3))))
     &        +DC(Mol,4))
 
            IF ( ALPHA.GT.0.0 ) THEN 
              TERM1 = EXP(-ABS(E-EE)/ALPHA)
              PROB = TERM1
              Hs = ALPHA
             ELSE
              PROB = 0.0D+00
              Hs = 1.0
            ENDIF
 
            
         ELSEIF ( ITYPE(Mol).EQ.12 ) THEN
            ALPHA1 = ABS( DC(Mol,1) + E*DC(Mol,2) + 
     &                DC(Mol,3)*EXP(-E/DC(Mol,4)) )

            IF ( ALPHA1.GT.0.0 ) THEN
               TERM1 = EXP(-ABS(E-EE)/ALPHA1)
               TERM11 = TERM1/ALPHA1
            ELSE
               TERM1 = 0.0D+00
               TERM11 = 0.0D+00
            ENDIF
 
            PROB = TERM1

            Hs = PROB/TERM11
            IF ( Hs.GT.ALPHA1 )  Hs = ALPHA1
c            IF (E .GT. Emax1+5.*Hs ) Hs = 5.*Hs    ! Big steps for exponential model at high E
 
            
         ELSEIF ( ITYPE(Mol).EQ.13 ) THEN
            ALPHA1 = ABS(DC(Mol,1) + E*DC(Mol,2))  ! Energy units
            ALPHA2 = ABS(DC(Mol,3) + E*DC(Mol,4))  ! Energy units
            BET = ALPHA1 + 0.5d+00*(ALPHA2 - ALPHA1)*(1.d+00 - 
     &             TANH((DC(Mol,5)-E)/DC(Mol,6)))                    ! Switching function
            TERM1 = EXP(-ABS(E-EE)/BET)
            TERM11 = TERM1/BET
            PROB = TERM1
            Hs = PROB/TERM11
            IF ( Hs.GT.BET )  Hs = BET
c            IF (E .GT. Emax1+5.*Hs ) Hs = 5.*Hs    ! Big steps for exponential model at high E

            
         ELSEIF ( ITYPE(Mol).EQ.14 ) THEN
            ALP = ABS( DC(Mol,1) )
            TERM1 = EXP(-(E-EE)/ALP )
            TERM11 = TERM1/ALP
            PROB = TERM1
            DENE = fxnexp(E, Mol, Dens) 
            DENEE = fxnexp(EE, Mol, Dens) 
            
            Statest = 0.5d+00/Egrain1
            IF ( DENE .GT. Statest .AND. DENEE .GT. Statest) THEN
              PROB = SQRT( DENEE*EXP(-(E-EE)*1.4388/TEMPI)/DENE )*TERM1
            ELSE
              PROB = 1.e-20
            ENDIF
            
            Hs = ALP

         ENDIF !  ITYPE

      ENDIF ! want to return titles or probs
 

      RETURN
  
      END SUBROUTINE PSTEP

!  ---------------------------------------------------------------------
      REAL(8) FUNCTION Estart(Iflag,Icalc,IR,Molinit,Einit,Tvib)

c            Selects initial energy for trial
c      Iflag      = 0      Initialize Pstart
c                 = 1      select Estart from Pstart
c      Icalc    : type of initial distribution
c                 = 1      Monoenergetic at energy Einit
c                 = 2      Thermal with energy offset Einit
c                 = 3      Chemical activation from 'product' #IR
c                 = 4      Read Pstart from file 'multiwell.pstart'
c      IR       : index of 'product' from which chemical activation originates
c      Molinit  : index of initial molecule (well)
c      Einit    : initial energy (relative to ZPE of Molinit; units: cm-1)
c      Tvib     : Temperature (vibrational)
c      Edels    : lowest energy (cm-1) of distribution (e.g. energy offset)
c      Hstart   : energy step-size (cm-1)
 
      USE declare_mod
      USE utili_mod
 
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: Iflag
      INTEGER, INTENT(IN) :: Icalc
      INTEGER, INTENT(IN) :: IR
      INTEGER, INTENT(IN) :: Molinit
      REAL(8), INTENT(IN) :: Tvib
      REAL(8), INTENT(IN) :: Einit
      INTEGER :: i , ichan , IE , IRUN , itest , Igg
      REAL(8) :: AA , RX
      REAL(8) :: E
      REAL(8) :: X , Y , Y1 , Y2 , Y3 , Term , Tlast ,
     &          Tiggp , Tiggm
     
      Statest = 0.1/Egrain1

      IF ( Iflag.EQ.0 ) THEN                    ! Initialize Pstart
         IF ( Icalc.EQ.1 ) THEN                             ! DELTA FUNCTION: No need for initialization
            RETURN
         ELSEIF ( Icalc.EQ.2 ) THEN                     ! THERMAL DISTRIBUTION
            Edels = 0.0d+00            ! Offset
         ELSEIF ( Icalc.EQ.3 ) THEN                     ! CHEMICAL ACTIVATION
            itest = 0
            DO 40 i = 1 , Nchan(Molinit)                ! Find reaction index
               IF ( Jto(Molinit,i).EQ.IR ) THEN
                  ichan = i
                  itest = 1
               ENDIF
 40         CONTINUE
            IF ( itest.EQ.0 ) THEN                  ! fatal error
               WRITE (6,99001) Molinit , IR
               STOP
            ENDIF
            if(itun( Molinit , ichan ) .eq. 0 ) then    ! NO TUNNELING
               Edels = Eor(Molinit,ichan)               !  (incorporates centrifugal correction)
            else
               Edels = Eor(Molinit,ichan) - Emax1      ! with tunneling, start at energy Emax1 below the barrier top
               IF ( Edels .LT. 0.0 ) Edels = 0.0
            endif

         ELSEIF ( Icalc.EQ.4 ) THEN                         ! READ FROM FILE (NORMALIZED CUMULATIVE DISTRIBUTION)
            OPEN (12,FILE='DensData/'//'multiwell.pstart',STATUS='OLD')
            READ (12,*) Jsize , Hstart , Edels     ! Edels and Hstart entered in cm-1
            IF ( Jsize .GT. Imax ) THEN
               write(*,*) '*** FATAL: Pstart array is too long ***'
               write(*,*) ' '
               STOP
            ENDIF

            DO 20 i = 1 , Jsize
               READ (12,*) AA , Pstart(i)       ! AA is dummy real variable
 20         CONTINUE
            CLOSE (12)
            RETURN
         ENDIF
 
         Hstart = Egrain1                                            ! ********* Set step size ************

         Tlast = 0.0D+00           ! previous term for trapezoidal rule
         DO 50 i = 1 , Jsize
            IF ( Icalc.EQ.2 ) THEN                         ! THERMAL, WITH OR WITHOUT OFFSET
               E = (i-1)*Hstart
               Term = fxnexp( E, Molinit, Dens )*exp(-hor*E/Tvib)
            ELSEIF ( Icalc.EQ.3 ) THEN                     ! CHEMICAL ACTIVATION, WITH OR WITHOUT OFFSET
               if(itun( Molinit , ichan ) .eq. 1 ) then    ! with tunneling
                 E = Edels + (i-1)*Hstart                       ! start at Emax1 below the barrier 
                 CALL NTERP2((E-Edels),IE,X)
                 Y1 = TRate(Molinit,ichan,IE-1)
                 Y2 = TRate(Molinit,ichan,IE)
                 Y3 = TRate(Molinit,ichan,IE+1)
                 Y  = EXP( XLINT(Y1,Y2,Y3,X) )                  ! with tunneling: Rate constant k(E) [stored as log]
                 Term = Y*fxnexp( E, Molinit, Dens )*exp(-hor*E/Tvib)
               else
                 E = Edels + (i-1)*Hstart                       ! Chemical activation WITHOUT TUNNELING
                 Y = rxnofe( E-Edels, Molinit , ichan , Rate )           ! Interpolated k(E)
                 Term = Y*fxnexp( E, Molinit, Dens )*exp(-hor*E/Tvib)    ! no tunneling
               end if
            ENDIF
            IF ( Hstart .EQ. Egrain1 ) THEN
               Pstart(i) = Hstart*Term                     ! Constant within each grain at low energies
            ELSE
               Pstart(i) = 0.5D+00*Hstart*(Term+Tlast)     ! Trapezoidal rule for larger steps
            ENDIF
            Tlast = Term           ! previous term for trapezoidal rule
 50      CONTINUE
 
         DO 100 i = 2 , Jsize
            Pstart(i) = Pstart(i) + Pstart(i-1)          ! Cumulative
 100     CONTINUE
 
         IF ( Pstart(Jsize) .LE. 0.0 ) THEN
           write(*,*) '**** FATAL: initial E distrib undefined ****'
           write(*,*) '    ** Probably because T is too low **'
           write(*,*) '  '
           write(KOUT,*) '**** FATAL: initial E distrib undefined ****'
           write(KOUT,*) '    ** Probably because T is too low **'
           write(*,*) '  '
           STOP
         ENDIF
 
         DO 150 i = 1 , Jsize
            Pstart(i) = Pstart(i)/Pstart(Jsize)          ! Normalize Cumulative distribution
 150     CONTINUE
 
c  ................................................................................

      ELSEIF ( Iflag.EQ.1 ) THEN                ! Select initial energy from Pstart
 
         IF ( Icalc.EQ.1 ) THEN
            Estart = Einit                              ! DELTA FUNCTION
            IF (Estart .LE. Emax1) THEN                 ! Align Estart with grains
                Igg = INT( Estart/Egrain1 + 0.5d+00 )
                Estart = Egrain1*Igg
                IF (Estart .LT. 0.0d+00) Estart = 0.0d+00
            ENDIF
            RETURN
         ELSE                                           ! THERMAL, CHEMACT, EXTERNAL
 
            RX = RAN1(IDUM)
            i = 0
            IRUN = 1
            DO WHILE ( IRUN .EQ. 1 )
               i = i + 1
               IF ( Pstart(i) .GT. RX ) IRUN = 0
            ENDDO
 
            Nstart(i) = Nstart(i) + 1           ! Record selected distribution for posterity

            IF ( i .GT. 1 ) THEN
               IF ( Hstart .EQ. Egrain1 ) THEN                          ! for Egrain1 steps
                  Estart = Edels + Einit + Hstart*(i-1)
                ELSE
                  Estart = Edels + Einit + Hstart*( (i-1) +             ! interpolate for other steps
     &               (RX - Pstart(i-1))/(Pstart(i) - Pstart(i-1)) )

                  Igg = NINT( Estart/Egrain1 )                    ! align with energy grains at low energy
                  Estart = Egrain1*Igg
                  IF ( fxnexp(E,Molinit,Dens ) .LT. Statest ) THEN  ! if a state is not present, move up or down one grain
                    Tiggp = fxnexp( Estart+Egrain1, Molinit, Dens)
                    Tiggm = fxnexp( Estart-Egrain1, Molinit, Dens)
                    IF ( Tiggp .GT. Tiggm ) THEN
                      Estart = Estart + Egrain1
                    ELSE
                      Estart = Estart - Egrain1
                    ENDIF  ! Tiggp
                  ENDIF  ! fxnofe
               ENDIF  ! Hstart
             ELSE
                 Estart = Edels + Einit
            ENDIF  ! i

         ENDIF  ! Icalc
      ENDIF     ! Iflag

      IF ( Estart.LT.0.0 ) Estart = 0.0d+00

c      write(*,*) Estart

      RETURN
 
99001 FORMAT (/'Data file error: CHEMACT rxn does not exist:',I3,' to',
     &        I3/'****** Hit RETURN to terminate run *******')
      END FUNCTION Estart

!  ---------------------------------------------------------------------
      REAL(8) FUNCTION Etherm(Mol,T)
c
c      Average thermal energy of Mol at temperature T
c      Uses numerical integration of densities of states.

      USE declare_mod
      USE utili_mod
      
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: Mol
      REAL(8), INTENT(IN) :: T
      REAL(8) :: Q , E , ET , Term , Test , H , 
     &                 Tlast , Elast
c
         ET = 0.0D+00
         Q = 0.0D+00
         Elast = 0.0D+00
         Tlast = 0.0D+00
         H = Egrain1
         E = -H
         Test = 1.0
         DO WHILE ( Test.GT.1.0E-08 )
            E = E + H
            Term = fxnexp(E, Mol, Dens)*exp(-hor*E/T)
            Q = Q + 0.5*H*(Tlast+Term)                  ! Trapezoidal rule
            ET = ET + 0.5*H*(Elast*Tlast+E*Term)        ! Trapezoidal rule
            Tlast = Term
            Elast = E
            IF ( E.GT.10.*T ) Test = E*Term/ET
         ENDDO
 
         Etherm = ET/Q
         RETURN
      END FUNCTION Etherm
      

!  ---------------------------------------------------------------------
      SUBROUTINE QKinf(Mol,nchann,tx,Kinf,Ainf,Einf,kcol,kosc)
c
c      Calculates high pressure limiting unimolecular rate constant Kinf
c         (See user Manual for details)
c
c      Mol      = index of molecule
c      nchann      = index of rxn channel
c      tx      = translational temperature
c      EC   = rxn "critical energy"; true Eor for classical, or lower energy bound for QM tunneling
c      Ainf      = A-factor for kinf
c      Einf      = Activation energy for kinf
c      rxnofe  = interpolated value of k(E)

      USE declare_mod
      USE utili_mod

      IMPLICIT NONE
      INTEGER, INTENT(IN) :: Mol
      INTEGER, INTENT(IN) :: nchann
      REAL(8), INTENT(IN) :: Tx
      REAL(8), INTENT(OUT) :: Kinf
      REAL(8), INTENT(OUT) :: Ainf
      REAL(8), INTENT(OUT) :: Einf
      REAL(8), INTENT(IN)  :: kcol
      REAL(8), INTENT(OUT) :: kosc

      INTEGER :: i , mm
      REAL(8) :: Q , RK , TT , E , D 
      REAL(8),DIMENSION(3) :: K
      REAL(8),DIMENSION(3) :: TS
      REAL(8) :: Et , arg , 
     &                 Estep, RR, QQ , Qlast , Rlast , TTSTEP , 
     &                 EC , rat , rrke
      REAL(8) :: Qts, QQts, Qtslast

      TTSTEP = 0.01d+00/tx                       ! Step size in 1/Temp; tx = input temperature
      TT = 1.0d+00/tx-2.d+00*TTSTEP               ! 1/Temp (starting value)
       
      DO mm = 1,3
        K(mm) = 0.0d+00
        RK = 0.0d+00
        Q  = 0.0d+00
        Qts = 0.0d+00
        Qlast = 0.0d+00
        Qtslast = 0.0d+00
        Rlast = 0.0d+00

        TT = TT + TTSTEP                           ! 1/T
        TS(mm) = TT                                ! 1/T

        IF ( NCENT(Mol,nchann) .NE. -2 ) THEN
          EC = Eor(Mol,nchann)                     ! NORMAL METHODS include centrifugal correction ( Eor = Eo for NOCENT )
        ELSE
          EC = Eo(Mol,nchann)                      ! LEGACY METHOD DOES NOT include centrifugal correction
        ENDIF
      
        Estep = Egrain1
        E = 0.0d+00                                ! active energy in thermal reactant, or TS
        arg = EC                                   ! active energy in excited reactant
        DO i = 1 , imax1                           ! summation continues until end of first part of double array
          QQ = fxnexp( E, Mol, Dens )*exp(-E*hor*TT)   ! reactant partition fxn; E=active energy in thermal reactant
          Q  = Q  + 0.5*Estep*(QQ + Qlast)             ! reactant partition fxn; Trapezoidal Rule
          Qlast = QQ
          D = fxnexp( arg, Mol, Dens )                ! density of states of EXCITED reactant
          QQts = D*exp(-arg*hor*TT)                    ! EXCITED reactant partition fxn for STRONG COLLIDER calculation
          RR = D*rxnofe( E,Mol,nchann,Rate )*exp(-arg*hor*TT)  ! thermal rate constant; arg=active energy in EXCITED reactant; rxnofe = interpolated k(E)
          RK = RK + 0.5*Estep*(RR + Rlast)             !  +++++++ rate constant k(inf) ++++++++++++
          Qts = Qts + 0.5*Estep*(QQts+Qtslast)         ! partition fxn for STRONG COLLIDER calculation
          Rlast = RR
          Qtslast = QQts
          E = E + Estep                                ! E = active energy in thermal reactant
          arg = EC + E                                 ! arg = active energy in TS and EXCITED reactant ( >EC ) 
        END DO  ! i
c
c       NOW CONTINUE TO HIGHER ENERGIES
c 
        Estep = 0.05*tx/hor                           ! New step size (cm-1)
        IF ( Estep. LT. Egrain1 ) Estep = Egrain1
        rrke = 0.0
        DO WHILE ( arg.LT.(Emax2-Estep) .AND. (rrke .GT. 100.) )     ! summation continues for second part of double array
          QQ = fxnexp( E, Mol, Dens )*exp(-E*hor*TT)                 ! reactant partition fxn; E=active energy in thermal reactant
          Q  = Q  + 0.5*Estep*(QQ + Qlast)                           ! reactant partition fxn; Trapezoidal Rule
          Qlast = QQ
          D = fxnexp( arg, Mol, Dens )                               ! density of states of EXCITED reactant
          QQts = D*exp(-arg*hor*TT)                                  ! EXCITED reactant partition fxn for STRONG COLLIDER calculation
          rrke = rxnofe( E,Mol,nchann,Rate )                         ! rxnofe = interpolated k(E)
          RR = D*rrke*exp(-arg*hor*TT)                               ! thermal rate constant; arg=active energy in EXCITED reactant
          RK = RK + 0.5*Estep*(RR + Rlast)                           !  +++++++ rate constant k(inf) ++++++++++++
          Qts = Qts + 0.5*Estep*(QQts+Qtslast)                       ! partition fxn for STRONG COLLIDER calculation
          Rlast = RR
          Qtslast = QQts
          E = E + Estep
          arg = EC + E                            ! arg=active energy in TS ( >EC )
        END DO  ! i

       K(mm) = RK/Q                               ! Rate constant, including any centrifugal corrections
                    
       IF ( mm .EQ. 2) THEN
         kosc = kcol*Qts/Q                        ! Strong-collider low pressure limit rate constant
       ENDIF

      END DO  ! m=1,3 for T-dependence
c
c     Now find Arrhenius parameters
c
      Kinf = K(2)
      IF ( K(1)/K(3) .GT. 0.d+00 ) THEN 
        Et = -LOG( K(1)/K(3) )/( TS(1) - TS(3) )              ! Et = Ea/R
        Einf = (1.9872065d-03)*Et                             ! kcal/mole
        Ainf = K(2)*exp( Et*TS(2) )                           ! 1/s
      ELSE
        Einf = 0.0
        Ainf = 0.0
      ENDIF
      
      IF ( NCENT(Mol,nchann) .EQ. -2 ) THEN      ! Legacy centrifugal correction (not recommended)
        IF (    (MolMom(Mol) .GT. 0.0 )            
     &   .AND. (TSmom(Mol,nchann) .GT. 0.0)  )  THEN
         rat = TSmom(Mol,nchann)/MolMom(Mol)
        ELSE
         rat = 1.0d+00
        ENDIF
        Kinf = rat*Kinf
        Ainf = rat*Ainf
        kosc = rat*kosc
      ENDIF
      
      RETURN
      END SUBROUTINE QKinf

!  ---------------------------------------------------------------------
      SUBROUTINE QKinftun(Mol,nchann,tx,Kinf,Ainf,Einf)
      USE declare_mod
      USE utili_mod

      IMPLICIT NONE
      INTEGER, INTENT(IN) :: Mol
      INTEGER, INTENT(IN) :: nchann
      REAL(8), INTENT(IN) :: tx
      REAL(8), INTENT(OUT) :: Kinf
      REAL(8), INTENT(OUT) :: Ainf
      REAL(8), INTENT(OUT) :: Einf
      INTEGER :: i , IE , Itop , mm, i2
      REAL(8) :: Q , RK , TT , E , D ,  
     &           Et , arg , Y1 , Y2 , Y3 , Y , 
     &           X , Estep, RR, QQ , Qlast , Rlast , TTSTEP , 
     &           EC, E2, Qcheck, Qts, QQts, Qtslast
      REAL(8),DIMENSION(3) :: K
      REAL(8),DIMENSION(3) :: TS
c
c      Calculates high pressure limiting unimolecular rate constant according to
c         Equation (4.30) of Robinson & Holbrook, "Unimolecular Reactions"
c         (Wiley-Interscience, London, 1972), p. 90.
c
c      Mol      = index of molecule
c      nchann      = index of rxn channel
c      t      = translational temperature
c      Kinf      = infinite pressure rate constant [Robinson & Holbrook Eq. 4.30]
c      Ainf      = A-factor for kinf
c      Einf      = Activation energy for kinf
c
c
c     Centrifugal correction at temperature = tx
c
      EC =  Eor(Mol,nchann)                 ! Includes centrifugal correction

      TTSTEP = 0.001d+00/tx                  ! Step size in 1/Temp
      TT = 1.0d+00/tx-2.d+00*TTSTEP          ! 1/Temp

      DO 1000 mm = 1,3
        K(mm) = 0.0d+00
        RK = 0.0d+00
        Q  = 0.0d+00
        RR = 0.0d+00
        QQ = 0.0d+00
        Qlast = 0.0d+00
        Rlast = 0.0d+00
        QQts=0.0d00
        Qtslast=0.0d0
        Qts=0.0d0

        TT = TT + TTSTEP                     ! 1/T
        TS(mm) = TT                          ! 1/T

        Estep = Egrain1
        E = -Estep*(imax1+1)        

      DO 100 i = 1 , 2*imax1                       ! start at imax1*Egrain1 below the reaction threshold
         E = E + Estep                             ! summation continues till end of first part of double array

        if(i.gt.imax1)then                         ! allows Q(Trans) to be calculated from 0-->inf
          QQ=exp( Dens(Mol,i-imax1)-(E*hor*TT) )   ! molecule partition function;correct E and Dens
                                                   ! both start at 0 Energy and index=0 i.e. (i-imax1)=0
          Q  = Q  + 0.5*Estep*(QQ + Qlast)         ! Trapezoidal Rule
          Qlast = QQ
          QQts=D*exp(-E*hor*TT)
          Qts=Qts+0.5*Estep*(QQts+Qtslast)
          Qtslast=QQts
        else
          continue
        end if

         arg = E + EC                              ! arg is E + barrier height; starts at -Egrain1*imax1
                                                   ! below the barrier

        if(arg.ge.0)then
          
          D = fxnexp( arg, Mol,Dens )               ! density of states of molecule evaluated at
                                                    ! the (top of the barrier - Egrain1*imax1)

        RR = D*exp(TRate(Mol,nchann,i)-(E*hor*TT))! this corresponds to Eq. 12 in MW pub. and is
                                                   ! the integral from (Eo-(Egrain1*imax1)) to inf
          RK = RK + 0.5*Estep*(RR + Rlast)
          Rlast = RR
        end if

 100  CONTINUE
c
c          NOW CONTINUE TO HIGHER ENERGIES
c 
      Estep = 10.*Egrain1                          ! Step size (cm-1)
      Itop = NINT( ( 25.0/(hor*TT))/Estep  )       ! Integral to be extended to 25 kT  
      E2=E
      DO i = 1 , Itop                              ! Continuing the integral to higher energies
         E = E + Estep                             ! last energy from above loop + new Estep is 1st energy evaluated
         D = fxnexp( E, Mol,Dens )
         QQ = D*exp(-E*hor*TT)
         Q  = Q  + 0.5*Estep*(QQ + Qlast)           ! Trapezoidal Rule
         Qlast = QQ
         Qcheck = dabs(QQ/Q)
         if(Qcheck.lt.Qthresh)then
           goto 17
         else
           continue
         end if

         if(E.gt.Emax2)then
           write(*,*) '*** Emax2 exceeded; extend Emax2 in Densum'
           write(KOUT,*) '*** Emax2 exceeded; extend Emax2 in Densum'
           STOP
          else
            continue
          end if
        end do

 17     continue

        arg=E2+EC
        do i2=1,Itop
           arg=arg+Estep
           E2=E2+Estep
           D = fxnexp(arg, Mol, Dens) 

           CALL NTERP2(E2,IE,X)                    ! NTERP2 determines IE and X for a k(E)
                                                   ! array with Isize+imax1 elements; of which
                                                   ! imax1 elements are below the reaction threshold
           Y1 = TRate(Mol,nchann,IE-1)
           Y2 = TRate(Mol,nchann,IE)
           Y3 = TRate(Mol,nchann,IE+1)
           Y = XLINT(Y1,Y2,Y3,X)                  ! Y = interpolated ln k(E)
           RR = D*exp(Y-E2*hor*TT)               ! k(E) stored as log
           RK = RK + 0.5*Estep*(RR + Rlast)
           Rlast = RR
           QQts=D*exp(-E2*hor*TT)              ! partition function evaluated to determine
           Qts=Qts+0.5*Estep*(QQts+Qtslast)    ! when the calculation of RK can be stopped
           Qtslast=QQts
           Qcheck=dabs(QQts/Qts)
           if(Qcheck.lt.Qthresh)then
             goto 18
           else
             continue
           end if

           if(E2.gt.Emax2)then
           write(*,*) '*** Emax2 exceeded; extend Emax2 in Densum'
           write(KOUT,*) '*** Emax2 exceeded; extend Emax2 in Densum'
           STOP
           else
             continue
           end if
         end do

 18      continue    
       RK = RK*exp(-Eo(Mol,nchann)*hor*TT)/Q
       IF ( NCENT(Mol,nchann).EQ.1 ) THEN               ! If Centrifugal Correction
          K(mm) = RK*(TSmom(Mol,nchann)/MolMom(Mol))    ! kinf with adiabatic Qrot ratio
       ELSE
          K(mm) = RK                                    ! Adiabatic Qrot ratio = 1.0
       ENDIF
1000  CONTINUE

      Kinf = K(2)
      Et = -(LOG(K(3)/K(1)))/(TS(3)-TS(1))                  ! Ea/R
      Ainf = K(2)*exp(Et*TS(2))                             ! 1/s
      Einf = (1.9872065d-03)*Et                             ! kcal/mole

      RETURN
      END SUBROUTINE QKinftun
      
!  ---------------------------------------------------------------------

        SUBROUTINE ECKART(Mol,i,EE,j,psum)
c
c This subroutine calculates the sum of states of the transition state
c modified to take into account the tunneling transmission coefficient.
c This modified sum of states is defined as psum below. psum is 
c calculated as the convolution of the density of states of the transition 
c state and the tunneling transmission coefficient. The program returns 
c psum to ratearray for calculation of k(E) corrected for tunneling

        USE declare_mod
        USE utili_mod
  
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: Mol        ! index number for TS
        INTEGER, INTENT(IN) :: i          ! index number for rxn channel
        REAL(8), INTENT(IN) :: EE         ! energy above the critical energy
        INTEGER, INTENT(IN) :: j          ! j is the energy bin from SUBROUTINE ratearray at energy EE
        REAL(8), INTENT(OUT) :: psum      ! sum of states, including tunneling

        INTEGER :: ij , ijmax
        REAL(8) :: pi , e2 , e3 , econstant , vh , pre
        REAL(8) :: E1 , E , rij, Egrain
        REAL(8) :: a1 , b1 , c1 , c2 , c3 , c4 
        REAL(8) :: tprob , tstep
        REAL(8) :: a5,b5,c5,d5
        REAL(8) :: term1,term2,term3,term4,Echeck
c        SAVE 
!      constants used in applying the tunneling correction

        pi=4.0D+00*DATAN(1.0D+00)
        e2=Eor(Mol,i)**(-0.5d0) 
        e3=V1(Mol,i)**(-0.5d0)
        econstant=(e2+e3)**(-1.0d0)
        vh=plank*avogadro*clight*jtocm/1000.0d0
        pre=(vimag(Mol,i)*vh)**(-1.0d0)

        c1=2.0d0*pi
        c2=(Eor(Mol,i)*V1(Mol,i))*pre**(2.0d0)
        c3=1.0d0/16.0d0
        c4=c1*(c2-c3)**(0.5d0)

        psum=0.0d0                 !zero out sum of states (ss) array

!      Egrain will be the array (lower or upper) dependent integration step size

        if(j.le.imax1*2)then       !these lines set the integration step size
          Egrain=Egrain1           !j is the index from ratearray.f which tells
          tstep=Egrain1            !us at which energy EE we are at; for j.le.imax1*2
        else                       !we use Egrain1 to perform the integration
          Egrain=Egrain2           
          tstep=Egrain2
        end if

        ijmax = 1 + NINT(EE/Egrain)             ! revised to NINT on 12/13/2021
        do ij = 1 , ijmax          !only calculate ss at grain energies
          rij = (ij-1)*Egrain

          E=rij                    !energy in the molecule's other degrees of freedom

          E1=EE-E-Eor(Mol,i)       !energy in the rxn coordinate

          if(EE.lt.0.0d0)goto 18   !EE is the energy from the zpe of the reactant it 
                                   !can never be negative

          if(Eor(Mol,i).ge.V1(Mol,i))then  !endothermic and barrier-less rxn condition
            Echeck=Eor(Mol,i)-V1(Mol,i)

            if((EE.ge.Echeck).and.(E1.ge.-V1(Mol,i)))then    !tunneling corrections are undefined if
                                                             ![EE < (Eor-V1)] for endothermic reactions
              continue
            else
              goto 18
            end if

          else                              !exothermic rxn condition

            if((EE.ge.0).and.(E1.ge.-Eor(Mol,i)))then
              continue         
            else
              goto 18
            end  if
                            
          end if
  
          a1=4.0d0*pi*pre*((E1+Eor(Mol,i))**(0.5d0))*econstant
          b1=4.0d0*pi*pre*((E1+V1(Mol,i))**(0.5d0))*econstant

!      approximate sinh(x) and cosh(x) by exp(x); this approximation
!      is valid since x is the value of a1, b1 and c4, and in this case
!      the values for these variables are large 

          a5=a1
          b5=b1
          c5=((a1+b1)/2.0d0)
          d5=c4

          term1=b5-d5
          term2=c5-a5
          term3=d5-a5
          term4=c5-d5

          tprob=exp(term1)/(exp(term2)*exp(term4)+exp(term3)) ! tunneling probability

          if(tprob.le.1.0d00)then  ! all tunneling probabilities should be 1 or less
            continue
          else
            tprob=1.0d00             ! this covers the case that the exponential terms
          end if                     ! in the above expression become too large; when this
                                     ! happens (and it sometimes happens for energies much
                                     ! greater then the reaction threshold and for small imaginary
                                     ! frequency) tprob=NAN when in fact it should be equal to unity.
                                     ! Setting tprob=1 takes care of this problem

        if(tprob.lt.(tunthresh))then ! if tprob is less than tunthresh the contribution to the sum
          return                     ! of states of the transition state is very small and we stop
        else                         ! peforming the convolution
          continue
        end if
          psum=psum+tstep*tprob*fxnexp( E, i, densTS )             ! psum = sum of states of TS for channel i is interpolated at energy E
 18    CONTINUE
       end do 

      RETURN
      END SUBROUTINE ECKART
      
!  ---------------------------------------------------------------------

      SUBROUTINE rotunits( VROTIN , X , VROTOUT )
      
      IMPLICIT NONE
      SAVE
      CHARACTER(len=4), INTENT(IN) :: VROTIN      ! declared rot units
      REAL(8), INTENT(INOUT)       :: X            ! numerical value      
      CHARACTER(len=4), INTENT(IN) :: VROTOUT      ! desired new rot units
     
c      CONVERT UNITS FROM VROTIN TO VROTOUT
c
c       VROTIN or VROTOUT
c               = 'AMUA' for moment of inertia in units of amu*Ang^2
c       or      = 'GMCM' for moment of inertia units of gram*cm^2
c       or      = 'CM-1' for rotational constant in units of cm^-1
c       or      = 'MHZ' for rotational constant in units of MHz
c       or      = 'GHZ' for rotational constant in units of GHz
c         
c         Convert VROTIN to AMUA
c
      IF ( VROTIN .EQ. 'AMUA' .OR. VROTIN .EQ. 'amua' ) THEN
      ELSEIF ( VROTIN .EQ. 'GMCM' .OR. VROTIN .EQ. 'gmcm' ) THEN
        X = X / 1.660538782d-040
      ELSEIF ( VROTIN .EQ. 'CM-1' .OR. VROTIN .EQ. 'cm-1' ) THEN
        X = 16.85763D+00 / X
      ELSEIF ( VROTIN .EQ. 'MHZ' .OR. VROTIN .EQ. 'MHz' 
     &      .OR. VROTIN .EQ. 'mhz' .OR. VROTIN .EQ. 'Mhz' ) THEN
        X = 5.05379D+005 / X
      ELSEIF ( VROTIN .EQ. 'GHZ' .OR. VROTIN .EQ. 'GHz' 
     &      .OR. VROTIN .EQ. 'ghz' .OR. VROTIN .EQ. 'Ghz' ) THEN
        X = 5.05379D+002 / X
      ELSE
        write (*,*) 'FATAL: units (VROTIN) for rotations not recongized'
        STOP
      ENDIF
c
c         Convert AMUA to VROTOUT and RETURN
c
      IF ( VROTOUT .EQ. 'AMUA' .OR. VROTOUT .EQ. 'amua' ) THEN
        RETURN
      ELSEIF ( VROTOUT .EQ. 'GMCM' .OR. VROTOUT .EQ. 'gmcm' ) THEN
        X = X * 1.660538782d-040
        RETURN
      ELSEIF ( VROTOUT .EQ. 'CM-1' .OR. VROTOUT .EQ. 'cm-1' ) THEN
        X = 16.85763D+00 / X
        RETURN
      ELSEIF ( VROTOUT .EQ. 'MHZ' .OR. VROTOUT .EQ. 'MHz' 
     &      .OR. VROTOUT .EQ. 'mhz' .OR. VROTOUT .EQ. 'Mhz' ) THEN
        X = 5.05379D+005 / X
        RETURN
      ELSEIF ( VROTOUT .EQ. 'GHZ' .OR. VROTOUT .EQ. 'GHz' 
     &      .OR. VROTOUT .EQ. 'ghz' .OR. VROTOUT .EQ. 'Ghz' ) THEN
        X = 5.05379D+002 / X
        RETURN
      ELSE
        write (*,*) 'FATAL: rotation units (VROTOUT) not recongized'
        STOP
      ENDIF

      END SUBROUTINE rotunits
      
!  ---------------------------------------------------------------------
      SUBROUTINE UnitTest ( DUM , Punits , Eunits , Rotatunits , KOUT )
c
c    Checks keywords ('DUM') for validity and then assigns units or STOPS
      
      IMPLICIT NONE
      SAVE

      INTEGER, INTENT(IN) :: KOUT
      CHARACTER(LEN=4), INTENT(IN), DIMENSION(3) :: DUM
      CHARACTER(LEN=3), INTENT(INOUT) :: Punits
      CHARACTER(LEN=4), INTENT(INOUT) :: Eunits
      CHARACTER(LEN=4), INTENT(INOUT) :: Rotatunits 

        IF (   ( DUM(1) .EQ. 'TORR' )
     &    .OR. ( DUM(1) .EQ. 'TOR' ) 
     &    .OR. ( DUM(1) .EQ. 'ATM' ) 
     &    .OR. ( DUM(1) .EQ. 'BAR' ) 
     &    .OR. ( DUM(1) .EQ. 'MCC' )  
     &                           ) THEN
           Punits = DUM(1)
         ELSE
           WRITE(KOUT,*) '**** FATAL: NO PRESSURE UNITS RECOGNIZED ****'
           WRITE(*,*)    '**** FATAL: NO PRESSURE UNITS RECOGNIZED ****'
           WRITE(*,*)    '          (in SUBROUTINE UnitTest)'
           STOP
        ENDIF

        IF (   ( DUM(2) .EQ. 'KJ' ) 
     &    .OR. ( DUM(2) .EQ. 'KJOU' )
     &    .OR. ( DUM(2) .EQ. 'CM-1' ) 
     &    .OR. ( DUM(2) .EQ. 'KCAL' )
     &                           ) THEN
           Eunits = DUM(2)
         ELSE
           WRITE(KOUT,*) '**** FATAL: NO ENERGY UNITS RECOGNIZED ****'
           WRITE(*,*)    '**** FATAL: NO ENERGY UNITS RECOGNIZED ****'
           WRITE(*,*)    '          (in SUBROUTINE UnitTest)'
           STOP
         ENDIF
      
        IF (   ( DUM(3) .EQ. 'AMUA' )
     &    .OR. ( DUM(3) .EQ. 'GMCM' )
     &    .OR. ( DUM(3) .EQ. 'CM-1' )
     &    .OR. ( DUM(3) .EQ. 'MHZ' )
     &    .OR. ( DUM(3) .EQ. 'GHZ' ) 
     &    .OR. ( DUM(3) .EQ. 'GHz' ) 
     &                           ) THEN
           Rotatunits = DUM(3)
         ELSE
           WRITE(KOUT,*) '**** FATAL: NO ROTATION UNITS RECOGNIZED ****'
           WRITE(*,*)    '**** FATAL: NO ROTATION UNITS RECOGNIZED ****'
           WRITE(*,*)    '          (in SUBROUTINE UnitTest)'
           STOP
        ENDIF
      
      RETURN
      END SUBROUTINE UnitTest

!  ---------------------------------------------------------------------
!  ---------------------------------------------------------------------

      END MODULE subs_mod