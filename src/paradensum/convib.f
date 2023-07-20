!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!    
!    LICENSE NOTICE
!
!    Copyright (C) 2019 Michele Ceotto, Chiara Aieta, Fabio Gabas, 
!                       Thanh Lam Nguyen, and John R. Barker
!
!    Contact:
!
!    Michele Ceotto  (email: michele.ceotto@unimi.it)
!    Dipartimento di Chimica
!    Università degli Studi di Milano
!    via Golgi 19, Milano 20133, ITALY
!
!    or:
!
!    John R. Barker   (email: jrbarker@umich.edu)
!    Department of Atmospheric, Oceanic, and Space Sciences
!    College of Engineering
!    University of Michigan
!    Ann Arbor, MI 48109
!
!    http://clasp-research.engin.umich.edu/multiwell
!    This program is free software; you can redistribute it and/or
!    modify it under the terms of the GNU General Public License (version 2)
!    as published by the Free Software Foundation.
!   
!    This program is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!    GNU General Public License for more details.
!   
!    See the 'ReadMe' file for a copy of the GNU General Public License,
!    or contact:
!   
!    Free Software Foundation, Inc.
!    59 Temple Place - Suite 330
!    Boston, MA 02111-1307, USA.
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!                                                                 
!                         PROGRAM paradensum    
!
!                               by      
!
!           by Michele Ceotto, Chiara Aieta, Fabio Gabas,
!               Thanh Lam Nguyen, and John R. Barker    
!                                                                  
!           ***PARallel Anharmonic DENsities and SUM of states***
!
!                             based on
!
!          The parallel implementation* of the Wang-Landau algorithms 
!          for densities of states**
!
!    Literature Citations:                                                    
!    *parallel implementation of anaharmonic density of states
!    C. Aieta, F. Gabas and M. Ceotto, J. Phys. Chem. A., 120(27), 
!    4853-4862 (2016).
!                                                                  
!    **Density of states algorithms
!    F. Wang and D. P. Landau, Phys. Rev. Letters 86, 2050-2053 (2001).           
!    M. Basire, P. Parneix, and F. Calvo, J. Chem.Phys. 129, 081101 (2008).    
!    T. L. Nguyen and J.R. Barker, J. Phys. Chem. A., 114, 3718–3730 (2010).
!    
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
c-----------------------------------------------------------------------------------------
c     From input as wa, w0, or wf frequencies, convert to the other forms:
c       wa(i): harmonic frequencies with zero of energy at the potential minimum
c          Evib = SUMi( wai*(vi+1/2) ) + SUMi( SUMj( Xij*(vi+1/2)*(vj+1/2) ) )
c       w0(i): frequencies with zero of energy at zpe
c          Evib = SUMi( w0i*vi ) + SUMi( SUMj( Xij*vi*vj ) )
c       wf(i) = fundamental frequencies (for 0-1 transitions with all vj = 0)
c
c     zpe = zero point energy, based on We and X matrix
c
c     Average Frequences:
c       ave = average wa(i)
c       av0 = average w0(i)
c       avf = average wf(i)
c
c     [Herzberg, Infrared and Raman Spectra (D. van Nostrand Co., 1945), p. 206ff]
c     Eq. numbers from Herzberg
c-----------------------------------------------------------------------------------------

      SUBROUTINE convib( ns , WW ) 
      
      USE decl_alloc

      IMPLICIT NONE
      INTEGER(4) ns , i , j
!      REAL(8) wa , w0 , wf , xa , ya , za , zpe, 
!      REAL(8)  zpe, 
!     &                 av0 , ave , avf , Viblo
      CHARACTER WW*2
!      COMMON/wxyz/ wa(100) , w0(100) , wf(100) , xa(100,100) ,
!      COMMON/wxyz/  
!     &   zpe , ave ,
!     &   av0 , avf , Viblo

!      SAVE
      
      IF ( WW .EQ. 'We' ) THEN

      DO i = 1 , ns                               ! start conversion
        w0(i) = wa(i) + xa(i,i)                   ! Eq. (II,273)
        wf(i) = wa(i) + 2.0d+00*xa(i,i)           ! fundamental
          DO j = 1, ns
            IF ( i .NE. j ) THEN
            w0(i) = w0(i) + 0.5d+00*xa(i,j)       ! Eq. (II,273)
            wf(i) = wf(i) + 0.5d+00*xa(i,j)       ! fundamental
            ENDIF
          END DO
      END DO                                      ! end conversion

      ELSEIF ( WW .EQ. 'W0' ) THEN

      DO i = 1 , ns                               ! start conversion
        wa(i) = w0(i) - xa(i,i)                   ! Eq. (II,273)
          DO j = 1, ns
            IF ( i .NE. j ) THEN
            wa(i) = wa(i) - 0.5d+00*xa(i,j)       ! Eq. (II,273)
            ENDIF
          END DO
      END DO                                      ! end conversion
      DO i = 1 , ns                               ! start conversion
        wf(i) = wa(i) + 2.0d+00*xa(i,i)           ! fundamental
          DO j = 1, ns
            IF ( i .NE. j ) THEN
            wf(i) = wf(i) + 0.5d+00*xa(i,j)       ! fundamental
            ENDIF
          END DO
      END DO                                      ! end conversion

      ELSEIF ( WW .EQ. 'Wf' ) THEN

      DO i = 1 , ns                               ! start conversion
        wa(i) = wf(i) - 2.d+00*xa(i,i)
          DO j = 1, ns
            IF ( i .NE. j ) THEN
            wa(i) = wa(i) - 0.5d+00*xa(i,j)       ! Eq. (II,273)
            ENDIF
          END DO
      END DO                                      ! end conversion
      DO i = 1 , ns                               ! start conversion
        w0(i) = wa(i) + xa(i,i)                   ! Eq. (II,273)
          DO j = 1, ns
            IF ( i .NE. j ) THEN
            w0(i) = w0(i) + 0.5d+00*xa(i,j)       ! Eq. (II,273)
            ENDIF
          END DO
      END DO                                      ! end conversion

      ENDIF
      
      zpe = 0.0d+00
      ave = 0.0d+00
      av0 = 0.0d+00
      avf = 0.0d+00
      Viblo = 10000.
      DO i = 1 , ns                           ! start ZPE
        zpe = zpe + 0.5d+00*wa(i)             ! Eq. (II,267)
        ave = ave + wa(i)                     ! summed frequency
        av0 = av0 + w0(i)                     ! summed frequency
        avf = avf + wf(i)                     ! summed frequency
        DO j = 1, i
          zpe = zpe + 0.25d+00*xa(j,i)        ! Eq. (II,267)
        END DO
        Viblo = MIN( Viblo , wa(i) )          ! find lowest frequency
      END DO                                  ! end 

      ave = ave/ns                            ! average frequency
      av0 = av0/ns                            ! average frequency
      avf = avf/ns                            ! average frequency


      RETURN
      END
