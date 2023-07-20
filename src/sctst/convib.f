c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c Copyright (C) 2009 John R. Barker and Thanh Lam Nguyen
c
c Authors: Thanh Lam Nguyen and John R. Barker
c          nguyenlt@umich.edu
c          September 2009
c
c John R. Barker
c jrbarker@umich.edu
c University of Michigan
c Ann Arbor, MI 48109-2143
c (734) 763 6239
c
c This program is free software; you can redistribute it and/or
c modify it under the terms of the GNU General Public License (version 2)
c as published by the Free Software Foundation.
c
c This program is distributed in the hope that it will be useful,
c but WITHOUT ANY WARRANTY; without even the implied warranty of
c MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
c GNU General Public License for more details.
c
c See the 'ReadMe' file for a copy of the GNU General Public License,
c or contact:
c
c Free Software Foundation, Inc.
c 59 Temple Place - Suite 330
c Boston, MA 02111-1307, USA.
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
      
      IMPLICIT NONE
      INTEGER(4) ns , i , j
      REAL(8) wa , w0 , wf , xa , ya , za , zpe, 
     &                 av0 , ave , avf , Viblo
      CHARACTER WW*2
      COMMON/wxyz/ wa(100) , w0(100) , wf(100) , xa(100,100) ,
     &   ya(100,100,100) , za(100,100,100,100) , zpe , ave ,
     &   av0 , avf , Viblo

      SAVE
      
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
