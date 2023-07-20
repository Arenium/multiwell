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
c----------------------------------------------------------------------
c
c      Read vibrational frequencies and anharmonicities
c
c       KIN = input UNIT for READ statements
c       ns  = number of vibrational frequencies
c       WW  = keyword ('we', 'w0', or 'wf') for the three types of frequencies:
c       wa(i): harmonic frequencies with zero of energy at the potential minimum
c          Evib = SUMi( wai*(vi+1/2) ) + SUMi( SUMj( Xij*(vi+1/2)*(vj+1/2) ) )
c       w0(i): frequencies with zero of energy at zpe
c          Evib = SUMi( w0i*vi ) + SUMi( SUMj( Xij*vi*vj ) )
c       wf(i) = fundamental frequencies (for 0-1 transitions with all vj = 0)
c----------------------------------------------------------------------
      SUBROUTINE READ_WXYZ( KIN , ns , NY, NZ, WW ) 
      
      IMPLICIT NONE
      CHARACTER WW*2 , KEYWORD*5
      INTEGER(4) KIN , ns , i , j , k , l
      INTEGER(4) NY , NZ , I1 , J1 , K1 , L1
      REAL(8) wa , w0 , wf , xa , ya , za
      REAL(8) ave , av0 , avf , zpe , Viblo
      COMMON/wxyz/ wa(100) , w0(100) , wf(100) , xa(100,100) ,
     &   ya(100,100,100) , za(100,100,100,100) , zpe , ave ,
     &   av0 , avf , Viblo
      SAVE

      IF ( (WW .EQ. 'We') .OR. (WW .EQ. 'we') .OR. (WW .EQ. 'WE') ) THEN
        READ(KIN,*)     (wa(i), i=1,ns)
        WW = 'We'
      ELSEIF ( (WW .EQ. 'W0') .OR. (WW .EQ. 'w0') ) THEN
        READ(KIN,*)     (w0(i), i=1,ns)
        WW = 'W0'
      ELSEIF ( (WW.EQ.'Wf') .OR. (WW.EQ.'wf') .OR. (WW.EQ.'WF') ) THEN
        READ(KIN,*)     (wf(i), i=1,ns)
        WW = 'Wf'
      ENDIF
c                                                        X matrix (vibrational anharmonicity)
      READ(KIN,9011) KEYWORD
      IF ( (KEYWORD .EQ. 'UPPER') .OR. (KEYWORD .EQ. 'upper') 
     &       .OR. (KEYWORD .EQ. 'Upper') ) THEN
        DO j = 1, ns
          READ(KIN,*) (xa(j,i), i=j, ns)
          IF ( j .LT. ns) THEN
            DO k = j+1 , ns
              xa(k,j) = xa(j,k)
            END DO
          ENDIF
        ENDDO
      ELSEIF ( (KEYWORD .EQ. 'LOWER') .OR. (KEYWORD .EQ. 'lower') 
     &           .OR. (KEYWORD .EQ. 'Lower') ) THEN
        DO j = 1, ns
          READ(KIN,*) (xa(j,i), i=1,j)
          IF ( j .GT. 1) THEN
            DO k = 1 , j-1
              xa(k,j) = xa(j,k)
            END DO
          ENDIF
        ENDDO
      ENDIF
c                                                        Y matrix (vibrational anharmonicity)
       DO i=1, ns
         DO j=1, ns
             DO k=1, ns
                ya(i,j,k)=0.0d0
             ENDDO
         ENDDO
       ENDDO

        DO i=1, NY
         READ(KIN,*) I1, J1, K1, ya(I1,J1,K1)
         ya(K1,J1,I1)=ya(I1,J1,K1)
         ya(K1,I1,J1)=ya(I1,J1,K1)
         ya(J1,I1,K1)=ya(I1,J1,K1)
         ya(J1,K1,I1)=ya(I1,J1,K1)
         ya(I1,K1,J1)=ya(I1,J1,K1)
        ENDDO 
c                                                        Z matrix (vibrational anharmonicity)
        DO i=1, ns
          DO j=1, ns
            DO k=1, ns    
              DO l=1, ns
                za(i,j,k,l)=0.0d0
              ENDDO
            ENDDO
          ENDDO
        ENDDO

        DO i=1, NZ
         READ(KIN,*) I1, J1, K1, L1, za(I1,J1,K1,L1)
         za(I1,J1,L1,K1)=za(I1,J1,K1,L1)
         za(I1,K1,J1,L1)=za(I1,J1,K1,L1)
         za(I1,K1,L1,J1)=za(I1,J1,K1,L1)
         za(I1,L1,K1,J1)=za(I1,J1,K1,L1)
         za(I1,L1,J1,K1)=za(I1,J1,K1,L1)
         za(J1,I1,K1,L1)=za(I1,J1,K1,L1)
         za(J1,I1,L1,K1)=za(I1,J1,K1,L1)
         za(J1,K1,I1,L1)=za(I1,J1,K1,L1)
         za(J1,K1,L1,I1)=za(I1,J1,K1,L1)
         za(J1,L1,I1,K1)=za(I1,J1,K1,L1)
         za(J1,L1,K1,I1)=za(I1,J1,K1,L1)
         za(K1,I1,J1,L1)=za(I1,J1,K1,L1)
         za(K1,I1,L1,J1)=za(I1,J1,K1,L1)
         za(K1,J1,I1,L1)=za(I1,J1,K1,L1)
         za(K1,J1,L1,I1)=za(I1,J1,K1,L1)
         za(K1,L1,I1,J1)=za(I1,J1,K1,L1)
         za(K1,L1,J1,I1)=za(I1,J1,K1,L1)
         za(L1,I1,J1,K1)=za(I1,J1,K1,L1)
         za(L1,I1,K1,J1)=za(I1,J1,K1,L1)
         za(L1,J1,I1,K1)=za(I1,J1,K1,L1)
         za(L1,J1,K1,I1)=za(I1,J1,K1,L1)
         za(L1,K1,I1,J1)=za(I1,J1,K1,L1)
         za(L1,K1,J1,I1)=za(I1,J1,K1,L1)
        ENDDO

9011  FORMAT(A5)

      RETURN
      END
      
      
