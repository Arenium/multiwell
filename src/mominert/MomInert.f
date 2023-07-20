 
      PROGRAM MOMINERT
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c MomInert: a code for moments of inertia.
c Copyright (C) 2010 John R. Barker
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
C     Designed to calculated moments of inertia for internal rotations.
C     These can either be connected to a larger mass object or stand-alone.
C
C     Input File Specifications:
C
C     1. Title
C             A line, up to 100 characters, describing the data.  This title is
C        reproduced in the output file.
C
C     2. Total number of atoms in the molecule, NATOMS.
C
C     3. Atom type, Atom index, X-coordinate, Y-coordinate, Z-coordinate
C        Atom type: Character*2, case sensitive atomic symbol, eg C, Br
C        Atom index: Identification number of the atom (1 to NATOMS)
C        X-coordinate, Y-coordinate, Z-coordinate in angstroms
C
C        Repeat for every atom.
C
C     4. Atom indices for the two atoms defining the axis of rotation. If the
C       atom indices are set equal to zero, then internal rotor is not calculated
C       and lines #5 and #6 can be omitted.
C
C     5. Number of atoms in one of the rotating moieties.
C
C     6. List of the atoms in that moiety.
C
C     7. If you wish to include more 1-D rotors, repeat steps 4-6 for each.
C        TO TERMINATE, ENTER A LINE CONTAINING TWO ZEROS: 0 , 0
C
C     8. BLANK LINE
C
C     Cartesian coordinates are in Angstroms, and the moments of inertia
C     are reported in amu*Ang*Ang.
 
      IMPLICIT NONE
 
      CHARACTER(len=6) AVERSION
      CHARACTER(len=8) ADATE
      PARAMETER ( AVERSION='2023', ADATE='Mar 2023' )
      INTEGER NATOMS , I , J , K , ROT1 , ROT2 , NUMROT , ROTING1(100)
      INTEGER ROTING2(100) , NUMROT2 , MOL , ATOM3 , N , NP , ierr , IV
      CHARACTER(len=1)   NFLAG
      CHARACTER(len=6)   ATYPE(100) , AT
      CHARACTER(len=100) TITLE
      DOUBLE PRECISION XCOORD(100) , YCOORD(100) , ZCOORD(100) , X1 , 
     &                 X2 , X3
      DOUBLE PRECISION Y1 , Y2 , Y3 , Z1 , Z2 , Z3 , R(100) , A1 , A2 , 
     &                 C1 , C2
      DOUBLE PRECISION D1 , D2 , A , C , DIST1(100) , DIST2(100) , X , 
     &                 Y , Z
      DOUBLE PRECISION MASS(100) , MOMENT1 , MOMENT2 , REDMOM , XC , 
     &                 YC , ZC
      DOUBLE PRECISION TOTMASS , XCM , YCM , ZCM , XPRIME(100) , 
     &                 YPRIME(100)
      DOUBLE PRECISION ZPRIME(100) , IXX , IYY , IZZ , IYZ , IXY , IXZ
      DOUBLE PRECISION MAT(3,3) , EIG(3) , fv1(3) , fv2(3) , ZZ(3,3)
      DOUBLE PRECISION TESTX12 , TESTX13 , TESTY12 , TESTY13 , TESTZ12 ,
     &                 TESTZ13 , CONV
      DOUBLE PRECISION XCM1 , YCM1 , ZCM1 , re
      INTEGER  FRGM1, FRGM2 , ncount
      DOUBLE PRECISION Hfxn
      CHARACTER(len=4) UNITS
      DOUBLE PRECISION scale
      CHARACTER(len=20) NameFile, NameOUT
 
c
c     FOR USER-DEFINED FILENAME ("NameFile")
c
      CALL get_command_argument( 1, NameFile )
      NameFile=NameFile(1:len_trim(NameFile))
      IF(NameFile.eq.'') NameFile="mominert.dat"              ! default filename
      NameOUT=NameFile(1:len_trim(NameFile)-4)//".out"
      OPEN (UNIT=2,FILE=NameFile,STATUS='OLD')
      OPEN (UNIT=4,FILE=NameOUT,STATUS='UNKNOWN')
 
      READ (2,99001) TITLE
      READ (2,*) UNITS
      IF ( (UNITS .EQ. 'angs') .OR. (UNITS .EQ. 'ANGS') ) THEN
         UNITS = 'ANGS'
         scale = 1.0d+00
      ELSEIF ( (UNITS .EQ. 'bohr') .OR. (UNITS .EQ. 'BOHR') ) THEN
         scale = 5.2917720859d-01
      ELSE
         WRITE (4,*)
         WRITE (4,*) '***FATAL: UNITS must be ''BOHR'' or ''ANGS'' ***'
         WRITE (*,*)
         WRITE (*,*) '***FATAL: UNITS must be ''BOHR'' or ''ANGS'' ***'
         STOP
      ENDIF
!-------------------------------------------------------------------
      READ (2,*,err=50) NATOMS, FRGM1, FRGM2
      READ (2,*,err=50) NATOMS
      IF(FRGM1+FRGM2.NE.NATOMS) then
       write(*,*) " Wrong fragments definition. Program stopped"
       stop
      ENDIF 
      goto 60
50    backspace(2)
      backspace(2)
      READ (2,*) NATOMS
!-------------------------------------------------------------------
 
60    WRITE (4,99002) AVERSION , ADATE , AVERSION , ADATE
      WRITE (4,*) TITLE
      IF ( UNITS .NE. 'ANGS' ) THEN
         WRITE (4,*)
         WRITE(4,*) '( Converted units from ',UNITS,' to ANGSTROMS )'
      ENDIF
      WRITE (4,*)
      WRITE (4,99003) 'NUMBER' , 'TYPE  ' , 'X' , 'Y' , 'Z'
  
C     Here zero all arrays to avoid external interference with calculations.
 
      DO 100 I = 1 , 100
         DIST1(I) = 0.D+00
         DIST2(I) = 0.D+00
         XCOORD(I) = 0.D+00
         YCOORD(I) = 0.D+00
         ZCOORD(I) = 0.D+00
         R(I) = 0.D+00
         ROTING1(I) = 0
         ROTING2(I) = 0
 100  CONTINUE
 
      DO 200 I = 1 , NATOMS
         READ (2,*) AT , MOL , XC , YC , ZC
         ATYPE(MOL) = AT
         XCOORD(MOL) = XC*scale
         YCOORD(MOL) = YC*scale
         ZCOORD(MOL) = ZC*scale
 200  CONTINUE
 
      DO 300 I = 1 , NATOMS
         WRITE (4,99004) I , ATYPE(I) , XCOORD(I) , YCOORD(I) , 
     &                   ZCOORD(I)
 300  CONTINUE
 
C     Define the mass of each atom.
C     This is done with simple character recognition, and the masses are
C     given in units of amu.
 
      DO 400 I = 1 , NATOMS
        CALL elemass( ATYPE(I), MASS(I) )  
 400  CONTINUE
 
       NFLAG = 'Y'
 
C     Read atoms which define rotation axis and report them to user.
 
      DO WHILE (NFLAG .EQ. 'Y' )
      
      WRITE (4,*)
      READ (2,*) ROT1 , ROT2

      IF ( ROT1 .EQ. 0 ) THEN
         NFLAG = 'N'
         
      ELSE
         WRITE (4,99005) ROT1 , ROT2 
         READ (2,*) NUMROT
 
C     Read atoms which make up the first rotating moiety and report them.
C     Reporting allows for ease in finding mistakes.
 
         WRITE (4,*)
         READ (2,*) (ROTING1(I),I=1,NUMROT)

99941  FORMAT('  MOIETY 1 CONTAINS ATOMS:',50I3)

         WRITE (4,99941) (ROTING1(I),I=1,NUMROT)
 
         WRITE (4,*)
 
C     All atoms which are not in the first moiety must be in the second
C     moiety.  Define it here.
 
         NUMROT2 = NATOMS - NUMROT
         K = 1
         DO 450 I = 1 , NATOMS
            DO 420 J = 1 , NUMROT
               IF ( I.EQ.ROTING1(J) ) GOTO 450
 420        CONTINUE
            ROTING2(K) = I
            K = K + 1
 450     CONTINUE
 
99942  FORMAT('  MOIETY 2 CONTAINS ATOMS:',50I3)

         WRITE (4,99942) (ROTING2(I),I=1,NUMROT2)
 
         WRITE (4,*)

c
c
c
C     Some quick idiot lights...
         
         IF ( ROT1.LT.0 ) THEN
            WRITE (4,*) 'ROTATION AXIS UNDEFINED'
            STOP
         ENDIF
 
         IF ( ROT1.GT.NATOMS ) THEN
            WRITE (4,*) 'ROTATION AXIS UNDEFINED'
            STOP
         ENDIF
 
         IF ( ROT2.LE.0 ) THEN
            WRITE (4,*) 'ROTATION AXIS UNDEFINED'
            STOP
         ENDIF
 
         IF ( ROT2.GT.NATOMS ) THEN
            WRITE (4,*) 'ROTATION AXIS UNDEFINED'
            STOP
         ENDIF
 
         IF ( NUMROT.GT.(NATOMS-1) ) THEN
            WRITE (4,*) 'ROTATING MOIETY IMPROPERLY DEFINED'
            STOP
         ENDIF
 
         DO 500 I = 1 , NUMROT
            IF ( ROTING1(I).LE.0 ) THEN
               WRITE (4,*) 'ROTATING MOIETY IMPROPERLY DEFINED'
               STOP
            ENDIF
 
            IF ( ROTING1(I).GT.NATOMS ) THEN
               WRITE (4,*) 'ROTATING MOIETY IMPROPERLY DEFINED'
               STOP
            ENDIF
 
 500     CONTINUE

C     Begin defining planes as used to determine radial distance to
C     axis of rotation.  Use the coord.s of each axis-defining atom to
C     start, but you need one more point. Take the first atom in the list
C     which is neither ROT1 nor ROT2 to be the other defining point for
C     the first plane.
 
         DO 550 I = 1 , NATOMS
            IF ( I.NE.ROT1 ) THEN
               IF ( I.NE.ROT2 ) GOTO 600
            ENDIF
 550     CONTINUE
 
 600     ATOM3 = I
 
C     If (X2-X1) or (X3-X1) or (Z2-Z1) or (Z3-Z1) is equal to zero, the
C     definition of the planes as here will not work.  So, if this is true,
C     send the coordinates to a subroutine which will rotate the molecule
C     while  leaving it unchanged.
c
c     If one rotation isn't effective, do it over, and over, and over...
 
         ncount = 0
         CONV = 1.0D-02
650      ncount = ncount + 1
         IF ( ncount .GT. 10 ) go to 651
         
         X1 = XCOORD(ROT1)
         Y1 = YCOORD(ROT1)
         Z1 = ZCOORD(ROT1)
 
         X2 = XCOORD(ROT2)
         Y2 = YCOORD(ROT2)
         Z2 = ZCOORD(ROT2)
 
         X3 = XCOORD(ATOM3)
         Y3 = YCOORD(ATOM3)
         Z3 = ZCOORD(ATOM3)
 
         TESTX12 = ABS(X1-X2)
         TESTX13 = ABS(X1-X3)
         TESTY12 = ABS(Y1-Y2)
         TESTY13 = ABS(Y1-Y3)
         TESTZ12 = ABS(Z1-Z2)
         TESTZ13 = ABS(Z1-Z3)
         
         IF ( (TESTX12.LE.CONV) .OR. (TESTX13.LE.CONV) .OR. 
     &        (TESTY13.LE.CONV) .OR. (TESTY12.LE.CONV) .OR.
     &        (TESTZ13.LE.CONV) .OR. (TESTZ12.LE.CONV) ) THEN
            CALL ROTATE(XCOORD,YCOORD,ZCOORD,NATOMS)
C            WRITE (*,*) 'Rotations: ', ncount
            GOTO 650
         ENDIF
         
651      CONTINUE
 
C     Define each plane according to 'Direction Numbers' A, B, C, and D
C     such that Ax+By+Cz+D=0.  Start by fitting the plane through the
C     three atoms selected.  That's why only one can have coordinates of
C     zero.  For simplicity, B=1.
 
         A = ((Y3-Y1)/(Z3-Z1)) - ((Y2-Y1)/(Z2-Z1))
         A1 = A/(((X2-X1)/(Z2-Z1))-((X3-X1)/(Z3-Z1)))
 
         C = ((Y3-Y1)/(X3-X1)) - ((Y2-Y1)/(X2-X1))
         C1 = C/(((Z2-Z1)/(X2-X1))-((Z3-Z1)/(X3-X1)))
 
         D1 = -(A1*X1) - Y1 - (C1*Z1)
 
C     Now, using the two axis-defining atoms and the other plane (which
C     is perpendicular to this one) define the second plane.  The
C     intersection of these two planes is the rotation axis in three
C     dimensions.  Again, B=1.
 
         A = (1/C1) - ((Y2-Y1)/(Z2-Z1))
         A2 = A/(((X2-X1)/(Z2-Z1))-(A1/C1))
 
         C = (1/A1) - ((Y2-Y1)/(X2-X1))
         C2 = C/(((Z2-Z1)/(X2-X1))-(C1/A1))
 
         D2 = -(A2*X1) - Y1 - (C2*Z1)
         
         
c         write(*,*) A1 , C1, D1, A2, C2, D2
         
C     Find the distance from each atom to each of the two planes.  Since
C     the planes are perpendicular, the distance to the intersection is
C     the hypoteneuse of the right triangle.  It's not necessary, but
C     to ensure positive radial distances take the absolute values of the
C     distance to each plane.
C     Then use the Pythagorean Theorem to calculate the radial distance of
C     each atom to the rotation axis.
 
         DO 700 I = 1 , NATOMS
            X = XCOORD(I)
            Y = YCOORD(I)
            Z = ZCOORD(I)
            DIST1(I) = ((A1*X)+Y+(C1*Z)+D1)/(SQRT((A1*A1)+1+(C1*C1)))
            DIST1(I) = ABS(DIST1(I))
            DIST2(I) = ((A2*X)+Y+(C2*Z)+D2)/(SQRT((A2*A2)+1+(C2*C2)))
            DIST2(I) = ABS(DIST2(I))
            R(I) = SQRT((DIST1(I)**2)+(DIST2(I)**2))
700     CONTINUE
 
C     Compute the moment of inertia for the first rotor,  using the
C     formula Im=sum(m*r*r)
 
         MOMENT1 = 0
         DO 750 I = 1 , NUMROT
            MOMENT1 = MOMENT1 + (MASS(ROTING1(I))*R(ROTING1(I))
     &                *R(ROTING1(I)))
 750     CONTINUE
 
         WRITE (4,99006) MOMENT1
         WRITE (4,99906) MOMENT1*1.66053878d-040
         WRITE (4,99996) 16.85763D+00/MOMENT1, 5.05379D+005/MOMENT1,
     &          5.05379D+002/MOMENT1
 
         WRITE (4,*)
 
C     Compute the moment of inertia for the second rotor,  using the
C     formula Im=sum(m*r*r)
 
         MOMENT2 = 0
         DO 800 I = 1 , NUMROT2
            MOMENT2 = MOMENT2 + (MASS(ROTING2(I))*R(ROTING2(I))
     &                *R(ROTING2(I)))
 800     CONTINUE
 
         WRITE (4,99007) MOMENT2
         WRITE (4,99906) MOMENT2*1.660538782d-040
         WRITE (4,99996) 16.85763D+00/MOMENT2, 5.05379D+005/MOMENT2,
     &          5.05379D+002/MOMENT2
 
         WRITE (4,*)
 
         REDMOM = (MOMENT1*MOMENT2)/(MOMENT1+MOMENT2)
         WRITE (4,99008) ROT1, ROT2, REDMOM
         WRITE (4,99906) REDMOM*1.660538782d-040
         WRITE (4,99996) 16.85763D+00/REDMOM, 5.05379D+005/REDMOM,
     &          5.05379D+002/REDMOM

      ENDIF
      
      enddo

999   TOTMASS = 0.D+00
      XCM = 0.D+00
      YCM = 0.D+00
      ZCM = 0.D+00
      IXX = 0.D+00
      IYY = 0.D+00
      IZZ = 0.D+00
      IXY = 0.D+00
      IYZ = 0.D+00
      IXZ = 0.D+00
 
C     Calculate center-of-mass of molecule
 
      DO 900 I = 1 , NATOMS
         TOTMASS = TOTMASS + MASS(I)
         XCM = XCM + (MASS(I)*XCOORD(I))
         YCM = YCM + (MASS(I)*YCOORD(I))
         ZCM = ZCM + (MASS(I)*ZCOORD(I))
 900  CONTINUE
 
      XCM = XCM/TOTMASS
      YCM = YCM/TOTMASS
      ZCM = ZCM/TOTMASS

C     Use existing coordinates and center-of-mass coordinates to find
C     new, 'internal' coordinates x', y', z'
 
      DO 1000 I = 1 , NATOMS
         XPRIME(I) = XCOORD(I) - XCM
         YPRIME(I) = YCOORD(I) - YCM
         ZPRIME(I) = ZCOORD(I) - ZCM
      
C    Use x', y', and z' to calculate moments and products of inertia

         IXX = IXX + (MASS(I)*(YPRIME(I)**2+ZPRIME(I)**2))
         IYY = IYY + (MASS(I)*(XPRIME(I)**2+ZPRIME(I)**2))
         IZZ = IZZ + (MASS(I)*(XPRIME(I)**2+YPRIME(I)**2))
         IXY = IXY + MASS(I)*XPRIME(I)*YPRIME(I)
         IXZ = IXZ + MASS(I)*XPRIME(I)*ZPRIME(I)
         IYZ = IYZ + MASS(I)*YPRIME(I)*ZPRIME(I)
 1000 CONTINUE
 
C     Compose matrix of moments and products of inertia
      MAT(1,1) = IXX
      MAT(1,2) = -IXY
      MAT(1,3) = -IXZ
      MAT(2,1) = -IXY
      MAT(2,2) = IYY
      MAT(2,3) = -IYZ
      MAT(3,1) = -IXZ
      MAT(3,2) = -IYZ
      MAT(3,3) = IZZ
 
      N = 3
      NP = 3
      IV = 0
 
C     Find eigenvalues
 
      CALL rs(N,NP,MAT,EIG,IV,ZZ,fv1,fv2,ierr)
 
      WRITE (4,99010) ( EIG(i) , i=1,3 ) ,
     &     ( EIG(i)*1.660538782d-040 , i=1,3 ) , 
     &     ( 16.85763D+00/EIG(i) , i=1,3 ) ,
     &     ( 5.05379D+005/EIG(i) , i=1,3 ) ,
     &     ( 5.05379D+002/EIG(i) , i=1,3 ) 


!----------------------------------------------------------------
c     Distance between center of mass

       IF(FRGM1.GT.0) then
        CALL CMdist(XCOORD,YCOORD,ZCOORD,MASS,FRGM1,FRGM2,NATOMS,re,
     & XCM1,XCM,YCM1,YCM,ZCM1,ZCM)
        WRITE(4,99012) re
        WRITE(4,99013)
        WRITE(4,99014) "   MOIETY  1 ( 1-",FRGM1,") :   ",XCM,YCM,ZCM
        WRITE(4,99015) "   MOIETY  2 (",FRGM1+1,"-",NATOMS,") :   ",
     & XCM1,YCM1,ZCM1
       ENDIF

!----------------------------------------------------------------
!
!      WRITE(4,99011) NATOMS-FRGM2, NATOMS-FRGM2+1, NATOMS 
!      WRITE(4,99012) re
!----------------------------------------------------------------

      CLOSE (2)
      CLOSE (4)
 
99001 FORMAT (A100)
99002 FORMAT (
     &'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     &%%%%%%%%%%%%%'//8x,
     &'                                John R. Barker'/8x,
     &'       Mominert-',A6,'          University of Michigan'/8x,
     &'                                Ann Arbor, MI 48109-2143'/8x,
     &'       ',A8,'                 jrbarker@umich.edu'/8x,
     &'                                (734) 763 6239'//8x,
     &'      http://clasp-research.engin.umich.edu/multiwell/'//,
     &'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     &%%%%%%%%%%%%%'//8x,
     &'Suggested Literature Citations:'//4x,
     &'a) J.R. Barker, T.L. Nguyen, J.F. Stanton, C. Aieta, M. Ceotto,',
     &/7x,'F. Gabas, T.J.D. Kumar, C.G.L.Li, L.L.Lohr, A. Maranzana,',
     &/7x,'N.F. Ortiz, J.M. Preses, J.M. Simmie, J.A. Sonk, and ',
     &'P.J. Stimac; ',
     &/7x,'MultiWell-',A7,' Software Suite; University of Michigan, ',
     &/7x,'Ann Arbor, Michigan, USA, ',A8,'; ',
     &/7x,'http://clasp-research.engin.umich.edu/multiwell/.',
     &//4x,'b) J.R. Barker, Int. J. Chem. Kinetics, 33, 232-45 (2001)',
     &//,
     &'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     &%%%%%%%%%%%%%'/)
99003 FORMAT (1X,A6,2X,A6,7X,A1,2(11X,A1))
99004 FORMAT (2X,I3,5X,A6,3(3X,F9.6))
99005 FORMAT (//,2X,'ATOMS ',I3,' AND ',I3,
     &         ' DEFINE INTERNAL ROTATION AXIS')
99006 FORMAT (2X,'MOMENT OF INERTIA ABOUT BOND FOR MOIETY #1: ',
     &            F9.3,'   amu*Ang^2')
99906 FORMAT (45X,1PE10.3,'   g*cm^2')
99996 FORMAT (14X,'Rotational constant:  B = ',F15.5,'   cm-1',/,
     &        14X,'                      B = ',F15.5,'   MHz',/ ,
     &        14X,'                      B = ',F15.5,'   GHz',/ )
99007 FORMAT (2X,'MOMENT OF INERTIA ABOUT BOND FOR MOIETY #2: ',
     &            F9.3,'   amu*Ang^2')
99008 FORMAT (2X,'REDUCED MOMENT OF INERTIA ABOUT ',I2,'--',I2,' BOND:',
     &         F9.3,'   amu*Ang^2')
99009 FORMAT (5X,A31,1X,I3,1X,A3,1X,I3)
99010 FORMAT (/,'PRINCIPAL MOMENTS OF INERTIA AND ROTATIONAL CONSTANTS',
     & /,
     & 2x,'Ia = ',F12.4,3x,'Ib = ',F12.4,3x,'Ic = ',F12.4,2x,
     & ' amu*ang^2',/,
     & 2x,'Ia = ',1pe12.4,3x,'Ib = ',1pe12.4,3x,'Ic = ',1pe12.4,2x,
     & ' g*cm^2',/,
     & 2x,'Ba = ',0pF12.4,3x,'Bb = ',0pF12.4, 3x,'Bc = ',0pF12.4,2x,
     & ' cm-1',/,
     & 2x,'Ba = ',0pF12.4,3x,'Bb = ',0pF12.4, 3x,'Bc = ',0pF12.4,2x,
     & ' MHz',/,
     & 2x,'Ba = ',0pF12.4,3x,'Bb = ',0pF12.4, 3x,'Bc = ',0pF12.4,2x,
     & ' GHz',/ )
99012 FORMAT (//,"DISTANCE BETWEEN CENTERS OF MASS = ", F11.6," Ang"/)
99013 FORMAT ("   CENTER OF MASS COORDINATES: ")
99014 FORMAT (A16,I2,A6,T23,3F12.6)
99015 FORMAT (A14,I2,A1,I2,A6,T23,3F12.6,//)
 
      END
