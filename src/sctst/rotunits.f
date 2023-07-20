c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c Thermo: a code for thermochemical calculations.
c Copyright (C) 2009 John R. Barker
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

      SUBROUTINE rotunits( VROTIN , X , VROTOUT )
      
      IMPLICIT NONE
      CHARACTER(len=4) VROTIN , VROTOUT     
      REAL(8) X
      SAVE
     
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

      END
      
      
