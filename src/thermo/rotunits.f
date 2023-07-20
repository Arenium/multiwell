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

      SUBROUTINE rotunits( VROT , X )
      
      IMPLICIT NONE
      CHARACTER(len=4) VROT      
      REAL(8) X
      SAVE
     
c      CONVERT ROTATIONAL UNITS TO amu*Ang^2
c
c      For qro and rot types [NOT HINDERED ROTORS]:
c          VROT = 'AMUA' for moment of inertia in units of amu*Ang^2
c       or      = 'GMCM' for moment of inertia units of gram*cm^2
c       or      = 'CM-1' for rotational constant in units of cm^-1
c       or      = 'MHZ' for rotational constant in units of MHz
c       or      = 'GHZ' for rotational constant in units of GHz
c         
      CALL ucase ( VROT )   ! convert to upper case
      IF ( VROT .EQ. 'AMUA' ) THEN
        RETURN
      ELSEIF ( VROT .EQ. 'GMCM' ) THEN
        X = X / 1.660538782d-040
        RETURN
      ELSEIF ( VROT .EQ. 'CM-1' ) THEN
        X = 16.85763D+00 / X
        RETURN
      ELSEIF ( VROT .EQ. 'MHZ' ) THEN
        X = 5.05379D+005 / X
        RETURN
      ELSEIF ( VROT .EQ. 'GHZ' ) THEN
        X = 5.05379D+002 / X
        RETURN
      ELSE
        write (*,*) 'FATAL: units (VROT) for rotations not recongized'
        STOP
      ENDIF

      END
      
      
