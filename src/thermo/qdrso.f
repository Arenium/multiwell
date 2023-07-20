c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c THERMO: a code for thermochemical calculations.
c Copyright (C) 2001-2019 John R. Barker
c
c SUBROUTINE qdro: a subroutine used by THERMO
c Copyright (C) 2019 John R. Barker
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
c See the "ReadMe" file for a copy of the GNU General Public License,
c or contact:
c
c Free Software Foundation, Inc.
c 59 Temple Place - Suite 330
c Boston, MA 02111-1307, USA.
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c    Combined Electronic + Rotation partition function 
c	
c	Implemented for DIATOMIC molecules with SPIN-ORBIT splitting
c	in 2Pi electronic states.
c
c	REFERENCES
c
c	E. L. Hill and J. H. Van Vleck, Phys. Rev. 32, 250 (1923)
c
c	as summarized in
c
c	G. Herzberg, Molecular Spectra and Molecular Structure. I. Spectra 
c	of Diatomic Molecules (Can Nostrand Reinhold Co.; 
c	copyright 1950 Litton Educationl Publishing, Inc.).
c
c	This approach is described mostly in Chapter V, Equation (V,28). 
c	Equation numbers here and in comments are from Herzberg.
c
c	***IMPORTANT***
c	When using this method, the electronic state degeneracy for this
c	species (entered on Line 9 in input file thermo.dat; see MultiWell 
c	User Manual) should be set equal to 1.
c
c	INPUT
c	T	= temperature
c	Aso	= spin orbit constant (cm-1)
c	B	= rotational constant (cm-1)
c	D	= centrifugal distortion constant (cm-1)
c	lamb	= Magnetic quantum number (Lambda) of the electronic state: 
c		  Sigma (Lamb=0), Pi (Lamb=1), Delta (Lamb=2), etc.
c	mult	= spin multiplicity of the electronic state.
c	sig  = rotation symmetry number
c   Nprnt = 1 (for printing rso levels), otherwise = 0
c
C	OUTPUT
c	Q0	= partition function for combined rotation + spin orbit
c	C	= constant p heat capacity (units of R)
c	S	= entropy (units of R)
c	H	= [H(T)-H(0)]/RT
c	Energy levels in file thermo.rso
c
c	*****  EXAMPLE PARAMETERS: OH(X 2Pi) *****
c	From Table 1 in:
c	Brooke et al., J Quant Spectros Rad Transfer 168 (2016) 142–157
c		http://dx.doi.org/10.1016/j.jqsrt.2015.07.021
c      Aso = -139.050877 cm-1
c      B = 18.53487308 cm-1
c      D = 1.9089352E-03 cm-1
c      mult = 2 
c      lamb = 1
c      sig = 1
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      SUBROUTINE qdrso(T,Aso,B,D,lamb,mult,sig,Q0,C,S,H,Nprnt) 
      IMPLICIT NONE
      REAL(8) Q0 , Q1 , Q2 , T , B , D , spin , omega , Aso ,
     & Y , aJ , aj2 , zero , con , test ,
     & C , S , H , Ered , E , hor , F , F1 , F2 , Fac 
      INTEGER(4) mult , lamb , sig , Nprnt
      PARAMETER (hor = 1.4387752D+00)
      SAVE 

c      IF (mult .NE. 2) THEN
c        write (*,*) 
c     &  'd.o.f. type "sor" implemented only for OH(X 2Pi)'
c        STOP
c      ENDIF 

      Y = Aso/B
      
      spin = ( mult - 1 )/2.0

c	Determine starting values for F1 and F2 loops
      IF ( Aso .GT. 0.0 ) THEN			! if "regular" case
        aj = lamb - spin
        aj2 = lamb + spin
      ELSE 							! if "inverted" case
        aj = lamb + spin
        aj2 = lamb - spin
      ENDIF

c	Zero of energy
      con = 0.5*SQRT( 4*(aj+0.5)**2 + Y*(Y-4)*lamb**2 )
      zero = B*( (aj+0.5)**2 - lamb**2 - con ) - D*aj**4    ! Expression for F1

      IF (Nprnt.EQ.1) THEN
        write(23,*) 'Zero of energy/cm-1 =', zero
      ENDIF 

c	Initialize
      Q0 = 0.0D+00
      Q1 = 0.0D+00
      Q2 = 0.0D+00
c
c     Now start F1 loop
      IF ( Nprnt .EQ. 1) THEN 
        write(23,*) '                F1'
        WRITE (23,*) ' J (total ang momentum)    Energy (cm-1)'
      ENDIF
      test = 1.0D+00
      DO WHILE ( test. GE. 1.0D-10 )
         con = 0.5*SQRT( 4*(aj+0.5)**2 + Y*(Y-4)*lamb**2 )
         F1 = B*( (aj+0.5)**2 - lamb**2 - con ) - D*aj**4 - zero
         Ered = hor*F1/T
         F = mult*( 2.d+00*aj + 1.d+00)
         Fac = F*EXP(-Ered)
         Q0 = Q0 + Fac
         Q1 = Q1 + Ered*Fac
         Q2 = Q2 + Ered*Ered*Fac
         IF ( aj .GT. 5 ) test = ABS(Fac/Q0)
         IF ( Nprnt .EQ. 1) write(23,9001) aj, F1
         aj = aj + 1.d+00
      ENDDO
c
c     Now start F2 loop
      IF ( Nprnt .EQ. 1) THEN 
        write(23,*) '                F2'
        WRITE (23,*) ' J (total ang momentum)    Energy (cm-1)'
      ENDIF
      aj = aj2
      test = 1.0D+00
      DO WHILE ( test. GE. 1.0D-10 )
         con = 0.5*SQRT( 4*(aj+0.5)**2 + Y*(Y-4)*lamb**2 )
         F2 = B*( (aj+0.5)**2 - lamb**2 + con ) - D*(aj+1)**4 - zero
         Ered = hor*F2/T
         F = mult*( 2.d+00*aj + 1.d+00)
         Fac = F*EXP(-Ered)
         Q0 = Q0 + Fac
         Q1 = Q1 + Ered*Fac
         Q2 = Q2 + Ered*Ered*Fac
        IF ( aj .GT. 5 ) test = ABS(Fac/Q0)
        IF ( Nprnt .EQ. 1) write(23,9001) aj, F2
        aj = aj + 1.d+00
      ENDDO
 
c      Thermo formulae from Pitzer and Gwinn, J. Chem. Phys. 10, 428 (1942)
 
      Q0 = Q0/sig
      Q1 = Q1/sig
      Q2 = Q2/sig
      C = Q2/Q0 - (Q1/Q0)**2
      H = Q1/Q0
      S = (H+log(Q0))
      
      IF (Nprnt.EQ.1) Nprnt = 0  			! now done with printing rso levels
      
9001  FORMAT(F5.1,20x,1pe11.4)
      
      RETURN
 
      END SUBROUTINE
