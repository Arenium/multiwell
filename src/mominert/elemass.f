c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c Element: a subroutine for calculating atomic masses and enthalpy functions.
c Copyright (C) 2011 John R. Barker
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

      SUBROUTINE elemass( ATYPE , MASS )
      
      DOUBLE PRECISION MASS
      CHARACTER(len=6) ATYPE
      SAVE
c
c      MASS (atomic mass units) from NIST / JANAF 1998
c
c      Updated from data on NIST web site (accessed Feb 24, 2011): 
c            Atomic Weights and Isotopic Compositions
c            http://www.nist.gov/pml/data/comp.cfm
c
         IF ( ATYPE.EQ.'H' ) THEN
            MASS = 1.00794D+00
         ELSEIF ( ATYPE.EQ.'H1' ) THEN
            MASS = 1.007825D+00
         ELSEIF ( ATYPE.EQ.'D' ) THEN
            MASS = 2.014102D+00
         ELSEIF ( ATYPE.EQ.'T' ) THEN
            MASS = 3.016049D+00
         ELSEIF ( ATYPE.EQ.'He' ) THEN
            MASS = 4.00260D+00
         ELSEIF ( ATYPE.EQ.'Li' ) THEN
            MASS = 6.941D+00
         ELSEIF ( ATYPE.EQ.'Li6' ) THEN
            MASS = 6.0151228D+00
         ELSEIF ( ATYPE.EQ.'Li7' ) THEN
            MASS = 7.0160046D+00
         ELSEIF ( ATYPE.EQ.'Be' ) THEN
            MASS = 9.012182D+00
         ELSEIF ( ATYPE.EQ.'B' ) THEN
            MASS = 10.81D+00
         ELSEIF ( ATYPE.EQ.'B10' ) THEN
            MASS = 10.0129D+00
         ELSEIF ( ATYPE.EQ.'B11' ) THEN
            MASS = 11.0093D+00
         ELSEIF ( ATYPE.EQ.'C' ) THEN
            MASS = 12.011D+00
         ELSEIF ( ATYPE.EQ.'C12' ) THEN
            MASS = 12.0D+00
         ELSEIF ( ATYPE.EQ.'C13' ) THEN
            MASS = 13.003355D+00
         ELSEIF ( ATYPE.EQ.'C14' ) THEN
            MASS = 14.003242D+00
         ELSEIF ( ATYPE.EQ.'C16' ) THEN
            MASS = 16.0147D+00
         ELSEIF ( ATYPE.EQ.'N' ) THEN
            MASS = 14.00674D+00
         ELSEIF ( ATYPE.EQ.'N14' ) THEN
            MASS = 14.003074D+00
         ELSEIF ( ATYPE.EQ.'N15' ) THEN
            MASS = 15.00010897D+00
         ELSEIF ( ATYPE.EQ.'O' ) THEN
            MASS = 15.9994D+00
         ELSEIF ( ATYPE.EQ.'O16' ) THEN
            MASS = 15.994915D+00
         ELSEIF ( ATYPE.EQ.'O17' ) THEN
            MASS = 16.9991315D+00
         ELSEIF ( ATYPE.EQ.'O18' ) THEN
            MASS = 17.999160D+00
         ELSEIF ( ATYPE.EQ.'F' ) THEN
            MASS = 18.998403D+00
         ELSEIF ( ATYPE.EQ.'Ne' ) THEN
            MASS = 20.1797D+00
         ELSEIF ( ATYPE.EQ.'Na' ) THEN
            MASS = 22.98977D+00
         ELSEIF ( ATYPE.EQ.'Mg' ) THEN
            MASS = 24.305D+00
         ELSEIF ( ATYPE.EQ.'Mg24' ) THEN
            MASS = 23.985D+00
         ELSEIF ( ATYPE.EQ.'Mg25' ) THEN
            MASS = 24.985D+00
         ELSEIF ( ATYPE.EQ.'Mg26' ) THEN
            MASS = 25.983D+00
         ELSEIF ( ATYPE.EQ.'Al' ) THEN
            MASS = 26.9815386D+00
         ELSEIF ( ATYPE.EQ.'Si' ) THEN
            MASS = 28.0855D+00
         ELSEIF ( ATYPE.EQ.'Si28' ) THEN
            MASS = 27.9769D+00
         ELSEIF ( ATYPE.EQ.'Si29' ) THEN
            MASS = 28.98D+00
         ELSEIF ( ATYPE.EQ.'Si30' ) THEN
            MASS = 29.97D+00
         ELSEIF ( ATYPE.EQ.'P' ) THEN
            MASS = 30.97376D+00
         ELSEIF ( ATYPE.EQ.'S' ) THEN
            MASS = 32.065D+00
         ELSEIF ( ATYPE.EQ.'S32' ) THEN
            MASS = 31.9720707D+00
         ELSEIF ( ATYPE.EQ.'S33' ) THEN
            MASS = 32.97145876D+00
         ELSEIF ( ATYPE.EQ.'S34' ) THEN
            MASS = 33.9678668D+00
         ELSEIF ( ATYPE.EQ.'Cl' ) THEN
            MASS = 35.453D+00
         ELSEIF ( ATYPE.EQ.'Cl35' ) THEN
            MASS = 34.96885271D+00
         ELSEIF ( ATYPE.EQ.'Cl37' ) THEN
            MASS = 36.9659D+00
         ELSEIF ( ATYPE.EQ.'Ar' ) THEN
            MASS = 39.948D+00
         ELSEIF ( ATYPE.EQ.'K' ) THEN
            MASS = 39.0983D+00
         ELSEIF ( ATYPE.EQ.'K39' ) THEN
            MASS = 38.96371D+00
         ELSEIF ( ATYPE.EQ.'K41' ) THEN
            MASS = 40.961826D+00
         ELSEIF ( ATYPE.EQ.'Ca' ) THEN
            MASS = 40.078D+00
         ELSEIF ( ATYPE.EQ.'Ti' ) THEN
            MASS = 47.867D+00
         ELSEIF ( ATYPE.EQ.'V' ) THEN
            MASS = 50.9415D+00
         ELSEIF ( ATYPE.EQ.'Cr' ) THEN
            MASS = 51.9961D+00
         ELSEIF ( ATYPE.EQ.'Fe' ) THEN
            MASS = 55.845D+00
         ELSEIF ( ATYPE.EQ.'Ni' ) THEN
            MASS = 58.6934D+00
         ELSEIF ( ATYPE.EQ.'Cu' ) THEN
            MASS = 63.546D+00
         ELSEIF ( ATYPE.EQ.'Zn' ) THEN
            MASS = 65.38D+00
         ELSEIF ( ATYPE.EQ.'Zn64' ) THEN
            MASS = 63.929D+00
         ELSEIF ( ATYPE.EQ.'Zn66' ) THEN
            MASS = 65.926D+00
         ELSEIF ( ATYPE.EQ.'Zn68' ) THEN
            MASS = 67.925D+00
         ELSEIF ( ATYPE.EQ.'Ga' ) THEN
            MASS = 69.723D+00
         ELSEIF ( ATYPE.EQ.'Ge' ) THEN
            MASS = 72.64D+00
         ELSEIF ( ATYPE.EQ.'Se' ) THEN
            MASS = 78.96D+00
         ELSEIF ( ATYPE.EQ.'Se78' ) THEN
            MASS = 77.917D+00
         ELSEIF ( ATYPE.EQ.'Se80' ) THEN
            MASS = 79.917D+00
         ELSEIF ( ATYPE.EQ.'Br' ) THEN
            MASS = 79.904D+00
         ELSEIF ( ATYPE.EQ.'Br79' ) THEN
            MASS = 78.918338D+00
         ELSEIF ( ATYPE.EQ.'Br81' ) THEN
            MASS = 80.916291D+00
         ELSEIF ( ATYPE.EQ.'Kr' ) THEN
            MASS = 83.80D+00
         ELSEIF ( ATYPE.EQ.'Rb' ) THEN
            MASS = 85.4678D+00
         ELSEIF ( ATYPE.EQ.'Sr' ) THEN
            MASS = 87.62D+00
         ELSEIF ( ATYPE.EQ.'Sn' ) THEN
            MASS = 118.710D+00
         ELSEIF ( ATYPE.EQ.'Te' ) THEN
            MASS = 127.60D+00
         ELSEIF ( ATYPE.EQ.'I ' ) THEN
            MASS = 126.904473D+00
         ELSEIF ( ATYPE.EQ.'Xe' ) THEN
            MASS = 131.29D+00
         ELSEIF ( ATYPE.EQ.'Cs' ) THEN
            MASS = 132.90545D+00
         ELSEIF ( ATYPE.EQ.'Ba' ) THEN
            MASS = 137.327D+00
         ELSEIF ( ATYPE.EQ.'Hg' ) THEN
            MASS = 200.59D+00
         ELSEIF ( ATYPE.EQ.'U' ) THEN
            MASS = 238.0289D+00
         ELSEIF ( ATYPE.EQ.'U235' ) THEN
            MASS = 235.04393D+00
         ELSEIF ( ATYPE.EQ.'U238' ) THEN
            MASS = 238.050788D+00
         ELSE
            STOP 'Element not recognized in subroutine elemmass.f'
         ENDIF

      RETURN
      END
      
      
