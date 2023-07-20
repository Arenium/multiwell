c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c Element: a subroutine for calculating atomic masses and enthalpy functions.
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

      SUBROUTINE element( ATYPE , MASS , Hfxn , Sfxn)
      
      REAL(8) MASS , Hfxn , R , Sf , Sfxn
      CHARACTER(len=6) ATYPE
      
      parameter (R=8.314472D+00)  ! J k-1 mol-1

      SAVE
c
c         MASS (atomic mass units) from NIST / JANAF 1998
c         Hfxn = [H(298.15) - H(0)] kJ mol-1 per atom from ref. state table for each element  [NIST / JANAF]
c         Sf   = entropy (per atom; J K-1 mol-1) of reference state at 298.15 K [NIST / JANAF]
c
         IF ( ATYPE.EQ.'Al' ) THEN
            MASS = 26.98154D+00
            Hfxn = 4.539
            Sf = 28.275
         ELSEIF ( ATYPE.EQ.'Ar' ) THEN
            MASS = 39.948D+00
            Hfxn = 6.197
            Sf = 154.845
         ELSEIF ( ATYPE.EQ.'B' ) THEN
            MASS = 10.81D+00
            Hfxn = 1.214
            Sf = 5.834
         ELSEIF ( ATYPE.EQ.'Ba' ) THEN
            MASS = 137.33D+00
            Hfxn = 6.912
            Sf = 62.475
         ELSEIF ( ATYPE.EQ.'Be' ) THEN
            MASS = 9.01218D+00
            Hfxn = 1.932
            Sf = 9.44
         ELSEIF ( ATYPE.EQ.'Br' ) THEN
            MASS = 79.904D+00
            Hfxn = 24.509 / 2.0
            Sf = 152.206 / 2.0
         ELSEIF ( ATYPE.EQ.'Br79' ) THEN
            MASS = 78.918338D+00
            Hfxn = 24.509 / 2.0    ! assumed to be the same as "Br"
            Sf = 152.206 / 2.0     ! assumed to be the same as "Br"
         ELSEIF ( ATYPE.EQ.'Br81' ) THEN
            MASS = 80.916291D+00
            Hfxn = 24.509 / 2.0    ! assumed to be the same as "Br"
            Sf = 152.206 / 2.0     ! assumed to be the same as "Br"
         ELSEIF ( ATYPE.EQ.'C' ) THEN
            MASS = 12.011D+00
            Hfxn = 1.051
            Sf = 5.740
         ELSEIF ( ATYPE.EQ.'Ca' ) THEN
            MASS = 40.08D+00
            Hfxn = 5.736
            Sf = 41.588
         ELSEIF ( ATYPE.EQ.'C12' ) THEN
            MASS = 12.0D+00
            Hfxn = 1.051                ! assumed to be the same as "C"
            Sf = 5.740
         ELSEIF ( ATYPE.EQ.'C13' ) THEN
            MASS = 13.003355D+00
            Hfxn = 1.051                ! assumed to be the same as "C"
            Sf = 5.740
         ELSEIF ( ATYPE.EQ.'C14' ) THEN
            MASS = 14.003242D+00
            Hfxn = 1.051                ! assumed to be the same as "C"
            Sf = 5.740
         ELSEIF ( ATYPE.EQ.'C16' ) THEN
            MASS = 16.0147D+00
            Hfxn = 1.051                ! assumed to be the same as "C"
            Sf = 5.740
         ELSEIF ( ATYPE.EQ.'Cl' ) THEN
            MASS = 35.453D+00
            Hfxn = 9.181 / 2.0
            Sf = 223.079 / 2.0
         ELSEIF ( ATYPE.EQ.'Cl35' ) THEN
            MASS = 34.96885271D+00
            Hfxn = 9.17 / 2.0      ! Calculated using rot and vib constants from Bermejo et al., J. Mol. Spectrosc. 212, 186-193 (2002)
            Sf = 222.69 / 2.0      ! Calculated using rot and vib constants from Bermejo et al., J. Mol. Spectrosc. 212, 186-193 (2002)
         ELSEIF ( ATYPE.EQ.'Cl37' ) THEN
            MASS = 36.9659D+00
            Hfxn = 9.16 / 2.0      ! Calculated using rot and vib constants from Bermejo et al., J. Mol. Spectrosc. 212, 186-193 (2002)
            Sf = 223.84 / 2.0      ! Calculated using rot and vib constants from Bermejo et al., J. Mol. Spectrosc. 212, 186-193 (2002)
         ELSEIF ( ATYPE.EQ.'Cr' ) THEN
            MASS = 51.996D+00
            Hfxn = 4.057
            Sf = 23.618
         ELSEIF ( ATYPE.EQ.'Cs' ) THEN
            MASS = 132.9054D+00
            Hfxn = 7.717
            Sf = 85.147
         ELSEIF ( ATYPE.EQ.'Cu' ) THEN
            MASS = 63.546D+00
            Hfxn = 5.007
            Sf = 33.164
         ELSEIF ( ATYPE.EQ.'Fe' ) THEN
            MASS = 55.847D+00
            Hfxn = 4.507
            Sf = 27.321
         ELSEIF ( ATYPE.EQ.'H' ) THEN
            MASS = 1.00794D+00
            Hfxn = 8.467 / 2.0
            Sf = 130.680 / 2.0
         ELSEIF ( ATYPE.EQ.'H1' ) THEN
            MASS = 1.007825D+00
            Hfxn = 8.467 / 2.0          ! assumed to be the same as "H"
            Sf = 130.680 / 2.0
         ELSEIF ( ATYPE.EQ.'D' ) THEN
            MASS = 2.014102D+00
            Hfxn = 8.569 / 2.0
            Sf = 144.960 / 2.0
         ELSEIF ( ATYPE.EQ.'T' ) THEN
            MASS = 3.016049D+00
            Hfxn = 8.60 / 2.0             ! calculated by using parameters found in NIST web-book
            Sf = 153.24 / 2.0             ! calculated by using parameters found in NIST web-book
         ELSEIF ( ATYPE.EQ.'F' ) THEN
            MASS = 18.998403D+00
            Hfxn = 8.825 / 2.0
            Sf = 202.789 / 2.0
         ELSEIF ( ATYPE.EQ.'Ga' ) THEN
            MASS = 69.72D+00
            Hfxn = 5.561
            Sf = 40.838
         ELSEIF ( ATYPE.EQ.'He' ) THEN
            MASS = 4.00260D+00
            Hfxn = 6.197
            Sf = 126.152
         ELSEIF ( ATYPE.EQ.'Hg' ) THEN
            MASS = 200.59D+00
            Hfxn = 9.343
            Sf = 76.028
         ELSEIF ( ATYPE.EQ.'I' ) THEN
            MASS = 126.9045D+00
            Hfxn = 13.198 / 2.0
            Sf = 116.142 / 2.0
         ELSEIF ( ATYPE.EQ.'K' ) THEN
            MASS = 39.0983D+00
            Hfxn = 7.082
            Sf = 64.670
         ELSEIF ( ATYPE.EQ.'Kr' ) THEN
            MASS = 83.80D+00
            Hfxn = 6.197
            Sf = 164.084
         ELSEIF ( ATYPE.EQ.'Li' ) THEN
            MASS = 6.941D+00
            Hfxn = 4.622
            Sf = 29.085
         ELSEIF ( ATYPE.EQ.'Mg' ) THEN
            MASS = 24.3057D+00
            Hfxn = 4.998
            Sf = 32.671
         ELSEIF ( ATYPE.EQ.'Mn' ) THEN
            MASS = 54.9380D+00
            Hfxn = 4.994
            Sf = 32.010
         ELSEIF ( ATYPE.EQ.'Mo' ) THEN
            MASS = 95.94D+00
            Hfxn = 4.585
            Sf = 28.605
         ELSEIF ( ATYPE.EQ.'N' ) THEN
            MASS = 14.00674D+00
            Hfxn = 8.670 / 2.0
            Sf = 191.609 / 2.0
         ELSEIF ( ATYPE.EQ.'N14' ) THEN
            MASS = 14.003074D+00
            Hfxn = 8.670 / 2.0    ! Calculated using rot and vib constants from Gilson et al., J. Raman Spectrosc. 9, 361-368 (1980) & Bendtsen, J. Raman Spectrosc. 2, 133-145 (1974)
            Sf = 191.60 / 2.0     ! Calculated using rot and vib constants from Gilson et al., J. Raman Spectrosc. 9, 361-368 (1980) & Bendtsen, J. Raman Spectrosc. 2, 133-145 (1974)
         ELSEIF ( ATYPE.EQ.'N15' ) THEN
            MASS = 15.00010897D+00
            Hfxn = 8.670 / 2.0    ! Calculated using rot and vib constants from Gilson et al., J. Raman Spectrosc. 9, 361-368 (1980) & Bendtsen, J. Raman Spectrosc. 2, 133-145 (1974)
            Sf = 193.03 / 2.0     ! Calculated using rot and vib constants from Gilson et al., J. Raman Spectrosc. 9, 361-368 (1980) & Bendtsen, J. Raman Spectrosc. 2, 133-145 (1974)
         ELSEIF ( ATYPE.EQ.'Na' ) THEN
            MASS = 22.98977D+00
            Hfxn = 6.447
            Sf = 51.455
         ELSEIF ( ATYPE.EQ.'Nb' ) THEN
            MASS = 92.9064D+00
            Hfxn = 5.241
            Sf = 36.464
         ELSEIF ( ATYPE.EQ.'Ne' ) THEN
            MASS = 20.179D+00
            Hfxn = 6.197
            Sf = 146.327
         ELSEIF ( ATYPE.EQ.'Ni' ) THEN
            MASS = 58.69D+00
            Hfxn = 4.786
            Sf = 29.870
         ELSEIF ( ATYPE.EQ.'O' ) THEN
            MASS = 15.9994D+00
            Hfxn = 8.683 / 2.0
            Sf = 205.147 / 2.0
         ELSEIF ( ATYPE.EQ.'O16' ) THEN
            MASS = 15.994915D+00
            Hfxn = 8.68 / 2.0    ! Calculated using rot and vib constants from Edwards et al., J. Raman Spectrosc., 10, 60-63 (1981)
            Sf = 205.13 / 2.0    ! Calculated using rot and vib constants from Edwards et al., J. Raman Spectrosc., 10, 60-63 (1981)
         ELSEIF ( ATYPE.EQ.'O17' ) THEN
            MASS = 16.9991315D+00
            Hfxn = 8.68 / 2.0    ! Calculated using rot and vib constants from Edwards et al., J. Raman Spectrosc., 10, 60-63 (1981)
            Sf = 206.40 / 2.0    ! Calculated using rot and vib constants from Edwards et al., J. Raman Spectrosc., 10, 60-63 (1981)
         ELSEIF ( ATYPE.EQ.'O18' ) THEN
            MASS = 17.999160D+00
            Hfxn = 8.69 / 2.0    ! Calculated using rot and vib constants from Edwards et al., J. Raman Spectrosc., 10, 60-63 (1981)
            Sf = 207.60 / 2.0    ! Calculated using rot and vib constants from Edwards et al., J. Raman Spectrosc., 10, 60-63 (1981)
         ELSEIF ( ATYPE.EQ.'P' ) THEN
            MASS = 30.97376D+00
            Hfxn = 5.360
            Sf = 41.077
         ELSEIF ( ATYPE.EQ.'Pb' ) THEN
            MASS = 207.2D+00
            Hfxn = 6.878
            Sf = 64.785
         ELSEIF ( ATYPE.EQ.'Rb' ) THEN
            MASS = 85.4678D+00
            Hfxn = 7.490
            Sf = 76.778
         ELSEIF ( ATYPE.EQ.'S' ) THEN
            MASS = 32.06D+00
            Hfxn = 4.412
            Sf = 32.056
         ELSEIF ( ATYPE.EQ.'S32' ) THEN
            MASS = 31.9720707D+00
            Hfxn = 4.412             ! assumed to be the same as "S"
            Sf = 32.056              ! assumed to be the same as "S"
         ELSEIF ( ATYPE.EQ.'S34' ) THEN
            MASS = 33.9678668D+00
            Hfxn = 4.412             ! assumed to be the same as "S"
            Sf = 32.056              ! assumed to be the same as "S"
         ELSEIF ( ATYPE.EQ.'Se' ) THEN
            MASS = 78.96D+00
            Hfxn = 5.519
            Sf =  42.442
         ELSEIF ( ATYPE.EQ.'Si' ) THEN
            MASS = 28.0855D+00
            Hfxn = 3.218
            Sf = 18.820
         ELSEIF ( ATYPE.EQ.'Si29' ) THEN
            MASS = 28.98D+00
            Hfxn = 3.218              ! assumed to be the same as Si
            Sf = 18.820               ! assumed to be the same as Si
         ELSEIF ( ATYPE.EQ.'Si30' ) THEN
            MASS = 29.97D+00
            Hfxn = 3.218              ! assumed to be the same as Si
            Sf = 18.820               ! assumed to be the same as Si
         ELSEIF ( ATYPE.EQ.'Sn' ) THEN
            MASS = 118.69D+00
            Hfxn = 6.297
            Sf = 51.18
         ELSEIF ( ATYPE.EQ.'Sr' ) THEN
            MASS = 87.62D+00
            Hfxn = 6.568
            Sf = 55.694
         ELSEIF ( ATYPE.EQ.'Ta' ) THEN
            MASS = 180.9479D+00
            Hfxn = 5.681
            Sf =  41.471
         ELSEIF ( ATYPE.EQ.'Te' ) THEN
            MASS = 127.60D+00
            Hfxn = 6.121
            Sf =  49.71
         ELSEIF ( ATYPE.EQ.'Ti' ) THEN
            MASS = 47.88D+00
            Hfxn = 4.830
            Sf =  30.759
         ELSEIF ( ATYPE.EQ.'U' ) THEN
            MASS = 238.0290D+00
            Hfxn = 6.364
            Sf =  50.21
         ELSEIF ( ATYPE.EQ.'V' ) THEN
            MASS = 50.9415D+00
            Hfxn = 4.640
            Sf =  28.936
         ELSEIF ( ATYPE.EQ.'W' ) THEN
            MASS = 183.85D+00
            Hfxn = 4.973
            Sf =  32.660
         ELSEIF ( ATYPE.EQ.'Xe' ) THEN
            MASS = 131.29D+00
            Hfxn = 6.197
            Sf = 169.684
         ELSEIF ( ATYPE.EQ.'Zn' ) THEN
            MASS = 65.38D+00
            Hfxn = 5.669
            Sf = 41.717
         ELSEIF ( ATYPE.EQ.'Zr' ) THEN
            MASS = 91.22D+00
            Hfxn = 5.497
            Sf = 38.869
         ELSE
            STOP 'Element not recognized in subroutine element.f'
         ENDIF

            Hfxn = Hfxn/R
            Sfxn = Sf/R

      RETURN
      END
      
      
