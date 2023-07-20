 
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c Thermo: a code for thermochemical calculations.
c Copyright (C) 2001 through 2022 John R. Barker
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
c
c
c      Equilibrium constants and thermodynamic parameters
c      based on molecular properties
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      PROGRAM THERMO

      include 'declare.inc' 

      REAL(8) rmax0 , I2dGorin0 , Veff0
      INTEGER(4) ix, mrso
      SAVE

      OPEN (UNIT=21,STATUS='UNKNOWN',FILE='thermo.partfxns') 
      OPEN (UNIT=22,STATUS='UNKNOWN',FILE='thermo.details') 

      call read_dat                          ! READ THERMO.DAT AND CALCULATE ZPE

      mrso = 0
      DO i = 1, Ns
       mrso = mrso + nso(i)
      END DO
      IF ( mrso .GT. 0 ) THEN
          OPEN (UNIT=23,STATUS='UNKNOWN',FILE='thermo.rso') 
      ENDIF

      
      T(Nt+1) = 298.15d+00                   ! always run at 298.15 K, but don't print this one
      

      IF ( Eunits.EQ.'KCAL' ) THEN
         WRITE (21,99004) UNITS 
         WRITE (22,99004) UNITS 
      ELSEIF ( Eunits.EQ.'KJOU' ) THEN
         WRITE (21,99005) UNITS 
         WRITE (22,99005) UNITS 
      ELSEIF ( Eunits.EQ.'CM-1' ) THEN
         WRITE (21,99055) UNITS 
         WRITE (22,99055) UNITS 
      ELSE
         WRITE(*,*) ' '
         WRITE(*,*) '*** FATAL: Enthalpy units not recognized ***'
         WRITE(*,*) ' '
         STOP
      ENDIF

      WRITE (21,9902)  
      WRITE (22,9905)  

      DO Temp=1,Nt+1  ! Nt  number of temperatures
       
        if(Gorin.gt.0 ) then 
         call MaxVeff(beta,POT,c,De,re,mu,T(Temp),rmax0,I2dGorin0,V0,
     &                Veff0,Eunits)                             ! I2d Gorin Subroutine

         V(Temp)=V0
         Veff(temp)=Veff0
         I2dGorin(Temp)=I2dGorin0
         bet(Temp)=beta
         rmax(Temp)=rmax0
        endif

        if(hindrance.gt.0 .AND. Temp.LE.Nt) then 
           call hindrance_fit                         !  Hindrance calculation for Gorin fit

        ELSE

           call thermochem                            !  Thermochemistry calculation at temperature Temp

        ENDIF

           T_AK(Temp)= AK(2)
           T_AA(Temp)= AA
           T_BB(Temp)=BB
           T_DSR(Temp)=DSR
           T_DHR(Temp)=DHR
           T_DCp(Temp)=DCp
           T_Sfac(Temp)=Sfac

           do ix=1,Ns
            T_S(Temp,ix)=S(ix)
            T_Cp(Temp,ix)=Cp(ix)
            T_H(Temp,ix)=H(ix)
            IF ( Temp .EQ. Nt+1 ) THEN
              delH298(ix) = delH(ix) + (H(ix) - Hfxn298(ix))/DIV      ! delH at 298.15 K; calculate enthalpy fxn for elements at 298.15 K
              delG298(ix) = delH298(ix) - 0.29815*(S(ix)-Sfel(ix))    ! delG at 298.15 K (in kcal or kJ)              
            ENDIF
           enddo
      ENDDO

      call write_out                                  ! Write thermo.out 
      
!   Close output files
      CLOSE(21)
      CLOSE(22)
      CLOSE(23)
      CLOSE(24)

!   FORMAT STATEMENTS
  
9902  FORMAT (' PARTITION FUNCTIONS (for some categories)',/, 
     &    3x, 'DISPLAYED: electronic, diatomic spin-orbit+rotation, ',
     &        'translation, vibration, rotation, hindered rotation.',/
     &    3x, 'NOTE: rotation includes free rotation types qro, rot, ',
     &         'and top.',/,
     &    3x, 'NOTE: transition state CRPs, particle in a box,'
     &      ' and some other types are not included in this table.',//,
     &     11x,'Species',30x,'T(K)',22x,'Q(elec)',19x,'Q(trans)',
     &     18x,'Q(vib)',20x,'Q(rot)',20x,'Q(hindro)',17x,'Q(rot+s.o.)' )
9905  FORMAT ( 
     &  'HEAT CAPACITY (Cp), ENTROPY (S), and ENTHALPY FUNCTION (H)',/,
     &    3x, 'DISPLAYED: symm.opt, electronic, translation, '
     &        'vibration, rotation, hindered rotation, '
     &        'rot+Spin.Orb., ',/
     &    3x, 'NOTE: symm.opt includes external symmetry number and '
     &            'number of optical isomers.',/,
     &    3x, 'NOTE: rotation includes free rotation types qro, rot, ',
     &            'and top.',/,
     &    3x, 'NOTE: transition state CRPs, particle in a box,'
     &      ' and some other types are not included in this table.',//,
     &      '     Species',31x,'T(K)  S(symm.opt)  ',
     &                      'Cp(elec)    S(elec)    H(elec)  ',
     &                      'Cp(trans)   S(trans)   H(trans)    ',
     &                      'Cp(vib)     S(vib)     H(vib)     ',
     &                      'Cp(rot)    S(rot)     H(rot) ',
     &                      'Cp(hindro)  S(hindro)  H(hindro)   ',
     &                      'Cp(r+s.o.)  S(r+s.o.)  H(r+s.o.)' )
99004 FORMAT (//'*  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  ',
     &'*  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *'//, 
     &'THERMODYNAMIC QUANTITIES',10x,'Std. State = ',A12,/,3x,
     &        'Entropy & Cp units: cal/K/mole'/3x,
     &        'DelS         units: cal/K/mole'/3x,
     &        '[H(T)-H(0)]  units: kcal/mole'/3x,
     &        'DelH & DelG  units: kcal/mole'/)
99005 FORMAT (//'*  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  ',
     &'*  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *'//, 
     &'THERMODYNAMIC QUANTITIES',10x,'Std. State = ',A12,/,3x,
     &        'Entropy & Cp units: J/K/mole'/3x,
     &        'DelS         units: J/K/mole'/3x,
     &        '[H(T)-H(0)]  units: kJ/mole'/3x,
     &        'DelH & DelG  units: kJ/mole'/)
99055 FORMAT (//'*  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  ',
     &'*  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *'//, 
     &'THERMODYNAMIC QUANTITIES',10x,'Std. State = ',A12,/,3x,
     &        'Entropy & Cp units: CM-1/K'/3x,
     &        'DelS         units: CM-1/K'/3x,
     &        '[H(T)-H(0)]  units: CM-1'/3x,
     &        'DelH & DelG  units: CM-1'/)
 
      END


