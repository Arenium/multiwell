c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c Copyright (C) 2009 Thanh Lam Nguyen
c
c Authors: Thanh Lam Nguyen 
c          nguyenlt@umich.edu
c          Dec. 23, 2009
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


C
C To write out the results
C


      SUBROUTINE WRITE_OUT(G44,ANG,M,C1,C2,C3,DELA,E)

      Double precision C, DELA, CC
      Parameter (C=16.85763D0)
      Integer MMAX, M, I, MM
      Parameter (MMAX=500)
      Double precision G44(MMAX), ANG(MMAX), E(MMAX)

      CHARACTER(len=100) C1 , C2 , C3 

      CHARACTER(len=6) AVERSION
      CHARACTER(len=8) ADATE 
      PARAMETER ( AVERSION='2022', ADATE='Jan 2022' )

      CC = acos(-1.0d0)/180.0d0

      OPEN(UNIT=4,STATUS='UNKNOWN',FILE='lamm.out')

      WRITE(4,9001) AVERSION , ADATE , AVERSION , ADATE

      write(4,*)
        
      WRITE(4,99) C1
        WRITE(4,99) C2

      write(4,*)
        WRITE(4,99) C3
        write(4,*)

      write(4,9999) "  INDEX    ANG(DEG)    ANG(RAD)      E(cm-1) 
     &I(amu.A**2)   B(cm-1)"

      IF(ANG(M).eq.(360.0d0-DELA)) THEN
            MM=M+1
            G44(MM)=G44(1)
            G44(MM-1)=G44(2)
            ANG(MM)=360.0d0
            E(MM)=E(1)
      ELSE
            MM=M
      ENDIF

      DO I=1, MM
        write(4,999) I, ANG(I), ANG(I)*CC, E(I), 1.0d0/G44(I), C*G44(I)
      ENDDO

        CLOSE (unit=4, status='keep')

99      FORMAT(A100)
999      FORMAT(1X,I5,7X,F5.1,7X,F8.6,5X,F7.1,4X,F8.4,6X,F8.4)   ! modified March 26, 2014 by JRB
9999      FORMAT(A73)

9001  FORMAT (
     &'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     &%%%%%%%%%%%%%'//,
     &'                      LAMM-',A6,//
     &'                       ',A8,//,
     &'                   Thanh Lam Nguyen',//
     &'                 Contact: jrbarker@umich.edu',/
     &'                   University of Michigan',/
     &'               Ann Arbor, Michigan 48109-2143',//
     &'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     &%%%%%%%%%%%%%',//,
     &7x,'J.R. Barker, T.L. Nguyen, J.F. Stanton, C. Aieta, M. Ceotto,',
     &/7x,'F. Gabas, T.J.D. Kumar, C.G.L.Li, L.L.Lohr, A. Maranzana,',
     &/7x,'N.F. Ortiz, J.M. Preses, J.M. Simmie, J.A. Sonk, and ',
     &'P.J. Stimac; ',
     &/7x,'MultiWell-',A7,' Software Suite; University of Michigan, ',
     &/7x,'Ann Arbor, Michigan, USA, ',A8,'; ',
     &/7x,'http://clasp-research.engin.umich.edu/multiwell/.', //,
     &'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     &%%%%%%%%%%%%%'/)

      RETURN
      END
