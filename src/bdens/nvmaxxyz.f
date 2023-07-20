
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c Copyright (C) 2009 John R. Barker
c
c Authors: Thanh Lam Nguyen and John R. Barker
c          nguyenlt@umich.edu   
c          Aug. 9, 2009
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

      INTEGER(4) FUNCTION nvmaxxyz( ns, w, x, y, z, nv, k, Eu )

c     Select highest quantum number to be sampled

      IMPLICIT NONE
      INTEGER(4) No
      PARAMETER (No=100)
c      INTEGER(4) ns, nv(No), i, j, k, nvd
      INTEGER(4) ns, nv(No), k
c      REAL(8) w(No), Eu, D , vd
      REAL(8) w(No), Eu, D
      REAL(8) x(No,No), y(No,No,No), z(No,No,No,No)
c      REAL(8) energy, ENER, EGS, energyxyz
      REAL(8) ENER, EGS, energyxyz
      REAL(8) Fo

      IF(Eu.GT.0.0d0) THEN
            D=0.0d0 
              nv(k)=0
              EGS=energyxyz( ns, w, x, y, z, nv )
11            nv(k)=nv(k)+1
              ENER=energyxyz( ns, w, x, y, z, nv ) - EGS
            IF(D.LT.ENER) THEN
                  D=ENER
                  IF(D.LE.Eu) THEN
                                   goto 11
                  ELSE
                        nv(k)=nv(k)-1
                        nvmaxxyz=nv(k)
                      ENDIF
            ELSE
                  nv(k)=nv(k)-1
                  CALL CALDERIXYZ(Fo,nv(k),ns,k,nv,w,x,y,z)
                  IF(Fo.GE.0.0d0) THEN
                        nvmaxxyz=nv(k)
                  ELSE
                        nvmaxxyz=nv(k)-1
                  ENDIF            
            ENDIF
      ELSEIF(Eu.EQ.0.0d0) THEN
            nvmaxxyz=0
      ELSE
            nvmaxxyz=-1
      ENDIF

      nv(k)=0

      RETURN
      END            
C==================================


