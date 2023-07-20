c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!    
!    LICENSE NOTICE
!
!    bdens: sums and densities of coupled anharmonic vibrations
!    Copyright (C) 2015 Collin Li, Thanh Lam Nguyen, and John R. Barker
!
!    This program is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    This program is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License: <http://www.gnu.org/licenses/>.
!
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c                                                                 
c    PROGRAM bdens
c    by John R. Barker, Thanh Lam Nguyen, and Collin G. L. Li    
c                                                                  
c  ***Direct Count and Wang-Landau Algorithm for Densities of States***
c
c    Contact:
c    John R. Barker   (email: jrbarker@umich.edu)
c    Department of Atmospheric, Oceanic, and Space Sciences
c    College of Engineering
c    University of Michigan
c    Ann Arbor, MI 48109
c
c    http://aoss.engin.umich.edu/multiwell
c
c@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
c
c     calein: a subroutine that energy levels using recursive loops method.
c	      It is used to replace calei1 to calei9 in doloops (which becomes newloops)
c	      A look at calei1.f and calei9.f should make this code easier to understand.
c
c     E = max energy
c     freq = frequency array
c     n = number of wave number
c
      SUBROUTINE CALEVN( E, freq, n, x, y, z, NY, NZ )

      Integer No
      Parameter (No=100)

      Integer n, NY, NZ, nvmaxdoloops, nvmaxxyzdoloops      
      Integer NN(No),Maxi(No)
      Real(8) freq(No), x(No,No), y(No,No,No), z(No,No,No,No)
      Integer lam, ntest
      Integer save Max
      Real(8) E, tp, Eu
      Integer c, Val
            
      c=1
      lam=0

      DO I=1, n
        NN(I)=0
      ENDDO
      Eu=E
      
      IF((NY.EQ.0).AND.(NZ.EQ.0)) THEN
        Maxi(1)= nvmaxdoloops( n, freq, x, NN, 1, Eu )
      ELSE
        Maxi(1)= nvmaxxyzdoloops( n, freq, x, y, z, NN, 1, Eu)
      ENDIF

      Max = Maxi(1)  
      DO WHILE(Maxi(1).GE.0)
        Val = Max-Maxi(1)
        CALL CALEVNSUB(c,Val,NY,NZ,E,freq,n,NN,x,y,z,Maxi,ntest,lam)
        Maxi(1) = Maxi(1) - 1
      ENDDO

      RETURN 
      END

      RECURSIVE SUBROUTINE CALEVNSUB(c,Val,NY,NZ,E,freq,n,NN,x,y,z,Maxi,
     &          ntest,lam)
      Parameter (No=100)
      Integer NY, NZ, n, Val, ntest, lam
      Integer save lic
      Integer save Max
      Integer c
      Real(8) E,Eu, freq(No),x(No,No),y(No,No,No)
      Real(8) z(No,No,No,No)
      Integer NN(No), Maxi(No)
      Real(8) tp,energydoloops,energyxyzdoloops
      
      IF(c.GE.n) THEN
         CALL CALEVNSUB2(c,Val,NN,NY,NZ,freq,n,x,y,z,ntest,lam,E)
      ELSE
        lic = c
        c = c+1
        DO I=lic, n
          NN(I)=0
        ENDDO
        NN(lic)=Val
        IF((NY.EQ.0).AND.(NZ.EQ.0)) THEN
          Eu = E - energydoloops(freq,n,NN,x)
          Maxi(lic+1) = nvmaxdoloops(n,freq,x,NN,lic+1,Eu)
        ELSE
          Eu = E - energyxyzdoloops(freq,n,NN,x,y,z)
          Maxi(lic+1) = nvmaxxyzdoloops(n,freq,x,y,z,NN,lic+1,Eu)
        ENDIF
      
        Max = Maxi(lic+1)

        DO WHILE(Maxi(lic+1).GE.0)
        Val = Max-Maxi(lic+1)
        CALL CALEVNSUB(c,Val,NY,NZ,E,freq,n,NN,x,y,z,Maxi,ntest,lam)
        Maxi(lic+1) = Maxi(lic+1) - 1
        ENDDO
        c=c-1
      ENDIF
      END

      SUBROUTINE CALEVNSUB2(c,Val,NN,NY,NZ,freq,n,x,y,z,ntest,lam,E)
      Parameter (No=100)
      Integer c,n,Val,NN(No),NY,NZ,ntest,lam
      Real(8) freq(No), x(No,No),y(No,No,No)
      Real(8) z(No,No,No,No)
      Real(8) tp,energydoloops,energyxyzdoloops,E
      REAL(8) Egrain1
      INTEGER nTmax, nT, nk
      COMMON/TBAR/ nTmax, nT(20003), Egrain1

      NN(c)=Val

      IF((NY.EQ.0).AND.(NZ.EQ.0)) THEN
        tp = energydoloops(freq,n,NN,x)
      ELSE
        tp = energyxyzdoloops(freq,n,NN,x,y,z)
      ENDIF

      IF((tp.LE.E).AND.(tp.GE.0.0d0)) THEN
        IF((NY.EQ.0).AND.(NZ.EQ.0)) THEN
                CALL ckderiv(n,NN,freq,x,ntest)
        ELSE
                CALL ckderivxyz(n,NN,freq,x,y,z,ntest)
        ENDIF
        IF(ntest.EQ.0) THEN
                lam=lam+1
                nk = CEILING(tp/Egrain1) + 1
                IF ( nk .LE. nTmax) nT(nk) = nT(nk) + 1
        ENDIF
      ENDIF
      END

