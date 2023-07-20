      subroutine calcRotor(Ix,Iy,Iz,Krot,ADrot)
      Double precision fxy,fxz,fyz
      Double precision Ix,Iy,Iz,Krot,ADrot
      
      ! identify K-rot  and calculate 2D- rotor 
      fxy=ABS(Ix-Iy)/(Ix+Iy)
      fxz=ABS(Ix-Iz)/(Ix+Iz)
      fyz=ABS(Iy-Iz)/(Iy+Iz)


      IF(fxy.lt.fxz.AND.fxy.lt.fyz) then
        Krot=Iz
        ADrot=sqrt(Ix*Iy)
      ELSEIF(fxz.lt.fxy.AND.fxz.lt.fyz) then
        Krot=Iy
        ADrot=sqrt(Ix*Iz)
      ELSEIF(fyz.lt.fxy.AND.fyz.lt.fxz) then
        Krot=Ix
        ADrot=sqrt(Iy*Iz)
      ELSE
        Krot=Ix                 ! Symmetryc top
        ADrot=sqrt(Iy*Iz)
      ENDIF
       
      end subroutine
