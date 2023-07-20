      subroutine extractDataMomInert(nameG03,countAtom,Ix,Iy,Iz)
      double precision Ix,Iy,Iz
      integer nul,countAtom
      character*52 Find
      character*72 nameG03
      countAtom=0
       
      Find=nameG03(1:lnblnk(nameG03)-4)//".coords.out" 
      open(unit=3,file=Find(1:lnblnk(Find)),status='old')

      do
      read(3,'(A7)') Find
      if (Find.eq." NUMBER") then
       do
        read(3,'(I8)') nul
          if (nul.eq.0) exit
        countAtom=countAtom+1
       enddo
      endif

100   if (Find.eq."  REDUC") write(*,*) " WARNING: Reduced moment of
     & inertia will NOT be automatically added to .vibs file!"
      if (Find.eq."PRINCIP") exit
      enddo

      read(3,*) Find,Find,Ix,Find,Find,Iy,Find,Find,Iz
      close(3)

      end subroutine
        
