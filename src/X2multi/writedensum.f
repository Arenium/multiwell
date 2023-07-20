      subroutine writedensum(nameG03,Find)
      integer nul,i, countAtom, Multiplicity,DoF,ESN
      double precision Ix,Iy,Iz,Krot,ADrot,Ezpe
      double precision Freq1,Freq2,Freq3,Freq(900)
      double precision MolecularMass
      character*4  rottype
      character*50 Stoichiometry
      character*72 Find
      character*72 FindFreq,nameG03

! read gaussian log (nameG03) : find cartesian coordinates, read X,Y,Z and
! write  Momintert.dat
! LIMIT: ca 300 atoms

      call extractDataMomInert(nameG03,countAtom,Ix,Iy,Iz)
      call extractDataOUT(nameG03,countAtom,Multiplicity,Stoichiometry,
     & Ezpe,MolecularMass,ESN,Freq,DoF)
      call calcRotor(Ix,Iy,Iz,Krot,ADrot)
      
!     open(unit=1,err=1000,file='gauss2multi.cfg',status='old')
!     do i=1,6
!     read(1,*)
!     enddo
!     read Egrain, Imax, isize from gauss2multi.cfg
!     read(1,'(A72)') Find 
!     close(1)

!write densum.dat
100   open(unit=2,file='densum.dat',status='unknown')
      write(2,*) "Title"
      write(2,*) nameG03(1:(lnblnk(nameG03)-4))
      if(Krot.eq.0) then
       write(2,'(I4,A16)')   DoF, "  0   HAR   AMUA"
      else
       write(2,'(I4,A16)') DoF+1, "  0   HAR   AMUA"
      endif
      write(2,*) Find
      call writeFreq(DoF,Freq,Krot,ADrot)
      close(2)
      

      open(unit=2,access='append',file='energies',status='unknown')
      write(2,*) nameG03(1:(lnblnk(nameG03)-4)), "  ", Ezpe
      close(2)
      return
1000  write(*,*) "File gauss2multi.cfg not found"
      write(*,*) "Default Egrain, Max1, Isize, Emax2 will be used in den
     &sum"
      Find="10    400     500     50000"
      goto 100
      end subroutine


