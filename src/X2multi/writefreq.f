      subroutine writeFreq(Dof,Freq,Krot,ADrot)
      integer Dof
      double precision Freq(900), Krot,ADrot
      character rottype*4
      
      do i=1, Dof
       write(2,'(I3,A9,F12.4,A10)') i, "   vib   ",Freq(i),"  0.0    1"
      enddo

      if(Krot.lt.11) then
       rottype="qrot "
      else
       rottype="rot  "
      endif

      if(Krot.ne.0) then
      write(2,'(I3,A2,A5,F10.4,A29)') i, "  ",rottype, Krot ,
     & "      1.0    1      ! K-rotor"
      else
      i=i-1
      endif

      if(ADrot.lt.11) then
       rottype="qrot "
      else
       rottype="rot  "
      endif
      if(ADrot.ne.0) then
        write(2,'(I3,A2,A5,F10.4,A40)') i+1, "  ",rottype, ADrot,
     & "     1.0    2      ! 2D adiabatic rotor"
        write(2,*)  
      endif
      
      end subroutine
      
