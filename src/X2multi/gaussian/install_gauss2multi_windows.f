      character*850  pathnew, pathold
      call system ("CD > pathnew")
      call system ("echo %PATH% > pathold")
      open (unit=1,file='pathnew',status='OLD')
      open (unit=2,file='pathold',status='OLD')
      open (unit=3,file='gauss2multi.bat',status='unknown')
      read(1,'(A850)') pathnew
      read(2,'(A850)') pathold
      pathnew= "PATH "//pathnew(1:lnblnk(pathnew))//";"//
     &pathold(1:lnblnk(pathold))
      write(3,*) pathnew(1:lnblnk(pathnew))
      write(3,*) "gaussian2multi %1 "
      close(1)
      close(2)
      close(3)
      call system ("del pathnew")
      call system ("del pathold")
      end
       
