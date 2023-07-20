      call system ("echo PATH=$PWD/../../bin:$PATH  > gauss2multi")
      open (unit=1,file='gauss2multi',status='OLD')
      read(1,*) 
      write(1,*) "gaussian2multi $1"
      close(1)
      call system ("chmod a+x gauss2multi")
      end
       
