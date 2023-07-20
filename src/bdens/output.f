       subroutine otput(sum,imax,range,den)
       character*80 output    !the name of output file
       integer sum(20000)      !the total number of states
       integer imax           !the number of energyGrains
       real*8 range(20000)    !the start number of each energy grain
       integer den(20000)      !the number of states in an energy grain

       output = 'result.txt'
       open(unit=1111,FILE=output)
       write(1111,*)'The result is: '
       do i = 1,imax+1
          write(1111,*)range(i),den(i),sum(i)
       end do

       end subroutine otput
