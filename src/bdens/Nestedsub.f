       subroutine NestedLoops(Emax,grain,s,w,limit,v,imax,sum,range,den)
       real*8 Emax            !upper limit of energy
       real*8 grain           !the length of energy grains
       integer s              !the number of frequencies
       real*8 w(20)           !the array store the frequencies
       integer limit(20)      !Emax / frequencies
       integer v(20)          !hold the number for calculation
       integer imax           !the number of energyGrains
       real*8 range(20000)      !the start number of each energy grain
       integer den(20000)       !the number of states in an energy grain
       real*8 product
       integer temp           !store the value of freqN for recursion
       integer sum            !the total number of states

       product = 0.0
       temp = s

       do i = 1,s
          limit(i) = int(Emax / w(i))
          v(i) = 0
       end do

       imax = int(Emax / grain)
       
       do i = 1,imax+1
          range(i) = (i-1)*grain
       end do 
       
       do i = 1,imax
          den(i) = 0
       end do

       call NL(w,v,Emax,s,product,limit,temp,sum,range,den,imax)

       end subroutine NestedLoops

       recursive subroutine NL(w,v,Emax,s,p,l,t,sum,range,den,imax)
       real*8 w(20)
       integer v(20)
       real*8 Emax
       integer s,l(20)    !l is limit array
       real*8 p               !p is product,need to compare with Emax
       integer t
       integer sum      !t store the value of freqN,r is result
       integer save lic
       real*8 range(20000)
       integer den(20000)
       integer imax

       !do i = 1,imax
       !   write(*,*)'in NL',den(i)
       !end do
       !write(*,*)'in NL',t
       p = 0
       if (t.EQ.1) then
           call NLH(w,v,Emax,p,l,s,sum,range,den,imax)
       else
           t = t - 1
           lic = t
          do while((v(lic+1).LE.l(lic+1)).AND.(p.LE.Emax))
             v(lic) = 0 
             t = lic
             !write(*,*)'in Loops',numberarr(1),numberarr(2),numberarr(3)
             product = 0
             call NL(w,v,Emax,s,p,l,t,sum,range,den,imax)
             v(lic+1) = v(lic+1) + 1             
          end do
       end if

       end subroutine NL

       subroutine NLH(w,v,Emax,product,l,s,sum,range,den,imax)
       real*8 w(20)
       integer v(20)
       real*8 Emax
       integer l(20),s,sum
       real*8 product
       real*8 range(20000)
       integer den(20000)
       integer imax
       
       !do i = 1,imax
       !   write(*,*)'in NLH',den(i)
       !end do
       !write(*,*)'in NLH',l(1),l(2),l(3)
       !write(*,*)'in NLH',numberarr(1),numberarr(2),numberarr(3)
       do while(v(1).LE.l(1))
          do i = 1,s
             product = product + w(i)*v(i)
          end do
          !write(*,*)'product',product
          !write(*,*)'Emax',Emax
          if (product.LE.Emax)then
             sum = sum + 1
             do i = 1,imax
                if(product.GE.range(i).AND.product.LT.range(i+1))then
                   den(i) = den(i) + 1
                end if
             end do 
          end if
          !write(*,*)'result',result
          v(1) = v(1) + 1
          product = 0
          if (product.GT.Emax)then   !if is not working here
             exit
          end if
       end do

       end subroutine NLH
