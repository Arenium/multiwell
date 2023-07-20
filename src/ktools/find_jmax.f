c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c
c copyright (c) 2017, 2020 john r. barker, jason a. sonk
c
c john r. barker
c jrbarker@umich.edu
c university of michigan
c ann arbor, mi 48109-2143
c (734) 763 6239
c
c this program is free software; you can redistribute it and/or
c modify it under the terms of the gnu general public license (version 2)
c as published by the free software foundation.
c
c this program is distributed in the hope that it will be useful,
c but without any warranty; without even the implied warranty of
c merchantability or fitness for a particular purpose. see the
c gnu general public license for more details.
c
c see the 'readme' file for a copy of the gnu general public license,
c or contact:
c
c free software foundation, inc.
c 59 temple place - suite 330
c boston, ma 02111-1307, usa.
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c	The find_jmax subroutine does 4 tests to try and find a maximum j-value.
c	Because we don't know how many trial transition states a user will pass, 
c	until they run the code, we check 
c	1) j corresponding to 20kbT
c	2) j corresponding to when the first trial transition state is higher in
c		energy than the second trial transition state (but not necessarily 
c		the reactants/products)
c	3) j corresponding to when the first trial transition state is higher in 
c		energy than all the other trial transition states (but not 
c		necessarily the reactants/products)
c	4) j corresponding to a totally repulsive effective potential
c
c	it stores them as "testj#" and then does a couple of checks to see which 
c	one to use when appropriate (depending on how many structures are passed 
c	to the code for example) and to make sure we're within the jtop value set 
c	by the program as the highest possible J.  Then it finally assigns "maxj"
c	as the maximum-j value the code will go up to for convergence.  
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine find_jmax(lunit,etopuser,jtopuser,jtop,rcnt,ntts,pcnt,
     $maxtemp,vmax,delh,epsil,barr,delhf,maxj,maxemax,testj3)

      implicit none
      real*8  v1,v2,b1,b2,jlow,temp
      real*8  test(rcnt+ntts+1),test2(rcnt+ntts+1)
      real*8  emax,minjmax,maxtest1,maxtest2
      real*8  maxemax							! maximum energy for allocating arrays
      real*8  testarr(4)
      integer lunit,etopuser,jtopuser,jtop,rcnt,ntts,pcnt,maxj
      integer i,j
      real*8  kb,maxtemp,vmax,delh(rcnt+ntts+1)
      real*8  epsil(rcnt+ntts+1),barr(rcnt+ntts+1)
      real*8  delhf(rcnt+ntts+1)
      parameter (kb=0.69503476d0)
      integer indarr(4)
      integer testj1,testj2,testj3,testj4,jlist(rcnt+ntts+1)
      character str*27,ans
      logical check
 
 100  format(a,2x,f10.2)
 101  format(2x,a10,2(4x,a10))
 102  format(2x,i10,2(4x,i20))
 112  format(a17,2x,i10)


      check=.false.
c      write(*,*)rcnt,ntts,pcnt
      do i=1,rcnt+ntts+1
c         write(*,*)delh(i),epsil(i),barr(i),delhf(i)
         test(i)=0.0d0
         test2(i)=0.0d0
      end do

      call sestamp('find_jmax',1)

c write individual jmax values for TTSs

      vmax=delh(rcnt+1)
      do i=rcnt+1,rcnt+ntts					! vmax is highest TTS pot energy relative to reactant
         if(delh(i).ge.vmax)vmax=delh(i)
      end do

      write(lunit,*)'For Emax: ',etopuser
      write(lunit,*)'following are jmax for each species'
      write(lunit,*)'designated "jtop"'
      write(lunit,*)
      write(lunit,*)
      write(lunit,101)'i','jtop'
      do i=1,rcnt+ntts+1					! epsil for reactant and all TTSs
         if(vmax-delh(i).lt.0)then
            epsil(i)=etopuser
         else
            epsil(i)=etopuser-delh(i)
         end if
         jlist(i)=nint(dsqrt(epsil(i)/barr(i)))
         write(lunit,102)i,jlist(i)					! jlist stores list of jtop
      end do

      write(lunit,*)
      write(lunit,*)

c  write table of kb*t values
      write(lunit,100)'for tmax of:',maxtemp
      write(lunit,*)'and increasing n multiples of kb*tmax'
      write(lunit,*)
      write(lunit,101)'n*kb*tmax','jmax'
      do i=1,20
         testj1=nint(dsqrt(((i*kb*maxtemp)/barr(rcnt+ntts))))
         write(lunit,102)i,testj1
      end do
      write(lunit,*)

c  test #1
c     find 20*kb*t=bin j*(j+1) --> 20kbt~=bj**2 --> j~=sqrt(20kbt/b)
c     using bmin as most reactant like trial ts (should be lowest 2d rotational constant)

      jlow=barr(rcnt+1)
      do i=1,ntts
         if(barr(rcnt+ntts).lt.jlow)jlow=barr(rcnt+ntts)
      end do

      testj1=nint(dsqrt(((20.0d0*kb*maxtemp)/jlow)))

c  test #2
c     if there are two or more points
c     and the potential is increasing find when
c     the first TS is higher than the second TS

      if((rcnt+ntts.ge.2).and.(delhf(rcnt+1).gt.0.0d0))then
        v1=delhf(rcnt+1)
        b1=barr(rcnt+1)

        v2=delh(rcnt+2)
        b2=barr(rcnt+2)

        if(b1.eq.b2)then
c           write(*,*)' '//achar(27)//'[31m '
c           write(*,*)'Rotational constants of the first'
c           write(*,*)'two Trial Transition States are identical.'
c           write(*,*)'This can lead to odd values of Jmax.'
c           write(*,*)'Continue? [Y/N]'//achar(27)//'[0m'
c           read(*,*)ans
c           call lowerc(ans)
c           if(ans.eq.'n')stop
c           b2=b2+1.0d-13
           b2=b2*(1.0e+00 + 1.0d-09)
c           check=.true.
        end if

        temp=(v2-v1)/(b1-b2)
        testj2=nint(dsqrt(dabs(temp)))

      else
      
        testj2=0

      end if

c  test #3
c     if there are more than 2 points find j where veff is entirely repulsive (first deriviative is negative at all points)

c  test #4
c     look for the first j where the first ts is higher in energy than everything else, there might still be well present


      if(ntts.ge.3)then
         if(testj1.le.testj2)then                                       ! start from lower of 2 previous jmax searches
            maxj=testj1
         else
            maxj=testj2
         end if
c         write(*,*)'maxj,t1,t2',maxj,testj1,testj2
         do j=maxj, 10*maxj                                             ! look over 10x previous jmax value
            do i=rcnt+1,rcnt+ntts                                       ! generate veff along TS potential
               test(i)=delh(i)+barr(i)*j*(j+1)
            end do

            do i=rcnt+1,rcnt+ntts                                       ! Midpoint forward differences derivative
               if(i.eq.rcnt+1)then
                  test2(i)=test(i+1)-test(i)
               elseif(i.eq.rcnt+ntts)then
                  test2(i)=test(i)-test(i-1)
               else
                  test2(i)=0.5d0*(test(i+1)-test(i-1))
               end if
            end do
            
            maxtest1=test(1)
            maxtest2=test2(1)
            do i=2,rcnt+ntts+1
               if(test(i).gt.maxtest1)maxtest1=test(i)
               if(test2(i).gt.maxtest2)maxtest2=test2(i)
            end do
            
c            write(*,*)j,maxtest1,maxtest2 

            if((maxtest1.eq.test(1)).and.(testj4.eq.0))then         ! Look for first j where ts1 > ts-rest
               testj4=j
            else
               testj4=j
            end if
      
            if(maxtest2.le.0)exit

         end do

         testj3=j

      else
        
         testj3=0

      end if

c  fixing tests depending on number of points on potential

      if(ntts.le.1)then                                                ! if there is 0 or 1 ts present then use maximum found jmax
         if(testj1.ge.testj2)then
            maxj=testj1
         else
            maxj=testj2
         end if
      else                                                              ! else use the minimum found jmax
c         if(testj1.le.testj2)then                           
c            if(testj1.ge.testj3)then
c               maxj=testj1
c            else
c               maxj=testj3
c            end if
c         else
c            if(testj2.ge.testj3)then
c               maxj=testj2
c            else
c              maxj=testj3
c            end if
c         end if
         maxj=testj4
      end if


      if(ntts.le.1)maxj=jtopuser
      if(jtopuser.ge.testj3)then
         write(lunit,*)'user j is higher than totally repulsive j'
      endif
      if((rcnt+ntts.gt.2).and.(maxj.gt.testj4))maxj=testj4
      if(maxj.gt.jtop-1)maxj=jtop-1
      if(maxj.gt.testj1)maxj=testj1
      testarr(1)=testj1*1.0d0
      testarr(2)=testj2*1.0d0
      testarr(3)=testj3*1.0d0
      testarr(4)=testj4*1.0d0

      call indexx(4,testarr,indarr)


c check to make sure j is not set to 0      

c      do i=1,4
c         if(testarr(indarr(i)).ne.0.0d0)then
c            maxj=int(testarr(indarr(i)))
c            exit
c         end if
c      end do



c      maxj=2                                                           ! to debug

c log print

      write(lunit,*)
      write(lunit,112)'looking for jmax'
      write(lunit,112)'20kbt method:'     ,testj1
      write(lunit,112)'pt1 > pt 2 method:',testj2
      write(lunit,112)'pt1 > everything:' ,testj4
      write(lunit,112)'tot.repul.method:' ,testj3
      str=adjustl('setting maxj to:')
      if(maxj.gt.jtop)then
         write(lunit,*)'maxj is larger than jtop'
         maxj=jtop
      elseif(ntts.eq.1)then
         write(lunit,*)'only one ts found, minimum jtop'
          maxj=jlist(1)
          do i=1,rcnt+ntts+1
            if(jlist(i).le.maxj)maxj=jlist(i)
          end do
         write(lunit,*)'only one ts found,using jtopuser'
         maxj=jtopuser
      elseif(check)then
         write(lunit,*)'identical rot.constants: using 20kbt method'
         maxj=testj1
      end if

c994   FORMAT (//,'*** FATAL: USER should increase Jmax to ≥',I4,
c     &       ' ***',/)
c      IF ( maxj .GT. jtopuser ) THEN
c        WRITE(*,994) maxj
c        WRITE (lunit,994) maxj
c        STOP 'Terminated in find_jmax.f'
c      ENDIF

      write(lunit,112)str,maxj
      write(lunit,*)
      write(lunit,*)

      maxemax =  vmax + etopuser + barr(1)*(maxj+1)*(maxj+2)			! maximum energy for allocating arrays 

      call sestamp('find_jmax',2)

      return

      end subroutine
