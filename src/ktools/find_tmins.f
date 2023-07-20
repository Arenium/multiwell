c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c
c copyright (c) 2017 john r. barker, jason a. sonk
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
      subroutine find_tmins(dat,mins,n,mincnt,hcheck)
      integer i,n,m,mincnt,mins(n),loopct
      real*8  dat(n),d
      integer fmax1,fmin1,fmax2,oldmin
      real*8  hcheck,lc,rc
      logical minskip

      mincnt=0
      do i=1,n
        mins(i)=0
      end do
      
      i=2
      do while(i.le.n)
         loopct=0
         minskip=.false.
         if(i.eq.2.and.dat(i).gt.dat(i-1))then                         ! Check to see if we're marching uphill
c            write(*,*)'checking i=2'
            minskip=.true.
            oldmin=1
         end if

 100     continue
         if(i.ge.n)exit
c         write(*,*)'100',loopct,i,n
         
         if(i.le.n)then
         do while(dat(i).ge.dat(i-1).and.i.le.n)                        ! march uphill to find max1
            i=i+1
            if(i.gt.n)exit
         end do
         end if

         fmax1=i-1

 200     continue
         if(i.ge.n)exit
c         write(*,*)'200',loopct,i,n

         if(loopct.ge.5)then
            i=i+1
         end if

         if(i.le.n)then
         do while(dat(i-1).ge.dat(i).and.i.le.n)                        ! march downhill to find min1
            i=i+1
            if(i.gt.n)exit
         end do
         end if
      
         fmin1=i-1
 
 300     continue
         if(i.ge.n)exit
c         write(*,*)'300',loopct,i,n

         if(minskip)then
c            write(*,*)'in minskip'
            if(dat(oldmin).lt.dat(fmin1))then
               fmin1=oldmin
            end if
         end if

         if(i.le.n)then
            do while((dat(i).ge.dat(i-1)))                    ! march uphill to find max2
               i=i+1
               if(i.gt.n)exit
            end do
         end if
         fmax2=i-1

         lc=dat(fmax1)/dat(fmin1)                                       ! check ratios
         rc=dat(fmax2)/dat(fmin1)
c         write(*,*)'Max @',fmax1,' Min @',fmin1,' Max @',fmax2,lc,rc

         if( ( (lc.lt.hcheck).or.(rc.lt.hcheck) ) .and. (i.lt.n) ) then
c            write(*,*)'min skipped'
            minskip=.true.
            oldmin=fmin1
            loopct=loopct+1
            go to 200
         end if

         if(mincnt.ge.1)then
            if(fmin1-mins(mincnt).eq.2)then
               if(dat(fmin1).lt.dat(mins(mincnt)))then
                  mins(mincnt)=fmin1
               end if
               go to 100
            end if
         end if

         if(fmax1.eq.fmin1.and.fmin1.eq.fmax2)then
            exit
         else
            mincnt=mincnt+1
            mins(mincnt)=fmin1
         end if
      end do



      if(mincnt.eq.0)then
c         write(*,*)'mincnt=0'
         fmin1=1
         do i=1,n
            if(dat(i).lt.dat(fmin1))fmin1=i
         end do
         mincnt=mincnt+1
         mins(mincnt)=fmin1
      end if
     
      return
      
      end subroutine
