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
      subroutine indexx(n,arrin,indx)
      implicit real*8 (a-h,o-z)
      dimension arrin(n),indx(n)

      do 11 j=1,n
         indx(j)=j
 11   continue
      nn = n/2 + 1
      ir=n
 10   continue
      if (nn.gt.1) then
         nn = nn-1
         indxt=indx(nn)
         q=arrin(indxt)
      else
         indxt=indx(ir)
         q=arrin(indxt)
         indx(ir)=indx(1)
         ir=ir-1
         if(ir.eq.1)then
            indx(1)=indxt            
            return
         endif
      endif
      i=nn
      j=2*nn
 20   if(j.le.ir)then
         if(j.lt.ir)then
            if(arrin(indx(j)).lt.arrin(indx(j+1)))j=j+1
         endif
         if(q.lt.arrin(indx(j)))then
            indx(i)=indx(j)
            i=j
            j=j+j
         else
            j=ir+1
         endif
         go to 20
      endif
      indx(i)=indxt
      go to 10
C 
      end
