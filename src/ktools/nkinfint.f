c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c
c copyright (c) 2014 john r. barker, jason a. sonk
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


      subroutine nkinfint(temps,nt,mat,jay,ee,bee,symm,de,logunit,q)
      implicit none
      integer nt,ee,jay,z,i,j,k,logunit,majorj
      real*8  mat(jay,ee),bee(jay),symm,q(nt)
      real*8  temps(nt),norm,term,inte,contest,pterm

      real*8  converge,de,kb,maxterm
      real*8  smallcon
      parameter ( converge = 1.0d-8, kb= 0.69503476d0 )
      parameter ( smallcon = 1.0d-8 )
      logical converged

 100  format (3(A6,2x),3(A15,2x))
 102  format (A6,2x,A10,2x,2(A6,2x),3(A15,2x),A5,2x)
 200  format (3(I6,2x),3(ES15.6,2x))
 201  format (4(I6,2x),3(ES15.6,2x))
 202  format (I6,2x,F10.2,2x,2(I6,2x),3(ES15.6,2x),I5,2x)


      contest=1.0d0
      
      write(*,*)'Integrating over E & J'
      do i=1,nt
         converged=.false.
         write(*,FMT="(A1,A,t21,F6.2,A)",ADVANCE="NO") achar(13),
     $ " Percent Complete: ", (real(i)/real(nt))*100.0, "%"
         norm=0.0d0
         j=1
         term=1.0d0
         do while((.not.converged).and.(j.le.jay))
            inte=0.0d0
            do k=1,ee
               term = -1.0d0*de*(k-1)
               term = term/(kb*temps(i))
               term = dexp(term)
               if(isnan(mat(j,k)*term))cycle
               pterm=inte
               inte = inte + (mat(j,k)*term)
               if(k.gt.ee/4.and.pterm.ne.0.0d0)then
                  if( ( (inte/pterm)-1.0d0 ) .le. converge )exit
               end if
            end do
            term = bee(j)
            term = term / (kb * temps(i))
            term = -1.0d0 * term
            term = dexp(term)
            term = term*inte*de
            if(isnan(term))then
               write(logunit,*)'hit nan for temps:',temps(i),'j:',j-1
            end if

            if(i.eq.1)then
               if(j.eq.1)then
                  write(logunit,102)'#','Temps','J','E','Norm','Term',
     $                                  'Maxterm','@J'
                  write(logunit,*)('-',z=1,92)
               end if
            end if

            pterm=norm
            norm = norm + term
            if(j.eq.1)contest=term
            if(term.gt.contest)then
               contest=norm
               majorj=j
            end if
            if( (j.gt.jay/4) .and. (pterm.ne.0.0d0) )then
               if( ( (norm/pterm)-1.0d0 ) .le. converge )then
                  converged=.true.
               end if
            end if
            j = j + 1
         end do

         write(logunit,202)i,temps(i),j-2,k,norm,term,contest,majorj
         q(i) = norm*(1.0d0/symm)

      end do
      write(*,*)

      end subroutine   
