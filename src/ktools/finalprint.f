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
      subroutine finalprint
      use DECLARE_MOD				! use module DECLARE_MOD
      include 'DECLARE_INC.inc'		! include 'DECLARE_MOD.inc'

 103  format(i3,2x,F9.2,2x,4(ES15.5,2x))
 104  format(a3,4x,a6,1x,3(a14,2x),a14,2x)
 105  format(3x,2x,6x,4x,a15,4x,a15,3x,a14x,2x)
 106  format(1x,'| ',a8,' rate and equilibrium constants |')
 107  format(a3,4x,a6,2x,3(a14,2x),a14,2x)
 108  format(3x,2x,6x,4x,a15,2x,a15,1x,a14x,2x)
 109  format(1x,'| ',a8,' rate constants |')

      if(whatdo.ne.'canon')then
         if(pcnt.ne.0)then
            if(.not.manymins)then
               write(canunit,*)('=',i=1,43)
               write(canunit,106)'forward'     
               write(canunit,*)('=',i=1,43)
               write(canunit,*)     

               write(canunit,105)'partition'
               write(canunit,105)'function','microcanonial'
               write(canunit,104)'#','T(K)','kinf','kinf','Keq'
               write(canunit,*)('-',i=1,67)
               do i=1,nt
                  write(canunit,103)i,temps(i),cratef(i),mrate(i),
     $                              queues(i)
               end do
               write(canunit,*)
               if(pcnt.ne.0)then
                  write(canunit,*)('=',i=1,43)
                  write(canunit,106)'reverse'     
                  write(canunit,*)('=',i=1,43)
                  write(canunit,*)     
                  write(canunit,105)'partition'
                  write(canunit,105)'function','microcanonial'
                  write(canunit,104)'#','T(K)','kinf','kinf','Keq'
                  write(canunit,*)('-',i=1,67)
                  do i=1,nt
                    write(canunit,103)i,temps(i),crater(i),
     $                                mrate(i)/queues(i),queues(i)
                  end do
               end if
            else
               write(canunit,*)('=',i=1,43)
               write(canunit,106)'forward'     
               write(canunit,*)('=',i=1,43)
               write(canunit,*)     

               write(canunit,108)'partition','','unified'
               write(canunit,105)'function','microcanonial',
     $                           'microcanonical'
               write(canunit,104)'#','T(K)','kinf','kinf','kinf','Keq'
               write(canunit,*)('-',i=1,82)
               do i=1,nt
                  write(canunit,103)i,temps(i),cratef(i),mrate(i),
     $                              umrate(i),queues(i)
               end do
               write(canunit,*)
               if(pcnt.ne.0)then
                  write(canunit,*)('=',i=1,43)
                  write(canunit,106)'reverse'     
                  write(canunit,*)('=',i=1,43)
                  write(canunit,*)     
                  write(canunit,108)'partition','','unified'
                  write(canunit,105)'function','microcanonial',
     $                              'microcanonial'
                  write(canunit,104)'#','T(K)','kinf','kinf','kinf',
     $                              'Keq'
                  write(canunit,*)('-',i=1,67)
                  do i=1,nt
                    write(canunit,103)i,temps(i),crater(i),
     $                                mrate(i)/queues(i),
     $                                umrate(i)/queues(i),queues(i)
                  end do
               end if
            end if
         else
            if(.not.manymins)then
               write(canunit,*)('=',i=1,27)
               write(canunit,109)'forward'     
               write(canunit,*)('=',i=1,27)
               write(canunit,*)     

               write(canunit,105)'partition'
               write(canunit,105)'function','microcanonial'
               write(canunit,104)'#','T(K)','kinf','kinf'
               write(canunit,*)('-',i=1,55)
               do i=1,nt
                  write(canunit,103)i,temps(i),cratef(i),mrate(i)
               end do
               write(canunit,*)
               if(pcnt.ne.0)then
                  write(canunit,*)('=',i=1,27)
                  write(canunit,109)'reverse'     
                  write(canunit,*)('=',i=1,27)
                  write(canunit,*)     
                  write(canunit,105)'partition'
                  write(canunit,105)'function','microcanonial'
                  write(canunit,104)'#','T(K)','kinf','kinf'
                  write(canunit,*)('-',i=1,55)
                  do i=1,nt
                    write(canunit,103)i,temps(i),crater(i),
     $                                mrate(i)/queues(i)
                  end do
               end if
            else
               write(canunit,*)('=',i=1,27)
               write(canunit,109)'forward'     
               write(canunit,*)('=',i=1,27)
               write(canunit,*)     

               write(canunit,108)'partition','','unified'
               write(canunit,105)'function','microcanonial',
     $                           'microcanonical'
               write(canunit,104)'#','T(K)','kinf','kinf','kinf'
               write(canunit,*)('-',i=1,82)
               do i=1,nt
                  write(canunit,103)i,temps(i),cratef(i),mrate(i),
     $                              umrate(i)
               end do
               write(canunit,*)
               if(pcnt.ne.0)then
                  write(canunit,*)
                  write(canunit,104)'#','T(K)','kinf(pf)R'
                  write(canunit,*)('-',i=1,82)
                  do i=1,nt
                    write(canunit,103)i,temps(i),crater(i)
                  end do
               end if
            end if  
         end if
      end if
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      write(canunit,*)
      write(canunit,*)
      write(canunit,*)
      if(pcnt.eq.0)then
         write(canunit,*)('%',i=1,57)
         write(canunit,*)'%%  FINAL RECOMMENDED REACTION RATE CONSTANTS 
     $   ',(' ',i=1,6),'%%'
         write(canunit,*)('%',i=1,57)

      else
         write(canunit,*)('%',i=1,76)
         write(canunit,*)'%%  FINAL RECOMMENDED REACTION RATE CONSTANTS 
     $AND EQUILIBRIUM CONSTANTS   %%'
         write(canunit,*)('%',i=1,76)

      end if

      write(canunit,*)
      write(canunit,992) de , etopuser , dj , jtopuser , maxj
      
992   FORMAT (15x,'    USER-SELECTED',/,
     & '   Energy grain:    ', F6.2,' cm-1',/,
     & '   Energy above TS: ', I6,  ' cm-1',/,
     & '   J grain:         ', I6,/,
     & '   Maximum J:       ', I6,6x,'(J ='I4,
     &                                   ': intrinsic upper limit)',/,
     &                   32x,'(used if smaller than user-selected)',/)

      if(manymins)then
      if(whatdo.eq.'canon')then
         if(pcnt.eq.0)then
            write(canunit,104)'','','kinf'
            write(canunit,107)'#','T(K)','forward'
            write(canunit,*)('-',i=1,32)
            do i=1,nt
               write(canunit,103)i,temps(i),ucratef(i)
            end do
         else
            write(canunit,104)'','','kinf','kinf'
            write(canunit,107)'#','T(K)','forward','reverse','Keq'
            write(canunit,*)('-',i=1,65)
            do i=1,nt
            write(canunit,103)i,temps(i),ucratef(i),ucrater(i),queues(i)
            end do
         end if
      else
         if(pcnt.eq.0)then
            write(canunit,104)'','','kinf'
            write(canunit,107)'#','T(K)','forward'
            write(canunit,*)('-',i=1,32)
            do i=1,nt
              write(canunit,103)i,temps(i),umrate(i)
            end do
         else
            write(canunit,104)'','','kinf','kinf'
            write(canunit,107)'#','T(K)','forward','reverse','Keq'
            write(canunit,*)('-',i=1,65)
            do i=1,nt
             write(canunit,103)i,temps(i),umrate(i),umrate(i)/queues(i),
     $                 queues(i)
            end do
         end if
      end if
      else
      if(whatdo.eq.'canon')then
         if(pcnt.eq.0)then
            write(canunit,104)'','','kinf'
            write(canunit,107)'#','T(K)','forward'
            write(canunit,*)('-',i=1,32)
            do i=1,nt
               write(canunit,103)i,temps(i),cratef(i)
            end do
         else
            write(canunit,104)'','','kinf','kinf'
            write(canunit,107)'#','T(K)','forward','reverse','Keq'
            write(canunit,*)('-',i=1,65)
            do i=1,nt
            write(canunit,103)i,temps(i),cratef(i),crater(i),queues(i)
            end do
         end if
      else
         if(pcnt.eq.0)then
            write(canunit,104)'','','kinf'
            write(canunit,107)'#','T(K)','forward'
            write(canunit,*)('-',i=1,32)
            do i=1,nt
              write(canunit,103)i,temps(i),mrate(i)
            end do
         else
            write(canunit,104)'','','kinf','kinf'
            write(canunit,107)'#','T(K)','forward','reverse','Keq'
            write(canunit,*)('-',i=1,65)
            do i=1,nt
             write(canunit,103)i,temps(i),mrate(i),mrate(i)/queues(i),
     $                 queues(i)
            end do
         end if
      end if
      end if 

            

      return

      end subroutine
