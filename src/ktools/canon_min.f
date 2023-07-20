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
      subroutine canon_min(kit,trate,utrate)
      use DECLARE_MOD				! use module DECLARE_MOD
      include 'DECLARE_INC.inc'		! include 'DECLARE_MOD.inc'
      real*8     kit(nt,ntts),minlist(nt,ntts)
      real*8     tdat(ntts),fminlist(nt,ntts),urate(nt)
      real*8     trate(nt),utrate(nt)
      integer    tmins(ntts),fndex(nt,ntts)
      integer    fct(nt),starlist(nt)
      integer    mincnt(nt),ndex(nt,ntts)
      character  frmt*80
 499  format('unified canonical rate constants and contributing minima')
 500  format(i3,1x,f9.2,4x,i3,6x,100(f4.2,2x))
 501  format(100(a5,1x))
 502  format(100(es10.4,2x))
 503  format(i3,1x,f9.2,"*",3x,i3,6x,100(f4.2,2x))
 1500 format(i3,1x,f9.2,2x,es10.4,2x,i3,6x,100(f4.2,2x))
 1501 format(100(a5,1x))
 1502 format(100(es10.4,2x))
 1503 format(i3,1x,f9.2,"*",1x,es10.4,2x,i3,6x,100(f4.2,2x))

c  begin search for all local minima

      do i=1,nt
         mincnt(i)=0
         do j=2,ntts-1
            if(j.eq.2.and.kit(i,j-1).lt.kit(i,j))then                   ! if starting point is minimum
                          mincnt(i)=mincnt(i)+1
                  ndex(i,mincnt(i))=j-1
               minlist(i,mincnt(i))=kit(i,j-1)
            elseif(j.eq.ntts-1.and.kit(i,j+1).lt.kit(i,j))then          ! if ending point is minimum
                          mincnt(i)=mincnt(i)+1
                  ndex(i,mincnt(i))=j+1
               minlist(i,mincnt(i))=kit(i,j+1)
            elseif( ( kit(i,j-1) .gt. kit(i,j) ) .and.                  ! search of minimum among all other points
     $              ( kit(i,j)   .lt. kit(i,j+1) ) )then
                          mincnt(i)=mincnt(i)+1
                  ndex(i,mincnt(i))=j
               minlist(i,mincnt(i))=kit(i,j)
            end if
         end do
      end do

c  begin search for 'true' minima

      do i=1,nt
         do j=1,ntts                                                    ! copy into 1d array
            tdat(j)=kit(i,j)              
         end do

         call find_tmins(tdat,tmins,ntts,fct(i),microh)                 ! send into find true minimums
         do j=1,fct(i)
            fminlist(i,j)=kit(i,tmins(j))
            fndex(i,j)=tmins(j)
         end do

         urate(i)=0.0d0
         trate(i)=fminlist(i,1)
         do j=1,fct(i)
            if(fminlist(i,j).le.trate(i))trate(i)=fminlist(i,j)
            urate(i)=urate(i)+(1.0d0/fminlist(i,j))
         end do
         utrate(i)=1.0d0/urate(i)
      end do

c  write out all found minima

      if(maxval(mincnt).ne.0)then
      write(canunit,*)'reporting all local minima'
      write(canunit,*)

c  setting variable formatting for writeout

      if(maxval(mincnt).gt.1)then
         write(frmt,'("(a3,1x,a9,1x,a8,1x,a8,1x,",i0,"x,a8,1x)")')6*
     $                                          (maxval(mincnt)-1)
         write(canunit,frmt)'#','temp','# min','r(min)','min k(t)'
         write(canunit,*)('-',i=1,53+12*(maxval(mincnt)-1))
      else
         write(frmt,'("(a3,1x,a9,1x,a8,1x,a8,1x,",i0,"x,a8,1x)")')6
         write(canunit,frmt)'#','temp','# min','r(min)','min k(t)'
         write(canunit,*)('-',i=1,50+12)
      end if

c  writing out all local minima info

      do i=1,nt
         if(mincnt(i).lt.maxval(mincnt))then
            write(canunit,500,advance='no')i,temps(i),mincnt(i),
     $                  (distl(ndex(i,j)+1),j=1,mincnt(i))
            write(canunit,501,advance='no')
     $                  ('',j=1,maxval(mincnt)-mincnt(i))
            write(canunit,502)(minlist(i,j),j=1,mincnt(i))

         else
            write(canunit,500,advance='no')i,temps(i),mincnt(i),
     $                  (distl(ndex(i,j)+1),j=1,mincnt(i))
            write(canunit,502)(minlist(i,j),j=1,mincnt(i))
         end if
      end do

      if(maxval(mincnt).gt.1)then
         write(canunit,*)('-',i=1,53+12*(maxval(mincnt)-1))
      else
         write(canunit,*)('-',i=1,50+12)
      end if

      write(canunit,*)
      write(canunit,*)
      end if
c  writing out 'significant' and unified rate constants

      write(canunit,499)
      write(canunit,*)
      if(maxval(fct).ne.1)then
        write(frmt,'("(a3,1x,a9,1x,a10,2x,a8,1x,a8,1x,",i0,
     $                                   "x,a8,1x)")')6*(maxval(fct)-1)
        write(canunit,frmt)'#','temp','uk(t)','# min','r(min)',
     $                                                        'min k(t)'
        write(canunit,*)('-',i=1,50+12*(maxval(fct)-1))
        do i=1,nt
           do j=1,fct(i)
          if((fndex(i,j).eq.rcnt+1).or.(fndex(i,j).eq.rcnt+ntts))then
                  starlist(i)=1
               else
                  starlist(i)=0
               end if
           end do
           if(fct(i).lt.maxval(fct))then
             if(starlist(i).eq.0)then
             write(canunit,1500,advance='no')i,temps(i),
     $              1.0d0/urate(i),fct(i),(distl(fndex(i,j)),j=1,fct(i))
             else
             write(canunit,1503,advance='no')i,temps(i),
     $              1.0d0/urate(i),fct(i),(distl(fndex(i,j)),j=1,fct(i))
             endif
             write(canunit,1501,advance='no')
     $                                       ('',j=1,maxval(fct)-fct(i))
             write(canunit,1502)(fminlist(i,j),j=1,fct(i))

           else
             if(starlist(i).eq.0)then
             write(canunit,1500,advance='no')i,temps(i),
     $              1.0d0/urate(i),fct(i),(distl(fndex(i,j)),j=1,fct(i))
             else
             write(canunit,1503,advance='no')i,temps(i),
     $              1.0d0/urate(i),fct(i),(distl(fndex(i,j)),j=1,fct(i))
             end if
             write(canunit,1502)(fminlist(i,j),j=1,fct(i))
           end if
        end do
      else
        write(frmt,'("(a3,1x,a9,1x,a8,1x,a8,1x,",i0,"x,a8,1x)")')
     $                                          (maxval(fct))
        write(canunit,frmt)'#','temp','# min','r(min)','min k(t)'
        write(canunit,*)('-',i=1,35+12*(maxval(fct)))
        do i=1,nt
           do j=1,fct(i)
          if((fndex(i,j).eq.rcnt+1).or.(fndex(i,j).eq.rcnt+ntts))then
                  starlist(i)=1
               else
                  starlist(i)=0
               end if
           end do
           if(fct(i).lt.maxval(fct))then
             if(starlist(i).eq.0)then
             write(canunit,500,advance='no')i,temps(i),
     $                  fct(i),(distl(fndex(i,j)),j=1,fct(i))
             else
             write(canunit,503,advance='no')i,temps(i),
     $                  fct(i),(distl(fndex(i,j)),j=1,fct(i))
             end if
             write(canunit,501,advance='no')
     $                                       ('',j=1,maxval(fct)-fct(i))
             write(canunit,502)(fminlist(i,j),j=1,fct(i))
           else
             if(starlist(i).eq.0)then
             write(canunit,500,advance='no')i,temps(i),
     $                  fct(i),(distl(fndex(i,j)),j=1,fct(i))
             else
             write(canunit,503,advance='no')i,temps(i),
     $                  fct(i),(distl(fndex(i,j)),j=1,fct(i))
             end if
             write(canunit,502)(fminlist(i,j),j=1,fct(i))
           end if
        end do
      end if

      if(maxval(fct).ne.1)then
         write(canunit,*)('-',i=1,50+12*(maxval(fct)-1))
      else
         write(canunit,*)('-',i=1,35+12*(maxval(fct)))
      end if
      write(canunit,*)
      do i=1,nt
         if(starlist(i).eq.1)then
            write(canunit,*)'Temps denoted with * have minimum rate cons
     $tants at the end of the given data array.'
            write(canunit,*)'More data points may be needed to ensureacc
     $urate rate constants have been calculated.'
            write(canunit,*)
            write(canunit,*)
            exit
         end if
      end do

      return
      
      end subroutine
