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
      subroutine canon_writeout(kif,kir,rqtot,tsqtot,pqtot,qvib,qrot,
     $                          qhind,qelec,qtrans)
      use DECLARE_MOD				! use module DECLARE_MOD
      include 'DECLARE_INC.inc'		! include 'DECLARE_MOD.inc'
      real*8   kif(nt,ntts)
      real*8   kir(nt,ntts)
      real*8   rqtot(nt)
      real*8   tsqtot(nt,ntts)
      real*8   pqtot(nt)
      real*8   qvib(nt,rcnt+ntts+pcnt)
      real*8   qrot(nt,rcnt+ntts+pcnt)
      real*8   qhind(nt,rcnt+ntts+pcnt)
      real*8   qelec(nt,rcnt+ntts+pcnt)
      real*8   qtrans(nt,rcnt+ntts+pcnt)
      integer  mde, line
      character namefile*80
 400  format(a9,2x,a3,2x,a8,2x,15(a12,2x))
 500  format(a9,2x,i3,2x,f8.2,2x,10(es12.5,2x))
 600  format(f10.2,2x,100(es13.5e3,2x))
 700  format(a10,2x,100(f13.2,2x))
 800  format('canonical rate constant data can be found in:',1x,a)
 900  format(a7,1x,'reaction canonical rate constants for all temperatur
     $es and transition states')
991   FORMAT (A180)

      mde = 0

c write partition functions to log file

      do i=1,nt

         write(lunit,400)'structure','#','T(K)','qele','qtrans','qvib',
     $                   'qrot','qhin','qtot','k(t)f','k(t)r'
         write(lunit,*)('-',j=1,138)

         do j=1,rcnt                                                    ! reactant partition functions
            if(j.eq.rcnt)then
               write(lunit,500)'reactant',j,temps(i),qelec(i,j),
     $         qtrans(i,j),qvib(i,j),qrot(i,j),qhind(i,j),rqtot(i)
            else
               write(lunit,500)'reactant',j,temps(i),qelec(i,j),
     $               qtrans(i,j),qvib(i,j),qrot(i,j),qhind(i,j)
            end if
         end do
         
         do j=1,ntts                                                    ! transition state partition functions
            k=rcnt+j
            write(lunit,500)'ts',j,temps(i),qelec(i,k),qtrans(i,k),
     $            qvib(i,k),qrot(i,k),qhind(i,k),tsqtot(i,j),kif(i,j),
     $            kir(i,j)
         end do
   
         do j=1,pcnt                                                    ! product partition functions
            k=rcnt+ntts+j
            if(j.eq.pcnt)then
               write(lunit,500)'product',j,temps(i),qelec(i,k),
     $               qtrans(i,k),qvib(i,k),qrot(i,k),qhind(i,k),pqtot(i)
            else
               write(lunit,500)'product',j,temps(i),qelec(i,k),
     $               qtrans(i,k),qvib(i,k),qrot(i,k),qhind(i,k)
            end if
         end do
         write(lunit,*)
         write(lunit,*)
      end do

c write data to canonical file
      if(ntts.ge.1)then
         if( (inputfile.eq.'') .or. (inputfile.eq.'ktools.dat') ) then
            namefile='kt-canon.canonical'
         else
            namefile=trim(fileroot)//".canonical"
         end if
         write(lunit,800)trim(namefile)
 
         open(unit=canunit,file=namefile)                               ! open file for writing
         call stamp(canunit,2)
         call dnt(canunit)
         write(canunit,*)'run id: ',tstmp  
         write(canunit,*)
         write(canunit,*) title
         write(canunit,*)
         do line = 1, nlines
           write(canunit,*) commentline(line)
         end do
         write(canunit,*)
1000     if(mde.eq.0)then
            write(canunit,900)'forward'
               write(canunit,700)'r',(temps(i),i=1,nt)
               do j=1,ntts
                  write(canunit,600)distl(j+rcnt),(kif(i,j),i=1,nt)
               end do
               write(canunit,*)
               write(canunit,*)
               call canon_min(kif,cratef,ucratef)
               mde=mde+1
               goto 1000
         elseif(pcnt.ge.1)then
            write(canunit,900)'reverse'
               write(canunit,700)'r',(temps(i),i=1,nt)
               do j=1,ntts
                  write(canunit,600)distl(j+rcnt),(kir(i,j),i=1,nt)
               end do
               write(canunit,*)
               write(canunit,*)
               call canon_min(kir,crater,ucrater)
               mde=mde+1
         end if
      end if

      return

      end subroutine


