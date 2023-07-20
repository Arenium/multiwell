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

      subroutine veff(jval,varray)
      use DECLARE_MOD				! use module DECLARE_MOD
      include 'DECLARE_INC.inc'		! include 'DECLARE_MOD.inc'
      integer jval
      real*8  varray(rcnt+ntts+1),temax
      character frmat*60,namefile*80
      real*8  x(3),y(3),fx,fy,fit(3),er1,vmin
      real*8  tx(3),ty(3),fa(3)

 100  format(1x,i3,2(2x,es13.6),3x,i7,2(3x,es13.6))
 201  format(1x,a3,2(2x,a13),3x,a7,2(3x,a13))
 200  format(2(3x,a),10x,a,9x,a,6x,a,11x,a)
 300  format(i5,2(2x,es15.5,2x,f5.2),2x,f7.2)
 400  format(i5,2(2x,es15.5,2x,f5.2),'*',1x,f7.2)
 500  format(a5,2(2x,a15,2x,a5),2x,a7)
 600  format('for j: ',I5,' vmax@',F10.2,'cm-1',1x,F10.1,'angstrom')
 700  format(A,ES10.2,'cm-1',2x,ES10.2,'angstrom')
 701  format(A,ES10.2,'cm-1')
 800  format(a,2x,i4,2x,F15.2)


c  calculated the effective potential for the passed value of j


      do i=1,rcnt+ntts+1
           varray(i)=delh(i)+barr(i)*(jval)*(jval+1)
      end do

c  find veff max among trial ts, do not include reactants or products

      if(ntts.gt.0)then
         vmax = varray(rcnt+1)
         do i=rcnt+1,rcnt+ntts
            if(varray(i).ge.vmax)then
               vmax=varray(i)
               vmaxi=i
            end if
         end do
      end if

c  find centrifugal maximum by inter/extrapolation

      if(ntts.ge.3)then
         if(vmaxi.eq.rcnt+ntts)then                                     ! if max is at the end of ntts
            do j=-2,0
               tx(3+j)=distl(vmaxi+j)
               ty(3+j)=varray(vmaxi+j)
            end do
         elseif(vmaxi.eq.rcnt+1)then                                    ! if max is at the start of ntts
            do j=0,2
               tx(1+j)=distl(vmaxi+j)
               ty(1+j)=varray(vmaxi+j)
            end do
         else
            do j=-1,1
               tx(2+j)=distl(vmaxi+j)
               ty(2+j)=varray(vmaxi+j)
            end do
         end if


         call parabfit(tx,ty,3,fa,er1)                                  ! test curvature
         x(1)=-fa(2)/(2.0d0*fa(3))
         y(1)=fa(1)+fa(2)*x(1)+fa(3)*x(1)*x(1)
         write(lunit,600)jval,vmax,distl(vmaxi)
         write(lunit,700)' centrifugal max @',y(1),x(1)
         write(lunit,701)' energy difference',y(1)-vmax

      end if

c  save the following
c     veff(1), the reactants in v1l
c     vmax in vefmax
c     vmax-veff(i) in vlist

         v1l(jval+1)=varray(1)
       vlist(jval+1)=vmax-varray(1)
      vefmax(jval+1)=vmax

c  find epsilon, the energy difference from the bottom of veff to vmax+etop

      do i=1,rcnt+ntts+1
         if(vmax-varray(i).lt.0)then
            epsil(i)=etopuser
         else
            epsil(i)=(vmax+etopuser)-varray(i)
         end if
      end do

c  fix reactant epsilons

      temax=epsil(1)
      do i=1,rcnt
         if(epsil(i).gt.temax)temax=epsil(i)
      end do

      do i=1,rcnt
         epsil(i)=temax
      end do
           
c  fix product epsilons

      temax=epsil(rcnt+ntts+1)
      do i=rcnt+ntts+1,rcnt+ntts+1
         if(epsil(i).gt.temax)temax=epsil(i)
      end do

      do i=rcnt+ntts+1,rcnt+ntts+1
         epsil(i)=temax
      end do
  

c  the largest value for epsilon will be for well based reactants at j=0
c     this gives the longest density of states vector

      if(jval.eq.0)then
         epsmax=nint(epsil(1)/de)+1
         if(nint(epsil(i)/de)+1.gt.epsmax)then                  ! make sure to find largest epsilon
            epsmax=nint(epsil(i)/de)+1
         end if
      end if

      
c  write veff out to veff.txt
      if(rcnt+ntts+1.gt.1)then
         if((inputfile.eq.'').or.(inputfile.eq.'ktools.dat'))then
            namefile='veff.veff'
         else
            namefile=trim(fileroot)//'.veff'
         end if
      
         if(jval.eq.0)then
            open(veffout,file=namefile)
            call stamp(veffout,2)
            call dnt(veffout)
            write(veffout,*)'run id: ',tstmp
            write(veffout,*)
            write(veffout,*)
         else
            open(veffout,file=namefile,position='append')
         end if

         write(frmat,'("(5x,",i0,"(es20.8))")')rcnt+ntts+1
         write(veffout,frmat)(vef(k),k=1,rcnt+ntts+1)
         close(veffout)
      end if

c  write to log file

      write(lunit,*)'start veff'
      write(lunit,*)
      write(lunit,*)'j:',jval
      write(lunit,*)
      write(lunit,201)'i','delh','b','j*(j+1)','veff','epsilon'
      do i=1, rcnt+ntts+1
         write(lunit,100)i, delh(i),barr(i),jval*(jval+1),
     $                   varray(i),epsil(i)
      end do
      write(lunit,800)'veff max is at point #',vmaxi,vmax
      write(lunit,*)
      call sestamp('veff',2)

      return
      
      end subroutine
