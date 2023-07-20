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
      subroutine canon_rates
      use DECLARE_MOD				! use module DECLARE_MOD
      include 'DECLARE_INC.inc'		! include 'DECLARE_INC.inc'
      real*8  tt,bj,bk,cstop,sstop,hstop
      real*8  he(nt,rcnt+ntts+pcnt)
      real*8   s(nt,rcnt+ntts+pcnt)
      real*8  cp(nt,rcnt+ntts+pcnt)
      real*8  qelec(nt,rcnt+ntts+pcnt)
      real*8  qqvib(nt,rcnt+ntts+pcnt)
      real*8  qrot(nt,rcnt+ntts+pcnt)
      real*8  qhind(nt,rcnt+ntts+pcnt)
      real*8  qtrans(nt,rcnt+ntts+pcnt)
      real*8  rqtot(nt)
      real*8  tsqtot(nt,ntts)
      real*8  pqtot(nt)
      real*8  kfor(nt,ntts)
      real*8  krev(nt,ntts)
      real*8  dsr(nt),dhr(nt),dcp(nt),dg(nt),ak(nt)
      real*8  an,cxx,ecrit,hxx,minq,qq,sxx,we,zpe
      real*8  qele,qvib,qroq,qhin,qstop,tqstop
      integer*4 ng,ss
      logical krof
 100  format(I3,1x,F10.2,1x,5(ES15.6,2x))
 200  format(I3,1x,10(ES10.4,2x))
111   format(I1,2x,I2,2x,10(ES15.6,2x))
222   format(i3,2x,F10.2,2x,1(ES10.4,2x),4(F10.2,2x))
223   format(i3,2x,F10.2,2x,1(ES10.4,2x),4(F10.2,2x))
333   format(A3,2x,A10,2x,1(A10,2x),4(A10,2x))
334   format(A3,2x,A10,2x,1(A10,2x),4(A10,2x))
99004 FORMAT (' THERMODYNAMIC QUANTITIES',/,
     $       ' Standard State = molecule/cc'/,3x,
     &       ' DelS         units: cal/K/mole'/3x,
     &       ' DelH         units: kcal/mole'/3x,
     &       ' DelCp        units: cal/K/mole'/3x,
     &       ' DelG         units: kcal/mole'/3x,
     $       ' Keq => Exp[(delS-delH/T)/R]'/3x,
     $       'vKeq => kforward/kreverse'//)

      call sestamp('canon_rates',1)

c  initialization 
      do i=1,nt
         do j=1,rcnt+ntts+pcnt
            he(i,j)=0.0d0
             s(i,j)=0.0d0
            cp(i,j)=0.0d0
      
            cxx = 0.0d0
            hxx = 0.0d0
            sxx = 0.0d0

            qelec(i,j)=1.0d0
            qqvib(i,j)=1.0d0
            qrot(i,j)=1.0d0
            qhind(i,j)=1.0d0
            qtrans(i,j)=1.0d0
         end do
         rqtot(i)=1.0d0
         pqtot(i)=1.0d0
         do j=1,ntts
            tsqtot(i,j)=1.0d0
         end do
      end do

c  calculate canonical rate constants

      write(*,*)'Calculating Canonical Rate Contants'
      do i=1,nt
         write(*,FMT="(A1,A,t21,F6.2,A)",ADVANCE="NO") achar(13),
     $ "  Percent Complete: ",(real(i)/real(nt))*100.0, "%"


         tt=temps(i)                                                    ! set temperature


         do j=1,rcnt+ntts+pcnt                                          ! calculate partition functions

            krof=.false.

            s(i,j) = rgascm1*dlog(DBLE(Sopt(j))/DBLE(Sym(j)))           ! entropy from optical isomers/symmetry

            qelec(i,j) = qele(tt,j,nele(j),elev,gele,datpts,cxx,sxx,hxx)       ! qelectronic
            he(i,j) = he(i,j) + rgascm1*hxx*tt        ! H(T)-H(0)
            cp(i,j) = cp(i,j) + rgascm1*cxx
             s(i,j) =  s(i,j) + rgascm1*sxx

            qtrans(i,j)=(DSqrt( ( amu(j)*tt )**3) )*1.879333D+20        ! qtranslation (std. state= 1 molecule/cc)
            he(i,j) = he(i,j) + 2.5d0*rgascm1*tt
            cp(i,j) = cp(i,j) + 2.5d0*rgascm1
             s(i,j) =  s(i,j) + rgascm1*(2.5d0 + dlog(qtrans(i,j)))

            if(ndof(j).ne.0)then                                        ! if polyatomic
               do k=1, ndof(j)                                          ! loop over all degrees of freedom
                  we=wel(j,k)
                  an=anhl(j,k)
                  ng=int(ngl(j,k))
                  if( idofl(j,k).eq.'vib' ) then                        ! qvibration
                     qq=qvib(tt,we,an,ng,cxx,sxx,hxx,zpe)
                     qqvib(i,j) = qqvib(i,j) * qq
                     cp(i,j) = cp(i,j) + cxx * rgascm1
                     he(i,j) = he(i,j) + hxx * rgascm1 * tt
                      s(i,j) =  s(i,j) + sxx * rgascm1

                  elseif( idofl(j,k).eq.'jro') then                     !  look for j-rotor
                     we=wel(j,k)
                     bj=we
                     do l=1,ndof(j)
                        if( idofl(j,l).eq.'kro')then                    !  chek for k-rotor
                           krof=.true.
                           an=wel(j,l)
                           bk=an
                           ng=int(ngl(j,l))
                        end if
                     end do

                     if(.not.krof)then                                  !  if there is no k-rotor then treat j-rotor as a qro dof
                        qq=qroq(tt,we,ng,int(an),cxx,sxx,hxx)
                        qrot(i,j) = qrot(i,j) * qq
                        cp(i,j) = cp(i,j) + cxx * rgascm1
                        he(i,j) = he(i,j) + hxx * rgascm1 * tt
                         s(i,j) =  s(i,j) + sxx * rgascm1
                     else
                        qq=tqstop(tt,we,an,ng,cxx,sxx,hxx)              ! otherwise use top
                        qrot(i,j) = qrot(i,j) * qq
                        cp(i,j) = cp(i,j) + cxx * rgascm1
                        he(i,j) = he(i,j) + hxx * rgascm1 * tt
                         s(i,j) =  s(i,j) + sxx * rgascm1
                     end if

                  elseif( idofl(j,k).eq.'qro')then                      ! qrotation 
                     qq=qroq(tt,we,ng,int(an),cxx,sxx,hxx)
                     qrot(i,j) = qrot(i,j) * qq
                     cp(i,j) = cp(i,j) + cxx * rgascm1
                     he(i,j) = he(i,j) + hxx * rgascm1 * tt
                      s(i,j) =  s(i,j) + sxx * rgascm1
                  elseif( idofl(j,k)(1:2).eq.'hr')then                  ! qhinderedrotation
                     do l=1,valmax(j,k)
                        evh(l)=evhl(j,k,l)
c                        write(*,*)l,evh(l)
                     end do
c               write(*,*)'qhin in ',j,tt,valmax(j,k),ng,cxx,sxx,hxx,zpe,
c     $ (evh(l),l=1,valmax(j,k))
                     qq=qhin(tt,evh,valmax(j,k),ng,cxx,sxx,hxx,zpe)
c               write(*,*)'qhin out',j,tt,valmax(j,k),ng,cxx,sxx,hxx,zpe,
c     $ (evh(l),l=1,valmax(j,k))
                     qhind(i,j)= qhind(i,j) * qq
                     cp(i,j) = cp(i,j) + cxx * rgascm1
                     he(i,j) = he(i,j) + hxx * rgascm1 * tt
                      s(i,j) =  s(i,j) + sxx * rgascm1
                  end if
               end do
            end if
         end do
      end do
      write(*,*)


c  multiply parition fuctions together to get total partition function


      do i=1,nt
         do j=1,rcnt                                                    ! reactant total partition function
            rqtot(i)=rqtot(i)*qtrans(i,j)                               ! qtrans
            rqtot(i)=rqtot(i)*qelec(i,j)                                ! qelec
            rqtot(i)=rqtot(i)*qqvib(i,j)                                ! qvib
            rqtot(i)=rqtot(i)*qrot(i,j)                                 ! qrot
            rqtot(i)=rqtot(i)*qhind(i,j)                                ! qhind
            rqtot(i)=rqtot(i)/sym(j)                                    ! symmetry
         end do

         do j=1,ntts                                                    ! ts total partition function
            tsqtot(i,j)=tsqtot(i,j)*qtrans(i,j+rcnt)                    ! qtrans
            tsqtot(i,j)=tsqtot(i,j)*qelec(i,j+rcnt)                     ! qelec
            tsqtot(i,j)=tsqtot(i,j)*qqvib(i,j+rcnt)                     ! qvib
            tsqtot(i,j)=tsqtot(i,j)*qrot(i,j+rcnt)                      ! qrot
            tsqtot(i,j)=tsqtot(i,j)*qhind(i,j+rcnt)                     ! qhind
            tsqtot(i,j)=tsqtot(i,j)/sym(rcnt+j)                         ! symmetry
         end do

c         do j=rcnt+ntts+1,rcnt+ntts+pcnt                                ! product total partition function
c            pqtot(i)=pqtot(i)*qtrans(i,j)                               ! qtrans
c            pqtot(i)=pqtot(i)*qelec(i,j)                                ! qelec
c            pqtot(i)=pqtot(i)*qqvib(i,j)                                ! qvib
c            pqtot(i)=pqtot(i)*qrot(i,j)                                 ! qrot
c            pqtot(i)=pqtot(i)*qhind(i,j)                                ! qhind
c            pqtot(i)=pqtot(i)/sym(j)                                    ! symmetry
c         end do

         do j=1,pcnt                                                    ! product total partition function
            pqtot(i)=pqtot(i)*qtrans(i,j+rcnt+ntts)                     ! qtrans
            pqtot(i)=pqtot(i)*qelec(i,j+rcnt+ntts)                      ! qelec
            pqtot(i)=pqtot(i)*qqvib(i,j+rcnt+ntts)                      ! qvib
            pqtot(i)=pqtot(i)*qrot(i,j+rcnt+ntts)                       ! qrot
            pqtot(i)=pqtot(i)*qhind(i,j+rcnt+ntts)                      ! qhind
            pqtot(i)=pqtot(i)/sym(j+rcnt+ntts)                          ! symmetry
         end do


      end do

c  find kforward and krev

      do i=1,nt
         do j=1,ntts
            ecrit=(-1.0d0*delhf(rcnt+j))/(kb*temps(i))                  ! forward critical energy
            kfor(i,j)=tsqtot(i,j)/rqtot(i)
            kfor(i,j)=kfor(i,j)*dexp(ecrit)
            kfor(i,j)=kfor(i,j)*((kbj*temps(i))/h)
            if(pcnt.ne.0)then
              ecrit=(-1.0d0*delhr(rcnt+j))/(kb*temps(i))                ! reverse critical energy
              krev(i,j)=tsqtot(i,j)/pqtot(i)
              krev(i,j)=krev(i,j)*dexp(ecrit)
              krev(i,j)=krev(i,j)*((kbj*temps(i))/h)
            end if
         end do
      end do
         
c  write out all k rates to canonical file

      call canon_writeout(kfor,krev,rqtot,tsqtot,pqtot,qqvib,qrot,qhind,
     $qelec,qtrans)

c      do i=1,nt
c         minq=tsqtot(i,1)
c         do j=1,ntts
c            if(tsqtot(i,j).lt.minq)minq=tsqtot(i,j)
c         end do
c         write(*,100)i,temps(i),rqtot(i),pqtot(i),minq,kfor(i,1),
c     $                 krev(i,1)
c      end do

c  if reactants and products are present, calculate Keq

      if( (rcnt.ge.1) .and. (pcnt.ge.1) ) then

         do i=1,nt
            dsr(i)=0.0d0
            dhr(i)=delh(rcnt+ntts+1)
            dcp(i)=0.0d0
            do j=1,rcnt+ntts+pcnt
            
               if(reprod(j).eq.'reac')then
                  dsr(i)=dsr(i)- s(i,j)
                  dhr(i)=dhr(i)-he(i,j)
                  dcp(i)=dcp(i)-cp(i,j)
               elseif(reprod(j).eq.'prod')then
                  dsr(i)=dsr(i)+ s(i,j)
                  dhr(i)=dhr(i)+he(i,j)
                  dcp(i)=dcp(i)+cp(i,j)
               end if
            end do

            dg(i)=(dsr(i)-(dhr(i)/temps(i)))
            dg(i)=dg(i)/rgascm1
            ak(i)=dexp(dg(i))
            queues(i)=ak(i)
            dsr(i)=(dsr(i)/kcal2nu)
            dhr(i)=(dhr(i)/kcal2nu)
            dcp(i)=(dcp(i)/kcal2nu)
            dg(i)=dhr(i)-(temps(i)*dsr(i))
         end do

      do j=lunit,canunit
         write(j,*)
         write(j,*)'Reactants and Products found reporting Kequil.'
         write(j,*)
         write(j,99004)

         if((ntts.gt.0).and.(pcnt.gt.0))then
            write(j,333)"#","T(K)",'Keq','delS(rxn)','delH(rxn)',
     $                'delCp(rxn)','delG(rxn)'
            write(j,*)('-',i=1,73)

            do i=1,nt
               write(j,222)i,temps(i),ak(i),
     $                  dsr(i)*1.0d+03,dhr(i),dcp(i)*1.0d+03,dg(i)
            end do
         elseif((ntts.eq.0).and.(pcnt.gt.0))then
            write(j,334)"#","T(K)",'Keq','delS(rxn)','delH(rxn)',
     $                'delCp(rxn)','delG(rxn)'
            write(j,*)('-',i=1,73)

            do i=1,nt
               write(j,223)i,temps(i),ak(i),
     $                  dsr(i)*1.0d+03,dhr(i),dcp(i)*1.0d+03,dg(i)
            end do
         elseif((ntts.gt.0).and.(pcnt.eq.0))then
            write(j,334)"#","T(K)",'Keq','delS(rxn)','delH(rxn)',
     $                'delCp(rxn)','delG(rxn)'
            write(j,*)('-',i=1,73)

            do i=1,nt
               write(j,223)i,temps(i),kfor(i,1)/krev(i,1),
     $                  dsr(i)*1.0d+03,dhr(i),dcp(i)*1.0d+03,dg(i)
            end do
         end if
         write(j,*)

      end do

      end if

      write(lunit,*)
      
      call sestamp('canon_rates',2)

      return

      end subroutine
