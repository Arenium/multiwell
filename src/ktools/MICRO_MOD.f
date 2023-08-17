c
c john r. barker
c jrbarker@umich.edu
c university of michigan
c ann arbor, mi 48109-2143
c (734) 763 6239
c
c Jason A. Sonk and john R. Barker
c University of Michigan
c Ann Arbor, MI 48109
c
c Contact: jsonk@umichedu and jrbarker@umich.edu
c
c Copyright 2020
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

      MODULE MICRO_MOD
!
!    Module for Microcanonical computations with large arrays in program KTOOLS.
!
      use DECLARE_MOD				! use module DECLARE_MOD

!	allocated/deallocated in subroutine micro_rates
      real (8), ALLOCATABLE, DIMENSION(:,:) :: tmat 	! allocated: micro_rates ; deallocated: micro_rates
      real (8), ALLOCATABLE, DIMENSION(:,:) :: gmat 	! allocated: micro_rates ; deallocated: micro_rates
      real (8), ALLOCATABLE, DIMENSION(:,:) :: rpmat	! allocated: micro_rates ; deallocated: calc_kinf
      real (8), ALLOCATABLE, DIMENSION(:,:) :: minpmat	! allocated: micro_rates ; deallocated: calc_kinf
      real (8), ALLOCATABLE, DIMENSION(:,:) :: mingmat	! allocated: micro_rates ; deallocated: calc_kinf
      real (8), ALLOCATABLE, DIMENSION(:,:) :: umingmat	! allocated: micro_rates ; deallocated: calc_kinf
      real (8), ALLOCATABLE, DIMENSION(:,:) :: kejmat	! allocated: micro_rates ; deallocated: calc_kinf
      real (8), ALLOCATABLE, DIMENSION(:,:) :: ukejmat	! allocated: micro_rates ; deallocated: calc_kinf

      integer :: nbinmax			! max number of energy bins for array allocations (subr. micro_rates)

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
!	THIS MODULE CONTAINS THE FOLLOWING SUBROUTINES
!
!	subroutine micro_rates
!	subroutine calc_kinf
!	subroutine calei
!	subroutine calvo
!	subroutine cmm
!	subroutine find_ming
!	subroutine fx
!	subroutine ghrlev
!	subroutine jaread_input
!	subroutine javerage
!	subroutine jthermavg
!	subroutine krotlev
!	subroutine micro_writeout
!	subroutine morlev
!	subroutine odqhr
!	subroutine prewrite
!	subroutine shrlev
!	subroutine sterab
!	subroutine sterabj
!	subroutine uhrlev
!	subroutine vx
!	subroutine write_mat
!	subroutine zpe_calc
!-----------------------------------------------------------------------
 
      CONTAINS

      subroutine micro_rates
      use DECLARE_MOD				! use module DECLARE_MOD
      include 'DECLARE_INC.inc'		! include 'DECLARE_INC.inc'

      ! Declare variable with explicit types
      integer :: jval, njmax, status
      real*8 :: viblo, maxemax
      character :: fname*80, namefile*80

      ! Format statements
 100  format(i7,2x,100(es15.4))
 200  format(3x,'veff',2x,100(es15.4))
 990  FORMAT ('allocation FAILED: gmat, tmat')
 991  FORMAT ('allocation FAILED: rpmat, minpmat,mingmat,umingmat...')
 992  FORMAT ('micro_rates: deallocated gmat and tmat',/)
 993  FORMAT ('micro_rates: deallocation FAILED: gmat and tmat',/)

      ! Timestamps
      call sestamp('micro_rates',1)

      ! Find jmax for system
      call find_jmax(lunit,etopuser,jtopuser,jtop,rcnt,ntts,pcnt,maxtemp
     $,vmax,delh,epsil,barr,delhf,maxj,maxemax,maxj3)
     
      write(lunit,*) 'max energy for allocating arrays: ', int(maxemax)

      ! if maxj is larger than what the user defined then limit to user
      ! defined otherwise use found maxj as njmax
      IF( maxj .ge. jtopuser ) then
        njmax = jtopuser + 1
        else
        njmax = maxj + 1
        end if
      write(lunit,*) 'max J for allocating arrays: ', njmax
      write(lunit,*)
      
      ! find number of bins for system depending on max emax and bin
      ! size
      nbinmax = nint( maxemax/de ) + 1

      ! allocate tmat/gmat dynamically based on structures and number of
      ! bins
      ALLOCATE( tmat(rcnt+ntts,nbinmax) , 
     &          gmat(rcnt+ntts,nbinmax) , STAT=status )

      ! check if allocation succeded
      allocate_ok: IF ( status .NE. 0 ) THEN   
         WRITE (lunit,990)
         STOP 'terminated in micro_rates'
      END IF allocate_ok
      
c      ALLOCATE(    rpmat(njmax+1,nbinmax+1) , 
c     &           minpmat(njmax+1,nbinmax+1) , 
c     &           mingmat(njmax+1,nbinmax+1) , 
c     &          umingmat(njmax+1,nbinmax+1) , 
c     &            kejmat(njmax+1,nbinmax+1) , 
c     &           ukejmat(njmax+1,nbinmax+1) , STAT=status )

      ! allocat matrixes to hold reactant DOS, SOS, kej, and unified variants
      ALLOCATE(    rpmat(maxj+1,nbinmax+1) , 
     &           minpmat(maxj+1,nbinmax+1) , 
     &           mingmat(maxj+1,nbinmax+1) , 
     &          umingmat(maxj+1,nbinmax+1) , 
     &            kejmat(maxj+1,nbinmax+1) , 
     &           ukejmat(maxj+1,nbinmax+1) , STAT=status )

      ! check if allocation succeded
      allocate2_ok: IF ( status .NE. 0 ) THEN   
         WRITE (lunit,991)
         STOP 'terminated in micro_rates'
      END IF allocate2_ok
      
      ! do an initial run to generate the t(r,i) matrix
      jval=0
      call veff(jval,vef)
      do i=1,rcnt+ntts
         call sterab(i,jval)
      end do

      manymins=.false.


      !run over all j values

      write(*,*)'Calculating Microcanonical Rate Constants'
      do jval=0,maxj,dj
         write(*,FMT="(A1,A,t21,F6.2,A)",ADVANCE="NO") achar(13),
     $ "  Percent Complete: ",
     $ ((real(jval)+1.0E-10)/(real(maxj)+1.0E-10))*100.0,
     $ "%"
         call veff(jval,vef)
         do i=1,rcnt+ntts
            call sterabj(i,jval)
         end do
      

      if(.false.)then
         if(
     $       (jval.eq.0)    .or.
     $       (jval.eq.maxj) .or.
     $       (mod(jval,jtopuser)).eq.0
     $     ) then
            if(fileroot.eq.'')then
               write(fname,"(a,i3.3)")'gmat-j-',jval
            else
               namefile=trim(fileroot)//"-j-"
               write(fname,"(a,i3.3)")trim(namefile),jval
            end if
            namefile=trim(fname)//".gmat"
            open(unit=12, file=namefile)
            call stamp(12,2)
            call dnt(12)
            write(12,*)
            write(12,*)
            write(12,200)(vef(j),j=rcnt+1,rcnt+ntts)
            do j=1, binmax
               write(12,100)j-1,(gmat(i,j),i=rcnt+1,rcnt+ntts)
            end do
            close(12)
 
         end if
      end if
         call find_ming(binmax,jval)
      end do
      write(*,*)

      call javerage(viblo)

      write_files:  IF ( whatdo .eq. 'savefiles' ) THEN
                       call micro_writeout(viblo)
                    END IF write_files

      call calc_kinf
      
c
c release dynamically allocated memory     ******** DEALLOCATE ************
c
      DEALLOCATE( gmat, tmat , STAT=status)
      deallocate_ok: IF ( status .NE. 0 ) THEN   
         WRITE (lunit,993)
         STOP 'terminated in subroutine micro_rates'
      END IF deallocate_ok      

      call sestamp('micro_rates',2)

      return

      end subroutine micro_rates

!-----------------------------------------------------------------------

      subroutine calc_kinf
      use DECLARE_MOD				! use module DECLARE_MOD
      include 'DECLARE_INC.inc'		! include 'DECLARE_INC.inc'
      real*8  qr(nt),qt(nt),qu(nt)
      real*8 t(maxj+1,maxval(nbinl))
      integer x,y
      integer status


 100  format(a8,2x,a8  ,2x,2(a12,2x))
 200  format(i8,2x,f8.2,2x,2(es12.5,2x))
 201  format(a8,2x,a8,2x,2(a12,2x))
 300  format(f8.2,2x,3(es12.5,2x))
 301  format(a8,2x,3(a12,2x))

      write(lunit,*)'convergence for reactant'
      x=maxj+1
      y=maxval(nbinl)
      do j=1,x
         do k=1,y
            t(j,k)=rpmat(j,k)
         end do
      end do
      call nkinfint(temps,nt,t,x,y,v1l,sym(1),de,lunit,qr)

      write(lunit,*)      
      write(lunit,*)'convergence for transition state'
      y=binmax
      do j=1,x
         do k=1,y
            t(j,k)=minpmat(j,k)*kejmat(j,k)
         end do
      end do
      call nkinfint(temps,nt,t,x,y,vefmax,sym(2),de,lunit,qt)

      if(manymins)then
         write(lunit,*)      
         write(lunit,*)'convergence for unified transition state'
         do j=1,x
            do k=1,y
               t(j,k)=minpmat(j,k)*ukejmat(j,k)
            end do
         end do
         call nkinfint(temps,nt,t,x,y,vefmax,sym(2),de,lunit,qu)
      end if

      write(lunit,*)
      write(lunit,*)'thermally averaged microcanonical partition functio
     $ns'
      write(lunit,*)
      if(manymins)then
      write(lunit,301)'temp','q(r)','q(ts)','q(uts)'
      write(lunit,*)('-',i=1,51)
      do i=1,nt
         write(lunit,300)temps(i),qr(i),qt(i)*(h/(kbj*temps(i))),qu(i)*
     $(h/(kbj*temps(i)))
      end do

      else

      write(lunit,301)'temp','q(r)','q(ts)'
      write(lunit,*)('-',i=1,51)
      do i=1,nt
         write(lunit,300)temps(i),qr(i),qt(i)*(h/(kbj*temps(i)))
      end do
      end if

      do i=1,nt
        ratekinf(1,i)=qt(i)/qr(i)
        mrate(i)=ratekinf(1,i)
        ratekinf(2,i)=qu(i)/qr(i)
       umrate(i)=ratekinf(2,i)
      end do

      if(manymins)then
         write(lunit,*)
         write(lunit,*)'thermal averaged microcanonical rate constants'
         write(lunit,*)
         write(lunit,201)'#','temp','<kv(e)>t','<ukv(e)>t'
         write(lunit,*)('-',i=1,47)
         do i=1,nt
            write(lunit,200)i,temps(i),ratekinf(1,i),ratekinf(2,i)
         end do
      else
         write(lunit,*)
         write(lunit,*)'thermal averaged microcanonical rate constants'
         write(lunit,*)
         write(lunit,201)'#','temp','<kv(e)>t'
         write(lunit,*)('-',i=1,47)
         do i=1,nt
            write(lunit,200)i,temps(i),ratekinf(1,i)
         end do
      end if

c
c release dynamically allocated memory     ******** DEALLOCATE ************
c
      DEALLOCATE ( rpmat , 
     &             minpmat , 
     &             mingmat , 
     &             umingmat , 
     &             kejmat , 
     &             ukejmat , STAT=status )
      deallocate_ok: IF ( status .NE. 0 ) THEN   
         WRITE (lunit,993)
         STOP 'terminated in subroutine calc_kinf'
      END IF deallocate_ok
993   FORMAT ('calc_kinf: deallocation FAILED: rpmat ... ukemat',/)

      return

      end subroutine calc_kinf

!-----------------------------------------------------------------------

      subroutine calei(ev,nn,nmax,emax,b)
 
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c calei: a code to compute/select eigenvalues of torsional motion based on emax 
c copyright (c) 2009 lam t. nguyen and john r. barker
c
c john r. barker
c jrbarker@umich.edu
c university of michigan
c ann arbor, mi 48109-2143
c (734) 763 6239
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c     input:      1) ev(nmax) : a vector of eigenvalues with nmax elements
c           2) nmax           : number of elements of vector ev
c           3) emax           : maximum energy in master equation
c           4) b        : average rotational constant, defined as:
c                 b = (integral of b(x)dx from 0 to 2*pi) / (2*pi)
c     output:     1) ev(nmax)
c           2) nmax     
c
c
      implicit double precision(a-h,o-z), integer(i-n)
      dimension ev(nn)

      i=nmax

c.......if the largest eigenvalue is larger than the maximum energy, 
c.......then searching for an eigenvalue, which is <= emax
11    e=ev(i)
      if(e.gt.emax) then
            i=i-1
            goto 11
      endif
      nmax=i

c.......looking for higher-lying eigenvalues with double degeneracy
111   ne1=nint(ev(i))
      ne2=nint(ev(i-1))
      if(ne1.gt.ne2) then
            i=i-1
            goto 111
      endif

c.......switching from a hindered rotor to a pitzer rotor
c.......pitzer rotor is defined as e = eo + b*(jp**2)
c.......two variables eo and jp will be computed as below
      tp=ev(i-1) - ev(i-2)
      jp=nint((tp/b - 1.0d0)/2.0d0)
      eo=ev(i-2)-b*jp*jp

1111  ev(i+1)=eo + b*(jp+2)*(jp+2)
      ev(i+2)=ev(i+1)
      if(ev(i+2).le.emax) then
            i=i+2
            jp=jp+1
            goto 1111
      endif
      nmax=i

      return
 
      end subroutine calei

!-----------------------------------------------------------------------

c
c calvo is used to compute vmax and vmin 
c of the torsional potential energy function 
c
        subroutine calvo(vo,vm,cv,ncv,nsiv,vhr,pha)

        implicit double precision(a-h,o-z), integer(i-n)
        dimension cv(ncv)
        character vhr*5

        pi=acos(-1.0d0)

        rmin=0.0d0
        rmax=pi*2
        zl=(rmax-rmin)
        dx=zl/3600.0d0
        call vx(vo,0.0d0+pha,cv,ncv,nsiv,vhr,0.0d0)
        vm=vo
        do i=1, 3600
                xa = rmin + dx*i
                call vx(vtp,xa+pha,cv,ncv,nsiv,vhr,0.0d0)
                if(vtp.gt.vo) vo=vtp
                if(vtp.lt.vm) vm=vtp
        enddo
        return
        end subroutine calvo
        
!-----------------------------------------------------------------------

c.....      cmm stands for constructing meyer's matrix, ar, for 1d-hindered internal rotor using meyer's method, 
c.....      r. meyer j. chem. phys. 52 (1970) 2053.
c.....      both rotational energy profile and moment inertia are properly treated as functions of 
c.....      internal rotation angle, i.e. non-rigid rotor
c.....      input:      1) ncf, ncv = number of coefficients drived from fitting rotational energy profile and moment inertia
c.....            2) cv(ncv) = vector of coefficients from rotational energy profile, with ncv elements
c.....            3) cf(ncf) = vector of coefficients from moment inertia, with ncf elements
c.....            4) nsig = rotational symmetry number, it can be either one or two or three or whatever 
c.....      output: 1) ar(nx,nx) = meyer's symmetrical matrix, with a size of nx--number of grid points   


      subroutine cmm(ar,ncv,ncf,cv,cf,nsiv,nsif,vhr,bhr,phav,phab,vm)
      implicit double precision(a-h,o-z), integer(i-n)
      parameter (nx=501)
c
c.....declare arrays.
      dimension ar(nx,nx),xa(nx) 
      dimension d(nx,nx), f(nx), v(nx), dt(nx,nx),
     & t(nx,nx), ft(nx,nx)
      dimension  cf(ncf), cv(ncv)
      character vhr*5, bhr*5

c
c.....variable data input.

      pi=acos(-1.d0)
c
c.....set up grid
      rmin=0.0d0
      rmax=pi*2
      zl=(rmax-rmin)
      dx=zl/nx

c.....now compute meyer's symmetrical matrix:

      do i=1, nx
            xa(i)=rmin+dx*i
            call vx(v(i),xa(i)+phav,cv,ncv,nsiv,vhr,vm)
            call fx(f(i),xa(i)+phab,cf,ncf,nsif,bhr)
            do j=1, nx
                  imj=i-j
                  if(imj.eq.0) then
                        d(i,j)=0.0d0
                  else
                        tp1=pi*imj/dble(nx)
                        tp2=(-1.0d0)**imj
                        d(i,j)=tp2/2.0d0/dsin(tp1)
                  endif
                  dt(i,j)=-d(i,j)
                  ft(i,j)=f(i)*delta(i,j)
            enddo
      enddo

      do i=1, nx
            do j=1, nx
                  su=0.0d0
                  do k=1, nx
                        su=su + ft(i,k)*d(k,j)
                  enddo
                  t(i,j)=su
            enddo
      enddo

      do i=1, nx
            do j=1, nx
                  su=v(i)*delta(i,j)
                  do k=1, nx
                        su=su + dt(i,k)*t(k,j)
                  enddo
                  ar(i,j)=su
            enddo
      enddo 

      return 
      end subroutine cmm

!-----------------------------------------------------------------------

      subroutine find_ming(n,jval)
      use DECLARE_MOD				! use module DECLARE_MOD
      include 'DECLARE_INC.inc'		! include 'DECLARE_INC.inc'
      integer jval,fend,minct,minil(rcnt+ntts)
c      real*8  gmat(rcnt+ntts,n),mingl(rcnt+ntts)
      real*8  mingl(rcnt+ntts)
      real*8  x(3),y(3),ming
      real*8  tdat(rcnt+ntts),urate(n)
      integer tmins(rcnt+ntts),fct(n)
      integer tempx,rx,lx
      character frmt*60,namefile*80
 100  format(i3,2x,f10.2,2x,i5)
 200  format(2x,a1,6x,a6,2x,a5,2x,4(a10,2x))
 300  format(i5,2x,4(es15.5,2x))

c
c   if there are 2 or more trial ts open nmins file
c      if(ntts.ge.2)then         
c         if(inputfile.eq.'')then
c            namefile='nmins.txt'
c         else
c            namefile=trim(fileroot)//"-nmins.txt"
c         end if
c
c         open(unit=nminsunit,file=namefile)
c         if(jval.eq.0)then
c            call stamp(nminsunit,2)
c            call dnt(nminsunit)
c            write(nminsunit,*)'run id: ',tstmp
c            write(nminsunit,*)
c            write(nminsunit,*)
c            write(nminsunit,200)'j','e cm-1','# min','found@','min g'
c         end if
c      end if

c
c   check for number of trial ts
      if(ntts.le.1)then                                  ! only 1 ts structure given
         do i=1,n                                        ! loop over all energy grains
            mingmat(jval+1,i)=gmat(rcnt+ntts,i)
           umingmat(jval+1,i)=gmat(rcnt+ntts,i)
         end do

      else                                               ! multiple ts given
c         write(*,*)'jval',jval 
         do i=1,n                                        ! loop over all energy grains

            fct(i)=0                                     ! found min count

            do j=rcnt+1,rcnt+ntts                        ! only over ntts
               tdat(j-rcnt)=gmat(j,i)
            end do
           
c            write(*,*)jval,i 
            call find_tmins(tdat,tmins,ntts,fct(i),microh)


            if(fct(i).gt.1)manymins=.true.            

            urate(i)=0.0d0
            ming=gmat(rcnt+tmins(1),i)
            do j=1,fct(i)
               minil(j)=rcnt+tmins(j)
               mingl(j)=gmat(rcnt+tmins(j),i)
               urate(i)=urate(i)+(1.0d0/mingl(j))
               if(mingl(j).lt.ming)ming=mingl(j)             
            end do
            urate(i)=1.0d0/urate(i)
           
c            if(ntts.ge.2)then 
c            write(frmt,'("(i3,2x,f10.2,2x,i3,2x,",i0,"(i3,2x),",i0,
c     $           "(es10.2,2x),i3)")')fct(i),fct(i)
c            write(nminsunit,frmt)jval,(i-1)*de,fct(i),
c     $           (minil(l)-rcnt,l=1,fct(i)),(mingl(l),l=1,fct(i))
c            end if

            mingmat(jval+1,i)=ming
           umingmat(jval+1,i)=urate(i)            


         end do

      end if

      return
      end subroutine find_ming
      
!-----------------------------------------------------------------------
      
      subroutine fx(f,x,cf,n,nsig,bhr)
      implicit double precision(a-h,o-z), integer(i-n)
      dimension cf(n)
      parameter (con=16.85763d0)
      character bhr*5

      if((bhr.eq.'bhrd1').or.(bhr.eq.'bhrd1')) then
            f=cf(1)
            do i=2, n
                  f=f+cf(i)*dcos(x*nsig*(i-1))
            enddo
        elseif((bhr.eq.'ihrd1').or.(bhr.eq.'ihrd1')) then
                f=cf(1)
                do i=2, n
                        f=f+cf(i)*dcos(x*nsig*(i-1))
                enddo
            f=con/f
      else
                write(*,*) "error at input for bhrd or ihrd"
            stop
      endif

      return
      end subroutine fx

!-----------------------------------------------------------------------

      subroutine ghrlev(ev,nn,nmax,emax,b,ncv,ncf,cv,cf,
     &            nsiv,nsif,vhr,bhr,phav,phab)
 
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c ghrlev: a code of general purposes to calculate eigenvalues of 1d-hindered internal rotation 
c copyright (c) 2009 lam t. nguyen
c
c     date: feb. 13, 2009
c
c lam t. nguyen
c nguyenlt@umich.edu
c university of michigan
c ann arbor, mi 48109-2143
c
c or contact
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
c see the "readme" file for a copy of the gnu general public license,
c or contact:
c
c free software foundation, inc.
c 59 temple place - suite 330
c boston, ma 02111-1307, usa.
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c     input:      
c           1) emax           : maximum energy in master equation
c           2) b        : average rotational constant, defined as:
c                 b = (integral of b(x)dx from 0 to 2*pi) / (2*pi)
c           3) ncv, ncf : number of elements in vectors of cv and cf, respectively
c           4) cv       : vector of coefficients with ncv elements in torsional 
c                 potential energy function, which is expressed as = cvo/2 + sum of cvi*(1-cos(i*x*nsig))/2
c           5) cf       : vector of coefficients with ncf elements in rotational constant,
c                 which is given by a fourier series = cfo + sum of cfi*cos(i*x*nsig)
c           6) nsiv, nsif     : symmetry numbers for tortional potential e function and rotational constant, respectively.
c           7) nn       : maximum number of elements in vector ev
c           8) vhr            : model for torsional potential energy function
c           9) bhr            : model for rotational constant
c           10) phav, phab    : phases of torsional angles in vhr and bhr, respectively.         
c
c     output:     1) ev(nmax) : vector of eigenvalues with nmax elements
c           2) nmax           : number of elements of vector ev
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


      implicit double precision(a-h,o-z), integer(i-n)
      parameter (nx=501)
c
c.....declare arrays.
      dimension ev(nn), cv(ncv), cf(ncf)
        character vhr*5, bhr*5
c
c.....calvo used to find vmax (vo) and vmin (vm) of torsional potential energy function
      call calvo(vo,vm,cv,ncv,nsiv,vhr,phav)

      call odqhr(ev,nn,ncv,ncf,cv,cf,nsiv,nsif,vhr,bhr,phav,phab,vm)
      nmax=nx
      call calei(ev,nn,nmax,emax,b)

        zpe = ev(1)
c        write(*,9991) vo, vm, ev(1)
9991    format(11x,'torsional potential: ',1x,'vmax (cm-1) = ',f8.1,1x,
     &  '; vmin (cm-1) = ',f8.1,1x,'; zpe (cm-1) = ',f5.1)
     
c*********************************************************************
c     the following lines are for output of hindered rotor eigenvalues
c*********************************************************************
c        write(12,*) '   '
c        write(12,9991) vo, vm, zpe

c        do i=1, imax+1
c          write(12, 99) i, ev(i), (ev(i)-zpe)
c        enddo
c99      format (3x, i5, 1x, f10.3 , 1x, f10.3 )
c*********************************************************************
     
      return 
      end subroutine ghrlev

!-----------------------------------------------------------------------

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
      subroutine jaread_input(inputfile,imax1,isize,temax2)

c	read output files, do J-summing, and create 1D ".dens' files for multiwell

      implicit none
      character(*) inputfile
      character filetype*4,fname*80
      real*8,allocatable,dimension(:,:)::mat
      real*8,allocatable,dimension(:)::khsg,sumdens
      real*8 emax,emax1,emax2,de,rotjay(1001),e,egrain2,egrain1,temax2
      real*8 viblo
      character(60000)line
      character(15)line2,junk,arg1,arg2,arg3
      character   sumtype
      logical there
      integer i,j,jdex,nbin,tnbin
      integer iunit,ounit,ccount,imax1,isize
      integer jays(1001),maxj,eol1,eol2,reason,lnct

 100  format('*** ',A,' not found ***')
 200  format('Unsupported file. Please use a .2dens or .2sums file')
 300  format(i10,1x,f10.1,2(1x,1pe12.5))
 400  format(5(F10.2,2x)) 
 500  format(1x,A,2x,F10.2)
 600  format(1x,A,2x,I10)

      iunit=111
      ounit=222

c
c  Get name of file to open

c      call get_command_argument(1,inputfile)
c      ccount=command_argument_count()
      ccount=10
      write(*,*)
      write(*,*)"File:",trim(inputfile)
      filetype=inputfile(len_trim(inputfile)-3:len_trim(inputfile))
      if(filetype.eq.'sums')then
         sumtype=inputfile(len_trim(inputfile)-5:len_trim(inputfile)-4)
      end if

c
c  Grab command line arguments for imax and isize
c      if(ccount.gt.1)then
c         call get_command_argument(2,arg1)
c         call get_command_argument(3,arg2)
c         read(arg1,*)imax1
c         read(arg2,*)isize
         write(*,600)'imax1:',imax1
         write(*,600)'isize:',isize
c         if(ccount.gt.3)then
c            call get_command_argument(4,arg3)
c            read(arg3,*)temax2
            write(*,500)'emax2:',temax2
c         end if
c      end if
      

c
c  Check if file is present
 111  inquire(file=inputfile,exist=there)
      if(.not.there)then
          write(*,100)trim(inputfile)
          stop
      end if

c
c  Check if valid file type

      if(trim(filetype).eq.'dens')then
         fname=inputfile(1:len_trim(inputfile)-6)//"-r."//trim(filetype)
      else if(trim(filetype).eq.'sums')then
            if(sumtype.eq.'u')then
               fname=inputfile(1:len_trim(inputfile)-7)//"-uts.dens"
            else
               fname=inputfile(1:len_trim(inputfile)-6)//"-ts.dens"
            end if
      else
           write(*,200)
           stop
      end if
      
c
c  Find max J value present

      open(unit=iunit,file=inputfile)
      lnct=0
      read(iunit,*)junk
      read(iunit,*)junk, maxj, emax, de
      read(iunit,*)viblo
      lnct=lnct+1

      read(iunit,"(A)",IOSTAT=eol1)line
      read(line,*,IOSTAT=eol2)jays

c  Read in rotational energies
      read(iunit,"(A)",IOSTAT=eol1)line
      read(line,*,IOSTAT=eol2)junk,rotjay
      lnct=lnct+1

c
c  Skip blank line

C      read(iunit,*)     
C      lnct=lnct+1
C      read(iunit,*,iostat=reason)emax1
c      lnct=lnct+1
c      read(iunit,*,iostat=reason)emax2
c      lnct=lnct+1
c      de=emax2-emax1
c      write(*,*)emax1,emax2,de

c
c  Find max E present
c      do
c        read(iunit,*,iostat=reason)emax
c        lnct=lnct+1
c        if(reason.lt.0)goto 700
c      end do
 700  close(iunit)

 
c
c  Set values for nbin, egrain2
      tnbin=nint(emax/de)+1
      if(temax2.ne.emax)then
         nbin=int(temax2/de)+1
         Egrain2=temax2/(isize-imax1-1)
      else
        nbin=nint(emax/de)+1
        Egrain2=emax/(isize-imax1-1)
      end if

c
c  Allocate array space
      allocate(mat(maxj,tnbin))
      allocate(khsg(tnbin))
      allocate(sumdens(tnbin))

c
c  Re-openfile

      open(unit=iunit,file=inputfile)

c
c  skip 4 lines

      do i=1,4
         read(iunit,*)
      end do

c
c  read in full matrix
      do i=1,tnbin
         read(iunit,*)line,(mat(j,i),j=1,maxj)
         khsg(i)=0.0D0
      end do
      close(iunit)

c
c  do j-summing placing all mat(j,i) into khsg(i+jdex)

      do j=1,maxj
         jdex=nint((rotjay(j)-rotjay(1))/de)
         do i=1,tnbin
            if((i+jdex).le.tnbin)then
                khsg(i+jdex)=khsg(i+jdex)+mat(j,i)
            end if
         end do
      end do
       

c      if(ccount.gt.1)then
         call chkdens(khsg,tnbin,dE,imax1,tnbin,de*(tnbin-1))
c      end if

c
c  calculate sums or dens
      if(trim(filetype).eq.'dens')then
         do i=1,tnbin
            if(i.eq.1)then
                sumdens(i)=khsg(i)*de
            else
                sumdens(i)=sumdens(i-1)+khsg(i)*de
            end if
         end do
      else
         do i=1,tnbin
            if(i.gt.1)then
                sumdens(i)=(khsg(i)-khsg(i-1))/de
            else
                sumdens(i)=khsg(i)/de
            end if
         end do
      end if

c
c write out .dens file
c data order is:
c Bin Number, Energy of Bin, Density of Bin, Sum of Bins.      
      open(unit=ounit,file=fname)
      call japrewrite(ounit,fname,de,imax1,temax2,isize,viblo)
      if(ccount.lt.2)then
         if(trim(filetype).eq.'dens')then
            do i=1,tnbin
               write(ounit,300)i,de*(i-1),khsg(i),sumdens(i)
            end do
         else
            do i=1,tnbin
               write(ounit,300)i,de*(i-1),sumdens(i),khsg(i)
            end do
         end if
      else
         if(trim(filetype).eq.'dens')then
            do i=1,imax1
               E = (i-1)*de
               if(i.gt.nbin)then
                  j=nbin
                  write(ounit,300)i,e,khsg(j),sumdens(j)
               else
                  write(ounit,300)i,e,khsg(i),sumdens(i)
               end if
            end do
   
c
c  do j-summing placing all mat(j,i) into khsg(i+jdex)
         do i=imax1+1, isize
               E = (i-imax1-1)*Egrain2
               j = int(e/de) + 1
               if(j.gt.tnbin)then
                  j=tnbin
                  write(ounit,300)i,e,khsg(j),sumdens(j)
               else
                  write(ounit,300)i,e,khsg(j),sumdens(j)
               end if
            end do
         else
            do i=1,imax1
               E = (i-1)*de
               if(i.gt.nbin)then
                  j=nbin
                  write(ounit,300)i,e,sumdens(j),khsg(j)
               else
                  write(ounit,300)i,e,sumdens(i),khsg(i)
               end if
            end do
            do i=imax1+1, isize
               E = (i-imax1-1)*Egrain2
               j = int(e/de) + 1
               if(j.gt.tnbin)then
                  j=tnbin
                  write(ounit,300)i,e,sumdens(j),khsg(j)
               else
                  write(ounit,300)i,e,sumdens(j),khsg(j)
               end if
            end do
         end if
      end if
 
      deallocate(sumdens)
      deallocate(khsg)
      deallocate(mat)

      end subroutine jaread_input

!-----------------------------------------------------------------------

      subroutine javerage(viblo)

      use DECLARE_MOD				! use module DECLARE_MOD
      include 'DECLARE_INC.inc'		! include 'DECLARE_INC.inc'
      integer jdex
c      real*8    gmat(rcnt+ntts,epsmax),khsp(5*etop),totp,viblo
c      real*8    khsp(5*etopuser),totp,viblo
c      real*8    khsg(5*etopuser),ikhsg(5*etopuser)
c      real*8    ukhsg(5*etopuser),uikhsg(5*etopuser)
      real*8    khsp(2*nbinmax),totp,viblo
      real*8    khsg(2*nbinmax),ikhsg(2*nbinmax)
      real*8    ukhsg(2*nbinmax),uikhsg(2*nbinmax)
 100  format(i7,2x,30(es15.4))
 200  format("**************input data summary**************")
 201  format(a,' data as generated from ',a)
 202  format(f10.2,2x,i10,2x,f10.2,2x,i10,2x,f10.2)
 203  format(a6,2x,a10,2x,5(a15))
 204  format(i6,2x,f10.2,2x,5(es15.6,2x))



c
c   calculate rate constants from sums and densities of states

      if(rcnt+ntts.eq.1)then                                            !  if only 1 species is passed write out .dens file for that species

         do j=1,maxj+1                                                  !  klippenstein, harding, smith, & gilbert j averaging 
            do i=1,maxval(nbinl)
               jdex=nint((barr(1)*j*(j-1))/de)
               khsp(i+jdex)=khsp(i+jdex)+rpmat(j,i)
            end do
         end do

         call get_viblo(rcnt+ntts,rcnt+ntts,viblo)

      else if(rcnt+ntts.ge.2)then                                       !  if r and ts passed write out .dens file for each species

         do j=1,maxj+1                                                  !  klippenstein, harding, smith, & gilbert j averaging 
            do i=1,maxval(nbinl)
               jdex=nint((barr(1)*j*(j-1))/de)
               khsp(i+jdex)=khsp(i+jdex)+rpmat(j,i)
               jdex=nint((vefmax(j)-vefmax(1))/de)
               khsg(i+jdex)=   khsg(i+jdex) +   mingmat(j,i)
              ukhsg(i+jdex)=  ukhsg(i+jdex) +  umingmat(j,i)
                 
            end do
         end do

         call get_viblo(rcnt,rcnt,viblo)


         do j=1,maxj+1,dj
            do i=1,binmax
                 kejmat(j,i)=(   mingmat(j,i)/minpmat(j,i))*c
                ukejmat(j,i)=(  umingmat(j,i)/minpmat(j,i))*c
            end do
         end do

         call jthermavg

      end if


      end subroutine javerage

!-----------------------------------------------------------------------

      subroutine jthermavg
      use DECLARE_MOD				! use module DECLARE_MOD
      include 'DECLARE_INC.inc'		! include 'DECLARE_INC.inc'
      real*8  nusum,desum,term
      real*8  jpmat(binmax,nt),jgmat(binmax,nt)
      real*8  jkmat(binmax,nt)
 100  format(i5,2x,f10.2,2x,100(ES15.6,2x))
 200  format(19x,100(F15.2,2x))

      do k=1,nt
         do i=1,binmax
            desum=0.0d0
            do j=1, maxj+1
               term=Barr(1)*j*(j+1)
               term=-1.0d0*term
               term=term/(kb*Temps(k))
               term=Dexp(term)
               term=minpmat(j,i)*term
               desum=desum+term
            end do
            jpmat(i,k)=desum
         end do

         do i=1,binmax
            nusum=0.0d0
            do j=1, maxj+1
               term=Barr(1)*j*(j+1)
               term=-1.0d0*term
               term=term/(kb*Temps(k))
               term=Dexp(term)
               term=mingmat(j,i)*term*c
               nusum=nusum+term
            end do
            jgmat(i,k)=nusum
         end do

         do i=1,binmax
            if(i.LT.(binmax-binmax))then
               jkmat(i,k)=0.0d0
            else
               jkmat(i,k)=(jgmat(i,k)/jpmat(i,k))
            end if
         end do
      end do

c j-averaged k(E) for each temperature
c      open(unit=junit,file="javge.txt")
c      write(junit,200)(temps(k),k=1,nt)
c      do i=1,binmax
c         write(junit,100)i,de*(i-1),(jkmat(i,k),k=1,nt)
c      end do
c      close(junit)

      end subroutine jthermavg
      
!-----------------------------------------------------------------------

      subroutine krotlev(at,t,nmax,dele,jmax,b,a,jk)
 
      implicit double precision(a-h,o-z)
      implicit integer(i-n)
      double precision at(nmax) , t(nmax)
      integer ixx
      save 
 
      rmax = sqrt(jmax*dele/b)           ! highest rotational quantum number consistent with energy nint
      imax = int(rmax) + 1               ! integer version
      if ( imax .gt. jk ) imax = jk
      if ( imax .gt. 20001 ) imax = 20001
             
      do j = 1, imax					! k = quantum number of 1-dimensional free rotor
         r = b*j*j+a*jk*(jk+1)			! eigenvalues
         r = b*j*j						! eigenvalues
         ir = nint( r/dele ) + 1		! nearest integer: number of grains
         
         if (j .eq. 0 ) then
             f = 1.0                       ! f is the multiplicity of the 1-d rotor (e.g. single-axis internal rotor) energy level (see j.l. mchale, molecular spectroscopy (prentice-hall, 1999), 216f.
         else
             f = 2.0                       ! f is the multiplicity of the 1-d rotor (e.g. single-axis internal rotor) energy level (see j.l. mchale, molecular spectroscopy (prentice-hall, 1999), 216f. 
         endif
c         write(12,*)j,r,ir,b,f

         do k = ir, jmax
            karg = k - ir + 1
            if (karg.le.jmax) at(k) = at(k) + f*t(karg)
         end do  ! k

      end do  ! j

      do j = 1 , jmax
          t(j) = at(j) + t(j)
          at(j) = 0.0d+00
      end do  ! j 

c      write(12,*)
 
      return
      end subroutine krotlev

!-----------------------------------------------------------------------

      subroutine micro_writeout(viblo)
      use DECLARE_MOD				! use module DECLARE_MOD
      include 'DECLARE_INC.inc'		! include 'DECLARE_INC.inc'
      real*8  viblo
      
      call prewrite('pmat',viblo)
      call prewrite('ming',viblo)
      if(rcnt+ntts.gt.1)then
         call prewrite('kej',viblo)
      end if

      if(manymins)then
         call prewrite('uming',viblo)
         call prewrite('ukej',viblo)
      end if


      end subroutine micro_writeout

!-----------------------------------------------------------------------

      subroutine morlev(at,t,nmax,dele,jmax,we,xe,ng,zpe)
      implicit none
      double precision at , t , dele
      dimension at(nmax) , t(nmax)
      integer i , ir , nmax , jmax , ng
      dimension ir(20001)
      double precision we , xe , zpe, w0 , rmax , vi , r
      integer imax , j , k , karg , l
      external sterab , gam , wrab , wrden , alngam
      save 
 
      zpe = 0.5d+00*we + 0.25d+00*xe                    ! zero point energy at v=0

      if ( xe .lt. 0.0 ) then
         w0 = we + xe
         rmax = -0.5d+00*w0/xe                          ! highest bound state, e relative to zpe
      else
         rmax = jmax*dele/we                            ! highest bound state
      endif
      imax = int(rmax) 
      if ( imax.gt.20001 ) imax = 20001

      do i = 1 , imax                               ! start at v=1
        vi = i + 0.5d+00
        r = (we + xe*vi)*vi - zpe                   ! state energy relative to zpe
c        ir(i) = 2 + int( r/dele )                   ! index of bin above zpe
        ir(i) = 1 + nint( r/dele )                ! nearest integer number of grains
        ir(i) = 1 + nint( r/dele )                ! nearest integer number of grains
      end do ! i

      do l = 1 , ng                                    ! loop for degenerate vibrations
         do j = 1 , imax                                  ! imax is the number of energy states
            do k = 1 , jmax                                  ! jmax is the number of energy grains
               karg = k - ir(j) + 1
               if ( karg.gt.0 .and. karg.le.jmax ) then
                 at(k) = at(k) + t(karg)               ! at(k) is the number of states in the kth grain
               endif
            end do  ! k
         end do  ! j
 
         do j = 1 , jmax
            t(j) = at(j) + t(j)
            at(j) = 0.0d+00
         end do
      end do   ! l
 
      return
      end subroutine morlev

!-----------------------------------------------------------------------

c****************************************************************************************************************
c     name: subroutine odqhr 
c     use:  to obtain a vector of quantum eigenvalues of 1d-hindered internal rotor, in which  
c           both torsional potential energy function and moment inertia are explicitly treated as 
c           functions of internal rotational angle, i.e. non-rigid rotor, based on meyer's algorithm:
c           r. meyer j. chem. phys. 52 (1970) 2053-2059.  
c
c     authors:    lam t. nguyen and john r. barker 
c
c     date: jan. 26, 2009
c
c****************************************************************************************************************
c     
c      nx  : number of grid points
c        xa  : position on a grid.
c        ar  : hamiltonian matrix, i.e. meyer's symmetrical matrix
c        er  : vector of eigenvalues 
c        zr  : eigenvectors (x,y); where
c                                      x : wavefunction.
c                                      y : energy level.       
c        npr      : 0 for eigenvalue only and 1 for both eigenvalue and eigenvector 
c        rmin     : starting point of grid
c        rmax     : end point of grid
c        zl       : grid length
c        dx       : grid spacings
c
c****************************************************************************************************************
c
c     input:      1) ncv, ncf = number of coefficients derived from fitting rotational energy profile and moment inertia
c           2) cv(ncv)  = vector of coefficients from rotational energy profile, with ncv elements
c           3) cf(ncf)  = vector of coefficients from moment inertia (rotational constant), with ncf elements
c           4) nsig           = rotational symmetry number
c     output:     1) er             = vector of eigenvalues of 1d-hindered internal rotor, with nx=501 elements 
c
c*****************************************************************************************************************

      subroutine odqhr(er,nn,ncv,ncf,cv,cf,nsiv,nsif,
     &      vhr,bhr,phav,phab,vm)
      implicit double precision(a-h,o-z), integer(i-n)
      parameter (nx=501)
c
c.....declare arrays.
      dimension ar(nx,nx),er(nn),zr(nx,nx),xa(nx) 
      dimension d(nx,nx), f(nx), v(nx), dt(nx,nx),
     & t(nx,nx), ft(nx,nx)
      dimension  cf(ncf), cv(ncv), fv1(nx), fv2(nx)
        character vhr*5, bhr*5

c.....      dimension work(nx*3)

c
c.....variable data input.

      data npr/0/

c.....now compute hamiltonian matrix:

      call cmm(ar,ncv,ncf,cv,cf,nsiv,nsif,vhr,bhr,phav,phab,vm)
                              
c.....now call eigenvalue solver using rs subroutine of eispack.

      call rs(nx,nx,ar,er,npr,zr,fv1,fv2,ierr) 

c...... alterative option: using subroutine dsyev from lapack 
c......     call dsyev( 'v', 'u', nx, ar, nx, er, work, nx*3, ierr )
c......

      return 
      end subroutine odqhr

!-----------------------------------------------------------------------

      subroutine prewrite(mat,viblo)
      use DECLARE_MOD				! use module DECLARE_MOD
      include 'DECLARE_INC.inc'		! include 'DECLARE_INC.inc'
      character(len=*) mat
      character fnme*15
      real*8   t(maxj+1,maxval(nbinl)),viblo
      integer  x,y


      select case(mat)
         case ('pmat')
            x=maxj+1
            y=maxval(nbinl)
            do i=1,x
               do j=1,y
                  t(i,j)=rpmat(i,j)
               end do
            end do
            fnme='2dens'
         case ('minp')
            x=maxj+1
            y=binmax
            do i=1,x
               do j=1,y
                  t(i,j)=minpmat(i,j)
               end do
            end do
            fnme='crit2dens'
         case ('ming')
            x=maxj+1
            y=binmax
            do i=1,x
               do j=1,y
                  t(i,j)=mingmat(i,j)
               end do
            end do
            fnme='2sums'
         case ('uming')
            x=maxj+1
            y=binmax
            do i=1,x
               do j=1,y
                  t(i,j)=umingmat(i,j)
               end do
            end do
            fnme='u2sums'
         case ('kej')
            x=maxj+1
            y=binmax
            do i=1,x
               do j=1,y
                  t(i,j)=kejmat(i,j)
               end do
            end do
            fnme='kej'
         case ('ukej')
            x=maxj+1
            y=binmax
            do i=1,x
               do j=1,y
                  t(i,j)=ukejmat(i,j)
               end do
            end do
            fnme='ukej'
      end select

      fnme=trim(fnme)
      numnames=numnames+1
      call write_mat(t,x,y,de,dj,v1l,vefmax,inputfile,fnme,rcnt+ntts,
     $               viblo,tstmp,names(numnames))

      end subroutine prewrite

!-----------------------------------------------------------------------

c
c shrlev: a code to calculate eigenvalues of symmetrical, rigid 1d-hindered internal rotation
c copyright (c) 2009 lam t. nguyen and john r. barker
c
c	date:	feb. 13, 2009
c
c john r. barker
c jrbarker@umich.edu
c university of michigan
c ann arbor, mi 48109-2143
c (734) 763 6239
c
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
c see the "readme" file for a copy of the gnu general public license,
c or contact:
c
c free software foundation, inc.
c 59 temple place - suite 330
c boston, ma 02111-1307, usa.
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        subroutine shrlev(t,at,dele,jmax,b,v,n1,n2,imax,zpe,
     &          vhr,bhr,phav,phab)

        implicit double precision(a-h,o-z), integer(i-n)
        parameter (nn=2000, no=1)
        dimension ev(nn), cv(no), cf(no), t(jmax), at(jmax)
        dimension ir(nn)
        character vhr*5, bhr*5

        emax=(jmax-1)*dele + 1000.              ! *****raise max energy to account for later subtraction of zpe******

        ncv=1
        ncf=1
        cv(1)=v
        cf(1)=b

        call ghrlev(ev,nn,imax,emax,b,ncv,ncf,cv,cf,n1,n2,
     &          vhr,bhr,phav,phab)

        zpe=ev(1)

c*********************************************************************
c     the following lines are for output of hindered rotor eigenvalues
c*********************************************************************

        do i=1, imax+1
          tp=ev(i)-zpe
        enddo
99    format (3x, i5, 1x, f10.3 , 1x, f10.3 )
c*********************************************************************

      do i=1, imax
         tp=ev(i)-zpe
         ir(i)=nint( tp/dele )
         ir(i)=nint( tp/dele )
      enddo 

      do j=1, imax
         ll=ir(j)
         do l=1, jmax - ll
            at(l+ll) = at(l+ll) + t(l)
         enddo
      enddo
      do j=1, jmax
         t(j)=at(j)/n1
        at(j)=0.0d0
      enddo 

      return
      end subroutine shrlev

!-----------------------------------------------------------------------

      subroutine sterab(z,jval)			!	CALLED BY micro_rates     
      use DECLARE_MOD				! use module DECLARE_MOD
      include 'DECLARE_INC.inc'		! include 'DECLARE_INC.inc'
      integer z,imax,jarg,nbin,jval,passno
      integer r(ndof(z)),vcnt,hracnt,nrg
      real*8  summ,rmax
      real*8  emax,b,v(ndof(z)),vv
      real*8,allocatable,dimension(:)::t,at
      integer,allocatable,dimension(:)::ir
      real*8  bb,vvv,ngg,zap
      real*8  hral(ndof(z)),zzpe(ndof(z))
c      real*8  tmat(rcnt+ntts,etop)
      real*8  wee,anhh
      integer nng,cnt,anh,we
      integer :: status						! Status: 0 for success
      logical tt,ff
 
 100  format(i3,2x,f10.2,2x,i5) 

!	DEFINITIONS
!	z		: index number for reactants and trial ts, looped in subroutine micro_rates
!	nbin	: number of bins, corresponding to epsil(z)
!	tmat	: matrix of Beyer-Swinehart tmat(z,i) with emax = epsil(z)
!	jval	: value of J 
!

c      
c setting data for initializing t and at arrays

      emax = epsil(z)                    ! settinng emax to epsilon for point z
      nbin = nint(emax/de) + 1           ! number of bins
      b=aarr(z)                          ! 1d external rotational constant for K-rotor
      rmax=dsqrt((nbin*de)/b)            ! max "k" for krotor 
      imax=nint(rmax)+1                  ! bins for krotor

      if(z.eq.1) nbinl(jval+1)=nbin     ! for a given a value of J = jval, nbinl = max number of energy bins for reactant 

c      
c dynamic allocation of memory      ******** ALLOCATE ************
c
      allocate( t(nbin) , at(nbin) , ir(imax) , STAT=status )

      allocate_ok: IF ( status .EQ. 0 ) THEN   
!         WRITE (lunit,990)
        ELSE
         WRITE (lunit,991)
         STOP 'terminated in sterab'
      END IF allocate_ok
990   FORMAT ('sterab: allocated t(nbin),at(nbin),ir(imax)',/)
991   FORMAT ('sterab: allocation FAILED: t(nbin),at(nbin),ir(imax)',/)
c       
c initalize t and at arrays
      do i=1,nbin
         t(i)=0.0d0
        at(i)=0.0d0
      end do  
      t(1)=1.0d0

c
c loop over internal degrees of freedom

      do i = 1 , (ndof(z))
       if    (idofl(z,i).eq.'vib')then                   ! ********  vibrations                 ********
          wee =wel(z,i)
          anhh=anhl(z,i)
          nng =ngl(z,i)
          call morlev(at,t,nbin,de,nbin,wee,anhh,nng,zap)
       elseif(idofl(z,i).eq.'hra')then                   ! ********  quantized hindered rotors  ********
          wee =wel(z,i)
          anhh=anhl(z,i)
          nng =ngl(z,i)
          bb  = amua2nu/anhl(z,i)                        ! rotational constant (cm-1)
          vv  = ((wel(z,i)/ngl(z,i))**2)/(bb)            ! hindrance barrier (cm-1)
          ngg = ngl(z,i)                                 ! potential energy symmetry (foldedness)
          call shrlev(t,at,de,nbin,bb,vv,int(ngg),1,imax,zap,'vhrd1',
     $                                             'bhrd1',0.0d0,0.0d0)
        elseif(idofl(z,i).eq.'hrb')then
          vv=anhl(z,i)                                   ! hindrance barrier (cm-1)
          ngg = ngl(z,i)                                 ! potential energy symmetry (foldedness)
          bb  = ((wel(z,i)/ngl(z,i))**2)/(vv)            ! rotational constant (cm-1)
          call shrlev(t,at,de,nbin,bb,vv,int(ngg),1,imax,zap,'vhrd1',
     $                                             'bhrd1',0.0d0,0.0d0)
        elseif(idofl(z,i).eq.'hrc')then
          bb  = amua2nu/wel(z,i)                        ! rotational constant (cm-1)
          vv  = anhl(z,i)                               ! hindrance barrier (cm-1)
          ngg = ngl(z,i)                                ! potential energy symmetry (foldedness)
          call shrlev(t,at,de,nbin,bb,vv,int(ngg),1,imax,zap,'vhrd1',
     $                                             'bhrd1',0.0d0,0.0d0)
        elseif(idofl(z,i).eq.'hrd')then

          do k=1, ncbl(z,i)
             cbb(k)=cbl(z,i,k)
          end do

          do k=1, ncvl(z,i)
             cvv(k)=cvl(z,i,k)
          end do
          call uhrlev(t,at,de,nbin,ncbl(z,i),ncvl(z,i),cbb,cvv,		! Hindered internal rotation
     $           nsvl(z,i),nsbl(z,i),valmax(z,i),zap,vhrl(z,i),
     $           bhrl(z,i),phavl(z,i),phabl(z,i),1,ngl(z,i))
c     $           bhrl(z,i),phavl(z,i),phabl(z,i),ngl(z,i),1)

        elseif(idofl(z,i).eq.'jro')then                   ! ignore 2D external rotation
          continue
        elseif(idofl(z,i).eq.'kro')then                   ! ignore 1D external rotation
          continue
       endif
      end do  

      do i=1,nbin
         tmat(z,i)=t(i)
      end do
c
c release dynamically allocated memory     ******** DEALLOCATE ************
c
      deallocate( ir , at , t , STAT=status)
      deallocate_ok: IF ( status .EQ. 0 ) THEN   
!         WRITE (lunit,992)
        ELSE
         WRITE (lunit,993)
         STOP 'terminated in sterab'
      END IF deallocate_ok
992   FORMAT ('sterab: deallocated t(nbin),at(nbin),ir(imax)',/)
993   FORMAT ('sterab: deallocation FAILED',/)

      return

      end subroutine sterab

!-----------------------------------------------------------------------

      subroutine sterabj(z,jval)
      use DECLARE_MOD				! use module DECLARE_MOD
      include 'DECLARE_INC.inc'		! include 'DECLARE_INC.inc'

!	DEFINITIONS
!	z		: index number for reactants and trial ts, looped in subroutine micro_rates
!	nbin	: number of bins, corresponding to epsil(z)
!	tmat	: matrix of Beyer-Swinehart tmat(z,i) calculated by subroutine sterab with emax = epsil(z)
!	jval	: value of J 
!

      integer z,imax,jarg,nbin,jval,passno
      integer r(ndof(z)),vcnt,hracnt,nrg
      real*8  summ,rmax
      real*8  emax,b,v(ndof(z)),vv
      real*8,allocatable,dimension(:)::t,at
      integer,allocatable,dimension(:)::ir
      real*8  bb,vvv,ngg,zap
      real*8  hral(ndof(z)),zzpe(ndof(z))
c      real*8  gmat(rcnt+ntts,etop)
c      real*8  tmat(rcnt+ntts,etop)
      real*8  wee,anhh,twee
      integer nng,cnt,anh,we
      integer :: status						! Status: 0 for success
      logical tt,ff
 
 100  format(i3,2x,f10.2,2x,i5) 

c      
c settinng data for initializingg t and at arrays

      emax = epsil(z)                    ! settinng emax to epsilon for point z
      nbin = nint(emax/de) + 1           ! number of bins
      b=aarr(z)                          ! 1d external rotational constant
      rmax=dsqrt((nbin*de)/b)            ! max "k" for krotor 
      imax=nint(rmax)+1                  ! bins for krotor

      if(z.eq.1) nbinl(jval+1)=nbin     ! for a given a value of J = jval, nbinl = max number of energy bins for reactant

c      write(*,*) z,' nbin: :', nbin, '    nbinmax: ', nbinmax
c      
c dynamic allocation of memory      ******** ALLOCATE ************

      allocate( t(nbin) , at(nbin) , ir(imax) , STAT=status )

      allocate_ok: IF ( status .EQ. 0 ) THEN   
!         WRITE (lunit,990)
        ELSE
         WRITE (lunit,991)
         STOP 'terminated in sterab'
      END IF allocate_ok
990   FORMAT ('sterabj: allocated t(nbin),at(nbin),ir(imax)',/)
991   FORMAT ('sterabj: allocation FAILED: t(nbin),at(nbin),ir(imax)',/)

c       
c initalize t and at arrays

       do i=1,nbin
         t(i)=tmat(z,i)
        at(i)=0.0d0
       end do  
c
c loop over degrees of freedom

111   continue

      do i = 1 , (ndof(z))
       if(idofl(z,i).eq.'kro')then                   ! ********  1-dim. free k-rotor        ********
          wee =wel(z,i)
          anhh=anhl(z,i)
          nng =ngl(z,i)
       elseif(idofl(z,i).eq.'jro')then
          twee=wel(z,i)
       end if
       
      end do 
      wee=wee-twee
      call krotlev(at,t,nbin,de,nbin,wee,twee,jval)

c
c include 2j+1 external rotation contribution
      do i=1, nbin
         t(i)=((2.0d0*jval)+1.0d0)*t(i)
      end do

c  
c sum up densities of states in t and store in vector at
      summ=0.0d0
      do i=1,nbin
         summ=summ+t(i)
         at(i)=summ
         t(i)=t(i)/de
      end do
c
c    find sums above veff max (ecrit) and store in gmat, store number of grains above vmax in smat, and densities in minpmat
      cnt=0
      if(z.eq.1)then                            ! store dos of reactant
         do i=1,nbin
            rpmat(jval+1,i)=t(i)                              ! store full dos
            nrg=(i-1)*de
            if(nrg.ge.((vmax-vef(z))))then                    ! store info above critical energy
               cnt = cnt + 1
               if(cnt.gt.binmax)exit
               minpmat(jval+1,cnt)=t(i)
               gmat(z,cnt)=at(i)
            end if
         end do
      else
         do i=1,nbin                                          ! store sos of trial ts
            nrg=(i-1)*de
            if(nrg.ge.((vmax-vef(z))))then
               cnt = cnt + 1
               if(cnt.gt.binmax)exit
               gmat(z,cnt)=at(i)
            end if
         end do
      end if

c
c release dynamically allocated memory     ******** DEALLOCATE ************
      deallocate( ir , at , t , STAT=status)
      deallocate_ok: IF ( status .EQ. 0 ) THEN   
!         WRITE (lunit,992)
        ELSE
         WRITE (lunit,993)
         STOP 'terminated in sterab'
      END IF deallocate_ok
992   FORMAT ('sterabj: deallocated t(nbin),at(nbin),ir(imax)',/)
993   FORMAT ('sterabj: deallocation FAILED',/)

      return
      end subroutine sterabj

!-----------------------------------------------------------------------

c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c uhrlev: a code to calculate eigenvalues of unsymmetrical, non-rigid 1d-hindered internal rotation
c copyright (c) 2009 lam t. nguyen and john r. barker
c
c     date: feb. 13, 2009
c
c john r. barker
c jrbarker@umich.edu
c university of michigan
c ann arbor, mi 48109-2143
c (734) 763 6239
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        subroutine uhrlev(t,at,dele,jmax,ncb,ncv,cb,cv,n1,n2,imax,
     &          zpe,vhr,bhr,phav,phab,control,nsig)

      implicit double precision(a-h,o-z), integer(i-n)
      parameter (nn=2000)
      integer control
      dimension ev(nn), cv(ncv), cb(ncb), t(jmax), at(jmax)
       dimension ir(nn)
        character vhr*5, bhr*5

        pi=acos(-1.0d0)
        emax=(jmax-1)*dele + 1000.              ! *****raise max energy to account for later subtraction of zpe******

        if((bhr.eq.'bhrd2').or.(bhr.eq.'bhrd2')) then   ! average rotational constant
                b=0.0d0
                do i=1, n2
                        b=b + (cb(i)/i)*((pi*2)**(i-1))
                enddo
        elseif((bhr.eq.'ihrd2').or.(bhr.eq.'ihrd2')) then       ! average moment of inertia
                b=0.0d0
                do i=1, n2
                        b=b + (cb(i)/i)*((pi*2)**(i-1))
                enddo
                b=16.85763d0 / b        ! convert from i to b
        elseif((bhr.eq.'ihrd1').or.(bhr.eq.'ihrd1')) then       ! average moment of inertia
                b=16.85763d0/cb(1)      ! convert from i to b
        elseif((bhr.eq.'bhrd1').or.(bhr.eq.'bhrd1')) then       ! average rotational constant
                b=cb(1)
        else
                write(*,*) "error at input for bhr or ihr"
        endif

        call ghrlev(ev,nn,imax,emax,b,ncv,ncb,cv,cb,n1,n2,
     &          vhr,bhr,phav,phab)

      zpe=ev(1)

      if(control.eq.0)return

c*********************************************************************
c     the following lines are for output of hindered rotor eigenvalues
c*********************************************************************
c        write(12,*) '    i     e(i)      e(i)-zpe     [uhrlev]'

        do i=1, imax+1
c          write(12, 99) i, ev(i), (ev(i)-zpe)
        enddo
99    format (3x, i5, 1x, f10.3 , 1x, f10.3 )
c*********************************************************************

      do i=1, imax
            tp=ev(i)-zpe
            ir(i)= nint( tp/dele )
            ir(i)= nint( tp/dele )
c......           write(*,*) i, tp
      enddo 

      do j=1, imax
            ll=ir(j)
            do l=1, jmax - ll
                  at(l+ll) = at(l+ll) + t(l)
            enddo
      enddo
      do j=1, jmax
            t(j)=at(j)/nsig
            at(j)=0.0d0
      enddo 

      return
      
      end subroutine uhrlev

!-----------------------------------------------------------------------

      subroutine vx(v,x,cv,n,nsig,vhr,vm)
      implicit double precision(a-h,o-z), integer(i-n)
      dimension cv(n)
      character vhr*5

      if((vhr.eq.'vhrd1').or.(vhr.eq.'vhrd1')) then
            v=0.0d0
            do i=1, n
                  v=v+cv(i)*(1.0d0-dcos(x*nsig*i))/2.0d0
            enddo
      elseif((vhr.eq.'vhrd2').or.(vhr.eq.'vhrd2')) then
            v=cv(1)
            do i=2, n
                  v=v+cv(i)*dcos(x*nsig*(i-1))
            enddo
        elseif((vhr.eq.'vhrd3').or.(vhr.eq.'vhrd3')) then
            v=cv(1)
            nc=(n+1)/2
            do i=2, n
                  if(i.le.nc) then
                        v=v+cv(i)*dcos(x*nsig*(i-1))
                  else
                        v=v+cv(i)*dsin(x*nsig*(i-nc))
                  endif
            enddo
      else
            write(*,*) "error at input for torsional function = vhrd"
            stop
      endif

      v=v-vm

      return
      end subroutine vx

!-----------------------------------------------------------------------

      subroutine write_mat(mat,dimmi,dimmj,de,dj,vr,vt,uname,fname,sw,
     $                     vl,tstmp,nname)
      
      integer i,j,dimmi,dimmj,tunit,dj,sw,one
      real*8 de,vl
      real*8 mat(dimmi,dimmj),vr(dimmi),vt(dimmi)
      character(len=*) fname,uname,nname
      character frmt*60,namefile*60,tstmp*10
      parameter (tunit=24,one=1)
     
      if(uname.eq.'')then
         namefile='dat.'//trim(fname)
      else
         namefile=uname(1:len_trim(uname)-4)//'.'//trim(fname)
      end if

      namefile=trim(namefile)
      nname=namefile
      open(tunit,file=namefile)
c      call stamp(tunit,2)
c      call dnt(tunit)
c      write(tunit,*)
c      write(tunit,*)
      write(tunit,*)'run id: ',tstmp
      write(tunit,*)'dimensions: ',(dimmi-1)/dj,dimmj,de
      write(tunit,*)vl
      if(dimmi.eq.1)then
            write(frmt,'("(14x,i9,3x,i15,3x)")')
      else
            write(frmt,'("(14x,i9,3x,",i0,"(i15,3x))")')dimmi-1
      end if
      write(tunit,frmt)(i-1,i=1,dimmi,dj)
      write(frmt,'("(a14,",i0,"(es15.6,3x))")')dimmi

      if(sw.gt.1)then
         if(fname.eq.'kej'.or.fname.eq.'ukej')then
            write(tunit,frmt)'Erot-Re',(vr(i),i=1,dimmi,dj)
            write(tunit,frmt)'Erot-TS',(vt(i),i=1,dimmi,dj)
         elseif(fname.eq.'2dens')then
            write(tunit,frmt)'Erot-Re',(vr(i),i=1,dimmi,dj)
         elseif(fname.eq.'2sums')then
            write(tunit,frmt)'Erot-TS',(vt(i),i=1,dimmi,dj)
         elseif(fname.eq.'u2sums')then
            write(tunit,frmt)'Erot-TS',(vt(i),i=1,dimmi,dj)
         end if
      else
         if(fname.eq.'2dens')then
            write(tunit,frmt)'Erot',(vr(i),i=1,dimmi,dj)
         else if(fname.eq.'2sums')then
c            write(tunit,frmt)'Erot',(vt(i),i=1,dimmi,dj)
            write(tunit,frmt)'Erot',(vr(i),i=1,dimmi,dj)
         end if
      end if
      
      write(tunit,*)
      write(frmt,'("(f10.1,4x,",i0,"(es15.6,3x))")')dimmi
      do j=1,dimmj
         write(tunit,frmt)(j-1)*de,(mat(i,j),i=1,dimmi,dj)
      end do
      write(tunit,*)
      close(tunit)

      end subroutine write_mat

!-----------------------------------------------------------------------

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
      subroutine zpe_calc
      use DECLARE_MOD				! use module DECLARE_MOD
      include 'DECLARE_INC.inc'		! include 'DECLARE_INC.inc'
      integer nn,no,ncv,ncf,n1,n2
      integer ncb,imax
      real*8  bt(datpts),bat(datpts)
      real time
      real*8 zpel(datpts)
      real*8 b,v,ngg,zap
      parameter (nn=2000,no=1)
      real*8 cv(no),cf(no),vo,vm
      real*8 phav,phab,emax,delmin
      character vhr*5,bhr*5
      character*80 efilename
      save

 111  format(5x,22(es20.8))
 222  format(1x,a3,2x,4(a12,2x))
 333  format(1x,i3,2x,4(f12.4,2x))
 444  format(a12,2x,a12,2x,10(a12,2x))
 555  format(9x,i3,9x,f5.2,2x,10(f12.4,2x))

      call sestamp('zpe_calc',1)

c
c create array zpel containing zpe data for 
c all structures rcnt + ntts

c loop over all structures rcnt + ntts
      do i=1,rcnt+ntts+pcnt
         zpel(i)=0.0d0
         do j=1,ndof(i)
            if(idofl(i,j).eq.'vib')then
               zpel(i)=zpel(i) + 0.50d0*(wel(i,j)*ngl(i,j))    ! zpe for vibrational dof is (0.5)*vibrational frequency
            else if(idofl(i,j).eq.'hra')then                   ! zpe for hindered rotors is calculated in ghrlev
               call nrgconvert(anhl(i,j),1,keyword(i,2),'cm-1')
               b  = anhl(i,j)                                  ! rotational constant (cm-1)
               v  = ((wel(i,j)/ngl(i,j))**2)/(b)               ! hindrance barrier (cm-1)
               ngg = ngl(i,j)                                  ! potential energy symmetry (foldedness)
               vhr='vhrd1'
               bhr='bhrd1'
               ncv=1
               ncf=1
               cv(1)=v
               cf(1)=b
               emax=20.0d0*(0.6950356D+00)*maxtemp
               n1=ngg
               n2=1
               phav=0.0d0
               phab=0.0d0
               vo=0.0d0
               vm=0.0d0
c               write(*,*)
c               write(*,*)'shr in ',emax,b,v,n1,n2,valmax(i,j),zap,vhr,
c     $bhr,phav,phab,vo,vm
               call tshrlev(emax,ev,b,v,n1,n2,valmax(i,j),zap,vhr,bhr,
     $                      phav,phab,vo,vm)
c               write(*,*)'shr out',emax,b,v,n1,n2,valmax(i,j),zap,vhr,
c     $bhr,phav,phab,vo,vm
c               write(*,*)
               do k=1,valmax(i,j)
                  evhl(i,j,k)=ev(k)
               end do
               zap = ev(1)
               zpel(i)=zpel(i) + zap
            else if(idofl(i,j).eq.'hrb')then                 ! to be implimented
               v  = anhl(i,j)                                ! hindrance barrier (cm-1)
               b  = (1.0d0/v)*((wel(i,j)/ngl(i,j))**2)       ! rotational constant (cm-1)
               ngg = ngl(i,j)                                ! potential energy symmetry (foldedness)
               vhr='vhrd1'
               bhr='bhrd1'
               ncv=1
               ncf=1
               cv(1)=v
               cf(1)=b
               emax=20.0d0*(0.6950356D+00)*maxtemp
               n1=ngg
               n2=1
               phav=0.0d0
               phab=0.0d0
               vo=0.0d0
               vm=0.0d0
               call tshrlev(emax,ev,b,v,n1,n2,valmax(i,j),zap,vhr,bhr,
     $                      phav,phab,vo,vm)
               do k=1,valmax(i,j)
                  evhl(i,j,k)=ev(k)
               end do
               zap = ev(1)
               zpel(i)=zpel(i) + zap
            else if(idofl(i,j).eq.'hrc')then
               v  = anhl(i,j)                                ! hindrance barrier (cm-1)
               b  = amua2nu/wel(i,j)                         ! rotational constant (cm-1)
               ngg = ngl(i,j)                                ! potential energy symmetry (foldedness)
               vhr='vhrd1'
               bhr='bhrd1'
               ncv=1
               ncf=1
               cv(1)=v
               cf(1)=b
               emax=20.0d0*(0.6950356D+00)*maxtemp
               n1=ngg
               n2=1
               phav=0.0d0
               phab=0.0d0
               vo=0.0d0
               vm=0.0d0
               call tshrlev(emax,ev,b,v,n1,n2,valmax(i,j),zap,vhr,bhr,
     $                      phav,phab,vo,vm)
               do k=1,valmax(i,j)
                  evhl(i,j,k)=ev(k)
               end do
               zap = ev(1)
               zpel(i)=zpel(i) + zap
            else if(idofl(i,j).eq.'hrd')then

               ncb=int(anhl(i,j))
               do k=1, ncb
                  cbb(k)=cbl(i,j,k)
               enddo

               ncv=int(wel(i,j))
               do k=1, ncv
                  cvv(k)=cvl(i,j,k)
               enddo

               call uhrlev(bt,bat,de,etopuser,ncb,ncv,cbb,cvv,
     $              nsvl(i,j),nsbl(i,j),valmax(i,j),zap,vhrl(i,j),
     $              bhrl(i,j),phavl(i,j),phabl(i,j),0, ngl(i,j))
c     $              bhrl(i,j),phavl(i,j),phabl(i,j),ngl(i,j),0,ev)

               do k=1, valmax(i,j)
                  evhl(i,j,k)=ev(k)
               end do
               zpel(i)=zpel(i) + zap

            end if
         end do
      end do

      write(lunit,*)'# of structures:',rcnt+ntts+pcnt
      write(lunit,*)
      write(lunit,*)'reporting v, zpe, and electronic (cm-1)'
      write(lunit,*)
      write(lunit,222)'#','v','zpe','electronic'
      write(lunit,*)('-',i=1,45)
      do i=1, rcnt+ntts+pcnt
         write(lunit,333)i,delh(i),zpel(i),(delh(i)-zpel(i))-
     $(delh(1)-zpel(1))
      end do
      write(lunit,*)

      if(rcnt+ntts+pcnt.ge.2)then
         if((inputfile.eq.'').or.(inputfile.eq.'ktools.dat'))then
            efilename='efile.txt'
         else
            efilename=trim(fileroot)//"-efile.txt"
         endif

         efilename=trim(efilename)

         open(unit=efunit,file=efilename)
         call stamp(efunit,2)
         call dnt(efunit)
         write(efunit,*)'run id: ',tstmp
         write(efunit,*)
         write(efunit,*)
         write(efunit,444)'struct#','rxncoord(a)','dh(cm-1)','zpe(cm-1)'
     $,'delect(cm-1)','k-rot(cm-1)','j-rot(cm-1)'
         do i=1,rcnt+ntts+pcnt
            write(efunit,555)i,distl(i),delh(i),zpel(i),
     $(delh(i)-zpel(i))-(delh(1)-zpel(1)),aarr(i),barr(i)
         end do
         close(efunit)
      end if

      call sestamp('zpe_calc',2)

      return

      end subroutine zpe_calc
      
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

      END MODULE MICRO_MOD
