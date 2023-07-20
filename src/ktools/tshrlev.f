c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c shrlev: a code to calculate eigenvalues of symmetrical, rigid 1d-hindered internal rotation
c copyright (c) 2009 lam t. nguyen and john r. barker
c
c      date: feb. 13, 2009
c
c john r. barker
c jrbarker@umich.edu
c university of michigan
c ann arbor, mi 48109-2143
c (734) 763 6239
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine tshrlev(emax,ev,b,v,n1,n2,nimax,zpe,
     &      vhr,bhr,phav,phab,vo,vm)

      implicit none
      integer(4) n1, n2, nimax, nn, no, ncv, ncf
      parameter (nn=2000, no=1)
      real(8) emax, ev, b,v, zpe, phav, phab, vo, vm
      real(8) cv, cf
      dimension ev(nn), cv(no), cf(no) 
      character(5) vhr , bhr
      save

      ncv=1
      ncf=1
      cv(1)=v
      cf(1)=b
      call tghrlev(ev,nn,nimax,emax,b,ncv,ncf,cv,cf,n1,n2,
     &      vhr,bhr,phav,phab,vo,vm)

      zpe=ev(1)	! zero-point energy (cm-1)


      return
      end


