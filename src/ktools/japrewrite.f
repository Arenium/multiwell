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
      subroutine japrewrite(iout,fname,egrain1,imax1,emax2,isize,viblo)
      USE DECLARE_MOD, ONLY: cut, aversion, adate
      real*8 emax2,egrain1,viblo
      integer imax1,isize,iout
      character(len=*) fname
99001 FORMAT (1x,f10.1,1x,I10,1x,f10.1,1x,I10,1x,f10.1)
99030 FORMAT (A46)
99016 FORMAT (
     &'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     &%%%%%%%%%%%%%'//8x,
     &'                              Jason A. Sonk & John R. Barker'/8x,
     &'        ktools-',A8,'       University of Michigan'/8x,
     &'                              Ann Arbor, MI 48109-2143'/8x,
     &'          ',A4,A4,'            Contact:',/8x,
     &'                                jrbarker@umich.edu'/8x,
     &'                                (734) 763 6239'//8x,
     &'      http://aoss.engin.umich.edu/multiwell/'//,
     &'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     &%%%%%%%%%%%%%'//8x,
     &'Suggested Literature Citation:'//8x,
     &'a) J.R. Barker, T.L. Nguyen, J.F. Stanton, C. Aieta, M. Ceotto,',
     &/7x,'F. Gabas, T.J.D. Kumar, C.G.L.Li, L.L.Lohr, A. Maranzana,',
     &/7x,'N.F. Ortiz, J.M. Preses, J.M. Simmie, J.A. Sonk, and ',
     &'P.J. Stimac; ',
     &/7x,'MultiWell-',A7,' Software Suite; University of Michigan, ',
     &/7x,'Ann Arbor, Michigan, USA, ',A8,'; ',
     &/7x,'http://clasp-research.engin.umich.edu/multiwell/.',//8x,
     &'b) J.R. Barker, Int. J. Chem. Kinetics, 33, 232-45 (2001)',//8x,
     &//,
     &'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     &%%%%%%%%%%%%%'/)
99010 FORMAT ('       No.    E(cm-1)    Density    Sum   ',
     &        '[Energy at TOP of energy grains]')


      write(iout,99030)cut
      write(iout,99016)aversion,adate,aversion,aversion,adate,aversion
      write(iout,99030)cut


      write(iout,*)trim(fname)
      write(iout,*)'J-summed Density and Sum of states.'
      write(iout,99001)egrain1,imax1,emax2,isize,viblo
      write(iout,99010)




      end subroutine
