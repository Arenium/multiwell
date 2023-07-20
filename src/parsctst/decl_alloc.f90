MODULE decl_alloc
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!    
!    LICENSE NOTICE
!
!    Copyright (C) 2019 Michele Ceotto, Chiara Aieta, Fabio Gabas, 
!                       Thanh Lam Nguyen, and John R. Barker
!
!    Contact:
!
!    Michele Ceotto  (email: michele.ceotto@unimi.it)
!    Dipartimento di Chimica
!    Università degli Studi di Milano
!    via Golgi 19, Milano 20133, ITALY
!
!    or:
!
!    John R. Barker   (email: jrbarker@umich.edu)
!    Department of Atmospheric, Oceanic, and Space Sciences
!    College of Engineering
!    University of Michigan
!    Ann Arbor, MI 48109
!
!    http://clasp-research.engin.umich.edu/multiwell
!    This program is free software; you can redistribute it and/or
!    modify it under the terms of the GNU General Public License (version 2)
!    as published by the Free Software Foundation.
!   
!    This program is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!    GNU General Public License for more details.
!   
!    See the 'ReadMe' file for a copy of the GNU General Public License,
!    or contact:
!   
!    Free Software Foundation, Inc.
!    59 Temple Place - Suite 330
!    Boston, MA 02111-1307, USA.
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!                                                                 
!                         PROGRAM parsctst    
!
!                               by      
!
!           by Michele Ceotto, Chiara Aieta, Fabio Gabas,
!               Thanh Lam Nguyen, and John R. Barker    
!                                                                  
!           ***PARalle Semi-Classical Transition State Theory***
!
!                             based on
!
!          The theory of W. H. Miller and coworkers* and the
!          parallel implementation** of the Wang-Landau algorithms 
!          for densities of states***
!
!    Literature Citations:                                                    
!    *Semi-Classical Transition State Theory
!    W. H. Miller, J. Chem. Phys. 62, 1899-1906 (1975).
!    W. H. Miller, Faraday Discuss. Chem. Soc. 62, 40-46 (1977).
!    W. H. Miller, R. Hernandez, N. C. Handy, D. Jayatilaka, and A. Willets,
!      Chem. Phys. Letters 172, 62-68 (1990).
!    R. Hernandez and W. H. Miller, Chem. Phys. Lett. 214, 129-136 (1993).
!    J. F. Stanton, J. Phys. Chem. Lett. 7, 2708-2713 (2016)
!
!    **Semi-Classical Transition State Theory parallel implementation
!    C. Aieta, F. Gabas and M. Ceotto, J. Chem. Theory Comput.,
!      15, 2142−2153 (2019).
!                                                                  
!    ***Density of states algorithms
!    F. Wang and D. P. Landau, Phys. Rev. Letters 86, 2050-2053 (2001).           
!    M. Basire, P. Parneix, and F. Calvo, J. Chem.Phys. 129, 081101 (2008).    
!    T. L. Nguyen and J.R. Barker, J. Phys. Chem. A., 114, 3718–3730 (2010).
!    C. Aieta, F. Gabas and M. Ceotto, J. Phys. Chem. A., 120(27), 4853-4862 (2016).
!                                                                  
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

    IMPLICIT NONE

    INTEGER(4) :: my_rank, num_procs, err
    INTEGER(4) nwind, nwalkers,idum
    INTEGER(4) AllocateStatus, DeAllocateStatus, ngrains
    INTEGER(4), DIMENSION(:), ALLOCATABLE :: ngrains_to_add,ngrains_per_wind
    INTEGER(4), DIMENSION(:), ALLOCATABLE :: ngrains_per_chunk
    REAL(8), DIMENSION(:), ALLOCATABLE :: lowbound,upbound,wind_time,wind_time_red
    REAL(8), DIMENSION(:), ALLOCATABLE :: Ttot,AT,DS,SS, RI, RJ
    REAL(8), DIMENSION(:), ALLOCATABLE :: wa, w0, wf, D, g_ave
    REAL(8), DIMENSION(:,:), ALLOCATABLE :: xa
    REAL(8), DIMENSION(:,:,:), ALLOCATABLE :: ya
    REAL(8), DIMENSION(:,:,:,:), ALLOCATABLE :: za
    REAL(8), DIMENSION(:), ALLOCATABLE :: chunk, half_perc_wind_overlap
    INTEGER(4), DIMENSION(:,:,:), ALLOCATABLE :: H, H_red
    REAL(8), DIMENSION(:,:,:), ALLOCATABLE :: g, g_red
    REAL(8), DIMENSION(:), ALLOCATABLE::HRfreq, HRVo, HRI, HRB,HRzpe,Phav,Phab
    INTEGER(4), DIMENSION(:), ALLOCATABLE:: NG,IMAX,NVV,NBB,NGV,NGB,NSIG
    INTEGER(4), DIMENSION(:), ALLOCATABLE:: IDIM,MODE,NRSN
    INTEGER(4), DIMENSION(:), ALLOCATABLE :: nv,nvold
    INTEGER(4), DIMENSION(:), ALLOCATABLE :: idum_cache, iy_cache, nold_cache
    INTEGER(4), DIMENSION(:,:), ALLOCATABLE :: iv_cache, nvold_cache
    REAL(8), DIMENSION(:,:), ALLOCATABLE :: Tw, Tw_red
    INTEGER(4), DIMENSION(:), ALLOCATABLE :: ref, ref_red
    REAL(8) :: zpe, ave, av0, avf, Viblo
    REAL(8), DIMENSION(:), ALLOCATABLE :: XI, TTFtot , Evtot
    INTEGER(8), DIMENSION(:), ALLOCATABLE :: H_check
    REAL(8), DIMENSION(:,:), ALLOCATABLE :: TTF_ave , Ev_ave, TTF_red, Ev_red
    REAL(8), DIMENSION(:,:,:), ALLOCATABLE :: TTF, Ev
    INTEGER(4) ngrains_corr
    REAL(8) :: ngrains_per_chunk_small, ngrains_per_chunk_big
    REAL(8) perc_wind_overlap 
    REAL(8) Emin, Emax, Emin_global, Emax_global 
    REAL(8) Egrain1 , Egrain2 , Emax1 , Emax2
    INTEGER(4) :: ngrains_per_wind_max, ngrains_to_add_max

CONTAINS

    SUBROUTINE alloc(nwind)
        INTEGER(4),INTENT(IN) :: nwind
               
        ALLOCATE ( ngrains_per_chunk(nwind), STAT = AllocateStatus)
        IF (AllocateStatus /= 0) STOP "*** Not enough memory ngrains_per_chunk***"

        ALLOCATE ( half_perc_wind_overlap(nwind), STAT = AllocateStatus)
        IF (AllocateStatus /= 0) STOP "*** Not enough memory half_perc_wind_overlap***"

        ALLOCATE ( chunk(nwind), STAT = AllocateStatus)
        IF (AllocateStatus /= 0) STOP "*** Not enough memory chunk***"

        ALLOCATE ( ngrains_to_add(nwind), STAT = AllocateStatus)
        IF (AllocateStatus /= 0) STOP "*** Not enough memory ngrains_to_add***"

        ALLOCATE ( ngrains_per_wind(nwind), STAT = AllocateStatus)
        IF (AllocateStatus /= 0) STOP "*** Not enough memory ngrains_per_wind***"

        ALLOCATE ( lowbound(nwind), STAT = AllocateStatus)
        IF (AllocateStatus /= 0) STOP "*** Not enough memory lowbound***"
        
        ALLOCATE ( upbound(nwind), STAT = AllocateStatus)
        IF (AllocateStatus /= 0) STOP "*** Not enough memory upbound ***"

        ALLOCATE ( wind_time(nwind), STAT = AllocateStatus)
        IF (AllocateStatus /= 0) STOP "*** Not enough memory wind_time ***"
        wind_time = 0.d0

        ALLOCATE ( wind_time_red(nwind), STAT = AllocateStatus)
        IF (AllocateStatus /= 0) STOP "*** Not enough memory wind_time ***"
        wind_time_red = 0.d0

    END SUBROUTINE alloc


    SUBROUTINE alloc_2(ngrains, ngrains_to_add_max,ngrains_per_wind_max, nwalkers,nwind,ns)
    
       INTEGER(4),INTENT(IN) :: ngrains, ngrains_to_add_max
       INTEGER(4),INTENT(IN) :: ngrains_per_wind_max, nwalkers, nwind, ns

       ALLOCATE ( nvold_cache(ns,nwalkers), STAT = AllocateStatus)
       IF (AllocateStatus /= 0) STOP &
       "*** Not enough memory nvold_cache***"
       nvold_cache = 0

       ALLOCATE ( g_ave(ngrains_per_wind_max), STAT = AllocateStatus)
       IF (AllocateStatus /= 0) STOP &
       "*** Not enough memory g_ave***"
       g_ave = 0.d0

       ALLOCATE ( ref(nwind), STAT = AllocateStatus)
       IF (AllocateStatus /= 0) STOP  &
       "*** Not enough memory ref***"
       ref = 0

       ALLOCATE ( Tw(ngrains_per_wind_max, nwind), STAT = AllocateStatus)
       IF (AllocateStatus /= 0) STOP  &
       "*** Not enough memory Tw***"
       Tw = 0

       ALLOCATE ( ref_red(nwind), STAT = AllocateStatus)
       IF (AllocateStatus /= 0) STOP  &
       "*** Not enough memory ref_red***"
       ref_red = 0

       ALLOCATE ( Tw_red(ngrains_per_wind_max, nwind), STAT = AllocateStatus)
       IF (AllocateStatus /= 0) STOP  &
       "*** Not enough memory Tw_red***"
       Tw_red = 0

       ALLOCATE ( Ttot(ngrains+ngrains_to_add_max), STAT = AllocateStatus)
       IF (AllocateStatus /= 0) STOP "*** Not enough memory Ttot***"
       Ttot(:) = 0.d0

       ALLOCATE ( AT(ngrains+ngrains_to_add_max), STAT = AllocateStatus)
       IF (AllocateStatus /= 0) STOP "*** Not enough memory AT***"
       AT(:) = 0.d0

       ALLOCATE ( DS(ngrains+ngrains_to_add_max), STAT = AllocateStatus)
       IF (AllocateStatus /= 0) STOP "*** Not enough memory DS***"
       DS(:) = 0.d0

       ALLOCATE ( SS(ngrains+ngrains_to_add_max), STAT = AllocateStatus)
       IF (AllocateStatus /= 0) STOP "*** Not enough memory SS***"
       SS(:) = 0.d0

       ALLOCATE ( H(ngrains_per_wind_max, nwalkers, nwind), STAT = AllocateStatus)
       IF (AllocateStatus /= 0) STOP "*** Not enough memory H ***"
       H(:,:,:) = 0

       ALLOCATE ( g(ngrains_per_wind_max, nwalkers, nwind), STAT = AllocateStatus)
       IF (AllocateStatus /= 0) STOP "*** Not enough memory g ***"
       g(:,:,:) = 0.d0

       ALLOCATE ( H_red(ngrains_per_wind_max, nwalkers, nwind), STAT =AllocateStatus)
       IF (AllocateStatus /= 0) STOP "*** Not enough memory H_red ***"
       H_red(:,:,:) = 0

       ALLOCATE ( g_red(ngrains_per_wind_max, nwalkers, nwind), STAT =AllocateStatus)
       IF (AllocateStatus /= 0) STOP "*** Not enough memory g_red ***"
       g_red(:,:,:) = 0.d0

       ALLOCATE ( idum_cache(nwalkers), STAT = AllocateStatus)
       IF (AllocateStatus /= 0) STOP  "*** Not enough memory idum_cache***"
       idum_cache = 0

       ALLOCATE ( iy_cache(nwalkers), STAT = AllocateStatus)
       IF (AllocateStatus /= 0) STOP  "*** Not enough memory iy_cache***"
       iy_cache = 0

       ALLOCATE ( iv_cache(32,nwalkers), STAT = AllocateStatus)
       IF (AllocateStatus /= 0) STOP  "*** Not enough memory iv_cache***"
       iv_cache = 0

       ALLOCATE ( nold_cache(nwalkers), STAT = AllocateStatus)
       IF (AllocateStatus /= 0) STOP "*** Not enough memory nold_cache***"
       nold_cache = 0

       ALLOCATE ( TTF(ngrains_per_wind_max,nwalkers,nwind), STAT = AllocateStatus)
       IF (AllocateStatus /= 0) STOP "*** Not enough memory TTF***"
       TTF = 0.d0

       ALLOCATE ( TTFtot(ngrains + ngrains_to_add_max), STAT = AllocateStatus)
       IF (AllocateStatus /= 0) STOP "*** Not enough memory TTFtot***"
       TTFtot = 0.d0

       ALLOCATE ( TTF_red(ngrains_per_wind_max,nwind), STAT = AllocateStatus)
       IF (AllocateStatus /= 0) STOP "*** Not enough memory TTF_red***"
       TTF_red = 0.d0

       ALLOCATE ( TTF_ave(ngrains_per_wind_max,nwind), STAT = AllocateStatus)
       IF (AllocateStatus /= 0) STOP "*** Not enough memory TTF_ave***"
       TTF_ave = 0.d0

       ALLOCATE ( Ev_red(ngrains_per_wind_max,nwind), STAT = AllocateStatus)
       IF (AllocateStatus /= 0) STOP "*** Not enough memory EV_red***"
       Ev_red = 0.d0

       ALLOCATE ( Ev(ngrains_per_wind_max,nwalkers,nwind), STAT = AllocateStatus)
       IF (AllocateStatus /= 0) STOP "*** Not enough memory Ev***"
       Ev = 0.d0

       ALLOCATE ( Ev_ave(ngrains_per_wind_max,nwind), STAT = AllocateStatus)
       IF (AllocateStatus /= 0) STOP "*** Not enough memory EV_ave***"
       Ev_ave = 0.d0

       ALLOCATE ( EVtot(ngrains + ngrains_to_add_max), STAT = AllocateStatus)
       IF (AllocateStatus /= 0) STOP "*** Not enough memory EVtot***"
       Evtot = 0.d0
  
       ALLOCATE ( H_check(ngrains_per_wind_max), STAT = AllocateStatus)
       IF (AllocateStatus /= 0) STOP "*** Not enough memory H_check***"
       H_check = 0.d0


    END SUBROUTINE alloc_2

    
    SUBROUTINE alloc_3(ns)
        INTEGER(4),INTENT(IN) :: ns

        ALLOCATE ( nv(ns), STAT = AllocateStatus)
        IF (AllocateStatus /= 0) STOP "*** Not enough memory nv***"
        nv(:) = 0

        ALLOCATE ( nvold(ns), STAT = AllocateStatus)
        IF (AllocateStatus /= 0) STOP "*** Not enough memory nvold***"
        nvold(:) = 0

        ALLOCATE ( wa(ns), STAT = AllocateStatus)
        IF (AllocateStatus /= 0) STOP "*** Not enough memory wa***"

        ALLOCATE ( w0(ns), STAT = AllocateStatus)
        IF (AllocateStatus /= 0) STOP "*** Not enough memory w0***"

        ALLOCATE ( wf(ns), STAT = AllocateStatus)
        IF (AllocateStatus /= 0) STOP "*** Not enough memory wf***"
        
        ALLOCATE ( xa(ns,ns), STAT = AllocateStatus)
        IF (AllocateStatus /= 0) STOP "*** Not enough memory xa***"

        ALLOCATE ( ya(ns,ns,ns), STAT = AllocateStatus)
        IF (AllocateStatus /= 0) STOP "*** Not enough memory ya***"

        ALLOCATE ( za(ns,ns,ns,ns), STAT = AllocateStatus)
        IF (AllocateStatus /= 0) STOP "*** Not enough memory za***"

        ALLOCATE ( D(ns), STAT = AllocateStatus)
        IF (AllocateStatus /= 0) STOP "*** Not enough memory D***"

        ALLOCATE ( XI(ns), STAT = AllocateStatus)
        IF (AllocateStatus /= 0) STOP "*** Not enough memory XI***"

     END SUBROUTINE alloc_3


     SUBROUTINE alloc_4(NSEP,ns)
        INTEGER(4),INTENT(IN) :: NSEP, ns

        ALLOCATE ( HRfreq(NSEP+ns), STAT = AllocateStatus)
        IF (AllocateStatus /= 0) STOP "*** Not enough memory HRfreq***"

        ALLOCATE ( HRVo(NSEP+ns), STAT = AllocateStatus)
        IF (AllocateStatus /= 0) STOP "*** Not enough memory HRVo***"

        ALLOCATE ( HRI(NSEP+ns), STAT = AllocateStatus)
        IF (AllocateStatus /= 0) STOP "*** Not enough memory HRI***"

        ALLOCATE ( HRB(NSEP+ns), STAT = AllocateStatus)
        IF (AllocateStatus /= 0) STOP "*** Not enough memory HRB***"

        ALLOCATE ( HRzpe(NSEP+ns), STAT = AllocateStatus)
        IF (AllocateStatus /= 0) STOP "*** Not enough memory HRzpe***"

        ALLOCATE ( Phav(NSEP+ns), STAT = AllocateStatus)
        IF (AllocateStatus /= 0) STOP "*** Not enough memory Phav***"
        Phav = 0

        ALLOCATE ( Phab(NSEP+ns), STAT = AllocateStatus)
        IF (AllocateStatus /= 0) STOP "*** Not enough memory Phab***"
        Phab = 0

        ALLOCATE ( RI(NSEP+ns), STAT = AllocateStatus)
        IF (AllocateStatus /= 0) STOP "*** Not enough memory RI***"

        ALLOCATE ( RJ(NSEP+ns), STAT = AllocateStatus)
        IF (AllocateStatus /= 0) STOP "*** Not enough memory RJ***"

        ALLOCATE ( NG(NSEP+ns), STAT = AllocateStatus)
        IF (AllocateStatus /= 0) STOP "*** Not enough memory NG***"

        ALLOCATE ( IMAX(NSEP+ns), STAT = AllocateStatus)
        IF (AllocateStatus /= 0) STOP "*** Not enough memory IMAX***"

        ALLOCATE ( NVV(NSEP+ns), STAT = AllocateStatus)
        IF (AllocateStatus /= 0) STOP "*** Not enough memory NVV***"

        ALLOCATE ( NBB(NSEP+ns), STAT = AllocateStatus)
        IF (AllocateStatus /= 0) STOP "*** Not enough memory NBB***"

        ALLOCATE ( NGV(NSEP+ns), STAT = AllocateStatus)
        IF (AllocateStatus /= 0) STOP "*** Not enough memory NGV***"

        ALLOCATE ( NGB(NSEP+ns), STAT = AllocateStatus)
        IF (AllocateStatus /= 0) STOP "*** Not enough memory NGB***"

        ALLOCATE ( NSIG(NSEP+ns), STAT = AllocateStatus)
        IF (AllocateStatus /= 0) STOP "*** Not enough memory NSIG***"

        ALLOCATE ( MODE(NSEP+ns), STAT = AllocateStatus)
        IF (AllocateStatus /= 0) STOP "*** Not enough memory MODE***"

        ALLOCATE ( NRSN(NSEP+ns), STAT = AllocateStatus)
        IF (AllocateStatus /= 0) STOP "*** Not enough memory NRSN***"

        ALLOCATE ( IDIM(NSEP+ns), STAT = AllocateStatus)
        IF (AllocateStatus /= 0) STOP "*** Not enough memory IDIM***"

    END SUBROUTINE alloc_4



END MODULE decl_alloc
