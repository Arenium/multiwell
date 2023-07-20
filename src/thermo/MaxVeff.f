
      subroutine MaxVeff(beta0,POT,c,De,re,mu,T,rmax,I2dGorin,V,
     &   Veff,Eunits)

      Character(len=7) POT
      CHARACTER(len=4) Eunits
      REAL(8) De, beta0 , re, f, fold , c, T, step , z
      REAL(8) I2dGorin, rmax, Veff, V, bet, mu
      INTEGER(4) J
      REAL(8) a, aold , x
      SAVE

! Secant method 
      J=0
      IF (POT.eq."VARSHNI") THEN                                   ! Varshni  V=De((1-re/r) exp(-beta(r**2-re**2)))**2   beta= 0.5(beta0-1/re)/re
        bet = beta0
        a = 1.5*re
        f = 2.d+0*De*re*(exp(-bet*(a**2-re**2)))*
     &    (1.d+0-(re/a)*exp(-bet*(a**2-re**2)))*(2*bet+1.d+0/(a**2))-
     &    2.d+0*(T/1.4388d+0)/a
      ELSEIF (POT.eq."sMORSE" .OR. POT.EQ.'MORSE') THEN           ! Stiff Morse  V=De(1-exp(-beta(r-re)))**2   beta=beta0(1+C(r-re)**2)
        IF(POT .EQ. 'MORSE') c = 0.d+00
        a = 1.5*re*(1. - 0.1*c)
        x = (a-re)
        bet = beta0*(1. + c*x*x)
        z = bet*x
        f = 2.*De*beta0*(1.+3.*c*x*x)*(exp(-z) - exp(-2.*z))
     &        -2.*(T/1.4388)/a
      ENDIF
      
c       write(*,*) J, a , step, f

      step = 0.01*a
      
      DO WHILE ( ABS(step) .GT. 1.e-7 )
        fold = f
        aold = a
        a = aold + step
      
        if(a .gt. 20.) then
         write(*,*) "**** FATAL: Unable max Veff(r): r -> INF ****"
         stop
        endif

        J=J+1
        if(J .gt. 1000) then                     ! Max iterations  1.0E3
         write(*,*) "**** FATAL: Unable max Veff(r): max iters ****"
         stop
        endif   

        IF (POT.eq."VARSHNI") THEN                                  ! Varshni  V=De((1-re/r) exp(-beta(r**2-re**2)))**2   beta= 0.5(beta0-1/re)/re
        f = 2.d+0*De*re*(exp(-bet*(a**2-re**2)))*
     &    (1.d+0-(re/a)*exp(-bet*(a**2-re**2)))*(2*bet+1.d+0/(a**2))-
     &    2.d+0*(T/1.4388d+0)/a
          V = De*(1-(re/a)*exp(-beta*(a**2-re**2)))**2 - De
        ELSEIF (POT.eq."sMORSE" .OR. POT.EQ.'MORSE') THEN           ! Stiff Morse  V=De(1-exp(-beta(r-re)))**2   beta=beta0(1+C(r-re)**2)
          IF(POT .EQ. 'MORSE') c = 0.d+00
          x = (a-re)
          bet = beta0*(1. + c*x*x)
          z = bet*x
          f = 2.*De*beta0*(1.+3.*c*x*x)*(exp(-z) - exp(-2.*z))
     &        -2.*(T/1.4388)/a
          V = De*( 1.-exp(-z) )**2 - De
      ENDIF

       step =  -f*( a - aold )/( f - fold )                       ! Secant method
       IF (ABS(step) .GT. 0.5) step = 0.5*step/ABS(step)          ! limit stepsize
       
c       write(*,*) J, a , step, f
      ENDDO
      
         rmax=a         
         I2dGorin=mu*(a**2)
         Veff=V+(T/1.4388)
         IF ( Eunits.EQ.'KCAL' .OR. Eunits.EQ.'kcal' ) THEN
            V=V/349.76
            Veff=Veff/349.76
         ELSEIF ( Eunits.EQ.'KJOU' .OR. Eunits.EQ.'kjou' ) THEN
            V=V/83.59
            Veff=Veff/83.59
         ENDIF

      RETURN
      end subroutine

