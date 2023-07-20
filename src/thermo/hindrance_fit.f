       subroutine hindrance_fit

       include 'declare.inc'
       REAL(8) a0
       SAVE

       x=1
       DO i=1,200             ! Max iterations

         call thermochem                         ! THERMOCHEMISTY CALCULATION
       
         Kcinold=Kcin
         Kcin=AK(2)                              ! AK(2) = rate constant (from thermochem) at middle T=Temp
         if(i.eq.1) then
           x=1.010
         elseif(i.gt.1) then 
           if(i.eq.2) xold=1.000
!          Approx:  Y=aX**m
           m=log(Kcin/Kcinold)/log(x/xold)
           a0=Kcin/(x**m)
           xold2=x

!          Y'=a*m*X**(m-1)
           x=x-((Kcin-Kexp(Temp))/(m*a0*x**(m-1)))   !Newton method with approx Y'
!          x=x-((Kcin-Kexp(temp))/((Kcin-Kcinold)/(x-xold)))   !Newton method
c          write(*,*) i, Kexp(temp), Kcin, x

           if(x.lt.0) x=ABS(xold/2)
           xold=xold2
           if( (ABS(x-xold).lt.1E-5) .AND. 
     &           (ABS(Kcin-Kexp(Temp))/Kexp(Temp) .lt. 1E-5) ) then    ! Convergence thresholds
              gamma(Temp)=x                   ! gamma = hindrance = sqrt(1-eta) 
c             write(*,*) i,' iterations'
              return
           endif
         endif
      ENDDO
      
      write(*,*) "FATAL: Unable to fit experimental rate constant ****"
      write(*,*) '      ',i,' iterations'
      stop
      
      end subroutine
 
