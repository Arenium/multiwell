      subroutine parabfit(x,y,n,a,chisq)
      integer n,i
      integer ia(n)
      real*8 x(n),y(n),a(3),covar(3,3),sig(n),chisq
      external funcs

      do i=1,n
         ia(i)=1
          a(i)=1.0d0
      end do

      do i=1,n
         sig(i)=1.0d0
      end do

      call lfit(x,y,sig,n,a,ia,3,covar,3,chisq,funcs)

      return
      end subroutine
