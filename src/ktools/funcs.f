      subroutine funcs(x,p,np)
      real*8 x,p
      integer np,j
      dimension p(np)

      p(1)=1.0d0
      do 11 j=2,np
            p(j)=p(j-1)*x
 11   continue
      return
      end
