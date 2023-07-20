      subroutine rank(n,indx,irank)
      integer n,indx(n),irank(n)
      integer j

      do j=1,n
         irank(indx(j))=j
      end do

      return

      end
