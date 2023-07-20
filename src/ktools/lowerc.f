      subroutine lowerc(word)
      character (len=*) , intent(in out) :: word
      integer                            :: i,ic,nlen

      nlen = len(word)

      do i=1,nlen
         ic = ichar(word(i:i))
         if (ic.ge. 65 .and. ic.le.90) word(i:i) = char(ic+32)
      end do 

      end subroutine lowerc

