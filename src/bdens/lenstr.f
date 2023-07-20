c
      function lenstr(string)
c
      implicit      none
c
      character*(*) string,  !input string
     2              null*1   !null character
c
      INTEGER(4)     i,       !string index
     2              iend,    !candidate for string length
     3              j,       !index used for stepping through string
     4              length,  !defined length of input
     5              lenstr   !length of string
c
      length=len(string)
c
c     Look for the existence of nulls; if there
c     are any, the first one marks one character past the
c     end of the defined string. If there are no nulls,
c     pick up the final nonspace character, and call that
c     the defined length of the string.
c
      null=char(0)
c
c     Look for nulls
c
      iend=index(string(1:length),null)
c
c     No nulls, look for the location of the last nonspace.
c
      if(iend.eq.0)then
          do i=1,length
              j=length+1-i
              if(string(j:j).ne.' ')then
                  lenstr=j
                  return
              endif
          enddo
c
c         all spaces
c
          lenstr=0
          return
      else
c
c         first null-1
c
          lenstr=iend-1
          return
      endif
      end
