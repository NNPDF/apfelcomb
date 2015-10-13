C
C
C For passing strings without making use of the "byte" type, see the
C following page
C http://dynopt.cheme.cmu.edu/roscoe/RTOp/doc/Misc/html/namespaceFortranTypes.html
C
C
C
C     interface to initPDFset that can be called from C. Inspired by 
C     http://www-mipl.jpl.nasa.gov/portguide/subsubsection3.11.4.2.html
C
C     Assumes a one-byte-per-character null-terminated C string as input
C     Makes no assumption about the structure of fortran strings.
C     Maybe this is overkill (cf Stefan Gieseke).
C
      subroutine initPDFsetC(name)
      implicit none
      byte    name(*)     ! intent(in)
      integer maxlen
      parameter(maxlen=200)
      character*(maxlen) f77name
      integer i,hitend
      hitend = 0
      do i = 1, maxlen
         if (name(i) .eq. 0) then
            hitend = i-1
            exit
         end if
         f77name(i:i) = achar(name(i))
      end do
      if (hitend == 0) then
         write(0,*) "Error in initPDFsetC: name is too long (max",
     $        maxlen," characters)"
         stop
      end if
      call initPDFset(f77name(1:hitend))
      end subroutine
