! AMAN PAREKH - 180073 - ME685 - Odd 2021

subroutine laginter(x, y, xp, yl, n, datap)

      implicit none
      integer :: n, datap, i, j, k
      real :: x(n), y(n), xp(datap), yl(datap), l
      open(unit=10, file = 'interpolated.dat', position = 'append')

      do k = 1,datap
        do i = 1,n
          l = 1.0
          do j = 1,n     
            if(i.ne.j) then 
              l = l*((xp(k) - x(j))/(x(i) - x(j)))
            end if
          end do
          yl(k) = yl(k) + l*y(i)
        end do
        write (10,*) xp(k), yl(k)
      end do

      close(10)

end subroutine laginter

program main

      implicit none
      real :: x(8), y(8), pi, maxerr, meanerr
      integer :: i, datap
      real :: xp(100), yl(100), ya(100)
      open(unit=11, file = 'analytical.dat', position = 'append')
      open(unit=12, file = 'given.dat', position = 'append')
      
      datap = 100
      pi = 3.14159265
      x(1) = 0
      y(1) = 0
      write(12,*) x(1), y(1)

      do i = 1,7
        x(i+1) = ((i*(2.0*pi))/7.0)
        y(i+1) = (sin(x(i+1)))**2
        write(12,*) x(i+1), y(i+1)
      end do

      xp(1) = 0
      do i = 1,(datap-1)
        xp(i+1) = ((i*(2.0*pi))/(datap-1))
      end do

      yl = 0.0
      call laginter(x, y, xp, yl, 8, datap)

      do i = 1,datap
        ya(i) = sin(xp(i))**2
        write(11,*) xp(i), ya(i)
      end do

      ! Error Calculation
      
      meanerr = 0.0
      maxerr = 0.0

      do i = 1,datap
        meanerr = meanerr + abs(yl(i) - ya(i))
        if(abs(yl(i) - ya(i)).gt.maxerr) then
          maxerr = abs(yl(i) - ya(i))
        end if
      end do

      meanerr = meanerr/datap

      print *, 'Mean Error= ', meanerr
      print *, 'Maximum Error= ', maxerr

      close(11)
      close(12)

end program main
