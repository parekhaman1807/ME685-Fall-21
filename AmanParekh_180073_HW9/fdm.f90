! AMAN PAREKH - 180073 - ME685 - Monsoon 2021

subroutine TDMA(theta,a,b,c,d,n)

      implicit none

      integer :: n, i
      real :: theta(n), a(n), b(n), c(n), d(n)

      ! Forward Elimination

      do i = 2,n
        b(i) = b(i) - (c(i-1)*(a(i)/b(i-1)))
        d(i) = d(i) - (d(i-1)*(a(i)/b(i-1)))
      end do

      ! Backward Substitution

      theta(n) = d(n)/b(n)

      do i = n-1, 1, -1
        theta(i) = (d(i) - (c(i)*theta(i+1)))/b(i)
      end do

end subroutine TDMA

program main

      implicit none

      integer :: m, i, n
      real :: dz, z(101), a(101), b(101), c(101), d(101)
      real :: theta(101)
      open(unit=10, file = 'fdm.dat', position = 'append')

      m = 50
      dz = 0.01
      n = (1.0/dz) + 1

      z(1) = 0.0

      do i = 2,n
        z(i) = ((i)*(1.0))/(n*1.0)
      end do

      ! Initializing the a and b arrays
      b(1) = 1.0
      b(n) = 1.0
      a(n) = 0.0

      do i = 2,n-1
        a(i) = 1.0
        b(i) = -2.0 - (m**2)*(dz**2)*(1.0)
      end do

      ! Initializing the c and d arrays
      c(1) = 0.0
      d(1) = 1.0
      d(n) = 0.0

      do i = 2,n-1
        c(i) = 1.0
        d(i) = 0.0
      end do

      call TDMA(theta,a,b,c,d,n)

      do i = 1,n
        write(10,*) z(i), theta(i)
      end do

end program main      
