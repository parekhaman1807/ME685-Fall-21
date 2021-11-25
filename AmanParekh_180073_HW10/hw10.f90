! AMAN PAREKH - 180073 - ME685 - Monsoon 2021

real function rms(u, u_old, nx)	! Function to Calculate the RMS Error

      implicit none
      integer :: i, nx
      real :: total, u(101), u_old(101)

      total = 0.0
      do i = 1,nx
        total = total + (u(i) - u_old(i))**2
      end do

      rms = sqrt(total)

end function rms

subroutine TDMA(a,b,c,d,u,nx)		

      implicit none

      integer :: nx, i
      real :: a(nx), b(nx), c(nx), d(nx), u(nx)

      c(1) = c(1)/b(1)
      d(1) = d(1)/b(1)

      ! Forward Elimination
      do i = 2,nx
        c(i) = (c(i))/(b(i) - a(i)*c(i-1))
        d(i) = (d(i) - a(i)*d(i-1))/(b(i) - a(i)*c(i-1))
      end do

      ! Backward Substitution
      u(nx) = d(nx)

      do i = nx-1,1,-1
        u(i) = d(i) - u(i+1)*c(i)
      end do

end subroutine TDMA

program main

      implicit none

      integer :: i, j, n, nx, m
      real :: x(101), u(101), u_old(101), a(101), b(101), c(101), d(101)
      real :: dx, dt, tol
      real, external :: rms
      open(unit=10, file = 'b.dat', position = 'append')

      nx = 101
      dx = 1.0/(nx-1)

      dt = 0.1*(dx**2)
      tol = 1e-8
      m = 100

      do i = 1,nx
        x(i) = (i-1)*dx
      end do
      
      write(10,*) x  

      ! Initializing u
      do i = 1,nx
          u(i) = 0.0
          u_old(i) = 1.0
      end do

      do while(rms(u,u_old,nx) > tol)
        u_old = u
        
       ! Calculating a array
        do j = 2,nx-1
          a(j) = ((-dt/(2.0*(dx**2))) + (dt/(4.0*dx*(x(j)-1))))
        end do

       ! Calculating b array
        do j = 2,nx-1
          b(j) = ((1 + (dt/(dx**2)) - ((m*dt)/(2*(x(j)-1)))))
        end do

       ! Calculating c array
        do j = 2,nx-1
          c(j) = ((-dt/(2.0*(dx**2))) - (dt/(4.0*dx*(x(j)-1))))
        end do

	! Applying the Boundary Conditions in Time
        b(1) = 1.0
        b(nx) = -3 - (c(nx-1)/a(nx-1))

        a(1) = 0.0
        a(nx) = -4 - (b(nx-1)/a(nx-1))

        c(1) = 0.0
        c(nx) = 0.0

        d(1) = 1.0
        d(nx) = 0.0
 
       ! Calculating d array
        do j = 2,nx-1
          d(j) = u_old(j-1)*((dt/(2.0*(dx**2))) - (dt/(4.0*dx*(x(j)-1)))) &
               + u_old(j)*(1 - (dt/(dx**2)) + ((m*dt)/(2*(x(j)-1)))) &
               + u_old(j+1)*((dt/(2.0*(dx**2))) + (dt/(4.0*dx*(x(j)-1))))
        end do

        call TDMA(a,b,c,d,u,nx)

        write(10,*) u
      end do

      close(10)

end program main

