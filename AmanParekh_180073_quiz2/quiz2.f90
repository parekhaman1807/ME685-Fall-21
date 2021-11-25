! AMAN PAREKH - Quiz 2 - ME685

real function rhs(j,y,u,w,n)

      implicit none

      integer :: j, n
      real :: y, u, w

      if (j.eq.1) then
          rhs = y  ! equation for y             
      else(j.eq.2) then
          rhs = u  ! equation for y prime
      else
          rhs = -((n+1)/2)*y*w + n*u**2 - n     ! equation for y double prime
      end if

end function rhs

subroutine rk4(y,u,w,dt,n)

      implicit none
        
      integer :: j, k, n
      real :: y(n), u(n), w(n), rk(3,4)
      real :: dt
      open(unit=10, file = 'rk4.dat', position = 'append')

      write (10,*) y(1), u(1), w(1)

      do j = 1,n        ! Iteration Loop

      ! RK-1 Step
        do k = 1,3
          rk(k,1) = rhs(k,y(j),u(j),w(j),n)
        end do


      ! RK-2 Step
        do k = 1,3
          rk(k,2) = rhs(k,y(j)+dt*rk(1,1)/2,u(j)+dt*rk(2,1)/2,w(j) + dt*rk(3,1)/2,n)
        end do

      ! RK-3 Step
        do k = 1,3
          rk(k,3) = rhs(k,y(j)+dt*rk(1,2)/2,u(j)+dt*rk(2,2)/2,w(j) + dt*rk(3,2)/2,n)
        end do

      ! RK-4 Step
        do k = 1,3
          rk(k,4) = rhs(k,y(j)+dt*rk(1,3),u(j)+dt*rk(2,3),w(j) + dt*rk(3,3),n)
        end do

        y(j+1) = y(j) + ((dt*(rk(1,1) + 2.0*rk(1,2) + 2.0*rk(1,3) + rk(1,4)))/6.0)
        u(j+1) = u(j) + ((dt*(rk(2,1) + 2.0*rk(2,2) + 2.0*rk(2,3) + rk(2,4)))/6.0)
        w(j+1) = w(j) + ((dt*(rk(3,1) + 2.0*rk(3,2) + 2.0*rk(3,3) + rk(3,4)))/6.0)

        write(10,*) y(j), u(j), w(j)

     end do
     
     close(10)

end subroutine SIR


program main

      implicit none

      n = -0.09
      dt = 0.1

      y(1) = 0.0
      u(1) = 0.0
      w(1) = 0.0 ! guess value

      ! Code for RK4 is written but shooting method could not be applied

      

end program main

