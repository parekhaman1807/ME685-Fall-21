! AMAN PAREKH - 180073 - ME685 - Monsoon 2021

real function rhs(j,s,i,beta,gamma)

      implicit none

      integer :: j
      real :: s, i, beta, gamma

      if (j.eq.1) then
          rhs = -beta*s*i               ! Equation on Change in S
      else
          rhs = beta*s*i - gamma*i      ! Equation on Change in I
      end if

end function rhs

subroutine SIR(s,i,r,n,beta,gamma,dt,maxi)

      implicit none
        
      integer :: j, k, n
      real :: s(n), i(n), r(n), rk(2,4), temp_s, temp_i
      real :: dt, beta, gamma, rhs, maxi
      open(unit=10, file = 'sirdata.dat', position = 'append')

      write(10,*) s(1), i(1), r(1)

      do j = 1,n        ! Iteration Loop

      ! RK-1 Step
        temp_s = s(j)  
        temp_i = i(j)
        do k = 1,2
          rk(k,1) = rhs(k,temp_s,temp_i,beta,gamma)
        end do

      ! RK-2 Step
        temp_s = s(j) + (dt*(rk(1,1)/2.0))
        temp_i = i(j) + (dt*(rk(2,1)/2.0))
        do k = 1,2
          rk(k,2) = rhs(k,temp_s,temp_i,beta,gamma)
        end do

      ! RK-3 Step
        temp_s = s(j) + (dt*(rk(1,2)/2.0))
        temp_i = i(j) + (dt*(rk(2,2)/2.0))
        do k = 1,2
          rk(k,3) = rhs(k,temp_s,temp_i,beta,gamma)
        end do

      ! RK-4 Step
        temp_s = s(j) + (dt*rk(1,3))
        temp_i = i(j) + (dt*rk(2,3))
        do k = 1,2
          rk(k,4) = rhs(k,temp_s,temp_i,beta,gamma)
        end do

        s(j+1) = s(j) + ((dt*(rk(1,1) + 2.0*rk(1,2) + 2.0*rk(1,3) + rk(1,4)))/6.0)
        i(j+1) = i(j) + ((dt*(rk(2,1) + 2.0*rk(2,2) + 2.0*rk(2,3) + rk(2,4)))/6.0)
        r(j+1) = 1 - s(j+1) - i(j+1)

        write(10,*) s(j+1), i(j+1), r(j+1)

      ! Calculating Maximum I
        if(i(j+1).gt.maxi) then
            maxi = i(j+1)
        end if

     end do

     close(10)

end subroutine SIR

program main

      implicit none

      integer :: j, ndr
      real :: beta, gamma, dt, maxi
      real :: s(100), i(100), r(100)
      real, external :: rhs

      beta = 1.66
      gamma = beta/3.65
      dt = 0.25

      ! Setting initial values
      i(1) = 1.0/763.0
      s(1) = 1.0 - i(1)
      r(1) = 0.0

      maxi = 0.0
      call SIR(s,i,r,100,beta,gamma,dt,maxi)

      ! Calculating the Days Required to 1% of Maximum
      do j = 100,1,-1
        if(i(j).gt.(0.01*maxi)) then
          ndr = j+1
          exit
        end if
      end do

      print '("Maximum Value of I = ",f10.8)', maxi
      print '("Days Required for Infections to Reach 1% of Maximum = ",f6.3)', ndr*dt

end program main
