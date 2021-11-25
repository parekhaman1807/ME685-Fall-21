! AMAN PAREKH - 180073 - QUIZ 1 - ME685 - Fall 2021

subroutine gaussel(a,b,n,ans)   ! Here a is the Jacobian matrix, b is the negative of function value function value, n is 7, and ans
!                                       is del(x)

        integer :: n
        real :: a(n,n), b(n,1)
        real :: ans(n,1)

! Beginining of Forward Elimination        
        do k = 1, n-1
            do i = k+1, n
                factor = a(i,k)/a(k,k)
                do j = k, n
                    a(i,j) = a(i,j) - a(k,j)*factor
                end do
                b(i,1) = b(i,1) - b(k,1)*factor
            end do
        end do
! End of Forward Elimination

! Beginning of Back Substituion
        ans(n,1) = b(n,1)/a(n,n)
        do i = n-1,1,-1
            sum = b(i,1)
            do j = i+1, n
                sum = sum - a(i,j)*ans(j,1)
            end do
            ans(i,1) = sum/a(i,i)
        end do

end subroutine gaussel

! Please note the code is code is commented because it is not functional, but all steps are highlighted

!program main

!        implicit none
        
!       ! All Definition of functions for the matrix entries are defined
!          J11 = 
!          J12 = 
!          J21 = 
!          etc...
!          f1 = 
!          f2 = 
!          etc..


!        guess =  will be a 7x1 matrix
!        xk =  will be a 7x1 vector
!        xkp =  will be a 7x1 vector

!        J = ! will be a 7x7 matrix with functions that calculate values


!       Next is the Newton-Raphson Loop
!        do while(convergence)  convergence is a function that checks for convergence of all variables
!                xk = xkp
!                 calculate J
!                 calculate f
!                 call gaussel(J,f,7,delx)      ! Call for Gaussian Elimination for solving for delx
!                 xkp = xk + delx
!        end do

!end program main
