! AMAN PAREKH - 180073 - ME685 - Fall 2021

real*4 function check_converge(ans, ans_old,n)          ! function to check if solution has converged

        implicit none        
        integer :: n, i
        real :: ans(n,1), ans_old(n,1), tolerance

        check_converge = 1.0                            ! Assuming it is converged
        tolerance = 1e-6

        do i = 1,n
            if (abs(ans(i,1) - ans_old(i,1)).gt.tolerance) then
                    check_converge = 0.0                    ! Setting check_converge as false if not converged
            endif
        end do

end function check_converge       

subroutine seidel(a,b,n,ans,maxiter)

        integer :: n, i, maxiter, k, j
        real :: a(n,n), b(n,1)
        real :: ans(n,1), ans_old(n,1)
        real :: r_sum

        i = 0
        ans_old =  reshape((/ 1,1,1,1 /), (/ 4,1 /))       ! random value of ans_old so that first check is passed

        do while(check_converge(ans, ans_old,n).eq.0.0 .AND. i.lt.maxiter)      
            ans_old = ans 
            do k = 1,n                                     ! iterating through the rows
                r_sum = 0
                do j = 1,n                                 ! ietrating through columns of the row
                if (j.lt.k) then                           ! if j is less than k, then adding updated values
                        r_sum = r_sum + (a(k,j)*ans(j,1))
                    else                                   ! otherwise adding previous iteration values
                        r_sum = r_sum + (a(k,j)*ans_old(j,1))
                    end if                        
                end do                    
                ans(k,1) = ans(k,1) + ((b(k,1) - r_sum)/a(k,k))
            end do
            i = i + 1
        end do

        print '("Number of Iterations Required= ", i3)', i

end subroutine seidel

program main

        implicit none

        real*4, external :: check_converge
        real :: A(4,4), B(4,1)
        real :: ans(4,1), guess(4,1)
        integer :: i, j, maxiter

        A = reshape((/ 4,1,1,1,1,4,1,1,1,1,4,1,1,1,1,4 /), (/ 4,4 /))
        B = reshape((/7,7,7,7/), (/4,1/))

        maxiter = 100
        guess = reshape((/ 0,0,0,0 /), (/ 4,1 /))       ! Initial Guess
        ans = guess

        print *, "A Matrix: "
        do i=1,4
            write(*,*) (A(i,j), j=1,4)
        end do

        print *, "B Matrix: "
        do i=1,4
            write(*,*) B(i,1)
        end do

        call seidel(A,B,4,ans,maxiter)

        print *, "Result Matrix: "
        do i=1,4
            write(*,*) ans(i,1)
        end do
        
end program main        
