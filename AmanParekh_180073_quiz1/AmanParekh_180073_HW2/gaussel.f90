! AMAN PAREKH - 180073 - ME685 - Summer 2021

subroutine gaussel(a,b,n,ans)					! Subrotuine that will find ans for any suppiled square matrix a and column vector b of size n
        
        integer :: n
        real :: a(n,n), b(n,1)
        real :: ans(n,1)

! Beginining of Forward Elimination        
        do k = 1, n-1						! Iterating through the rows
            do i = k+1, n					! Selecting the row we are modifying
                factor = a(i,k)/a(k,k)			! Calculating the factor which depends on pivot element and leading element of the row
                do j = k, n					! Iterating through the columns of the selected row
                    a(i,j) = a(i,j) - a(k,j)*factor		! Modifying the column entries
                end do
                b(i,1) = b(i,1) - b(k,1)*factor		! Modifying the entry in the B matrix
            end do
        end do
! End of Forward Elimination

! Beginning of Back Substituion
        ans(n,1) = b(n,1)/a(n,n)				! Calculating from the last row which has only one unknown
        do i = n-1,1,-1					! Iterating through rows in reverse manner as number of unknowns increase
            sum = b(i,1)					
            do j = i+1, n					
                sum = sum - a(i,j)*ans(j,1)			! Calculating the sum with the help of solved unknowns
            end do
            ans(i,1) = sum/a(i,i)				! Finding the unknown  of the i'th row
        end do
! End of Back Substitution        

end subroutine gaussel      

program main							! This main function is written for the current problem

        implicit none

        real :: A(4,4), B(4,1)
        real :: ans(4,1)
        integer :: i, j

        A = reshape((/ 4,1,1,1,1,4,1,1,1,1,4,1,1,1,1,4 /), (/ 4,4 /))	! PLease note fortran takes entries in column-wise format
        B = reshape((/7,7,7,7/), (/4,1/))
        
        print *, "A Matrix: "
        do i=1,4
            write(*,*) (A(i,j), j=1,4)			! Printing the value original A matrix
        end do

        print *, "B Matrix: "
        do i=1,4
            write(*,*) B(i,1)					! Printing the value original B matrix
        end do
 
        call gaussel(A,B,4,ans)         			! Calling the Gauss Elimination Subroutine

        print *, "Modified A Matrix: "
        do i=1,4
            write(*,*) (A(i,j), j=1,4)			! Printing the value modified A matrix
        end do

        print *, "Modified B Matrix: "
        do i=1,4
            write(*,*) B(i,1)					! Printing the value modified B matrix
        end do

        print *, "X (Solution) Matrix: "
        do i=1,4
            write(*,*) ans(i,1)				! Printing the solution matrix
        end do

end program main        
