! AMAN PAREKH - 180073 - ME685 - Summer 2021

subroutine matrixmul(a, b, c, i, j, k)     ! this subroutine can work for any two matrices with sizes i*k and k*j

        integer :: i, j, k              ! i is the number of rows in left-matrix
                                        ! j is the number of columns in right-matrix
                                        ! k is the number of columns in left and number of rows in right matrix
        real :: a(i,k), b(k,j), c(i,j)

        do l = 1,i                      ! iterating through rows in left matrix
            do m = 1,j                  ! iterating through columns in right matrix
                c(l,m) = 0
                do  n = 1,k             ! iterating through elements of row in left and elements of column in right matrix
                    c(l,m) = c(l,m) + a(l,n)*b(n,m)
                end do
            end do
        end do
        
end subroutine matrixmul


program main                            ! this main program is specific to the supplied matrices

        implicit none

        integer :: i,j
        real :: A(3,3), C(3,3), D(3,3), E(3,3)

        read(*,*) ((A(i,j), j=1,3),i=1,3)

        call matrixmul(A,A,C,3,3,3)        ! multiplying C = A*A
        
        print *, "Result(A*A): "
        do i=1,3
            write(*,*) (C(i,j), j=1,3)
        end do

        call matrixmul(C,A,D,3,3,3)        ! multiplying D = A^2 * A

        print *, "Result(A^2 * A): "
        do i=1,3
            write(*,*) (D(i,j), j=1,3)
        end do

        call matrixmul(D,A,E,3,3,3)        ! multiplying E = A^3 * A

        print *, "Result(A^3 * A): "
        do i=1,3
            write(*,*) (E(i,j), j=1,3)
        end do

end program main        
