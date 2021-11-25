! AMAN PAREKH - 180073 - ME685 - Fall 2021

real*16 function f(x)                   ! Function Definition for f
        implicit none
        real*16 :: x
        f = tan(x) - x
end function f

real*16 function f_prime(x)             ! Function Definition for f_prime
        implicit none
        real*16 :: x
        f_prime = (1/(cos(x)*cos(x))) - 1
end function f_prime

subroutine newrap(f, f_prime, x0, tolerance, maxiter)

        integer :: i, maxiter
        real*16 :: root, x0, tolerance
        real*16 :: f, f_prime

        i = 0
        root = x0                       ! Initially setting Root as Guess

        print '("Value of Root After",i2," iteration= ",f10.8,", and Function Value= ",es14.7)',i,root,abs(f(root))

        do while(abs(f(root)).gt.tolerance .AND. i.lt.maxiter)  ! Checking if Tolerance has been reached
            root = root - (f(root)/f_prime(root))               ! Newton-Raphson Update Step
            i = i + 1                                           ! Counting Number of Iterations
            print '("Value of Root After",i2," iteration= ",f10.8,", and Function Value= ",es14.7)',i,root,abs(f(root))
        end do
       
        print '("Number of Iterations Required= ",i4)',i
        print '("Final Root= ",f10.8,", and Final Function Value= ",es14.7)', root, abs(f(root))

end subroutine newrap
        
program main

        implicit none

        real*16, external :: f, f_prime
        real*16 :: guess, tolerance

        guess = 4.5
        tolerance = 1e-6
        
        print '("Initial Guess= ",f7.5)', guess

        call newrap(f, f_prime, guess, tolerance, 100) 

end program main
