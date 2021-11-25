! AMAN PAREKH - 180073 - ME685 - Fall 2021

subroutine gaussel(a,b,n,ans)
        
        integer :: n, idmax, idx
        real :: a(n,n), b(n,1)
        real :: ans(n,1), temp, rmax, maxi

        idmax = -1
        idx = -1
        rmax = -1.0
        maxi = 0.0

! Beginining of Forward Elimination        
        do k = 1, n-1
          
            ! Scaling 
            do i = k,n
                if(rmax.lt.abs(a(k,i))) then
                    idmax = i
                end if
            end do

            do j = k,n
                if(maxi.lt.abs(a(j,k)/(abs(a(j,idmax))))) then
                    maxi = a(j,k)
                    idx = j
                end if
            end do
        
            ! Pivoting
            if(idx.ne.k) then
                do i = k,n
                    temp = a(k,i)
                    a(k,i) = a(idx,i)
                    a(idx,i) = temp
                end do
                temp = b(k,1)
                b(k,1) = b(idx,1)
                b(idx,1) = temp
            end if       
                
            ! Elimination
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

real*4 function residual(n, phi, phi_old)

        implicit none
        
        integer :: n, i, j
        real :: phi(n,1), phi_old(n,1), s

        s = 0
        do i = 1,n
          s = s + (abs(phi(i,1)) - abs(phi_old(i,1)))**2
        end do

        residual = sqrt(s)

end function residual

real*4 function norm(A,v,n)

        implicit none
        integer :: n, i, j
        real :: A(n,n), v(n,n), s1, s2

        s2 = 0
        do i = 1,n
          s1 = 0
          do j = 1,n
            s1 = s1 + A(i,j)*v(j,1)
          end do
          s2 = s2 + s1**2
        end do

        norm = sqrt(s2)

end function norm        

subroutine matvecprod(A,v,n,w)

        implicit none
        integer :: n, i, j
        real :: A(n,n), v(n,1), s1, w(n,1)

        do i = 1,n
          s1 = 0
          do j = 1,n
            s1 = s1 + A(i,j)*v(j,1)
          end do
          w(i,1) = s1
        end do

end subroutine matvecprod

subroutine power(A,n,v,l)

        implicit none
        integer :: n, i
        real :: A(n,n), v(n,1), vo(n,1), l, no, w(n,1)
        real*4 :: norm, residual

        do while(residual(n,v,vo).gt.1e-5)
            vo = v
            no = norm(A,vo,n)
            call matvecprod(A,vo,n,w)

            do i = 1,n
              v(i,1) = w(i,1)/no
            end do
        end do

        l = norm(A,v,n)

end subroutine power      

subroutine inverse_power(A,n,v,l)

        implicit none
        integer :: n, i, j
        real :: A(n,n), Atemp(n,n), v(n,1), vo(n,1), votemp(n,1), l, s
        real*4 :: norm, residual 

        vo = 0.0
        votemp = 0.0
        Atemp = A
        do while(residual(n,v,votemp).gt.1e-6)
            vo = v

            ! Have to reallocate the data as subroutine gaussel modifies it
            votemp = v
            A = Atemp

            call gaussel(A,vo,n,v)

            ! Normalizing the Vector
            s = 0
            do i = 1,n
                s = s + v(i,1)**2
            end do

            s = sqrt(s)

            do i = 1,n
                v(i,1) = v(i,1)/s
            end do
        end do

        A = Atemp        
        l = norm(A,v,n)

end subroutine inverse_power

program main

        implicit none

        real*4, external :: check_converge, norm, residual
        real*16, external :: det33
        real :: A(3,3)
        real :: v1(3,1), l1, v2(3,1), l2, ans(3,1)
        integer :: i, j, maxiter

        A = reshape((/ 0.0,0.6,0.0,2.3,0.0,0.3,0.4,0.0,0.0  /), (/ 3,3 /))
        v1 = reshape((/ 1.0,0.0,0.0 /), (/ 3,1 /))
        v2 = v1

        call power(A,3,v1,l1)

        print *, "Maximum EigenValue= ", l1
        print *, "Corresponding EigenVector= ", v1

        call inverse_power(A,3,v2,l2)

        print *, "Minimum EigenValue= ", l2
        print *, "Corresponding EigenVector= ", v2

end program main        
