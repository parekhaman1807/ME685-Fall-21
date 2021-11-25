! AMAN PAREKH - 180073 - ME685 - Odd 2021

real*8 function den(x,y)

      implicit none
      real :: x,y

      den = abs(x+y) + 1

end function den

subroutine trapz(meshx,meshy,n,auc,hy,flag)

      implicit none
      integer :: n, i, j, flag
      real :: meshx(n+1,n+1), meshy(n+1,n+1), auc, hy(n+1)
      real*8 :: den

      if(flag.eq.1) then                ! For calculating Volume
          do i = 1,n+1                  ! Outer Loop of integral (x-variable)
            auc = auc + ((1.0*hy(i))/2.0)
            do j = 2,n                  ! Inner Loop of Integral (y-variable)
              auc = auc + (1.0*hy(i))
            end do
            auc = auc + ((1.0*hy(i))/2.0)
          end do
      else if(flag.eq.2) then           ! For calculating Mass
          do i = 1,n+1
            auc = auc + ((den(meshx(i,1),meshy(i,1))*hy(i))/2.0)
            do j = 2,n
              auc = auc + (den(meshx(i,j),meshy(i,j))*hy(i))
            end do
            auc = auc + ((den(meshx(i,n+1),meshy(i,n+1))*hy(i))/2.0)
          end do
      else if(flag.eq.3) then           ! For calculating X-COM
           do i = 1,n+1
            auc = auc + ((den(meshx(i,1),meshy(i,1))*meshx(i,1)*hy(i))/2.0)
            do j = 2,n
              auc = auc + (den(meshx(i,j),meshy(i,j))*meshx(i,j)*hy(i))
            end do
            auc = auc + ((den(meshx(i,n+1),meshy(i,n+1))*meshx(i,n+1)*hy(i))/2.0)
          end do
      else                              ! For calculating Y-COM
            do i = 1,n+1
            auc = auc + ((den(meshx(i,1),meshy(i,1))*meshy(i,1)*hy(i))/2.0)
            do j = 2,n
              auc = auc + (den(meshx(i,j),meshy(i,j))*meshy(i,j)*hy(i))
            end do
            auc = auc + ((den(meshx(i,n+1),meshy(i,n+1))*meshy(i,n+1)*hy(i))/2.0)
          end do
      end if
        
end subroutine trapz      

program main

      implicit none
      real*8, external :: den
      integer :: n, i, j
      real :: meshx(401,401), meshy(401,401), x, y, hy(401)
      real :: auc, vol, mass

      ! Mesh Generation
      do i = 1,401
        x = -1 + (0.005*(i-1))
        hy(i) = (sqrt(1 - x**2)+sqrt(1 - x**2))/400.0
        do j = 1,401
          y = -sqrt(1 - x **2) + hy(i)*(j-1)
          meshx(i,j) = x
          meshy(i,j) = y
        end do
      end do

      auc = 0.0
      call trapz(meshx,meshy,400,auc,hy,1)

      vol = auc

      auc = 0.0
      call trapz(meshx,meshy,400,auc,hy,2)

      mass = auc

      print *,"Average Density = ", mass/vol

      auc = 0.0
      call trapz(meshx,meshy,400,auc,hy,3)

      print *,"X Coordinate = ", auc/mass

      auc = 0.0
      call trapz(meshx,meshy,400,auc,hy,4)

      print *,"Y Coordinate = ", auc/mass

end program main      
