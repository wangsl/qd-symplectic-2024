c /* $Id$ */

      subroutine rvadd(c, a, b, n)
      implicit none
      integer n, i
      Real*8 c(n), a(n), b(n)
      do i = 1, n
         c(i) = a(i) + b(i)
      end do
      end

      subroutine rvsubt(c, a, b, n)
      implicit none
      integer n, i
      Real*8 c(n), a(n), b(n)
      do i = 1, n
         c(i) = a(i) - b(i)
      end do
      end
      
      subroutine rsvmult(c, a, b, n)
      implicit none
      integer n, i
      Real*8 c(n), a, b(n)
      do i = 1, n
         c(i) = a * b(i)
      end do
      end

      subroutine rvvmult(c, a, b, n)
      implicit none
      integer n, i
      Real*8 c, a(n), b(n)
      c = 0.0d0
      do i = 1, n
         c = c + a(i) * b(i)
      end do
      end
      
      subroutine rmadd(c, a, b, n, m)
      implicit none
      integer n, m, i, j
      Real*8 c(n,m), a(n,m), b(n,m)
      do i = 1, n
         do j = 1, m
            c(i,j) = a(i,j) + b(i,j)
         end do
      end do
      end

      subroutine rmsubt(c, a, b, n, m)
      implicit none
      integer n, m, i, j
      Real*8 c(n,m), a(n,m), b(n,m)
      do i = 1, n
         do j = 1, m
            c(i,j) = a(i,j) - b(i,j)
         end do
      end do
      end      

      subroutine rsmmult(c, a, b, n, m)
      implicit none
      integer n, m, i, j
      Real*8 c(n,m), a, b(n,m)
      do i = 1, n
         do j = 1, m
            c(i,j) = a * b(i,j)
         end do
      end do
      end

      subroutine rvmmult(c, a, b, n, m)
      implicit none
      integer n, m, i, j
      Real*8 c(m), a(n), b(n,m), x
      do j = 1, m
         x = 0.0d0
         do i = 1, n
            x = x + a(i) * b(i,j)
         end do
         c(j) = x
      end do
      end

      subroutine rmvmult(c, a, b, n, m)
      implicit none
      integer n, m, i, j
      Real*8 c(n), a(n,m), b(m), x
      do i = 1, n
         x = 0.0d0
         do j = 1, m
            x = x + a(i,j) * b(j)
         end do
         c(i) = x
      end do
      end
      
      subroutine rmmmult(c, a, b, n, m, l)
      implicit none
      integer n, m, l, i, j, k
      Real*8 c(n,m), a(n,l), b(l,m), x
      do i = 1, n
         do j = 1, m
            x = 0.0d0
            do k = 1, l
               x = x + a(i,k) * b(k,j)
            end do
            c(i,j) = x
         end do
      end do
      end      
      
c      subroutine rvmvmult(c, a, b, n)
c      implicit none
c      integer n, i, j
c      double precision c, a(n,n), b(n)
c      c = 0.0d0
c      do i = 1, n
c         c = c + a(i,i) * b(i)**2
c         do j = i+1, n
c            c = c + 2 * a(i,j) * b(i) * b(j)
c         end do
c      end do
c      end

c modified by Shenglong Wang 1-10-2005

      subroutine rvmvmult(c, a, b, n)
      implicit none
      integer n, i, j
      Real*8 c, a(n,n), b(n)
      Real*8 d
      c = 0.0d0
      do i = 1, n
         c = c + a(i,i) * b(i)**2
         d = 0.0d0
         do j = i+1, n
            d = d + a(i,j)*b(j)
         end do
         d = 2.0d0*b(i)*d
         c = c + d
      end do
      end
      
