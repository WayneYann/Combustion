      subroutine ctoprim(lo,hi,  uin, q, 
                         courno,dx,dy,dz,ng,gamma,cv)
      !
      !     Will give primitive variables on lo-ngp:hi+ngp, and flatn on lo-ngf:hi+ngf
      !     if iflaten=1.  Declared dimensions of q,c,gamc,csml,flatn are given
      !     by DIMS(q).  This declared region is assumed to encompass lo-ngp:hi+ngp.
      !     Also, uflaten call assumes ngp>=ngf+3 (ie, primitve data is used by the
      !     routine that computes flatn).  
      !
      implicit none

      integer lo(3), hi(3), ng

      double precision :: uin(lo(1)-ng:hi(1)+ng,lo(2)-ng:hi(2)+ng,lo(3)-ng:hi(3)+ng,5)
      double precision :: q(lo(1)-ng:hi(1)+ng,lo(2)-ng:hi(2)+ng,lo(3)-ng:hi(3)+ng,6)
      double precision :: dx, dy, dz, dt, courno, gamma, cv

      integer          :: i, j, k
      double precision :: courx, coury, courz, courmx, courmy, courmz

      double precision, parameter :: gamma = 1.4d0
      double precision, parameter :: cv = 8.3333333333d6



      !$OMP PARALLEL DO PRIVATE(i,j,k)
      do k = lo(3)-ng,hi(3)+ng
         do j = lo(2)-ng,hi(2)+ng
            do i = lo(1)-ng,hi(1)+ng

               q(i,j,k,1) = uin(i,j,k,1)
               q(i,j,k,2) = uin(i,j,k,2)/uin(i,j,k,1)
               q(i,j,k,3) = uin(i,j,k,3)/uin(i,j,k,1)
               q(i,j,k,4) = uin(i,j,k,4)/uin(i,j,k,1)
            enddo
         enddo
      enddo
      !$OMP END PARALLEL DO

      ! Get gamc, p, T, c, csml using q state
      !$OMP PARALLEL DO PRIVATE(i,j,k,pt_index)
      do k = lo(3)-ng,hi(3)+ng
         do j = lo(2)-ng,hi(2)+ng
            do i = lo(1)-ng,hi(1)+ng

!  should set up call to eos
!  simplify to simple things

            eint =  uin(i,j,k,5)/uin(i,j,k,1)      & 
               - 0.5d0*( q(i,j,k,2)**2 + q(i,j,k,3)**2 + q(i,j,k,4) **2)
            q(i,j,k,5) = (gamma-1.d0)*eint*uin(i,j,k,1)
            q(i,j,k,6) = eint/cv

            end do
         end do
      end do
      !$OMP END PARALLEL DO

      ! Compute running max of Courant number over grids
      courmx = courno
      courmy = courno
      courmz = courno

      !$OMP PARALLEL DO PRIVATE(i,j,k,courx,coury,courz) REDUCTION(max:courmx,courmy,courmz)
      do k = lo(3),hi(3)
         do j = lo(2),hi(2)
            do i = lo(1),hi(1)

               c = dsqrt(gamma*q(i,j,k,5)/q(i,j,k,1)
               courx = ( c+abs(q(i,j,k,2)) ) * dx
               coury = ( c+abs(q(i,j,k,3)) ) * dy
               courz = ( c+abs(q(i,j,k,4)) ) * dz

               courmx = max( courmx, courx )
               courmy = max( courmy, coury )
               courmz = max( courmz, courz )


            enddo
         enddo
      enddo
      !$OMP END PARALLEL DO

      courno = max( courmx, courmy, courmz )


      end subroutine ctoprim
