      subroutine update_rhoh(scal_old,scal_new,aofs,alpha,
     &                       dRhs,Rhs,dt,lo,hi,bc)
      implicit none
      include 'spec.h'
      real*8 scal_old(-2:nfine+1,nscal)
      real*8 scal_new(-2:nfine+1,nscal)
      real*8     aofs(0 :nfine-1,nscal)
      real*8    alpha(0 :nfine-1)
      real*8     beta(-1:nfine  ,nscal)
      real*8     dRhs(0 :nfine-1)
      real*8      Rhs(0 :nfine-1)
      real*8 dx,dt,be_cn_theta
      integer lo,hi,bc(2)

      integer i

      do i=lo,hi
         scal_new(i,RhoH) = scal_old(i,RhoH) + dt*aofs(i,RhoH)
         Rhs(i) = dRhs(i) + scal_new(i,RhoH)
         alpha(i) = scal_new(i,Density)
      enddo
      
      call set_bc_s(scal_new,lo,hi,bc)

      end
