      subroutine update_rhoh(scal_old,scal_new,aofs,alpha,
     &                       beta,dRhs,Rhs,dx,dt,be_cn_theta,lo,hi,bc)
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

      real*8  visc(-1:nfine)
      real*8  visc_term
      integer i

      call get_rhoh_visc_terms(scal_old,beta,visc,dx,lo,hi)
      
      do i=lo,hi
         visc_term = dt*(1.d0 - be_cn_theta)*visc(i)
         
         Rhs(i) = dRhs(i) + scal_old(i,RhoH) + visc_term! + dt*aofs(i,RhoH)
         alpha(i) = scal_new(i,Density)
      enddo
      
      call set_bc_s(scal_new,lo,hi,bc)

      end
