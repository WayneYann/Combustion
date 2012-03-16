      subroutine update_rhoh(scal_old,scal_new,aofs,alpha,
     &                       beta,dRhs,Rhs,dx,dt,be_cn_theta)
      implicit none
      include 'spec.h'
      real*8 scal_old(-2:nx+1,nscal)
      real*8 scal_new(-2:nx+1,nscal)
      real*8     aofs(0 :nx-1,nscal)
      real*8    alpha(0 :nx-1)
      real*8     beta(-1:nx  ,nscal)
      real*8     dRhs(0 :nx-1)
      real*8      Rhs(0 :nx-1)
      real*8 dx,dt,be_cn_theta

      real*8  visc(-1:nx)
      real*8  visc_term
      integer i

      call get_rhoh_visc_terms(scal_old,beta,visc,dx)
      
      do i = 0,nx-1
         visc_term = dt*(1.d0 - be_cn_theta)*visc(i)
         
         scal_new(i,RhoH) = scal_old(i,RhoH) + dt*aofs(i,RhoH)
         Rhs(i) = dRhs(i) + scal_new(i,RhoH) + visc_term
         alpha(i) = scal_new(i,Density)
      enddo
      
      call set_bc_s(scal_new)

      end
