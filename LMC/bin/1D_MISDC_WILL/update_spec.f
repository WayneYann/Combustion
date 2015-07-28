      subroutine update_spec(scal_old,scal_new,aofs,alpha,beta,
     &                       dRhs,Rhs,dx,dt,be_cn_theta,lo,hi,bc)
      implicit none
      include 'spec.h'
      real*8 scal_old(-2:nfine+1,nscal)
      real*8 scal_new(-2:nfine+1,nscal)
      real*8     aofs(0 :nfine-1,nscal)
      real*8    alpha(0 :nfine-1)
      real*8     beta(-1:nfine  ,nscal)
      real*8      Rhs(0:nfine-1,Nspec)
      real*8     dRhs(0:nfine-1,Nspec)
      real*8 dx,dt,be_cn_theta
      integer lo,hi,bc(2)
      
      real*8 visc(-1:nfine,Nspec)
      real*8  visc_term
      integer i,n,is
      
      real*8 gamma_lo(0:nfine-1,Nspec)
      real*8 gamma_hi(0:nfine-1,Nspec)

      call get_spec_visc_terms(scal_old,beta,visc,
     $                         gamma_lo,gamma_hi,dx,lo,hi)
      do i=lo,hi
         do n=1,Nspec
            is = FirstSpec + n - 1
            visc_term = dt*(1.d0 - be_cn_theta)*visc(i,n)

            Rhs(i,n) = dRhs(i,n) + visc_term + scal_old(i,is)! + dt*aofs(i,is)
            alpha(i) = scal_new(i,Density)
         enddo
      enddo

      call set_bc_s(scal_new,lo,hi,bc)

      end


