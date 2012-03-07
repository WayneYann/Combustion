      subroutine update_spec(scal_old,scal_new,aofs,
     &     alpha,beta,dRhs,Rhs,dx,dt,be_cn_theta,time)
      implicit none
      include 'spec.h'
      real*8 scal_old(-2:nx+1,nscal)
      real*8 scal_new(-2:nx+1,nscal)
      real*8     aofs(0 :nx-1,nscal)
      real*8    alpha(0 :nx-1)
      real*8     beta(-1:nx  ,nscal)
      real*8      Rhs(0 :nx-1,*)
      real*8     dRhs(0:nx-1,1:maxspec)
      real*8 dx,dt,be_cn_theta,time
      
      real*8 visc(0 :nx-1,Nspec)
      real*8  visc_term
      integer i,n,is
      
      real*8 spec_flux_lo(0:nx-1,maxspec)
      real*8 spec_flux_hi(0:nx-1,maxspec)

      call get_spec_visc_terms(scal_old,beta,visc,
     $                         spec_flux_lo,spec_flux_hi,dx)
      do i = 0,nx-1
         do n=1,Nspec
            is = FirstSpec + n - 1
            visc_term = dt*(1.d0 - be_cn_theta)*visc(i,n)

            scal_new(i,is) = scal_old(i,is) + dt*aofs(i,is)
            Rhs(i,n) = dRhs(i,n) + visc_term + scal_new(i,is)
            alpha(i) = scal_new(i,Density)
         enddo
      enddo

      call set_bc_s(scal_new,dx,time)

      end


