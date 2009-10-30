        subroutine update_rhoh(nx,scal_old,scal_new,
     $                         beta_old,beta_new,
     $                         aofs,tforce,alpha,Rhs,dx,dt,be_cn_theta)

        implicit none

        include 'nums.fi'
        include 'sndata.fi'

        integer nx

c       Quantities passed in
        real*8  scal_old(-1:nx  ,nscal)
        real*8  scal_new(-1:nx  ,nscal)
        real*8      aofs(0 :nx-1,nscal)
        real*8    tforce(0 :nx-1,nscal)
        real*8     alpha(0 :nx-1)
        real*8  beta_old(-1:nx,nscal)
        real*8  beta_new(-1:nx,nscal)
        real*8       Rhs(0 :nx-1)
        real*8 delta_rhs(0 :nx-1)
        real*8 dx
        real*8 dt
        real*8 be_cn_theta

c       Local variables
        real*8 dth
        real*8 beta_lo,beta_hi
        real*8 visc_term
        real*8 dxsqinv
        real*8 spec_grad_terms(0:nx)
        real*8 rhohalf
        integer i,n

        dth = 0.5d0 * dt
        dxsqinv = 1.d0/(dx*dx)

        n        = RhoH

        be_cn_theta = 1.0d0
c       be_cn_theta = 0.5d0

c*************************************************************************
c       Create delta_rhs.
c*************************************************************************

        call get_spec_grad_terms(nx,scal_old,beta_old,
     $                           spec_grad_terms,dx)

        do i = 0,nx-1
          delta_rhs(i) = -(1.d0-be_cn_theta)*dt*spec_grad_terms(i)
        enddo

c       We assume here that scal_new is already holding good values
c          of (rho Y)_i, rho and Temp.
        call get_spec_grad_terms(nx,scal_new,beta_new,
     $                           spec_grad_terms,dx)

        do i = 0,nx-1
          delta_rhs(i) = delta_rhs(i)-be_cn_theta*dt*spec_grad_terms(i)
        enddo

c*************************************************************************
c       Create Snew = Sold + dt*aofs
c*************************************************************************
        do i = 0,nx-1
          scal_new(i,n) = scal_old(i,n) + dt * aofs(i,n)
     $                                  + dt * tforce(i,n)
        enddo

c*************************************************************************
c       Create RHS = time n diffusive term.
c*************************************************************************

        do i = -1,nx
          scal_old(i,RhoH) = scal_old(i,RhoH) / scal_old(i,Density)
        enddo

        do i = 0,nx-1

          beta_lo = 2.d0 / 
     $      (1.d0/beta_old(i,RhoH)+1.d0/beta_old(i-1,RhoH))
          beta_hi = 2.d0 / 
     $      (1.d0/beta_old(i,RhoH)+1.d0/beta_old(i+1,RhoH))

          visc_term =  
     $       beta_hi*(scal_old(i+1,RhoH)-scal_old(i  ,RhoH)) -
     $       beta_lo*(scal_old(i  ,RhoH)-scal_old(i-1,RhoH))
          Rhs(i) = (1.d0-be_cn_theta) * dt * visc_term * dxsqinv

        enddo

        do i = -1,nx
          scal_old(i,RhoH) = scal_old(i,RhoH) * scal_old(i,Density)
        enddo

c*************************************************************************
c       Construct alpha.
c*************************************************************************

        do i = 0,nx-1
          alpha(i) = scal_new(i,Density)
        enddo

c*************************************************************************
c       Update Rhs.
c*************************************************************************

        do i = 0,nx-1
          Rhs(i) = Rhs(i) + delta_rhs(i)
        enddo

        do i = 0,nx-1
          Rhs(i) = Rhs(i) + scal_new(i,RhoH)
        enddo

        end
