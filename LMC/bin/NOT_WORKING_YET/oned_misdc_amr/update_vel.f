      subroutine update_vel(vel_old,vel_new,gp,rhohalf,
     &                      macvel,veledge,alpha,beta,
     &                      Rhs,dx,dt,be_cn_theta,lo,hi)
      implicit none
      include 'spec.h'
      real*8 vel_old(-2:nfine+1)
      real*8 vel_new(-2:nfine+1)
      real*8      gp(-1:nfine)
      real*8 rhohalf(0 :nfine-1)
      real*8  macvel(0 :nfine)
      real*8 veledge(0 :nfine)
      real*8   alpha(0 :nfine-1)
      real*8    beta(-1:nfine  )
      real*8     Rhs(0 :nfine-1)
      real*8 dx,dt,be_cn_theta
      integer lo,hi
      
! local
      real*8  visc(-1:nfine)
      real*8  aofs,visc_term
      integer i


      call get_vel_visc_terms(vel_old,beta,visc,dx,lo,hi)

c     rho.DU/Dt + G(pi) = D(tau), here D(tau) = d/dx ( a . du/dx ), a=4.mu/3
      do i=lo,hi
         visc_term = dt*(1.d0 - be_cn_theta)*visc(i)

         aofs = ( (macvel(i+1)*veledge(i+1) - macvel(i)*veledge(i))  -
     $        0.5d0*(macvel(i+1)-macvel(i))*(veledge(i)+veledge(i+1)) )/dx

         vel_new(i) = vel_old(i) - dt * ( aofs + gp(i)/rhohalf(i) )
         alpha(i) = rhohalf(i)
         Rhs(i) = vel_new(i)*alpha(i) + visc_term
      enddo

      call set_bc_v(vel_new)

      end
      
