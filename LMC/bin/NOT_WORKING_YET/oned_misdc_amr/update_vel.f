      subroutine update_vel(vel_old,vel_new,gp,rhohalf,
     &                      macvel,sedge,alpha,beta,
     &                      Rhs,dx,dt,be_cn_theta,time)
      implicit none
      include 'spec.h'
      real*8 vel_old(-1:nx  )
      real*8 vel_new(-1:nx  )
      real*8      gp(0 :nx-1)
      real*8   alpha(0 :nx-1)
      real*8    beta(-1:nx  )
      real*8     Rhs(0 :nx-1)
      real*8  rhohalf(0 :nx-1)
      real*8   macvel(0 :nx)
      real*8    sedge(0 :nx)
      real*8 dx,dt,be_cn_theta,time
      
      real*8  visc(0:nx-1)
      real*8  aofs,visc_term
      integer i


      call get_vel_visc_terms(vel_old,beta,visc,dx,time)

c     rho.DU/Dt + G(pi) = D(tau), here D(tau) = d/dx ( a . du/dx ), a=4.mu/3
      do i = 0,nx-1
         visc_term = dt*(1.d0 - be_cn_theta)*visc(i)

         aofs = ( (macvel(i+1)*sedge(i+1) - macvel(i)*sedge(i))  -
     $        0.5d0*(macvel(i+1)-macvel(i))*(sedge(i)+sedge(i+1)) )/dx

         vel_new(i) = vel_old(i) - dt * ( aofs + gp(i)/rhohalf(i) )
         alpha(i) = rhohalf(i)
         Rhs(i) = vel_new(i)*alpha(i) + visc_term
      enddo

      call set_bc_v(vel_new,dx,time)

      end
      
