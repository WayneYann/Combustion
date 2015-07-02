      subroutine update_rho(scal_old,scal_new,aofs,dt,lo,hi,bc)
      implicit none
      include 'spec.h'
      real*8 scal_old(-2:nfine+1,nscal)
      real*8 scal_new(-2:nfine+1,nscal)
      real*8     aofs(0 :nfine-1,nscal)
      real*8 dt
      integer lo,hi,bc(2)

      integer i

      do i=lo,hi
         scal_new(i,Density) = scal_old(i,Density) + dt*aofs(i,Density)
      enddo
      
      call set_bc_s(scal_new,lo,hi,bc)

      end
