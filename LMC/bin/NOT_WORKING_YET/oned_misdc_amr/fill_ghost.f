      subroutine fill_ghost(scal_old,scal_new,finelev,lo,hi,
     $                      time,curtime,dx,dt)
      implicit none
      include 'spec.h'
      real*8 scal_old(0:nlevs-1,-2:nfine+1)
      real*8 scal_new(0:nlevs-1,-2:nfine+1)
      integer finelev
      integer lo(0:nlevs-1), hi(0:nlevs-1)
      real*8 time, curtime, dx(0:nlevs-1), dt(0:nlevs-1)
      
!     fill fine ghost cells in scal_old by interpolating coarse 
!     data in space and time

!     local
      real*8 scal_c(-2:nfine+1)
      real*8 factor

      integer i,i_c
      real*8 dc,dp,dm,slim,sflag,slpe,x_c,x_f

      factor = (curtime-time) / dt(finelev-1)

!     scal_c is interpolated in time to be consistent with curtime
      scal_c(:) =         factor *scal_new(finelev-1,:) 
     $            + (1.d0-factor)*scal_old(finelev-1,:)

!     ghost cells on low side of the fine grid
      do i=lo(finelev)-2,lo(finelev-1)
!     compute limited second-order slope in underlying coarse cell
         i_c = i / rr
         dc = 0.5d0*( scal_c(i_c+1) - scal_c(i_c-1) )
         dp = 2.d0* ( scal_c(i_c+1) - scal_c(i_c  ) )
         dm = 2.d0* ( scal_c(i_c  ) - scal_c(i_c-1) )
         slim = min(abs(dp), abs(dm))
         sflag = sign(1.d0,dc)
         slpe = sflag*min(slim,abs(dc))
!     fill ghost cell value using conservative linear interpolation
         x_c = (dble(i_c)+0.5d0)*dx(finelev-1)
         x_f = (dble(i  )+0.5d0)*dx(finelev  )
         scal_old(finelev,i) = scal_c(i_c) + slpe*(x_f-x_c)/dx(finelev-1)
      end do

!     ghost cells on the hi side of the fine grid
      do i=hi(finelev)+1,hi(finelev)+2
!     compute limited second-order slope in underlying coarse cell
         i_c = i / rr
         dc = 0.5d0*( scal_c(i_c+1) - scal_c(i_c-1) )
         dp = 2.d0* ( scal_c(i_c+1) - scal_c(i_c  ) )
         dm = 2.d0* ( scal_c(i_c  ) - scal_c(i_c-1) )
         slim = min(abs(dp), abs(dm))
         sflag = sign(1.d0,dc)
         slpe = sflag*min(slim,abs(dc))
!     fill ghost cell value using conservative linear interpolation
         x_c = (dble(i_c)+0.5d0)*dx(1)
         x_f = (dble(i  )+0.5d0)*dx(2)
         scal_old(finelev,i) = scal_c(i_c) + slpe*(x_f-x_c)/dx(1)
      end do

      end
