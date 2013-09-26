subroutine rns_makeplotvar(lo, hi, dx, &
     prim, pm_l1, pm_l2, pm_l3, pm_h1, pm_h2, pm_h3, &
     plot, pt_l1, pt_l2, pt_l3, pt_h1, pt_h2, pt_h3, &
     npv, icomp_magvel, icomp_Mach, icomp_divu, icomp_magvort, &
     icomp_X, icomp_omegadot, icomp_dYdt, icomp_heatRelease, icomp_fuelConsumption, &
     fuelID)
  use meth_params_module, only : NVAR, QRHO, QU, QV, QW, QTEMP, QFY, NSPEC
  use eos_module, only : eos_get_c
  implicit none

  double precision, intent(in) :: dx(3)
  integer, intent(in) :: lo(3), hi(3)
  integer, intent(in) :: pm_l1, pm_h1, pm_l2, pm_h2, pm_l3, pm_h3
  integer, intent(in) :: pt_l1, pt_h1, pt_l2, pt_h2, pt_l3, pt_h3
  integer, intent(in) :: npv, icomp_magvel, icomp_Mach, icomp_divu, icomp_magvort, &
       icomp_X, icomp_omegadot,	icomp_dYdt, icomp_heatRelease, icomp_fuelConsumption, &
       fuelID
  double precision,intent(in   )::prim(pm_l1:pm_h1,pm_l2:pm_h2,pm_l3:pm_h3,NVAR) 
  double precision,intent(inout)::plot(pt_l1:pt_h1,pt_l2:pt_h2,pt_l3:pt_h3,0:npv-1) 
    
  integer :: i, j, k
  double precision :: Yt(NSPEC), v, c, dxinv(3), uy,uz,vx,vz,wx,wy

  !$omp parallel private(i,j,k,Yt,v,c,dxinv,uy,uz,vx,vz,wx,wy)

  if (icomp_magvel .ge. 0) then
     !$omp do collapse(2)
     do k=lo(3), hi(3)
        do j=lo(2), hi(2)
           do i=lo(1), hi(1)
              plot(i,j,k,icomp_magvel) = sqrt(prim(i,j,k,QU)**2 &
                   +                          prim(i,j,k,QV)**2 &
                   +                          prim(i,j,k,QW)**2 )
           end do
        end do
     end do
     !$omp end do nowait
  end if

  if (icomp_Mach .ge. 0) then
     !$omp do collapse(2)
     do k=lo(3), hi(3)
        do j=lo(2), hi(2)
           do i=lo(1), hi(1)
              Yt = prim(i,j,k,QFY:QFY+NSPEC-1)
              call eos_get_c(c,prim(i,j,k,QRHO),prim(i,j,k,QTEMP),Yt)
              if (icomp_magvel .ge. 0) then
                 v = plot(i,j,k,icomp_magvel)
              else
                 v = sqrt(prim(i,j,k,QU)**2+prim(i,j,k,QV)**2+prim(i,j,k,QW)**2)
              end if
              plot(i,j,k,icomp_Mach) = v / c
           end do
        end do
     end do
     !$omp end do nowait
  end if

  if (icomp_divu .ge. 0) then
     dxinv = 1.d0/dx
     !$omp do collapse(2)
     do k=lo(3), hi(3)
        do j=lo(2), hi(2)
           do i=lo(1),hi(1)
              plot(i,j,k,icomp_divu) = 0.5d0 * &
                   ( (prim(i+1,j,k,QU)-prim(i-1,j,k,QU))*dxinv(1) &
                   + (prim(i,j+1,k,QV)-prim(i,j-1,k,QV))*dxinv(2) &
                   + (prim(i,j,k+1,QW)-prim(i,j,k-1,QW))*dxinv(3) )
           end do
        end do
     end do
     !$omp end do nowait
  end if

  if (icomp_magvort .ge. 0) then
     dxinv = 1.d0/dx
     !$omp do collapse(2)
     do k=lo(3), hi(3)
        do j=lo(2), hi(2)
           do i=lo(1),hi(1)
              uy = 0.5d0*dxinv(2)*(-prim(i,j-1,k,QU)+prim(i,j+1,k,QU))
              uz = 0.5d0*dxinv(3)*(-prim(i,j,k-1,QU)+prim(i,j,k+1,QU))
              vx = 0.5d0*dxinv(1)*(-prim(i-1,j,k,QV)+prim(i+1,j,k,QV))
              vz = 0.5d0*dxinv(3)*(-prim(i,j,k-1,QV)+prim(i,j,k+1,QV))
              wx = 0.5d0*dxinv(1)*(-prim(i-1,j,k,QW)+prim(i+1,j,k,QW))
              wy = 0.5d0*dxinv(2)*(-prim(i,j-1,k,QW)+prim(i,j+1,k,QW))
              plot(i,j,k,icomp_magvort) = sqrt((wy-vz)**2+(uz-wx)**2+(vx-uy)**2)
           end do
        end do
     end do
     !$omp end do
  end if

  !$omp end parallel

end subroutine rns_makeplotvar
