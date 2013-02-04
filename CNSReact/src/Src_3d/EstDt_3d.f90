subroutine cns_estdt(u,u_l1,u_l2,u_l3,u_h1,u_h2,u_h3,lo,hi,dx,dt)

  use chemistry_module, only : nspecies
  use eos_module
  use meth_params_module, only : NVAR, URHO, UMX, UMY, UMZ, UEINT, UTEMP, UFS

  implicit none

  integer, intent(in) :: u_l1,u_l2,u_l3,u_h1,u_h2,u_h3
  integer, intent(in) :: lo(3), hi(3)
  double precision, intent(in) :: u(u_l1:u_h1,u_l2:u_h2,u_l3:u_h3,NVAR)
  double precision, intent(in) :: dx(3) 
  double precision, intent(inout) :: dt

  double precision :: c, T, xn(nspecies)
  double precision :: rhoInv,ux,uy,uz,dt1,dt2,dt3
  integer          :: i,j,k
  integer          :: pt_index(3)

  !$OMP PARALLEL DO PRIVATE(i,j,k,rhoInv,ux,uy,uz,T,xn,pt_index,c,dt1,dt2,dt3) REDUCTION(min:dt)
  do k = lo(3),hi(3)
     do j = lo(2),hi(2)
        do i = lo(1),hi(1)

           rhoInv = 1.d0/u(i,j,k,URHO)

           ux = u(i,j,k,UMX)*rhoInv
           uy = u(i,j,k,UMY)*rhoInv
           uz = u(i,j,k,UMZ)*rhoInv
           T  = u(i,j,k,UTEMP)

           xn(1:nspecies)=u(i,j,k,UFS:UFS+nspecies-1)*rhoInv

           pt_index(1) = i
           pt_index(2) = j
           pt_index(3) = k
           call eos_get_c(c, u(i,j,k,URHO), T, xn, pt_index)

           dt1 = dx(1)/(c + abs(ux))
           dt2 = dx(2)/(c + abs(uy))
           dt3 = dx(3)/(c + abs(uz))
           dt = min(dt,dt1,dt2,dt3)

        enddo
     enddo
  enddo
  !$OMP END PARALLEL DO

end subroutine cns_estdt
