subroutine rns_estdt_diff(u,u_l1,u_l2,u_l3,u_h1,u_h2,u_h3,lo,hi,dx,dt)
  use meth_params_module
  use egz_module, only : egzini, egzpar, egzvr1
  use chemistry_module, only : inv_mwt
  integer, intent(in) :: u_l1,u_l2,u_l3,u_h1,u_h2,u_h3
  integer, intent(in) :: lo(3), hi(3)
  double precision, intent(in) :: u(u_l1:u_h1,u_l2:u_h2,u_l3:u_h3,NVAR)
  double precision, intent(in) :: dx(3)
  double precision, intent(inout) :: dt

  integer :: i, j, k, n, np, iwrk
  double precision :: rwrk, Yt(NSPEC), Xt(NSPEC), rhoInv, Wbar, D, dx2
  double precision, allocatable :: DZ(:,:), XZ(:,:), CPZ(:,:)
  
  np = hi(1)-lo(1)+1
  
  dx2 = 0.5*(min(dx(1),dx(2),dx(3)))**2
  
  !$omp parallel private(i,j,k,n,iwrk,rwrk,Yt,Xt,rhoInv,Wbar,D, DZ,XZ,CPZ) &
  !$omp reduction(min:dt)
  
  allocate(DZ (lo(1):hi(1),nspec))
  allocate(XZ (lo(1):hi(1),nspec))
  allocate(CPZ(lo(1):hi(1),nspec))
  
  CPZ = 0.d0
  
  call egzini(np)
  
  !$omp do collapse(2)
  do k=lo(3),hi(3)
  do j=lo(2),hi(2)
     
     do i=lo(1),hi(1)
        rhoInv = 1.d0/u(i,j,k,URHO)
        Yt = u(i,j,k,UFS:UFS+nspec-1)*rhoInv
        call ckytx(Yt, iwrk, rwrk, Xt)
        XZ(i,:) = Xt
     end do

     call egzpar(u(lo(1):hi(1),j,k,UTEMP), XZ, CPZ)
     
     call egzvr1(u(lo(1):hi(1),j,k,UTEMP), DZ)
     
     do i=lo(1),hi(1)
        Xt = XZ(i,:)
        call ckmmwx(Xt, iwrk, rwrk, Wbar)
        D = 0.d0
        do n=1,nspec
           D = max(D, DZ(i,n)*Wbar*inv_mwt(n))
        end do
        dt = min(dt, dx2/D*u(i,j,k,URHO))
     end do
     
  end do
  end do
  !$omp end do
  
  deallocate(DZ,XZ,CPZ)
  
  !$omp end parallel
  
end subroutine rns_estdt_diff

