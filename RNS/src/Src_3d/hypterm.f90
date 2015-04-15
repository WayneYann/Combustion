module hypterm_module

  implicit none

  private

  public :: hypterm

contains

  subroutine hypterm(lo,hi,U,Ulo,Uhi,fx,fxlo,fxhi,fy,fylo,fyhi,fz,fzlo,fzhi,dx)

    use meth_params_module, only : NVAR, difmag

    integer, intent(in) :: lo(3), hi(3), Ulo(3), Uhi(3), fxlo(3), fxhi(3), &
         fylo(3), fyhi(3), fzlo(3), fzhi(3)
    double precision,intent(in   )::dx(3)
    double precision,intent(in   ):: U( Ulo(1): Uhi(1), Ulo(2): Uhi(2), Ulo(3): Uhi(3),NVAR)
    double precision,intent(inout)::fx(fxlo(1):fxhi(1),fxlo(2):fxhi(2),fxlo(3):fxhi(3),NVAR)
    double precision,intent(inout)::fy(fylo(1):fyhi(1),fylo(2):fyhi(2),fylo(3):fyhi(3),NVAR)
    double precision,intent(inout)::fz(fzlo(1):fzhi(1),fzlo(2):fzhi(2),fzlo(3):fzhi(3),NVAR)

    call hypterm_x(lo,hi,U,Ulo,Uhi,fx,fxlo,fxhi)
    call hypterm_y(lo,hi,U,Ulo,Uhi,fy,fylo,fyhi)
    call hypterm_z(lo,hi,U,Ulo,Uhi,fz,fzlo,fzhi)
    
    if (difmag .gt. 0.0d0) then
       call add_artifical_viscocity(lo,hi,U,Ulo,Uhi,fx,fxlo,fxhi,fy,fylo,fyhi,fz,fzlo,fzhi,dx)
    end if

  end subroutine hypterm


  subroutine hypterm_x(lo,hi,U,Ulo,Uhi,fx,fxlo,fxhi)
    use meth_params_module, only : NVAR, URHO, UMX, UMY, UMZ, UTEMP, UFS, UEDEN, NSPEC, NCHARV, CFS, do_mdcd
    use eigen_module, only : get_eigen_matrices
    use renorm_module, only : floor_species
    use eos_module, only : eos_get_eref
    use mdcd_module, only : mdcd
    use weno_module, only : weno4_gauss, weno5_face
    use riemann_module, only : riemann

    integer, intent(in) :: lo(3), hi(3), Ulo(3), Uhi(3), fxlo(3), fxhi(3)
    double precision, intent(in   ) ::  U( Ulo(1): Uhi(1), Ulo(2): Uhi(2), Ulo(3): Uhi(3),NVAR)
    double precision, intent(inout) :: fx(fxlo(1):fxhi(1),fxlo(2):fxhi(2),fxlo(3):fxhi(3),NVAR)

    integer :: i, j, k, n, ii, jj, kk, ivar, m, gy, gz
    double precision :: RoeWl, RoeWr, rhoInvl, rhoInvr
    double precision, dimension(:,:,:,:), allocatable, target :: UZ1,UZ2
    double precision, dimension(:,:,:)  , allocatable, target :: UY1,UY2
    double precision, dimension(:,:)    , allocatable, target :: UL,UR
    double precision, dimension(:,:)    , allocatable         :: flux
    double precision, dimension(:,:,:,:), pointer             :: UZ
    double precision, dimension(:,:,:)  , pointer             :: UY
    double precision :: rhoInv, rho0, Y0(nspec), T0, v0(3), fac, eref
    double precision :: egv1(NCHARV,NCHARV), egv2(NCHARV,NCHARV), Uch(NCHARV), charv(-3:2,NCHARV)
    double precision :: cvl(NCHARV), cvr(NCHARV)

    allocate(UZ1(lo(1)-3:hi(1)+3,lo(2)-2:hi(2)+2,lo(3):hi(3),NVAR))
    allocate(UZ2(lo(1)-3:hi(1)+3,lo(2)-2:hi(2)+2,lo(3):hi(3),NVAR))

    allocate(UY1(lo(1)-3:hi(1)+3,lo(2):hi(2),NVAR))
    allocate(UY2(lo(1)-3:hi(1)+3,lo(2):hi(2),NVAR))

    allocate(UL(lo(1):hi(1)+1,NVAR))
    allocate(UR(lo(1):hi(1)+1,NVAR))
    allocate(flux(lo(1):hi(1)+1,NVAR))

    UZ1 = 0.d0
    UZ2 = 0.d0

    do k = lo(3), hi(3)
       do j = lo(2)-2, hi(2)+2
          do i = lo(1)-3, hi(1)+3
             rho0 = U(i,j,k,URHO)
             rhoInv = 1.d0/rho0
             do n=1,nspec
                Y0(n) = U(i,j,k,UFS+n-1)*rhoInv
             end do
             v0(1) = U(i,j,k,UMZ)*rhoInv
             v0(2) = U(i,j,k,UMX)*rhoInv
             v0(3) = U(i,j,k,UMY)*rhoInv
             T0 = U(i,j,k,UTEMP)

             ! egv1: left matrix;  egv2: right matrix
             call get_eigen_matrices(rho0, Y0, T0, v0, egv1, egv2)

             do kk=-2,2
                Uch(1) = U(i,j,k+kk,UMZ)
                Uch(2) = U(i,j,k+kk,UMX)
                Uch(3) = U(i,j,k+kk,UMY)
                Uch(4) = U(i,j,k+kk,UEDEN)
                rhoInv = 1.d0/U(i,j,k+kk,URHO)
                do n=1,nspec
                   Uch(CFS+n-1) = U(i,j,k+kk,UFS+n-1)
                   Y0(n) = Uch(CFS+n-1)*rhoInv
                end do
                eref = eos_get_eref(Y0)
                Uch(4) = Uch(4) - U(i,j,k+kk,URHO)*eref

                do n=1,NCHARV
                   charv(kk,n) = dot_product(egv1(:,n),Uch)
                end do
             end do

             do ivar=1,NCHARV
                call weno4_gauss(charv(-2:2,ivar), cvl(ivar), cvr(ivar))
             end do

             egv1 = transpose(egv2)  ! egv1 now holds transposed right matrix

             do n=1,NCHARV
                UZ1(i,j,k,UMZ  ) = UZ1(i,j,k,UMZ  ) + cvl(n)*egv2(1,n)
                UZ2(i,j,k,UMZ  ) = UZ2(i,j,k,UMZ  ) + cvr(n)*egv2(1,n)
                UZ1(i,j,k,UMX  ) = UZ1(i,j,k,UMX  ) + cvl(n)*egv2(2,n)
                UZ2(i,j,k,UMX  ) = UZ2(i,j,k,UMX  ) + cvr(n)*egv2(2,n)
                UZ1(i,j,k,UMY  ) = UZ1(i,j,k,UMY  ) + cvl(n)*egv2(3,n)
                UZ2(i,j,k,UMY  ) = UZ2(i,j,k,UMY  ) + cvr(n)*egv2(3,n)
                UZ1(i,j,k,UEDEN) = UZ1(i,j,k,UEDEN) + cvl(n)*egv2(4,n)
                UZ2(i,j,k,UEDEN) = UZ2(i,j,k,UEDEN) + cvr(n)*egv2(4,n)
             end do

             do m=1,nspec
                UZ1(i,j,k,UFS+m-1) = UZ1(i,j,k,UFS+m-1) + dot_product(cvl, egv1(:,CFS+m-1))
                UZ2(i,j,k,UFS+m-1) = UZ2(i,j,k,UFS+m-1) + dot_product(cvr, egv1(:,CFS+m-1))
                UZ1(i,j,k,URHO   ) = UZ1(i,j,k,URHO   ) + UZ1(i,j,k,UFS+m-1)
                UZ2(i,j,k,URHO   ) = UZ2(i,j,k,URHO   ) + UZ2(i,j,k,UFS+m-1)
             end do
             
             rhoInv = 1.d0/UZ1(i,j,k,URHO)
             do n=1,nspec
                Y0(n) = UZ1(i,j,k,UFS+n-1) * rhoInv
             end do
             call floor_species(nspec, Y0)
             eref = eos_get_eref(Y0)
             UZ1(i,j,k,UEDEN) = UZ1(i,j,k,UEDEN) + UZ1(i,j,k,URHO) * eref
             do n=1,nspec
                UZ1(i,j,k,UFS+n-1) = UZ1(i,j,k,URHO)*Y0(n)
             end do
             
             rhoInv = 1.d0/UZ2(i,j,k,URHO)
             do n=1,nspec
                Y0(n) = UZ2(i,j,k,UFS+n-1) * rhoInv
             end do
             call floor_species(nspec, Y0)
             eref = eos_get_eref(Y0)
             UZ2(i,j,k,UEDEN) = UZ2(i,j,k,UEDEN) + UZ2(i,j,k,URHO) * eref
             do n=1,nspec
                UZ2(i,j,k,UFS+n-1) = UZ2(i,j,k,URHO)*Y0(n)
             end do
             
             UZ1(i,j,k,UTEMP) = U(i,j,k,UTEMP)
             UZ2(i,j,k,UTEMP) = U(i,j,k,UTEMP)
          end do
       end do
    end do

    do gz = 1, 2
       if (gz.eq.1) then
          UZ => UZ1
       else
          UZ => UZ2
       end if

       do k=lo(3),hi(3)
          UY1 = 0.d0
          UY2 = 0.d0

          do j = lo(2), hi(2)
             do i = lo(1)-3, hi(1)+3
                rho0 = UZ(i,j,k,URHO)
                rhoInv = 1.d0/rho0
                do n=1,nspec
                   Y0(n) = UZ(i,j,k,UFS+n-1)*rhoInv
                end do
                v0(1) = UZ(i,j,k,UMY) * rhoInv
                v0(2) = UZ(i,j,k,UMZ) * rhoInv
                v0(3) = UZ(i,j,k,UMX) * rhoInv
                T0 = UZ(i,j,k,UTEMP)
                
                ! egv1: left matrix;  egv2: right matrix
                call get_eigen_matrices(rho0, Y0, T0, v0, egv1, egv2)

                do jj=-2,2
                   Uch(1) = UZ(i,j+jj,k,UMY)
                   Uch(2) = UZ(i,j+jj,k,UMZ)
                   Uch(3) = UZ(i,j+jj,k,UMX)
                   Uch(4) = UZ(i,j+jj,k,UEDEN)
                   rhoinV = 1.d0/UZ(i,j+jj,k,URHO)
                   do n=1,nspec
                      Uch(CFS+n-1) = UZ(i,j+jj,k,UFS+n-1)
                      Y0(n) = Uch(CFS+n-1)*rhoInv
                   end do
                   eref = eos_get_eref(Y0)
                   Uch(4) = Uch(4) - UZ(i,j+jj,k,URHO)*eref                   
                   do n=1,NCHARV
                      charv(jj,n) = dot_product(egv1(:,n),Uch)
                   end do
                end do

                do ivar=1,NCHARV
                   call weno4_gauss(charv(-2:2,ivar), cvl(ivar), cvr(ivar))
                end do
                
                egv1 = transpose(egv2)  ! egv1 now holds transposed right matrix
                
                do n=1,NCHARV
                   UY1(i,j,UMY  ) = UY1(i,j,UMY  ) + cvl(n)*egv2(1,n)
                   UY2(i,j,UMY  ) = UY2(i,j,UMY  ) + cvr(n)*egv2(1,n)
                   UY1(i,j,UMZ  ) = UY1(i,j,UMZ  ) + cvl(n)*egv2(2,n)
                   UY2(i,j,UMZ  ) = UY2(i,j,UMZ  ) + cvr(n)*egv2(2,n)
                   UY1(i,j,UMX  ) = UY1(i,j,UMX  ) + cvl(n)*egv2(3,n)
                   UY2(i,j,UMX  ) = UY2(i,j,UMX  ) + cvr(n)*egv2(3,n)
                   UY1(i,j,UEDEN) = UY1(i,j,UEDEN) + cvl(n)*egv2(4,n)
                   UY2(i,j,UEDEN) = UY2(i,j,UEDEN) + cvr(n)*egv2(4,n)
                end do
          
                do m=1,nspec
                   UY1(i,j,UFS+m-1) = UY1(i,j,UFS+m-1) + dot_product(cvl, egv1(:,CFS+m-1))
                   UY2(i,j,UFS+m-1) = UY2(i,j,UFS+m-1) + dot_product(cvr, egv1(:,CFS+m-1))
                   UY1(i,j,URHO   ) = UY1(i,j,URHO   ) + UY1(i,j,UFS+m-1)
                   UY2(i,j,URHO   ) = UY2(i,j,URHO   ) + UY2(i,j,UFS+m-1)
                end do

                rhoInv = 1.d0/UY1(i,j,URHO)
                do n=1,nspec
                   Y0(n) = UY1(i,j,UFS+n-1) * rhoInv
                end do
                call floor_species(nspec, Y0)
                eref = eos_get_eref(Y0)
                UY1(i,j,UEDEN) = UY1(i,j,UEDEN) + UY1(i,j,URHO) * eref
                do n=1,nspec
                   UY1(i,j,UFS+n-1) = UY1(i,j,URHO)*Y0(n)
                end do
                
                rhoInv = 1.d0/UY2(i,j,URHO)
                do n=1,nspec
                   Y0(n) = UY2(i,j,UFS+n-1) * rhoInv
                end do
                call floor_species(nspec, Y0)
                eref = eos_get_eref(Y0)
                UY2(i,j,UEDEN) = UY2(i,j,UEDEN) + UY2(i,j,URHO) * eref
                do n=1,nspec
                   UY2(i,j,UFS+n-1) = UY2(i,j,URHO)*Y0(n)
                end do
                
                UY1(i,j,UTEMP) = UZ(i,j,k,UTEMP)
                UY2(i,j,UTEMP) = UZ(i,j,k,UTEMP)
             end do
          end do

          do gy = 1, 2
             if (gy.eq.1) then
                UY => UY1
             else
                UY => UY2
             end if

             do j = lo(2), hi(2)
                UL = 0.d0
                UR = 0.d0

                do i = lo(1),hi(1)+1
                   rhoInvl = 1.d0/UY(i-1,j,URHO)
                   rhoInvr = 1.d0/UY(i  ,j,URHO)
                   RoeWl = sqrt(UY(i-1,j,URHO))
                   RoeWr = sqrt(UY(i  ,j,URHO))
                   rho0 = RoeWl*RoeWr
                   fac = 1.d0/(RoeWl+RoeWr)
                   do n=1,nspec
                      Y0(n) = (UY(i-1,j,UFS+n-1)*rhoInvl*RoeWl+UY(i,j,UFS+n-1)*rhoInvr*RoeWr)*fac
                   end do
                   T0 = (UY(i-1,j,UTEMP)*RoeWl+UY(i,j,UTEMP)*RoeWr)*fac
                   v0(1) = (UY(i-1,j,UMX)*RoeWl+UY(i,j,UMX)*RoeWr)*fac
                   v0(2) = (UY(i-1,j,UMY)*RoeWl+UY(i,j,UMY)*RoeWr)*fac
                   v0(3) = (UY(i-1,j,UMZ)*RoeWl+UY(i,j,UMZ)*RoeWr)*fac

                   call floor_species(nspec,Y0)

                   ! egv1: left matrix;  egv2: right matrix
                   call get_eigen_matrices(rho0, Y0, T0, v0, egv1, egv2)

                   do ii=-3,2
                      Uch(1) = UY(i+ii,j,UMX)
                      Uch(2) = UY(i+ii,j,UMY)
                      Uch(3) = UY(i+ii,j,UMZ)
                      Uch(4) = UY(i+ii,j,UEDEN)
                      rhoInv = 1.d0/UY(i+ii,j,URHO)
                      do n=1,nspec
                         Uch(CFS+n-1) = UY(i+ii,j,UFS+n-1)
                         Y0(n) = Uch(CFS+n-1)*rhoInv
                      end do
                      eref = eos_get_eref(Y0)
                      Uch(4) = Uch(4) - UY(i+ii,j,URHO) * eref
                      do n=1,NCHARV
                         charv(ii,n) = dot_product(egv1(:,n),Uch)
                      end do
                   end do
                   
                   if (do_mdcd) then
                      do ivar=1,NCHARV
                         call mdcd(charv(:,ivar), cvl(ivar), cvr(ivar))
                      end do
                   else
                      do ivar=1,NCHARV
                         call weno5_face(charv(:,ivar), cvl(ivar), cvr(ivar))
                      end do
                   end if
                   
                   egv1 = transpose(egv2)  ! egv1 now holds transposed right matrix

                   do n=1,NCHARV
                      UL(i,UMX  ) = UL(i,UMX  ) + cvl(n)*egv2(1,n)
                      UR(i,UMX  ) = UR(i,UMX  ) + cvr(n)*egv2(1,n)
                      UL(i,UMY  ) = UL(i,UMY  ) + cvl(n)*egv2(2,n)
                      UR(i,UMY  ) = UR(i,UMY  ) + cvr(n)*egv2(2,n)
                      UL(i,UMZ  ) = UL(i,UMZ  ) + cvl(n)*egv2(3,n)
                      UR(i,UMZ  ) = UR(i,UMZ  ) + cvr(n)*egv2(3,n)
                      UL(i,UEDEN) = UL(i,UEDEN) + cvl(n)*egv2(4,n)
                      UR(i,UEDEN) = UR(i,UEDEN) + cvr(n)*egv2(4,n)
                   end do
             
                   do m=1,nspec
                      UL(i,UFS+m-1) = UL(i,UFS+m-1) + dot_product(cvl, egv1(:,CFS+m-1))
                      UR(i,UFS+m-1) = UR(i,UFS+m-1) + dot_product(cvr, egv1(:,CFS+m-1))
                      UL(i,URHO   ) = UL(i,URHO) + UL(i,UFS+m-1)
                      UR(i,URHO   ) = UR(i,URHO) + UR(i,UFS+m-1)
                   end do
             
                   rhoInv = 1.d0/UL(i,URHO)
                   do n=1,nspec
                      Y0(n) = UL(i,UFS+n-1) * rhoInv
                   end do
                   call floor_species(nspec, Y0)
                   eref = eos_get_eref(Y0)
                   UL(i,UEDEN) = UL(i,UEDEN) + UL(i,URHO) * eref
                   do n=1,nspec
                      UL(i,UFS+n-1) = UL(i,URHO)*Y0(n)
                   end do
             
                   rhoInv = 1.d0/UR(i,URHO)
                   do n=1,nspec
                      Y0(n) = UR(i,UFS+n-1) * rhoInv
                   end do
                   call floor_species(nspec, Y0)
                   eref = eos_get_eref(Y0)
                   UR(i,UEDEN) = UR(i,UEDEN) + UR(i,URHO) * eref
                   do n=1,nspec
                      UR(i,UFS+n-1) = UR(i,URHO)*Y0(n)
                   end do
                   
                   UL(i,UTEMP) = T0
                   UR(i,UTEMP) = T0
                end do
             
                call riemann(lo(1),hi(1),UL,UR,lo(1),hi(1)+1,flux,lo(1),hi(1)+1,dir=1)
                do n=1,NVAR
                   do i=lo(1),hi(1)+1
                      fx(i,j,k,n) = fx(i,j,k,n) + 0.25d0*flux(i,n)
                   end do
                end do
             end do
             Nullify(UY)
          end do
       end do
       Nullify(UZ)
    end do

    deallocate(UZ1,UZ2,UY1,UY2,UL,UR,flux)

  end subroutine hypterm_x


  subroutine hypterm_y(lo,hi,U,Ulo,Uhi,fy,fylo,fyhi)
    use meth_params_module, only : NVAR, URHO, UMX, UMY, UMZ, UTEMP, UFS, UEDEN, NSPEC, NCHARV, CFS, do_mdcd
    use eigen_module, only : get_eigen_matrices
    use renorm_module, only : floor_species
    use eos_module, only : eos_get_eref
    use mdcd_module, only : mdcd
    use weno_module, only : weno4_gauss, weno5_face
    use riemann_module, only : riemann

    integer, intent(in) :: lo(3), hi(3), Ulo(3), Uhi(3), fylo(3), fyhi(3)
    double precision, intent(in   ) ::  U( Ulo(1): Uhi(1), Ulo(2): Uhi(2), Ulo(3): Uhi(3),NVAR)
    double precision, intent(inout) :: fy(fylo(1):fyhi(1),fylo(2):fyhi(2),fylo(3):fyhi(3),NVAR)

    integer :: i, j, k, n, ii, jj, kk, ivar, m, gx, gz
    double precision :: RoeWl, RoeWr, rhoInvl, rhoInvr
    double precision, dimension(:,:,:,:), allocatable, target :: UX1,UX2
    double precision, dimension(:,:,:)  , allocatable, target :: UZ1,UZ2
    double precision, dimension(:,:)    , allocatable, target :: UL,UR
    double precision, dimension(:,:)    , allocatable         :: flux
    double precision, dimension(:,:,:,:), pointer             :: UX
    double precision, dimension(:,:,:)  , pointer             :: UZ
    double precision :: rhoInv, rho0, Y0(nspec), T0, v0(3), fac, eref
    double precision :: egv1(NCHARV,NCHARV), egv2(NCHARV,NCHARV), Uch(NCHARV), charv(-3:2,NCHARV)
    double precision :: cvl(NCHARV), cvr(NCHARV)

    allocate(UX1(lo(1):hi(1),lo(2)-3:hi(2)+3,lo(3)-2:hi(3)+2,NVAR))
    allocate(UX2(lo(1):hi(1),lo(2)-3:hi(2)+3,lo(3)-2:hi(3)+2,NVAR))

    allocate(UZ1(lo(2)-3:hi(2)+3,lo(3):hi(3),NVAR))
    allocate(UZ2(lo(2)-3:hi(2)+3,lo(3):hi(3),NVAR))

    allocate(UL(lo(2):hi(2)+1,NVAR))
    allocate(UR(lo(2):hi(2)+1,NVAR))
    allocate(flux(lo(2):hi(2)+1,NVAR))

    UX1 = 0.d0
    UX2 = 0.d0

    do k = lo(3)-2, hi(3)-2
       do j = lo(2)-3, hi(2)+3
          do i = lo(1), hi(1)
             rho0 = U(i,j,k,URHO)
             rhoInv = 1.d0/rho0
             do n=1,nspec
                Y0(n) = U(i,j,k,UFS+n-1)*rhoInv
             end do
             v0(1) = U(i,j,k,UMX)*rhoInv
             v0(2) = U(i,j,k,UMY)*rhoInv
             v0(3) = U(i,j,k,UMZ)*rhoInv
             T0 = U(i,j,k,UTEMP)

             ! egv1: left matrix;  egv2: right matrix
             call get_eigen_matrices(rho0, Y0, T0, v0, egv1, egv2)

             do ii=-2,2
                Uch(1) = U(i+ii,j,k,UMX)
                Uch(2) = U(i+ii,j,k,UMY)
                Uch(3) = U(i+ii,j,k,UMZ)
                Uch(4) = U(i+ii,j,k,UEDEN)
                rhoInv = 1.d0/U(i+ii,j,k,URHO)
                do n=1,nspec
                   Uch(CFS+n-1) = U(i+ii,j,k,UFS+n-1)
                   Y0(n) = Uch(CFS+n-1)*rhoInv
                end do
                eref = eos_get_eref(Y0)
                Uch(4) = Uch(4) - U(i+ii,j,k,URHO)*eref
                do n=1,NCHARV
                   charv(ii,n) = dot_product(egv1(:,n),Uch)
                end do
             end do

             do ivar=1,NCHARV
                call weno4_gauss(charv(-2:2,ivar), cvl(ivar), cvr(ivar))
             end do

             egv1 = transpose(egv2)  ! egv1 now holds transposed right matrix

             do n=1,NCHARV
                UX1(i,j,k,UMX  ) = UX1(i,j,k,UMX  ) + cvl(n)*egv2(1,n)
                UX2(i,j,k,UMX  ) = UX2(i,j,k,UMX  ) + cvr(n)*egv2(1,n)
                UX1(i,j,k,UMY  ) = UX1(i,j,k,UMY  ) + cvl(n)*egv2(2,n)
                UX2(i,j,k,UMY  ) = UX2(i,j,k,UMY  ) + cvr(n)*egv2(2,n)
                UX1(i,j,k,UMZ  ) = UX1(i,j,k,UMZ  ) + cvl(n)*egv2(3,n)
                UX2(i,j,k,UMZ  ) = UX2(i,j,k,UMZ  ) + cvr(n)*egv2(3,n)
                UX1(i,j,k,UEDEN) = UX1(i,j,k,UEDEN) + cvl(n)*egv2(4,n)
                UX2(i,j,k,UEDEN) = UX2(i,j,k,UEDEN) + cvr(n)*egv2(4,n)
             end do

             do m=1,nspec
                UX1(i,j,k,UFS+m-1) = UX1(i,j,k,UFS+m-1) + dot_product(cvl, egv1(:,CFS+m-1))
                UX2(i,j,k,UFS+m-1) = UX2(i,j,k,UFS+m-1) + dot_product(cvr, egv1(:,CFS+m-1))
                UX1(i,j,k,URHO   ) = UX1(i,j,k,URHO   ) + UX1(i,j,k,UFS+m-1)
                UX2(i,j,k,URHO   ) = UX2(i,j,k,URHO   ) + UX2(i,j,k,UFS+m-1)
             end do
             
             rhoInv = 1.d0/UX1(i,j,k,URHO)
             do n=1,nspec
                Y0(n) = UX1(i,j,k,UFS+n-1) * rhoInv
             end do
             call floor_species(nspec, Y0)
             eref = eos_get_eref(Y0)
             UX1(i,j,k,UEDEN) = UX1(i,j,k,UEDEN) + UX1(i,j,k,URHO) * eref
             do n=1,nspec
                UX1(i,j,k,UFS+n-1) = UX1(i,j,k,URHO)*Y0(n)
             end do
             
             rhoInv = 1.d0/UX2(i,j,k,URHO)
             do n=1,nspec
                Y0(n) = UX2(i,j,k,UFS+n-1) * rhoInv
             end do
             call floor_species(nspec, Y0)
             eref = eos_get_eref(Y0)
             UX2(i,j,k,UEDEN) = UX2(i,j,k,UEDEN) + UX2(i,j,k,URHO) * eref
             do n=1,nspec
                UX2(i,j,k,UFS+n-1) = UX2(i,j,k,URHO)*Y0(n)
             end do
             
             UX1(i,j,k,UTEMP) = U(i,j,k,UTEMP)
             UX2(i,j,k,UTEMP) = U(i,j,k,UTEMP)
          end do
       end do
    end do

    do gx = 1, 2
       if (gx.eq.1) then
          UX => UX1
       else
          UX => UX2
       end if

       do i=lo(1),hi(1)
          UZ1 = 0.d0
          UZ2 = 0.d0

          do k=lo(3), hi(3)
             do j = lo(2)-3, hi(2)+3
                rho0 = UX(i,j,k,URHO)
                rhoInv = 1.d0/rho0
                do n=1,nspec
                   Y0(n) = UX(i,j,k,UFS+n-1)*rhoInv
                end do
                v0(1) = UX(i,j,k,UMZ) * rhoInv
                v0(2) = UX(i,j,k,UMX) * rhoInv
                v0(3) = UX(i,j,k,UMY) * rhoInv
                T0 = UX(i,j,k,UTEMP)

                ! egv1: left matrix;  egv2: right matrix
                call get_eigen_matrices(rho0, Y0, T0, v0, egv1, egv2)

                do kk=-2,2
                   Uch(1) = UX(i,j,k+kk,UMZ)
                   Uch(2) = UX(i,j,k+kk,UMX)
                   Uch(3) = UX(i,j,k+kk,UMY)
                   Uch(4) = UX(i,j,k+kk,UEDEN)
                   rhoinV = 1.d0/UX(i,j,k+kk,URHO)
                   do n=1,nspec
                      Uch(CFS+n-1) = UX(i,j,k+kk,UFS+n-1)
                      Y0(n) = Uch(CFS+n-1)*rhoInv
                   end do
                   eref = eos_get_eref(Y0)
                   Uch(4) = Uch(4) - UX(i,j,k+kk,URHO)*eref
                   do n=1,NCHARV
                      charv(kk,n) = dot_product(egv1(:,n),Uch)
                   end do
                end do

                do ivar=1,NCHARV
                   call weno4_gauss(charv(-2:2,ivar), cvl(ivar), cvr(ivar))
                end do
                
                egv1 = transpose(egv2)  ! egv1 now holds transposed right matrix

                do n=1,NCHARV
                   UZ1(j,k,UMZ  ) = UZ1(j,k,UMZ  ) + cvl(n)*egv2(1,n)
                   UZ2(j,k,UMZ  ) = UZ2(j,k,UMZ  ) + cvr(n)*egv2(1,n)
                   UZ1(j,k,UMX  ) = UZ1(j,k,UMX  ) + cvl(n)*egv2(2,n)
                   UZ2(j,k,UMX  ) = UZ2(j,k,UMX  ) + cvr(n)*egv2(2,n)
                   UZ1(j,k,UMY  ) = UZ1(j,k,UMY  ) + cvl(n)*egv2(3,n)
                   UZ2(j,k,UMY  ) = UZ2(j,k,UMY  ) + cvr(n)*egv2(3,n)
                   UZ1(j,k,UEDEN) = UZ1(j,k,UEDEN) + cvl(n)*egv2(4,n)
                   UZ2(j,k,UEDEN) = UZ2(j,k,UEDEN) + cvr(n)*egv2(4,n)
                end do
          
                do m=1,nspec
                   UZ1(j,k,UFS+m-1) = UZ1(j,k,UFS+m-1) + dot_product(cvl, egv1(:,CFS+m-1))
                   UZ2(j,k,UFS+m-1) = UZ2(j,k,UFS+m-1) + dot_product(cvr, egv1(:,CFS+m-1))
                   UZ1(j,k,URHO   ) = UZ1(j,k,URHO   ) + UZ1(j,k,UFS+m-1)
                   UZ2(j,k,URHO   ) = UZ2(j,k,URHO   ) + UZ2(j,k,UFS+m-1)
                end do

                rhoInv = 1.d0/UZ1(j,k,URHO)
                do n=1,nspec
                   Y0(n) = UZ1(j,k,UFS+n-1) * rhoInv
                end do
                call floor_species(nspec, Y0)
                eref = eos_get_eref(Y0)
                UZ1(j,k,UEDEN) = UZ1(j,k,UEDEN) + UZ1(j,k,URHO) * eref
                do n=1,nspec
                   UZ1(j,k,UFS+n-1) = UZ1(j,k,URHO)*Y0(n)
                end do
                
                rhoInv = 1.d0/UZ2(j,k,URHO)
                do n=1,nspec
                   Y0(n) = UZ2(j,k,UFS+n-1) * rhoInv
                end do
                call floor_species(nspec, Y0)
                eref = eos_get_eref(Y0)
                UZ2(j,k,UEDEN) = UZ2(j,k,UEDEN) + UZ2(j,k,URHO) * eref
                do n=1,nspec
                   UZ2(j,k,UFS+n-1) = UZ2(j,k,URHO)*Y0(n)
                end do
                
                UZ1(j,k,UTEMP) = UX(i,j,k,UTEMP)
                UZ2(j,k,UTEMP) = UX(i,j,k,UTEMP)
             end do
          end do

          do gz = 1, 2
             if (gz .eq. 1) then
                UZ => UZ1
             else
                UZ => UZ2
             end if

             do k=lo(3), hi(3)
                UL = 0.d0
                UR = 0.d0

                do j = lo(2),hi(2)+1
                   rhoInvl = 1.d0/UZ(j-1,k,URHO)
                   rhoInvr = 1.d0/UZ(j  ,k,URHO)
                   RoeWl = sqrt(UZ(j-1,k,URHO))
                   RoeWr = sqrt(UZ(j  ,k,URHO))
                   rho0 = RoeWl*RoeWr
                   fac = 1.d0/(RoeWl+RoeWr)
                   do n=1,nspec
                      Y0(n) = (UZ(j-1,k,UFS+n-1)*rhoInvl*RoeWl+UZ(j,k,UFS+n-1)*rhoInvr*RoeWr)*fac
                   end do
                   T0 = (UZ(j-1,k,UTEMP)*RoeWl+UZ(j,k,UTEMP)*RoeWr)*fac
                   v0(1) = (UZ(j-1,k,UMY)*RoeWl+UZ(j,k,UMY)*RoeWr)*fac
                   v0(2) = (UZ(j-1,k,UMZ)*RoeWl+UZ(j,k,UMZ)*RoeWr)*fac
                   v0(3) = (UZ(j-1,k,UMX)*RoeWl+UZ(j,k,UMX)*RoeWr)*fac

                   call floor_species(nspec,Y0)

                   ! egv1: left matrix;  egv2: right matrix
                   call get_eigen_matrices(rho0, Y0, T0, v0, egv1, egv2)

                   do jj=-3,2
                      Uch(1) = UZ(j+jj,k,UMY)
                      Uch(2) = UZ(j+jj,k,UMZ)
                      Uch(3) = UZ(j+jj,k,UMX)
                      Uch(4) = UZ(j+jj,k,UEDEN)
                      rhoInv = 1.d0/UZ(j+jj,k,URHO)
                      do n=1,nspec
                         Uch(CFS+n-1) = UZ(j+jj,k,UFS+n-1)
                         Y0(n) = Uch(CFS+n-1)*rhoInv
                      end do
                      eref = eos_get_eref(Y0)
                      Uch(4) = Uch(4) - UZ(j+jj,k,URHO) * eref
                      do n=1,NCHARV
                         charv(jj,n) = dot_product(egv1(:,n),Uch)
                      end do
                   end do
                   
                   if (do_mdcd) then
                      do ivar=1,NCHARV
                         call mdcd(charv(:,ivar), cvl(ivar), cvr(ivar))
                      end do
                   else
                      do ivar=1,NCHARV
                         call weno5_face(charv(:,ivar), cvl(ivar), cvr(ivar))
                      end do
                   end if
                   
                   egv1 = transpose(egv2)  ! egv1 now holds transposed right matrix

                   do n=1,NCHARV
                      UL(j,UMY  ) = UL(j,UMY  ) + cvl(n)*egv2(1,n)
                      UR(j,UMY  ) = UR(j,UMY  ) + cvr(n)*egv2(1,n)
                      UL(j,UMZ  ) = UL(j,UMZ  ) + cvl(n)*egv2(2,n)
                      UR(j,UMZ  ) = UR(j,UMZ  ) + cvr(n)*egv2(2,n)
                      UL(j,UMX  ) = UL(j,UMX  ) + cvl(n)*egv2(3,n)
                      UR(j,UMX  ) = UR(j,UMX  ) + cvr(n)*egv2(3,n)
                      UL(j,UEDEN) = UL(j,UEDEN) + cvl(n)*egv2(4,n)
                      UR(j,UEDEN) = UR(j,UEDEN) + cvr(n)*egv2(4,n)
                   end do
             
                   do m=1,nspec
                      UL(j,UFS+m-1) = UL(j,UFS+m-1) + dot_product(cvl, egv1(:,CFS+m-1))
                      UR(j,UFS+m-1) = UR(j,UFS+m-1) + dot_product(cvr, egv1(:,CFS+m-1))
                      UL(j,URHO   ) = UL(j,URHO) + UL(j,UFS+m-1)
                      UR(j,URHO   ) = UR(j,URHO) + UR(j,UFS+m-1)
                   end do
             
                   rhoInv = 1.d0/UL(j,URHO)
                   do n=1,nspec
                      Y0(n) = UL(j,UFS+n-1) * rhoInv
                   end do
                   call floor_species(nspec, Y0)
                   eref = eos_get_eref(Y0)
                   UL(j,UEDEN) = UL(j,UEDEN) + UL(j,URHO) * eref
                   do n=1,nspec
                      UL(j,UFS+n-1) = UL(j,URHO)*Y0(n)
                   end do
             
                   rhoInv = 1.d0/UR(j,URHO)
                   do n=1,nspec
                      Y0(n) = UR(j,UFS+n-1) * rhoInv
                   end do
                   call floor_species(nspec, Y0)
                   eref = eos_get_eref(Y0)
                   UR(j,UEDEN) = UR(j,UEDEN) + UR(j,URHO) * eref
                   do n=1,nspec
                      UR(j,UFS+n-1) = UR(j,URHO)*Y0(n)
                   end do
                   
                   UL(j,UTEMP) = T0
                   UR(j,UTEMP) = T0
                end do
             
                call riemann(lo(2),hi(2),UL,UR,lo(2),hi(2)+1,flux,lo(2),hi(2)+1,dir=2)
                do n=1,NVAR
                   do j=lo(2),hi(2)+1
                      fy(i,j,k,n) = fy(i,j,k,n) + 0.25d0*flux(j,n)
                   end do
                end do
             end do
             Nullify(UZ)
          end do
       end do
       Nullify(UX)
    end do

    deallocate(UX1,UX2,UZ1,UZ2,UL,UR,flux)

  end subroutine hypterm_y


  subroutine hypterm_z(lo,hi,U,Ulo,Uhi,fz,fzlo,fzhi)
    use meth_params_module, only : NVAR, URHO, UMX, UMY, UMZ, UTEMP, UFS, UEDEN, NSPEC, NCHARV, CFS, do_mdcd
    use eigen_module, only : get_eigen_matrices
    use renorm_module, only : floor_species
    use eos_module, only : eos_get_eref
    use mdcd_module, only : mdcd
    use weno_module, only : weno4_gauss, weno5_face
    use riemann_module, only : riemann

    integer, intent(in) :: lo(3), hi(3), Ulo(3), Uhi(3), fzlo(3), fzhi(3)
    double precision, intent(in   ) ::  U( Ulo(1): Uhi(1), Ulo(2): Uhi(2), Ulo(3): Uhi(3),NVAR)
    double precision, intent(inout) :: fz(fzlo(1):fzhi(1),fzlo(2):fzhi(2),fzlo(3):fzhi(3),NVAR)

    integer :: i, j, k, n, ii, jj, kk, ivar, m, gx, gy
    double precision :: RoeWl, RoeWr, rhoInvl, rhoInvr
    double precision, dimension(:,:,:,:), allocatable, target :: UY1,UY2
    double precision, dimension(:,:,:)  , allocatable, target :: UX1,UX2
    double precision, dimension(:,:)    , allocatable, target :: UL,UR
    double precision, dimension(:,:)    , allocatable         :: flux
    double precision, dimension(:,:,:,:), pointer             :: UY
    double precision, dimension(:,:,:)  , pointer             :: UX
    double precision :: rhoInv, rho0, Y0(nspec), T0, v0(3), fac, eref
    double precision :: egv1(NCHARV,NCHARV), egv2(NCHARV,NCHARV), Uch(NCHARV), charv(-3:2,NCHARV)
    double precision :: cvl(NCHARV), cvr(NCHARV)

    allocate(UY1(lo(1)-2:hi(1)+2,lo(2):hi(2),lo(3)-3:hi(3)+3,NVAR))
    allocate(UY2(lo(1)-2:hi(1)+2,lo(2):hi(2),lo(3)-3:hi(3)+3,NVAR))

    allocate(UX1(lo(1):hi(1),lo(3)-3:hi(3)+3,NVAR))
    allocate(UX2(lo(1):hi(1),lo(3)-3:hi(3)+3,NVAR))

    allocate(UL(lo(3):hi(3)+1,NVAR))
    allocate(UR(lo(3):hi(3)+1,NVAR))
    allocate(flux(lo(3):hi(3)+1,NVAR))

    UY1 = 0.d0
    UY2 = 0.d0

    do k = lo(3)-3, hi(3)+3
       do j = lo(2), hi(2)
          do i = lo(1)-2, hi(1)+2
             rho0 = U(i,j,k,URHO)
             rhoInv = 1.d0/rho0
             do n=1,nspec
                Y0(n) = U(i,j,k,UFS+n-1)*rhoInv
             end do
             v0(1) = U(i,j,k,UMY)*rhoInv
             v0(2) = U(i,j,k,UMZ)*rhoInv
             v0(3) = U(i,j,k,UMX)*rhoInv
             T0 = U(i,j,k,UTEMP)

             ! egv1: left matrix;  egv2: right matrix
             call get_eigen_matrices(rho0, Y0, T0, v0, egv1, egv2)

             do jj=-2,2
                Uch(1) = U(i,j+jj,k,UMY)
                Uch(2) = U(i,j+jj,k,UMZ)
                Uch(3) = U(i,j+jj,k,UMX)
                Uch(4) = U(i,j+jj,k,UEDEN)
                rhoInv = 1.d0/U(i,j+jj,k,URHO)
                do n=1,nspec
                   Uch(CFS+n-1) = U(i,j+jj,k,UFS+n-1)
                   Y0(n) = Uch(CFS+n-1)*rhoInv
                end do
                eref = eos_get_eref(Y0)
                Uch(4) = Uch(4) - U(i,j+jj,k,URHO)*eref
                do n=1,NCHARV
                   charv(jj,n) = dot_product(egv1(:,n),Uch)
                end do
             end do

             do ivar=1,NCHARV
                call weno4_gauss(charv(-2:2,ivar), cvl(ivar), cvr(ivar))
             end do

             egv1 = transpose(egv2)  ! egv1 now holds transposed right matrix

             do n=1,NCHARV
                UY1(i,j,k,UMY  ) = UY1(i,j,k,UMY  ) + cvl(n)*egv2(1,n)
                UY2(i,j,k,UMY  ) = UY2(i,j,k,UMY  ) + cvr(n)*egv2(1,n)
                UY1(i,j,k,UMZ  ) = UY1(i,j,k,UMZ  ) + cvl(n)*egv2(2,n)
                UY2(i,j,k,UMZ  ) = UY2(i,j,k,UMZ  ) + cvr(n)*egv2(2,n)
                UY1(i,j,k,UMX  ) = UY1(i,j,k,UMX  ) + cvl(n)*egv2(3,n)
                UY2(i,j,k,UMX  ) = UY2(i,j,k,UMX  ) + cvr(n)*egv2(3,n)
                UY1(i,j,k,UEDEN) = UY1(i,j,k,UEDEN) + cvl(n)*egv2(4,n)
                UY2(i,j,k,UEDEN) = UY2(i,j,k,UEDEN) + cvr(n)*egv2(4,n)
             end do

             do m=1,nspec
                UY1(i,j,k,UFS+m-1) = UY1(i,j,k,UFS+m-1) + dot_product(cvl, egv1(:,CFS+m-1))
                UY2(i,j,k,UFS+m-1) = UY2(i,j,k,UFS+m-1) + dot_product(cvr, egv1(:,CFS+m-1))
                UY1(i,j,k,URHO   ) = UY1(i,j,k,URHO   ) + UY1(i,j,k,UFS+m-1)
                UY2(i,j,k,URHO   ) = UY2(i,j,k,URHO   ) + UY2(i,j,k,UFS+m-1)
             end do
             
             rhoInv = 1.d0/UY1(i,j,k,URHO)
             do n=1,nspec
                Y0(n) = UY1(i,j,k,UFS+n-1) * rhoInv
             end do
             call floor_species(nspec, Y0)
             eref = eos_get_eref(Y0)
             UY1(i,j,k,UEDEN) = UY1(i,j,k,UEDEN) + UY1(i,j,k,URHO) * eref
             do n=1,nspec
                UY1(i,j,k,UFS+n-1) = UY1(i,j,k,URHO)*Y0(n)
             end do
             
             rhoInv = 1.d0/UY2(i,j,k,URHO)
             do n=1,nspec
                Y0(n) = UY2(i,j,k,UFS+n-1) * rhoInv
             end do
             call floor_species(nspec, Y0)
             eref = eos_get_eref(Y0)
             UY2(i,j,k,UEDEN) = UY2(i,j,k,UEDEN) + UY2(i,j,k,URHO) * eref
             do n=1,nspec
                UY2(i,j,k,UFS+n-1) = UY2(i,j,k,URHO)*Y0(n)
             end do
             
             UY1(i,j,k,UTEMP) = U(i,j,k,UTEMP)
             UY2(i,j,k,UTEMP) = U(i,j,k,UTEMP)
          end do
       end do
    end do

    do gy = 1, 2
       if(gy .eq. 1) then
          UY => UY1
       else
          UY => UY2
       end if

       do j=lo(2),hi(2)
          UX1 = 0.d0
          UX2 = 0.d0

          do k = lo(3)-3, hi(3)+3
             do i = lo(1), hi(1)
                rho0 = UY(i,j,k,URHO)
                rhoInv = 1.d0/rho0
                do n=1,nspec
                   Y0(n) = UY(i,j,k,UFS+n-1)*rhoInv
                end do
                v0(1) = UY(i,j,k,UMX) * rhoInv
                v0(2) = UY(i,j,k,UMY) * rhoInv
                v0(3) = UY(i,j,k,UMZ) * rhoInv
                T0 = UY(i,j,k,UTEMP)

                ! egv1: left matrix;  egv2: right matrix
                call get_eigen_matrices(rho0, Y0, T0, v0, egv1, egv2)

                do ii=-2,2
                   Uch(1) = UY(i+ii,j,k,UMX)
                   Uch(2) = UY(i+ii,j,k,UMY)
                   Uch(3) = UY(i+ii,j,k,UMZ)
                   Uch(4) = UY(i+ii,j,k,UEDEN)
                   rhoinV = 1.d0/UY(i+ii,j,k,URHO)
                   do n=1,nspec
                      Uch(CFS+n-1) = UY(i+ii,j,k,UFS+n-1)
                      Y0(n) = Uch(CFS+n-1)*rhoInv
                   end do
                   eref = eos_get_eref(Y0)
                   Uch(4) = Uch(4) - UY(i+ii,j,k,URHO)*eref
                   do n=1,NCHARV
                      charv(ii,n) = dot_product(egv1(:,n),Uch)
                   end do
                end do

                do ivar=1,NCHARV
                   call weno4_gauss(charv(-2:2,ivar), cvl(ivar), cvr(ivar))
                end do
                
                egv1 = transpose(egv2)  ! egv1 now holds transposed right matrix

                do n=1,NCHARV
                   UX1(i,k,UMX  ) = UX1(i,k,UMX  ) + cvl(n)*egv2(1,n)
                   UX2(i,k,UMX  ) = UX2(i,k,UMX  ) + cvr(n)*egv2(1,n)
                   UX1(i,k,UMY  ) = UX1(i,k,UMY  ) + cvl(n)*egv2(2,n)
                   UX2(i,k,UMY  ) = UX2(i,k,UMY  ) + cvr(n)*egv2(2,n)
                   UX1(i,k,UMZ  ) = UX1(i,k,UMZ  ) + cvl(n)*egv2(3,n)
                   UX2(i,k,UMZ  ) = UX2(i,k,UMZ  ) + cvr(n)*egv2(3,n)
                   UX1(i,k,UEDEN) = UX1(i,k,UEDEN) + cvl(n)*egv2(4,n)
                   UX2(i,k,UEDEN) = UX2(i,k,UEDEN) + cvr(n)*egv2(4,n)
                end do
          
                do m=1,nspec
                   UX1(i,k,UFS+m-1) = UX1(i,k,UFS+m-1) + dot_product(cvl, egv1(:,CFS+m-1))
                   UX2(i,k,UFS+m-1) = UX2(i,k,UFS+m-1) + dot_product(cvr, egv1(:,CFS+m-1))
                   UX1(i,k,URHO   ) = UX1(i,k,URHO   ) + UX1(i,k,UFS+m-1)
                   UX2(i,k,URHO   ) = UX2(i,k,URHO   ) + UX2(i,k,UFS+m-1)
                end do

                rhoInv = 1.d0/UX1(i,k,URHO)
                do n=1,nspec
                   Y0(n) = UX1(i,k,UFS+n-1) * rhoInv
                end do
                call floor_species(nspec, Y0)
                eref = eos_get_eref(Y0)
                UX1(i,k,UEDEN) = UX1(i,k,UEDEN) + UX1(i,k,URHO) * eref
                do n=1,nspec
                   UX1(i,k,UFS+n-1) = UX1(i,k,URHO)*Y0(n)
                end do
                
                rhoInv = 1.d0/UX2(i,k,URHO)
                do n=1,nspec
                   Y0(n) = UX2(i,k,UFS+n-1) * rhoInv
                end do
                call floor_species(nspec, Y0)
                eref = eos_get_eref(Y0)
                UX2(i,k,UEDEN) = UX2(i,k,UEDEN) + UX2(i,k,URHO) * eref
                do n=1,nspec
                   UX2(i,k,UFS+n-1) = UX2(i,k,URHO)*Y0(n)
                end do
                
                UX1(i,k,UTEMP) = UY(i,j,k,UTEMP)
                UX2(i,k,UTEMP) = UY(i,j,k,UTEMP)
             end do
          end do

          do gx = 1, 2
             if (gx .eq. 1) then
                UX => UX1
             else
                UX => UX2
             end if

             do i = lo(1), hi(1)
                UL = 0.d0
                UR = 0.d0

                do k = lo(3), hi(3)+1
                   rhoInvl = 1.d0/UX(i,k-1,URHO)
                   rhoInvr = 1.d0/UX(i,k  ,URHO)
                   RoeWl = sqrt(UX(i,k-1,URHO))
                   RoeWr = sqrt(UX(i,k  ,URHO))
                   rho0 = RoeWl*RoeWr
                   fac = 1.d0/(RoeWl+RoeWr)
                   do n=1,nspec
                      Y0(n) = (UX(i,k-1,UFS+n-1)*rhoInvl*RoeWl+UX(i,k,UFS+n-1)*rhoInvr*RoeWr)*fac
                   end do
                   T0 = (UX(i,k-1,UTEMP)*RoeWl+UX(i,k,UTEMP)*RoeWr)*fac
                   v0(1) = (UX(i,k-1,UMZ)*RoeWl+UX(i,k,UMZ)*RoeWr)*fac
                   v0(2) = (UX(i,k-1,UMX)*RoeWl+UX(i,k,UMX)*RoeWr)*fac
                   v0(3) = (UX(i,k-1,UMY)*RoeWl+UX(i,k,UMY)*RoeWr)*fac

                   call floor_species(nspec,Y0)

                   ! egv1: left matrix;  egv2: right matrix
                   call get_eigen_matrices(rho0, Y0, T0, v0, egv1, egv2)

                   do kk=-3,2
                      Uch(1) = UX(i,k+kk,UMZ)
                      Uch(2) = UX(i,k+kk,UMX)
                      Uch(3) = UX(i,k+kk,UMY)
                      Uch(4) = UX(i,k+kk,UEDEN)
                      rhoInv = 1.d0/UX(i,k+kk,URHO)
                      do n=1,nspec
                         Uch(CFS+n-1) = UX(i,k+kk,UFS+n-1)
                         Y0(n) = Uch(CFS+n-1)*rhoInv
                      end do
                      eref = eos_get_eref(Y0)
                      Uch(4) = Uch(4) - UX(i,k+kk,URHO) * eref
                      do n=1,NCHARV
                         charv(kk,n) = dot_product(egv1(:,n),Uch)
                      end do
                   end do
                   
                   if (do_mdcd) then
                      do ivar=1,NCHARV
                         call mdcd(charv(:,ivar), cvl(ivar), cvr(ivar))
                      end do
                   else
                      do ivar=1,NCHARV
                         call weno5_face(charv(:,ivar), cvl(ivar), cvr(ivar))
                      end do
                   end if
                   
                   egv1 = transpose(egv2)  ! egv1 now holds transposed right matrix

                   do n=1,NCHARV
                      UL(k,UMZ  ) = UL(k,UMZ  ) + cvl(n)*egv2(1,n)
                      UR(k,UMZ  ) = UR(k,UMZ  ) + cvr(n)*egv2(1,n)
                      UL(k,UMX  ) = UL(k,UMX  ) + cvl(n)*egv2(2,n)
                      UR(k,UMX  ) = UR(k,UMX  ) + cvr(n)*egv2(2,n)
                      UL(k,UMY  ) = UL(k,UMY  ) + cvl(n)*egv2(3,n)
                      UR(k,UMY  ) = UR(k,UMY  ) + cvr(n)*egv2(3,n)
                      UL(k,UEDEN) = UL(k,UEDEN) + cvl(n)*egv2(4,n)
                      UR(k,UEDEN) = UR(k,UEDEN) + cvr(n)*egv2(4,n)
                   end do
             
                   do m=1,nspec
                      UL(k,UFS+m-1) = UL(k,UFS+m-1) + dot_product(cvl, egv1(:,CFS+m-1))
                      UR(k,UFS+m-1) = UR(k,UFS+m-1) + dot_product(cvr, egv1(:,CFS+m-1))
                      UL(k,URHO   ) = UL(k,URHO) + UL(k,UFS+m-1)
                      UR(k,URHO   ) = UR(k,URHO) + UR(k,UFS+m-1)
                   end do
             
                   rhoInv = 1.d0/UL(k,URHO)
                   do n=1,nspec
                      Y0(n) = UL(k,UFS+n-1) * rhoInv
                   end do
                   call floor_species(nspec, Y0)
                   eref = eos_get_eref(Y0)
                   UL(k,UEDEN) = UL(k,UEDEN) + UL(k,URHO) * eref
                   do n=1,nspec
                      UL(k,UFS+n-1) = UL(k,URHO)*Y0(n)
                   end do
             
                   rhoInv = 1.d0/UR(k,URHO)
                   do n=1,nspec
                      Y0(n) = UR(k,UFS+n-1) * rhoInv
                   end do
                   call floor_species(nspec, Y0)
                   eref = eos_get_eref(Y0)
                   UR(k,UEDEN) = UR(k,UEDEN) + UR(k,URHO) * eref
                   do n=1,nspec
                      UR(k,UFS+n-1) = UR(k,URHO)*Y0(n)
                   end do
                   
                   UL(k,UTEMP) = T0
                   UR(k,UTEMP) = T0
                end do

                call riemann(lo(3),hi(3),UL,UR,lo(3),hi(3)+1,flux,lo(3),hi(3)+1,dir=3)
                do n=1,NVAR
                   do k=lo(3),hi(3)+1
                      fz(i,j,k,n) = fz(i,j,k,n) + 0.25d0*flux(k,n)
                   end do
                end do
             end do
             Nullify(UX)
          end do
       end do
       Nullify(UY)
    end do

    deallocate(UY1,UY2,UX1,UX2,UL,UR,flux)

  end subroutine hypterm_z


  subroutine add_artifical_viscocity(lo,hi,U,Ulo,Uhi,fx,fxlo,fxhi,fy,fylo,fyhi,fz,fzlo,fzhi,dx)

    use meth_params_module, only : NVAR, URHO, UMX, UMY, UMZ, UTEMP, difmag

    integer, intent(in) :: lo(3), hi(3), Ulo(3), Uhi(3), fxlo(3), fxhi(3), &
         fylo(3), fyhi(3), fzlo(3), fzhi(3)
    double precision,intent(in   )::dx(3)
    double precision,intent(in   ):: U( Ulo(1): Uhi(1), Ulo(2): Uhi(2), Ulo(3): Uhi(3),NVAR)
    double precision,intent(inout)::fx(fxlo(1):fxhi(1),fxlo(2):fxhi(2),fxlo(3):fxhi(3),NVAR)
    double precision,intent(inout)::fy(fylo(1):fyhi(1),fylo(2):fyhi(2),fylo(3):fyhi(3),NVAR)
    double precision,intent(inout)::fz(fzlo(1):fzhi(1),fzlo(2):fzhi(2),fzlo(3):fzhi(3),NVAR)

    double precision, allocatable :: divv(:,:,:), vx(:,:,:), vy(:,:,:), vz(:,:,:)
    double precision :: rhoInv, dvdx, dvdy, dvdz, div1, dxinv(3)
    integer :: i, j, k, n

    dxinv = 1.d0/dx

    allocate(vx  (lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1))
    allocate(vy  (lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1))
    allocate(vz  (lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1))
    allocate(divv(lo(1)  :hi(1)+1,lo(2)  :hi(2)+1,lo(3)  :hi(3)+1))

    do       k=lo(3)-1,hi(3)+1
       do    j=lo(2)-1,hi(2)+1
          do i=lo(1)-1,hi(1)+1
             rhoInv = 1.d0/U(i,j,k,URHO)
             vx(i,j,k) = U(i,j,k,UMX)*rhoInv
             vy(i,j,k) = U(i,j,k,UMY)*rhoInv
             vz(i,j,k) = U(i,j,k,UMZ)*rhoInv
          end do
       end do
    end do

    do k=lo(3),hi(3)+1
       do j=lo(2),hi(2)+1
          do i=lo(1),hi(1)+1

             dvdx = .25d0*dxinv(1)*( &
                    + vx(i  ,j  ,k  ) - vx(i-1,j  ,k  ) &
                    + vx(i  ,j  ,k-1) - vx(i-1,j  ,k-1) &
                    + vx(i  ,j-1,k  ) - vx(i-1,j-1,k  ) &
                    + vx(i  ,j-1,k-1) - vx(i-1,j-1,k-1) )

             dvdy = .25d0*dxinv(2)*( &
                    + vy(i  ,j  ,k  ) - vy(i  ,j-1,k  ) &
                    + vy(i  ,j  ,k-1) - vy(i  ,j-1,k-1) &
                    + vy(i-1,j  ,k  ) - vy(i-1,j-1,k  ) &
                    + vy(i-1,j  ,k-1) - vy(i-1,j-1,k-1) )

             dvdz = .25d0*dxinv(3)*( &
                    + vz(i  ,j  ,k  ) - vz(i  ,j  ,k-1) &
                    + vz(i  ,j-1,k  ) - vz(i  ,j-1,k-1) &
                    + vz(i-1,j  ,k  ) - vz(i-1,j  ,k-1) &
                    + vz(i-1,j-1,k  ) - vz(i-1,j-1,k-1) )

             divv(i,j,k) = dvdx + dvdy + dvdz

          enddo
       enddo
    enddo

    do n = 1, NVAR
       if ( n.ne.UTEMP ) then
          
          do k = lo(3),hi(3)
             do j = lo(2),hi(2)
                do i = lo(1),hi(1)+1
                   div1 = .25d0*(divv(i,j,k) + divv(i,j+1,k) + divv(i,j,k+1) + divv(i,j+1,k+1))
                   div1 = difmag*min(0.d0,div1)
                   fx(i,j,k,n) = fx(i,j,k,n) + dx(1)*div1*(U(i,j,k,n)-U(i-1,j,k,n))
                enddo
             enddo
          enddo

          do k = lo(3),hi(3)
             do j = lo(2),hi(2)+1
                do i = lo(1),hi(1)
                   div1 = .25d0*(divv(i,j,k) + divv(i+1,j,k) + divv(i,j,k+1) + divv(i+1,j,k+1))
                   div1 = difmag*min(0.d0,div1)
                   fy(i,j,k,n) = fy(i,j,k,n) + dx(2)*div1*(U(i,j,k,n)-U(i,j-1,k,n))
                enddo
             enddo
          enddo

          do k = lo(3),hi(3)+1
             do j = lo(2),hi(2)
                do i = lo(1),hi(1)
                   div1 = .25d0*(divv(i,j,k) + divv(i+1,j,k) + divv(i,j+1,k) + divv(i+1,j+1,k))
                   div1 = difmag*min(0.d0,div1)
                   fz(i,j,k,n) = fz(i,j,k,n) + dx(3)*div1*(U(i,j,k,n)-U(i,j,k-1,n))
                enddo
             enddo
          enddo
          
       endif
    enddo

  end subroutine add_artifical_viscocity

end module hypterm_module
