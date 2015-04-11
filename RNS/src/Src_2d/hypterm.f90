module hypterm_module

  implicit none

  private

  public :: hypterm

contains

  subroutine hypterm(lo,hi,U,Ulo,Uhi,fx,fxlo,fxhi,fy,fylo,fyhi,dx)

    use meth_params_module, only : NVAR, difmag
    use hypterm_xy_module, only : hypterm_xy

    integer, intent(in) :: lo(2), hi(2), Ulo(2), Uhi(2), fxlo(2), fxhi(2), fylo(2), fyhi(2)
    double precision, intent(in   ) :: dx(2)
    double precision, intent(in   ) ::  U( Ulo(1): Uhi(1), Ulo(2): Uhi(2),NVAR)
    double precision, intent(inout) :: fx(fxlo(1):fxhi(1),fxlo(2):fxhi(2),NVAR)
    double precision, intent(inout) :: fy(fylo(1):fyhi(1),fylo(2):fyhi(2),NVAR)

    integer :: tlo(3), thi(3), tUlo(3), tUhi(3), tfxlo(3), tfxhi(3), tfylo(3), tfyhi(3)
    double precision :: tdx(3)

    tlo(1:2) = lo
    tlo(3) = 1
    thi(1:2) = hi
    thi(3) = 1
    
    tUlo(1:2) = Ulo
    tUlo(3) = 1
    tUhi(1:2) = Uhi
    tUhi(3) = 1
    
    tfxlo(1:2) = fxlo
    tfxlo(3) = 1
    tfxhi(1:2) = fxhi
    tfxhi(3) = 1
    
    tfylo(1:2) = fylo
    tfylo(3) = 1
    tfyhi(1:2) = fyhi
    tfyhi(3) = 1
    
    tdx(1:2) = dx
    tdx(3) = 0.d0
    
    call hypterm_xy(1.d0,tlo,thi,U,tUlo,tUhi,fx,tfxlo,tfxhi,fy,tfylo,tfyhi,tdx,tlo,thi)

!    fx = 0.d0
!    call hypterm_x(lo,hi,U,Ulo,Uhi,fx,fxlo,fxhi,dx)

    if (difmag .gt. 0.0d0) then
       call add_artifical_viscocity(lo,hi,U,Ulo,Uhi,fx,fxlo,fxhi,fy,fylo,fyhi,dx)
    end if

  end subroutine hypterm


  subroutine hypterm_x(lo,hi,U,Ulo,Uhi,fx,fxlo,fxhi,dx)

    use meth_params_module, only : NVAR, URHO, UMX, UMY, UTEMP, UFS, UEDEN, NSPEC, NCHARV, CFS
    use reconstruct_module, only : get_eigen_matrices_q
    use renorm_module, only : floor_species
    use eos_module, only : eos_get_eref
    use mdcd_module, only : mdcd
    use weno_module, only : weno5
    use riemann_module, only : riemann

    integer, intent(in) :: lo(2), hi(2), Ulo(2), Uhi(2), fxlo(2), fxhi(2)
    double precision, intent(in   ) :: dx(2)
    double precision, intent(in   ) ::  U( Ulo(1): Uhi(1), Ulo(2): Uhi(2),NVAR)
    double precision, intent(inout) :: fx(fxlo(1):fxhi(1),fxlo(2):fxhi(2),NVAR)

    integer :: i, j, n, ii, jj, ivar, m, ilr, ierr
    double precision, allocatable :: Y(:,:), RoeW(:), v(:,:), flux(:,:)
    double precision, dimension(:,:,:), allocatable, target :: UL0, UR0, UL1, UR1, UL2, UR2
    double precision, dimension(:,:,:), pointer :: U0, U1, U2
    double precision :: rhoInv, rho0, Y0(nspec), T0, v0(3), fac, eref
    double precision :: p, c, gamc, dpdr(nspec), dpde, e
    double precision :: egv1(NCHARV,NCHARV), egv2(NCHARV,NCHARV), Uii(NCHARV), charv(-3:2,NCHARV)
    double precision :: cvl(NCHARV), cvr(NCHARV)

    allocate(Y   (lo(1)-3:hi(1)+3,nspec))
    allocate(RoeW(lo(1)-3:hi(1)+3))
    allocate(v   (lo(1)-3:hi(1)+3,2))
    
    allocate(flux(lo(1):hi(1)+1,NVAR))

    ! x-face, y-average
    allocate(UL0(lo(1):hi(1)+1,lo(2)-3:hi(2)+3,NVAR))
    allocate(UR0(lo(1):hi(1)+1,lo(2)-3:hi(2)+3,NVAR))

    ! x-face, y-Gauss
    allocate(UL1(lo(1):hi(1)+1,lo(2)  :hi(2)  ,NVAR))
    allocate(UR1(lo(1):hi(1)+1,lo(2)  :hi(2)  ,NVAR))
    allocate(UL2(lo(1):hi(1)+1,lo(2)  :hi(2)  ,NVAR))
    allocate(UR2(lo(1):hi(1)+1,lo(2)  :hi(2)  ,NVAR))

    UL0 = 0.d0
    UR0 = 0.d0
    UL1 = 0.d0
    UR1 = 0.d0
    UL2 = 0.d0
    UR2 = 0.d0

    do j = lo(2)-3, hi(2)+3
       do i = lo(1)-3, hi(1)+3
          rhoInv = 1.d0/U(i,j,URHO)
          RoeW(i) = sqrt(U(i,j,URHO))
          do n=1,nspec
             Y(i,n) = U(i,j,UFS+n-1)*rhoInv
          end do
          v(i,1) = U(i,j,UMX)*rhoInv
          v(i,2) = U(i,j,UMY)*rhoInv
       end do

       do i = lo(1), hi(1)+1
          rho0 = RoeW(i-1)*RoeW(i)
          fac = 1.d0/(RoeW(i-1)+RoeW(i))
          Y0 = (Y(i-1,:)*RoeW(i-1)+Y(i,:)*RoeW(i))*fac
          T0 = (U(i-1,j,UTEMP)*RoeW(i-1)+U(i,j,UTEMP)*RoeW(i))*fac
          v0(1) = (v(i-1,1)*RoeW(i-1)+v(i,1)*RoeW(i))*fac
          v0(2) = (v(i-1,2)*RoeW(i-1)+v(i,2)*RoeW(i))*fac
          v0(3) = 0.d0

          call floor_species(nspec,Y0)

          ! egv1: left matrix;  egv2: right matrix
          call get_eigen_matrices_q(rho0, Y0, T0, v0, egv1, egv2)

          do ii=-3,2
             Uii(1) = U(i+ii,j,UMX)
             Uii(2) = U(i+ii,j,UMY)
             Uii(3) = 0.d0  ! because of 2d
             Uii(4) = U(i+ii,j,UEDEN)
          
             eref = eos_get_eref(Y(i+ii,:))
             Uii(4) = Uii(4) - U(i+ii,j,URHO) * eref

             do n=1,nspec
                Uii(CFS+n-1) = U(i+ii,j,UFS+n-1)
             end do

             do n=1,NCHARV
                charv(ii,n) = dot_product(egv1(:,n),Uii)
             end do
          end do

          do ivar=1,NCHARV
             call mdcd(charv(:,ivar), cvl(ivar), cvr(ivar))
          end do

          egv1 = transpose(egv2)  ! egv1 now holds transposed right matrix
          
          do n=1,NCHARV
             UL0(i,j,UMX  ) = UL0(i,j,UMX  ) + cvl(n)*egv2(1,n)
             UR0(i,j,UMX  ) = UR0(i,j,UMX  ) + cvr(n)*egv2(1,n)
             UL0(i,j,UMY  ) = UL0(i,j,UMY  ) + cvl(n)*egv2(2,n)
             UR0(i,j,UMY  ) = UR0(i,j,UMY  ) + cvr(n)*egv2(2,n)
             UL0(i,j,UEDEN) = UL0(i,j,UEDEN) + cvl(n)*egv2(4,n)
             UR0(i,j,UEDEN) = UR0(i,j,UEDEN) + cvr(n)*egv2(4,n)
          end do

          do m=1,nspec
             UL0(i,j,UFS+m-1) = UL0(i,j,UFS+m-1) + dot_product(cvl, egv1(:,CFS+m-1))
             UR0(i,j,UFS+m-1) = UR0(i,j,UFS+m-1) + dot_product(cvr, egv1(:,CFS+m-1))
             UL0(i,j,URHO   ) = UL0(i,j,URHO) + UL0(i,j,UFS+m-1)
             UR0(i,j,URHO   ) = UR0(i,j,URHO) + UR0(i,j,UFS+m-1)
          end do
          
          rhoInv = 1.d0/UL0(i,j,URHO)
          do n=1,nspec
             Y0(n) = UL0(i,j,UFS+n-1) * rhoInv
          end do
          call floor_species(nspec, Y0)
          eref = eos_get_eref(Y0)
          UL0(i,j,UEDEN) = UL0(i,j,UEDEN) + UL0(i,j,URHO) * eref
          do n=1,nspec
             UL0(i,j,UFS+n-1) = UL0(i,j,URHO)*Y0(n)
          end do
          
          rhoInv = 1.d0/UR0(i,j,URHO)
          do n=1,nspec
             Y0(n) = UR0(i,j,UFS+n-1) * rhoInv
          end do
          call floor_species(nspec, Y0)
          eref = eos_get_eref(Y0)
          UR0(i,j,UEDEN) = UR0(i,j,UEDEN) + UR0(i,j,URHO) * eref
          do n=1,nspec
             UR0(i,j,UFS+n-1) = UR0(i,j,URHO)*Y0(n)
          end do

          call mdcd(U(i-3:i+2,j,UTEMP), UL0(i,j,UTEMP), UR0(i,j,UTEMP))
       enddo
    end do

    do ilr = 1, 2
       if (ilr .eq. 1) then
          U0 => UL0
          U1 => UL1
          U2 => UL2
       else
          U0 => UR0
          U1 => UR1
          U2 => UR2
       end if

       do j = lo(2), hi(2)
          do i = lo(1), hi(1)+1

             rho0 = U0(i,j,URHO)
             rhoInv = 1.d0/rho0
             do n=1,nspec
                Y0(n) = U0(i,j,UFS+n-1) * rhoInv
             end do
             v0(1) = U0(i,j,UMY) * rhoInv
             v0(2) = U0(i,j,UMX) * rhoInv
             v0(3) = 0.d0
             T0 = U0(i,j,UTEMP)

             ! egv1: left matrix;  egv2: right matrix
             call get_eigen_matrices_q(rho0, Y0, T0, v0, egv1, egv2)

             do jj=-2,2
                Uii(1) = U0(i,j+jj,UMY)
                Uii(2) = U0(i,j+jj,UMX)
                Uii(3) = 0.d0   ! becuase of 2d 
                Uii(4) = U0(i,j+jj,UEDEN)

                do n=1,nspec
                   Uii(CFS+n-1) = U(i,j+jj,UFS+n-1)
                end do

                rhoInv = 1.d0/U0(i,j+jj,URHO)
                Y0 = Uii(CFS:CFS+nspec-1)*rhoInv
                eref = eos_get_eref(Y0)
                Uii(4) = Uii(4) - U(i,j+jj,URHO)*eref

                do n=1,NCHARV
                   charv(jj,n) = dot_product(egv1(:,n),Uii)
                enddo
             enddo

             do ivar=1,NCHARV
                call weno5(charv(-2:2,ivar), vg1=cvl(ivar), vg2=cvr(ivar))
             enddo

             egv1 = transpose(egv2)  ! egv1 now holds transposed right matrix

             do n=1,NCHARV
                U1(i,j,UMY  ) = U1(i,j,UMY  ) + cvl(n)*egv2(1,n)
                U2(i,j,UMY  ) = U2(i,j,UMY  ) + cvr(n)*egv2(1,n)
                U1(i,j,UMX  ) = U1(i,j,UMX  ) + cvl(n)*egv2(2,n)
                U2(i,j,UMX  ) = U2(i,j,UMX  ) + cvr(n)*egv2(2,n)
                U1(i,j,UEDEN) = U1(i,j,UEDEN) + cvl(n)*egv2(4,n)
                U2(i,j,UEDEN) = U2(i,j,UEDEN) + cvr(n)*egv2(4,n)
             end do

             do m=1,nspec
                U1(i,j,UFS+m-1) = U1(i,j,UFS+m-1) + dot_product(cvl, egv1(:,CFS+m-1))
                U2(i,j,UFS+m-1) = U2(i,j,UFS+m-1) + dot_product(cvr, egv1(:,CFS+m-1))
                U1(i,j,URHO   ) = U1(i,j,URHO   ) + U1(i,j,UFS+m-1)
                U2(i,j,URHO   ) = U2(i,j,URHO   ) + U2(i,j,UFS+m-1)
             end do

             rhoInv = 1.d0/U1(i,j,URHO)
             do n=1,nspec
                Y0(n) = U1(i,j,UFS+n-1) * rhoInv
             end do
             call floor_species(nspec, Y0)
             eref = eos_get_eref(Y0)
             U1(i,j,UEDEN) = U1(i,j,UEDEN) + U1(i,j,URHO) * eref
             do n=1,nspec
                U1(i,j,UFS+n-1) = U1(i,j,URHO)*Y0(n)
             end do

             rhoInv = 1.d0/U2(i,j,URHO)
             do n=1,nspec
                Y0(n) = U2(i,j,UFS+n-1) * rhoInv
             end do
             call floor_species(nspec, Y0)
             eref = eos_get_eref(Y0)
             U2(i,j,UEDEN) = U2(i,j,UEDEN) + U2(i,j,URHO) * eref
             do n=1,nspec
                U2(i,j,UFS+n-1) = U2(i,j,URHO)*Y0(n)
             end do

             U1(i,j,UTEMP) = U0(i,j,UTEMP)
             U2(i,j,UTEMP) = U0(i,j,UTEMP)
          end do
       end do

       Nullify(U0,U1,U2)
    end do

    do j=lo(2),hi(2)
       call riemann(lo(1),hi(1),UL1(:,j,:),UR1(:,j,:),lo(1),hi(1)+1,flux,lo(1),hi(1)+1)
       do n=1,NVAR
          do i=lo(1),hi(1)+1
             fx(i,j,n) = fx(i,j,n) + 0.5d0*flux(i,n)
          end do
       end do
       call riemann(lo(1),hi(1),UL2(:,j,:),UR2(:,j,:),lo(1),hi(1)+1,flux,lo(1),hi(1)+1)
       do n=1,NVAR
          do i=lo(1),hi(1)+1
             fx(i,j,n) = fx(i,j,n) + 0.5d0*flux(i,n)
          end do
       end do
    end do

    deallocate(Y, RoeW, v, flux)
    deallocate(UL0, UR0, UL1, UR1, UL2, UR2)

  end subroutine hypterm_x


  subroutine add_artifical_viscocity(lo,hi,U,Ulo,Uhi,fx,fxlo,fxhi,fy,fylo,fyhi,dx)

    use meth_params_module, only : NVAR, URHO, UMX, UMY, UTEMP, UEDEN, UFS, NSPEC, difmag
    use eos_module, only : eos_get_T, eos_get_c

    integer, intent(in) :: lo(2), hi(2), Ulo(2), Uhi(2), fxlo(2), fxhi(2), fylo(2), fyhi(2)
    double precision, intent(in   ) :: dx(2)
    double precision, intent(in   ) ::  U( Ulo(1): Uhi(1), Ulo(2): Uhi(2),NVAR)
    double precision, intent(inout) :: fx(fxlo(1):fxhi(1),fxlo(2):fxhi(2),NVAR)
    double precision, intent(inout) :: fy(fylo(1):fyhi(1),fylo(2):fyhi(2),NVAR)
    
    double precision, allocatable :: vx(:,:), vy(:,:), cs(:,:)
    double precision :: rhoInv, hdvdx, hdvdy, hdivv, dxinv(2), e, T, Y(nspec)
    double precision :: nu, fac, cmin
    integer :: i, j, n

    dxinv = 1.d0/dx
    
    allocate(vx(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1))
    allocate(vy(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1))
    allocate(cs(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1))
    
    do    j=lo(2)-1,hi(2)+1
       do i=lo(1)-1,hi(1)+1
          rhoInv = 1.d0/U(i,j,URHO)
          vx(i,j) = U(i,j,UMX)*rhoInv
          vy(i,j) = U(i,j,UMY)*rhoInv
          e  = U(i,j,UEDEN)*rhoInv - 0.5d0*(vx(i,j)**2+vy(i,j)**2)
          Y = U(i,j,UFS:UFS+NSPEC-1)*rhoInv
          T = U(i,j,UTEMP)
          call eos_get_T(T, e, Y)
          call eos_get_c(cs(i,j), U(i,j,URHO), T, Y)
       end do
    end do

    ! x-direction

    fac = 0.25d0*dx(1)*dxinv(2)
    
    do j=lo(2),hi(2)
       do i=lo(1),hi(1)+1
          hdvdx = vx(i,j)-vx(i-1,j)
          hdvdy = fac*(-vy(i-1,j-1)-vy(i,j-1)+vy(i-1,j+1)+vy(i,j+1))
          hdivv = hdvdx + hdvdy
          if (hdivv .lt. 0.d0) then
             cmin = min(cs(i-1,j),cs(i,j))
             nu = difmag * hdivv * min(hdivv*hdivv/(cmin*cmin*difmag), 1.d0)
             do n=1,NVAR
                if (n.ne.UTEMP) then
                   fx(i,j,n) = fx(i,j,n) - nu*(U(i-1,j,n)-U(i,j,n))
                end if
             end do
          end if
       end do
    end do

    ! y-direction

    fac = 0.25d0*dxinv(1)*dx(2)
    
    do j=lo(2),hi(2)+1
       do i=lo(1),hi(1)
          hdvdx = fac*(-vx(i-1,j-1)-vx(i-1,j)+vx(i+1,j-1)+vx(i+1,j))
          hdvdy = vy(i,j)-vy(i,j-1)
          hdivv = hdvdx + hdvdy
          if (hdivv .lt. 0.d0) then
             cmin = min(cs(i,j-1),cs(i,j))
             nu = difmag * hdivv * min(hdivv*hdivv/(cmin*cmin*difmag), 1.d0)
             do n=1,NVAR
                if (n.ne.UTEMP) then
                   fy(i,j,n) = fy(i,j,n) - nu*(U(i,j-1,n)-U(i,j,n))
                end if
             end do
          end if
       end do
    end do

    deallocate(vx,vy,cs)

  end subroutine add_artifical_viscocity

end module hypterm_module

