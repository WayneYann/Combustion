module kernels_module
  use chemistry_module, only : nspecies, molecular_weight
  use derivative_stencil_module, only : stencil_ng, first_deriv_8, M8 
  use variables_module
  implicit none

  private

  public :: hypterm_3d, compact_diffterm_3d, chemterm_3d, comp_courno_3d, &
       S3D_diffterm_1, S3D_diffterm_2

contains

  subroutine hypterm_3d (lo,hi,ng,dx,cons,q,rhs)

    integer,          intent(in ) :: lo(3),hi(3),ng
    double precision, intent(in ) :: dx(3)
    double precision, intent(in ) :: cons(-ng+lo(1):hi(1)+ng,-ng+lo(2):hi(2)+ng,-ng+lo(3):hi(3)+ng,ncons)
    double precision, intent(in ) ::    q(-ng+lo(1):hi(1)+ng,-ng+lo(2):hi(2)+ng,-ng+lo(3):hi(3)+ng,nprim)
    double precision, intent(out) ::  rhs(    lo(1):hi(1)   ,    lo(2):hi(2)   ,    lo(3):hi(3)   ,ncons)

    integer          :: i,j,k,n
    double precision :: un(-4:4)
    double precision :: dxinv(3)

    do i=1,3
       dxinv(i) = 1.0d0 / dx(i)
    end do

    rhs = 0.d0
    
    !$omp parallel private(i,j,k,n,un)
    !$omp do 
    do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)

             un = q(i-4:i+4,j,k,qu)

             rhs(i,j,k,irho) = rhs(i,j,k,irho) - dxinv(1) * &
                  first_deriv_8( cons(i-4:i+4,j,k,imx) ) 

             rhs(i,j,k,imx) = rhs(i,j,k,imx) - dxinv(1) * &
                  first_deriv_8( cons(i-4:i+4,j,k,imx)*un+q(i-4:i+4,j,k,qpres) )

             rhs(i,j,k,imy) = rhs(i,j,k,imy) - dxinv(1) * &
                  first_deriv_8( cons(i-4:i+4,j,k,imy)*un ) 

             rhs(i,j,k,imz) = rhs(i,j,k,imz) - dxinv(1) * &
                  first_deriv_8( cons(i-4:i+4,j,k,imz)*un ) 

             rhs(i,j,k,iene) = rhs(i,j,k,iene) - dxinv(1) * &
                  first_deriv_8( (cons(i-4:i+4,j,k,iene)+q(i-4:i+4,j,k,qpres))*un )

             do n = iry1, iry1+nspecies-1
                rhs(i,j,k,n) = rhs(i,j,k,n) - dxinv(1) * &
                     first_deriv_8( cons(i-4:i+4,j,k,n)*un )
             end do

          enddo
       enddo
    enddo
    !$omp end do

    !$omp do
    do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)

             un = q(i,j-4:j+4,k,qv)

             rhs(i,j,k,irho)=rhs(i,j,k,irho) - dxinv(2) * &
                  first_deriv_8( cons(i,j-4:j+4,k,imy) )

             rhs(i,j,k,imx)=rhs(i,j,k,imx) - dxinv(2) * &
                  first_deriv_8( cons(i,j-4:j+4,k,imx)*un )

             rhs(i,j,k,imy)=rhs(i,j,k,imy) - dxinv(2) * &
                  first_deriv_8( cons(i,j-4:j+4,k,imy)*un+q(i,j-4:j+4,k,qpres) )

             rhs(i,j,k,imz)=rhs(i,j,k,imz) - dxinv(2) * &
                  first_deriv_8( cons(i,j-4:j+4,k,imz)*un )

             rhs(i,j,k,iene)=rhs(i,j,k,iene) - dxinv(2) * &
                  first_deriv_8( (cons(i,j-4:j+4,k,iene)+q(i,j-4:j+4,k,qpres))*un )

             do n = iry1, iry1+nspecies-1
                rhs(i,j,k,n) = rhs(i,j,k,n) - dxinv(2) * &
                     first_deriv_8( cons(i,j-4:j+4,k,n)*un )
             end do

          enddo
       enddo
    enddo
    !$omp end do

    !$omp do
    do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)

             un = q(i,j,k-4:k+4,qw)

             rhs(i,j,k,irho)=rhs(i,j,k,irho) - dxinv(3) * &
                  first_deriv_8( cons(i,j,k-4:k+4,imz) )

             rhs(i,j,k,imx)=rhs(i,j,k,imx) - dxinv(3) * &
                  first_deriv_8( cons(i,j,k-4:k+4,imx)*un )

             rhs(i,j,k,imy)=rhs(i,j,k,imy) - dxinv(3) * &
                  first_deriv_8( cons(i,j,k-4:k+4,imy)*un )

             rhs(i,j,k,imz)=rhs(i,j,k,imz) - dxinv(3) * &
                  first_deriv_8( cons(i,j,k-4:k+4,imz)*un+q(i,j,k-4:k+4,qpres) )

             rhs(i,j,k,iene)=rhs(i,j,k,iene) - dxinv(3) * &
                  first_deriv_8( (cons(i,j,k-4:k+4,iene)+q(i,j,k-4:k+4,qpres))*un )

             do n = iry1, iry1+nspecies-1
                rhs(i,j,k,n) = rhs(i,j,k,n) - dxinv(3) * &
                     first_deriv_8(cons(i,j,k-4:k+4,n)*un)
             end do

          enddo
       enddo
    enddo
    !$omp end do
    !$omp end parallel

  end subroutine hypterm_3d


  subroutine compact_diffterm_3d (lo,hi,ng,dx,q,rhs,mu,xi,lam,dxy)

    integer,          intent(in ) :: lo(3),hi(3),ng
    double precision, intent(in ) :: dx(3)
    double precision, intent(in ) :: q  (-ng+lo(1):hi(1)+ng,-ng+lo(2):hi(2)+ng,-ng+lo(3):hi(3)+ng,nprim)
    double precision, intent(in ) :: mu (-ng+lo(1):hi(1)+ng,-ng+lo(2):hi(2)+ng,-ng+lo(3):hi(3)+ng)
    double precision, intent(in ) :: xi (-ng+lo(1):hi(1)+ng,-ng+lo(2):hi(2)+ng,-ng+lo(3):hi(3)+ng)
    double precision, intent(in ) :: lam(-ng+lo(1):hi(1)+ng,-ng+lo(2):hi(2)+ng,-ng+lo(3):hi(3)+ng)
    double precision, intent(in ) :: dxy(-ng+lo(1):hi(1)+ng,-ng+lo(2):hi(2)+ng,-ng+lo(3):hi(3)+ng,nspecies)
    double precision, intent(out) :: rhs(    lo(1):hi(1)   ,    lo(2):hi(2)   ,    lo(3):hi(3)   ,ncons)

    double precision, allocatable, dimension(:,:,:) :: ux,uy,uz,vx,vy,vz,wx,wy,wz
    double precision, allocatable, dimension(:,:,:) :: vsp,vsm, dpe
    double precision, allocatable, dimension(:,:,:,:) :: Hg, dpy, dxe
    ! dxy: diffusion coefficient of X in equation for Y
    ! dpy: diffusion coefficient of p in equation for Y
    ! dxe: diffusion coefficient of X in equation for energy
    ! dpe: diffusion coefficient of p in equation for energy

    double precision :: dxinv(3), dx2inv(3), divu
    double precision :: dmvxdy,dmwxdz,dmvywzdx
    double precision :: dmuydx,dmwydz,dmuxwzdy
    double precision :: dmuzdx,dmvzdy,dmuxvydz
    double precision :: tauxx,tauyy,tauzz 
    double precision :: Htot, Htmp(nspecies), Ytmp(nspecies), hhalf
    integer          :: i,j,k,n, qxn, qyn, qhn
    integer :: dlo(3), dhi(3)

    double precision :: muM8(8), M8p(8), M8X(8), tmp(8)

    do i = 1,3
       dxinv(i) = 1.0d0 / dx(i)
       dx2inv(i) = dxinv(i)**2
    end do

    dlo = lo - ng
    dhi = hi + ng

    allocate(ux( lo(1): hi(1),dlo(2):dhi(2),dlo(3):dhi(3)))
    allocate(vx( lo(1): hi(1),dlo(2):dhi(2),dlo(3):dhi(3)))
    allocate(wx( lo(1): hi(1),dlo(2):dhi(2),dlo(3):dhi(3)))

    allocate(uy(dlo(1):dhi(1), lo(2): hi(2),dlo(3):dhi(3)))
    allocate(vy(dlo(1):dhi(1), lo(2): hi(2),dlo(3):dhi(3)))
    allocate(wy(dlo(1):dhi(1), lo(2): hi(2),dlo(3):dhi(3)))

    allocate(uz(dlo(1):dhi(1),dlo(2):dhi(2), lo(3): hi(3)))
    allocate(vz(dlo(1):dhi(1),dlo(2):dhi(2), lo(3): hi(3)))
    allocate(wz(dlo(1):dhi(1),dlo(2):dhi(2), lo(3): hi(3)))

    allocate(vsp(dlo(1):dhi(1),dlo(2):dhi(2),dlo(3):dhi(3)))
    allocate(vsm(dlo(1):dhi(1),dlo(2):dhi(2),dlo(3):dhi(3)))

    !$omp parallel &
    !$omp private(i,j,k,tauxx,tauyy,tauzz,divu) &
    !$omp private(dmvxdy,dmwxdz,dmvywzdx,dmuydx,dmwydz,dmuxwzdy,dmuzdx,dmvzdy,dmuxvydz)

    !$omp workshare
    rhs = 0.d0
    !$omp end workshare

    !$OMP DO
    do k=dlo(3),dhi(3)
       do j=dlo(2),dhi(2)
          do i=dlo(1),dhi(1)
             vsp(i,j,k) = xi(i,j,k) + FourThirds*mu(i,j,k)
             vsm(i,j,k) = xi(i,j,k) -  TwoThirds*mu(i,j,k)
          enddo
       enddo
    enddo
    !$OMP END DO NOWAIT

    !$omp do
    do k=dlo(3),dhi(3)
       do j=dlo(2),dhi(2)
          do i=lo(1),hi(1)
             ux(i,j,k) = dxinv(1)*first_deriv_8(q(i-4:i+4,j,k,qu))
             vx(i,j,k) = dxinv(1)*first_deriv_8(q(i-4:i+4,j,k,qv))
             wx(i,j,k) = dxinv(1)*first_deriv_8(q(i-4:i+4,j,k,qw))
          enddo
       enddo
    enddo
    !$OMP END DO NOWAIT

    !$OMP DO
    do k=dlo(3),dhi(3)
       do j=lo(2),hi(2)   
          do i=dlo(1),dhi(1)
             uy(i,j,k) = dxinv(2)*first_deriv_8(q(i,j-4:j+4,k,qu))
             vy(i,j,k) = dxinv(2)*first_deriv_8(q(i,j-4:j+4,k,qv))
             wy(i,j,k) = dxinv(2)*first_deriv_8(q(i,j-4:j+4,k,qw))
          enddo
       enddo
    enddo
    !$OMP END DO NOWAIT

    !$OMP DO
    do k=lo(3),hi(3)
       do j=dlo(2),dhi(2)
          do i=dlo(1),dhi(1)
             uz(i,j,k) = dxinv(3)*first_deriv_8(q(i,j,k-4:k+4,qu))
             vz(i,j,k) = dxinv(3)*first_deriv_8(q(i,j,k-4:k+4,qv))
             wz(i,j,k) = dxinv(3)*first_deriv_8(q(i,j,k-4:k+4,qw))
          enddo
       enddo
    enddo
    !$OMP END DO NOWAIT

    !$omp barrier

    !$omp do
    do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)

             divu = (ux(i,j,k)+vy(i,j,k)+wz(i,j,k))*vsm(i,j,k)
             tauxx = 2.d0*mu(i,j,k)*ux(i,j,k) + divu
             tauyy = 2.d0*mu(i,j,k)*vy(i,j,k) + divu
             tauzz = 2.d0*mu(i,j,k)*wz(i,j,k) + divu
             
             ! change in internal energy
             rhs(i,j,k,iene) = rhs(i,j,k,iene) + &
                  tauxx*ux(i,j,k) + tauyy*vy(i,j,k) + tauzz*wz(i,j,k) &
                  + mu(i,j,k)*((uy(i,j,k)+vx(i,j,k))**2 &
                  &          + (wx(i,j,k)+uz(i,j,k))**2 &
                  &          + (vz(i,j,k)+wy(i,j,k))**2 )

          end do
       end do
    end do
    !$omp end do nowait

    ! d()/dx
    !$omp do
    do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             ! d((xi-2/3*mu)*(vy+wz))/dx
             dmvywzdx = dxinv(1) * &
                  first_deriv_8( vsm(i-4:i+4,j,k)*(vy(i-4:i+4,j,k)+wz(i-4:i+4,j,k)) )
             ! d(mu*du/dy)/dx
             dmuydx = dxinv(1) * first_deriv_8( mu(i-4:i+4,j,k)*uy(i-4:i+4,j,k) )
             ! d(mu*du/dz)/dx
             dmuzdx = dxinv(1) * first_deriv_8( mu(i-4:i+4,j,k)*uz(i-4:i+4,j,k) )
             rhs(i,j,k,imx) = rhs(i,j,k,imx) + dmvywzdx
             rhs(i,j,k,imy) = rhs(i,j,k,imy) + dmuydx
             rhs(i,j,k,imz) = rhs(i,j,k,imz) + dmuzdx
          end do
       end do
    end do
    !$omp end do

    ! d()/dy
    !$omp do
    do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             ! d(mu*dv/dx)/dy
             dmvxdy = dxinv(2) * first_deriv_8( mu(i,j-4:j+4,k)*vx(i,j-4:j+4,k) )
             ! d((xi-2/3*mu)*(ux+wz))/dy
             dmuxwzdy = dxinv(2) * &
                  first_deriv_8( vsm(i,j-4:j+4,k)*(ux(i,j-4:j+4,k)+wz(i,j-4:j+4,k)) )
             ! d(mu*dv/dz)/dy
             dmvzdy = dxinv(2) * first_deriv_8( mu(i,j-4:j+4,k)*vz(i,j-4:j+4,k) )
             rhs(i,j,k,imx) = rhs(i,j,k,imx) + dmvxdy
             rhs(i,j,k,imy) = rhs(i,j,k,imy) + dmuxwzdy
             rhs(i,j,k,imz) = rhs(i,j,k,imz) + dmvzdy
          end do
       end do
    end do
    !$omp end do 
    
    ! d()/dz
    !$omp do
    do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             ! d(mu*dw/dx)/dz
             dmwxdz = dxinv(3) * first_deriv_8( mu(i,j,k-4:k+4)*wx(i,j,k-4:k+4) )
             ! d(mu*dw/dy)/dz
             dmwydz = dxinv(3) * first_deriv_8( mu(i,j,k-4:k+4)*wy(i,j,k-4:k+4) )
             ! d((xi-2/3*mu)*(ux+vy))/dz
             dmuxvydz = dxinv(3) * &
                  first_deriv_8( vsm(i,j,k-4:k+4)*(ux(i,j,k-4:k+4)+vy(i,j,k-4:k+4)) )
             rhs(i,j,k,imx) = rhs(i,j,k,imx) + dmwxdz
             rhs(i,j,k,imy) = rhs(i,j,k,imy) + dmwydz
             rhs(i,j,k,imz) = rhs(i,j,k,imz) + dmuxvydz
          end do
       end do
    end do
    !$omp end do nowait

    !$omp end parallel

    deallocate(ux,uy,uz,vx,vy,vz,wx,wy,wz)

    allocate(dpy(dlo(1):dhi(1),dlo(2):dhi(2),dlo(3):dhi(3),nspecies))
    allocate(dxe(dlo(1):dhi(1),dlo(2):dhi(2),dlo(3):dhi(3),nspecies))
    allocate(dpe(dlo(1):dhi(1),dlo(2):dhi(2),dlo(3):dhi(3)))

    allocate(Hg(lo(1):hi(1)+1,lo(2):hi(2)+1,lo(3):hi(3)+1,2:ncons))

    !$omp parallel &
    !$omp private(i,j,k,n,qxn,qyn,qhn,Htot,Htmp,Ytmp,hhalf,muM8,M8p,M8X,tmp) 

    !$omp workshare
    dpe = 0.d0
    !$omp end workshare

    do n=1,nspecies
       qxn = qx1+n-1
       qyn = qy1+n-1
       qhn = qh1+n-1
       !$OMP DO
       do k=dlo(3),dhi(3)
          do j=dlo(2),dhi(2)
             do i=dlo(1),dhi(1)
                dpy(i,j,k,n) = dxy(i,j,k,n)/q(i,j,k,qpres)*(q(i,j,k,qxn)-q(i,j,k,qyn))
                dxe(i,j,k,n) = dxy(i,j,k,n)*q(i,j,k,qhn)
                dpe(i,j,k) = dpe(i,j,k) + dpy(i,j,k,n)*q(i,j,k,qhn)
             end do
          end do
       end do
       !$omp end do nowait
    end do

    !$omp barrier

    ! ------- BEGIN x-direction -------
    !$omp do
    do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)+1

             muM8 = matmul(   mu(i-4:i+3,j,k)      , M8)
             M8p  = matmul(M8, q(i-4:i+3,j,k,qpres))

             tmp = matmul(M8, q(i-4:i+3,j,k,qu))
             Hg(i,j,k,imx) = dot_product(vsp(i-4:i+3,j,k), tmp)

             Hg(i,j,k,imy) = dot_product(muM8, q(i-4:i+3,j,k,qv))
             Hg(i,j,k,imz) = dot_product(muM8, q(i-4:i+3,j,k,qw))

             tmp = matmul(M8, q(i-4:i+3,j,k,qtemp))
             Hg(i,j,k,iene) = dot_product(lam(i-4:i+3,j,k), tmp) &
                  &         + dot_product(dpe(i-4:i+3,j,k), M8p)

             Htot = 0.d0
             do n = 1, nspecies
                qxn = qx1+n-1
                qyn = qy1+n-1

                M8X = matmul(M8, q(i-4:i+3,j,k,qxn))
                
                Htmp(n) = dot_product(dpy(i-4:i+3,j,k,n), M8p) &
                     +    dot_product(dxy(i-4:i+3,j,k,n), M8X)

                Hg(i,j,k,iene) = Hg(i,j,k,iene) &
                     +    dot_product(dxe(i-4:i+3,j,k,n), M8X)

                Htot = Htot + Htmp(n)
                Ytmp(n) = (q(i-1,j,k,qyn) + q(i,j,k,qyn)) / 2.d0
             end do

             do n = 1, nspecies
                Hg(i,j,k,iry1+n-1) = Htmp(n) - Ytmp(n)*Htot
             end do

             do n = 1, nspecies
                qhn = qh1+n-1
                hhalf = (q(i-1,j,k,qhn) + q(i,j,k,qhn)) / 2.d0
                Hg(i,j,k,iene) =  Hg(i,j,k,iene) - Ytmp(n) * hhalf * Htot
             end do

          end do
       end do
    end do
    !$omp end do

    ! add x-direction rhs
    do n=2,ncons
       !$omp do
       do k=lo(3),hi(3)
          do j=lo(2),hi(2)
             do i=lo(1),hi(1)
                rhs(i,j,k,n) = rhs(i,j,k,n) + (Hg(i+1,j,k,n) - Hg(i,j,k,n)) * dx2inv(1)
             end do
          end do
       end do
       !$omp end do nowait
    end do
    ! ------- END x-direction -------

    !$omp barrier

    ! ------- BEGIN y-direction -------
    !$omp do
    do k=lo(3),hi(3)
       do j=lo(2),hi(2)+1
          do i=lo(1),hi(1)
             
             muM8 = matmul(   mu(i,j-4:j+3,k)      , M8)
             M8p  = matmul(M8, q(i,j-4:j+3,k,qpres))

             Hg(i,j,k,imx) = dot_product(muM8, q(i,j-4:j+3,k,qu))
             Hg(i,j,k,imz) = dot_product(muM8, q(i,j-4:j+3,k,qw))

             tmp = matmul(M8, q(i,j-4:j+3,k,qv))
             Hg(i,j,k,imy) = dot_product(vsp(i,j-4:j+3,k), tmp)

             tmp = matmul(M8, q(i,j-4:j+3,k,qtemp))
             Hg(i,j,k,iene) = dot_product(lam(i,j-4:j+3,k), tmp)      &
                  +           dot_product(dpe(i,j-4:j+3,k), M8p)

             Htot = 0.d0
             do n = 1, nspecies
                qxn = qx1+n-1
                qyn = qy1+n-1

                M8X = matmul(M8, q(i,j-4:j+3,k,qxn))

                Htmp(n) = dot_product(dpy(i,j-4:j+3,k,n), M8P) &
                     +    dot_product(dxy(i,j-4:j+3,k,n), M8X)

                Hg(i,j,k,iene) = Hg(i,j,k,iene) &
                     +    dot_product(dxe(i,j-4:j+3,k,n), M8X)

                Htot = Htot + Htmp(n)
                Ytmp(n) = (q(i,j-1,k,qyn) + q(i,j,k,qyn)) / 2.d0
             end do

             do n = 1, nspecies
                Hg(i,j,k,iry1+n-1) = Htmp(n) - Ytmp(n)*Htot
             end do

             do n = 1, nspecies
                qhn = qh1+n-1
                hhalf = (q(i,j-1,k,qhn) + q(i,j,k,qhn)) / 2.d0
                Hg(i,j,k,iene) =  Hg(i,j,k,iene) - Ytmp(n) * hhalf * Htot
             end do

          end do
       end do
    end do
    !$omp end do

    ! add y-direction rhs
    do n=2,ncons
       !$omp do
       do k=lo(3),hi(3)
          do j=lo(2),hi(2)
             do i=lo(1),hi(1)
                rhs(i,j,k,n) = rhs(i,j,k,n) + (Hg(i,j+1,k,n) - Hg(i,j,k,n)) * dx2inv(2)
             end do
          end do
       end do
       !$omp end do nowait
    end do
    ! ------- END y-direction -------

    !$omp barrier

    ! ------- BEGIN z-direction -------
    !$omp do
    do k=lo(3),hi(3)+1
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)

             muM8 = matmul(   mu(i,j,k-4:k+3)      , M8)
             M8p  = matmul(M8, q(i,j,k-4:k+3,qpres))

             Hg(i,j,k,imx) = dot_product(muM8, q(i,j,k-4:k+3,qu))
             Hg(i,j,k,imy) = dot_product(muM8, q(i,j,k-4:k+3,qv))

             tmp = matmul(M8, q(i,j,k-4:k+3,qw))
             Hg(i,j,k,imz) = dot_product(vsp(i,j,k-4:k+3), tmp)

             tmp = matmul(M8, q(i,j,k-4:k+3,qtemp))
             Hg(i,j,k,iene) = dot_product(lam(i,j,k-4:k+3), tmp)      &
                  +           dot_product(dpe(i,j,k-4:k+3), M8p)

             Htot = 0.d0
             do n = 1, nspecies
                qxn = qx1+n-1
                qyn = qy1+n-1

                M8X = matmul(M8, q(i,j,k-4:k+3,qxn))

                Htmp(n) = dot_product(dpy(i,j,k-4:k+3,n), M8P) &
                     +    dot_product(dxy(i,j,k-4:k+3,n), M8X)

                Hg(i,j,k,iene) = Hg(i,j,k,iene) &
                     +    dot_product(dxe(i,j,k-4:k+3,n), M8X)

                Htot = Htot + Htmp(n)
                Ytmp(n) = (q(i,j,k-1,qyn) + q(i,j,k,qyn)) / 2.d0
             end do

             do n = 1, nspecies
                Hg(i,j,k,iry1+n-1) = Htmp(n) - Ytmp(n)*Htot
             end do

             do n = 1, nspecies
                qhn = qh1+n-1
                hhalf = (q(i,j,k-1,qhn) + q(i,j,k,qhn)) / 2.d0
                Hg(i,j,k,iene) =  Hg(i,j,k,iene) - Ytmp(n) * hhalf * Htot
             end do

          end do
       end do
    end do
    !$omp end do

    ! add z-direction rhs
    do n=2,ncons
       !$omp do
       do k=lo(3),hi(3)
          do j=lo(2),hi(2)
             do i=lo(1),hi(1)
                rhs(i,j,k,n) = rhs(i,j,k,n) + (Hg(i,j,k+1,n) - Hg(i,j,k,n)) * dx2inv(3)
             end do
          end do
       end do
       !$omp end do nowait
    end do
    ! ------- END z-direction -------
    
    !$omp barrier

    ! add kinetic energy
    !$omp do
    do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             rhs(i,j,k,iene) = rhs(i,j,k,iene) &
                  + rhs(i,j,k,imx)*q(i,j,k,qu) &
                  + rhs(i,j,k,imy)*q(i,j,k,qv) &
                  + rhs(i,j,k,imz)*q(i,j,k,qw)
          end do
       end do
    end do
    !$omp end do 
    
    !$omp end parallel

    deallocate(Hg,dpy,dxe,dpe,vsp,vsm)

  end subroutine compact_diffterm_3d


  subroutine chemterm_3d(lo,hi,ng,q,up) ! up is UPrime that has no ghost cells
    integer,          intent(in ) :: lo(3),hi(3),ng
    double precision, intent(in )   :: q (-ng+lo(1):hi(1)+ng,-ng+lo(2):hi(2)+ng,-ng+lo(3):hi(3)+ng,nprim)
    double precision, intent(inout) :: up(    lo(1):hi(1)   ,    lo(2):hi(2)   ,    lo(3):hi(3)   ,ncons)

    integer :: iwrk, i,j,k
    double precision :: Yt(nspecies), wdot(nspecies), rwrk

    !$omp parallel do private(i,j,k,iwrk,rwrk,Yt,wdot)
    do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)

             Yt = q(i,j,k,qy1:qy1+nspecies-1)
             call ckwyr(q(i,j,k,qrho), q(i,j,k,qtemp), Yt, iwrk, rwrk, wdot)
             up(i,j,k,iry1:) = up(i,j,k,iry1:) + wdot * molecular_weight
             
          end do
       end do
    end do
    !$omp end parallel do 

  end subroutine chemterm_3d

  subroutine comp_courno_3d(lo,hi,ng,dx,Q,courno)
    integer, intent(in) :: lo(3), hi(3), ng
    double precision, intent(in) :: dx(3)
    double precision, intent(in) :: q(-ng+lo(1):hi(1)+ng,-ng+lo(2):hi(2)+ng,-ng+lo(3):hi(3)+ng,nprim)
    double precision, intent(inout) :: courno

    integer :: i,j,k, iwrk
    double precision :: dxinv(3), c, rwrk, Cv, Cp
    double precision :: Tt, X(nspecies), gamma
    double precision :: courx, coury, courz

    double precision, parameter :: Ru = 8.31451d7

    do i=1,3
       dxinv(i) = 1.0d0 / dx(i)
    end do

    !$omp parallel do private(i,j,k,iwrk,rwrk,Tt,X,gamma,Cv,Cp,c) &
    !$omp private(courx,coury,courz) &
    !$omp reduction(max:courno)
    do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)

             Tt = q(i,j,k,qtemp)
             X  = q(i,j,k,qx1:qx1+nspecies-1)
             call ckcvbl(Tt, X, iwrk, rwrk, Cv)
             Cp = Cv + Ru
             gamma = Cp / Cv
             c = sqrt(gamma*q(i,j,k,qpres)/q(i,j,k,qrho))

             courx = (c+abs(q(i,j,k,qu))) * dxinv(1)
             coury = (c+abs(q(i,j,k,qv))) * dxinv(2)
             courz = (c+abs(q(i,j,k,qw))) * dxinv(3)
             
             courno = max( courx, coury, courz , courno )

          end do
       end do
    end do
    !$omp end parallel do

  end subroutine comp_courno_3d

  subroutine S3D_diffterm_1(lo,hi,ng,ndq,dx,q,rhs,mu,xi,qx,qy,qz)
 
    integer,          intent(in ) :: lo(3),hi(3),ng,ndq
    double precision, intent(in ) :: dx(3)
    double precision, intent(in ) :: q  (-ng+lo(1):hi(1)+ng,-ng+lo(2):hi(2)+ng,-ng+lo(3):hi(3)+ng,nprim)
    double precision, intent(in ) :: mu (-ng+lo(1):hi(1)+ng,-ng+lo(2):hi(2)+ng,-ng+lo(3):hi(3)+ng)
    double precision, intent(in ) :: xi (-ng+lo(1):hi(1)+ng,-ng+lo(2):hi(2)+ng,-ng+lo(3):hi(3)+ng)
    double precision, intent(out) :: qx (-ng+lo(1):hi(1)+ng,-ng+lo(2):hi(2)+ng,-ng+lo(3):hi(3)+ng,ndq)
    double precision, intent(out) :: qy (-ng+lo(1):hi(1)+ng,-ng+lo(2):hi(2)+ng,-ng+lo(3):hi(3)+ng,ndq)
    double precision, intent(out) :: qz (-ng+lo(1):hi(1)+ng,-ng+lo(2):hi(2)+ng,-ng+lo(3):hi(3)+ng,ndq)
    double precision, intent(out) :: rhs(    lo(1):hi(1)   ,    lo(2):hi(2)   ,    lo(3):hi(3)   ,ncons)

    double precision, allocatable, dimension(:,:,:) :: vsm

    double precision :: dxinv(3), divu
    double precision :: dmvxdy,dmwxdz,dmvywzdx
    double precision :: dmuydx,dmwydz,dmuxwzdy
    double precision :: dmuzdx,dmvzdy,dmuxvydz
    double precision :: tauxx,tauyy,tauzz 
    integer :: i,j,k,n, qxn, qdxn

    allocate(vsm(-ng+lo(1):hi(1)+ng,-ng+lo(2):hi(2)+ng,-ng+lo(3):hi(3)+ng))

    do i = 1,3
       dxinv(i) = 1.0d0 / dx(i)
    end do

    !$omp parallel private(i,j,k,n,qxn,qdxn,divu,tauxx,tauyy,tauzz) &
    !$omp   private(dmvxdy,dmwxdz,dmvywzdx,dmuydx,dmwydz,dmuxwzdy,dmuzdx,dmvzdy,dmuxvydz)

    !$omp workshare
    rhs(:,:,:,irho) = 0.d0
    !$omp end workshare

    !$OMP DO
    do k=lo(3)-ng,hi(3)+ng
       do j=lo(2)-ng,hi(2)+ng
          do i=lo(1)-ng,hi(1)+ng
             vsm(i,j,k) = xi(i,j,k) -  TwoThirds*mu(i,j,k)
          enddo
       enddo
    enddo
    !$OMP END DO NOWAIT

    !$omp do
    do k=lo(3)-ng,hi(3)+ng
       do j=lo(2)-ng,hi(2)+ng
          do i=lo(1),hi(1)
             qx(i,j,k,idu) = dxinv(1) * first_deriv_8( q(i-4:i+4,j,k,qu) )
             qx(i,j,k,idv) = dxinv(1) * first_deriv_8( q(i-4:i+4,j,k,qv) )
             qx(i,j,k,idw) = dxinv(1) * first_deriv_8( q(i-4:i+4,j,k,qw) )
          enddo
       enddo
    enddo
    !$OMP END DO NOWAIT

    !$OMP DO
    do k=lo(3)-ng,hi(3)+ng
       do j=lo(2),hi(2)   
          do i=lo(1)-ng,hi(1)+ng
             qy(i,j,k,idu) = dxinv(2) * first_deriv_8( q(i,j-4:j+4,k,qu) )
             qy(i,j,k,idv) = dxinv(2) * first_deriv_8( q(i,j-4:j+4,k,qv) )
             qy(i,j,k,idw) = dxinv(2) * first_deriv_8( q(i,j-4:j+4,k,qw) )
          enddo
       enddo
    enddo
    !$OMP END DO NOWAIT

    !$OMP DO
    do k=lo(3),hi(3)
       do j=lo(2)-ng,hi(2)+ng
          do i=lo(1)-ng,hi(1)+ng
             qz(i,j,k,idu) = dxinv(3) * first_deriv_8( q(i,j,k-4:k+4,qu) )
             qz(i,j,k,idv) = dxinv(3) * first_deriv_8( q(i,j,k-4:k+4,qv) )
             qz(i,j,k,idw) = dxinv(3) * first_deriv_8( q(i,j,k-4:k+4,qw) )
          enddo
       enddo
    enddo
    !$OMP END DO


    !$omp do
    do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)

             ! d(mu*dv/dx)/dy
             dmvxdy = dxinv(2) * first_deriv_8( mu(i,j-4:j+4,k)*qx(i,j-4:j+4,k,idv) )

             ! d(mu*dw/dx)/dz
             dmwxdz = dxinv(3) * first_deriv_8( mu(i,j,k-4:k+4)*qx(i,j,k-4:k+4,idw) )

             ! d((xi-2/3*mu)*(vy+wz))/dx
             dmvywzdx = dxinv(1) * &
                  first_deriv_8( vsm(i-4:i+4,j,k)*(qy(i-4:i+4,j,k,idv)+qz(i-4:i+4,j,k,idw)) )

             ! d(mu*du/dy)/dx
             dmuydx = dxinv(1) * first_deriv_8( mu(i-4:i+4,j,k)*qy(i-4:i+4,j,k,idu) )

             ! d(mu*dw/dy)/dz
             dmwydz = dxinv(3) * first_deriv_8( mu(i,j,k-4:k+4)*qy(i,j,k-4:k+4,idw) )

             ! d((xi-2/3*mu)*(ux+wz))/dy
             dmuxwzdy = dxinv(2) * &
                  first_deriv_8( vsm(i,j-4:j+4,k)*(qx(i,j-4:j+4,k,idu)+qz(i,j-4:j+4,k,idw)) )

             ! d(mu*du/dz)/dx
             dmuzdx = dxinv(1) * first_deriv_8( mu(i-4:i+4,j,k)*qz(i-4:i+4,j,k,idu) )

             ! d(mu*dv/dz)/dy
             dmvzdy = dxinv(2) * first_deriv_8( mu(i,j-4:j+4,k)*qz(i,j-4:j+4,k,idv) )

             ! d((xi-2/3*mu)*(ux+vy))/dz
             dmuxvydz = dxinv(3) * &
                  first_deriv_8( vsm(i,j,k-4:k+4)*(qx(i,j,k-4:k+4,idu)+qy(i,j,k-4:k+4,idv)) )

             rhs(i,j,k,imx) = dmvxdy + dmwxdz + dmvywzdx
             rhs(i,j,k,imy) = dmuydx + dmwydz + dmuxwzdy
             rhs(i,j,k,imz) = dmuzdx + dmvzdy + dmuxvydz

             divu = (qx(i,j,k,idu)+qy(i,j,k,idv)+qz(i,j,k,idw))*vsm(i,j,k)
             tauxx = 2.d0*mu(i,j,k)*qx(i,j,k,idu) + divu
             tauyy = 2.d0*mu(i,j,k)*qy(i,j,k,idv) + divu
             tauzz = 2.d0*mu(i,j,k)*qz(i,j,k,idw) + divu
             
             ! change in internal energy
             rhs(i,j,k,iene) = tauxx*qx(i,j,k,idu) + tauyy*qy(i,j,k,idv) + tauzz*qz(i,j,k,idw) &
                  + mu(i,j,k)*((qy(i,j,k,idu)+qx(i,j,k,idv))**2 &
                  &          + (qx(i,j,k,idw)+qz(i,j,k,idu))**2 &
                  &          + (qz(i,j,k,idv)+qy(i,j,k,idw))**2 )

          end do
       end do
    end do
    !$omp end do nowait

    !$omp workshare
    rhs(:,:,:,iry1:) = 0.d0
    !$omp end workshare

    !$omp do
    do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)

             qx(i,j,k,idT) = dxinv(1) * first_deriv_8( q(i-4:i+4,j,k,qtemp) )
             qx(i,j,k,idp) = dxinv(1) * first_deriv_8( q(i-4:i+4,j,k,qpres) )

             qy(i,j,k,idT) = dxinv(2) * first_deriv_8( q(i,j-4:j+4,k,qtemp) )
             qy(i,j,k,idp) = dxinv(2) * first_deriv_8( q(i,j-4:j+4,k,qpres) )

             qz(i,j,k,idT) = dxinv(3) * first_deriv_8( q(i,j,k-4:k+4,qtemp) )
             qz(i,j,k,idp) = dxinv(3) * first_deriv_8( q(i,j,k-4:k+4,qpres) )

          enddo
       enddo
    enddo
    !$omp end do nowait

    do n=1,nspecies
       qxn = qx1 + n - 1
       qdxn = idX1 + n -1
       !$omp do
       do k=lo(3),hi(3)
          do j=lo(2),hi(2)
             do i=lo(1),hi(1)
                qx(i,j,k,qdxn) = dxinv(1) * first_deriv_8( q(i-4:i+4,j,k,qxn) )
                qy(i,j,k,qdxn) = dxinv(2) * first_deriv_8( q(i,j-4:j+4,k,qxn) )
                qz(i,j,k,qdxn) = dxinv(3) * first_deriv_8( q(i,j,k-4:k+4,qxn) )
             enddo
          enddo
       enddo
       !$omp end do nowait
    enddo

    !$omp end parallel

    deallocate(vsm)

  end subroutine S3D_diffterm_1


  subroutine S3D_diffterm_2(lo,hi,ng,ndq,dx,q,rhs,mu,xi,lam,dxy,qx,qy,qz)

    integer,          intent(in )  :: lo(3),hi(3),ng,ndq
    double precision, intent(in )  :: dx(3)
    double precision, intent(in )  :: q  (-ng+lo(1):hi(1)+ng,-ng+lo(2):hi(2)+ng,-ng+lo(3):hi(3)+ng,nprim)
    double precision, intent(in )  :: mu (-ng+lo(1):hi(1)+ng,-ng+lo(2):hi(2)+ng,-ng+lo(3):hi(3)+ng)
    double precision, intent(in )  :: xi (-ng+lo(1):hi(1)+ng,-ng+lo(2):hi(2)+ng,-ng+lo(3):hi(3)+ng)
    double precision, intent(in )  :: lam(-ng+lo(1):hi(1)+ng,-ng+lo(2):hi(2)+ng,-ng+lo(3):hi(3)+ng)
    double precision, intent(in )  :: dxy(-ng+lo(1):hi(1)+ng,-ng+lo(2):hi(2)+ng,-ng+lo(3):hi(3)+ng,nspecies)
    double precision, intent(in)   :: qx (-ng+lo(1):hi(1)+ng,-ng+lo(2):hi(2)+ng,-ng+lo(3):hi(3)+ng,ndq)
    double precision, intent(in)   :: qy (-ng+lo(1):hi(1)+ng,-ng+lo(2):hi(2)+ng,-ng+lo(3):hi(3)+ng,ndq)
    double precision, intent(in)   :: qz (-ng+lo(1):hi(1)+ng,-ng+lo(2):hi(2)+ng,-ng+lo(3):hi(3)+ng,ndq)
    double precision, intent(inout):: rhs(    lo(1):hi(1)   ,    lo(2):hi(2)   ,    lo(3):hi(3)   ,ncons)
 
    double precision, allocatable, dimension(:,:,:) :: vp, dpe, FE
    double precision, allocatable, dimension(:,:,:,:) :: dpy, FY
    ! dxy: diffusion coefficient of X in equation for Y
    ! dpy: diffusion coefficient of p in equation for Y
    ! NOT USING ! dxe: diffusion coefficient of X in equation for energy
    ! dpe: diffusion coefficient of p in equation for energy

    double precision :: dxinv(3), rhoVc
    integer          :: i,j,k,n, qxn, qyn, qhn, idXn, iryn

    allocate(vp(-ng+lo(1):hi(1)+ng,-ng+lo(2):hi(2)+ng,-ng+lo(3):hi(3)+ng))

    allocate(dpy(-ng+lo(1):hi(1)+ng,-ng+lo(2):hi(2)+ng,-ng+lo(3):hi(3)+ng,nspecies))
    allocate(dpe(-ng+lo(1):hi(1)+ng,-ng+lo(2):hi(2)+ng,-ng+lo(3):hi(3)+ng))

    allocate(FY(-ng+lo(1):hi(1)+ng,-ng+lo(2):hi(2)+ng,-ng+lo(3):hi(3)+ng,nspecies))
    allocate(FE(-ng+lo(1):hi(1)+ng,-ng+lo(2):hi(2)+ng,-ng+lo(3):hi(3)+ng))

    do i = 1,3
       dxinv(i) = 1.0d0 / dx(i)
    end do

    !$omp parallel private(i,j,k,n,qxn,qyn,qhn,idXn,iryn,rhoVc)

    !$omp do
    do k=lo(3)-ng,hi(3)+ng
       do j=lo(2)-ng,hi(2)+ng
          do i=lo(1)-ng,hi(1)+ng
             vp(i,j,k) = xi(i,j,k) + FourThirds*mu(i,j,k)
          enddo
       enddo
    enddo
    !$omp end do nowait

    !$omp workshare
    dpe = 0.d0
    !$omp end workshare

    do n=1,nspecies
       qxn = qx1+n-1
       qyn = qy1+n-1
       qhn = qh1+n-1
       !$OMP DO
       do k=lo(3)-ng,hi(3)+ng
          do j=lo(2)-ng,hi(2)+ng
             do i=lo(1)-ng,hi(1)+ng
                dpy(i,j,k,n) = dxy(i,j,k,n)/q(i,j,k,qpres)*(q(i,j,k,qxn)-q(i,j,k,qyn))
                dpe(i,j,k) = dpe(i,j,k) + dpy(i,j,k,n)*q(i,j,k,qhn)
             end do
          end do
       end do
       !$omp end do nowait
    end do

    ! ===== mx =====
    !$omp do
    do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             rhs(i,j,k,imx) = rhs(i,j,k,imx) &
                  + dxinv(1) * first_deriv_8( vp(i-4:i+4,j,k)*qx(i-4:i+4,j,k,idu)) &
                  + dxinv(2) * first_deriv_8( mu(i,j-4:j+4,k)*qy(i,j-4:j+4,k,idu)) &
                  + dxinv(3) * first_deriv_8( mu(i,j,k-4:k+4)*qz(i,j,k-4:k+4,idu))
          end do
       end do
    end do
    !$omp end do nowait
    
    ! ===== my =====
    !$omp do
    do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             rhs(i,j,k,imy) = rhs(i,j,k,imy) &
                  + dxinv(1) * first_deriv_8( mu(i-4:i+4,j,k)*qx(i-4:i+4,j,k,idv) ) &
                  + dxinv(2) * first_deriv_8( vp(i,j-4:j+4,k)*qy(i,j-4:j+4,k,idv) ) &
                  + dxinv(3) * first_deriv_8( mu(i,j,k-4:k+4)*qz(i,j,k-4:k+4,idv) ) 
          end do
       end do
    end do
    !$omp end do nowait
    
    ! ===== mz =====
    !$omp do
    do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             rhs(i,j,k,imz) = rhs(i,j,k,imz) &
                  + dxinv(1) * first_deriv_8( mu(i-4:i+4,j,k)*qx(i-4:i+4,j,k,idw) ) &
                  + dxinv(2) * first_deriv_8( mu(i,j-4:j+4,k)*qy(i,j-4:j+4,k,idw) ) &
                  + dxinv(3) * first_deriv_8( vp(i,j,k-4:k+4)*qz(i,j,k-4:k+4,idw) )
          end do
       end do
    end do
    !$omp end do
    
    ! add kinetic energy
    !$omp do
    do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             rhs(i,j,k,iene) = rhs(i,j,k,iene) &
                  + rhs(i,j,k,imx)*q(i,j,k,qu) &
                  + rhs(i,j,k,imy)*q(i,j,k,qv) &
                  + rhs(i,j,k,imz)*q(i,j,k,qw)
          end do
       end do
    end do
    !$omp end do

    ! thermal conduction
    !$omp do
    do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             rhs(i,j,k,iene) = rhs(i,j,k,iene) &
                  + dxinv(1) * first_deriv_8( lam(i-4:i+4,j,k)*qx(i-4:i+4,j,k,idT) ) &
                  + dxinv(2) * first_deriv_8( lam(i,j-4:j+4,k)*qy(i,j-4:j+4,k,idT) ) &
                  + dxinv(3) * first_deriv_8( lam(i,j,k-4:k+4)*qz(i,j,k-4:k+4,idT) )
          end do
       end do
    end do
    !$omp end do nowait

    ! x-direction
    !$omp do
    do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          do i=lo(1)-ng,hi(1)+ng

             rhoVc = 0.d0
             FE(i,j,k) = dpe(i,j,k) * qx(i,j,k,idp)

             do n=1,nspecies
                idXn = idX1+n-1
                qhn = qh1+n-1
                FY(i,j,k,n) = dxy(i,j,k,n)*qx(i,j,k,idXn) + dpy(i,j,k,n)*qx(i,j,k,idp)
                FE(i,j,k) = FE(i,j,k) + dxy(i,j,k,n)*qx(i,j,k,idXn)*q(i,j,k,qhn)
                rhoVc = rhoVc + FY(i,j,k,n)
             end do

             do n=1,nspecies
                qyn = qy1+n-1
                qhn = qh1+n-1
                FY(i,j,k,n) = FY(i,j,k,n) - rhoVc*q(i,j,k,qyn)
                FE(i,j,k) = FE(i,j,k) - rhoVc*q(i,j,k,qyn)*q(i,j,k,qhn)
             end do
          end do
       end do
    end do
    !$omp end do

    do n=1,nspecies    
       iryn = iry1+n-1
       !$omp do
       do k=lo(3),hi(3)
          do j=lo(2),hi(2)
             do i=lo(1),hi(1)
                rhs(i,j,k,iryn) = rhs(i,j,k,iryn) + &
                     dxinv(1) * first_deriv_8( FY(i-4:i+4,j,k,n) )
             end do
          end do
       end do
       !$omp end do nowait
    end do
    
    !$omp do
    do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             rhs(i,j,k,iene) = rhs(i,j,k,iene) + &
                  dxinv(1) * first_deriv_8( FE(i-4:i+4,j,k) )
          end do
       end do
    end do
    !$omp end do

    ! y-direction
    !$omp do
    do k=lo(3),hi(3)
       do j=lo(2)-ng,hi(2)+ng
          do i=lo(1),hi(1)

             rhoVc = 0.d0
             FE(i,j,k) = dpe(i,j,k) * qy(i,j,k,idp)

             do n=1,nspecies
                idXn = idX1+n-1
                qhn = qh1+n-1
                FY(i,j,k,n) = dxy(i,j,k,n)*qy(i,j,k,idXn) + dpy(i,j,k,n)*qy(i,j,k,idp)
                FE(i,j,k) = FE(i,j,k) + dxy(i,j,k,n)*qy(i,j,k,idXn)*q(i,j,k,qhn)
                rhoVc = rhoVc + FY(i,j,k,n)
             end do

             do n=1,nspecies
                qyn = qy1+n-1
                qhn = qh1+n-1
                FY(i,j,k,n) = FY(i,j,k,n) - rhoVc*q(i,j,k,qyn)
                FE(i,j,k) = FE(i,j,k) - rhoVc*q(i,j,k,qyn)*q(i,j,k,qhn)
             end do
          end do
       end do
    end do
    !$omp end do

    do n=1,nspecies    
       iryn = iry1+n-1
       !$omp do
       do k=lo(3),hi(3)
          do j=lo(2),hi(2)
             do i=lo(1),hi(1)
                rhs(i,j,k,iryn) = rhs(i,j,k,iryn) + &
                     dxinv(2) * first_deriv_8( FY(i,j-4:j+4,k,n) )
             end do
          end do
       end do
       !$omp end do nowait
    end do
    
    !$omp do
    do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             rhs(i,j,k,iene) = rhs(i,j,k,iene) + &
                  dxinv(2) * first_deriv_8( FE(i,j-4:j+4,k) )
          end do
       end do
    end do
    !$omp end do

    ! z-direction
    !$omp do
    do k=lo(3)-ng,hi(3)+ng
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)

             rhoVc = 0.d0
             FE(i,j,k) = dpe(i,j,k) * qz(i,j,k,idp)

             do n=1,nspecies
                idXn = idX1+n-1
                qhn = qh1+n-1
                FY(i,j,k,n) = dxy(i,j,k,n)*qz(i,j,k,idXn) + dpy(i,j,k,n)*qz(i,j,k,idp)
                FE(i,j,k) = FE(i,j,k) + dxy(i,j,k,n)*qz(i,j,k,idXn)*q(i,j,k,qhn)
                rhoVc = rhoVc + FY(i,j,k,n)
             end do

             do n=1,nspecies
                qyn = qy1+n-1
                qhn = qh1+n-1
                FY(i,j,k,n) = FY(i,j,k,n) - rhoVc*q(i,j,k,qyn)
                FE(i,j,k) = FE(i,j,k) - rhoVc*q(i,j,k,qyn)*q(i,j,k,qhn)
             end do
          end do
       end do
    end do
    !$omp end do

    do n=1,nspecies    
       iryn = iry1+n-1
       !$omp do
       do k=lo(3),hi(3)
          do j=lo(2),hi(2)
             do i=lo(1),hi(1)
                rhs(i,j,k,iryn) = rhs(i,j,k,iryn) + &
                     dxinv(3) * first_deriv_8( FY(i,j,k-4:k+4,n) )
             end do
          end do
       end do
       !$omp end do nowait
    end do

    !$omp do
    do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             rhs(i,j,k,iene) = rhs(i,j,k,iene) + &
                  dxinv(3) * first_deriv_8( FE(i,j,k-4:k+4) )
          end do
       end do
    end do
    !$omp end do

    !$omp end parallel

    deallocate(vp,dpy,dpe,FY,FE)

  end subroutine S3D_diffterm_2

end module kernels_module
