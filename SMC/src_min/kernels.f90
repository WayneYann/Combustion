module kernels_module
  use chemistry_module, only : nspecies, molecular_weight
  use derivative_stencil_module, only : stencil_ng, first_deriv_8, M8, D8 
  use variables_module
  implicit none

  private

  public :: hypterm_3d, narrow_diffterm_3d, chemterm_3d, comp_courno_3d, &
       S3D_diffterm_1, S3D_diffterm_2

contains

  subroutine hypterm_3d (lo,hi,ng,dx,cons,q,rhs)

    integer,          intent(in ) :: lo(3),hi(3),ng
    double precision, intent(in ) :: dx(3)
    double precision, intent(in ) :: cons(-ng+lo(1):hi(1)+ng,-ng+lo(2):hi(2)+ng,-ng+lo(3):hi(3)+ng,ncons)
    double precision, intent(in ) ::    q(-ng+lo(1):hi(1)+ng,-ng+lo(2):hi(2)+ng,-ng+lo(3):hi(3)+ng,nprim)
    double precision, intent(out) ::  rhs(    lo(1):hi(1)   ,    lo(2):hi(2)   ,    lo(3):hi(3)   ,ncons)

    integer          :: i,j,k,n
    double precision :: dxinv(3)
    double precision, allocatable :: tmpx(:), tmpy(:,:),tmpz(:,:,:)

    integer :: jj
    integer, parameter :: jblocksize=8

    do i=1,3
       dxinv(i) = 1.0d0 / dx(i)
    end do

    allocate(tmpz(lo(1)  :hi(1)  ,lo(2)  :hi(2)  ,lo(3)-4:hi(3)+4))

    !$omp parallel private(i,j,k,n,jj,tmpx,tmpy)

    allocate(tmpx(lo(1)-4:hi(1)+4))
    allocate(tmpy(lo(1)  :hi(1)  ,lo(2)-4:hi(2)+4))
    
    !$omp workshare
    rhs = 0.d0
    !$omp end workshare

    ! ------- BEGIN x-direction -------

    !$omp do
    do k=lo(3),hi(3)
       do j=lo(2),hi(2)

          do i=lo(1),hi(1)
             rhs(i,j,k,irho) = rhs(i,j,k,irho) - dxinv(1) * &
                  first_deriv_8( cons(i-4:i+4,j,k,imx) ) 
          end do

          do i=lo(1)-4,hi(1)+4
             tmpx(i) = cons(i,j,k,imx)*q(i,j,k,qu)+q(i,j,k,qpres)
          end do
          do i=lo(1),hi(1)
             rhs(i,j,k,imx) = rhs(i,j,k,imx) - dxinv(1) * first_deriv_8(tmpx(i-4:i+4))
          end do

          do i=lo(1)-4,hi(1)+4
             tmpx(i) = cons(i,j,k,imy)*q(i,j,k,qu)
          end do
          do i=lo(1),hi(1)
             rhs(i,j,k,imy) = rhs(i,j,k,imy) - dxinv(1) * first_deriv_8(tmpx(i-4:i+4))
          end do

          do i=lo(1)-4,hi(1)+4
             tmpx(i) = cons(i,j,k,imz)*q(i,j,k,qu)
          end do
          do i=lo(1),hi(1)
             rhs(i,j,k,imz) = rhs(i,j,k,imz) - dxinv(1) * first_deriv_8(tmpx(i-4:i+4))
          end do

          do i=lo(1)-4,hi(1)+4
             tmpx(i) = (cons(i,j,k,iene)+q(i,j,k,qpres))*q(i,j,k,qu)
          end do
          do i=lo(1),hi(1)
             rhs(i,j,k,iene) = rhs(i,j,k,iene) - dxinv(1) * first_deriv_8(tmpx(i-4:i+4))
          end do

       end do
    end do
    !$omp end do

    do n = iry1, iry1+nspecies-1
       !$omp do 
       do k=lo(3),hi(3)
          do j=lo(2),hi(2)
    
             do i=lo(1)-4,hi(1)+4
                tmpx(i) = cons(i,j,k,n)*q(i,j,k,qu)
             end do
             do i=lo(1),hi(1)
                rhs(i,j,k,n) = rhs(i,j,k,n) - dxinv(1) * first_deriv_8(tmpx(i-4:i+4))
             end do

          enddo
       enddo
       !$omp end do
    enddo

    ! ------- END x-direction -------

    ! ------- BEGIN y-direction -------

    !$omp do
    do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             rhs(i,j,k,irho)=rhs(i,j,k,irho) - dxinv(2) * &
                  first_deriv_8( cons(i,j-4:j+4,k,imy) )
          enddo
       enddo

       do j=lo(2)-4,hi(2)+4
          do i=lo(1),hi(1)
             tmpy(i,j) = cons(i,j,k,imx)*q(i,j,k,qv)
          end do
       end do
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             rhs(i,j,k,imx)=rhs(i,j,k,imx) - dxinv(2) * first_deriv_8(tmpy(i,j-4:j+4))
          enddo
       enddo

       do j=lo(2)-4,hi(2)+4
          do i=lo(1),hi(1)
             tmpy(i,j) = cons(i,j,k,imy)*q(i,j,k,qv)+q(i,j,k,qpres)
          end do
       end do
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             rhs(i,j,k,imy)=rhs(i,j,k,imy) - dxinv(2) * first_deriv_8(tmpy(i,j-4:j+4))
          enddo
       enddo

       do j=lo(2)-4,hi(2)+4
          do i=lo(1),hi(1)
             tmpy(i,j) = cons(i,j,k,imz)*q(i,j,k,qv)
          end do
       end do
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             rhs(i,j,k,imz)=rhs(i,j,k,imz) - dxinv(2) * first_deriv_8(tmpy(i,j-4:j+4))
          enddo
       enddo

       do j=lo(2)-4,hi(2)+4
          do i=lo(1),hi(1)
             tmpy(i,j) = (cons(i,j,k,iene)+q(i,j,k,qpres))*q(i,j,k,qv)
          end do
       end do
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             rhs(i,j,k,iene)=rhs(i,j,k,iene) - dxinv(2) * first_deriv_8(tmpy(i,j-4:j+4))
          enddo
       enddo
    enddo
    !$omp end do

    do n = iry1, iry1+nspecies-1
       !$omp do
       do k=lo(3),hi(3)
          do j=lo(2)-4,hi(2)+4
             do i=lo(1),hi(1)
                tmpy(i,j) = cons(i,j,k,n)*q(i,j,k,qv)
             end do
          end do
          do j=lo(2),hi(2)
             do i=lo(1),hi(1)
                rhs(i,j,k,n) = rhs(i,j,k,n) - dxinv(2) * first_deriv_8(tmpy(i,j-4:j+4))
             end do
          enddo
       enddo
       !$omp end do
    enddo

    ! ------- END y-direction -------

    ! ------- BEGIN z-direction -------

    !$omp do
    do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             rhs(i,j,k,irho)=rhs(i,j,k,irho) - dxinv(3) * &
                  first_deriv_8( cons(i,j,k-4:k+4,imz) )
          end do
       end do
    end do
    !$omp end do nowait

    !$omp do
    do k=lo(3)-4,hi(3)+4
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             tmpz(i,j,k) = cons(i,j,k,imx) * q(i,j,k,qw)
          end do
       end do
    end do
    !$omp end do
    !$omp do
    do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             rhs(i,j,k,imx)=rhs(i,j,k,imx) - dxinv(3) * first_deriv_8(tmpz(i,j,k-4:k+4))
          end do
       end do
    end do
    !$omp end do

    !$omp do
    do k=lo(3)-4,hi(3)+4
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             tmpz(i,j,k) = cons(i,j,k,imy) * q(i,j,k,qw)
          end do
       end do
    end do
    !$omp end do
    !$omp do
    do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             rhs(i,j,k,imy)=rhs(i,j,k,imy) - dxinv(3) * first_deriv_8(tmpz(i,j,k-4:k+4))
          end do
       end do
    end do
    !$omp end do

    !$omp do
    do k=lo(3)-4,hi(3)+4
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             tmpz(i,j,k) = cons(i,j,k,imz)*q(i,j,k,qw) + q(i,j,k,qpres)
          end do
       end do
    end do
    !$omp end do
    !$omp do
    do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             rhs(i,j,k,imz)=rhs(i,j,k,imz) - dxinv(3) * first_deriv_8(tmpz(i,j,k-4:k+4))
          end do
       end do
    end do
    !$omp end do

    !$omp do
    do k=lo(3)-4,hi(3)+4
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             tmpz(i,j,k) = (cons(i,j,k,iene)+q(i,j,k,qpres))*q(i,j,k,qw)
          end do
       end do
    end do
    !$omp end do
    !$omp do
    do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             rhs(i,j,k,iene)=rhs(i,j,k,iene) - dxinv(3) * first_deriv_8(tmpz(i,j,k-4:k+4))
          end do
       end do
    end do
    !$omp end do

    do n = iry1, iry1+nspecies-1
       do jj=lo(2),hi(2),jblocksize
       !$omp do
       do k=lo(3)-4,hi(3)+4
!          do j=lo(2),hi(2)
          do j=jj,min(jj+jblocksize-1,hi(2))
             do i=lo(1),hi(1)
                tmpz(i,j,k) = cons(i,j,k,n)*q(i,j,k,qw)
             end do
          end do
       end do
       !$omp end do
       end do ! jj

       do jj=lo(2),hi(2),jblocksize
       !$omp do 
       do k=lo(3),hi(3)
!          do j=lo(2),hi(2)
          do j=jj,min(jj+jblocksize-1,hi(2))
             do i=lo(1),hi(1)
                rhs(i,j,k,n) = rhs(i,j,k,n) - dxinv(3) * first_deriv_8(tmpz(i,j,k-4:k+4))
             end do
          enddo
       enddo
       !$omp end do
       end do ! jj
    enddo

    deallocate(tmpx,tmpy)

    !$omp end parallel

    deallocate(tmpz)

  end subroutine hypterm_3d


  subroutine narrow_diffterm_3d (lo,hi,ng,dx,q,rhs,mu,xi,lam,dxy)

    integer,          intent(in ) :: lo(3),hi(3),ng
    double precision, intent(in ) :: dx(3)
    double precision, intent(in ) :: q  (-ng+lo(1):hi(1)+ng,-ng+lo(2):hi(2)+ng,-ng+lo(3):hi(3)+ng,nprim)
    double precision, intent(in ) :: mu (-ng+lo(1):hi(1)+ng,-ng+lo(2):hi(2)+ng,-ng+lo(3):hi(3)+ng)
    double precision, intent(in ) :: xi (-ng+lo(1):hi(1)+ng,-ng+lo(2):hi(2)+ng,-ng+lo(3):hi(3)+ng)
    double precision, intent(in ) :: lam(-ng+lo(1):hi(1)+ng,-ng+lo(2):hi(2)+ng,-ng+lo(3):hi(3)+ng)
    double precision, intent(in ) :: dxy(-ng+lo(1):hi(1)+ng,-ng+lo(2):hi(2)+ng,-ng+lo(3):hi(3)+ng,nspecies)
    double precision, intent(out) :: rhs(    lo(1):hi(1)   ,    lo(2):hi(2)   ,    lo(3):hi(3)   ,ncons)

    double precision, allocatable, dimension(:,:,:) :: ux,uy,uz,vx,vy,vz,wx,wy,wz
    double precision, allocatable :: tmpx(:), tmpy(:,:),tmpz(:,:,:)
    double precision, allocatable, dimension(:,:,:) :: vsp,vsm, dpe
    double precision, allocatable, dimension(:,:,:,:) :: Hg, dpy, dxe
    ! dxy: diffusion coefficient of X in equation for Y
    ! dpy: diffusion coefficient of p in equation for Y
    ! dxe: diffusion coefficient of X in equation for energy
    ! dpe: diffusion coefficient of p in equation for energy

    double precision :: dxinv(3), dx2inv(3), divu
    double precision :: tauxx,tauyy,tauzz 
    integer          :: i,j,k,n, qxn, qyn, qhn
    integer :: dlo(3), dhi(3)

    double precision :: mmtmp(8), Yhalf, hhalf
    double precision, allocatable, dimension(:,:,:,:) :: M8p
    double precision, allocatable, dimension(:,:,:) :: Hry

    integer :: jj 
    integer, parameter :: jblocksize=8

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

    allocate(tmpz( lo(1): hi(1), lo(2): hi(2),dlo(3):dhi(3)))

    !$omp parallel &
    !$omp private(i,j,k,tauxx,tauyy,tauzz,divu,tmpx,tmpy)

    allocate(tmpx(dlo(1):dhi(1)))
    allocate(tmpy( lo(1): hi(1),dlo(2):dhi(2)))

    !$omp workshare
    rhs = 0.d0
    !$omp end workshare

    !$omp do
    do k=dlo(3),dhi(3)
       do j=dlo(2),dhi(2)
          do i=dlo(1),dhi(1)
             vsp(i,j,k) = xi(i,j,k) + FourThirds*mu(i,j,k)
             vsm(i,j,k) = xi(i,j,k) -  TwoThirds*mu(i,j,k)
          enddo
       enddo
    enddo
    !$omp end do nowait

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
    !$omp end do nowait

    !$omp do
    do k=dlo(3),dhi(3)
       do j=lo(2),hi(2)   

          do i=dlo(1),dhi(1)
             uy(i,j,k) = dxinv(2)*first_deriv_8(q(i,j-4:j+4,k,qu))
          enddo

          do i=dlo(1),dhi(1)
             vy(i,j,k) = dxinv(2)*first_deriv_8(q(i,j-4:j+4,k,qv))
          enddo

          do i=dlo(1),dhi(1)
             wy(i,j,k) = dxinv(2)*first_deriv_8(q(i,j-4:j+4,k,qw))
          enddo

       enddo
    enddo
    !$omp end do nowait

    !$omp do
    do k=lo(3),hi(3)
       do j=dlo(2),dhi(2)    
          do i=dlo(1),dhi(1)
             uz(i,j,k) = dxinv(3)*first_deriv_8(q(i,j,k-4:k+4,qu))
             vz(i,j,k) = dxinv(3)*first_deriv_8(q(i,j,k-4:k+4,qv))
             wz(i,j,k) = dxinv(3)*first_deriv_8(q(i,j,k-4:k+4,qw))
          end do      
       end do
    end do
    !$omp end do

    !----- mx -----

    !----- mx : d()/dx -----
    !$omp do
    do k=lo(3),hi(3)
       do j=lo(2),hi(2)

          do i=lo(1)-4,hi(1)+4
             tmpx(i) = vsm(i,j,k)*(vy(i,j,k)+wz(i,j,k))
          end do

          do i=lo(1),hi(1)
             rhs(i,j,k,imx) = rhs(i,j,k,imx) + dxinv(1) * first_deriv_8(tmpx(i-4:i+4)) 
          end do

       end do
    end do
    !$omp end do

    !----- mx : d()/dy -----
    !$omp do
    do k=lo(3),hi(3)

       do j=lo(2)-4,hi(2)+4
          do i=lo(1),hi(1)
             tmpy(i,j) = mu(i,j,k)*vx(i,j,k)
          end do
       end do

       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             rhs(i,j,k,imx) = rhs(i,j,k,imx) + dxinv(2) * first_deriv_8(tmpy(i,j-4:j+4)) 
          end do
       end do

    end do
    !$omp end do

    !----- mx : d()/dz -----
    !$omp do
    do k=lo(3)-4,hi(3)+4
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             tmpz(i,j,k) = mu(i,j,k)*wx(i,j,k)
          end do
       end do
    end do
    !$omp end do
    !$omp do
    do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             rhs(i,j,k,imx) = rhs(i,j,k,imx) + dxinv(3) * first_deriv_8(tmpz(i,j,k-4:k+4)) 
          end do
       end do
    end do
    !$omp end do

    !----- my -----

    !$omp do
    do k=lo(3),hi(3)
       do j=lo(2),hi(2)

          do i=lo(1)-4,hi(1)+4
             tmpx(i) = mu(i,j,k)*uy(i,j,k)
          end do

          do i=lo(1),hi(1)
             rhs(i,j,k,imy) = rhs(i,j,k,imy) + dxinv(1) * first_deriv_8(tmpx(i-4:i+4)) 
          end do

       end do
    end do
    !$omp end do

    !$omp do
    do k=lo(3),hi(3)

       do j=lo(2)-4,hi(2)+4
          do i=lo(1),hi(1)
             tmpy(i,j) = vsm(i,j,k)*(ux(i,j,k)+wz(i,j,k))
          end do
       end do

       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             rhs(i,j,k,imy) = rhs(i,j,k,imy) + dxinv(2) * first_deriv_8(tmpy(i,j-4:j+4)) 
          end do
       end do

    end do
    !$omp end do

    !$omp do
    do k=lo(3)-4,hi(3)+4
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             tmpz(i,j,k) = mu(i,j,k)*wy(i,j,k)
          end do
       end do
    end do
    !$omp end do
    !$omp do
    do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             rhs(i,j,k,imy) = rhs(i,j,k,imy) + dxinv(3) * first_deriv_8(tmpz(i,j,k-4:k+4)) 
          end do
       end do
    end do
    !$omp end do

    !----- mz -----

    !$omp do
    do k=lo(3),hi(3)
       do j=lo(2),hi(2)

          do i=lo(1)-4,hi(1)+4
             tmpx(i) = mu(i,j,k)*uz(i,j,k)
          end do

          do i=lo(1),hi(1)
             rhs(i,j,k,imz) = rhs(i,j,k,imz) + dxinv(1) * first_deriv_8(tmpx(i-4:i+4))
          end do

       end do
    end do
    !$omp end do

    !$omp do
    do k=lo(3),hi(3)

       do j=lo(2)-4,hi(2)+4
          do i=lo(1),hi(1)
             tmpy(i,j) = mu(i,j,k)*vz(i,j,k)
          end do
       end do

       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             rhs(i,j,k,imz) = rhs(i,j,k,imz) + dxinv(2) * first_deriv_8(tmpy(i,j-4:j+4))
          end do
       end do

    end do
    !$omp end do 

    !$omp do
    do k=lo(3)-4,hi(3)+4
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             tmpz(i,j,k) = vsm(i,j,k)*(ux(i,j,k)+vy(i,j,k))
          end do
       end do
    end do
    !$omp end do
    !$omp do
    do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             rhs(i,j,k,imz) = rhs(i,j,k,imz) + dxinv(3) * first_deriv_8(tmpz(i,j,k-4:k+4))
          end do
       end do
    end do
    !$omp end do nowait

    !----- energy -----

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
    !$omp end do

    deallocate(tmpx,tmpy)

    !$omp end parallel

    deallocate(ux,uy,uz,vx,vy,vz,wx,wy,wz,tmpz)

    allocate(dpy(dlo(1):dhi(1),dlo(2):dhi(2),dlo(3):dhi(3),nspecies))
    allocate(dxe(dlo(1):dhi(1),dlo(2):dhi(2),dlo(3):dhi(3),nspecies))
    allocate(dpe(dlo(1):dhi(1),dlo(2):dhi(2),dlo(3):dhi(3)))

    allocate(Hg(lo(1):hi(1)+1,lo(2):hi(2)+1,lo(3):hi(3)+1,2:ncons))

    allocate(M8p(8,lo(1):hi(1)+1,lo(2):hi(2)+1,lo(3):hi(3)+1))
    allocate(Hry(  lo(1):hi(1)+1,lo(2):hi(2)+1,lo(3):hi(3)+1))

    !$omp parallel &
    !$omp private(i,j,k,n,qxn,qyn,qhn,Yhalf,hhalf,mmtmp,jj)

    !$omp workshare
    dpe = 0.d0
    !$omp end workshare

    do n=1,nspecies
       qxn = qx1+n-1
       qyn = qy1+n-1
       qhn = qh1+n-1
       !$omp do
       do k=dlo(3),dhi(3)
          do j=dlo(2),dhi(2)
             do i=dlo(1),dhi(1)
                dpy(i,j,k,n) = dxy(i,j,k,n)/q(i,j,k,qpres)*(q(i,j,k,qxn)-q(i,j,k,qyn))
                dxe(i,j,k,n) = dxy(i,j,k,n)*q(i,j,k,qhn)
                dpe(i,j,k) = dpe(i,j,k) + dpy(i,j,k,n)*q(i,j,k,qhn)
             end do
          end do
       end do
       !$omp end do
    end do

    ! ------- BEGIN x-direction -------

    !$omp do
    do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)+1
             mmtmp = matmul(vsp(i-4:i+3,j,k), M8)
             Hg(i,j,k,imx) = dot_product(mmtmp, q(i-4:i+3,j,k,qu))

             mmtmp = matmul(mu(i-4:i+3,j,k), M8)
             Hg(i,j,k,imy) = dot_product(mmtmp, q(i-4:i+3,j,k,qv))
             Hg(i,j,k,imz) = dot_product(mmtmp, q(i-4:i+3,j,k,qw))

             mmtmp = matmul(lam(i-4:i+3,j,k), M8)
             Hg(i,j,k,iene) = dot_product(mmtmp, q(i-4:i+3,j,k,qtemp))

             mmtmp = matmul(M8, q(i-4:i+3,j,k,qpres))
             Hg(i,j,k,iene) = Hg(i,j,k,iene) + dot_product(dpe(i-4:i+3,j,k), mmtmp)
             M8p(:,i,j,k) = mmtmp
          end do
       end do
    end do
    !$omp end do

    do n = 1, nspecies
       qxn = qx1+n-1
       !$omp do
       do k=lo(3),hi(3)
          do j=lo(2),hi(2)
   
             do i=lo(1),hi(1)+1
                Hg(i,j,k,iry1+n-1) = dot_product(dpy(i-4:i+3,j,k,n), M8p(:,i,j,k))
             end do
    
             do i=lo(1),hi(1)+1
                mmtmp = matmul(M8, q(i-4:i+3,j,k,qxn))
                Hg(i,j,k,iene) = Hg(i,j,k,iene) + dot_product(dxe(i-4:i+3,j,k,n), mmtmp)
                Hg(i,j,k,iry1+n-1) = Hg(i,j,k,iry1+n-1) &
                     + dot_product(dxy(i-4:i+3,j,k,n), mmtmp)
             end do

          end do
       end do
       !$omp end do
    end do

    ! correction
    
    !$omp workshare
    Hry = 0.d0
    !$omp end workshare

    do n = 1, nspecies
       !$omp do
       do k=lo(3),hi(3)
          do j=lo(2),hi(2)
             do i=lo(1),hi(1)+1
                Hry(i,j,k) = Hry(i,j,k) + Hg(i,j,k,iry1+n-1)
             end do
          end do
       end do
       !$omp end do
    end do
    
    do n = 1, nspecies
       qyn = qy1+n-1
       qhn = qh1+n-1
       !$omp do
       do k=lo(3),hi(3)
          do j=lo(2),hi(2)
             do i=lo(1),hi(1)+1
                Yhalf = 0.5d0*(q(i-1,j,k,qyn) + q(i,j,k,qyn))
                hhalf = 0.5d0*(q(i-1,j,k,qhn) + q(i,j,k,qhn))
                Hg(i,j,k,iry1+n-1) = Hg(i,j,k,iry1+n-1)- (Yhalf*Hry(i,j,k))
                Hg(i,j,k,iene) = Hg(i,j,k,iene) - (Yhalf*Hry(i,j,k))*hhalf
             end do
          end do
       end do
       !$omp end do
    end do

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
             mmtmp = matmul(mu(i,j-4:j+3,k), M8)
             Hg(i,j,k,imx) = dot_product(mmtmp, q(i,j-4:j+3,k,qu))
             Hg(i,j,k,imz) = dot_product(mmtmp, q(i,j-4:j+3,k,qw))
          end do

          do i=lo(1),hi(1)
             mmtmp = matmul(vsp(i,j-4:j+3,k), M8)
             Hg(i,j,k,imy) = dot_product(mmtmp, q(i,j-4:j+3,k,qv))
          end do

          do i=lo(1),hi(1)
             mmtmp = matmul(lam(i,j-4:j+3,k), M8)
             Hg(i,j,k,iene) = dot_product(mmtmp, q(i,j-4:j+3,k,qtemp))
          end do

          do i=lo(1),hi(1)
             mmtmp = matmul(M8, q(i,j-4:j+3,k,qpres))
             M8p(:,i,j,k) = mmtmp
             Hg(i,j,k,iene) = Hg(i,j,k,iene) + dot_product(dpe(i,j-4:j+3,k), mmtmp)
          end do

       end do
    end do
    !$omp end do

    do n = 1, nspecies
       qxn = qx1+n-1

      !$omp do
       do k=lo(3),hi(3)
          do j=lo(2),hi(2)+1

             do i=lo(1),hi(1)
                Hg(i,j,k,iry1+n-1) = dot_product(dpy(i,j-4:j+3,k,n), M8p(:,i,j,k))
             end do

             do i=lo(1),hi(1)
                mmtmp = matmul(M8, q(i,j-4:j+3,k,qxn))
                Hg(i,j,k,iene) = Hg(i,j,k,iene) + dot_product(dxe(i,j-4:j+3,k,n), mmtmp)
                Hg(i,j,k,iry1+n-1) = Hg(i,j,k,iry1+n-1) &
                     + dot_product(dxy(i,j-4:j+3,k,n), mmtmp)
             end do

          end do
       end do

    end do
       
    ! correction

    !$omp workshare
    Hry = 0.d0
    !$omp end workshare

    do n = 1, nspecies
       !$omp do
       do k=lo(3),hi(3)
          do j=lo(2),hi(2)+1
             do i=lo(1),hi(1)
                Hry(i,j,k) = Hry(i,j,k) + Hg(i,j,k,iry1+n-1)
             end do
          end do
       end do
       !$omp end do
    end do

    do n = 1, nspecies
       qyn = qy1+n-1
       qhn = qh1+n-1
       !$omp do
       do k=lo(3),hi(3)
          do j=lo(2),hi(2)+1
             do i=lo(1),hi(1)
                Yhalf = 0.5d0*(q(i,j-1,k,qyn) + q(i,j,k,qyn))
                hhalf = 0.5d0*(q(i,j-1,k,qhn) + q(i,j,k,qhn))
                Hg(i,j,k,iry1+n-1) = Hg(i,j,k,iry1+n-1)- (Yhalf*Hry(i,j,k))
                Hg(i,j,k,iene) = Hg(i,j,k,iene) - (Yhalf*Hry(i,j,k))*hhalf
             end do
          end do
       end do
       !$omp end do
    end do

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
             mmtmp = matmul(mu(i,j,k-4:k+3), M8)
             Hg(i,j,k,imx) = dot_product(mmtmp, q(i,j,k-4:k+3,qu))
             Hg(i,j,k,imy) = dot_product(mmtmp, q(i,j,k-4:k+3,qv))
          end do
       end do
    end do
    !$omp end do

    !$omp do
    do k=lo(3),hi(3)+1
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             mmtmp = matmul(vsp(i,j,k-4:k+3), M8)
             Hg(i,j,k,imz) = dot_product(mmtmp, q(i,j,k-4:k+3,qw))
          end do
       end do
    end do
    !$omp end do

    !$omp do
    do k=lo(3),hi(3)+1
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             mmtmp = matmul(lam(i,j,k-4:k+3), M8)
             Hg(i,j,k,iene) = dot_product(mmtmp, q(i,j,k-4:k+3,qtemp))
          end do
       end do
    end do
    !$omp end do

    !$omp do
    do k=lo(3),hi(3)+1
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             mmtmp = matmul(M8, q(i,j,k-4:k+3,qpres))
             M8p(:,i,j,k) = mmtmp
             Hg(i,j,k,iene) = Hg(i,j,k,iene) + dot_product(dpe(i,j,k-4:k+3), mmtmp)
          end do
       end do
    end do
    !$omp end do

    do n = 1, nspecies
       qxn = qx1+n-1

       do jj=lo(2),hi(2),jblocksize
       !$omp do
       do k=lo(3),hi(3)+1
!          do j=lo(2),hi(2)
          do j=jj,min(jj+jblocksize-1,hi(2))
             do i=lo(1),hi(1)
                Hg(i,j,k,iry1+n-1) = dot_product(dpy(i,j,k-4:k+3,n), M8p(:,i,j,k))
             end do
          end do
       end do
       !$omp end do
       enddo

       do jj=lo(2),hi(2),jblocksize
       !$omp do
       do k=lo(3),hi(3)+1
!          do j=lo(2),hi(2)
          do j=jj,min(jj+jblocksize-1,hi(2))
             do i=lo(1),hi(1)
                mmtmp = matmul(M8, q(i,j,k-4:k+3,qxn))
                Hg(i,j,k,iene) = Hg(i,j,k,iene) + dot_product(dxe(i,j,k-4:k+3,n), mmtmp)
                Hg(i,j,k,iry1+n-1) = Hg(i,j,k,iry1+n-1) &
                     + dot_product(dxy(i,j,k-4:k+3,n), mmtmp)
             end do
          end do
       end do
       !$omp end do
       enddo
    end do

    ! correction

    !$omp workshare
    Hry = 0.d0
    !$omp end workshare

    do n = 1, nspecies
       !$omp do
       do k=lo(3),hi(3)+1
          do j=lo(2),hi(2)
             do i=lo(1),hi(1)
                Hry(i,j,k) = Hry(i,j,k) + Hg(i,j,k,iry1+n-1)
             end do
          end do
       end do
       !$omp end do
    end do

    do n = 1, nspecies
       qyn = qy1+n-1
       qhn = qh1+n-1
       !$omp do
       do k=lo(3),hi(3)+1
          do j=lo(2),hi(2)
             do i=lo(1),hi(1)
                Yhalf = 0.5d0*(q(i,j,k-1,qyn) + q(i,j,k,qyn))
                hhalf = 0.5d0*(q(i,j,k-1,qhn) + q(i,j,k,qhn))
                Hg(i,j,k,iry1+n-1) = Hg(i,j,k,iry1+n-1)- (Yhalf*Hry(i,j,k))
                Hg(i,j,k,iene) = Hg(i,j,k,iene) - (Yhalf*Hry(i,j,k))*hhalf
             end do
          end do
       end do
       !$omp end do
    end do
    
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

    deallocate(Hg,dpy,dxe,dpe,vsp,vsm,M8p,Hry)

  end subroutine narrow_diffterm_3d


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
             up(i,j,k,iry1:) = wdot * molecular_weight
             
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
    double precision, allocatable :: tmpx(:), tmpy(:,:),tmpz(:,:,:)

    double precision :: dxinv(3), divu
    double precision :: tauxx,tauyy,tauzz 
    integer :: i,j,k,n, qxn, qdxn

    allocate(vsm (-ng+lo(1):hi(1)+ng,-ng+lo(2):hi(2)+ng,-ng+lo(3):hi(3)+ng))
    allocate(tmpz(    lo(1):hi(1)   ,    lo(2):hi(2)   ,-ng+lo(3):hi(3)+ng))

    do i = 1,3
       dxinv(i) = 1.0d0 / dx(i)
    end do

    !$omp parallel private(i,j,k,n,qxn,qdxn,divu,tauxx,tauyy,tauzz,tmpx,tmpy)

    allocate(tmpx(-ng+lo(1):hi(1)+ng))
    allocate(tmpy(    lo(1):hi(1)   ,-ng+lo(2):hi(2)+ng))

    !$omp workshare
    rhs = 0.d0
    !$omp end workshare

    !$omp do
    do k=lo(3)-ng,hi(3)+ng
       do j=lo(2)-ng,hi(2)+ng
          do i=lo(1)-ng,hi(1)+ng
             vsm(i,j,k) = xi(i,j,k) -  TwoThirds*mu(i,j,k)
          enddo
       enddo
    enddo
    !$omp end do nowait

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
    !$omp end do nowait

    !$omp do
    do k=lo(3)-ng,hi(3)+ng
       do j=lo(2),hi(2)   
          do i=lo(1)-ng,hi(1)+ng
             qy(i,j,k,idu) = dxinv(2) * first_deriv_8( q(i,j-4:j+4,k,qu) )
          enddo
       enddo
    enddo
    !$omp end do nowait

    !$omp do
    do k=lo(3)-ng,hi(3)+ng
       do j=lo(2),hi(2)   
          do i=lo(1)-ng,hi(1)+ng
             qy(i,j,k,idv) = dxinv(2) * first_deriv_8( q(i,j-4:j+4,k,qv) )
          enddo
       enddo
    enddo
    !$omp end do nowait

    !$omp do
    do k=lo(3)-ng,hi(3)+ng
       do j=lo(2),hi(2)   
          do i=lo(1)-ng,hi(1)+ng
             qy(i,j,k,idw) = dxinv(2) * first_deriv_8( q(i,j-4:j+4,k,qw) )
          enddo
       enddo
    enddo
    !$omp end do nowait

    !$omp do
    do k=lo(3),hi(3)
       do j=lo(2)-ng,hi(2)+ng
          do i=lo(1)-ng,hi(1)+ng
             qz(i,j,k,idu) = dxinv(3) * first_deriv_8( q(i,j,k-4:k+4,qu) )
             qz(i,j,k,idv) = dxinv(3) * first_deriv_8( q(i,j,k-4:k+4,qv) )
             qz(i,j,k,idw) = dxinv(3) * first_deriv_8( q(i,j,k-4:k+4,qw) )
          enddo
       enddo
    enddo
    !$omp end do

    !----- mx -----

    !$omp do
    do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          do i=lo(1)-4,hi(1)+4
             tmpx(i) = vsm(i,j,k)*(qy(i,j,k,idv)+qz(i,j,k,idw))
          end do
          do i=lo(1),hi(1)
             rhs(i,j,k,imx) = rhs(i,j,k,imx) + dxinv(1) * first_deriv_8(tmpx(i-4:i+4))
          end do
       end do
    end do
    !$omp end do

    !$omp do
    do k=lo(3),hi(3)
       do j=lo(2)-4,hi(2)+4
          do i=lo(1),hi(1)
             tmpy(i,j) = mu(i,j,k)*qx(i,j,k,idv)
          end do
       end do
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             rhs(i,j,k,imx) = rhs(i,j,k,imx) + dxinv(2) * first_deriv_8(tmpy(i,j-4:j+4))
          end do
       end do
    end do
    !$omp end do

    !$omp do
    do k=lo(3)-4,hi(3)+4
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             tmpz(i,j,k) = mu(i,j,k)*qx(i,j,k,idw)
          end do
       end do
    end do
    !$omp end do
    !$omp do
    do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             rhs(i,j,k,imx) = rhs(i,j,k,imx) + dxinv(3) * first_deriv_8(tmpz(i,j,k-4:k+4))
          end do
       end do
    end do
    !$omp end do

    !----- my -----

    !$omp do
    do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          do i=lo(1)-4,hi(1)+4
             tmpx(i) = mu(i,j,k)*qy(i,j,k,idu)
          end do
          do i=lo(1),hi(1)
             rhs(i,j,k,imy) = rhs(i,j,k,imy) + dxinv(1) * first_deriv_8(tmpx(i-4:i+4))
          end do
       end do
    end do
    !$omp end do

    !$omp do
    do k=lo(3),hi(3)
       do j=lo(2)-4,hi(2)+4
          do i=lo(1),hi(1)
             tmpy(i,j) = vsm(i,j,k)*(qx(i,j,k,idu)+qz(i,j,k,idw))
          end do
       end do
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             rhs(i,j,k,imy) = rhs(i,j,k,imy) + dxinv(2) * first_deriv_8(tmpy(i,j-4:j+4))
          end do
       end do
    end do
    !$omp end do

    !$omp do
    do k=lo(3)-4,hi(3)+4
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             tmpz(i,j,k) = mu(i,j,k)*qy(i,j,k,idw)
          end do
       end do
    end do
    !$omp end do
    !$omp do
    do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             rhs(i,j,k,imy) = rhs(i,j,k,imy) + dxinv(3) * first_deriv_8(tmpz(i,j,k-4:k+4))
          end do
       end do
    end do
    !$omp end do

    !----- mz -----

    !$omp do
    do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          do i=lo(1)-4,hi(1)+4
             tmpx(i) = mu(i,j,k)*qz(i,j,k,idu)
          end do
          do i=lo(1),hi(1)
             rhs(i,j,k,imz) = rhs(i,j,k,imz) + dxinv(1) * first_deriv_8(tmpx(i-4:i+4))
          end do
       end do
    end do
    !$omp end do

    !$omp do
    do k=lo(3),hi(3)
       do j=lo(2)-4,hi(2)+4
          do i=lo(1),hi(1)
             tmpy(i,j) = mu(i,j,k)*qz(i,j,k,idv)
          end do
       end do
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             rhs(i,j,k,imz) = rhs(i,j,k,imz) + dxinv(2) * first_deriv_8(tmpy(i,j-4:j+4))
          end do
       end do
    end do
    !$omp end do

    !$omp do
    do k=lo(3)-4,hi(3)+4
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             tmpz(i,j,k) = vsm(i,j,k)*(qx(i,j,k,idu)+qy(i,j,k,idv))
          end do
       end do
    end do
    !$omp end do
    !$omp do
    do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             rhs(i,j,k,imz) = rhs(i,j,k,imz) + dxinv(3) * first_deriv_8(tmpz(i,j,k-4:k+4))
          end do
       end do
    end do
    !$omp end do nowait

    !----- energy -----

    !$omp do
    do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             divu = (qx(i,j,k,idu)+qy(i,j,k,idv)+qz(i,j,k,idw))*vsm(i,j,k)
             tauxx = 2.d0*mu(i,j,k)*qx(i,j,k,idu) + divu
             tauyy = 2.d0*mu(i,j,k)*qy(i,j,k,idv) + divu
             tauzz = 2.d0*mu(i,j,k)*qz(i,j,k,idw) + divu

             ! change in internal energy
             rhs(i,j,k,iene) = tauxx*qx(i,j,k,idu) + tauyy*qy(i,j,k,idv) + tauzz*qz(i,j,k,idw) &
                  + mu(i,j,k)*((qy(i,j,k,idu)+qx(i,j,k,idv))**2 &
                  &          + (qx(i,j,k,idw)+qz(i,j,k,idu))**2 &
                  &          + (qz(i,j,k,idv)+qy(i,j,k,idw))**2 )
          enddo
       enddo
    enddo
    !$omp end do nowait

    !$omp do
    do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             qx(i,j,k,idT) = dxinv(1) * first_deriv_8( q(i-4:i+4,j,k,qtemp) )
          enddo
       enddo
    enddo
    !$omp end do nowait

    !$omp do
    do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             qy(i,j,k,idT) = dxinv(2) * first_deriv_8( q(i,j-4:j+4,k,qtemp) )
          enddo
       enddo
    enddo
    !$omp end do nowait

    !$omp do
    do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             qz(i,j,k,idT) = dxinv(3) * first_deriv_8( q(i,j,k-4:k+4,qtemp) )
          enddo
       enddo
    enddo
    !$omp end do nowait

    !$omp do
    do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             qx(i,j,k,idp) = dxinv(1) * first_deriv_8( q(i-4:i+4,j,k,qpres) )
          enddo
       enddo
    enddo
    !$omp end do nowait

    !$omp do
    do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             qy(i,j,k,idp) = dxinv(2) * first_deriv_8( q(i,j-4:j+4,k,qpres) )
          enddo
       enddo
    enddo
    !$omp end do nowait

    !$omp do
    do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
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
             enddo
          enddo
       enddo
       !$omp end do nowait
       !$omp do
       do k=lo(3),hi(3)
          do j=lo(2),hi(2)
             do i=lo(1),hi(1)
                qy(i,j,k,qdxn) = dxinv(2) * first_deriv_8( q(i,j-4:j+4,k,qxn) )
             enddo
          enddo
       enddo
       !$omp end do nowait
       !$omp do
       do k=lo(3),hi(3)
          do j=lo(2),hi(2)
             do i=lo(1),hi(1)
                qz(i,j,k,qdxn) = dxinv(3) * first_deriv_8( q(i,j,k-4:k+4,qxn) )
             enddo
          enddo
       enddo
       !$omp end do nowait
    enddo

    deallocate(tmpx,tmpy)

    !$omp end parallel

    deallocate(vsm,tmpz)

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

    double precision :: dxinv(3)
    integer          :: i,j,k,n, qxn, qyn, qhn, idXn, iryn
    double precision, allocatable, dimension(:,:,:) :: rvc
    double precision, allocatable :: tmpx(:), tmpy(:,:),tmpz(:,:,:)

    allocate(vp(-ng+lo(1):hi(1)+ng,-ng+lo(2):hi(2)+ng,-ng+lo(3):hi(3)+ng))

    allocate(dpy(-ng+lo(1):hi(1)+ng,-ng+lo(2):hi(2)+ng,-ng+lo(3):hi(3)+ng,nspecies))
    allocate(dpe(-ng+lo(1):hi(1)+ng,-ng+lo(2):hi(2)+ng,-ng+lo(3):hi(3)+ng))

    allocate(FY(-ng+lo(1):hi(1)+ng,-ng+lo(2):hi(2)+ng,-ng+lo(3):hi(3)+ng,nspecies))
    allocate(FE(-ng+lo(1):hi(1)+ng,-ng+lo(2):hi(2)+ng,-ng+lo(3):hi(3)+ng))

    allocate(rvc (-ng+lo(1):hi(1)+ng,-ng+lo(2):hi(2)+ng,-ng+lo(3):hi(3)+ng))
    allocate(tmpz(    lo(1):hi(1)   ,    lo(2):hi(2)   ,-ng+lo(3):hi(3)+ng))

    do i = 1,3
       dxinv(i) = 1.0d0 / dx(i)
    end do

    !$omp parallel private(i,j,k,n,qxn,qyn,qhn,idXn,iryn,tmpx,tmpy)

    allocate(tmpx(-ng+lo(1):hi(1)+ng))
    allocate(tmpy(    lo(1):hi(1)   ,-ng+lo(2):hi(2)+ng))

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
       !$omp do
       do k=lo(3)-ng,hi(3)+ng
          do j=lo(2)-ng,hi(2)+ng
             do i=lo(1)-ng,hi(1)+ng
                dpy(i,j,k,n) = dxy(i,j,k,n)/q(i,j,k,qpres)*(q(i,j,k,qxn)-q(i,j,k,qyn))
                dpe(i,j,k) = dpe(i,j,k) + dpy(i,j,k,n)*q(i,j,k,qhn)
             end do
          end do
       end do
       !$omp end do
    end do

    ! x-direction

    !$omp do
    do k=lo(3),hi(3)
       do j=lo(2),hi(2)

          do i=lo(1)-4,hi(1)+4
             tmpx(i) = vp(i,j,k)*qx(i,j,k,idu) 
          end do
          do i=lo(1),hi(1)
             rhs(i,j,k,imx) = rhs(i,j,k,imx) + dxinv(1) * first_deriv_8(tmpx(i-4:i+4))
          end do

          do i=lo(1)-4,hi(1)+4
             tmpx(i) = mu(i,j,k)*qx(i,j,k,idv)
          end do
          do i=lo(1),hi(1)
             rhs(i,j,k,imy) = rhs(i,j,k,imy) + dxinv(1) * first_deriv_8(tmpx(i-4:i+4))
          end do

          do i=lo(1)-4,hi(1)+4
             tmpx(i) = mu(i,j,k)*qx(i,j,k,idw)
          end do
          do i=lo(1),hi(1)
             rhs(i,j,k,imz) = rhs(i,j,k,imz) + dxinv(1) * first_deriv_8(tmpx(i-4:i+4))
          end do

          do i=lo(1)-4,hi(1)+4
             tmpx(i) = lam(i,j,k)*qx(i,j,k,idT)
          end do
          do i=lo(1),hi(1)
             rhs(i,j,k,iene) = rhs(i,j,k,iene) + dxinv(1) * first_deriv_8(tmpx(i-4:i+4))
          end do

       end do
    end do
    !$omp end do

    ! y-direction

    !$omp do
    do k=lo(3),hi(3)

       do j=lo(2)-4,hi(2)+4
          do i=lo(1),hi(1)
             tmpy(i,j) = mu(i,j,k)*qy(i,j,k,idu)
          end do
       end do
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             rhs(i,j,k,imx) = rhs(i,j,k,imx) + dxinv(2) * first_deriv_8(tmpy(i,j-4:j+4))
          end do
       end do

       do j=lo(2)-4,hi(2)+4
          do i=lo(1),hi(1)
             tmpy(i,j) = vp(i,j,k)*qy(i,j,k,idv)
          end do
       end do
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             rhs(i,j,k,imy) = rhs(i,j,k,imy) + dxinv(2) * first_deriv_8(tmpy(i,j-4:j+4))
          end do
       end do

       do j=lo(2)-4,hi(2)+4
          do i=lo(1),hi(1)
             tmpy(i,j) = mu(i,j,k)*qy(i,j,k,idw)
          end do
       end do
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             rhs(i,j,k,imz) = rhs(i,j,k,imz) + dxinv(2) * first_deriv_8(tmpy(i,j-4:j+4))
          end do
       end do

       do j=lo(2)-4,hi(2)+4
          do i=lo(1),hi(1)
             tmpy(i,j) = lam(i,j,k)*qy(i,j,k,idT)
          end do
       end do
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             rhs(i,j,k,iene) = rhs(i,j,k,iene) + dxinv(2) * first_deriv_8(tmpy(i,j-4:j+4))
          end do
       end do

    end do
    !$omp end do

    !$omp do
    do k=lo(3)-4,hi(3)+4
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             tmpz(i,j,k) = mu(i,j,k)*qz(i,j,k,idu)
          end do
       end do
    end do
    !$omp end do
    !$omp do
    do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             rhs(i,j,k,imx) = rhs(i,j,k,imx) + dxinv(3) * first_deriv_8(tmpz(i,j,k-4:k+4))
          end do
       end do
    end do
    !$omp end do
    
    !$omp do
    do k=lo(3)-4,hi(3)+4
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             tmpz(i,j,k) = mu(i,j,k)*qz(i,j,k,idv)
          end do
       end do
    end do
    !$omp end do
    !$omp do
    do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             rhs(i,j,k,imy) = rhs(i,j,k,imy) + dxinv(3) * first_deriv_8(tmpz(i,j,k-4:k+4))
          end do
       end do
    end do
    !$omp end do
    
    !$omp do
    do k=lo(3)-4,hi(3)+4
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             tmpz(i,j,k) = vp(i,j,k)*qz(i,j,k,idw)
          end do
       end do
    end do
    !$omp end do
    !$omp do
    do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             rhs(i,j,k,imz) = rhs(i,j,k,imz) + dxinv(3) * first_deriv_8(tmpz(i,j,k-4:k+4))
          end do
       end do
    end do
    !$omp end do
    
    !$omp do
    do k=lo(3)-4,hi(3)+4
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             tmpz(i,j,k) = lam(i,j,k)*qz(i,j,k,idT)
          end do
       end do
    end do
    !$omp end do
    !$omp do
    do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             rhs(i,j,k,iene) = rhs(i,j,k,iene) + dxinv(3) * first_deriv_8(tmpz(i,j,k-4:k+4))
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

    ! x-direction
    !$omp do
    do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          do i=lo(1)-ng,hi(1)+ng
             FE(i,j,k) = dpe(i,j,k) * qx(i,j,k,idp)
          end do
       end do
    end do
    !$omp end do nowait

    !$omp workshare
    rvc = 0.d0
    !$omp end workshare

    do n=1,nspecies
       idXn = idX1+n-1
       qhn = qh1+n-1
       !$omp do
       do k=lo(3),hi(3)
          do j=lo(2),hi(2)
             do i=lo(1)-ng,hi(1)+ng
                FE(i,j,k) = FE(i,j,k) + dxy(i,j,k,n)*qx(i,j,k,idXn)*q(i,j,k,qhn)
                FY(i,j,k,n) = dxy(i,j,k,n)*qx(i,j,k,idXn) + dpy(i,j,k,n)*qx(i,j,k,idp)
                rvc(i,j,k) = rvc(i,j,k) + FY(i,j,k,n)
             end do
          end do
       end do
       !$omp end do
    end do

    do n=1,nspecies
       qyn = qy1+n-1
       qhn = qh1+n-1
       !$omp do
       do k=lo(3),hi(3)
          do j=lo(2),hi(2)
             do i=lo(1)-ng,hi(1)+ng
                FY(i,j,k,n) = FY(i,j,k,n) - rvc(i,j,k)*q(i,j,k,qyn)
                FE(i,j,k) = FE(i,j,k) - rvc(i,j,k)*q(i,j,k,qyn)*q(i,j,k,qhn)
             end do
          end do
       end do
       !$omp end do
    end do

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
             FE(i,j,k) = dpe(i,j,k) * qy(i,j,k,idp)
          end do
       end do
    end do
    !$omp end do nowait

    !$omp workshare
    rvc = 0.d0
    !$omp end workshare

    do n=1,nspecies
       idXn = idX1+n-1
       qhn = qh1+n-1
       !$omp do    
       do k=lo(3),hi(3)
          do j=lo(2)-ng,hi(2)+ng
             do i=lo(1),hi(1)
                FE(i,j,k) = FE(i,j,k) + dxy(i,j,k,n)*qy(i,j,k,idXn)*q(i,j,k,qhn)
                FY(i,j,k,n) = dxy(i,j,k,n)*qy(i,j,k,idXn) + dpy(i,j,k,n)*qy(i,j,k,idp)
                rvc(i,j,k) = rvc(i,j,k) + FY(i,j,k,n)                
             end do
          end do
       end do
       !$omp end do    
    end do

    do n=1,nspecies
       qyn = qy1+n-1
       qhn = qh1+n-1
       !$omp do    
       do k=lo(3),hi(3)
          do j=lo(2)-ng,hi(2)+ng
             do i=lo(1),hi(1)
                FY(i,j,k,n) = FY(i,j,k,n) - rvc(i,j,k)*q(i,j,k,qyn)
                FE(i,j,k) = FE(i,j,k) - rvc(i,j,k)*q(i,j,k,qyn)*q(i,j,k,qhn)                
             end do
          end do
       end do
       !$omp end do
    end do
    
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
             FE(i,j,k) = dpe(i,j,k) * qz(i,j,k,idp)
          end do
       end do
    end do
    !$omp end do

    !$omp workshare
    rvc = 0.d0
    !$omp end workshare

    do n=1,nspecies
       idXn = idX1+n-1
       qhn = qh1+n-1
       !$omp do
       do k=lo(3)-ng,hi(3)+ng
          do j=lo(2),hi(2)
             do i=lo(1),hi(1)
                FE(i,j,k) = FE(i,j,k) + dxy(i,j,k,n)*qz(i,j,k,idXn)*q(i,j,k,qhn)
                FY(i,j,k,n) = dxy(i,j,k,n)*qz(i,j,k,idXn) + dpy(i,j,k,n)*qz(i,j,k,idp)
                rvc(i,j,k) = rvc(i,j,k) + FY(i,j,k,n)                
             end do
          end do
       end do
       !$omp end do
    end do

    do n=1,nspecies
       qyn = qy1+n-1
       qhn = qh1+n-1
       !$omp do
       do k=lo(3)-ng,hi(3)+ng
          do j=lo(2),hi(2)
             do i=lo(1),hi(1)
                FY(i,j,k,n) = FY(i,j,k,n) - rvc(i,j,k)*q(i,j,k,qyn)
                FE(i,j,k) = FE(i,j,k) - rvc(i,j,k)*q(i,j,k,qyn)*q(i,j,k,qhn)
             end do
          end do
       end do
       !$omp end do
    end do

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
       !$omp end do
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

    deallocate(tmpx,tmpy)

    !$omp end parallel

    deallocate(vp,dpy,dpe,FY,FE,rvc,tmpz)

  end subroutine S3D_diffterm_2

end module kernels_module
