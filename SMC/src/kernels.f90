module kernels_module
  use chemistry_module, only : nspecies, molecular_weight
  use derivative_stencil_module
  use variables_module
  implicit none

contains

  subroutine hypterm_3d (lo,hi,ng,dx,cons,q,flux)

    integer,          intent(in ) :: lo(3),hi(3),ng
    double precision, intent(in ) :: dx(3)
    double precision, intent(in ) :: cons(-ng+lo(1):hi(1)+ng,-ng+lo(2):hi(2)+ng,-ng+lo(3):hi(3)+ng,ncons)
    double precision, intent(in ) ::    q(-ng+lo(1):hi(1)+ng,-ng+lo(2):hi(2)+ng,-ng+lo(3):hi(3)+ng,nprim)
    double precision, intent(out) :: flux(    lo(1):hi(1)   ,    lo(2):hi(2)   ,    lo(3):hi(3)   ,ncons)

    integer          :: i,j,k,n
    double precision :: unp1,unp2,unp3,unp4,unm1,unm2,unm3,unm4
    double precision :: dxinv(3)

    do i=1,3
       dxinv(i) = 1.0d0 / dx(i)
    end do

    !$omp parallel private(i,j,k,n,unp1,unp2,unp3,unp4,unm1,unm2,unm3,unm4)
    !$omp do 
    do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)

             unp1 = q(i+1,j,k,qu)
             unp2 = q(i+2,j,k,qu)
             unp3 = q(i+3,j,k,qu)
             unp4 = q(i+4,j,k,qu)

             unm1 = q(i-1,j,k,qu)
             unm2 = q(i-2,j,k,qu)
             unm3 = q(i-3,j,k,qu)
             unm4 = q(i-4,j,k,qu)

             flux(i,j,k,irho)= - &
                   (ALP*(cons(i+1,j,k,imx)-cons(i-1,j,k,imx)) &
                  + BET*(cons(i+2,j,k,imx)-cons(i-2,j,k,imx)) &
                  + GAM*(cons(i+3,j,k,imx)-cons(i-3,j,k,imx)) &
                  + DEL*(cons(i+4,j,k,imx)-cons(i-4,j,k,imx)))*dxinv(1)

             flux(i,j,k,imx)= - &
                   (ALP*(cons(i+1,j,k,imx)*unp1-cons(i-1,j,k,imx)*unm1 &
                  +        (q(i+1,j,k,qpres)   -   q(i-1,j,k,qpres)))  &
                  + BET*(cons(i+2,j,k,imx)*unp2-cons(i-2,j,k,imx)*unm2 &
                  +        (q(i+2,j,k,qpres)   -   q(i-2,j,k,qpres)))  &
                  + GAM*(cons(i+3,j,k,imx)*unp3-cons(i-3,j,k,imx)*unm3 &
                  +        (q(i+3,j,k,qpres)   -   q(i-3,j,k,qpres)))  &
                  + DEL*(cons(i+4,j,k,imx)*unp4-cons(i-4,j,k,imx)*unm4 &
                  +        (q(i+4,j,k,qpres)   -   q(i-4,j,k,qpres))))*dxinv(1)

             flux(i,j,k,imy)= - &
                   (ALP*(cons(i+1,j,k,imy)*unp1-cons(i-1,j,k,imy)*unm1) &
                  + BET*(cons(i+2,j,k,imy)*unp2-cons(i-2,j,k,imy)*unm2) &
                  + GAM*(cons(i+3,j,k,imy)*unp3-cons(i-3,j,k,imy)*unm3) &
                  + DEL*(cons(i+4,j,k,imy)*unp4-cons(i-4,j,k,imy)*unm4))*dxinv(1)

             flux(i,j,k,imz)= - &
                   (ALP*(cons(i+1,j,k,imz)*unp1-cons(i-1,j,k,imz)*unm1) &
                  + BET*(cons(i+2,j,k,imz)*unp2-cons(i-2,j,k,imz)*unm2) &
                  + GAM*(cons(i+3,j,k,imz)*unp3-cons(i-3,j,k,imz)*unm3) &
                  + DEL*(cons(i+4,j,k,imz)*unp4-cons(i-4,j,k,imz)*unm4))*dxinv(1)

             flux(i,j,k,iene)= - &
                   (ALP*(cons(i+1,j,k,iene )*unp1-cons(i-1,j,k,iene )*unm1 &
                  +        (q(i+1,j,k,qpres)*unp1-   q(i-1,j,k,qpres)*unm1)) &
                  + BET*(cons(i+2,j,k,iene )*unp2-cons(i-2,j,k,iene )*unm2 &
                  +        (q(i+2,j,k,qpres)*unp2-   q(i-2,j,k,qpres)*unm2)) &
                  + GAM*(cons(i+3,j,k,iene )*unp3-cons(i-3,j,k,iene )*unm3 &
                  +        (q(i+3,j,k,qpres)*unp3-   q(i-3,j,k,qpres)*unm3)) &
                  + DEL*(cons(i+4,j,k,iene )*unp4-cons(i-4,j,k,iene )*unm4 &
                  +        (q(i+4,j,k,qpres)*unp4-   q(i-4,j,k,qpres)*unm4)))*dxinv(1) 

             do n = iry1, iry1+nspecies-1
                flux(i,j,k,n) = - &
                     ( ALP*(cons(i+1,j,k,n)*unp1-cons(i-1,j,k,n)*unm1) &
                     + BET*(cons(i+2,j,k,n)*unp2-cons(i-2,j,k,n)*unm2) &
                     + GAM*(cons(i+3,j,k,n)*unp3-cons(i-3,j,k,n)*unm3) &
                     + DEL*(cons(i+4,j,k,n)*unp4-cons(i-4,j,k,n)*unm4))*dxinv(1)
             end do

          enddo
       enddo
    enddo
    !$omp end do

    !$omp do
    do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)

             unp1 = q(i,j+1,k,qv)
             unp2 = q(i,j+2,k,qv)
             unp3 = q(i,j+3,k,qv)
             unp4 = q(i,j+4,k,qv)

             unm1 = q(i,j-1,k,qv)
             unm2 = q(i,j-2,k,qv)
             unm3 = q(i,j-3,k,qv)
             unm4 = q(i,j-4,k,qv)

             flux(i,j,k,irho)=flux(i,j,k,irho) - &
                   (ALP*(cons(i,j+1,k,imy)-cons(i,j-1,k,imy)) &
                  + BET*(cons(i,j+2,k,imy)-cons(i,j-2,k,imy)) &
                  + GAM*(cons(i,j+3,k,imy)-cons(i,j-3,k,imy)) &
                  + DEL*(cons(i,j+4,k,imy)-cons(i,j-4,k,imy)))*dxinv(2)

             flux(i,j,k,imx)=flux(i,j,k,imx) - &
                   (ALP*(cons(i,j+1,k,imx)*unp1-cons(i,j-1,k,imx)*unm1) &
                  + BET*(cons(i,j+2,k,imx)*unp2-cons(i,j-2,k,imx)*unm2) &
                  + GAM*(cons(i,j+3,k,imx)*unp3-cons(i,j-3,k,imx)*unm3) &
                  + DEL*(cons(i,j+4,k,imx)*unp4-cons(i,j-4,k,imx)*unm4))*dxinv(2)

             flux(i,j,k,imy)=flux(i,j,k,imy) - &
                   (ALP*(cons(i,j+1,k,imy)*unp1-cons(i,j-1,k,imy)*unm1 &
                  +        (q(i,j+1,k,qpres)   -   q(i,j-1,k,qpres)))  &
                  + BET*(cons(i,j+2,k,imy)*unp2-cons(i,j-2,k,imy)*unm2 &
                  +        (q(i,j+2,k,qpres)   -   q(i,j-2,k,qpres)))  &
                  + GAM*(cons(i,j+3,k,imy)*unp3-cons(i,j-3,k,imy)*unm3 &
                  +        (q(i,j+3,k,qpres)   -   q(i,j-3,k,qpres)))  &
                  + DEL*(cons(i,j+4,k,imy)*unp4-cons(i,j-4,k,imy)*unm4 &
                  +        (q(i,j+4,k,qpres)   -   q(i,j-4,k,qpres))))*dxinv(2)

             flux(i,j,k,imz)=flux(i,j,k,imz) - &
                   (ALP*(cons(i,j+1,k,imz)*unp1-cons(i,j-1,k,imz)*unm1) &
                  + BET*(cons(i,j+2,k,imz)*unp2-cons(i,j-2,k,imz)*unm2) &
                  + GAM*(cons(i,j+3,k,imz)*unp3-cons(i,j-3,k,imz)*unm3) &
                  + DEL*(cons(i,j+4,k,imz)*unp4-cons(i,j-4,k,imz)*unm4))*dxinv(2)

             flux(i,j,k,iene)=flux(i,j,k,iene) - &
                   (ALP*(cons(i,j+1,k,iene )*unp1-cons(i,j-1,k,iene )*unm1 &
                  +        (q(i,j+1,k,qpres)*unp1-   q(i,j-1,k,qpres)*unm1)) &
                  + BET*(cons(i,j+2,k,iene )*unp2-cons(i,j-2,k,iene )*unm2 &
                  +        (q(i,j+2,k,qpres)*unp2-   q(i,j-2,k,qpres)*unm2)) &
                  + GAM*(cons(i,j+3,k,iene )*unp3-cons(i,j-3,k,iene )*unm3 &
                  +        (q(i,j+3,k,qpres)*unp3-   q(i,j-3,k,qpres)*unm3)) &
                  + DEL*(cons(i,j+4,k,iene )*unp4-cons(i,j-4,k,iene )*unm4 &
                  +        (q(i,j+4,k,qpres)*unp4-   q(i,j-4,k,qpres)*unm4)))*dxinv(2)

             do n = iry1, iry1+nspecies-1
                flux(i,j,k,n) = flux(i,j,k,n) - &
                     ( ALP*(cons(i,j+1,k,n)*unp1-cons(i,j-1,k,n)*unm1) &
                     + BET*(cons(i,j+2,k,n)*unp2-cons(i,j-2,k,n)*unm2) &
                     + GAM*(cons(i,j+3,k,n)*unp3-cons(i,j-3,k,n)*unm3) &
                     + DEL*(cons(i,j+4,k,n)*unp4-cons(i,j-4,k,n)*unm4))*dxinv(2)
             end do

          enddo
       enddo
    enddo
    !$omp end do

    !$omp do
    do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)

             unp1 = q(i,j,k+1,qw)
             unp2 = q(i,j,k+2,qw)
             unp3 = q(i,j,k+3,qw)
             unp4 = q(i,j,k+4,qw)

             unm1 = q(i,j,k-1,qw)
             unm2 = q(i,j,k-2,qw)
             unm3 = q(i,j,k-3,qw)
             unm4 = q(i,j,k-4,qw)

             flux(i,j,k,irho)=flux(i,j,k,irho) - &
                   (ALP*(cons(i,j,k+1,imz)-cons(i,j,k-1,imz)) &
                  + BET*(cons(i,j,k+2,imz)-cons(i,j,k-2,imz)) &
                  + GAM*(cons(i,j,k+3,imz)-cons(i,j,k-3,imz)) &
                  + DEL*(cons(i,j,k+4,imz)-cons(i,j,k-4,imz)))*dxinv(3)

             flux(i,j,k,imx)=flux(i,j,k,imx) - &
                   (ALP*(cons(i,j,k+1,imx)*unp1-cons(i,j,k-1,imx)*unm1) &
                  + BET*(cons(i,j,k+2,imx)*unp2-cons(i,j,k-2,imx)*unm2) &
                  + GAM*(cons(i,j,k+3,imx)*unp3-cons(i,j,k-3,imx)*unm3) &
                  + DEL*(cons(i,j,k+4,imx)*unp4-cons(i,j,k-4,imx)*unm4))*dxinv(3)

             flux(i,j,k,imy)=flux(i,j,k,imy) - &
                   (ALP*(cons(i,j,k+1,imy)*unp1-cons(i,j,k-1,imy)*unm1) &
                  + BET*(cons(i,j,k+2,imy)*unp2-cons(i,j,k-2,imy)*unm2) &
                  + GAM*(cons(i,j,k+3,imy)*unp3-cons(i,j,k-3,imy)*unm3) &
                  + DEL*(cons(i,j,k+4,imy)*unp4-cons(i,j,k-4,imy)*unm4))*dxinv(3)

             flux(i,j,k,imz)=flux(i,j,k,imz) - &
                   (ALP*(cons(i,j,k+1,imz)*unp1-cons(i,j,k-1,imz)*unm1 &
                  +        (q(i,j,k+1,qpres)   -   q(i,j,k-1,qpres)))  &
                  + BET*(cons(i,j,k+2,imz)*unp2-cons(i,j,k-2,imz)*unm2 &
                  +        (q(i,j,k+2,qpres)   -   q(i,j,k-2,qpres)))  &
                  + GAM*(cons(i,j,k+3,imz)*unp3-cons(i,j,k-3,imz)*unm3 &
                  +        (q(i,j,k+3,qpres)   -   q(i,j,k-3,qpres)))  &
                  + DEL*(cons(i,j,k+4,imz)*unp4-cons(i,j,k-4,imz)*unm4 &
                  +        (q(i,j,k+4,qpres)   -   q(i,j,k-4,qpres))))*dxinv(3)

             flux(i,j,k,iene)=flux(i,j,k,iene) - &
                   (ALP*(cons(i,j,k+1,iene )*unp1-cons(i,j,k-1,iene )*unm1 &
                  +        (q(i,j,k+1,qpres)*unp1-   q(i,j,k-1,qpres)*unm1)) &
                  + BET*(cons(i,j,k+2,iene )*unp2-cons(i,j,k-2,iene )*unm2 &
                  +        (q(i,j,k+2,qpres)*unp2-   q(i,j,k-2,qpres)*unm2)) &
                  + GAM*(cons(i,j,k+3,iene )*unp3-cons(i,j,k-3,iene )*unm3 &
                  +        (q(i,j,k+3,qpres)*unp3-   q(i,j,k-3,qpres)*unm3)) &
                  + DEL*(cons(i,j,k+4,iene )*unp4-cons(i,j,k-4,iene )*unm4 &
                  +        (q(i,j,k+4,qpres)*unp4-   q(i,j,k-4,qpres)*unm4)))*dxinv(3)

             do n = iry1, iry1+nspecies-1
                flux(i,j,k,n) = flux(i,j,k,n) - &
                     ( ALP*(cons(i,j,k+1,n)*unp1-cons(i,j,k-1,n)*unm1) &
                     + BET*(cons(i,j,k+2,n)*unp2-cons(i,j,k-2,n)*unm2) &
                     + GAM*(cons(i,j,k+3,n)*unp3-cons(i,j,k-3,n)*unm3) &
                     + DEL*(cons(i,j,k+4,n)*unp4-cons(i,j,k-4,n)*unm4))*dxinv(3)
             end do

          enddo
       enddo
    enddo
    !$omp end do
    !$omp end parallel

  end subroutine hypterm_3d


  subroutine compact_diffterm_3d (lo,hi,ng,dx,q,flx,mu,xi,lam,dxy)

    integer,          intent(in ) :: lo(3),hi(3),ng
    double precision, intent(in ) :: dx(3)
    double precision, intent(in ) :: q  (-ng+lo(1):hi(1)+ng,-ng+lo(2):hi(2)+ng,-ng+lo(3):hi(3)+ng,nprim)
    double precision, intent(in ) :: mu (-ng+lo(1):hi(1)+ng,-ng+lo(2):hi(2)+ng,-ng+lo(3):hi(3)+ng)
    double precision, intent(in ) :: xi (-ng+lo(1):hi(1)+ng,-ng+lo(2):hi(2)+ng,-ng+lo(3):hi(3)+ng)
    double precision, intent(in ) :: lam(-ng+lo(1):hi(1)+ng,-ng+lo(2):hi(2)+ng,-ng+lo(3):hi(3)+ng)
    double precision, intent(in ) :: dxy(-ng+lo(1):hi(1)+ng,-ng+lo(2):hi(2)+ng,-ng+lo(3):hi(3)+ng,nspecies)
    double precision, intent(out) :: flx(    lo(1):hi(1)   ,    lo(2):hi(2)   ,    lo(3):hi(3)   ,ncons)
 
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

    allocate(ux(    lo(1):hi(1)   ,-ng+lo(2):hi(2)+ng,-ng+lo(3):hi(3)+ng))
    allocate(vx(    lo(1):hi(1)   ,-ng+lo(2):hi(2)+ng,-ng+lo(3):hi(3)+ng))
    allocate(wx(    lo(1):hi(1)   ,-ng+lo(2):hi(2)+ng,-ng+lo(3):hi(3)+ng))
    allocate(uy(-ng+lo(1):hi(1)+ng,    lo(2):hi(2)   ,-ng+lo(3):hi(3)+ng))
    allocate(vy(-ng+lo(1):hi(1)+ng,    lo(2):hi(2)   ,-ng+lo(3):hi(3)+ng))
    allocate(wy(-ng+lo(1):hi(1)+ng,    lo(2):hi(2)   ,-ng+lo(3):hi(3)+ng))
    allocate(uz(-ng+lo(1):hi(1)+ng,-ng+lo(2):hi(2)+ng,    lo(3):hi(3)   ))
    allocate(vz(-ng+lo(1):hi(1)+ng,-ng+lo(2):hi(2)+ng,    lo(3):hi(3)   ))
    allocate(wz(-ng+lo(1):hi(1)+ng,-ng+lo(2):hi(2)+ng,    lo(3):hi(3)   ))

    allocate(vsp(-ng+lo(1):hi(1)+ng,-ng+lo(2):hi(2)+ng,-ng+lo(3):hi(3)+ng))
    allocate(vsm(-ng+lo(1):hi(1)+ng,-ng+lo(2):hi(2)+ng,-ng+lo(3):hi(3)+ng))

    allocate(Hg(lo(1):hi(1)+1,lo(2):hi(2)+1,lo(3):hi(3)+1,2:ncons))

    allocate(dpy(-ng+lo(1):hi(1)+ng,-ng+lo(2):hi(2)+ng,-ng+lo(3):hi(3)+ng,nspecies))
    allocate(dxe(-ng+lo(1):hi(1)+ng,-ng+lo(2):hi(2)+ng,-ng+lo(3):hi(3)+ng,nspecies))
    allocate(dpe(-ng+lo(1):hi(1)+ng,-ng+lo(2):hi(2)+ng,-ng+lo(3):hi(3)+ng))

    do i = 1,3
       dxinv(i) = 1.0d0 / dx(i)
       dx2inv(i) = dxinv(i)**2
    end do

    flx(:,:,:,irho) = 0.d0

    !$omp parallel private(i,j,k,n,qxn,qyn,qhn,Htot,Htmp,Ytmp,hhalf) &
    !$omp private(tauxx,tauyy,tauzz,dmuzdx,dmvzdy,dmuxvydz,dmuydx,dmwydz,dmuxwzdy) &
    !$omp private(dmvxdy,dmwxdz,dmvywzdx,divu)

    !$OMP DO
    do k=lo(3)-ng,hi(3)+ng
       do j=lo(2)-ng,hi(2)+ng
          do i=lo(1)-ng,hi(1)+ng
             vsp(i,j,k) = xi(i,j,k) + FourThirds*mu(i,j,k)
             vsm(i,j,k) = xi(i,j,k) -  TwoThirds*mu(i,j,k)
          enddo
       enddo
    enddo
    !$OMP END DO NOWAIT

    !$omp do
    do k=lo(3)-ng,hi(3)+ng
       do j=lo(2)-ng,hi(2)+ng
          do i=lo(1),hi(1)

             ux(i,j,k)= &
                   (ALP*(q(i+1,j,k,qu)-q(i-1,j,k,qu)) &
                  + BET*(q(i+2,j,k,qu)-q(i-2,j,k,qu)) &
                  + GAM*(q(i+3,j,k,qu)-q(i-3,j,k,qu)) &
                  + DEL*(q(i+4,j,k,qu)-q(i-4,j,k,qu)))*dxinv(1)

             vx(i,j,k)= &
                   (ALP*(q(i+1,j,k,qv)-q(i-1,j,k,qv)) &
                  + BET*(q(i+2,j,k,qv)-q(i-2,j,k,qv)) &
                  + GAM*(q(i+3,j,k,qv)-q(i-3,j,k,qv)) &
                  + DEL*(q(i+4,j,k,qv)-q(i-4,j,k,qv)))*dxinv(1)

             wx(i,j,k)= &
                   (ALP*(q(i+1,j,k,qw)-q(i-1,j,k,qw)) &
                  + BET*(q(i+2,j,k,qw)-q(i-2,j,k,qw)) &
                  + GAM*(q(i+3,j,k,qw)-q(i-3,j,k,qw)) &
                  + DEL*(q(i+4,j,k,qw)-q(i-4,j,k,qw)))*dxinv(1)
          enddo
       enddo
    enddo
    !$OMP END DO NOWAIT

    !$OMP DO
    do k=lo(3)-ng,hi(3)+ng
       do j=lo(2),hi(2)   
          do i=lo(1)-ng,hi(1)+ng

             uy(i,j,k)= &
                   (ALP*(q(i,j+1,k,qu)-q(i,j-1,k,qu)) &
                  + BET*(q(i,j+2,k,qu)-q(i,j-2,k,qu)) &
                  + GAM*(q(i,j+3,k,qu)-q(i,j-3,k,qu)) &
                  + DEL*(q(i,j+4,k,qu)-q(i,j-4,k,qu)))*dxinv(2)

             vy(i,j,k)= &
                   (ALP*(q(i,j+1,k,qv)-q(i,j-1,k,qv)) &
                  + BET*(q(i,j+2,k,qv)-q(i,j-2,k,qv)) &
                  + GAM*(q(i,j+3,k,qv)-q(i,j-3,k,qv)) &
                  + DEL*(q(i,j+4,k,qv)-q(i,j-4,k,qv)))*dxinv(2)

             wy(i,j,k)= &
                   (ALP*(q(i,j+1,k,qw)-q(i,j-1,k,qw)) &
                  + BET*(q(i,j+2,k,qw)-q(i,j-2,k,qw)) &
                  + GAM*(q(i,j+3,k,qw)-q(i,j-3,k,qw)) &
                  + DEL*(q(i,j+4,k,qw)-q(i,j-4,k,qw)))*dxinv(2)
          enddo
       enddo
    enddo
    !$OMP END DO NOWAIT

    !$OMP DO
    do k=lo(3),hi(3)
       do j=lo(2)-ng,hi(2)+ng
          do i=lo(1)-ng,hi(1)+ng

             uz(i,j,k)= &
                   (ALP*(q(i,j,k+1,qu)-q(i,j,k-1,qu)) &
                  + BET*(q(i,j,k+2,qu)-q(i,j,k-2,qu)) &
                  + GAM*(q(i,j,k+3,qu)-q(i,j,k-3,qu)) &
                  + DEL*(q(i,j,k+4,qu)-q(i,j,k-4,qu)))*dxinv(3)

             vz(i,j,k)= &
                   (ALP*(q(i,j,k+1,qv)-q(i,j,k-1,qv)) &
                  + BET*(q(i,j,k+2,qv)-q(i,j,k-2,qv)) &
                  + GAM*(q(i,j,k+3,qv)-q(i,j,k-3,qv)) &
                  + DEL*(q(i,j,k+4,qv)-q(i,j,k-4,qv)))*dxinv(3)

             wz(i,j,k)= &
                   (ALP*(q(i,j,k+1,qw)-q(i,j,k-1,qw)) &
                  + BET*(q(i,j,k+2,qw)-q(i,j,k-2,qw)) &
                  + GAM*(q(i,j,k+3,qw)-q(i,j,k-3,qw)) &
                  + DEL*(q(i,j,k+4,qw)-q(i,j,k-4,qw)))*dxinv(3)
          enddo
       enddo
    enddo
    !$OMP END DO NOWAIT

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
                dxe(i,j,k,n) = dxy(i,j,k,n)*q(i,j,k,qhn)
                dpe(i,j,k) = dpe(i,j,k) + dpy(i,j,k,n)*q(i,j,k,qhn)
             end do
          end do
       end do
       !$omp end do nowait
    end do

    !$omp do
    do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)

             ! d(mu*dv/dx)/dy
             dmvxdy = (ALP*(mu(i,j+1,k)*vx(i,j+1,k)-mu(i,j-1,k)*vx(i,j-1,k)) &
                  +    BET*(mu(i,j+2,k)*vx(i,j+2,k)-mu(i,j-2,k)*vx(i,j-2,k)) &
                  +    GAM*(mu(i,j+3,k)*vx(i,j+3,k)-mu(i,j-3,k)*vx(i,j-3,k)) &
                  +    DEL*(mu(i,j+4,k)*vx(i,j+4,k)-mu(i,j-4,k)*vx(i,j-4,k)))*dxinv(2) 

             ! d(mu*dw/dx)/dz
             dmwxdz = (ALP*(mu(i,j,k+1)*wx(i,j,k+1)-mu(i,j,k-1)*wx(i,j,k-1)) &
                  +    BET*(mu(i,j,k+2)*wx(i,j,k+2)-mu(i,j,k-2)*wx(i,j,k-2)) &
                  +    GAM*(mu(i,j,k+3)*wx(i,j,k+3)-mu(i,j,k-3)*wx(i,j,k-3)) &
                  +    DEL*(mu(i,j,k+4)*wx(i,j,k+4)-mu(i,j,k-4)*wx(i,j,k-4)))*dxinv(3) 

             ! d((xi-2/3*mu)*(vy+wz))/dx
             dmvywzdx = (ALP*(vsm(i+1,j,k)*(vy(i+1,j,k)+wz(i+1,j,k))-vsm(i-1,j,k)*(vy(i-1,j,k)+wz(i-1,j,k))) &
                  +      BET*(vsm(i+2,j,k)*(vy(i+2,j,k)+wz(i+2,j,k))-vsm(i-2,j,k)*(vy(i-2,j,k)+wz(i-2,j,k))) &
                  +      GAM*(vsm(i+3,j,k)*(vy(i+3,j,k)+wz(i+3,j,k))-vsm(i-3,j,k)*(vy(i-3,j,k)+wz(i-3,j,k))) &
                  +      DEL*(vsm(i+4,j,k)*(vy(i+4,j,k)+wz(i+4,j,k))-vsm(i-4,j,k)*(vy(i-4,j,k)+wz(i-4,j,k))) &
                  ) * dxinv(1)

             ! d(mu*du/dy)/dx
             dmuydx = (ALP*(mu(i+1,j,k)*uy(i+1,j,k)-mu(i-1,j,k)*uy(i-1,j,k)) &
                  +    BET*(mu(i+2,j,k)*uy(i+2,j,k)-mu(i-2,j,k)*uy(i-2,j,k)) &
                  +    GAM*(mu(i+3,j,k)*uy(i+3,j,k)-mu(i-3,j,k)*uy(i-3,j,k)) &
                  +    DEL*(mu(i+4,j,k)*uy(i+4,j,k)-mu(i-4,j,k)*uy(i-4,j,k)))*dxinv(1) 

             ! d(mu*dw/dy)/dz
             dmwydz = (ALP*(mu(i,j,k+1)*wy(i,j,k+1)-mu(i,j,k-1)*wy(i,j,k-1)) &
                  +    BET*(mu(i,j,k+2)*wy(i,j,k+2)-mu(i,j,k-2)*wy(i,j,k-2)) &
                  +    GAM*(mu(i,j,k+3)*wy(i,j,k+3)-mu(i,j,k-3)*wy(i,j,k-3)) &
                  +    DEL*(mu(i,j,k+4)*wy(i,j,k+4)-mu(i,j,k-4)*wy(i,j,k-4)))*dxinv(3) 

             ! d((xi-2/3*mu)*(ux+wz))/dy
             dmuxwzdy = (ALP*(vsm(i,j+1,k)*(ux(i,j+1,k)+wz(i,j+1,k))-vsm(i,j-1,k)*(ux(i,j-1,k)+wz(i,j-1,k))) &
                  +      BET*(vsm(i,j+2,k)*(ux(i,j+2,k)+wz(i,j+2,k))-vsm(i,j-2,k)*(ux(i,j-2,k)+wz(i,j-2,k))) &
                  +      GAM*(vsm(i,j+3,k)*(ux(i,j+3,k)+wz(i,j+3,k))-vsm(i,j-3,k)*(ux(i,j-3,k)+wz(i,j-3,k))) &
                  +      DEL*(vsm(i,j+4,k)*(ux(i,j+4,k)+wz(i,j+4,k))-vsm(i,j-4,k)*(ux(i,j-4,k)+wz(i,j-4,k))) &
                  ) * dxinv(2)

             ! d(mu*du/dz)/dx
             dmuzdx = (ALP*(mu(i+1,j,k)*uz(i+1,j,k)-mu(i-1,j,k)*uz(i-1,j,k)) &
                  +    BET*(mu(i+2,j,k)*uz(i+2,j,k)-mu(i-2,j,k)*uz(i-2,j,k)) &
                  +    GAM*(mu(i+3,j,k)*uz(i+3,j,k)-mu(i-3,j,k)*uz(i-3,j,k)) &
                  +    DEL*(mu(i+4,j,k)*uz(i+4,j,k)-mu(i-4,j,k)*uz(i-4,j,k)))*dxinv(1) 

             ! d(mu*dv/dz)/dy
             dmvzdy = (ALP*(mu(i,j+1,k)*vz(i,j+1,k)-mu(i,j-1,k)*vz(i,j-1,k)) &
                  +    BET*(mu(i,j+2,k)*vz(i,j+2,k)-mu(i,j-2,k)*vz(i,j-2,k)) &
                  +    GAM*(mu(i,j+3,k)*vz(i,j+3,k)-mu(i,j-3,k)*vz(i,j-3,k)) &
                  +    DEL*(mu(i,j+4,k)*vz(i,j+4,k)-mu(i,j-4,k)*vz(i,j-4,k)))*dxinv(2) 

             ! d((xi-2/3*mu)*(ux+vy))/dz
             dmuxvydz = (ALP*(vsm(i,j,k+1)*(ux(i,j,k+1)+vy(i,j,k+1))-vsm(i,j,k-1)*(ux(i,j,k-1)+vy(i,j,k-1))) &
                  +      BET*(vsm(i,j,k+2)*(ux(i,j,k+2)+vy(i,j,k+2))-vsm(i,j,k-2)*(ux(i,j,k-2)+vy(i,j,k-2))) &
                  +      GAM*(vsm(i,j,k+3)*(ux(i,j,k+3)+vy(i,j,k+3))-vsm(i,j,k-3)*(ux(i,j,k-3)+vy(i,j,k-3))) &
                  +      DEL*(vsm(i,j,k+4)*(ux(i,j,k+4)+vy(i,j,k+4))-vsm(i,j,k-4)*(ux(i,j,k-4)+vy(i,j,k-4))) &
                  ) * dxinv(3)

             flx(i,j,k,imx) = dmvxdy + dmwxdz + dmvywzdx
             flx(i,j,k,imy) = dmuydx + dmwydz + dmuxwzdy
             flx(i,j,k,imz) = dmuzdx + dmvzdy + dmuxvydz

             divu = (ux(i,j,k)+vy(i,j,k)+wz(i,j,k))*vsm(i,j,k)
             tauxx = 2.d0*mu(i,j,k)*ux(i,j,k) + divu
             tauyy = 2.d0*mu(i,j,k)*vy(i,j,k) + divu
             tauzz = 2.d0*mu(i,j,k)*wz(i,j,k) + divu
             
             ! change in internal energy
             flx(i,j,k,iene) = tauxx*ux(i,j,k) + tauyy*vy(i,j,k) + tauzz*wz(i,j,k) &
                  + mu(i,j,k)*((uy(i,j,k)+vx(i,j,k))**2 &
                  &          + (wx(i,j,k)+uz(i,j,k))**2 &
                  &          + (vz(i,j,k)+wy(i,j,k))**2 )

          end do
       end do
    end do
    !$omp end do 

    do n=1,nspecies
       !$OMP DO
       do k=lo(3),hi(3)
          do j=lo(2),hi(2)
             do i=lo(1),hi(1)
                flx(i,j,k,iry1+n-1) = 0.d0
             end do
          end do
       end do
       !$omp end do nowait
    end do

    ! ------- BEGIN x-direction -------
    !$omp do
    do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)+1
             Hg(i,j,k,imx) = m11*(vsp(i-4,j,k)*q(i-4,j,k,qu)-vsp(i+3,j,k)*q(i+3,j,k,qu)) &
                  +          m12*(vsp(i-4,j,k)*q(i-3,j,k,qu)-vsp(i+3,j,k)*q(i+2,j,k,qu)) &
                  +          m13*(vsp(i-4,j,k)*q(i-2,j,k,qu)-vsp(i+3,j,k)*q(i+1,j,k,qu)) &
                  +          m14*(vsp(i-4,j,k)*q(i-1,j,k,qu)-vsp(i+3,j,k)*q(i  ,j,k,qu)) &
                  +          m15*(vsp(i-4,j,k)*q(i  ,j,k,qu)-vsp(i+3,j,k)*q(i-1,j,k,qu)) &
                  &        + m21*(vsp(i-3,j,k)*q(i-4,j,k,qu)-vsp(i+2,j,k)*q(i+3,j,k,qu)) &
                  +          m22*(vsp(i-3,j,k)*q(i-3,j,k,qu)-vsp(i+2,j,k)*q(i+2,j,k,qu)) &
                  +          m23*(vsp(i-3,j,k)*q(i-2,j,k,qu)-vsp(i+2,j,k)*q(i+1,j,k,qu)) &
                  +          m24*(vsp(i-3,j,k)*q(i-1,j,k,qu)-vsp(i+2,j,k)*q(i  ,j,k,qu)) &
                  +          m25*(vsp(i-3,j,k)*q(i  ,j,k,qu)-vsp(i+2,j,k)*q(i-1,j,k,qu)) &
                  +          m26*(vsp(i-3,j,k)*q(i+1,j,k,qu)-vsp(i+2,j,k)*q(i-2,j,k,qu)) &
                  &        + m31*(vsp(i-2,j,k)*q(i-4,j,k,qu)-vsp(i+1,j,k)*q(i+3,j,k,qu)) &
                  +          m32*(vsp(i-2,j,k)*q(i-3,j,k,qu)-vsp(i+1,j,k)*q(i+2,j,k,qu)) &
                  +          m33*(vsp(i-2,j,k)*q(i-2,j,k,qu)-vsp(i+1,j,k)*q(i+1,j,k,qu)) &
                  +          m34*(vsp(i-2,j,k)*q(i-1,j,k,qu)-vsp(i+1,j,k)*q(i  ,j,k,qu)) &
                  +          m35*(vsp(i-2,j,k)*q(i  ,j,k,qu)-vsp(i+1,j,k)*q(i-1,j,k,qu)) &
                  +          m36*(vsp(i-2,j,k)*q(i+1,j,k,qu)-vsp(i+1,j,k)*q(i-2,j,k,qu)) &
                  +          m37*(vsp(i-2,j,k)*q(i+2,j,k,qu)-vsp(i+1,j,k)*q(i-3,j,k,qu)) &
                  &        + m41*(vsp(i-1,j,k)*q(i-4,j,k,qu)-vsp(i  ,j,k)*q(i+3,j,k,qu)) &
                  +          m42*(vsp(i-1,j,k)*q(i-3,j,k,qu)-vsp(i  ,j,k)*q(i+2,j,k,qu)) &
                  +          m43*(vsp(i-1,j,k)*q(i-2,j,k,qu)-vsp(i  ,j,k)*q(i+1,j,k,qu)) &
                  +          m44*(vsp(i-1,j,k)*q(i-1,j,k,qu)-vsp(i  ,j,k)*q(i  ,j,k,qu)) &
                  +          m45*(vsp(i-1,j,k)*q(i  ,j,k,qu)-vsp(i  ,j,k)*q(i-1,j,k,qu)) &
                  +          m46*(vsp(i-1,j,k)*q(i+1,j,k,qu)-vsp(i  ,j,k)*q(i-2,j,k,qu)) &
                  +          m47*(vsp(i-1,j,k)*q(i+2,j,k,qu)-vsp(i  ,j,k)*q(i-3,j,k,qu)) &
                  +          m48*(vsp(i-1,j,k)*q(i+3,j,k,qu)-vsp(i  ,j,k)*q(i-4,j,k,qu))
             
             Hg(i,j,k,imy) = m11*(mu(i-4,j,k)*q(i-4,j,k,qv)-mu(i+3,j,k)*q(i+3,j,k,qv)) &
                  +          m12*(mu(i-4,j,k)*q(i-3,j,k,qv)-mu(i+3,j,k)*q(i+2,j,k,qv)) &
                  +          m13*(mu(i-4,j,k)*q(i-2,j,k,qv)-mu(i+3,j,k)*q(i+1,j,k,qv)) &
                  +          m14*(mu(i-4,j,k)*q(i-1,j,k,qv)-mu(i+3,j,k)*q(i  ,j,k,qv)) &
                  +          m15*(mu(i-4,j,k)*q(i  ,j,k,qv)-mu(i+3,j,k)*q(i-1,j,k,qv)) &
                  &        + m21*(mu(i-3,j,k)*q(i-4,j,k,qv)-mu(i+2,j,k)*q(i+3,j,k,qv)) &
                  +          m22*(mu(i-3,j,k)*q(i-3,j,k,qv)-mu(i+2,j,k)*q(i+2,j,k,qv)) &
                  +          m23*(mu(i-3,j,k)*q(i-2,j,k,qv)-mu(i+2,j,k)*q(i+1,j,k,qv)) &
                  +          m24*(mu(i-3,j,k)*q(i-1,j,k,qv)-mu(i+2,j,k)*q(i  ,j,k,qv)) &
                  +          m25*(mu(i-3,j,k)*q(i  ,j,k,qv)-mu(i+2,j,k)*q(i-1,j,k,qv)) &
                  +          m26*(mu(i-3,j,k)*q(i+1,j,k,qv)-mu(i+2,j,k)*q(i-2,j,k,qv)) &
                  &        + m31*(mu(i-2,j,k)*q(i-4,j,k,qv)-mu(i+1,j,k)*q(i+3,j,k,qv)) &
                  +          m32*(mu(i-2,j,k)*q(i-3,j,k,qv)-mu(i+1,j,k)*q(i+2,j,k,qv)) &
                  +          m33*(mu(i-2,j,k)*q(i-2,j,k,qv)-mu(i+1,j,k)*q(i+1,j,k,qv)) &
                  +          m34*(mu(i-2,j,k)*q(i-1,j,k,qv)-mu(i+1,j,k)*q(i  ,j,k,qv)) &
                  +          m35*(mu(i-2,j,k)*q(i  ,j,k,qv)-mu(i+1,j,k)*q(i-1,j,k,qv)) &
                  +          m36*(mu(i-2,j,k)*q(i+1,j,k,qv)-mu(i+1,j,k)*q(i-2,j,k,qv)) &
                  +          m37*(mu(i-2,j,k)*q(i+2,j,k,qv)-mu(i+1,j,k)*q(i-3,j,k,qv)) &
                  &        + m41*(mu(i-1,j,k)*q(i-4,j,k,qv)-mu(i  ,j,k)*q(i+3,j,k,qv)) &
                  +          m42*(mu(i-1,j,k)*q(i-3,j,k,qv)-mu(i  ,j,k)*q(i+2,j,k,qv)) &
                  +          m43*(mu(i-1,j,k)*q(i-2,j,k,qv)-mu(i  ,j,k)*q(i+1,j,k,qv)) &
                  +          m44*(mu(i-1,j,k)*q(i-1,j,k,qv)-mu(i  ,j,k)*q(i  ,j,k,qv)) &
                  +          m45*(mu(i-1,j,k)*q(i  ,j,k,qv)-mu(i  ,j,k)*q(i-1,j,k,qv)) &
                  +          m46*(mu(i-1,j,k)*q(i+1,j,k,qv)-mu(i  ,j,k)*q(i-2,j,k,qv)) &
                  +          m47*(mu(i-1,j,k)*q(i+2,j,k,qv)-mu(i  ,j,k)*q(i-3,j,k,qv)) &
                  +          m48*(mu(i-1,j,k)*q(i+3,j,k,qv)-mu(i  ,j,k)*q(i-4,j,k,qv))

             Hg(i,j,k,imz) = m11*(mu(i-4,j,k)*q(i-4,j,k,qw)-mu(i+3,j,k)*q(i+3,j,k,qw)) &
                  +          m12*(mu(i-4,j,k)*q(i-3,j,k,qw)-mu(i+3,j,k)*q(i+2,j,k,qw)) &
                  +          m13*(mu(i-4,j,k)*q(i-2,j,k,qw)-mu(i+3,j,k)*q(i+1,j,k,qw)) &
                  +          m14*(mu(i-4,j,k)*q(i-1,j,k,qw)-mu(i+3,j,k)*q(i  ,j,k,qw)) &
                  +          m15*(mu(i-4,j,k)*q(i  ,j,k,qw)-mu(i+3,j,k)*q(i-1,j,k,qw)) &
                  &        + m21*(mu(i-3,j,k)*q(i-4,j,k,qw)-mu(i+2,j,k)*q(i+3,j,k,qw)) &
                  +          m22*(mu(i-3,j,k)*q(i-3,j,k,qw)-mu(i+2,j,k)*q(i+2,j,k,qw)) &
                  +          m23*(mu(i-3,j,k)*q(i-2,j,k,qw)-mu(i+2,j,k)*q(i+1,j,k,qw)) &
                  +          m24*(mu(i-3,j,k)*q(i-1,j,k,qw)-mu(i+2,j,k)*q(i  ,j,k,qw)) &
                  +          m25*(mu(i-3,j,k)*q(i  ,j,k,qw)-mu(i+2,j,k)*q(i-1,j,k,qw)) &
                  +          m26*(mu(i-3,j,k)*q(i+1,j,k,qw)-mu(i+2,j,k)*q(i-2,j,k,qw)) &
                  &        + m31*(mu(i-2,j,k)*q(i-4,j,k,qw)-mu(i+1,j,k)*q(i+3,j,k,qw)) &
                  +          m32*(mu(i-2,j,k)*q(i-3,j,k,qw)-mu(i+1,j,k)*q(i+2,j,k,qw)) &
                  +          m33*(mu(i-2,j,k)*q(i-2,j,k,qw)-mu(i+1,j,k)*q(i+1,j,k,qw)) &
                  +          m34*(mu(i-2,j,k)*q(i-1,j,k,qw)-mu(i+1,j,k)*q(i  ,j,k,qw)) &
                  +          m35*(mu(i-2,j,k)*q(i  ,j,k,qw)-mu(i+1,j,k)*q(i-1,j,k,qw)) &
                  +          m36*(mu(i-2,j,k)*q(i+1,j,k,qw)-mu(i+1,j,k)*q(i-2,j,k,qw)) &
                  +          m37*(mu(i-2,j,k)*q(i+2,j,k,qw)-mu(i+1,j,k)*q(i-3,j,k,qw)) &
                  &        + m41*(mu(i-1,j,k)*q(i-4,j,k,qw)-mu(i  ,j,k)*q(i+3,j,k,qw)) &
                  +          m42*(mu(i-1,j,k)*q(i-3,j,k,qw)-mu(i  ,j,k)*q(i+2,j,k,qw)) &
                  +          m43*(mu(i-1,j,k)*q(i-2,j,k,qw)-mu(i  ,j,k)*q(i+1,j,k,qw)) &
                  +          m44*(mu(i-1,j,k)*q(i-1,j,k,qw)-mu(i  ,j,k)*q(i  ,j,k,qw)) &
                  +          m45*(mu(i-1,j,k)*q(i  ,j,k,qw)-mu(i  ,j,k)*q(i-1,j,k,qw)) &
                  +          m46*(mu(i-1,j,k)*q(i+1,j,k,qw)-mu(i  ,j,k)*q(i-2,j,k,qw)) &
                  +          m47*(mu(i-1,j,k)*q(i+2,j,k,qw)-mu(i  ,j,k)*q(i-3,j,k,qw)) &
                  +          m48*(mu(i-1,j,k)*q(i+3,j,k,qw)-mu(i  ,j,k)*q(i-4,j,k,qw))

             Hg(i,j,k,iene) = m11*(lam(i-4,j,k)*q(i-4,j,k,qtemp)-lam(i+3,j,k)*q(i+3,j,k,qtemp)) &
                  +           m12*(lam(i-4,j,k)*q(i-3,j,k,qtemp)-lam(i+3,j,k)*q(i+2,j,k,qtemp)) &
                  +           m13*(lam(i-4,j,k)*q(i-2,j,k,qtemp)-lam(i+3,j,k)*q(i+1,j,k,qtemp)) &
                  +           m14*(lam(i-4,j,k)*q(i-1,j,k,qtemp)-lam(i+3,j,k)*q(i  ,j,k,qtemp)) &
                  +           m15*(lam(i-4,j,k)*q(i  ,j,k,qtemp)-lam(i+3,j,k)*q(i-1,j,k,qtemp)) &
                  &         + m21*(lam(i-3,j,k)*q(i-4,j,k,qtemp)-lam(i+2,j,k)*q(i+3,j,k,qtemp)) &
                  +           m22*(lam(i-3,j,k)*q(i-3,j,k,qtemp)-lam(i+2,j,k)*q(i+2,j,k,qtemp)) &
                  +           m23*(lam(i-3,j,k)*q(i-2,j,k,qtemp)-lam(i+2,j,k)*q(i+1,j,k,qtemp)) &
                  +           m24*(lam(i-3,j,k)*q(i-1,j,k,qtemp)-lam(i+2,j,k)*q(i  ,j,k,qtemp)) &
                  +           m25*(lam(i-3,j,k)*q(i  ,j,k,qtemp)-lam(i+2,j,k)*q(i-1,j,k,qtemp)) &
                  +           m26*(lam(i-3,j,k)*q(i+1,j,k,qtemp)-lam(i+2,j,k)*q(i-2,j,k,qtemp)) &
                  &         + m31*(lam(i-2,j,k)*q(i-4,j,k,qtemp)-lam(i+1,j,k)*q(i+3,j,k,qtemp)) &
                  +           m32*(lam(i-2,j,k)*q(i-3,j,k,qtemp)-lam(i+1,j,k)*q(i+2,j,k,qtemp)) &
                  +           m33*(lam(i-2,j,k)*q(i-2,j,k,qtemp)-lam(i+1,j,k)*q(i+1,j,k,qtemp)) &
                  +           m34*(lam(i-2,j,k)*q(i-1,j,k,qtemp)-lam(i+1,j,k)*q(i  ,j,k,qtemp)) &
                  +           m35*(lam(i-2,j,k)*q(i  ,j,k,qtemp)-lam(i+1,j,k)*q(i-1,j,k,qtemp)) &
                  +           m36*(lam(i-2,j,k)*q(i+1,j,k,qtemp)-lam(i+1,j,k)*q(i-2,j,k,qtemp)) &
                  +           m37*(lam(i-2,j,k)*q(i+2,j,k,qtemp)-lam(i+1,j,k)*q(i-3,j,k,qtemp)) &
                  &         + m41*(lam(i-1,j,k)*q(i-4,j,k,qtemp)-lam(i  ,j,k)*q(i+3,j,k,qtemp)) &
                  +           m42*(lam(i-1,j,k)*q(i-3,j,k,qtemp)-lam(i  ,j,k)*q(i+2,j,k,qtemp)) &
                  +           m43*(lam(i-1,j,k)*q(i-2,j,k,qtemp)-lam(i  ,j,k)*q(i+1,j,k,qtemp)) &
                  +           m44*(lam(i-1,j,k)*q(i-1,j,k,qtemp)-lam(i  ,j,k)*q(i  ,j,k,qtemp)) &
                  +           m45*(lam(i-1,j,k)*q(i  ,j,k,qtemp)-lam(i  ,j,k)*q(i-1,j,k,qtemp)) &
                  +           m46*(lam(i-1,j,k)*q(i+1,j,k,qtemp)-lam(i  ,j,k)*q(i-2,j,k,qtemp)) &
                  +           m47*(lam(i-1,j,k)*q(i+2,j,k,qtemp)-lam(i  ,j,k)*q(i-3,j,k,qtemp)) &
                  +           m48*(lam(i-1,j,k)*q(i+3,j,k,qtemp)-lam(i  ,j,k)*q(i-4,j,k,qtemp))

             Htot = 0.d0
             do n = 1, nspecies
                qxn = qx1+n-1
                qyn = qy1+n-1
                Htmp(n) = m11*(dxy(i-4,j,k,n)*q(i-4,j,k,qxn)-dxy(i+3,j,k,n)*q(i+3,j,k,qxn)) &
                  +       m12*(dxy(i-4,j,k,n)*q(i-3,j,k,qxn)-dxy(i+3,j,k,n)*q(i+2,j,k,qxn)) &
                  +       m13*(dxy(i-4,j,k,n)*q(i-2,j,k,qxn)-dxy(i+3,j,k,n)*q(i+1,j,k,qxn)) &
                  +       m14*(dxy(i-4,j,k,n)*q(i-1,j,k,qxn)-dxy(i+3,j,k,n)*q(i  ,j,k,qxn)) &
                  +       m15*(dxy(i-4,j,k,n)*q(i  ,j,k,qxn)-dxy(i+3,j,k,n)*q(i-1,j,k,qxn)) &
                  &     + m21*(dxy(i-3,j,k,n)*q(i-4,j,k,qxn)-dxy(i+2,j,k,n)*q(i+3,j,k,qxn)) &
                  +       m22*(dxy(i-3,j,k,n)*q(i-3,j,k,qxn)-dxy(i+2,j,k,n)*q(i+2,j,k,qxn)) &
                  +       m23*(dxy(i-3,j,k,n)*q(i-2,j,k,qxn)-dxy(i+2,j,k,n)*q(i+1,j,k,qxn)) &
                  +       m24*(dxy(i-3,j,k,n)*q(i-1,j,k,qxn)-dxy(i+2,j,k,n)*q(i  ,j,k,qxn)) &
                  +       m25*(dxy(i-3,j,k,n)*q(i  ,j,k,qxn)-dxy(i+2,j,k,n)*q(i-1,j,k,qxn)) &
                  +       m26*(dxy(i-3,j,k,n)*q(i+1,j,k,qxn)-dxy(i+2,j,k,n)*q(i-2,j,k,qxn)) &
                  &     + m31*(dxy(i-2,j,k,n)*q(i-4,j,k,qxn)-dxy(i+1,j,k,n)*q(i+3,j,k,qxn)) &
                  +       m32*(dxy(i-2,j,k,n)*q(i-3,j,k,qxn)-dxy(i+1,j,k,n)*q(i+2,j,k,qxn)) &
                  +       m33*(dxy(i-2,j,k,n)*q(i-2,j,k,qxn)-dxy(i+1,j,k,n)*q(i+1,j,k,qxn)) &
                  +       m34*(dxy(i-2,j,k,n)*q(i-1,j,k,qxn)-dxy(i+1,j,k,n)*q(i  ,j,k,qxn)) &
                  +       m35*(dxy(i-2,j,k,n)*q(i  ,j,k,qxn)-dxy(i+1,j,k,n)*q(i-1,j,k,qxn)) &
                  +       m36*(dxy(i-2,j,k,n)*q(i+1,j,k,qxn)-dxy(i+1,j,k,n)*q(i-2,j,k,qxn)) &
                  +       m37*(dxy(i-2,j,k,n)*q(i+2,j,k,qxn)-dxy(i+1,j,k,n)*q(i-3,j,k,qxn)) &
                  &     + m41*(dxy(i-1,j,k,n)*q(i-4,j,k,qxn)-dxy(i  ,j,k,n)*q(i+3,j,k,qxn)) &
                  +       m42*(dxy(i-1,j,k,n)*q(i-3,j,k,qxn)-dxy(i  ,j,k,n)*q(i+2,j,k,qxn)) &
                  +       m43*(dxy(i-1,j,k,n)*q(i-2,j,k,qxn)-dxy(i  ,j,k,n)*q(i+1,j,k,qxn)) &
                  +       m44*(dxy(i-1,j,k,n)*q(i-1,j,k,qxn)-dxy(i  ,j,k,n)*q(i  ,j,k,qxn)) &
                  +       m45*(dxy(i-1,j,k,n)*q(i  ,j,k,qxn)-dxy(i  ,j,k,n)*q(i-1,j,k,qxn)) &
                  +       m46*(dxy(i-1,j,k,n)*q(i+1,j,k,qxn)-dxy(i  ,j,k,n)*q(i-2,j,k,qxn)) &
                  +       m47*(dxy(i-1,j,k,n)*q(i+2,j,k,qxn)-dxy(i  ,j,k,n)*q(i-3,j,k,qxn)) &
                  +       m48*(dxy(i-1,j,k,n)*q(i+3,j,k,qxn)-dxy(i  ,j,k,n)*q(i-4,j,k,qxn))
                Htmp(n) = Htmp(n)  &
                  +       m11*(dpy(i-4,j,k,n)*q(i-4,j,k,qpres)-dpy(i+3,j,k,n)*q(i+3,j,k,qpres)) &
                  +       m12*(dpy(i-4,j,k,n)*q(i-3,j,k,qpres)-dpy(i+3,j,k,n)*q(i+2,j,k,qpres)) &
                  +       m13*(dpy(i-4,j,k,n)*q(i-2,j,k,qpres)-dpy(i+3,j,k,n)*q(i+1,j,k,qpres)) &
                  +       m14*(dpy(i-4,j,k,n)*q(i-1,j,k,qpres)-dpy(i+3,j,k,n)*q(i  ,j,k,qpres)) &
                  +       m15*(dpy(i-4,j,k,n)*q(i  ,j,k,qpres)-dpy(i+3,j,k,n)*q(i-1,j,k,qpres)) &
                  &     + m21*(dpy(i-3,j,k,n)*q(i-4,j,k,qpres)-dpy(i+2,j,k,n)*q(i+3,j,k,qpres)) &
                  +       m22*(dpy(i-3,j,k,n)*q(i-3,j,k,qpres)-dpy(i+2,j,k,n)*q(i+2,j,k,qpres)) &
                  +       m23*(dpy(i-3,j,k,n)*q(i-2,j,k,qpres)-dpy(i+2,j,k,n)*q(i+1,j,k,qpres)) &
                  +       m24*(dpy(i-3,j,k,n)*q(i-1,j,k,qpres)-dpy(i+2,j,k,n)*q(i  ,j,k,qpres)) &
                  +       m25*(dpy(i-3,j,k,n)*q(i  ,j,k,qpres)-dpy(i+2,j,k,n)*q(i-1,j,k,qpres)) &
                  +       m26*(dpy(i-3,j,k,n)*q(i+1,j,k,qpres)-dpy(i+2,j,k,n)*q(i-2,j,k,qpres)) &
                  &     + m31*(dpy(i-2,j,k,n)*q(i-4,j,k,qpres)-dpy(i+1,j,k,n)*q(i+3,j,k,qpres)) &
                  +       m32*(dpy(i-2,j,k,n)*q(i-3,j,k,qpres)-dpy(i+1,j,k,n)*q(i+2,j,k,qpres)) &
                  +       m33*(dpy(i-2,j,k,n)*q(i-2,j,k,qpres)-dpy(i+1,j,k,n)*q(i+1,j,k,qpres)) &
                  +       m34*(dpy(i-2,j,k,n)*q(i-1,j,k,qpres)-dpy(i+1,j,k,n)*q(i  ,j,k,qpres)) &
                  +       m35*(dpy(i-2,j,k,n)*q(i  ,j,k,qpres)-dpy(i+1,j,k,n)*q(i-1,j,k,qpres)) &
                  +       m36*(dpy(i-2,j,k,n)*q(i+1,j,k,qpres)-dpy(i+1,j,k,n)*q(i-2,j,k,qpres)) &
                  +       m37*(dpy(i-2,j,k,n)*q(i+2,j,k,qpres)-dpy(i+1,j,k,n)*q(i-3,j,k,qpres)) &
                  &     + m41*(dpy(i-1,j,k,n)*q(i-4,j,k,qpres)-dpy(i  ,j,k,n)*q(i+3,j,k,qpres)) &
                  +       m42*(dpy(i-1,j,k,n)*q(i-3,j,k,qpres)-dpy(i  ,j,k,n)*q(i+2,j,k,qpres)) &
                  +       m43*(dpy(i-1,j,k,n)*q(i-2,j,k,qpres)-dpy(i  ,j,k,n)*q(i+1,j,k,qpres)) &
                  +       m44*(dpy(i-1,j,k,n)*q(i-1,j,k,qpres)-dpy(i  ,j,k,n)*q(i  ,j,k,qpres)) &
                  +       m45*(dpy(i-1,j,k,n)*q(i  ,j,k,qpres)-dpy(i  ,j,k,n)*q(i-1,j,k,qpres)) &
                  +       m46*(dpy(i-1,j,k,n)*q(i+1,j,k,qpres)-dpy(i  ,j,k,n)*q(i-2,j,k,qpres)) &
                  +       m47*(dpy(i-1,j,k,n)*q(i+2,j,k,qpres)-dpy(i  ,j,k,n)*q(i-3,j,k,qpres)) &
                  +       m48*(dpy(i-1,j,k,n)*q(i+3,j,k,qpres)-dpy(i  ,j,k,n)*q(i-4,j,k,qpres))
                Htot = Htot + Htmp(n)
                Ytmp(n) = (q(i-1,j,k,qyn) + q(i,j,k,qyn)) / 2.d0
             end do

             do n = 1, nspecies
                Hg(i,j,k,iry1+n-1) = Htmp(n) - Ytmp(n)*Htot
             end do

             do n = 1, nspecies
                qxn = qx1+n-1
                Hg(i,j,k,iene) =  Hg(i,j,k,iene) &
                  + m11*(dxe(i-4,j,k,n)*q(i-4,j,k,qxn)-dxe(i+3,j,k,n)*q(i+3,j,k,qxn)) &
                  + m12*(dxe(i-4,j,k,n)*q(i-3,j,k,qxn)-dxe(i+3,j,k,n)*q(i+2,j,k,qxn)) &
                  + m13*(dxe(i-4,j,k,n)*q(i-2,j,k,qxn)-dxe(i+3,j,k,n)*q(i+1,j,k,qxn)) &
                  + m14*(dxe(i-4,j,k,n)*q(i-1,j,k,qxn)-dxe(i+3,j,k,n)*q(i  ,j,k,qxn)) &
                  + m15*(dxe(i-4,j,k,n)*q(i  ,j,k,qxn)-dxe(i+3,j,k,n)*q(i-1,j,k,qxn)) &
                  + m21*(dxe(i-3,j,k,n)*q(i-4,j,k,qxn)-dxe(i+2,j,k,n)*q(i+3,j,k,qxn)) &
                  + m22*(dxe(i-3,j,k,n)*q(i-3,j,k,qxn)-dxe(i+2,j,k,n)*q(i+2,j,k,qxn)) &
                  + m23*(dxe(i-3,j,k,n)*q(i-2,j,k,qxn)-dxe(i+2,j,k,n)*q(i+1,j,k,qxn)) &
                  + m24*(dxe(i-3,j,k,n)*q(i-1,j,k,qxn)-dxe(i+2,j,k,n)*q(i  ,j,k,qxn)) &
                  + m25*(dxe(i-3,j,k,n)*q(i  ,j,k,qxn)-dxe(i+2,j,k,n)*q(i-1,j,k,qxn)) &
                  + m26*(dxe(i-3,j,k,n)*q(i+1,j,k,qxn)-dxe(i+2,j,k,n)*q(i-2,j,k,qxn)) &
                  + m31*(dxe(i-2,j,k,n)*q(i-4,j,k,qxn)-dxe(i+1,j,k,n)*q(i+3,j,k,qxn)) &
                  + m32*(dxe(i-2,j,k,n)*q(i-3,j,k,qxn)-dxe(i+1,j,k,n)*q(i+2,j,k,qxn)) &
                  + m33*(dxe(i-2,j,k,n)*q(i-2,j,k,qxn)-dxe(i+1,j,k,n)*q(i+1,j,k,qxn)) &
                  + m34*(dxe(i-2,j,k,n)*q(i-1,j,k,qxn)-dxe(i+1,j,k,n)*q(i  ,j,k,qxn)) &
                  + m35*(dxe(i-2,j,k,n)*q(i  ,j,k,qxn)-dxe(i+1,j,k,n)*q(i-1,j,k,qxn)) &
                  + m36*(dxe(i-2,j,k,n)*q(i+1,j,k,qxn)-dxe(i+1,j,k,n)*q(i-2,j,k,qxn)) &
                  + m37*(dxe(i-2,j,k,n)*q(i+2,j,k,qxn)-dxe(i+1,j,k,n)*q(i-3,j,k,qxn)) &
                  + m41*(dxe(i-1,j,k,n)*q(i-4,j,k,qxn)-dxe(i  ,j,k,n)*q(i+3,j,k,qxn)) &
                  + m42*(dxe(i-1,j,k,n)*q(i-3,j,k,qxn)-dxe(i  ,j,k,n)*q(i+2,j,k,qxn)) &
                  + m43*(dxe(i-1,j,k,n)*q(i-2,j,k,qxn)-dxe(i  ,j,k,n)*q(i+1,j,k,qxn)) &
                  + m44*(dxe(i-1,j,k,n)*q(i-1,j,k,qxn)-dxe(i  ,j,k,n)*q(i  ,j,k,qxn)) &
                  + m45*(dxe(i-1,j,k,n)*q(i  ,j,k,qxn)-dxe(i  ,j,k,n)*q(i-1,j,k,qxn)) &
                  + m46*(dxe(i-1,j,k,n)*q(i+1,j,k,qxn)-dxe(i  ,j,k,n)*q(i-2,j,k,qxn)) &
                  + m47*(dxe(i-1,j,k,n)*q(i+2,j,k,qxn)-dxe(i  ,j,k,n)*q(i-3,j,k,qxn)) &
                  + m48*(dxe(i-1,j,k,n)*q(i+3,j,k,qxn)-dxe(i  ,j,k,n)*q(i-4,j,k,qxn))
             end do

             Hg(i,j,k,iene) =  Hg(i,j,k,iene) &
                  + m11*(dpe(i-4,j,k)*q(i-4,j,k,qpres)-dpe(i+3,j,k)*q(i+3,j,k,qpres)) &
                  + m12*(dpe(i-4,j,k)*q(i-3,j,k,qpres)-dpe(i+3,j,k)*q(i+2,j,k,qpres)) &
                  + m13*(dpe(i-4,j,k)*q(i-2,j,k,qpres)-dpe(i+3,j,k)*q(i+1,j,k,qpres)) &
                  + m14*(dpe(i-4,j,k)*q(i-1,j,k,qpres)-dpe(i+3,j,k)*q(i  ,j,k,qpres)) &
                  + m15*(dpe(i-4,j,k)*q(i  ,j,k,qpres)-dpe(i+3,j,k)*q(i-1,j,k,qpres)) &
                  + m21*(dpe(i-3,j,k)*q(i-4,j,k,qpres)-dpe(i+2,j,k)*q(i+3,j,k,qpres)) &
                  + m22*(dpe(i-3,j,k)*q(i-3,j,k,qpres)-dpe(i+2,j,k)*q(i+2,j,k,qpres)) &
                  + m23*(dpe(i-3,j,k)*q(i-2,j,k,qpres)-dpe(i+2,j,k)*q(i+1,j,k,qpres)) &
                  + m24*(dpe(i-3,j,k)*q(i-1,j,k,qpres)-dpe(i+2,j,k)*q(i  ,j,k,qpres)) &
                  + m25*(dpe(i-3,j,k)*q(i  ,j,k,qpres)-dpe(i+2,j,k)*q(i-1,j,k,qpres)) &
                  + m26*(dpe(i-3,j,k)*q(i+1,j,k,qpres)-dpe(i+2,j,k)*q(i-2,j,k,qpres)) &
                  + m31*(dpe(i-2,j,k)*q(i-4,j,k,qpres)-dpe(i+1,j,k)*q(i+3,j,k,qpres)) &
                  + m32*(dpe(i-2,j,k)*q(i-3,j,k,qpres)-dpe(i+1,j,k)*q(i+2,j,k,qpres)) &
                  + m33*(dpe(i-2,j,k)*q(i-2,j,k,qpres)-dpe(i+1,j,k)*q(i+1,j,k,qpres)) &
                  + m34*(dpe(i-2,j,k)*q(i-1,j,k,qpres)-dpe(i+1,j,k)*q(i  ,j,k,qpres)) &
                  + m35*(dpe(i-2,j,k)*q(i  ,j,k,qpres)-dpe(i+1,j,k)*q(i-1,j,k,qpres)) &
                  + m36*(dpe(i-2,j,k)*q(i+1,j,k,qpres)-dpe(i+1,j,k)*q(i-2,j,k,qpres)) &
                  + m37*(dpe(i-2,j,k)*q(i+2,j,k,qpres)-dpe(i+1,j,k)*q(i-3,j,k,qpres)) &
                  + m41*(dpe(i-1,j,k)*q(i-4,j,k,qpres)-dpe(i  ,j,k)*q(i+3,j,k,qpres)) &
                  + m42*(dpe(i-1,j,k)*q(i-3,j,k,qpres)-dpe(i  ,j,k)*q(i+2,j,k,qpres)) &
                  + m43*(dpe(i-1,j,k)*q(i-2,j,k,qpres)-dpe(i  ,j,k)*q(i+1,j,k,qpres)) &
                  + m44*(dpe(i-1,j,k)*q(i-1,j,k,qpres)-dpe(i  ,j,k)*q(i  ,j,k,qpres)) &
                  + m45*(dpe(i-1,j,k)*q(i  ,j,k,qpres)-dpe(i  ,j,k)*q(i-1,j,k,qpres)) &
                  + m46*(dpe(i-1,j,k)*q(i+1,j,k,qpres)-dpe(i  ,j,k)*q(i-2,j,k,qpres)) &
                  + m47*(dpe(i-1,j,k)*q(i+2,j,k,qpres)-dpe(i  ,j,k)*q(i-3,j,k,qpres)) &
                  + m48*(dpe(i-1,j,k)*q(i+3,j,k,qpres)-dpe(i  ,j,k)*q(i-4,j,k,qpres))

             do n = 1, nspecies
                qhn = qh1+n-1
                hhalf = (q(i-1,j,k,qhn) + q(i,j,k,qhn)) / 2.d0
                Hg(i,j,k,iene) =  Hg(i,j,k,iene) - Ytmp(n) * hhalf * Htot
             end do

          end do
       end do
    end do
    !$omp end do

    ! add x-direction flux
    do n=2,ncons
       !$omp do
       do k=lo(3),hi(3)
          do j=lo(2),hi(2)
             do i=lo(1),hi(1)
                flx(i,j,k,n) = flx(i,j,k,n) + (Hg(i+1,j,k,n) - Hg(i,j,k,n)) * dx2inv(1)
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
             Hg(i,j,k,imx) = m11*(mu(i,j-4,k)*q(i,j-4,k,qu)-mu(i,j+3,k)*q(i,j+3,k,qu)) &
                  +          m12*(mu(i,j-4,k)*q(i,j-3,k,qu)-mu(i,j+3,k)*q(i,j+2,k,qu)) &
                  +          m13*(mu(i,j-4,k)*q(i,j-2,k,qu)-mu(i,j+3,k)*q(i,j+1,k,qu)) &
                  +          m14*(mu(i,j-4,k)*q(i,j-1,k,qu)-mu(i,j+3,k)*q(i,j  ,k,qu)) &
                  +          m15*(mu(i,j-4,k)*q(i,j  ,k,qu)-mu(i,j+3,k)*q(i,j-1,k,qu)) &
                  &        + m21*(mu(i,j-3,k)*q(i,j-4,k,qu)-mu(i,j+2,k)*q(i,j+3,k,qu)) &
                  +          m22*(mu(i,j-3,k)*q(i,j-3,k,qu)-mu(i,j+2,k)*q(i,j+2,k,qu)) &
                  +          m23*(mu(i,j-3,k)*q(i,j-2,k,qu)-mu(i,j+2,k)*q(i,j+1,k,qu)) &
                  +          m24*(mu(i,j-3,k)*q(i,j-1,k,qu)-mu(i,j+2,k)*q(i,j  ,k,qu)) &
                  +          m25*(mu(i,j-3,k)*q(i,j  ,k,qu)-mu(i,j+2,k)*q(i,j-1,k,qu)) &
                  +          m26*(mu(i,j-3,k)*q(i,j+1,k,qu)-mu(i,j+2,k)*q(i,j-2,k,qu)) &
                  &        + m31*(mu(i,j-2,k)*q(i,j-4,k,qu)-mu(i,j+1,k)*q(i,j+3,k,qu)) &
                  +          m32*(mu(i,j-2,k)*q(i,j-3,k,qu)-mu(i,j+1,k)*q(i,j+2,k,qu)) &
                  +          m33*(mu(i,j-2,k)*q(i,j-2,k,qu)-mu(i,j+1,k)*q(i,j+1,k,qu)) &
                  +          m34*(mu(i,j-2,k)*q(i,j-1,k,qu)-mu(i,j+1,k)*q(i,j  ,k,qu)) &
                  +          m35*(mu(i,j-2,k)*q(i,j  ,k,qu)-mu(i,j+1,k)*q(i,j-1,k,qu)) &
                  +          m36*(mu(i,j-2,k)*q(i,j+1,k,qu)-mu(i,j+1,k)*q(i,j-2,k,qu)) &
                  +          m37*(mu(i,j-2,k)*q(i,j+2,k,qu)-mu(i,j+1,k)*q(i,j-3,k,qu)) &
                  &        + m41*(mu(i,j-1,k)*q(i,j-4,k,qu)-mu(i,j  ,k)*q(i,j+3,k,qu)) &
                  +          m42*(mu(i,j-1,k)*q(i,j-3,k,qu)-mu(i,j  ,k)*q(i,j+2,k,qu)) &
                  +          m43*(mu(i,j-1,k)*q(i,j-2,k,qu)-mu(i,j  ,k)*q(i,j+1,k,qu)) &
                  +          m44*(mu(i,j-1,k)*q(i,j-1,k,qu)-mu(i,j  ,k)*q(i,j  ,k,qu)) &
                  +          m45*(mu(i,j-1,k)*q(i,j  ,k,qu)-mu(i,j  ,k)*q(i,j-1,k,qu)) &
                  +          m46*(mu(i,j-1,k)*q(i,j+1,k,qu)-mu(i,j  ,k)*q(i,j-2,k,qu)) &
                  +          m47*(mu(i,j-1,k)*q(i,j+2,k,qu)-mu(i,j  ,k)*q(i,j-3,k,qu)) &
                  +          m48*(mu(i,j-1,k)*q(i,j+3,k,qu)-mu(i,j  ,k)*q(i,j-4,k,qu))

             Hg(i,j,k,imy) = m11*(vsp(i,j-4,k)*q(i,j-4,k,qv)-vsp(i,j+3,k)*q(i,j+3,k,qv)) &
                  +          m12*(vsp(i,j-4,k)*q(i,j-3,k,qv)-vsp(i,j+3,k)*q(i,j+2,k,qv)) &
                  +          m13*(vsp(i,j-4,k)*q(i,j-2,k,qv)-vsp(i,j+3,k)*q(i,j+1,k,qv)) &
                  +          m14*(vsp(i,j-4,k)*q(i,j-1,k,qv)-vsp(i,j+3,k)*q(i,j  ,k,qv)) &
                  +          m15*(vsp(i,j-4,k)*q(i,j  ,k,qv)-vsp(i,j+3,k)*q(i,j-1,k,qv)) &
                  &        + m21*(vsp(i,j-3,k)*q(i,j-4,k,qv)-vsp(i,j+2,k)*q(i,j+3,k,qv)) &
                  +          m22*(vsp(i,j-3,k)*q(i,j-3,k,qv)-vsp(i,j+2,k)*q(i,j+2,k,qv)) &
                  +          m23*(vsp(i,j-3,k)*q(i,j-2,k,qv)-vsp(i,j+2,k)*q(i,j+1,k,qv)) &
                  +          m24*(vsp(i,j-3,k)*q(i,j-1,k,qv)-vsp(i,j+2,k)*q(i,j  ,k,qv)) &
                  +          m25*(vsp(i,j-3,k)*q(i,j  ,k,qv)-vsp(i,j+2,k)*q(i,j-1,k,qv)) &
                  +          m26*(vsp(i,j-3,k)*q(i,j+1,k,qv)-vsp(i,j+2,k)*q(i,j-2,k,qv)) &
                  &        + m31*(vsp(i,j-2,k)*q(i,j-4,k,qv)-vsp(i,j+1,k)*q(i,j+3,k,qv)) &
                  +          m32*(vsp(i,j-2,k)*q(i,j-3,k,qv)-vsp(i,j+1,k)*q(i,j+2,k,qv)) &
                  +          m33*(vsp(i,j-2,k)*q(i,j-2,k,qv)-vsp(i,j+1,k)*q(i,j+1,k,qv)) &
                  +          m34*(vsp(i,j-2,k)*q(i,j-1,k,qv)-vsp(i,j+1,k)*q(i,j  ,k,qv)) &
                  +          m35*(vsp(i,j-2,k)*q(i,j  ,k,qv)-vsp(i,j+1,k)*q(i,j-1,k,qv)) &
                  +          m36*(vsp(i,j-2,k)*q(i,j+1,k,qv)-vsp(i,j+1,k)*q(i,j-2,k,qv)) &
                  +          m37*(vsp(i,j-2,k)*q(i,j+2,k,qv)-vsp(i,j+1,k)*q(i,j-3,k,qv)) &
                  &        + m41*(vsp(i,j-1,k)*q(i,j-4,k,qv)-vsp(i,j  ,k)*q(i,j+3,k,qv)) &
                  +          m42*(vsp(i,j-1,k)*q(i,j-3,k,qv)-vsp(i,j  ,k)*q(i,j+2,k,qv)) &
                  +          m43*(vsp(i,j-1,k)*q(i,j-2,k,qv)-vsp(i,j  ,k)*q(i,j+1,k,qv)) &
                  +          m44*(vsp(i,j-1,k)*q(i,j-1,k,qv)-vsp(i,j  ,k)*q(i,j  ,k,qv)) &
                  +          m45*(vsp(i,j-1,k)*q(i,j  ,k,qv)-vsp(i,j  ,k)*q(i,j-1,k,qv)) &
                  +          m46*(vsp(i,j-1,k)*q(i,j+1,k,qv)-vsp(i,j  ,k)*q(i,j-2,k,qv)) &
                  +          m47*(vsp(i,j-1,k)*q(i,j+2,k,qv)-vsp(i,j  ,k)*q(i,j-3,k,qv)) &
                  +          m48*(vsp(i,j-1,k)*q(i,j+3,k,qv)-vsp(i,j  ,k)*q(i,j-4,k,qv))

             Hg(i,j,k,imz) = m11*(mu(i,j-4,k)*q(i,j-4,k,qw)-mu(i,j+3,k)*q(i,j+3,k,qw)) &
                  +          m12*(mu(i,j-4,k)*q(i,j-3,k,qw)-mu(i,j+3,k)*q(i,j+2,k,qw)) &
                  +          m13*(mu(i,j-4,k)*q(i,j-2,k,qw)-mu(i,j+3,k)*q(i,j+1,k,qw)) &
                  +          m14*(mu(i,j-4,k)*q(i,j-1,k,qw)-mu(i,j+3,k)*q(i,j  ,k,qw)) &
                  +          m15*(mu(i,j-4,k)*q(i,j  ,k,qw)-mu(i,j+3,k)*q(i,j-1,k,qw)) &
                  &        + m21*(mu(i,j-3,k)*q(i,j-4,k,qw)-mu(i,j+2,k)*q(i,j+3,k,qw)) &
                  +          m22*(mu(i,j-3,k)*q(i,j-3,k,qw)-mu(i,j+2,k)*q(i,j+2,k,qw)) &
                  +          m23*(mu(i,j-3,k)*q(i,j-2,k,qw)-mu(i,j+2,k)*q(i,j+1,k,qw)) &
                  +          m24*(mu(i,j-3,k)*q(i,j-1,k,qw)-mu(i,j+2,k)*q(i,j  ,k,qw)) &
                  +          m25*(mu(i,j-3,k)*q(i,j  ,k,qw)-mu(i,j+2,k)*q(i,j-1,k,qw)) &
                  +          m26*(mu(i,j-3,k)*q(i,j+1,k,qw)-mu(i,j+2,k)*q(i,j-2,k,qw)) &
                  &        + m31*(mu(i,j-2,k)*q(i,j-4,k,qw)-mu(i,j+1,k)*q(i,j+3,k,qw)) &
                  +          m32*(mu(i,j-2,k)*q(i,j-3,k,qw)-mu(i,j+1,k)*q(i,j+2,k,qw)) &
                  +          m33*(mu(i,j-2,k)*q(i,j-2,k,qw)-mu(i,j+1,k)*q(i,j+1,k,qw)) &
                  +          m34*(mu(i,j-2,k)*q(i,j-1,k,qw)-mu(i,j+1,k)*q(i,j  ,k,qw)) &
                  +          m35*(mu(i,j-2,k)*q(i,j  ,k,qw)-mu(i,j+1,k)*q(i,j-1,k,qw)) &
                  +          m36*(mu(i,j-2,k)*q(i,j+1,k,qw)-mu(i,j+1,k)*q(i,j-2,k,qw)) &
                  +          m37*(mu(i,j-2,k)*q(i,j+2,k,qw)-mu(i,j+1,k)*q(i,j-3,k,qw)) &
                  &        + m41*(mu(i,j-1,k)*q(i,j-4,k,qw)-mu(i,j  ,k)*q(i,j+3,k,qw)) &
                  +          m42*(mu(i,j-1,k)*q(i,j-3,k,qw)-mu(i,j  ,k)*q(i,j+2,k,qw)) &
                  +          m43*(mu(i,j-1,k)*q(i,j-2,k,qw)-mu(i,j  ,k)*q(i,j+1,k,qw)) &
                  +          m44*(mu(i,j-1,k)*q(i,j-1,k,qw)-mu(i,j  ,k)*q(i,j  ,k,qw)) &
                  +          m45*(mu(i,j-1,k)*q(i,j  ,k,qw)-mu(i,j  ,k)*q(i,j-1,k,qw)) &
                  +          m46*(mu(i,j-1,k)*q(i,j+1,k,qw)-mu(i,j  ,k)*q(i,j-2,k,qw)) &
                  +          m47*(mu(i,j-1,k)*q(i,j+2,k,qw)-mu(i,j  ,k)*q(i,j-3,k,qw)) &
                  +          m48*(mu(i,j-1,k)*q(i,j+3,k,qw)-mu(i,j  ,k)*q(i,j-4,k,qw))

             Hg(i,j,k,iene) = m11*(lam(i,j-4,k)*q(i,j-4,k,qtemp)-lam(i,j+3,k)*q(i,j+3,k,qtemp)) &
                  +           m12*(lam(i,j-4,k)*q(i,j-3,k,qtemp)-lam(i,j+3,k)*q(i,j+2,k,qtemp)) &
                  +           m13*(lam(i,j-4,k)*q(i,j-2,k,qtemp)-lam(i,j+3,k)*q(i,j+1,k,qtemp)) &
                  +           m14*(lam(i,j-4,k)*q(i,j-1,k,qtemp)-lam(i,j+3,k)*q(i,j  ,k,qtemp)) &
                  +           m15*(lam(i,j-4,k)*q(i,j  ,k,qtemp)-lam(i,j+3,k)*q(i,j-1,k,qtemp)) &
                  &         + m21*(lam(i,j-3,k)*q(i,j-4,k,qtemp)-lam(i,j+2,k)*q(i,j+3,k,qtemp)) &
                  +           m22*(lam(i,j-3,k)*q(i,j-3,k,qtemp)-lam(i,j+2,k)*q(i,j+2,k,qtemp)) &
                  +           m23*(lam(i,j-3,k)*q(i,j-2,k,qtemp)-lam(i,j+2,k)*q(i,j+1,k,qtemp)) &
                  +           m24*(lam(i,j-3,k)*q(i,j-1,k,qtemp)-lam(i,j+2,k)*q(i,j  ,k,qtemp)) &
                  +           m25*(lam(i,j-3,k)*q(i,j  ,k,qtemp)-lam(i,j+2,k)*q(i,j-1,k,qtemp)) &
                  +           m26*(lam(i,j-3,k)*q(i,j+1,k,qtemp)-lam(i,j+2,k)*q(i,j-2,k,qtemp)) &
                  &         + m31*(lam(i,j-2,k)*q(i,j-4,k,qtemp)-lam(i,j+1,k)*q(i,j+3,k,qtemp)) &
                  +           m32*(lam(i,j-2,k)*q(i,j-3,k,qtemp)-lam(i,j+1,k)*q(i,j+2,k,qtemp)) &
                  +           m33*(lam(i,j-2,k)*q(i,j-2,k,qtemp)-lam(i,j+1,k)*q(i,j+1,k,qtemp)) &
                  +           m34*(lam(i,j-2,k)*q(i,j-1,k,qtemp)-lam(i,j+1,k)*q(i,j  ,k,qtemp)) &
                  +           m35*(lam(i,j-2,k)*q(i,j  ,k,qtemp)-lam(i,j+1,k)*q(i,j-1,k,qtemp)) &
                  +           m36*(lam(i,j-2,k)*q(i,j+1,k,qtemp)-lam(i,j+1,k)*q(i,j-2,k,qtemp)) &
                  +           m37*(lam(i,j-2,k)*q(i,j+2,k,qtemp)-lam(i,j+1,k)*q(i,j-3,k,qtemp)) &
                  &         + m41*(lam(i,j-1,k)*q(i,j-4,k,qtemp)-lam(i,j  ,k)*q(i,j+3,k,qtemp)) &
                  +           m42*(lam(i,j-1,k)*q(i,j-3,k,qtemp)-lam(i,j  ,k)*q(i,j+2,k,qtemp)) &
                  +           m43*(lam(i,j-1,k)*q(i,j-2,k,qtemp)-lam(i,j  ,k)*q(i,j+1,k,qtemp)) &
                  +           m44*(lam(i,j-1,k)*q(i,j-1,k,qtemp)-lam(i,j  ,k)*q(i,j  ,k,qtemp)) &
                  +           m45*(lam(i,j-1,k)*q(i,j  ,k,qtemp)-lam(i,j  ,k)*q(i,j-1,k,qtemp)) &
                  +           m46*(lam(i,j-1,k)*q(i,j+1,k,qtemp)-lam(i,j  ,k)*q(i,j-2,k,qtemp)) &
                  +           m47*(lam(i,j-1,k)*q(i,j+2,k,qtemp)-lam(i,j  ,k)*q(i,j-3,k,qtemp)) &
                  +           m48*(lam(i,j-1,k)*q(i,j+3,k,qtemp)-lam(i,j  ,k)*q(i,j-4,k,qtemp))

             Htot = 0.d0
             do n = 1, nspecies
                qxn = qx1+n-1
                qyn = qy1+n-1
                Htmp(n) = m11*(dxy(i,j-4,k,n)*q(i,j-4,k,qxn)-dxy(i,j+3,k,n)*q(i,j+3,k,qxn)) &
                  +       m12*(dxy(i,j-4,k,n)*q(i,j-3,k,qxn)-dxy(i,j+3,k,n)*q(i,j+2,k,qxn)) &
                  +       m13*(dxy(i,j-4,k,n)*q(i,j-2,k,qxn)-dxy(i,j+3,k,n)*q(i,j+1,k,qxn)) &
                  +       m14*(dxy(i,j-4,k,n)*q(i,j-1,k,qxn)-dxy(i,j+3,k,n)*q(i,j  ,k,qxn)) &
                  +       m15*(dxy(i,j-4,k,n)*q(i,j  ,k,qxn)-dxy(i,j+3,k,n)*q(i,j-1,k,qxn)) &
                  &     + m21*(dxy(i,j-3,k,n)*q(i,j-4,k,qxn)-dxy(i,j+2,k,n)*q(i,j+3,k,qxn)) &
                  +       m22*(dxy(i,j-3,k,n)*q(i,j-3,k,qxn)-dxy(i,j+2,k,n)*q(i,j+2,k,qxn)) &
                  +       m23*(dxy(i,j-3,k,n)*q(i,j-2,k,qxn)-dxy(i,j+2,k,n)*q(i,j+1,k,qxn)) &
                  +       m24*(dxy(i,j-3,k,n)*q(i,j-1,k,qxn)-dxy(i,j+2,k,n)*q(i,j  ,k,qxn)) &
                  +       m25*(dxy(i,j-3,k,n)*q(i,j  ,k,qxn)-dxy(i,j+2,k,n)*q(i,j-1,k,qxn)) &
                  +       m26*(dxy(i,j-3,k,n)*q(i,j+1,k,qxn)-dxy(i,j+2,k,n)*q(i,j-2,k,qxn)) &
                  &     + m31*(dxy(i,j-2,k,n)*q(i,j-4,k,qxn)-dxy(i,j+1,k,n)*q(i,j+3,k,qxn)) &
                  +       m32*(dxy(i,j-2,k,n)*q(i,j-3,k,qxn)-dxy(i,j+1,k,n)*q(i,j+2,k,qxn)) &
                  +       m33*(dxy(i,j-2,k,n)*q(i,j-2,k,qxn)-dxy(i,j+1,k,n)*q(i,j+1,k,qxn)) &
                  +       m34*(dxy(i,j-2,k,n)*q(i,j-1,k,qxn)-dxy(i,j+1,k,n)*q(i,j  ,k,qxn)) &
                  +       m35*(dxy(i,j-2,k,n)*q(i,j  ,k,qxn)-dxy(i,j+1,k,n)*q(i,j-1,k,qxn)) &
                  +       m36*(dxy(i,j-2,k,n)*q(i,j+1,k,qxn)-dxy(i,j+1,k,n)*q(i,j-2,k,qxn)) &
                  +       m37*(dxy(i,j-2,k,n)*q(i,j+2,k,qxn)-dxy(i,j+1,k,n)*q(i,j-3,k,qxn)) &
                  &     + m41*(dxy(i,j-1,k,n)*q(i,j-4,k,qxn)-dxy(i,j  ,k,n)*q(i,j+3,k,qxn)) &
                  +       m42*(dxy(i,j-1,k,n)*q(i,j-3,k,qxn)-dxy(i,j  ,k,n)*q(i,j+2,k,qxn)) &
                  +       m43*(dxy(i,j-1,k,n)*q(i,j-2,k,qxn)-dxy(i,j  ,k,n)*q(i,j+1,k,qxn)) &
                  +       m44*(dxy(i,j-1,k,n)*q(i,j-1,k,qxn)-dxy(i,j  ,k,n)*q(i,j  ,k,qxn)) &
                  +       m45*(dxy(i,j-1,k,n)*q(i,j  ,k,qxn)-dxy(i,j  ,k,n)*q(i,j-1,k,qxn)) &
                  +       m46*(dxy(i,j-1,k,n)*q(i,j+1,k,qxn)-dxy(i,j  ,k,n)*q(i,j-2,k,qxn)) &
                  +       m47*(dxy(i,j-1,k,n)*q(i,j+2,k,qxn)-dxy(i,j  ,k,n)*q(i,j-3,k,qxn)) &
                  +       m48*(dxy(i,j-1,k,n)*q(i,j+3,k,qxn)-dxy(i,j  ,k,n)*q(i,j-4,k,qxn))
                Htmp(n) = Htmp(n)  &
                  +       m11*(dpy(i,j-4,k,n)*q(i,j-4,k,qpres)-dpy(i,j+3,k,n)*q(i,j+3,k,qpres)) &
                  +       m12*(dpy(i,j-4,k,n)*q(i,j-3,k,qpres)-dpy(i,j+3,k,n)*q(i,j+2,k,qpres)) &
                  +       m13*(dpy(i,j-4,k,n)*q(i,j-2,k,qpres)-dpy(i,j+3,k,n)*q(i,j+1,k,qpres)) &
                  +       m14*(dpy(i,j-4,k,n)*q(i,j-1,k,qpres)-dpy(i,j+3,k,n)*q(i,j  ,k,qpres)) &
                  +       m15*(dpy(i,j-4,k,n)*q(i,j  ,k,qpres)-dpy(i,j+3,k,n)*q(i,j-1,k,qpres)) &
                  &     + m21*(dpy(i,j-3,k,n)*q(i,j-4,k,qpres)-dpy(i,j+2,k,n)*q(i,j+3,k,qpres)) &
                  +       m22*(dpy(i,j-3,k,n)*q(i,j-3,k,qpres)-dpy(i,j+2,k,n)*q(i,j+2,k,qpres)) &
                  +       m23*(dpy(i,j-3,k,n)*q(i,j-2,k,qpres)-dpy(i,j+2,k,n)*q(i,j+1,k,qpres)) &
                  +       m24*(dpy(i,j-3,k,n)*q(i,j-1,k,qpres)-dpy(i,j+2,k,n)*q(i,j  ,k,qpres)) &
                  +       m25*(dpy(i,j-3,k,n)*q(i,j  ,k,qpres)-dpy(i,j+2,k,n)*q(i,j-1,k,qpres)) &
                  +       m26*(dpy(i,j-3,k,n)*q(i,j+1,k,qpres)-dpy(i,j+2,k,n)*q(i,j-2,k,qpres)) &
                  &     + m31*(dpy(i,j-2,k,n)*q(i,j-4,k,qpres)-dpy(i,j+1,k,n)*q(i,j+3,k,qpres)) &
                  +       m32*(dpy(i,j-2,k,n)*q(i,j-3,k,qpres)-dpy(i,j+1,k,n)*q(i,j+2,k,qpres)) &
                  +       m33*(dpy(i,j-2,k,n)*q(i,j-2,k,qpres)-dpy(i,j+1,k,n)*q(i,j+1,k,qpres)) &
                  +       m34*(dpy(i,j-2,k,n)*q(i,j-1,k,qpres)-dpy(i,j+1,k,n)*q(i,j  ,k,qpres)) &
                  +       m35*(dpy(i,j-2,k,n)*q(i,j  ,k,qpres)-dpy(i,j+1,k,n)*q(i,j-1,k,qpres)) &
                  +       m36*(dpy(i,j-2,k,n)*q(i,j+1,k,qpres)-dpy(i,j+1,k,n)*q(i,j-2,k,qpres)) &
                  +       m37*(dpy(i,j-2,k,n)*q(i,j+2,k,qpres)-dpy(i,j+1,k,n)*q(i,j-3,k,qpres)) &
                  &     + m41*(dpy(i,j-1,k,n)*q(i,j-4,k,qpres)-dpy(i,j  ,k,n)*q(i,j+3,k,qpres)) &
                  +       m42*(dpy(i,j-1,k,n)*q(i,j-3,k,qpres)-dpy(i,j  ,k,n)*q(i,j+2,k,qpres)) &
                  +       m43*(dpy(i,j-1,k,n)*q(i,j-2,k,qpres)-dpy(i,j  ,k,n)*q(i,j+1,k,qpres)) &
                  +       m44*(dpy(i,j-1,k,n)*q(i,j-1,k,qpres)-dpy(i,j  ,k,n)*q(i,j  ,k,qpres)) &
                  +       m45*(dpy(i,j-1,k,n)*q(i,j  ,k,qpres)-dpy(i,j  ,k,n)*q(i,j-1,k,qpres)) &
                  +       m46*(dpy(i,j-1,k,n)*q(i,j+1,k,qpres)-dpy(i,j  ,k,n)*q(i,j-2,k,qpres)) &
                  +       m47*(dpy(i,j-1,k,n)*q(i,j+2,k,qpres)-dpy(i,j  ,k,n)*q(i,j-3,k,qpres)) &
                  +       m48*(dpy(i,j-1,k,n)*q(i,j+3,k,qpres)-dpy(i,j  ,k,n)*q(i,j-4,k,qpres))
                Htot = Htot + Htmp(n)
                Ytmp(n) = (q(i,j-1,k,qyn) + q(i,j,k,qyn)) / 2.d0
             end do

             do n = 1, nspecies
                Hg(i,j,k,iry1+n-1) = Htmp(n) - Ytmp(n)*Htot
             end do

             do n = 1, nspecies
                qxn = qx1+n-1
                Hg(i,j,k,iene) =  Hg(i,j,k,iene) &
                  + m11*(dxe(i,j-4,k,n)*q(i,j-4,k,qxn)-dxe(i,j+3,k,n)*q(i,j+3,k,qxn)) &
                  + m12*(dxe(i,j-4,k,n)*q(i,j-3,k,qxn)-dxe(i,j+3,k,n)*q(i,j+2,k,qxn)) &
                  + m13*(dxe(i,j-4,k,n)*q(i,j-2,k,qxn)-dxe(i,j+3,k,n)*q(i,j+1,k,qxn)) &
                  + m14*(dxe(i,j-4,k,n)*q(i,j-1,k,qxn)-dxe(i,j+3,k,n)*q(i,j  ,k,qxn)) &
                  + m15*(dxe(i,j-4,k,n)*q(i,j  ,k,qxn)-dxe(i,j+3,k,n)*q(i,j-1,k,qxn)) &
                  + m21*(dxe(i,j-3,k,n)*q(i,j-4,k,qxn)-dxe(i,j+2,k,n)*q(i,j+3,k,qxn)) &
                  + m22*(dxe(i,j-3,k,n)*q(i,j-3,k,qxn)-dxe(i,j+2,k,n)*q(i,j+2,k,qxn)) &
                  + m23*(dxe(i,j-3,k,n)*q(i,j-2,k,qxn)-dxe(i,j+2,k,n)*q(i,j+1,k,qxn)) &
                  + m24*(dxe(i,j-3,k,n)*q(i,j-1,k,qxn)-dxe(i,j+2,k,n)*q(i,j  ,k,qxn)) &
                  + m25*(dxe(i,j-3,k,n)*q(i,j  ,k,qxn)-dxe(i,j+2,k,n)*q(i,j-1,k,qxn)) &
                  + m26*(dxe(i,j-3,k,n)*q(i,j+1,k,qxn)-dxe(i,j+2,k,n)*q(i,j-2,k,qxn)) &
                  + m31*(dxe(i,j-2,k,n)*q(i,j-4,k,qxn)-dxe(i,j+1,k,n)*q(i,j+3,k,qxn)) &
                  + m32*(dxe(i,j-2,k,n)*q(i,j-3,k,qxn)-dxe(i,j+1,k,n)*q(i,j+2,k,qxn)) &
                  + m33*(dxe(i,j-2,k,n)*q(i,j-2,k,qxn)-dxe(i,j+1,k,n)*q(i,j+1,k,qxn)) &
                  + m34*(dxe(i,j-2,k,n)*q(i,j-1,k,qxn)-dxe(i,j+1,k,n)*q(i,j  ,k,qxn)) &
                  + m35*(dxe(i,j-2,k,n)*q(i,j  ,k,qxn)-dxe(i,j+1,k,n)*q(i,j-1,k,qxn)) &
                  + m36*(dxe(i,j-2,k,n)*q(i,j+1,k,qxn)-dxe(i,j+1,k,n)*q(i,j-2,k,qxn)) &
                  + m37*(dxe(i,j-2,k,n)*q(i,j+2,k,qxn)-dxe(i,j+1,k,n)*q(i,j-3,k,qxn)) &
                  + m41*(dxe(i,j-1,k,n)*q(i,j-4,k,qxn)-dxe(i,j  ,k,n)*q(i,j+3,k,qxn)) &
                  + m42*(dxe(i,j-1,k,n)*q(i,j-3,k,qxn)-dxe(i,j  ,k,n)*q(i,j+2,k,qxn)) &
                  + m43*(dxe(i,j-1,k,n)*q(i,j-2,k,qxn)-dxe(i,j  ,k,n)*q(i,j+1,k,qxn)) &
                  + m44*(dxe(i,j-1,k,n)*q(i,j-1,k,qxn)-dxe(i,j  ,k,n)*q(i,j  ,k,qxn)) &
                  + m45*(dxe(i,j-1,k,n)*q(i,j  ,k,qxn)-dxe(i,j  ,k,n)*q(i,j-1,k,qxn)) &
                  + m46*(dxe(i,j-1,k,n)*q(i,j+1,k,qxn)-dxe(i,j  ,k,n)*q(i,j-2,k,qxn)) &
                  + m47*(dxe(i,j-1,k,n)*q(i,j+2,k,qxn)-dxe(i,j  ,k,n)*q(i,j-3,k,qxn)) &
                  + m48*(dxe(i,j-1,k,n)*q(i,j+3,k,qxn)-dxe(i,j  ,k,n)*q(i,j-4,k,qxn))
             end do

             Hg(i,j,k,iene) =  Hg(i,j,k,iene) &
                  + m11*(dpe(i,j-4,k)*q(i,j-4,k,qpres)-dpe(i,j+3,k)*q(i,j+3,k,qpres)) &
                  + m12*(dpe(i,j-4,k)*q(i,j-3,k,qpres)-dpe(i,j+3,k)*q(i,j+2,k,qpres)) &
                  + m13*(dpe(i,j-4,k)*q(i,j-2,k,qpres)-dpe(i,j+3,k)*q(i,j+1,k,qpres)) &
                  + m14*(dpe(i,j-4,k)*q(i,j-1,k,qpres)-dpe(i,j+3,k)*q(i,j  ,k,qpres)) &
                  + m15*(dpe(i,j-4,k)*q(i,j  ,k,qpres)-dpe(i,j+3,k)*q(i,j-1,k,qpres)) &
                  + m21*(dpe(i,j-3,k)*q(i,j-4,k,qpres)-dpe(i,j+2,k)*q(i,j+3,k,qpres)) &
                  + m22*(dpe(i,j-3,k)*q(i,j-3,k,qpres)-dpe(i,j+2,k)*q(i,j+2,k,qpres)) &
                  + m23*(dpe(i,j-3,k)*q(i,j-2,k,qpres)-dpe(i,j+2,k)*q(i,j+1,k,qpres)) &
                  + m24*(dpe(i,j-3,k)*q(i,j-1,k,qpres)-dpe(i,j+2,k)*q(i,j  ,k,qpres)) &
                  + m25*(dpe(i,j-3,k)*q(i,j  ,k,qpres)-dpe(i,j+2,k)*q(i,j-1,k,qpres)) &
                  + m26*(dpe(i,j-3,k)*q(i,j+1,k,qpres)-dpe(i,j+2,k)*q(i,j-2,k,qpres)) &
                  + m31*(dpe(i,j-2,k)*q(i,j-4,k,qpres)-dpe(i,j+1,k)*q(i,j+3,k,qpres)) &
                  + m32*(dpe(i,j-2,k)*q(i,j-3,k,qpres)-dpe(i,j+1,k)*q(i,j+2,k,qpres)) &
                  + m33*(dpe(i,j-2,k)*q(i,j-2,k,qpres)-dpe(i,j+1,k)*q(i,j+1,k,qpres)) &
                  + m34*(dpe(i,j-2,k)*q(i,j-1,k,qpres)-dpe(i,j+1,k)*q(i,j  ,k,qpres)) &
                  + m35*(dpe(i,j-2,k)*q(i,j  ,k,qpres)-dpe(i,j+1,k)*q(i,j-1,k,qpres)) &
                  + m36*(dpe(i,j-2,k)*q(i,j+1,k,qpres)-dpe(i,j+1,k)*q(i,j-2,k,qpres)) &
                  + m37*(dpe(i,j-2,k)*q(i,j+2,k,qpres)-dpe(i,j+1,k)*q(i,j-3,k,qpres)) &
                  + m41*(dpe(i,j-1,k)*q(i,j-4,k,qpres)-dpe(i,j  ,k)*q(i,j+3,k,qpres)) &
                  + m42*(dpe(i,j-1,k)*q(i,j-3,k,qpres)-dpe(i,j  ,k)*q(i,j+2,k,qpres)) &
                  + m43*(dpe(i,j-1,k)*q(i,j-2,k,qpres)-dpe(i,j  ,k)*q(i,j+1,k,qpres)) &
                  + m44*(dpe(i,j-1,k)*q(i,j-1,k,qpres)-dpe(i,j  ,k)*q(i,j  ,k,qpres)) &
                  + m45*(dpe(i,j-1,k)*q(i,j  ,k,qpres)-dpe(i,j  ,k)*q(i,j-1,k,qpres)) &
                  + m46*(dpe(i,j-1,k)*q(i,j+1,k,qpres)-dpe(i,j  ,k)*q(i,j-2,k,qpres)) &
                  + m47*(dpe(i,j-1,k)*q(i,j+2,k,qpres)-dpe(i,j  ,k)*q(i,j-3,k,qpres)) &
                  + m48*(dpe(i,j-1,k)*q(i,j+3,k,qpres)-dpe(i,j  ,k)*q(i,j-4,k,qpres))

             do n = 1, nspecies
                qhn = qh1+n-1
                hhalf = (q(i,j-1,k,qhn) + q(i,j,k,qhn)) / 2.d0
                Hg(i,j,k,iene) =  Hg(i,j,k,iene) - Ytmp(n) * hhalf * Htot
             end do

          end do
       end do
    end do
    !$omp end do

    ! add y-direction flux
    do n=2,ncons
       !$omp do
       do k=lo(3),hi(3)
          do j=lo(2),hi(2)
             do i=lo(1),hi(1)
                flx(i,j,k,n) = flx(i,j,k,n) + (Hg(i,j+1,k,n) - Hg(i,j,k,n)) * dx2inv(2)
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
             Hg(i,j,k,imx) = m11*(mu(i,j,k-4)*q(i,j,k-4,qu)-mu(i,j,k+3)*q(i,j,k+3,qu)) &
                  +          m12*(mu(i,j,k-4)*q(i,j,k-3,qu)-mu(i,j,k+3)*q(i,j,k+2,qu)) &
                  +          m13*(mu(i,j,k-4)*q(i,j,k-2,qu)-mu(i,j,k+3)*q(i,j,k+1,qu)) &
                  +          m14*(mu(i,j,k-4)*q(i,j,k-1,qu)-mu(i,j,k+3)*q(i,j,k  ,qu)) &
                  +          m15*(mu(i,j,k-4)*q(i,j,k  ,qu)-mu(i,j,k+3)*q(i,j,k-1,qu)) &
                  &        + m21*(mu(i,j,k-3)*q(i,j,k-4,qu)-mu(i,j,k+2)*q(i,j,k+3,qu)) &
                  +          m22*(mu(i,j,k-3)*q(i,j,k-3,qu)-mu(i,j,k+2)*q(i,j,k+2,qu)) &
                  +          m23*(mu(i,j,k-3)*q(i,j,k-2,qu)-mu(i,j,k+2)*q(i,j,k+1,qu)) &
                  +          m24*(mu(i,j,k-3)*q(i,j,k-1,qu)-mu(i,j,k+2)*q(i,j,k  ,qu)) &
                  +          m25*(mu(i,j,k-3)*q(i,j,k  ,qu)-mu(i,j,k+2)*q(i,j,k-1,qu)) &
                  +          m26*(mu(i,j,k-3)*q(i,j,k+1,qu)-mu(i,j,k+2)*q(i,j,k-2,qu)) &
                  &        + m31*(mu(i,j,k-2)*q(i,j,k-4,qu)-mu(i,j,k+1)*q(i,j,k+3,qu)) &
                  +          m32*(mu(i,j,k-2)*q(i,j,k-3,qu)-mu(i,j,k+1)*q(i,j,k+2,qu)) &
                  +          m33*(mu(i,j,k-2)*q(i,j,k-2,qu)-mu(i,j,k+1)*q(i,j,k+1,qu)) &
                  +          m34*(mu(i,j,k-2)*q(i,j,k-1,qu)-mu(i,j,k+1)*q(i,j,k  ,qu)) &
                  +          m35*(mu(i,j,k-2)*q(i,j,k  ,qu)-mu(i,j,k+1)*q(i,j,k-1,qu)) &
                  +          m36*(mu(i,j,k-2)*q(i,j,k+1,qu)-mu(i,j,k+1)*q(i,j,k-2,qu)) &
                  +          m37*(mu(i,j,k-2)*q(i,j,k+2,qu)-mu(i,j,k+1)*q(i,j,k-3,qu)) &
                  &        + m41*(mu(i,j,k-1)*q(i,j,k-4,qu)-mu(i,j,k  )*q(i,j,k+3,qu)) &
                  +          m42*(mu(i,j,k-1)*q(i,j,k-3,qu)-mu(i,j,k  )*q(i,j,k+2,qu)) &
                  +          m43*(mu(i,j,k-1)*q(i,j,k-2,qu)-mu(i,j,k  )*q(i,j,k+1,qu)) &
                  +          m44*(mu(i,j,k-1)*q(i,j,k-1,qu)-mu(i,j,k  )*q(i,j,k  ,qu)) &
                  +          m45*(mu(i,j,k-1)*q(i,j,k  ,qu)-mu(i,j,k  )*q(i,j,k-1,qu)) &
                  +          m46*(mu(i,j,k-1)*q(i,j,k+1,qu)-mu(i,j,k  )*q(i,j,k-2,qu)) &
                  +          m47*(mu(i,j,k-1)*q(i,j,k+2,qu)-mu(i,j,k  )*q(i,j,k-3,qu)) &
                  +          m48*(mu(i,j,k-1)*q(i,j,k+3,qu)-mu(i,j,k  )*q(i,j,k-4,qu))

             Hg(i,j,k,imy) = m11*(mu(i,j,k-4)*q(i,j,k-4,qv)-mu(i,j,k+3)*q(i,j,k+3,qv)) &
                  +          m12*(mu(i,j,k-4)*q(i,j,k-3,qv)-mu(i,j,k+3)*q(i,j,k+2,qv)) &
                  +          m13*(mu(i,j,k-4)*q(i,j,k-2,qv)-mu(i,j,k+3)*q(i,j,k+1,qv)) &
                  +          m14*(mu(i,j,k-4)*q(i,j,k-1,qv)-mu(i,j,k+3)*q(i,j,k  ,qv)) &
                  +          m15*(mu(i,j,k-4)*q(i,j,k  ,qv)-mu(i,j,k+3)*q(i,j,k-1,qv)) &
                  &        + m21*(mu(i,j,k-3)*q(i,j,k-4,qv)-mu(i,j,k+2)*q(i,j,k+3,qv)) &
                  +          m22*(mu(i,j,k-3)*q(i,j,k-3,qv)-mu(i,j,k+2)*q(i,j,k+2,qv)) &
                  +          m23*(mu(i,j,k-3)*q(i,j,k-2,qv)-mu(i,j,k+2)*q(i,j,k+1,qv)) &
                  +          m24*(mu(i,j,k-3)*q(i,j,k-1,qv)-mu(i,j,k+2)*q(i,j,k  ,qv)) &
                  +          m25*(mu(i,j,k-3)*q(i,j,k  ,qv)-mu(i,j,k+2)*q(i,j,k-1,qv)) &
                  +          m26*(mu(i,j,k-3)*q(i,j,k+1,qv)-mu(i,j,k+2)*q(i,j,k-2,qv)) &
                  &        + m31*(mu(i,j,k-2)*q(i,j,k-4,qv)-mu(i,j,k+1)*q(i,j,k+3,qv)) &
                  +          m32*(mu(i,j,k-2)*q(i,j,k-3,qv)-mu(i,j,k+1)*q(i,j,k+2,qv)) &
                  +          m33*(mu(i,j,k-2)*q(i,j,k-2,qv)-mu(i,j,k+1)*q(i,j,k+1,qv)) &
                  +          m34*(mu(i,j,k-2)*q(i,j,k-1,qv)-mu(i,j,k+1)*q(i,j,k  ,qv)) &
                  +          m35*(mu(i,j,k-2)*q(i,j,k  ,qv)-mu(i,j,k+1)*q(i,j,k-1,qv)) &
                  +          m36*(mu(i,j,k-2)*q(i,j,k+1,qv)-mu(i,j,k+1)*q(i,j,k-2,qv)) &
                  +          m37*(mu(i,j,k-2)*q(i,j,k+2,qv)-mu(i,j,k+1)*q(i,j,k-3,qv)) &
                  &        + m41*(mu(i,j,k-1)*q(i,j,k-4,qv)-mu(i,j,k  )*q(i,j,k+3,qv)) &
                  +          m42*(mu(i,j,k-1)*q(i,j,k-3,qv)-mu(i,j,k  )*q(i,j,k+2,qv)) &
                  +          m43*(mu(i,j,k-1)*q(i,j,k-2,qv)-mu(i,j,k  )*q(i,j,k+1,qv)) &
                  +          m44*(mu(i,j,k-1)*q(i,j,k-1,qv)-mu(i,j,k  )*q(i,j,k  ,qv)) &
                  +          m45*(mu(i,j,k-1)*q(i,j,k  ,qv)-mu(i,j,k  )*q(i,j,k-1,qv)) &
                  +          m46*(mu(i,j,k-1)*q(i,j,k+1,qv)-mu(i,j,k  )*q(i,j,k-2,qv)) &
                  +          m47*(mu(i,j,k-1)*q(i,j,k+2,qv)-mu(i,j,k  )*q(i,j,k-3,qv)) &
                  +          m48*(mu(i,j,k-1)*q(i,j,k+3,qv)-mu(i,j,k  )*q(i,j,k-4,qv))

             Hg(i,j,k,imz) = m11*(vsp(i,j,k-4)*q(i,j,k-4,qw)-vsp(i,j,k+3)*q(i,j,k+3,qw)) &
                  +          m12*(vsp(i,j,k-4)*q(i,j,k-3,qw)-vsp(i,j,k+3)*q(i,j,k+2,qw)) &
                  +          m13*(vsp(i,j,k-4)*q(i,j,k-2,qw)-vsp(i,j,k+3)*q(i,j,k+1,qw)) &
                  +          m14*(vsp(i,j,k-4)*q(i,j,k-1,qw)-vsp(i,j,k+3)*q(i,j,k  ,qw)) &
                  +          m15*(vsp(i,j,k-4)*q(i,j,k  ,qw)-vsp(i,j,k+3)*q(i,j,k-1,qw)) &
                  &        + m21*(vsp(i,j,k-3)*q(i,j,k-4,qw)-vsp(i,j,k+2)*q(i,j,k+3,qw)) &
                  +          m22*(vsp(i,j,k-3)*q(i,j,k-3,qw)-vsp(i,j,k+2)*q(i,j,k+2,qw)) &
                  +          m23*(vsp(i,j,k-3)*q(i,j,k-2,qw)-vsp(i,j,k+2)*q(i,j,k+1,qw)) &
                  +          m24*(vsp(i,j,k-3)*q(i,j,k-1,qw)-vsp(i,j,k+2)*q(i,j,k  ,qw)) &
                  +          m25*(vsp(i,j,k-3)*q(i,j,k  ,qw)-vsp(i,j,k+2)*q(i,j,k-1,qw)) &
                  +          m26*(vsp(i,j,k-3)*q(i,j,k+1,qw)-vsp(i,j,k+2)*q(i,j,k-2,qw)) &
                  &        + m31*(vsp(i,j,k-2)*q(i,j,k-4,qw)-vsp(i,j,k+1)*q(i,j,k+3,qw)) &
                  +          m32*(vsp(i,j,k-2)*q(i,j,k-3,qw)-vsp(i,j,k+1)*q(i,j,k+2,qw)) &
                  +          m33*(vsp(i,j,k-2)*q(i,j,k-2,qw)-vsp(i,j,k+1)*q(i,j,k+1,qw)) &
                  +          m34*(vsp(i,j,k-2)*q(i,j,k-1,qw)-vsp(i,j,k+1)*q(i,j,k  ,qw)) &
                  +          m35*(vsp(i,j,k-2)*q(i,j,k  ,qw)-vsp(i,j,k+1)*q(i,j,k-1,qw)) &
                  +          m36*(vsp(i,j,k-2)*q(i,j,k+1,qw)-vsp(i,j,k+1)*q(i,j,k-2,qw)) &
                  +          m37*(vsp(i,j,k-2)*q(i,j,k+2,qw)-vsp(i,j,k+1)*q(i,j,k-3,qw)) &
                  &        + m41*(vsp(i,j,k-1)*q(i,j,k-4,qw)-vsp(i,j,k  )*q(i,j,k+3,qw)) &
                  +          m42*(vsp(i,j,k-1)*q(i,j,k-3,qw)-vsp(i,j,k  )*q(i,j,k+2,qw)) &
                  +          m43*(vsp(i,j,k-1)*q(i,j,k-2,qw)-vsp(i,j,k  )*q(i,j,k+1,qw)) &
                  +          m44*(vsp(i,j,k-1)*q(i,j,k-1,qw)-vsp(i,j,k  )*q(i,j,k  ,qw)) &
                  +          m45*(vsp(i,j,k-1)*q(i,j,k  ,qw)-vsp(i,j,k  )*q(i,j,k-1,qw)) &
                  +          m46*(vsp(i,j,k-1)*q(i,j,k+1,qw)-vsp(i,j,k  )*q(i,j,k-2,qw)) &
                  +          m47*(vsp(i,j,k-1)*q(i,j,k+2,qw)-vsp(i,j,k  )*q(i,j,k-3,qw)) &
                  +          m48*(vsp(i,j,k-1)*q(i,j,k+3,qw)-vsp(i,j,k  )*q(i,j,k-4,qw))

             Hg(i,j,k,iene) = m11*(lam(i,j,k-4)*q(i,j,k-4,qtemp)-lam(i,j,k+3)*q(i,j,k+3,qtemp)) &
                  +           m12*(lam(i,j,k-4)*q(i,j,k-3,qtemp)-lam(i,j,k+3)*q(i,j,k+2,qtemp)) &
                  +           m13*(lam(i,j,k-4)*q(i,j,k-2,qtemp)-lam(i,j,k+3)*q(i,j,k+1,qtemp)) &
                  +           m14*(lam(i,j,k-4)*q(i,j,k-1,qtemp)-lam(i,j,k+3)*q(i,j,k  ,qtemp)) &
                  +           m15*(lam(i,j,k-4)*q(i,j,k  ,qtemp)-lam(i,j,k+3)*q(i,j,k-1,qtemp)) &
                  &         + m21*(lam(i,j,k-3)*q(i,j,k-4,qtemp)-lam(i,j,k+2)*q(i,j,k+3,qtemp)) &
                  +           m22*(lam(i,j,k-3)*q(i,j,k-3,qtemp)-lam(i,j,k+2)*q(i,j,k+2,qtemp)) &
                  +           m23*(lam(i,j,k-3)*q(i,j,k-2,qtemp)-lam(i,j,k+2)*q(i,j,k+1,qtemp)) &
                  +           m24*(lam(i,j,k-3)*q(i,j,k-1,qtemp)-lam(i,j,k+2)*q(i,j,k  ,qtemp)) &
                  +           m25*(lam(i,j,k-3)*q(i,j,k  ,qtemp)-lam(i,j,k+2)*q(i,j,k-1,qtemp)) &
                  +           m26*(lam(i,j,k-3)*q(i,j,k+1,qtemp)-lam(i,j,k+2)*q(i,j,k-2,qtemp)) &
                  &         + m31*(lam(i,j,k-2)*q(i,j,k-4,qtemp)-lam(i,j,k+1)*q(i,j,k+3,qtemp)) &
                  +           m32*(lam(i,j,k-2)*q(i,j,k-3,qtemp)-lam(i,j,k+1)*q(i,j,k+2,qtemp)) &
                  +           m33*(lam(i,j,k-2)*q(i,j,k-2,qtemp)-lam(i,j,k+1)*q(i,j,k+1,qtemp)) &
                  +           m34*(lam(i,j,k-2)*q(i,j,k-1,qtemp)-lam(i,j,k+1)*q(i,j,k  ,qtemp)) &
                  +           m35*(lam(i,j,k-2)*q(i,j,k  ,qtemp)-lam(i,j,k+1)*q(i,j,k-1,qtemp)) &
                  +           m36*(lam(i,j,k-2)*q(i,j,k+1,qtemp)-lam(i,j,k+1)*q(i,j,k-2,qtemp)) &
                  +           m37*(lam(i,j,k-2)*q(i,j,k+2,qtemp)-lam(i,j,k+1)*q(i,j,k-3,qtemp)) &
                  &         + m41*(lam(i,j,k-1)*q(i,j,k-4,qtemp)-lam(i,j,k  )*q(i,j,k+3,qtemp)) &
                  +           m42*(lam(i,j,k-1)*q(i,j,k-3,qtemp)-lam(i,j,k  )*q(i,j,k+2,qtemp)) &
                  +           m43*(lam(i,j,k-1)*q(i,j,k-2,qtemp)-lam(i,j,k  )*q(i,j,k+1,qtemp)) &
                  +           m44*(lam(i,j,k-1)*q(i,j,k-1,qtemp)-lam(i,j,k  )*q(i,j,k  ,qtemp)) &
                  +           m45*(lam(i,j,k-1)*q(i,j,k  ,qtemp)-lam(i,j,k  )*q(i,j,k-1,qtemp)) &
                  +           m46*(lam(i,j,k-1)*q(i,j,k+1,qtemp)-lam(i,j,k  )*q(i,j,k-2,qtemp)) &
                  +           m47*(lam(i,j,k-1)*q(i,j,k+2,qtemp)-lam(i,j,k  )*q(i,j,k-3,qtemp)) &
                  +           m48*(lam(i,j,k-1)*q(i,j,k+3,qtemp)-lam(i,j,k  )*q(i,j,k-4,qtemp))

             Htot = 0.d0
             do n = 1, nspecies
                qxn = qx1+n-1
                qyn = qy1+n-1
                Htmp(n) = m11*(dxy(i,j,k-4,n)*q(i,j,k-4,qxn)-dxy(i,j,k+3,n)*q(i,j,k+3,qxn)) &
                  +       m12*(dxy(i,j,k-4,n)*q(i,j,k-3,qxn)-dxy(i,j,k+3,n)*q(i,j,k+2,qxn)) &
                  +       m13*(dxy(i,j,k-4,n)*q(i,j,k-2,qxn)-dxy(i,j,k+3,n)*q(i,j,k+1,qxn)) &
                  +       m14*(dxy(i,j,k-4,n)*q(i,j,k-1,qxn)-dxy(i,j,k+3,n)*q(i,j,k  ,qxn)) &
                  +       m15*(dxy(i,j,k-4,n)*q(i,j,k  ,qxn)-dxy(i,j,k+3,n)*q(i,j,k-1,qxn)) &
                  &     + m21*(dxy(i,j,k-3,n)*q(i,j,k-4,qxn)-dxy(i,j,k+2,n)*q(i,j,k+3,qxn)) &
                  +       m22*(dxy(i,j,k-3,n)*q(i,j,k-3,qxn)-dxy(i,j,k+2,n)*q(i,j,k+2,qxn)) &
                  +       m23*(dxy(i,j,k-3,n)*q(i,j,k-2,qxn)-dxy(i,j,k+2,n)*q(i,j,k+1,qxn)) &
                  +       m24*(dxy(i,j,k-3,n)*q(i,j,k-1,qxn)-dxy(i,j,k+2,n)*q(i,j,k  ,qxn)) &
                  +       m25*(dxy(i,j,k-3,n)*q(i,j,k  ,qxn)-dxy(i,j,k+2,n)*q(i,j,k-1,qxn)) &
                  +       m26*(dxy(i,j,k-3,n)*q(i,j,k+1,qxn)-dxy(i,j,k+2,n)*q(i,j,k-2,qxn)) &
                  &     + m31*(dxy(i,j,k-2,n)*q(i,j,k-4,qxn)-dxy(i,j,k+1,n)*q(i,j,k+3,qxn)) &
                  +       m32*(dxy(i,j,k-2,n)*q(i,j,k-3,qxn)-dxy(i,j,k+1,n)*q(i,j,k+2,qxn)) &
                  +       m33*(dxy(i,j,k-2,n)*q(i,j,k-2,qxn)-dxy(i,j,k+1,n)*q(i,j,k+1,qxn)) &
                  +       m34*(dxy(i,j,k-2,n)*q(i,j,k-1,qxn)-dxy(i,j,k+1,n)*q(i,j,k  ,qxn)) &
                  +       m35*(dxy(i,j,k-2,n)*q(i,j,k  ,qxn)-dxy(i,j,k+1,n)*q(i,j,k-1,qxn)) &
                  +       m36*(dxy(i,j,k-2,n)*q(i,j,k+1,qxn)-dxy(i,j,k+1,n)*q(i,j,k-2,qxn)) &
                  +       m37*(dxy(i,j,k-2,n)*q(i,j,k+2,qxn)-dxy(i,j,k+1,n)*q(i,j,k-3,qxn)) &
                  &     + m41*(dxy(i,j,k-1,n)*q(i,j,k-4,qxn)-dxy(i,j,k  ,n)*q(i,j,k+3,qxn)) &
                  +       m42*(dxy(i,j,k-1,n)*q(i,j,k-3,qxn)-dxy(i,j,k  ,n)*q(i,j,k+2,qxn)) &
                  +       m43*(dxy(i,j,k-1,n)*q(i,j,k-2,qxn)-dxy(i,j,k  ,n)*q(i,j,k+1,qxn)) &
                  +       m44*(dxy(i,j,k-1,n)*q(i,j,k-1,qxn)-dxy(i,j,k  ,n)*q(i,j,k  ,qxn)) &
                  +       m45*(dxy(i,j,k-1,n)*q(i,j,k  ,qxn)-dxy(i,j,k  ,n)*q(i,j,k-1,qxn)) &
                  +       m46*(dxy(i,j,k-1,n)*q(i,j,k+1,qxn)-dxy(i,j,k  ,n)*q(i,j,k-2,qxn)) &
                  +       m47*(dxy(i,j,k-1,n)*q(i,j,k+2,qxn)-dxy(i,j,k  ,n)*q(i,j,k-3,qxn)) &
                  +       m48*(dxy(i,j,k-1,n)*q(i,j,k+3,qxn)-dxy(i,j,k  ,n)*q(i,j,k-4,qxn))
                Htmp(n) = Htmp(n)  &                   
                  +       m11*(dpy(i,j,k-4,n)*q(i,j,k-4,qpres)-dpy(i,j,k+3,n)*q(i,j,k+3,qpres)) &
                  +       m12*(dpy(i,j,k-4,n)*q(i,j,k-3,qpres)-dpy(i,j,k+3,n)*q(i,j,k+2,qpres)) &
                  +       m13*(dpy(i,j,k-4,n)*q(i,j,k-2,qpres)-dpy(i,j,k+3,n)*q(i,j,k+1,qpres)) &
                  +       m14*(dpy(i,j,k-4,n)*q(i,j,k-1,qpres)-dpy(i,j,k+3,n)*q(i,j,k  ,qpres)) &
                  +       m15*(dpy(i,j,k-4,n)*q(i,j,k  ,qpres)-dpy(i,j,k+3,n)*q(i,j,k-1,qpres)) &
                  &     + m21*(dpy(i,j,k-3,n)*q(i,j,k-4,qpres)-dpy(i,j,k+2,n)*q(i,j,k+3,qpres)) &
                  +       m22*(dpy(i,j,k-3,n)*q(i,j,k-3,qpres)-dpy(i,j,k+2,n)*q(i,j,k+2,qpres)) &
                  +       m23*(dpy(i,j,k-3,n)*q(i,j,k-2,qpres)-dpy(i,j,k+2,n)*q(i,j,k+1,qpres)) &
                  +       m24*(dpy(i,j,k-3,n)*q(i,j,k-1,qpres)-dpy(i,j,k+2,n)*q(i,j,k  ,qpres)) &
                  +       m25*(dpy(i,j,k-3,n)*q(i,j,k  ,qpres)-dpy(i,j,k+2,n)*q(i,j,k-1,qpres)) &
                  +       m26*(dpy(i,j,k-3,n)*q(i,j,k+1,qpres)-dpy(i,j,k+2,n)*q(i,j,k-2,qpres)) &
                  &     + m31*(dpy(i,j,k-2,n)*q(i,j,k-4,qpres)-dpy(i,j,k+1,n)*q(i,j,k+3,qpres)) &
                  +       m32*(dpy(i,j,k-2,n)*q(i,j,k-3,qpres)-dpy(i,j,k+1,n)*q(i,j,k+2,qpres)) &
                  +       m33*(dpy(i,j,k-2,n)*q(i,j,k-2,qpres)-dpy(i,j,k+1,n)*q(i,j,k+1,qpres)) &
                  +       m34*(dpy(i,j,k-2,n)*q(i,j,k-1,qpres)-dpy(i,j,k+1,n)*q(i,j,k  ,qpres)) &
                  +       m35*(dpy(i,j,k-2,n)*q(i,j,k  ,qpres)-dpy(i,j,k+1,n)*q(i,j,k-1,qpres)) &
                  +       m36*(dpy(i,j,k-2,n)*q(i,j,k+1,qpres)-dpy(i,j,k+1,n)*q(i,j,k-2,qpres)) &
                  +       m37*(dpy(i,j,k-2,n)*q(i,j,k+2,qpres)-dpy(i,j,k+1,n)*q(i,j,k-3,qpres)) &
                  &     + m41*(dpy(i,j,k-1,n)*q(i,j,k-4,qpres)-dpy(i,j,k  ,n)*q(i,j,k+3,qpres)) &
                  +       m42*(dpy(i,j,k-1,n)*q(i,j,k-3,qpres)-dpy(i,j,k  ,n)*q(i,j,k+2,qpres)) &
                  +       m43*(dpy(i,j,k-1,n)*q(i,j,k-2,qpres)-dpy(i,j,k  ,n)*q(i,j,k+1,qpres)) &
                  +       m44*(dpy(i,j,k-1,n)*q(i,j,k-1,qpres)-dpy(i,j,k  ,n)*q(i,j,k  ,qpres)) &
                  +       m45*(dpy(i,j,k-1,n)*q(i,j,k  ,qpres)-dpy(i,j,k  ,n)*q(i,j,k-1,qpres)) &
                  +       m46*(dpy(i,j,k-1,n)*q(i,j,k+1,qpres)-dpy(i,j,k  ,n)*q(i,j,k-2,qpres)) &
                  +       m47*(dpy(i,j,k-1,n)*q(i,j,k+2,qpres)-dpy(i,j,k  ,n)*q(i,j,k-3,qpres)) &
                  +       m48*(dpy(i,j,k-1,n)*q(i,j,k+3,qpres)-dpy(i,j,k  ,n)*q(i,j,k-4,qpres))
                Htot = Htot + Htmp(n)
                Ytmp(n) = (q(i,j,k-1,qyn) + q(i,j,k,qyn)) / 2.d0
             end do

             do n = 1, nspecies
                Hg(i,j,k,iry1+n-1) = Htmp(n) - Ytmp(n)*Htot
             end do

             do n = 1, nspecies
                qxn = qx1+n-1
                Hg(i,j,k,iene) =  Hg(i,j,k,iene) &
                  + m11*(dxe(i,j,k-4,n)*q(i,j,k-4,qxn)-dxe(i,j,k+3,n)*q(i,j,k+3,qxn)) &
                  + m12*(dxe(i,j,k-4,n)*q(i,j,k-3,qxn)-dxe(i,j,k+3,n)*q(i,j,k+2,qxn)) &
                  + m13*(dxe(i,j,k-4,n)*q(i,j,k-2,qxn)-dxe(i,j,k+3,n)*q(i,j,k+1,qxn)) &
                  + m14*(dxe(i,j,k-4,n)*q(i,j,k-1,qxn)-dxe(i,j,k+3,n)*q(i,j,k  ,qxn)) &
                  + m15*(dxe(i,j,k-4,n)*q(i,j,k  ,qxn)-dxe(i,j,k+3,n)*q(i,j,k-1,qxn)) &
                  + m21*(dxe(i,j,k-3,n)*q(i,j,k-4,qxn)-dxe(i,j,k+2,n)*q(i,j,k+3,qxn)) &
                  + m22*(dxe(i,j,k-3,n)*q(i,j,k-3,qxn)-dxe(i,j,k+2,n)*q(i,j,k+2,qxn)) &
                  + m23*(dxe(i,j,k-3,n)*q(i,j,k-2,qxn)-dxe(i,j,k+2,n)*q(i,j,k+1,qxn)) &
                  + m24*(dxe(i,j,k-3,n)*q(i,j,k-1,qxn)-dxe(i,j,k+2,n)*q(i,j,k  ,qxn)) &
                  + m25*(dxe(i,j,k-3,n)*q(i,j,k  ,qxn)-dxe(i,j,k+2,n)*q(i,j,k-1,qxn)) &
                  + m26*(dxe(i,j,k-3,n)*q(i,j,k+1,qxn)-dxe(i,j,k+2,n)*q(i,j,k-2,qxn)) &
                  + m31*(dxe(i,j,k-2,n)*q(i,j,k-4,qxn)-dxe(i,j,k+1,n)*q(i,j,k+3,qxn)) &
                  + m32*(dxe(i,j,k-2,n)*q(i,j,k-3,qxn)-dxe(i,j,k+1,n)*q(i,j,k+2,qxn)) &
                  + m33*(dxe(i,j,k-2,n)*q(i,j,k-2,qxn)-dxe(i,j,k+1,n)*q(i,j,k+1,qxn)) &
                  + m34*(dxe(i,j,k-2,n)*q(i,j,k-1,qxn)-dxe(i,j,k+1,n)*q(i,j,k  ,qxn)) &
                  + m35*(dxe(i,j,k-2,n)*q(i,j,k  ,qxn)-dxe(i,j,k+1,n)*q(i,j,k-1,qxn)) &
                  + m36*(dxe(i,j,k-2,n)*q(i,j,k+1,qxn)-dxe(i,j,k+1,n)*q(i,j,k-2,qxn)) &
                  + m37*(dxe(i,j,k-2,n)*q(i,j,k+2,qxn)-dxe(i,j,k+1,n)*q(i,j,k-3,qxn)) &
                  + m41*(dxe(i,j,k-1,n)*q(i,j,k-4,qxn)-dxe(i,j,k  ,n)*q(i,j,k+3,qxn)) &
                  + m42*(dxe(i,j,k-1,n)*q(i,j,k-3,qxn)-dxe(i,j,k  ,n)*q(i,j,k+2,qxn)) &
                  + m43*(dxe(i,j,k-1,n)*q(i,j,k-2,qxn)-dxe(i,j,k  ,n)*q(i,j,k+1,qxn)) &
                  + m44*(dxe(i,j,k-1,n)*q(i,j,k-1,qxn)-dxe(i,j,k  ,n)*q(i,j,k  ,qxn)) &
                  + m45*(dxe(i,j,k-1,n)*q(i,j,k  ,qxn)-dxe(i,j,k  ,n)*q(i,j,k-1,qxn)) &
                  + m46*(dxe(i,j,k-1,n)*q(i,j,k+1,qxn)-dxe(i,j,k  ,n)*q(i,j,k-2,qxn)) &
                  + m47*(dxe(i,j,k-1,n)*q(i,j,k+2,qxn)-dxe(i,j,k  ,n)*q(i,j,k-3,qxn)) &
                  + m48*(dxe(i,j,k-1,n)*q(i,j,k+3,qxn)-dxe(i,j,k  ,n)*q(i,j,k-4,qxn))
             end do

             Hg(i,j,k,iene) =  Hg(i,j,k,iene) &
                  + m11*(dpe(i,j,k-4)*q(i,j,k-4,qpres)-dpe(i,j,k+3)*q(i,j,k+3,qpres)) &
                  + m12*(dpe(i,j,k-4)*q(i,j,k-3,qpres)-dpe(i,j,k+3)*q(i,j,k+2,qpres)) &
                  + m13*(dpe(i,j,k-4)*q(i,j,k-2,qpres)-dpe(i,j,k+3)*q(i,j,k+1,qpres)) &
                  + m14*(dpe(i,j,k-4)*q(i,j,k-1,qpres)-dpe(i,j,k+3)*q(i,j,k  ,qpres)) &
                  + m15*(dpe(i,j,k-4)*q(i,j,k  ,qpres)-dpe(i,j,k+3)*q(i,j,k-1,qpres)) &
                  + m21*(dpe(i,j,k-3)*q(i,j,k-4,qpres)-dpe(i,j,k+2)*q(i,j,k+3,qpres)) &
                  + m22*(dpe(i,j,k-3)*q(i,j,k-3,qpres)-dpe(i,j,k+2)*q(i,j,k+2,qpres)) &
                  + m23*(dpe(i,j,k-3)*q(i,j,k-2,qpres)-dpe(i,j,k+2)*q(i,j,k+1,qpres)) &
                  + m24*(dpe(i,j,k-3)*q(i,j,k-1,qpres)-dpe(i,j,k+2)*q(i,j,k  ,qpres)) &
                  + m25*(dpe(i,j,k-3)*q(i,j,k  ,qpres)-dpe(i,j,k+2)*q(i,j,k-1,qpres)) &
                  + m26*(dpe(i,j,k-3)*q(i,j,k+1,qpres)-dpe(i,j,k+2)*q(i,j,k-2,qpres)) &
                  + m31*(dpe(i,j,k-2)*q(i,j,k-4,qpres)-dpe(i,j,k+1)*q(i,j,k+3,qpres)) &
                  + m32*(dpe(i,j,k-2)*q(i,j,k-3,qpres)-dpe(i,j,k+1)*q(i,j,k+2,qpres)) &
                  + m33*(dpe(i,j,k-2)*q(i,j,k-2,qpres)-dpe(i,j,k+1)*q(i,j,k+1,qpres)) &
                  + m34*(dpe(i,j,k-2)*q(i,j,k-1,qpres)-dpe(i,j,k+1)*q(i,j,k  ,qpres)) &
                  + m35*(dpe(i,j,k-2)*q(i,j,k  ,qpres)-dpe(i,j,k+1)*q(i,j,k-1,qpres)) &
                  + m36*(dpe(i,j,k-2)*q(i,j,k+1,qpres)-dpe(i,j,k+1)*q(i,j,k-2,qpres)) &
                  + m37*(dpe(i,j,k-2)*q(i,j,k+2,qpres)-dpe(i,j,k+1)*q(i,j,k-3,qpres)) &
                  + m41*(dpe(i,j,k-1)*q(i,j,k-4,qpres)-dpe(i,j,k  )*q(i,j,k+3,qpres)) &
                  + m42*(dpe(i,j,k-1)*q(i,j,k-3,qpres)-dpe(i,j,k  )*q(i,j,k+2,qpres)) &
                  + m43*(dpe(i,j,k-1)*q(i,j,k-2,qpres)-dpe(i,j,k  )*q(i,j,k+1,qpres)) &
                  + m44*(dpe(i,j,k-1)*q(i,j,k-1,qpres)-dpe(i,j,k  )*q(i,j,k  ,qpres)) &
                  + m45*(dpe(i,j,k-1)*q(i,j,k  ,qpres)-dpe(i,j,k  )*q(i,j,k-1,qpres)) &
                  + m46*(dpe(i,j,k-1)*q(i,j,k+1,qpres)-dpe(i,j,k  )*q(i,j,k-2,qpres)) &
                  + m47*(dpe(i,j,k-1)*q(i,j,k+2,qpres)-dpe(i,j,k  )*q(i,j,k-3,qpres)) &
                  + m48*(dpe(i,j,k-1)*q(i,j,k+3,qpres)-dpe(i,j,k  )*q(i,j,k-4,qpres))

             do n = 1, nspecies
                qhn = qh1+n-1
                hhalf = (q(i,j,k-1,qhn) + q(i,j,k,qhn)) / 2.d0
                Hg(i,j,k,iene) =  Hg(i,j,k,iene) - Ytmp(n) * hhalf * Htot
             end do

          end do
       end do
    end do
    !$omp end do

    ! add z-direction flux
    do n=2,ncons
       !$omp do
       do k=lo(3),hi(3)
          do j=lo(2),hi(2)
             do i=lo(1),hi(1)
                flx(i,j,k,n) = flx(i,j,k,n) + (Hg(i,j,k+1,n) - Hg(i,j,k,n)) * dx2inv(3)
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
             flx(i,j,k,iene) = flx(i,j,k,iene) &
                  + flx(i,j,k,imx)*q(i,j,k,qu) &
                  + flx(i,j,k,imy)*q(i,j,k,qv) &
                  + flx(i,j,k,imz)*q(i,j,k,qw)
          end do
       end do
    end do
    !$omp end do 
    
    !$omp end parallel

    deallocate(ux,uy,uz,vx,vy,vz,wx,wy,wz,vsp,vsm,Hg,dpy,dxe,dpe)

  end subroutine compact_diffterm_3d


  subroutine chemterm_3d(lo,hi,ng,q,up) ! up is UPrime that has no ghost cells
    integer,          intent(in ) :: lo(3),hi(3),ng
    double precision, intent(in )   :: q (-ng+lo(1):hi(1)+ng,-ng+lo(2):hi(2)+ng,-ng+lo(3):hi(3)+ng,nprim)
    double precision, intent(inout) :: up(    lo(1):hi(1)   ,    lo(2):hi(2)   ,    lo(3):hi(3)   ,ncons)

    integer :: iwrk, i,j,k
    double precision :: Xt(nspecies), wdot(nspecies), rwrk

    !$omp parallel do private(i,j,k,iwrk,rwrk,Xt,wdot)
    do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)

             Xt = q(i,j,k,qx1:qx1+nspecies-1)
             call ckwxr(q(i,j,k,qrho), q(i,j,k,qtemp), Xt, iwrk, rwrk, wdot)
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
    double precision :: dxinv(3), c, rwrk, Ru, Ruc, Pa, Cv, Cp
    double precision :: Tt, X(nspecies), gamma
    double precision :: courx, coury, courz

    do i=1,3
       dxinv(i) = 1.0d0 / dx(i)
    end do

    call ckrp(iwrk, rwrk, Ru, Ruc, Pa)

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

  subroutine S3D_diffterm_1(lo,hi,ng,ndq,dx,q,flx,mu,xi,qx,qy,qz)
 
    integer,          intent(in ) :: lo(3),hi(3),ng,ndq
    double precision, intent(in ) :: dx(3)
    double precision, intent(in ) :: q  (-ng+lo(1):hi(1)+ng,-ng+lo(2):hi(2)+ng,-ng+lo(3):hi(3)+ng,nprim)
    double precision, intent(in ) :: mu (-ng+lo(1):hi(1)+ng,-ng+lo(2):hi(2)+ng,-ng+lo(3):hi(3)+ng)
    double precision, intent(in ) :: xi (-ng+lo(1):hi(1)+ng,-ng+lo(2):hi(2)+ng,-ng+lo(3):hi(3)+ng)
    double precision, intent(out) :: qx (-ng+lo(1):hi(1)+ng,-ng+lo(2):hi(2)+ng,-ng+lo(3):hi(3)+ng,ndq)
    double precision, intent(out) :: qy (-ng+lo(1):hi(1)+ng,-ng+lo(2):hi(2)+ng,-ng+lo(3):hi(3)+ng,ndq)
    double precision, intent(out) :: qz (-ng+lo(1):hi(1)+ng,-ng+lo(2):hi(2)+ng,-ng+lo(3):hi(3)+ng,ndq)
    double precision, intent(out) :: flx(    lo(1):hi(1)   ,    lo(2):hi(2)   ,    lo(3):hi(3)   ,ncons)

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
    flx(:,:,:,irho) = 0.d0
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

             qx(i,j,k,idu)= &
                   (ALP*(q(i+1,j,k,qu)-q(i-1,j,k,qu)) &
                  + BET*(q(i+2,j,k,qu)-q(i-2,j,k,qu)) &
                  + GAM*(q(i+3,j,k,qu)-q(i-3,j,k,qu)) &
                  + DEL*(q(i+4,j,k,qu)-q(i-4,j,k,qu)))*dxinv(1)

             qx(i,j,k,idv)= &
                   (ALP*(q(i+1,j,k,qv)-q(i-1,j,k,qv)) &
                  + BET*(q(i+2,j,k,qv)-q(i-2,j,k,qv)) &
                  + GAM*(q(i+3,j,k,qv)-q(i-3,j,k,qv)) &
                  + DEL*(q(i+4,j,k,qv)-q(i-4,j,k,qv)))*dxinv(1)

             qx(i,j,k,idw)= &
                   (ALP*(q(i+1,j,k,qw)-q(i-1,j,k,qw)) &
                  + BET*(q(i+2,j,k,qw)-q(i-2,j,k,qw)) &
                  + GAM*(q(i+3,j,k,qw)-q(i-3,j,k,qw)) &
                  + DEL*(q(i+4,j,k,qw)-q(i-4,j,k,qw)))*dxinv(1)

          enddo
       enddo
    enddo
    !$OMP END DO NOWAIT

    !$OMP DO
    do k=lo(3)-ng,hi(3)+ng
       do j=lo(2),hi(2)   
          do i=lo(1)-ng,hi(1)+ng

             qy(i,j,k,idu)= &
                   (ALP*(q(i,j+1,k,qu)-q(i,j-1,k,qu)) &
                  + BET*(q(i,j+2,k,qu)-q(i,j-2,k,qu)) &
                  + GAM*(q(i,j+3,k,qu)-q(i,j-3,k,qu)) &
                  + DEL*(q(i,j+4,k,qu)-q(i,j-4,k,qu)))*dxinv(2)

             qy(i,j,k,idv)= &
                   (ALP*(q(i,j+1,k,qv)-q(i,j-1,k,qv)) &
                  + BET*(q(i,j+2,k,qv)-q(i,j-2,k,qv)) &
                  + GAM*(q(i,j+3,k,qv)-q(i,j-3,k,qv)) &
                  + DEL*(q(i,j+4,k,qv)-q(i,j-4,k,qv)))*dxinv(2)

             qy(i,j,k,idw)= &
                   (ALP*(q(i,j+1,k,qw)-q(i,j-1,k,qw)) &
                  + BET*(q(i,j+2,k,qw)-q(i,j-2,k,qw)) &
                  + GAM*(q(i,j+3,k,qw)-q(i,j-3,k,qw)) &
                  + DEL*(q(i,j+4,k,qw)-q(i,j-4,k,qw)))*dxinv(2)
          enddo
       enddo
    enddo
    !$OMP END DO NOWAIT

    !$OMP DO
    do k=lo(3),hi(3)
       do j=lo(2)-ng,hi(2)+ng
          do i=lo(1)-ng,hi(1)+ng

             qz(i,j,k,idu)= &
                   (ALP*(q(i,j,k+1,qu)-q(i,j,k-1,qu)) &
                  + BET*(q(i,j,k+2,qu)-q(i,j,k-2,qu)) &
                  + GAM*(q(i,j,k+3,qu)-q(i,j,k-3,qu)) &
                  + DEL*(q(i,j,k+4,qu)-q(i,j,k-4,qu)))*dxinv(3)

             qz(i,j,k,idv)= &
                   (ALP*(q(i,j,k+1,qv)-q(i,j,k-1,qv)) &
                  + BET*(q(i,j,k+2,qv)-q(i,j,k-2,qv)) &
                  + GAM*(q(i,j,k+3,qv)-q(i,j,k-3,qv)) &
                  + DEL*(q(i,j,k+4,qv)-q(i,j,k-4,qv)))*dxinv(3)

             qz(i,j,k,idw)= &
                   (ALP*(q(i,j,k+1,qw)-q(i,j,k-1,qw)) &
                  + BET*(q(i,j,k+2,qw)-q(i,j,k-2,qw)) &
                  + GAM*(q(i,j,k+3,qw)-q(i,j,k-3,qw)) &
                  + DEL*(q(i,j,k+4,qw)-q(i,j,k-4,qw)))*dxinv(3)
          enddo
       enddo
    enddo
    !$OMP END DO


    !$omp do
    do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)

             ! d(mu*dv/dx)/dy
             dmvxdy = (ALP*(mu(i,j+1,k)*qx(i,j+1,k,idv)-mu(i,j-1,k)*qx(i,j-1,k,idv)) &
                  +    BET*(mu(i,j+2,k)*qx(i,j+2,k,idv)-mu(i,j-2,k)*qx(i,j-2,k,idv)) &
                  +    GAM*(mu(i,j+3,k)*qx(i,j+3,k,idv)-mu(i,j-3,k)*qx(i,j-3,k,idv)) &
                  +    DEL*(mu(i,j+4,k)*qx(i,j+4,k,idv)-mu(i,j-4,k)*qx(i,j-4,k,idv)))*dxinv(2) 

             ! d(mu*dw/dx)/dz
             dmwxdz = (ALP*(mu(i,j,k+1)*qx(i,j,k+1,idw)-mu(i,j,k-1)*qx(i,j,k-1,idw)) &
                  +    BET*(mu(i,j,k+2)*qx(i,j,k+2,idw)-mu(i,j,k-2)*qx(i,j,k-2,idw)) &
                  +    GAM*(mu(i,j,k+3)*qx(i,j,k+3,idw)-mu(i,j,k-3)*qx(i,j,k-3,idw)) &
                  +    DEL*(mu(i,j,k+4)*qx(i,j,k+4,idw)-mu(i,j,k-4)*qx(i,j,k-4,idw)))*dxinv(3) 

             ! d((xi-2/3*mu)*(vy+wz))/dx
             dmvywzdx = (ALP*(vsm(i+1,j,k)*(qy(i+1,j,k,idv)+qz(i+1,j,k,idw))-vsm(i-1,j,k)*(qy(i-1,j,k,idv)+qz(i-1,j,k,idw))) &
                  +      BET*(vsm(i+2,j,k)*(qy(i+2,j,k,idv)+qz(i+2,j,k,idw))-vsm(i-2,j,k)*(qy(i-2,j,k,idv)+qz(i-2,j,k,idw))) &
                  +      GAM*(vsm(i+3,j,k)*(qy(i+3,j,k,idv)+qz(i+3,j,k,idw))-vsm(i-3,j,k)*(qy(i-3,j,k,idv)+qz(i-3,j,k,idw))) &
                  +      DEL*(vsm(i+4,j,k)*(qy(i+4,j,k,idv)+qz(i+4,j,k,idw))-vsm(i-4,j,k)*(qy(i-4,j,k,idv)+qz(i-4,j,k,idw))) &
                  ) * dxinv(1)

             ! d(mu*du/dy)/dx
             dmuydx = (ALP*(mu(i+1,j,k)*qy(i+1,j,k,idu)-mu(i-1,j,k)*qy(i-1,j,k,idu)) &
                  +    BET*(mu(i+2,j,k)*qy(i+2,j,k,idu)-mu(i-2,j,k)*qy(i-2,j,k,idu)) &
                  +    GAM*(mu(i+3,j,k)*qy(i+3,j,k,idu)-mu(i-3,j,k)*qy(i-3,j,k,idu)) &
                  +    DEL*(mu(i+4,j,k)*qy(i+4,j,k,idu)-mu(i-4,j,k)*qy(i-4,j,k,idu)))*dxinv(1) 

             ! d(mu*dw/dy)/dz
             dmwydz = (ALP*(mu(i,j,k+1)*qy(i,j,k+1,idw)-mu(i,j,k-1)*qy(i,j,k-1,idw)) &
                  +    BET*(mu(i,j,k+2)*qy(i,j,k+2,idw)-mu(i,j,k-2)*qy(i,j,k-2,idw)) &
                  +    GAM*(mu(i,j,k+3)*qy(i,j,k+3,idw)-mu(i,j,k-3)*qy(i,j,k-3,idw)) &
                  +    DEL*(mu(i,j,k+4)*qy(i,j,k+4,idw)-mu(i,j,k-4)*qy(i,j,k-4,idw)))*dxinv(3) 

             ! d((xi-2/3*mu)*(ux+wz))/dy
             dmuxwzdy = (ALP*(vsm(i,j+1,k)*(qx(i,j+1,k,idu)+qz(i,j+1,k,idw))-vsm(i,j-1,k)*(qx(i,j-1,k,idu)+qz(i,j-1,k,idw))) &
                  +      BET*(vsm(i,j+2,k)*(qx(i,j+2,k,idu)+qz(i,j+2,k,idw))-vsm(i,j-2,k)*(qx(i,j-2,k,idu)+qz(i,j-2,k,idw))) &
                  +      GAM*(vsm(i,j+3,k)*(qx(i,j+3,k,idu)+qz(i,j+3,k,idw))-vsm(i,j-3,k)*(qx(i,j-3,k,idu)+qz(i,j-3,k,idw))) &
                  +      DEL*(vsm(i,j+4,k)*(qx(i,j+4,k,idu)+qz(i,j+4,k,idw))-vsm(i,j-4,k)*(qx(i,j-4,k,idu)+qz(i,j-4,k,idw))) &
                  ) * dxinv(2)

             ! d(mu*du/dz)/dx
             dmuzdx = (ALP*(mu(i+1,j,k)*qz(i+1,j,k,idu)-mu(i-1,j,k)*qz(i-1,j,k,idu)) &
                  +    BET*(mu(i+2,j,k)*qz(i+2,j,k,idu)-mu(i-2,j,k)*qz(i-2,j,k,idu)) &
                  +    GAM*(mu(i+3,j,k)*qz(i+3,j,k,idu)-mu(i-3,j,k)*qz(i-3,j,k,idu)) &
                  +    DEL*(mu(i+4,j,k)*qz(i+4,j,k,idu)-mu(i-4,j,k)*qz(i-4,j,k,idu)))*dxinv(1) 

             ! d(mu*dv/dz)/dy
             dmvzdy = (ALP*(mu(i,j+1,k)*qz(i,j+1,k,idv)-mu(i,j-1,k)*qz(i,j-1,k,idv)) &
                  +    BET*(mu(i,j+2,k)*qz(i,j+2,k,idv)-mu(i,j-2,k)*qz(i,j-2,k,idv)) &
                  +    GAM*(mu(i,j+3,k)*qz(i,j+3,k,idv)-mu(i,j-3,k)*qz(i,j-3,k,idv)) &
                  +    DEL*(mu(i,j+4,k)*qz(i,j+4,k,idv)-mu(i,j-4,k)*qz(i,j-4,k,idv)))*dxinv(2) 

             ! d((xi-2/3*mu)*(ux+vy))/dz
             dmuxvydz = (ALP*(vsm(i,j,k+1)*(qx(i,j,k+1,idu)+qy(i,j,k+1,idv))-vsm(i,j,k-1)*(qx(i,j,k-1,idu)+qy(i,j,k-1,idv))) &
                  +      BET*(vsm(i,j,k+2)*(qx(i,j,k+2,idu)+qy(i,j,k+2,idv))-vsm(i,j,k-2)*(qx(i,j,k-2,idu)+qy(i,j,k-2,idv))) &
                  +      GAM*(vsm(i,j,k+3)*(qx(i,j,k+3,idu)+qy(i,j,k+3,idv))-vsm(i,j,k-3)*(qx(i,j,k-3,idu)+qy(i,j,k-3,idv))) &
                  +      DEL*(vsm(i,j,k+4)*(qx(i,j,k+4,idu)+qy(i,j,k+4,idv))-vsm(i,j,k-4)*(qx(i,j,k-4,idu)+qy(i,j,k-4,idv))) &
                  ) * dxinv(3)

             flx(i,j,k,imx) = dmvxdy + dmwxdz + dmvywzdx
             flx(i,j,k,imy) = dmuydx + dmwydz + dmuxwzdy
             flx(i,j,k,imz) = dmuzdx + dmvzdy + dmuxvydz

             divu = (qx(i,j,k,idu)+qy(i,j,k,idv)+qz(i,j,k,idw))*vsm(i,j,k)
             tauxx = 2.d0*mu(i,j,k)*qx(i,j,k,idu) + divu
             tauyy = 2.d0*mu(i,j,k)*qy(i,j,k,idv) + divu
             tauzz = 2.d0*mu(i,j,k)*qz(i,j,k,idw) + divu
             
             ! change in internal energy
             flx(i,j,k,iene) = tauxx*qx(i,j,k,idu) + tauyy*qy(i,j,k,idv) + tauzz*qz(i,j,k,idw) &
                  + mu(i,j,k)*((qy(i,j,k,idu)+qx(i,j,k,idv))**2 &
                  &          + (qx(i,j,k,idw)+qz(i,j,k,idu))**2 &
                  &          + (qz(i,j,k,idv)+qy(i,j,k,idw))**2 )

          end do
       end do
    end do
    !$omp end do nowait

    !$omp workshare
    flx(:,:,:,iry1:) = 0.d0
    !$omp end workshare

    !$omp do
    do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             qx(i,j,k,idT) = (ALP*(q(i+1,j,k,qtemp)-q(i-1,j,k,qtemp)) &
                  +           BET*(q(i+2,j,k,qtemp)-q(i-2,j,k,qtemp)) &
                  +           GAM*(q(i+3,j,k,qtemp)-q(i-3,j,k,qtemp)) &
                  +           DEL*(q(i+4,j,k,qtemp)-q(i-4,j,k,qtemp)))*dxinv(1)

             qx(i,j,k,idp) = (ALP*(q(i+1,j,k,qpres)-q(i-1,j,k,qpres)) &
                  +           BET*(q(i+2,j,k,qpres)-q(i-2,j,k,qpres)) &
                  +           GAM*(q(i+3,j,k,qpres)-q(i-3,j,k,qpres)) &
                  +           DEL*(q(i+4,j,k,qpres)-q(i-4,j,k,qpres)))*dxinv(1)

             qy(i,j,k,idT) = (ALP*(q(i,j+1,k,qtemp)-q(i,j-1,k,qtemp)) &
                  +           BET*(q(i,j+2,k,qtemp)-q(i,j-2,k,qtemp)) &
                  +           GAM*(q(i,j+3,k,qtemp)-q(i,j-3,k,qtemp)) &
                  +           DEL*(q(i,j+4,k,qtemp)-q(i,j-4,k,qtemp)))*dxinv(2)

             qy(i,j,k,idp) = (ALP*(q(i,j+1,k,qpres)-q(i,j-1,k,qpres)) &
                  +           BET*(q(i,j+2,k,qpres)-q(i,j-2,k,qpres)) &
                  +           GAM*(q(i,j+3,k,qpres)-q(i,j-3,k,qpres)) &
                  +           DEL*(q(i,j+4,k,qpres)-q(i,j-4,k,qpres)))*dxinv(2)

             qz(i,j,k,idT) = (ALP*(q(i,j,k+1,qtemp)-q(i,j,k-1,qtemp)) &
                  +           BET*(q(i,j,k+2,qtemp)-q(i,j,k-2,qtemp)) &
                  +           GAM*(q(i,j,k+3,qtemp)-q(i,j,k-3,qtemp)) &
                  +           DEL*(q(i,j,k+4,qtemp)-q(i,j,k-4,qtemp)))*dxinv(3)

             qz(i,j,k,idp) = (ALP*(q(i,j,k+1,qpres)-q(i,j,k-1,qpres)) &
                  +           BET*(q(i,j,k+2,qpres)-q(i,j,k-2,qpres)) &
                  +           GAM*(q(i,j,k+3,qpres)-q(i,j,k-3,qpres)) &
                  +           DEL*(q(i,j,k+4,qpres)-q(i,j,k-4,qpres)))*dxinv(3)
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
                qx(i,j,k,qdxn) = (ALP*(q(i+1,j,k,qxn)-q(i-1,j,k,qxn)) &
                     +            BET*(q(i+2,j,k,qxn)-q(i-2,j,k,qxn)) &
                     +            GAM*(q(i+3,j,k,qxn)-q(i-3,j,k,qxn)) &
                     +            DEL*(q(i+4,j,k,qxn)-q(i-4,j,k,qxn)))*dxinv(1)

                qy(i,j,k,qdxn) = (ALP*(q(i,j+1,k,qxn)-q(i,j-1,k,qxn)) &
                     +            BET*(q(i,j+2,k,qxn)-q(i,j-2,k,qxn)) &
                     +            GAM*(q(i,j+3,k,qxn)-q(i,j-3,k,qxn)) &
                     +            DEL*(q(i,j+4,k,qxn)-q(i,j-4,k,qxn)))*dxinv(2)

                qz(i,j,k,qdxn) = (ALP*(q(i,j,k+1,qxn)-q(i,j,k-1,qxn)) &
                     +            BET*(q(i,j,k+2,qxn)-q(i,j,k-2,qxn)) &
                     +            GAM*(q(i,j,k+3,qxn)-q(i,j,k-3,qxn)) &
                     +            DEL*(q(i,j,k+4,qxn)-q(i,j,k-4,qxn)))*dxinv(3)
             enddo
          enddo
       enddo
       !$omp end do nowait
    enddo

    !$omp end parallel

    deallocate(vsm)

  end subroutine S3D_diffterm_1


  subroutine S3D_diffterm_2(lo,hi,ng,ndq,dx,q,flx,mu,xi,lam,dxy,qx,qy,qz)

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
    double precision, intent(inout):: flx(    lo(1):hi(1)   ,    lo(2):hi(2)   ,    lo(3):hi(3)   ,ncons)
 
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
             flx(i,j,k,imx) = flx(i,j,k,imx) &
                  + (ALP*(vp(i+1,j,k)*qx(i+1,j,k,idu)-vp(i-1,j,k)*qx(i-1,j,k,idu)) &
                  +  BET*(vp(i+2,j,k)*qx(i+2,j,k,idu)-vp(i-2,j,k)*qx(i-2,j,k,idu)) &
                  +  GAM*(vp(i+3,j,k)*qx(i+3,j,k,idu)-vp(i-3,j,k)*qx(i-3,j,k,idu)) &
                  +  DEL*(vp(i+4,j,k)*qx(i+4,j,k,idu)-vp(i-4,j,k)*qx(i-4,j,k,idu)))*dxinv(1)&
                  + (ALP*(mu(i,j+1,k)*qy(i,j+1,k,idu)-mu(i,j-1,k)*qy(i,j-1,k,idu)) &
                  +  BET*(mu(i,j+2,k)*qy(i,j+2,k,idu)-mu(i,j-2,k)*qy(i,j-2,k,idu)) &
                  +  GAM*(mu(i,j+3,k)*qy(i,j+3,k,idu)-mu(i,j-3,k)*qy(i,j-3,k,idu)) &
                  +  DEL*(mu(i,j+4,k)*qy(i,j+4,k,idu)-mu(i,j-4,k)*qy(i,j-4,k,idu)))*dxinv(2)&
                  + (ALP*(mu(i,j,k+1)*qz(i,j,k+1,idu)-mu(i,j,k-1)*qz(i,j,k-1,idu)) &
                  +  BET*(mu(i,j,k+2)*qz(i,j,k+2,idu)-mu(i,j,k-2)*qz(i,j,k-2,idu)) &
                  +  GAM*(mu(i,j,k+3)*qz(i,j,k+3,idu)-mu(i,j,k-3)*qz(i,j,k-3,idu)) &
                  +  DEL*(mu(i,j,k+4)*qz(i,j,k+4,idu)-mu(i,j,k-4)*qz(i,j,k-4,idu)))*dxinv(3)
          end do
       end do
    end do
    !$omp end do nowait
    
    ! ===== my =====
    !$omp do
    do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             flx(i,j,k,imy) = flx(i,j,k,imy) &
                  + (ALP*(mu(i+1,j,k)*qx(i+1,j,k,idv)-mu(i-1,j,k)*qx(i-1,j,k,idv)) &
                  +  BET*(mu(i+2,j,k)*qx(i+2,j,k,idv)-mu(i-2,j,k)*qx(i-2,j,k,idv)) &
                  +  GAM*(mu(i+3,j,k)*qx(i+3,j,k,idv)-mu(i-3,j,k)*qx(i-3,j,k,idv)) &
                  +  DEL*(mu(i+4,j,k)*qx(i+4,j,k,idv)-mu(i-4,j,k)*qx(i-4,j,k,idv)))*dxinv(1)&
                  + (ALP*(vp(i,j+1,k)*qy(i,j+1,k,idv)-vp(i,j-1,k)*qy(i,j-1,k,idv)) &
                  +  BET*(vp(i,j+2,k)*qy(i,j+2,k,idv)-vp(i,j-2,k)*qy(i,j-2,k,idv)) &
                  +  GAM*(vp(i,j+3,k)*qy(i,j+3,k,idv)-vp(i,j-3,k)*qy(i,j-3,k,idv)) &
                  +  DEL*(vp(i,j+4,k)*qy(i,j+4,k,idv)-vp(i,j-4,k)*qy(i,j-4,k,idv)))*dxinv(2)&
                  + (ALP*(mu(i,j,k+1)*qz(i,j,k+1,idv)-mu(i,j,k-1)*qz(i,j,k-1,idv)) &
                  +  BET*(mu(i,j,k+2)*qz(i,j,k+2,idv)-mu(i,j,k-2)*qz(i,j,k-2,idv)) &
                  +  GAM*(mu(i,j,k+3)*qz(i,j,k+3,idv)-mu(i,j,k-3)*qz(i,j,k-3,idv)) &
                  +  DEL*(mu(i,j,k+4)*qz(i,j,k+4,idv)-mu(i,j,k-4)*qz(i,j,k-4,idv)))*dxinv(3)
          end do
       end do
    end do
    !$omp end do nowait
    
    ! ===== mz =====
    !$omp do
    do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             flx(i,j,k,imz) = flx(i,j,k,imz) &
                  + (ALP*(mu(i+1,j,k)*qx(i+1,j,k,idw)-mu(i-1,j,k)*qx(i-1,j,k,idw)) &
                  +  BET*(mu(i+2,j,k)*qx(i+2,j,k,idw)-mu(i-2,j,k)*qx(i-2,j,k,idw)) &
                  +  GAM*(mu(i+3,j,k)*qx(i+3,j,k,idw)-mu(i-3,j,k)*qx(i-3,j,k,idw)) &
                  +  DEL*(mu(i+4,j,k)*qx(i+4,j,k,idw)-mu(i-4,j,k)*qx(i-4,j,k,idw)))*dxinv(1)&
                  + (ALP*(mu(i,j+1,k)*qy(i,j+1,k,idw)-mu(i,j-1,k)*qy(i,j-1,k,idw)) &
                  +  BET*(mu(i,j+2,k)*qy(i,j+2,k,idw)-mu(i,j-2,k)*qy(i,j-2,k,idw)) &
                  +  GAM*(mu(i,j+3,k)*qy(i,j+3,k,idw)-mu(i,j-3,k)*qy(i,j-3,k,idw)) &
                  +  DEL*(mu(i,j+4,k)*qy(i,j+4,k,idw)-mu(i,j-4,k)*qy(i,j-4,k,idw)))*dxinv(2)&
                  + (ALP*(vp(i,j,k+1)*qz(i,j,k+1,idw)-vp(i,j,k-1)*qz(i,j,k-1,idw)) &
                  +  BET*(vp(i,j,k+2)*qz(i,j,k+2,idw)-vp(i,j,k-2)*qz(i,j,k-2,idw)) &
                  +  GAM*(vp(i,j,k+3)*qz(i,j,k+3,idw)-vp(i,j,k-3)*qz(i,j,k-3,idw)) &
                  +  DEL*(vp(i,j,k+4)*qz(i,j,k+4,idw)-vp(i,j,k-4)*qz(i,j,k-4,idw)))*dxinv(3)
          end do
       end do
    end do
    !$omp end do
    
    ! add kinetic energy
    !$omp do
    do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             flx(i,j,k,iene) = flx(i,j,k,iene) &
                  + flx(i,j,k,imx)*q(i,j,k,qu) &
                  + flx(i,j,k,imy)*q(i,j,k,qv) &
                  + flx(i,j,k,imz)*q(i,j,k,qw)
          end do
       end do
    end do
    !$omp end do

    ! thermal conduction
    !$omp do
    do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             flx(i,j,k,iene) = flx(i,j,k,iene) &
                  + (ALP*(lam(i+1,j,k)*qx(i+1,j,k,idT)-lam(i-1,j,k)*qx(i-1,j,k,idT)) &
                  +  BET*(lam(i+2,j,k)*qx(i+2,j,k,idT)-lam(i-2,j,k)*qx(i-2,j,k,idT)) &
                  +  GAM*(lam(i+3,j,k)*qx(i+3,j,k,idT)-lam(i-3,j,k)*qx(i-3,j,k,idT)) &
                  +  DEL*(lam(i+4,j,k)*qx(i+4,j,k,idT)-lam(i-4,j,k)*qx(i-4,j,k,idT)))*dxinv(1)&
                  + (ALP*(lam(i,j+1,k)*qy(i,j+1,k,idT)-lam(i,j-1,k)*qy(i,j-1,k,idT)) &
                  +  BET*(lam(i,j+2,k)*qy(i,j+2,k,idT)-lam(i,j-2,k)*qy(i,j-2,k,idT)) &
                  +  GAM*(lam(i,j+3,k)*qy(i,j+3,k,idT)-lam(i,j-3,k)*qy(i,j-3,k,idT)) &
                  +  DEL*(lam(i,j+4,k)*qy(i,j+4,k,idT)-lam(i,j-4,k)*qy(i,j-4,k,idT)))*dxinv(2)&
                  + (ALP*(lam(i,j,k+1)*qz(i,j,k+1,idT)-lam(i,j,k-1)*qz(i,j,k-1,idT)) &
                  +  BET*(lam(i,j,k+2)*qz(i,j,k+2,idT)-lam(i,j,k-2)*qz(i,j,k-2,idT)) &
                  +  GAM*(lam(i,j,k+3)*qz(i,j,k+3,idT)-lam(i,j,k-3)*qz(i,j,k-3,idT)) &
                  +  DEL*(lam(i,j,k+4)*qz(i,j,k+4,idT)-lam(i,j,k-4)*qz(i,j,k-4,idT)))*dxinv(3)
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
                flx(i,j,k,iryn) = flx(i,j,k,iryn) + &
                     ( ALP*(FY(i+1,j,k,n)-FY(i-1,j,k,n)) &
                     + BET*(FY(i+2,j,k,n)-FY(i-2,j,k,n)) &
                     + GAM*(FY(i+3,j,k,n)-FY(i-3,j,k,n)) &
                     + DEL*(FY(i+4,j,k,n)-FY(i-4,j,k,n)))*dxinv(1)
             end do
          end do
       end do
       !$omp end do nowait
    end do
    
    !$omp do
    do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             flx(i,j,k,iene) = flx(i,j,k,iene) + &
                  ( ALP*(FE(i+1,j,k)-FE(i-1,j,k)) &
                  + BET*(FE(i+2,j,k)-FE(i-2,j,k)) &
                  + GAM*(FE(i+3,j,k)-FE(i-3,j,k)) &
                  + DEL*(FE(i+4,j,k)-FE(i-4,j,k)))*dxinv(1)
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
                flx(i,j,k,iryn) = flx(i,j,k,iryn) + &
                     ( ALP*(FY(i,j+1,k,n)-FY(i,j-1,k,n)) &
                     + BET*(FY(i,j+2,k,n)-FY(i,j-2,k,n)) &
                     + GAM*(FY(i,j+3,k,n)-FY(i,j-3,k,n)) &
                     + DEL*(FY(i,j+4,k,n)-FY(i,j-4,k,n)))*dxinv(2)
             end do
          end do
       end do
       !$omp end do nowait
    end do
    
    !$omp do
    do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             flx(i,j,k,iene) = flx(i,j,k,iene) + &
                  ( ALP*(FE(i,j+1,k)-FE(i,j-1,k)) &
                  + BET*(FE(i,j+2,k)-FE(i,j-2,k)) &
                  + GAM*(FE(i,j+3,k)-FE(i,j-3,k)) &
                  + DEL*(FE(i,j+4,k)-FE(i,j-4,k)))*dxinv(2)
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
                flx(i,j,k,iryn) = flx(i,j,k,iryn) + &
                     ( ALP*(FY(i,j,k+1,n)-FY(i,j,k-1,n)) &
                     + BET*(FY(i,j,k+2,n)-FY(i,j,k-2,n)) &
                     + GAM*(FY(i,j,k+3,n)-FY(i,j,k-3,n)) &
                     + DEL*(FY(i,j,k+4,n)-FY(i,j,k-4,n)))*dxinv(3)
             end do
          end do
       end do
       !$omp end do nowait
    end do

    !$omp do
    do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             flx(i,j,k,iene) = flx(i,j,k,iene) + &
                  ( ALP*(FE(i,j,k+1)-FE(i,j,k-1)) &
                  + BET*(FE(i,j,k+2)-FE(i,j,k-2)) &
                  + GAM*(FE(i,j,k+3)-FE(i,j,k-3)) &
                  + DEL*(FE(i,j,k+4)-FE(i,j,k-4)))*dxinv(3)
          end do
       end do
    end do
    !$omp end do

    !$omp end parallel

    deallocate(vp,dpy,dpe,FY,FE)

  end subroutine S3D_diffterm_2

end module kernels_module
