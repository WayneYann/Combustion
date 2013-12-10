module difterm_module
  implicit none
  double precision, parameter :: twoThirds = 2.d0/3.d0
  private
  public :: difterm
contains

  subroutine difterm(lo,hi,U,Ulo,Uhi,fx,fxlo,fxhi,fy,fylo,fyhi,fz,fzlo,fzhi,dxinv)

    use meth_params_module, only : NVAR, URHO, UMX, UMY, UMZ, NSPEC, QCVAR
    use weno_module, only : cellavg2dergausspt_1d, cellavg2gausspt_1d, &
         cellavg2dergausspt_2d, cellavg2gausspt_2d
    use polyinterp_module, only : cc2zgauss_3d, cc2zface_3d
    use convert_3d_module, only : cellavg2cc_3d
    use variables_module, only : ctoprim
    use transport_properties, only : get_transport_properties

    integer, intent(in) :: lo(3), hi(3), Ulo(3), Uhi(3), fxlo(3), fxhi(3), &
         fylo(3), fyhi(3), fzlo(3), fzhi(3)
    double precision,intent(in   )::dxinv(3)
    double precision,intent(in   ):: U( Ulo(1): Uhi(1), Ulo(2): Uhi(2), Ulo(3): Uhi(3),NVAR)
    double precision,intent(inout)::fx(fxlo(1):fxhi(1),fxlo(2):fxhi(2),fxlo(3):fxhi(3),NVAR)
    double precision,intent(inout)::fy(fylo(1):fyhi(1),fylo(2):fyhi(2),fylo(3):fyhi(3),NVAR)
    double precision,intent(inout)::fz(fzlo(1):fzhi(1),fzlo(2):fzhi(2),fzlo(3):fzhi(3),NVAR)

    integer :: i, j, k, n 
    integer :: g2lo(3), g2hi(3), tlo(3), thi(3)
    integer :: Uzlo(3), Uzhi(3), dUdzlo(3), dUdzhi(3)
    integer :: Uxylo(3), Uxyhi(3), dUdxylo(3), dUdxyhi(3)
    double precision, allocatable, dimension(:,:,:,:) :: UZ1, UZ2, dUdz1, dUdz2
    double precision, allocatable, dimension(:,:,:,:) :: Uxy1, Uxy2, Uxy3, Uxy4
    double precision, allocatable, dimension(:,:,:,:) :: dUdx1, dUdx2, dUdx3, dUdx4
    double precision, allocatable, dimension(:,:,:,:) :: dUdy1, dUdy2, dUdy3, dUdy4
    double precision, allocatable :: Ucc(:,:,:,:)
    double precision, allocatable :: mucc(:,:,:),xicc(:,:,:),lamcc(:,:,:),Ddiacc(:,:,:,:)
    double precision, allocatable, dimension(:,:,:) :: muz1, muz2, xiz1, xiz2, lamz1, lamz2
    double precision, allocatable, dimension(:,:,:,:) :: Ddiaz1, Ddiaz2
    double precision, allocatable, dimension(:,:,:) ::   muxy1,  muxy2,  muxy3,  muxy4
    double precision, allocatable, dimension(:,:,:) ::   xixy1,  xixy2,  xixy3,  xixy4
    double precision, allocatable, dimension(:,:,:) ::  lamxy1, lamxy2, lamxy3, lamxy4
    double precision, allocatable, dimension(:,:,:,:) :: Ddiaxy1,Ddiaxy2,Ddiaxy3,Ddiaxy4

    g2lo = lo-2 ;    g2hi = hi+2

    allocate(Ucc   (g2lo(1):g2hi(1),g2lo(2):g2hi(2),g2lo(3):g2hi(3),QCVAR))
    allocate(mucc  (g2lo(1):g2hi(1),g2lo(2):g2hi(2),g2lo(3):g2hi(3)))
    allocate(xicc  (g2lo(1):g2hi(1),g2lo(2):g2hi(2),g2lo(3):g2hi(3)))
    allocate(lamcc (g2lo(1):g2hi(1),g2lo(2):g2hi(2),g2lo(3):g2hi(3)))
    allocate(Ddiacc(g2lo(1):g2hi(1),g2lo(2):g2hi(2),g2lo(3):g2hi(3),NSPEC))

    do n=1,NVAR
       call cellavg2cc_3d(g2lo,g2hi, U(:,:,:,n), Ulo,Uhi, Ucc(:,:,:,n), g2lo,g2hi)
    end do

    call ctoprim(g2lo, g2hi, Ucc, g2lo, g2hi, QCVAR)
    ! Ucc now contains Q at cell centers

    ! transport coefficients at cell centers
    call get_transport_properties(g2lo,g2hi, Ucc,g2lo,g2hi,QCVAR, &
         mucc,xicc,lamcc,Ddiacc, g2lo,g2hi)

    deallocate(Ucc)

    Uzlo(1) = lo(1)-3
    Uzlo(2) = lo(2)-3
    Uzlo(3) = lo(3)
    Uzhi(1) = hi(1)+3
    Uzhi(2) = hi(2)+3
    Uzhi(3) = hi(3)

    dUdzlo(1) = lo(1)-2
    dUdzlo(2) = lo(2)-2
    dUdzlo(3) = lo(3)
    dUdzhi(1) = hi(1)+2
    dUdzhi(2) = hi(2)+2
    dUdzhi(3) = hi(3)

    allocate(  UZ1(  Uzlo(1):  Uzhi(1),  Uzlo(2):  Uzhi(2),  Uzlo(3):  Uzhi(3),NVAR))
    allocate(  UZ2(  Uzlo(1):  Uzhi(1),  Uzlo(2):  Uzhi(2),  Uzlo(3):  Uzhi(3),NVAR))
    allocate(dUdz1(dUdzlo(1):dUdzhi(1),dUdzlo(2):dUdzhi(2),dUdzlo(3):dUdzhi(3),4))
    allocate(dUdz2(dUdzlo(1):dUdzhi(1),dUdzlo(2):dUdzhi(2),dUdzlo(3):dUdzhi(3),4))

    do n=1,4  ! rho, mx, my, & mz
       do j=lo(2)-2,hi(2)+2
          do i=lo(1)-2,hi(1)+2
             call cellavg2dergausspt_1d(lo(3),hi(3), U(i,j,:,n), Ulo(3), Uhi(3), &
                  dUdz1(i,j,:,n), dUdz2(i,j,:,n), lo(3), hi(3))
          end do
       end do       
    end do

    do n=1,NVAR
       do j=lo(2)-3,hi(2)+3
          do i=lo(1)-3,hi(1)+3
             call cellavg2gausspt_1d(lo(3),hi(3), U(i,j,:,n), Ulo(3), Uhi(3), &
                  UZ1(i,j,:,n), UZ2(i,j,:,n), lo(3), hi(3))
          end do
       end do
    end do

    tlo(1:2) = g2lo(1:2)
    tlo(3) = lo(3)
    thi(1:2) = g2hi(1:2)
    thi(3) = hi(3)

    allocate(  muz1(tlo(1):thi(1),tlo(2):thi(2),tlo(3):thi(3)))
    allocate(  xiz1(tlo(1):thi(1),tlo(2):thi(2),tlo(3):thi(3)))
    allocate( lamz1(tlo(1):thi(1),tlo(2):thi(2),tlo(3):thi(3)))
    allocate(Ddiaz1(tlo(1):thi(1),tlo(2):thi(2),tlo(3):thi(3),NSPEC))

    allocate(  muz2(tlo(1):thi(1),tlo(2):thi(2),tlo(3):thi(3)))
    allocate(  xiz2(tlo(1):thi(1),tlo(2):thi(2),tlo(3):thi(3)))
    allocate( lamz2(tlo(1):thi(1),tlo(2):thi(2),tlo(3):thi(3)))
    allocate(Ddiaz2(tlo(1):thi(1),tlo(2):thi(2),tlo(3):thi(3),NSPEC))

    call cc2zgauss_3d(tlo,thi,  mucc, g2lo,g2hi,  muz1,  muz2, tlo,thi)
    call cc2zgauss_3d(tlo,thi,  xicc, g2lo,g2hi,  xiz1,  xiz2, tlo,thi)
    call cc2zgauss_3d(tlo,thi, lamcc, g2lo,g2hi, lamz1, lamz2, tlo,thi)
    do n=1,NSPEC
       call cc2zgauss_3d(tlo,thi, Ddiacc(:,:,:,n), g2lo,g2hi, &
            Ddiaz1(:,:,:,n), Ddiaz2(:,:,:,n), tlo,thi)
    end do

    call diff_xy(lo,hi, UZ1,Uzlo,Uzhi, dUdz1,dUdzlo,dUdzhi, &
         muz1, xiz1, lamz1, Ddiaz1, tlo, thi, &
         fx, fxlo, fxhi, fy, fylo, fyhi, dxinv)

    call diff_xy(lo,hi, UZ2,Uzlo,Uzhi, dUdz2,dUdzlo,dUdzhi, &
         muz2, xiz2, lamz2, Ddiaz2, tlo, thi, &
         fx, fxlo, fxhi, fy, fylo, fyhi, dxinv)

    deallocate(UZ1,UZ2,dUdz1,dUdz2)
    deallocate(muz1,muz2,xiz1,xiz2,lamz1,lamz2,Ddiaz1,Ddiaz2)

    ! work on z-direction

    Uxylo(1:2) = lo(1:2)
    Uxylo(3) = lo(3)-3
    Uxyhi(1:2) = hi(1:2)
    Uxyhi(3) = hi(3)+3
    
    dUdxylo(1:2) = lo(1:2)
    dUdxylo(3) = lo(3)-2
    dUdxyhi(1:2) = hi(1:2)
    dUdxyhi(3) = hi(3)+2

    tlo = lo
    thi(1:2) = hi(1:2)
    thi(3) = hi(3)+1

    allocate( Uxy1(  Uxylo(1):  Uxyhi(1),  Uxylo(2):  Uxyhi(2),  Uxylo(3):  Uxyhi(3),NVAR))
    allocate( Uxy2(  Uxylo(1):  Uxyhi(1),  Uxylo(2):  Uxyhi(2),  Uxylo(3):  Uxyhi(3),NVAR))
    allocate( Uxy3(  Uxylo(1):  Uxyhi(1),  Uxylo(2):  Uxyhi(2),  Uxylo(3):  Uxyhi(3),NVAR))
    allocate( Uxy4(  Uxylo(1):  Uxyhi(1),  Uxylo(2):  Uxyhi(2),  Uxylo(3):  Uxyhi(3),NVAR))
    allocate(dUdx1(dUdxylo(1):dUdxyhi(1),dUdxylo(2):dUdxyhi(2),dUdxylo(3):dUdxyhi(3),3))
    allocate(dUdx2(dUdxylo(1):dUdxyhi(1),dUdxylo(2):dUdxyhi(2),dUdxylo(3):dUdxyhi(3),3))
    allocate(dUdx3(dUdxylo(1):dUdxyhi(1),dUdxylo(2):dUdxyhi(2),dUdxylo(3):dUdxyhi(3),3))
    allocate(dUdx4(dUdxylo(1):dUdxyhi(1),dUdxylo(2):dUdxyhi(2),dUdxylo(3):dUdxyhi(3),3))
    allocate(dUdy1(dUdxylo(1):dUdxyhi(1),dUdxylo(2):dUdxyhi(2),dUdxylo(3):dUdxyhi(3),3))
    allocate(dUdy2(dUdxylo(1):dUdxyhi(1),dUdxylo(2):dUdxyhi(2),dUdxylo(3):dUdxyhi(3),3))
    allocate(dUdy3(dUdxylo(1):dUdxyhi(1),dUdxylo(2):dUdxyhi(2),dUdxylo(3):dUdxyhi(3),3))
    allocate(dUdy4(dUdxylo(1):dUdxyhi(1),dUdxylo(2):dUdxyhi(2),dUdxylo(3):dUdxyhi(3),3))

    allocate(muxy1(tlo(1):thi(1),tlo(2):thi(2),tlo(3):thi(3)))
    allocate(muxy2(tlo(1):thi(1),tlo(2):thi(2),tlo(3):thi(3)))
    allocate(muxy3(tlo(1):thi(1),tlo(2):thi(2),tlo(3):thi(3)))
    allocate(muxy4(tlo(1):thi(1),tlo(2):thi(2),tlo(3):thi(3)))

    allocate(xixy1(tlo(1):thi(1),tlo(2):thi(2),tlo(3):thi(3)))
    allocate(xixy2(tlo(1):thi(1),tlo(2):thi(2),tlo(3):thi(3)))
    allocate(xixy3(tlo(1):thi(1),tlo(2):thi(2),tlo(3):thi(3)))
    allocate(xixy4(tlo(1):thi(1),tlo(2):thi(2),tlo(3):thi(3)))

    allocate(lamxy1(tlo(1):thi(1),tlo(2):thi(2),tlo(3):thi(3)))
    allocate(lamxy2(tlo(1):thi(1),tlo(2):thi(2),tlo(3):thi(3)))
    allocate(lamxy3(tlo(1):thi(1),tlo(2):thi(2),tlo(3):thi(3)))
    allocate(lamxy4(tlo(1):thi(1),tlo(2):thi(2),tlo(3):thi(3)))

    allocate(Ddiaxy1(tlo(1):thi(1),tlo(2):thi(2),tlo(3):thi(3),NSPEC))
    allocate(Ddiaxy2(tlo(1):thi(1),tlo(2):thi(2),tlo(3):thi(3),NSPEC))
    allocate(Ddiaxy3(tlo(1):thi(1),tlo(2):thi(2),tlo(3):thi(3),NSPEC))
    allocate(Ddiaxy4(tlo(1):thi(1),tlo(2):thi(2),tlo(3):thi(3),NSPEC))

    ! d./dx
    n = 1  ! rho
    do k=lo(3)-2,hi(3)+2
       call cellavg2dergausspt_2d(lo(1:2),hi(1:2), U(:,:,k,URHO), Ulo(1:2), Uhi(1:2), &
            dUdx1(:,:,k,n), dUdx2(:,:,k,n), dUdx3(:,:,k,n), dUdx4(:,:,k,n), &
            dUdxylo(1:2), dUdxyhi(1:2), 1)
    end do
    n = 2  ! my
    do k=lo(3)-2,hi(3)+2
       call cellavg2dergausspt_2d(lo(1:2),hi(1:2), U(:,:,k,UMY), Ulo(1:2), Uhi(1:2), &
            dUdx1(:,:,k,n), dUdx2(:,:,k,n), dUdx3(:,:,k,n), dUdx4(:,:,k,n), &
            dUdxylo(1:2), dUdxyhi(1:2), 1)
    end do
    n = 3  ! mz
    do k=lo(3)-2,hi(3)+2
       call cellavg2dergausspt_2d(lo(1:2),hi(1:2), U(:,:,k,UMZ), Ulo(1:2), Uhi(1:2), &
            dUdx1(:,:,k,n), dUdx2(:,:,k,n), dUdx3(:,:,k,n), dUdx4(:,:,k,n), &
            dUdxylo(1:2), dUdxyhi(1:2), 1)
    end do
 
    ! d./dy
    n = 1  ! rho
    do k=lo(3)-2,hi(3)+2
       call cellavg2dergausspt_2d(lo(1:2),hi(1:2), U(:,:,k,URHO), Ulo(1:2), Uhi(1:2), &
            dUdy1(:,:,k,n), dUdy2(:,:,k,n), dUdy3(:,:,k,n), dUdy4(:,:,k,n), &
            dUdxylo(1:2), dUdxyhi(1:2), 2)
    end do
    n = 2  ! mx
    do k=lo(3)-2,hi(3)+2
       call cellavg2dergausspt_2d(lo(1:2),hi(1:2), U(:,:,k,UMX), Ulo(1:2), Uhi(1:2), &
            dUdy1(:,:,k,n), dUdy2(:,:,k,n), dUdy3(:,:,k,n), dUdy4(:,:,k,n), &
            dUdxylo(1:2), dUdxyhi(1:2), 2)
    end do
    n = 3  ! mz
    do k=lo(3)-2,hi(3)+2
       call cellavg2dergausspt_2d(lo(1:2),hi(1:2), U(:,:,k,UMZ), Ulo(1:2), Uhi(1:2), &
            dUdy1(:,:,k,n), dUdy2(:,:,k,n), dUdy3(:,:,k,n), dUdy4(:,:,k,n), &
            dUdxylo(1:2), dUdxyhi(1:2), 2)
    end do
   
    do n=1,NVAR
       do k=lo(3)-3,hi(3)+3
          call cellavg2gausspt_2d(lo(1:2),hi(1:2), U(:,:,k,n), Ulo(1:2), Uhi(1:2), &
               Uxy1(:,:,k,n), Uxy2(:,:,k,n), Uxy3(:,:,k,n), Uxy4(:,:,k,n), &
               Uxylo(1:2), Uxyhi(1:2))
       end do
    end do

    call cc2zface_3d(lo,hi, mucc,g2lo,g2hi, muxy1, muxy2, muxy3, muxy4,tlo,thi)
    call cc2zface_3d(lo,hi, xicc,g2lo,g2hi, xixy1, xixy2, xixy3, xixy4,tlo,thi)
    call cc2zface_3d(lo,hi,lamcc,g2lo,g2hi,lamxy1,lamxy2,lamxy3,lamxy4,tlo,thi)
    do n=1,NSPEC
       call cc2zface_3d(lo,hi,Ddiacc(:,:,:,n),g2lo,g2hi, &
            Ddiaxy1(:,:,:,n),Ddiaxy2(:,:,:,n),Ddiaxy3(:,:,:,n),Ddiaxy4(:,:,:,n),tlo,thi)
    end do

    do j=lo(2),hi(2)
       do i=lo(1),hi(1)
          call diff_z(lo(3),hi(3), Uxy1(i,j,:,:), Uxylo(3), Uxyhi(3), &
               dUdx1(i,j,:,:), dUdy1(i,j,:,:), dUdxylo(3), dUdxyhi(3), &
               muxy1(i,j,:), xixy1(i,j,:), lamxy1(i,j,:), Ddiaxy1(i,j,:,:), tlo(3),thi(3), &
               fz(i,j,:,:), fzlo(3), fzhi(3), dxinv)
       end do
    end do
    !
    do j=lo(2),hi(2)
       do i=lo(1),hi(1)
          call diff_z(lo(3),hi(3), Uxy2(i,j,:,:), Uxylo(3), Uxyhi(3), &
               dUdx2(i,j,:,:), dUdy2(i,j,:,:), dUdxylo(3), dUdxyhi(3), &
               muxy2(i,j,:), xixy2(i,j,:), lamxy2(i,j,:), Ddiaxy2(i,j,:,:), tlo(3),thi(3), &
               fz(i,j,:,:), fzlo(3), fzhi(3), dxinv)
       end do
    end do
    !
    do j=lo(2),hi(2)
       do i=lo(1),hi(1)
          call diff_z(lo(3),hi(3), Uxy3(i,j,:,:), Uxylo(3), Uxyhi(3), &
               dUdx3(i,j,:,:), dUdy3(i,j,:,:), dUdxylo(3), dUdxyhi(3), &
               muxy3(i,j,:), xixy3(i,j,:), lamxy3(i,j,:), Ddiaxy3(i,j,:,:), tlo(3),thi(3), &
               fz(i,j,:,:), fzlo(3), fzhi(3), dxinv)
       end do
    end do
    !
    do j=lo(2),hi(2)
       do i=lo(1),hi(1)
          call diff_z(lo(3),hi(3), Uxy4(i,j,:,:), Uxylo(3), Uxyhi(3), &
               dUdx4(i,j,:,:), dUdy4(i,j,:,:), dUdxylo(3), dUdxyhi(3), &
               muxy4(i,j,:), xixy4(i,j,:), lamxy4(i,j,:), Ddiaxy4(i,j,:,:), tlo(3),thi(3), &
               fz(i,j,:,:), fzlo(3), fzhi(3), dxinv)
       end do
    end do

    deallocate(Uxy1,Uxy2,Uxy3,Uxy4,dUdx1,dUdx2,dUdx3,dUdx4,dUdy1,dUdy2,dUdy3,dUdy4)
    deallocate(muxy1,xixy1,lamxy1,Ddiaxy1)
    deallocate(muxy2,xixy2,lamxy2,Ddiaxy2)
    deallocate(muxy3,xixy3,lamxy3,Ddiaxy3)
    deallocate(muxy4,xixy4,lamxy4,Ddiaxy4)
    deallocate(mucc,xicc,lamcc,Ddiacc)

    return
  end subroutine difterm


  subroutine diff_xy(lo, hi, U, Ulo, Uhi, dUdz, dUdzlo, dUdzhi, &
       mu, xi, lam, Ddia, dfclo, dfchi, &
       fx, fxlo, fxhi, fy, fylo, fyhi, dxinv)
    
    use meth_params_module, only : NVAR, URHO, UMX, UMY, UMZ, NSPEC, QCVAR, QFVAR
    use weno_module, only : cellavg2gausspt_1d, cellavg2face_1d, cellavg2dergausspt_1d
    use convert_3d_module, only : cellavg2cc_3d
    use polyinterp_module, only : cc2xface_2d, cc2yface_2d
    use variables_module, only : ctoprim

    integer, intent(in) :: lo(3), hi(3), Ulo(3), Uhi(3), dUdzlo(3), dUdzhi(3), &
         dfclo(3), dfchi(3), fxlo(3), fxhi(3), fylo(3), fyhi(3)
    double precision,intent(in   )::dxinv(3)
    double precision,intent(in   )::   U(   Ulo(1):   Uhi(1),   Ulo(2):   Uhi(2),   Ulo(3):   Uhi(3),NVAR)
    double precision,intent(in   )::dUdz(dUdzlo(1):dUdzhi(1),dUdzlo(2):dUdzhi(2),dUdzlo(3):dUdzhi(3),4)
    double precision,intent(in   )::  mu( dfclo(1): dfchi(1), dfclo(2): dfchi(2), dfclo(3): dfchi(3)) 
    double precision,intent(in   )::  xi( dfclo(1): dfchi(1), dfclo(2): dfchi(2), dfclo(3): dfchi(3)) 
    double precision,intent(in   ):: lam( dfclo(1): dfchi(1), dfclo(2): dfchi(2), dfclo(3): dfchi(3)) 
    double precision,intent(in   )::Ddia( dfclo(1): dfchi(1), dfclo(2): dfchi(2), dfclo(3): dfchi(3),NSPEC) 
    double precision,intent(inout)::  fx(  fxlo(1):  fxhi(1),  fxlo(2):  fxhi(2),  fxlo(3):  fxhi(3),NVAR)
    double precision,intent(inout)::  fy(  fylo(1):  fyhi(1),  fylo(2):  fyhi(2),  fylo(3):  fyhi(3),NVAR)

    double precision, dimension(:,:,:), pointer ::  Uag, dUag, dUdzag, Ddia0
    double precision, dimension(:,:), pointer :: mu0, xi0, lam0
    double precision, allocatable, target ::    U1(:,:,:),    U2(:,:,:)
    double precision, allocatable, target ::   dU1(:,:,:),   dU2(:,:,:)
    double precision, allocatable, target :: dUdz1(:,:,:), dUdz2(:,:,:)
    double precision, allocatable :: Qc(:,:,:), Qf(:,:,:), dmom(:,:,:), dmdz(:,:,:)
    double precision, dimension(:,:), allocatable, target :: mu1, mu2, xi1, xi2, lam1, lam2
    double precision, dimension(:,:,:), allocatable, target :: Ddia1, Ddia2
    integer :: i, j, k, n, g
    integer :: g2lo(3), g2hi(3), g3lo(3), g3hi(3)
    integer :: tlo(3), thi(3), Qclo(3), Qchi(3), Qflo(3), Qfhi(3)
    double precision, parameter :: fac = 0.25d0 ! due to Gauss quadrature

    g2lo=1; g2hi=1; g3lo=1; g3hi=1; tlo=1; thi=1; Qclo=1; Qchi=1; Qflo=1; Qfhi=1;

    g3lo(1:2) = lo(1:2)-3  
    g3hi(1:2) = hi(1:2)+3
    allocate(U1(g3lo(1):g3hi(1),g3lo(2):g3hi(2),NVAR))
    allocate(U2(g3lo(1):g3hi(1),g3lo(2):g3hi(2),NVAR))

    g2lo(1:2) = lo(1:2)-2  
    g2hi(1:2) = hi(1:2)+2
    allocate(dU1(g2lo(1):g2hi(1),g2lo(2):g2hi(2),3))
    allocate(dU2(g2lo(1):g2hi(1),g2lo(2):g2hi(2),3))    

    allocate(dUdz1(g2lo(1):g2hi(1),g2lo(2):g2hi(2),3))
    allocate(dUdz2(g2lo(1):g2hi(1),g2lo(2):g2hi(2),3))    

    Qclo(1:2) = lo(1:2)-2
    Qchi(1:2) = hi(1:2)+2

    Qflo(1:2) = lo(1:2)
    Qfhi(1:2) = hi(1:2)+1

    Qflo(1:2) = lo(1:2)
    Qfhi(1:2) = hi(1:2)+1

    allocate(Qc  (Qclo(1):Qchi(1),Qclo(2):Qchi(2),QCVAR))
    allocate(Qf  (Qflo(1):Qfhi(1),Qflo(2):Qfhi(2),QFVAR))
    allocate(dmom(Qflo(1):Qfhi(1),Qflo(2):Qfhi(2),3))
    allocate(dmdz(Qflo(1):Qfhi(1),Qflo(2):Qfhi(2),3))

    allocate(  mu1(Qflo(1):Qfhi(1),Qflo(2):Qfhi(2)))
    allocate(  xi1(Qflo(1):Qfhi(1),Qflo(2):Qfhi(2)))
    allocate( lam1(Qflo(1):Qfhi(1),Qflo(2):Qfhi(2)))
    allocate(Ddia1(Qflo(1):Qfhi(1),Qflo(2):Qfhi(2),NSPEC))

    allocate(  mu2(Qflo(1):Qfhi(1),Qflo(2):Qfhi(2)))
    allocate(  xi2(Qflo(1):Qfhi(1),Qflo(2):Qfhi(2)))
    allocate( lam2(Qflo(1):Qfhi(1),Qflo(2):Qfhi(2)))
    allocate(Ddia2(Qflo(1):Qfhi(1),Qflo(2):Qfhi(2),NSPEC))

    do k = lo(3), hi(3)

       ! ----- compute x-direction flux first -----

       ! cell center => Gauss points on x-face
       call cc2xface_2d(lo,hi, mu(:,:,k),dfclo(1:2),dfchi(1:2), mu1, mu2,Qflo(1:2),Qfhi(1:2))
       call cc2xface_2d(lo,hi, xi(:,:,k),dfclo(1:2),dfchi(1:2), xi1, xi2,Qflo(1:2),Qfhi(1:2))
       call cc2xface_2d(lo,hi,lam(:,:,k),dfclo(1:2),dfchi(1:2),lam1,lam2,Qflo(1:2),Qfhi(1:2))
       do n=1,NSPEC
          call cc2xface_2d(lo,hi,Ddia(:,:,k,n),dfclo(1:2),dfchi(1:2), &
               Ddia1(:,:,n),Ddia2(:,:,n),Qflo(1:2),Qfhi(1:2))
       end do

       ! cell-average => cell-avg-in-x and Gauss-point-in-y of d?/dy
       do n=1,3
          do i=lo(1)-2,hi(1)+2
             call cellavg2dergausspt_1d(lo(2),hi(2), U(i,:,k,n), Ulo(2), Uhi(2), &
                  dU1(i,:,n), dU2(i,:,n), g2lo(2), g2hi(2))
          end do
       end do

       ! cell-average => cell-avg-in-x and Gauss-point-in-y
       do n=1,NVAR
          do i=lo(1)-3,hi(1)+3
             call cellavg2gausspt_1d(lo(2),hi(2), U(i,:,k,n), Ulo(2), Uhi(2), &
                  U1(i,:,n), U2(i,:,n), g3lo(2), g3hi(2))
          end do
       end do

       n = 1  ! rho
       do i=lo(1)-2,hi(1)+2
          call cellavg2gausspt_1d(lo(2),hi(2), dUdz(i,:,k,URHO), dUdzlo(2), dUdzhi(2), &
               dUdz1(i,:,n), dUdz2(i,:,n), g2lo(2), g2hi(2))
       end do
       !
       n = 2  ! mx
       do i=lo(1)-2,hi(1)+2
          call cellavg2gausspt_1d(lo(2),hi(2), dUdz(i,:,k,UMX), dUdzlo(2), dUdzhi(2), &
               dUdz1(i,:,n), dUdz2(i,:,n), g2lo(2), g2hi(2))
       end do
       !
       n = 3  ! mz 
       do i=lo(1)-2,hi(1)+2
          call cellavg2gausspt_1d(lo(2),hi(2), dUdz(i,:,k,UMZ), dUdzlo(2), dUdzhi(2), &
               dUdz1(i,:,n), dUdz2(i,:,n), g2lo(2), g2hi(2))
       end do
       
       do g=1,2
          
          if (g .eq. 1) then
             Uag    =>    U1
             dUag   =>   dU1
             dUdzag => dUdz1
             mu0    =>   mu1
             xi0    =>   xi1
             lam0   =>  lam1
             Ddia0  => Ddia1
          else
             Uag    =>    U2
             dUag   =>   dU2
             dUdzag => dUdz2
             mu0    =>   mu2
             xi0    =>   xi2
             lam0   =>  lam2
             Ddia0  => Ddia2
          end if

          do n=1,3
             do j=lo(2),hi(2)
                ! cell-avg-in-x and Gauss-point-in-y => xface and Gauss-point-in-y
                call cellavg2face_1d(lo(1),hi(1)+1, dUag(:,j,n),g2lo(1),g2hi(1), &
                     dmom(:,j,n),Qflo(1),Qfhi(1))
             end do
          end do

          do n=1,3
             do j=lo(2),hi(2)
                ! cell-avg-in-x and Gauss-point-in-y => xface and Gauss-point-in-y
                call cellavg2face_1d(lo(1),hi(1)+1, dUdzag(:,j,n),g2lo(1),g2hi(1), &
                     dmdz(:,j,n),Qflo(1),Qfhi(1))
             end do
          end do

          do n=1,NVAR
             do j=lo(2),hi(2)
                ! cell-avg-in-x and Gauss-point-in-y => xface and Gauss-point-in-y
                call cellavg2face_1d(lo(1),hi(1)+1, Uag(:,j,n),g3lo(1),g3hi(1), &
                     Qf(:,j,n),Qflo(1),Qfhi(1))
             end do
          end do

          tlo(1) = lo(1)-2
          tlo(2) = lo(2)
          thi(1) = hi(1)+2
          thi(2) = hi(2)
          do n=1,NVAR
             ! cell-avg-in-x and Gauss-point-in-y => cell-center-in-x and Gauss-point-in-y
             call cellavg2cc_3d(tlo,thi, Uag(:,:,n),g3lo,g3hi, &
                  Qc(:,:,n),Qclo,Qchi,idir=1)
          end do

          tlo(1:2) = lo(1:2)
          thi(1) = hi(1)+1
          thi(2) = hi(2) 
          call ctoprim(tlo,thi, Qf, Qflo,Qfhi,QFVAR)

          tlo(1) = lo(1)-2
          thi(1) = hi(1)+2
          call ctoprim(tlo,thi, Qc, Qclo,Qchi,QCVAR)

          tlo(1:2) = lo(1:2)
          thi(1) = hi(1)+1
          thi(2) = hi(2)
          call comp_diff_flux_x(tlo(1:2), thi(1:2), k, fx, fxlo, fxhi, &
               Qf, mu0, xi0, lam0, Ddia0, dmom, dmdz, Qflo, Qfhi, &
               Qc, Qclo, Qchi, dxinv, fac)

          Nullify(Uag,dUag,dUdzag,mu0,xi0,lam0,Ddia0)
       end do

       ! ----- compute y-direction flux -----

       ! cell center => Gauss points on y-face
       call cc2yface_2d(lo,hi, mu(:,:,k),dfclo(1:2),dfchi(1:2), mu1, mu2,Qflo(1:2),Qfhi(1:2))
       call cc2yface_2d(lo,hi, xi(:,:,k),dfclo(1:2),dfchi(1:2), xi1, xi2,Qflo(1:2),Qfhi(1:2))
       call cc2yface_2d(lo,hi,lam(:,:,k),dfclo(1:2),dfchi(1:2),lam1,lam2,Qflo(1:2),Qfhi(1:2))
       do n=1,NSPEC
          call cc2yface_2d(lo,hi,Ddia(:,:,k,n),dfclo(1:2),dfchi(1:2), &
               Ddia1(:,:,n),Ddia2(:,:,n),Qflo(1:2),Qfhi(1:2))
       end do
       
       ! cell-average => cell-avg-in-y and Gauss-point-in-x of d?/dx
       do n=1,3
          do j=lo(2)-2,hi(2)+2
             call cellavg2dergausspt_1d(lo(1),hi(1), U(:,j,k,n), Ulo(1), Uhi(1), &
                  dU1(:,j,n), dU2(:,j,n), g2lo(1), g2hi(1))
          end do
       end do

       ! cell-average => cell-avg-in-y and Gauss-point-in-x
       do n=1,NVAR
          do j=lo(2)-3,hi(2)+3
             call cellavg2gausspt_1d(lo(1),hi(1), U(:,j,k,n), Ulo(1), Uhi(1), &
                  U1(:,j,n), U2(:,j,n), g3lo(1), g3hi(1))
          end do
       end do

       n = 1 ! rho
       do j=lo(2)-2,hi(2)+2
          call cellavg2gausspt_1d(lo(1),hi(1), dUdz(:,j,k,URHO), dUdzlo(1), dUdzhi(1), &
               dUdz1(:,j,n), dUdz2(:,j,n), g2lo(1), g2hi(1))
       end do
       !
       n = 2 ! my 
       do j=lo(2)-2,hi(2)+2
          call cellavg2gausspt_1d(lo(1),hi(1), dUdz(:,j,k,UMY), dUdzlo(1), dUdzhi(1), &
               dUdz1(:,j,n), dUdz2(:,j,n), g2lo(1), g2hi(1))
       end do
       n = 3 ! mz
       do j=lo(2)-2,hi(2)+2
          call cellavg2gausspt_1d(lo(1),hi(1), dUdz(:,j,k,UMZ), dUdzlo(1), dUdzhi(1), &
               dUdz1(:,j,n), dUdz2(:,j,n), g2lo(1), g2hi(1))
       end do

       do g=1,2
          
          if (g .eq. 1) then
             Uag    =>    U1
             dUag   =>   dU1
             dUdzag => dUdz1
             mu0    =>   mu1
             xi0    =>   xi1
             lam0   =>  lam1
             Ddia0  => Ddia1
          else
             Uag    =>    U2
             dUag   =>   dU2
             dUdzag => dUdz2
             mu0    =>   mu2
             xi0    =>   xi2
             lam0   =>  lam2
             Ddia0  => Ddia2
          end if

          do n=1,3
             do i=lo(1),hi(1)
                ! cell-avg-in-y and Gauss-point-in-x => yface and Gauss-point-in-x
                call cellavg2face_1d(lo(2),hi(2)+1, dUag(i,:,n),g2lo(2),g2hi(2), &
                     dmom(i,:,n),Qflo(2),Qfhi(2))
             end do
          end do

          do n=1,3
             do i=lo(1),hi(1)
                ! cell-avg-in-y and Gauss-point-in-x => yface and Gauss-point-in-x
                call cellavg2face_1d(lo(2),hi(2)+1, dUdzag(i,:,n),g2lo(2),g2hi(2), &
                     dmdz(i,:,n),Qflo(2),Qfhi(2))
             end do
          end do

          do n=1,NVAR
             do i=lo(1),hi(1)
                ! cell-avg-in-y and Gauss-point-in-x => yface and Gauss-point-in-x
                call cellavg2face_1d(lo(2),hi(2)+1, Uag(i,:,n),g3lo(2),g3hi(2), &
                     Qf(i,:,n),Qflo(2),Qfhi(2))
             end do
          end do

          tlo(1) = lo(1)
          tlo(2) = lo(2)-2
          thi(1) = hi(1)
          thi(2) = hi(2)+2
          do n=1,NVAR
             ! cell-avg-in-y and Gauss-point-in-x => cell-center-in-y and Gauss-point-in-x
             call cellavg2cc_3d(tlo,thi, Uag(:,:,n),g3lo,g3hi, &
                  Qc(:,:,n),Qclo,Qchi,idir=2)
          end do

          tlo(1:2) = lo(1:2)
          thi(1) = hi(1)
          thi(2) = hi(2)+1
          call ctoprim(tlo,thi, Qf, Qflo,Qfhi,QFVAR)

          tlo(2) = lo(2)-2
          thi(2) = hi(2)+2
          call ctoprim(tlo,thi, Qc, Qclo,Qchi,QCVAR)

          tlo(1:2) = lo(1:2)
          thi(1) = hi(1)
          thi(2) = hi(2)+1
          call comp_diff_flux_y(tlo(1:2), thi(1:2), k, fy, fylo, fyhi, &
               Qf, mu0, xi0, lam0, Ddia0, dmom, dmdz, Qflo, Qfhi, &
               Qc, Qclo, Qchi, dxinv, fac)
          
          Nullify(Uag,dUag,dUdzag,mu0,xi0,lam0,Ddia0)
       end do

    end do

    deallocate(U1,U2,dU1,dU2,dUdz1,dUdz2,Qc,Qf,dmom,dmdz)
    deallocate(mu1,xi1,lam1,Ddia1)
    deallocate(mu2,xi2,lam2,Ddia2)

  end subroutine diff_xy


  subroutine comp_diff_flux_x(lo, hi, k, flx, flo, fhi, &
       Qf, mu, xi, lam, Ddia, dmdy, dmdz, Qflo, Qfhi, &
       Qc, Qclo, Qchi, dxinv, fac)

    use meth_params_module
    use derivative_stencil_module, only : FD4

    integer, intent(in) :: lo(2), hi(2), k, flo(3), fhi(3), Qflo(2), Qfhi(2), Qclo(2), Qchi(2)
    double precision, intent(in) :: dxinv(3), fac
    double precision, intent(inout) ::  flx( flo(1): fhi(1), flo(2): fhi(2),flo(3):fhi(3),NVAR)
    double precision, intent(in   ) ::   Qf(Qflo(1):Qfhi(1),Qflo(2):Qfhi(2),QFVAR)
    double precision, intent(in   ) ::   mu(Qflo(1):Qfhi(1),Qflo(2):Qfhi(2))
    double precision, intent(in   ) ::   xi(Qflo(1):Qfhi(1),Qflo(2):Qfhi(2))
    double precision, intent(in   ) ::  lam(Qflo(1):Qfhi(1),Qflo(2):Qfhi(2))
    double precision, intent(in   ) :: Ddia(Qflo(1):Qfhi(1),Qflo(2):Qfhi(2),NSPEC)
    double precision, intent(in   ) :: dmdy(Qflo(1):Qfhi(1),Qflo(2):Qfhi(2),3)
    double precision, intent(in   ) :: dmdz(Qflo(1):Qfhi(1),Qflo(2):Qfhi(2),3)
    double precision, intent(in   ) ::   Qc(Qclo(1):Qchi(1),Qclo(2):Qchi(2),QCVAR)

    integer :: i, j, n, UYN, QYN, QXN, QHN
    double precision :: tauxx, tauxy, tauxz
    double precision :: dudx, dudy, dudz, dvdx, dvdy, dwdx, dwdz, divu, rhoinv
    double precision :: dTdx, dXdx, Vd
    double precision, dimension(lo(1):hi(1)) :: dlnpdx, Vc

    do j=lo(2),hi(2)
       do i=lo(1),hi(1)
          ! viscous stress
          dudx = dxinv(1)*(FD4(-2)*Qc(i-2,j,QU) + FD4(-1)*Qc(i-1,j,QU) &
               + FD4(0)*Qc(i,j,QU) + FD4(1)*Qc(i+1,j,QU))
          dvdx = dxinv(1)*(FD4(-2)*Qc(i-2,j,QV) + FD4(-1)*Qc(i-1,j,QV) &
               + FD4(0)*Qc(i,j,QV) + FD4(1)*Qc(i+1,j,QV))
          dwdx = dxinv(1)*(FD4(-2)*Qc(i-2,j,QW) + FD4(-1)*Qc(i-1,j,QW) &
               + FD4(0)*Qc(i,j,QW) + FD4(1)*Qc(i+1,j,QW))
          rhoinv = 1.d0/Qf(i,j,QRHO)
          dudy = dxinv(2)*rhoinv*(dmdy(i,j,2)-Qf(i,j,QU)*dmdy(i,j,1))
          dvdy = dxinv(2)*rhoinv*(dmdy(i,j,3)-Qf(i,j,QV)*dmdy(i,j,1))
          dudz = dxinv(3)*rhoinv*(dmdz(i,j,2)-Qf(i,j,QU)*dmdz(i,j,1))
          dwdz = dxinv(3)*rhoinv*(dmdz(i,j,3)-Qf(i,j,QW)*dmdz(i,j,1))
          divu = dudx + dvdy + dwdz
          tauxx = mu(i,j)*(2.d0*dudx-twoThirds*divu) + xi(i,j)*divu
          tauxy = mu(i,j)*(dudy+dvdx)
          tauxz = mu(i,j)*(dudz+dwdx)
          flx(i,j,k,UMX)   = flx(i,j,k,UMX)   - fac*tauxx
          flx(i,j,k,UMY)   = flx(i,j,k,UMY)   - fac*tauxy
          flx(i,j,k,UMZ)   = flx(i,j,k,UMZ)   - fac*tauxz
          flx(i,j,k,UEDEN) = flx(i,j,k,UEDEN) - fac* &
               (tauxx*Qf(i,j,QU) + tauxy*Qf(i,j,QV) + tauxz*Qf(i,j,QW))

          ! thermal conduction
          dTdx = dxinv(1) * (FD4(-2)*Qc(i-2,j,QTEMP) + FD4(-1)*Qc(i-1,j,QTEMP) &
               + FD4(0)*Qc(i,j,QTEMP) + FD4(1)*Qc(i+1,j,QTEMP))
          flx(i,j,k,UEDEN) = flx(i,j,k,UEDEN) - fac*lam(i,j)*dTdx

          ! compute dpdx
          dlnpdx(i) = dxinv(1) * (FD4(-2)*Qc(i-2,j,QPRES) + FD4(-1)*Qc(i-1,j,QPRES) &
               + FD4(0)*Qc(i,j,QPRES) + FD4(1)*Qc(i+1,j,QPRES)) / Qf(i,j,QPRES)
          Vc(i) = 0.d0
       end do

       do n=1,NSPEC
          UYN = UFS+n-1
          QYN = QFY+n-1
          QXN = QFX+n-1
          QHN = QFH+n-1
          do i = lo(1), hi(1)
             dXdx = dxinv(1) * (FD4(-2)*Qc(i-2,j,QXN) + FD4(-1)*Qc(i-1,j,QXN) &
                  + FD4(0)*Qc(i,j,QXN) + FD4(1)*Qc(i+1,j,QXN))
             Vd = -Ddia(i,j,n)*(dXdx + (Qf(i,j,QXN)-Qf(i,j,QYN))*dlnpdx(i))
             
             flx(i,j,k,UYN) = flx(i,j,k,UYN) + fac*Vd
             Vc(i) = Vc(i) + Vd
             flx(i,j,k,UEDEN) = flx(i,j,k,UEDEN) + fac*Vd*Qf(i,j,QHN)
          end do
       end do

       do n=1,NSPEC
          UYN = UFS+n-1
          QYN = QFY+n-1
          QHN = QFH+n-1
          do i = lo(1), hi(1)
             flx(i,j,k,UYN )  = flx(i,j,k,UYN  ) - (fac*Qf(i,j,QYN)*Vc(i))
             flx(i,j,k,UEDEN) = flx(i,j,k,UEDEN) - (fac*Qf(i,j,QYN)*Vc(i))*Qf(i,j,QHN)
          end do
       end do
    end do

  end subroutine comp_diff_flux_x


  subroutine comp_diff_flux_y(lo, hi, k, flx, flo, fhi, &
       Qf, mu, xi, lam, Ddia, dmdx, dmdz, Qflo, Qfhi, &
       Qc, Qclo, Qchi, dxinv, fac)

    use meth_params_module
    use derivative_stencil_module, only : FD4

    integer, intent(in) :: lo(2), hi(2), k, flo(3), fhi(3), Qflo(2), Qfhi(2), Qclo(2), Qchi(2)
    double precision, intent(in) :: dxinv(3), fac
    double precision, intent(inout) ::  flx( flo(1): fhi(1), flo(2): fhi(2),flo(3):fhi(3),NVAR)
    double precision, intent(in   ) ::   Qf(Qflo(1):Qfhi(1),Qflo(2):Qfhi(2),QFVAR)
    double precision, intent(in   ) ::   mu(Qflo(1):Qfhi(1),Qflo(2):Qfhi(2))
    double precision, intent(in   ) ::   xi(Qflo(1):Qfhi(1),Qflo(2):Qfhi(2))
    double precision, intent(in   ) ::  lam(Qflo(1):Qfhi(1),Qflo(2):Qfhi(2))
    double precision, intent(in   ) :: Ddia(Qflo(1):Qfhi(1),Qflo(2):Qfhi(2),NSPEC)
    double precision, intent(in   ) :: dmdx(Qflo(1):Qfhi(1),Qflo(2):Qfhi(2),3)
    double precision, intent(in   ) :: dmdz(Qflo(1):Qfhi(1),Qflo(2):Qfhi(2),3)
    double precision, intent(in   ) ::   Qc(Qclo(1):Qchi(1),Qclo(2):Qchi(2),QCVAR)

    integer :: i, j, n, UYN, QYN, QXN, QHN
    double precision :: tauyy, tauxy, tauyz
    double precision :: dudx, dudy, dvdx, dvdy, dvdz, dwdy, dwdz, divu, rhoinv
    double precision :: dTdy, dXdy, Vd
    double precision, allocatable :: dlnpdy(:,:), Vc(:,:)

    allocate(dlnpdy(lo(1):hi(1),lo(2):hi(2)))
    allocate(    Vc(lo(1):hi(1),lo(2):hi(2)))

    do j=lo(2),hi(2)
       do i=lo(1),hi(1)
          ! viscous stress
          dudy = dxinv(2)*(FD4(-2)*Qc(i,j-2,QU) + FD4(-1)*Qc(i,j-1,QU) &
               + FD4(0)*Qc(i,j,QU) + FD4(1)*Qc(i,j+1,QU))
          dvdy = dxinv(2)*(FD4(-2)*Qc(i,j-2,QV) + FD4(-1)*Qc(i,j-1,QV) &
               + FD4(0)*Qc(i,j,QV) + FD4(1)*Qc(i,j+1,QV))
          dwdy = dxinv(2)*(FD4(-2)*Qc(i,j-2,QW) + FD4(-1)*Qc(i,j-1,QW) &
               + FD4(0)*Qc(i,j,QW) + FD4(1)*Qc(i,j+1,QW))
          rhoinv = 1.d0/Qf(i,j,QRHO)
          dudx = dxinv(1)*rhoinv*(dmdx(i,j,2)-Qf(i,j,QU)*dmdx(i,j,1))
          dvdx = dxinv(1)*rhoinv*(dmdx(i,j,3)-Qf(i,j,QV)*dmdx(i,j,1))
          dvdz = dxinv(3)*rhoinv*(dmdz(i,j,2)-Qf(i,j,QV)*dmdz(i,j,1))
          dwdz = dxinv(3)*rhoinv*(dmdz(i,j,3)-Qf(i,j,QW)*dmdz(i,j,1))
          divu = dudx + dvdy + dwdz
          tauyy = mu(i,j)*(2.d0*dvdy-twoThirds*divu) + xi(i,j)*divu
          tauxy = mu(i,j)*(dudy+dvdx)
          tauyz = mu(i,j)*(dwdy+dvdz)
          flx(i,j,k,UMX)   = flx(i,j,k,UMX)   - fac*tauxy
          flx(i,j,k,UMY)   = flx(i,j,k,UMY)   - fac*tauyy
          flx(i,j,k,UMZ)   = flx(i,j,k,UMZ)   - fac*tauyz
          flx(i,j,k,UEDEN) = flx(i,j,k,UEDEN) - fac* &
               (tauxy*Qf(i,j,QU) + tauyy*Qf(i,j,QV) + tauyz*Qf(i,j,QW))

          ! thermal conduction
          dTdy = dxinv(2) * (FD4(-2)*Qc(i,j-2,QTEMP) + FD4(-1)*Qc(i,j-1,QTEMP) &
               + FD4(0)*Qc(i,j,QTEMP) + FD4(1)*Qc(i,j+1,QTEMP))
          flx(i,j,k,UEDEN) = flx(i,j,k,UEDEN) - fac*lam(i,j)*dTdy

          ! compute dpdy
          dlnpdy(i,j) = dxinv(2) * (FD4(-2)*Qc(i,j-2,QPRES) + FD4(-1)*Qc(i,j-1,QPRES) &
               + FD4(0)*Qc(i,j,QPRES) + FD4(1)*Qc(i,j+1,QPRES)) / Qf(i,j,QPRES)
          Vc(i,j) = 0.d0
       end do
    end do

    do n=1,NSPEC
       UYN = UFS+n-1
       QYN = QFY+n-1
       QXN = QFX+n-1
       QHN = QFH+n-1
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)
             dXdy = dxinv(2) * (FD4(-2)*Qc(i,j-2,QXN) + FD4(-1)*Qc(i,j-1,QXN) &
                  + FD4(0)*Qc(i,j,QXN) + FD4(1)*Qc(i,j+1,QXN))
             Vd = -Ddia(i,j,n)*(dXdy + (Qf(i,j,QXN)-Qf(i,j,QYN))*dlnpdy(i,j))
             
             flx(i,j,k,UYN) = flx(i,j,k,UYN) + fac*Vd
             Vc(i,j) = Vc(i,j) + Vd
             flx(i,j,k,UEDEN) = flx(i,j,k,UEDEN) + fac*Vd*Qf(i,j,QHN)
          end do
       end do
    end do

    do n=1,NSPEC
       UYN = UFS+n-1
       QYN = QFY+n-1
       QHN = QFH+n-1
       do j= lo(2), hi(2)
          do i = lo(1), hi(1)
             flx(i,j,k,UYN )  = flx(i,j,k,UYN  ) - (fac*Qf(i,j,QYN)*Vc(i,j))
             flx(i,j,k,UEDEN) = flx(i,j,k,UEDEN) - (fac*Qf(i,j,QYN)*Vc(i,j))*Qf(i,j,QHN)
          end do
       end do
    end do

    deallocate(dlnpdy,Vc)

  end subroutine comp_diff_flux_y


  subroutine diff_z(lo, hi, U, Ulo, Uhi, dUdx, dUdy, dlo, dhi, &
       mu, xi, lam, Ddia, clo, chi, f, flo, fhi, dxinv)

    use meth_params_module, only : NVAR, URHO, UMX, UMY, UMZ, NSPEC, QCVAR, QFVAR
    use convert_3d_module, only : cellavg2cc_3d
    use weno_module, only : cellavg2face_1d
    use variables_module, only : ctoprim

    integer, intent(in) :: lo, hi, Ulo, Uhi, dlo, dhi, clo, chi, flo, fhi
    double precision, intent(in   ) ::    U(Ulo:Uhi,NVAR)
    double precision, intent(in   ) :: dUdx(dlo:dhi,3)
    double precision, intent(in   ) :: dUdy(dlo:dhi,3)
    double precision, intent(in   ) ::   mu(clo:chi)
    double precision, intent(in   ) ::   xi(clo:chi)
    double precision, intent(in   ) ::  lam(clo:chi)
    double precision, intent(in   ) :: Ddia(clo:chi,NSPEC)
    double precision, intent(inout) ::    f(flo:fhi,NVAR)
    double precision, intent(in) :: dxinv(3)

    integer :: n
    integer :: Qclo(3), Qchi(3), Qflo(3), Qfhi(3), tlo(3), thi(3)
    double precision, allocatable :: Qc(:,:), Qf(:,:)  ! cell-center and face
    double precision, allocatable :: dmdx(:,:), dmdy(:,:)
    double precision, parameter :: fac = 0.25d0 ! due to Gauss quadrature

    Qclo = 1;  Qchi = 1;  Qflo = 1;  Qfhi = 1; tlo = 1; thi = 1;

    Qclo(3) = lo-2   ! need 2 ghost cells for fourth-order derivatives on face
    Qchi(3) = hi+2

    Qflo(3) = lo
    Qfhi(3) = hi+1

    allocate(  Qc(Qclo(3):Qchi(3), QCVAR))
    allocate(  Qf(Qflo(3):Qfhi(3), QFVAR))
    allocate(dmdx(Qflo(3):Qfhi(3), 3))
    allocate(dmdy(Qflo(3):Qfhi(3), 3))

    ! cell-centered variables
    tlo(3) = Ulo
    thi(3) = Uhi
    do n=1,NVAR
       call cellavg2cc_3d(Qclo, Qchi, U(:,n), tlo, thi, Qc(:,n), Qclo, Qchi, 3)
    end do
    !
    call ctoprim(Qclo,Qchi, Qc, Qclo,Qchi,QCVAR)

    ! face variables
    do n=1,NVAR
       call cellavg2face_1d(Qflo(3),Qfhi(3), U(:,n), Ulo,Uhi, Qf(:,n), Qflo(3), Qfhi(3))
    end do
    !
    call ctoprim(Qflo,Qfhi, Qf, Qflo,Qfhi,QFVAR)

    ! d./dx
    do n=1,3
       call cellavg2face_1d(Qflo(3),Qfhi(3), dUdx(:,n), dlo,dhi, dmdx(:,n), Qflo(3), Qfhi(3))
    end do
    ! d./dy
    do n=1,3
       call cellavg2face_1d(Qflo(3),Qfhi(3), dUdy(:,n), dlo,dhi, dmdy(:,n), Qflo(3), Qfhi(3))
    end do

    ! It is assumed that Qflo(3) == clo & Qfhi(3) == chi

    call comp_diff_flux_z(lo, hi+1, f, flo, fhi, &
         Qf, mu, xi, lam, Ddia, dmdx, dmdy, Qflo(3), Qfhi(3), &
         Qc, Qclo(3), Qchi(3), dxinv, fac)

    deallocate(Qc,Qf,dmdx,dmdy)

  end subroutine diff_z


  subroutine comp_diff_flux_z(lo, hi, flx, flo, fhi, &
       Qf, mu, xi, lam, Ddia, dmdx, dmdy, Qflo, Qfhi, &
       Qc, Qclo, Qchi, dxinv, fac)

    use meth_params_module
    use derivative_stencil_module, only : FD4

    integer, intent(in) :: lo, hi, flo, fhi, Qflo, Qfhi, Qclo, Qchi
    double precision, intent(in) :: dxinv(3), fac
    double precision, intent(inout) ::  flx( flo: fhi,NVAR)
    double precision, intent(in   ) ::   Qf(Qflo:Qfhi,QFVAR)
    double precision, intent(in   ) ::   mu(Qflo:Qfhi)
    double precision, intent(in   ) ::   xi(Qflo:Qfhi)
    double precision, intent(in   ) ::  lam(Qflo:Qfhi)
    double precision, intent(in   ) :: Ddia(Qflo:Qfhi,NSPEC)
    double precision, intent(in   ) :: dmdx(Qflo:Qfhi,3)
    double precision, intent(in   ) :: dmdy(Qflo:Qfhi,3)
    double precision, intent(in   ) ::   Qc(Qclo:Qchi,QCVAR)

    integer :: k, n, UYN, QYN, QXN, QHN
    double precision :: tauzz, tauxz, tauyz
    double precision :: dudx, dudz, dvdy, dvdz, dwdx, dwdy, dwdz, divu, rhoinv
    double precision :: dTdz, dXdz, Vd
    double precision, dimension(lo:hi) :: dlnpdz, Vc

    do k=lo,hi
       ! viscous stress
       dudz = dxinv(3)*(FD4(-2)*Qc(k-2,QU) + FD4(-1)*Qc(k-1,QU) &
            + FD4(0)*Qc(k,QU) + FD4(1)*Qc(k+1,QU))
       dvdz = dxinv(3)*(FD4(-2)*Qc(k-2,QV) + FD4(-1)*Qc(k-1,QV) &
            + FD4(0)*Qc(k,QV) + FD4(1)*Qc(k+1,QV))
       dwdz = dxinv(3)*(FD4(-2)*Qc(k-2,QW) + FD4(-1)*Qc(k-1,QW) &
            + FD4(0)*Qc(k,QW) + FD4(1)*Qc(k+1,QW))
       rhoinv = 1.d0/Qf(k,QRHO)
       dudx = dxinv(1)*rhoinv*(dmdx(k,2)-Qf(k,QU)*dmdx(k,1))
       dwdx = dxinv(1)*rhoinv*(dmdx(k,3)-Qf(k,QW)*dmdx(k,1))
       dvdy = dxinv(2)*rhoinv*(dmdy(k,2)-Qf(k,QV)*dmdy(k,1))
       dwdy = dxinv(2)*rhoinv*(dmdy(k,3)-Qf(k,QW)*dmdy(k,1))
       divu = dudx + dvdy + dwdz
       tauxz = mu(k)*(dudz+dwdx)
       tauyz = mu(k)*(dvdz+dwdy)
       tauzz = mu(k)*(2.d0*dwdz-twoThirds*divu) + xi(k)*divu
       flx(k,UMX)   = flx(k,UMX)   - fac*tauxz
       flx(k,UMY)   = flx(k,UMY)   - fac*tauyz
       flx(k,UMZ)   = flx(k,UMZ)   - fac*tauzz
       flx(k,UEDEN) = flx(k,UEDEN) - fac* &
            (tauxz*Qf(k,QU) + tauyz*Qf(k,QV) + tauzz*Qf(k,QW))

       !thermal conduction
       dTdz = dxinv(3) * (FD4(-2)*Qc(k-2,QTEMP) + FD4(-1)*Qc(k-1,QTEMP) &
            + FD4(0)*Qc(k,QTEMP) + FD4(1)*Qc(k+1,QTEMP))
       flx(k,UEDEN) = flx(k,UEDEN) - fac*lam(k)*dTdz

       ! compute dpdz
       dlnpdz(k) = dxinv(3) * (FD4(-2)*Qc(k-2,QPRES) + FD4(-1)*Qc(k-1,QPRES) &
            + FD4(0)*Qc(k,QPRES) + FD4(1)*Qc(k+1,QPRES)) / Qf(k,QPRES)
       Vc(k) = 0.d0
    end do

    do n=1,NSPEC
       UYN = UFS+n-1
       QYN = QFY+n-1
       QXN = QFX+n-1
       QHN = QFH+n-1
       do k = lo, hi
          dXdz = dxinv(3) * (FD4(-2)*Qc(k-2,QXN) + FD4(-1)*Qc(k-1,QXN) &
               + FD4(0)*Qc(k,QXN) + FD4(1)*Qc(k+1,QXN))
          Vd = -Ddia(k,n)*(dXdz + (Qf(k,QXN)-Qf(k,QYN))*dlnpdz(k))
             
          flx(k,UYN) = flx(k,UYN) + fac*Vd
          Vc(k) = Vc(k) + Vd
          flx(k,UEDEN) = flx(k,UEDEN) + fac*Vd*Qf(k,QHN)
       end do
    end do
    
    do n=1,NSPEC
       UYN = UFS+n-1
       QYN = QFY+n-1
       QHN = QFH+n-1
       do k = lo, hi
          flx(k,UYN )  = flx(k,UYN  ) - (fac*Qf(k,QYN)*Vc(k))
          flx(k,UEDEN) = flx(k,UEDEN) - (fac*Qf(k,QYN)*Vc(k))*Qf(k,QHN)
       end do
    end do

  end subroutine comp_diff_flux_z

end module difterm_module
