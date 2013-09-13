module difterm_module

  use meth_params_module, only : NVAR, NSPEC, QCVAR, QFVAR
  use weno_module, only : cellavg2gausspt_1d, cellavg2face_1d
  use convert_module, only : cellavg2cc
  use variables_module, only : ctoprim
  use transport_properties, only : get_transport_properties

  implicit none

  private

  public :: difterm

contains

  subroutine difterm(lo,hi,U,Ulo,Uhi,fx,fxlo,fxhi,fy,fylo,fyhi,dxinv)
    integer, intent(in) :: lo(2), hi(2), Ulo(2), Uhi(2), fxlo(2), fxhi(2), fylo(2), fyhi(2)
    double precision, intent(in   ) :: dxinv(2)
    double precision, intent(in   ) ::  U( Ulo(1): Uhi(1), Ulo(2): Uhi(2),NVAR)
    double precision, intent(inout) :: fx(fxlo(1):fxhi(1),fxlo(2):fxhi(2),NVAR)
    double precision, intent(inout) :: fy(fylo(1):fyhi(1),fylo(2):fyhi(2),NVAR)

    double precision, dimension(:,:,:), pointer :: Uag => Null()
    double precision, allocatable, target :: U1(:,:,:), U2(:,:,:)
    double precision, allocatable :: Qc(:,:,:), Qf(:,:,:)
    double precision, allocatable :: mu(:,:), xi(:,:), lam(:,:), Ddia(:,:,:)
    integer :: i, j, n, g
    integer :: tlo(3), thi(3), Qclo(3), Qchi(3), Qflo(3), Qfhi(3), g3lo(2), g3hi(2)

    g3lo = lo-3;  g3hi = hi+3
    allocate(U1(g3lo(1):g3hi(1),g3lo(2):g3hi(2),NVAR))
    allocate(U2(g3lo(1):g3hi(1),g3lo(2):g3hi(2),NVAR))

    tlo = 1; thi = 1; Qclo = 1;  Qchi = 1;  Qflo = 1;  Qfhi = 1

    Qclo(1:2) = lo(1:2)-2
    Qchi(1:2) = hi(1:2)+2

    Qflo(1:2) = lo(1:2)
    Qfhi(1:2) = hi(1:2)+1

    allocate(Qc  (Qclo(1):Qchi(1),Qclo(2):Qchi(2),QCVAR))
    allocate(Qf  (Qflo(1):Qfhi(1),Qflo(2):Qfhi(2),QFVAR))
    allocate(mu  (Qflo(1):Qfhi(1),Qflo(2):Qfhi(2)))
    allocate(xi  (Qflo(1):Qfhi(1),Qflo(2):Qfhi(2)))
    allocate(lam (Qflo(1):Qfhi(1),Qflo(2):Qfhi(2)))
    allocate(Ddia(Qflo(1):Qfhi(1),Qflo(2):Qfhi(2),NSPEC))

    ! ----- compute x-direction flux first -----

    ! cell-average => cell-avg-in-x and Gauss-point-in-y
    do n=1,NVAR
       do i=lo(1)-3,hi(1)+3
          call cellavg2gausspt_1d(lo(2),hi(2), U(i,:,n), Ulo(2), Uhi(2), &
               U1(i,:,n), U2(i,:,n), g3lo(2), g3hi(2))
       end do
    end do
    
    do g=1,2
       
       if (g .eq. 1) then
          Uag => U1
       else
          Uag => U2
       end if

       do n=1,NVAR
          do j=lo(2),hi(2)
             ! cell-avg-in-x and Gauss-point-in-y => xface and Gauss-point-in-y
             call cellavg2face_1d(lo(1),hi(1)+1, Uag(:,j,n),g3lo(1),g3hi(1), &
                  Qf(:,j,n),Qflo(1),Qfhi(1))

             ! cell-avg-in-x and Gauss-point-in-y => cell-center-in-x and Gauss-point-in-y
             tlo(1) = lo(1)-2
             tlo(2) = j
             thi(1) = hi(1)+2
             thi(2) = j
             call cellavg2cc(tlo(1:2),thi(1:2), Uag(:,:,n),g3lo,g3hi, &
                  Qc(:,:,n),Qclo(1:2),Qchi(1:2))
          end do
       end do

       tlo(1:2) = lo(1:2)
       thi(1) = hi(1)+1
       thi(2) = hi(2)
       call ctoprim(tlo,thi, Qf, Qflo,Qfhi,QFVAR)

       tlo(1) = lo(1)-2
       thi(1) = hi(1)+2
       call ctoprim(tlo,thi, Qc, Qclo,Qchi,QCVAR)

       ! transport coefficients on face
       call get_transport_properties(Qflo,Qfhi, Qf,Qflo,Qfhi,QFVAR, mu,xi,lam,Ddia,Qflo,Qfhi)

!       call comp_diff_flux_x()

       Nullify(Uag)
    end do

    ! ----- compute y-direction flux -----

  end subroutine difterm

end module difterm_module
