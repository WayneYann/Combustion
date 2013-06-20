module eos_module 

  use chemistry_module, only : nspecies, Ru, inv_mwt

  implicit none

  double precision, save, public :: Tref = 298.d0
  double precision, allocatable, save, public :: eref(:)

  double precision, save, private :: smalld = 1.d-50
  double precision, save, private :: smallt = 1.d-50
  double precision, save, private :: smallp = 1.d-50
  double precision, save, private :: smallc = 1.d-50

  logical, save, private :: initialized = .false.

contains

  subroutine eos_init(small_temp, small_dens, gamma_in, Tref_in)
    double precision, intent(in), optional :: small_temp
    double precision, intent(in), optional :: small_dens
    double precision, intent(in), optional :: gamma_in
    double precision, intent(in), optional :: Tref_in
    
    integer :: n, iwrk
    double precision :: Y(nspecies), rwrk

    if (present(small_temp)) then
       if (small_temp > 0.d0) then
          smallt = small_temp
       end if
    endif
 
    if (present(small_dens)) then
       if (small_dens > 0.d0) then
          smalld = small_dens
       end if
    endif

    if (present(Tref_in)) then
       Tref = Tref_in
    endif

    allocate(eref(nspecies))

    if (Tref .gt. 0.0) then
       do n=1,nspecies
          Y = 0.d0
          Y(n) = 1.d0
          call ckubms(Tref,Y,iwrk,rwrk,eref(n))
       end do
    else
       eref = 0.d0
    end if

    initialized = .true.
    
  end subroutine eos_init


  subroutine eos_get_small_temp(small_temp_out)
    
    double precision, intent(out) :: small_temp_out
 
    small_temp_out = smallt
 
  end subroutine eos_get_small_temp
 
  subroutine eos_get_small_dens(small_dens_out)
    
    double precision, intent(out) :: small_dens_out
    
    small_dens_out = smalld
    
  end subroutine eos_get_small_dens


  subroutine eos_get_c(c, rho, T, Y, pt_index)
    double precision, intent(out) :: c
    double precision, intent(in) :: rho, T, Y(nspecies)
    integer, optional, intent(in) :: pt_index(:)
    integer :: iwrk
    double precision :: rwrk, Cv, gamma, p, X(nspecies)

    call ckytx(Y, iwrk, rwrk, X)
    call ckcvbl(T, X, iwrk, rwrk, Cv)

    gamma = (Cv + Ru) / Cv

    call ckpy(rho, T, Y, iwrk, rwrk, p)

    c = sqrt(gamma*p/rho)
    c = max(c, smallc)
  end subroutine eos_get_c


  subroutine eos_get_T(T, e, Y, pt_index)
    double precision, intent(out) :: T
    double precision, intent(in) :: e, Y(nspecies)
    integer, optional, intent(in) :: pt_index(:)

    integer :: iwrk, ierr
    double precision :: rwrk
    call get_T_given_eY(e, Y, iwrk, rwrk, T, ierr)
    if (ierr .ne. 0) then
       print *, 'EOS: get_T failed, T, e, Y = ', T, e, Y
       call flush(6)
       call bl_error("Error: eos_get_T")
    end if
    T = max(T, smallt)
  end subroutine eos_get_T


  ! subroutine eos_get_p(p, rho, T, Y, pt_index)
  !   double precision, intent(out) :: p
  !   double precision, intent(in) :: rho, T, Y(nspecies)
  !   integer, optional, intent(in) :: pt_index(:)

  !   integer :: iwrk
  !   double precision :: rwrk

  !   call ckpy(rho, T, Y, iwrk, rwrk, p)
  !   p = max(p,smallp)

  ! end subroutine eos_get_p


  ! subroutine eos_get_TP(T, P, rho, e, Y, pt_index)
  !   double precision, intent(out) :: T, P
  !   double precision, intent(in ) :: rho, e, Y(nspecies)
  !   integer, optional, intent(in   ) :: pt_index(:)
    
  !   integer :: iwrk, ierr
  !   double precision :: rwrk

  !   call get_T_given_eY(e, Y, iwrk, rwrk, T, ierr)
  !   T = max(T, smallt)

  !   call ckpy(rho, T, Y, iwrk, rwrk, p)
  !   p = max(p,smallp)    

  ! end subroutine eos_get_TP


  ! subroutine eos_get_eP(e, P, rho, T, Y, pt_index)
  !   double precision, intent(  out) :: e, P
  !   double precision, intent(in   ) :: rho, T, Y(nspecies)
  !   integer, optional, intent(in   ) :: pt_index(:)
    
  !   integer :: iwrk
  !   double precision :: rwrk

  !   call ckpy(rho, T, Y, iwrk, rwrk, P)
  !   P = max(P,smallp)

  !   call ckubms(T,Y,iwrk,rwrk,e)

  ! end subroutine eos_get_eP


  subroutine eos_given_RTY(e, P, C, dpdr, dpde, rho, T, Y, pt_index)
    double precision, intent(  out) :: e, P, C, dpdr(nspecies), dpde
    double precision, intent(in   ) :: rho, T, Y(nspecies)
    integer, optional, intent(in   ) :: pt_index(:)
    
    integer :: iwrk
    double precision :: rwrk, Cv, G, X(nspecies)

    call ckpy(rho, T, Y, iwrk, rwrk, P)
    P = max(P,smallp)

    call ckubms(T,Y,iwrk,rwrk,e)

    call ckytx(Y, iwrk, rwrk, X)
    call ckcvbl(T, Y, iwrk, rwrk, Cv)
    G = (Cv + Ru) / Cv

    C = sqrt(G*P/rho)
    C = max(C, smallc)

    dpdr = Ru*T*inv_mwt
    dpde = (G-1.d0)*rho

  end subroutine eos_given_RTY


  subroutine eos_given_ReY(P, C, T, dpdr, dpde, rho, e, Y, pt_index)
    double precision, intent(  out) :: P, C, T, dpdr(nspecies), dpde
    double precision, intent(in   ) :: rho, e, Y(nspecies)
    integer, optional, intent(in  ) :: pt_index(:)

    integer :: iwrk, ierr
    double precision :: rwrk, Cv, G, X(nspecies)

    call get_T_given_eY(e, Y, iwrk, rwrk, T, ierr)

    if (ierr .ne. 0) then
       print *, 'EOS: get_T failed, T, e, Y = ', T, e, Y
       call flush(6)
       call bl_error("Error: eos_given_ReY")
    end if
    T = max(T, smallt)
    
    call ckytx(Y, iwrk, rwrk, X)
    call ckcvbl(T, X, iwrk, rwrk, Cv)
    G = (Cv + Ru) / Cv

    call ckpy(rho, T, Y, iwrk, rwrk, P)

    C = sqrt(G*P/rho)
    C = max(C, smallc)

    dpdr = Ru*T*inv_mwt
    dpde = (G-1.d0)*rho

  end subroutine eos_given_ReY


  subroutine eos_given_PTY(rho, e, Y, P, T, X, pt_index)
    double precision, intent(  out) :: rho, e, Y(nspecies)
    double precision, intent(in   ) :: P, T, X(nspecies)
    integer, optional, intent(in  ) :: pt_index(:)
    integer :: iwrk
    double precision :: rwrk
    call ckxty (X, iwrk, rwrk, Y)
    call ckrhoy(P,T,Y,iwrk,rwrk,rho)
    call ckubms(T   ,Y,iwrk,rwrk,e)
  end subroutine eos_given_PTY


  pure function eos_get_eref(Y) result(r)
    double precision, intent(in) :: Y(nspecies)
    double precision :: r
    r = dot_product(eref, Y)
  end function eos_get_eref

end module eos_module 
