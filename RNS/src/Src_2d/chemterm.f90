module chemterm_module

  use meth_params_module
  use burner_module, only : burn, compute_rhodYdt, splitburn, beburn
  use eos_module, only : eos_get_T
  use renorm_module, only : floor_species
  use passinfo_module, only : level

  implicit none

  private

  public :: chemterm, dUdt_chem

contains

  subroutine chemterm(lo, hi, U, Ulo, Uhi, st, stlo, sthi, dt, Up)
    integer, intent(in) :: lo(2), hi(2), Ulo(2), Uhi(2), stlo(2), sthi(2) 
    double precision, intent(inout) :: U(Ulo(1):Uhi(1),Ulo(2):Uhi(2),NVAR)
    double precision, intent(inout) :: st(stlo(1):sthi(1),stlo(2):sthi(2))
    double precision, intent(in) :: dt
    double precision, intent(in), optional :: Up(lo(1):hi(1),lo(2):hi(2),NVAR)

    select case (chem_solver)
       case (cc_burning)
          call chemterm_cellcenter(lo, hi, U, Ulo, Uhi, st, stlo, sthi, dt)
       case (Gauss_burning)
          call chemterm_gauss(lo, hi, U, Ulo, Uhi, dt)
       case (split_burning)
          call chemterm_split(lo, hi, U, Ulo, Uhi, dt)
       case (BEcc_burning)
          call chemterm_becc(lo, hi, U, Ulo, Uhi, st, stlo, sthi, dt, Up)
       case (BEGp_burning)
          call chemterm_begp(lo, hi, U, Ulo, Uhi, dt, Up)
       case default
          call bl_error("unknown chem_solver")
       end select

  end subroutine chemterm


  subroutine chemterm_cellcenter(lo, hi, U, Ulo, Uhi, dt)
    use convert_module, only : cellavg2cc_2d, cc2cellavg_2d
    integer, intent(in) :: lo(2), hi(2), Ulo(2), Uhi(2)
    double precision, intent(inout) :: U(Ulo(1):Uhi(1),Ulo(2):Uhi(2),NVAR)
    double precision, intent(in) :: dt

    integer :: i, j, n
    logical :: force_new_J
    double precision :: rhot(1), rhoinv, ei
    double precision :: Yt(nspec+1)
    double precision, allocatable :: Ucc(:,:,:)

    allocate(Ucc(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,NVAR))

    !$omp parallel private(i,j,n,rhot,rhoinv,ei,Yt,force_new_J)

    !$omp do
    do n=1,NVAR
       call cellavg2cc_2d(lo-1,hi+1, U(:,:,n), Ulo,Uhi, Ucc(:,:,n), lo-1,hi+1)
    end do
    !$omp end do

    force_new_J = .true.  ! always recompute Jacobina when a new FAB starts

    !$omp do
    do j=lo(2)-1,hi(2)+1
       do i=lo(1)-1,hi(1)+1

          rhot = 0.d0
          do n=1,NSPEC
             Yt(n) = Ucc(i,j,UFS+n-1)
             rhot(1) = rhot(1) + Yt(n)
          end do
          rhoinv = 1.d0/rhot(1)

          Yt(1:nspec) = Yt(1:nspec) * rhoinv
          Yt(nspec+1) = Ucc(i,j,UTEMP)

          ei = rhoinv*( Ucc(i,j,UEDEN) - 0.5d0*rhoinv*(Ucc(i,j,UMX)**2 &
               + Ucc(i,j,UMY)**2) )

          call eos_get_T(Yt(nspec+1), ei, Yt(1:nspec))

          call burn(1, rhot, Yt, dt, force_new_J)

          force_new_J = new_J_cell

          do n=1,nspec
             Ucc(i,j,UFS+n-1) = rhot(1)*Yt(n)
          end do

       end do
    end do
    !$omp end do

    !$omp do
    do n=UFS,UFS+nspec-1
       call cc2cellavg_2d(lo,hi, Ucc(:,:,n), lo-1,hi+1, U(:,:,n), Ulo,Uhi)
    end do
    !$omp end do

    !$omp end parallel

    deallocate(Ucc)

  end subroutine chemterm_cellcenter


  subroutine chemterm_gauss(lo, hi, U, Ulo, Uhi, dt)
    use weno_module, only : cellavg2gausspt_2d
    integer, intent(in) :: lo(2), hi(2), Ulo(2), Uhi(2)
    double precision, intent(inout) :: U(Ulo(1):Uhi(1),Ulo(2):Uhi(2),NVAR)
    double precision, intent(in) :: dt

    integer :: i, j, n, g
    logical :: force_new_J
    double precision :: rhot(4), rhoinv, ei
    double precision :: Yt(nspec+1,4)
    double precision, allocatable :: UG(:,:,:,:)

    allocate(UG(lo(1):hi(1),lo(2):hi(2),4,NVAR))

    !$omp parallel private(i,j,n,g,rhot,rhoinv,ei,Yt,force_new_J)

    !$omp do
    do n=1,NVAR
       call cellavg2gausspt_2d(lo,hi, U(:,:,n), Ulo,Uhi, UG(:,:,:,n), lo,hi)
    end do
    !$omp end do

    force_new_J = .true.  ! always recompute Jacobian when a new FAB starts

    !$omp do
    do j=lo(2),hi(2)
       do i=lo(1),hi(1)

          do g=1,4

             rhot(g) = 0.d0
             do n=1,NSPEC
                Yt(n,g) = UG(i,j,g,UFS+n-1)
                rhot(g) = rhot(g) + Yt(n,g)
             end do
             rhoinv = 1.d0/rhot(g)

             Yt(1:nspec,g) = Yt(1:nspec,g) * rhoinv
             Yt(nspec+1,g) = UG(i,j,g,UTEMP)

             ei = rhoinv*( UG(i,j,g,UEDEN) - 0.5d0*rhoinv*(UG(i,j,g,UMX)**2 &
                  + UG(i,j,g,UMY)**2) )

             call eos_get_T(Yt(nspec+1,g), ei, Yt(1:nspec,g))

          end do

          call burn(4, rhot, Yt, dt, force_new_J)

          force_new_J = new_J_cell

          U(i,j,UFS:UFS+nspec-1) = 0.d0 
          do g=1,4
             do n=1,nspec
                U(i,j,UFS+n-1) = U(i,j,UFS+n-1) + 0.25d0*rhot(g)*Yt(n,g)
             end do
          end do

       end do
    end do
    !$omp end do
    !$omp end parallel

    deallocate(UG)

  end subroutine chemterm_gauss


  subroutine chemterm_split(lo, hi, U, Ulo, Uhi, dt)
    use weno_module, only : cellavg2gausspt_2d
    integer, intent(in) :: lo(2), hi(2), Ulo(2), Uhi(2)
    double precision, intent(inout) :: U(Ulo(1):Uhi(1),Ulo(2):Uhi(2),NVAR)
    double precision, intent(in) :: dt

    integer :: i, j, n, g
    logical :: force_new_J
    double precision :: rhot(4), rhoinv, ei, rho0(1)
    double precision :: Yt(nspec+1,4), Y0(nspec+1)
    double precision, allocatable :: UG(:,:,:,:)

    allocate(UG(lo(1):hi(1),lo(2):hi(2),4,NVAR))

    !$omp parallel private(i,j,n,g,rhot,rhoinv,ei,Yt,force_new_J,rho0,Y0)

    !$omp do
    do n=1,NVAR
       call cellavg2gausspt_2d(lo,hi, U(:,:,n), Ulo,Uhi, UG(:,:,:,n), lo,hi)
    end do
    !$omp end do

    force_new_J = .true.  ! always recompute Jacobina when a new FAB starts

    !$omp do
    do j=lo(2),hi(2)
       do i=lo(1),hi(1)

          Y0 = 0.d0
          rho0(1) = 0.d0

          do g=1,4

             rhot(g) = 0.d0
             do n=1,NSPEC
                Yt(n,g) = UG(i,j,g,UFS+n-1)
                rhot(g) = rhot(g) + Yt(n,g)
             end do
             rhoinv = 1.d0/rhot(g)

             Yt(1:nspec,g) = Yt(1:nspec,g) * rhoinv
             Yt(nspec+1,g) = UG(i,j,g,UTEMP)

             ei = rhoinv*( UG(i,j,g,UEDEN) - 0.5d0*rhoinv*(UG(i,j,g,UMX)**2 &
                  + UG(i,j,g,UMY)**2) )

             call eos_get_T(Yt(nspec+1,g), ei, Yt(1:nspec,g))

             Y0 = Y0 + 0.25d0*Yt(:,g)
             rho0(1) = rho0(1) + 0.25d0*rhot(g)

          end do

          call burn(1, rho0(1), Y0, dt, force_new_J)

          force_new_J = new_J_cell

          call splitburn(4, rho0(1), Y0, rhot, Yt, dt) 
          ! Now Yt is \Delta Y and T

          U(i,j,UFS:UFS+nspec-1) = 0.d0 
          do g=1,4
             do n=1,nspec
                U(i,j,UFS+n-1) = U(i,j,UFS+n-1) +  &
                     0.25d0*(UG(i,j,g,UFS+n-1) + rhot(g)*Yt(n,g))
             end do
          end do

       end do
    end do
    !$omp end do
    !$omp end parallel

    deallocate(UG)

  end subroutine chemterm_split


  subroutine chemterm_begp(lo, hi, U, Ulo, Uhi, dt, Up)
    use weno_module, only : cellavg2gausspt_2d
    integer, intent(in) :: lo(2), hi(2), Ulo(2), Uhi(2)
    double precision, intent(inout) :: U(Ulo(1):Uhi(1),Ulo(2):Uhi(2),NVAR)
    double precision, intent(in) :: dt
    double precision, intent(in), optional :: Up(lo(1):hi(1),lo(2):hi(2),NVAR)

    integer :: i, j, n, g, ierr
    logical :: force_new_J
    double precision :: rhot(4), rhoinv, rho0(1)
    double precision :: Yt(nspec+1,4), Y0(nspec+1), rhoY(nspec)
    double precision, allocatable :: UG(:,:,:,:)

    allocate(UG(lo(1):hi(1),lo(2):hi(2),4,NVAR))

    if (chem_do_weno) then
       call chem_weno(lo, hi, U, Ulo, Uhi, UG)
    else
       !$omp parallel do private(n)
       do n=1,NVAR
          call cellavg2gausspt_2d(lo,hi, U(:,:,n), Ulo,Uhi, UG(:,:,:,n), lo,hi)
       end do
       !$omp end parallel do
    end if

    !$omp parallel private(i,j,n,g,ierr,rhot,rhoinv,Yt,force_new_J,rho0,Y0,rhoY)

    force_new_J = .true.  ! always recompute Jacobina when a new FAB starts

    !$omp do schedule(dynamic)
    do j=lo(2),hi(2)
       do i=lo(1),hi(1)

          do g=1,4
             call get_rhoYT(UG(i,j,g,:), rhot(g), YT(1:nspec,g), YT(nspec+1,g), ierr)
             if (ierr .ne. 0) then
                print *, 'chemterm_be: eos_get_T failed for UG at ', i,j,g,UG(i,j,g,:)
                call bl_error("chemterm_be failed at eos_get_T")
             end if
          end do

          if (present(Up)) then

             rho0(1) = Up(i,j,URHO)
             rhoinv = 1.d0/rho0(1)
             Y0(1:nspec) = Up(i,j,UFS:UFS+nspec-1)*rhoinv
             Y0(nspec+1) = Up(i,j,UTEMP)

          else

             call get_rhoYT(U(i,j,:), rho0(1), Y0(1:nspec), Y0(nspec+1), ierr)
             if (ierr .ne. 0) then
                print *, 'chemterm_be: eos_get_T failed for U at ', i,j,U(i,j,:)
                call bl_error("chemterm_be failed at eos_get_T for U")
             end if

             call burn(1, rho0(1), Y0, dt, force_new_J, ierr)
             force_new_J = new_J_cell
             if (ierr .ne. 0) then
                print *, 'chemterm_be: burn failed at ', i,j,U(i,j,:)
                print *, '   rho0, Y0 =', rho0(1), Y0
                call bl_error("chemterm_be failed at burn")
             end if

          end if

          call renorm(nspec, Y0(1:nspec), ierr)
          if (ierr .ne. 0) then
             call bl_error("chemterm_be failed at renormalizing Y0")
          end if

          rhoY = 0.d0
          do g=1,4
             call beburn(rho0(1), Y0, rhot(g), Yt(:,g), dt, g, ierr)
             if (ierr .ne. 0) then ! beburn failed
                print *, 'chemterm_be: beburn failed at ',i,j,g,UG(i,j,g,:)
                call bl_error("chemterm_be: beburn failed at g")
             end if
             do n=1,nspec
                rhoY(n) = rhoY(n) + 0.25d0*rhot(g)*Yt(n,g)
             end do
          end do

          U(i,j,UFS:UFS+nspec-1) = rhoY

       end do
    end do
    !$omp end do
    !$omp end parallel

    deallocate(UG)

  end subroutine chemterm_begp


  subroutine dUdt_chem(lo, hi, U, Ulo, Uhi, Ut, Utlo, Uthi)
    use weno_module, only : cellavg2gausspt_2d
    integer, intent(in) :: lo(2), hi(2), Ulo(2), Uhi(2), Utlo(2), Uthi(2)
    double precision, intent(in ) ::  U( Ulo(1): Uhi(1), Ulo(2): Uhi(2),NVAR)
    double precision, intent(out) :: Ut(Utlo(1):Uthi(1),Utlo(2):Uthi(2),NVAR)

    integer :: i, j, n, g, np, ierr
    double precision :: rho(lo(1):hi(1)), T(lo(1):hi(1))
    double precision :: Y(lo(1):hi(1),nspec), rdYdt(lo(1):hi(1),nspec)
    double precision, allocatable :: UG(:,:,:,:)

    np = hi(1)-lo(1)+1

    allocate(UG(lo(1):hi(1),lo(2):hi(2),4,NVAR))

    if (chem_do_weno) then
       call chem_weno(lo, hi, U, Ulo, Uhi, UG)
    else
       !$omp parallel do private(n)
       do n=1,NVAR
          call cellavg2gausspt_2d(lo,hi, U(:,:,n), Ulo,Uhi, UG(:,:,:,n), lo,hi)
       end do
       !$omp end parallel do
    end if

    !$omp parallel private(i,j,n,g,rho,T,Y,rdYdt)

    !$omp do
    do n=1,NVAR
       do j=lo(2),hi(2)
       do i=lo(1),hi(1)
          Ut(i,j,n) = 0.d0
       end do
       end do
    end do
    !$omp end do

    !$omp do
    do j=lo(2),hi(2)

       do g=1,4
          
          do i=lo(1),hi(1)
             call get_rhoYT(UG(i,j,g,:), rho(i), Y(i,:), T(i), ierr)
             if (ierr .ne. 0) then
                print *, 'dUdt_chem: eos_get_T failed for UG at ', i,j,g,UG(i,j,g,:)
                call bl_error("dUdt_chem failed at eos_get_T for UG")
             end if
          end do

          call compute_rhodYdt(np,rho,T,Y,rdYdt)

          do n=1,nspec
             do i=lo(1),hi(1)
                Ut(i,j,UFS+n-1) = Ut(i,j,UFS+n-1) + 0.25d0*rdYdt(i,n)
             end do
          end do

       end do
    end do
    !$omp end do

    !$omp end parallel

    deallocate(UG)

  end subroutine dUdt_chem


  subroutine get_rhoYT(U, rho, Y, T, ierr)
    double precision, intent(in) :: U(NVAR)
    double precision, intent(out) :: rho, Y(nspec), T
    integer, intent(out) :: ierr

    integer :: n
    double precision :: rhoinv, ei
    
    rho = 0.d0
    do n=1,NSPEC
       Y(n) = U(UFS+n-1)
       rho = rho + Y(n)
    end do
    rhoinv = 1.d0/rho
    
    Y = Y * rhoinv
    T = U(UTEMP)
    
    ei = rhoinv*( U(UEDEN) - 0.5d0*rhoinv*(U(UMX)**2 + U(UMY)**2) )
    
    call eos_get_T(T, ei, Y, ierr=ierr)
  end subroutine get_rhoYT


  subroutine chem_weno(lo, hi, U, Ulo, Uhi, UG)
    use reconstruct_module, only : reconstruct_comp
    integer, intent(in) :: lo(2), hi(2), Ulo(2), Uhi(2)
    double precision, intent(in) :: U(Ulo(1):Uhi(1),Ulo(2):Uhi(2),NVAR)
    double precision, intent(out):: UG(lo(1):hi(1),lo(2):hi(2),4,NVAR)

    integer :: i,j
    double precision, allocatable :: UG1y(:,:,:), UG2y(:,:,:)

    allocate(UG1y(lo(1)-2:hi(1)+2, lo(2):hi(2), NVAR))
    allocate(UG2y(lo(1)-2:hi(1)+2, lo(2):hi(2), NVAR))

    !$omp parallel private(i,j)

    !$omp do
    do i=lo(1)-2,hi(1)+2
       call reconstruct_comp(lo(2),hi(2), &
            Ulo(2),Uhi(2),     &  ! for U
            0, 0,              &  ! for UL & UR, not present
            lo(2)  ,hi(2),     &  ! for UG1 & UG2
            U(i,:,:), &
            UG1=UG1y(i,:,:), UG2=UG2y(i,:,:) )
    end do
    !$omp end do
    
    !$omp do
    do j=lo(2), hi(2)
       call reconstruct_comp(lo(1),hi(1), &
            lo(1)-2,hi(1)+2,   &  ! for U
            0, 0,              &  ! for UL & UR, not present
            lo(1)  ,hi(1),     &  ! for UG1 & UG2
            UG1y(:,j,:), &
            UG1=UG(:,j,1,:), UG2=UG(:,j,2,:) )
       call reconstruct_comp(lo(1),hi(1), &
            lo(1)-2,hi(1)+2,   &  ! for U
            0, 0,              &  ! for UL & UR, not present
            lo(1)  ,hi(1),     &  ! for UG1 & UG2
            UG2y(:,j,:), &
            UG1=UG(:,j,3,:), UG2=UG(:,j,4,:) )
    end do
    !$omp end do 

    !$omp end parallel

    deallocate(UG1y,UG2y)

  end subroutine chem_weno

end module chemterm_module
