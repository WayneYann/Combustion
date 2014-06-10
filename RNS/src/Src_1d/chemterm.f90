module chemterm_module

  use meth_params_module
  use burner_module, only : burn, compute_rhodYdt, splitburn, beburn
  use eos_module, only : eos_get_T
  use weno_module, only : cellavg2gausspt_1d
  use convert_module, only : cellavg2cc_1d, cc2cellavg_1d
  use renorm_module, only : floor_species
  use passinfo_module, only : level

  implicit none

  private

  public :: chemterm, dUdt_chem

contains

  subroutine chemterm(lo, hi, U, Ulo, Uhi, st, stlo, sthi, dt, Up)
    integer, intent(in) :: lo(1), hi(1), Ulo(1), Uhi(1), stlo(1), sthi(1)
    double precision, intent(inout) :: U(Ulo(1):Uhi(1),NVAR)
    double precision, intent(inout) :: st(stlo(1):sthi(1))
    double precision, intent(in) :: dt
    double precision, intent(in), optional :: Up(lo(1):hi(1),NVAR)

    select case (chem_solver)
       case (cc_burning)
          call chemterm_cellcenter(lo, hi, U, Ulo, Uhi, st, stlo, sthi, dt)
       case (Gauss_burning)
          call chemterm_gauss(lo, hi, U, Ulo, Uhi, st, stlo, sthi, dt)
       case (split_burning)
          call chemterm_split(lo, hi, U, Ulo, Uhi, st, stlo, sthi, dt)
       case (BEcc_burning)
          call chemterm_becc(lo, hi, U, Ulo, Uhi, st, stlo, sthi, dt, Up)
       case (BEGp_burning)
          call chemterm_begp(lo, hi, U, Ulo, Uhi, st, stlo, sthi, dt, Up)
       case default
          call bl_error("unknown chem_solver")
       end select

  end subroutine chemterm


  subroutine dUdt_chem(lo, hi, U, Ulo, Uhi, Ut, Utlo, Uthi, st, stlo, sthi)
    integer, intent(in) :: lo(1), hi(1), Ulo(1), Uhi(1), Utlo(1), Uthi(1), stlo(1), sthi(1)
    double precision, intent(in)  ::  U( Ulo(1): Uhi(1),NVAR)
    double precision, intent(out) :: Ut(Utlo(1):Uthi(1),NVAR)
    double precision, intent(inout) :: st(stlo(1):sthi(1))
    if (chem_solver .eq. cc_burning .or. chem_solver .eq. BEcc_burning) then
       call dUdt_chem_cellcenter(lo, hi, U, Ulo, Uhi, Ut, Utlo, Uthi, st, stlo, sthi)
    else
       call dUdt_chem_gauss(lo, hi, U, Ulo, Uhi, Ut, Utlo, Uthi, st, stlo, sthi)
    end if
  end subroutine dUdt_chem


  subroutine chemterm_cellcenter(lo, hi, U, Ulo, Uhi, st, stlo, sthi, dt)
    integer, intent(in) :: lo(1), hi(1), Ulo(1), Uhi(1), stlo(1), sthi(1)
    double precision, intent(inout) :: U(Ulo(1):Uhi(1),NVAR)
    double precision, intent(inout) :: st(stlo(1):sthi(1))
    double precision, intent(in) :: dt

    integer :: i, n, ierr
    logical :: force_new_J
    double precision :: rhot(1), Yt(nspec+1)
    double precision, allocatable :: Ucc(:,:)

    allocate(Ucc(lo(1)-1:hi(1)+1,NVAR))

    do n=1,NVAR
       call cellavg2cc_1d(lo(1)-1,hi(1)+1, U(:,n), Ulo(1),Uhi(1), Ucc(:,n), lo(1)-1,hi(1)+1)
    end do

    force_new_J = .true.  ! always recompute Jacobina when a new FAB starts

    do i=lo(1)-1,hi(1)+1

       if (st(i) .eq. 0.d0) then
          call get_rhoYT(Ucc(i,:), rhot(1), YT(1:nspec), YT(nspec+1), ierr)

          if (ierr .ne. 0) then
             st(i) = -1.d0
             force_new_J = .true.
          else
             call burn(1, rhot, YT, dt, force_new_J, ierr)
             if (ierr .ne. 0) then
                st(i) = -1.d0
                force_new_J = .true.
             else
                force_new_J = new_J_cell
             end if
          end if
       end if
          
       if (st(i) .ne. 0.d0) then ! burn cell average instead
          call get_rhoYT(U(i,:), rhot(1), YT(1:nspec), YT(nspec+1), ierr)
          if (ierr .ne. 0) then
             print *, 'chemterm_cellcenter: eos_get_T failed for U at ', &
                  level, i,U(i,:)
             call bl_error("chemterm_cellcenter failed at eos_get_T")
          end if

          force_new_J = .true.
          call burn(1, rhot, YT, dt, force_new_J, ierr)
          force_new_J = new_J_cell
          if (ierr .ne. 0) then
             print *, 'chemterm_cellcenter: bdf burn failed for U at ', &
                  level,i,U(i,:)
             call bl_error("chemterm_cellcenter failed at bdf burn for U")
          end if
       end if

       do n=1,nspec
          Ucc(i,UFS+n-1) = rhot(1)*YT(n)
       end do

    end do

    do n=UFS,UFS+nspec-1
       call cc2cellavg_1d(lo(1),hi(1), Ucc(:,n), lo(1)-1,hi(1)+1, U(:,n), Ulo(1),Uhi(1))
    end do

    deallocate(Ucc)

  end subroutine chemterm_cellcenter


  subroutine chemterm_gauss(lo, hi, U, Ulo, Uhi, st, stlo, sthi, dt)
    integer, intent(in) :: lo(1), hi(1), Ulo(1), Uhi(1), stlo(1), sthi(1)
    double precision, intent(inout) :: U(Ulo(1):Uhi(1),NVAR)
    double precision, intent(inout) :: st(stlo(1):sthi(1))
    double precision, intent(in) :: dt

    integer :: i, n, g, ierr
    logical :: force_new_J
    double precision :: rhot(2), YT(nspec+1,2)
    double precision, allocatable :: UG(:,:,:)

    allocate(UG(lo(1):hi(1),NVAR,2))

    if (chem_do_weno) then
       call chem_weno(lo(1), hi(1), U, Ulo(1), Uhi(1), UG)
    else
       do n=1,NVAR
          call cellavg2gausspt_1d(lo(1),hi(1), U(:,n), Ulo(1),Uhi(1), &
               UG(:,n,1), UG(:,n,2), lo(1),hi(1))
       end do
    end if
       
    force_new_J = .true.  ! always recompute Jacobian when a new FAB starts
    
    do i=lo(1),hi(1)
       if (st(i) .eq. 0.d0) then
          do g=1,2
             call get_rhoYT(UG(i,:,g), rhot(g), YT(1:nspec,g), YT(nspec+1,g), ierr)
             if (ierr .ne. 0) then
                st(i) = -1.d0
                force_new_J = .true.
                exit
             end if
          end do
       end if

       if (st(i) .eq. 0.d0) then
          call burn(2, rhot, Yt, dt, force_new_J)
          if (ierr .ne. 0) then
             st(i) = -1.d0
             force_new_J = .true.
          else
             force_new_J = new_J_cell
          end if
       end if
       
       if (st(i) .eq. 0.d0) then
          
          U(i,UFS:UFS+nspec-1) = 0.d0 
          do g=1,2
             do n=1,nspec
                U(i,UFS+n-1) = U(i,UFS+n-1) + 0.5d0*rhot(g)*Yt(n,g)
             end do
          end do
             
       else ! burn cell average instead

          call get_rhoYT(U(i,:), rhot(1), YT(1:nspec,1), YT(nspec+1,1), ierr)
          if (ierr .ne. 0) then
             print *, 'chemterm_gauss: eos_get_T failed for U at ', level,i,U(i,:)
             call bl_error("chemterm_gauss failed at eos_get_T")
          end if
          
          force_new_J = .true.
          call burn(1, rhot(1:1), YT(:,1), dt, force_new_J, ierr)
          force_new_J = new_J_cell
          if (ierr .ne. 0) then
             print *, "chemterm_gauss: bdf burn failed for U at ", level,i,U(i,:)
             call bl_error("chemterm_gauss failed at bdf burn for U")
          end if
             
          do n=1,nspec
             U(i,UFS+n-1) = rhot(1)*YT(n,1)
          end do

       end if

    end do

    deallocate(UG)

  end subroutine chemterm_gauss


  subroutine chemterm_split(lo, hi, U, Ulo, Uhi, st, stlo, sthi, dt)
    integer, intent(in) :: lo(1), hi(1), Ulo(1), Uhi(1), stlo(1), sthi(1)
    double precision, intent(inout) :: U(Ulo(1):Uhi(1),NVAR)
    double precision, intent(inout) :: st(stlo(1):sthi(1))
    double precision, intent(in) :: dt

    integer :: i, n, g, ierr
    logical :: force_new_J
    double precision :: rhot(2), rho0(1)
    double precision :: Yt(nspec+1,2), YT0(nspec+1)
    double precision, allocatable :: UG(:,:,:)

    allocate(UG(lo(1):hi(1),NVAR,2))

    if (chem_do_weno) then
       call chem_weno(lo(1),hi(1),U,Ulo(1),Uhi(1),UG)
    else
       do n=1,NVAR
          call cellavg2gausspt_1d(lo(1),hi(1), U(:,n), Ulo(1),Uhi(1), &
               UG(:,n,1), UG(:,n,2), lo(1),hi(1))
       end do
    end if

    force_new_J = .true.  ! always recompute Jacobina when a new FAB starts

    do i=lo(1),hi(1)
       
       if (st(i) .eq. 0.d0) then
          YT0 = 0.d0
          rho0(1) = 0.d0
          do g=1,2
             call get_rhoYT(UG(i,:,g), rhot(g), YT(1:nspec,g), YT(nspec+1,g), ierr)
             if (ierr .ne. 0) then
                st(i) = -1.d0
                force_new_J = .true.
                exit
             end if
             YT0 = YT0 + 0.5d0*YT(:,g)
             rho0(1) = rho0(1) + 0.5d0*rhot(g)
          end do
       end if

       if (st(i) .eq. 0.d0) then
          call burn(1, rho0, YT0, dt, force_new_J, ierr)
          if (ierr .ne. 0) then
             st(i) = -1.d0
             force_new_J = .true.
          else
             force_new_J = new_J_cell
          end if
       end if

       if (st(i) .eq. 0.d0) then

          call splitburn(2, rho0(1), YT0, rhot, Yt, dt) 
          ! Now Yt is \Delta Y and T

          U(i,UFS:UFS+nspec-1) = 0.d0 
          do g=1,2
             do n=1,nspec
                U(i,UFS+n-1) = U(i,UFS+n-1) +  &
                     0.5d0*(UG(i,g,UFS+n-1) + rhot(g)*Yt(n,g))
             end do
          end do

       else ! burn cell average instead
          
          call get_rhoYT(U(i,:), rhot(1), YT(1:nspec,1), YT(nspec+1,1), ierr)
          if (ierr .ne. 0) then
             print *, "chemterm_split: eos_get_T faile for U at ", level,i,U(i,:)
             call bl_error("chemterm_split failed at eos_get_T")
          end if

          call burn(1, rhot(1:1), YT(:,1), dt, force_new_J, ierr)
          force_new_J = new_J_cell
          if (ierr .ne. 0) then
             print *, "chemterm_split: bdf burn failed for U at ", level, i,U(i,:)
             call bl_error("chemterm_split failed at bdf burn for U")
          end if
          
          do n=1,nspec
             U(i,UFS+n-1) = rhot(1)*YT(n,1)
          end do

       end if

    end do

    deallocate(UG)

  end subroutine chemterm_split


  subroutine chemterm_becc(lo, hi, U, Ulo, Uhi, st, stlo, sthi, dt, Up)
    integer, intent(in) :: lo(1), hi(1), Ulo(1), Uhi(1), stlo(1), sthi(1)
    double precision, intent(inout) :: U(Ulo(1):Uhi(1),NVAR)
    double precision, intent(inout) :: st(stlo(1):sthi(1))
    double precision, intent(in) :: dt
    double precision, intent(in), optional :: Up(lo(1):hi(1),NVAR)

    integer :: i, n, ierr
    logical :: force_new_J
    double precision :: rho0(1), rhot(1), YT0(nspec+1), YT(nspec+1)
    double precision, allocatable :: Ucc(:,:)

    allocate(Ucc(lo(1)-1:hi(1)+1,NVAR))

    do n=1,NVAR
       call cellavg2cc_1d(lo(1)-1,hi(1)+1, U(:,n), Ulo(1),Uhi(1), Ucc(:,n), lo(1)-1,hi(1)+1)
    end do

    force_new_J = .true.

    do i=lo(1)-1,hi(1)+1

       if (st(i) .eq. 0.d0) then
             
          call get_rhoYT(Ucc(i,:), rhot(1), YT(1:nspec), YT(nspec+1), ierr)
             
          if (ierr .ne. 0) then
             st(i) = -1.d0
             force_new_J = .true.
          else
             
             if (present(Up) .and. i.ge.lo(1) .and. i.le.hi(1)) then

                rho0(1) = Up(i,URHO)
                YT0(1:nspec) = Up(i,UFS:UFS+nspec-1)/rho0(1)
                YT0(nspec+1) = Up(i,UTEMP)

             else

                rho0(1) = rhot(1)
                YT0 = YT

                call burn(1,rho0, YT0, dt, force_new_J, ierr)
                if (ierr .ne. 0) then
                   st(i) = -1.d0
                   force_new_J = .true.
                else
                   force_new_J = new_J_cell
                end if
                
             end if
             
          end if

       end if

       if (st(i) .ne. 0.d0) then ! burn cell average
             
          call get_rhoYT(U(i,:), rhot(1), YT(1:nspec), YT(nspec+1), ierr)
          if (ierr .ne. 0) then
             print *, "chemterm_becc: eos_get_T failed for U at ", level, i,U(i,:)
             call bl_error("chemterm_becc failed at eos_get_T")
          end if
             
          if (present(Up) .and. i.ge.lo(1) .and. i.le.hi(1)) then
                
             rho0(1) = Up(i,URHO)
             YT0(1:nspec) = Up(i,UFS:UFS+nspec-1)/rho0(1)
             YT0(nspec+1) = Up(i,UTEMP)
                
          else
                
             rho0(1) = rhot(1)
             YT0 = YT
                
             force_new_J = .true.
             call burn(1, rho0, YT0, dt, force_new_J, ierr)
             force_new_J = new_J_cell
             if (ierr .ne. 0) then
                print *, "chemterm_becc: bdf burn failed for U at ", level,i,U(i,:)
                call bl_error("chemterm_becc failed at bdf burn for U")
             end if
                
          end if
             
       end if
          
       call floor_species(nspec, YT0(1:nspec))
          
       call beburn(rho0(1), YT0, rhot(1), YT, dt, 1, ierr)
       if (ierr .ne. 0) then 
          print *, "chemterm_becc: beburn failed at ", level,i,U(i,:)
          call bl_error("chemterm_becc: beburn failed")
       end if
          
       do n=1,nspec
          Ucc(i,UFS+n-1) = rhot(1)*YT(n)
       end do
          
    end do

    do n=UFS,UFS+nspec-1
       call cc2cellavg_1d(lo(1),hi(1), Ucc(:,n), lo(1)-1,hi(1)+1, U(:,n), Ulo(1),Uhi(1))
    end do

    deallocate(Ucc)

  end subroutine chemterm_becc


  subroutine chemterm_begp(lo, hi, U, Ulo, Uhi, st, stlo, sthi, dt, Up)
    integer, intent(in) :: lo(1), hi(1), Ulo(1), Uhi(1), stlo(1), sthi(1)
    double precision, intent(inout) :: U(Ulo(1):Uhi(1),NVAR)
    double precision, intent(inout) :: st(stlo(1):sthi(1))
    double precision, intent(in) :: dt
    double precision, intent(in), optional :: Up(lo(1):hi(1),NVAR)

    integer :: i, n, g, ierr
    logical :: force_new_J
    double precision :: rhot(2), rho0(1)
    double precision :: Yt(nspec+1,2), YT0(nspec+1), rhoY(nspec)
    double precision, allocatable :: UG(:,:,:)

    allocate(UG(lo(1):hi(1),NVAR,2))

    if (chem_do_weno) then
       call chem_weno(lo(1), hi(1), U, Ulo(1), Uhi(1), UG)
    else
       do n=1,NVAR
          call cellavg2gausspt_1d(lo(1),hi(1), U(:,n), Ulo(1),Uhi(1), &
               UG(:,n,1), UG(:,n,2), lo(1),hi(1))
       end do
    end if

    force_new_J = .true.  ! always recompute Jacobina when a new FAB starts

    do i=lo(1),hi(1)

       if (st(i).eq.0.d0) then
          do g=1,2
             call get_rhoYT(UG(i,:,g), rhot(g), YT(1:nspec,g), YT(nspec+1,g), ierr)
             if (ierr .ne. 0) then
                force_new_J = .true.
                st(i) = -1.d0
                exit
             end if
          end do
       end if
       
       if (st(i) .eq. 0.d0) then
          
          if (present(Up)) then
             rho0(1) = Up(i,URHO)
             YT0(1:nspec) = Up(i,UFS:UFS+nspec-1)/rho0(1)
             YT0(nspec+1) = Up(i,UTEMP)
          else
             call get_rhoYT(U(i,:), rho0(1), YT0(1:nspec), YT0(nspec+1), ierr)
             if (ierr .ne. 0) then
                print *, 'chemterm_begp: eos_get_T failed for U at ', level,i,U(i,:)
                call bl_error("chemterm_begp failed at eos_get_T for U")
             end if
             
             call burn(1, rho0(1), YT0, dt, force_new_J, ierr)
             force_new_J = new_J_cell
             if (ierr .ne. 0) then
                print *, 'chemterm_begp: bdf burn failed at ', i,U(i,:)
                print *, '   rho0, YT0 =', rho0(1), YT0
                call bl_error("chemterm_begp failed at bdf burn")
             end if
          end if

          call floor_species(nspec, YT0(1:nspec))

          rhoY = 0.d0
          do g=1,2
             call beburn(rho0(1), YT0, rhot(g), Yt(:,g), dt, g, ierr)
             if (ierr .ne. 0) then ! beburn failed
                st(i) = -1.d0
                force_new_J = .true.
                exit
             end if
             do n=1,nspec
                rhoY(n) = rhoY(n) + 0.5d0*rhot(g)*Yt(n,g)
             end do
          end do
          
       end if

       if (st(i) .eq. 0.d0) then

          U(i,UFS:UFS+nspec-1) = rhoY

       else

          call get_rhoYT(U(i,:), rhot(1), YT(1:nspec,1), YT(nspec+1,1), ierr)
          if (ierr .ne. 0) then
             print *, "chemterm_begp: eos_get_T failed for U at ", level,i,U(i,:)
             call bl_error("chemterm_begp failed at eos_get_T")
          end if
          
          if (present(Up)) then
             
             rho0(1) = Up(i,URHO)
             YT0(1:nspec) = Up(i,UFS:UFS+nspec-1)/rho0(1)
             YT0(nspec+1) = Up(i,UTEMP)
             
          else
             
             rho0(1) = rhot(1)
             YT0 = YT(:,1)
             
             force_new_J = .true.
             call burn(1, rho0, YT0, dt, force_new_J, ierr)
             force_new_J = new_J_cell
             if (ierr .ne. 0) then
                print *, "chemterm_begp: bdf burn failed for U at ", level, i,U(i,:)
                call bl_error("chemterm_begp failed at bdf burn for U")
             end if
             
          end if
          
          call floor_species(nspec, YT0(1:nspec))
          
          call beburn(rho0(1), YT0, rhot(1), YT(:,1), dt, 1, ierr)
          if (ierr .ne. 0) then
             print *, "chemterm_begp: beburn failed for U at ", level,i,U(i,:)
             call bl_error("chemterm_begp: beburn failed")
          end if
          
          do n=1,nspec
             U(i,UFS+n-1) = rhot(1)*YT(n,1)
          end do
          
       end if

    end do

    deallocate(UG)

  end subroutine chemterm_begp


  subroutine dUdt_chem_gauss(lo, hi, U, Ulo, Uhi, Ut, Utlo, Uthi, st, stlo, sthi)
    integer, intent(in) :: lo(1), hi(1), Ulo(1), Uhi(1), Utlo(1), Uthi(1), stlo(1), sthi(1)
    double precision, intent(in ) ::  U( Ulo(1): Uhi(1),NVAR)
    double precision, intent(out) :: Ut(Utlo(1):Uthi(1),NVAR)
    double precision, intent(inout) :: st(stlo(1):sthi(1))

    integer :: i, n, g, np, ierr
    double precision :: rho(lo(1):hi(1)), T(lo(1):hi(1))
    double precision :: Y(lo(1):hi(1),nspec), rdYdt(lo(1):hi(1),nspec)
    double precision, allocatable :: UG(:,:,:)

    np = hi(1)-lo(1)+1

    allocate(UG(lo(1):hi(1),NVAR,2))

    if (chem_do_weno) then
       call chem_weno(lo(1), hi(1), U, Ulo(1), Uhi(1), UG)
    else
       do n=1,NVAR
          call cellavg2gausspt_1d(lo(1),hi(1), U(:,n), Ulo(1),Uhi(1), UG(:,n,1), UG(:,n,2), lo(1),hi(1))
       end do
    end if

    do n=1,NVAR
       do i=lo(1),hi(1)
          Ut(i,n) = 0.d0
       end do
    end do
    
    do g=1,2
          
       do i=lo(1),hi(1)
          call get_rhoYT(UG(i,:,g), rho(i), Y(i,:), T(i), ierr)
          if (ierr .ne. 0) then
             st(i) = -1.d0
          end if
       end do

       call compute_rhodYdt(np,rho,T,Y,rdYdt)

       do n=1,nspec
          do i=lo(1),hi(1)
             Ut(i,UFS+n-1) = Ut(i,UFS+n-1) + 0.5d0*rdYdt(i,n)
          end do
       end do
       
    end do

    do i=lo(1),hi(1)
       if (st(i) .ne. 0.d0) then
          call get_rhoYT(U(i,:), rho(i), Y(i,:), T(i), ierr)
          if (ierr .ne. 0) then
             print *, "dUdt_chem_gauss: eos_get_T failed for U at ", level,i,U(i,:)
             call bl_error("dUdt_chem_gauss failed at eos_get_T for U")
          end if
          call compute_rhodYdt(1,rho(i:i),T(i:i),Y(i:i,:),Ut(i:i,UFS:UFS+nspec-1))
       end if
    end do

    deallocate(UG)

  end subroutine dUdt_chem_gauss


  subroutine dUdt_chem_cellcenter(lo, hi, U, Ulo, Uhi, Ut, Utlo, Uthi, st, stlo, sthi)
    integer, intent(in) :: lo(1), hi(1), Ulo(1), Uhi(1), Utlo(1), Uthi(1), stlo(1), sthi(1)
    double precision, intent(in ) ::  U( Ulo(1): Uhi(1),NVAR)
    double precision, intent(out) :: Ut(Utlo(1):Uthi(1),NVAR)
    double precision, intent(inout) :: st(stlo(1):sthi(1))

    integer :: i, n, np, ierr
    double precision :: rho(lo(1)-1:hi(1)+1), Y(lo(1)-1:hi(1)+1,nspec), &
         T(lo(1)-1:hi(1)+1)
    double precision, allocatable :: Ucc(:,:), Utcc(:,:)

    np = hi(1)-lo(1)+3

    allocate( Ucc(lo(1)-1:hi(1)+1,NVAR))
    allocate(Utcc(lo(1)-1:hi(1)+1,nspec))

    do n=1,NVAR
       call cellavg2cc_1d(lo(1)-1,hi(1)+1, U(:,n), Ulo(1),Uhi(1), Ucc(:,n), lo(1)-1,hi(1)+1)
    end do
       
    do i = lo(1)-1, hi(1)+1
       if (st(i) .eq. 0.d0) then
          call get_rhoYT(Ucc(i,:), rho(i), Y(i,:), T(i), ierr)
          if (ierr .ne. 0) then
             st(i) = -1.d0
          end if
       end if
       if (st(i) .ne. 0.d0) then
          call get_rhoYT(U(i,:), rho(i), Y(i,:), T(i), ierr)
          if (ierr .ne. 0) then
             print *, 'dUdt_chem_cellcenter: eos_get_T failed for U at ', level,i,Ucc(i,:)
             call bl_error("dUdt_chem_cellcenter failed at eos_get_T for U")
          end if
       end if
    end do

    call compute_rhodYdt(np,rho,T,Y,Utcc)

    do n=1,UFS-1
       do i=lo(1),hi(1)
          Ut(i,n) = 0.d0
       end do
    end do

    do n=1,nspec
       call cc2cellavg_1d(lo(1),hi(1), Utcc(:,n), lo(1)-1,hi(1)+1, &
            Ut(:,UFS+n-1), Utlo(1),Uthi(1))
    end do

    deallocate(Ucc,Utcc)

  end subroutine dUdt_chem_cellcenter


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
    
    ei = rhoinv*( U(UEDEN) - 0.5d0*rhoinv*U(UMX)**2 )
    
    call floor_species(nspec, Y)

    call eos_get_T(T, ei, Y, ierr=ierr)
  end subroutine get_rhoYT


  subroutine chem_weno(lo, hi, U, Ulo, Uhi, UG)
    use reconstruct_module, only : reconstruct_comp
    integer, intent(in) :: lo, hi, Ulo, Uhi
    double precision, intent(in) :: U(Ulo:Uhi,NVAR)
    double precision, intent(out):: UG(lo:hi,NVAR,2)
    call reconstruct_comp(lo,hi, &
         Ulo,Uhi,           &  ! for U
         0, 0,              &  ! for UL & UR, not present
         lo,hi,             &  ! for UG1 & UG2
         U,     &
         UG1=UG(:,:,1), UG2=UG(:,:,2) )
  end subroutine chem_weno

end module chemterm_module
