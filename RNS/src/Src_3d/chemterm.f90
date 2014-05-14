module chemterm_module

  use meth_params_module
  use burner_module, only : burn, compute_rhodYdt, splitburn, beburn
  use eos_module, only : eos_get_T
  use renorm_module, only : renorm

  implicit none

  private

  public :: chemterm, dUdt_chem

contains

  subroutine chemterm(lo, hi, U, Ulo, Uhi, dt, Up)
    integer, intent(in) :: lo(3), hi(3), Ulo(3), Uhi(3)
    double precision, intent(inout) :: U(Ulo(1):Uhi(1),Ulo(2):Uhi(2),Ulo(3):Uhi(3),NVAR)
    double precision, intent(in) :: dt
    double precision, intent(in), optional :: Up(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),NVAR)

    select case (chem_solver)
       case (cc_burning)
          call chemterm_cellcenter(lo, hi, U, Ulo, Uhi, dt)
       case (Gauss_burning)
          call chemterm_gauss(lo, hi, U, Ulo, Uhi, dt)
       case (split_burning)
          call chemterm_split(lo, hi, U, Ulo, Uhi, dt)
       case (BE_burning)
          call chemterm_be(lo, hi, U, Ulo, Uhi, dt, Up)
       case default
          call bl_error("unknown chem_solver")
       end select

  end subroutine chemterm


  subroutine chemterm_cellcenter(lo, hi, U, Ulo, Uhi, dt)
    use convert_module, only : cellavg2cc_3d, cc2cellavg_3d
    integer, intent(in) :: lo(3), hi(3), Ulo(3), Uhi(3)
    double precision, intent(inout) :: U(Ulo(1):Uhi(1),Ulo(2):Uhi(2),Ulo(3):Uhi(3),NVAR)
    double precision, intent(in) :: dt

    integer :: i, j, k, n
    logical :: force_new_J
    double precision :: rhot(1), rhoinv, ei
    double precision :: Yt(nspec+1)
    double precision, allocatable :: Ucc(:,:,:,:)

    allocate(Ucc(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1,NVAR))

    !$omp parallel private(i,j,k,n,rhot,rhoinv,ei,Yt,force_new_J)

    !$omp do
    do n=1,NVAR
       call cellavg2cc_3d(lo-1,hi+1, U(:,:,:,n), Ulo,Uhi, Ucc(:,:,:,n), lo-1,hi+1)
    end do
    !$omp end do

    force_new_J = .true.  ! always recompute Jacobian when a new FAB starts

    !$omp do collapse(2)
    do k=lo(3)-1,hi(3)+1
       do j=lo(2)-1,hi(2)+1
          do i=lo(1)-1,hi(1)+1

             rhot = 0.d0
             do n=1,NSPEC
                Yt(n) = Ucc(i,j,k,UFS+n-1)
                rhot(1) = rhot(1) + Yt(n)
             end do
             rhoinv = 1.d0/rhot(1)

             Yt(1:nspec) = Yt(1:nspec) * rhoinv
             Yt(nspec+1) = Ucc(i,j,k,UTEMP)

             ei = rhoinv*( Ucc(i,j,k,UEDEN) - 0.5d0*rhoinv*(Ucc(i,j,k,UMX)**2 &
                  + Ucc(i,j,k,UMY)**2 + Ucc(i,j,k,UMZ)**2) )

             call eos_get_T(Yt(nspec+1), ei, Yt(1:nspec))

             call burn(1, rhot, Yt, dt, force_new_J)

             force_new_J = new_J_cell

             do n=1,nspec
                Ucc(i,j,k,UFS+n-1) = rhot(1)*Yt(n)
             end do

          end do
       end do
    end do
    !$omp end do

    !$omp do
    do n=UFS,UFS+nspec-1
       call cc2cellavg_3d(lo,hi, Ucc(:,:,:,n), lo-1,hi+1, U(:,:,:,n), Ulo,Uhi)
    end do
    !$omp end do

    !$omp end parallel

    deallocate(Ucc)

  end subroutine chemterm_cellcenter


  subroutine chemterm_gauss(lo, hi, U, Ulo, Uhi, dt)
    use weno_module, only : cellavg2gausspt_3d
    integer, intent(in) :: lo(3), hi(3), Ulo(3), Uhi(3)
    double precision, intent(inout) :: U(Ulo(1):Uhi(1),Ulo(2):Uhi(2),Ulo(3):Uhi(3),NVAR)
    double precision, intent(in) :: dt

    integer :: i, j, k, n, g
    logical :: force_new_J
    double precision :: rhot(8), rhoinv, ei
    double precision :: Yt(nspec+1,8)
    double precision, allocatable :: UG(:,:,:,:,:)

    allocate(UG(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),8,NVAR))

    !$omp parallel private(i,j,k,n,g,rhot,rhoinv,ei,Yt,force_new_J)

    !$omp do
    do n=1,NVAR
       call cellavg2gausspt_3d(lo,hi, U(:,:,:,n), Ulo,Uhi, UG(:,:,:,:,n), lo,hi)
    end do
    !$omp end do

    force_new_J = .true.  ! always recompute Jacobian when a new FAB starts

    !$omp do collapse(2)
    do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)

             do g=1,8

                rhot(g) = 0.d0
                do n=1,NSPEC
                   Yt(n,g) = UG(i,j,k,g,UFS+n-1)
                   rhot(g) = rhot(g) + Yt(n,g)
                end do
                rhoinv = 1.d0/rhot(g)

                Yt(1:nspec,g) = Yt(1:nspec,g) * rhoinv
                Yt(nspec+1,g) = UG(i,j,k,g,UTEMP)

                ei = rhoinv*( UG(i,j,k,g,UEDEN) - 0.5d0*rhoinv*(UG(i,j,k,g,UMX)**2 &
                     + UG(i,j,k,g,UMY)**2 + UG(i,j,k,g,UMZ)**2) )

                call eos_get_T(Yt(nspec+1,g), ei, Yt(1:nspec,g))

             end do

             call burn(8, rhot, Yt, dt, force_new_J)

             force_new_J = new_J_cell

             U(i,j,k,UFS:UFS+nspec-1) = 0.d0 
             do g=1,8
                do n=1,nspec
                   U(i,j,k,UFS+n-1) = U(i,j,k,UFS+n-1) + 0.125d0*rhot(g)*Yt(n,g)
                end do
             end do

          end do
       end do
    end do
    !$omp end do
    !$omp end parallel

    deallocate(UG)

  end subroutine chemterm_gauss


  subroutine chemterm_split(lo, hi, U, Ulo, Uhi, dt)
    use weno_module, only : cellavg2gausspt_3d
    integer, intent(in) :: lo(3), hi(3), Ulo(3), Uhi(3)
    double precision, intent(inout) :: U(Ulo(1):Uhi(1),Ulo(2):Uhi(2),Ulo(3):Uhi(3),NVAR)
    double precision, intent(in) :: dt

    integer :: i, j, k, n, g
    logical :: force_new_J
    double precision :: rhot(8), rhoinv, ei, rho0(1)
    double precision :: Yt(nspec+1,8), Y0(nspec+1)
    double precision, allocatable :: UG(:,:,:,:,:)

    allocate(UG(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),8,NVAR))

    !$omp parallel private(i,j,k,n,g,rhot,rhoinv,ei,Yt,force_new_J,rho0,Y0)

    !$omp do
    do n=1,NVAR
       call cellavg2gausspt_3d(lo,hi, U(:,:,:,n), Ulo,Uhi, UG(:,:,:,:,n), lo,hi)
    end do
    !$omp end do

    force_new_J = .true.  ! always recompute Jacobina when a new FAB starts

    !$omp do collapse(2)
    do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)

             Y0 = 0.d0
             rho0(1) = 0.d0

             do g=1,8

                rhot(g) = 0.d0
                do n=1,NSPEC
                   Yt(n,g) = UG(i,j,k,g,UFS+n-1)
                   rhot(g) = rhot(g) + Yt(n,g)
                end do
                rhoinv = 1.d0/rhot(g)

                Yt(1:nspec,g) = Yt(1:nspec,g) * rhoinv
                Yt(nspec+1,g) = UG(i,j,k,g,UTEMP)

                ei = rhoinv*( UG(i,j,k,g,UEDEN) - 0.5d0*rhoinv*(UG(i,j,k,g,UMX)**2 &
                     + UG(i,j,k,g,UMY)**2 + UG(i,j,k,g,UMZ)**2) )

                call eos_get_T(Yt(nspec+1,g), ei, Yt(1:nspec,g))

                Y0 = Y0 + 0.125d0*Yt(:,g)
                rho0(1) = rho0(1) + 0.125d0*rhot(g)

             end do

             call burn(1, rho0, Y0, dt, force_new_J)

             force_new_J = new_J_cell

             call splitburn(8, rho0(1), Y0, rhot, Yt, dt)
             ! Now Yt is \Delta Y and T

             U(i,j,k,UFS:UFS+nspec-1) = 0.d0 
             do g=1,8
                do n=1,nspec
                   U(i,j,k,UFS+n-1) = U(i,j,k,UFS+n-1) + &
                        0.125d0*(UG(i,j,k,g,UFS+n-1) + rhot(g)*Yt(n,g))
                end do
             end do

          end do
       end do
    end do
    !$omp end do
    !$omp end parallel

    deallocate(UG)

  end subroutine chemterm_split


  subroutine chemterm_be(lo, hi, U, Ulo, Uhi, dt, Up)
    use weno_module, only : cellavg2gausspt_3d
    integer, intent(in) :: lo(3), hi(3), Ulo(3), Uhi(3)
    double precision, intent(inout) :: U(Ulo(1):Uhi(1),Ulo(2):Uhi(2),Ulo(3):Uhi(3),NVAR)
    double precision, intent(in) :: dt
    double precision, intent(in), optional :: Up(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),NVAR)

    integer :: i, j, k, n, g, ierr
    logical :: force_new_J
    double precision :: rhot(8), rhoinv, rho0(1)
    double precision :: Yt(nspec+1,8), Y0(nspec+1), rhoY(nspec)
    double precision, allocatable :: UG(:,:,:,:,:)

    allocate(UG(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),8,NVAR))

    if (chem_do_weno) then
       call chem_weno(lo, hi, U, Ulo, Uhi, UG)
    else
       !$omp parallel do private(n)
       do n=1,NVAR
          call cellavg2gausspt_3d(lo,hi, U(:,:,:,n), Ulo,Uhi, UG(:,:,:,:,n), lo,hi)
       end do
       !$omp end parallel do
    end if

    !$omp parallel private(i,j,k,n,g,ierr,rhot,rhoinv,Yt,force_new_J,rho0,Y0,rhoY)

    force_new_J = .true.  ! always recompute Jacobina when a new FAB starts

    !$omp do schedule(dynamic) collapse(2)
     do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)

             do g=1,8
                call get_rhoYT(UG(i,j,k,g,:), rhot(g), YT(:,g), ierr)
                if (ierr .ne. 0) then
                   print *, 'chemterm_be: eos_get_T failed at ', i,j,k,g,UG(i,j,k,g,:)
                   call flush(6)
                   call bl_error("chemterm_be failed at eos_get_T")
                end if
             end do

             if (present(Up)) then

                rho0(1) = Up(i,j,k,URHO)
                rhoinv = 1.d0/rho0(1)
                Y0(1:nspec) = Up(i,j,k,UFS:UFS+nspec-1)*rhoinv
                Y0(nspec+1) = Up(i,j,k,UTEMP)

             else

                call get_rhoYT(U(i,j,k,:), rho0(1), Y0, ierr)
                if (ierr .ne. 0) then
                   print *, 'chemterm_be: eos_get_T failed for U at ', i,j,k,U(i,j,k,:)
                   call flush(6)
                   call bl_error("chemterm_be failed at eos_get_T for U")
                end if
                
                call burn(1, rho0, Y0, dt, force_new_J, ierr)
                force_new_J = new_J_cell
                if (ierr .ne. 0) then
                   print *, 'chemterm_be: burn failed at ', i,j,k,U(i,j,k,:)
                   print *, '   rho0, Y0 =', rho0(1), Y0
                   call flush(6)
                   call bl_error("chemterm_be failed at burn")
                end if

             end if

             call renorm(nspec, Y0(1:nspec), ierr)
             if (ierr .ne. 0) then
                call bl_error("chemterm_be failed at renormalizing Y0")
             end if

             rhoY = 0.d0
             do g=1,8
                call beburn(rho0(1), Y0, rhot(g), Yt(:,g), dt, g, ierr)
                if (ierr .ne. 0) then ! beburn failed
                   print *, 'chemterm_be: beburn failed at ',i,j,k,g,UG(i,j,k,g,:)
                   call flush(6)
                   call bl_error("chemterm_be: beburn failed at g")
                end if
                do n=1,nspec
                   rhoY(n) = rhoY(n) + 0.125d0*rhot(g)*Yt(n,g)
                end do
             end do

             U(i,j,k,UFS:UFS+nspec-1) = rhoY

          end do
       end do
    end do
    !$omp end do
    !$omp end parallel

    deallocate(UG)

  end subroutine chemterm_be


  subroutine dUdt_chem(lo, hi, U, Ulo, Uhi, Ut, Utlo, Uthi)
    use weno_module, only : cellavg2gausspt_3d
    integer, intent(in) :: lo(3), hi(3), Ulo(3), Uhi(3), Utlo(3), Uthi(3)
    double precision, intent(in ) ::  U( Ulo(1): Uhi(1), Ulo(2): Uhi(2), Ulo(3): Uhi(3),NVAR)
    double precision, intent(out) :: Ut(Utlo(1):Uthi(1),Utlo(2):Uthi(2),Utlo(3):Uthi(3),NVAR)

    integer :: i, j, k, n, g, np, ierr
    double precision :: rhoinv, ei
    double precision :: rho(lo(1):hi(1)), T(lo(1):hi(1))
    double precision :: Ytmp(nspec)
    double precision :: Y(lo(1):hi(1),nspec), rdYdt(lo(1):hi(1),nspec)
    double precision, allocatable :: UG(:,:,:,:,:)

    np = hi(1)-lo(1)+1

    allocate(UG(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),8,NVAR))

    if (chem_do_weno) then
       call chem_weno(lo, hi, U, Ulo, Uhi, UG)
    else
       !$omp parallel do private(n)
       do n=1,NVAR
          call cellavg2gausspt_3d(lo,hi, U(:,:,:,n), Ulo,Uhi, UG(:,:,:,:,n), lo,hi)
       end do
       !$omp end parallel do
    end if

    !$omp parallel private(i,j,k,n,g,ierr,rhoinv,ei,rho,T,Ytmp,Y,rdYdt)

    !$omp do
    do n=1,NVAR
       do k=lo(3),hi(3)
       do j=lo(2),hi(2)
       do i=lo(1),hi(1)
          Ut(i,j,k,n) = 0.d0
       end do
       end do
       end do
    end do
    !$omp end do

    !$omp do collapse(2)
    do k=lo(3),hi(3)
       do j=lo(2),hi(2)

          do g=1,8

             do i=lo(1),hi(1)
                rho(i) = 0.d0
                do n=1,nspec
                   Y(i,n) = UG(i,j,k,g,UFS+n-1)
                   rho(i) = rho(i) + Y(i,n)
                end do
                rhoinv = 1.d0/rho(i)
                
                do n=1,nspec
                   Y(i,n) = Y(i,n) * rhoinv
                   Ytmp(n) = Y(i,n)
                end do
                
                ei = rhoinv*( UG(i,j,k,g,UEDEN) - 0.5d0*rhoinv*(UG(i,j,k,g,UMX)**2 &
                     + UG(i,j,k,g,UMY)**2 + UG(i,j,k,g,UMZ)**2) )
                
                T(i) = UG(i,j,k,g,UTEMP)
                call eos_get_T(T(i), ei, Ytmp, ierr=ierr)

                if (ierr .ne. 0) then
                   print *, 'dUdt_chem failed at ', i,j,k,g,UG(i,j,k,g,:)
                   call flush(6)
                   call bl_error("dUdt_chem failed")
                end if
             end do
             
             call compute_rhodYdt(np,rho,T,Y,rdYdt)
             
             do n=1,nspec
                do i=lo(1),hi(1)
                   Ut(i,j,k,UFS+n-1) = Ut(i,j,k,UFS+n-1) + 0.125d0*rdYdt(i,n)
                end do
             end do
             
          end do

       end do
    end do
    !$omp end do 

    !$omp end parallel

    deallocate(UG)

  end subroutine dUdt_chem

  subroutine get_rhoYT(U,rho,YT,ierr)
    double precision, intent(in) :: U(NVAR)
    double precision, intent(out) :: rho, YT(nspec+1)
    integer, intent(out) :: ierr

    integer :: n
    double precision :: rhoinv, ei
    
    rho = 0.d0
    do n=1,NSPEC
       Yt(n) = U(UFS+n-1)
       rho = rho + Yt(n)
    end do
    rhoinv = 1.d0/rho
    
    Yt(1:nspec) = Yt(1:nspec) * rhoinv
    Yt(nspec+1) = U(UTEMP)
    
    ei = rhoinv*( U(UEDEN) - 0.5d0*rhoinv*(U(UMX)**2 &
         + U(UMY)**2 + U(UMZ)**2) )
    
    call eos_get_T(Yt(nspec+1), ei, Yt(1:nspec), ierr=ierr)
  end subroutine get_rhoYT


  subroutine chem_weno(lo, hi, U, Ulo, Uhi, UG)
    use reconstruct_module, only : reconstruct_comp
    integer, intent(in) :: lo(3), hi(3), Ulo(3), Uhi(3)
    double precision, intent(in) :: U(Ulo(1):Uhi(1),Ulo(2):Uhi(2),Ulo(3):Uhi(3),NVAR)
    double precision, intent(out):: UG(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),8,NVAR)

    integer :: i,j,k
    double precision, allocatable :: UG1z(:,:,:,:), UG2z(:,:,:,:)
    double precision, allocatable :: UG1y(:,:,:,:), UG2y(:,:,:,:)

    allocate(UG1z(lo(1)-2:hi(1)+2,lo(2)-2:hi(2)+2,lo(3):hi(3),NVAR))
    allocate(UG2z(lo(1)-2:hi(1)+2,lo(2)-2:hi(2)+2,lo(3):hi(3),NVAR))
    allocate(UG1y(lo(1)-2:hi(1)+2,lo(2)  :hi(2)  ,lo(3):hi(3),NVAR))
    allocate(UG2y(lo(1)-2:hi(1)+2,lo(2)  :hi(2)  ,lo(3):hi(3),NVAR))

    !$omp parallel private(i,j,k)

    !$omp do collapse(2)
    do j   =lo(2)-2,hi(2)+2
       do i=lo(1)-2,hi(1)+2
          call reconstruct_comp(lo(3),hi(3), &
               Ulo(3),Uhi(3),   &  ! for U
               0, 0,            &  ! for UL & UR, not present
               lo (3), hi(3),   &  ! for UG1 & UG2
               U(i,j,:,:), &
               UG1=UG1z(i,j,:,:), UG2=UG2z(i,j,:,:) )
       end do
    end do
    !$omp end do

    !$omp do collapse(2)
    do k   =lo(3)  ,hi(3)
       do i=lo(1)-2,hi(1)+2
          call reconstruct_comp(lo(2),hi(2), &
               lo(2)-2,hi(2)+2,   &  ! for U
               0, 0,              &  ! for UL & UR, not present
               lo(2)  ,hi(2),     &  ! for UG1 & UG2
               UG1Z(i,:,k,:), &
               UG1=UG1y(i,:,k,:), UG2=UG2y(i,:,k,:) )
       end do
    end do
    !$omp end do
    
    !$omp do collapse(2)
    do k   =lo(3), hi(3)
       do j=lo(2), hi(2)
          call reconstruct_comp(lo(1),hi(1), &
               lo(1)-2,hi(1)+2,   &  ! for U
               0, 0,              &  ! for UL & UR, not present
               lo(1)  ,hi(1),     &  ! for UG1 & UG2
               UG1y(:,j,k,:), &
               UG1=UG(:,j,k,1,:), UG2=UG(:,j,k,2,:) )
          call reconstruct_comp(lo(1),hi(1), &
               lo(1)-2,hi(1)+2,   &  ! for U
               0, 0,              &  ! for UL & UR, not present
               lo(1)  ,hi(1),     &  ! for UG1 & UG2
               UG2y(:,j,k,:), &
               UG1=UG(:,j,k,3,:), UG2=UG(:,j,k,4,:) )
       end do
    end do
    !$omp end do 
    
    !$omp do collapse(2)
    do k   =lo(3)  ,hi(3)
       do i=lo(1)-2,hi(1)+2
          call reconstruct_comp(lo(2),hi(2), &
               lo(2)-2,hi(2)+2,   &  ! for U
               0, 0,              &  ! for UL & UR, not present
               lo(2)  ,hi(2),     &  ! for UG1 & UG2
               UG2Z(i,:,k,:), &
               UG1=UG1y(i,:,k,:), UG2=UG2y(i,:,k,:) )
       end do
    end do
    !$omp end do
    
    !$omp do collapse(2)
    do k   =lo(3), hi(3)
       do j=lo(2), hi(2)
          call reconstruct_comp(lo(1),hi(1), &
               lo(1)-2,hi(1)+2,   &  ! for U
               0, 0,              &  ! for UL & UR, not present
               lo(1)  ,hi(1),     &  ! for UG1 & UG2
               UG1y(:,j,k,:), &
               UG1=UG(:,j,k,5,:), UG2=UG(:,j,k,6,:) )
          call reconstruct_comp(lo(1),hi(1), &
               lo(1)-2,hi(1)+2,   &  ! for U
               0, 0,              &  ! for UL & UR, not present
               lo(1)  ,hi(1),     &  ! for UG1 & UG2
               UG2y(:,j,k,:), &
               UG1=UG(:,j,k,7,:), UG2=UG(:,j,k,8,:) )
       end do
    end do
    !$omp end do 

    !$omp end parallel

    deallocate(UG1z,UG2z,UG1y,UG2y)

  end subroutine chem_weno

end module chemterm_module
