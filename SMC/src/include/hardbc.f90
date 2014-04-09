  subroutine store_inflow(qin, U)
    type(physbndry_reg), intent(inout) :: qin
    type(multifab), intent(in) :: U

    integer :: dm,n,ng
    integer ::  lo(U%dim),  hi(U%dim)
    integer :: qlo(U%dim), qhi(U%dim)
    double precision, pointer, dimension(:,:,:,:) :: qp, up
    
    dm = U%dim
    ng = nghost(U)

    do n=1,nfabs(U)
       if (isValid(qin,n)) then
          lo = lwb(get_box(U,n))
          hi = upb(get_box(U,n))

          qlo = lwb(get_box(qin%data,n))
          qhi = upb(get_box(qin%data,n))

          qp => dataptr(qin%data, n)
          up => dataptr(U,n)

          if (dm .eq. 1) then
             call ctoprim_inflow_1d(lo,hi,ng,up,qlo,qhi,qp)
          else if (dm .eq. 2) then
             call ctoprim_inflow_2d(lo,hi,ng,up,qlo,qhi,qp)
          else
             call ctoprim_inflow_3d(lo,hi,ng,up,qlo,qhi,qp)
          end if
       end if
    end do

  end subroutine store_inflow

  subroutine ctoprim_inflow_1d(lo,hi,ng,cons,qlo,qhi,qin)
    integer, intent(in) :: lo(1), hi(1), ng, qlo(1), qhi(1)
    double precision, intent(in) :: cons(-ng+lo(1):hi(1)+ng,ncons)
    double precision, intent(out) :: qin(nqin,qlo(1):qhi(1))

    integer :: i, iwrk, ierr
    double precision :: rhoinv, ei, rwrk

    do i=qlo(1),qhi(1)
       rhoinv = 1.d0/cons(i,irho)
       qin(iuin,i) = cons(i,imx)*rhoinv
       qin(iYin1:,i) = cons(i,iry1:iry1+nspecies-1)*rhoinv
       
       ei = rhoinv*cons(i,iene) - 0.5d0*qin(iuin,i)**2
       qin(iTin,i) = 0.d0
       call get_t_given_ey(ei, qin(iYin1:,i), iwrk, rwrk, qin(iTin,i), ierr)
    end do
  end subroutine ctoprim_inflow_1d

  subroutine ctoprim_inflow_2d(lo,hi,ng,cons,qlo,qhi,qin)
    integer, intent(in) :: lo(2), hi(2), ng, qlo(2), qhi(2)
    double precision, intent(in) :: cons(-ng+lo(1):hi(1)+ng,-ng+lo(2):hi(2)+ng,ncons)
    double precision, intent(out) :: qin(nqin,qlo(1):qhi(1),qlo(2):qhi(2))

    integer :: i, j, iwrk, ierr
    double precision :: rhoinv, ei, rwrk

    !$omp parallel private(i,j,iwrk,rhoinv, ei, rwrk, ierr)
    !$omp do collapse(2)
    do j=qlo(2),qhi(2)
    do i=qlo(1),qhi(1)
       rhoinv = 1.d0/cons(i,j,irho)
       qin(iuin,i,j) = cons(i,j,imx)*rhoinv
       qin(ivin,i,j) = cons(i,j,imy)*rhoinv
       qin(iYin1:,i,j) = cons(i,j,iry1:iry1+nspecies-1)*rhoinv
       
       ei = rhoinv*cons(i,j,iene) - 0.5d0*(qin(iuin,i,j)**2  &
            + qin(ivin,i,j)**2)
       qin(iTin,i,j) = 0.d0
       call get_t_given_ey(ei, qin(iYin1:,i,j), iwrk, rwrk, qin(iTin,i,j), ierr)
    end do
    end do
    !$omp end do
    !$omp end parallel
  end subroutine ctoprim_inflow_2d

  subroutine ctoprim_inflow_3d(lo,hi,ng,cons,qlo,qhi,qin)
    integer, intent(in) :: lo(3), hi(3), ng, qlo(3), qhi(3)
    double precision, intent(in) :: cons(-ng+lo(1):hi(1)+ng,-ng+lo(2):hi(2)+ng,-ng+lo(3):hi(3)+ng,ncons)
    double precision, intent(out) :: qin(nqin,qlo(1):qhi(1),qlo(2):qhi(2),qlo(3):qhi(3))

    integer :: i, j, k, iwrk, ierr
    double precision :: rhoinv, ei, rwrk

    !$omp parallel private(i,j,k,iwrk,rhoinv, ei, rwrk, ierr)
    !$omp do collapse(2)
    do k=qlo(3),qhi(3)
    do j=qlo(2),qhi(2)
    do i=qlo(1),qhi(1)
       rhoinv = 1.d0/cons(i,j,k,irho)
       qin(iuin,i,j,k) = cons(i,j,k,imx)*rhoinv
       qin(ivin,i,j,k) = cons(i,j,k,imy)*rhoinv
       qin(iwin,i,j,k) = cons(i,j,k,imz)*rhoinv
       qin(iYin1:,i,j,k) = cons(i,j,k,iry1:iry1+nspecies-1)*rhoinv
       
       ei = rhoinv*cons(i,j,k,iene) - 0.5d0*(qin(iuin,i,j,k)**2  &
            + qin(ivin,i,j,k)**2 + qin(iwin,i,j,k)**2)
       qin(iTin,i,j,k) = 0.d0
       call get_t_given_ey(ei, qin(iYin1:,i,j,k), iwrk, rwrk, qin(iTin,i,j,k), ierr)
    end do
    end do
    end do
    !$omp end do
    !$omp end parallel
  end subroutine ctoprim_inflow_3d

  subroutine impose_hard_bc(U, t, dx)
    double precision, intent(in) :: t, dx(3)
    type(multifab), intent(inout) :: U

    integer :: n, ng, dm
    integer ::  lo(U%dim),  hi(U%dim)
    integer :: qlo(U%dim), qhi(U%dim)
    double precision, pointer, dimension(:,:,:,:) :: qp, up

    dm = U%dim
    ng = nghost(U)

    do n=1,nfabs(U)

       lo = lwb(get_box(U,n))
       hi = upb(get_box(U,n))
       
       up => dataptr(U,n)

       if (isValid(qin_xlo,n)) then
          qlo = lwb(get_box(qin_xlo%data,n))
          qhi = upb(get_box(qin_xlo%data,n))
          qp  => dataptr(   qin_xlo%data,n)
          if (dm .eq. 1) then
             call update_inlet_xlo_1d(lo,hi,qp(:,lo(1),1,1),t,dx)
             call impose_inflow_1d(lo,hi,ng,up,qlo,qhi,qp,1)
          else if (dm .eq. 2) then
             call update_inlet_xlo_2d(lo,hi,qp(:,lo(1),:,1),t,dx)
             call impose_inflow_2d(lo,hi,ng,up,qlo,qhi,qp,1)
          else
             call update_inlet_xlo_3d(lo,hi,qp(:,lo(1),:,:),t,dx)
             call impose_inflow_3d(lo,hi,ng,up,qlo,qhi,qp,1)
          end if
       end if

       if (isValid(qin_xhi,n)) then
          qlo = lwb(get_box(qin_xhi%data,n))
          qhi = upb(get_box(qin_xhi%data,n))
          qp  => dataptr(   qin_xhi%data,n)
          call bl_error("inflow for xhi not implemented")
          if (dm .eq. 1) then
             call impose_inflow_1d(lo,hi,ng,up,qlo,qhi,qp,1)
          else if (dm .eq. 2) then
             call impose_inflow_2d(lo,hi,ng,up,qlo,qhi,qp,1)
          else
             call impose_inflow_3d(lo,hi,ng,up,qlo,qhi,qp,1)
          end if
       end if

       if (dm.ge.2) then
          if (isValid(qin_ylo,n)) then
             qlo = lwb(get_box(qin_ylo%data,n))
             qhi = upb(get_box(qin_ylo%data,n))
             qp  => dataptr(   qin_ylo%data,n)
             if (dm .eq. 2) then
                call update_inlet_ylo_2d(lo,hi,qp(:,:,lo(2),1),t,dx)
                call impose_inflow_2d(lo,hi,ng,up,qlo,qhi,qp,2)
             else
                call update_inlet_ylo_3d(lo,hi,qp(:,:,lo(2),:),t,dx)
                call impose_inflow_3d(lo,hi,ng,up,qlo,qhi,qp,2)
             end if
          end if
          
          if (isValid(qin_yhi,n)) then
             qlo = lwb(get_box(qin_yhi%data,n))
             qhi = upb(get_box(qin_yhi%data,n))
             qp  => dataptr(   qin_yhi%data,n)
             call bl_error("inflow for yhi not implemented")
             if (dm .eq. 2) then
                call impose_inflow_3d(lo,hi,ng,up,qlo,qhi,qp,2)
             else
                call impose_inflow_3d(lo,hi,ng,up,qlo,qhi,qp,2)
             end if
          end if
       end if

       if (dm .eq. 3) then
          if (isValid(qin_zlo,n)) then
             qlo = lwb(get_box(qin_zlo%data,n))
             qhi = upb(get_box(qin_zlo%data,n))
             qp  => dataptr(   qin_zlo%data,n)
             call bl_error("inflow for zlo not implemented")
             call impose_inflow_3d(lo,hi,ng,up,qlo,qhi,qp,3)
          end if
          
          if (isValid(qin_zhi,n)) then
             qlo = lwb(get_box(qin_zhi%data,n))
             qhi = upb(get_box(qin_zhi%data,n))
             qp  => dataptr(   qin_zhi%data,n)
             call bl_error("inflow for zhi not implemented")
             call impose_inflow_3d(lo,hi,ng,up,qlo,qhi,qp,3)
          end if
       end if

    end do

  end subroutine impose_hard_bc


  subroutine impose_inflow_1d(lo,hi,ng,cons,qlo,qhi,qin,idim)
    integer, intent(in) :: lo(1), hi(1), ng, qlo(1), qhi(1), idim
    double precision, intent(inout) :: cons(-ng+lo(1):hi(1)+ng,ncons)
    double precision, intent(in) :: qin(nqin,qlo(1):qhi(1))

    integer :: i, iwrk
    double precision :: ei, rwrk

    do i=qlo(1),qhi(1)
       cons(i,iry1:iry1+nspecies-1) = cons(i,irho) * qin(iYin1:,i)

       call CKUBMS(qin(iTin,i),qin(iYin1:,i),iwrk,rwrk,ei)

       cons(i,iene) = cons(i,irho)*ei + cons(i,imx)**2 / (2.d0*cons(i,irho))
    end do
  end subroutine impose_inflow_1d

  subroutine impose_inflow_2d(lo,hi,ng,cons,qlo,qhi,qin,idim)
    integer, intent(in) :: lo(2), hi(2), ng, qlo(2), qhi(2), idim
    double precision, intent(inout) :: cons(-ng+lo(1):hi(1)+ng,-ng+lo(2):hi(2)+ng,ncons)
    double precision, intent(in) :: qin(nqin,qlo(1):qhi(1),qlo(2):qhi(2))

    integer :: i, j, iwrk
    double precision :: ei, rwrk

    !$omp parallel private(i,j,iwrk, ei, rwrk)
    !$omp do collapse(2)
    do j=qlo(2),qhi(2)
    do i=qlo(1),qhi(1)
       if (idim.eq.1) then
          cons(i,j,imy) = cons(i,j,irho) * qin(ivin,i,j)
       else
          cons(i,j,imx) = cons(i,j,irho) * qin(iuin,i,j)
       end if

       cons(i,j,iry1:iry1+nspecies-1) = cons(i,j,irho) * qin(iYin1:,i,j)

       call CKUBMS(qin(iTin,i,j),qin(iYin1:,i,j),iwrk,rwrk,ei)

       cons(i,j,iene) = cons(i,j,irho)*ei + (cons(i,j,imx)**2 &
            + cons(i,j,imy)**2) / (2.d0*cons(i,j,irho))
    end do
    end do
    !$omp end do
    !$omp end parallel
  end subroutine impose_inflow_2d

  subroutine impose_inflow_3d(lo,hi,ng,cons,qlo,qhi,qin,idim)
    integer, intent(in) :: lo(3), hi(3), ng, qlo(3), qhi(3), idim
    double precision, intent(inout) :: cons(-ng+lo(1):hi(1)+ng,-ng+lo(2):hi(2)+ng,-ng+lo(3):hi(3)+ng,ncons)
    double precision, intent(in) :: qin(nqin,qlo(1):qhi(1),qlo(2):qhi(2),qlo(3):qhi(3))

    integer :: i, j, k, iwrk
    double precision :: ei, rwrk

    !$omp parallel private(i,j,k,iwrk, ei, rwrk)
    !$omp do collapse(2)
    do k=qlo(3),qhi(3)
    do j=qlo(2),qhi(2)
    do i=qlo(1),qhi(1)
       if (idim.eq.1) then
          cons(i,j,k,imy) = cons(i,j,k,irho) * qin(ivin,i,j,k)
          cons(i,j,k,imz) = cons(i,j,k,irho) * qin(iwin,i,j,k)
       else if (idim.eq.2) then
          cons(i,j,k,imx) = cons(i,j,k,irho) * qin(iuin,i,j,k)
          cons(i,j,k,imz) = cons(i,j,k,irho) * qin(iwin,i,j,k)
       else
          cons(i,j,k,imx) = cons(i,j,k,irho) * qin(iuin,i,j,k)
          cons(i,j,k,imy) = cons(i,j,k,irho) * qin(ivin,i,j,k)
       end if

       cons(i,j,k,iry1:iry1+nspecies-1) = cons(i,j,k,irho) * qin(iYin1:,i,j,k)

       call CKUBMS(qin(iTin,i,j,k),qin(iYin1:,i,j,k),iwrk,rwrk,ei)

       cons(i,j,k,iene) = cons(i,j,k,irho)*ei + (cons(i,j,k,imx)**2 &
            + cons(i,j,k,imy)**2 + cons(i,j,k,imz)**2) / (2.d0*cons(i,j,k,irho))
    end do
    end do
    end do
    !$omp end do
    !$omp end parallel
  end subroutine impose_inflow_3d
