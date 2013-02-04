  subroutine store_inflow(qin, U)
    type(physbndry_reg), intent(inout) :: qin
    type(multifab), intent(in) :: U

    integer :: n,ng
    integer ::  lo(U%dim),  hi(U%dim)
    integer :: qlo(U%dim), qhi(U%dim)
    double precision, pointer, dimension(:,:,:,:) :: qp, up
    
    ng = nghost(U)

    do n=1,nfabs(U)
       if (isValid(qin,n)) then
          lo = lwb(get_box(U,n))
          hi = upb(get_box(U,n))

          qlo = lwb(get_box(qin%data,n))
          qhi = upb(get_box(qin%data,n))

          qp => dataptr(qin%data, n)
          up => dataptr(U,n)

          call ctoprim_inflow(lo,hi,ng,up,qlo,qhi,qp)
       end if
    end do

  end subroutine store_inflow

  subroutine ctoprim_inflow(lo,hi,ng,cons,qlo,qhi,qin)
    integer, intent(in) :: lo(3), hi(3), ng, qlo(3), qhi(3)
    double precision, intent(in) :: cons(-ng+lo(1):hi(1)+ng,-ng+lo(2):hi(2)+ng,-ng+lo(3):hi(3)+ng,ncons)
    double precision, intent(out) :: qin(nqin,qlo(1):qhi(1),qlo(2):qhi(2),qlo(3):qhi(3))

    integer :: i, j, k, iwrk
    double precision :: rhoinv, ei, rwrk

    !$omp parallel private(i,j,k,iwrk,rhoinv, ei, rwrk)
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
       call feeytt(ei, qin(iYin1:,i,j,k), iwrk, rwrk, qin(iTin,i,j,k))
    end do
    end do
    end do
    !$omp end do
    !$omp end parallel
  end subroutine ctoprim_inflow

  subroutine impose_hard_bc(U)
    type(multifab), intent(inout) :: U

    integer :: n, ng
    integer ::  lo(U%dim),  hi(U%dim)
    integer :: qlo(U%dim), qhi(U%dim)
    double precision, pointer, dimension(:,:,:,:) :: qp, up

    ng = nghost(U)

    do n=1,nfabs(U)

       lo = lwb(get_box(U,n))
       hi = upb(get_box(U,n))
       
       up => dataptr(U,n)

       if (isValid(qin_xlo,n)) then
          qlo = lwb(get_box(qin_xlo%data,n))
          qhi = upb(get_box(qin_xlo%data,n))
          qp  => dataptr(   qin_xlo%data,n)
          call impose_inflow(lo,hi,ng,up,qlo,qhi,qp,1)
       end if

       if (isValid(qin_xhi,n)) then
          qlo = lwb(get_box(qin_xhi%data,n))
          qhi = upb(get_box(qin_xhi%data,n))
          qp  => dataptr(   qin_xhi%data,n)
          call impose_inflow(lo,hi,ng,up,qlo,qhi,qp,1)
       end if

       if (isValid(qin_ylo,n)) then
          qlo = lwb(get_box(qin_ylo%data,n))
          qhi = upb(get_box(qin_ylo%data,n))
          qp  => dataptr(   qin_ylo%data,n)
          call impose_inflow(lo,hi,ng,up,qlo,qhi,qp,2)
       end if

       if (isValid(qin_yhi,n)) then
          qlo = lwb(get_box(qin_yhi%data,n))
          qhi = upb(get_box(qin_yhi%data,n))
          qp  => dataptr(   qin_yhi%data,n)
          call impose_inflow(lo,hi,ng,up,qlo,qhi,qp,2)
       end if

       if (isValid(qin_zlo,n)) then
          qlo = lwb(get_box(qin_zlo%data,n))
          qhi = upb(get_box(qin_zlo%data,n))
          qp  => dataptr(   qin_zlo%data,n)
          call impose_inflow(lo,hi,ng,up,qlo,qhi,qp,3)
       end if

       if (isValid(qin_zhi,n)) then
          qlo = lwb(get_box(qin_zhi%data,n))
          qhi = upb(get_box(qin_zhi%data,n))
          qp  => dataptr(   qin_zhi%data,n)
          call impose_inflow(lo,hi,ng,up,qlo,qhi,qp,3)
       end if

    end do

  end subroutine impose_hard_bc


  subroutine impose_inflow(lo,hi,ng,cons,qlo,qhi,qin,idim)
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
  end subroutine impose_inflow
