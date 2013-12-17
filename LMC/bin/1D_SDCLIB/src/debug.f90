module debug
  implicit none
  integer, parameter :: gpun = 66, datun = 67
  logical, save      :: connected = .false.
contains

  subroutine plotX(q, window, options, wait)
    real(8),      intent(in   ) :: q(:,:)
    integer,      intent(in   ) :: window
    character(*), intent(in   ) :: options
    logical,      intent(in   ) :: wait

    character(128) :: fname

    integer :: i

#ifdef GNUPLOT
    if (.not. connected) then
       open(unit=gpun, file="/tmp/gp")
       connected = .true.
    end if

    write(fname,"(a,i1)") "/tmp/debug.dat", window
    open(unit=datun, file=fname, action="write")
    do i = 1, size(q, dim=2)
       write(datun,*) i, q(:,i)
    end do
    close(datun)

    write(gpun,*) "set term wxt ", window
    write(gpun,*) "plot '" // trim(fname) // "' " // options

    if (wait) then
       write (*,*) '==> paused'
       read  (*,*)
    end if
#endif
  end subroutine plotX

  subroutine plot1(q1, window, options, wait)
    real(8),      intent(in   ) :: q1(:)
    integer,      intent(in   ) :: window
    character(*), intent(in   ) :: options
    logical,      intent(in   ) :: wait
    real(8) :: q(1,size(q1))
    q(1,:) = q1
    call plotX(q, window, options, wait)
  end subroutine plot1

  subroutine plot2(q1, q2, window, options, wait)
    real(8),      intent(in   ) :: q1(:), q2(:)
    integer,      intent(in   ) :: window
    character(*), intent(in   ) :: options
    logical,      intent(in   ) :: wait
    real(8) :: q(2,size(q1))
    q(1,:) = q1
    q(2,:) = q2
    call plotX(q, window, options, wait)
  end subroutine plot2

end module debug
