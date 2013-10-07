module debug_module
  use iso_c_binding, only: c_ptr, c_int, c_null_ptr, c_double, c_associated
  implicit none

  logical, save :: debug = .false.
  type(c_ptr), save :: default_zmq_context = c_null_ptr

#ifdef ZMQ
  interface
     type(c_ptr) function dzmq_connect() bind(c)
       import :: c_ptr
     end function dzmq_connect

     subroutine dzmq_send(ptr, q, nx) bind(c)
       import :: c_ptr, c_int, c_double
       type(c_ptr), intent(in), value :: ptr
       integer(c_int), intent(in), value :: nx
       real(c_double), intent(in) :: q(nx)
     end subroutine dzmq_send

     subroutine dzmq_send_size(ptr, nx, ny) bind(c)
       import :: c_ptr, c_int
       type(c_ptr), intent(in), value :: ptr
       integer(c_int), intent(in), value :: nx, ny
     end subroutine dzmq_send_size

     subroutine dzmq_close(ptr) bind(c)
       import :: c_ptr
       type(c_ptr), intent(in), value :: ptr
     end subroutine dzmq_close

  end interface
#endif

contains

  subroutine dsend(q, wait)
    real(8), intent(in) :: q(:)
    logical, intent(in) :: wait

    if (.not. c_associated(default_zmq_context)) then
       default_zmq_context = dzmq_connect()
    end if

#ifdef ZMQ    
    call dzmq_send(default_zmq_context, q/maxval(abs(q)), size(q))

    if (wait) then
       write (*,*) '==> paused'
       read  (*,*)
    end if
#endif

  end subroutine dsend

end module debug_module
