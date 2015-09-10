program chemh_converge6

  implicit none

  character in0*(32)
  character in1*(32)
  character in2*(32)

  integer i,j,k,nsteps_0,nsteps_1,nsteps_2
  integer nx_0,nx_1,nx_2
  real*8 time_0,time_1,time_2,sum
  
  real*8 data0  (0:4096,16)
  real*8 data1  (0:4096,16)
  real*8 data2  (0:4096,16)
  real*8 data2_c(0:4096,16)

  real*8 L0_0(16)
  real*8 L1_0(16)
  real*8 L2_0(16)
  real*8 L0_1(16)
  real*8 L1_1(16)
  real*8 L2_1(16)
  real*8 L0_2(16)
  real*8 L1_2(16)
  real*8 L2_2(16)

  read(*,*) in0
  read(*,*) in1
  read(*,*) in2

  open(99,file=in0,form='formatted')
  read(99,*) nsteps_0
  read(99,*) nx_0
  read(99,*) time_0

  do i=0,nx_0-1
     read(99,*) data0(i,1:16)
  end do

  open(10,file=in1,form='formatted')
  read(10,*) nsteps_1
  read(10,*) nx_1
  read(10,*) time_1

  do i=0,nx_1-1
     read(10,*) data1(i,1:16)
  end do

  open(20,file=in2,form='formatted')
  read(20,*) nsteps_2
  read(20,*) nx_2
  read(20,*) time_2

  do i=0,nx_2-1
     read(20,*) data2(i,1:16)
  end do

  if (nx_0 .ne. nx_1) then
     print*,'ERROR: nx needs to match in the first two inputs'
     stop
  end if

  !  coarsen fine data
  do i=0,nx_0-1
     do j=1,16
        sum = 0.d0
        do k=0,1
           sum = sum + data2(2*i+k,j)
        end do
        data2_c(i,j) = sum / 2.d0
     end do
  end do
  
  j=14
  do i=0,nx_0-1
     data2_c(i,j) = data2(i*2,j)
  end do
  

  !!!!!!!!!!!!!!!!!!!
  ! error for input1
  !!!!!!!!!!!!!!!!!!!

  L0_0 = 0.d0
  L1_0 = 0.d0
  L2_0 = 0.d0

  do i=0,nx_0-1
     do j=1,16
        L0_0(j) = max(L0_0(j), abs(data0(i,j)-data2_c(i,j)))
        L1_0(j) = L1_0(j) + abs(data0(i,j)-data2_c(i,j))
        L2_0(j) = L2_0(j) + (data0(i,j)-data2_c(i,j))**2
     end do
  end do
  L1_0 = L1_0 / dble(nx_0)
  L2_0 = sqrt(L2_0/nx_0)

  !  coarsen fine data
  do i=0,nx_1-1
     do j=1,16
        sum = 0.d0
        do k=0,1
           sum = sum + data2(2*i+k,j)
        end do
        data2_c(i,j) = sum / 2.d0
     end do
  end do

  j=14  
  do i=0,nx_2-1
     data2_c(i,j) = data2(i*2,j)
  end do
  

  !!!!!!!!!!!!!!!!!!!
  ! error for input2
  !!!!!!!!!!!!!!!!!!!

  L0_1 = 0.d0
  L1_1 = 0.d0
  L2_1 = 0.d0

  do i=0,nx_1-1
     do j=1,16
        L0_1(j) = max(L0_1(j), abs(data1(i,j)-data2_c(i,j)))
        L1_1(j) = L1_1(j) + abs(data1(i,j)-data2_c(i,j))
        L2_1(j) = L2_1(j) + (data1(i,j)-data2_c(i,j))**2
     end do
  end do
  L1_1 = L1_1 / dble(nx_1)
  L2_1 = sqrt(L2_1/nx_1)

1000 format(a,es12.4,es12.4,es12.4)
1001 format(a,es9.2,a,f6.3,a,es9.2,a,f6.3,a,es9.2,a)

  print*,"nsteps =",nsteps_0,nsteps_1,nsteps_2
  print*,"nx     =",nx_0,nx_1,nx_2
  write(*,1000) "time   =",time_0,time_1,time_2
  print*,""
  write(*,1001) "$Y({\rm H}_2)$          &",L1_0( 2)," & ",log(L1_0( 2)/L1_1( 2))/log(2.d0)," &",L1_1( 2)
  write(*,1001) "$Y({\rm O}_2)$          &",L1_0( 5)," & ",log(L1_0( 5)/L1_1( 5))/log(2.d0)," &",L1_1( 5)
  write(*,1001) "$Y({\rm H}_2{\rm O})$   &",L1_0( 7)," & ",log(L1_0( 7)/L1_1( 7))/log(2.d0)," &",L1_1( 7)
  write(*,1001) "$Y({\rm H})$            &",L1_0( 3)," & ",log(L1_0( 3)/L1_1( 3))/log(2.d0)," &",L1_1( 3)
  write(*,1001) "$Y({\rm O})$            &",L1_0( 4)," & ",log(L1_0( 4)/L1_1( 4))/log(2.d0)," &",L1_1( 4)
  write(*,1001) "$Y({\rm OH})$           &",L1_0( 6)," & ",log(L1_0( 6)/L1_1( 6))/log(2.d0)," &",L1_1( 6)
  write(*,1001) "$Y({\rm HO}_2)$         &",L1_0( 8)," & ",log(L1_0( 8)/L1_1( 8))/log(2.d0)," &",L1_1( 8)
  write(*,1001) "$Y({\rm H}_2{\rm O}_2)$ &",L1_0( 9)," & ",log(L1_0( 9)/L1_1( 9))/log(2.d0)," &",L1_1( 9)
  write(*,1001) "$Y({\rm N}_2)$          &",L1_0(10)," & ",log(L1_0(10)/L1_1(10))/log(2.d0)," &",L1_1(10)
  write(*,1001) "$\rho$                  &",L1_0(11)," & ",log(L1_0(11)/L1_1(11))/log(2.d0)," &",L1_1(11)
  write(*,1001) "$T$                     &",L1_0(12)," & ",log(L1_0(12)/L1_1(12))/log(2.d0)," &",L1_1(12)
  write(*,1001) "$\rho h$                &",L1_0(13)," & ",log(L1_0(13)/L1_1(13))/log(2.d0)," &",L1_1(13)
  write(*,1001) "$U$                     &",L1_0(14)," & ",log(L1_0(14)/L1_1(14))/log(2.d0)," &",L1_1(14)
  write(*,1001) "$\nabla \cdot U$        &",L1_0(16)," & ",log(L1_0(16)/L1_1(16))/log(2.d0)," &",L1_1(16)

end program chemh_converge6
