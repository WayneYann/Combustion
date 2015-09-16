program dme_converge

  implicit none

  character in0*(32)
  character in1*(32)
  character in2*(32)
  character in3*(32)

  integer i,j,k,nsteps_0,nsteps_1,nsteps_2,nsteps_3
  integer nx_0,nx_1,nx_2,nx_3
  real*8 time_0,time_1,time_2,time_3,sum
  
  real*8 data0  (0:4096,56)
  real*8 data1  (0:4096,56)
  real*8 data2  (0:4096,56)
  real*8 data3  (0:4096,56)
  real*8 data1_0(0:4096,56)
  real*8 data2_1(0:4096,46)
  real*8 data3_2(0:4096,46)

  real*8 L0_10(46)
  real*8 L1_10(46)
  real*8 L2_10(46)
  real*8 L0_21(46)
  real*8 L1_21(46)
  real*8 L2_21(46)
  real*8 L0_32(46)
  real*8 L1_32(46)
  real*8 L2_32(46)

  read(*,*) in0
  read(*,*) in1
  read(*,*) in2
  read(*,*) in3

  open(99,file=in0,form='formatted')
  read(99,*) nsteps_0
  read(99,*) nx_0
  read(99,*) time_0

  do i=0,nx_0-1
     read(99,*) data0(i,1:46)
  end do

  open(10,file=in1,form='formatted')
  read(10,*) nsteps_1
  read(10,*) nx_1
  read(10,*) time_1

  do i=0,nx_1-1
     read(10,*) data1(i,1:46)
  end do

  open(20,file=in2,form='formatted')
  read(20,*) nsteps_2
  read(20,*) nx_2
  read(20,*) time_2

  do i=0,nx_2-1
     read(20,*) data2(i,1:46)
  end do

  open(30,file=in3,form='formatted')
  read(30,*) nsteps_3
  read(30,*) nx_3
  read(30,*) time_3

  do i=0,nx_3-1
     read(30,*) data3(i,1:46)
  end do

  !!!!!!!!!!!!!!!!!!!
  ! error for input0
  !!!!!!!!!!!!!!!!!!!

  !  coarsen fine data
  do i=0,nx_0-1
     do j=1,46
        sum = 0.d0
        do k=0,1
           sum = sum + data1(2*i+k,j)
        end do
        data1_0(i,j) = sum / 2.d0
     end do
  end do

  j=44  
  do i=0,nx_1-1
     data1_0(i,j) = data1(i*2,j)
  end do
  
  L0_10 = 0.d0
  L1_10 = 0.d0
  L2_10 = 0.d0

  do i=0,nx_0-1
     do j=1,46
        L0_10(j) = max(L0_10(j), abs(data0(i,j)-data1_0(i,j)))
        L1_10(j) = L1_10(j) + abs(data0(i,j)-data1_0(i,j))
        L2_10(j) = L2_10(j) + (data0(i,j)-data1_0(i,j))**2
     end do
  end do
  L1_10 = L1_10 / dble(nx_0)
  L2_10 = sqrt(L2_10/nx_0)

  !!!!!!!!!!!!!!!!!!!
  ! error for input1
  !!!!!!!!!!!!!!!!!!!

  !  coarsen fine data
  do i=0,nx_1-1
     do j=1,46
        sum = 0.d0
        do k=0,1
           sum = sum + data2(2*i+k,j)
        end do
        data2_1(i,j) = sum / 2.d0
     end do
  end do

  j=44  
  do i=0,nx_2-1
     data2_1(i,j) = data2(i*2,j)
  end do
  
  L0_21 = 0.d0
  L1_21 = 0.d0
  L2_21 = 0.d0

  do i=0,nx_1-1
     do j=1,46
        L0_21(j) = max(L0_21(j), abs(data1(i,j)-data2_1(i,j)))
        L1_21(j) = L1_21(j) + abs(data1(i,j)-data2_1(i,j))
        L2_21(j) = L2_21(j) + (data1(i,j)-data2_1(i,j))**2
     end do
  end do
  L1_21 = L1_21 / dble(nx_1)
  L2_21 = sqrt(L2_21/nx_1)

  !!!!!!!!!!!!!!!!!!!
  ! error for input2
  !!!!!!!!!!!!!!!!!!!

  !  coarsen fine data
  do i=0,nx_2-1
     do j=1,46
        sum = 0.d0
        do k=0,1
           sum = sum + data3(2*i+k,j)
        end do
        data3_2(i,j) = sum / 2.d0
     end do
  end do

  j=44  
  do i=0,nx_3-1
     data3_2(i,j) = data3(i*2,j)
  end do
  
  L0_32 = 0.d0
  L1_32 = 0.d0
  L2_32 = 0.d0

  do i=0,nx_2-1
     do j=1,46
        L0_32(j) = max(L0_32(j), abs(data2(i,j)-data3_2(i,j)))
        L1_32(j) = L1_32(j) + abs(data2(i,j)-data3_2(i,j))
        L2_32(j) = L2_32(j) + (data2(i,j)-data3_2(i,j))**2
     end do
  end do
  L1_32 = L1_32 / dble(nx_2)
  L2_32 = sqrt(L2_32/nx_2)

1000 format(a,es12.4,es12.4,es12.4,es12.4)
1001 format(a,es9.2,a,f4.2,a,es9.2,a,f4.2,a,es9.2,a)

  print*,"nsteps =",nsteps_0,nsteps_1,nsteps_2,nsteps_3
  print*,"nx     =",nx_0,nx_1,nx_2,nx_3
  write(*,1000) "time   =",time_0,time_1,time_2,time_3
  print*,""
  write(*,1001) "$Y({\rm CH}_3{\rm OCH}_3)$           &",L1_10(28)," & ",log(L1_10(28)/L1_21(28))/log(2.d0)," &",L1_21(28)," & ", &
       log(L1_21(28)/L1_32(28))/log(2.d0)," &",L1_32(28)," \\"
  write(*,1001) "$Y({\rm O}_2)$                       &",L1_10(21)," & ",log(L1_10(21)/L1_21(21))/log(2.d0)," &",L1_21(21)," & ", &
       log(L1_21(21)/L1_32(21))/log(2.d0)," &",L1_32(21)," \\"
  write(*,1001) "$Y({\rm CO}_2)$                      &",L1_10(24)," & ",log(L1_10(24)/L1_21(24))/log(2.d0)," &",L1_21(24)," & ", &
       log(L1_21(24)/L1_32(24))/log(2.d0)," &",L1_32(24)," \\"
  write(*,1001) "$Y({\rm H}_2{\rm O})$                &",L1_10(10)," & ",log(L1_10(10)/L1_21(10))/log(2.d0)," &",L1_21(10)," & ", &
       log(L1_21(10)/L1_32(10))/log(2.d0)," &",L1_32(10)," \\"
  write(*,1001) "$Y({\rm CH}_3{\rm OCH}_2{\rm O}_2)$  &",L1_10(34)," & ",log(L1_10(34)/L1_21(34))/log(2.d0)," &",L1_21(34)," & ", &
       log(L1_21(34)/L1_32(34))/log(2.d0)," &",L1_32(34)," \\"
  write(*,1001) "$Y({\rm OH})$                        &",L1_10( 9)," & ",log(L1_10( 9)/L1_21( 9))/log(2.d0)," &",L1_21( 9)," & ", &
       log(L1_21( 9)/L1_32( 9))/log(2.d0)," &",L1_32( 9)," \\"
  write(*,1001) "$Y({\rm HO}_2)$                      &",L1_10(22)," & ",log(L1_10(22)/L1_21(22))/log(2.d0)," &",L1_21(22)," & ", &
       log(L1_21(22)/L1_32(22))/log(2.d0)," &",L1_32(22)," \\"
  write(*,1001) "$Y({\rm O})$                         &",L1_10( 7)," & ",log(L1_10( 7)/L1_21( 7))/log(2.d0)," &",L1_21( 7)," & ", &
       log(L1_21( 7)/L1_32( 7))/log(2.d0)," &",L1_32( 7)," \\"
  write(*,1001) "$Y({\rm H})$                         &",L1_10( 2)," & ",log(L1_10( 2)/L1_21( 2))/log(2.d0)," &",L1_21( 2)," & ", &
       log(L1_21( 2)/L1_32( 2))/log(2.d0)," &",L1_32( 2)," \\"
  write(*,1001) "$Y({\rm N}_2)$                       &",L1_10(40)," & ",log(L1_10(40)/L1_21(40))/log(2.d0)," &",L1_21(40)," & ", &
       log(L1_21(40)/L1_32(40))/log(2.d0)," &",L1_32(40)," \\"
  write(*,1001) "$\rho$                               &",L1_10(41)," & ",log(L1_10(41)/L1_21(41))/log(2.d0)," &",L1_21(41)," & ", &
       log(L1_21(41)/L1_32(41))/log(2.d0)," &",L1_32(41)," \\"
  write(*,1001) "$T$                                  &",L1_10(42)," & ",log(L1_10(42)/L1_21(42))/log(2.d0)," &",L1_21(42)," & ", &
       log(L1_21(42)/L1_32(42))/log(2.d0)," &",L1_32(42)," \\"
  write(*,1001) "$\rho h$                             &",L1_10(43)," & ",log(L1_10(43)/L1_21(43))/log(2.d0)," &",L1_21(43)," & ", &
       log(L1_21(43)/L1_32(43))/log(2.d0)," &",L1_32(43)," \\"
  write(*,1001) "$U$                                  &",L1_10(44)," & ",log(L1_10(44)/L1_21(44))/log(2.d0)," &",L1_21(44)," & ", &
       log(L1_21(44)/L1_32(44))/log(2.d0)," &",L1_32(44)," \\"
  write(*,1001) "$\nabla \cdot U$                     &",L1_10(46)," & ",log(L1_10(46)/L1_21(46))/log(2.d0)," &",L1_21(46)," & ", &
       log(L1_21(46)/L1_32(46))/log(2.d0)," &",L1_32(46)," \\"

end program dme_converge
