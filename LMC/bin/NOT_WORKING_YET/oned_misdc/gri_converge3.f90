program gri_converge3

  implicit none

  character in0*(32)
  character in1*(32)
  character in2*(32)
  character in3*(32)

  integer i,j,k,nsteps_0,nsteps_1,nsteps_2,nsteps_3
  integer nx_0,nx_1,nx_2,nx_3,rr_03,rr_13,rr_23
  real*8 time_0,time_1,time_2,time_3,sum
  
  real*8 data0  (4096,114)
  real*8 data1  (4096,114)
  real*8 data2  (4096,114)
  real*8 data3  (4096,114)
  real*8 data3_0(4096,114)
  real*8 data3_1(4096,114)
  real*8 data3_2(4096,114)

  real*8 L0_03(114)
  real*8 L1_03(114)
  real*8 L2_03(114)
  real*8 L0_13(114)
  real*8 L1_13(114)
  real*8 L2_13(114)
  real*8 L0_23(114)
  real*8 L1_23(114)
  real*8 L2_23(114)

  read(*,*) in0
  read(*,*) in1
  read(*,*) in2
  read(*,*) in3

  open(99,file=in0,form='formatted')
  read(99,*) nsteps_0
  read(99,*) nx_0
  read(99,*) time_0

  do i=0,nx_0-1
     read(99,*) data0(i,1:114)
  end do

  open(10,file=in1,form='formatted')
  read(10,*) nsteps_1
  read(10,*) nx_1
  read(10,*) time_1

  do i=0,nx_1-1
     read(10,*) data1(i,1:114)
  end do

  open(20,file=in2,form='formatted')
  read(20,*) nsteps_2
  read(20,*) nx_2
  read(20,*) time_2

  do i=0,nx_2-1
     read(20,*) data2(i,1:114)
  end do

  open(30,file=in3,form='formatted')
  read(30,*) nsteps_3
  read(30,*) nx_3
  read(30,*) time_3

  do i=0,nx_3-1
     read(30,*) data3(i,1:114)
  end do

  !!!!!!!!!!!!!!!!!!!
  ! error for input0
  !!!!!!!!!!!!!!!!!!!

  rr_03 = nx_3 / nx_0

  !  coarsen fine data
  do i=0,nx_0-1
     do j=1,114
        sum = 0.d0
        do k=0,rr_03-1
           sum = sum + data3(rr_03*i+k,j)
        end do
        data3_0(i,j) = sum / dble(rr_03)
     end do
  end do
  
  L0_03 = 0.d0
  L1_03 = 0.d0
  L2_03 = 0.d0

  do i=0,nx_0-1
     do j=1,114
        L0_03(j) = max(L0_03(j), abs(data0(i,j)-data3_0(i,j)))
        L1_03(j) = L1_03(j) + abs(data0(i,j)-data3_0(i,j))
        L2_03(j) = L2_03(j) + (data0(i,j)-data3_0(i,j))**2
     end do
  end do
  L1_03 = L1_03 / dble(nx_0)
  L2_03 = sqrt(L2_03/nx_0)

  !!!!!!!!!!!!!!!!!!!
  ! error for input1
  !!!!!!!!!!!!!!!!!!!

  rr_13 = nx_3 / nx_1

  !  coarsen fine data
  do i=0,nx_1-1
     do j=1,114
        sum = 0.d0
        do k=0,rr_13-1
           sum = sum + data3(rr_13*i+k,j)
        end do
        data3_1(i,j) = sum / dble(rr_13)
     end do
  end do
  
  L0_13 = 0.d0
  L1_13 = 0.d0
  L2_13 = 0.d0

  do i=0,nx_1-1
     do j=1,114
        L0_13(j) = max(L0_13(j), abs(data1(i,j)-data3_1(i,j)))
        L1_13(j) = L1_13(j) + abs(data1(i,j)-data3_1(i,j))
        L2_13(j) = L2_13(j) + (data1(i,j)-data3_1(i,j))**2
     end do
  end do
  L1_13 = L1_13 / dble(nx_1)
  L2_13 = sqrt(L2_13/nx_1)

  !!!!!!!!!!!!!!!!!!!
  ! error for input2
  !!!!!!!!!!!!!!!!!!!

  rr_23 = nx_3 / nx_2

  !  coarsen fine data
  do i=0,nx_2-1
     do j=1,114
        sum = 0.d0
        do k=0,rr_23-1
           sum = sum + data3(rr_23*i+k,j)
        end do
        data3_1(i,j) = sum / dble(rr_23)
     end do
  end do
  
  L0_23 = 0.d0
  L1_23 = 0.d0
  L2_23 = 0.d0

  do i=0,nx_2-1
     do j=1,114
        L0_23(j) = max(L0_23(j), abs(data2(i,j)-data3_1(i,j)))
        L1_23(j) = L1_23(j) + abs(data2(i,j)-data3_1(i,j))
        L2_23(j) = L2_23(j) + (data2(i,j)-data3_1(i,j))**2
     end do
  end do
  L1_23 = L1_23 / dble(nx_2)
  L2_23 = sqrt(L2_23/nx_2)

1000 format(a,es12.4,es12.4,es12.4,es12.4)
1001 format(a,es9.2,a,f4.2,a,es9.2,a,f4.2,a,es9.2,a)

  print*,"nsteps =",nsteps_0,nsteps_1,nsteps_2,nsteps_3
  print*,"nx     =",nx_0,nx_1,nx_2,nx_3
  write(*,1000) "time   =",time_0,time_1,time_2,time_3
  print*,""
  write(*,1001) "$Y({\rm O}_2)$        &",L1_03( 5)," & ",log(L1_03( 5)/L1_13( 5))/log(2.d0)," &",L1_13( 5)," & ", &
       log(L1_13( 5)/L1_23( 5))/log(2.d0)," &",L1_23( 5)," \\"
  write(*,1001) "$Y({\rm OH})$         &",L1_03( 6)," & ",log(L1_03( 6)/L1_13( 6))/log(2.d0)," &",L1_13( 6)," & ", &
       log(L1_13( 6)/L1_23( 6))/log(2.d0)," &",L1_23( 6)," \\"
  write(*,1001) "$Y({\rm H}_2{\rm O})$ &",L1_03( 7)," & ",log(L1_03( 7)/L1_13( 7))/log(2.d0)," &",L1_13( 7)," & ", &
       log(L1_13( 7)/L1_23( 7))/log(2.d0)," &",L1_23( 7)," \\"
  write(*,1001) "$Y({\rm CH}_4)$       &",L1_03(15)," & ",log(L1_03(15)/L1_13(15))/log(2.d0)," &",L1_13(15)," & ", &
       log(L1_13(15)/L1_23(15))/log(2.d0)," &",L1_23(15)," \\"
  write(*,1001) "$Y({\rm CO})$         &",L1_03(16)," & ",log(L1_03(16)/L1_13(16))/log(2.d0)," &",L1_13(16)," & ", &
       log(L1_13(16)/L1_23(16))/log(2.d0)," &",L1_23(16)," \\"
  write(*,1001) "$Y({\rm CO}_2)$       &",L1_03(17)," & ",log(L1_03(17)/L1_13(17))/log(2.d0)," &",L1_13(17)," & ", &
       log(L1_13(17)/L1_23(17))/log(2.d0)," &",L1_23(17)," \\"
  write(*,1001) "$Y({\rm N}_2)$        &",L1_03(49)," & ",log(L1_03(49)/L1_13(49))/log(2.d0)," &",L1_13(49)," & ", &
       log(L1_13(49)/L1_23(49))/log(2.d0)," &",L1_23(49)," \\"
  write(*,1001) "$\rho$                &",L1_03(55)," & ",log(L1_03(55)/L1_13(55))/log(2.d0)," &",L1_13(55)," & ", &
       log(L1_13(55)/L1_23(55))/log(2.d0)," &",L1_23(55)," \\"
  write(*,1001) "$T$                   &",L1_03(56)," & ",log(L1_03(56)/L1_13(56))/log(2.d0)," &",L1_13(56)," & ", &
       log(L1_13(56)/L1_23(56))/log(2.d0)," &",L1_23(56)," \\"
  write(*,1001) "$\rho h$              &",L1_03(57)," & ",log(L1_03(57)/L1_13(57))/log(2.d0)," &",L1_13(57)," & ", &
       log(L1_13(57)/L1_23(57))/log(2.d0)," &",L1_23(57)," \\"
  write(*,1001) "$U$                   &",L1_03(58)," & ",log(L1_03(58)/L1_13(58))/log(2.d0)," &",L1_13(58)," & ", &
       log(L1_13(58)/L1_23(58))/log(2.d0)," &",L1_23(58)," \\"

end program gri_converge3
