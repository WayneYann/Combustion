program chemh_converge4

  implicit none

  character in0*(32)
  character in1*(32)
  character in2*(32)
  character in3*(32)

  integer i,j,k,nsteps_0,nsteps_1,nsteps_2,nsteps_3
  integer nx_0,nx_1,nx_2,nx_3,rr_03,rr_13,rr_23
  real*8 time_0,time_1,time_2,time_3,sum
  
  real*8 data0  (4096,26)
  real*8 data1  (4096,26)
  real*8 data2  (4096,26)
  real*8 data3  (4096,26)
  real*8 data3_0(4096,26)
  real*8 data3_1(4096,26)
  real*8 data3_2(4096,26)

  real*8 L0_03(26)
  real*8 L1_03(26)
  real*8 L2_03(26)
  real*8 L0_13(26)
  real*8 L1_13(26)
  real*8 L2_13(26)
  real*8 L0_23(26)
  real*8 L1_23(26)
  real*8 L2_23(26)

  read(*,*) in0
  read(*,*) in1
  read(*,*) in2
  read(*,*) in3

  open(99,file=in0,form='formatted')
  read(99,*) nsteps_0
  read(99,*) nx_0
  read(99,*) time_0

  do i=0,nx_0-1
     read(99,*) data0(i,1:26)
  end do

  open(10,file=in1,form='formatted')
  read(10,*) nsteps_1
  read(10,*) nx_1
  read(10,*) time_1

  do i=0,nx_1-1
     read(10,*) data1(i,1:26)
  end do

  open(20,file=in2,form='formatted')
  read(20,*) nsteps_2
  read(20,*) nx_2
  read(20,*) time_2

  do i=0,nx_2-1
     read(20,*) data2(i,1:26)
  end do

  open(30,file=in3,form='formatted')
  read(30,*) nsteps_3
  read(30,*) nx_3
  read(30,*) time_3

  do i=0,nx_3-1
     read(30,*) data3(i,1:26)
  end do

  !!!!!!!!!!!!!!!!!!!
  ! error for input0
  !!!!!!!!!!!!!!!!!!!

  rr_03 = nx_3 / nx_0

  !  coarsen fine data
  do i=0,nx_0-1
     do j=1,26
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
     do j=1,26
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
     do j=1,26
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
     do j=1,26
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
     do j=1,26
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
     do j=1,26
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
  write(*,1001) "$Y({\rm H}_2)$          &",L1_03( 2)," & ",log(L1_03( 2)/L1_13( 2))/log(2.d0)," &",L1_13( 2)," & ", &
       log(L1_13( 2)/L1_23( 2))/log(2.d0)," &",L1_23( 2)," \\"
  write(*,1001) "$Y({\rm H})$            &",L1_03( 3)," & ",log(L1_03( 3)/L1_13( 3))/log(2.d0)," &",L1_13( 3)," & ", &
       log(L1_13( 3)/L1_23( 3))/log(2.d0)," &",L1_23( 3)," \\"
  write(*,1001) "$Y({\rm H})$            &",L1_03( 4)," & ",log(L1_03( 4)/L1_13( 4))/log(2.d0)," &",L1_13( 4)," & ", &
       log(L1_13( 4)/L1_23( 4))/log(2.d0)," &",L1_23( 4)," \\"
  write(*,1001) "$Y({\rm O}_2)$          &",L1_03( 5)," & ",log(L1_03( 5)/L1_13( 5))/log(2.d0)," &",L1_13( 5)," & ", &
       log(L1_13( 5)/L1_23( 5))/log(2.d0)," &",L1_23( 5)," \\"
  write(*,1001) "$Y({\rm OH})$           &",L1_03( 6)," & ",log(L1_03( 6)/L1_13( 6))/log(2.d0)," &",L1_13( 6)," & ", &
       log(L1_13( 6)/L1_23( 6))/log(2.d0)," &",L1_23( 6)," \\"
  write(*,1001) "$Y({\rm H}_2{\rm O})$   &",L1_03( 7)," & ",log(L1_03( 7)/L1_13( 7))/log(2.d0)," &",L1_13( 7)," & ", &
       log(L1_13( 7)/L1_23( 7))/log(2.d0)," &",L1_23( 7)," \\"
  write(*,1001) "$Y({\rm HO}_2)$         &",L1_03( 8)," & ",log(L1_03( 8)/L1_13( 8))/log(2.d0)," &",L1_13( 8)," & ", &
       log(L1_13( 8)/L1_23( 8))/log(2.d0)," &",L1_23( 8)," \\"
  write(*,1001) "$Y({\rm H}_2{\rm O}_2)$ &",L1_03( 9)," & ",log(L1_03( 9)/L1_13( 9))/log(2.d0)," &",L1_13( 9)," & ", &
       log(L1_13( 9)/L1_23( 9))/log(2.d0)," &",L1_23( 9)," \\"
  write(*,1001) "$Y({\rm N}_2)$          &",L1_03(10)," & ",log(L1_03(10)/L1_13(10))/log(2.d0)," &",L1_13(10)," & ", &
       log(L1_13(10)/L1_23(10))/log(2.d0)," &",L1_23(10)," \\"
  write(*,1001) "$\rho$                  &",L1_03(11)," & ",log(L1_03(11)/L1_13(11))/log(2.d0)," &",L1_13(11)," & ", &
       log(L1_13(11)/L1_23(11))/log(2.d0)," &",L1_23(11)," \\"
  write(*,1001) "$T$                     &",L1_03(12)," & ",log(L1_03(12)/L1_13(12))/log(2.d0)," &",L1_13(12)," & ", &
       log(L1_13(12)/L1_23(12))/log(2.d0)," &",L1_23(12)," \\"
  write(*,1001) "$\rho h$                &",L1_03(13)," & ",log(L1_03(13)/L1_13(13))/log(2.d0)," &",L1_13(13)," & ", &
       log(L1_13(13)/L1_23(13))/log(2.d0)," &",L1_23(13)," \\"
  write(*,1001) "$U$                     &",L1_03(14)," & ",log(L1_03(14)/L1_13(14))/log(2.d0)," &",L1_13(14)," & ", &
       log(L1_13(14)/L1_23(14))/log(2.d0)," &",L1_23(14)," \\"

end program chemh_converge4
