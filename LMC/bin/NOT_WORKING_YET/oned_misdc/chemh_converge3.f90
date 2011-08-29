program chemh_converge2

  character in1*(32)
  character in2*(32)
  character in3*(32)

  integer i,j,k,nsteps_1,nsteps_2,nsteps_3,nx_1,nx_2,nx_3,rr_13,rr_23
  real*8 time_1,time_3,sum
  
  real*8 data1  (4096,26)
  real*8 data2  (4096,26)
  real*8 data3_f(4096,26)
  real*8 data3_c(4096,26)

  real*8 L0_13(26)
  real*8 L1_13(26)
  real*8 L2_13(26)
  real*8 L0_23(26)
  real*8 L1_23(26)
  real*8 L2_23(26)

  read(*,*) in1
  read(*,*) in2
  read(*,*) in3

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
     read(30,*) data3_f(i,1:26)
  end do

  !!!!!!!!!!!!!!!!!!!
  ! error for input1
  !!!!!!!!!!!!!!!!!!!

  rr_13 = nx_3 / nx_1

  !  coarsen fine data
  do i=0,nx_1-1
     do j=1,26
        sum = 0.d0
        do k=0,rr_13-1
           sum = sum + data3_f(rr_13*i+k,j)
        end do
        data3_c(i,j) = sum / dble(rr_13)
     end do
  end do
  
  L0_13 = 0.d0
  L1_13 = 0.d0
  L2_13 = 0.d0

  do i=0,nx_1-1
     do j=1,26
        L0_13(j) = max(L0_13(j), abs(data1(i,j)-data3_c(i,j)))
        L1_13(j) = L1_13(j) + abs(data1(i,j)-data3_c(i,j))
        L2_13(j) = L2_13(j) + (data1(i,j)-data3_c(i,j))**2
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
           sum = sum + data3_f(rr_23*i+k,j)
        end do
        data3_c(i,j) = sum / dble(rr_23)
     end do
  end do
  
  L0_23 = 0.d0
  L1_23 = 0.d0
  L2_23 = 0.d0

  do i=0,nx_2-1
     do j=1,26
        L0_23(j) = max(L0_23(j), abs(data2(i,j)-data3_c(i,j)))
        L1_23(j) = L1_23(j) + abs(data2(i,j)-data3_c(i,j))
        L2_23(j) = L2_23(j) + (data2(i,j)-data3_c(i,j))**2
     end do
  end do
  L1_23 = L1_23 / dble(nx_2)
  L2_23 = sqrt(L2_23/nx_2)

1000 format(a,e12.5,e12.5,e12.5)
1001 format(a,e12.3,e12.3,f7.2)

  print*,"nsteps =",nsteps_1,nsteps_2,nsteps_3
  print*,"nx     =",nx_1,nx_2,nx_3
  write(*,1000) "time   =",time_1,time_3,time_3
  print*,""
  write(*,1001) "L1 NORM: H2   =",L1_13( 2),L1_23( 2),log(L1_13( 2)/L1_23( 2))/log(2.d0)
  write(*,1001) "L1 NORM: H    =",L1_13( 3),L1_23( 3),log(L1_13( 3)/L1_23( 3))/log(2.d0)
  write(*,1001) "L1 NORM: O    =",L1_13( 4),L1_23( 4),log(L1_13( 4)/L1_23( 4))/log(2.d0)
  write(*,1001) "L1 NORM: O2   =",L1_13( 5),L1_23( 5),log(L1_13( 5)/L1_23( 5))/log(2.d0)
  write(*,1001) "L1 NORM: OH   =",L1_13( 6),L1_23( 6),log(L1_13( 6)/L1_23( 6))/log(2.d0)
  write(*,1001) "L1 NORM: H2O  =",L1_13( 7),L1_23( 7),log(L1_13( 7)/L1_23( 7))/log(2.d0)
  write(*,1001) "L1 NORM: HO2  =",L1_13( 8),L1_23( 8),log(L1_13( 8)/L1_23( 8))/log(2.d0)
  write(*,1001) "L1 NORM: H2O2 =",L1_13( 9),L1_23( 9),log(L1_13( 9)/L1_23( 9))/log(2.d0)
  write(*,1001) "L1 NORM: N2   =",L1_13(10),L1_23(10),log(L1_13(10)/L1_23(10))/log(2.d0)
  write(*,1001) "L1 NORM: rho  =",L1_13(11),L1_23(11),log(L1_13(11)/L1_23(11))/log(2.d0)
  write(*,1001) "L1 NORM: temp =",L1_13(12),L1_23(12),log(L1_13(12)/L1_23(12))/log(2.d0)
  write(*,1001) "L1 NORM: rhoh =",L1_13(13),L1_23(13),log(L1_13(13)/L1_23(13))/log(2.d0)
  write(*,1001) "L1 NORM: vel  =",L1_13(14),L1_23(14),log(L1_13(14)/L1_23(14))/log(2.d0)
  write(*,1001) "L1 NORM: pi   =",L1_13(15),L1_23(15),log(L1_13(15)/L1_23(15))/log(2.d0)
  write(*,1001) "L1 NORM: pthm =",L1_13(16),L1_23(16),log(L1_13(16)/L1_23(16))/log(2.d0)
  write(*,1001) "L1 NORM: divu =",L1_13(17),L1_23(17),log(L1_13(17)/L1_23(17))/log(2.d0)
  write(*,1001) "L1 NORM:IH2   =",L1_13(18),L1_23(18),log(L1_13(18)/L1_23(18))/log(2.d0)
  write(*,1001) "L1 NORM:IH    =",L1_13(19),L1_23(19),log(L1_13(19)/L1_23(19))/log(2.d0)
  write(*,1001) "L1 NORM:IO    =",L1_13(20),L1_23(20),log(L1_13(20)/L1_23(20))/log(2.d0)
  write(*,1001) "L1 NORM:IO2   =",L1_13(21),L1_23(21),log(L1_13(21)/L1_23(21))/log(2.d0)
  write(*,1001) "L1 NORM:IOH   =",L1_13(22),L1_23(22),log(L1_13(22)/L1_23(22))/log(2.d0)
  write(*,1001) "L1 NORM:IH2O  =",L1_13(23),L1_23(23),log(L1_13(23)/L1_23(23))/log(2.d0)
  write(*,1001) "L1 NORM:IHO2  =",L1_13(24),L1_23(24),log(L1_13(24)/L1_23(24))/log(2.d0)
  write(*,1001) "L1 NORM:IH2O2 =",L1_13(25),L1_23(25),log(L1_13(25)/L1_23(25))/log(2.d0)
  write(*,1001) "L1 NORM:IN2   =",L1_13(26),L1_23(26),log(L1_13(26)/L1_23(26))/log(2.d0)

end program chemh_converge2
