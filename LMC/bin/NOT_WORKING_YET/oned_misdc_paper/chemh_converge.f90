program chemh_converge

  character in1*(32)
  character in2*(32)

  integer i,j,k,nsteps_c,nsteps_f,nx_c,nx_f,rr
  real*8 time_c,time_f,sum
  
  real*8 data1  (4096,26)
  real*8 data2_f(4096,26)
  real*8 data2_c(4096,26)

  real*8 L0(26)
  real*8 L1(26)
  real*8 L2(26)

  read(*,*) in1
  read(*,*) in2

  open(10,file=in1,form='formatted')
  read(10,*) nsteps_c
  read(10,*) nx_c
  read(10,*) time_c

  do i=0,nx_c-1
     read(10,*) data1(i,1:26)
  end do

  open(20,file=in2,form='formatted')
  read(20,*) nsteps_f
  read(20,*) nx_f
  read(20,*) time_f

  do i=0,nx_f-1
     read(20,*) data2_f(i,1:26)
  end do

  rr = nx_f / nx_c

  !  coarsen fine data
  do i=0,nx_c-1
     do j=1,26
        sum = 0.d0
        do k=0,rr-1
           sum = sum + data2_f(rr*i+k,j)
        end do
        data2_c(i,j) = sum / dble(rr)
     end do
  end do
  
  L0 = 0.d0
  L1 = 0.d0
  L2 = 0.d0

  do i=0,nx_c-1
     do j=1,26
        L0(j) = max(L0(j), abs(data1(i,j)-data2_c(i,j)))
        L1(j) = L1(j) + abs(data1(i,j)-data2_c(i,j))
        L2(j) = L2(j) + (data1(i,j)-data2_c(i,j))**2
     end do
  end do
  L1 = L1 / dble(nx_c)
  L2 = sqrt(L2/nx_c)

  write(*,*) "nsteps =",nsteps_c,nsteps_f
  write(*,*) "nx     =",nx_c,nx_f
  write(*,*) "time   =",time_c,time_f
  write(*,*)
  write(*,*) "L0 NORM: H2   =",L0( 2)
  write(*,*) "L0 NORM: H    =",L0( 3)
  write(*,*) "L0 NORM: O    =",L0( 4)
  write(*,*) "L0 NORM: O2   =",L0( 5)
  write(*,*) "L0 NORM: OH   =",L0( 6)
  write(*,*) "L0 NORM: H2O  =",L0( 7)
  write(*,*) "L0 NORM: HO2  =",L0( 8)
  write(*,*) "L0 NORM: H2O2 =",L0( 9)
  write(*,*) "L0 NORM: N2   =",L0(10)
  write(*,*) "L0 NORM: rho  =",L0(11)
  write(*,*) "L0 NORM: temp =",L0(12)
  write(*,*) "L0 NORM: rhoh =",L0(13)
  write(*,*) "L0 NORM: vel  =",L0(14)
  write(*,*)
  write(*,*) "L1 NORM: H2   =",L1( 2)
  write(*,*) "L1 NORM: H    =",L1( 3)
  write(*,*) "L1 NORM: O    =",L1( 4)
  write(*,*) "L1 NORM: O2   =",L1( 5)
  write(*,*) "L1 NORM: OH   =",L1( 6)
  write(*,*) "L1 NORM: H2O  =",L1( 7)
  write(*,*) "L1 NORM: HO2  =",L1( 8)
  write(*,*) "L1 NORM: H2O2 =",L1( 9)
  write(*,*) "L1 NORM: N2   =",L1(10)
  write(*,*) "L1 NORM: rho  =",L1(11)
  write(*,*) "L1 NORM: temp =",L1(12)
  write(*,*) "L1 NORM: rhoh =",L1(13)
  write(*,*) "L1 NORM: vel  =",L1(14)
  write(*,*)
  write(*,*) "L2 NORM: H2   =",L2( 2)
  write(*,*) "L2 NORM: H    =",L2( 3)
  write(*,*) "L2 NORM: O    =",L2( 4)
  write(*,*) "L2 NORM: O2   =",L2( 5)
  write(*,*) "L2 NORM: OH   =",L2( 6)
  write(*,*) "L2 NORM: H2O  =",L2( 7)
  write(*,*) "L2 NORM: HO2  =",L2( 8)
  write(*,*) "L2 NORM: H2O2 =",L2( 9)
  write(*,*) "L2 NORM: N2   =",L2(10)
  write(*,*) "L2 NORM: rho  =",L2(11)
  write(*,*) "L2 NORM: temp =",L2(12)
  write(*,*) "L2 NORM: rhoh =",L2(13)
  write(*,*) "L2 NORM: vel  =",L2(14)

end program chemh_converge
