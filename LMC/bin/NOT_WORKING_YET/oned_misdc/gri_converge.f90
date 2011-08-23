program gri_converge

  character in1*(32)
  character in2*(32)

  integer i,j,k,nsteps_c,nsteps_f,nx_c,nx_f,rr
  real*8 time_c,time_f,sum
  
  real*8 data1  (1024,114)
  real*8 data2_f(4096,114)
  real*8 data2_c(1024,114)

  real*8 L0(114)
  real*8 L1(114)
  real*8 L2(114)

  read(*,*) in1
  read(*,*) in2

  open(10,file=in1,form='formatted')
  read(10,*) nsteps_c
  read(10,*) nx_c
  read(10,*) time_c

  do i=0,nx_c-1
     read(10,*) data1(i,1:114)
  end do

  open(20,file=in2,form='formatted')
  read(20,*) nsteps_f
  read(20,*) nx_f
  read(20,*) time_f

  do i=0,nx_f-1
     read(20,*) data2_f(i,1:114)
  end do

  rr = nx_f / nx_c

  !  coarsen fine data
  do i=0,nx_c-1
     do j=1,114
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
     do j=1,114
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
  write(*,*) "L1 NORM: H2     =",L1( 2)
  write(*,*) "L1 NORM: H      =",L1( 3)
  write(*,*) "L1 NORM: O      =",L1( 4)
  write(*,*) "L1 NORM: O2     =",L1( 5)
  write(*,*) "L1 NORM: OH     =",L1( 6)
  write(*,*) "L1 NORM: H2O    =",L1( 7)
  write(*,*) "L1 NORM: HO2    =",L1( 8)
  write(*,*) "L1 NORM: H2O2   =",L1( 9)
  write(*,*) "L1 NORM: C      =",L1(10)
  write(*,*) "L1 NORM: CH     =",L1(11)
  write(*,*) "L1 NORM: CH2    =",L1(12)
  write(*,*) "L1 NORM: CH2(S) =",L1(13)
  write(*,*) "L1 NORM: CH3    =",L1(14)
  write(*,*) "L1 NORM: CH4    =",L1(15)
  write(*,*) "L1 NORM: CO     =",L1(16)
  write(*,*) "L1 NORM: CO2    =",L1(17)
  write(*,*) "L1 NORM: HCO    =",L1(18)
  write(*,*) "L1 NORM: CH2O   =",L1(19)
  write(*,*) "L1 NORM: CH2OH  =",L1(20)
  write(*,*) "L1 NORM: CH3O   =",L1(21)
  write(*,*) "L1 NORM: CH3OH  =",L1(22)
  write(*,*) "L1 NORM: C2H    =",L1(23)
  write(*,*) "L1 NORM: C2H2   =",L1(24)
  write(*,*) "L1 NORM: C2H3   =",L1(25)
  write(*,*) "L1 NORM: C2H4   =",L1(26)
  write(*,*) "L1 NORM: C2H5   =",L1(27)
  write(*,*) "L1 NORM: C2H6   =",L1(28)
  write(*,*) "L1 NORM: HCCO   =",L1(29)
  write(*,*) "L1 NORM: CH2CO  =",L1(30)
  write(*,*) "L1 NORM: HCCOH  =",L1(31)
  write(*,*) "L1 NORM: N      =",L1(32)
  write(*,*) "L1 NORM: NH     =",L1(33)
  write(*,*) "L1 NORM: NH2    =",L1(34)
  write(*,*) "L1 NORM: NH3    =",L1(35)
  write(*,*) "L1 NORM: NH4    =",L1(36)
  write(*,*) "L1 NORM: NO     =",L1(37)
  write(*,*) "L1 NORM: NO2    =",L1(38)
  write(*,*) "L1 NORM: N2O    =",L1(39)
  write(*,*) "L1 NORM: HNO    =",L1(40)
  write(*,*) "L1 NORM: CN     =",L1(41)
  write(*,*) "L1 NORM: HCN    =",L1(42)
  write(*,*) "L1 NORM: H2CN   =",L1(43)
  write(*,*) "L1 NORM: HCNN   =",L1(44)
  write(*,*) "L1 NORM: HCNO   =",L1(45)
  write(*,*) "L1 NORM: HOCN   =",L1(46)
  write(*,*) "L1 NORM: HNCO   =",L1(47)
  write(*,*) "L1 NORM: NCO    =",L1(48)
  write(*,*) "L1 NORM: N2     =",L1(49)
  write(*,*) "L1 NORM: AR     =",L1(50)
  write(*,*) "L1 NORM: C3H7   =",L1(51)
  write(*,*) "L1 NORM: C3H8   =",L1(52)
  write(*,*) "L1 NORM: CH2CHO =",L1(53)
  write(*,*) "L1 NORM: CH3CHO =",L1(54)
  write(*,*) "L1 NORM: rho    =",L1(55)
  write(*,*) "L1 NORM: temp   =",L1(56)
  write(*,*) "L1 NORM: rhoh   =",L1(57)
  write(*,*) "L1 NORM: vel    =",L1(58)

end program gri_converge
