program gri_analysis

  implicit none

  character in1*(32)

  integer i,j,k,nsteps_1,nx_1
  real*8 time_1
  
  real*8 data1  (4096,114)

  real*8 sum(114)

  real*8 L0_13(114)
  real*8 L1_13(114)
  real*8 L2_13(114)
  real*8 L0_23(114)
  real*8 L1_23(114)
  real*8 L2_23(114)

  read(*,*) in1

  open(10,file=in1,form='formatted')
  read(10,*) nsteps_1
  read(10,*) nx_1
  read(10,*) time_1

  do i=0,nx_1-1
     read(10,*) data1(i,1:114)
  end do

  !!!!!!!!!!!!!!!!!!!
  ! error for input1
  !!!!!!!!!!!!!!!!!!!

  sum = 0.d0

  do i=0,nx_1-1
     do j=1,114
        sum(j) = sum(j) + data1(i,j)
     end do
  end do

1000 format(a,e12.5,e12.5,e12.5)
1001 format(a,e12.3)

  print*,"nsteps =",nsteps_1
  print*,"nx     =",nx_1
  write(*,1000) "time   =",time_1
  print*,""
  write(*,1001) "SUM H2     =",sum(2)
  write(*,1001) "SUM H      =",sum(3)
  write(*,1001) "SUM O      =",sum(4)
  write(*,1001) "SUM O2     =",sum(5)
  write(*,1001) "SUM OH     =",sum(6)
  write(*,1001) "SUM H2O    =",sum(7)
  write(*,1001) "SUM HO2    =",sum(8)
  write(*,1001) "SUM H2O2   =",sum(9)
  write(*,1001) "SUM C      =",sum(10)
  write(*,1001) "SUM CH     =",sum(11)
  write(*,1001) "SUM CH2    =",sum(12)
  write(*,1001) "SUM CH2(S) =",sum(13)
  write(*,1001) "SUM CH3    =",sum(14)
  write(*,1001) "SUM CH4    =",sum(15)
  write(*,1001) "SUM CO     =",sum(16)
  write(*,1001) "SUM CO2    =",sum(17)
  write(*,1001) "SUM HCO    =",sum(18)
  write(*,1001) "SUM CH2O   =",sum(19)
  write(*,1001) "SUM CH2OH  =",sum(20)
  write(*,1001) "SUM CH3O   =",sum(21)
  write(*,1001) "SUM CH3OH  =",sum(22)
  write(*,1001) "SUM C2H    =",sum(23)
  write(*,1001) "SUM C2H2   =",sum(24)
  write(*,1001) "SUM C2H3   =",sum(25)
  write(*,1001) "SUM C2H4   =",sum(26)
  write(*,1001) "SUM C2H5   =",sum(27)
  write(*,1001) "SUM C2H6   =",sum(28)
  write(*,1001) "SUM HCCO   =",sum(29)
  write(*,1001) "SUM CH2CO  =",sum(30)
  write(*,1001) "SUM HCCOH  =",sum(31)
  write(*,1001) "SUM N      =",sum(32)
  write(*,1001) "SUM NH     =",sum(33)
  write(*,1001) "SUM NH2    =",sum(34)
  write(*,1001) "SUM NH3    =",sum(35)
  write(*,1001) "SUM NH4    =",sum(36)
  write(*,1001) "SUM NO     =",sum(37)
  write(*,1001) "SUM NO2    =",sum(38)
  write(*,1001) "SUM N2O    =",sum(39)
  write(*,1001) "SUM HNO    =",sum(40)
  write(*,1001) "SUM CN     =",sum(41)
  write(*,1001) "SUM HCN    =",sum(42)
  write(*,1001) "SUM H2CN   =",sum(43)
  write(*,1001) "SUM HCNN   =",sum(44)
  write(*,1001) "SUM HCNO   =",sum(45)
  write(*,1001) "SUM HOCN   =",sum(46)
  write(*,1001) "SUM HNCO   =",sum(47)
  write(*,1001) "SUM NCO    =",sum(48)
  write(*,1001) "SUM N2     =",sum(49)
  write(*,1001) "SUM AR     =",sum(50)
  write(*,1001) "SUM C3H7   =",sum(51)
  write(*,1001) "SUM C3H8   =",sum(52)
  write(*,1001) "SUM CH2CHO =",sum(53)
  write(*,1001) "SUM CH3CHO =",sum(54)
  write(*,1001) "SUM IH2    =",sum(60+2)
  write(*,1001) "SUM IH     =",sum(60+3)
  write(*,1001) "SUM IO     =",sum(60+4)
  write(*,1001) "SUM IO2    =",sum(60+5)
  write(*,1001) "SUM IOH    =",sum(60+6)
  write(*,1001) "SUM IH2O   =",sum(60+7)
  write(*,1001) "SUM IHO2   =",sum(60+8)
  write(*,1001) "SUM IH2O2  =",sum(60+9)
  write(*,1001) "SUM IC     =",sum(60+10)
  write(*,1001) "SUM ICH    =",sum(60+11)
  write(*,1001) "SUM ICH2   =",sum(60+12)
  write(*,1001) "SUM ICH2(S)=",sum(60+13)
  write(*,1001) "SUM ICH3   =",sum(60+14)
  write(*,1001) "SUM ICH4   =",sum(60+15)
  write(*,1001) "SUM ICO    =",sum(60+16)
  write(*,1001) "SUM ICO2   =",sum(60+17)
  write(*,1001) "SUM IHCO   =",sum(60+18)
  write(*,1001) "SUM ICH2O  =",sum(60+19)
  write(*,1001) "SUM ICH2OH =",sum(60+20)
  write(*,1001) "SUM ICH3O  =",sum(60+21)
  write(*,1001) "SUM ICH3OH =",sum(60+22)
  write(*,1001) "SUM IC2H   =",sum(60+23)
  write(*,1001) "SUM IC2H2  =",sum(60+24)
  write(*,1001) "SUM IC2H3  =",sum(60+25)
  write(*,1001) "SUM IC2H4  =",sum(60+26)
  write(*,1001) "SUM IC2H5  =",sum(60+27)
  write(*,1001) "SUM IC2H6  =",sum(60+28)
  write(*,1001) "SUM IHCCO  =",sum(60+29)
  write(*,1001) "SUM ICH2CO =",sum(60+30)
  write(*,1001) "SUM IHCCOH =",sum(60+31)
  write(*,1001) "SUM IN     =",sum(60+32)
  write(*,1001) "SUM INH    =",sum(60+33)
  write(*,1001) "SUM INH2   =",sum(60+34)
  write(*,1001) "SUM INH3   =",sum(60+35)
  write(*,1001) "SUM INH4   =",sum(60+36)
  write(*,1001) "SUM INO    =",sum(60+37)
  write(*,1001) "SUM INO2   =",sum(60+38)
  write(*,1001) "SUM IN2O   =",sum(60+39)
  write(*,1001) "SUM IHNO   =",sum(60+40)
  write(*,1001) "SUM ICN    =",sum(60+41)
  write(*,1001) "SUM IHCN   =",sum(60+42)
  write(*,1001) "SUM IH2CN  =",sum(60+43)
  write(*,1001) "SUM IHCNN  =",sum(60+44)
  write(*,1001) "SUM IHCNO  =",sum(60+45)
  write(*,1001) "SUM IHOCN  =",sum(60+46)
  write(*,1001) "SUM IHNCO  =",sum(60+47)
  write(*,1001) "SUM INCO   =",sum(60+48)
  write(*,1001) "SUM IN2    =",sum(60+49)
  write(*,1001) "SUM IAR    =",sum(60+50)
  write(*,1001) "SUM IC3H7  =",sum(60+51)
  write(*,1001) "SUM IC3H8  =",sum(60+52)
  write(*,1001) "SUM ICH2CHO=",sum(60+53)
  write(*,1001) "SUM ICH3CHO=",sum(60+54)

end program gri_analysis
