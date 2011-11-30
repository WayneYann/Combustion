program chemh_analysis

  implicit none

  character in1*(32)

  integer i,j,k,nsteps_1,nx_1
  real*8 time_1
  
  real*8 data1  (4096,26)

  real*8 sum(26)

  real*8 L0_13(26)
  real*8 L1_13(26)
  real*8 L2_13(26)
  real*8 L0_23(26)
  real*8 L1_23(26)
  real*8 L2_23(26)

  read(*,*) in1

  open(10,file=in1,form='formatted')
  read(10,*) nsteps_1
  read(10,*) nx_1
  read(10,*) time_1

  do i=0,nx_1-1
     read(10,*) data1(i,1:26)
  end do

  !!!!!!!!!!!!!!!!!!!
  ! error for input1
  !!!!!!!!!!!!!!!!!!!

  sum = 0.d0

  do i=0,nx_1-1
     do j=1,26
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
  write(*,1001) "SUM N2     =",sum(10)
  write(*,1001) "SUM IH2    =",sum(2+16)
  write(*,1001) "SUM IH     =",sum(3+16)
  write(*,1001) "SUM IO     =",sum(4+16)
  write(*,1001) "SUM IO2    =",sum(5+16)
  write(*,1001) "SUM IOH    =",sum(6+16)
  write(*,1001) "SUM IH2O   =",sum(7+16)
  write(*,1001) "SUM IHO2   =",sum(8+16)
  write(*,1001) "SUM IH2O2  =",sum(9+16)
  write(*,1001) "SUM IN2    =",sum(10+16)

end program chemh_analysis
