program gri_converge

  character in1*(32)
  character in2*(32)

  integer i,j,nsteps,nx
  real*8 time
  
  real*8 data1(1024,114)
  real*8 data2(1024,114)

  real*8 L0(114)
  real*8 L1(114)
  real*8 L2(114)

  read(*,*) in1
  read(*,*) in2

  open(10,file=in1,form='formatted')
  read(10,*) nsteps
  read(10,*) nx
  read(10,*) time

  do i=0,nx-1
     read(10,*) data1(i,1), &
                data1(i,2:54), &
                data1(i,55), &
                data1(i,56), &
                data1(i,57), &
                data1(i,58), &
                data1(i,59), &
                data1(i,60), &
                data1(i,61), &
                data1(i,62:114)
  end do

  open(20,file=in2,form='formatted')
  read(20,*) nsteps
  read(20,*) nx
  read(20,*) time

  do i=0,nx-1
     read(20,*) data2(i,1), &
                data2(i,2:54), &
                data2(i,55), &
                data2(i,56), &
                data2(i,57), &
                data2(i,58), &
                data2(i,59), &
                data2(i,60), &
                data2(i,61), &
                data2(i,62:114)
  end do

  L0 = 0.d0
  L1 = 0.d0
  L2 = 0.d0

  do i=0,nx-1
     do j=1,114
        L0(j) = max(L0(j), abs(data1(i,j)-data2(i,j)))
        L1(j) = L1(j) + abs(data1(i,j)-data2(i,j))
        L2(j) = L2(j) + (data1(i,j)-data2(i,j))**2
     end do
  end do
  L1 = L1 / dble(nx)
  L2 = sqrt(L2/nx)

  write(*,*) "nsteps =",nsteps
  write(*,*) "nx     =",nx
  write(*,*) "time   =",time
  write(*,*)
  write(*,*) "L0 NORM: H2          =",L0(2)
  write(*,*) "L0 NORM: Density     =",L0(55)
  write(*,*) "L0 NORM: Temperature =",L0(56)
  write(*,*) "L0 NORM: RhoH        =",L0(57)
  write(*,*) "L0 NORM: velocity    =",L0(58)
  write(*,*)
  write(*,*) "L1 NORM: H2          =",L1(2)
  write(*,*) "L1 NORM: Density     =",L1(55)
  write(*,*) "L1 NORM: Temperature =",L1(56)
  write(*,*) "L1 NORM: RhoH        =",L1(57)
  write(*,*) "L1 NORM: velocity    =",L1(58)
  write(*,*)
  write(*,*) "L2 NORM: H2          =",L2(2)
  write(*,*) "L2 NORM: Density     =",L2(55)
  write(*,*) "L2 NORM: Temperature =",L2(56)
  write(*,*) "L2 NORM: RhoH        =",L2(57)
  write(*,*) "L2 NORM: velocity    =",L2(58)

end program gri_converge
