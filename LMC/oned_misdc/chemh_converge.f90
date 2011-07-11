program chemh_converge

  character in1*(32)
  character in2*(32)

  integer i,j,nsteps,nx
  real*8 time
  
  real*8 data1(1024,26)
  real*8 data2(1024,26)

  real*8 L0(26)
  real*8 L1(26)
  real*8 L2(26)

  read(*,*) in1
  read(*,*) in2

  open(10,file=in1,form='formatted')
  read(10,*) nsteps
  read(10,*) nx
  read(10,*) time

  do i=0,nx-1
     read(10,*) data1(i,1), &
                data1(i,2:10), &
                data1(i,11), &
                data1(i,12), &
                data1(i,13), &
                data1(i,14), &
                data1(i,15), &
                data1(i,16), &
                data1(i,17), &
                data1(i,18:26)
  end do

  open(20,file=in2,form='formatted')
  read(20,*) nsteps
  read(20,*) nx
  read(20,*) time

  do i=0,nx-1
     read(20,*) data2(i,1), &
                data2(i,2:10), &
                data2(i,11), &
                data2(i,12), &
                data2(i,13), &
                data2(i,14), &
                data2(i,15), &
                data2(i,16), &
                data2(i,17), &
                data2(i,18:26)
  end do

  L0 = 0.d0
  L1 = 0.d0
  L2 = 0.d0

  do i=0,nx-1
     do j=1,26
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
  write(*,*) "L0 NORM: O2          =",L0(21)
  write(*,*) "L0 NORM: Density     =",L0(11)
  write(*,*) "L0 NORM: Temperature =",L0(12)
  write(*,*) "L0 NORM: RhoH        =",L0(13)
  write(*,*) "L0 NORM: velocity    =",L0(14)
  write(*,*)
  write(*,*) "L1 NORM: O2          =",L1(21)
  write(*,*) "L1 NORM: Density     =",L1(11)
  write(*,*) "L1 NORM: Temperature =",L1(12)
  write(*,*) "L1 NORM: RhoH        =",L1(13)
  write(*,*) "L1 NORM: velocity    =",L1(14)
  write(*,*)
  write(*,*) "L2 NORM: O2          =",L2(21)
  write(*,*) "L2 NORM: Density     =",L2(11)
  write(*,*) "L2 NORM: Temperature =",L2(12)
  write(*,*) "L2 NORM: RhoH        =",L2(13)
  write(*,*) "L2 NORM: velocity    =",L2(14)

end program chemh_converge
