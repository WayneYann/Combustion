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
     read(10,*) data1(i,1), &
                data1(i,2:10), &
                data1(i,11), &
                data1(i,12), &
                data1(i,13), &
                data1(i,14), &
                data1(i,15), &
                data1(i,16), &
                data1(i,17), &
                data1(i,18:114)
  end do

  open(20,file=in2,form='formatted')
  read(20,*) nsteps_f
  read(20,*) nx_f
  read(20,*) time_f

  do i=0,nx_f-1
     read(20,*) data2_f(i,1), &
                data2_f(i,2:10), &
                data2_f(i,11), &
                data2_f(i,12), &
                data2_f(i,13), &
                data2_f(i,14), &
                data2_f(i,15), &
                data2_f(i,16), &
                data2_f(i,17), &
                data2_f(i,18:114)
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
  write(*,*) "rr     =",rr
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
