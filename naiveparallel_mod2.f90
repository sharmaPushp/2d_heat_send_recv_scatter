program heat2d
  use mpi
  implicit none
  integer :: nprocs, myrank, ierr
  integer, parameter :: source=0
  integer, parameter :: tag1=99
  integer, parameter :: tag2=199
  integer :: status(mpi_status_size)

  integer,parameter :: nn=512
  integer :: row, col,nx,ny,nxx,nyy,ntimesteps,t
  real :: t0(nn,nn)
  real :: tn(nn,nn)
  real, allocatable :: t0slice(:,:)
  real, allocatable :: tnslice(:,:)
  real :: convergence
  integer :: slice,il,rad

  character(1024):: filename
  nx=nn
  ny=nn
  ntimesteps=1000
  rad=100
  
  call mpi_init(ierr)
  call mpi_comm_size(mpi_comm_world, nprocs, ierr)
  call mpi_comm_rank(mpi_comm_world, myrank, ierr)

  slice = nx/nprocs
  nxx = slice+2
  nyy = ny
  allocate(t0slice(ny,nxx))
  allocate(tnslice(ny,nxx))

  t0=0.0          !!!!!  t0 intialised as zero everywhere
  if (myrank==0) then
    do col = 1,nx
      do row = 1,ny
          if((( row-ny/2 )*( row-ny/2 )+( col-nx/2)*( col - nx/2 ))&
            & .le.( rad * rad )) then
            t0(row,col) = 100.0
          endif
        enddo
      enddo
  endif
  
  call mpi_scatter(t0(1,1), ny*slice, mpi_real, t0slice(1,2), ny*slice, mpi_real, source, mpi_comm_world, ierr)

write(filename,FMT='(A9,I6.6,A4)') "send_recv",myrank,".dat"
filename=trim(filename)
open(unit=12,file=filename) !,status="replace")

  do t=1, ntimesteps

      if (t.eq.1) write(12,*)'proc', myrank

    if (myrank.gt.0) then

      !send to left
      call mpi_send(t0slice(1,2),nyy,mpi_real, myrank-1, tag1, mpi_comm_world, ierr)
      if (t.eq.1) write(12,*)'proc_except0', myrank,'send to', myrank-1
      !receive from left
      call mpi_recv(t0slice(1,1),nyy,mpi_real, myrank-1, tag2, mpi_comm_world,status, ierr)
      if (t.eq.1) write(12,*),'proc_except0', myrank,'receive from', myrank-1

    endif

    if (myrank.lt.nprocs-1) then
      !send to right
      call mpi_send(t0slice(1,nxx-1),nyy,mpi_real, myrank+1, tag2, mpi_comm_world, ierr)
      if (t.eq.1) write(12,*),'proc_except3', myrank,'send to', myrank+1
      !receive from right
      call mpi_recv(t0slice(1,nxx),nyy,mpi_real, myrank+1, tag1, mpi_comm_world,status, ierr)
      if (t.eq.1) write(12,*),'proc_except3', myrank,'receive from', myrank+1

    endif
    if (myrank==0) then
      t0slice(:,1)=t0slice(:,2)
    endif
    if (myrank==nprocs-1) then
      t0slice(:,nxx)=t0slice(:,nxx-1)
    endif

    do col = 2,nxx-1
      do row = 2,nyy-1
        tnslice(row,col)   = 0.25*(t0slice(row,col + 1) + t0slice(row,col-1) +  t0slice(row+1,col) + t0slice(row-1,col))
      enddo
    enddo
    tnslice(1,:)=tnslice(2,:)
    tnslice(nyy,:)=tnslice(nyy-1,:)

    t0slice=tnslice
 

    if (t == 1 .or. mod(t,1000)==0) then
      write(filename,FMT='(A5,I6.6,A1,I3.3,A4)') "temp_",t,"_",myrank,".dat"
      filename=trim(filename)

!      xstart = rank*128 + 1
!      xend = gg
!!!       if(myrank == 2) then
      open(unit=52,file=filename,status="replace")
          write(52,*) 'Variables = "X" "Y" "T"'
!          write(52,*) 'Zone I=,',nn, 'j=',nn/nprocs !'F=Block'
          write(52,*) 'Zone I=,',nn, 'j=',nn/nprocs !'F=Block'
        do col=2,nxx-1
         do row=1,nyy
!          write(52,'(500F13.6)')(t0slice(row,col),col=2,nxx-1)
!          write(52,*)(row, col-1, t0slice(row,col),col=2,nxx-1)
!           write(52,*) row, col-1,t0slice(row,col)
           write(52,*) (myrank*128+col-1), row, t0slice(row,col)
!            write(52,'(500F13.6)')(t0slice(row,col),row=1,nyy)
         enddo
        enddo

      close(52)
      print*,'file written at',t,'timestep'  
    endif
  enddo
close(12)
  call mpi_finalize(ierr)
 
end program heat2d
