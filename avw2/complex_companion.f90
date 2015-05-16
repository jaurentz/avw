
program main 
   
  implicit none 
   
  complex(kind(1d0)), allocatable :: poly(:),roots(:),coeffs(:,:),allroots(:,:)
  integer, allocatable :: iterations(:) 
  double precision, allocatable :: res(:,:)
  double precision :: time
   
  integer :: ii,jj,kk,ll,mm,N,zero,flag
  integer :: clock_start,clock_end,clock_rate 
  
  open (unit=8, file='poly.txt', status='unknown')
  open (unit=9, file='roots.txt', status='unknown')
  open (unit=10, file='errors.txt', status='unknown')


    	N = 10

	allocate(poly(N+1),roots(N),coeffs(N,3),allroots(N,2),iterations(N),res(N,6))

	call cnormalpoly(N,poly(2:N+1))
	poly(1) = complex(1d0,0d0)

	do kk=1,N+1
		write(8,*) dble(poly(kk)),dimag(poly(kk))
	end do

	CALL SYSTEM_CLOCK(COUNT_RATE=clock_rate) 
	CALL SYSTEM_CLOCK(COUNT=clock_start)

	call ZAVW2(0,N,1,POLY,COEFFS,ROOTS,ITERATIONS,FLAG)
	call RESCHECK(0,N,1,1,POLY(2:N+1)/poly(1),COEFFS,ROOTS,ALLROOTS,RES)	

	CALL SYSTEM_CLOCK(COUNT=clock_end)  
	time = dble(clock_end - clock_start)/dble(clock_rate) 
	write (*,*) "Run time =", time

	do kk=1,N
		write(10,*) res(kk,:)
		write(9,*) dble(allroots(kk,2)),dimag(allroots(kk,2))
	end do

   	deallocate(poly,roots,coeffs,allroots,iterations,res)

  close(8)
  close(9)
  close(10)


end program main
