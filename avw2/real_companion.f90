
program main 
   
  implicit none 
   
  complex(kind(1d0)), allocatable :: wpoly(:),wroots(:),wcoeffs(:,:),allroots(:,:)
  integer, allocatable :: iterations(:) 
  double precision, allocatable :: res(:,:),poly(:),rroots(:),iroots(:),coeffs(:,:)
  double precision :: time
   
  integer :: ii,jj,kk,ll,mm,N,zero,flag
  integer :: clock_start,clock_end,clock_rate 
  
  open (unit=8, file='poly.txt', status='unknown')
  open (unit=9, file='roots.txt', status='unknown')
  open (unit=10, file='errors.txt', status='unknown')


    	N = 10

	allocate(poly(N+1),rroots(N),iroots(N),coeffs(N,3),allroots(N,2),iterations(N),res(N,6))
	allocate(wpoly(N+1),wroots(N),wcoeffs(N,3))

	call dnormalpoly(N,poly(2:N+1))
	poly(1) = 1d0

	do kk=1,N+1
		wpoly(kk) = complex(poly(kk),0d0)
		write(8,*) poly(kk)
	end do

	CALL SYSTEM_CLOCK(COUNT_RATE=clock_rate) 
	CALL SYSTEM_CLOCK(COUNT=clock_start)

	call DAVW2(0,N,1,POLY,COEFFS,RROOTS,IROOTS,ITERATIONS,FLAG)
	
	do kk=1,N
		wroots(kk) = complex(rroots(kk),iroots(kk))
	end do

	call RESCHECK(0,N,1,1,wPOLY(2:N+1)/wpoly(1),wCOEFFS,wROOTS,ALLROOTS,RES)	

	CALL SYSTEM_CLOCK(COUNT=clock_end)  
	time = dble(clock_end - clock_start)/dble(clock_rate) 
	write (*,*) "Run time =", time

	do kk=1,N
		write(10,*) res(kk,:)
		write(9,*) dble(allroots(kk,2)),dimag(allroots(kk,2))
	end do

	deallocate(poly,rroots,iroots,coeffs,allroots,iterations,res)
	deallocate(wpoly,wroots,wcoeffs)

  close(8)
  close(9)
  close(10)


end program main
