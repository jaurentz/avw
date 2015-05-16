! ********************************************** 
!   1/11/2012 
! ********************************************** 
! 
! Main program for computing roots of a polynomial 
! using a factorization of the companion matrix 
! and a modification of Francis' single shift algorithm
! that preserves the factored structure. This method was
! developed by David Watkins and Raf Vandebril.
! 
 
 
program main 
   
  implicit none 
   
  complex(kind(1d0)), allocatable :: poly(:),roots(:)
  double precision, allocatable :: residuals(:,:)
  integer, allocatable :: iterations(:) 
  double precision :: our_time,res_time,preal,pimag
   
  integer :: i,j,degree,zero,flag 
  integer :: clock_start,clock_end,clock_rate 
  
degree = 1200
flag = 0

  allocate(poly(degree),roots(degree),iterations(degree),residuals(degree,6)) 

  call init_random_seed()
  call cnormalpoly(degree,poly)

  CALL SYSTEM_CLOCK(COUNT_RATE=clock_rate) 
  CALL SYSTEM_CLOCK(COUNT=clock_start) 

  call avw(poly,degree,roots,iterations,flag)

  CALL SYSTEM_CLOCK(COUNT=clock_end)  
  our_time = dble(clock_end - clock_start)/dble(clock_rate) 
  write(*,*)"AVW time =",our_time 

  CALL SYSTEM_CLOCK(COUNT_RATE=clock_rate) 
  CALL SYSTEM_CLOCK(COUNT=clock_start) 

  call rescheck(degree,poly,roots,1,residuals)

  CALL SYSTEM_CLOCK(COUNT=clock_end)  
  res_time = dble(clock_end - clock_start)/dble(clock_rate) 
  write(*,*)"RESCHECK time =",res_time 
  write(*,*)"RESCHECK RATIO =",res_time/our_time

  open (unit=7, file='roots.txt', status='unknown')
  do j=1,degree
  	write (7,*) dble(roots(j)),dimag(roots(j))
  end do
  close(7)

  open (unit=8, file='poly.txt', status='unknown')
  do j=1,degree
  	write (8,*) dble(poly(j)),dimag(poly(j))
  end do
  close(8)

  open (unit=9, file='iterations.txt', status='unknown')
  do j=1,degree
  	write (9,*) iterations(j)
  end do
  close(9)

  open (unit=10, file='residuals.txt', status='unknown')
  do j=1,degree
  	write (10,*) residuals(j,:)
  end do
  close(10)

  deallocate(poly,roots,iterations,residuals)

end program main
