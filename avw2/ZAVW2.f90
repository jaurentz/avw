! ********************************************** 
!  Last updated: December 10, 2012
! ********************************************** 
! 
! This subroutine computes the roots of a complex polynomial of degree N
! in any recursively defined polynomial basis. The general 
! algorithm is outlined in: 
!
! Fast Computation of Eigenvalues of Companion, Comrade and Related Matrices, BIT Numerical Mathematics 2012
!
! **************************************************************
!
! Input variables:
!
! BANDSWITCH   0 indicates no COEFFS are input, 1 indicates COEFFS are given
!
! N            degree of polynomial
!
! K            maximum length of the recursion
!
! POLY         array of dimension N+1 containing the coefficients 
!              of the polynomial in terms of the basis elements
!
! COEFFS       array of dimension (N x K) containing the recursion
!	       coefficients of the basis elements, ignored if BANDSWITCH = 0 
!
! Output variables:
!
! ROOTS         contains the computed roots, dimension N
!
! ITERATIONS   integer array storing the number of iterations per root calculation
!
! FLAG         reports the position of any instance where the algorithm 
!              tried to divide by 0
!
! ***************************************************************


subroutine ZAVW2(BANDSWITCH,N,K,POLY,COEFFS,ROOTS,ITERATIONS,FLAG)


   
	implicit none 

	integer, intent(in) :: BANDSWITCH,N,K
	integer, intent(inout) :: ITERATIONS(N),FLAG
	complex(kind(1d0)), intent(in) :: POLY(N+1),COEFFS(N,K) 
	complex(kind(1d0)), intent(inout) :: ROOTS(N)

	FLAG = 0

	if(BANDSWITCH == 0)then
		
		call ZLR(N,POLY,ROOTS,ITERATIONS,FLAG)

	else if(BANDSWITCH == 1)then

		call ZBLR(N,K,POLY,COEFFS,ROOTS,ITERATIONS,FLAG)

	else

		write(*,*) "Not a valid argument for BANDSWITCH"
		write(*,*) ""
		return
	end if

end subroutine


! ************************************************************** 
! 
! This subroutine calculates the roots using the companion matrix
!
! **************************************************************

subroutine ZLR(N,POLY,ROOTS,ITERATIONS,FLAG)

   
  implicit none 

  integer, intent(in) :: N
  complex(kind(1d0)), intent(in) :: poly(N+1) 
  complex(kind(1d0)), intent(inout) :: roots(N)
  integer, intent(out) :: iterations(N),flag
  integer :: stop_index,start_index,zero_index,ii,jj,kk,mm,it_count=0,it_max
  real(kind(1d0)) :: tolerance=1d-16,zero_point,diagonals,diag1,diag2,norm
  complex(kind(1d0)) :: trace,detm,disc,DET,shift1,shift2,rho(1)
  complex(kind(1d0)) :: Blocks(2,2,N),L

! Initialize and normalize blocks

	if(abs(poly(1)) == 0)then
		write(*,*) "Polynomial is of degree less than N!"
		flag = 1
		return
	end if

  do kk=2,N
	Blocks(1,1,kk) = -poly(kk)/poly(1)
	Blocks(1,2,kk) = complex(1d0,0d0)
	Blocks(2,1,kk) = complex(1d0,0d0)
	Blocks(2,2,kk) = complex(0d0,0d0)
  end do

  Blocks(1,2,N) = -poly(N+1)/poly(1)

  Blocks(1,1,1) = complex(0d0,0d0)
  Blocks(1,2,1) = complex(0d0,0d0)
  Blocks(2,1,1) = complex(0d0,0d0)
  Blocks(2,2,1) = complex(1d0,0d0)
	
! Main loop

  stop_index = N
  start_index = 2
  zero_index = 1
  it_max = 20*N
  it_count = 0
  flag = 0

  do mm=1,it_max

	! Zero check
	do ii=1,stop_index
		zero_point = abs(Blocks(2,1,stop_index+1-ii))

		if(ii>1)then
			diag2 = abs(Blocks(1,1,stop_index+2-ii)*Blocks(2,2,stop_index+1-ii))
		else
			diag2 = abs(Blocks(2,2,stop_index+1-ii))
		end if

		if(ii<stop_index)then
			diag1 = abs(Blocks(1,1,stop_index+1-ii)*Blocks(2,2,stop_index-ii))
		else
			diag1 = abs(Blocks(1,1,stop_index+1-ii))
		end if

		diagonals = diag1 + diag2

		if(zero_point < tolerance*diagonals)then
			zero_index = stop_index+1-ii
			Blocks(2,1,zero_index) = complex(0d0,0d0)

			! Deflate down if not at bottom
			if(zero_index<stop_index)then
				Blocks(1,1,zero_index+1) = Blocks(1,1,zero_index+1)*Blocks(2,2,zero_index)
				Blocks(2,1,zero_index+1) = Blocks(2,1,zero_index+1)*Blocks(2,2,zero_index)
				Blocks(2,2,zero_index) = complex(1d0,0d0)
				start_index = zero_index+1
			end if
			
			exit

		end if
	end do

	! zero at bottom
	if(zero_index==stop_index)then
		roots(zero_index) = Blocks(2,2,zero_index)
		
		iterations(N-zero_index+1) = it_count
		it_count = 0

		if(stop_index<=1)then
			exit
		end if

		! deflate up
		Blocks(2,1,zero_index-1) = Blocks(2,1,zero_index-1)*Blocks(1,1,zero_index)
		Blocks(2,2,zero_index-1) = Blocks(2,2,zero_index-1)*Blocks(1,1,zero_index)
		Blocks(1,1,zero_index) = complex(1d0,0d0)

		stop_index = stop_index-1

	! 2x2 case
	else if(stop_index==start_index)then
		! calculate eigenvalues using modified quadratic
		trace = Blocks(1,1,stop_index) + Blocks(2,2,stop_index)
		detm = Blocks(1,1,stop_index)*Blocks(2,2,stop_index) - Blocks(1,2,stop_index)*Blocks(2,1,stop_index)

		disc = sqrt(trace*trace - 4*detm)

		if(abs(trace-disc)<abs(trace+disc))then
			if(abs(trace+disc)==0)then
				shift1 = complex(0d0,0d0)
				shift2 = complex(0d0,0d0)
			else
				shift1 = (trace+disc)/2
				shift2 = 2*detm/(trace+disc)
			end if
		else
			if(abs(trace-disc)==0)then
				shift1 = complex(0d0,0d0)
				shift2 = complex(0d0,0d0)
			else
				shift1 = (trace-disc)/2
				shift2 = 2*detm/(trace-disc)
			end if
		end if

		! store eigenvalues
		roots(stop_index) = shift1
		Blocks(2,2,stop_index-1) = shift2

		iterations(N-stop_index) = 0
		it_count = 0

		stop_index = stop_index-1

	! more than one block
	else
		it_count = it_count + 1

		! shift calc
		trace = Blocks(2,2,stop_index-1)*Blocks(1,1,stop_index) + Blocks(2,2,stop_index)
		detm = Blocks(1,1,stop_index)*Blocks(2,2,stop_index) - Blocks(1,2,stop_index)*Blocks(2,1,stop_index)
		detm = Blocks(2,2,stop_index-1)*detm

		disc = sqrt(trace*trace - 4*detm)

		if(abs(trace-disc)<abs(trace+disc))then
			if(abs(trace+disc)==0)then
				shift1 = complex(0d0,0d0)
				shift2 = complex(0d0,0d0)
			else
				shift1 = (trace+disc)/complex(2d0,0d0)
				shift2 = complex(2d0,0d0)*detm/(trace+disc)
			end if
		else
			if(abs(trace-disc)==0)then
				shift1 = complex(0d0,0d0)
				shift2 = complex(0d0,0d0)
			else
				shift1 = (trace-disc)/complex(2d0,0d0)
				shift2 = complex(2d0,0d0)*detm/(trace-disc)
			end if
		end if

		! pick best shift
		if(mm == 0)then
			rho(1) = complex(dsqrt(2d0),dsqrt(2d0))
!		else if(mm == 2)then
!			rho(1) = complex(dsqrt(2d0),-dsqrt(2d0))
!		else if(mm == 3)then
!			rho(1) = complex(-dsqrt(2d0),dsqrt(2d0))
!		else if(mm == 4)then
!			rho(1) = complex(-dsqrt(2d0),-dsqrt(2d0))
		else if(mod(it_count,15) == 0)then
			call cnormalpoly(1,rho(1))
			rho(1) = rho(1)
			write(*,*) "Exceptional shift in LR!"
		else if(abs(Blocks(2,2,stop_index)-shift1)<abs(Blocks(2,2,stop_index)-shift2))then
			rho(1) = shift1
		else
			rho(1) = shift2
		end if

		! Build L
		if(abs(Blocks(1,1,start_index)-rho(1)) == 0)then
			if(abs(rho(1)) == 0)then
				call cnormalpoly(1,rho(1))
				rho(1) = rho(1)/abs(rho(1))
			else
				rho(1) = -rho(1)
			end if
		end if

		L = Blocks(2,1,start_index)/(Blocks(1,1,start_index)-rho(1))

		! Initial transform
		Blocks(2,:,start_index) = Blocks(2,:,start_index) - L*Blocks(1,:,start_index)	

		! Turnover iterations
		do ii=(start_index+1),(stop_index)
			

			call ZGATO(Blocks(:,:,ii-1),Blocks(:,:,ii),L,flag)
			if(flag>0)then
				write(*,*) "Flag set at outer iteration ",mm," and inner iteration ",ii
				return
			end if

			! Dump the Bulge
			if(ii == stop_index)then

				Blocks(:,1,ii) = Blocks(:,1,ii) + L*Blocks(:,2,ii)

				exit
			end if

		end do

	end if

  end do

  ! set flag for maxing out loop
  if(mm == it_max)then
  	flag = 1
  end if

end subroutine


! ************************************************************** 
! 
! This subroutine calculates the roots using any congenial matrix
!
! **************************************************************

subroutine ZBLR(N,K,POLY,COEFFS,ROOTS,ITERATIONS,FLAG)

	implicit none 

	integer, intent(in) :: N,K
	complex(kind(1d0)), intent(in) :: POLY(N+1),COEFFS(N,K)
	complex(kind(1d0)), intent(inout) :: ROOTS(N)
	integer, intent(out) :: iterations(N),flag

	integer :: stop_index,start_index,zero_index,ii,jj,hh,mm,it_count=0,it_max
	integer :: length
	real(kind(1d0)) :: tolerance=1d-16,zero_point,diagonals,diag1,diag2
	complex(kind(1d0)) :: trace,detm,disc,DET,shift1,shift2,rho(1)
	complex(kind(1d0)) :: L,W(2,2),bulge
	complex(kind(1d0)), allocatable :: spike(:),Blocks(:,:,:),Bands(:,:) 


! Initialize Bands and Blocks
	
	allocate(spike(N),Blocks(2,2,N),Bands(N,K))

	if(abs(poly(1)) == 0)then
		write(*,*) "Polynomial is of degree less than N!"
		flag = 1
		return
	end if

	spike = -Poly(2:(N+1))/Poly(1)*Coeffs(1,1)
	spike(1:(K-1)) = spike(1:(K-1)) + Coeffs(1,2:K)

	Bands = complex(0d0,0d0)
	Bands(1:(N-1),1:K) = Coeffs(2:N,1:K)

	Blocks(:,:,1) = complex(0d0,0d0)
	Blocks(2,2,1) = complex(1d0,0d0)

	do ii=1,(N-1)
		if(abs(Bands(ii,1)) == 0)then
			write(*,*) "Matrix is not properly upper-hessenberg!"
			flag = 2
			return
		end if

		L = spike(ii)/Bands(ii,1)

		length = min(N+1-ii,K)
		spike(ii:(ii+length-1)) = spike(ii:(ii+length-1)) - L*Bands(ii,1:length)

		Blocks(1,1,ii+1) = L
		Blocks(1,2,ii+1) = complex(1d0,0d0)	
		Blocks(2,1,ii+1) = complex(1d0,0d0)
		Blocks(2,2,ii+1) = complex(0d0,0d0)
	end do

	Bands(N,1) = spike(N)



! Main loop

  stop_index = N
  start_index = 2
  zero_index = 1
  it_max = 20*N
  it_count = 0
  flag = 0

  do mm=1,it_max

	! Zero check
	do ii=1,stop_index
		zero_point = abs(Blocks(2,1,stop_index+1-ii))

		if(ii>1)then
			diag2 = abs(Blocks(1,1,stop_index+2-ii)*Blocks(2,2,stop_index+1-ii)*Bands(stop_index+1-ii,1))
		else
			diag2 = abs(Blocks(2,2,stop_index+1-ii)*Bands(stop_index+1-ii,1))
		end if

		if(ii<stop_index)then
			diag1 = abs(Blocks(1,1,stop_index+1-ii)*Blocks(2,2,stop_index-ii)*Bands(stop_index-ii,1))
		else
			diag1 = abs(Blocks(1,1,stop_index+1-ii)*Bands(stop_index-ii,1))
		end if

		diagonals = diag1 + diag2

		if(zero_point < tolerance*diagonals)then
			zero_index = stop_index+1-ii
			Blocks(2,1,zero_index) = complex(0d0,0d0)

			! Deflate down if not at bottom
			if(zero_index<stop_index)then
				Blocks(1,1,zero_index+1) = Blocks(1,1,zero_index+1)*Blocks(2,2,zero_index)
				Blocks(1,2,zero_index+1) = Blocks(1,2,zero_index+1)*Blocks(2,2,zero_index)
				Blocks(1,2,zero_index) = Blocks(1,2,zero_index)/Blocks(2,2,zero_index)
				Blocks(2,2,zero_index) = complex(1d0,0d0)
				start_index = zero_index+1
			end if

			exit

		end if
	end do

	! zero at bottom
	if(zero_index==stop_index)then
		ROOTS(zero_index) = Blocks(2,2,zero_index)*Bands(zero_index,1)
		
		iterations(N-zero_index+1) = it_count
		it_count = 0

		if(stop_index<=1)then
			exit
		end if

		! deflate up
		Blocks(:,2,zero_index-1) = Blocks(1,1,zero_index)*Blocks(:,2,zero_index-1)
		Blocks(1,1,zero_index) = complex(1d0,0d0)

		stop_index = stop_index-1

	! 2x2 case
	else if(stop_index==start_index)then

		! calculate eigenvalues using modified quadratic
		if(K == 1)then
		W(1,1) = Blocks(1,1,stop_index)*Bands(stop_index-1,1)
		W(1,2) = Blocks(1,2,stop_index)*Bands(stop_index,1)
		W(2,1) = Blocks(2,1,stop_index)*Bands(stop_index-1,1)
		W(2,2) = Blocks(2,2,stop_index)*Bands(stop_index,1)
		else
		W(1,1) = Blocks(1,1,stop_index)*Bands(stop_index-1,1)
		W(1,2) = Blocks(1,1,stop_index)*Bands(stop_index-1,2) + Blocks(1,2,stop_index)*Bands(stop_index,1)
		W(2,1) = Blocks(2,1,stop_index)*Bands(stop_index-1,1)
		W(2,2) = Blocks(2,1,stop_index)*Bands(stop_index-1,2) + Blocks(2,2,stop_index)*Bands(stop_index,1)
		end if

		trace = W(1,1) + W(2,2)
		detm = W(1,1)*W(2,2) - W(1,2)*W(2,1)

		disc = sqrt(trace*trace - 4*detm)

		if(abs(trace-disc)<abs(trace+disc))then
			if(abs(trace+disc)==0)then
				shift1 = complex(0d0,0d0)
				shift2 = complex(0d0,0d0)
			else
				shift1 = (trace+disc)/2
				shift2 = 2*detm/(trace+disc)
			end if
		else
			if(abs(trace-disc)==0)then
				shift1 = complex(0d0,0d0)
				shift2 = complex(0d0,0d0)
			else
				shift1 = (trace-disc)/2
				shift2 = 2*detm/(trace-disc)
			end if
		end if

		! store eigenvalues
		ROOTS(stop_index) = shift1
		ROOTS(stop_index-1) = shift2

		iterations(N-stop_index) = 0
		it_count = 0

		if(stop_index <= 2)then
			exit
		end if

		! deflate up
		Blocks(:,2,stop_index-2) = Blocks(1,1,stop_index-1)*Blocks(:,2,stop_index-2)
		Blocks(1,1,stop_index-1) = complex(1d0,0d0)

		stop_index = stop_index-2

	! more than one block
	else

		it_count = it_count + 1

		! shift calc
		if(K == 1)then
		W(1,1) = Blocks(1,1,stop_index)*Bands(stop_index-1,1)*Blocks(2,2,stop_index-1)
		W(1,2) = Blocks(1,2,stop_index)*Bands(stop_index,1)
		W(1,2) = W(1,2)*Blocks(2,2,stop_index-1)
		W(2,1) = Blocks(2,1,stop_index)*Bands(stop_index-1,1)
		W(2,2) = Blocks(2,2,stop_index)*Bands(stop_index,1)
		else
		W(1,1) = Blocks(1,1,stop_index)*Bands(stop_index-1,1)*Blocks(2,2,stop_index-1)
		W(1,2) = Blocks(1,1,stop_index)*Bands(stop_index-1,2) + Blocks(1,2,stop_index)*Bands(stop_index,1)
		W(1,2) = W(1,2)*Blocks(2,2,stop_index-1)
		W(2,1) = Blocks(2,1,stop_index)*Bands(stop_index-1,1)
		W(2,2) = Blocks(2,1,stop_index)*Bands(stop_index-1,2) + Blocks(2,2,stop_index)*Bands(stop_index,1)
		end if


		trace = W(1,1) + W(2,2)
		detm = W(1,1)*W(2,2) - W(1,2)*W(2,1)

		disc = sqrt(trace*trace - 4*detm)

		if(abs(trace-disc)<abs(trace+disc))then
			if(abs(trace+disc)==0)then
				shift1 = complex(0d0,0d0)
				shift2 = complex(0d0,0d0)
			else
				shift1 = (trace+disc)/complex(2d0,0d0)
				shift2 = complex(2d0,0d0)*detm/(trace+disc)
			end if
		else
			if(abs(trace-disc)==0)then
				shift1 = complex(0d0,0d0)
				shift2 = complex(0d0,0d0)
			else
				shift1 = (trace-disc)/complex(2d0,0d0)
				shift2 = complex(2d0,0d0)*detm/(trace-disc)
			end if
		end if

		if(mm == 0)then
			rho(1) = complex(dsqrt(2d0),dsqrt(2d0))
		else if(mod(it_count,15) == 0)then
			call cnormalpoly(1,rho(1))
			rho(1) = rho(1)
			write(*,*) "Exceptional shift in BLR!"
		else if(abs(Blocks(2,2,stop_index)*Bands(stop_index,1)-shift1)<abs(Blocks(2,2,stop_index)*Bands(stop_index,1)-shift2))then
			rho(1) = shift1
		else
			rho(1) = shift2
		end if

		! Build L
		if(abs(Blocks(1,1,start_index)*Bands(start_index-1,1)-rho(1)) == 0)then
			if(abs(rho(1)) == 0)then
				call cnormalpoly(1,rho(1))
				rho(1) = rho(1)/abs(rho(1))
			else
				rho(1) = -rho(1)
			end if
		end if

		L = Blocks(2,1,start_index)/(Blocks(1,1,start_index)*Bands(start_index-1,1)-rho(1))
		L = L*Bands(start_index-1,1)

		! Initial transform
		Blocks(2,:,start_index) = Blocks(2,:,start_index) - L*Blocks(1,:,start_index)	

		! Turnover iterations
		do hh=(start_index),(stop_index-1)
			! Pass through Bands
			length = min(hh-start_index+1,K-1)			
			do jj=1,length
				Bands(hh-jj,jj) = Bands(hh-jj,jj) + L*Bands(hh-jj,jj+1)
			end do
			
			bulge = L*Bands(hh,1)

			L = bulge/Bands(hh-1,1)

			length = min(K-1,stop_index+1-hh)
			do jj=1,length
				Bands(hh,jj) = Bands(hh,jj) - L*Bands(hh-1,jj+1)
			end do

			! Pass through Blocks
			call ZGATO(Blocks(:,:,hh),Blocks(:,:,hh+1),L,flag)
			if(flag>0)then
				write(*,*) "Flag set at outer iteration ",mm," and inner iteration ",hh
				return
			end if
		end do

		! Dump the Bulge			
			! Pass through Bands
			length = min(stop_index-start_index+1,K-1)			
			do jj=1,length
				Bands(stop_index-jj,jj) = Bands(stop_index-jj,jj) + L*Bands(stop_index-jj,jj+1)
			end do
			
			bulge = L*Bands(stop_index,1)

			L = bulge/Bands(stop_index-1,1)

			if(K == 1)then
				Bands(stop_index,1) = Bands(stop_index,1)
			else
				Bands(stop_index,1) = Bands(stop_index,1) - L*Bands(stop_index-1,2)
			end if

		Blocks(:,1,stop_index) = Blocks(:,1,stop_index) + L*Blocks(:,2,stop_index)

	end if

  end do

  ! set flag for maxing out loop
  if(mm == it_max)then
  	flag = 1
  end if

	! free memory
	deallocate(spike,Blocks,Bands)

end subroutine

! ****************************************************************************
!
! This subroutine passes a Gauss transform through two adjacent core 
! transformations.
!
! ****************************************************************************

subroutine ZGATO(A,B,L,flag)


	implicit none 

	integer, intent(inout) :: flag
	complex(kind(1d0)), intent(inout) :: A(2,2),B(2,2),L

	! Update A
	A(:,1) = A(:,1) + L*B(1,1)*A(:,2)

	! pivot decision
	if(abs(A(2,1)) == 0)then
		flag = 1
		return
	end if

	L = L*B(2,1)/A(2,1)
	
	! Update B
	B(2,:) = B(2,:) - L*A(2,2)*B(1,:)
	
end subroutine

! ****************************************************************************
!
! This subroutine generates a one dimensional complex array whose entries are 
! normally distributed with mean 0 and variance 1 in both the real and imaginary parts
!
! ****************************************************************************

subroutine cnormalpoly(degree,poly)

	implicit none
	
	integer, intent(in) :: degree
	complex(kind(1d0)), intent(inout) :: poly(degree) 
	
	double precision :: u,v,s,pi = 3.141592653589793239d0
	integer :: i,j

	do i=1,degree
		do j=1,20
			
			call random_number(u)
			call random_number(v)
	
			s = u**2 + v**2
	
			if(s > 0 .and. s < 1)then				
				poly(i) = complex(cos(2.d0*pi*v)*dsqrt(-2.d0*log(u)),sin(2.d0*pi*v)*dsqrt(-2.d0*log(u)))
				exit
			end if
		end do

	end do


end subroutine

! ****************************************************************************
!
! This subroutine initializes the random seed using the cpu time
!
! ****************************************************************************

SUBROUTINE init_random_seed()

	implicit none

        INTEGER :: i, n, clock
        INTEGER, DIMENSION(:), ALLOCATABLE :: seed
          
        CALL RANDOM_SEED(size = n)
        ALLOCATE(seed(n))
          
        CALL SYSTEM_CLOCK(COUNT=clock)
          
        seed = clock + 37 * (/ (i - 1, i = 1, n) /)
        CALL RANDOM_SEED(PUT = seed)
          
	DEALLOCATE(seed)

END SUBROUTINE
