! ********************************************** 
!   July 31, 2012
! ********************************************** 
! 
! This subroutine performs the algorithm described in:
!
!	Aurentz, J. L., Vandebril, R., Watkins, D. S. (2012). Fast Computation of the Zeros of a Polynomial via
!	Factorization of the Companion Matrix. SIAM Journal on Scientific Computing.
!
!
! **************************************************************
!
! Input variables:
!
! POLY         complex array of length DEGREE, containing the nonleading 
!              coefficients of P (ordered with decreasing degree)
!
! DEGREE       degree of polynomial 
!
! Output variables:
!
! ROOTS        complex array of length DEGREE, containing the roots
!              of Poly
!
! ITERATIONS   integer array of length DEGREE the stores the number of iterations 
!	       per deflation
!
! FLAG         reports the position of any instance where the algorithm 
!              tried to divide by 0
!
! ***************************************************************


subroutine avw(POLY,DEGREE,ROOTS,ITERATIONS,FLAG)

  implicit none 

  integer, intent(in) :: degree
  complex(kind(1d0)), intent(in) :: poly(degree) 
  complex(kind(1d0)), intent(out) :: roots(degree)
  integer, intent(out) :: iterations(degree),flag
  integer :: stop_index,start_index,zero_index,i,j,k,n,m,it_count=0,it_max
  double precision :: tolerance=1d-16,zero_point,diagonals,diag1,diag2,norm,nr,nc,nm,bgst,x1,x2,x3,x4
  complex(kind(1d0)) :: trace,detm,disc,DET,shift1,shift2,rho(1),fctr
  complex(kind(1d0)) :: Q(2,2),Q_star(2,2),Temp(2,2),A(2,2),B(2,2),C(2,2),H(3,3),TempH(2,3),rr(1,2),cc(2,1),Blocks(2,2,degree)

  double precision :: DZNRM2

! Initialize and normalize blocks

  do k=2,degree
	Blocks(1,1,k) = complex(-dble(poly(k-1)),-dimag(poly(k-1)))
	Blocks(1,2,k) = complex(1d0,0d0)
	Blocks(2,1,k) = complex(1d0,0d0)
	Blocks(2,2,k) = complex(0d0,0d0)
  end do

  Blocks(1,2,degree) = complex(-dble(poly(degree)),-dimag(poly(degree)))

  Blocks(1,1,1) = complex(0d0,0d0)
  Blocks(1,2,1) = complex(0d0,0d0)
  Blocks(2,1,1) = complex(0d0,0d0)
  Blocks(2,2,1) = complex(1d0,0d0)

! Main loop

  stop_index = degree
  start_index = 2
  zero_index = 1
  it_max = 20*degree
  it_count = 0

  do m=1,it_max

	! Zero check
	do i=1,stop_index
		zero_point = abs(Blocks(2,1,stop_index+1-i))

		if(i>1)then
			diag2 = abs(Blocks(1,1,stop_index+2-i)*Blocks(2,2,stop_index+1-i))
		else
			diag2 = abs(Blocks(2,2,stop_index+1-i))
		end if

		if(i<stop_index)then
			diag1 = abs(Blocks(1,1,stop_index+1-i)*Blocks(2,2,stop_index-i))
		else
			diag1 = abs(Blocks(1,1,stop_index+1-i))
		end if

		diagonals = diag1 + diag2

		if(zero_point < tolerance*diagonals)then
			zero_index = stop_index+1-i
			Blocks(2,1,zero_index) = complex(0d0,0d0)

			! Deflate down if not at bottom
			if(zero_index<stop_index)then
				Blocks(1,1,zero_index+1) = Blocks(1,1,zero_index+1)*Blocks(2,2,zero_index)
				Blocks(1,2,zero_index+1) = Blocks(1,2,zero_index+1)*Blocks(2,2,zero_index)
				Blocks(2,2,zero_index) = complex(1d0,0d0)
				start_index = zero_index+1
			end if
			
			exit

		end if
	end do

	! zero at bottom
	if(zero_index==stop_index)then
		roots(zero_index) = Blocks(2,2,zero_index)
		
		iterations(degree-zero_index+1) = it_count
		it_count = 0

		if(stop_index<=1)then
			exit
		end if

		! deflate up
		Blocks(1,2,zero_index-1) = Blocks(1,2,zero_index-1)*Blocks(1,1,zero_index)
		Blocks(2,2,zero_index-1) = Blocks(2,2,zero_index-1)*Blocks(1,1,zero_index)
		Blocks(1,1,zero_index) = complex(1d0,0d0)

		stop_index = stop_index-1

	! 2x2 case
	else if(stop_index==start_index)then
		! calculate eigenvalues using modified quadratic
		trace = Blocks(1,1,stop_index) + Blocks(2,2,stop_index)
		detm = Blocks(1,1,stop_index)*Blocks(2,2,stop_index) - Blocks(1,2,stop_index)*Blocks(2,1,stop_index)

		disc = sqrt(trace*trace - complex(4d0,0d0)*detm)

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

		! store eigenvalues
		roots(stop_index) = shift1
		Blocks(2,2,stop_index-1) = shift2

		iterations(degree-stop_index) = 0
		it_count = 0

		stop_index = stop_index-1

	! more than one block
	else
		it_count = it_count + 1

		! shift calc
		trace = Blocks(2,2,stop_index-1)*Blocks(1,1,stop_index) + Blocks(2,2,stop_index)
		detm = Blocks(1,1,stop_index)*Blocks(2,2,stop_index) - Blocks(1,2,stop_index)*Blocks(2,1,stop_index)
		detm = Blocks(2,2,stop_index-1)*detm

		disc = sqrt(trace*trace - complex(4d0,0d0)*detm)

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
		if(m == 1)then
			rho(1) = complex(sqrt(2d0),sqrt(2d0))
		else if(m == 2)then
			rho(1) = complex(sqrt(2d0),-sqrt(2d0))
		else if(m == 3)then
			rho(1) = complex(-sqrt(2d0),sqrt(2d0))
		else if(m == 4)then
			rho(1) = complex(-sqrt(2d0),-sqrt(2d0))
		else if(mod(it_count,15)==0)then			
			call cnormalpoly(1,rho(1))
			rho(1) = rho(1)/complex(abs(rho(1)),0d0)
			write(*,*) "Exceptional shift"
		else if(abs(Blocks(2,2,stop_index)-shift1)<abs(Blocks(2,2,stop_index)-shift2))then
			rho(1) = shift1
		else
			rho(1) = shift2
		end if

		! Build Q and Q_star
		Q(1,1) = Blocks(1,1,start_index)-rho(1)
		Q(2,1) = Blocks(2,1,start_index)
		norm = dznrm2(2,Q(:,1),1)

		if(norm==0)then
			flag = 2
		end if

		Q(1,1) = Q(1,1)/norm
		Q(2,1) = Q(2,1)/norm
		Q(1,2) = -conjg(Q(2,1))
		Q(2,2) = conjg(Q(1,1))

		Q_star(1,1) = conjg(Q(1,1))
		Q_star(2,1) = conjg(Q(1,2))
		Q_star(1,2) = conjg(Q(2,1))
		Q_star(2,2) = conjg(Q(2,2))

		! initialize A,B and C
		A(1,1) = (Q_star(1,1)*Blocks(1,1,start_index) + Q_star(1,2)*Blocks(2,1,start_index))
		A(1,2) = (Q_star(1,1)*Blocks(1,2,start_index) + Q_star(1,2)*Blocks(2,2,start_index))
		A(2,1) = (Q_star(2,1)*Blocks(1,1,start_index) + Q_star(2,2)*Blocks(2,1,start_index))
		A(2,2) = (Q_star(2,1)*Blocks(1,2,start_index) + Q_star(2,2)*Blocks(2,2,start_index))

		B(1,1) = Blocks(1,1,start_index+1)
		B(1,2) = Blocks(1,2,start_index+1)
		B(2,1) = Blocks(2,1,start_index+1)
		B(2,2) = Blocks(2,2,start_index+1)

		C(1,1) = Q(1,1)
		C(1,2) = Q(1,2)
		C(2,1) = Q(2,1)
		C(2,2) = Q(2,2)		

		! Turnover iterations
		do n=(start_index+1),(stop_index)
			

			call zto(A,B,C,flag)
			if(flag>0)then
				write(*,*) "Flag set at outer iteration ",m," and inner iteration ",n
				return
			end if

			! Update Blocks and C
			Blocks(1,1,n-1) = B(1,1)
			Blocks(1,2,n-1) = B(1,2)
			Blocks(2,1,n-1) = B(2,1)
			Blocks(2,2,n-1) = B(2,2)

			Blocks(1,1,n) = C(1,1)
			Blocks(1,2,n) = C(1,2)
			Blocks(2,1,n) = C(2,1)
			Blocks(2,2,n) = C(2,2)

			! Dump the Bulge
			if(n == stop_index)then
				Temp(1,1) = Blocks(1,1,n)*A(1,1) + Blocks(1,2,n)*A(2,1)
				Temp(1,2) = Blocks(1,1,n)*A(1,2) + Blocks(1,2,n)*A(2,2)
				Temp(2,1) = Blocks(2,1,n)*A(1,1) + Blocks(2,2,n)*A(2,1)
				Temp(2,2) = Blocks(2,1,n)*A(1,2) + Blocks(2,2,n)*A(2,2)

				Blocks(1,1,n) = Temp(1,1)
				Blocks(1,2,n) = Temp(1,2)
				Blocks(2,1,n) = Temp(2,1)
				Blocks(2,2,n) = Temp(2,2)

				exit
			end if

			! Update A, B and C
			Temp(1,1) = A(1,1)
			Temp(1,2) = A(1,2)
			Temp(2,1) = A(2,1)
			Temp(2,2) = A(2,2)

			A(1,1) = C(1,1)
			A(1,2) = C(1,2)
			A(2,1) = C(2,1)
			A(2,2) = C(2,2)

			C(1,1) = Temp(1,1)
			C(1,2) = Temp(1,2)
			C(2,1) = Temp(2,1)
			C(2,2) = Temp(2,2)

			B(1,1) = Blocks(1,1,n+1)
			B(1,2) = Blocks(1,2,n+1)
			B(2,1) = Blocks(2,1,n+1)
			B(2,2) = Blocks(2,2,n+1)

		end do

	end if

  end do

  ! set flag for maxing out loop
  if(m == it_max)then
  	flag = 1
  end if

end subroutine


! ************************************************************** 
! 
! This subroutine performs the turnover step in the AVW algorithm.
!
! **************************************************************

subroutine zto(A,B,C,flag)


	implicit none 

	integer, intent(inout) :: flag
	complex(kind(1d0)), intent(inout) :: A(2,2),B(2,2),C(2,2)
	double precision :: norm,nr,nc,nm,bgst,x1,x2,x3,x4
	complex(kind(1d0)) :: Temp(2,2),H(3,3),TempH(2,3),rr(1,2),cc(2,1),fctr,DET,QA(2,2),QB(2,2),QC(2,2)

	double precision :: DZNRM2

			! Build H
			H(1,1) = A(1,1)*C(1,1) + A(1,2)*B(1,1)*C(2,1)
			H(1,2) = A(1,1)*C(1,2) + A(1,2)*B(1,1)*C(2,2)
			H(1,3) = A(1,2)*B(1,2)
			H(2,1) = A(2,1)*C(1,1) + A(2,2)*B(1,1)*C(2,1)
			H(2,2) = A(2,1)*C(1,2) + A(2,2)*B(1,1)*C(2,2)
			H(2,3) = A(2,2)*B(1,2)
			H(3,1) = B(2,1)*C(2,1)
			H(3,2) = B(2,1)*C(2,2)
			H(3,3) = B(2,2)

			! Build first inverse
			Temp(1,1) = H(1,2)*H(3,3) - H(1,3)*H(3,2)
			Temp(1,2) = -H(1,2)*H(2,3) + H(1,3)*H(2,2)
			norm = dznrm2(2,Temp(1,:),1)

			if(norm == 0)then
				flag = 3
				write(*,*) "Flag set in zto!"
				return
			end if

			Temp(1,1) = Temp(1,1)/norm
			Temp(1,2) = Temp(1,2)/norm

			Temp(2,1) = -H(3,1)
			Temp(2,2) = H(2,1)
			norm = dznrm2(2,Temp(2,:),1)

			if(norm == 0)then
				flag = 4
				write(*,*) "Flag set in zto!"
				return
			end if

			Temp(2,1) = Temp(2,1)/norm
			Temp(2,2) = Temp(2,2)/norm

			! invert and store in A
			DET = Temp(1,1)*Temp(2,2) - Temp(1,2)*Temp(2,1)

			if(DET == 0)then
				flag = 5
				write(*,*) "Flag set in zto!"
				return
			end if
			
			A(1,1) = Temp(2,2)/DET
			A(2,2) = Temp(1,1)/DET
			A(2,1) = -Temp(2,1)/DET
			A(1,2) = -Temp(1,2)/DET

			! Update H
			TempH(1,1) = Temp(1,1)*H(2,1) + Temp(1,2)*H(3,1)
			TempH(1,2) = Temp(1,1)*H(2,2) + Temp(1,2)*H(3,2)
			TempH(1,3) = Temp(1,1)*H(2,3) + Temp(1,2)*H(3,3)
			TempH(2,1) = Temp(2,1)*H(2,1) + Temp(2,2)*H(3,1)
			TempH(2,2) = Temp(2,1)*H(2,2) + Temp(2,2)*H(3,2)
			TempH(2,3) = Temp(2,1)*H(2,3) + Temp(2,2)*H(3,3)

			
			H(2,1) = TempH(1,1)
			H(2,2) = TempH(1,2)
			H(2,3) = TempH(1,3)
			H(3,1) = TempH(2,1)
			H(3,2) = TempH(2,2)
			H(3,3) = TempH(2,3)

			! Prepare to split updated H
			x1 = abs(H(1,2))
			x2 = abs(H(1,3))
			x3 = abs(H(2,2))
			x4 = abs(H(2,3))

			bgst = max(x1,x2,x3,x4)

			if(bgst == 0)then
				flag = 6
				write(*,*) "Flag set in zto!"
				return
			end if

			if(x1 == bgst)then
				fctr = sqrt(H(1,2))
				rr(1,1) = H(1,2)/fctr
				rr(1,2) = H(1,3)/fctr
				cc(1,1) = H(1,2)/fctr
				cc(2,1) = H(2,2)/fctr
			else if(x2 == bgst)then
				fctr = sqrt(H(1,3))
				rr(1,1) = H(1,2)/fctr
				rr(1,2) = H(1,3)/fctr
				cc(1,1) = H(1,3)/fctr
				cc(2,1) = H(2,3)/fctr
			else if(x3 == bgst)then
				fctr = sqrt(H(2,2))
				rr(1,1) = H(2,2)/fctr
				rr(1,2) = H(2,3)/fctr
				cc(1,1) = H(1,2)/fctr
				cc(2,1) = H(2,2)/fctr
			else
				fctr = sqrt(H(2,3))
				rr(1,1) = H(2,2)/fctr
				rr(1,2) = H(2,3)/fctr
				cc(1,1) = H(1,3)/fctr
				cc(2,1) = H(2,3)/fctr
			end if

			nr = dznrm2(2,rr(1,:),1)
			nc = dznrm2(2,cc(:,1),1)

			nm = sqrt(nc*nr)

			rr(1,1) = rr(1,1)*nm/nr
			rr(1,2) = rr(1,2)*nm/nr

			cc(1,1) = cc(1,1)*nm/nc
			cc(2,1) = cc(2,1)*nm/nc

			! Split updated H into B and C
			B(1,1) = H(1,1)
			B(2,1) = H(2,1)
			B(1,2) = cc(1,1)
			B(2,2) = cc(2,1)

			C(1,1) = rr(1,1)
			C(1,2) = rr(1,2)
			C(2,1) = H(3,2)
			C(2,2) = H(3,3)


end subroutine


! ************************************************************** 
! 
! This subroutine initializes the random seed using the CPU time.
!
! **************************************************************

subroutine init_random_seed()

	implicit none

        integer :: i, n, clock
        integer, dimension(:), allocatable :: seed
          
        call RANDOM_SEED(size = n)
        allocate(seed(n))
          
        call SYSTEM_CLOCK(COUNT=clock)
          
        seed = clock + 37 * (/ (i - 1, i = 1, n) /)
        call RANDOM_SEED(PUT = seed)
          
	deallocate(seed)
end subroutine


! ************************************************************** 
! 
! This subroutine generates nonleading coefficients of a complex,
! monic polynomial where both the real and imaginary parts are normally
! distributed with 0 mean and a variance of 1.
!
! **************************************************************


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
				poly(i) = complex(cos(2d0*pi*v)*sqrt(-2d0*log(u)),sin(2d0*pi*v)*sqrt(-2d0*log(u)))
				exit
			end if
		end do

	end do


end subroutine

