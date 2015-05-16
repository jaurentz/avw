! ********************************************** 
!  August 1, 2012
! ********************************************** 
! 
! This subroutine computes three residuals for each computed root, lambda, 
! of a polynomial P(x). It is also capable of applying an arbitrary number
! of Newton iterations to each computed root. 
!
! The residuals are |P(lambda)/P'(lambda)|, |P(lambda)/P'(lambda)/lambda|, and 
! ||Cv-lambda v||/||C||/||v||, in the inifinity norm where C is the Companion Matrix
! and v is the eigenvectro associated with lambda.
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
!
! ROOTS        complex array of length DEGREE, containing the roots
!              of Poly
!
! NEWTNUM      a non-negative integer specifying the number of Newton iterations
!              zero is acceptable.
!
! Output variables:
!
! RESIDUALS    double precision array of dimension (DEGREE,3*(NEWTNUM+1)).
!	       each row of RESIDUALS corresponds to one root of POLY.
!	       columns are in sets of three, with the columns 1, 2 and 3 corresponding
!	       to the three residuals mentioned above respectively. Every set of three
!              columns corresponds to a Newton iteration except for the first one.
!
! ***************************************************************





subroutine rescheck(degree,poly,roots,newtnum,residuals)

	implicit none

	integer, intent(in) :: degree,newtnum
	complex(kind(1d0)), intent(in) :: poly(degree)
	double precision, intent(inout) :: residuals(degree,3*(newtnum+1))
	complex(kind(1d0)), intent(inout) :: roots(degree)

	integer ii,jj,kk
	double precision :: Cnorm
	complex(kind(1d0)) :: f, fprime, lambda

	! Matrix infinity norms
	Cnorm = 0d0
	do ii=1,degree
		Cnorm = Cnorm + abs(poly(ii))
	end do

	Cnorm = max(1d0,Cnorm)

	! Function evaluations and Newton Corrections
	do ii=1,degree

		! Roots inside or on the unit circle
		if(abs(roots(ii)) <= 1d0)then
		do jj=1,(newtnum+1)
			! function evals
			lambda = roots(ii)
			f = lambda + poly(1)
			fprime = complex(dble(degree),0d0)*f - poly(1)
			do kk=2,(degree-1)
				f = lambda*f + poly(kk)
				fprime = lambda*fprime + complex(dble(degree-kk),0d0)*poly(kk)
			end do
			f = f*lambda + poly(degree)
	
			! Store residuals
			residuals(ii,3*(jj-1)+1) = abs(f/fprime)
			residuals(ii,3*(jj-1)+2) = abs(f/fprime/lambda)
			residuals(ii,3*(jj-1)+3) = abs(f)/Cnorm

			! Newton correction
			if((newtnum+1-jj) > 0)then
				lambda = lambda - f/fprime
				roots(ii) = lambda
			end if
		end do
		

		! Roots outside the unit circle
		else
		do jj=1,(newtnum+1)
			! function evals
			lambda = complex(1d0,0d0)/roots(ii)
			f = poly(degree)*lambda + poly(degree-1)
			fprime = complex(dble(degree),0d0)*f - poly(degree-1)
			do kk=2,(degree-1)
				f = lambda*f + poly(degree-kk)
				fprime = lambda*fprime + complex(dble(degree-kk),0d0)*poly(degree-kk)
			end do
			f = f*lambda + complex(1d0,0d0)
	
			! Store residuals
			residuals(ii,3*(jj-1)+1) = abs(f/fprime*roots(ii)*roots(ii))
			residuals(ii,3*(jj-1)+2) = abs(f/fprime*roots(ii))
			residuals(ii,3*(jj-1)+3) = abs(f)*roots(ii)/Cnorm

			! Newton correction
			if((newtnum+1-jj) > 0)then
				lambda = lambda - f/fprime
				roots(ii) = complex(1d0,0d0)/lambda
			end if
		end do
		end if

	end do


end subroutine
