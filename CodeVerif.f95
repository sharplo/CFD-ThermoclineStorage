MODULE Code_Verification
	
	USE Input_Output
	USE Semi_Discretized_Model

	PRIVATE

	PUBLIC CalErrorNorms, OrderVerification

CONTAINS

	SUBROUTINE CalErrorNorms(nCells, arr_f, arr_i, Rel, errL1, errL2, errInf, errInfLoc)
		
		IMPLICIT NONE

		INTEGER, INTENT(IN) :: nCells
		REAL, INTENT(IN) :: arr_f(nCells), arr_i(nCells)
		LOGICAL, INTENT(IN) :: Rel
		REAL, INTENT(OUT) :: errL1, errL2, errInf
		INTEGER, INTENT(OUT) :: errInfLoc

		INTEGER :: i
		REAL :: err_i
		
		errL1 = 0; errL2 = 0; errInf = -1.

		DO i = 1,nCells
			
			IF (Rel .EQV. .FALSE.) THEN
				err_i = abs(arr_f(i) - arr_i(i))
			ELSE IF ( (Rel .EQV. .TRUE.) .AND. (arr_i(i) /= 0.) ) THEN
				err_i = abs( (arr_f(i) - arr_i(i)) / arr_i(i) )
			ELSE
				WRITE(*,*) "Error: devided by 0!"
				STOP
			END IF
			
			errL1 = errL1 + err_i
			errL2 = errL2 + err_i*err_i
			IF (err_i > errInf) THEN
				errInf = err_i
				errInfLoc = i
			END IF

		END DO
		
		errL1 = errL1/nCells
		errL2 = sqrt(errL2/nCells)

	END SUBROUTINE CalErrorNorms
	
	SUBROUTINE OrderVerification(waveNum, Pi, ErrThd, MaxTStep)

		IMPLICIT NONE
		
		REAL, INTENT(IN) :: Pi, ErrThd
		INTEGER, INTENT(IN) :: waveNum, MaxTStep
		
		REAL :: height, sigma, d_f, u_f, Pe, Temp_in, k, dx, dt, &
			discL1, discL2, discInf, stdyL1, stdyL2, stdyInf
		INTEGER :: nCells, errorFlag, discInfLoc, stdyInfLoc, pt, tStep, i, j
		REAL, ALLOCATABLE :: manSol(:), Temp_f(:), Temp_pre(:), log_dx(:), err_string(:,:)

!		height = 16.; d_f = 0.4; dt = 1E-4
		height = 64.; u_f = 1E-2; dt = 1E-5; Pe = 1000.
!		height = 1.; sigma = 0.6; d_f = 0.1; dt = 1E-4
		pt = 15; nCells = 8
		k = 2*Pi*waveNum/height
			
		ALLOCATE(log_dx(pt), err_string(4,pt), STAT=errorFlag)
		IF (errorFlag /= 0) THEN
			WRITE(*,*) "Error: Could not allocate log_dx, err_string!"
			STOP
		END IF
		
		DO i = 1,pt
			log_dx(i) = LOG10(height/nCells) - (i-1)*LOG10(2.)
		END DO

		i = 0

		DO WHILE (i < pt) ! number of points in graphs for OV study
			
			ALLOCATE(manSol(nCells), Temp_f(nCells), Temp_pre(nCells), STAT=errorFlag)
			IF (errorFlag /= 0) THEN
				WRITE(*,*) "Error: Could not allocate manSol, Temp_f, Temp_pre!"
				STOP
			END IF
			
!			dx = height/nCells; sigma = d_f/nCells
			dx = height/nCells; sigma = u_f*dt/dx; d_f = nCells*sigma/Pe
!			dx = height/nCells

			DO j = 1,nCells
				manSol(j) = cos(k*dx*(j-1./2.)) ! manufactured solution T = cos(kx)
			END DO
			
			DO j = 1,nCells
				Temp_f(j) = cos(k*dx*(j-1./2.)) ! manSol taken as initial condition
!				Temp_f(j) = 2./(k*dx) * cos(k*dx*(j-1)) * sin(k*dx/2.) ! cell average of T(x)=cos(kx)
			END DO
			
			DO tStep = 1,MaxTStep ! evolve manSol over time to confirm steadiness
				
				DO j = 1,nCells
					Temp_pre(j) = Temp_f(j)
				END DO

				CALL EvolveTempFluid(sigma, d_f, manSol(1), nCells, dx, .TRUE., k, Temp_f)
				
!				WRITE(*,*) "before", stdyL1, stdyL2, stdyInf, stdyInfLoc
				CALL CalErrorNorms(nCells, Temp_f, Temp_pre, .FALSE., stdyL1, stdyL2, stdyInf, stdyInfLoc)
!				WRITE(*,*) "after", stdyL1, stdyL2, stdyInf, stdyInfLoc

				IF (stdyL2/dt < ErrThd) THEN	
					WRITE(*,"(A, 1X, F3.0, 1X, A, 1X, I5)") "For 2^", LOG(REAL(nCells))/LOG(2.), "cells, &
						manufactured solution converges after step", tStep
					WRITE(*,*) stdyL1, stdyL2, stdyInf, stdyInfLoc, "stdy"
					WRITE(*,*) dt, stdyL2/dt
					EXIT
				ELSE IF (tStep == MaxTStep) THEN
					WRITE(*,"(A, 1X, F3.0, 1X, A, 1X, I5)") "For 2^", LOG(REAL(nCells))/LOG(2.), "cells, &
						manufactured solution CANNOT converge."
				END IF

			END DO
			
			CALL CalErrorNorms(nCells, Temp_f, manSol, .TRUE., discL1, discL2, discInf, discInfLoc)
			WRITE(*,*) discL1, discL2, discInf, discInfLoc, "disc"
			err_string(1, i+1) = LOG10(discL1)
			err_string(2, i+1) = LOG10(discL2)
			err_string(3, i+1) = LOG10(discInf)
			err_string(4, i+1) = REAL(discInfLoc)
			
			DEALLOCATE(manSol, Temp_f, Temp_pre, STAT=errorFlag)
			IF (errorFlag /= 0) THEN
				WRITE(*,*) "Error: Could not deallocate manSol, Temp_f, Temp_pre!"
				STOP
			END IF
			
			nCells = nCells*2
			i = i+1

		END DO
		
		CALL PlotFigure(3, "discL1.dat", "log10(dx)", log_dx, "log10(L1)", err_string(1,:), pt)
		CALL PlotFigure(4, "discL2.dat", "log10(dx)", log_dx, "log10(L2)", err_string(2,:), pt)
		CALL PlotFigure(7, "discInf.dat", "log10(dx)", log_dx, "log10(L_Inf)", err_string(3,:), pt)
		CALL PlotFigure(8, "discInfLoc.dat", "log10(dx)", log_dx, "cell", err_string(4,:), pt)

	END SUBROUTINE OrderVerification

END MODULE Code_Verification
