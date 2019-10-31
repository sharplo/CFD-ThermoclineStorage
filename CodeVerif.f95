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
			
			IF ( (Rel .EQV. .TRUE.) .AND. arr_i(i) /= 0) THEN
				err_i = abs( (arr_f(i) - arr_i(i)) / arr_i(i) )
			ELSE IF (Rel .EQV. .FALSE.) THEN
				err_i = abs( (arr_f(i) - arr_i(i)) )
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
		
		REAL :: height, sigma, d_f, u_f, Temp_in, k, dx, dt, &
			discL1, discL2, discInf, stdyL1, stdyL2, stdyInf
		INTEGER :: nCells, errorFlag, discInfLoc, stdyInfLoc, pt, tStep, i, j
		REAL, ALLOCATABLE :: manSol(:), Temp_f(:), Temp_pre(:), log_dx(:), err_string(:,:)

		height = 1.; d_f = 0.1; dt = 1E-3; Temp_in = 1.
		pt = 11; nCells = 8
		k = 2*Pi*waveNum/height
			
		ALLOCATE(log_dx(pt), err_string(4,pt), STAT=errorFlag)
		IF (errorFlag /= 0) THEN
			WRITE(*,*) "Error: Could not allocate log_dx, err_string!"
			STOP
		END IF
		
		i = 0

		DO WHILE (i < pt) ! number of points in graphs for OV study
			
			ALLOCATE(manSol(nCells), Temp_f(nCells), Temp_pre(nCells), STAT=errorFlag)
			IF (errorFlag /= 0) THEN
				WRITE(*,*) "Error: Could not allocate manSol, Temp_f, Temp_pre!"
				STOP
			END IF
			
			dx = height/nCells
			sigma = d_f/nCells

			DO j = 1,nCells
				manSol(j) = cos(k*dx*(j-1./2.)) ! manufactured solution T = cos(kx)
!				manSol(j) = 2./(k*dx) * cos(k*dx*(j-1)) * sin(k*dx/2.) ! cell average of T(x)=cos(kx)
			END DO
			
			DO j = 1,nCells
				Temp_f(j) = cos(k*dx*(j-1./2.)) ! manSol taken as initial condition
!				Temp_f(j) = 2./(k*dx) * cos(k*dx*(j-1)) * sin(k*dx/2.) ! cell average of T(x)=cos(kx)
			END DO
			
			DO tStep = 1,MaxTStep ! evolve manSol over time to confirm steadiness
				
				DO j = 1,nCells
					Temp_pre(j) = Temp_f(j)
				END DO

				CALL EvolveTempFluid(sigma, d_f, Temp_in, nCells, dx, .TRUE., k, Temp_f)
				IF (tStep > 1) THEN ! otherwise divided by 0
!					WRITE(*,*) "before", stdyL1, stdyL2, stdyInf, stdyInfLoc
					CALL CalErrorNorms(nCells, Temp_f, Temp_pre, .FALSE., stdyL1, stdyL2, stdyInf, stdyInfLoc)
!					WRITE(*,*) "after", stdyL1, stdyL2, stdyInf, stdyInfLoc
					IF (stdyL2/dt < ErrThd) THEN
						WRITE(*,"(A, 1X, I5, 1X, A, 1X, I5)") "For", nCells, "nCells, manufactured solution converges after step", tStep
						WRITE(*,*) "stdy", stdyL1, stdyL2, stdyInf, stdyInfLoc
						WRITE(*,*) dt, stdyL2/dt
						EXIT
					ELSE IF (tStep == MaxTStep) THEN
						WRITE(*,"(A, 1X, I5, 1X, A, 1X, I5)") "For", nCells, "nCells, manufactured solution CANNOT converge."
						WRITE(*,*) "stdy", stdyL1, stdyL2, stdyInf, stdyInfLoc
						WRITE(*,*) dt, stdyL2/dt
					END IF

				END IF
			END DO
			
			CALL CalErrorNorms(nCells, Temp_f, manSol, .TRUE., discL1, discL2, discInf, discInfLoc)
			WRITE(*,*) "disc", discL1, discL2, discInf, discInfLoc
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
		
		DO i = 1,pt
			log_dx(i) = LOG10(dx) - (i-1)*LOG10(2.)
		END DO

		CALL PlotFigure(3, "discL1.dat", "log10(dx)", log_dx, "log10(L1)", err_string(1,:), pt)
		CALL PlotFigure(4, "discL2.dat", "log10(dx)", log_dx, "log10(L2)", err_string(2,:), pt)
		CALL PlotFigure(7, "discInf.dat", "log10(dx)", log_dx, "log10(L_Inf)", err_string(3,:), pt)
		CALL PlotFigure(8, "discInfLoc.dat", "log10(dx)", log_dx, "cell", err_string(4,:), pt)

	END SUBROUTINE OrderVerification

END MODULE Code_Verification
