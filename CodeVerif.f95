MODULE Code_Verification
	
	USE Input_Output
	USE Semi_Discretized_Model

	PRIVATE

	PUBLIC CalErrorNorms, OrderVerification
	
	REAL, PARAMETER :: Pi = 4*ATAN(1.), ErrThd = 1E-4
	INTEGER, PARAMETER :: waveNum = 1, MaxTStep = 1E6
	
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
				err_i = abs( (arr_f(i) - arr_i(i) ) / arr_i(i) )
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
	
	SUBROUTINE OrderVerification()

		IMPLICIT NONE
		
		REAL :: height, sigma, d_f, u_f, alpha_f, Pe, Temp_in, k, dx, dt, &
			discL1, discL2, discInf, stdyL1, stdyL2, stdyInf
		INTEGER :: nCells, errorFlag, discInfLoc, stdyInfLoc, pt, tStep, i, j
		REAL, ALLOCATABLE :: manSol(:), Temp_f(:), Temp_pre(:), err_arr(:,:), loc_arr(:,:), sol_arr(:,:)
		CHARACTER(LEN=6) :: label(5), label2(3), label3(2)

		height = 3.2E6; u_f = 1E0; alpha_f = 3.2E9

		k = 2*Pi*waveNum/height; nCells = 2**3
		pt = 10

		ALLOCATE(err_arr(4,pt), loc_arr(2,pt), STAT=errorFlag)
		IF (errorFlag /= 0) THEN
			WRITE(*,*) "Error: Could not allocate err_arr, loc_arr!"
			STOP
		END IF
		
		i = 0

		DO WHILE (i < pt) ! number of points in graphs for OVS
			
			ALLOCATE(manSol(nCells), Temp_f(nCells), Temp_pre(nCells), STAT=errorFlag)
			IF (errorFlag /= 0) THEN
				WRITE(*,*) "Error: Could not allocate manSol, Temp_f, Temp_pre!"
				STOP
			END IF
			
			dx = height/nCells;
			! dt = 0.47*dx*dx/alpha_f ! Pe = 1
			dt = 0.3/(u_f/dx + 2*alpha_f/dx/dx)
			sigma = u_f*dt/dx; d_f = alpha_f*dt/(dx*dx)

			DO j = 1,nCells
				manSol(j) = cos(k*dx*(j-1./2.)) ! manufactured solution T = cos(kx)
				Temp_f(j) = manSol(j) ! manSol taken as initial condition
			END DO
			
			DO tStep = 1,MaxTStep ! evolve manSol over time to confirm steadiness
				
				DO j = 1,nCells
					Temp_pre(j) = Temp_f(j)
				END DO

				CALL EvolveTempFluid(sigma, d_f, 1., nCells, dx, .TRUE., k, Temp_f)
				
				CALL CalErrorNorms(nCells, Temp_f, Temp_pre, .FALSE., stdyL1, stdyL2, stdyInf, stdyInfLoc)
				IF (mod(tStep,MaxTStep/10) == 0) THEN
					WRITE(*,*) "after step", tStep, "dT/dt =", stdyL1/dt, stdyL2/dt, stdyInfLoc
				END IF

				IF (stdyL2/dt < ErrThd) THEN
					
					WRITE(*,"(A, 1X, F3.0, 1X, A, 1X, I7)") "For 2^", LOG(REAL(nCells))/LOG(2.), "cells, &
						manufactured solution converges after step", tStep
					WRITE(*,*) "Temp_f =", Temp_f(1), Temp_f(nCells/2), Temp_f(nCells)
					WRITE(*,*) "dT/dt =", stdyL2/dt, "stdy error =", stdyL1, stdyL2, stdyInf, stdyInfLoc
					
					CALL CalErrorNorms(nCells, Temp_f, manSol, .TRUE., discL1, discL2, discInf, discInfLoc)
					WRITE(*,*) "disc error =", discL1, discL2, discInf, discInfLoc
					err_arr(1, i+1) = LOG10(dx)
					err_arr(2, i+1) = LOG10(discL1)
					err_arr(3, i+1) = LOG10(discL2)
					err_arr(4, i+1) = LOG10(discInf)
					loc_arr(1, i+1) = LOG10(dx)
					loc_arr(2, i+1) = REAL(discInfLoc)/REAL(nCells)
			
					EXIT

				ELSE IF (tStep == MaxTStep) THEN
					
					WRITE(*,"(A, 1X, F3.0, 1X, A, 1X, I7)") "For 2^", LOG(REAL(nCells))/LOG(2.), "cells, &
						manufactured solution CANNOT converge."
					WRITE(*,*) "Temp_f =", Temp_f(1), Temp_f(nCells/2), Temp_f(nCells)
					WRITE(*,*) "dT/dt =", stdyL2/dt, "stdy error =", stdyL1, stdyL2, stdyInf, stdyInfLoc
					
					CALL CalErrorNorms(nCells, Temp_f, manSol, .TRUE., discL1, discL2, discInf, discInfLoc)
					WRITE(*,*) "disc error =", discL1, discL2, discInf, discInfLoc

				END IF

			END DO
			
			IF (nCells == 8) THEN
				ALLOCATE(sol_arr(3,nCells))
				DO j = 1,nCells
					sol_arr(1,j) = (j-1)
					sol_arr(2,j) = manSol(j)
					sol_arr(3,j) = Temp_f(j)
				END DO
				label2(1) = "x_i"; label2(2) = "manSol"; label2(3) = "appSol"
				CALL PlotFigure(4, "Temp_f.dat", "x_i", "temperature", label2, 3, nCells, sol_arr)
				DEALLOCATE(sol_arr)
			END IF

			DEALLOCATE(manSol, Temp_f, Temp_pre, STAT=errorFlag)
			IF (errorFlag /= 0) THEN
				WRITE(*,*) "Error: Could not deallocate manSol, Temp_f, Temp_pre!"
				STOP
			END IF
			
			nCells = nCells*2
			i = i+1

		END DO
	
		label(1) = "dx"; label(2) = "L1"; label(3) = "L2"; label(4) = "Inf"
		label3(1) = "dx"; label3(2) = "Inf"
		CALL PlotFigure(2, "discErr.dat", "log10(dx)", "log10(error)", label, 4, pt, err_arr)
		CALL PlotFigure(3, "discLoc.dat", "log10(dx)", "relPos", label3, 2, pt, loc_arr)
		
	END SUBROUTINE OrderVerification

END MODULE Code_Verification
