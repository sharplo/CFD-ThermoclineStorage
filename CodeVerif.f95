MODULE Code_Verification
	
	USE Input_Output
	USE Semi_Discretized_Model

	PRIVATE

	PUBLIC OrderVerificationFluid
	
	REAL, PARAMETER :: Pi = 4*ATAN(1.), ErrThd = 8E-3, &
		height = 2E3, u_f = 1E-3, alpha_f = 2E3, dt = 3.69E-6
	INTEGER, PARAMETER :: MaxTStep = 1E7, waveNum = 4, pt = 12

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
	
	SUBROUTINE OrderVerificationFluid()

		IMPLICIT NONE
		
		REAL :: k, dx, sigma, d_f, discL1, discL2, discInf, stdyL1, stdyL2, stdyInf
		REAL, ALLOCATABLE :: manSol(:), Temp_f(:), Temp_pre(:), err_arr(:,:), loc_arr(:,:), sol_arr(:,:)
		INTEGER :: nCells, errorFlag, discInfLoc, stdyInfLoc, tStep, i, j
		CHARACTER(LEN=6) :: label(5), label2(2), label3(3)
		LOGICAL :: Rel

		k = 2*Pi*waveNum/height
		IF (waveNum == 4) THEN
			Rel = .FALSE. ! because cos(kx) is 0 at grid centres for nCells = 8
		ELSE
			Rel = .TRUE.
		END IF

		ALLOCATE(err_arr(4,pt), loc_arr(2,pt), STAT=errorFlag)
		IF (errorFlag /= 0) THEN
			WRITE(*,*) "Error: Could not allocate err_arr, loc_arr!"
			STOP
		END IF
		
		i = 0; nCells = 2**3

		DO WHILE (i < pt) ! number of points in graphs for OVS

			ALLOCATE(manSol(nCells), Temp_f(nCells), Temp_pre(nCells), STAT=errorFlag)
			IF (errorFlag /= 0) THEN
				WRITE(*,*) "Error: Could not allocate manSol, Temp_f, Temp_pre!"
				STOP
			END IF
			
			dx = height/nCells; sigma = u_f*dt/dx; d_f = alpha_f*dt/(dx*dx)

			DO j = 1,nCells
				manSol(j) = cos(k*dx*(j-1./2.)) ! manufactured solution T = cos(kx)
				Temp_f(j) = manSol(j) ! manSol taken as initial condition
			END DO
			
			! Evolve Temp_f over time to confirm steadiness
			DO tStep = 1,MaxTStep

				DO j = 1,nCells
					Temp_pre(j) = Temp_f(j) ! for calculating dT/dt
				END DO

				CALL EvolveTempFluid(sigma, d_f, 1., 1., nCells, dx, .TRUE., k, Temp_f)
				CALL CalErrorNorms(nCells, Temp_f, Temp_pre, .FALSE., stdyL1, stdyL2, stdyInf, stdyInfLoc)
				
				! Addressing the progress to user
				IF (mod(tStep,MaxTStep/10) == 0) THEN
					WRITE(*,*) "after step", tStep, "dT/dt =", stdyL2/dt, "stdyInfLoc =", stdyInfLoc
				END IF

				IF (stdyL2/dt < ErrThd) THEN
					
					CALL CalErrorNorms(nCells, Temp_f, manSol, Rel, discL1, discL2, discInf, discInfLoc)
					
					! Write on file only if the approximate solution is steady
					err_arr(1, i+1) = LOG10(dx)
					err_arr(2, i+1) = LOG10(discL1)
					err_arr(3, i+1) = LOG10(discL2)
					err_arr(4, i+1) = LOG10(discInf)
					loc_arr(1, i+1) = LOG10(dx)
					loc_arr(2, i+1) = REAL(discInfLoc)/REAL(nCells)
					
					WRITE(*,"(A, 1X, F3.0, 1X, A, 1X, I8)") "For 2^", LOG(REAL(nCells))/LOG(2.), "cells, &
						manufactured solution converges after step", tStep
					WRITE(*,*) "dT/dt =", stdyL2/dt, "and stdyInfLoc =", stdyInfLoc
					WRITE(*,*) "discL2 =", discL2, "dicsInf =", discInf, "at", discInfLoc
					WRITE(*,"(A)")
			
					EXIT

				ELSE IF (tStep == MaxTStep) THEN
					
					CALL CalErrorNorms(nCells, Temp_f, manSol, Rel, discL1, discL2, discInf, discInfLoc)
					
					WRITE(*,"(A, 1X, F3.0, 1X, A)") &
						"For 2^", LOG(REAL(nCells))/LOG(2.), "cells, manufactured solution CANNOT converge!"
					WRITE(*,*) "dT/dt =", stdyL2/dt, "and stdyInfLoc =", stdyInfLoc
					WRITE(*,*) "discL2 =", discL2, "dicsInf =", discInf, "at", discInfLoc
					WRITE(*,"(A)")

				END IF

			END DO
			
			! Visualize approximate solution and manufactured solution
			IF (nCells == 8) THEN
				
				ALLOCATE(sol_arr(3,nCells), STAT=errorFlag)
				IF (errorFlag /= 0) THEN
					WRITE(*,*) "Error: Could not allocate sol_arr!"
					STOP
				END IF
				
				DO j = 1,nCells
					sol_arr(1,j) = (j-1)
					sol_arr(2,j) = manSol(j)
					sol_arr(3,j) = Temp_f(j)
				END DO
				label3(1) = "x_i"; label3(2) = "manSol"; label3(3) = "appSol"
				CALL PlotFigure(4, "Temp_f.dat", "x_i", "temperature", label3, 3, nCells, sol_arr)
				
				DEALLOCATE(sol_arr, STAT=errorFlag)
				IF (errorFlag /= 0) THEN
					WRITE(*,*) "Error: Could not allocate sol_arr!"
					STOP
				END IF

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
		label2(1) = "dx"; label2(2) = "Inf"
		CALL PlotFigure(2, "discErr.dat", "log10(dx)", "log10(error)", label, 4, pt, err_arr)
		CALL PlotFigure(3, "discLoc.dat", "log10(dx)", "RelPos", label2, 2, pt, loc_arr)
		
	END SUBROUTINE OrderVerificationFluid

END MODULE Code_Verification
