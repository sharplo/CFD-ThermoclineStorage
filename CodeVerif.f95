MODULE Code_Verification
	
	USE Input_Output
	USE Semi_Discretized_Model

	PRIVATE

	PUBLIC OrderVerificationFluid, OrderVerificationSolid, ErrThd, MaxTStep
	
	REAL, PARAMETER :: Pi = 4*ATAN(1.), ErrThd = 8E-3, &
		height = 2E3, u_f = 0, alpha_f = 2E3, alpha_s = 2E3, dt = 5.9E-5
	INTEGER, PARAMETER :: MaxTStep = 1E7, waveNum = 1, pt = 10

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

	FUNCTION SteadyState(Temp, sigma, d, Temp_in, nCells, dx, dt, MMS, k)

		IMPLICIT NONE

		REAL :: SteadyState(nCells)
		REAL, INTENT(IN)  :: Temp(nCells), sigma, d, Temp_in, dx, dt, k
		INTEGER, INTENT(IN) :: nCells
		LOGICAL, INTENT(IN) :: MMS

		REAL :: Temp_pre(nCells), stdyL1, stdyL2, stdyInf
		INTEGER :: tStep, stdyInfLoc, i

		SteadyState(:) = Temp(:)

		DO tStep = 1,MaxTStep

			Temp_pre(:) = SteadyState(:) ! for calculating dT/dt

			CALL EvolveTemp(SteadyState, sigma, d, Temp_in, nCells, dx, MMS, k)
			CALL CalErrorNorms(nCells, SteadyState, Temp_pre, .FALSE., stdyL1, stdyL2, stdyInf, stdyInfLoc)
				
			! Addressing the progress to user
			IF (mod(tStep,MaxTStep/10) == 0) THEN
				WRITE(*,*) "after step", tStep, "dT/dt =", stdyL2/dt, "stdyInfLoc =", stdyInfLoc
			END IF

			IF (stdyL2/dt < ErrThd) THEN		
				WRITE(*,"(A, 1X, F3.0, 1X, A, 1X, I8)") "For 2^", LOG(REAL(nCells))/LOG(2.), "cells, &
					manufactured solution converges after step", tStep
				WRITE(*,*) "dT/dt =", stdyL2/dt, "and stdyInfLoc =", stdyInfLoc
				EXIT
			ELSE IF (tStep == MaxTStep) THEN
				WRITE(*,"(A, 1X, F3.0, 1X, A)") &
					"For 2^", LOG(REAL(nCells))/LOG(2.), "cells, manufactured solution CANNOT converge!"
				WRITE(*,*) "dT/dt =", stdyL2/dt, "and stdyInfLoc =", stdyInfLoc
				STOP
			END IF

		END DO

	END FUNCTION SteadyState

	SUBROUTINE OrderVerificationFluid()

		IMPLICIT NONE
		
		REAL :: k, dx, sigma, d_f, discL1, discL2, discInf
		REAL, ALLOCATABLE :: manSol(:), Temp_f(:), err_arr(:,:), loc_arr(:,:), sol_arr(:,:)
		INTEGER :: nCells, errorFlag, discInfLoc, i, j
		CHARACTER(LEN=6) :: label(5), label2(2), label3(3)

		k = 2*Pi*waveNum/height

		ALLOCATE(err_arr(4,pt), loc_arr(2,pt), STAT=errorFlag)
		IF (errorFlag /= 0) THEN
			WRITE(*,*) "Error: Could not allocate err_arr, loc_arr!"
			STOP
		END IF
		
		i = 0; nCells = 2**3

		DO WHILE (i < pt) ! number of points in graphs for OVS

			ALLOCATE(manSol(nCells), Temp_f(nCells), STAT=errorFlag)
			IF (errorFlag /= 0) THEN
				WRITE(*,*) "Error: Could not allocate manSol, Temp_f!"
				STOP
			END IF
			
			dx = height/nCells; sigma = u_f*dt/dx; d_f = alpha_f*dt/(dx*dx)

			DO j = 1,nCells
				manSol(j) = cos(k*dx*(j-1./2.)) ! manufactured solution T = cos(kx)
			END DO
			
			! Take manSol as initial condition and evolves it over time to until steady
			Temp_f(:) = SteadyState(manSol, sigma, d_f, 1., nCells, dx, dt, .TRUE., k)
			
			CALL CalErrorNorms(nCells, Temp_f, manSol, .TRUE., discL1, discL2, discInf, discInfLoc)		
			WRITE(*,*) "discL2 =", discL2, "dicsInf =", discInf, "at", discInfLoc
			WRITE(*,"(A)")
			
			! Write on file only if the approximate solution is steady
			err_arr(1, i+1) = LOG10(dx)
			err_arr(2, i+1) = LOG10(discL1)
			err_arr(3, i+1) = LOG10(discL2)
			err_arr(4, i+1) = LOG10(discInf)
			loc_arr(1, i+1) = LOG10(dx)
			loc_arr(2, i+1) = REAL(discInfLoc)/REAL(nCells)
			
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

			DEALLOCATE(manSol, Temp_f, STAT=errorFlag)
			IF (errorFlag /= 0) THEN
				WRITE(*,*) "Error: Could not deallocate manSol, Temp_f!"
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

	SUBROUTINE OrderVerificationSolid()

		IMPLICIT NONE
		
		REAL :: k, dx, d_s, discL1, discL2, discInf
		REAL, ALLOCATABLE :: manSol(:), Temp_s(:), err_arr(:,:), loc_arr(:,:), sol_arr(:,:)
		INTEGER :: nCells, errorFlag, discInfLoc, i, j
		CHARACTER(LEN=6) :: label(5), label2(2), label3(3)

		k = 4*Pi*waveNum/height

		ALLOCATE(err_arr(4,pt), loc_arr(2,pt), STAT=errorFlag)
		IF (errorFlag /= 0) THEN
			WRITE(*,*) "Error: Could not allocate err_arr, loc_arr!"
			STOP
		END IF
		
		i = 0; nCells = 2**3

		DO WHILE (i < pt) ! number of points in graphs for OVS

			ALLOCATE(manSol(nCells), Temp_s(nCells), STAT=errorFlag)
			IF (errorFlag /= 0) THEN
				WRITE(*,*) "Error: Could not allocate manSol, Temp_s!"
				STOP
			END IF
			
			dx = height/nCells; d_s = alpha_s*dt/(dx*dx)

			DO j = 1,nCells
				manSol(j) = cos(k*dx*(j-1./2.)) ! manufactured solution T = cos(kx)
			END DO
			
			! Take manSol as initial condition and evolves it over time to until steady
			Temp_s(:) = SteadyState(manSol, 0., d_s, 0., nCells, dx, dt, .TRUE., k)
			
			CALL CalErrorNorms(nCells, Temp_s, manSol, .TRUE., discL1, discL2, discInf, discInfLoc)		
			WRITE(*,*) "discL2 =", discL2, "dicsInf =", discInf, "at", discInfLoc
			WRITE(*,"(A)")
			
			! Write on file only if the approximate solution is steady
			err_arr(1, i+1) = LOG10(dx)
			err_arr(2, i+1) = LOG10(discL1)
			err_arr(3, i+1) = LOG10(discL2)
			err_arr(4, i+1) = LOG10(discInf)
			loc_arr(1, i+1) = LOG10(dx)
			loc_arr(2, i+1) = REAL(discInfLoc)/REAL(nCells)
			
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
					sol_arr(3,j) = Temp_s(j)
				END DO
				label3(1) = "x_i"; label3(2) = "manSol"; label3(3) = "appSol"
				CALL PlotFigure(4, "Temp_s.dat", "x_i", "temperature", label3, 3, nCells, sol_arr)
				
				DEALLOCATE(sol_arr, STAT=errorFlag)
				IF (errorFlag /= 0) THEN
					WRITE(*,*) "Error: Could not allocate sol_arr!"
					STOP
				END IF

			END IF

			DEALLOCATE(manSol, Temp_s, STAT=errorFlag)
			IF (errorFlag /= 0) THEN
				WRITE(*,*) "Error: Could not deallocate manSol, Temp_s!"
				STOP
			END IF
			
			nCells = nCells*2
			i = i+1

		END DO
	
		label(1) = "dx"; label(2) = "L1"; label(3) = "L2"; label(4) = "Inf"
		label2(1) = "dx"; label2(2) = "Inf"
		CALL PlotFigure(2, "discErr.dat", "log10(dx)", "log10(error)", label, 4, pt, err_arr)
		CALL PlotFigure(3, "discLoc.dat", "log10(dx)", "RelPos", label2, 2, pt, loc_arr)
		
	END SUBROUTINE OrderVerificationSolid

END MODULE Code_Verification
