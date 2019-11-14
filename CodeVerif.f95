MODULE Code_Verification
	
	USE Input_Output
	USE Semi_Discretized_Model

	PRIVATE

	PUBLIC OrderVerification
	
	REAL, PARAMETER :: Pi = 4*ATAN(1.), ErrThd = 8E-3, &
		height = 2E3, u_f = 1E-3, alpha_f = 2E3, alpha_s = 2E3, dt = 3.69E-6
	INTEGER, PARAMETER :: MaxTStep = 1E7, waveNum = 1, pt = 12

CONTAINS

	FUNCTION ErrorNorms(nCells, arr_f, arr_i, Rel)
		
		IMPLICIT NONE

		REAL :: ErrorNorms(4)
		INTEGER, INTENT(IN) :: nCells
		REAL, INTENT(IN) :: arr_f(nCells), arr_i(nCells)
		LOGICAL, INTENT(IN) :: Rel

		REAL :: errL1, errL2, errInf, err_i
		INTEGER :: i, errInfLoc

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

		ErrorNorms(1) = errL1; ErrorNorms(2) = errL2;
		ErrorNorms(3) = errInf; ErrorNorms(4) = REAL(errInfLoc)

	END FUNCTION ErrorNorms

	FUNCTION SteadyState(Temp, sigma, d, Temp_in, nCells, dx, dt, MMS, k)

		IMPLICIT NONE

		REAL :: SteadyState(nCells)
		REAL, INTENT(IN)  :: Temp(nCells), sigma, d, Temp_in, dx, dt, k
		INTEGER, INTENT(IN) :: nCells
		LOGICAL, INTENT(IN) :: MMS

		REAL :: Temp_pre(nCells), stdyErr(4)
		INTEGER :: tStep, i

		SteadyState(:) = Temp(:)

		DO tStep = 1,MaxTStep

			Temp_pre(:) = SteadyState(:) ! for calculating dT/dt

			CALL EvolveTemp(SteadyState, sigma, d, Temp_in, nCells, dx, MMS, k)
			stdyErr(:) = ErrorNorms(nCells, SteadyState, Temp_pre, .FALSE.)
				
			! Address the progress to user
			IF (mod(tStep,MaxTStep/10) == 0) THEN
				WRITE(*,*) "after step", tStep, "dT/dt =", stdyErr(2)/dt, "stdyInfLoc =", INT(stdyErr(4))
			END IF

			IF (stdyErr(2)/dt < ErrThd) THEN ! basically dT/dt
				WRITE(*,*) "Manufactured solution converges after step", tStep
				WRITE(*,*) "dT/dt =", stdyErr(2)/dt, "and stdyInfLoc =", INT(stdyErr(4))
				EXIT
			ELSE IF (tStep == MaxTStep) THEN
				WRITE(*,*) "Manufactured solution CANNOT converge!"
				WRITE(*,*) "dT/dt =", stdyErr(2)/dt, "and stdyInfLoc =", INT(stdyErr(4))
				STOP
			END IF

		END DO

	END FUNCTION SteadyState

	SUBROUTINE OrderVerification()

		IMPLICIT NONE
		
		REAL :: k_f, k_s, dx, sigma, d_f, d_s, discErr_f(4), discErr_s(4)
		REAL, ALLOCATABLE :: manSol_f(:), manSol_s(:), Temp_f(:), Temp_s(:), &
			errArr_f(:,:), errArr_s(:,:), locArr_f(:,:), locArr_s(:,:), solArr(:,:)
		INTEGER :: nCells, errorFlag, i, j
		CHARACTER(LEN=6) :: label(3), label2(4), label3(2)

		k_f = 2*Pi*waveNum/height; k_s = 2*k_f

		ALLOCATE(errArr_f(4,pt), errArr_s(4,pt), locArr_f(2,pt), locArr_s(2,pt), STAT=errorFlag)
		IF (errorFlag /= 0) THEN
			WRITE(*,*) "Error: Could not allocate errArr, locArr!"
			STOP
		END IF
		
		i = 0; nCells = 2**3

		DO WHILE (i < pt) ! number of points in graphs for OVS

			ALLOCATE(manSol_f(nCells), manSol_s(nCells), Temp_f(nCells), Temp_s(nCells), STAT=errorFlag)
			IF (errorFlag /= 0) THEN
				WRITE(*,*) "Error: Could not allocate manSol, Temp!"
				STOP
			END IF
			
			dx = height/nCells
			sigma = u_f*dt/dx
			d_f = alpha_f*dt/(dx*dx); d_s = alpha_s*dt/(dx*dx)

			DO j = 1,nCells
				manSol_f(j) = cos(k_f*dx*(j-1./2.)) ! manufactured solution cos(kx) for fluid phase
				manSol_s(j) = cos(k_s*dx*(j-1./2.)) ! manufactured solution cos(2kx) for solid phase
			END DO

			WRITE(*,"(A, 1X, F3.0, 1X, A)") "For 2^", LOG(REAL(nCells))/LOG(2.), "cells,"
			
			! Take manSol as initial condition and evolves it over time to reach steady state
			WRITE(*,*) "Fluid phase:"
			Temp_f(:) = SteadyState(manSol_f, sigma, d_f, 1., nCells, dx, dt, .TRUE., k_f)
			discErr_f(:) = ErrorNorms(nCells, Temp_f, manSol_f, .TRUE.)
			WRITE(*,*) "discL2 =", discErr_f(2), "dicsInf =", discErr_f(3), "at", INT(discErr_f(4))
			
			WRITE(*,*) "Solid phase:"
			Temp_s(:) = SteadyState(manSol_s, 0., d_s, 0., nCells, dx, dt, .TRUE., k_s)			
			discErr_s(:) = ErrorNorms(nCells, Temp_s, manSol_s, .TRUE.)
			WRITE(*,*) "discL2 =", discErr_s(2), "dicsInf =", discErr_s(3), "at", INT(discErr_s(4))
			WRITE(*,"(A)")
			
			! Store the error norms for plotting later
			errArr_f(1, i+1) = LOG10(dx)
			errArr_f(2, i+1) = LOG10(discErr_f(1))
			errArr_f(3, i+1) = LOG10(discErr_f(2))
			errArr_f(4, i+1) = LOG10(discErr_f(3))
			locArr_f(1, i+1) = LOG10(dx)
			locArr_f(2, i+1) = discErr_f(4)/REAL(nCells)
			
			errArr_s(1, i+1) = LOG10(dx)
			errArr_s(2, i+1) = LOG10(discErr_s(1))
			errArr_s(3, i+1) = LOG10(discErr_s(2))
			errArr_s(4, i+1) = LOG10(discErr_s(3))
			locArr_s(1, i+1) = LOG10(dx)
			locArr_s(2, i+1) = discErr_s(4)/REAL(nCells)

			! Visualize approximate solutions and manufactured solutions
			IF (nCells == 8) THEN
				
				ALLOCATE(solArr(3,nCells), STAT=errorFlag)
				IF (errorFlag /= 0) THEN
					WRITE(*,*) "Error: Could not allocate solArr!"
					STOP
				END IF
				
				DO j = 1,nCells
					solArr(1,j) = j
					solArr(2,j) = manSol_f(j)
					solArr(3,j) = Temp_f(j)
				END DO
				label(1) = "x_i"; label(2) = "manSol_f"; label(3) = "appSol"
				CALL PlotFigure(1, "Temp_f.dat", "x_i", "temperature", label, 3, nCells, solArr)
				
				DO j = 1,nCells ! reuse solArr
					solArr(2,j) = manSol_s(j)
					solArr(3,j) = Temp_s(j)
				END DO
				label(2) = "manSol_s"
				CALL PlotFigure(2, "Temp_s.dat", "x_i", "temperature", label, 3, nCells, solArr)
				
				DEALLOCATE(solArr, STAT=errorFlag)
				IF (errorFlag /= 0) THEN
					WRITE(*,*) "Error: Could not allocate solArr!"
					STOP
				END IF

			END IF

			DEALLOCATE(manSol_f, manSol_s, Temp_f, Temp_s, STAT=errorFlag)
			IF (errorFlag /= 0) THEN
				WRITE(*,*) "Error: Could not deallocate manSol, Temp!"
				STOP
			END IF
			
			nCells = nCells*2
			i = i+1

		END DO
	
		label2(1) = "dx"; label2(2) = "L1"; label2(3) = "L2"; label2(4) = "Inf"
		label3(1) = "dx"; label3(2) = "Inf"
		CALL PlotFigure(3, "discErr_f.dat", "log10(dx)", "log10(error)", label2, 4, pt, errArr_f)
		CALL PlotFigure(4, "discLoc_f.dat", "log10(dx)", "RelPos", label3, 2, pt, locArr_f)
		
		CALL PlotFigure(7, "discErr_s.dat", "log10(dx)", "log10(error)", label2, 4, pt, errArr_s)
		CALL PlotFigure(8, "discLoc_s.dat", "log10(dx)", "RelPos", label3, 2, pt, locArr_s)

	END SUBROUTINE OrderVerification

END MODULE Code_Verification
