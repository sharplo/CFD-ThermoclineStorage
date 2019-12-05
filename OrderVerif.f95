MODULE Order_Verification
	
	USE Input_Output
	USE Tools
	USE Dynamics
	USE OMP_LIB

	PRIVATE

	PUBLIC GridRefinement

CONTAINS

	SUBROUTINE SteadyState(Temp_f, Temp_s, nCells, Temp_in, charge, MMS, manSol_f, manSol_s)

		IMPLICIT NONE

		REAL, INTENT(INOUT) :: Temp_f(nCells), Temp_s(nCells)
		INTEGER, INTENT(IN) :: nCells
		REAL, INTENT(IN) :: Temp_in, manSol_f(nCells), manSol_s(nCells)
		CHARACTER(LEN=6), INTENT(IN) :: charge
		LOGICAL, INTENT(IN) :: MMS

		REAL :: preT_f(nCells), preT_s(nCells), stdyErr_f(4), stdyErr_s(4), discErr_f(4), discErr_s(4)
		INTEGER*8 :: tStep

		DO tStep = 1,MaxTStep

			! Record previous Temp for calculating dT/dt
			preT_f(:) = Temp_f(:)
			preT_s(:) = Temp_s(:)

			! Iteratate both phases for one time step
			CALL EvolveTemp(Temp_f, nCells, Temp_in, "fluid", charge, MMS)
			CALL EvolveTemp(Temp_s, nCells, 0., "solid", charge, MMS)
			
			! Couple the two phases usiing point-implicit method
			CALL PointImplicitMethod(Temp_f, Temp_s, nCells)

			! Calculate error norms for calculating dT/dt
			stdyErr_f(:) = ErrorNorms(nCells, Temp_f, preT_f, .FALSE.)
			stdyErr_s(:) = ErrorNorms(nCells, Temp_s, preT_s, .FALSE.)
				
			! Calculate the difference between approximate solutions and manufactured solutions
			discErr_f(:) = ErrorNorms(nCells, Temp_f, manSol_f, .TRUE.)
			discErr_s(:) = ErrorNorms(nCells, Temp_s, manSol_s, .TRUE.)
			
			! Address the progress to user
			IF (mod(tStep,MaxTStep/10) == 0) THEN
				WRITE(*,"(A, 1X, F3.0, 1X, A)") "For 2^", LOG(REAL(nCells))/LOG(2.), "cells,"
				WRITE(*,*) "After step", tStep
				WRITE(*,*) "Fluid phase:"
				WRITE(*,*) "dT/dt =", stdyErr_f(2)/dt, "stdyInfLoc =", INT(stdyErr_f(4))
				WRITE(*,*) "discL2 =", discErr_f(2), "dicsInfLoc =", INT(discErr_f(4))
				WRITE(*,*) "Solid phase:"
				WRITE(*,*) "dT/dt =", stdyErr_s(2)/dt, "stdyInfLoc =", INT(stdyErr_s(4))
				WRITE(*,*) "discL2 =", discErr_s(2), "dicsInfLoc =", INT(discErr_s(4))
			END IF

			IF (stdyErr_f(2)/dt < ErrThd*discErr_f(2) &
				.AND. stdyErr_s(2)/dt < ErrThd*discErr_s(2)) THEN ! basically dT/dt << discetization error
				WRITE(*,"(A, 1X, F3.0, 1X, A)") "For 2^", LOG(REAL(nCells))/LOG(2.), "cells,"
				WRITE(*,*) "Approximate solution converges after step", tStep
				WRITE(*,*) "Fluid phase:"
				WRITE(*,*) "dT/dt =", stdyErr_f(2)/dt, "stdyInfLoc =", INT(stdyErr_f(4))
				WRITE(*,*) "discL2 =", discErr_f(2), "dicsInfLoc =", INT(discErr_f(4))
				WRITE(*,*) "Solid phase:"
				WRITE(*,*) "dT/dt =", stdyErr_s(2)/dt, "stdyInfLoc =", INT(stdyErr_s(4))
				WRITE(*,*) "discL2 =", discErr_s(2), "dicsInfLoc =", INT(discErr_s(4))
				WRITE(*,"(A)")
				EXIT
			ELSE IF (tStep == MaxTStep) THEN
				WRITE(*,"(A, 1X, F3.0, 1X, A)") "For 2^", LOG(REAL(nCells))/LOG(2.), "cells,"
				WRITE(*,*) "Approximate solution CANNOT converge!"
				STOP
			END IF

		END DO

	END SUBROUTINE SteadyState

	SUBROUTINE GridRefinement()

		IMPLICIT NONE
		
		REAL :: k, dx, discErr_f(4), discErr_s(4)
		REAL, ALLOCATABLE :: manSol_f(:), manSol_s(:), Temp_f(:), Temp_s(:), &
			errArr_f(:,:), errArr_s(:,:), locArr_f(:,:), locArr_s(:,:), solArr(:,:)
		INTEGER :: nCells, errorFlag, pt, i, j
		CHARACTER(LEN=6) :: label(3), label2(4), label3(2)

		k = 2*Pi*waveNum/height	
		pt = pt_f - pt_i + 1 ! number of points in graphs for OVS

		ALLOCATE(errArr_f(4,pt), errArr_s(4,pt), locArr_f(2,pt), locArr_s(2,pt), STAT=errorFlag)
		IF (errorFlag /= 0) THEN
			WRITE(*,*) "Error: Could not allocate errArr, locArr!"
			STOP
		END IF

		!$OMP PARALLEL DO &
		!$OMP PRIVATE(i, j, nCells, dx, manSol_f, manSol_s, Temp_f, Temp_s, discErr_f, discErr_s)
		DO i = 1, pt

			nCells = 2**(pt_i + i - 1)	
			dx = height/nCells

			ALLOCATE(manSol_f(nCells), manSol_s(nCells), &
				Temp_f(nCells), Temp_s(nCells), STAT=errorFlag)
			IF (errorFlag /= 0) THEN
				WRITE(*,*) "Error: Could not allocate manSol, Temp!"
				STOP
			END IF

			! Manufactured solutions
			DO j = 1,nCells
				manSol_f(j) = 2 + cos(k*dx*(j-1./2)) ! T(x) = 2 + cos(kx)
				manSol_s(j) = 2 + cos(2*k*dx*(j-1./2)) ! T(x) = 2 + cos(2kx)
			END DO
			
			! Take manufactured solutions as initial conditions
			Temp_f(:) = manSol_f(:)
			Temp_s(:) = manSol_s(:)

			! Evolve the approximate solutions with coupling to reach steady states
			CALL SteadyState(Temp_f, Temp_s, nCells, 3., "charge", .TRUE., manSol_f, manSol_s)
			
			! Calculate the difference between approximate solutions and manufactured solutions
			discErr_f(:) = ErrorNorms(nCells, Temp_f, manSol_f, .TRUE.)
			discErr_s(:) = ErrorNorms(nCells, Temp_s, manSol_s, .TRUE.)
			
			! Visualize approximate solutions and manufactured solutions
			IF (nCells == 2**3) THEN
				CALL VisualTemp(nCells, Temp_f, Temp_s, .TRUE., manSol_f, manSol_s)
			END IF

			! Store the error norms for plotting later
			errArr_f(1, i) = LOG10(dx)
			errArr_f(2, i) = LOG10(discErr_f(1))
			errArr_f(3, i) = LOG10(discErr_f(2))
			errArr_f(4, i) = LOG10(discErr_f(3))
			locArr_f(1, i) = LOG10(dx)
			locArr_f(2, i) = (discErr_f(4)-1./2)*dx
			
			errArr_s(1, i) = LOG10(dx)
			errArr_s(2, i) = LOG10(discErr_s(1))
			errArr_s(3, i) = LOG10(discErr_s(2))
			errArr_s(4, i) = LOG10(discErr_s(3))
			locArr_s(1, i) = LOG10(dx)
			locArr_s(2, i) = (discErr_s(4)-1./2)*dx

			DEALLOCATE(manSol_f, manSol_s, Temp_f, Temp_s, STAT=errorFlag)
			IF (errorFlag /= 0) THEN
				WRITE(*,*) "Error: Could not deallocate manSol, Temp!"
				STOP
			END IF

		END DO
		!$OMP END PARALLEL DO
	
		label2(1) = "dx"; label2(2) = "L1"; label2(3) = "L2"; label2(4) = "Inf"
		label3(1) = "dx"; label3(2) = "Inf"
		CALL PlotFigure(3, "discErr_f.dat", "log10(dx)", "log10(error)", label2, 4, pt, errArr_f)
		CALL PlotFigure(4, "discLoc_f.dat", "log10(dx)", "height", label3, 2, pt, locArr_f)
		
		CALL PlotFigure(7, "discErr_s.dat", "log10(dx)", "log10(error)", label2, 4, pt, errArr_s)
		CALL PlotFigure(8, "discLoc_s.dat", "log10(dx)", "height", label3, 2, pt, locArr_s)

	END SUBROUTINE GridRefinement

END MODULE Order_Verification
