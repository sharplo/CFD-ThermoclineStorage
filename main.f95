PROGRAM main
	
	USE Input_Output
	USE Dynamics
	USE Order_Verification

	IMPLICIT NONE
	
	REAL :: dx, sigma, d_f, d_s, state(2,100)
	REAL, ALLOCATABLE :: initT_f(:), initT_s(:), Temp_f(:), Temp_s(:), stgArr_f(:,:), stgArr_s(:,:)
	INTEGER :: nCells, nCycles, nTSteps, errorFlag, i
	CHARACTER(LEN=6) :: label(2), label2(3)

	! CALL Input(height, diameter, nCells, Temp_i, u_f, alpha_f, alpha_s, duration, Temp_f, Temp_s, nCycles, nTSteps)

	label(1) = "time"; label(2) = "state"
	DO i = 1,100
		state(1,i) = 60.*i
		state(2,i) = REAL(StorageState(state(1,i), (/2500.,500.,2000.,1000./)))
	END DO
	CALL PlotFigure(1, "storstate.dat", "t(s)", "State", label, 2, 100, state)
	
!	CALL GridRefinement()

	CALL StorageCycle(2**8, 288., 773., 5000, 3000)

CONTAINS

	INTEGER FUNCTION StorageState(time, duration)
		
		IMPLICIT NONE

		REAL, INTENT(IN) :: time
		REAL, DIMENSION(4), INTENT(IN) :: duration
		
		REAL :: equivTime
		
		equivTime = time - sum(duration)*real(floor(time/sum(duration)))
		
		IF (equivTime < duration(1)) THEN
			StorageState = 1
		ELSE IF (duration(1) <= equivTime .AND. equivTime < duration(1)+duration(2)) THEN
			StorageState = 2
		ELSE IF (duration(1)+duration(2) <= equivTime .AND. equivTime < duration(1)+duration(2)+duration(3)) THEN
			StorageState = 3
		ELSE
			StorageState = 4
		END IF

	END FUNCTION StorageState

	SUBROUTINE StorageCycle(nCells, Temp_0, Temp_in, TStp_chg, TStp_dis)
		
		IMPLICIT NONE
	
		INTEGER, INTENT(IN) :: nCells, TStp_chg, TStp_dis
		REAL, INTENT(IN) :: Temp_0, Temp_in

		REAL :: Temp_f(nCells), Temp_s(nCells), stgArr_f(3,nCells), stgArr_s(3,nCells)
		INTEGER :: i
		CHARACTER(LEN=6) :: label(3)

		Temp_f(:) = Temp_0; Temp_s(:) = Temp_0
	
		DO i = 1, TStp_chg
			CALL EvolveTemp(Temp_f, nCells, Temp_in, "fluid", "charge", .FALSE.)
			CALL EvolveTemp(Temp_s, nCells, 0., "solid", "charge", .FALSE.)	
			CALL PointImplicitMethod(Temp_f, Temp_s, nCells)
		END DO
		DO i = 1,nCells
			stgArr_f(1,i) = i
			stgArr_f(2,i) = Temp_f(i)
			stgArr_s(1,i) = i
			stgArr_s(2,i) = Temp_s(i)
		END DO
		DO i = 1, TStp_dis
			CALL EvolveTemp(Temp_f, nCells, Temp_0, "fluid", "dischg", .FALSE.)
			CALL EvolveTemp(Temp_s, nCells, 0., "solid", "dischg", .FALSE.)	
			CALL PointImplicitMethod(Temp_f, Temp_s, nCells)
		END DO
		DO i = 1,nCells
			stgArr_f(3,i) = Temp_f(i)
			stgArr_s(3,i) = Temp_s(i)
		END DO
		label(1) = "x_i"; label(2) = "charge"; label(3) = "dischg"
		CALL PlotFigure(1, "Cycle_f.dat", "x_i", "temperature", label, 3, nCells, stgArr_f)
		CALL PlotFigure(2, "Cycle_s.dat", "x_i", "temperature", label, 3, nCells, stgArr_s)

	END SUBROUTINE StorageCycle

END PROGRAM main
