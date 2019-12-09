MODULE Thermocline_Cycle
	
	USE Input_Output
	USE Dynamics
	USE Order_Verification

	PRIVATE

	PUBLIC StorageState, VisualMotion, CompareExSol

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

	SUBROUTINE VisualMotion(nCells, Temp_c, Temp_d, duration)
		
		IMPLICIT NONE
	
		INTEGER, INTENT(IN) :: nCells
		REAL, INTENT(IN) :: Temp_c, Temp_d, duration(4)

		REAL :: dx, time, Temp_f(nCells), Temp_s(nCells), stgArr_f(5,nCells), stgArr_s(5,nCells)
		INTEGER :: i
		CHARACTER(LEN=6) :: label(5)

		dx = height/nCells
		Temp_f(:) = Temp_d; Temp_s(:) = Temp_d

		! Charging Phase
		time = 0
		DO WHILE (time < duration(1))
			CALL EvolveTemp(Temp_f, nCells, Temp_c, "fluid", "charge", .FALSE.)
			CALL EvolveTemp(Temp_s, nCells, 0., "solid", "charge", .FALSE.)
			CALL PointImplicitMethod(Temp_f, Temp_s, nCells)
			time = time + dt
		END DO
		DO i = 1,nCells
			stgArr_f(1,i) = (i-1./2)*dx
			stgArr_f(2,i) = Temp_f(i)
			stgArr_s(1,i) = (i-1./2)*dx
			stgArr_s(2,i) = Temp_s(i)
		END DO
		
		! Idling phase
		time = 0.
		DO WHILE (time < duration(2))
			CALL EvolveTemp(Temp_f, nCells, 0., "fluid", "idling", .FALSE.)
			CALL EvolveTemp(Temp_s, nCells, 0., "solid", "idling", .FALSE.)
			CALL PointImplicitMethod(Temp_f, Temp_s, nCells)
			time = time + dt
		END DO
		DO i = 1,nCells
			stgArr_f(3,i) = Temp_f(i)
			stgArr_s(3,i) = Temp_s(i)
		END DO
		
		! Discharging phase
		time = 0.
		DO WHILE (time < duration(3))
			CALL EvolveTemp(Temp_f, nCells, Temp_d, "fluid", "dischg", .FALSE.)
			CALL EvolveTemp(Temp_s, nCells, 0., "solid", "dischg", .FALSE.)
			CALL PointImplicitMethod(Temp_f, Temp_s, nCells)
			time = time + dt
		END DO
		DO i = 1,nCells
			stgArr_f(4,i) = Temp_f(i)
			stgArr_s(4,i) = Temp_s(i)
		END DO
		
		! Idling phase
		time = 0.
		DO WHILE (time < duration(4))
			CALL EvolveTemp(Temp_f, nCells, 0., "fluid", "idling", .FALSE.)
			CALL EvolveTemp(Temp_s, nCells, 0., "solid", "idling", .FALSE.)
			CALL PointImplicitMethod(Temp_f, Temp_s, nCells)
			time = time + dt
		END DO
		DO i = 1,nCells
			stgArr_f(5,i) = Temp_f(i)
			stgArr_s(5,i) = Temp_s(i)
		END DO

		label(1) = "height"; label(2) = "charge"; label(3) = "idle_1"; label(4) = "dischg"; label(5) = "idle_2"
		CALL PlotFigure(1, "Motion_f.dat", "height", "temperature", label, 5, nCells, stgArr_f)
		CALL PlotFigure(2, "Motion_s.dat", "height", "temperature", label, 5, nCells, stgArr_s)

	END SUBROUTINE VisualMotion

	SUBROUTINE CompareExSol(nCells, Temp_c, Temp_d, duration)
		
		IMPLICIT NONE
	
		INTEGER, INTENT(IN) :: nCells
		REAL, INTENT(IN) :: Temp_c, Temp_d, duration(4)

		REAL :: dx, time, Temp_f(nCells), Temp_s(nCells), solArr(3,nCells)
		INTEGER :: i
		CHARACTER(LEN=6) :: label(3)

		dx = height/nCells
		Temp_f(:) = Temp_d; Temp_s(:) = Temp_d
		
		label(1) = "height"; label(2) = "fluid"; label(3) = "solid"
		
		! Charging phase
		time = 0.
		DO WHILE (time < duration(1))
			CALL EvolveTemp(Temp_f, nCells, Temp_c, "fluid", "charge", .FALSE.)
			CALL EvolveTemp(Temp_s, nCells, 0., "solid", "charge", .FALSE.)
			CALL PointImplicitMethod(Temp_f, Temp_s, nCells)
			time = time + dt
		END DO
		DO i = 1,nCells
			solArr(1,i) = (i-1./2)*dx
			solArr(2,i) = Temp_f(i)
			solArr(3,i) = Temp_s(i)
		END DO
		CALL PlotFigure(1, "charge.dat", "height", "temperature", label, 3, nCells, solArr)
		
	END SUBROUTINE CompareExSol

END MODULE Thermocline_Cycle
