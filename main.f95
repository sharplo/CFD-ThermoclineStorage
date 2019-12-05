PROGRAM main
	
	USE Input_Output
	USE Dynamics
	USE Order_Verification

	IMPLICIT NONE
	
	REAL :: state(2,100)
	INTEGER :: i, tStps(4)
	CHARACTER(LEN=6) :: label(2)

	! CALL Input(height, diameter, nCells, Temp_i, u_f, alpha_f, alpha_s, duration, Temp_f, Temp_s, nCycles, nTSteps)

	label(1) = "time"; label(2) = "state"
	DO i = 1,100
		state(1,i) = 60.*i
		state(2,i) = REAL(StorageState(state(1,i), (/2500.,500.,2000.,1000./)))
	END DO
	CALL PlotFigure(1, "storstate.dat", "t(s)", "State", label, 2, 100, state)

	CALL GridRefinement()

	! Part 4
!	CALL StorageCycle(2**10, 288., 773., 50000, 10000, 30000, 10000)
	
	! Part 5
!	CALL RealSimulation(2**10, 293., 873., 5000.)

CONTAINS

	FUNCTION PhaseTStps(duration)
		
		IMPLICIT NONE
		
		INTEGER :: PhaseTStps(4)
		REAL, INTENT(IN) :: duration(4)
	
		PhaseTStps(1) = NINT(duration(1)/dt)
		PhaseTStps(2) = NINT(duration(2)/dt)
		PhaseTStps(3) = NINT(duration(3)/dt)
		PhaseTStps(4) = NINT(duration(4)/dt)

	END FUNCTION PhaseTStps

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

	SUBROUTINE StorageCycle(nCells, Temp_0, Temp_in, TStp_chg, TStp_idl1, TStp_dis, TStp_idl2)
		
		IMPLICIT NONE
	
		INTEGER, INTENT(IN) :: nCells, TStp_chg, TStp_dis, TStp_idl1, TStp_idl2
		REAL, INTENT(IN) :: Temp_0, Temp_in

		REAL :: dx, Temp_f(nCells), Temp_s(nCells), stgArr_f(5,nCells), stgArr_s(5,nCells)
		INTEGER :: i
		CHARACTER(LEN=6) :: label(5)

		dx = height/nCells
		Temp_f(:) = Temp_0; Temp_s(:) = Temp_0

		DO i = 1, TStp_chg
			CALL EvolveTemp(Temp_f, nCells, Temp_in, "fluid", "charge", .FALSE.)
			CALL EvolveTemp(Temp_s, nCells, 0., "solid", "charge", .FALSE.)
			CALL PointImplicitMethod(Temp_f, Temp_s, nCells)
		END DO
		DO i = 1,nCells
			stgArr_f(1,i) = (i-1./2)*dx
			stgArr_f(2,i) = Temp_f(i)
			stgArr_s(1,i) = (i-1./2)*dx
			stgArr_s(2,i) = Temp_s(i)
		END DO
		DO i = 1, TStp_idl1
			CALL EvolveTemp(Temp_f, nCells, 0., "fluid", "idling", .FALSE.)
			CALL EvolveTemp(Temp_s, nCells, 0., "solid", "idling", .FALSE.)
			CALL PointImplicitMethod(Temp_f, Temp_s, nCells)
		END DO
		DO i = 1,nCells
			stgArr_f(3,i) = Temp_f(i)
			stgArr_s(3,i) = Temp_s(i)
		END DO
		DO i = 1, TStp_dis
			CALL EvolveTemp(Temp_f, nCells, Temp_0, "fluid", "dischg", .FALSE.)
			CALL EvolveTemp(Temp_s, nCells, 0., "solid", "dischg", .FALSE.)
			CALL PointImplicitMethod(Temp_f, Temp_s, nCells)
		END DO
		DO i = 1,nCells
			stgArr_f(4,i) = Temp_f(i)
			stgArr_s(4,i) = Temp_s(i)
		END DO
		DO i = 1, TStp_idl2
			CALL EvolveTemp(Temp_f, nCells, 0., "fluid", "idling", .FALSE.)
			CALL EvolveTemp(Temp_s, nCells, 0., "solid", "idling", .FALSE.)
			CALL PointImplicitMethod(Temp_f, Temp_s, nCells)
		END DO
		DO i = 1,nCells
			stgArr_f(5,i) = Temp_f(i)
			stgArr_s(5,i) = Temp_s(i)
		END DO
		label(1) = "height"; label(2) = "charge"; label(3) = "idle_1"; label(4) = "dischg"; label(5) = "idle_2"
		CALL PlotFigure(1, "Cycle_f.dat", "height", "temperature", label, 5, nCells, stgArr_f)
		CALL PlotFigure(2, "Cycle_s.dat", "height", "temperature", label, 5, nCells, stgArr_s)

	END SUBROUTINE StorageCycle

	SUBROUTINE RealSimulation(nCells, Temp_0, Temp_in, duration)
		
		IMPLICIT NONE
	
		INTEGER, INTENT(IN) :: nCells
		REAL, INTENT(IN) :: Temp_0, Temp_in, duration(4)

		REAL :: dx, time, Temp_f(nCells), Temp_s(nCells), solArr(3,nCells)
		INTEGER :: i
		CHARACTER(LEN=6) :: label(3)

		dx = height/nCells
		Temp_f(:) = Temp_0; Temp_s(:) = Temp_0
		
		label(1) = "height"; label(2) = "fluid"; label(3) = "solid"
		
		! Charging phase
		time = 0.
		DO WHILE (time < duration(1))
			CALL EvolveTemp(Temp_f, nCells, Temp_in, "fluid", "charge", .FALSE.)
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
		
		! Idling phase
		time = 0.
		DO WHILE (time < duration(2))
			CALL EvolveTemp(Temp_f, nCells, 0., "fluid", "idling", .FALSE.)
			CALL EvolveTemp(Temp_s, nCells, 0., "solid", "idling", .FALSE.)
			CALL PointImplicitMethod(Temp_f, Temp_s, nCells)
			time = time + dt
		END DO
		DO i = 1,nCells
			solArr(2,i) = Temp_f(i)
			solArr(3,i) = Temp_s(i)
		END DO
		CALL PlotFigure(2, "idle_1.dat", "height", "temperature", label, 3, nCells, solArr)
		
		! Discharging phase
		time = 0.
		DO WHILE (time < duration(3))
			CALL EvolveTemp(Temp_f, nCells, Temp_0, "fluid", "dischg", .FALSE.)
			CALL EvolveTemp(Temp_s, nCells, 0., "solid", "dischg", .FALSE.)
			CALL PointImplicitMethod(Temp_f, Temp_s, nCells)
			time = time + dt
		END DO
		DO i = 1,nCells
			solArr(2,i) = Temp_f(i)
			solArr(3,i) = Temp_s(i)
		END DO
		CALL PlotFigure(3, "dischg.dat", "height", "temperature", label, 3, nCells, solArr)
		
		! Idling phase
		time = 0.
		DO WHILE (time < duration(4))
			CALL EvolveTemp(Temp_f, nCells, 0., "fluid", "idling", .FALSE.)
			CALL EvolveTemp(Temp_s, nCells, 0., "solid", "idling", .FALSE.)
			CALL PointImplicitMethod(Temp_f, Temp_s, nCells)
			time = time + dt
		END DO
		DO i = 1,nCells
			solArr(2,i) = Temp_f(i)
			solArr(3,i) = Temp_s(i)
		END DO
		CALL PlotFigure(4, "idle_2.dat", "height", "temperature", label, 3, nCells, solArr)

	END SUBROUTINE RealSimulation

END PROGRAM main
