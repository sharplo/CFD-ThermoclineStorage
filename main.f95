PROGRAM main

	USE Input_Output
	USE Measurements
	USE Dynamics
	USE Order_Verification
	USE Thermocline_Cycle

	! Order verification study with/without source term
!	CALL GridRefinement()

	! Visualize thermocline motions
!	CALL VisualMotion(2**12, 873., 293., (/21600., 21600., 21600., 21600./))

	! Compare with exact solutions
!	CALL CompareExSol(2**12, 873., 293., (/5000., 0., 0., 0./))

	! Real simulation
	CALL SteadyCycle(2**12, 873., 293., (/21600., 21600., 21600., 21600./))


CONTAINS

	SUBROUTINE SteadyCycle(nCells, Temp_c, Temp_d, duration)

		IMPLICIT NONE

		REAL, INTENT(IN)  :: Temp_d, Temp_c, duration(4)
		INTEGER, INTENT(IN) :: nCells

		REAL :: Temp_f(nCells), Temp_s(nCells), dx, time, &
			enEff_c, enEff_d, preEff, cycEff, storEn, maxEn, capFact
		INTEGER :: nCycle

		dx = height/nCells
		Temp_f(:) = Temp_d; Temp_s(:) = Temp_d
		
		enEff_c = 0.
		enEff_d = 0.
		cycEff = 0.
		preEff = 0.

		nCycle = 0
		DO WHILE (ABS(cycEff - preEff) >= ErrThd*preEff)
		
			preEff = cycEff

			! Charging phase
			time = 0.
			! Trapezoid rule: boundary term
			enEff_c = enEff_c + Temp_f(nCells) - Temp_f(1) - Temp_ref*LOG(Temp_f(nCells)/Temp_f(1))
			DO WHILE (time < duration(1))
				CALL EvolveTemp(Temp_f, nCells, Temp_c, "fluid", "charge", .FALSE.)
				CALL EvolveTemp(Temp_s, nCells, 0., "solid", "charge", .FALSE.)
				CALL PointImplicitMethod(Temp_f, Temp_s, nCells)
				time = time + dt
				! Trapezoid rule: interior terms
				enEff_c = enEff_c + 2*( Temp_f(nCells) - Temp_f(1) - Temp_ref*LOG(Temp_f(nCells)/Temp_f(1)) )
			END DO
			! Trapezoid rule: boundary term was double counted
			enEff_c = enEff_c - ( Temp_f(nCells) - Temp_f(1) - Temp_ref*LOG(Temp_f(nCells)/Temp_f(1)) )
		
			! Idling phase
			time = 0.
			DO WHILE (time < duration(2))
				CALL EvolveTemp(Temp_f, nCells, 0., "fluid", "idling", .FALSE.)
				CALL EvolveTemp(Temp_s, nCells, 0., "solid", "idling", .FALSE.)
				CALL PointImplicitMethod(Temp_f, Temp_s, nCells)
				time = time + dt
			END DO

			storEn = ThermalEnergy(nCells, Temp_f, Temp_s, Temp_d)

			! Discharging phase
			time = 0.
			enEff_d = enEff_d + Temp_f(nCells) - Temp_f(1) - Temp_ref*LOG(Temp_f(nCells)/Temp_f(1))
			DO WHILE (time < duration(3))
				CALL EvolveTemp(Temp_f, nCells, Temp_d, "fluid", "dischg", .FALSE.)
				CALL EvolveTemp(Temp_s, nCells, 0., "solid", "dischg", .FALSE.)
				CALL PointImplicitMethod(Temp_f, Temp_s, nCells)
				time = time + dt
				enEff_d = enEff_d + 2*( Temp_f(nCells) - Temp_f(1) - Temp_ref*LOG(Temp_f(nCells)/Temp_f(1)) )
			END DO
			enEff_d = enEff_d - ( Temp_f(nCells) - Temp_f(1) - Temp_ref*LOG(Temp_f(nCells)/Temp_f(1)) )
		
			! Idling phase
			time = 0.
			DO WHILE (time < duration(4))
				CALL EvolveTemp(Temp_f, nCells, 0., "fluid", "idling", .FALSE.)
				CALL EvolveTemp(Temp_s, nCells, 0., "solid", "idling", .FALSE.)
				CALL PointImplicitMethod(Temp_f, Temp_s, nCells)
				time = time + dt
			END DO

			cycEff = enEff_d/enEff_c
			nCycle = nCycle + 1
			
			IF (enEff_c == 0) THEN
				WRITE(*,*) "Error: Division of 0!"
				STOP
			ELSE IF (MOD(nCycle,100) == 0) THEN
				WRITE(*,*) "nCycle", nCycle, "cycEff", cycEff
			END IF

		END DO
		PRINT*, "Total nCycle", nCycle

		storEn = storEn - ThermalEnergy(nCells, Temp_f, Temp_s, Temp_d)
		maxEn = MaxEnergyStored(Temp_c, Temp_d)
		capFact = storEn/maxEn
			
		! Undergo the charging phase once more to get Temperature outflow
		time = 0.
		DO WHILE (time < duration(1))
			CALL EvolveTemp(Temp_f, nCells, Temp_c, "fluid", "charge", .FALSE.)
			CALL EvolveTemp(Temp_s, nCells, 0., "solid", "charge", .FALSE.)
			CALL PointImplicitMethod(Temp_f, Temp_s, nCells)
			time = time + dt
		END DO
		
		WRITE(*, "(A)")
		WRITE(*,*) "Temperature increase at outflow", Temp_f(nCells) - Temp_d
		WRITE(*,*) "Cycle energy efficiency", cycEff
		WRITE(*,*) "Capacity factor", capFact

	END SUBROUTINE SteadyCycle

END PROGRAM main
