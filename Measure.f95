MODULE Measurements

	USE Input_Output

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

	REAL FUNCTION EnergyFlux(nTStps, Temp, Temp_ref)

		IMPLICIT NONE

		INTEGER, INTENT(IN) :: nTStps
		REAL, INTENT(IN) :: Temp(nTStps), Temp_ref

		INTEGER :: i
		REAL :: var
		REAL :: dm=1 ! need to be deleted
		EnergyFlux = 0

		! Trapezoid rule
		var = Temp(1) - Temp_ref - Temp_ref*LOG(Temp(1)/Temp_ref)
		EnergyFlux = EnergyFlux + var
		DO i = 2,nTStps-1
			var = Temp(i) - Temp_ref - Temp_ref*LOG(Temp(i)/Temp_ref)
			EnergyFlux = EnergyFlux + 2*var
		END DO
		var = Temp(nTStps) - Temp_ref - Temp_ref*LOG(Temp(nTStps)/Temp_ref)
		EnergyFlux = EnergyFlux + var

		EnergyFlux = EnergyFlux*dm*C_f*dt/2

	END FUNCTION EnergyFlux

	REAL FUNCTION ThermalEnergy(nCells, Temp_f, Temp_s, Temp_d)
		
		IMPLICIT NONE

		INTEGER, INTENT(IN) :: nCells
		REAL, INTENT(IN) :: Temp_f(nCells), Temp_s(nCells), Temp_d

		REAL :: dx, sum_f, sum_s

		ThermalEnergy = 0

		dx = height/nCells
		sum_f = SUM(Temp_f)*dx - Temp_d*height
		sum_s = SUM(Temp_s)*dx - Temp_d*height
		
		ThermalEnergy = ThermalEnergy + eps*rho_f*C_f*sum_f + (1-eps)*rho_s*C_s*sum_s
		ThermalEnergy = ThermalEnergy*Pi/4*diameter**2

	END FUNCTION ThermalEnergy

	REAL FUNCTION MaxEnergyStored(Temp_c, Temp_d)

		IMPLICIT NONE

		REAL, INTENT(IN) :: Temp_c, Temp_d

		MaxEnergyStored = eps*rho_f*C_f + (1-eps)*rho_s*C_s
		MaxEnergyStored = MaxEnergyStored*Pi/4*diameter**2
		MaxEnergyStored = MaxEnergyStored*height*(Temp_c - Temp_d)

	END FUNCTION MaxEnergyStored

END MODULE Measurements
