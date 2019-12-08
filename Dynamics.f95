MODULE Dynamics

	USE Input_Output
	USE Measurements

	PRIVATE

	PUBLIC EvolveTemp, PointImplicitMethod

CONTAINS

	SUBROUTINE EvolveTemp(Temp, nCells, Temp_in, phase, charge, MMS)

		IMPLICIT NONE

		REAL, INTENT(INOUT) :: Temp(nCells)
		INTEGER, INTENT(IN) :: nCells
		REAL, INTENT(IN) :: Temp_in
		CHARACTER(LEN=5), INTENT(IN) :: phase
		CHARACTER(LEN=6), INTENT(IN) :: charge
		LOGICAL, INTENT(IN) :: MMS

		INTEGER :: i
		REAL :: dx, sigma, d, k, k_f, k_s, h_c, flux, manSol, allFlux(nCells)
		
		dx = height/nCells
		k_f = 2*Pi*waveNum/height; k_s = 2*k_f ! for MMS only
		IF (phase == "fluid") THEN
			sigma = u_f*dt/dx
			d = alpha_f*dt/(dx*dx)
			k = k_f ! for MMS only
			h_c = h_vf*dt/dx ! for MMS only
		ELSE IF (phase == "solid") THEN
			sigma = 0
			d = alpha_s*dt/(dx*dx)
			k = k_s ! for MMS only
			h_c = - h_vs*dt/dx ! for MMS only
		ELSE
			WRITE(*,*) "Error: phase can be only either 'fluid' or 'solid'!"
			STOP
		END IF

		! Calculate all fluxes to be added to Temp
		allFlux(:) = 0
		IF (MMS .AND. charge == "charge") THEN
			! Method of manufactured solution with T = cos(kx)

			! Boundary face where flux comes in
			flux = -sigma*(-Temp(2)/2. + Temp(1)/2. + Temp_in) ! assume no conductive flux
			manSol = sigma
			allFlux(1) = allFlux(1) - flux - manSol
			
			! Interior faces
			DO i = 1,nCells-1
				flux = -sigma*Temp(i) + d*(Temp(i+1) - Temp(i))
				manSol = sigma*cos(k*dx*i) + d*dx*k*sin(k*dx*i) &
					+ h_c*(sin(k_f*dx*i)/k_f - sin(k_s*dx*i)/k_s)
				allFlux(i) = allFlux(i) + flux + manSol
				allFlux(i+1) = allFlux(i+1) - flux - manSol
			END DO
			
			! Boundary face where flux goes out
			flux = -sigma*Temp(nCells) ! assume no conductive flux
			manSol = sigma
			allFlux(nCells) = allFlux(nCells) + flux + manSol

		ELSE IF (charge == "charge") THEN ! MMS .EQV. .FALSE.
			! Real simulation

			! Boundary face where flux comes in
			flux = -sigma*(-Temp(2)/2. + Temp(1)/2. + Temp_in) ! assume no conductive flux
			allFlux(1) = allFlux(1) - flux

			! Interior faces
			DO i = 1,nCells-1
				flux = -sigma*Temp(i) + d*(Temp(i+1) - Temp(i))
				allFlux(i) = allFlux(i) + flux
				allFlux(i+1) = allFlux(i+1) - flux
			END DO
			
			! Boundary face where flux goes out
			flux = -sigma*Temp(nCells) ! assume no conductive flux
			allFlux(nCells) = allFlux(nCells) + flux

		ELSE IF (charge == "dischg") THEN ! MMS .EQV. .FALSE.
			! Reflection of above case to implement upwind scheme

			flux = -sigma*(-Temp(nCells-1)/2. + Temp(nCells)/2. + Temp_in)
			allFlux(nCells) = allFlux(nCells) - flux
			
			DO i = 1,nCells-1
				flux = -sigma*Temp(nCells-i+1) + d*(Temp(nCells-i) - Temp(nCells-i+1))
				allFlux(nCells-i+1) = allFlux(nCells-i+1) + flux
				allFlux(nCells-i) = allFlux(nCells-i) - flux
			END DO
			
			flux = -sigma*Temp(1)
			allFlux(1) = allFlux(1) + flux

		ELSE IF(charge == "idling") THEN ! MMS .EQV. .FALSE.

			! During idle phase, only diffusion takes place
			DO i = 1,nCells-1
				flux = d*(Temp(i+1) - Temp(i))
				allFlux(i) = allFlux(i) + flux
				allFlux(i+1) = allFlux(i+1) - flux
			END DO
			
		ELSE
			WRITE(*,*) "Error: should not fall in this case!"
			STOP
		END IF

		! Update Temp
		DO i = 1,nCells
			Temp(i) = Temp(i) + allFlux(i)
		END DO

	END SUBROUTINE EvolveTemp

	SUBROUTINE PointImplicitMethod(Temp_f, Temp_s, nCells)

		IMPLICIT NONE

		REAL, INTENT(INOUT) :: Temp_f(nCells), Temp_s(nCells)
		INTEGER, INTENT(IN) :: nCells

		REAL :: det, matrix(2,2), var1, var2
		INTEGER :: i
		
		det = 1. + (h_vf+h_vs)*dt

		! Below order is same as the order of the matrix
		! but indices are adapted for data structure of fortran arrays
		matrix(1,1) = (1.+h_vs*dt)/det; matrix(2,1) = h_vf*dt/det
		matrix(1,2) = h_vs*dt/det; matrix(2,2) = (1.+h_vf*dt)/det

		DO i = 1,nCells
			var1 = Temp_f(i)
			var2 = Temp_s(i)
			Temp_f(i) = matrix(1,1)*var1 + matrix(2,1)*var2
			Temp_s(i) = matrix(1,2)*var1 + matrix(2,2)*var2
		END DO

	END SUBROUTINE PointImplicitMethod

END MODULE Dynamics
