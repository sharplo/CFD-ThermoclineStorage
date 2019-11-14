MODULE Semi_Discretized_Model

	PRIVATE

	PUBLIC EvolveTemp

CONTAINS

	SUBROUTINE EvolveTemp(Temp, sigma, d, Temp_in, nCells, dx, MMS, k)

		IMPLICIT NONE

		REAL, INTENT(INOUT) :: Temp(nCells)
		REAL, INTENT(IN) :: sigma, d, Temp_in, dx, k
		INTEGER, INTENT(IN) :: nCells
		LOGICAL, INTENT(IN) :: MMS

		INTEGER :: i
		REAL :: flux, manSol
		
		IF (MMS .EQV. .TRUE.) THEN ! using method of manufactured solution with T = cos(kx)
			
			! Boundary face where flux comes in
			flux = -sigma*(-Temp(2)/2. + Temp(1)/2. + Temp_in) ! assume no conductive flux
			manSol = sigma
			Temp(1) = Temp(1) - flux - manSol
			
			! Interior faces
			DO i = 1,nCells-1
				flux = -sigma*Temp(i) + d*(Temp(i+1) - Temp(i))
				manSol = sigma*cos(k*dx*i) + d*dx*k*sin(k*dx*i)
				Temp(i) = Temp(i) + flux + manSol
				Temp(i+1) = Temp(i+1) - flux - manSol
			END DO
			
			! Boundary face where flux goes out
			flux = -sigma*Temp(nCells) ! assume no conductive flux
			manSol = sigma
			Temp(nCells) = Temp(nCells) + flux + manSol

		ELSE ! real simulation
			WRITE(*,*) "Error: nothing here yet!"
			STOP
		END IF

	END SUBROUTINE EvolveTemp
	
END MODULE Semi_Discretized_Model
