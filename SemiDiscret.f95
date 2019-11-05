MODULE Semi_Discretized_Model

	PRIVATE

	PUBLIC EvolveTempFluid

CONTAINS

	SUBROUTINE EvolveTempFluid(sigma, d_f, Temp_in, nCells, dx, MMS, k, Temp_f)

		IMPLICIT NONE

		REAL, INTENT(IN) :: sigma, d_f, Temp_in, dx, k
		INTEGER, INTENT(IN) :: nCells
		LOGICAL, INTENT(IN) :: MMS
		REAL, INTENT(INOUT) :: Temp_f(nCells)

		INTEGER :: i
		REAL :: flux, manSol
		
		IF (MMS .EQV. .TRUE.) THEN ! using method of manufactured solution with T = cos(kx)
			
			! Boundary face where flux comes in
			flux = -sigma*Temp_in ! assume no conductive flux
			manSol = sigma*cos(k*dx*(-1./2.)) + d_f*dx*k*sin(k*dx*(-1./2.))
			Temp_f(1) = Temp_f(1) - flux - manSol
			
			! Interior faces
			DO i = 1,nCells-1
				flux = -sigma*Temp_f(i) + d_f*(Temp_f(i+1) - Temp_f(i))
				manSol = sigma*cos(k*dx*(i-1./2.)) + d_f*dx*k*sin(k*dx*(i-1./2.))
				Temp_f(i) = Temp_f(i) + flux + manSol
				Temp_f(i+1) = Temp_f(i+1) - flux - manSol
			END DO
			
			! Boundary face where flux goes out
			flux = -sigma*Temp_f(nCells) ! assume no conductive flux
			manSol = sigma*cos(k*dx*(nCells-1./2.)) + d_f*dx*k*sin(k*dx*(nCells-1./2.))
			Temp_f(nCells) = Temp_f(nCells) + flux + manSol

		ELSE ! real simulation
			WRITE(*,*) "Error: nothing here yet!"
			STOP
		END IF

	END SUBROUTINE EvolveTempFluid
	
END MODULE Semi_Discretized_Model
