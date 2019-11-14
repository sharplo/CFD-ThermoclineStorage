PROGRAM main
	
	USE Input_Output
	USE Semi_Discretized_Model
	USE Code_Verification

	IMPLICIT NONE
	
!	REAL :: height, diameter, Temp_i, u_f, alpha_f, alpha_s
	REAL, DIMENSION(4) :: duration
	INTEGER :: nCells, nCycles, nTSteps, i
	REAL, DIMENSION(2,100) :: state
	CHARACTER(LEN=6) :: label(2)

	! CALL Input(height, diameter, nCells, Temp_i, u_f, alpha_f, alpha_s, duration, Temp_f, Temp_s, nCycles, nTSteps)

	label(1) = "time"; label(2) = "state"
	DO i = 1,100
		state(1,i) = 60.*i
		state(2,i) = REAL(StorageState(state(1,i), (/2500.,500.,2000.,1000./)))
	END DO
	CALL PlotFigure(1, "storstate.dat", "t(s)", "State", label, 2, 100, state)
	
	CALL OrderVerification

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

END PROGRAM main
