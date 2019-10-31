PROGRAM main
	
	USE Input_Output
	USE Semi_Discretized_Model
	USE Code_Verification

	IMPLICIT NONE
	
	REAL :: height, diameter, Temp_i, u_f, alpha_f, alpha_s, dx, k, sigma, d_f, d_s, &
		err_L1, err_L2, err_inf, Temp_in
	REAL, DIMENSION(4) :: duration
	REAL, ALLOCATABLE :: Temp_f(:), Temp_s(:), err_f(:)
	INTEGER :: nCells, nCycles, nTSteps, err_inf_loc, waveNum, i
	REAL, DIMENSION(100) :: time, state
	
	REAL, PARAMETER :: Pi = 3.141592654, ErrThd = 1E-4
	INTEGER, PARAMETER :: MaxTStep = 1E5
	
	waveNum = 1
	
	! CALL Input(height, diameter, nCells, Temp_i, u_f, alpha_f, alpha_s, duration, Temp_f, Temp_s, nCycles, nTSteps)

	DO i = 1,100
		time(i) = 60.*i
		state(i) = REAL(StorageState(time(i), (/2500.,500.,2000.,1000./)))
	END DO
	CALL PlotFigure(1, "storstate.dat", "t(s)", time, "State", state, 100)
	
	CALL OrderVerification(waveNum, Pi, ErrThd, MaxTStep)

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
		
		RETURN

	END FUNCTION StorageState

END PROGRAM main
