MODULE Input_Output

	PRIVATE

	PUBLIC Input, PlotFigure

CONTAINS

	SUBROUTINE Input(height, diameter, nCells, Temp_i, u_f, alpha_f, alpha_s, duration, Temp_f, Temp_s, nCycles, nTSteps)
	
	IMPLICIT NONE

	REAL, INTENT(OUT) :: height, diameter, Temp_i, u_f, alpha_f, alpha_s
	REAL, DIMENSION(4), INTENT(OUT) :: duration
	REAL, ALLOCATABLE, INTENT(OUT) :: Temp_f(:), Temp_s(:)
	INTEGER, INTENT(OUT) :: nCells, nCycles, nTSteps
	
	INTEGER :: errorFlag

	WRITE(*,*) "Please input the height and diameter of the cylindrical storage:"
	READ(*,*) height, diameter
	WRITE(*,*) "Please input the number of cells inside the storage:"
	READ(*,*) nCells
	WRITE(*,*) "Please input the initial temperarutre of the fluid and solid phases (equal):"
	READ(*,*) Temp_i
	
	ALLOCATE(Temp_f(nCells), Temp_s(nCells), STAT=errorFlag)
	IF (errorFlag /= 0) THEN
		WRITE(*,*) "Error: could not allocate Temp_f, Temp_s!"
		STOP
	END IF
	
	Temp_f(:) = Temp_i
	Temp_s(:) = Temp_i

	WRITE(*,*) "Please input the durations of the four states:"
	READ(*,*) duration(1), duration(2), duration(3), duration(4)
	WRITE(*,*) "Please input the number of cycles:"
	READ(*,*) nCycles
	WRITE(*,*) "Please input the number of time steps per cycle:"
	READ(*,*) nTSteps

	WRITE(*,*) "Please input the velocity u_f:"
	READ(*,*) u_f
	WRITE(*,*) "Please input the diffusivities for fluid and solid phases respectively:"
	READ(*,*) alpha_f, alpha_s

	END SUBROUTINE Input

	SUBROUTINE PlotFigure(fileUnit, fileName, xLabel, x, yLabel, y, n)
		
		IMPLICIT NONE
		
		INTEGER, INTENT(IN) :: fileUnit, n
		CHARACTER(LEN=*), INTENT(IN) :: fileName, xLabel, yLabel
		REAL, INTENT(IN) :: x(n), y(n)

		INTEGER :: errorFlag, i
		
		IF (fileUnit == 0 .OR. fileUnit == 5 .OR. fileUnit == 6) THEN
			WRITE(*,*) "ERROR: fileUnit = 0,5,6 are reserved!"
			STOP
		END IF

		OPEN(UNIT=fileUnit, FILE=TRIM(fileName), IOSTAT=errorFlag)
		IF (errorFlag /= 0) THEN
			WRITE(*,*) "ERROR: Could not open file!"
			STOP
		END IF
		
		WRITE(fileUnit,"(A, 1X, A)") TRIM(xLabel), TRIM(yLabel)
		DO i=1,n
			WRITE(fileUnit,"(E13.6, 1X, E13.6)") x(i), y(i)
		END DO

		CLOSE(UNIT=fileUnit, IOSTAT=errorFlag)
		IF (errorFlag /= 0) THEN
			WRITE(*,*) "ERROR: Could not close file!"
			STOP
		END IF

	END SUBROUTINE PlotFigure

END MODULE	Input_Output

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
			flux = -sigma*Temp_in ! assume no conductive flux
			manSol = sigma*cos(k*dx*(-1./2.)) + d_f*dx*k*sin(k*dx*(-1./2.))
			Temp_f(1) = Temp_f(1) - flux - manSol
			DO i = 1,nCells-1
				flux = -sigma*Temp_f(i) + d_f*(Temp_f(i+1) - Temp_f(i))
				manSol = sigma*cos(k*dx*(i-1./2.)) + d_f*dx*k*sin(k*dx*(i-1./2.))
				Temp_f(i) = Temp_f(i) + flux + manSol
				Temp_f(i+1) = Temp_f(i+1) - flux - manSol
			END DO
			flux = -sigma*Temp_f(nCells) ! assume no conductive flux
			manSol = sigma*cos(k*dx*(nCells-1./2.)) + d_f*dx*k*sin(k*dx*(nCells-1./2.))
			Temp_f(nCells) = Temp_f(nCells) + flux + manSol
		ELSE ! real simulation
			flux = -sigma*Temp_in ! assume no conductive flux
			Temp_f(1) = Temp_f(1) - flux
			DO i = 1,nCells-1
				flux = -sigma*Temp_f(i) + d_f*(Temp_f(i+1) - Temp_f(i))
				Temp_f(i) = Temp_f(i) + flux
				Temp_f(i+1) = Temp_f(i+1) - flux
			END DO
			flux = -sigma*Temp_f(nCells) ! assume no conductive flux
			Temp_f(nCells) = Temp_f(nCells) + flux
		END IF

	END SUBROUTINE EvolveTempFluid

END MODULE Semi_Discretized_Model

MODULE Code_Verification
	
	USE Input_Output
	USE Semi_Discretized_Model
	
	PRIVATE

	PUBLIC CalErrorFluid, OrderVerification

CONTAINS

	SUBROUTINE CalErrorFluid(nCells, dx, Temp_a, MMS, k, err_f, err_L1, err_L2, err_inf, err_inf_loc)
		
		IMPLICIT NONE

		INTEGER, INTENT(IN) :: nCells
		REAL, INTENT(IN) :: dx, Temp_a(nCells), k
		LOGICAL, INTENT(IN) :: MMS
		REAL, INTENT(OUT) :: err_f(nCells), err_L1, err_L2, err_inf
		INTEGER, INTENT(OUT) :: err_inf_loc

		INTEGER :: i
		REAL :: exSol, err_i
		
		err_L1 = 0; err_L2 = 0; err_inf = -1.

		IF (MMS .EQV. .TRUE.) THEN
			DO i = 1,nCells
				exSol = 2./(k*dx) * cos(k*dx*(i-1)) * sin(k*dx/2.) ! cell average of T(x)=cos(kx)
				err_i = abs( (Temp_a(i)-exSol)/(2./(k*dx)) ) ! typical value
				err_L1 = err_L1 + err_i
				err_L2 = err_L2 + err_i*err_i
				IF (err_i > err_inf) THEN
					err_inf = err_i
					err_inf_loc = i
				END IF
!				WRITE(*,*) i, "error is", err_i, "Temp_a", Temp_a(i), "exSol is", exSol
			END DO
		ELSE
			WRITE(*,*) "Error: No exact solution to compare!"
			STOP
		END IF

		err_L1 = err_L1/nCells
		err_L2 = sqrt(err_L2/nCells)

	END SUBROUTINE CalErrorFluid
	
	SUBROUTINE OrderVerification(waveNum, Pi, ErrThd, MaxTStep)

		IMPLICIT NONE
		
		REAL, INTENT(IN) :: Pi, ErrThd
		INTEGER, INTENT(IN) :: waveNum, MaxTStep
		
		REAL :: height, sigma, d_f, Temp_in, k, dx, err_L2_pre, err_L1, err_L2, err_inf, log_dx(1000), err_string(4,1000)
		INTEGER :: nCells, errorFlag, i, err_inf_loc, n, s
		REAL, ALLOCATABLE :: Temp_f(:), err_f(:)

		height = 16.; sigma = 0.00005; d_f = 0.4; Temp_in = 0.
		n = 260; s = 8
		k = 2*Pi*waveNum/height
		
		DO i = 1,n
			log_dx(i) = LOG10(height/(s*i))
		END DO

		DO nCells = s, 8*n, s
			
			ALLOCATE(Temp_f(nCells), err_f(nCells), STAT=errorFlag)
			IF (errorFlag /= 0) THEN
				WRITE(*,*) "Error: Could not allocate Temp_f, err_f!"
				STOP
			END IF
			
			dx = height/nCells

			DO i = 1,nCells
				! Temp_f(i) = cos(k*dx*(i-1./2.)) ! manufactured solution T = cos(kx)
				Temp_f(i) = 2./(k*dx) * cos(k*dx*(i-1)) * sin(k*dx/2.) ! cell average of T(x)=cos(kx)
			END DO
			
			err_L2_pre = 0.

			DO i = 1,MaxTStep
				
				CALL EvolveTempFluid(sigma, d_f, Temp_in, nCells, dx, .TRUE., k, Temp_f)
				CALL CalErrorFluid(nCells, dx, Temp_f, .TRUE., k, err_f, err_L1, err_L2, err_inf, err_inf_loc)
				
				IF (i > 80 .AND. & ! prevent fluctuations in the beginning
					
					abs(err_L2 - err_L2_pre) < ErrThd * abs(err_L2_pre)) THEN
						
						IF (mod(nCells,8) == 0) THEN ! to reduce computational cost
							WRITE(*,"(A, 1X, I4, 1X, A, 1X, I4)") "For", nCells, "nCells, manufactured solution converges after step", i
							WRITE(*,*) err_L2, err_inf, err_inf_loc
							WRITE(*,*) abs(err_L2 - err_L2_pre), ErrThd * abs(err_L2_pre)
						END IF

					EXIT
				
				ELSE IF (i == MaxTStep) THEN
					WRITE(*,"(A, 1X, I4, 1X, A, 1X, I4)") "For", nCells, "nCells, manufactured solution CANNOT converge."
				END IF
				
				err_L2_pre = err_L2

			END DO
	
			err_string(1, nCells/s) = LOG10(err_L1)
			err_string(2, nCells/s) = LOG10(err_L2)
			err_string(3, nCells/s) = LOG10(err_inf)
			err_string(4, nCells/s) = REAL(err_inf_loc)
			
			DEALLOCATE(Temp_f, err_f, STAT=errorFlag)
			IF (errorFlag /= 0) THEN
				WRITE(*,*) "Error: Could not deallocate Temp_f, err_f!"
				STOP
			END IF

		END DO
		
		CALL PlotFigure(3, "err_L1.dat", "log10(dx)", log_dx, "log10(L1)", err_string(1,:), 8*n/s)
		CALL PlotFigure(4, "err_L2.dat", "log10(dx)", log_dx, "log10(L2)", err_string(2,:), 8*n/s)
		CALL PlotFigure(7, "err_inf.dat", "log10(dx)", log_dx, "log10(L_inf)", err_string(3,:), 8*n/s)
		CALL PlotFigure(8, "err_inf_loc.dat", "log10(dx)", log_dx, "cell", err_string(4,:), 8*n/s)

	END SUBROUTINE OrderVerification

END MODULE Code_Verification

!---------------------------------

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
	
	REAL, PARAMETER :: dt = 2 ! unit is 0.1 second
	REAL, PARAMETER :: Pi = 3.141592654, ErrThd = 1E-3
	INTEGER, PARAMETER :: MaxTStep = 1000
	
	waveNum = 1
	
	! CALL Input(height, diameter, nCells, Temp_i, u_f, alpha_f, alpha_s, duration, Temp_f, Temp_s, nCycles, nTSteps)

	CALL PlotFigure(1, "figdata.dat", "X_Axis", (/1.,2.,3./), "Y_Axix", (/10.,20.,30./), 3)
	
	DO i = 1,100
		time(i) = 60.*i
		state(i) = REAL(StorageState(time(i), (/2500.,500.,2000.,1000./)))
	END DO
	CALL PlotFigure(2, "storstate.dat", "t(s)", time, "State", state, 100)
	
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
