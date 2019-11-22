MODULE Input_Output
	
	IMPLICIT NONE
	
	REAL, PARAMETER :: Pi = 4.*ATAN(1.), ErrThd = 8E-3, &
		height = 10.5, u_f = 1E-1, alpha_f = 2E-7, alpha_s = 9E-7, dt = 1.51E-2, &
		h_v = 1000, rho_f = 1899, C_f = 1495, rho_s = 2600, C_s = 900, eps = 0.4
	INTEGER, PARAMETER :: MaxTStep = 1E7, waveNum = 1, pt = 10

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

	SUBROUTINE PlotFigure(fileUnit, fileName, xLabel, yLabel, lineLabel, nx, ny, arr)
		
		IMPLICIT NONE
		
		INTEGER, INTENT(IN) :: fileUnit, nx, ny
		CHARACTER(LEN=*), INTENT(IN) :: fileName, xLabel, yLabel, lineLabel(nx)
		REAL, INTENT(IN) :: arr(nx,ny)

		INTEGER :: errorFlag, i, j
		
		IF (fileUnit == 0 .OR. fileUnit == 5 .OR. fileUnit == 6) THEN
			WRITE(*,*) "ERROR: fileUnit = 0,5,6 are reserved!"
			STOP
		END IF

		OPEN(UNIT=fileUnit, FILE=TRIM(fileName), IOSTAT=errorFlag)
		IF (errorFlag /= 0) THEN
			WRITE(*,*) "ERROR: Could not open file!"
			STOP
		END IF
		
		! Write labels
		WRITE(fileUnit,"(A, 1X, A)") TRIM(xLabel), TRIM(yLabel)
		DO i = 1,nx
			WRITE(fileUnit,"(1X, A, $)") TRIM(lineLabel(i))
		END DO
		WRITE(fileUnit, "(A)")
		
		! Write data
		DO j = 1,ny 
			WRITE(fileUnit,"(E17.10, $)") arr(1,j)
			DO i = 2,nx
				WRITE(fileUnit, "(1X, E17.10, $)") arr(i,j)
			END DO
			WRITE(fileUnit, "(A)") 
		END DO

		CLOSE(UNIT=fileUnit, IOSTAT=errorFlag)
		IF (errorFlag /= 0) THEN
			WRITE(*,*) "ERROR: Could not close file!"
			STOP
		END IF

	END SUBROUTINE PlotFigure

	SUBROUTINE VisualTemp(nCells, Temp_f, Temp_s, MMS, manSol_f, manSol_s)

		IMPLICIT NONE

		INTEGER, INTENT(IN) :: nCells
		REAL, INTENT(IN) :: Temp_f(nCells), Temp_s(nCells)
		REAL, INTENT(IN), OPTIONAL :: manSol_f(nCells), manSol_s(nCells)
		LOGICAL, INTENT(IN) :: MMS

		REAL, ALLOCATABLE :: solArr(:,:)
		CHARACTER(LEN=6), ALLOCATABLE :: label(:)
		INTEGER :: errorFlag, i

		IF (MMS) THEN
			
			ALLOCATE(solArr(3,nCells), label(3), STAT=errorFlag)
			IF (errorFlag /= 0) THEN
				WRITE(*,*) "Error: could not allocate solArr, label!"
				STOP
			END IF

			DO i = 1,nCells
				solArr(1,i) = i
				solArr(2,i) = manSol_f(i)
				solArr(3,i) = Temp_f(i)
			END DO
			label(1) = "x_i"; label(2) = "manSol"; label(3) = "appSol"
			CALL PlotFigure(1, "Temp_f.dat", "x_i", "temperature", label, 3, nCells, solArr)
				
			DO i = 1,nCells ! reuse solArr
				solArr(2,i) = manSol_s(i)
				solArr(3,i) = Temp_s(i)
			END DO
			CALL PlotFigure(2, "Temp_s.dat", "x_i", "temperature", label, 3, nCells, solArr)
		
		ELSE

			ALLOCATE(solArr(2,nCells), label(2), STAT=errorFlag)
			IF (errorFlag /= 0) THEN
				WRITE(*,*) "Error: could not allocate solArr, label!"
				STOP
			END IF

			DO i = 1,nCells
				solArr(1,i) = i
				solArr(2,i) = Temp_f(i)
			END DO
			label(1) = "x_i"; label(2) = "appSol"
			CALL PlotFigure(1, "Temp_f.dat", "x_i", "temperature", label, 2, nCells, solArr)
				
			DO i = 1,nCells ! reuse solArr
				solArr(2,i) = Temp_s(i)
			END DO
			CALL PlotFigure(2, "Temp_s.dat", "x_i", "temperature", label, 2, nCells, solArr)
		
		END IF

	END SUBROUTINE VisualTemp

END MODULE Input_Output
