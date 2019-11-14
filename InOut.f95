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

END MODULE Input_Output
