MODULE Input_Output

	IMPLICIT NONE
	
	REAL, PARAMETER :: Pi = 4.*ATAN(1.), Temp_ref = 288.15
	REAL :: ErrThd, volume, diameter, dt, &
		rho_f, C_f, rho_s, C_s, eps, d_s, k_s, k_f, mu_f, dm, &
		height, u_f, alpha_f, alpha_s, Pr, Re, Nu, h_fs, h, h_v, h_vf, h_vs
	INTEGER :: waveNum, pt_i, pt_f
	INTEGER*8 :: MaxTStep
	CHARACTER(LEN=7) :: study

CONTAINS

	SUBROUTINE PlotFigure(fileUnit, fileName, xLabel, yLabel, lineLabel, nx, ny, arr)
		
		IMPLICIT NONE
		
		INTEGER, INTENT(IN) :: fileUnit, nx, ny
		CHARACTER(LEN=*), INTENT(IN) :: fileName, xLabel, yLabel, lineLabel(nx)
		REAL, INTENT(IN) :: arr(nx,ny)

		INTEGER :: errorFlag, i, j
		
		IF (fileUnit == 0 .OR. fileUnit == 5 .OR. fileUnit == 6) THEN
			WRITE(*,*) "Error: fileUnit = 0,5,6 are reserved!"
			STOP
		END IF

		OPEN(UNIT=fileUnit, FILE=TRIM(fileName), IOSTAT=errorFlag)
		IF (errorFlag /= 0) THEN
			WRITE(*,*) "Error: Could not open file!"
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
			WRITE(*,*) "Error: Could not close file!"
			STOP
		END IF

	END SUBROUTINE PlotFigure

	SUBROUTINE VisualTemp(nCells, Temp_f, Temp_s, MMS, manSol_f, manSol_s)

		IMPLICIT NONE

		INTEGER, INTENT(IN) :: nCells
		REAL, INTENT(IN) :: Temp_f(nCells), Temp_s(nCells)
		REAL, INTENT(IN), OPTIONAL :: manSol_f(nCells), manSol_s(nCells)
		LOGICAL, INTENT(IN) :: MMS

		REAL :: dx
		REAL, ALLOCATABLE :: solArr(:,:)
		CHARACTER(LEN=6), ALLOCATABLE :: label(:)
		INTEGER :: errorFlag, i

		dx = height/nCells

		IF (MMS) THEN
			
			ALLOCATE(solArr(3,nCells), label(3), STAT=errorFlag)
			IF (errorFlag /= 0) THEN
				WRITE(*,*) "Error: could not allocate solArr, label!"
				STOP
			END IF

			DO i = 1,nCells
				solArr(1,i) = (i-1./2)*dx
				solArr(2,i) = manSol_f(i)
				solArr(3,i) = Temp_f(i)
			END DO
			label(1) = "x_i"; label(2) = "manSol"; label(3) = "appSol"
			CALL PlotFigure(1, "Temp_f.dat", "height", "temperature", label, 3, nCells, solArr)
				
			DO i = 1,nCells ! reuse solArr
				solArr(2,i) = manSol_s(i)
				solArr(3,i) = Temp_s(i)
			END DO
			CALL PlotFigure(2, "Temp_s.dat", "height", "temperature", label, 3, nCells, solArr)
		
		ELSE

			ALLOCATE(solArr(2,nCells), label(2), STAT=errorFlag)
			IF (errorFlag /= 0) THEN
				WRITE(*,*) "Error: could not allocate solArr, label!"
				STOP
			END IF

			DO i = 1,nCells
				solArr(1,i) = (i-1./2)*dx
				solArr(2,i) = Temp_f(i)
			END DO
			label(1) = "x_i"; label(2) = "appSol"
			CALL PlotFigure(1, "Temp_f.dat", "height", "temperature", label, 2, nCells, solArr)
				
			DO i = 1,nCells ! reuse solArr
				solArr(2,i) = Temp_s(i)
			END DO
			CALL PlotFigure(2, "Temp_s.dat", "height", "temperature", label, 2, nCells, solArr)
		
		END IF

	END SUBROUTINE VisualTemp

END MODULE Input_Output
