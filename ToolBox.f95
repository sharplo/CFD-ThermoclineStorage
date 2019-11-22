MODULE Tools

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

END MODULE Tools
