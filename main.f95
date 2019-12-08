PROGRAM main
	
	USE Order_Verification
	USE Thermocline_Cycle

	IMPLICIT NONE
	
!	CALL GridRefinement()

	! Part 4
	CALL VisualMotion(2**8, 288., 773., 400, 0, 300, 0)
	
	! Part 5
	CALL StorageCycle(2**8, 288., 773., (/80., 0., 60., 0./))

END PROGRAM main
