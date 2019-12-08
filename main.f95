PROGRAM main
	
	USE Order_Verification
	USE Thermocline_Cycle

	IMPLICIT NONE

	! Order verification study with/without source term
!	CALL GridRefinement()

	! Visualize thermocline motions
	CALL VisualMotion(2**8, 293., 873., 5E6, 1E6, 3E6, 1E6)

	! Compare with exact solutions
!	CALL StorageCycle(2**12, 293., 873., (/5000., 0., 0., 0./))

	! Real simulation
	CALL StorageCycle(2**12, 293., 873., (/5000., 0., 0., 0./))

END PROGRAM main
