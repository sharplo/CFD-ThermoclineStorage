PROGRAM main

	USE Input_Output
	USE Order_Verification
	USE Thermocline_Cycle

	IMPLICIT NONE

	INTEGER :: errorFlag, nCells

	NAMELIST /inputs/ study, nCells, ErrThd, volume, diameter, dt, &
		rho_f, C_f, rho_s, C_s, eps, d_s, k_s, k_f, mu_f, dm, &
		height, u_f, alpha_f, alpha_s, h_v, &
		waveNum, pt_i, pt_f, MaxTStep

	! Parameters to be used when not properly specified in "parameters.dat"
	ErrThd = 1E-5; nCells = 1024; MaxTStep = 1E7; waveNum = 1; pt_i = 3; pt_f = 6
	volume = 300.; diameter = 4; dt = 10.
	rho_f = 1835.6; C_f = 1511.8; rho_s = 2600.; C_s = 900; eps = 0.4
	d_s = 0.03; k_s = 2.; k_f = 0.52; mu_f = 2.63; dm = 10.

	height = volume*4/(Pi*diameter**2); u_f = dm/(rho_f*eps*Pi*diameter**2/4)
	alpha_f = k_f/(eps*rho_f*C_f); alpha_s = k_s/((1.-eps)*rho_s*C_s)
	Pr = mu_f*C_f/k_f; Re = eps*rho_f*u_f*d_s/mu_f
	Nu = 0.255/eps*Pr**(1./3)*Re**(2./3)
	h_fs = Nu*k_f/d_s; h = 1/(1/h_fs + d_s/(10*k_s)); h_v = 6*(1-eps)*h/d_s
	h_vf = h_v/(eps*rho_f*C_f); h_vs = h_v/((1.-eps)*rho_s*C_s)
		
	OPEN(1, FILE="parameters.dat", STATUS="old", IOSTAT=errorFlag)
	IF (errorFlag /= 0) THEN
		WRITE(*,*) "Error: Could not open file!"
		STOP
	END IF
	READ(1, inputs)
	CLOSE(1, IOSTAT=errorFlag)
	IF (errorFlag /= 0) THEN
		WRITE(*,*) "Error: Could not close file!"
		STOP
	END IF
	WRITE(*,inputs)

	IF (TRIM(study) == "OVS") THEN

		h_vf = h_v/(eps*rho_f*C_f); h_vs = h_v/((1.-eps)*rho_s*C_s)
		
		! Order verification study with source term
		CALL GridRefinement()

	ELSE IF (TRIM(study) == "exSol") THEN

		u_f = dm/(rho_f*eps*Pi*diameter**2/4)
		alpha_f = k_f/(eps*rho_f*C_f); alpha_s = k_s/((1.-eps)*rho_s*C_s)
		Pr = mu_f*C_f/k_f; Re = eps*rho_f*u_f*d_s/mu_f
		Nu = 0.255/eps*Pr**(1./3)*Re**(2./3)
		h_fs = Nu*k_f/d_s; h = 1/(1/h_fs + d_s/(10*k_s)); h_v = 6*(1-eps)*h/d_s
		h_vf = h_v/(eps*rho_f*C_f); h_vs = h_v/((1.-eps)*rho_s*C_s)

		! Visualize thermocline motions of one cycle
!		CALL VisualMotion(2**10, 873., 293., (/21600., 21600., 21600., 21600./))

		! Compare with exact solutions
		CALL CompareExSol(nCells, 873., 293., (/5000., 0., 0., 0./))

	ELSE IF (TRIM(study) == "realSim") THEN
		
		height = volume*4/(Pi*diameter**2); u_f = dm/(rho_f*eps*Pi*diameter**2/4)
		alpha_f = k_f/(eps*rho_f*C_f); alpha_s = k_s/((1.-eps)*rho_s*C_s)
		Pr = mu_f*C_f/k_f; Re = eps*rho_f*u_f*d_s/mu_f
		Nu = 0.255/eps*Pr**(1./3)*Re**(2./3)
		h_fs = Nu*k_f/d_s; h = 1/(1/h_fs + d_s/(10*k_s)); h_v = 6*(1-eps)*h/d_s
		h_vf = h_v/(eps*rho_f*C_f); h_vs = h_v/((1.-eps)*rho_s*C_s)

		! Visualize thermocline motions of the first cycle
!		CALL VisualMotion(2**10, 873., 293., (/21600., 21600., 21600., 21600./))
		
		! Real simulation
		CALL SteadyCycle(nCells, 873., 293., (/21600., 21600., 21600., 21600./))

	ELSE
		WRITE(*,*) "Error: 'study' does not have this option!"
		STOP
	END IF

END PROGRAM main
