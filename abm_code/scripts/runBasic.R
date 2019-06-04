

writeSimControl.reps(
	reps = 100,
	NumDays = 365,
	InitialProportionInfected_DENV3 = 5e-4)

system('./DengueSim SimControl.csv')