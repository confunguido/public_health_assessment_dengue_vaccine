source('functions.R')

reps = 1000

HumanImmunityDays = rep(686,2)
EmergenceFactor = c(.01,1.99)
MosquitoLifespan = c(2,6)
MosquitoInfectiousness = c(.01,.4)
MosquitoMoveProbability = rep(0.3,2)
MosquitoRestDays = rep(2/3,2)
ForceOfImportation = c(-log(.01)/(9*365+182),-log(.1)/(9*365+182))
Seed = c(1,as.integer(Sys.time()))
SeedInf = Seed

scalars = sobol(vars = list(
  HumanImmunityDays = HumanImmunityDays,
  EmergenceFactor = EmergenceFactor,
  MosquitoLifespan = MosquitoLifespan,
  MosquitoInfectiousness = MosquitoInfectiousness,
  MosquitoMoveProbability = MosquitoMoveProbability,
  MosquitoRestDays = MosquitoRestDays,
  ForceOfImportation = ForceOfImportation,
  Seed = Seed,
  SeedInf = SeedInf),
  reps)

for(ii in 1 : nrow(scalars)){
  writeSimControl.sobol(
    simControlNum = ii,
    reps = 1,
    NumDays = 365 * 100,
    Vaccination = 0,
    OutputPath = '.',
    LocationFile = 'locations_20150801.csv',
    VaccineProfileFile = 'vaccine_profile_placebo.csv',
    DemographyFile = 'demo_rates_peru.csv',
    scalars = scalars[ii,])
  
  writeSimControl.sobol(
    simControlNum = nrow(scalars) + ii,
    reps = 1,
    NumDays = 365 * 100,
    VaccineDay = 365 * 70,
    Vaccination = 1,
    VaccineCoverage = .8,
    VaccineAge = 9,
    Catchup = 0,
    OutputPath = '.',
    LocationFile = 'locations_20150801.csv',
    VaccineProfileFile = 'vaccine_profile_placebo.csv',
    DemographyFile = 'demo_rates_peru.csv',
    scalars = scalars[ii,])
  
  writeSimControl.sobol(
    simControlNum = 2 * nrow(scalars) + ii,
    reps = 1,
    NumDays = 365 * 100,
    VaccineDay = 365 * 70,
    Vaccination = 1,
    VaccineCoverage = .8,
    VaccineAge = 9,
    Catchup = 1,
    OutputPath = '.',
    LocationFile = 'locations_20150801.csv',
    VaccineProfileFile = 'vaccine_profile_placebo.csv',
    DemographyFile = 'demo_rates_peru.csv',
    scalars = scalars[ii,])
  
  writeSimControl.sobol(
    simControlNum = 3 * nrow(scalars) + ii,
    reps = 1,
    NumDays = 365 * 100,
    VaccineDay = 365 * 70,
    Vaccination = 1,
    VaccineCoverage = .5,
    VaccineAge = 9,
    Catchup = 0,
    OutputPath = '.',
    LocationFile = 'locations_20150801.csv',
    VaccineProfileFile = 'vaccine_profile_placebo.csv',
    DemographyFile = 'demo_rates_peru.csv',
    scalars = scalars[ii,])
  
  writeSimControl.sobol(
    simControlNum = 4 * nrow(scalars) + ii,
    reps = 1,
    NumDays = 365 * 100,
    VaccineDay = 365 * 70,
    Vaccination = 1,
    VaccineCoverage = .8,
    VaccineAge = 16,
    Catchup = 0,
    OutputPath = '.',
    LocationFile = 'locations_20150801.csv',
    VaccineProfileFile = 'vaccine_profile_placebo.csv',
    DemographyFile = 'demo_rates_peru.csv',
    scalars = scalars[ii,])
}
