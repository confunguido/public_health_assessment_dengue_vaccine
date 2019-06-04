source('functions_importation.R')

reps = 100

HumanImmunityDays = rep(686,2)
EmergenceFactor = c(0,0)
MosquitoLifespan = c(.1,.1)
MosquitoInfectiousness = c(0,.001)
MosquitoMoveProbability = rep(0.3,2)
MosquitoRestDays = rep(2/3,2)
Seed = c(1,as.integer(Sys.time()))
SeedInf = Seed

scalars = sobol(vars = list(
  HumanImmunityDays = HumanImmunityDays,
  EmergenceFactor = EmergenceFactor,
  MosquitoLifespan = MosquitoLifespan,
  MosquitoInfectiousness = MosquitoInfectiousness,
  MosquitoMoveProbability = MosquitoMoveProbability,
  MosquitoRestDays = MosquitoRestDays,
  Seed = Seed,
  SeedInf = SeedInf),
  reps)

scalars = rbind(scalars,scalars,scalars,scalars,scalars)
scalars = cbind(scalars,ForceOfImportation = c(
  rep(-log(.1)/(9*365+182),reps),
  rep(-log(.3)/(9*365+182),reps),
  rep(-log(.5)/(9*365+182),reps),
  rep(-log(.7)/(9*365+182),reps),
  rep(-log(.9)/(9*365+182),reps)))

for(ii in 1 : nrow(scalars)){
  writeSimControl.sobol(
    simControlNum = ii,
    reps = 1,
    NumDays = 365 * 50,
    Vaccination = 0,
    OutputPath = '.',
    LocationFile = 'locations_20150801.csv',
    VaccineProfileFile = 'vaccine_profile_placebo.csv',
    DemographyFile = 'demo_rates_peru.csv',
    scalars = scalars[ii,])
  
  writeSimControl.sobol(
    simControlNum = nrow(scalars) + ii,
    reps = 1,
    NumDays = 365 * 50,
    VaccineDay = 365 * 20,
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
    NumDays = 365 * 50,
    VaccineDay = 365 * 20,
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
    NumDays = 365 * 50,
    VaccineDay = 365 * 20,
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
    NumDays = 365 * 50,
    VaccineDay = 365 * 20,
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
