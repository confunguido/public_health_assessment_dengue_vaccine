#ifndef SIMULATION_H
#define SIMULATION_H

#include <map>
#include <set>
#include <memory>
#include <fstream>
#include <queue>
#include <stdexcept>
#include <array>
// 
#include "defines.h"
#include "Location.h"
#include "Mosquito.h"
#include "RandomNumGenerator.h"
#include "Report.h"
#include "Vaccine.h"
#include "Recruitment.h"
#include "Human.h"
#include "OutbreakResponse.h"

using std::string;
using std::vector;
using std::ifstream;
using std::make_pair;
using traject_t = std::array<vpath_t, N_TRAJECTORY >;


class Simulation {
public:
    RandomNumGenerator rGen;
    RandomNumGenerator rGenInf;
    RandomNumGenerator rGenControl;

    Simulation(string);
    Simulation() = delete;
    Simulation(const Simulation& orig) = delete;
    string readInputs();
    void readTrajectoryFile(string);
    void readBirthsFile(string);
    void readInitialInfectionsFile(string);
    void readSimControlFile(string);
    void readLocationFile(string);
    void readDiseaseRatesFile();
    void readVaccineProfilesFile();
    void readVaccinationGroupsFile();
    void readVaccineSettingsFile();
    void readAegyptiFile(string);
    void readInitialFOI(string);
    void readDailyFOI(string);
    void readDailyLocalFOI(string);
    void setLocNeighborhood(double);
    void setRadiusLocations(double);
    void simEngine();
    void humanDynamics();
    void tests();
    void mosquitoDynamics();
    void generateMosquitoes();
    unsigned setInitialInfection(double, unsigned);
    void simulate();
    void updatePop();
    void selectEligibleTrialParticipants();
    bool checkAgeToVaccinate(int age_);
    //virtual ~Simulation();
private:
    vector<string>getParamsLine(string);
    int parseInteger(string);
    string parseString(string);
    double parseDouble(string);
    void parseVector(string line, vector<int> *);
    void parseVector(string line, vector<double> *);

    map <string,std::unique_ptr<Location>> locations;    
    map<unsigned,double> halflife;
    map<unsigned,double> disRates;
    map<unsigned,double> hospRates;
    map<int,int> ageGroups;
    map<unsigned, Vaccine> vaccines;
    std::multimap<string,std::unique_ptr<Mosquito>> mosquitoes;
    std::multimap<string,sp_human_t> humans;
    std::multimap<int,sp_human_t> future_humans;
    map<string, sp_human_t> total_humans_by_id;

    std::set<std::string> zones;
    std::map<std::string,unsigned> zones_counts;
    unsigned currentDay;
    unsigned numDays;
    string trajectoryFile;
    string birthsFile;
    string configLine;
    string locationFile;
    string vaccineProfilesFile;
    string vaccineSettingsFile;
    string trialSettingsFile;
    string vaccinationStrategy;
    string reportsFile;
    string diseaseRatesFile;
    string vaccinationGroupsFile;
    string aegyptiRatesFile;

    bool routineVaccination;
    bool catchupFlag;
    bool trialVaccination;
    bool randomTrial;
    bool testBeforeVaccine;
    double routineTestSpecificity;
    double routineTestSensitivity;
    int vaccineID;
    unsigned vaccineDay;
    unsigned vaccineAge;
    double vaccineCoverage;
    Report outputReport;
    Recruitment recruitmentTrial;
    string outputFile;
    string outputPopFile;
    string outputPrevacFile;
    string simName;
    unsigned rSeed;
    unsigned rSeedInf;
    string outputPath;
    unsigned humanInfectionDays;
    unsigned huImm;
    double mlife;
    double mbite;
    double emergeFactor;
    double biteProbablity;
    double mozInfectiousness;   
    double mozMoveProbability;
    double attractShape;
    double selfReportProb;
    unsigned year;

    OutbreakResponse outbreakIntervention;
    string outbreakFile;
    
    int deathMoz;
    int lifeMoz;
    int humanDeaths;
    int visitorsCounter;
    vector<double> meanDailyEIP;
    vector<double> firstBiteRate;
    vector<double> secondBiteRate;
    vector<double> mozDailyDeathRate;
    vector<double> dailyEmergenceFactor;
    vector<double> InitialConditionsFOI;
    vector<map<unsigned, double>> dailyForceOfImportation;
    vector<map<unsigned, double>> dailyForceOfLocalImportation;

    std::ofstream out;
    std::ofstream outpop;
    std::ofstream outprevac;
};

#endif  /* SIMULATION_H */
