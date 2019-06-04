#ifndef REPORT_H
#define	REPORT_H

#include <fstream>
#include <string>
#include <vector>
#include <set>
#include <unordered_map>
#include <iostream>
#include <utility>
#include <memory>
#include <cmath>
#include "Human.h"
#include "Mosquito.h"
#include "RandomNumGenerator.h"

using std::string;
using std::vector;
using std::set;
using std::unordered_map;

class Report {    
 public:    
    struct rangeStruct{
	int min;
	int max;
    };
    
    struct eventStats{
	int events[5];
	int nonevents[5];
    };
    
    struct reportStats{
	eventStats status[4];
	eventStats total;
    };    
    
 private:
    
    string outputPopFile;   
    string outputCohortFile;
    string outputAgesFile;
    string outputFOIFile;

    std::ofstream outCohort;
    std::ofstream outAges;
    std::ofstream outPop;
    std::ofstream outFOI;
    std::ofstream outSpatial;

    bool reportHouseBites;
    bool reportCohort;
    bool reportAges;
    bool reportGroups;
    bool reportFOI;
    bool reportFOIOutbreakSymptomatics;
    bool printR0;
    bool reportSpatial;
    
    bool printAgesPop;
    bool printGroupsPop;
    bool printGroupsAgeFirst;
    bool printGroupsTotalAges;
    bool printZonesFOI;
    
    vector<string> events;
    vector<string> status;

    unordered_map<string, string> parameters;
    set<string> zonesToPrint;
    unordered_map<int,unordered_map<unsigned,unordered_map<string,int>>> secondaryCases;
    unordered_map<int,unordered_map<unsigned,unordered_map<string,int>>> dailyFOI;
    unordered_map<int,int> human_counts;
    unordered_map<int,int> dailyOutbreakResponseSymp;
    unordered_map<int,int> dailyOutbreakResponseLocations;
    unordered_map<int,unordered_map<string,unordered_map<string,unordered_map<string, unordered_map<string,int>>>>> cohortData;
    // <Cohort # <tested status < vaccine status < serostatus at vac. <output name, output value>>>>    
    unordered_map<string,unordered_map<string,unordered_map<string, unordered_map<string,int>>>> popData; // < vaccine status < serostatus at vac. <<age group<output name, output value>>>
    
    int cohortReportPeriod[3];

    int groupsEvents[5];
    int groupsStatus[2];
    int groupsReportPeriod[3];
    int groupsMaxIndex;
    int groupsAvgAgeFirst;
    int groupsTotalFirstInf;
    vector<rangeStruct> groupsAges;
    vector<reportStats> groupsStats;
    reportStats groupsTotalAgeStats;

    int falsePositives;
    int falseNegatives;
    int truePositives;
    int trueNegatives;
    int hospitalizations;
    int deathsFromInfection;
    int discountedLifeExpectancy;
    int cumAgeOfHospitalization;
    int totalVaccinated;
    int ageEvents[5];
    int ageStatus[4];
    rangeStruct discreteAges;
    int ageReportPeriod[3];
    int ageMaxIndex;
    vector<reportStats> ageStats;

    int foiReportPeriod[3];
    int foiTypes[4];
    int newInfections[4];
    int susceptibles[4];
    int susceptibles_temp[4];
    int importations[4];
    int mozSusceptibles[4];
    int mozExposed[4];
    int mozInfectious[4];    
    
    int spatialReportPeriod[3];
    bool spatialMosquitoes;
    bool spatialOutbreakReport;
    bool spatialSymptomatics;
    vector<string> spatialData;
	
    void parseEvents(string line, int *, int);
    void parseGroupsAges(string line, vector<rangeStruct> *);
    rangeStruct parseDiscreteAges(string);
    void parsePeriod(string line, int *);
    void addParameter(string);

    void readParameter(string,string,int *);
    void readParameter(string,string, vector<rangeStruct> *);
    rangeStruct readParameter(string,string, rangeStruct);
    bool readParameter(string, bool);
    bool parseBoolean(string line);
    
    int getGroup(int, vector<rangeStruct>);
    int calculateExpectedLifeInDays(int);

public:


    Report(string);
    Report();
    Report(const Report& orig);
    //virtual ~Report();

    void setupReport(string, string, string);
    void setupZones(set<string>);
    void updateReport(int, Human *);
    void join(const vector<string>&,char, string &);
    void printReport(int, int);
    void printHeaders();
    void printPopHeader();
    void printCohortHeader();
    void printAgesHeader();
    void printSpatialHeader();
    void printPopReport(int);
    void printCohortReport(int, int);
    void printAgesReport(int);
    void printFOIReport(int);
    void printSpatialReport(int);
    void updateCohortReport(int, Human *);
    void updatePopReport(int, Human *);
    void updatePossibleVaccineeStatus(bool, bool);
    void addVaccinated();
    void updateAgesReport(int, Human *);
    void updateSecondaryCases(int, Human *);
    void updateFOI(int, Human *);
    void updateOutbreakSymptomatics(int, int);
    void addImportation(int, int, Human *);
    void updateSpatialReport(int, Human *);
    void updateMosquitoReport(int, Mosquito *);
    void addOutbreakSymptomatic(int , Human *, bool );
    void addOutbreakResponseLocations(int , Location *);
    void resetReports();
    void resetGroupStats();
    void resetCohortStats();
    void resetAgeStats();
    void resetFOIStats();
    void resetSpatialStats();    
    void finalizeReport(int);
    bool isReportHouseBitesEnabled(){return reportHouseBites;}

};

#endif	/* REPORT_H */
