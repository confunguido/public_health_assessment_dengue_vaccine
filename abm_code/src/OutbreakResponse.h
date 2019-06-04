#ifndef OUTBREAKRESPONSE_H
#define	OUTBREAKRESPONSE_H

#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <queue>
#include <unordered_set>
#include <iostream>
#include <utility>
#include <algorithm>
#include <memory>
#include <cmath>
#include "RandomNumGenerator.h"
#include "Human.h"
#include "Location.h"
#include "Report.h"

using std::string;
using std::queue;
using std::map;
using std::unordered_set;
using std::vector;

using symptomatic_t = vector<Human *>;

class OutbreakResponse {
 public:
    OutbreakResponse();
    OutbreakResponse(const OutbreakResponse& orig);

    void setup(string file);
    void setLocations(map <string,std::unique_ptr<Location>> *);
    void update(unsigned currDay, RandomNumGenerator * rGen, Report * r);
    void addSymptomatic(unsigned currDay, Human * h, RandomNumGenerator * rGen, Report * r);

    bool isOutbreakResponseActive(unsigned currDay);
    bool isOutbreakResponseEnabled();

    double getResponseRadius(){return spatialRadius;}
    unsigned getTodaySymptomatics(){return todaySympCounts;}
    
    void finalizeOutbreakResponse();
    
 private:
    
    double surveillanceEffort;        
    double responseThreshold;
    double aggressiveness;
    double thoroughness;
    double spatialRadius;
    double residuality;
    double compliance;
    
    unsigned surveillanceDelay;
    unsigned symptomaticCases;
    unsigned maxHouses;
    unsigned maxHousesPerDay;
    unsigned startDay;
    unsigned endDay;
    unsigned totalHousesSprayed;
    unsigned todaySympCounts;
    
    string responseStrategy;
    bool activeOutbreakResponse;
        
    symptomatic_t todaySymptomatics;
    symptomatic_t delayedSymptomatics;

    vector<symptomatic_t> symptomaticRecords;
    queue<string> housesToSpray;
    unordered_set<string> housesEnlisted;   
    map < string,std::unique_ptr<Location> > * locations_ptr;    

    map<string, string> parameters;

    void sprayTodayLocations(unsigned currDay, Report * r);
    
    // Routines to read parameters    
    void addParameter(string);
    int parseInteger(string);
    int readParameter(string, int );
    bool parseBoolean(string);
    double parseDouble(string);
    double readParameter(string, double);
    bool readParameter(string , bool);
    void readParameter(string , string *);
    string parseString(string);
};

#endif	/* OUTBREAKRESPONSE_H */
