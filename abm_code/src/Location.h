#ifndef LOCATION_H
#define	LOCATION_H

#include <string>
#include <vector>
#include <set>
#include <forward_list>
#include <memory>
#include <unordered_map>
#include <fstream>
#include "RandomNumGenerator.h"
#include "Human.h"

using std::string;
using std::vector;

class Location {
private:
    string locID;
    string neighborhoodID;
    string zoneID;
    double xCor;
    double yCor;
    double emergenceRate;
    bool infectedVisitor;
    string locType;
    bool bitesCounterEnabled;
    double insecticideDecayRate;
    double insecticideEfficacy;
    unsigned insecticideStartDay;
    bool insecticideSprayed;
    
    std::unique_ptr<vector<string>> closeLocs;
    std::unique_ptr<vector<string>> radiusLocs;
    std::set<sp_human_t,Human::sortid> humans;
    std::unordered_map<string,int>visitorBites;    
public:
    string getRandomCloseLoc(RandomNumGenerator&);
    size_t getNumberRadiusLocations(){return radiusLocs->size();}
    void addHuman(sp_human_t);
    void increaseBites(string);
    void enableBitesCounter(){bitesCounterEnabled = true;}
    void disableBitesCounter(){bitesCounterEnabled = false;}
    void sprayAdultInsecticide(unsigned,double,double);
    
    double getIncreasedMortalityInsecticide(unsigned, double);
    double calculateGiniIndex();
    void removeHuman(sp_human_t);
    std::set<sp_human_t,Human::sortid> & getHumans(){return humans;}
    void addCloseLoc(string);
    void addRadiusLoc(string);
    void printHumans();
    double getDistanceFromLoc(Location &) const;
    double getLocX() const;
    double getLocY() const;
    double getEmergenceRate() const;
    bool getInfectedVisitor(){return infectedVisitor;}
    string getLocID() const;
    string getLocType() const;
    string getNeighID() const;
    string getZoneID() const;
    vector<string> getRadiusLocations();
    
    Location(string, string, double, double, double);
    Location(string, string, string, string, double, double, double);
    Location();
    Location(const Location& orig);
    //virtual ~Location();
    void updateInfectedVisitor();

private:

};

#endif	/* LOCATION_H */

