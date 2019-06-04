#include "Location.h"
#include <sstream>
#include <iostream>
#include <cmath>

string Location::getLocID() const {
    return locID;
}

string Location::getLocType() const {
    return locType;
}

string Location::getNeighID() const {
    return neighborhoodID;
}

string Location::getZoneID() const{
  return zoneID;
}

double Location::getLocX() const {
    return xCor;
}

double Location::getLocY() const {
    return yCor;
}

double Location::getEmergenceRate() const {
    return emergenceRate;
}

Location::Location(string lID, string lType, double x, double y, double e) {
    locID = lID;
    locType = lType;
    xCor = x;
    yCor = y;
    emergenceRate = e;
    closeLocs.reset(new vector<string>());
    radiusLocs.reset(new vector<string>());
    infectedVisitor = false;
    visitorBites.clear();
    bitesCounterEnabled = false;
    insecticideDecayRate = 1.0;
    insecticideEfficacy = 0.0;
    insecticideStartDay = 0;
    insecticideSprayed = false;
}

Location::Location(string lID, string lType, string NID, string ZID, double x, double y, double e) {
    locID = lID;
    locType = lType;
    neighborhoodID = NID;
    zoneID = ZID;
    xCor = x;
    yCor = y;
    emergenceRate = e;
    closeLocs.reset(new vector<string>());
    radiusLocs.reset(new vector<string>());
    infectedVisitor = false;
    visitorBites.clear();
    bitesCounterEnabled = false;
    insecticideDecayRate = 1.0;
    insecticideEfficacy = 0.0;
    insecticideStartDay = 0;
    insecticideSprayed = false;
}

string Location::getRandomCloseLoc(RandomNumGenerator& rGen) {
    int i = closeLocs->size();
    if (i > 0)
        return (*closeLocs)[rGen.getMozNextLoc(i)];
    else return "TOO_FAR_FROM_ANYWHERE";
}

double Location::getDistanceFromLoc(Location& loc) const {
    return sqrt((xCor - loc.getLocX()) * (xCor - loc.getLocX()) + (yCor - loc.getLocY()) * (yCor - loc.getLocY()));
}

vector<string> Location::getRadiusLocations(){
    return(*(radiusLocs.get()));
}

void Location::addCloseLoc(string loc) {
    closeLocs->push_back(loc);
}

void Location::addRadiusLoc(string loc) {
    radiusLocs->push_back(loc);
}

void Location::sprayAdultInsecticide(unsigned currDay,double efficacy_,double residuality_){
    insecticideStartDay = currDay;
    insecticideSprayed = true;
    insecticideEfficacy = efficacy_;
    insecticideDecayRate = residuality_;
}

double Location::getIncreasedMortalityInsecticide(unsigned currDay, double mozMortality_){
    if(insecticideSprayed){
	double mortality_rate = -log(1 - mozMortality_);
	mortality_rate += insecticideEfficacy * exp(-insecticideDecayRate * (currDay - insecticideStartDay));
	return(1-exp(-mortality_rate));
    }else{
	return(mozMortality_);
    }
}

void Location::addHuman(sp_human_t h) {
    humans.insert(h);
}

void Location::removeHuman(sp_human_t h){
    // set.erase(val) requires no check
    humans.erase(h);
}

//Location::Location() {
//}

//Location::Location(const Location& orig) {
//}

//Location::~Location() {
//}



void Location::updateInfectedVisitor(){
    infectedVisitor = false;
    for(auto itHum = humans.begin(); itHum != humans.end(); itHum++){
        if((*itHum)->infection != nullptr){
            infectedVisitor = true;
            return;
        }
    }
}

void Location::printHumans(){
    for(auto itHum = humans.begin(); itHum != humans.end(); itHum++){
	printf("Human %s in Location %s\n", (*itHum)->getPersonID().c_str(), locID.c_str());
    }
}



void Location::increaseBites(string personID){
    if(bitesCounterEnabled == true){
	visitorBites[personID]++;
    }
}

double Location::calculateGiniIndex(){
    double num = 0;
    double den = 0;
    for(auto it = visitorBites.begin(); it != visitorBites.end();){
	den += (*it).second;
	for(auto jt = visitorBites.begin(); jt != visitorBites.end();){
	    num += std::abs((*it).second - (*jt).second);
	    ++jt;
	}
	++it;
    }
    if(den > 0){
	return(num / (2.0 * visitorBites.size() * den));
    }else{
	return(-1.0);
    }
}
