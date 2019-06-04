#include <fstream>
#include <string>
#include <iostream>
#include <utility>
#include <cstdlib>
#include <vector>
#include <algorithm>
#include <sstream>
#include "OutbreakResponse.h"
#include "Report.h"

OutbreakResponse::OutbreakResponse(){
    surveillanceEffort = 0.0;        
    responseThreshold = 0.0;
    aggressiveness = 0.0;
    thoroughness = 0.0;
    spatialRadius = 0.0;
    residuality = 0.0;
    compliance = 0.0;
    
    symptomaticCases = 0;
    totalHousesSprayed = 0;
    surveillanceDelay = 0;
    maxHouses = 0;
    maxHousesPerDay = 0;
    startDay = 0;
    endDay = 0;

    responseStrategy = "none";
    activeOutbreakResponse = false;
}

void OutbreakResponse::setup(string file){
    if (file.length() == 0) {
	exit(1);
    }
    string line;
    std::ifstream infile(file);
    if(!infile.good()){
	exit(1);
    }
    while(getline(infile,line,'\n')){
	this->addParameter(line);
    }
    infile.close();
    
    surveillanceEffort = this->readParameter("outbreak_surveillance_effort",0.0);
    responseThreshold = this->readParameter("outbreak_response_threshold", 0.0);
    aggressiveness = this->readParameter("outbreak_aggressiveness", 0.0);    
    thoroughness = this->readParameter("outbreak_thoroughness",0.0);
    spatialRadius = this->readParameter("outbreak_spatial_radius",0.0);
    residuality = this->readParameter("outbreak_residuality",0.0);
    compliance = this->readParameter("outbreak_compliance",0.0);
    
    surveillanceDelay = this->readParameter("outbreak_surveillance_delay",0);
    maxHouses = this->readParameter("outbreak_max_houses",0);
    maxHousesPerDay = this->readParameter("outbreak_max_houses_per_day",0);
    startDay = this->readParameter("outbreak_start_day",0);
    endDay = this->readParameter("outbreak_end_day",0);

    this->readParameter("outbreak_response_strategy",&responseStrategy);
    
    printf("surveillance effort: %.2f | response threshold: %.2f | aggressiveness: %.2f | thoroughness: %.2f | spatial radius: %.2f\n", surveillanceEffort, responseThreshold, aggressiveness, thoroughness, spatialRadius);
    printf("residuality: %.2f | compliance: %.2f | surveillance delay: %u | max houses: %u | max houses per day: %u\n", residuality, compliance, surveillanceDelay, maxHouses,maxHousesPerDay);
    printf("start day: %u | end day: %u\n", startDay, endDay);
    
}

void OutbreakResponse::setLocations(map < string,std::unique_ptr<Location> > * locations_in){
    locations_ptr = locations_in;
}

bool OutbreakResponse::isOutbreakResponseActive(unsigned currDay){
    if(responseStrategy == "ring" && currDay >= startDay && currDay <= endDay && totalHousesSprayed < maxHouses){
	return(true);
    }else{
	return(false);
    }
}

bool OutbreakResponse::isOutbreakResponseEnabled(){
    if(responseStrategy == "ring"){
	return(true);
    }else{
	return(false);
    }
}

void OutbreakResponse::addSymptomatic(unsigned currDay, Human * h, RandomNumGenerator * rGen, Report * r){
    if(rGen->getEventProbability() < surveillanceEffort){	
	symptomaticCases++;
	todaySympCounts++;
	if(symptomaticCases >= responseThreshold){	    
	    if(rGen->getEventProbability() < aggressiveness){
		r->addOutbreakSymptomatic(currDay,h,true);
		if(this->isOutbreakResponseActive(currDay)){
		    todaySymptomatics.push_back(h);
		}
	    }else{
		r->addOutbreakSymptomatic(currDay,h,false);
	    }
	}else{
	    r->addOutbreakSymptomatic(currDay,h,false);
	}
    }
}

void OutbreakResponse::update(unsigned currDay, RandomNumGenerator * rGen, Report * r){
    todaySympCounts = 0;
    if(this->isOutbreakResponseActive(currDay)){
	printf("Updating Outbreakresponse is active day %d records: %lu symptomatics today: %lu\n", currDay, symptomaticRecords.size(), todaySymptomatics.size());
	// Check symptomatic records and get the delayed records
	if(symptomaticRecords.size() == surveillanceDelay){
	    delayedSymptomatics = symptomaticRecords.front();
	    symptomaticRecords.erase(symptomaticRecords.begin());
	}
	symptomaticRecords.push_back(todaySymptomatics);
	todaySymptomatics.clear();

	if(!delayedSymptomatics.empty()){
	    printf("Day: %d There are %lu symptomatics delayed\n",currDay, delayedSymptomatics.size());
	    for(auto it_h:delayedSymptomatics){
		string house_id_tmp = (*it_h).getHouseID();
		Location * house_tmp = (*locations_ptr)[house_id_tmp].get();
		size_t n_nearby_locs = house_tmp->getNumberRadiusLocations();
		vector<string> nearby_locs = house_tmp->getRadiusLocations();
		printf("Day: %d house: %s No Locations: %lu\n", currDay,house_id_tmp.c_str(), n_nearby_locs);
		if(rGen->getEventProbability() < compliance){
		    if(housesEnlisted.empty()){		    
			housesEnlisted.insert(house_id_tmp);
			housesToSpray.push(house_id_tmp);
		    }else if(housesEnlisted.find(house_id_tmp) == housesEnlisted.end()){
			housesEnlisted.insert(house_id_tmp);
			housesToSpray.push(house_id_tmp);
		    }
		}
		for(auto it_loc: nearby_locs){
		    if(rGen->getEventProbability() < compliance){
			if(housesEnlisted.find(it_loc) == housesEnlisted.end()){
			    housesEnlisted.insert(it_loc);
			    housesToSpray.push(it_loc);
			}
		    }
		}		
	    }
	}
	delayedSymptomatics.clear();
	printf("Day %u size of houses enlisted %lu\n", currDay, housesEnlisted.size());
	printf("Day %u size of houses to spray %lu\n", currDay, housesToSpray.size());
	this->sprayTodayLocations(currDay, r);
    }
    if(totalHousesSprayed >= maxHouses || currDay == endDay){
	symptomaticRecords.clear();
	housesEnlisted.clear();
	while(!housesToSpray.empty()){
	    housesToSpray.pop();
	}
    }
}

void OutbreakResponse::sprayTodayLocations(unsigned currDay, Report * r){
    unsigned houses_sprayed_today = 0;
    while(houses_sprayed_today < maxHousesPerDay && totalHousesSprayed < maxHouses && !housesToSpray.empty()){
	string loc_id_tmp = housesToSpray.front();
	housesToSpray.pop();
	Location * loc_tmp = (*locations_ptr)[loc_id_tmp].get();	
	loc_tmp->sprayAdultInsecticide(currDay, thoroughness, residuality);
	houses_sprayed_today++;
	totalHousesSprayed++;
	r->addOutbreakResponseLocations(currDay, loc_tmp);
    }
    printf("Day %u houses sprayed %u total sprayed %u houses to spray %lu\n", currDay, houses_sprayed_today, totalHousesSprayed,housesToSpray.size());
}

// Parameters parsing
void OutbreakResponse::addParameter(string line){
    if(line.size() > 0 && line[0] != '#' && line[0] != ' '){
	string param_name, param_value;
	size_t pos_equal = line.find_first_of('=');
	if(pos_equal != string::npos){
	    param_name = line.substr(0,pos_equal);
	    param_value = line.substr(pos_equal + 1);	    
	    // trim trailing spaces and weird stuff for param_name
	    pos_equal = param_name.find_first_of(" \t");
	    if(pos_equal != string::npos){
		param_name = param_name.substr(0,pos_equal);
	    }
	    // trim trailing and leading spaces and weird stuff from param_value
	    pos_equal = param_value.find_first_not_of(" \t");
	    if(pos_equal != string::npos){
		param_value = param_value.substr(pos_equal);
	    }
	    pos_equal = param_value.find_first_of("#");
	    if(pos_equal != string::npos){
		param_value = param_value.substr(0,pos_equal);
	    }
	    // Add the parameter name and value to the map
	    parameters.insert(make_pair(param_name,param_value));
	}
    }
}

int OutbreakResponse::readParameter(string param_name, int vtemp){
    map<string, string>::iterator it;
    int values_ = vtemp;
    it = parameters.find(param_name);
    if(it != parameters.end()){
	values_ = this->parseInteger(it->second);
    }
    return values_;
}
double OutbreakResponse::readParameter(string param_name, double vtemp){
    map<string, string>::iterator it;
    double values_ = vtemp;
    it = parameters.find(param_name);
    if(it != parameters.end()){
	values_ = this->parseDouble(it->second);
    }
    return values_;
}

bool OutbreakResponse::readParameter(string param_name, bool vtemp){
    map<string, string>::iterator it;
    bool values_ = vtemp;
    it = parameters.find(param_name);
    if(it != parameters.end()){
	values_ = this-parseBoolean(it->second);
    }
    return values_;
}

bool OutbreakResponse::parseBoolean(string line){
    bool flag_temp = (std::stoi(line.c_str(), NULL, 10) == 0 ? false : true);
    return flag_temp;
}

void OutbreakResponse::readParameter(string param_name, string * param_var){
   map<string, string>::iterator it;
    it = parameters.find(param_name);
    if(it != parameters.end()){
	*param_var = this->parseString(it->second);
    }
}

string OutbreakResponse::parseString(string line){
    size_t first_ = line.find_first_not_of(' ');
    size_t last_ = line.find_last_not_of(' ');
    return line.substr(first_,(last_ - first_ + 1));
}


int OutbreakResponse::parseInteger(string line){
    return strtol(line.c_str(), NULL, 10);
}

double OutbreakResponse::parseDouble(string line){
    return strtod(line.c_str(), NULL);
}
