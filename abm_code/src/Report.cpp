#include <fstream>
#include <string>
#include <iostream>
#include <utility>
#include <cstdlib>
#include <vector>
#include <algorithm>
#include <sstream>
#include "Report.h"
#include <cmath>

using std::stringstream;

Report::Report(){
    reportHouseBites = false;
    reportCohort = false;
    reportAges = false;
    reportGroups = false;
    reportFOI = false;
    reportFOIOutbreakSymptomatics = false;    
    printR0 = false;
    reportSpatial = false;    
    printAgesPop = false;
    printGroupsPop = false;
    printZonesFOI = false;    
    printGroupsAgeFirst = false;
    printGroupsTotalAges = false;
    spatialMosquitoes = false;
    spatialOutbreakReport = false;
    spatialSymptomatics = false;
    groupsMaxIndex = 0;
    ageMaxIndex = 0;
    groupsAges.clear();
    groupsStats.clear();
    groupsAvgAgeFirst = 0;
    secondaryCases.clear();
    dailyFOI.clear();	
    discreteAges.min = 0;
    discreteAges.max = 0;
    truePositives = 0;
    falsePositives = 0;
    falseNegatives = 0;
    trueNegatives = 0;
    totalVaccinated = 0;
    hospitalizations = 0;
    deathsFromInfection = 0;
    discountedLifeExpectancy = 0;
    cumAgeOfHospitalization = 0;
    ageStats.clear();
    parameters.clear();
    spatialData.clear();
    cohortData.clear();
    popData.clear();
    
    for(int i = 0;i < 5;i++){
	groupsEvents[i] = 0;
	ageEvents[i] = 0;
    }    
    for(int i = 0; i < 2; i++){
	groupsStatus[i] = 0;
	ageStatus[i] = 0;
    }    
    for(int i = 0; i < 3; i++){
	groupsReportPeriod[i] = 0;
	cohortReportPeriod[i] = 0;
	ageReportPeriod[i] = 0;
	foiReportPeriod[i] = 0;
	spatialReportPeriod[i] = 0;
    }    
    for(int i = 0; i < 4; i++){
	foiTypes[i] = 0;
	newInfections[i] = 0;
	susceptibles[i] = 0;
	susceptibles_temp[i] = 0;
	mozSusceptibles[i] = 0;
	mozExposed[i] = 0;
	mozInfectious[i] = 0;
	importations[i] = 0;
    }    
    for(int i = 0; i < 5; i++){
	// There are four status VacSero+ VacSero- PlacSero+ PlacSero-
	for(int j = 0; j < 4; j++){
	    groupsTotalAgeStats.status[j].events[i] = 0;
	    groupsTotalAgeStats.status[j].nonevents[i] = 0;
	}
	groupsTotalAgeStats.total.events[i] = 0;
	groupsTotalAgeStats.total.nonevents[i] = 0;
    }
    
    outCohort.close();
    outAges.close();
    outPop.close();
    events = {"inf", "dis", "hosp", "serop", "vac"};
    status = {"vac", "plac", "serop", "seron"};
}

void Report::setupZones(set<string> zonesIn){
    zonesToPrint = zonesIn;
    for(auto locIt = zonesToPrint.begin(); locIt != zonesToPrint.end();){
	std::string tmpstr = (*locIt);
	printf("Zone %s in Report\n", tmpstr.c_str());
	++locIt;
    }
}

void Report::setupReport(string file, string outputPath_, string simName_) {    
    //test the life expectancy
    int age_tmp = calculateExpectedLifeInDays(41*365);
    int age_child = calculateExpectedLifeInDays(10*365);
    int age_old = calculateExpectedLifeInDays(80*365);
    printf("LIFE EXPECTANCY TEST:: at 41: %d, at 10: %d, at 80: %d\n", age_tmp, age_child,age_old);
    
    //       Read the reporting file and assign variable values
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

    this->readParameter("groups_events", "events",groupsEvents);
    this->readParameter("age_events","events", ageEvents);
    this->readParameter("groups_status","status",groupsStatus);
    this->readParameter("age_status","status",ageStatus);
    this->readParameter("groups_report_period","period",groupsReportPeriod);
    this->readParameter("cohort_report_period","period",cohortReportPeriod);
    this->readParameter("age_report_period","period",ageReportPeriod);
    this->readParameter("foi_report_period","period",foiReportPeriod);
    this->readParameter("foi_serotypes","serotypes",foiTypes);
    this->readParameter("spatial_report_period", "period", spatialReportPeriod);

    this->readParameter("groups_ages","ages",&groupsAges);

    discreteAges = this->readParameter("age_ages","ages", discreteAges);
    printGroupsPop = this->readParameter("groups_complement", printGroupsPop);
    printAgesPop = this->readParameter("age_complement", printAgesPop);
    printGroupsAgeFirst = this->readParameter("groups_avg_first", printGroupsAgeFirst);
    printGroupsTotalAges = this->readParameter("groups_print_total_ages", printGroupsTotalAges);
    printZonesFOI = this->readParameter("foi_print_zones", printZonesFOI);
    reportGroups = this->readParameter("groups_print", reportGroups);
    reportCohort = this->readParameter("cohort_print", reportCohort);
    reportHouseBites = this->readParameter("housebites_print", reportHouseBites);
    reportAges = this->readParameter("age_print", reportAges);
    reportFOI = this->readParameter("foi_print", reportFOI);
    reportFOIOutbreakSymptomatics = this->readParameter("foi_outbreak_symptomatics", reportFOIOutbreakSymptomatics);
    printR0 = this->readParameter("foi_print_R0", printR0);
    reportSpatial = this->readParameter("spatial_print", reportSpatial);
    spatialMosquitoes = this->readParameter("spatial_mosquitoes", spatialMosquitoes);
    spatialSymptomatics = this->readParameter("spatial_symptomatics", spatialSymptomatics);
    spatialOutbreakReport = this->readParameter("spatial_outbreak_report", spatialOutbreakReport);

    if(reportFOI == true){
	outputFOIFile = outputPath_ + "/" + simName_ + "_foi.csv";
    }
    if(reportSpatial == true){
	string outputSpatialFile = outputPath_ + "/" + simName_ + "_spatial.csv";
	outSpatial.open(outputSpatialFile);
	if (!outSpatial.good()) {
	    exit(1);
	}
    }
    if(reportGroups == true){
	outputPopFile = outputPath_ + "/" + simName_ + "_pop.csv";
	outPop.open(outputPopFile);
	if (!outPop.good()) {
	    exit(1);
	}
    }
    if(reportCohort == true){
	outputCohortFile = outputPath_ + "/" + simName_ + "_cohort.csv";
	outCohort.open(outputCohortFile);
	if (!outCohort.good()) {
	    exit(1);
	}
    }
    if(reportAges == true){
	outputAgesFile = outputPath_ + "/" + simName_ + "_ages.csv";
	outAges.open(outputAgesFile);
	if (!outAges.good()) {
	    exit(1);
	}
    }
    resetReports();
    printHeaders();
}

void Report::readParameter(string param_name,string param_type, vector<rangeStruct> * values_){
    unordered_map<string, string>::iterator it;
    it = parameters.find(param_name);
    if(it != parameters.end()){
	if(param_type == "ages"){
	    if(param_name != "age_ages"){
		parseGroupsAges(it->second,values_);
	    }
	}
    }
}

Report::rangeStruct Report::readParameter(string param_name,string param_type, rangeStruct vtemp){
   unordered_map<string, string>::iterator it;
   rangeStruct values_ = vtemp;
    it = parameters.find(param_name);
    if(it != parameters.end()){
	if(param_type == "ages" && param_name == "age_ages"){
	    values_ = parseDiscreteAges(it->second);
	}
    }
    return values_;
}

bool Report::readParameter(string param_name, bool vtemp){
    unordered_map<string, string>::iterator it;
    bool values_ = vtemp;
    it = parameters.find(param_name);
    if(it != parameters.end()){
	values_ = parseBoolean(it->second);
    }
    return values_;
}

void Report::readParameter(string param_name, string param_type, int * values_){
    unordered_map<string, string>::iterator it;
    it = parameters.find(param_name);
    if(it != parameters.end()){
	if(param_type == "events"){
	    parseEvents(it->second,values_,5);
	}
	if(param_type == "status"){
	    parseEvents(it->second,values_,2);
	}
	if(param_type == "period"){
	    parsePeriod(it->second,values_);
	}
	if(param_type == "serotypes"){
	    parseEvents(it->second,values_,4);
	}
    }
}

void Report::addParameter(string line){
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

bool Report::parseBoolean(string line){
    bool flag_temp = (std::stoi(line.c_str(), NULL, 10) == 0 ? false : true);
    return flag_temp;
}

void Report::parsePeriod(string line, int * period_temp){
    stringstream linetemp;
    string line2;
    linetemp.clear();
    int count =0;
    linetemp << line;
    while(getline(linetemp,line2,',')){
		if(count > 2){
		    break;
		}
		period_temp[count] = strtol(line2.c_str(), NULL, 10);
		if(period_temp[count] < 0){
		    period_temp[count] = 0;
		}
		count++;
    }
    // If there are less than 3 values in the period, or the the start:increase:end don't make sense
    if(count < 3 || !(period_temp[0] < period_temp[1] + period_temp[2] && period_temp[1] < period_temp[2])){
	printf("Reporting issue: The period has not been setup properly\n");
	exit(1);
    }
}

void Report::parseGroupsAges(string line, vector<rangeStruct> * ages_temp){
    stringstream linetemp;
    string line2;
    linetemp.clear();
    linetemp << line;
    ages_temp->clear();
    while(getline(linetemp,line2,';')){
	stringstream lTemp; lTemp << line2;
	string line3;
	rangeStruct rangeTemp;
	getline(lTemp,line3,',');
	rangeTemp.min = strtol(line3.c_str(), NULL, 10);
	getline(lTemp,line3,',');
	rangeTemp.max = strtol(line3.c_str(), NULL, 10);
	if(rangeTemp.min + rangeTemp.max > 0){
	    ages_temp->push_back(rangeTemp);
	}
    }
    if(ages_temp->empty()){
	printf("Report::ParseGroupsAges ages_temp is empty\n");
	exit(1);
    }
}

Report::rangeStruct Report::parseDiscreteAges(string line){
    stringstream linetemp;
    string line2;
    linetemp.clear();
    linetemp << line;
    rangeStruct rangeTemp;
    getline(linetemp,line2,',');
    rangeTemp.min = strtol(line2.c_str(), NULL, 10);
    getline(linetemp,line2,',');
    rangeTemp.max = strtol(line2.c_str(), NULL, 10);
    if(rangeTemp.min + rangeTemp.max > 0 && rangeTemp.min < rangeTemp.max && rangeTemp.min >= 0){
	return rangeTemp;
    }else{
	printf("Report::parseDiscreteAges minimum and maximum ages do not make sense\n");
	exit(1);
    }
}

void Report::parseEvents(string line, int * Events_, int len){
    stringstream linetemp;
    string line2;
    linetemp.clear();
    int count =0;
    linetemp << line;
    while(getline(linetemp,line2,',')){
	if(count >= len){
	    break;
	}
	Events_[count] = strtol(line2.c_str(), NULL, 10);
	if(Events_[count] > 1){
	    Events_[count] = 1;
	}
	if(Events_[count] < 0){
	    Events_[count] = 0;
	}
	count++;
    }
}

void Report::updateMosquitoReport(int currDay, Mosquito * m){
    if(reportFOI == true){
	if(m->infection == nullptr){
	    for(unsigned i = 1; i < 5; i++){
		dailyFOI[currDay][i]["mozsus"]++;
	    }
	}else{
	    unsigned sero = m->infection->getInfectionType();
	    if(m->infection->getInfectiousness() > 0.0){
		dailyFOI[currDay][sero]["mozinf"]++;
	    }else{
		dailyFOI[currDay][sero]["mozexp"]++;
	    }
	}
    }
    if(reportSpatial == true && spatialMosquitoes == true){
	if(currDay >= spatialReportPeriod[0] && currDay <= spatialReportPeriod[2] && (currDay - spatialReportPeriod[0]) % spatialReportPeriod[1] == 0){
	    string tmp_str = "Mosquitoes," + m->getLocationID()  + "," + std::to_string(currDay) + "," + "NaN";
	    spatialData.push_back(tmp_str);
	}
    }
}

void Report::updateReport(int currDay, Human * h){
    int reportNum = 0;
    if(reportGroups == true){
	if(currDay >= groupsReportPeriod[0] && currDay <= groupsReportPeriod[2]){	    
	    updatePopReport(currDay, h);	    
	    if((currDay - groupsReportPeriod[0]) % groupsReportPeriod[1] == 0){
		h->setPopCount(false);
		h->setSeroprevalenceCount(false);
	    }
	}
    }
    if(reportCohort == true){
	if(currDay >= cohortReportPeriod[0] && currDay <= cohortReportPeriod[2]){	    
	    updateCohortReport(currDay, h);
	    if((currDay - cohortReportPeriod[0]) % cohortReportPeriod[1] == 0){
		h->setCohortCount(false);
	    }
	}
    }
    if(reportAges == true){
	if(currDay >= ageReportPeriod[0] && currDay <= ageReportPeriod[2] && (currDay - ageReportPeriod[0]) % ageReportPeriod[1] == 0){
	    updateAgesReport(currDay, h);
	    reportNum++;
	}
    }
    if(reportFOI == true){
	if(printR0 == true){
	    updateSecondaryCases(currDay,h);
	}
	updateFOI(currDay, h);
    }
    if(reportSpatial == true){
	if(currDay >= spatialReportPeriod[0] && currDay <= spatialReportPeriod[2] && (currDay - spatialReportPeriod[0]) % spatialReportPeriod[1] == 0){
	    updateSpatialReport(currDay, h);
	}
    }
    if(reportNum > 0){
	h->resetRecent();
    }
}

void Report::printReport(int currDay, int vDay_){
    if(reportGroups == true){
	if(currDay >= groupsReportPeriod[0] && currDay <= groupsReportPeriod[2] && (currDay - groupsReportPeriod[0]) % groupsReportPeriod[1] == 0){
	    printPopReport(currDay);
	    resetGroupStats();
	}
    }
    if(reportCohort == true){

	if(currDay >= cohortReportPeriod[0] && currDay <= cohortReportPeriod[2] && (currDay - cohortReportPeriod[0]) % cohortReportPeriod[1] == 0){
	    printCohortReport(currDay, vDay_);
	    resetCohortStats();	    
	}
    }
    if(reportAges == true){
	if(currDay >= ageReportPeriod[0] && currDay <= ageReportPeriod[2] && (currDay - ageReportPeriod[0]) % ageReportPeriod[1] == 0){
	    printAgesReport(currDay);
	    resetAgeStats();
	}
    }
    if(reportSpatial == true){
	if(currDay >= spatialReportPeriod[0] && currDay <= spatialReportPeriod[2] && (currDay - spatialReportPeriod[0]) % spatialReportPeriod[1] == 0){
	    printSpatialReport(currDay);
	    resetSpatialStats();
	}
    }
}

void Report::updateSpatialReport(int currDay, Human * h){
    if(h->infection != nullptr){
	unsigned sero = h->infection->getInfectionType();
	if(currDay == h->infection->getStartDay()){
	    string tmp_str = "InfectiousHum," + h->getHouseID() + "," + std::to_string(currDay) + "," + std::to_string(sero);
	    spatialData.push_back(tmp_str);
	}
	if(h->isSymptomatic() && currDay == floor(h->infection->getSymptomOnset())){
	    string tmp_str = "SymptomaticHum," + h->getHouseID() + "," + std::to_string(currDay) + "," + std::to_string(sero);
	    spatialData.push_back(tmp_str);
	}
    }
}
void Report::addOutbreakResponseLocations(int currDay, Location * locNow){
    dailyOutbreakResponseLocations[currDay]++;
    if(reportSpatial == true && spatialOutbreakReport){
	string tmp_str = "OutbreakResponseLocation," + locNow->getLocID() + "," + std::to_string(currDay) + ",NaN";
	spatialData.push_back(tmp_str);
    }
}
void Report::addOutbreakSymptomatic(int currDay, Human * h, bool respond_){
    dailyOutbreakResponseSymp[currDay]++;
    if(reportSpatial == true && spatialOutbreakReport){
	string tmp_str = "OutbreakReportedHum," + h->getHouseID() + "," + std::to_string(currDay) + ",NaN";
	spatialData.push_back(tmp_str);
	tmp_str = "OutbreakReportedHumInfLoc," + h->getInfLocID() + "," + std::to_string(currDay) + ",NaN";
	spatialData.push_back(tmp_str);
	if(respond_){
	    tmp_str = "OutbreakRespondToHum," +  h->getHouseID() + "," + std::to_string(currDay) + ",NaN";
	    spatialData.push_back(tmp_str);
	    tmp_str = "OutbreakRespondToHumInfLoc," + h->getInfLocID() + "," + std::to_string(currDay) + ",NaN";
	    spatialData.push_back(tmp_str);
	}
    }
}

void Report::updateAgesReport(int currDay, Human * h){

    int age_temp = floor((double) h->getAgeDays(currDay) / 365.0);
    if(age_temp > discreteAges.max || age_temp < discreteAges.min){
	return;
    }
    int group = age_temp - discreteAges.min;

    if(h->getRecentInf() > 0){
	ageStats[group].total.events[0]++;
	if(h->isVaccinated() == true){
	    if(h->getSeroStatusAtVaccination() == true){
		ageStats[group].status[0].events[0]++;
	    }else{
		ageStats[group].status[1].events[0]++;
	    }
	}else{
	    if(h->getSeroStatusAtVaccination() == true){
		ageStats[group].status[2].events[0]++;
	    }else{
		ageStats[group].status[3].events[0]++;
	    }
	}
	if(h->getRecentDis() > 0){
	    ageStats[group].total.events[1]++;
	    if(h->isVaccinated() == true){
		if(h->getSeroStatusAtVaccination() == true){
		    ageStats[group].status[0].events[1]++;
		}else{
		    ageStats[group].status[1].events[1]++;
		}
	    }else{
		if(h->getSeroStatusAtVaccination() == true){
		    ageStats[group].status[2].events[1]++;
		}else{
		    ageStats[group].status[3].events[1]++;
		}
	    }
	    if(h->getRecentHosp() > 0){
		ageStats[group].total.events[2]++;
		if(h->isVaccinated() == true){
		    if(h->getSeroStatusAtVaccination() == true){
			ageStats[group].status[0].events[2]++;
		    }else{
			ageStats[group].status[1].events[2]++;
		    }
		}else{
		    if(h->getSeroStatusAtVaccination() == true){
			ageStats[group].status[2].events[2]++;
		    }else{
			ageStats[group].status[3].events[2]++;
		    }
		}
	    } else {
		ageStats[group].total.nonevents[2]++;
		if(h->isVaccinated() == true){
	    if(h->getSeroStatusAtVaccination() == true){
		ageStats[group].status[0].nonevents[2]++;
	    }else{
		ageStats[group].status[1].nonevents[2]++;
	    }
		}else{
		    if(h->getSeroStatusAtVaccination() == true){
			ageStats[group].status[2].nonevents[2]++;
		    }else{
			ageStats[group].status[3].nonevents[2]++;
		    }
		}
	    }
	} else {
	    ageStats[group].total.nonevents[1]++;
	    if(h->isVaccinated() == true){
		if(h->getSeroStatusAtVaccination() == true){
		    ageStats[group].status[0].nonevents[1]++;
		}else{
		    ageStats[group].status[1].nonevents[1]++;
		}
	    }else{
		if(h->getSeroStatusAtVaccination() == true){
		    ageStats[group].status[2].nonevents[1]++;
		}else{
		    ageStats[group].status[3].nonevents[1]++;
		}
	    }
	}
    } else {
	ageStats[group].total.nonevents[0]++;
	if(h->isVaccinated() == true){
	    if(h->getSeroStatusAtVaccination() == true){
		ageStats[group].status[0].nonevents[0]++;
	    }else{
		ageStats[group].status[1].nonevents[0]++;
	    }
	}else{
	    if(h->getSeroStatusAtVaccination() == true){
		ageStats[group].status[2].nonevents[0]++;
	    }else{
		ageStats[group].status[3].nonevents[0]++;
	    }
	}
    }

    if(h->getPreviousInfections() > 0){
	ageStats[group].total.events[3]++;
	if(h->isVaccinated() == true){
	    if(h->getSeroStatusAtVaccination() == true){
		ageStats[group].status[0].events[3]++;
	    }else{
		ageStats[group].status[1].events[3]++;
	    }
	}else{
	    if(h->getSeroStatusAtVaccination() == true){
		ageStats[group].status[2].events[3]++;
	    }else{
		ageStats[group].status[3].events[3]++;
	    }
	}
    } else {
	ageStats[group].total.nonevents[3]++;
	if(h->isVaccinated() == true){
	    if(h->getSeroStatusAtVaccination() == true){
		ageStats[group].status[0].nonevents[3]++;
	    }else{
		ageStats[group].status[1].nonevents[3]++;
	    }
	}else{
	    if(h->getSeroStatusAtVaccination() == true){
		ageStats[group].status[2].nonevents[3]++;
	    }else{
		ageStats[group].status[3].nonevents[3]++;
	    }
	}
    }
    if(h->isVaccinated() == true){
	ageStats[group].total.events[4]++;
	if(h->isVaccinated() == true){
	    if(h->getSeroStatusAtVaccination() == true){
		ageStats[group].status[0].events[4]++;
	    }else{
		ageStats[group].status[1].events[4]++;
	    }
	}else{
	    if(h->getSeroStatusAtVaccination() == true){
		ageStats[group].status[2].events[4]++;
	    }else{
		ageStats[group].status[3].events[4]++;
	    }
	}
    } else {
	ageStats[group].total.nonevents[4]++;
	if(h->isVaccinated() == true){
	    if(h->getSeroStatusAtVaccination() == true){
		ageStats[group].status[0].nonevents[4]++;
	    }else{
		ageStats[group].status[1].nonevents[4]++;
	    }
	}else{
	    if(h->getSeroStatusAtVaccination() == true){
		ageStats[group].status[2].nonevents[4]++;
	    }else{
		ageStats[group].status[3].nonevents[4]++;
	    }
	}
    }
}

void Report::updatePossibleVaccineeStatus(bool test_, bool serostatus_){
    if(test_ == true){
	if(serostatus_ == true){
	    popData["All"]["All"]["All"]["truePositives"]++;
	    truePositives++;
	}else{
	    popData["All"]["All"]["All"]["falsePositives"]++;
	    falsePositives++;
	}
    }else{
	if(serostatus_ == false){
	    popData["All"]["All"]["All"]["trueNegatives"]++;
	    trueNegatives++;
	}else{
	    popData["All"]["All"]["All"]["falseNegatives"]++;
	    falseNegatives++;
	}
    }
}

void Report::addVaccinated(){
    totalVaccinated++;
    popData["All"]["All"]["All"]["totalVaccinated"]++;
}

int Report::calculateExpectedLifeInDays(int current_age_in){
    double expected_life_days = 0.0;
    double S_previous = 1.0;
    for(int aa = 1; aa < 36501; aa++){
	double S_t = S_previous * (1.0 - (double)((6.9e-9) * (aa + current_age_in)));
	expected_life_days += aa * (S_previous - S_t);
	S_previous = S_t;
    }
    return(ceil(expected_life_days));
}

void Report::updatePopReport(int currDay, Human * h){    
    if(popData.empty()){
	vector<string> health_outcomes; health_outcomes.clear();
	health_outcomes.push_back("Size");
	health_outcomes.push_back("Inf");
	health_outcomes.push_back("Dis");
	health_outcomes.push_back("Hosp");
	health_outcomes.push_back("Seroprevalence");
	health_outcomes.push_back("deathsFromInfection");
	health_outcomes.push_back("firstInfections");
	health_outcomes.push_back("cumAgeOfInfection");
	health_outcomes.push_back("cumAgeOfHospitalization");
	health_outcomes.push_back("discountedLifeExpectancy");	
	health_outcomes.push_back("truePositives");	
	health_outcomes.push_back("trueNegatives");	
	health_outcomes.push_back("falsePositives");	
	health_outcomes.push_back("falseNegatives");	
	health_outcomes.push_back("totalVaccinated");	
	for(unsigned int i = 0; i < health_outcomes.size(); i++){
	    string health_str = health_outcomes[i];
	    for(unsigned k = 0; k < (groupsAges.size() + 1); k++){
		string ageGroupStr;
		if(k < groupsAges.size()){
		    ageGroupStr = std::to_string(groupsAges[k].min) + "_" + std::to_string(groupsAges[k].max);
		}else{
		    ageGroupStr = "All";
		}	    
		popData["Vaccine"]["Seropositive"][ageGroupStr][health_str] = 0;
		popData["Vaccine"]["Seronegative"][ageGroupStr][health_str] = 0;
		popData["Placebo"]["Seropositive"][ageGroupStr][health_str] = 0;
		popData["Placebo"]["Seronegative"][ageGroupStr][health_str] = 0;
		popData["All"]["All"][ageGroupStr][health_str] = 0;
		popData["All"]["Seropositive"][ageGroupStr][health_str] = 0;
		popData["All"]["Seronegative"][ageGroupStr][health_str] = 0;
	    }
	}
    }
    string vac_status;
    string sero_status;
    if(h->getSeroStatusAtVaccination() == true){
	sero_status = "Seropositive";
    }else{
	sero_status = "Seronegative";
    }
    if(h->isVaccinated() == true){
	vac_status = "Vaccine";
    }else{
	vac_status = "Placebo";
    }
    
    int group = getGroup(h->getAgeDays(currDay), groupsAges);    
    string ageGroupString;
    if(group >= 0){
	ageGroupString = std::to_string(groupsAges[group].min) + "_" + std::to_string(groupsAges[group].max);
    }
    
    if(h->getPopCount() == false){
	popData["All"]["All"]["All"]["Size"]++;
	popData[vac_status][sero_status]["All"]["Size"]++;    	
	popData["All"][sero_status]["All"]["Size"]++;    	
		
	if(group >= 0){
	    popData[vac_status][sero_status][ageGroupString]["Size"]++;
	    popData["All"]["All"][ageGroupString]["Size"]++;
	    popData["All"][sero_status][ageGroupString]["Size"]++;
	}	
	h->setPopCount(true);
    }
    if(h->getSeroprevalenceCount() == false){
	if(h->getPreviousInfections() > 0){
	    popData["All"]["All"]["All"]["Seroprevalence"]++;
	    popData[vac_status][sero_status]["All"]["Seroprevalence"]++;    	
	    popData["All"][sero_status]["All"]["Seroprevalence"]++;    	
	    if(group >= 0){
		popData[vac_status][sero_status][ageGroupString]["Seroprevalence"]++;
		popData["All"]["All"][ageGroupString]["Seroprevalence"]++;
		popData["All"][sero_status][ageGroupString]["Seroprevalence"]++;
	    }
	    h->setSeroprevalenceCount(true);
	}
    }
    
    
    if(h->infection != nullptr){
	if(floor(h->infection->getSymptomOnset()) == currDay){
	    popData["All"]["All"]["All"]["Inf"]++;
	    popData["All"][sero_status]["All"]["Inf"]++;
	    popData[vac_status][sero_status]["All"]["Inf"]++;
	    
	    if(h->getPreviousInfections() == 1){
		popData["All"]["All"]["All"]["cumAgeOfInfection"]+= h->getAgeDays(currDay);
		popData["All"][sero_status]["All"]["cumAgeOfInfection"]+= h->getAgeDays(currDay);
		popData[vac_status][sero_status]["All"]["cumAgeOfInfection"]+= h->getAgeDays(currDay);
		popData["All"]["All"]["All"]["firstInfections"]++;
		popData["All"][sero_status]["All"]["firstInfections"]++;
		popData[vac_status][sero_status]["All"]["firstInfections"]++;
	    }
	    
	    if(group >=0){
		popData["All"]["All"][ageGroupString]["Inf"]++;
		popData["All"][sero_status][ageGroupString]["Inf"]++;
		popData[vac_status][sero_status][ageGroupString]["Inf"]++;
	    }
	    if(h->isSymptomatic() == true){
		popData["All"]["All"]["All"]["Dis"]++;
		popData["All"][sero_status]["All"]["Dis"]++;
		popData[vac_status][sero_status]["All"]["Dis"]++;
		if(group >=0){
		    popData["All"]["All"][ageGroupString]["Dis"]++;
		    popData["All"][sero_status][ageGroupString]["Dis"]++;
		    popData[vac_status][sero_status][ageGroupString]["Dis"]++;
		}
		if(h->isHospitalized() == true){
		    popData["All"]["All"]["All"]["Hosp"]++;
		    popData["All"][sero_status]["All"]["Hosp"]++;
		    popData[vac_status][sero_status]["All"]["Hosp"]++;
		    
		    popData["All"]["All"]["All"]["cumAgeOfHospitalization"] += h->getAgeDays(currDay);
		    popData["All"][sero_status]["All"]["cumAgeOfHospitalization"] += h->getAgeDays(currDay);
		    popData[vac_status][sero_status]["All"]["cumAgeOfHospitalization"] += h->getAgeDays(currDay);
		    if(group >=0){
			popData["All"]["All"][ageGroupString]["Hosp"]++;
			popData["All"][sero_status][ageGroupString]["Hosp"]++;
			popData[vac_status][sero_status][ageGroupString]["Hosp"]++;
		    }
		    if(h->isDeadFromInfection()){
			int le = h->getDeathday();
			if(le == -1){
			    le = calculateExpectedLifeInDays(h->getAgeDays(currDay)) + currDay; 
			}
			if(le < currDay){
			    le = currDay;
			}
			discountedLifeExpectancy += le - currDay;
			
			popData["All"]["All"]["All"]["deathsFromInfection"]++;
			popData["All"][sero_status]["All"]["deathsFromInfection"]++; 
			popData[vac_status][sero_status]["All"]["deathsFromInfection"]++;

			popData["All"]["All"]["All"]["discountedLifeExpectancy"]+= le - currDay;
			popData["All"][sero_status]["All"]["discountedLifeExpectancy"] += le - currDay; 
			popData[vac_status][sero_status]["All"]["discountedLifeExpectancy"] += le - currDay;

			if(group >=0){
			    popData["All"]["All"][ageGroupString]["deathsFromInfection"]++;
			    popData["All"][sero_status][ageGroupString]["deathsFromInfection"]++; 
			    popData[vac_status][sero_status][ageGroupString]["deathsFromInfection"]++; 
			    
			    popData["All"]["All"][ageGroupString]["discountedLifeExpectancy"]+= le - currDay;
			    popData["All"][sero_status][ageGroupString]["discountedLifeExpectancy"]+= le - currDay;
			    popData[vac_status][sero_status][ageGroupString]["discountedLifeExpectancy"]+= le - currDay;
			}
		    }
		}
	    }
	}
    }
}

void Report::updateCohortReport(int currDay, Human * h){
    int cohortNum = h->getCohort();
    if(cohortNum < 1){
	return;
    }
    if(cohortData[cohortNum].empty()){
      vector<string> health_outcomes; health_outcomes.clear();
      vector<string> test_status_vec; test_status_vec.clear();
      vector<string> vax_status_vec; vax_status_vec.clear();
      vector<string> sero_status_vec; sero_status_vec.clear();
      health_outcomes.push_back("Size");
      health_outcomes.push_back("Inf");
      health_outcomes.push_back("Dis");
      health_outcomes.push_back("Hosp");
      health_outcomes.push_back("cumAgeOfHospitalization");
      health_outcomes.push_back("discountedLifeExpectancy");
      test_status_vec.push_back("Tested");
      test_status_vec.push_back("NonTested");
      vax_status_vec.push_back("Vaccine");
      vax_status_vec.push_back("Placebo");
      sero_status_vec.push_back("Seropositive");
      sero_status_vec.push_back("Seronegative");
      for(unsigned i = 0; i < health_outcomes.size(); i++){
	string health_str = health_outcomes[i];
	for(unsigned t = 0; t < test_status_vec.size(); t++){
	  string test_str = test_status_vec[t];
	  for(unsigned v = 0; v < vax_status_vec.size(); v++){
	    string vax_str = vax_status_vec[v];
	    for(unsigned s = 0; s < sero_status_vec.size(); s++){
	      string sero_str = sero_status_vec[s];		  
	      cohortData[cohortNum][test_str][vax_str][sero_str][health_str] = 0;
	    }
	  }
	}
      }
    }
    string tested_status;
    string vac_status;
    string sero_status;
    if(h->isTestedPreviousExposure()){
      tested_status = "Tested";
    }else{
      tested_status = "NonTested";
    }
    if(h->getSeroStatusAtVaccination() == true){
      sero_status = "Seropositive";
    }else{
      sero_status = "Seronegative";
    }
    if(h->isVaccinated() == true){
      vac_status = "Vaccine";
    }else{
      vac_status = "Placebo";
    }
    
    if(h->getCohortCount() == false){
	cohortData[cohortNum][tested_status][vac_status][sero_status]["Size"]++;    
	h->setCohortCount(true);
    }

    if(h->infection != nullptr){
	if(floor(h->infection->getSymptomOnset()) == currDay){
	    cohortData[cohortNum][tested_status][vac_status][sero_status]["Inf"]++;
	    if(h->isSymptomatic() == true){
		cohortData[cohortNum][tested_status][vac_status][sero_status]["Dis"]++;		
		if(h->isHospitalized()){
		    cohortData[cohortNum][tested_status][vac_status][sero_status]["Hosp"]++;
		    cohortData[cohortNum][tested_status][vac_status][sero_status]["cumAgeOfHospitalization"] += h->getAgeDays(currDay);
		    if(h->isDeadFromInfection()){
			deathsFromInfection++;
			int le = h->getDeathday();
			if(le == -1){
			    le = calculateExpectedLifeInDays(h->getAgeDays(currDay)) + currDay; 
			}
			if(le < currDay){
			    le = currDay;
			}
			cohortData[cohortNum][tested_status][vac_status][sero_status]["discountedLifeExpectancy"] += le - currDay;
		    }
		}
	    }
	}
    }   
}



int Report::getGroup(int age_, vector<rangeStruct> groups_temp){
    vector<rangeStruct>::iterator itAge = groups_temp.begin();
    int count = 0;
    for(; itAge != groups_temp.end(); itAge++){
	if((double )age_ / 365.0 >= (*itAge).min && (double) age_ / 365.0 < (*itAge).max){
	    return count;
	}
	count++;
    }
    return -1;
}

void Report::printSpatialReport(int currDay){
    if(!spatialData.empty()){
	for(auto it = spatialData.begin(); it != spatialData.end(); ++it){
	    outSpatial << (*it) << "\n";
	}
    }
}

void Report::printAgesReport(int currDay){
    outAges << currDay << ",";
    for(int i = 0; i < 5 ; i++){
	if(ageEvents[i] == 1){
	    if(ageStatus[0] > 0 || ageStatus[1] > 0){
		if(ageStatus[0] > 0 && ageStatus[1] > 0){
		    for(int j = 0; j < 4; j++){
			for(int k = 0; k <= discreteAges.max - discreteAges.min; k++){
			    if(printAgesPop == true){
				outAges << ageStats[k].status[j].nonevents[i]<<",";
			    }
			    outAges << ageStats[k].status[j].events[i];
			    if(i == ageMaxIndex && k == discreteAges.max - discreteAges.min && j == 3){
				outAges << "\n";
			    }else{
				outAges << ",";
			    }
			}
		    }
		}else{
		    int inc_ = 0;
		    int sum_ = 0;
		    if(ageStatus[0] > 0){
			inc_ = 2;
			sum_ = 1;
		    }else if(ageStatus[1] > 0){ 
			inc_ = 1;
			sum_ = 2;
		    }
		    for(int j = 0;j < 1 +  inc_;j += inc_){
			for(int k = 0; k <= discreteAges.max - discreteAges.min; k++){
			    if(printAgesPop == true){
				outAges << ageStats[k].status[j].nonevents[i] + ageStats[k].status[j + sum_].nonevents[i]<<",";
			    }
			    outAges << ageStats[k].status[j].events[i] + ageStats[k].status[j + sum_].events[i];
			    if(i == ageMaxIndex && k == discreteAges.max - discreteAges.min && j == inc_){
				outAges << "\n";
			    }else{
				outAges << ",";
			    }
			}
		    }
		}

	    }else{
		for(int k = 0; k <= discreteAges.max - discreteAges.min; k++){
		    if(printAgesPop == true){
			outAges << ageStats[k].total.nonevents[i]<<",";
		    }
		    outAges << ageStats[k].total.events[i];
		    if(i == ageMaxIndex && k == discreteAges.max - discreteAges.min){
			outAges << "\n";
		    }else{
			outAges << ",";
		    }
		}
	    }
	}
    }
}      
void Report::printCohortReport(int currDay, int vDay_){
  if(cohortData.size() > 0){
    //unordered_map<int,unordered_map<string,unordered_map<string,unordered_map<string, unordered_map<string,int>>>>>
    // <Cohort # < tested_status <vaccine status < serostatus at vac. <output name, output value>>>>
    for(unsigned i = 0; i < cohortData.size(); i++){
      if( ((double)(currDay - vDay_) / 365.0 ) >= i){
	for(auto it = cohortData[i].begin(); it != cohortData[i].end(); ++it){
	  for(auto vIt = (*it).second.begin(); vIt != (*it).second.end(); ++vIt){
	    for(auto sIt = (*vIt).second.begin(); sIt != (*vIt).second.end(); ++sIt){
	      vector<string> cohort_values; cohort_values.clear();
	      string outstring;
	      cohort_values.push_back(std::to_string(currDay));
	      cohort_values.push_back(std::to_string(i));
	      cohort_values.push_back((*it).first);
	      cohort_values.push_back((*vIt).first);
	      cohort_values.push_back((*sIt).first);
	      cohort_values.push_back(std::to_string((*sIt).second["Size"]));
	      cohort_values.push_back(std::to_string((*sIt).second["Inf"]));
	      cohort_values.push_back(std::to_string((*sIt).second["Dis"]));
	      cohort_values.push_back(std::to_string((*sIt).second["Hosp"]));
	      cohort_values.push_back(std::to_string((*sIt).second["cumAgeOfHospitalization"]));
	      cohort_values.push_back(std::to_string((*sIt).second["discountedLifeExpectancy"]));
	      Report::join(cohort_values,',',outstring);
	      outCohort << outstring;
	    }
	  }
	}
      }
    }    
  }
}

void Report::printPopReport(int currDay){    
    if(popData.size() > 0){	
	for(auto it = popData.begin(); it != popData.end(); ++it){
	    for(auto sIt = (*it).second.begin(); sIt != (*it).second.end(); ++sIt){
		for(auto aIt = (*sIt).second.begin(); aIt != (*sIt).second.end(); ++aIt){
		    vector<string> pop_values; pop_values.clear();
		    string outstring;
		    pop_values.push_back(std::to_string(currDay));
		    pop_values.push_back((*it).first);
		    pop_values.push_back((*sIt).first);
		    pop_values.push_back((*aIt).first);
		    pop_values.push_back(std::to_string((*aIt).second["Size"]));
		    pop_values.push_back(std::to_string((*aIt).second["Inf"]));
		    pop_values.push_back(std::to_string((*aIt).second["Dis"]));
		    pop_values.push_back(std::to_string((*aIt).second["Hosp"]));
		    pop_values.push_back(std::to_string((*aIt).second["Seroprevalence"]));
		    pop_values.push_back(std::to_string((*aIt).second["deathsFromInfection"]));
		    pop_values.push_back(std::to_string((*aIt).second["firstInfections"]));
		    pop_values.push_back(std::to_string((*aIt).second["cumAgeOfInfection"]));
		    pop_values.push_back(std::to_string((*aIt).second["cumAgeOfHospitalization"]));
		    pop_values.push_back(std::to_string((*aIt).second["discountedLifeExpectancy"]));
		    pop_values.push_back(std::to_string((*aIt).second["truePositives"]));
		    pop_values.push_back(std::to_string((*aIt).second["trueNegatives"]));
		    pop_values.push_back(std::to_string((*aIt).second["falsePositives"]));
		    pop_values.push_back(std::to_string((*aIt).second["falseNegatives"]));
		    pop_values.push_back(std::to_string((*aIt).second["totalVaccinated"]));
		    Report::join(pop_values,',',outstring);
		    outPop << outstring;
		}
	    }
	}
    }
}



void Report::printHeaders(){
    if(reportGroups == true){
	printPopHeader();
    }
    if(reportCohort == true){
	printCohortHeader();
    }
    if(reportAges == true){
	printAgesHeader();
    }
    if(reportSpatial == true){
	printSpatialHeader();
    }
}

void Report::join(const vector<string>& v, char c, string& s) {

    s.clear();

    for (vector<string>::const_iterator p = v.begin(); p != v.end(); ++p) {
	s += *p;
	if (p != v.end() - 1){
	    s += c;
	}else{
	    s +="\n";
	}
    }

}

void Report::printSpatialHeader(){
    outSpatial << "Variable,code,Day,Serotype\n";
}


void Report::printAgesHeader(){
    outAges << "day,";
    for(int i = 0; i < 5 ; i++){
	if(ageEvents[i] == 1){
	    if(ageStatus[0] > 0 || ageStatus[1] > 0){
		if(ageStatus[0] > 0 && ageStatus[1] > 0){
		    for(int j = 0;j < 2 ;j++){
			for(int jj = 2; jj < 4; jj++){
			    for(int k = 0; k <= discreteAges.max - discreteAges.min; k++){
				if(printAgesPop == true){
				    outAges << status[j].c_str() <<"_"<< status[jj].c_str() << "_age_" << k + discreteAges.min << "_";
				    outAges << "no"<<events[i].c_str()<<",";
				}
				outAges << status[j].c_str() << "_" << status[jj].c_str() << "_age_" << k + discreteAges.min << "_";
				outAges << events[i].c_str();
				if(i == ageMaxIndex && k == discreteAges.max - discreteAges.min && j == 1 && jj == 3){
				    outAges << "\n";
				}else{
				    outAges << ",";
				}
			    }
			}
		    }
		}else{
		    int ind_ = 0;
		    if(ageStatus[0] > 0){
			ind_ = 0;
		    }else if(ageStatus[1] > 0){ 
			ind_ = 2;
		    }
		    for(int j = ind_;j < ind_ + 2;j++){
			for(int k = 0; k <= discreteAges.max - discreteAges.min; k++){
			    if(printAgesPop == true){
				outAges << status[j].c_str() << "_age_" << k + discreteAges.min << "_";
				outAges << "no"<<events[i].c_str()<<",";
			    }
			    outAges << status[j].c_str() << "_age_" << k + discreteAges.min << "_";
			    outAges << events[i].c_str();
			    if(i == ageMaxIndex && k == discreteAges.max - discreteAges.min && j == ind_ + 1){
				outAges << "\n";
			    }else{
				outAges << ",";
			    }
			}
		    }
		}

	    }else{
		for(int k = 0; k <= discreteAges.max - discreteAges.min; k++){
		    if(printAgesPop == true){
			outAges << "age_" << k + discreteAges.min << "_";
			outAges << "no"<<events[i].c_str()<<",";
		    }
		    outAges << "age_" << k + discreteAges.min << "_";
		    outAges << events[i].c_str();
		    if(i == ageMaxIndex && k == discreteAges.max - discreteAges.min){
			outAges << "\n";
		    }else{
			outAges << ",";
		    }
		}
	    }
	}
    }
}
void Report::printPopHeader(){
    outPop << "Day,Vaccine,Serostatus,AgeGroup,Size,Infectious,Dis,Hosp,Seroprevalence,deathsFromInfection,firstInfections,cumAgeOfInfection,cumAgeOfHospitalization,discountedLifeExpectancy,truePositives,trueNegatives,falsePositives,falseNegatives,totalVaccinated\n";
}

void Report::printCohortHeader(){
    outCohort << "Day,Cohort,Tested,Vaccine,Serostatus,PopSize,Infectious,Dis,Hosp,cumAgeOfHospitalization,discountedLifeExpectancy\n";
}

void Report::resetReports(){
    if(reportGroups == true){
	resetGroupStats();
    }
    if(reportCohort == true){
	resetCohortStats();
    }
    if(reportAges == true){
	resetAgeStats();
    }
    if(reportSpatial == true){
	resetSpatialStats();
    }
}

void Report::resetSpatialStats(){
    spatialData.clear();
}


void Report::resetAgeStats(){
    ageStats.clear();
    if(discreteAges.min + discreteAges.max <= 0 || discreteAges.min >= discreteAges.max || discreteAges.min < 0){
	exit(1);
    }

    for(int k = 0; k <= discreteAges.max - discreteAges.min; k++){
	reportStats tempStats;
	for(int i = 0; i < 5; i++){
	    // There are four status VacSero+ VacSero- PlacSero+ PlacSero-
	    for(int j = 0; j < 4; j++){
		tempStats.status[j].events[i] = 0;
		tempStats.status[j].nonevents[i] = 0;
	    }
	    tempStats.total.events[i] = 0;
	    tempStats.total.nonevents[i] = 0;
	}
	ageStats.push_back(tempStats);
    }

    for(int i = 0; i < 5; i++){
	if(ageEvents[i] == 1){
	    ageMaxIndex = i;
	}
    }
}

void Report::resetCohortStats(){
    cohortData.clear();
}

void Report::resetGroupStats(){
    popData.clear();
}

void Report::updateFOI(int currDay, Human * h){
    human_counts[currDay]++;
    string tmpzone = h->getZoneID();
    if(h->isInfected()){
	if(h->infection != nullptr){
	    if(h->infection->getStartDay() == currDay){
		int sero = h->infection->getInfectionType();
		if(sero > 0){
		    dailyFOI[currDay][sero]["newinf"]++;
		    dailyFOI[currDay][sero][tmpzone + "newinf"]++;
		}
	    }
	    if(h->isSymptomatic() && floor(h->infection->getSymptomOnset()) == currDay){
		int sero = h->infection->getInfectionType();
		if(sero > 0){
		    dailyFOI[currDay][sero]["newsymp"]++;
		}	       		
	    }
	}
    }    
    // Check for immunity against all the serotypes
    for(unsigned i = 1; i < 5; i++){
	if(!h->isPermImmune(i)){
	    dailyFOI[currDay][i]["sus"]++;
	    dailyFOI[currDay][i][tmpzone + "sus"]++;
	}
	if(!h->isImmune(i)){
	    dailyFOI[currDay][i]["sustmp"]++;
	}	    
    }
}

void Report::updateSecondaryCases(int currDay, Human * h){
    if(h != nullptr){
	if(h->isInfected()){
	    if(h->infection != nullptr){
		unsigned sero = h->infection->getInfectionType();
		int infDay = h->infection->getStartDay();
		string id = h->getPersonID();
		secondaryCases[infDay][sero][id] = h->getR0(sero);
	    }
	}
	if(h->isImported()){
	    if(h->infection != nullptr){
		unsigned sero = h->infection->getInfectionType();
		int infDay = h->infection->getStartDay();
		string id = h->getVisitorID();
		secondaryCases[infDay][sero][id] = h->getVisitorR0(sero);
	    }
	}
    }
}

void Report::addImportation(int currDay, int sero_in, Human * h){
    if(sero_in > 0 && sero_in < 5){
	dailyFOI[currDay][sero_in]["imports"]++;
    }
}

void Report::printFOIReport(int lastDay){    
    outFOI.open(outputFOIFile);
    if (!outFOI.good()) {
	exit(1);
    }
    std::vector<string> headerstr; string headout;
    headerstr.clear();
    for(int s = 1; s < 5; s++){
	if(foiTypes[s - 1]){
	    if(printR0 == true){
		headerstr.push_back("R0_Denv"+std::to_string(s));
	    }
	    headerstr.push_back("Denv"+std::to_string(s));
	    headerstr.push_back("Susceptible_" + std::to_string(s));
	    headerstr.push_back("Susceptible_temp_" + std::to_string(s));
	    headerstr.push_back("Infectious_" + std::to_string(s));
	    headerstr.push_back("Symptomatic_" + std::to_string(s));	    
	    headerstr.push_back("MozSusceptible_" + std::to_string(s));
	    headerstr.push_back("MozExposed_" + std::to_string(s));
	    headerstr.push_back("MozInfectious_" + std::to_string(s));
	    headerstr.push_back("Importations_" + std::to_string(s));
	}
	if(printZonesFOI == true){
	    for(auto locIt = zonesToPrint.begin(); locIt != zonesToPrint.end();){
		headerstr.push_back("FOI_" + *locIt + std::to_string(s));
		++locIt;
	    }
	}
    }
    if(reportFOIOutbreakSymptomatics){
	headerstr.push_back("ReportedSymptomatics");
	headerstr.push_back("OutbreakResponseLocations");
    }
    headerstr.push_back("Humans");
    
    if(!headerstr.empty()){
	outFOI << "day,";
	Report::join(headerstr,',',headout);
	outFOI << headout;	
    }
    
    for(int d = 0; d < lastDay; d++){
	string outstring;
	vector<string> foi_values; foi_values.clear();
	for(unsigned s = 1; s < 5; s++){
	    if(foiTypes[s - 1] != 1){
		continue;
	    }
	    if(printR0 == true){
		string R0str = "0.000";
		if(secondaryCases.find(d) != secondaryCases.end()){		
		    if(secondaryCases[d].find(s) != secondaryCases[d].end()){
			int dailyR0 = 0;
			int sumInfectors = 0;
			for(auto idt = secondaryCases[d][s].begin(); idt != secondaryCases[d][s].end();idt++){
			    if(idt->second >= 0){
				dailyR0 += idt->second;
				sumInfectors++;
			    }
			}
			R0str = sumInfectors > 0 ? std::to_string((double) dailyR0 / (double) sumInfectors) : 0;
		    }
		}
		foi_values.push_back(R0str);
	    }

	    string foistr = "0.000";	    
	    if(dailyFOI.find(d) != dailyFOI.end()){
		if(dailyFOI[d].find(s) != dailyFOI[d].end()){
		    double foi_temp = dailyFOI[d][s]["sustmp"] > 0 ? (double) dailyFOI[d][s]["newinf"] / (double) dailyFOI[d][s]["sustmp"] : 0;
		    foistr = std::to_string(foi_temp);
		}
	    }
	    foi_values.push_back(foistr);

	    foistr = "0.000";	    
	    if(dailyFOI.find(d) != dailyFOI.end()){
		if(dailyFOI[d].find(s) != dailyFOI[d].end()){
		    foistr = std::to_string(dailyFOI[d][s]["sus"]);
		}
	    }
	    foi_values.push_back(foistr);

	    foistr = "0.000";	    
	    if(dailyFOI.find(d) != dailyFOI.end()){
		if(dailyFOI[d].find(s) != dailyFOI[d].end()){
		    foistr = std::to_string(dailyFOI[d][s]["sustmp"]);
		}
	    }
	    foi_values.push_back(foistr);

	    foistr = "0.000";	    
	    if(dailyFOI.find(d) != dailyFOI.end()){
		if(dailyFOI[d].find(s) != dailyFOI[d].end()){
		    foistr = std::to_string(dailyFOI[d][s]["newinf"]);
		}
	    }
	    foi_values.push_back(foistr);

	    foistr = "0.000";	    
	    if(dailyFOI.find(d) != dailyFOI.end()){
		if(dailyFOI[d].find(s) != dailyFOI[d].end()){
		    foistr = std::to_string(dailyFOI[d][s]["newsymp"]);
		}
	    }
	    foi_values.push_back(foistr);
	    
	    foistr = "0.000";	    
	    if(dailyFOI.find(d) != dailyFOI.end()){
		if(dailyFOI[d].find(s) != dailyFOI[d].end()){
		    foistr = std::to_string(dailyFOI[d][s]["mozsus"]);
		}
	    }
	    foi_values.push_back(foistr);

	    foistr = "0.000";	    
	    if(dailyFOI.find(d) != dailyFOI.end()){
		if(dailyFOI[d].find(s) != dailyFOI[d].end()){
		    foistr = std::to_string(dailyFOI[d][s]["mozexp"]);
		}
	    }
	    foi_values.push_back(foistr);
	    
	    foistr = "0.000";	    
	    if(dailyFOI.find(d) != dailyFOI.end()){
		if(dailyFOI[d].find(s) != dailyFOI[d].end()){
		    foistr = std::to_string(dailyFOI[d][s]["mozinf"]);
		}
	    }
	    foi_values.push_back(foistr);

	    foistr = "0.000";	    
	    if(dailyFOI.find(d) != dailyFOI.end()){
		if(dailyFOI[d].find(s) != dailyFOI[d].end()){
		    foistr = std::to_string(dailyFOI[d][s]["imports"]);
		}
	    }
	    foi_values.push_back(foistr);
	    
	    if(printZonesFOI == true){		
		for(auto locIt = zonesToPrint.begin(); locIt != zonesToPrint.end();){
		    foistr = "0.000";
		    if(dailyFOI[d].find(s) != dailyFOI[d].end()){
			if(dailyFOI[d][s].find(*locIt + "newinf") != dailyFOI[d][s].end()){
			    double foi_temp = dailyFOI[d][s][*locIt + "sus"] > 0 ? (double) dailyFOI[d][s][*locIt + "newinf"] / (double) dailyFOI[d][s][*locIt + "sus"] : 0;
			    foistr = std::to_string(foi_temp);			    
			}
		    }
		    ++locIt;
		    foi_values.push_back(foistr);
		}
	    }	    
	}
	if(reportFOIOutbreakSymptomatics){
	    string outbreakstr= "0";
	    if(dailyOutbreakResponseSymp.find(d) != dailyOutbreakResponseSymp.end()){
		outbreakstr = std::to_string(dailyOutbreakResponseSymp[d]);
	    }
	    foi_values.push_back(outbreakstr);
	    outbreakstr= "0";
	    
	    if(dailyOutbreakResponseLocations.find(d) != dailyOutbreakResponseLocations.end()){
		outbreakstr = std::to_string(dailyOutbreakResponseLocations[d]);
	    }
	    foi_values.push_back(outbreakstr);	    
	}

	string humstr = "0";
	if(human_counts.find(d) != human_counts.end()){
	    humstr = std::to_string(human_counts[d]);
	}
	foi_values.push_back(humstr);
	
	if(!foi_values.empty()){
	    Report::join(foi_values,',',outstring);
	    outFOI << d << ",";
	    outFOI << outstring;
	}    	
    }
    outFOI.close();
}


void Report::finalizeReport(int currDay){
    outCohort.close();
    outAges.close();
    outPop.close();
    outSpatial.close();
    if(reportFOI == true){
	this->printFOIReport(currDay);
    }
}



//Report::~Report() {
//}

