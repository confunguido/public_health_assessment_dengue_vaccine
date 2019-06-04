#ifndef _CITY_H
#define _CITY_H



#include "house.h"
#include "scales.h"

extern gsl_rng *rGsl;



using namespace std;



class City
{
	map<string,Lot> locations;
	map<string,House> houses;
	void getHouseList();
	void loadHouseComp(
		vector<vector<pair<char,int> > > *);
	void loadLocations();
	void loadParmIn(
	  map<string,map<string,float> > *,
	  map<string,map<string,float> > *);
	void processParmFD(
	  map<string,float> *,
	  map<string,map<string,float> > *);
public:
	City();
	map<string,Lot>* getLocations(){return &locations;}
	map<string,House>* getHouses(){return &houses;}
}; // class City



#endif // _CITY_H
