#include "house.h"



void House::setDist(map<string,Lot> *locations)
{
  float x_house, y_house;
  x_house = (*locations)[code].xcor;
  y_house = (*locations)[code].ycor;

  for(map<string,Lot>::iterator itr = locations->begin();
      itr != locations->end();
      ++itr)
  {
    itr->second.dist = sqrt(
      pow(itr->second.xcor - x_house, 2) +
      pow(itr->second.ycor - y_house, 2));
  } // end for
} // setDist()



void House::setHouse(
  string codeIn,
  vector<pair<char,int> > * houseComp,
  map<string,map<string,float> > * parmIn,
  map<string,map<string,float> > * parmFD,
  map<string,Lot> * locations,
  map<string,vector<double> > * probLoc,
  map<string,vector<vector<string> > > * distLocs)
{
  code = codeIn;
  setDist(locations);

  setProbLoc(locations, &(*parmIn)["bins.txt"], probLoc, distLocs);

  int residentNumber = 1;
  for(vector<pair<char,int> >::iterator itr = houseComp->begin(); itr != houseComp->end(); ++itr)
  {
    Human * house_resident =
      new Human(
        code,
        itr->first,
        itr->second,
        parmIn,
        parmFD,
        locations,
        probLoc,
        distLocs);
    residents.push_back(*house_resident);

    house_resident->writeTrajectories(residentNumber++);
  } // for each resident of the house
} // setHouse()



void House::setProbLoc(
  map<string,Lot> *locations,
  map<string,float> *parmIn,
  map<string,vector<double> > *probLoc,
  map<string,vector<vector<string> > > *distLocs)
{
  (*probLoc)["ALOJAMIENTO"] = vector<double>((*parmIn)["max"] / (*parmIn)["interval"], 0);
  (*probLoc)["CEMENTERIOS"] = vector<double>((*parmIn)["max"] / (*parmIn)["interval"], 0);
  (*probLoc)["EDUCACION"] = vector<double>((*parmIn)["max"] / (*parmIn)["interval"], 0);
  (*probLoc)["HOUSE"] = vector<double>((*parmIn)["max"] / (*parmIn)["interval"], 0);
  (*probLoc)["IGLESIAS"] = vector<double>((*parmIn)["max"] / (*parmIn)["interval"], 0);
  (*probLoc)["INSTITUCIONES"] = vector<double>((*parmIn)["max"] / (*parmIn)["interval"], 0);
  (*probLoc)["OTROS"] = vector<double>((*parmIn)["max"] / (*parmIn)["interval"], 0);
  (*probLoc)["RECREACION"] = vector<double>((*parmIn)["max"] / (*parmIn)["interval"], 0);
  (*probLoc)["SALUD"] = vector<double>((*parmIn)["max"] / (*parmIn)["interval"], 0);
  (*probLoc)["TIENDAS"] = vector<double>((*parmIn)["max"] / (*parmIn)["interval"], 0);

  vector<vector<string> > initDistLocs;
  vector<string> initLocs;
  for(int i = 0; i < (*parmIn)["max"] / (*parmIn)["interval"]; ++i)
    initDistLocs.push_back(initLocs);
  (*distLocs)["ALOJAMIENTO"] = initDistLocs;
  (*distLocs)["CEMENTERIOS"] = initDistLocs;
  (*distLocs)["EDUCACION"] = initDistLocs;
  (*distLocs)["HOUSE"] = initDistLocs;
  (*distLocs)["IGLESIAS"] = initDistLocs;
  (*distLocs)["INSTITUCIONES"] = initDistLocs;
  (*distLocs)["OTROS"] = initDistLocs;
  (*distLocs)["RECREACION"] = initDistLocs;
  (*distLocs)["SALUD"] = initDistLocs;
  (*distLocs)["TIENDAS"] = initDistLocs;

  for(map<string,Lot>::iterator itr = locations->begin();
      itr != locations->end(); ++itr)
  {
    ++(*probLoc)[itr->second.type][floor(itr->second.dist / (*parmIn)["interval"])];
    (*distLocs)[itr->second.type][floor(itr->second.dist / (*parmIn)["interval"])].push_back(itr->first);
  } // end for
} // setProbLoc()
