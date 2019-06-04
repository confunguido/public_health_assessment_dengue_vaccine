#ifndef _HOUSE_H
#define _HOUSE_H



#include "human.h"
#include "scales.h"

extern gsl_rng *rGsl;



using namespace std;



class House
{
  string code;
  vector<Human> residents;
  void setDist(
    map<string,Lot> *);
  void setProbLoc(
    map<string,Lot> *,
    map<string,float> *,
    map<string,vector<double> > *,
    map<string,vector<vector<string> > > *);
public:
  House(){;};
  string getCode(){return code;}
  unsigned getNumRes(){return residents.size();}
  vector<Human>* getRes(){return &residents;}
  void setHouse(
    string codeIn,
    vector<pair<char,int> > *,
    map<string,map<string,float> > *,
    map<string,map<string,float> > *,
    map<string,Lot> *,
    map<string,vector<double> > *,
    map<string,vector<vector<string> > > *);
}; // class House



#endif // _HOUSE_H
