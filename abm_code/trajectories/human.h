#ifndef _HUMAN_H
#define _HUMAN_H



#include "scales.h"

extern gsl_rng *rGsl;



using namespace std;



class Human;
struct Lot{
  string type;
  float xcor;
  float ycor;
  float dist;
  map<Human *,float> time;
}; // struct Lot



class Human
{
  string homeCode;
  char sex;
  int age;
  map<string,float> durnRep;
  map<string,float> freqRep;
  unsigned setRepSize(
    map<string,float> parm);
  void setFreqDurn(
    map<string,string> *,
    map<string,float> *,
    map<string,float> *,
    map<string,map<string,float> > *,
    map<string,Lot> *);
  void setType(
    unsigned repSize,
    map<string,float> *,
    multiset<string> *);
  void setWhere(
    multiset<string> *,
    map<string,vector<double> > *,
    map<string,vector<vector<string> > > *,
    map<string,Lot> *,
    map<string,string> *);
public:
  Human(
    string codeIn,
    char sexIn,
    int ageIn,
    map<string,map<string,float> > *,
    map<string,map<string,float> > *,
    map<string,Lot> *,
    map<string,vector<double> > *,
    map<string,vector<vector<string> > > *);
  string getHomeCode(){return homeCode;}
  map<string,float>* getDurnRep(){return &durnRep;}
  map<string,float>* getFreqRep(){return &freqRep;}
  void writeTrajectories(int);
}; // class Human



#endif // _HUMAN_H
