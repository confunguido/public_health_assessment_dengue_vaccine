#include "human.h"



Human::Human(
  string codeIn,
  char sexIn,
  int ageIn,
  map<string,map<string,float> > *parmIn,
  map<string,map<string,float> > *parmFD,
  map<string,Lot> *locations,
  map<string,vector<double> > *probLoc,
  map<string,vector<vector<string> > > *distLocs)
{
  homeCode = codeIn;
  sex = sexIn;
  age = ageIn;

  unsigned repSize = setRepSize((*parmIn)["repSize.txt"]);

  multiset<string> repType;
  setType(
    repSize,
    &(*parmIn)["propLoc.txt"],
    &repType);

  map<string,string> repWhere;
  setWhere(
    &repType,
    probLoc,
    distLocs,
    locations,
    &repWhere);

  setFreqDurn(
    &repWhere,
    &(*parmIn)["timeHome.txt"],
    &(*parmIn)["timeFI.txt"],
    parmFD,
    locations);
} // Human constructor



unsigned Human::setRepSize(
  map<string,float> parm)
{
  unsigned repSize;

  do
  {
    repSize = gsl_ran_negative_binomial(rGsl, parm["p"], parm["n"]);
  }
  while(repSize < 2);

  return repSize;
} // setRepSize()



void Human::setFreqDurn(
  map<string,string> *repWhere,
  map<string,float> *parmHome,
  map<string,float> *parmFI,
  map<string,map<string,float> > *parmFD,
  map<string,Lot> *locations)
{
  durnRep.insert(pair<string,float>
    (homeCode,
    exp((*parmHome)["durnMean"] + gsl_ran_gaussian(rGsl, (*parmHome)["durnStd"]))));
  freqRep.insert(pair<string,float>
    (homeCode,
    exp((*parmHome)["freqMean"] + gsl_ran_gaussian(rGsl, (*parmHome)["freqStd"]))));

  double freq, durn, dist;

  for(map<string,string>::iterator itr = ++repWhere->begin();
      itr != repWhere->end();
      ++itr)
  {
    if(itr->second == "FI")
    {
      durnRep.insert(pair<string,float>
        (itr->first,
        exp((*parmFI)["FI_durnMeanA"] + gsl_ran_gaussian(rGsl, (*parmFI)["FI_durnStd"]))));
      freqRep.insert(pair<string,float>
        (itr->first,
        exp((*parmFI)["FI_freqMeanA"] + gsl_ran_gaussian(rGsl, (*parmFI)["FI_freqStd"]))));
    }
    else
    {
      dist = (*locations)[itr->first].dist;

      gsl_ran_bivariate_gaussian(rGsl,
        (*parmFD)[itr->second]["freqStd"],
        (*parmFD)[itr->second]["durnStd"],
        (*parmFD)[itr->second]["corr"],
        &freq, &durn);

      freq += .5 * (*parmFD)[itr->second]["freqMeanA"] +
        (*parmFD)[itr->second]["freqMeanA"] * (*parmFD)[itr->second]["freqMeanB"] *
        exp((*parmFD)[itr->second]["freqMeanC"] * dist) /
        (1 + exp((*parmFD)[itr->second]["freqMeanC"] * dist));
      durn += .5 * (*parmFD)[itr->second]["durnMeanA"] +
        (*parmFD)[itr->second]["durnMeanA"] * (*parmFD)[itr->second]["durnMeanB"] *
        exp((*parmFD)[itr->second]["durnMeanC"] * dist) /
        (1 + exp((*parmFD)[itr->second]["durnMeanC"] * dist));

      durnRep.insert(pair<string,float>(itr->first, exp(durn)));
      freqRep.insert(pair<string,float>(itr->first, exp(freq)));
    } // if FI else not
  } // for each repertoire location
} // setTime()



void Human::setType(
  unsigned repSize,
  map<string,float> *parm,
  multiset<string> *types)
{
  double p[parm->size()];
  map<string,float>::iterator itr = parm->begin();
  for(int i = 0; i < parm->size(); ++i)
    p[i] = (itr++)->second;

  unsigned n[parm->size()];
  gsl_ran_multinomial(rGsl, parm->size(), repSize, p, n);

  itr = parm->begin();
  for(int i = 0; i < parm->size(); ++i)
  {
    for(int j = 0; j < n[i]; ++j)
      types->insert(itr->first);
    ++itr;
  } // end for i  
} // setType()



void Human::setWhere(
  multiset<string> *repType,
  map<string,vector<double> > *probLoc,
  map<string,vector<vector<string> > > *distLocs,
  map<string,Lot> *locations,
  map<string,string> *repWhere)
{
  repWhere->insert(pair<string,string>(homeCode,"HOME"));
  (*locations)[homeCode].time.insert(make_pair(this,0));

  // find number of each location type in the repertoire
  map<string,unsigned> types;
  for(multiset<string>::iterator itr = repType->begin();
      itr != repType->end();
      ++itr)
    if(types.find(*itr) == types.end())
      types.insert(pair<string,unsigned>(*itr, repType->count(*itr)));

  // select locations of each type
  char countFI = '1';
  for(map<string,unsigned>::iterator itrType = types.begin();
      itrType != types.end();
      ++itrType)
  {
    unsigned count = 0;

    if(itrType->first == "FI")
      while(count++ < itrType->second)
      {
        string nameFI = "FI";
        nameFI.append(1, countFI++);
        (*repWhere)[nameFI] = "FI";
      }
    else
    {
      vector<double> *vecPtr = &((probLoc->find(itrType->first))->second);
      double * p = &(*vecPtr)[0];
      unsigned n[vecPtr->size()];

      while(count < itrType->second)
      {
        gsl_ran_multinomial(rGsl, vecPtr->size(), 1, p, n);

        unsigned i = 0, dest;
        for(; n[i] != 1; ++i){;}

        unsigned numLocs = (*distLocs)[itrType->first][i].size();        
        unsigned src[numLocs], k = 0;
        for(; k < numLocs; ++k){src[k] = k;}

        gsl_ran_choose(rGsl, &dest, 1, src, numLocs, sizeof(unsigned));

        if(repWhere->find((*distLocs)[itrType->first][i][dest]) == repWhere->end())
        {
          (*repWhere)[(*distLocs)[itrType->first][i][dest]] = itrType->first;
          (*locations)[(*distLocs)[itrType->first][i][dest]].time.insert(make_pair(this,0));
          ++count;
        } // end if new location
      } // end while there are still some repertoire locations unchosen
    } // end if not FI
  } // end for each location type
} // setWhere()



void Human::writeTrajectories(int residentNumber)
{
  float freqSum, runif, rcum, timeBout, timeSum;
  int next, i, n;
  string nextCode, nowCode;
  map<string,float>::iterator mapItr;

  freqSum = 0;
  for(mapItr = freqRep.begin(); mapItr != freqRep.end(); mapItr++)
    freqSum += mapItr->second;

  for(n = 0; n < 5; n++)
  {
    timeSum = 0;

    cout <<
      homeCode << "," <<
      residentNumber << "," <<
      sex << "," << age;

    next = gsl_rng_uniform_int(rGsl, durnRep.size());
    i = 0, mapItr = durnRep.begin();
    for(; i < next; i++){mapItr++;}
    nowCode = mapItr->first;

    // pick location
    do
    {
      do
      {
        runif = gsl_rng_uniform(rGsl) * freqSum; rcum = 0;
        mapItr = freqRep.begin();
        do
        {
          rcum += mapItr->second;
          mapItr++;
        }
        while(rcum < runif);
        nextCode = (--mapItr)->first;            
      }
      while(nextCode == nowCode);
      nowCode = nextCode;

      // draw time at location
      timeBout = min(1.0, gsl_ran_exponential(rGsl, 1 / durnRep.at(nowCode)) / 17);
      timeSum += timeBout;

      if(timeSum < 1)
        // continue;
        cout << "," << nowCode << "," << timeBout;
      else
      {
        cout << "," << nowCode << "," << 1 - timeSum + timeBout;
        break;        
      }
    }
    while(timeSum < 1);

    cout << endl;
  }

} // writeTrajectories()

