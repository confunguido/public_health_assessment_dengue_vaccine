#include "city.h"



City::City()
{
  vector<vector<pair<char,int> > > houseComp;
  loadHouseComp(&houseComp);
  loadLocations();
  getHouseList();

  map<string,map<string,float> > parmIn;
  map<string,map<string,float> > parmFD;
  loadParmIn(&parmIn, &parmFD);

  map<string,vector<double> > probLoc;
  map<string,vector<vector<string> > > distLocs;

  vector<vector<pair<char,int> > >::iterator itrComp = houseComp.begin();
  map<string,House>::iterator itrHouse = houses.begin();

  std::random_device rd;
  std::mt19937 g(rd());
  shuffle(houseComp.begin(), houseComp.end(), g);

  for(; itrHouse != houses.end() && itrComp != houseComp.end();)
  {
    itrHouse->second.setHouse(
      itrHouse->first,
      &(*itrComp),
      &parmIn,
      &parmFD,
      &locations,
      &probLoc,
      &distLocs);
    itrHouse++;
    itrComp++;
  }
} // City constructor



void City::getHouseList()
{
  for(map<string,Lot>::iterator itr = locations.begin();
      itr != locations.end();
      ++itr)
    if(itr->second.type == "HOUSE")
      houses.insert(pair<string,House>(itr->first,House()));
} // getHouseList()



void City::loadHouseComp(
  vector<vector<pair<char,int> > > * houseCompPtr)
{
  string filename = "wholePop.txt";
  ifstream NewInfile(filename.c_str());

  char sexFile;
  int ageFile;
  vector<pair<char,int> > compFile;

  while(!NewInfile.eof())
  {
    NewInfile >> sexFile;
    while(sexFile != 'N')
    {
      NewInfile >> ageFile;
      compFile.push_back(pair<char,int>(sexFile,ageFile));
      NewInfile >> sexFile;
    }
    NewInfile.ignore(200, '\n');
// cout << compFile.size() << ",";
    if(!compFile.empty())
      houseCompPtr->push_back(compFile);

    compFile.clear();
  } // end while
} // loadHouseComp()



void City::loadLocations()
{
  string filename = "locations.txt";
  ifstream NewInfile(filename.c_str());

  string codeFile, typeFile;
  unsigned blockFile;
  float xcorFile, ycorFile;
  Lot temp;

  while(!NewInfile.eof())
  {
    NewInfile >> codeFile >> typeFile >> blockFile >> xcorFile >> ycorFile;
    NewInfile.ignore(200, '\n');

    temp.type = typeFile;
    temp.xcor = xcorFile;
    temp.ycor = ycorFile;
    temp.dist = -1;
    locations.insert(pair<string,Lot>(codeFile, temp));
  } // end while
} // loadLocations()



void City::loadParmIn(
  map<string,map<string,float> > *parmIn,
  map<string,map<string,float> > *parmFD)
{
  vector<string> filename;
  filename.push_back("bins.txt"),
  filename.push_back("biting.txt"),
  filename.push_back("freqDurn.txt"),
  filename.push_back("propLoc.txt"),
  filename.push_back("repSize.txt"),
  filename.push_back("timeFI.txt"),
  filename.push_back("timeHome.txt"),
  filename.push_back("where.txt");

  string parFile;
  float valFile;
  map<string,float> temp;

  for(vector<string>::iterator itr = filename.begin();
      itr != filename.end();
      ++itr)
  {
    ifstream NewInfile(itr->c_str());
    while(!NewInfile.eof())
    {
      NewInfile >> parFile >> valFile;
      NewInfile.ignore(200, '\n');

      temp.insert(pair<string,float>(parFile, valFile));
    } // end while

    parmIn->insert(pair<string,map<string,float> >(*itr, temp));
    temp.clear();
  } // end for

  processParmFD(&(*parmIn)["freqDurn.txt"], parmFD);

  parmIn->erase("freqDurn.txt");
} // loadParmIn()



void City::processParmFD(
  map<string,float> *parmIn,
  map<string,map<string,float> > *parmOut)
{
  unsigned pos;
  string type, parm;

  for(map<string,float>::iterator itr = parmIn->begin();
      itr != parmIn->end();
      ++itr)
  {
    pos = itr->first.find("_");
    type = itr->first.substr(0, pos);
    parm = itr->first.substr(pos + 1, itr->first.length() - pos + 1);

    (*parmOut)[type].insert(pair<string,float>(parm, itr->second));
  }
} // processParmFD
