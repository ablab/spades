/***************************************************************************
 * Title:          IntersectProfiles.cpp 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  11/26/2009
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#include <iostream>
#include <set>
#include <vector>
#include <string>
#include <algorithm>
#include <fstream>


class Profile {
public:
  std::string title, value;
  Profile &operator=(const Profile &p) {
		if (this != &p) {
			title = p.title;
			value = p.value;
		}
    return *this;
  }
  Profile(std::string &t, std::string &v) {
    title = t;
    value = v;
  }
  Profile() {
  }
  Profile(const Profile &p) {
    *this = p;
  }
  bool operator<(const Profile &p) const {
    return title < p.title;
  }
};


class CmpProfile {
public:
  ssize_t operator()(const Profile &a, const Profile &b) {
    return a.title < b.title;
  }
};

class CmpTitle {
public:
  ssize_t operator()(Profile &a, const std::string &title) {
    return a.title < title;
  }
  ssize_t operator()(const std::string &title, Profile &a ) {
    return title < a.title;
  }

};

int main(int argc, char *argv[]) {

  std::string listA, listB, intersectionA, intersectionB;
  if (argc < 5) {
    std::cout << "usage: intersectProfiles listA listB intersectionA, intersectionB" << std::endl;
    exit(1);
  }

  std::vector<std::string> aTitles, bTitles, aValues, bValues;
  std::set<std::string> aTitleSet, bTitleSet;
  
  std::ifstream astrm, bstrm;
  astrm.open(argv[1]);
  bstrm.open(argv[2]);
  
  if (!astrm or !bstrm) {
    std::cout << "problem opening input " << std::endl;
    exit(1);
  }

  std::ofstream aistrm, bistrm;
  aistrm.open(argv[3]);
  bistrm.open(argv[4]);
  std::string title;
  std::string value;
  //UNUSED// ssize_t naRead = 0;
  std::vector<Profile> aProf, bProf;
  while(astrm) {
    astrm >> title;
    std::getline(astrm, value);
    aProf.push_back(Profile(title, value));
  }
  astrm.close();
  while(bstrm) {
    bstrm >> title;
    std::getline(bstrm, value);
    bProf.push_back(Profile(title, value));
  }
  bstrm.close();
  CmpProfile cmp;
  std::cout << "sorting " << argv[1] << std::endl;
  std::sort(aProf.begin(), aProf.end(), cmp);
  std::cout << "sorting " << argv[2] << std::endl;
  std::sort(bProf.begin(), bProf.end(), cmp);
  ssize_t i;
  
  std::set<std::string> titleIsect;
  std::vector<Profile> profIsect;

  std::back_insert_iterator<std::vector<Profile> > profIsectInsert(profIsect);
  std::cout << "interseting " << argv[1] << " and " << argv[2] << std::endl;
  std::set_intersection(aProf.begin(), aProf.end(),
			bProf.begin(), bProf.end(),
			profIsectInsert);

  std::cout << "done with intersect, got " << profIsect.size() << " common elements " << std::endl;
  // Now output the results

  ssize_t foundA, foundB;
  CmpTitle titleCmp;
  foundA = foundB = 0;
  for (i = 0; i < bProf.size(); i++ ){ 
    if (std::binary_search(profIsect.begin(), 
			   profIsect.end(), bProf[i].title, titleCmp)) {
      bistrm << bProf[i].title << " " << bProf[i].value << std::endl;
      ++foundB;
    }
    else {
      //      std::cout << "didn't find " << bProf[i].title << std::endl;
    }
  }  
  std::cout << "matchedB " << foundB << " / " << bProf.size() << std::endl;
  for (i = 0; i < aProf.size(); i++ ){ 
    if (std::binary_search(profIsect.begin(), 
			   profIsect.end(), aProf[i].title, titleCmp)) {
      aistrm << aProf[i].title << " " << aProf[i].value << std::endl;
      ++foundA;
    }
  }
  std::cout << "matchedA " << foundA << " / " << aProf.size() << std::endl;
  return 0;
}
  
