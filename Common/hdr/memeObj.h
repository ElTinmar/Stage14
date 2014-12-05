#ifndef MEMEOBJ_H
#define MEMEOBJ_H
#include <iostream>
#include <list>
#include <vector>
#include <map>
#include <set>
#include <string>
#include <cstdlib>
#include <fstream>
#include "motif.h"

using namespace std;

class MemeObj
{
	friend class HitList;
	
  public:
  void      					afficher            (ostream&) const;
  void								set_version					(string);
  void								set_alphabet				(string);
  void								add_strand					(string);
  void								set_bg_frq					(map<string, double>);
  void								set_motifs					(vector<Motif>);
  const set<char>& 		get_alphabet	      () const;
  
  void								entropy							(ostream&);
  string							ref_to_altname			(string);
  int									ref_to_index				(string);
  map<string, int>		name_index					();
  set<string>					motifs_altname			();
  
  static MemeObj			load_meme						(string);
  
  private:
  string 								version;
  set<char>							alphabet;
  set<string>						strands;
  map<string, double> 	bg_frequencies;
  vector<Motif>					motifs;
};

ostream &operator<<(ostream &flux, MemeObj const& d);

#endif
