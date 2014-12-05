#ifndef HITLIST_H
#define HITLIST_H
#include <iostream>
#include <list>
#include <vector>
#include <map>
#include <set>
#include <tuple>
#include <string>
#include <cstdlib>
#include <fstream>
#include "cluster.h"
#include "affichage.h"
#include "memeObj.h"

using namespace std;

class HitList
{
	friend class Affichage;
	
	public:
  
  					 HitList								(int, MemeObj, Dictionnaire);
  void       afficher            		(ostream&) const;
  void       to_hitlist            	(ostream&) const;
  void			 afficher_annuaire			(ostream&);
  void			 afficher_annuaire2			(ostream&);
  void       add_element         		(Cluster);
  void 			 construire_voisinages	();
  void 			 clear_overlap					(ostream&); 
  void			 to_svg									();
  
  static HitList load_hitFile			  (string,int, MemeObj, Dictionnaire);
  //comptage motifs individuels
  map<string, set<string> >			motifs_par_seq					();
  void													output_motifs_par_seq		(map<string, set<string> >, ostream&);
  
  //comptage couples
  map<string, map<string, int> >								compter_couples							();
  map<string, set<pair<string, string> > >      presence_couples_par_seq		();
  void			 																		output_couples							(map<string, map<string, int> >, ostream&);
  void			 																		output_presence_couples_seq	(map<string, set<pair<string, string> > >, ostream&);
  
  //comptage pseudo palindromes
  map<string, int>			compter_pseudo_pal					();
  void			 						output_pseudo_pal		  			(map<string, int>, ostream&);

  private:
  
  int 						taille_voisinage;
  Affichage 			affichage;
  Dictionnaire		dictionnaire;
  MemeObj					meme;
  set<string>			motifs;
  list<Cluster> 	collection_clusters;
  
  // sequence > brin > motif > (start,stop) > iterateur de la liste
  map<string, map<string, map<string, map<pair<int,int>, list<Cluster>::iterator> > > > annuaire;
  
  // sequence > (start,stop) > brin > motif > iterateur de la liste
  map<string, multimap<pair<int,int>, map<string, map<string, list<Cluster>::iterator> > > > annuaire2;
  
  map<string, vector<list<Cluster>::iterator> > mtf_index;
};

ostream &operator<<(ostream &flux, HitList const& d);

#endif
