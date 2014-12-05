#ifndef CLUSTER_H
#define CLUSTER_H
#include <iostream>
#include <list>
#include <vector>
#include <map>
#include <tuple>
#include <string>
#include <cstdlib>
#include <fstream>

using namespace std;

class Cluster
{
	friend class Dom;
	friend class Affichage;
	
  public:
  
  							Cluster				(string,string,string,int,string,int,int);
  void      		afficher      (ostream&) const;
  bool          egal          (Cluster) const;
  void					add_vois			(list<Cluster>::iterator);
  void					to_svg				(ostream&);
  
  private:
  
  string 				nom_sequence;
  string				ref_motif;
  string 				nom_motif;
  int						id_motif;
  string				strand;
  int 					start;
  int 					stop;
  
  vector<list<Cluster>::iterator > voisinage;
  vector<list<Cluster>::iterator > doublons;
  vector<list<Cluster>::iterator > chevauchants;
};

ostream &operator<<(ostream &flux, Cluster const& d);
bool operator==(Cluster const&, Cluster const&);

#endif
