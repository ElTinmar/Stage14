#ifndef Affichage_H
#define Affichage_H
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
#include "dictionnaire.h" 

// definitions pour la sortie SVG
#define MARGIN_TOP 30
#define MARGIN_BOTTOM 30
#define MARGIN_LEFT 30
#define MARGIN_RIGHT 50
#define V_SEQUENCE 15
#define X_SEQUENCE 50
#define V_MOTIF	10
#define V_STRING 5
#define H_STRING 80
#define TRAIT_SEQ 2

#define DEBUG 0

using namespace std;

struct _sequence {
  string name;
  string sequence;
  int x;
  int y;
  int l;
  int v;
} ;

struct _motif {
  int c;
  int x;
  int y;
  int l;
  int v;
  list<_sequence>::iterator seq;
} ;

class Affichage
{
  public:
  
  void       									svg 	           		();
  void       									debug           		(ostream&);
  void												add_l_sequences			(map<string, int>);
  void												add_name_index			(map<string, int>);
  void			 									add_motif						(Cluster, Dictionnaire);
  void												resize							(int);
  static bool 								superposition       (_motif const&, _motif const&);
  bool												en_superposition		(_motif);
  list<_sequence>::iterator 	add_sequence				(string,string);
  list<_sequence>::iterator		push_sequence				(_sequence);

  private:
  
  map<string, int>	l_sequences;
  map<string, int>  name_index;
  list<_sequence> 	sequences;
  list<_motif> 			motifs;
  int								hauteur;
  int								largeur;
  
};

#endif
