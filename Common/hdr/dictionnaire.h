#ifndef DICTIONNAIRE_H
#define DICTIONNAIRE_H
#include <iostream>
#include <list>
#include <vector>
#include <map>
#include <tuple>
#include <string>
#include <cstdlib>
#include <fstream>
#include "element.h"

using namespace std;

class Dictionnaire
{
	friend class Affichage;
	
  public:
  							Dictionnaire				();
  void      		afficher            (ostream&) const;
  void      		to_fasta    				(ostream&) const;
  void      		add_element         (Element);
  void					statistics					(ostream&) const;
  void					sequence_length_distrib (ostream&) const;
  void					cat									(ostream&) const;
  
  string				get_sequence				(string);
  
  map<string,int>	sequence_length_distrib ();
  vector<string> list_ids						() const;
  
  static string valid_filename  (string);
  static Dictionnaire load_fasta(string);
  
  //public attributes for option gestion
  int output_format; //1:classique par défaut, 2:vista
  
  private:
  // (gene ID, expression domain, organism, chromosom, position on the chromosom, sequence);
  list<Element> dictionnaire;
};

ostream &operator<<(ostream &flux, Dictionnaire const& d);

#endif
