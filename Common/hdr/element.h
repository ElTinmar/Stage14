#ifndef ELEMENT_H
#define ELEMENT_H
#include <iostream>
#include <list>
#include <vector>
#include <deque>
#include <map>
#include <tuple>
#include <string>
#include <cstdlib>
#include <fstream>

using namespace std;

class Element
{
	friend class Dictionnaire;
	
	public:
					Element						(string,string);
					Element						(string);
	void		add_seq						(string);
	void    afficher 					(ostream&) const; 
	void    to_fasta  				(ostream&) const; // écrire au format FASTA 80
	void		print_seq60 			(ostream&) const; // rajouter des sauts de ligne tous les 60 caractères
	void		print_seq80 			(ostream&) const; // ------------------------------------ 80 ----------
	void		print_seq120			(ostream&) const; // ------------------------------------ 120 ---------
  void 		valid_id					();
  int 		seq_length				() const;
	
	private:

	string  				id; //nom de la sequence
	string  				sequence;
};

ostream &operator<<(ostream&, Element const&);

#endif
