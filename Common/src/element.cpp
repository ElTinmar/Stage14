#include "element.h"

Element::Element(string i,string seq) :
	id(i),sequence(seq) {}
	
Element::Element(string i) :
	id(i){}
	
void Element::add_seq(string seq) {

	sequence = seq;
}
				
ostream &operator<<(ostream &flux, Element const& e) {

  e.afficher(flux);
  return flux;
}

void Element::print_seq60 (ostream& flux) const {

	int compteur = 0;
  for (auto it=sequence.begin(); it!=sequence.end(); it++)
  {
  	flux << *it;
  	compteur++;
  	if(compteur == 60)
  	{
  		compteur = 0;
  		flux << endl;
  	}
  }
  if(compteur != 0) flux << endl;
}

void Element::print_seq80 (ostream& flux) const {

	int compteur = 0;
  for (auto it=sequence.begin(); it!=sequence.end(); it++)
  {
  	flux << *it;
  	compteur++;
  	if(compteur == 80)
  	{
  		compteur = 0;
  		flux << endl;
  	}
  }
  if(compteur != 0) flux << endl;
}

void Element::print_seq120 (ostream& flux) const {

	int compteur = 0;
  for (auto it=sequence.begin(); it!=sequence.end(); it++)
  {
  	flux << *it;
  	compteur++;
  	if(compteur == 120)
  	{
  		compteur = 0;
  		flux << endl;
  	}
  }
  if(compteur != 0) flux << endl;
}

void Element::valid_id() {

	for(auto it=id.begin(); it!=id.end(); it++)
	{
		if(*it == ' ') 
		{
			id.erase(it);
			it--;
		}
	}
}

int Element::seq_length() const {

	return sequence.length();
}

void Element::afficher(ostream& flux) const {

	flux << "Identifiant : " << id << endl
			 << "Sequence : " << endl;
  
  //Séquence, affichage sur 60 colonnes
  print_seq80(flux);
}

void Element::to_fasta(ostream& flux) const {

  flux << ">" << id << endl;
  //Séquence
  print_seq60(flux);
}

