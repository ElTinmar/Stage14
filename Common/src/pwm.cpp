#include "pwm.h"

ostream &operator<<(ostream &flux, Pwm const& pwm) {

  pwm.afficher(flux);
  return flux;
}

vector<double>& Pwm::operator[](char letter) {

	return matrice[letter];
}

void Pwm::afficher(ostream& flux) const {

	flux << "Longueur : " << w << endl
			 << "Card Alphabet : " << alength << endl;
			 
  for(auto L_it = matrice.begin(); L_it != matrice.end(); L_it++) {
  	flux << L_it->first << " | ";
  	for(auto P_it = L_it->second.begin(); P_it != L_it->second.end(); P_it++) {
  		flux << *P_it << "\t";
  	}
  	flux << endl;
  }
}

void Pwm::set_w(int width) {

	w = width;
}

void Pwm::set_alength(int card) {

	alength = card;
}

void Pwm::set_mat(map<char, vector<double> > mat) {

	matrice = mat;
}

int Pwm::longueur() {

	return w;
}

double Pwm::entropy() {
	
	double h = 0;	
	double h_i = 0;
	for(int i=0; i<w; i++) {
		h_i = 2;
		for(auto it=matrice.begin(); it!=matrice.end(); it++) {
			double p_i = it->second[i];
			if(p_i != 0) h_i = h_i + (p_i*log2(p_i));
		}
		h = h + h_i;
	}
	return h/w;
}
