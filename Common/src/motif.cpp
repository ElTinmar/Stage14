#include "motif.h"

Motif::Motif(string id, Pwm p) {

	identifier = id;
	alt_name = "";
	matrice = p;
	url = "";
}

Motif::Motif(string id, Pwm p, string u) {

	identifier = id;
	alt_name = "";
	matrice = p;
	url = u;
}

Motif::Motif(string id, string altn, Pwm p) {

	identifier = id;
	alt_name = altn;
	matrice = p;
	url = "";
}

Motif::Motif(string id, string altn, Pwm p, string u) {

	identifier = id;
	alt_name = altn;
	matrice = p;
	url = u;
}
	
ostream &operator<<(ostream &flux, Motif const& m) {

  m.afficher(flux);
  return flux;
}

void Motif::afficher(ostream& flux) const {

	flux << "Identifiant : " << identifier << endl
			 << "Nom alternatif : " << alt_name << endl
			 << "Matrice : " << endl << matrice << endl
			 << "URL : " << url << endl;
}

double Motif::entropy() {

	return matrice.entropy();
}

int	Motif::longueur() {

	return matrice.longueur();
}
