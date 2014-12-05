#include "cluster.h"

Cluster::Cluster(string n_seq,string r_mtf,string n_mtf,int i_mtf,string strd,
	int strt,int stp) :
	
	nom_sequence(n_seq),
	ref_motif(r_mtf),
	nom_motif(n_mtf),
	id_motif(i_mtf),
	strand(strd),
	start(strt),
	stop(stp) 
	{}
	
ostream &operator<<(ostream &flux, Cluster const& c) {

  c.afficher(flux);
  return flux;
}

bool operator==(Cluster const& c1, Cluster const& c2) {

	return c1.egal(c2);
}

bool Cluster::egal(Cluster c) const {

	return ( (nom_sequence == c.nom_sequence) &&
					 (nom_motif == c.nom_motif) &&
					 (start == c.start) &&
					 (stop == c.stop) &&
					 (strand == c.strand) );
}

void Cluster::afficher(ostream& flux) const {

	flux << "Sequence : " << nom_sequence << endl
			 << "Motif : " << nom_motif << endl
			 << "Brin : " << strand << endl
			 << "Start : " << start << endl
			 << "Stop : " << stop << endl
			 << "Voisinage : " << endl;
			 
	for(auto v_it=voisinage.begin(); v_it!=voisinage.end(); v_it++)
	{
		flux << "\t" << (*v_it)->nom_motif << "\t" << (*v_it)->strand
				 << "\t" << (*v_it)->start << "\t" << (*v_it)->stop
				 << endl;
	}
	
	flux << "Doublons : " << endl;
			 
	for(auto do_it=doublons.begin(); do_it!=doublons.end(); do_it++)
	{
		flux << "\t" << (*do_it)->nom_motif << "\t" << (*do_it)->strand
				 << "\t" << (*do_it)->start << "\t" << (*do_it)->stop
				 << endl;
	}
	
	flux << "Chevauchants : " << endl;
			 
	for(auto ch_it=chevauchants.begin(); ch_it!=chevauchants.end(); ch_it++)
	{
		flux << "\t" << (*ch_it)->nom_motif << "\t" << (*ch_it)->strand
				 << "\t" << (*ch_it)->start << "\t" << (*ch_it)->stop
				 << endl;
	}
	
	flux << endl;
}

// Si il y a déjà un exemplaire du motif dans le voisinage, on n'en rajoute pas
// d'autre
void Cluster::add_vois(list<Cluster>::iterator l_it) {

	bool deja_present=false;
	for(auto v_it=voisinage.begin(); v_it!=voisinage.end(); v_it++) {
		
		if((*v_it)->nom_motif == l_it->nom_motif) deja_present=true;
	}
	if(!deja_present) voisinage.push_back(l_it);
	else doublons.push_back(l_it);
}
