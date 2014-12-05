#include "cluster.h"

Cluster::Cluster() {}

Cluster::Cluster(string n_seq,string r_mtf,string n_mtf,int i_mtf,string strd,
	int strt,int stp,double scr, double q_val, string mseq) :
	
	nom_sequence(n_seq),
	ref_motif(r_mtf),
	nom_motif(n_mtf),
	id_motif(i_mtf),
	strand(strd),
	start(strt),
	stop(stp),
	score(scr),
	q_value(q_val),
	matched_seq(mseq)
	{}
	
ostream &operator<<(ostream &flux, Cluster const& c) {

  c.afficher(flux);
  return flux;
}

bool operator==(Cluster const& c1, Cluster const& c2) {

	return c1.egal(c2);
}

bool operator<(Cluster const& c1, Cluster const& c2) {

	return c1.inf(c2);
}

bool Cluster::egal(Cluster c) const {

	return ( (nom_sequence == c.nom_sequence) &&
					 (nom_motif == c.nom_motif) &&
					 (start == c.start) &&
					 (stop == c.stop) &&
					 (strand == c.strand) );
}

bool Cluster::inf(Cluster c) const {

	return ( ( (nom_sequence == c.nom_sequence) && (start < c.start) ) ||
					 ( (nom_sequence == c.nom_sequence) && 
					 	 (start == c.start) &&
					 	 lexicographical_compare(
					 			nom_motif.begin(),nom_motif.end(),
					 			c.nom_motif.begin(),c.nom_motif.end() ) 
					 	) ||
					 ( lexicographical_compare(
					 			nom_sequence.begin(),nom_sequence.end(),
					 			c.nom_sequence.begin(),c.nom_sequence.end() ) )
				);
}

void Cluster::afficher(ostream& flux) const {

	flux << "Sequence : " << nom_sequence << endl
			 << "Motif : " << nom_motif << endl
			 << "Brin : " << strand << endl
			 << "Start : " << start << endl
			 << "Stop : " << stop << endl
			 << "Score : " << score << endl
			 << "q-value : " << q_value << endl
			 << "Matched Sequence : " << matched_seq << endl 
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

void Cluster::to_hit(ostream& flux) const {

	flux << ref_motif << "\t"
			 << nom_sequence << "\t"
			 << start << "\t"
			 << stop << "\t"
			 << strand << "\t"
			 << score << "\t"
			 << q_value << "\t"
			 << matched_seq << endl;
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
