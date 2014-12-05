#include "affichage.h"

void Affichage::add_l_sequences(map<string, int> l) {

	l_sequences = l;
}

void Affichage::add_name_index(map<string, int> n){

	name_index = n;
}
  
bool Affichage::superposition(_motif const& m1,_motif const& m2) {
	
	if( m1.seq == m2.seq && m1.v > 0 && m2.v < 0 ) return false;
	if( m1.seq == m2.seq && m1.v < 0 && m2.v > 0 ) return false;
	if( m1.seq == m2.seq && m1.v > 0 && m2.v > 0 ) {
		return (
			!(
					(m1.x > m2.x + m2.l) ||
					(m1.y > m2.y + m2.v) ||
					(m2.x > m1.x + m1.l) ||
					(m2.y > m1.y + m1.v)
			) 
		);
	}
	if( m1.seq == m2.seq && m1.v < 0 && m2.v < 0 ) {
		return (
			!(
					(m1.x > m2.x + m2.l) ||
					(m1.y > m2.y - m2.v) ||
					(m2.x > m1.x + m1.l) ||
					(m2.y > m1.y - m1.v)
			) 
		);
	}
}

list<_sequence>::iterator	Affichage::push_sequence(_sequence s) {

	for(auto it = sequences.begin(); it != sequences.end(); it++) {
		if(it->name == s.name) return it;
	}
	sequences.push_back(s);
	auto res_it = sequences.end();
	res_it--;
	return res_it;
}

list<_sequence>::iterator Affichage::add_sequence(string seq_name, string sequence) {

	_sequence new_seq;
	new_seq.name = seq_name;
	new_seq.sequence = sequence;
	new_seq.x = X_SEQUENCE;
  new_seq.y = 0;
  new_seq.l = l_sequences[seq_name];
  new_seq.v = V_SEQUENCE;
  return push_sequence(new_seq);
}

bool Affichage::en_superposition(_motif m) {

	for(auto it = motifs.begin(); it != motifs.end(); it++)
	{
		if( superposition(m,*it) ) { 
			return true;
		}
	}
	return false;
}

void Affichage::add_motif(Cluster c, Dictionnaire d) {

	auto seq = add_sequence(c.nom_sequence, d.get_sequence(c.nom_sequence));
	
	_motif new_mtf;
	new_mtf.c = c.id_motif;
	new_mtf.x = c.start;
	new_mtf.l = c.stop - c.start;
	if (c.strand == "+") new_mtf.v = -V_MOTIF;
	else new_mtf.v = V_MOTIF;
	new_mtf.seq = seq;
	
	new_mtf.y = 0;
	int decalage = 0;
	while( en_superposition(new_mtf) ) {
		if(new_mtf.v < 0) new_mtf.y -= V_MOTIF + 1;
		else new_mtf.y += V_MOTIF + 1;
		decalage += V_MOTIF + 1;
	}
	
	(new_mtf.seq)->v = max((new_mtf.seq)->v, V_SEQUENCE + decalage);
	
	motifs.push_back(new_mtf);
}

void Affichage::svg() {
	
	ofstream svg_file;
	int ncol_legende = 0;
	int nrow_legende = 0;
	int current_pos = 0;
	
	for(auto s_it = sequences.begin(); s_it != sequences.end(); s_it++) {
	
		svg_file.open(string(s_it->name) + ".svg");
		
		if(!svg_file.is_open()) {
			cerr << "Une erreur est survenue lors de l'ouverture du fichier " << string(s_it->name) + ".svg";
			exit(EXIT_FAILURE);
		}
		
		largeur = MARGIN_LEFT + s_it->l + MARGIN_RIGHT;
		
		// initialisation
		ncol_legende = ( largeur - MARGIN_RIGHT - MARGIN_LEFT ) / H_STRING;
		nrow_legende = ( name_index.size() / ncol_legende ) + 1;
	
		// Commencer le tracé de la figure en bas de la légende
		current_pos = V_STRING * nrow_legende;
		
		current_pos = V_STRING * nrow_legende;
		s_it->y = current_pos + s_it->v;
		current_pos += 2*s_it->v; 
	
		hauteur = MARGIN_TOP + 
							MARGIN_BOTTOM + 
							current_pos;
	
		// En tête SVG
		svg_file << "<?xml version=\"1.0\" standalone=\"no\"?>\n"
		     << "\t<!DOCTYPE svg PUBLIC \"-//W3C//DTD SVG 1.1//EN\" \"http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd\">\n"
		     << "\t\t<svg width=\"" << largeur << "\" height=\"" << hauteur << "\" version=\"1.1\" xmlns=\"http://www.w3.org/2000/svg\" >\n";
		     
		// Definition des motifs
		svg_file << "\t\t\t<defs>\n"
		     << "\t\t\t\t<rect id=\"Motif\" height=\"" << V_MOTIF << "\" stroke=\"none\" fill=\"black\"/>\n"
		     << "\t\t\t</defs>\n";

		// Tracé de la légende
		int mtf_count = 0;
		int ncol_legende = ( largeur - MARGIN_RIGHT - MARGIN_LEFT ) / H_STRING;
		int nrow_legende = ( name_index.size() / ncol_legende ) + 1;
	
		svg_file << "\t\t\t<g id=\"Legende\">\n";
		for(auto l_it = name_index.begin(); l_it != name_index.end(); l_it++)
		{
		  svg_file << "\t\t\t\t<text text-anchor=\"begin\" font-family=\"Verdana\" font-size=\"" << V_STRING << "\" fill=\"black\" y=\""
		       << V_STRING * ( (mtf_count % nrow_legende) + 1 )
		       << "\" x=\""
		       << H_STRING * (mtf_count / nrow_legende) + MARGIN_LEFT
		       << "\" >"
		       << l_it->second << " : " << l_it->first
		       << "</text>\n";
		  mtf_count++;
		}
		svg_file << "\t\t\t</g>\n";
		
		// Tracé des séquences
  
  	//Commencer le groupe et ecrire le nom de la séquence
    svg_file << "\t\t\t<g>\n"
         << "\t\t\t\t<text x=\"" << MARGIN_LEFT << "\" y=\"" 
         << s_it->y 
         << "\" text-anchor=\"end\" font-family=\"Verdana\" font-size=\"" << V_STRING << "\" fill=\"black\" >" 
         << s_it->name
         << "</text>\n";
    
    //tracer le trait représentant la séquence
    svg_file << "\t\t\t\t<path d=\"M " << s_it->x << " " 
         << s_it->y 
         << " H "
         << s_it->x + l_sequences[s_it->name]
         << "\" fill=\"black\" stroke=\"black\" stroke-width=\"" << TRAIT_SEQ << "\" />\n"; 
         
    //écrire la séquence
    svg_file << "\t\t\t\t<text x=\"" << s_it->x << "\" y=\"" 
         << s_it->y+0.5
         << "\" fill=\"white\" font-size=\"1.66235\" font-family=\"Monospace\">"
         << s_it->sequence
         << "</text>\n" 
         << "\t\t\t</g>\n";
  
		// Tracé des motifs
		for(auto m_it = motifs.begin(); m_it != motifs.end(); m_it++) {
			if(m_it->seq->name == s_it->name) {
		
				// Rectangle
				svg_file << "\t\t\t<rect transform=\"translate(" 
						 << m_it->seq->x + m_it->x 
						 << ",";		
						 
				if(m_it->v < 0) {
					svg_file << m_it->seq->y - 2 + m_it->y + m_it->v;
				}
				else {
					svg_file << m_it->seq->y + 2 + m_it->y;
				}
			
				svg_file << ")\" height=\"" << V_MOTIF << "\" stroke=\"none\" fill=\"black\"  x=\"0\" y=\"0\" width=\""
						 << m_it->l
						 << "\" />\n";
						 
				// Index du motif
				svg_file << "\t\t\t<text x=\""
				     << m_it->seq->x + m_it->x + m_it->l/2
				     << "\" y=\""
				     << m_it->seq->y + 2 + m_it->y + m_it->v/2
				 		 << "\" text-anchor=\"middle\" font-family=\"Verdana\" font-size=\"" << V_STRING << "\" fill=\"white\" >" 
				     << m_it->c
				     << "</text>\n";
			}
		}
		
		svg_file << "\t\t</svg>\n\n";
		svg_file.close();
	}
}

void Affichage::resize(int largeur) {
 	//TODO
}

void Affichage::debug(ostream& flux) {

	//init();
	flux << "longueurs : " << endl << endl;
	for(auto l_it = l_sequences.begin(); l_it != l_sequences.end(); l_it++) {
		flux << l_it->first << "\t" << l_it->second << endl;
	}
	
	flux << "Sequences : " << endl << endl;
	for(auto s_it = sequences.begin(); s_it != sequences.end(); s_it++) {
		flux << s_it->name << "\t" << s_it->x << "\t" << s_it->y << "\t" << s_it->l << "\t" << s_it->v << endl;
	}
	flux << "Motifs : " << endl << endl;
	for(auto m_it = motifs.begin(); m_it != motifs.end(); m_it++) {
		flux << m_it->c << "\t" << m_it->x << "\t" << m_it->y << "\t" << m_it->l << "\t" << m_it->v << "\t" << (m_it->seq)->name << endl;
	}
}
