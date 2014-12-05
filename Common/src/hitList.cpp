#include "hitList.h"

HitList::HitList(int t, MemeObj m, Dictionnaire d) { 
	
	taille_voisinage = t; 
	dictionnaire = d; 
	meme = m; 
	motifs = m.motifs_altname();
	affichage.add_l_sequences(dictionnaire.sequence_length_distrib());
	affichage.add_name_index(meme.name_index());
}

ostream &operator<<(ostream &flux, HitList const& d) {

  d.afficher(flux);
  return flux;
}

void HitList::afficher(ostream& flux) const {

	for(auto it=collection_clusters.begin();
					 it!=collection_clusters.end();
					 it++) {
					 
		flux << *it << "\n";
	}
}

void HitList::to_hitlist(ostream& flux) const {

	for(auto it=collection_clusters.begin();
					 it!=collection_clusters.end();
					 it++) {
					 
		it->to_hit(flux);
	}
}

void HitList::afficher_annuaire(ostream& flux) {

	for (auto annuaire_it = annuaire.begin(); 
					 	annuaire_it != annuaire.end();
					  annuaire_it++) {
					 
		for (auto sequence_it = annuaire_it->second.begin();
							sequence_it != annuaire_it->second.end();
							sequence_it++) {
							
			for (auto brin_it = sequence_it->second.begin();
								brin_it != sequence_it->second.end();
								brin_it++) {
								
				for (auto motif_it = brin_it->second.begin();
									motif_it != brin_it->second.end();
									motif_it++) {

					flux << *motif_it->second;
				
				}
			}
		}
	}
}

void HitList::afficher_annuaire2(ostream& flux) {

	for (auto annuaire_it = annuaire2.begin(); 
					 	annuaire_it != annuaire2.end();
					  annuaire_it++) {
					 
		for (auto sequence_it = annuaire_it->second.begin();
							sequence_it != annuaire_it->second.end();
							sequence_it++) {
							
			for (auto paire_it = sequence_it->second.begin();
								paire_it != sequence_it->second.end();
								paire_it++) {
								
				for (auto brin_it = paire_it->second.begin();
									brin_it != paire_it->second.end();
									brin_it++) {

					flux << brin_it->second->nom_motif << "\t"
							 << brin_it->second->start << "\t"
							 << brin_it->second->stop - brin_it->second->start<< "\n";
				
				}
			}
		}
	}
}

void HitList::clear_overlap(ostream& flux) {

	set<Cluster> clear_cluster_collection;
	Cluster min_cluster;
	
	for (auto annuaire_it = annuaire2.begin(); 
					 	annuaire_it != annuaire2.end();
					  annuaire_it++) {
					 
		for (auto sequence_it = annuaire_it->second.begin();
							sequence_it != annuaire_it->second.end();
							sequence_it++) {
							
			for (auto paire_it = sequence_it->second.begin();
								paire_it != sequence_it->second.end();
								paire_it++) {
								
				for (auto brin_it = paire_it->second.begin();
									brin_it != paire_it->second.end();
									brin_it++) {
					
					min_cluster = *(brin_it->second);
					
					for (auto chev_it = brin_it->second->chevauchants.begin();
										chev_it != brin_it->second->chevauchants.end();
										chev_it++) {
						
						if( (*chev_it)->nom_motif == brin_it->second->nom_motif ) {
							
							if( (*chev_it)->q_value <= brin_it->second->q_value ) {
							
								min_cluster = **chev_it;
							}
						}
					}

					clear_cluster_collection.insert(min_cluster);
				}
			}
		}
	}
	
	flux << "#pattern name	sequence name	start	stop	strand	score	p-value	q-value	matched sequence\n";
	for(auto it = clear_cluster_collection.begin();
					 it != clear_cluster_collection.end();
					 it++) {
					 
		 it->to_hit(flux);				 
	}
}

void HitList::add_element(Cluster c) {

  collection_clusters.push_back(c);
  auto it=collection_clusters.end();
  it--;
  
  annuaire[c.nom_sequence][c.strand][c.nom_motif][make_pair(c.start,c.stop)]=it;
  
  map<string, map<string, list<Cluster>::iterator> > temp_map;
  temp_map[c.strand][c.nom_motif]=it;
  annuaire2[c.nom_sequence].insert(make_pair(make_pair(c.start,c.stop), temp_map));
  mtf_index[c.nom_motif].push_back(it);
  
  affichage.add_motif(c,dictionnaire);
}

void HitList::construire_voisinages() {

	for(auto c_it = collection_clusters.begin();
					 c_it != collection_clusters.end();
					 c_it++) {
					 
		for(auto aSeq_it = annuaire[c_it->nom_sequence].begin();
		 				 aSeq_it != annuaire[c_it->nom_sequence].end();
		 				 aSeq_it++) {
		 				 
		 	for(auto aStd_it = aSeq_it->second.begin();
		 				 aStd_it != aSeq_it->second.end();
		 				 aStd_it++) {
		 				 	 
		 		for(auto aStart_it = aStd_it->second.begin();
		 				 	 	 aStart_it != aStd_it->second.end();
		 				 	   aStart_it++) {
			 	 	
			 	 	//Ne pas rajouter le motif courant dans ses propres voisins
			 		if( aStart_it->second != c_it ) {
			 		 
						if( aStart_it->first.first < c_it->start ) {
					
							//Pas de chevauchement : on vérifie la distance
							if( aStart_it->first.second < c_it->start ) {
						
								if( c_it->start - aStart_it->first.second <= taille_voisinage ) {
								
									c_it->add_vois(aStart_it->second);
								}
							}
							//Chevauchement : on ne prend pas en compte
							else {
					
								c_it->chevauchants.push_back(aStart_it->second);
							}
						}
					
						else if ( aStart_it->first.first >= c_it->start ) {
				
							//Pas de chevauchement
							if( aStart_it->first.first > c_it->stop ) {
						
								if( aStart_it->first.first - c_it->stop  <= taille_voisinage ) {
								
									c_it->add_vois(aStart_it->second);
								}
							}
							//Chevauchement : on ne prend pas en compte
							else {
						
								c_it->chevauchants.push_back(aStart_it->second);
							}
						}
					}		
				}
			}
		}
	}
}

map<string, set<string> >	HitList::motifs_par_seq() {
	
	map<string, set<string> > results;
	for(auto it=collection_clusters.begin(); it!=collection_clusters.end(); it++){
		results[it->nom_sequence].insert(it->nom_motif);
	}
	return(results);
}

void HitList::output_motifs_par_seq(map<string, set<string> > table, ostream& flux) {
	
	for(auto it=table.begin(); it!=table.end(); it++){
		flux << it->first << "\t";
		for(auto mtf_it=it->second.begin(); mtf_it!=it->second.end(); mtf_it++) {
		
			flux << *mtf_it << "\t";
		}
		flux << endl;
	}
}

map<string, map<string, int> > HitList::compter_couples() {

	map<string, map<string, int> > results;
	for(auto mtf_it1=motifs.begin(); mtf_it1!=motifs.end(); mtf_it1++) {
	
		for(auto mtf_it2=motifs.begin(); mtf_it2!=motifs.end(); mtf_it2++) {
		
			for(auto lst_it = mtf_index[*mtf_it1].begin(); 
							 lst_it != mtf_index[*mtf_it1].end();
							 lst_it++) {
							 
				for(auto vsn_it = (*lst_it)->voisinage.begin(); 
							 	 vsn_it != (*lst_it)->voisinage.end();
							 	 vsn_it++) {
					
					//Si on trouve motif2 dans le voisinage de motif1		 
					if((*vsn_it)->nom_motif == *mtf_it2) {
					
						results[*mtf_it1][*mtf_it2]++; 
					}
				}		
			}
		}
	}
	return results;
}

map<string, set<pair<string, string> >	> HitList::presence_couples_par_seq() {

	map<string, set<pair<string, string> >	> results;
	for(auto a_it=annuaire.begin(); a_it!=annuaire.end(); a_it++) {
	
		for(auto mtf_it1=motifs.begin(); mtf_it1!=motifs.end(); mtf_it1++) {
	
			for(auto mtf_it2=motifs.begin(); mtf_it2!=motifs.end(); mtf_it2++) {
		
				for(auto lst_it = mtf_index[*mtf_it1].begin(); 
								 lst_it != mtf_index[*mtf_it1].end();
								 lst_it++) {
					
					//si m1 est dans la bonne séquence
					if((*lst_it)->nom_sequence == a_it->first) {
								 
						for(auto vsn_it = (*lst_it)->voisinage.begin(); 
									 	 vsn_it != (*lst_it)->voisinage.end();
									 	 vsn_it++) {
					
							//Si on trouve motif2 dans le voisinage de motif1		 
							if((*vsn_it)->nom_motif == *mtf_it2) {
					
								results[(*lst_it)->nom_sequence].insert(make_pair(*mtf_it1,*mtf_it2)); 
							}
						}		
					}
				}
			}
		}
	}
	return results;
}

void HitList::output_presence_couples_seq(map<string, set<pair<string, string> > > table, ostream& flux) {

	for(auto seq_it=table.begin(); seq_it != table.end(); seq_it++) {
		
		flux << seq_it->first << "\t";
		for(auto pair_it=seq_it->second.begin(); pair_it!=seq_it->second.end(); pair_it++) {
		
			flux << pair_it->first << "-" << pair_it->second << "\t";
		}
		flux << endl;
	}
}

void HitList::output_couples(map<string, map<string, int> > table, ostream& flux) {

	for(auto it1=motifs.begin(); it1!=motifs.end(); it1++) {
  
  	 flux << *it1 << "\t";
  }
  flux << endl;
  
  for(auto it1=motifs.begin(); it1!=motifs.end(); it1++) {
  
  	flux << *it1 << "\t";
  	for(auto it2=motifs.begin(); it2!=motifs.end(); it2++) {
  	 
  		flux << table[*it1][*it2] << "\t";
  	}
  	flux << endl;
  }
}

// Compter pour chaque motif le nombre d'occurences qui chevauche entièrement ou
// partiellement une ou plusieurs occurences du même motif sur le brin 
// complémentaire
map<string, int> HitList::compter_pseudo_pal() {
	
	map<string, int> results;
	for(auto c_it=collection_clusters.begin();
	 				 c_it!=collection_clusters.end();
	 				 c_it++) {
	 				 
		for(auto ch_it=c_it->chevauchants.begin(); 
						 ch_it!=c_it->chevauchants.end(); 
						 ch_it++) {
						 
			if( ((*ch_it)->nom_motif == c_it->nom_motif) && 
					((*ch_it)->strand != c_it->strand) ) {
					
				results[c_it->nom_motif]++;
				break; // au moins un 
			}
		}
	}
	return results;
}

void HitList::output_pseudo_pal(map<string, int> table, ostream& flux) {
  
  flux << "\t" << "compte" << endl;
  for(auto it1=motifs.begin(); it1!=motifs.end(); it1++) {
  
  	flux << *it1 << "\t" << table[*it1] << endl;
  }
}

void HitList::to_svg() {

  affichage.svg();
}

HitList HitList::load_hitFile(string filename,int t, MemeObj m, Dictionnaire d) {

	HitList hitList(t,m,d);
	
	ifstream hit_file(filename, std::ifstream::in);
  if (!hit_file.is_open())
  {
    cerr << "erreur à l'ouverture de " << filename << " en lecture" << endl;
    exit(EXIT_FAILURE);
  }
  
  string junk;
  string motif_ref;
	string sequence_name;
  int start;
  int stop;
  string strand;
  double score;
  double q_value;
  string matched_seq;
  
  while(hit_file >> motif_ref)
  {
    if(motif_ref.find("#")!=string::npos)
    {
    	getline(hit_file, junk);
    	continue;
    }
    
    hit_file >> sequence_name >> start >> stop >> strand >> score >> q_value >> matched_seq;
    getline(hit_file, junk);
    
    Cluster c(
    	sequence_name, 
    	motif_ref, 
    	m.ref_to_altname(motif_ref), 
    	m.ref_to_index(motif_ref), 
    	strand, 
    	start, 
    	stop,
    	score,
    	q_value,
    	matched_seq
    );
    
    hitList.add_element(c);
  }
  
  hitList.construire_voisinages();
  
  hit_file.close();
  
  return hitList;
}
