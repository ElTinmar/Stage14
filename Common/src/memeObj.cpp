#include "memeObj.h"

extern int mmparse(memeObj& meme);
extern FILE* mmin;

ostream &operator<<(ostream &flux, memeObj const& mo) {

  mo.afficher(flux);
  return flux;
}

void memeObj::afficher(ostream& flux) const {

	flux << "Version : " << version << endl
			 << "Alphabet : ";
			 
	for(auto A_it=alphabet.begin(); A_it!=alphabet.end(); A_it++) {
		flux << *A_it << " ";
	}
	
	flux << endl
			 << "Strands : ";
			 
	for(auto S_it=strands.begin(); S_it!=strands.end(); S_it++) {
		flux << *S_it << " ";
	}
	
	flux << endl
			 << "Background Frequencies : " ;
		
	for(auto B_it=bg_frequencies.begin(); B_it!=bg_frequencies.end(); B_it++) {
		flux << B_it->first << " " << B_it->second << "\t";
	}		 
			 
	flux << endl
			 << "Motifs : ";
			 
	for(auto M_it=motifs.begin(); M_it!=motifs.end(); M_it++) {
		flux << *M_it << endl;
	}
}

void memeObj::set_version(string v) {

	version = v;
}

void memeObj::set_alphabet(string a) {

	for(auto it=a.begin(); it!=a.end(); it++) {
		alphabet.insert(*it);
	}
}

void memeObj::add_strand(string s) {

	strands.insert(s);
}

void memeObj::set_bg_frq(map<string, double> m) {
	
	bg_frequencies=m;
}

void memeObj::set_motifs(vector<Motif> v) {

	motifs = v;
}

memeObj memeObj::load_meme(string filename) {

	memeObj meme;
	FILE *input = fopen(filename.c_str(),"r");
	
  if (input == NULL)
  {
    cerr << "impossible d'ouvrir le fichier %s en lecture\n" << filename;
    exit(EXIT_FAILURE);
  }
  
  mmin = input;  // fichier en entree du parser
  int parser_result;
  
  //envoyer probleme comme parametre du parser
  parser_result = mmparse(meme); 
  fclose(input);
  
  if (parser_result != 0)
  {
    cerr << "parser error : " << filename << "\n";
    exit(EXIT_FAILURE);
  }
 
	return meme;
}

const set<char>& memeObj::get_alphabet() const {
	
	return alphabet;
}

void memeObj::entropy(ostream& flux) {
	
	//sort by alt_name
	map<string, pair<double, int> > resultat;
	for(auto it=motifs.begin(); it!=motifs.end(); it++) {
		resultat[it->alt_name]=make_pair(it->entropy(), it->longueur());
	}
	
	flux << "\t" << "entropy" << "\t" << "longueur" << endl;
	for(auto it=resultat.begin(); it!=resultat.end(); it++) {
		flux << it->first << "\t" 
				 << it->second.first << "\t" 
				 << it->second.second << endl;
	}
}

string memeObj::ref_to_altname(string ref) {

	for(auto it=motifs.begin(); it != motifs.end(); it++) {
		if(it->identifier == ref) return it->alt_name;
	}
	return("");
}

int memeObj::ref_to_index(string ref) {

	for(auto it=motifs.begin(); it != motifs.end(); it++) {
		if(it->identifier == ref) return (distance(motifs.begin(), it)+1);
	}
	return(-1);
}

map<string, int> memeObj::name_index() {

	map<string, int> result;
	for(auto it=motifs.begin(); it != motifs.end(); it++) {
		result[it->alt_name]=(distance(motifs.begin(), it)+1);
	}
	return result;
}

set<string> memeObj::motifs_altname() {

	set<string> result;
	for(auto it=motifs.begin(); it != motifs.end(); it++) {
		result.insert(it->alt_name);
	}
	return result;
}

