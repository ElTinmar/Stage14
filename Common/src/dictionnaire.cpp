#include "dictionnaire.h"

extern int yyparse(Dictionnaire& dico);
extern FILE* yyin;


ostream &operator<<(ostream &flux, Dictionnaire const& d) {

  d.afficher(flux);
  return flux;
}

Dictionnaire::Dictionnaire() {

	output_format = 0;
}

void Dictionnaire::to_fasta(ostream& flux) const {

  for(auto it=dictionnaire.begin(); it != dictionnaire.end(); it++)
  {
    it->to_fasta(flux);
  }
}

string Dictionnaire::valid_filename(string filename) {

  //Supprimer ","
  size_t pos = filename.find(",");
  while(pos!=string::npos)
  {
    filename.replace(pos,1,"_");
    pos = filename.find(",");
  }
  
  //Supprimer "/"
  pos = filename.find("/");
  while(pos!=string::npos)
  {
    filename.replace(pos,1,"_");
    pos = filename.find("/");
  }
  return filename;
}

void Dictionnaire::afficher(ostream& flux) const {

  for(auto it=dictionnaire.begin(); it != dictionnaire.end(); it++)
  {
    flux << *it;
  }
}

void Dictionnaire::cat(ostream& flux) const {

	string concatenate_seq;
	for(auto it=dictionnaire.begin(); it != dictionnaire.end(); it++)
  {
    concatenate_seq.append(it->sequence);
  }
  
  flux << ">cat\n"
  		 << concatenate_seq
  		 << "\n";
}

void Dictionnaire::statistics(ostream& flux) const {

  flux << "Nombre de séquences total : " << dictionnaire.size() << endl;
}

void Dictionnaire::sequence_length_distrib(ostream& flux) const {

	flux << "sequence\tsize" << endl;
	for (auto it=dictionnaire.begin(); it!=dictionnaire.end(); it++)
	{
		flux << (*it).id << "\t" << (*it).seq_length() << endl;
	}
}

map<string, int> Dictionnaire::sequence_length_distrib() {
	
	map<string, int> result;
	for (auto it=dictionnaire.begin(); it!=dictionnaire.end(); it++)
	{
		result[it->id]=it->seq_length();
	}
	return result; 
}

void Dictionnaire::add_element(Element e) {

	e.valid_id();
  dictionnaire.push_back(e);
}

Dictionnaire Dictionnaire::load_fasta(string filename) {

	Dictionnaire dico;
	FILE *input = fopen(filename.c_str(),"r");
	
  if (input == NULL)
  {
    cerr << "impossible d'ouvrir le fichier %s en lecture\n" << filename;
    exit(EXIT_FAILURE);
  }
  
  yyin = input;  // fichier en entree du parser
  int parser_result;
  
  //envoyer probleme comme parametre du parser
  parser_result = yyparse(dico); 
  fclose(input);
  
  if (parser_result != 0)
  {
    cerr << "parser error :\n";
  }
 
	return dico;
}

vector<string> Dictionnaire::list_ids() const {

	vector<string> result;
	for(auto it=dictionnaire.begin(); it != dictionnaire.end(); it++)
  {
    result.push_back(it->id);
  }
  return result;
}

string Dictionnaire::get_sequence	(string id_sequence) {
	
	for(auto it=dictionnaire.begin(); it != dictionnaire.end(); it++) {
		if(it->id == id_sequence) {
			return it->sequence;
		}
	}
}

