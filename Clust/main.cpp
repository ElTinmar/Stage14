#include <iostream>
#include <fstream>
#include <vector>
#include <deque>
#include <string>
#include <map>
#include <cstdlib>
#include <string>
#include <cstdio>
#include <limits>
#include <set>
#include "dictionnaire.h"
#include "memeObj.h"
#include "dom.h" 

using namespace std;

void usage(char* nom)
{
  cout << endl
       << "nombre de paramètres incorrects\n"
       << "\t usage: " << nom << " database.meme sequence.fna hit_list.txt taille_voisinage [option = l|c|p|s|i|d]\n\n"
       << "\t database.meme : base de donnees contenant les motifs\n"
       << "\t sequence.fna : sequences au format FASTA\n"
       << "\t hit_list.txt : alignement de motifs par FIMO\n"
       << "\t taille_voisinage : fenêtre dans laquelle rechercher des voisins (en pb)\n"
       << "\t options :\n"
       << "\t\t c compter les couples\n"
       << "\t\t p compter les palindromes dégénérés\n"
       << "\t\t l afficher le voisinage\n"
       << "\t\t s couples présents par séquence\n"
       << "\t\t i motifs individuels présents par séquence\n"
       << "\t\t d afficher les motifs au format SVG\n"
       << "\t\t v afficher dans l'ordre des sequences\n"
       << endl << endl;
       
  exit(EXIT_FAILURE);
}

int option = 0; //0:l , 1:c, 2:p

int main(int argc, char** argv)
{
  //gestion des paramètres
  if (argc != 6) usage(argv[0]);
  
  //Parser les options
	memeObj	meme = memeObj::load_meme(argv[1]);  
  Dictionnaire dico = Dictionnaire::load_fasta(argv[2]);
  
  ifstream hit_file(argv[3], std::ifstream::in);
  if (!hit_file.is_open())
  {
    cerr << "erreur à l'ouverture de " << argv[3] << " en lecture" << endl;
    exit(EXIT_FAILURE);
  }
  
  int l = atoi(argv[4]);
  
  switch(argv[5][0]) {
  
  	case 'l':
  		option=0;
  		break;
  	case 'c':
  		option=1;
  		break;
  	case 'p':
  		option=2;
  		break;
  	case 's':
  		option=3;
  		break;
  	case 'i':
  		option=4;
  		break;
  	case 'd':
  		option=5;
  		break;
  	case 'v':
  		option=6;
  		break;
  	case 'k':
  		option=7;
  		break;
  	default:
  		usage(argv[0]);
  		break;
  }
  
  if(option == 7) {
  	dico.cat(cout);
  	exit(EXIT_SUCCESS);
  }
  
  string junk;

  //remplissage de la table
  Dom table(l, dico, meme);
  string sequence_name;
  string motif_ref;
  string strand;
  int start;
  int stop;
  while(hit_file >> motif_ref)
  {
    if(motif_ref.find("#")!=string::npos)
    {
    	getline(hit_file, junk);
    	continue;
    }
    
    hit_file >> sequence_name >> start >> stop >> strand;
    getline(hit_file, junk);
    
    Cluster c(
    	sequence_name, 
    	motif_ref, 
    	meme.ref_to_altname(motif_ref), 
    	meme.ref_to_index(motif_ref), 
    	strand, 
    	start, 
    	stop
    );
    
    table.add_element(c);
  }
  
  table.construire_voisinages();
  
  if(option == 0) cout << table;
  if(option == 1) table.output_couples(table.compter_couples(),cout);
	if(option == 2) table.output_pseudo_pal(table.compter_pseudo_pal(),cout);
	if(option == 3) table.output_presence_couples_seq(table.presence_couples_par_seq(),cout);
  if(option == 4) table.output_motifs_par_seq(table.motifs_par_seq(),cout);
  if(option == 5) table.to_svg();
  if(option == 6) table.afficher_annuaire2(cout);
  if(option == 7) dico.cat(cout);
  
  //Fermer les fichiers avant de quitter
  hit_file.close();
  
  return 0;
}
