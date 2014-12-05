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
#include "hitList.h"

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
	MemeObj	meme = MemeObj::load_meme(argv[1]);  
  Dictionnaire dico = Dictionnaire::load_fasta(argv[2]);
  int l = atoi(argv[4]);
  HitList hitList = HitList::load_hitFile(argv[3],l,meme,dico);
  
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
  
  if(option == 0) cout << hitList;
  if(option == 1) hitList.output_couples(hitList.compter_couples(),cout);
	if(option == 2) hitList.output_pseudo_pal(hitList.compter_pseudo_pal(),cout);
	if(option == 3) hitList.output_presence_couples_seq(hitList.presence_couples_par_seq(),cout);
  if(option == 4) hitList.output_motifs_par_seq(hitList.motifs_par_seq(),cout);
  if(option == 5) hitList.to_svg();
  if(option == 6) hitList.afficher_annuaire2(cout);
  if(option == 7) dico.cat(cout);
  
  return 0;
}
