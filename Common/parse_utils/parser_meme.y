%{

#include <vector>
#include <deque>
#include <iostream>
#include <map>
#include <fstream>
#include "memeObj.h"
#include "motif.h"
#include "pwm.h"

using namespace std;

int mmlex(void);
int mmerror(MemeObj& meme, const char* s);

%}

/*definit le nom des fichiers c et h crées par bison*/
%output  "parser_meme.cpp"
%defines "parser_meme.hpp"

/* permet d'inclure du code dans le fichier h généré par bison. */
%code requires {
#include <vector>
#include <iostream>
#include <map>
#include <deque>
#include <fstream>
#include "memeObj.h"
#include "motif.h"
#include "pwm.h"
}

/*par défaut yyparse ne prends pas d'argument. %parse-param permet
de changer ce comportement et d'éviter de déclarer des variables globales
dans le cas où le main est séparé*/
%parse-param {MemeObj& meme} 
%error-verbose
%name-prefix "mm"
%union {
	char* c;
	double d;
	int i;
	map<string, double>* m;
	vector<string>*	s;
	deque<double>* dd;
	map<char, vector<double> >* mvd;
	vector<Motif>* vmtf;
	Pwm* p;
}

%type <m> frq_list
%type <s> name_line
%type <p> pwm
%type <mvd> mat
%type <dd> position
%type <c> url
%type <vmtf> motifs

%token <i> NOMBRE
%token <d> DOUBLE
%token <c> STRING
%token MEME_VERSION ALPHABET STRANDS BG_LETTER_FRQ PWM ALEN MOTIF_W NSITES E MOTIF URL EOL

%start meme
%%

eol:
	EOL
	| EOL eol
	;
	
junk:
	STRING eol {free($1);}
	| STRING junk {free($1);}
	;
	
meme:
  version alphabet strands bg_frq motifs {
  	meme.set_motifs(*$5);
  	delete($5);
  }
  | version alphabet strands motifs {
			meme.set_motifs(*$4);
			delete($4);
		}
  | version alphabet bg_frq motifs {
			meme.set_motifs(*$4);
			delete($4);
		}
  | version alphabet motifs {
			meme.set_motifs(*$3);
			delete($3);
		}
  | version strands bg_frq motifs {
			meme.set_motifs(*$4);
			delete($4);
		}
  | version strands motifs {
			meme.set_motifs(*$3);
			delete($3);
		}
  | version bg_frq motifs {
			meme.set_motifs(*$3);
			delete($3);
		}
  | version motifs {
			meme.set_motifs(*$2);
			delete($2);
		}
  ;

version:
	MEME_VERSION STRING eol {
	
		meme.set_version($2); 
		free($2);
	}
	| MEME_VERSION DOUBLE eol {
	
		meme.set_version(to_string($2));
	}
	;

alphabet:
	ALPHABET STRING eol {
	
		meme.set_alphabet($2);
		free($2);
	}
	;

strands:
	STRANDS STRING eol {
		meme.add_strand($2);
		free($2);
	}
	| STRANDS STRING STRING eol {
			meme.add_strand($2); 
			meme.add_strand($3);
			free($2); free($3);
		}
	;
	
bg_frq:
	BG_LETTER_FRQ junk frq_list eol {
	
		map<string, double>* m = $3;
		meme.set_bg_frq(*m);
		delete($3);
	}
	;

frq_list:
	STRING DOUBLE {
	
		map<string, double>* m = new map<string, double>;
		(*m)[string($1)]=$2;
		free($1);
		$$=m;
	}
	| STRING DOUBLE frq_list {
	
			map<string, double>* m = $3;
			(*m)[string($1)]=$2;
			free($1);
			$$=m;
		}
	;
	
motifs:
	name_line pwm url {
		vector<string>* names = $1;
		vector<Motif>* v = new vector<Motif>;
		if(names->size() == 1) {
			Motif m((*names)[0], *$2, $3);
			v->push_back(m);
			delete($1);
			delete($2);
			free($3);
			$$ = v;
		}
		if(names->size() == 2) {
			Motif m((*names)[0], (*names)[1], *$2, $3);
			v->push_back(m);
			delete($1);
			delete($2);
			free($3);
			$$ = v;
		}
	}
	| name_line pwm {
			vector<string>* names = $1;
			vector<Motif>* v = new vector<Motif>;
			if(names->size() == 1) {
				Motif m((*names)[0], *$2);
				v->push_back(m);
				delete($1);
				delete($2);
				$$ = v;
			}
			if(names->size() == 2) {
				Motif m((*names)[0], (*names)[1], *$2);
				v->push_back(m);
				delete($1);
				delete($2);
				$$ = v;
			}
		}
	| name_line pwm url motifs {
			vector<string>* names = $1;
			vector<Motif>* v = $4;
			if(names->size() == 1) {
				Motif m((*names)[0], *$2, $3);
				v->push_back(m);
				delete($1);
				delete($2);
				free($3);
				$$ = v;
			}
			if(names->size() == 2) {
				Motif m((*names)[0], (*names)[1], *$2, $3);
				v->push_back(m);
				delete($1);
				delete($2);
				free($3);
				$$ = v;
			}
		}
	| name_line pwm motifs {
			vector<string>* names = $1;
			vector<Motif>* v = $3;
			if(names->size() == 1) {
				Motif m((*names)[0], *$2);
				v->push_back(m);
				delete($1);
				delete($2);
				$$ = v;
			}
			if(names->size() == 2) {
				Motif m((*names)[0], (*names)[1], *$2);
				v->push_back(m);
				delete($1);
				delete($2);
				$$ = v;
			}
		}
	;
	
name_line:
	MOTIF STRING eol {
	
		vector<string>* v = new vector<string>;
		v->push_back($2);
		free($2);
		$$=v;
	}
	| MOTIF STRING STRING eol {
	
			vector<string>* v = new vector<string>;
			v->push_back($2);
			v->push_back($3);
			free($2);
			free($3);
			$$=v;
		}
	;
	
pwm:
  PWM ALEN NOMBRE MOTIF_W NOMBRE NSITES NOMBRE E NOMBRE eol mat {
	
		Pwm* pwm = new Pwm;
		pwm->set_alength($3);
		pwm->set_w($5);
		pwm->set_mat(*$11); 
		$$=pwm;
	}
	;
	
mat:
	position eol {
	
		deque<double>* pos = $1;
		map<char, vector<double> >* mat = new map<char, vector<double> >;
		const set<char> alphabet = meme.get_alphabet();
		for(auto it=alphabet.begin(); it!=alphabet.end(); it++) {
			(*mat)[*it].push_back(pos->front());
			pos->pop_front();
		}
		delete($1);
		$$=mat;
	}
	| mat position eol {
	
			deque<double>* pos = $2;
			map<char, vector<double> >* mat = $1;
			const set<char> alphabet = meme.get_alphabet();
			for(auto it=alphabet.begin(); it!=alphabet.end(); it++) {
				(*mat)[*it].push_back(pos->front());
				pos->pop_front();
			}
			delete($2);
			$$=mat;
		}
	;

position:
	DOUBLE {
	
		deque<double>* v = new deque<double>;
		v->push_front($1);
		$$=v;
	}
	| DOUBLE position {
	
			deque<double>* v = $2;
			v->push_front($1);
			$$=v;
		}
	;
	
url:
	URL STRING eol {
	
		$$=$2;
	}
	;

%%

int mmerror(MemeObj& meme, const char *s) 
{
  cerr << "Erreur rencontree : " << s << endl;
  cerr << "Dernier token : " << mmlval.c << endl;
 	cerr << "Objet meme : " << endl << meme;
  return 0;
}

