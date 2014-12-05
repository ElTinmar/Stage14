%{

#include <iostream>
#include "dictionnaire.h"
#include <string>
#include <cstring>
#include <cstdio>

using namespace std;
int yylex(void);
int yyerror(Dictionnaire& dico, const char* s);

%}

/*definit le nom des fichiers c et h crées par bison*/
%output  "parser.cpp"
%defines "parser.hpp"

/* permet d'inclure du code dans le fichier h généré par bison. */
%code requires {
#include "dictionnaire.h"
}

/*par défaut yyparse ne prends pas d'argument. %parse-param permet
de changer ce comportement et d'éviter de déclarer des variables globales
dans le cas où le main est séparé*/
%parse-param {Dictionnaire& dico}
%error-verbose
%union {
	char* c;
	Element* e;
}

%token CHEV PIPE COL EOL
%token <c> STRING 

%type <e> Cpremiere_ligne
%type <c> sequence

%start fichier_comp
%%

fichier_comp:
  comp_sequence
  | fichier_comp comp_sequence
  ;
  
sequence:
	STRING EOL {$$ = $1;}
	| STRING EOL sequence 
	{
		char* seq = (char*) malloc(strlen($3)+121*sizeof(char));
		strcpy(seq, $1);
		strcat(seq, $3);
		free($1);free($3);
		$$ = seq;
	}
	;
	
comp_sequence:
  Cpremiere_ligne sequence EOL
  {
  	$1->add_seq($2);
    dico.add_element(*$1);
    delete($1); free($2);
  }
  | Cpremiere_ligne sequence
  {
  	$1->add_seq($2);
    dico.add_element(*$1);
    delete($1); free($2);
  }
  ;
  
junk:
	| PIPE junk
	| COL junk
	| STRING junk {free($1);}
	;

Cpremiere_ligne:
	CHEV STRING junk EOL
	{
		Element* e = new Element ($2);
		free($2);
		$$ = e;
	}
	;

%%

int yyerror(Dictionnaire& dico, const char *s) 
{
  cout << s << endl;
  (void)dico;
  return 0;
}

