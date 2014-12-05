%{

#include <iostream>
#include "parser_meme.hpp"
#define mmterminate() mmlex_destroy(); return mm_NULL

using namespace std;

%}

/* nom du fichier c généré par flex */
%option outfile="lexer_meme.cpp"
%option header-file="lexer_meme.hpp"
%option nounput
%option noinput
%option noyywrap prefix="mm"

blancs    [ \t]+
entier    [0-9]+
dseq_opt  [0-9]*
frac      (({dseq_opt}"."{entier})|{entier}".")
string	  [A-Za-z0-9_=\[\]\(\)\-"+",/".""?""&"":"]+

%%

{blancs}  { /* On ignore */}

"\n"	{ return(EOL); }
"MEME version"	{ return(MEME_VERSION); }
"ALPHABET="	{ return(ALPHABET); }
"strands:"	{ return(STRANDS); }
"Background letter frequencies"	{ return(BG_LETTER_FRQ); }
"letter-probability matrix:"	{ return(PWM); }
"alength="	{ return(ALEN); }
"w="	{ return(MOTIF_W); }
"nsites="	{ return(NSITES); }
"E="	{ return(E); }
"MOTIF"	{ return(MOTIF); }
"URL"	{ return(URL); }
{entier}  { mmlval.i=atoi(mmtext); return(NOMBRE); }
{frac} { mmlval.d=atof(mmtext); return(DOUBLE); }
{string}  { mmlval.c=strdup(mmtext); return(STRING); }


