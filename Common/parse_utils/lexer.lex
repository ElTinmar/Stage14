%{

#include <iostream>
#include <cstdlib>
#include "parser.hpp"
#define yyterminate() yylex_destroy(); return YY_NULL

using namespace std;

%}

/* nom du fichier c généré par flex */
%option outfile="lexer.cpp"
%option header-file="lexer.hpp"
%option nounput
%option noinput
%option noyywrap


string	  [ \tA-Za-z0-9_\[\]\(\)\-,/]+

%%

">"     { return(CHEV); }
"|"     { return(PIPE); }
":"			{ return(COL); }
"\n"		{ return(EOL); }
{string}  { yylval.c = strdup(yytext); return(STRING); }



