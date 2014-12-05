/* A Bison parser, made by GNU Bison 3.0.2.  */

/* Bison interface for Yacc-like parsers in C

   Copyright (C) 1984, 1989-1990, 2000-2013 Free Software Foundation, Inc.

   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.  */

/* As a special exception, you may create a larger work that contains
   part or all of the Bison parser skeleton and distribute that work
   under terms of your choice, so long as that work isn't itself a
   parser generator using the skeleton or a modified version thereof
   as a parser skeleton.  Alternatively, if you modify or redistribute
   the parser skeleton itself, you may (at your option) remove this
   special exception, which will cause the skeleton and the resulting
   Bison output files to be licensed under the GNU General Public
   License without this special exception.

   This special exception was added by the Free Software Foundation in
   version 2.2 of Bison.  */

#ifndef YY_MM_PARSER_MEME_HPP_INCLUDED
# define YY_MM_PARSER_MEME_HPP_INCLUDED
/* Debug traces.  */
#ifndef YYDEBUG
# define YYDEBUG 0
#endif
#if YYDEBUG
extern int mmdebug;
#endif
/* "%code requires" blocks.  */
#line 24 "parser_meme.y" /* yacc.c:1909  */

#include <vector>
#include <iostream>
#include <map>
#include <deque>
#include <fstream>
#include "memeObj.h"
#include "motif.h"
#include "pwm.h"

#line 55 "parser_meme.hpp" /* yacc.c:1909  */

/* Token type.  */
#ifndef YYTOKENTYPE
# define YYTOKENTYPE
  enum yytokentype
  {
    NOMBRE = 258,
    DOUBLE = 259,
    STRING = 260,
    MEME_VERSION = 261,
    ALPHABET = 262,
    STRANDS = 263,
    BG_LETTER_FRQ = 264,
    PWM = 265,
    ALEN = 266,
    MOTIF_W = 267,
    NSITES = 268,
    E = 269,
    MOTIF = 270,
    URL = 271,
    EOL = 272
  };
#endif

/* Value type.  */
#if ! defined YYSTYPE && ! defined YYSTYPE_IS_DECLARED
typedef union YYSTYPE YYSTYPE;
union YYSTYPE
{
#line 41 "parser_meme.y" /* yacc.c:1909  */

	char* c;
	double d;
	int i;
	map<string, double>* m;
	vector<string>*	s;
	deque<double>* dd;
	map<char, vector<double> >* mvd;
	vector<Motif>* vmtf;
	Pwm* p;

#line 97 "parser_meme.hpp" /* yacc.c:1909  */
};
# define YYSTYPE_IS_TRIVIAL 1
# define YYSTYPE_IS_DECLARED 1
#endif


extern YYSTYPE mmlval;

int mmparse (memeObj& meme);

#endif /* !YY_MM_PARSER_MEME_HPP_INCLUDED  */
