DEBUG=yes
CXX=g++
OBJDIR = ../obj
VPATH = ../src:../hdr:../obj
ifeq ($(DEBUG),yes)
	CFLAGS=-W -Wall -ansi -pedantic -g -std=c++0x -I ../hdr
else
	CFLAGS=-O3 -g -std=c++0x -I ../hdr
endif

all: parser.o lexer.o parser_meme.o lexer_meme.o

meme: parser_meme.o lexer_meme.o

fasta: parser.o lexer.o

$(OBJDIR)/parser.o: parser.cpp parser.hpp lexer.hpp lexer.cpp 
	$(CXX) -c $(CFLAGS) $< -o $@
  
$(OBJDIR)/lexer.o: lexer.cpp lexer.hpp parser.cpp
	$(CXX) -c $(CFLAGS) $< -o $@

parser.cpp: parser.y lexer.cpp
	bison -d -r all $<

lexer.cpp: lexer.lex
	flex $<
	
$(OBJDIR)/parser_meme.o: parser_meme.cpp parser_meme.hpp lexer_meme.hpp lexer_meme.cpp 
	$(CXX) -c $(CFLAGS) $< -o $@
  
$(OBJDIR)/lexer_meme.o: lexer_meme.cpp lexer_meme.hpp parser_meme.cpp
	$(CXX) -c $(CFLAGS) $< -o $@

parser_meme.cpp: parser_meme.y lexer_meme.cpp
	bison -d -r all $<

lexer_meme.cpp: lexer_meme.lex
	flex $<
	
clean:
	rm -f $(OBJDIR)/parser*.o $(OBJDIR)/lexer*.o *.cpp *.hpp *.output
	
superclean:
	rm -f $(OBJDIR)/parser*.o $(OBJDIR)/lexer*.o *.cpp *.hpp *.output

