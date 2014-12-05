DEBUG=yes
CXX=g++
COMMON_DIR = ../Common

VPATH = $(COMMON_DIR)/src:$(COMMON_DIR)/hdr:$(COMMON_DIR)/obj:$(COMMON_DIR)/parse_utils
OBJDIR = $(COMMON_DIR)/obj
PARSETOOLS_DIR = $(COMMON_DIR)/parse_utils
SRCDIR = $(COMMON_DIR)/src
HDRDIR = $(COMMON_DIR)/hdr
INCLUDES = -I $(COMMON_DIR)/hdr -I $(COMMON_DIR)/parse_utils

ifeq ($(DEBUG),yes)
	CFLAGS=-W -Wall -ansi -pedantic -g -std=c++0x $(INCLUDES)
	LDFLAGS=
else
	CFLAGS=-O3 -g -std=c++0x $(INCLUDES)
	LDFLAGS=
endif

all: parse_tools main

main: main.o dictionnaire.o element.o lexer.o parser.o lexer_meme.o parser_meme.o cluster.o hitList.o memeObj.o motif.o pwm.o affichage.o
	$(CXX) -o clust $^

main.o: main.cpp parser.cpp lexer.cpp
	$(CXX) -c $(CFLAGS) $< -o $@
	
$(OBJDIR)/hitList.o: hitList.cpp cluster.cpp hitList.h cluster.h
	$(CXX) -c $(CFLAGS) $< -o $@
	
$(OBJDIR)/cluster.o: cluster.cpp cluster.h
	$(CXX) -c $(CFLAGS) $< -o $@

$(OBJDIR)/dictionnaire.o: dictionnaire.cpp dictionnaire.h
	$(CXX) -c $(CFLAGS) $< -o $@
	
$(OBJDIR)/element.o: element.cpp element.h
	$(CXX) -c $(CFLAGS) $< -o $@

$(OBJDIR)/affichage.o: affichage.cpp affichage.h
	$(CXX) -c $(CFLAGS) $< -o $@
	
$(OBJDIR)/memeObj.o: memeObj.cpp memeObj.h
	$(CXX) -c $(CFLAGS) $< -o $@
	
$(OBJDIR)/motif.o: motif.cpp motif.h
	$(CXX) -c $(CFLAGS) $< -o $@
	
$(OBJDIR)/pwm.o: pwm.cpp pwm.h
	$(CXX) -c $(CFLAGS) $< -o $@
	
parse_tools: 
	$(MAKE) -C $(PARSETOOLS_DIR) all
	
$(PARSETOOLS_DIR)/parser.cpp:

$(PARSETOOLS_DIR)/lexer.cpp:

$(PARSETOOLS_DIR)/parser.hpp:

$(PARSETOOLS_DIR)/lexer.hpp:

$(OBJDIR)/lexer.o:

$(OBJDIR)/parser.o:

$(PARSETOOLS_DIR)/parser_meme.cpp:

$(PARSETOOLS_DIR)/lexer_meme.cpp:

$(PARSETOOLS_DIR)/parser_meme.hpp:

$(PARSETOOLS_DIR)/lexer_meme.hpp:

$(OBJDIR)/lexer_meme.o:

$(OBJDIR)/parser_meme.o:
	
clean:
	$(MAKE) -C $(PARSETOOLS_DIR) clean
	rm -f $(OBJDIR)/element.o $(OBJDIR)/dictionnaire.o $(OBJDIR)/memeObj.o $(OBJDIR)/motif.o $(OBJDIR)/pwm.o $(OBJDIR)/hitList.o $(OBJDIR)/cluster.o $(OBJDIR)/affichage.o *.o

superclean:
	$(MAKE) -C $(PARSETOOLS_DIR) clean
	rm -f clust $(OBJDIR)/element.o $(OBJDIR)/dictionnaire.o $(OBJDIR)/memeObj.o $(OBJDIR)/motif.o $(OBJDIR)/pwm.o $(OBJDIR)/hitList.o $(OBJDIR)/cluster.o $(OBJDIR)/affichage.o *.o

