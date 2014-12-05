#ifndef MOTIF_H
#define MOTIF_H
#include <iostream>
#include <list>
#include <vector>
#include <deque>
#include <map>
#include <tuple>
#include <string>
#include <cstdlib>
#include <fstream>
#include "pwm.h"

using namespace std;

class Motif
{
	friend class MemeObj;
	
	public:
									Motif							(string,Pwm);
									Motif							(string,Pwm,string);
									Motif							(string,string,Pwm);
									Motif							(string,string,Pwm,string);
	void    				afficher 					(ostream&) const;
	double					entropy						();
	int							longueur					();
	
	private:
	string					identifier;
	string					alt_name;
	Pwm							matrice;
	string					url;
};

ostream &operator<<(ostream&, Motif const&);

#endif
