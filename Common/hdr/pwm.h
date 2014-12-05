#ifndef PWM_H
#define PWM_H
#include <iostream>
#include <list>
#include <vector>
#include <deque>
#include <map>
#include <tuple>
#include <string>
#include <cstdlib>
#include <fstream>
#include <cmath>

using namespace std;

class Pwm
{
	public:
	void    				afficher 					(ostream&) const; 
	void						set_w							(int);
	void						set_alength				(int);
	void						set_mat						(map<char, vector<double> >);
	double					entropy						();
	int							longueur					();
	vector<double>& operator[](char letter);
	
	private:
	int							w;
	int							alength;
	map<char, vector<double> >  matrice;
};

ostream& operator<<(ostream&, Pwm const&);

#endif
