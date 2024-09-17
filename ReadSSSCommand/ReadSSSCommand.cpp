// ReadSSSCommand.cpp : Parse helpers of SSS command line.
//

#ifndef READSSSCOMMAND_H
#define READSSSCOMMAND_H

#include<iostream>
#include<sstream>
#include<string>
#include<vector>
#include<Eigen/Dense>
using namespace std;
using namespace Eigen;

vector<string> split(const string& s, char delim)
{
	vector<string> result;
	stringstream ss(s);
	string item;
	while(getline(ss, item, delim))
		result.push_back(item);
	return result;
}

bool isv3(string word)
{
	if (word[0] != '(')
		return false;
	if (word[word.length() - 1] != ')')
		return false;
	if (split(word.substr(1, word.length() - 2), ',').size() != 3)
		return false;
	return true;
}

// Turn string "(x,y,z)" into Vector3d x y z.
Vector3d stov3(string coord)
{
	vector<string> coords = split(coord.substr(1, coord.length() - 2), ',');
	return Vector3d(stod(coords[0]), stod(coords[1]), stod(coords[2]));
}

bool isv7(string word)
{
	if (word[0] != '(')
		return false;
	if (word[word.length() - 1] != ')')
		return false;
	if (split(word.substr(1, word.length() - 2), ',').size() != 7)
		return false;
	return true;
}

// Turn string "(x,y,z)" into Vector3d x y z.
VectorXd stov7(string coord)
{
	vector<string> coords = split(coord.substr(1, coord.length() - 2), ',');
	VectorXd vec(7);
	vec << stod(coords[0]), stod(coords[1]), stod(coords[2]), stod(coords[3]), stod(coords[4]), stod(coords[5]), stod(coords[6]);
	return vec;
}

enum fpred { DEF, DEL, SHOW, SET, RUN, EXIT, ERROR };

fpred tofpred(string word)
{
	if (word == "def")
		return DEF;
	if (word == "del")
		return DEL;
	if (word == "run")
		return RUN;
	if (word == "show")
		return SHOW;
	if (word == "set")
		return SET;
	if (word == "exit")
		return EXIT;
	return ERROR;
}

enum ftype { TPOINT, TEDGE, TFACE, TMESH, TWRONG };

ftype toftype(string word)
{
	if (word == "point")
		return TPOINT;
	if (word == "edge")
		return TEDGE;
	if (word == "face")
		return TFACE;
	if (word == "mesh")
		return TMESH;
	return TWRONG;
}

enum config { ALPHA, BETA, RANGE, HEURISTIC, LIMIT, EPSILON, OTHER };

config toconfig(string word)
{
	if (word == "alpha" || word == "aph")
		return ALPHA;
	if (word == "beta" || word == "bet")
		return BETA;
	if (word == "range" || word == "rag")
		return RANGE;
	if (word == "heuristic" || word == "heu")
		return HEURISTIC;
	if (word == "limit" || word == "lim")
		return LIMIT;
	if (word == "epsilon" || word == "varepsilon" || word == "eps" || word == "var")
		return EPSILON;
	return OTHER;
}

//************* Heuristic type functions *************//

enum heutype { RAND, BYID, WIDTH, TARGET, GBF, DIS, GBFDIS, WIDIS, RECUR };

heutype toheutype(string word)
{
	if (word == "id")
		return BYID;
	if (word == "width")
		return WIDTH;
	if (word == "target")
		return TARGET;
	if (word == "gbf")
		return GBF;
	if (word == "dis")
		return DIS;
	if (word == "gbfdis")
		return GBFDIS;
	if (word == "widis" || word == "width_dis" || word == "widthdis")
		return WIDIS;
	if (word == "recursive" || word == "recur" || word == "rec")
		return RECUR;
	return RAND;
}

ostream& operator<<(ostream& os, heutype heu)
{
	switch (heu)
	{
	case BYID: os << "id"; break;
	case WIDTH: os << "width"; break;
	case TARGET: os << "target"; break;
	case GBF: os << "greedy"; break;
	case DIS: os << "distance"; break;
	case GBFDIS: os << "greedy-distance"; break;
	case WIDIS:os << "width-distance"; break;
	case RECUR:os << "recursive"; break;
	default:os << "random";
	}
	return os;
}

//*****************************************************// 

#endif