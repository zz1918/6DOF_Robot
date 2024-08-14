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

enum config { ALPHA, BETA, RANGE, OTHER };

config toconfig(string word)
{
	if (word == "alpha")
		return ALPHA;
	if (word == "beta")
		return BETA;
	if (word == "range")
		return RANGE;
	return OTHER;
}

#endif