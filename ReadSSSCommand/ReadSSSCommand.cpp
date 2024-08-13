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
	return word[0] == '(';
}

// Turn string "(x,y,z)" into Vector3d x y z.
Vector3d stov3(string coord)
{
	vector<string> coords = split(coord.substr(1, coord.length() - 2), ',');
	return Vector3d(stod(coords[0]), stod(coords[1]), stod(coords[2]));
}

enum fpred { DEF, DEL, SHOW, RUN, EXIT, ERROR };

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
	if (word == "exit")
		return EXIT;
	return ERROR;
}

enum ftype { POINT, EDGE, FACE, MESH, WRONG };

ftype toftype(string word)
{
	if (word == "point")
		return POINT;
	if (word == "edge")
		return EDGE;
	if (word == "face")
		return FACE;
	if (word == "mesh")
		return MESH;
	return WRONG;
}

#endif