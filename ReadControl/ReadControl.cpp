// ReadControl.cpp: This file controls the input of SSS main program.

#ifndef READCONTROL_H
#define READCONTROL_H

#include<iostream>
#include<string>
#include<Eigen/Dense>
#include<FeatureSet.h>
#include<ReadSSSCommand.h>
#include<FileCheck.h>
using namespace std;

// JSON parser library (https://github.com/nlohmann/json)
#include "json.hpp"
using json = nlohmann::json;

int extern ArgcSize;
int extern ModeSize;

string extern fileplace;
string extern fileformat;
string extern jsonformat;

EnvironmentFeature extern env;
VectorXd extern SSSalpha, SSSbeta;
double extern r0;
double extern envrange;
double extern varepsilon;
int extern ExpandLimit;
heutype extern SSSheu;
bool extern SSSshow;
int extern ExpandShow;
bool extern box_draw_strategy;

// Read by json.
void read_json(string filename)
{
	// Load json data from scene file
	json data;
	ifstream in(filename);

	in >> data;

	// Helper function to read a Vector3d from a json array
	auto read_vec3 = [](const json& x) {
		return Vector3d(x[0], x[1], x[2]);
		};

	// Read required arguments.
	double alphaOx = data["startOx"], alphaOy = data["startOy"], alphaOz = data["startOz"];
	double alphaAphi = data["startAphi"], alphaAtheta = data["startAtheta"], alphaBtheta = data["startBtheta"];
	double betaOx = data["goalOx"], betaOy = data["goalOy"], betaOz = data["goalOz"];
	double betaAphi = data["goalAphi"], betaAtheta = data["goalAtheta"], betaBtheta = data["goalBtheta"];
	ExpandLimit = data["ExpandLimit"];
	r0 = data["edgeL"];
	envrange = data["envrange"];
	varepsilon = data["epsilon"];
	SSSheu = toheutype(data["Qtype"]);
	SSSshow = data["show"];
	ExpandShow = data["ExpandShow"];
	box_draw_strategy = data["DrawBox"];
	Vector3d alphaO(alphaOx, alphaOy, alphaOz);
	Vector4d alphaQ = SO3(alphaAphi, alphaAtheta, alphaBtheta).Q();
	Vector3d betaO(betaOx, betaOy, betaOz);
	Vector4d betaQ = SO3(betaAphi, betaAtheta, betaBtheta).Q();
	SSSalpha << alphaO(0), alphaO(1), alphaO(2), alphaQ(0), alphaQ(1), alphaQ(2), alphaQ(3);
	SSSbeta << betaO(0), betaO(1), betaO(2), betaQ(0), betaQ(1), betaQ(2), betaQ(3);

	// Read obstacles.
	for (const auto& entryc : data["obstacle"]["corner"]) {
		string name = entryc["name"];
		Vector3d coord = read_vec3(entryc["position"]);
		env.add_point(coord, name);
		cout << "Define point " << name << " as " << coord.transpose() << endl;
	}

	for (const auto& entrye : data["obstacle"]["edge"]) {
		string name = entrye["name"];
		string name1 = entrye["p1"];
		string name2 = entrye["p2"];
		if (env.add_edge(name1, name2, name) == NULL)
			cout << "Some point not found" << endl;
		else
			cout << "Define edge " << name << " as the edge connecting point " << name1 << " and point " << name2 << endl;
	}

	for (const auto& entryf : data["obstacle"]["face"]) {
		string name = entryf["name"];
		string name1 = entryf["p1"];
		string name2 = entryf["p2"];
		string name3 = entryf["p3"];
		if (env.add_face(name1, name2, name3, name) == NULL)
			cout << "Some point not found" << endl;
		else
			cout << "Define face " << name << " as the triangle connecting point " << name1 << " and point " << name2 << " and point " << name3 << endl;
	}

	for (const auto& entrym : data["obstacle"]["mesh"]) {
		if (string(entrym).empty())
			continue;
		string offfilename = fileplace + string(entrym) + fileformat;
		if (!exist(offfilename))
		{
			cout << "bash: file " << offfilename << " not found" << endl;
			continue;
		}
		env.add_mesh(offfilename, string(entrym));
	}
}

// Read by arguments.
void read_argv(int argc, char* argv[])
{
	double alphaOx = -envrange / 2.0 - 0.5, alphaOy = -envrange / 2.0 - 0.5, alphaOz = -envrange / 2.0 - 0.5;
	double alphaAphi = 0, alphaAtheta = 0, alphaBtheta = 0;
	double betaOx = envrange / 2.0 + 0.5, betaOy = envrange / 2.0 + 0.5, betaOz = envrange / 2.0 + 0.5;
	double betaAphi = 0, betaAtheta = 0, betaBtheta = 0;
	if (argc > ModeSize + ArgcSize)
	{
		for (int i = ModeSize + ArgcSize; i < argc; ++i)
		{
			if (string(argv[i]).empty())
				continue;
			string filename = fileplace + string(argv[i]) + fileformat;
			if (!exist(filename))
			{
				cout << "bash: file " << filename << " not found" << endl;
				continue;
			}
			env.add_mesh(filename, string(argv[i]));
		}
	}
	switch (argc - ModeSize)
	{
	default:
	case 17: ExpandLimit = stoi(argv[ModeSize + 16]);
	case 16: betaBtheta = stod(argv[ModeSize + 15]);
	case 15: betaAtheta = stod(argv[ModeSize + 14]);
	case 14: betaAphi = stod(argv[ModeSize + 13]);
	case 13: betaOz = stod(argv[ModeSize + 12]);
	case 12: betaOy = stod(argv[ModeSize + 11]);
	case 11: betaOx = stod(argv[ModeSize + 10]);
	case 10: alphaBtheta = stod(argv[ModeSize + 9]);
	case 9: alphaAtheta = stod(argv[ModeSize + 8]);
	case 8: alphaAphi = stod(argv[ModeSize + 7]);
	case 7: alphaOz = stod(argv[ModeSize + 6]);
	case 6: alphaOy = stod(argv[ModeSize + 5]);
	case 5: alphaOx = stod(argv[ModeSize + 4]);
	case 4: r0 = stod(argv[ModeSize + 3]);
	case 3: envrange = stod(argv[ModeSize + 2]);
	case 2: varepsilon = stod(argv[ModeSize + 1]);
	case 1: SSSheu = toheutype(argv[ModeSize]);
	case 0: break;
	}
	Vector3d alphaO(alphaOx, alphaOy, alphaOz);
	Vector4d alphaQ = SO3(alphaAphi, alphaAtheta, alphaBtheta).Q();
	Vector3d betaO(betaOx, betaOy, betaOz);
	Vector4d betaQ = SO3(betaAphi, betaAtheta, betaBtheta).Q();
	SSSalpha << alphaO(0), alphaO(1), alphaO(2), alphaQ(0), alphaQ(1), alphaQ(2), alphaQ(3);
	SSSbeta << betaO(0), betaO(1), betaO(2), betaQ(0), betaQ(1), betaQ(2), betaQ(3);
}

// Deal arguments for json.
bool json_argv(int argc, char* argv[], string& filename)
{
	if (argc < 4)
	{
		cout << "bash: not enough input arguments!" << endl;
		cout << "MISSING: json filename" << endl;
		return false;
	}
	filename = fileplace + string(argv[3]) + jsonformat;
	if (!exist(filename))
	{
		cout << "bash: file " << filename << " not found" << endl;
		return false;
	}
	return true;
}

#endif