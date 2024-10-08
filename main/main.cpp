// main.cpp: This file contains the main function.
#pragma warning (disable : 4819)
#include<iostream>
#include<string>
#include<Eigen/Dense>
#include<vector>
#include<FeatureSet.h>
#include<ReadSSSCommand.h>
#include<ReadControl.h>
#include<ViewerControl.h>
#include<RunControl.h>
#include<InteractiveControl.h>
#include<DefMesh.h>
using namespace std;
using namespace Eigen;

string fileplace = "Input/";
string outputplace = "Output/";
string modelplace = "Input/Model/";
string envplace = "Input/Environment/";
string fileformat = ".off";
string outputformat = ".log";
string jsonformat = ".json";
string version = "1.6";

double Point_Size = 15.0;
double Line_Width = 1.0;

int ArgcSize = 17;
int ModeSize = 3;

EnvironmentFeature env;
VectorXd SSSalpha(7), SSSbeta(7);
double r0 = 1.0;
double envrange = 4.0;
double varepsilon = 0.05;
int ExpandLimit = 262144;
heutype SSSheu = WIDTH;
bool SSSshow = false;
bool SSSCheckYellow = false;
int ExpandShow = 200;
bool box_draw_strategy = false;
vector<Vector3d> SSShints;
string SSSfilename;
string EnvName;
Vector3d ViewPoint;

// Any possible tests.
void test()
{
	cout << "Hello World!" << endl;
}

int main(int argc,char* argv[])
{
	int mode = -1, arg = 0;
	string filename;
	if (argc > 1)
		mode = stoi(argv[1]);
	// SSS framework main program.
	if (abs(mode) < 2)
	{
		if (argc > 2)
			arg = stoi(argv[2]);
		if (arg)
		{
			if (!json_argv(argc, argv, filename))
				return 1;
			read_json(filename);
		}
		else
			read_argv(argc, argv);
		switch (mode)
		{
		case 1: interactive(); break;
		case 0: non_interactive(SSSshow); break;
		case -1: show_only(); break;
		default: test(); break;
		}
	}
	// Assistant functions.
	else
	{
		switch (mode)
		{
		case 2:set_mesh(argc, argv); break;
		case 3:merge_mesh(argc, argv); break;
		case 4:split_mesh(argc, argv); break;
		default:test(); break;
		}
	}
	return 0;
}