// main.cpp: This file contains the main function.
#pragma warning (disable : 4819)
#include<iostream>
#include<Eigen/Dense>
#include<FeatureSet.h>
#include<ReadSSSCommand.h>
#include<ReadControl.h>
#include<ViewerControl.h>
#include<RunControl.h>
#include<InteractiveControl.h>
using namespace std;
using namespace Eigen;

string fileplace = "Input/";
string fileformat = ".off";
string jsonformat = ".json";
string version = "1.2";

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
	if (argc > 2)
		arg = stoi(argv[2]);
	if (arg)
	{
		if (!json_argv(argc, argv, filename))
			return -1;
		read_json(filename);
	}
	else
		read_argv(argc, argv);
	switch (mode)
	{
	case 1: interactive(); break;
	case 0: non_interactive(); break;
	case -1: show_only(); break;
	default: test(); break;
	}
	return 1;
}