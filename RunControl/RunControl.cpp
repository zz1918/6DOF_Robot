// RunControl.cpp: This file controls the running of SSS main program.

#ifndef RUNCONTROL_H
#define RUNCONTROL_H

#pragma warning (disable : 4819)
#include<iostream>
#include<fstream>
#include<string>
#include<vector>
#include<Eigen/Dense>
#include<FeatureSet.h>
#include<ReadSSSCommand.h>
#include<SSS.h>
#include<IntroControl.h>
#include<ViewerControl.h>
#include <chrono>
using namespace std;
using namespace Eigen;

EnvironmentFeature extern env;
VectorXd extern SSSalpha, SSSbeta;
double extern envrange;
double extern r0;
double extern varepsilon;
int extern ExpandLimit;
heutype extern SSSheu;
bool extern SSSshow;
SSSViewer extern viewer;
string extern SSSfilename;
string extern outputplace;
string extern outputformat;

void write_in_file(string filename, vector<VectorXd>& path, double time_in_second)
{
	ofstream os(filename);

	if (path.empty())
		os << "No-Path!" << endl;
	else
	{
		os << "The path is the line segments connecting the following points: " << endl;
		os << SSSalpha.transpose() << endl;
		for (int i = 0; i < path.size(); ++i)
			os << path[i].transpose() << endl;
		os << SSSbeta.transpose() << endl;
	}

	os << endl;
	os << "Find path total time: " << time_in_second << "s." << endl;
}

// Run SSS main program.
void SSS_run(bool show = false)
{
	MatrixId EnvRange(Vector3d(-envrange, -envrange, -envrange), Vector3d(envrange, envrange, envrange));
	SE3Tree* T = new SE3Tree(EnvRange);
	DeltaPredicate* C = new DeltaPredicate();
	DeltaFeature* Omega = env.make_feature();
	SSS<VectorXd, SE3Box, SE3Tree, DeltaPredicate, DeltaFeature> SE3SSS(T, C);
	viewer.set_env(env, envrange, SSSalpha, SSSbeta);
	cout << "Simple find path algorithm with epsilon = " << varepsilon << endl;
	auto start_time = std::chrono::high_resolution_clock::now();
	vector<VectorXd> path = SE3SSS.Find_Path(SSSalpha, SSSbeta, *Omega, varepsilon, SSSheu, SSSshow);
	auto end_time = std::chrono::high_resolution_clock::now();
	double time_in_second = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time).count() / 1000000.0;
	cout << "Find path total time: " << time_in_second << "s." << endl;
	if (SSSshow)
		viewer.view();
	if (!path.empty())
		output_path(SSSalpha, SSSbeta, path);
	else
		cout << "Algorithm terminates. Probably no-path, or the maximum expansion limit is too small." << endl;
	viewer.clear();
	viewer.set_env(env, envrange, SSSalpha, SSSbeta);
	viewer.set_path(path);
	viewer.view();
	write_in_file(outputplace + SSSfilename + outputformat, path, time_in_second);
	delete T;
	delete C;
	delete Omega;
}

// Non-interactive mode.
void non_interactive(bool show=false)
{
	SSS_intro();
	SSS_run(show);
}

#endif