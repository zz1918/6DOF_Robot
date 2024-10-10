// RunControl.cpp: This file controls the running of SSS main program.

#ifndef RUNCONTROL_H
#define RUNCONTROL_H

#pragma warning (disable : 4819)
#include<iostream>
#include<fstream>
#include<string>
#include<iosfwd>
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
string extern EnvName;
int extern expanded;
int extern stat_mixed;
int extern stat_free;
int extern stat_stuck;
int extern stat_small;
double extern MaxSize;
double extern MinSize;
int extern Noisity;

string VecXd(VectorXd v)
{
	stringstream vi;
	vi << fixed << setprecision(2);
	vi << "(" << v(0);
	for (int i = 1; i < v.size(); ++i)
		vi << ", " << v(i);
	vi << ")";
	return vi.str();
}

void output_summarize(ostream &os, vector<VectorXd>& path, double time_in_second)
{
	os << endl;
	os << "Find path total time: " << time_in_second << "s." << endl;

	if (path.empty())
		os << "No-Path!" << endl;
	else
		os << "Path!" << endl;

	os << expanded << " cubes have been expanded." << endl;
	os << "There are " << stat_mixed << " mixed boxes, " << stat_free << " free boxes and " << stat_stuck << " stuck boxes and " << stat_small << " veps-small boxes." << endl;
	os << "Maximum leaf size in the end: " << MaxSize << endl;
	os << "Minimum leaf size in the end: " << MinSize << endl;

	/*
	os << "#################### Latex Statistic Table ####################" << endl;

	os << "\\begin{tabular}{ | c | c | c | c | c | c | c | c | c | c | c | c | }" << endl;
	os << "\\hline" << endl;
	os << " Env & Start Conf. & Goal Conf. & $\\veps$ & $l_0$ & Qtype & time(s) & Path & Free / Stuck / Mixed / $\\veps$-small & \\#Boxes & Max size & Min size \\\\" << endl;
	os << "\\hline" << endl;
	os << " " << EnvName << " & " << VecXd(SSSalpha) << " & " << VecXd(SSSbeta) << " & " << varepsilon << " & " << r0 << " & " << SSSheu << " & " << time_in_second << " & " << (path.empty() ? "N" : "Y");
	os << " & " << stat_free << "/" << stat_stuck << "/" << stat_mixed << "/" << stat_small << " & " << expanded << " & " << MaxSize << " & " << MinSize << " \\\\" << endl;
	os << "\\hline" << endl;
	os << "\\end{tabular}" << endl;

	os << "###############################################################" << endl;*/

	if (!path.empty() && Noisity >= 1)
	{
		os << "The path is the line segments connecting the following points: " << endl;
		os << SSSalpha.transpose() << endl;
		for (int i = 0; i < path.size(); ++i)
			os << path[i].transpose() << endl;
		os << SSSbeta.transpose() << endl;
	}
}

void write_in_file(string filename, vector<VectorXd>& path, double time_in_second)
{
	ofstream os(filename);

	os << endl;
	os << "Find path total time: " << time_in_second << "s." << endl;

	if (path.empty())
		os << "No-Path!" << endl;
	else
		os << "Path!" << endl;

	os << expanded << " cubes have been expanded." << endl;
	os << "There are " << stat_mixed << " mixed boxes, " << stat_free << " free boxes and " << stat_stuck << " stuck boxes and " << stat_small << " veps-small boxes." << endl;

	os << "#################### Latex Statistic Table ####################" << endl;

	os << "\\begin{tabular}{ | c | c | c | c | c | c | c | c | c | c | c | c | }" << endl;
	os << "\\hline" << endl;
	os << " Env & Start Conf. & Goal Conf. & $\\veps$ & $l_0$ & Qtype & time(s) & Path & Free / Stuck / Mixed / $\\veps$-small & \\#Boxes & Max size & Min size \\\\" << endl;
	os << "\\hline" << endl;
	os << " " << EnvName << " & " << VecXd(SSSalpha) << " & " << VecXd(SSSbeta) << " & " << varepsilon << " & " << r0 << " & " << SSSheu << " & " << time_in_second << " & " << (path.empty() ? "N" : "Y");
	os << " & " << stat_free << "/" << stat_stuck << "/" << stat_mixed << "/" << stat_small << " & " << expanded << " & " << MaxSize << " & " << MinSize << " \\\\" << endl;
	os << "\\hline" << endl;
	os << "\\end{tabular}" << endl;

	os << "###############################################################" << endl;

	if (!path.empty())
	{
		os << "The path is the line segments connecting the following points: " << endl;
		os << SSSalpha.transpose() << endl;
		for (int i = 0; i < path.size(); ++i)
			os << path[i].transpose() << endl;
		os << SSSbeta.transpose() << endl;
	}
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
	double time_in_second = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time).count() / (long double)(1000000.0);
	if (Noisity >= 5)
		viewer.view();
	output_summarize(cout, path, time_in_second);
	write_in_file(outputplace + SSSfilename + outputformat, path, time_in_second);
	if (Noisity >= 0)
	{
		viewer.clear();
		viewer.set_env(env, envrange, SSSalpha, SSSbeta);
		viewer.set_path(path);
		viewer.view();
	}
	delete T;
	delete C;
	delete Omega;
}

// Non-interactive mode.
void non_interactive(bool show = false)
{
	SSS_intro();
	SSS_run(show);
}

#endif