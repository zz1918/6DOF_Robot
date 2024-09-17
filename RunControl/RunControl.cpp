// RunControl.cpp: This file controls the running of SSS main program.

#ifndef RUNCONTROL_H
#define RUNCONTROL_H

#pragma warning (disable : 4819)
#include<iostream>
#include<Eigen/Dense>
#include<FeatureSet.h>
#include<ReadSSSCommand.h>
#include<SSS.h>
#include<IntroControl.h>
#include<ViewerControl.h>

EnvironmentFeature extern env;
VectorXd extern SSSalpha, SSSbeta;
double extern envrange;
double extern r0;
double extern varepsilon;
int extern ExpandLimit;
heutype extern SSSheu;
bool extern SSSshow;
SSSViewer extern viewer;

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
	vector<VectorXd> path = SE3SSS.Find_Path(SSSalpha, SSSbeta, *Omega, varepsilon, SSSheu, SSSshow);
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