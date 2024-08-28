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
double extern r0;
double extern envrange;
double extern varepsilon;
int extern ExpandLimit;
heutype extern SSSheu;
SSSViewer extern viewer;

// Run SSS main program.
void SSS_run()
{
	MatrixId EnvRange(Vector3d(-envrange, -envrange, -envrange), Vector3d(envrange, envrange, envrange));
	SE3Tree* T = new SE3Tree(EnvRange);
	DeltaPredicate* C = new DeltaPredicate();
	DeltaFeature* Omega = env.make_feature();
	SSS<VectorXd, SE3Box, SE3Tree, DeltaPredicate, DeltaFeature> SE3SSS(T, C);
	cout << "Simple find path algorithm with epsilon = " << varepsilon << endl;
	vector<VectorXd> path = SE3SSS.Find_Path(SSSalpha, SSSbeta, *Omega, varepsilon, SSSheu);
	if (!path.empty())
		output_path(SSSalpha, SSSbeta, path);
	else
		cout << "Algorithm terminates. Probably no-path, or the maximum expansion limit is too small." << endl;
	viewer.view_path(env, envrange, SSSalpha, SSSbeta, path);
	delete T;
	delete C;
	delete Omega;
}

// Non-interactive mode.
void non_interactive()
{
	SSS_intro();
	SSS_run();
}

#endif