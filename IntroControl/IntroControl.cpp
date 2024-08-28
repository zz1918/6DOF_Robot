// IntroControl.cpp: This file controls the introduction of SSS main program.

#ifndef INTROCONTROL_H
#define INTROCONTROL_H

#include<iostream>
#include<Eigen/Dense>
#include<FeatureSet.h>
#include<ReadSSSCommand.h>
using namespace std;

string extern version;
EnvironmentFeature extern env;
VectorXd extern SSSalpha, SSSbeta;
double extern r0;
double extern envrange;
double extern varepsilon;
int extern ExpandLimit;
heutype extern SSSheu;

// Introduce default settings.
void SSS_intro()
{
	cout << "Delta robot find path algorithm by SSS framework." << endl;
	cout << "Demo -- Version " << version << endl;
	cout << "Default environment range: " << "[-" << envrange << "," << envrange << "] * [-" << envrange << "," << envrange << "] * [-" << envrange << "," << envrange << "]" << endl;
	cout << "Default alpha: (" << SSSalpha(0) << "," << SSSalpha(1) << "," << SSSalpha(2) << "," << SSSalpha(3) << "," << SSSalpha(4) << "," << SSSalpha(5) << "," << SSSalpha(6) << ")" << endl;
	cout << "Default beta: (" << SSSbeta(0) << "," << SSSbeta(1) << "," << SSSbeta(2) << "," << SSSbeta(3) << "," << SSSbeta(4) << "," << SSSbeta(5) << "," << SSSbeta(6) << ")" << endl;
	cout << "Default epsilon: " << varepsilon << endl;
	cout << "Default heuristic type: " << SSSheu << "." << endl;
	cout << "Default maximum expansion: " << ExpandLimit << " cubes." << endl;
	env.show_mesh();
}

#endif