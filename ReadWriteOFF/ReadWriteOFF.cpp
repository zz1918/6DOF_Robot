// ReadWriteOFF.cpp: this file defines the read and write from an .off file.

#include <iostream>
#include <string>
#include <fstream>
#include <Eigen/Dense>

using namespace std;
using namespace Eigen;

// Read a triangle mesh from an off file
void read_OFF(const string& filename, MatrixXd& V, MatrixXi& F) {
	ifstream in(filename);
	string token;
	in >> token;
	int nv, nf, ne;
	in >> nv >> nf >> ne;
	MatrixXd _V(nv, 3);
	MatrixXi _F(nf, 3);
	for (int i = 0; i < nv; ++i) {
		in >> _V(i, 0) >> _V(i, 1) >> _V(i, 2);
	}
	for (int i = 0; i < nf; ++i) {
		int s;
		in >> s >> _F(i, 0) >> _F(i, 1) >> _F(i, 2);
		assert(s == 3);
	}
	V = _V;
	F = _F;
}

// Write a triangle mesh to an off file
void write_OFF(const string& filename, MatrixXd& V, MatrixXi& F) {
	ofstream out(filename);
	string token = "OFF";
	out << token << endl;
	int nv = V.rows(), nf = F.rows(), ne = 0;
	out << nv << " " << nf << " " << ne << endl;
	for (int i = 0; i < nv; ++i) {
		out << V(i, 0) << " " << V(i, 1) << " " << V(i, 2) << endl;
	}
	for (int i = 0; i < nf; ++i) {
		int s = 3;
		out << s << " " << F(i, 0) << " " << F(i, 1) << " " << F(i, 2) << endl;
	}
}
