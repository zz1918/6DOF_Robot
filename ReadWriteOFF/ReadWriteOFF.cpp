// ReadWriteOFF.cpp: this file defines the read and write from an .off file.

#ifndef READWRITEOFF_H
#define READWRITEOFF_H

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

// Read a cube mesh from a cube file
void read_CUBE(const string& filename, MatrixXd& V, MatrixXi& F)
{
	ifstream in(filename);
	string token;
	in >> token;
	Vector3d rc, rx, ry, rz;
	in >> rc(0) >> rc(1) >> rc(2);
	in >> rx(0) >> rx(1) >> rx(2);
	in >> ry(0) >> ry(1) >> ry(2);
	in >> rz(0) >> rz(1) >> rz(2);
	V.resize(8, 3);
	F.resize(12, 3);
	V.row(0) = rc;
	V.row(1) = rc + rx;
	V.row(2) = rc + ry;
	V.row(3) = rc + rx + ry;
	V.row(4) = rc + rz;
	V.row(5) = rc + rx + rz;
	V.row(6) = rc + ry + rz;
	V.row(7) = rc + rx + ry + rz;
	F << 0, 1, 5,
		0, 5, 4,
		0, 4, 6,
		0, 6, 2,
		0, 2, 3,
		0, 3, 1,
		4, 5, 6,
		6, 5, 7,
		7, 5, 3,
		3, 5, 1,
		6, 7, 3,
		6, 3, 2;
}

#endif