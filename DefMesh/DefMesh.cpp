// DefMesh.cpp: This file can rotate/translate/scale given meshes by linear transformations.

#ifndef DEFMESH_H
#define DEFMESH_H

#include<string>
#include<Eigen/Dense>
#include<SO3.h>
#include<ReadWriteOFF.h>
using namespace std;
using namespace Eigen;

class OffLinearTool
{
	Matrix4d T;			// An SE(3) matrix that deal with all linear transformations.
public:
	OffLinearTool()
	{
		T = Matrix4d::Identity();
	}
	void add_translation(Vector3d v)
	{
		Matrix4d _T = Matrix4d::Identity();
		_T.block<3, 1>(0, 3) = v;
		T = _T * T;
	}
	void add_rotation(double phi, double theta, double xi)
	{
		Matrix4d _T = Matrix4d::Identity();
		_T.block<3, 3>(0, 0) = SO3(phi, theta, xi).R();
		T = _T * T;
	}
	void add_rotation(Vector3d v, double tau)
	{
		Matrix4d _T = Matrix4d::Identity();
		_T.block<3, 3>(0, 0) = SO3(v, tau).R();
		T = _T * T;
	}
	void add_scale(double t)
	{
		Matrix4d _T = Matrix4d::Identity();
		_T.block<3, 3>(0, 0) *= t;
		T = _T * T;
	}
	void make_OFF(string input, string output)
	{
		MatrixXd V;
		MatrixXi F;
		read_OFF(input, V, F);
		V.conservativeResize(V.rows(), V.cols() + 1);
		V.col(V.cols() - 1) = VectorXd::Ones(V.rows());
		V = V * T.transpose();
		V.conservativeResize(V.rows(), V.cols() - 1);
		write_OFF(output, V, F);
	}
};

extern string fileplace;
extern string fileformat;

// Translate, rotate or scale meshes.
void set_mesh(int argc, char* argv[])
{
	string input = fileplace + string(argv[2]) + fileformat;
	string output = fileplace + string(argv[3]) + fileformat;
	Vector3d rotA(stod(argv[4]), stod(argv[5]), stod(argv[6]));
	double rotT(stod(argv[7]));
	Vector3d rotC(stod(argv[8]), stod(argv[9]), stod(argv[10]));
	Vector3d trans(stod(argv[11]), stod(argv[12]), stod(argv[13]));
	double scale = stod(argv[14]);
	OffLinearTool OLT;
	OLT.add_translation(-rotC);
	OLT.add_scale(scale);
	OLT.add_rotation(rotA, rotT);
	OLT.add_translation(rotC);
	OLT.add_translation(trans);
	OLT.make_OFF(input, output);
}

#endif