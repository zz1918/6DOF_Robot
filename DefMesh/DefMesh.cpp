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
	void add_scale(double t,int c)
	{
		Matrix4d _T = Matrix4d::Identity();
		_T.block<1, 1>(c, c) *= t;
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

class OffMergeTool
{
	MatrixXd V;
	MatrixXi F;
public:
	OffMergeTool()
	{
		V.resize(0, 3);
		F.resize(0, 3);
	}
	void add_mesh(string input)
	{
		MatrixXd _V;
		MatrixXi _F;
		read_OFF(input, _V, _F);
		F.conservativeResize(F.rows() + _F.rows(), F.cols());
		for (int i = 0; i < _F.rows(); ++i)
			for (int j = 0; j < F.cols(); ++j)
				F(F.rows() - _F.rows() + i, j) = _F(i, j) + V.rows();
		V.conservativeResize(V.rows() + _V.rows(), V.cols());
		for (int i = 0; i < _V.rows(); ++i)
			V.row(V.rows() - _V.rows() + i) = _V.row(i);
	}
	void make_OFF(string output)
	{
		write_OFF(output, V, F);
	}
};

class OffSplitTool
{
	MatrixXd V, V1, V2;
	MatrixXi F, F1, F2;
public:
	OffSplitTool()
	{
		V.resize(0, 3);
		V1.resize(0, 3);
		V2.resize(0, 3);
		F.resize(0, 3);
		F1.resize(0, 3);
		F2.resize(0, 3);
	}
	void def_mesh(string input)
	{
		read_OFF(input, V, F);
	}
	void split(int V1size)
	{
		V1.resize(V1size, 3);
		V2.resize(V.rows() - V1size, 3);
		for (int i = 0; i < V1size; ++i)
			V1.row(i) = V.row(i);
		for (int i = V1size; i < V.rows(); ++i)
			V2.row(i - V1size) = V.row(i);
		int F1size = 0;
		for (int i = 0; i < F.rows(); ++i)
			if (F(i, 0) < V1size)
				F1size++;
		F1.resize(F1size, 3);
		F2.resize(F.rows() - F1size, 3);
		int i1 = 0, i2 = 0;
		for (int i = 0; i < F.rows(); ++i)
			if (F(i, 0) < V1size)
				F1.row(i1++) = F.row(i);
			else
			{
				for (int j = 0; j < F.cols(); ++j)
					F2(i2, j) = F(i, j) - V1size;
				i2++;
			}
	}
	void make_OFF(string output1, string output2)
	{
		write_OFF(output1, V1, F1);
		write_OFF(output2, V2, F2);
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
	double scaleX = stod(argv[15]);
	double scaleY = stod(argv[16]);
	double scaleZ = stod(argv[17]);
	OffLinearTool OLT;
	OLT.add_translation(-rotC);
	OLT.add_scale(scale);
	OLT.add_scale(scaleX, 0);
	OLT.add_scale(scaleY, 1);
	OLT.add_scale(scaleZ, 2);
	OLT.add_rotation(rotA, rotT);
	OLT.add_translation(rotC);
	OLT.add_translation(trans);
	OLT.make_OFF(input, output);
}

// Merge two meshes.
void merge_mesh(int argc, char* argv[])
{
	string input1 = fileplace + string(argv[2]) + fileformat;
	string input2 = fileplace + string(argv[3]) + fileformat;
	string output = fileplace + string(argv[4]) + fileformat;
	OffMergeTool OMT;
	OMT.add_mesh(input1);
	OMT.add_mesh(input2);
	OMT.make_OFF(output);
}

// Split a mesh into two parts by the number of vertices in the first mesh.
void split_mesh(int argc, char* argv[])
{
	string input = fileplace + string(argv[2]) + fileformat;
	string output1 = fileplace + string(argv[3]) + fileformat;
	string output2 = fileplace + string(argv[4]) + fileformat;
	int V1size = stoi(argv[5]);
	OffSplitTool OST;
	OST.def_mesh(input);
	OST.split(V1size);
	OST.make_OFF(output1, output2);
}

#endif