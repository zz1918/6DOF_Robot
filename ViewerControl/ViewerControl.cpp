// ViewerControl.cpp: This file controls the viewer of SSS main program.

#ifndef VIEWERCONTROL_H
#define VIEWERCONTROL_H

#pragma warning (disable : 4819)
#include<vector>
#include<map>
#include<Eigen/Dense>
#include<FeatureSet.h>
#include<IntroControl.h>
#include<igl/opengl/glfw/Viewer.h>
using namespace Eigen;

Vector3d red(1, 0, 0), green(0, 1, 0), blue(0, 0, 1), yellow(1, 1, 0), purple(1, 0, 1), cyan(0, 1, 1), white(1, 1, 1), black(0, 0, 0);
Vector4f ared(1, 0, 0, 0), agreen(0, 1, 0, 0), ablue(0, 0, 1, 0), ayellow(1, 1, 0, 0), apurple(1, 0, 1, 0), acyan(0, 1, 1, 0), awhite(1, 1, 1, 0), ablack(0, 0, 0, 0);

double extern Point_Size;
double extern Line_Width;

class SSSViewer
{
public:
	// Edges matrix 1.
	MatrixXd E1;
	// Edges matrix 2.
	MatrixXd E2;
	// Vertices matrix.
	MatrixXd Vs;
	// Edge color matrix.
	MatrixXd EC;
	// Vertex color matrix.
	MatrixXd VC;
	// Mesh V matrix.
	MatrixXd V;
	// Mesh F matrix.
	MatrixXi F;
	// Mesh color matrix.
	MatrixXd C;
	// Path edge matrix 1.
	MatrixXd PE1;
	// Path edge matrix 2.
	MatrixXd PE2;
	// Path edge color matrix.
	MatrixXd PEC;
	// Constructor.
	SSSViewer()
	{
		clear();
	}
	// Remove all environments.
	void clear()
	{
		E1.resize(0, 3);
		E2.resize(0, 3);
		Vs.resize(0, 3);
		EC.resize(0, 3);
		VC.resize(0, 3);
		V.resize(0, 3);
		F.resize(0, 3);
		C.resize(0, 3);
		PE1.resize(0, 3);
		PE2.resize(0, 3);
		PEC.resize(0, 3);
	}
	// Set the range box.
	void set_range(double range)
	{
		Vs.resize(8, 3);
		VC.resize(8, 3);
		for (int i = 0; i < 8; ++i)
			for (int j = 0; j < 3; ++j)
				Vs(i, j) = ((i >> j & 1) ? range : -range);
		for (int i = 0; i < 8; ++i)
			VC.row(i) = red;
		E1.resize(24, 3);
		E2.resize(24, 3);
		EC.resize(24, 3);
		for (int i = 0; i < 8; ++i)
		{
			E1.row(3 * i + 0) = Vs.row(i);
			E1.row(3 * i + 1) = Vs.row(i);
			E1.row(3 * i + 2) = Vs.row(i);
			E2.row(3 * i + 0) = Vs.row(i) - RowVector3d(Vs(i, 0) / range, 0, 0);
			E2.row(3 * i + 1) = Vs.row(i) - RowVector3d(0, Vs(i, 1) / range, 0);
			E2.row(3 * i + 2) = Vs.row(i) - RowVector3d(0, 0, Vs(i, 2) / range);
			EC.row(3 * i + 0) = blue;
			EC.row(3 * i + 1) = blue;
			EC.row(3 * i + 2) = blue;
		}
	}
	// Merge a mesh from V-F matrix.
	void merge(MatrixXd _V, MatrixXi _F)
	{
		V.conservativeResize(V.rows() + _V.rows(), V.cols());
		for (int i = 0; i < _V.rows(); ++i)
			V.row(V.rows() - _V.rows() + i) = _V.row(i);
		F.conservativeResize(F.rows() + _F.rows(), F.cols());
		for (int i = 0; i < _F.rows(); ++i)
			for (int j = 0; j < F.cols(); ++j)
				F(F.rows() - _F.rows() + i, j) = _F(i, j) + V.rows() - _V.rows();
		C.conservativeResize(C.rows() + _F.rows(), C.cols());
		for (int i = 0; i < _F.rows(); ++i)
			C.row(C.rows() - _F.rows() + i) = yellow;
	}
	// Add a mesh.
	void add(Mesh* M)
	{
		merge(M->Vertices(), M->Faces());
	}
	// Add a triangle.
	void add(Triangle* T)
	{
		Matrix3d _V;
		RowVector3i _F;
		for (int i = 0; i < 3; ++i)
			_V.row(i) = T->P(i)->p;
		_F = RowVector3i(0, 1, 2);
		merge(_V, _F);
	}
	// Add an edge.
	void add(Edge* E)
	{
		E1.conservativeResize(E1.rows() + 1, E1.cols());
		E1.row(E1.rows() - 1) = E->P(0)->p;
		E2.conservativeResize(E2.rows() + 1, E2.cols());
		E2.row(E2.rows() - 1) = E->P(1)->p;
		EC.conservativeResize(EC.rows() + 1, EC.cols());
		EC.row(EC.rows() - 1) = blue;
	}
	// Add a point.
	void add(Point* V)
	{
		Vs.conservativeResize(Vs.rows() + 1, Vs.cols());
		Vs.row(Vs.rows() - 1) = V->p;
		VC.conservativeResize(VC.rows() + 1, VC.cols());
		VC.row(VC.rows() - 1) = red;
	}
	// Add a footprint for time t.
	void add(VectorXd gamma, double time)
	{
		Matrix3d FpGamma = DeltaExFp(gamma).Fp;
		V.conservativeResize(V.rows() + 3, V.cols());
		for (int i = 0; i < 3; ++i)
			V.row(V.rows() - 3 + i) = FpGamma.row(i);
		F.conservativeResize(F.rows() + 1, F.cols());
		F.row(F.rows() - 1) = RowVector3i(V.rows() - 3, V.rows() - 2, V.rows() - 1);
		C.conservativeResize(C.rows() + 1, C.cols());
		C.row(C.rows() - 1) = (1 - time) * red + time * blue;
	}
	// Add the connection edge for consecutive footprints.
	void add(VectorXd gamma1, VectorXd gamma2)
	{
		PE1.conservativeResize(PE1.rows() + 1, PE1.cols());
		PE1.row(PE1.rows() - 1) = Vector3d(gamma1(0), gamma1(1), gamma1(2));
		PE2.conservativeResize(PE2.rows() + 1, PE2.cols());
		PE2.row(PE2.rows() - 1) = Vector3d(gamma2(0), gamma2(1), gamma2(2));
		PEC.conservativeResize(PEC.rows() + 1, PEC.cols());
		PEC.row(PEC.rows() - 1) = black;
	}
	// View environment.
	void view_env(EnvironmentFeature _env, double _envrange, VectorXd _alpha, VectorXd _beta)
	{
		igl::opengl::glfw::Viewer _viewer;
		clear();
		set_range(_envrange);
		for (map<string, Mesh*>::iterator it = _env.Mlist.begin(); it != _env.Mlist.end(); it++)
			add(it->second);
		for (map<string, Triangle*>::iterator it = _env.Tlist.begin(); it != _env.Tlist.end(); it++)
			add(it->second);
		for (map<string, Edge*>::iterator it = _env.Elist.begin(); it != _env.Elist.end(); it++)
			add(it->second);
		for (map<string, Point*>::iterator it = _env.Vlist.begin(); it != _env.Vlist.end(); it++)
			add(it->second);
		add(_alpha, 0);
		add(_beta, 1);
		_viewer.data().point_size = Point_Size;
		_viewer.data().line_width = Line_Width;
		_viewer.data().set_face_based(true);
		_viewer.data().double_sided = true;
		_viewer.data().show_lines = false;
		_viewer.data().set_mesh(V, F);
		_viewer.data().set_colors(C);
		_viewer.data().add_edges(E1, E2, EC);
		_viewer.data().add_points(Vs, VC);
		_viewer.core().align_camera_center(Vector3d(-_envrange - 1, 0, 0));
		_viewer.core().background_color = agreen;
		_viewer.launch();
	}
	// View path in environment.
	void view_path(EnvironmentFeature _env, double _envrange, VectorXd _alpha, VectorXd _beta, vector<VectorXd> _path)
	{
		igl::opengl::glfw::Viewer _viewer;
		clear();
		set_range(_envrange);
		for (map<string, Mesh*>::iterator it = _env.Mlist.begin(); it != _env.Mlist.end(); it++)
			add(it->second);
		for (map<string, Triangle*>::iterator it = _env.Tlist.begin(); it != _env.Tlist.end(); it++)
			add(it->second);
		for (map<string, Edge*>::iterator it = _env.Elist.begin(); it != _env.Elist.end(); it++)
			add(it->second);
		for (map<string, Point*>::iterator it = _env.Vlist.begin(); it != _env.Vlist.end(); it++)
			add(it->second);
		add(_alpha, 0);
		if (!_path.empty())
			add(_alpha, _path[0]);
		for (int i = 0; i < _path.size(); ++i)
		{
			add(_path[i], (i + 1) * 1.0 / (_path.size() + 1));
			if (i < _path.size() - 1)
				add(_path[i], _path[i + 1]);
			else
				add(_path[i], _beta);
		}
		add(_beta, 1);
		_viewer.data().point_size = Point_Size;
		_viewer.data().line_width = Line_Width;
		_viewer.data().set_face_based(true);
		_viewer.data().double_sided = true;
		_viewer.data().show_lines = false;
		_viewer.data().set_mesh(V, F);
		_viewer.data().set_colors(C);
		_viewer.data().add_edges(E1, E2, EC);
		_viewer.data().add_edges(PE1, PE2, PEC);
		_viewer.data().add_points(Vs, VC);
		_viewer.core().align_camera_center(Vector3d(-_envrange - 1, 0, 0));
		_viewer.core().background_color = agreen;
		_viewer.launch();
	}
};

SSSViewer viewer;

// Show environment.
void SSS_show()
{
	viewer.view_env(env, envrange, SSSalpha, SSSbeta);
}

// Show environment mode.
void show_only()
{
	SSS_intro();
	SSS_show();
}

#endif