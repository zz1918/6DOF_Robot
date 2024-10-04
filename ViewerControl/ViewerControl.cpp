// ViewerControl.cpp: This file controls the viewer of SSS main program.

#ifndef VIEWERCONTROL_H
#define VIEWERCONTROL_H

#pragma warning (disable : 4819)
#include<vector>
#include<map>
#include<set>
#include<SE3Box.h>
#include<Graph.h>
#include<Eigen/Dense>
#include<FeatureSet.h>
#include<IntroControl.h>
#include<igl/opengl/glfw/Viewer.h>
using namespace Eigen;

Vector3d red(1, 0, 0), green(0, 1, 0), blue(0, 0, 1), yellow(1, 1, 0), purple(1, 0, 1), cyan(0, 1, 1), white(1, 1, 1), black(0, 0, 0), grey(0.5, 0.5, 0.5);
Vector4f ared(1, 0, 0, 0), agreen(0, 1, 0, 0), ablue(0, 0, 1, 0), ayellow(1, 1, 0, 0), apurple(1, 0, 1, 0), acyan(0, 1, 1, 0), awhite(1, 1, 1, 0), ablack(0, 0, 0, 0), agrey(0.5, 0.5, 0.5, 0);

double extern Point_Size;
double extern Line_Width;
bool extern box_draw_strategy;

Vector3d to_color(Vcolor c)
{
	switch (c)
	{
	case GREEN:return green;
	case YELLOW:return yellow;
	case RED:return red;
	case BLACK:return black;
	case GREY:return grey;
	default:return blue;
	}
}

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
	// Environment.
	EnvironmentFeature env;
	// Environment range.
	double envrange;
	// Initial configuration.
	VectorXd alpha;
	// Target configuration.
	VectorXd beta;
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
	void set_range(double range, Vector3d bColor = blue)
	{
		MatrixXd Corner(8, 3);
		for (int i = 0; i < 8; ++i)
			for (int j = 0; j < 3; ++j)
				Corner(i, j) = ((i >> j & 1) ? range : -range);
		for (int i = 0; i < 8; ++i)
			add(new Point(Corner.row(i)), red);
		add(new Edge(Corner.row(0), Corner.row(1)), bColor);
		add(new Edge(Corner.row(1), Corner.row(3)), bColor);
		add(new Edge(Corner.row(3), Corner.row(2)), bColor);
		add(new Edge(Corner.row(2), Corner.row(0)), bColor);
		add(new Edge(Corner.row(4), Corner.row(5)), bColor);
		add(new Edge(Corner.row(5), Corner.row(7)), bColor);
		add(new Edge(Corner.row(7), Corner.row(6)), bColor);
		add(new Edge(Corner.row(6), Corner.row(4)), bColor);
		add(new Edge(Corner.row(0), Corner.row(4)), bColor);
		add(new Edge(Corner.row(5), Corner.row(1)), bColor);
		add(new Edge(Corner.row(2), Corner.row(6)), bColor);
		add(new Edge(Corner.row(7), Corner.row(3)), bColor);
		/*
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
		}*/
	}
	// Merge a mesh from V-F matrix.
	void merge(MatrixXd _V, MatrixXi _F, Vector3d _C = white)
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
			C.row(C.rows() - _F.rows() + i) = _C;
	}
	// Add a mesh.
	void add(Mesh* M, Vector3d _C = white)
	{
		merge(M->Vertices(), M->Faces(), _C);
	}
	// Add a triangle.
	void add(Triangle* T, Vector3d _C = white)
	{
		Matrix3d _V;
		RowVector3i _F;
		for (int i = 0; i < 3; ++i)
			_V.row(i) = T->P(i)->p;
		_F = RowVector3i(0, 1, 2);
		merge(_V, _F, _C);
	}
	// Add an edge.
	void add(Edge* E, Vector3d _C = blue)
	{
		E1.conservativeResize(E1.rows() + 1, E1.cols());
		E1.row(E1.rows() - 1) = E->P(0)->p;
		E2.conservativeResize(E2.rows() + 1, E2.cols());
		E2.row(E2.rows() - 1) = E->P(1)->p;
		EC.conservativeResize(EC.rows() + 1, EC.cols());
		EC.row(EC.rows() - 1) = _C;
	}
	// Add a point.
	void add(Point* P, Vector3d _C = red)
	{
		Vs.conservativeResize(Vs.rows() + 1, Vs.cols());
		Vs.row(Vs.rows() - 1) = P->p;
		VC.conservativeResize(VC.rows() + 1, VC.cols());
		VC.row(VC.rows() - 1) = _C;
	}
	// Add a footprint.
	void add(VectorXd gamma, Vector3d _C = green)
	{
		Matrix3d FpGamma = DeltaExFp(gamma).Fp;
		V.conservativeResize(V.rows() + 3, V.cols());
		for (int i = 0; i < 3; ++i)
			V.row(V.rows() - 3 + i) = FpGamma.row(i);
		F.conservativeResize(F.rows() + 1, F.cols());
		F.row(F.rows() - 1) = RowVector3i(V.rows() - 3, V.rows() - 2, V.rows() - 1);
		C.conservativeResize(C.rows() + 1, C.cols());
		C.row(C.rows() - 1) = _C;
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
	// Set environment.
	void set_env(EnvironmentFeature _env, double _envrange, VectorXd _alpha, VectorXd _beta)
	{
		clear();
		env = _env;
		envrange = _envrange;
		alpha = _alpha;
		beta = _beta;
		set_range(envrange);
		for (map<string, Mesh*>::iterator it = env.Mlist.begin(); it != env.Mlist.end(); it++)
			add(it->second, white);
		for (map<string, Triangle*>::iterator it = env.Tlist.begin(); it != env.Tlist.end(); it++)
			add(it->second, white);
		for (map<string, Edge*>::iterator it = env.Elist.begin(); it != env.Elist.end(); it++)
			add(it->second, white);
		for (map<string, Point*>::iterator it = env.Vlist.begin(); it != env.Vlist.end(); it++)
			add(it->second, white);
		add(alpha, 0.0);
		add(beta, 1.0);
	}
	// View everyting.
	void view(Vector4f bColor = agrey)
	{
		// The viewer.
		igl::opengl::glfw::Viewer _viewer;
		_viewer.data().point_size = Point_Size;
		_viewer.data().line_width = Line_Width;
		_viewer.data().set_face_based(true);
		_viewer.data().double_sided = true;
		_viewer.data().show_lines = false;
		_viewer.data().set_mesh(V, F);
		_viewer.data().set_colors(C);
		_viewer.data().add_edges(E1, E2, EC);
		_viewer.data().add_points(Vs, VC);
		_viewer.core().align_camera_center(Vector3d(-envrange - 1, 0, 0));
		_viewer.core().background_color = bColor;
		/*
		_viewer.callback_key_down = [&](decltype(_viewer)&, unsigned int k, int m)
			{
				switch (char(k))
				{
				case 'c':_viewer.close(); return true;
				default:return true;
				}
			};*/
		_viewer.launch();
	}
	// View environment.
	void view_env(Vector4f bColor = agrey)
	{
		view(bColor);
	}
	// Add a FREE box.
	void add_box(SE3Box* b, Vector3d bColor = green)
	{
		if (box_draw_strategy)
		{
			R3Box* bt = b->BT();
			MatrixXd Corner(8, 3);
			for (int i = 0; i < 8; ++i)
				for (int j = 0; j < 3; ++j)
					Corner(i, j) = ((i >> j & 1) ? bt->Range().max()(j) : bt->Range().min()(j));
			add(new Edge(Corner.row(0), Corner.row(1)), bColor);
			add(new Edge(Corner.row(1), Corner.row(3)), bColor);
			add(new Edge(Corner.row(3), Corner.row(2)), bColor);
			add(new Edge(Corner.row(2), Corner.row(0)), bColor);
			add(new Edge(Corner.row(4), Corner.row(5)), bColor);
			add(new Edge(Corner.row(5), Corner.row(7)), bColor);
			add(new Edge(Corner.row(7), Corner.row(6)), bColor);
			add(new Edge(Corner.row(6), Corner.row(4)), bColor);
			add(new Edge(Corner.row(0), Corner.row(4)), bColor);
			add(new Edge(Corner.row(5), Corner.row(1)), bColor);
			add(new Edge(Corner.row(2), Corner.row(6)), bColor);
			add(new Edge(Corner.row(7), Corner.row(3)), bColor);
		}
		else
		{
			VectorXd gamma = b->center();
			add(gamma, bColor);
		}
	}
	// Set the global graph.
	void set_graph(Graph<SE3Box> G)
	{
		if (G.Vlist.empty())
			return;
		for (auto it = G.Vlist.begin(); it != G.Vlist.end(); ++it)
			add_box(G.node(*it)->content, to_color(G.node(*it)->colors[0]));
	}
	// Set the global graph.
	void set_graph(Graph<SE3Box> G, set<int> forbid)
	{
		if (G.Vlist.empty())
			return;
		for (auto it = G.Vlist.begin(); it != G.Vlist.end(); ++it)
			if (forbid.find(*it) == forbid.end())
				add_box(G.node(*it)->content, to_color(G.node(*it)->colors[0]));
	}
	// Set the FREE graph.
	void set_graph(vector<SE3Box*> G, Vector3d gColor = green)
	{
		for (auto it = G.begin(); it != G.end(); ++it)
			add_box(*it, gColor);
	}
	// View FREE graph in environment.
	void view_graph(vector<SE3Box*> G, Vector3d gColor = green, Vector4f bColor = agrey)
	{
		set_graph(G, gColor);
		view(bColor);
	}
	// Set the path.
	void set_path(vector<VectorXd> _path)
	{
		if (!_path.empty())
			add(alpha, _path[0]);
		for (int i = 0; i < _path.size(); ++i)
		{
			add(_path[i], (i + 1) * 1.0 / (_path.size() + 1));
			if (i < _path.size() - 1)
				add(_path[i], _path[i + 1]);
			else
				add(_path[i], beta);
		}
	}
	// View path in environment.
	void view_path(vector<VectorXd> _path, Vector4f bColor = agrey)
	{
		set_path(_path);
		view(bColor);
	}
};

SSSViewer viewer;

// Show environment.
void SSS_show()
{
	viewer.set_env(env, envrange, SSSalpha, SSSbeta);
	viewer.view_env();
}

// Show environment mode.
void show_only()
{
	SSS_intro();
	SSS_show();
}

#endif