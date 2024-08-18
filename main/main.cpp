// main.cpp: This file contains the main function.
#pragma warning (disable : 4819)
#define Point_Size 15
#define Line_Width 1
#include<iostream>
#include<sys/stat.h>
#include<string>
#include<Eigen/Dense>
#include<SSS.h>
#include<ReadSSSCommand.h>
#include<igl/opengl/glfw/Viewer.h>
using namespace std;
using namespace Eigen;

Vector3d red(1, 0, 0), green(0, 1, 0), blue(0, 0, 1), yellow(1, 1, 0), purple(1, 0, 1), cyan(0, 1, 1), white(1, 1, 1), black(0, 0, 0);
Vector4f ared(1, 0, 0, 0), agreen(0, 1, 0, 0), ablue(0, 0, 1, 0), ayellow(1, 1, 0, 0), apurple(1, 0, 1, 0), acyan(0, 1, 1, 0), awhite(1, 1, 1, 0), ablack(0, 0, 0, 0);

string fileplace = "../Input/";

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

EnvironmentFeature env;
VectorXd SSSalpha(7), SSSbeta(7);
double envrange = 4;
heutype SSSheu = WIDTH;
SSSViewer viewer;

heutype toheutype(string word)
{
	if (word == "id")
		return BYID;
	if (word == "width")
		return WIDTH;
	if (word == "target")
		return TARGET;
	if (word == "gbf")
		return GBF;
	if (word == "dis")
		return DIS;
	if (word == "gbfdis")
		return GBFDIS;
	if (word == "widis" || word == "width_dis" || word == "widthdis")
		return WIDIS;
	return RAND;
}

bool exist(const string& filename)
{
	struct stat buffer;
	return (stat(filename.c_str(), &buffer) == 0);
}

// Decode enviroment command.
bool decode(string command)
{
	vector<string> words = split(command, ' ');
	if (words.size() < 1)
	{
		cout << "bash: not enough argument" << endl;
		return true;
	}
	switch (tofpred(words[0]))
	{
	case DEF:
	{
		if (words.size() < 2)
		{
			cout << "bash: not enough argument" << endl;
			return true;
		}
		switch (toftype(words[1]))
		{
		case TPOINT:
		{
			if (words.size() < 4)
			{
				cout << "bash: not enough arguments" << endl;
				return true;
			}
			string name = words[2];
			if (!isv3(words[3]))
			{
				cout << "bash: invalid vector type" << endl;
				return true;
			}
			Vector3d coord = stov3(words[3]);
			env.add_point(coord, name);
			cout << "Define point " << name << " as " << coord.transpose() << endl;
			return true;
		}
		case TEDGE:
		{
			if (words.size() < 5)
			{
				cout << "bash: not enough arguments" << endl;
				return true;
			}
			string name = words[2];
			if (isv3(words[3]))
			{
				Vector3d coord1 = stov3(words[3]);
				if (!isv3(words[4]))
				{
					cout << "bash: invalid vector type" << endl;
					return true;
				}
				Vector3d coord2 = stov3(words[4]);
				env.add_edge(coord1, coord2, name);
				cout << "Define edge " << name << " as the edge connecting " << coord1.transpose() << " and " << coord2.transpose() << endl;
			}
			else
			{
				string name1 = words[3];
				string name2 = words[4];
				if (env.add_edge(name1, name2, name) == NULL)
					cout << "Some point not found" << endl;
				else
					cout << "Define edge " << name << " as the edge connecting point " << name1 << " and point " << name2 << endl;
			}
			return true;
		}
		case TFACE:
		{
			if (words.size() < 6)
			{
				cout << "bash: not enough argument" << endl;
				return true;
			}
			string name = words[2];
			if (isv3(words[3]))
			{
				Vector3d coord1 = stov3(words[3]);
				if (!isv3(words[4]))
				{
					cout << "bash: invalid vector type" << endl;
					return true;
				}
				Vector3d coord2 = stov3(words[4]);
				if (!isv3(words[5]))
				{
					cout << "bash: invalid vector type" << endl;
					return true;
				}
				Vector3d coord3 = stov3(words[5]);
				env.add_face(coord1, coord2, coord3, name);
				cout << "Define face " << name << " as the triangle connecting " << coord1.transpose() << " and " << coord2.transpose() << " and " << coord3.transpose() << endl;
			}
			else
			{
				string name1 = words[3];
				string name2 = words[4];
				string name3 = words[5];
				if (env.add_face(name1, name2, name3, name) == NULL)
					cout << "Some point not found" << endl;
				else
					cout << "Define face " << name << " as the triangle connecting point " << name1 << " and point " << name2 << " and point " << name3 << endl;
			}
			return true;
		}
		case TMESH:
		{
			if (words.size() < 4)
			{
				cout << "bash: not enough argument" << endl;
				return true;
			}
			string name = words[2];
			string filename = words[3];
			if (filename.find("/") == string::npos && filename.find("\\") == string::npos)
				filename = fileplace + filename;
			if (!exist(filename))
			{
				cout << "bash: file " << filename << " not found" << endl;
				return true;
			}
			env.add_mesh(filename, name);
			cout << "Define mesh " << name << " from file " << filename << endl;
			return true;
		}
		default:cout << "bash: " << words[1] << ": command not found" << endl; return true;
		}
	}
	case DEL:
	{
		if (words.size() < 2)
		{
			cout << "bash: not enough argument" << endl;
			return true;
		}
		switch (toftype(words[1]))
		{
		case TPOINT:
		{
			if (words.size() < 3)
			{
				cout << "bash: not enough argument" << endl;
				return true;
			}
			string name = words[2];
			env.remove_point(name);
			cout << "Delete point " << name << endl;
			return true;
		}
		case TEDGE:
		{
			if (words.size() < 3)
			{
				cout << "bash: not enough argument" << endl;
				return true;
			}
			string name = words[2];
			env.remove_edge(name);
			cout << "Delete edge " << name << endl;
			return true;
		}
		case TFACE:
		{
			if (words.size() < 3)
			{
				cout << "bash: not enough argument" << endl;
				return true;
			}
			string name = words[2];
			env.remove_face(name);
			cout << "Delete face " << name << endl;
			return true;
		}
		case TMESH:
		{
			if (words.size() < 3)
			{
				cout << "bash: not enough argument" << endl;
				return true;
			}
			string name = words[2];
			env.remove_mesh(name);
			cout << "Delete mesh " << name << endl;
			return true;
		}
		default:cout << "bash: " << words[1] << ": command not found" << endl; return true;
		}
	}
	case SHOW:
	{
		if (words.size() < 2)
		{
			cout << "Show the environment." << endl;
			viewer.view_env(env, envrange, SSSalpha, SSSbeta);
			return true;
		}
		switch (toftype(words[1]))
		{
		case TPOINT:
		{
			env.show_point();
			return true;
		}
		case TEDGE:
		{
			env.show_edge();
			return true;
		}
		case TFACE:
		{
			env.show_face();
			return true;
		}
		case TMESH:
		{
			env.show_mesh();
			return true;
		}
		default:cout << "bash: " << words[1] << ": command not found" << endl; return true;
		}
	}
	case SET:
	{
		if (words.size() < 2)
		{
			cout << "bash: not enough argument" << endl;
			return true;
		}
		switch (toconfig(words[1]))
		{
		case ALPHA:
		{
			if (words.size() < 3)
			{
				cout << "bash: not enough argument" << endl;
				return true;
			}
			SSSalpha = stov7(words[2]);
			cout << "Set configuration alpha as " << SSSalpha.transpose() << endl;
			return true;
		}
		case BETA:
		{
			if (words.size() < 3)
			{
				cout << "bash: not enough argument" << endl;
				return true;
			}
			SSSbeta = stov7(words[2]);
			cout << "Set configuration beta as " << SSSbeta.transpose() << endl;
			return true;
		}
		case RANGE:
		{
			if (words.size() < 3)
			{
				cout << "bash: not enough argument" << endl;
				return true;
			}
			envrange = stod(words[2]);
			cout << "Set the environment range as " << "[-" << envrange << "," << envrange << "] * [-" << envrange << "," << envrange << "] * [-" << envrange << "," << envrange << "]" << endl;
			return true;
		}
		case HEURISTIC:
		{
			if (words.size() < 3)
			{
				cout << "bash: not enough argument" << endl;
				return true;
			}
			switch (toheutype(words[2]))
			{
			case RAND:SSSheu = RAND; cout << "Set heuristic by random." << endl; return true;
			case BYID:SSSheu = BYID; cout << "Set heuristic by the construction time of boxes." << endl; return true;
			case WIDTH:SSSheu = WIDTH; cout << "Set heuristic by the width of boxes." << endl; return true;
			case TARGET:SSSheu = TARGET; cout << "Set heuristic by the distance from boxes to target configuration." << endl; return true;
			case GBF:SSSheu = GBF; cout << "Set heuristic by greedy best first." << endl; return true;
			case DIS:SSSheu = DIS; cout << "Set heuristic by distance from features." << endl; return true;
			case GBFDIS:SSSheu = GBFDIS; cout << "Set heuristic by the product of GBF and DIS." << endl; return true;
			case WIDIS: SSSheu = WIDIS; cout << "Set heuristic by the mixed of WIDTH and DIS." << endl; return true;
			default:cout << "bash: " << words[2] << ": command not found" << endl; return true;
			}
		}
		default:cout << "bash: " << words[1] << ": command not found" << endl; return true;
		}
	}
	case RUN:
	{
		SE3Tree* T = new SE3Tree(MatrixId(Vector3d(-envrange, -envrange, -envrange), Vector3d(envrange, envrange, envrange)));
		DeltaPredicate* C = new DeltaPredicate();
		DeltaFeature* Omega = env.make_feature();
		double varepsilon = 0.05;
		if (words.size() > 1)
			varepsilon = stod(words[1]);
		SSS<VectorXd, SE3Box, SE3Tree, DeltaPredicate, DeltaFeature> SE3SSS(T, C);
		cout << "Simple find path algorithm with epsilon = " << varepsilon << endl;
		vector<VectorXd> path = SE3SSS.Find_Path(SSSalpha, SSSbeta, *Omega, varepsilon, SSSheu);
		if (path.empty())
			cout << "No Path!" << endl;
		else
			output_path(SSSalpha, SSSbeta, path);
		viewer.view_path(env, envrange, SSSalpha, SSSbeta, path);
		delete T;
		delete C;
		delete Omega;
		return true; 
	}
	case EXIT:
		return false;
	default:cout << "bash: " << words[0] << ": command not found" << endl; return true;
	}
}

// Command test.
void test(int argc, char* argv[])
{
	char command[1000];
	cout << "Delta robot find path algorithm by SSS framework." << endl;
	cout << "Demo -- Version 1.0" << endl;
	cout << "Default environment range: " << "[-" << envrange << "," << envrange << "] * [-" << envrange << "," << envrange << "] * [-" << envrange << "," << envrange << "]" << endl;
	cout << "Default alpha: (" << -envrange / 2.0 - 0.5 << "," << -envrange / 2.0 - 0.5 << "," << -envrange / 2.0 - 0.5 << "," << 1 << "," << 0 << "," << 0 << "," << 0 << ")" << endl;
	cout << "Default beta: (" << envrange / 2.0 + 0.5 << "," << envrange / 2.0 + 0.5 << "," << envrange / 2.0 + 0.5 << "," << 1 << "," << 0 << "," << 0 << "," << 0 << ")" << endl;
	cout << "Default epsilon: " << 0.05 << endl;
	SSSalpha << -envrange / 2.0 - 0.5, -envrange / 2.0 - 0.5, -envrange / 2.0 - 0.5, 1, 0, 0, 0;
	SSSbeta << envrange / 2.0 + 0.5, envrange / 2.0 + 0.5, envrange / 2.0 + 0.5, 1, 0, 0, 0;
	while (true)
	{
		cout << ">> ";
		cin.getline(command, 1000);
		if (!decode(string(command)))
			break;
	}
	cout << "Find path algorithm terminates." << endl;
}

// SSS test.
void test0(int argc, char* argv[])
{
	SE3Tree* T = new SE3Tree(MatrixId(Vector3d(-envrange, -envrange, -envrange), Vector3d(envrange, envrange, envrange)));
	DeltaPredicate* C = new DeltaPredicate();
	DeltaFeature Omega;
	/*
	// Do something with Omega.
	Point* ob1 = new Point(Vector3d(0, 0, 0));
	Point* ob2 = new Point(Vector3d(1, 0, 0));
	Point* ob3 = new Point(Vector3d(0, 0, 1));
	Edge* ob4 = new Edge(ob1, ob2);
	Edge* ob5 = new Edge(ob2, ob3);
	Edge* ob6 = new Edge(ob3, ob1);
	Triangle* ob = new Triangle(ob4, ob5, ob6);
	Omega.Vlist.push_back(ob1);
	Omega.Vlist.push_back(ob2);
	Omega.Vlist.push_back(ob3);
	Omega.Elist.push_back(ob4);
	Omega.Elist.push_back(ob5);
	Omega.Elist.push_back(ob6);
	Omega.Tlist.push_back(ob);
	// ******************* //
	*/
	double varepsilon = 0.1;
	SSS<VectorXd, SE3Box, SE3Tree, DeltaPredicate, DeltaFeature> SE3SSS(T, C);
	VectorXd alpha(7), beta(7);
	alpha << -127, -127, -127, 1, 0, 0, 0;
	beta << 127, 127, 127, 1, 0, 0, 0;
	/*C->set_feature(Omega);
	SE3SSS.Expand(SE3SSS.B->Root(), true);
	SE3SSS.B->Root()->child(1)->out();
	SE3SSS.Expand(SE3SSS.B->Root()->child(0));
	cout << endl;
	SE3SSS.B->Root()->child(0)->child(7)->out();
	cout << endl;
	SE3SSS.C->AppFp(SE3SSS.B->Root()->child(0)->child(7))->out();
	cout << endl << SE3SSS.C->classify(SE3SSS.B->Root()->child(0)->child(7)) << endl;*/
	output_path(alpha, beta, SE3SSS.Find_Path(alpha, beta, Omega, varepsilon));
}

// Viewer test.
void test1(int argc, char* argv[])
{
	// Inline mesh of a cube
	const Eigen::MatrixXd V = (Eigen::MatrixXd(8, 3) <<
		0.0, 0.0, 0.0,
		0.0, 0.0, 1.0,
		0.0, 1.0, 0.0,
		0.0, 1.0, 1.0,
		1.0, 0.0, 0.0,
		1.0, 0.0, 1.0,
		1.0, 1.0, 0.0,
		1.0, 1.0, 1.0).finished();
	const Eigen::MatrixXi F = (Eigen::MatrixXi(12, 3) <<
		0, 6, 4,
		0, 2, 6,
		0, 3, 2,
		0, 1, 3,
		2, 7, 6,
		2, 3, 7,
		4, 6, 7,
		4, 7, 5,
		0, 4, 5,
		0, 5, 1,
		1, 5, 7,
		1, 7, 3).finished();

	// Plot the mesh
	igl::opengl::glfw::Viewer viewer;
	viewer.data().set_mesh(V, F);
	viewer.data().set_face_based(true);
	viewer.launch();
}

// Load OFF test.
void test2(int argc, char* argv[])
{
	MatrixXd V;
	MatrixXi F;
	string filename;
	cout << "Input filename: " << endl;
	cin >> filename;
	read_OFF(filename, V, F);
	Mesh* new_mesh = new Mesh(V, F, true);
	cout << "Mesh built successfully!" << endl;
}

// Other tests.
void test3(int argc, char* argv[])
{
	while (true)
	{
		cout << "Input .off filename >> " << endl;
		string filename;
		cin >> filename;
		if (filename.find("/") == string::npos && filename.find("\\") == string::npos)
			filename = fileplace + filename;
		if (!exist(filename))
		{
			cout << "bash: file " << filename << " not found" << endl;
			break;
		}
		MatrixXd V;
		MatrixXi F;
		read_OFF(filename, V, F);
		Mesh* M = new Mesh(V, F);
		Vector3d p;
		cout << "Input point coordinate >> " << endl;
		cin >> p(0) >> p(1) >> p(2);
		if (M->inside(p))
			cout << p.transpose() << " is in the mesh " << filename << endl;
		else
			cout << p.transpose() << " is not in the mesh " << filename << endl;
	}
}

int main(int argc,char* argv[])
{
	int mode = -1;
	if (argc > 1)
		mode = stoi(argv[1]);
	switch (mode)
	{
	case 0:test0(argc, argv); break;
	case 1:test1(argc, argv); break;
	case 2:test2(argc, argv); break;
	case 3:test3(argc, argv); break;
	default:test(argc, argv);
	}
	return 0;
}