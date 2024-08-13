// main.cpp: This file contains the main function.

#include<iostream>
#include<sys/stat.h>
#include<string>
#include<Eigen/Dense>
#include<SSS.h>
#include<ReadSSSCommand.h>
using namespace std;
using namespace Eigen;

EnvironmentFeature env;
VectorXd alpha(7), beta(7);
double envrange = 256;

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
		case POINT:
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
		case EDGE:
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
		case FACE:
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
		case MESH:
		{
			if (words.size() < 4)
			{
				cout << "bash: not enough argument" << endl;
				return true;
			}
			string name = words[2];
			string filename = words[3];
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
		case POINT:
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
		case EDGE:
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
		case FACE:
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
		case MESH:
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
			cout << "bash: not enough argument" << endl;
			return true;
		}
		switch (toftype(words[1]))
		{
		case POINT:
		{
			env.show_point();
			return true;
		}
		case EDGE:
		{
			env.show_edge();
			return true;
		}
		case FACE:
		{
			env.show_face();
			return true;
		}
		case MESH:
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
			alpha = stov7(words[2]);
			cout << "Set configuration alpha as " << alpha.transpose() << endl;
			return true;
		}
		case BETA:
		{
			if (words.size() < 3)
			{
				cout << "bash: not enough argument" << endl;
				return true;
			}
			beta = stov7(words[2]);
			cout << "Set configuration beta as " << beta.transpose() << endl;
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
		output_path(alpha, beta, SE3SSS.Find_Path(alpha, beta, *Omega, varepsilon, true));
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
	cout << "Default alpha: (-127,-127,-127,1,0,0,0)" << endl;
	cout << "Default beta: (127,127,127,1,0,0,0)" << endl;
	cout << "Default epsilon: " << 0.05 << endl;
	alpha << -127, -127, -127, 1, 0, 0, 0;
	beta << 127, 127, 127, 1, 0, 0, 0;
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
	output_path(alpha, beta, SE3SSS.Find_Path(alpha, beta, Omega, varepsilon, true));
}

void test1(int argc, char* argv[])
{}

void test2(int argc, char* argv[])
{}

void test3(int argc, char* argv[])
{}

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