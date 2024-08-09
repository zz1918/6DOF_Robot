#include<iostream>
#include<string>
#include<Eigen/Dense>
#include<SSS.h>
using namespace std;
using namespace Eigen;

// SSS test.
void test(int argc, char* argv[])
{
	SE3Tree* T = new SE3Tree(MatrixId(Vector3d(-256, -256, -256), Vector3d(256, 256, 256)));
	DeltaPredicate* C = new DeltaPredicate();
	DeltaFeature Omega;

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

void test0(int argc, char* argv[])
{}

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