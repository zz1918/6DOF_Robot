#include<iostream>
#include<string>
#include<Eigen/Dense>
//#include<SO3.h>
//#include<interval.h>
//#include<BoxTree.h>
//#include<Solid.h>
#include<WtFp.h>
using namespace std;

void test0(int argc, char* argv[])
{
	MatrixXd V(4, 3);
	MatrixXi F(4, 3);
	V << -2, 0, 0,
		-1.5, 0, 0,
		-1.65377, -0.88236, 0,
		-0.20404, 0.22943, 0.10952;
	F << 0, 1, 2,
		2, 1, 3,
		3, 1, 0,
		0, 2, 3;
	Mesh* M = new Mesh(V, F);
	cout << V << endl;
	cout << F << endl;
	double w_cell = 0.1;
	MatrixId Bt(Vector3d(-w_cell / 2, -w_cell / 2, -w_cell / 2), Vector3d(w_cell / 2, w_cell / 2, w_cell / 2));
	MatrixId Br = Bt;
	WtFp footprint(Bt, Br, 0);
	cout << footprint.classify(M) << endl;
}

void test1(int argc, char* argv[])
{
	//Pyramid(Vector3d p, Vector3d q, Vector3d r, Vector3d dir, double l1, double l2, double l3)
	Pyramid* T = new Pyramid(Vector3d(0, 0, 0), Vector3d(0, 4, 0), Vector3d(4, 0, 0), Vector3d(0, 0, 2), 1, 1, 1);
	Point* P = new Point(Vector3d(0, 0, 6));
	cout << "Point" << endl;
	cout << T->Sep(P) << endl;
	//******************
	/*Edge* edge = new Edge(new Point(Vector3d(10, 10, 10)), new Point(Vector3d(11, 11, 11)));
	cout << "Edge" << endl;
	cout << pyramid->Sep(edge) << endl;*/
	//******************
	for (int i = 0; i < 3; ++i)
	{
		cout << "T->V(" << i << ")->Sep(P): " << T->V(i)->Sep(P, 0) << endl;
		for (int j = 0; j < 4; ++j)
			cout << "T->V(" << i << ")->E(" << j << ")->Sep(P) : " << T->V(i)->E(j)->Sep(P, 0) << endl;
	}
	cout << T->V(0)->E(0)->P(0)->p.transpose() << endl;
	cout << T->V(0)->E(0)->P(1)->p.transpose() << endl;
	cout << T->V(0)->E(0)->int_Sep(P, 0) << endl;
	cout << T->V(0)->E(0)->projects_on(P) << endl;
}

void test2(int argc, char* argv[])
{
	Point* P1 = new Point(Vector3d(-2, -2, -2));
	Point* P2 = new Point(Vector3d(-2, -2, 2));
	Point* P3 = new Point(Vector3d(-2, 2, 2));
	Edge* E1 = new Edge(P1, P2);
	Edge* E2 = new Edge(P2, P3);
	Edge* E3 = new Edge(P3, P1);
	Triangle* T = new Triangle(E1, E2, E3);
	Point* Q1 = new Point(Vector3d(1, 0, 1));
	Point* Q2 = new Point(Vector3d(1, 0, -1));
	Edge* F = new Edge(Q1, Q2);
	cout << F->Sep(T, 0) << endl;
	for (int i = 0; i < 3; ++i)
		cout << T->E(i)->Sep(F, 0) << endl;
	for (int i = 0; i < 2; ++i)
		cout << T->Sep(F->P(i), 0) << endl;
	cout << E3->int_Sep(F, 0) << endl;
	cout<<E3->line_distance(F)<<endl;
}

void test3(int argc, char* argv[])
{
	double w_cell = 0.1;
	MatrixId Bt(Vector3d(-w_cell / 2, -w_cell / 2, -w_cell / 2), Vector3d(w_cell / 2, w_cell / 2, w_cell / 2));
	MatrixId Br = Bt;
	WtFp footprint(Bt, Br, 0);
	Point* P1 = new Point(Vector3d(-10, -10, -10));
	Point* P2 = new Point(Vector3d(-10, -10, 10));
	Point* P3 = new Point(Vector3d(-10, 10, 10));
	Edge* E1 = new Edge(P1, P2);
	Edge* E2 = new Edge(P2, P3);
	Edge* E3 = new Edge(P3, P1);
	Triangle* T = new Triangle(E1, E2, E3);
	cout << "dB: " << footprint.dB << endl;
	cout << "rB: " << footprint.rB << endl;
	cout << "R: " << footprint.dB + footprint.rB << endl;
	cout << "r: " << footprint.rB << endl;
	cout << "Is singular? " << footprint.singular << endl;
	cout << "Cyl^+: " << footprint.SegAB->Sep(T, footprint.dB + footprint.rB) << endl;
	cout << "H^+: " << footprint.H->Sep(T, 0) << endl;
	cout << "IccA: " << footprint.IccA->Sep(T, footprint.rB) << endl;
	cout << "IccB: " << footprint.IccB->Sep(T, footprint.rB) << endl;
}

int main(int argc,char* argv[])
{
	int mode = 0;
	if (argc > 1)
		mode = stoi(argv[1]);
	switch (mode)
	{
	case 0:test0(argc, argv); break;
	case 1:test1(argc, argv); break;
	case 2:test2(argc, argv); break;
	case 3:test3(argc, argv); break;
	}
	return 0;
}