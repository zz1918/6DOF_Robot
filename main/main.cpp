#include<iostream>
#include<string>
#include<Eigen/Dense>
#include<SO3.h>
#include<BoxTree.h>
using namespace std;

class SE3
{
	
};

void test0(int argc, char* argv[])
{
	if (argc < 5)
		return;
	Vector4d a(stod(argv[1]), stod(argv[2]), stod(argv[3]), stod(argv[4]));
	SO3 A(a);
	cout << A.Q().transpose() << endl;
	cout << A.R() << endl;
	Vector3d v_A;
	double theta_A;
	A.Axis_Angle(v_A, theta_A);
	cout << v_A.transpose() << endl;
	cout << theta_A << endl;
	if (argc < 9)
		return;
	Vector4d b(stod(argv[5]), stod(argv[6]), stod(argv[7]), stod(argv[8]));
	SO3 B(b);
	cout << B.Q().transpose() << endl;
	cout << B.R() << endl;
	Vector3d v_B;
	double theta_B;
	B.Axis_Angle(v_B, theta_B);
	cout << v_B.transpose() << endl;
	cout << theta_B << endl;
	cout << (A == B) << endl;
}

void test1(int argc, char* argv[])
{
	if (argc < 4)
		return;
	int a[3];
	a[0] = stoi(argv[1]);
	a[1] = stoi(argv[2]);
	a[2] = stoi(argv[3]);
	SymG<3> t(a);
	cout << t << endl;
	cout << t.inverse() << endl;
	cout << t.act(0) << endl;
	cout << t.act(1) << endl;
	cout << t.act(2) << endl;
	cout << t.inverse().act(0) << endl;
	cout << t.inverse().act(1) << endl;
	cout << t.inverse().act(2) << endl;
	cout << t.binAct(0) << endl;
	cout << t.binAct(1) << endl;
	cout << t.binAct(2) << endl;
	cout << t.binAct(3) << endl;
	cout << t.binAct(4) << endl;
	cout << t.binAct(5) << endl;
	cout << t.binAct(6) << endl;
	cout << t.binAct(7) << endl;
}

void test2(int argc, char* argv[])
{
	SO3Tree T;
	T.cell(0)->subdivide();
	T.cell(1)->subdivide();
	T.cell(2)->subdivide();
	T.cell(3)->subdivide();
	T.cell(0)->child(0)->subdivide();
	for (int i = 0; i < 4; ++i)
		for (int j = 0; j < 8; ++j)
		{
			cout << "Cell " << i << ", child " << j << ":" << endl;
			if (T.cell(i)->child(j) == NULL)
			{
				cout << "No this child." << endl;
				continue;
			}
			cout << *(T.cell(i)->child(j));
			T.cell(i)->child(j)->show_neighbor();
			cout << endl;
		}
}

void test3(int argc, char* argv[])
{
	;
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