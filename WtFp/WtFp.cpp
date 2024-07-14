// WtFp.cpp: This file implements the special Sigma_2 decomposition 
//			of the approximate footprint of the Delta robot.

#ifndef WTFP_H
#define WTFP_H

#include<vector>
#include<Eigen/Dense>
#include<Interval.h>
#include<SO3.h>
#include<Solid.h>
#include<bimap.h>
using namespace std;
using namespace Eigen;

enum pvalue { MIXED, FREE, STUCK, UNKNOWN };

#define pA Vector3d(1, 0, 0)
#define pB Vector3d(0, 1, 0)

// *********************Help functions********************* //

// Min of 3 numbers.
double minof(double a, double b, double c)
{
	return min(min(a, b), c);
}
// Max of 3 numbers.
double maxof(double a, double b, double c)
{
	return max(max(a, b), c);
}
// Bounding box of an edge.
MatrixId B_Box(Edge* f)
{
	return MatrixId(f->P(0)->p, f->P(1)->p);
}
// Bounding box of a triangle.
MatrixId B_Box(Triangle* f)
{
	Vector3d m, M;
	m(0) = minof(f->P(0)->p(0), f->P(1)->p(0), f->P(2)->p(0));
	m(1) = minof(f->P(0)->p(1), f->P(1)->p(1), f->P(2)->p(1));
	m(2) = minof(f->P(0)->p(2), f->P(1)->p(2), f->P(2)->p(2));
	M(0) = maxof(f->P(0)->p(0), f->P(1)->p(0), f->P(2)->p(0));
	M(1) = maxof(f->P(0)->p(1), f->P(1)->p(1), f->P(2)->p(1));
	M(2) = maxof(f->P(0)->p(2), f->P(1)->p(2), f->P(2)->p(2));
	return MatrixId(m, M);
}


// Psuedo feature Mesh, which can only be detected by WtFp but cannot generally do parametric queries.
class Mesh
{
	MatrixXd V;
	MatrixXi F;
	bimap<int, int, int> VE;
	MatrixXi EV;

	void edge_topology(MatrixXi F)
	{
		VE.clear();
		EV.resize(edges.size(), 2);
		int edge_size = 0;
		for (int i = 0; i < F.rows(); ++i)
			for (int j = 0; j < 3; ++j)
			{
				VE.insert(F(i, j), F(i, (j + 1) % 3), edge_size + 1);
				EV(edge_size, 0) = F(i, j);
				EV(edge_size, 1) = F(i, (j + 1) % 3);
				edge_size++;
			}
	}
	int Eid(int Vid1, int Vid2)
	{
		return VE.coeff(Vid1, Vid2) - 1;
	}

	// Bounding box of mesh for quick exclusion.
	MatrixId bbox;
public:
	vector<Point*> corners;
	vector<Edge*> edges;
	vector<Triangle*> faces;
	Mesh(MatrixXd _V, MatrixXi _F, bool show = false)
	{
		V = _V;
		F = _F;
		// Construct point features.
		for (int i = 0; i < V.rows(); ++i)
		{
			Point* P = new Point(V.row(i));
			corners.push_back(P);
		}
		if (show)
			cout << "Points constructed." << endl;
		// Construct edge and triangle features.
		for (int i = 0; i < F.rows(); ++i)
		{
			Edge* E0 = new Edge(corners[F(i, 0)], corners[F(i, 1)]);
			Edge* E1 = new Edge(corners[F(i, 1)], corners[F(i, 2)]);
			Edge* E2 = new Edge(corners[F(i, 2)], corners[F(i, 0)]);
			Triangle* T = new Triangle(E0, E1, E2);
			edges.push_back(E0);
			edges.push_back(E1);
			edges.push_back(E2);

			faces.push_back(T);
		}
		if (show)
			cout << "Faces constructed." << endl;
		// Set edge and inv-edge relations.
		edge_topology(F);
		if (show)
			cout << "Topology built." << endl;
		for (int i = 0; i < edges.size(); ++i)
			if (Eid(EV(i, 1), EV(i, 0)))
				edges[i]->set_inv(edges[Eid(EV(i, 1), EV(i, 0))]);
		if (show)
			cout << "Topology constructed." << endl;
		bbox = MatrixId(Vector3d(V.col(0).minCoeff(), V.col(1).minCoeff(), V.col(2).minCoeff()), Vector3d(V.col(0).maxCoeff(), V.col(1).maxCoeff(), V.col(2).maxCoeff()));
		if (show)
			cout << "Bounding box constructed." << endl;
	}
	// If a point p is positive to the normals of the triangle mesh.
	bool inside(Vector3d p)
	{
		for (int i = 0; i < faces.size(); ++i)
			if (faces[i]->N().dot(faces[i]->P(0)->D(p)) < 0)
				return false;
		return true;
	}
	// Return the bbox.
	MatrixId BBox()
	{
		return bbox;
	}
};

// Approximate footprint for Delta robot.
class DeltaWtFp
{
public:
	// \icc_\pA and \icc_\pB.
	IceCream* IccA, * IccB;
	// \seg_{\pA\pB}.
	Edge* SegAB;
	// Pyramid^+.
	Pyramid* H;
	// r(B) and d(B) for box B.
	double rB, dB;
	// Center of box Bt.
	Vector3d mB;
	// If the approximate footprint is reduced (OO is inside ball SA+).
	bool singular;
	// Centers of SA and SB (pure rotation part).
	Vector3d opA, opB;
	// Bounding box of approximate footprint for quick exclusion.
	MatrixId bbox;

	// Finds the 4D corners of Br.
	vector<Vector4d> make_corners(MatrixId Br, int wxyz)
	{
		Vector3d LCorners[8];
		vector<Vector4d> Corners;
		Vector3d MIN = Br.min();
		Vector3d MAX = Br.max();
		LCorners[0] = Vector3d(MIN(0), MIN(1), MIN(2));
		LCorners[1] = Vector3d(MAX(0), MIN(1), MIN(2));
		LCorners[2] = Vector3d(MIN(0), MAX(1), MIN(2));
		LCorners[3] = Vector3d(MAX(0), MAX(1), MIN(2));
		LCorners[4] = Vector3d(MIN(0), MIN(1), MAX(2));
		LCorners[5] = Vector3d(MAX(0), MIN(1), MAX(2));
		LCorners[6] = Vector3d(MIN(0), MAX(1), MAX(2));
		LCorners[7] = Vector3d(MAX(0), MAX(1), MAX(2));
		switch (wxyz)
		{
		case 0:for (int i = 0; i < 8; ++i)
			Corners.push_back(Vector4d(1, LCorners[i](0), LCorners[i](1), LCorners[i](2)));
			break;
		case 1:for (int i = 0; i < 8; ++i)
			Corners.push_back(Vector4d(LCorners[i](0), 1, LCorners[i](1), LCorners[i](2)));
			break;
		case 2:for (int i = 0; i < 8; ++i)
			Corners.push_back(Vector4d(LCorners[i](0), LCorners[i](1), 1, LCorners[i](2)));
			break;
		case 3:for (int i = 0; i < 8; ++i)
			Corners.push_back(Vector4d(LCorners[i](0), LCorners[i](1), LCorners[i](2), 1));
			break;
		default:return Corners;
		}
		return Corners;
	}
	// Footprint of A under Q.
	Vector3d FpA(Vector4d Q)
	{
		return SO3(Q).act_on(pA);
	}
	// Footprint of B under Q.
	Vector3d FpB(Vector4d Q)
	{
		return SO3(Q).act_on(pB);
	}
	// Center of gravity of ps.
	Vector3d CenGra(vector<Vector3d> ps)
	{
		Vector3d sum = Vector3d::Zero();
		for (int i = 0; i < ps.size(); ++i)
			sum += ps[i];
		return sum / ps.size();
	}
	// Maximum distance from points ps to point q.
	double MaxDis(vector<Vector3d> ps, Vector3d q)
	{
		double MD = 0;
		for (int i = 0; i < ps.size(); ++i)
		{
			double dis = (ps[i] - q).norm();
			if (MD < dis)
				MD = dis;
		}
		return MD;
	}
	// Construction of the pyramid given by center of balls and radiuses.
	Pyramid* Hull()
	{
		// Computations for the pyramid.
		Vector3d uAB = (opB - opA).normalized();
		Vector3d uOA = opA.normalized();
		Vector3d uOB = opB.normalized();
		double R = dB + rB, r = rB;
		double R2 = R * R, r2 = r * r;

		double dA = opA.norm() * R / (R - r);
		double dOA = opA.norm() * r / (R - r);
		double dB = opB.norm() * R / (R - r);
		double dOB = opB.norm() * r / (R - r);
		double dA2 = dA * dA;
		double dB2 = dB * dB;
		double dOA2 = dOA * dOA;
		double dOB2 = dOB * dOB;
		Vector3d dir = uAB.cross(uOA);
		double dd = dir.norm();
		double dd2 = dd * dd;

		Vector3d cA = opA + mB + ((uAB.dot(uOA) * R2) / (dd2 * dA)) * uAB - (R2 / (dd2 * dA)) * uOA;
		Vector3d cB = opB + mB + ((uAB.dot(uOB) * R2) / (dd2 * dB)) * uAB - (R2 / (dd2 * dB)) * uOB;
		Vector3d cO = mB + (r2 / dd2) * (uOA.dot(uOB) / dOB - 1 / dOA) * uOA + (r2 / dd2) * (uOA.dot(uOB) / dOA - 1 / dOB) * uOB;

		double lA = R * sqrt(dd2 * dA2 - R2) / (dd2 * dA);
		double lB = R * sqrt(dd2 * dB2 - R2) / (dd2 * dB);
		double lO = (r2 / dd2) * sqrt((dd2 / r2) - (1 / dOA2 - 2 * uOA.dot(uOB) / (dOA * dOB) + 1 / dOB2));

		return new Pyramid(cA - lA * dir, cB - lB * dir, cO - lO * dir, dir, 2 * lA, 2 * lB, 2 * lO);
	}
	// Construct the bounding box of the WtFp
	MatrixId BBox()
	{
		if (singular)
		{
			// If the approximate footprint is singular, then the bounding box is given by the two balls.
			double R = dB + rB;
			Vector3d bBA = mB + opA + Vector3d(R, R, R), bbA = mB + opA - Vector3d(R, R, R);
			Vector3d bBB = mB + opB + Vector3d(R, R, R), bbB = mB + opB - Vector3d(R, R, R);
			Vector3d posmin = Vector3d(min(bbA(0), bbB(0)), min(bbA(1), bbB(1)), min(bbA(2), bbB(2)));
			Vector3d posmax = Vector3d(max(bBA(0), bBB(0)), max(bBA(1), bBB(1)), max(bBA(2), bBB(2)));
			return MatrixId(posmin, posmax);
		}
		else
		{
			// If the approximate footprint is not singular, then the bounding box is given by the three balls.
			double R = dB + rB, r = rB;
			Vector3d bBA = mB + opA + Vector3d(R, R, R), bbA = mB + opA - Vector3d(R, R, R);
			Vector3d bBB = mB + opB + Vector3d(R, R, R), bbB = mB + opB - Vector3d(R, R, R);
			Vector3d bBO = mB + Vector3d(r, r, r), bbO = mB - Vector3d(r, r, r);
			Vector3d posmin = Vector3d(minof(bbA(0), bbB(0), bbO(0)), minof(bbA(1), bbB(1), bbO(1)), minof(bbA(2), bbB(2), bbO(2)));
			Vector3d posmax = Vector3d(maxof(bBA(0), bBB(0), bBO(0)), maxof(bBA(1), bBB(1), bBO(1)), maxof(bBA(2), bBB(2), bBO(2)));
			return MatrixId(posmin, posmax);
		}
	}
	// Constructor by a box B in \intbox W.
	DeltaWtFp(MatrixId Bt, MatrixId Br, int wxyz=0)
	{
		// Bt part.
		mB = (Bt.min() + Bt.max()) / 2;
		rB = (Bt.max() - Bt.min()).norm() / 2;

		// Br part.
		vector<Vector4d> Corners = make_corners(Br, wxyz);
		vector<Vector3d> FpAs, FpBs;
		for (int i = 0; i < 8; ++i)
		{
			FpAs.push_back(FpA(Corners[i]));
			FpBs.push_back(FpB(Corners[i]));
		}
		opA = CenGra(FpAs);
		opB = CenGra(FpBs);
		dB = max(MaxDis(FpAs, opA), MaxDis(FpBs, opB));
		Point* OpA = new Point(opA + mB), * OpB = new Point(opB + mB), * OO = new Point(mB);

		if ((dB + rB) > opA.norm() || (dB + rB) > opB.norm())
		{
			singular = true;
			SegAB = new Edge(OpA, OpB);
			IccA = NULL;
			IccB = NULL;
			H = NULL;
			bbox = BBox();
		}
		else
		{
			singular = false;
			SegAB = new Edge(OpA, OpB);
			IccA = new IceCream(OO, OpA, dB);
			IccB = new IceCream(OO, OpB, dB);
			H = Hull();
			bbox = BBox();
		}
	}
	// If the WtFp does not intersect point feature f.
	bool free_from(Point* f)
	{
		if (bbox.ncontains(f->p))
			return true;
		if (singular)
			return SegAB->Sep(f, dB + rB);
		return SegAB->Sep(f, dB + rB) && H->Sep(f, 0) && IccA->Sep(f, rB) && IccB->Sep(f, rB);
	}
	// If the WtFp does not intersect edge feature f.
	bool free_from(Edge* f)
	{
		if (!bbox.intersects(B_Box(f)))
			return true;
		if (singular)
			return SegAB->Sep(f, dB + rB);
		return SegAB->Sep(f, dB + rB) && H->Sep(f, 0) && IccA->Sep(f, rB) && IccB->Sep(f, rB);
	}
	// If the WtFp does not intersect triangle feature f.
	bool free_from(Triangle* f)
	{
		if (!bbox.intersects(B_Box(f)))
			return true;
		if (singular)
			return SegAB->Sep(f, dB + rB);
		return SegAB->Sep(f, dB + rB) && H->Sep(f, 0) && IccA->Sep(f, rB) && IccB->Sep(f, rB);
	}
	// Classify the relation between WtFp with a point feature.
	int classify(Point* f)
	{
		if (free_from(f))
			return FREE;
		else
			return MIXED;
	}
	// Classify the relation between WtFp with a point feature.
	int classify(Edge* f)
	{
		if (free_from(f))
			return FREE;
		else
			return MIXED;
	}
	// Classify the relation between WtFp with a point feature.
	int classify(Triangle* f)
	{
		if (free_from(f))
			return FREE;
		else
			return MIXED;
	}
	// Classify the relation between WtFp with a mesh.
	int classify(Mesh* f)
	{
		if (!bbox.intersects(f->BBox()))
			return FREE;
		for (int i = 0; i < f->corners.size(); ++i)
			if (!free_from(f->corners[i]))
				return MIXED;
		for (int i = 0; i < f->edges.size(); ++i)
			if (!free_from(f->edges[i]))
				return MIXED;
		for (int i = 0; i < f->faces.size(); ++i)
			if (!free_from(f->faces[i]))
				return MIXED;
		if (f->inside(mB))
			return STUCK;
		else
			return FREE;
	}
};

#endif