// FeatureSet.cpp: This file defines the set of features and give methods of adding features.

#ifndef FEATURESET_H
#define FEATURESET_H

#include<iostream>
#include<sstream> 
#include<string>
#include<vector>
#include<map>
#include<ReadWriteOFF.h>
#include<Eigen/Dense>
#include<WtFp.h>
using namespace std;
using namespace Eigen;

class DeltaFeature
{
public:
	vector<Point*> Vlist;
	vector<Edge*> Elist;
	vector<Triangle*> Tlist;
	vector<Mesh*> Mlist;
	bool empty()
	{
		return Vlist.empty() && Elist.empty() && Tlist.empty() && Mlist.empty();
	}
	void clear()
	{
		Vlist.clear();
		Elist.clear();
		Tlist.clear();
		Mlist.clear();
	}
	void out()
	{
		if (empty())
		{
			cout << "The feature set is empty." << endl;
			return;
		}
		cout << "The feature set contains: " << endl;
		for (int i = 0; i < Vlist.size(); ++i)
			cout << "    Point (" << Vlist[i]->p.transpose() << ")" << endl;
		for (int i = 0; i < Elist.size(); ++i)
			cout << "    Edge (" << Elist[i]->P(0)->p.transpose() << ")-(" << Elist[i]->P(1)->p.transpose() << ")" << endl;
		for (int i = 0; i < Tlist.size(); ++i)
			cout << "    Triangle (" << Tlist[i]->P(0)->p.transpose() << ")-(" << Tlist[i]->P(1)->p.transpose() << ")-(" << Tlist[i]->P(2)->p.transpose() << ")" << endl;
	}
};

class EnvironmentFeature
{
public:
	map<string, Point*> Vlist;
	map<string, Edge*> Elist;
	map<string, Triangle*> Tlist;
	map<string, Mesh*> Mlist;
	Point* find_point(string name)
	{
		map<string, Point*>::iterator it = Vlist.find(name);
		if (it == Vlist.end())
			return NULL;
		else
			return it->second;
	}
	Edge* find_edge(string name)
	{
		map<string, Edge*>::iterator it = Elist.find(name);
		if (it == Elist.end())
			return NULL;
		else
			return it->second;
	}
	Triangle* find_face(string name)
	{
		 map<string, Triangle*>::iterator it = Tlist.find(name);
		 if (it == Tlist.end())
			 return NULL;
		 else
			 return it->second;
	}
	Mesh* find_mesh(string name)
	{
		map<string, Mesh*>::iterator it = Mlist.find(name);
		if (it == Mlist.end())
			return NULL;
		else
			return it->second;
	}
	// Construct a new point from a 3d point.
	Point* add_point(Vector3d p, string name)
	{
		Point* new_point = new Point(p);
		Vlist.insert(make_pair(name, new_point));
		return new_point;
	}
	// Construct a new edge from two 3d points.
	Edge* add_edge(Vector3d p, Vector3d q, string name)
	{
		Point* P = add_point(p, name + "_P");
		Point* Q = add_point(q, name + "_Q");
		Edge* E = new Edge(P, Q);
		Elist.insert(make_pair(name, E));
		return E;
	}
	// Construct a new edge from a solid edge.
	Edge* add_edge(Edge* E, string name)
	{
		if (E == NULL)
			return NULL;
		Elist.insert(make_pair(name, E));
		return E;
	}
	// Construct a new edge from two points with names.
	Edge* add_edge(string name1, string name2, string name)
	{
		Point* P = find_point(name1);
		if (P == NULL)
			return NULL;
		Point* Q = find_point(name2);
		if (Q == NULL)
			return NULL;
		Edge* E = new Edge(P, Q);
		Elist.insert(make_pair(name, E));
		return E;
	}
	// Construct a new face from three 3d points.
	Triangle* add_face(Vector3d p, Vector3d q, Vector3d r, string name)
	{
		Point* P = add_point(p, name + "_P");
		Point* Q = add_point(q, name + "_Q");
		Point* R = add_point(r, name + "_R");
		Edge* E = new Edge(P, Q);
		Edge* F = new Edge(Q, R);
		Edge* G = new Edge(R, P);
		Elist.insert(make_pair(name + "_E", E));
		Elist.insert(make_pair(name + "_F", F));
		Elist.insert(make_pair(name + "_G", G));
		Triangle* T = new Triangle(E, F, G);
		Tlist.insert(make_pair(name, T));
		return T;
	}
	// Construct a new face from three points with names.
	Triangle* add_face(string name1, string name2, string name3, string name)
	{
		Point* P = find_point(name1);
		if (P == NULL)
			return NULL;
		Point* Q = find_point(name2);
		if (Q == NULL)
			return NULL;
		Point* R = find_point(name3);
		if (R == NULL)
			return NULL;
		Edge* E = new Edge(P, Q);
		Edge* F = new Edge(Q, R);
		Edge* G = new Edge(R, P);
		Elist.insert(make_pair(name1 + "_" + name2, E));
		Elist.insert(make_pair(name2 + "_" + name3, F));
		Elist.insert(make_pair(name3 + "_" + name1, G));
		Triangle* T = new Triangle(E, F, G);
		Tlist.insert(make_pair(name, T));
		return T;
	}
	// Construct a new mesh from V-F matrix.
	Mesh* add_mesh(MatrixXd V, MatrixXi F, string name)
	{
		Mesh* new_mesh = new Mesh(V, F);
		/*for (int i = 0; i < new_mesh->corners.size(); ++i)
			Vlist.insert(make_pair(name + "_P_" + to_string(i), new_mesh->corners[i]));
		for (int i = 0; i < new_mesh->edges.size(); ++i)
			Elist.insert(make_pair(name + "_E_" + to_string(i), new_mesh->edges[i]));
		for (int i = 0; i < new_mesh->faces.size(); ++i)
			Tlist.insert(make_pair(name + "_F_" + to_string(i), new_mesh->faces[i]));*/
		Mlist.insert(make_pair(name, new_mesh));
		return new_mesh;
	}
	// Construct a new mesh from .OFF address.
	Mesh* add_mesh(string filename, string name)
	{
		MatrixXd V;
		MatrixXi F;
		read_OFF(filename, V, F);
		return add_mesh(V, F, name);
	}
	// Remove a point by name.
	int remove_point(string name)
	{
		return Vlist.erase(name);
	}
	// Remove an edge by name.
	int remove_edge(string name)
	{
		return Elist.erase(name);
	}
	// Remove a point by name.
	int remove_face(string name)
	{
		return Tlist.erase(name);
	}
	// Remove an edge by name.
	int remove_mesh(string name)
	{
		return Mlist.erase(name);
	}
	// Construct the delta feature from the environment.
	DeltaFeature* make_feature()
	{
		DeltaFeature* Phi = new DeltaFeature();
		Phi->clear();
		for (map<string, Point*>::iterator it = Vlist.begin(); it != Vlist.end(); it++)
			Phi->Vlist.push_back(it->second);
		for (map<string, Edge*>::iterator it = Elist.begin(); it != Elist.end(); it++)
			Phi->Elist.push_back(it->second);
		for (map<string, Triangle*>::iterator it = Tlist.begin(); it != Tlist.end(); it++)
			Phi->Tlist.push_back(it->second);
		for (map<string, Mesh*>::iterator it = Mlist.begin(); it != Mlist.end(); it++)
			Phi->Mlist.push_back(it->second);
		return Phi;
	}
	// Show Vlist.
	void show_point()
	{
		cout << "The environment contains:" << endl;
		for (map<string, Point*>::iterator it = Vlist.begin(); it != Vlist.end(); it++)
			cout << "point " << it->first << ": " << it->second->p.transpose() << endl;
	}
	// Show Elist.
	void show_edge()
	{
		cout << "The environment contains:" << endl;
		for (map<string, Edge*>::iterator it = Elist.begin(); it != Elist.end(); it++)
			cout << "edge " << it->first << endl;
	}
	// Show Tlist.
	void show_face()
	{
		cout << "The environment contains:" << endl;
		for (map<string, Triangle*>::iterator it = Tlist.begin(); it != Tlist.end(); it++)
			cout << "face " << it->first << endl;
	}
	// Show Mlist.
	void show_mesh()
	{
		cout << "The environment contains:" << endl;
		for (map<string, Mesh*>::iterator it = Mlist.begin(); it != Mlist.end(); it++)
			cout << "mesh " << it->first << endl;
	}
};

#endif