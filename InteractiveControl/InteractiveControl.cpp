// InteractiveControl.cpp: This file controls the interactive of SSS main program.

#ifndef INTERACTIVECONTROL_H
#define INTERACTIVECONTROL_H

#include<iostream>
#include<string>
#include<Eigen/Dense>
#include<FeatureSet.h>
#include<ReadSSSCommand.h>
#include<IntroControl.h>
#include<RunControl.h>
#include<FileCheck.h>

string extern fileplace;
string extern fileformat;

EnvironmentFeature extern env;
VectorXd extern SSSalpha, SSSbeta;
double extern r0;
double extern envrange;
double extern varepsilon;
int extern ExpandLimit;
heutype extern SSSheu;
SSSViewer extern viewer;

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
			if (filename.find("/") == string::npos && filename.find("\\") == string::npos && filename.find(".") == string::npos)
				filename = fileplace + filename + fileformat;
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
			SSS_show();
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
		case LIMIT:
		{
			if (words.size() < 3)
			{
				cout << "bash: not enough argument" << endl;
				return true;
			}
			ExpandLimit = stoi(words[2]);
			cout << "Set the maximum expand limit to be " << ExpandLimit << "." << endl;
			return true;
		}
		case EPSILON:
		{
			if (words.size() < 3)
			{
				cout << "bash: not enough argument" << endl;
				return true;
			}
			varepsilon = stod(words[2]);
			cout << "Set epsilon to be " << varepsilon << "." << endl;
			return true;
		}
		default:cout << "bash: " << words[1] << ": command not found" << endl; return true;
		}
	}
	case RUN:
	{
		SSS_run();
		return true;
	}
	case EXIT:
		return false;
	default:cout << "bash: " << words[0] << ": command not found" << endl; return true;
	}
}

// Interactive mode.
void interactive()
{
	SSS_intro();
	cout << "Type \"exit\" to exit." << endl;
	char command[1000];
	while (true)
	{
		cout << ">> ";
		cin.getline(command, 1000);
		if (!decode(string(command)))
			break;
	}
	cout << "Find path algorithm terminates." << endl;
}

#endif