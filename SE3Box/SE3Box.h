// SE3Box.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include <vector>
#include <bimap.h>
#include "SE3Box.cpp"
using namespace Eigen;
using namespace std;

class R3Box;
class SO3Box;
class SE3Box;

// Binary string for the tree code, 0 for \bar{1}, 1 for 1. We add a "1" at the beginning to avoid starting with "0".
class bit;
// Output a bit.
ostream& operator<<(ostream& os, bit t);

// Get the t-th binary bit of n.
bool getBin(int n, int t);