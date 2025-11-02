// BinStr.h: This is the head file for BinStr.cpp.
//

#include"BinStr.cpp"

// Binary string for the tree code, 0 for \bar{1}, 1 for 1. We add a "1" at the beginning to avoid starting with "0".
class bit;

// Output a bit.
ostream& operator<<(ostream& os, bit t);

// Get the t-th binary bit of n.
bool getBin(int n, int t);