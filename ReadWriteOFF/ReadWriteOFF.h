// ReadWriteOFF.h

#include "ReadWriteOFF.cpp"

void read_OFF(const string& filename, MatrixXd& V, MatrixXi& F);
void write_OFF(const string& filename, MatrixXd& V, MatrixXi& F);