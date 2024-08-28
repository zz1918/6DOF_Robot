// FileCheck.cpp: This file checks if a file address is valid or not.

#ifndef FILECHECK_H
#define FILECHECK_H

#include<string>
#include<sys/stat.h>

bool exist(const std::string& filename)
{
	struct stat buffer;
	return (stat(filename.c_str(), &buffer) == 0);
}

#endif