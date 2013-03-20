#include <iostream>
#include <string>
#include <fstream>
#include <stdlib.h>

using namespace std;

#define DEBUG 0

void show_version()
{
	string version = "2010-01-28";
	cout << "Version: PSLDoc " << version << endl;
}

int check_pg_exist(string pg_name)
{
	string check_pg="which " + pg_name + " >/dev/null 2>/dev/null";
	int i_return = system(check_pg.c_str());
	if(i_return != 0)
	{
		cout << "(ERROR) program \"" << pg_name << "\" does not exist" << endl;
		return 0;
	}
	return 1;
}

bool fexists(string filename)
{
  ifstream ifile(filename.c_str());
  return ifile;
}

void my_system(string cmd)
{
	if(DEBUG)
	{
		cout << "system cmd = " << cmd << endl;
	}

	if(system(cmd.c_str()) != 0)
	{
		cout << "ERROR in " << cmd << endl;
		exit(1);
	}
}

bool plsaMatrixExist(string data_path)
{
	string wt_file = data_path+".wt";
	string td_file = data_path+".td";

	return (fexists(wt_file)&&fexists(td_file));
}

