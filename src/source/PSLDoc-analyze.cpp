//Command
//mpirun -machinefile hosts_file -np #_of_Nodes PSLDocOO inputfile
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <fcntl.h>
#include <dirent.h>
#include <string.h>

#include "config.h"
#include "PSLDoc.h"
#include "term.h"
#include "data.h"
#include "train_data.h"
#include "plsa.h"
#include "util.h"

using namespace std;

void parse_command_line(int argc, char **argv, char *input_file_name);
void show_param(char *input_file_name);
void initial();

struct PSLDoc_PARAMETER param; // set by parse_command_line
string str_signature_file;
Term gapped_dipeptide_term;

void exit_with_help()
{
	printf(
		"Usage: PSLDoc-analyze [options] data_set_file\n"
		"options:\n"
		"\t-s signature : set signature file (NULL)\n"
		"\t-d distance : set gapped-dipeptide distance (default 13)\n"
		"\nFor instance\n"
		"PSLDoc-analyze data/simple_train\n"
		"Output:\n"
		"\tTopic vs localization class = tc.csv\n"
		"\tAnalyze result = data/simple_train.analysis\n"
	);
	exit(0);
}

int main(int argc, char *argv[])
{
	char input_file_name[1024];

	parse_command_line(argc, argv, input_file_name);

cout << endl << "[1.INITIAL]" << endl;	
	show_param(input_file_name);
	initial();

cout << endl << "[2.BEGIN]" << endl;
	cout << "(STRAT) Preform PLSA Analysis..." << endl;
	Train_data inputVec(input_file_name, &param);
	vector<string> vClass;
	inputVec.getClass(vClass);

	cout << "# of proteins =\t" << inputVec.size() << endl
	<< "# of localization class =\t" << (int)vClass.size() << endl;

	if(plsaMatrixExist((string)input_file_name))
	{
		string wt_file = (string)input_file_name+".wt";
		string td_file = (string)input_file_name+".td";
		PLSA plsa_train(&inputVec, wt_file.c_str(), td_file.c_str(), &param, gapped_dipeptide_term.size(), inputVec.size());
		plsa_train.analyze((string)input_file_name+".analysis", &inputVec, vClass);
	}
	else
	{
		cout << "[ERROR] There is no wt and td binary files" << endl
		<< "PLSA analyzation is based on the wt and td binary files." << endl;
		exit(1);
	}

cout << endl << "[3.END]" << endl;
	return 1;
}

void parse_command_line(int argc, char **argv, char *input_file_name)
{
	int i;
	// default values
	param.data_path = (getenv ("PSLDOC_PATH")!= NULL?string(getenv ("PSLDOC_PATH")): DATA_PATH);
	param.gapped_distance = 13;
	param.feature_size = 80;

	// parse options
	for(i=1;i<argc;i++)
	{
		if(argv[i][0] != '-') break;
		if(++i>=argc)
			exit_with_help();
		switch(argv[i-1][1])
		{
			case 'd':
				param.gapped_distance = atoi(argv[i]);
				break;
			case 'f':
				param.feature_size = atoi(argv[i]);
				break;
			case 's':
				str_signature_file = argv[i];
				break;
			default:
				fprintf(stderr,"unknown option\n");
				exit_with_help();
		}
	}
	// determine filenames
	if(i>=argc)
		exit_with_help();

	strcpy(input_file_name, argv[i]);	
}

void show_param(char *input_file_name)
{
	cout << "Analyzed file =\t" << input_file_name << endl
	<< "What to do =\tPerform PLSA analysis" << endl;

	if((int)str_signature_file.length() > 0)
		cout << "\tGapped-Dipeptide Initialized by :\t" << str_signature_file << endl;
	cout << endl;

	cout << "----------------------PSLDoc PARAMETERS----------------------" << endl
		<< "data_path=\t" << param.data_path << endl
		<< "gapped_distance=\t" << param.gapped_distance << endl;
	cout << "-------------------------------------------------------------" << endl;
}

void initial()
{
	if((int)str_signature_file.length() > 0)
		gapped_dipeptide_term.initial(str_signature_file);
	else
		gapped_dipeptide_term.initial(param.gapped_distance);
}
