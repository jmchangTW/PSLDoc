//Command
//mpirun -machinefile hosts_file -np #_of_Nodes PSLDocOO inputfile
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <fcntl.h>
#include <dirent.h>
#include <string.h>
#include <omp.h>

#include "config.h"
#include "PSLDoc.h"
#include "term.h"
#include "data.h"
#include "util.h"

using namespace std;

void parse_command_line(int argc, char **argv, char *input_file_name, char *model_file_name);
void show_param(char *input_file_name, char *model_file_name);
void initial();

struct PSLDoc_PARAMETER param; // set by parse_command_line
int nr_fold, ner_num;
string str_signature_file;
Term gapped_dipeptide_term;

void exit_with_help()
{
	printf(
		"Usage: PSLDoc-prepare [options] data_set_file\n"
		"options:\n"
		"\t-o overwrite_TFPSSM \toverwrite TFPSSM data (default 0)\n"
		"\t\t0 -- No\n"
		"\t\t1 -- Yes\n"
		"\t-s signature \tset signature file (NULL)\n"
		"\t-d distance \tset gapped-dipeptide distance (default 13)\n"
		"\t-j loop \tset the number of the Loop in PSIBLAST (default 2)\n"		
		"\t-e value \tset the e-value of PSIBLAST (default 0.01)\n"
		"\nFor instance:\n"
		"PSLDoc-prepare -e 0.001 data/simple_train\n"
		"Output: *.pssm, *.tfpssm in the path specified by param.data_path\n"
	);
	exit(0);
}

void test(int index)
{
	cout << "process index = " << index << " in " << omp_get_thread_num() << endl;
}

int main(int argc, char *argv[])
{
	char input_file_name[1024];
	char test_file_name[1024];
	int run_index = (getenv("RUN_INDEX")!= NULL?atoi(getenv("RUN_INDEX")):-1);
	
	parse_command_line(argc, argv, input_file_name, test_file_name);
	cout << endl << "[1.INITIAL]" << endl;
	show_param(input_file_name, test_file_name);

	initial();
	cout << endl << "[2.BEGIN]" << endl;
	Data all(input_file_name, &param);
	
	cout << "protein size = " << all.size() << endl;
// 	#pragma omp parallel for
	for(int i = 0; i< all.size(); i++)
	{
 		if( (run_index == -1)||(i == run_index) )
		{
		   all.protein2feature(i, &param, gapped_dipeptide_term);
		}
	}

	cout << endl << "[3.END]" << endl;
	return 0;
}

void parse_command_line(int argc, char **argv, char *input_file_name, char *test_file_name)
{
	int i;
	// default values
	param.data_path = (getenv ("PSLDOC_PATH")!= NULL?string(getenv ("PSLDOC_PATH")): DATA_PATH);
	param.overwrite_tfpssm = 0;
	param.gapped_distance = 13;
	param.PSIBLAST_iteration = "2";
	param.PSIBLAST_e_value = "0.001";
	
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
			case 'j':
				param.PSIBLAST_iteration = argv[i];
				break;
			case 'e':
				param.PSIBLAST_e_value = argv[i];
				break;
			case 'o':
				param.overwrite_tfpssm = atoi(argv[i]);
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

void show_param(char *input_file_name, char *test_file_name)
{
	cout << "Input File=\t" << input_file_name << endl;	
	cout << "What to do=\tFEATURIZED : generate PSSM and TFPSSM" << endl;
	cout << endl;

	if((int)str_signature_file.length() > 0)
		cout << "\tGapped-Dipeptide Initialized by File :\t" << str_signature_file << endl;
	
	cout << "----------------------PSLDoc PARAMETERS----------------------" << endl
		<< "overwrite_tfpssm\t\t" << (param.overwrite_tfpssm?"YES":"NO") << endl
		<< "PSLDOC_PATH\t\t" << param.data_path << endl
		<< "gapped_distance\t\t" << param.gapped_distance << endl
		<< "PSIBLAST_iteration\t\t" << param.PSIBLAST_iteration << endl //Set PSIBLAST parameter
		<< "PSIBLAST_e_value\t\t" << param.PSIBLAST_e_value << endl;
	cout << "-------------------------------------------------------------" << endl;
}

void initial()
{
	if((int)str_signature_file.length() > 0)
		gapped_dipeptide_term.initial(str_signature_file);
	else
		gapped_dipeptide_term.initial(param.gapped_distance);
}
