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
#include "test_data.h"
#include "plsa.h"
#include "util.h"

using namespace std;

void parse_command_line(int argc, char **argv, char *input_file_name, string& SVM_train_file);
void show_param(char *input_file_name);
void initial();
void move_svm_output(string ori_name, string des_name);

struct PSLDoc_PARAMETER param; // set by parse_command_line
bool per_PLSA, use_exist_data;
string str_signature_file;
Term gapped_dipeptide_term;

void exit_with_help()
{
	printf(
		"Usage: PSLDoc-train [options] train_file\n"
		"options:\n"
		"\t-r reduction : feature reduction or not (default 1)\n"
		"\t\t0 -- Do not perform feature reduction\n"
		"\t\t1 -- Perform feature reduction by PLSA\n"
		"\t-s signature : set signature file (NULL)\n"
		"\t-d distance : set gapped-dipeptide distance (default 13)\n"
		"\t-f size : set the reduced feature size of PLSA (default 80)\n"
		"\t-i iteration : set the iteration of PLSA (default 300)\n"	
		"\t-e exist_data : use existing data (default 1)\n"
		"\t\t0 -- No\n"
		"\t\t1 -- Yes\n"
		"\nFor instance:\n"
		"PSLDoc-train data/simple_train\n"
		"Output:\n"
		"\tPLSA wt matrix (binary format) = data/simple_train.wt\n"
		"\tPLSA td matrix (binary format) = data/simple_train.td\n"
		"\tPLSA wt matrix (csv format) = data/simple_train_wt.csv\n"
		"\tPLSA td matrix (csv format) = data/simple_train_td.csv\n"
		"\tSVM input = data/simple_train.svm_input\n"
		"\tSVM model = data/simple_train.svm_model\n"
		"\tSVM range = data/simple_train.svm_range\n"
	);
	exit(0);
}

int main(int argc, char *argv[])
{
	char input_file_name[1024];
	string SVM_train_file;

	parse_command_line(argc, argv, input_file_name, SVM_train_file);

cout << endl << "[1.INITIAL]" << endl;
	show_param(input_file_name);
	initial();

cout << endl << "[2.BEGIN]" << endl;
	Train_data trainVec(input_file_name, &param);
	cout << "# of proteins in train data = " << trainVec.size() << endl;

	if(per_PLSA)
	{
		if(use_exist_data && plsaMatrixExist((string)input_file_name))
		{
			string wt_file = (string)input_file_name+".wt";
			string td_file = (string)input_file_name+".td";
			cout << "(NOTICE) PLSA will use existing files" << endl
			<< "wt_file = " << wt_file << endl
			<< "td_file = " << td_file << endl;

			PLSA plsa_train(&trainVec, wt_file.c_str(), td_file.c_str(), &param, gapped_dipeptide_term.size(), trainVec.size());
			trainVec.update_vec(plsa_train.td);
		}
		else
		{
			cout << "(START)PLSA-train..." << endl;
			PLSA plsa_train(&trainVec, &param, gapped_dipeptide_term.size(), trainVec.size());
			plsa_train.EM(true);
			plsa_train.write_matrix((string)input_file_name, "binary");
			plsa_train.write_matrix((string)input_file_name, "csv");
			trainVec.update_vec(plsa_train.td);
		}
	}

//	Output Reduced PLSA feature
	if(!fexists(SVM_train_file))
		trainVec.gen_svm_file(SVM_train_file.c_str());
	else
		cout << "(NOTICE) SVM will use existing svm_train file = " << SVM_train_file << endl;
	
	cout << "(START) libsvm..." << endl;
	check_pg_exist(LIBSVM_EASYPY_CMD);
	string cmd = (string)LIBSVM_EASYPY_CMD + " " + (string)SVM_train_file;
	my_system(cmd);
	move_svm_output(SVM_train_file, (string)input_file_name);

cout << endl << "[3.END]" << endl;
	return 0;
}

void parse_command_line(int argc, char **argv, char *input_file_name, string& SVM_train_file)
{
	int i;

	// default values
	param.data_path = (getenv ("PSLDOC_PATH")!= NULL?string(getenv ("PSLDOC_PATH")): DATA_PATH);
	param.gapped_distance = 13;
	param.feature_size = 80;
	param.iteration = 300; //EM interation
	
	per_PLSA = use_exist_data = true;
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
			case 'i':
				param.iteration = atoi(argv[i]);
				break;
			case 'r':
				per_PLSA = (atoi(argv[i])==1);
				break;
			case 'e':
				use_exist_data = (atoi(argv[i])==1);
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
	SVM_train_file = (string)input_file_name +".svm_input";
}

void show_param(char *input_file_name)
{
	cout << "Train file =\t" << input_file_name << endl
	<< "What to do =\tperform PSLDoc training based on the train file" << endl << endl;

	if((int)str_signature_file.length() > 0)
		cout << "\tGapped-Dipeptide Initialized by :\t" << str_signature_file << endl;

	cout << "----------------------PSLDoc PARAMETERS----------------------" << endl
		<< "data_path =\t" << param.data_path << endl
		<< "gapped_distance =\t" << param.gapped_distance << endl
		<< "PLSA_feature_size =\t" << param.feature_size << endl  	//PLSA parameter
		<< "PLSA_iteration =\t" << param.iteration << endl;
	cout << "-------------------------------------------------------------" << endl;
}

void initial()
{
	if((int)str_signature_file.length() > 0)
		gapped_dipeptide_term.initial(str_signature_file);
	else
		gapped_dipeptide_term.initial(param.gapped_distance);
}

void move_svm_output(string ori_name, string des_name)
{
	string ori_range_name = ori_name.substr(ori_name.find_last_of("/")+1) + ".range";
	string ori_model_name = ori_name.substr(ori_name.find_last_of("/")+1) + ".model";
	string des_range_name = des_name + ".svm_range";
	string des_model_name = des_name + ".svm_model";
	string cmd;

	if(fexists(ori_range_name))
	{
		cmd = "mv " + ori_range_name + " " + des_range_name;
		my_system(cmd);
	}
	else
	{
		cout << "(ERROR) svm range file = " << ori_range_name << " does not exist";
		exit(1);
	}

	if(fexists(ori_model_name))
	{
		cmd = "mv " + ori_model_name + " " + des_model_name;
		my_system(cmd);
	}
	else
	{
		cout << "(ERROR) svm model file = " << ori_model_name << " does not exist";
		exit(1);
	}

	cmd = "rm " + ori_name.substr(ori_name.find_last_of("/")+1) + ".scale*";
	my_system(cmd);
}
