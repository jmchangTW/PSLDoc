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

void parse_command_line(int argc, char **argv, char *train_file_name, char *model_file_name);
void show_param(char *train_file_name, char *model_file_name);
void initial();
void perform_svm(Test_data *testVec, const char *train_file_name,const char *test_file_name);

struct PSLDoc_PARAMETER param; // set by parse_command_line
bool per_PLSA;
bool SVM_OUTPUT_PRO;
int ner_num;
string str_signature_file;
Term gapped_dipeptide_term;

void exit_with_help()
{
	show_version();
	printf(
		"Usage: PSLDoc-test [options] model_file test_data_file\n"
		"options:\n"
		"\t-r reduction : feature reduction or not (default 1)\n"
		"\t\t0 -- Do not perform feature reduction\n"
		"\t\t1 -- Perform feature reduction by PLSA\n"
		"\t-p predict_method : set prediction method (default 1)\n"
		"\t\t0 -- k Nearest Neighbor\n"
		"\t\t1 -- SVM\n"
		"\t-s signature : set signature file (NULL)\n"
		"\t-k # of neighgors : set the number of Nearest Neighbors (default 1)\n"		
		"\t-d distance : set gapped-dipeptide distance (default 13)\n"
		"\t-f size : set the reduced feature size of PLSA (default 80)\n"
		"\t-i iteration : set the iteration of PLSA (default 300)\n"
		"\t-b probability : whether to predict probability estimates inside SVM (default 1)\n"
		"\t\t0 -- Do not output\n"
		"\t\t1 -- Output\n"
		"\nFor instance:\n"
		"PSLDoc-test data/simple_train data/simple_test\n"
		"Output:\n"
		"\tSVM input = data/simple_test.svm_input\n"
		"\tprediction result = data/simple_test.predict\n"
		"\tprediction result (csv format) = data/simple_test.csv\n"
	);
	exit(0);
}

int main(int argc, char *argv[])
{
	char train_file_name[1024];
	char test_file_name[1024];

	parse_command_line(argc, argv, train_file_name, test_file_name);	

cout << endl << "[1.INITIAL]" << endl;
	show_param(train_file_name, test_file_name);	
	initial();		

cout << endl << "[2.BEGIN]" << endl;
	string train_wt_file = (string)train_file_name + ".wt";
	string csv_result_file = (string)test_file_name + ".csv";

	Test_data testVec(test_file_name, &param);
	Train_data trainVec(train_file_name, &param);
	cout << "# of proteins in the test data = " << testVec.size() << endl;

	//PLSA
	if(per_PLSA)
	{
		cout << "(START) PLSA-foldIn..." << endl;
		PLSA plsa_test(&testVec, train_wt_file.c_str(), &param, gapped_dipeptide_term.size(), testVec.size());
		plsa_test.EM(false);
		plsa_test.write_matrix((string)test_file_name, "binary");
		testVec.update_vec(plsa_test.td);
	}
	//K-NN
	if(param.pred_type == 0)
	{
		cout << "(START) Nearest neighboring prediction : " << ner_num << "-NN" << endl;
		testVec.pred_kNN(trainVec, ner_num);
	}
	//SVM
	else if(param.pred_type == 1)
		perform_svm(&testVec, train_file_name, test_file_name);
	
	cout << "PREDICTION RESULT (csv format) = " << csv_result_file << endl; 
	testVec.show_pred(csv_result_file.c_str());

cout << endl << "[3.END]" << endl;
	return 0;
}

void parse_command_line(int argc, char **argv, char *train_file_name, char *test_file_name)
{
	int i;
	// default values
	param.data_path = (getenv ("PSLDOC_PATH")!= NULL?string(getenv ("PSLDOC_PATH")): DATA_PATH);
	param.gapped_distance = 13;
	param.feature_size = 80;
	param.iteration = 300; //EM interation
	param.pred_type = 1;
	
	per_PLSA = SVM_OUTPUT_PRO = true;
	ner_num = 1;
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
			case 'k':
				ner_num = atoi(argv[i]);
				break;
			case 'p':
				param.pred_type = atoi(argv[i]);
				break;				
			case 'r':
				per_PLSA = (atoi(argv[i])==1);
				break;
			case 'b':
				SVM_OUTPUT_PRO = (atoi(argv[i])==1);
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
	if((i>=argc)||(i+1>=argc))
		exit_with_help();

	strcpy(train_file_name, argv[i]);
	strcpy(test_file_name, argv[i+1]);
}

void show_param(char *train_file_name, char *test_file_name)
{

	cout << "Model file =\t" << train_file_name << endl
	<< "Test file =\t" << test_file_name << endl
	<< "What to do =\tperform PSLDoc prediction for test file based on model file" << endl;

	if((int)str_signature_file.length() > 0)
		cout << "\t\tGapped-Dipeptide Initialized by :\t" << str_signature_file << endl;

	if(param.pred_type == 0)
		cout	<< "\t\t# Nearest Neighbors\t" << ner_num << endl;
	cout << endl;

	cout << "----------------------PSLDoc PARAMETERS----------------------" << endl
		<< "data_path =\t" << param.data_path << endl
		<< "gapped_distance =\t" << param.gapped_distance << endl
		<< "PLSA_feature_size =\t" << param.feature_size << endl  	//PLSA parameter
		<< "PLSA_iteration =\t" << param.iteration << endl
		<< "pred_type =\t" << (param.pred_type==0?"Nearest Neighbors":"libsvm") << endl;
	cout << "-------------------------------------------------------------" << endl;
}

void initial()
{
	if((int)str_signature_file.length() > 0)
		gapped_dipeptide_term.initial(str_signature_file);
	else
		gapped_dipeptide_term.initial(param.gapped_distance);
}

void perform_svm(Test_data *testVec, const char *train_file_name,const char *test_file_name)
{
	string cmd;
	string SVM_test_file = (string)test_file_name + ".svm_input";
	string SVM_model_file = (string)train_file_name + ".svm_model";
	string SVM_range_file = (string)train_file_name + ".svm_range";
	string SVM_result_file = (string)test_file_name + ".predict";
	string SVM_output_prob_cmd = " -b 1 ";

	//Output Reduced PLSA feature
	if(!fexists(SVM_test_file))
		testVec->gen_svm_file(SVM_test_file.c_str());
	else
		cout << "(NOTICE) SVM will use the existing svm file = " << SVM_test_file << endl;

	cout << "(START) libsvm..." << endl;		

//scale svm input file according to original range file
	cmd  = (string)SVM_SCALE_CMD + " -r " + SVM_range_file + " " + SVM_test_file + " > tmp.txt";
	cout << "scaling input : " << cmd << endl;
	my_system(cmd);

//perform SVM predict
	if(SVM_OUTPUT_PRO)
		cmd  = (string)SVM_PREDICT_CMD + SVM_output_prob_cmd + " tmp.txt " + SVM_model_file + " " + SVM_result_file + " > /dev/null";
	else
		cmd  = (string)SVM_PREDICT_CMD + " tmp.txt " + SVM_model_file + " " + SVM_result_file + " > /dev/null";
	cout << "predicting : " << cmd << endl;
	my_system(cmd);

	testVec->parse_svm_result(SVM_result_file.c_str());
//	cmd = "rm tmp.txt; rm " + SVM_result_file;
//	my_system(cmd);
}
