#include <iostream>
#include <fstream>
#include <algorithm>
#include <math.h>
#include <float.h>
#include <sstream>
#include "term.h"
#include "protein.h"
#include "test_data.h"
#include "train_data.h"

Test_data::Test_data(const char *file_name, bool bSelectAllData, const struct PSLDoc_PARAMETER *param, int index, const int *nr_fold):Data(file_name, param)
{
	int i;
	vector<Protein>::iterator iter;
		
	for(iter = vec.begin(), i = 0; iter != vec.end(); iter++, i++) 
	{
		if(bSelectAllData ||((i%*nr_fold) == index)) //bSelectAllData: select all
			iter->construct_fea_vec(data_path.c_str());
		else
			vec.erase(iter--);
	}
	class_num = 0;
}

Test_data::Test_data(const char *file_name, const struct PSLDoc_PARAMETER *param):Data(file_name, param)
{
	int i;
	vector<Protein>::iterator iter;
		
	for(iter = vec.begin(), i = 0; iter != vec.end(); iter++, i++) 
		iter->construct_fea_vec(data_path.c_str());
	class_num = 0;
}

float CosineSimilarity(const vector<float> &testVec, const vector<float> &trainVec)
{
	float sum, testSum, trainSum;
	sum = testSum = trainSum = 0;
	
	if((int)testVec.size() != (int)trainVec.size())
	{
		cout << "ERROR in C osineSimilarity size not match " << (int)testVec.size() << " " <<  (int)trainVec.size() << endl;
		exit(0);
	}
	
	for(int i = 0; i < (int)testVec.size(); i++)
	{
		sum += (testVec[i]*trainVec[i]);
		trainSum += (trainVec[i]*trainVec[i]);
		testSum += (testVec[i]*testVec[i]);				
	}
	return sum/(float)(sqrt(testSum)*sqrt(trainSum));
}

//For kNN
struct sim_pair
{ 
	int index;
	float sim;
};

bool cmp(const sim_pair& a, const sim_pair& b)
{
    return a.sim > b.sim;
}

void Test_data::pred_kNN(Train_data train_data, int kNN)
{
	float sim, tmpSim;
	int index;
	vector<sim_pair> sim_vec;
	sim_pair tmp_sim_pair;

	if(kNN > (int)train_data.size())
		kNN = (int)train_data.size();
		
	for(int i = 0; i < (int)vec.size(); i++)
	{
		sim = tmpSim = DBL_MIN; index = 0;
		sim_vec.clear();
		for(int j = 0; j < (int)train_data.size(); j++)
		{
			tmp_sim_pair.index = j;
			tmp_sim_pair.sim = CosineSimilarity(vec[i].get_fea_vec(), train_data[j].get_fea_vec());
			sim_vec.push_back(tmp_sim_pair);
		}
		
		sort(sim_vec.begin(), sim_vec.end(), cmp);
		
		for(int k = 0; k < kNN; k++)
			vec[i].set_pred(train_data[sim_vec[k].index]);
	}
}

void Test_data::show_pred(const char *file_name) const
{
	int correct_count = 0;
	fstream file;
	file.open(file_name, ios::out);
	if(!file)
	{
		perror(file_name);
		exit(0);
	}
	
	file << "ProteinName,Original,Predicted,Correct,";
	for(int i = 0 ; i < class_num; i++)
		file << "prob_class-" << i+1 << ",";
	file << endl;

	for(int i = 0; i < (int)vec.size(); i++)
	{
		correct_count += vec[i].pred_correct();
		file << vec[i].get_name() << ",";
		vec[i].show_site(file);
		file << ",";
		vec[i].show_pred(file);
		file << "," << vec[i].pred_correct() << ",";
		vec[i].show_probability(file);
		file << endl;
	}
	file << "#accuracy=" << 100*correct_count/(float)vec.size() << endl;
	file.close();
}

PRE_ELEMENT str2pre_ele(string str_input)
{
	istringstream tmp_sstr(str_input);
	float tmp_flt;
	PRE_ELEMENT tmp_pre;

	tmp_sstr >> tmp_pre.pre_label;
	while (tmp_sstr >> tmp_flt)
	{
		tmp_pre.confidence.push_back(tmp_flt);
	}
	
	return tmp_pre;
}

int count_label(string label_str)
{
	istringstream tmp_sstr(label_str);
	int count = 0;
	string tmp_str;
	while (tmp_sstr >> tmp_str)
		count++;
	return --count;
}

void Test_data::parse_svm_result(const char * filename)
{
	fstream file;
	vector<string> file_content;
	string tmp_string;
	PRE_ELEMENT tmp_element;
	int index = 0;

	file.open(filename, ios::in);
	if(!file)
	{
		perror(filename);
		exit(0);
	}
	while(getline(file, tmp_string))
	{	
		if(tmp_string.find("labels") == string::npos)
			file_content.push_back(tmp_string);
		else	
			class_num = count_label(tmp_string);
	}	
	file.close();

	for(int i = 0; i < (int)file_content.size(); i++)
	{
		tmp_element = str2pre_ele(file_content[i]);
		if(index >= (int)vec.size())
		{
			cout << "[ERROR] # of predictions in SVM result = " << index 
			<< " > # of test set = " << (int)vec.size() << endl;
			exit(0);
		}
		vec[index].push_pre(tmp_element);
		index++;
	}
}

/*
void Test_data::parse_svm_result(const char * filename)
{
	fstream file;
	string pre_site;
	int index = 0;

	cout << "Parse the svm result file = " << filename << endl;
	file.open(filename, ios::in);
	if(!file)
	{
		perror(filename);
		exit(0);
	}
	while(file >> pre_site)
	{	
		if(index >= (int)vec.size())
		{
			cout << "[ERROR] # of predictions in SVM result = " << index 
			<< " > # of test set = " << (int)vec.size() << endl;
			exit(0);
		}
		vec[index].push_pre(pre_site);
		index++;			
	}
}
*/
//void Test_data::Evaluation()
//{
//	
//}
