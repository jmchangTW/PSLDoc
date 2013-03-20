#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <algorithm>
#include <omp.h>

#include <stdio.h>
#include <string.h>
#include <errno.h>

#include "term.h"
#include "protein.h"
#include "plsa.h"
#include "data.h"
#include "PSLDoc.h"
#include "util.h"

#include <unistd.h>
#include <fcntl.h>
#include <cerrno> // for errno 

extern float CosineSimilarity(const vector<float> &testVec, const vector<float> &trainVec);

Data::Data()
{
}

Data::Data(const char *file_name, const struct PSLDoc_PARAMETER *param)
{
	fstream file;
	string str, name, seq, tmp;
	vector<string> site;

	file.open(file_name, ios_base::in);
	if(!file)
	{
		system("hostname");
		perror(file_name);
		exit(0);
	}
	while(!file.eof())
	{
		getline(file, str);
		if(str[0] == '>')
		{
		  if(seq.length() > 0)
		  {
		    vec.push_back(Protein(name, seq, site));
		    seq="";
		    site.clear();
		  }  
		  name = str.substr(1, str.find(" ")-1);
		  tmp = str.substr(str.find("=")+1, str.find(";")-str.find("=")-1);
		  site.push_back(tmp);
		}
		else
		{
		  seq += str;
		}
	}
	vec.push_back(Protein(name, seq, site));
	file.close();
	file.clear();

	data_path = param->data_path;
/*mk dir for self name
	data_path = param->data_path + (strrchr(file_name,'/')?strrchr(file_name,'/')+1:file_name) + "/";
	string str_make_dir = "mkdir -p " + data_path;
	if(system(str_make_dir.c_str()) != 0)
	{
		cout << "(ERROR) system cmd : " << str_make_dir << endl;
		exit(1);
	}
*/
}

Data::~Data()
{
}

Protein Data::operator[](int i) const
{
	if(i < (int)vec.size())
		return vec[i];
	else
		cout << "ERROR in acess vec " << i << " > size " << (int)vec.size() << endl;
	return vec[0];
}

void Data::update_vec(Matrix matrix)
{
	if(matrix.get_size2() != (int)vec.size())
	{
		cout << "Error: Size Not Match " << matrix.get_size2() << " " <<  (int)vec.size() << endl;
		exit(0);
	}
	
	gsl_vector *gsl_vec = gsl_vector_alloc(matrix.get_size1());	
	for(int i = 0; i < (int)matrix.get_size2(); i++)
	{
		gsl_matrix_get_col(gsl_vec, matrix.matrix, i);
		vec[i].construct_fea_vec(gsl_vec);
	}
	gsl_vector_free(gsl_vec);
}

void Data::gen_svm_file(const char* filename) const
{
	fstream file;

	cout << "Output svm file : " << filename << endl;
	file.open(filename, ios::out);
	if(!file)
	{
		perror(filename);
		exit(0);
	}
	for(int i = 0; i < (int)vec.size(); i++)
	{
		vec[i].show_site(file);
		vec[i].show_fea(file);
		file << endl;
	}
	file.close();
}
// int check_run(string fileName)
// {
//    int done = 0;
//    
//    int fd = open(fileName.c_str(), O_WRONLY|O_CREAT|O_EXCL, 0644);
//    if (fd != -1){     done = 1;   }
//    if (errno != EEXIST)
//    {
//       /* do some sort of error handling here... */
//       cout << "[WARNNING] in lock access " << fileName << endl
//       << "message : " << strerror(errno) << endl;
//    }
//    close(fd);
//    return done;
// }

void Data::protein2feature(int i, const struct PSLDoc_PARAMETER *param, Term term) const
{
// 	string run, cmd;
// 	run = vec[i].get_name() + ".run";
	
// 	if(i < (int)vec.size() && !fexists(run) && check_run(run))
	if(i < (int)vec.size())
	{
	  vec[i].gen_PSSM(param->data_path.c_str(), param);
	  vec[i].gen_TFPSSM(param->data_path.c_str(), term, param);
	}
}

void Data::gen_distance_matrix()
{
	fstream file;
	
	file.open("matrix.txt", ios::out);
	if(!file)
	{
		cout << "ERROR: matrix file can not be opened" << endl;
		exit(0);
	}	
	file << (int)vec.size() << endl;
	for(int i = 0; i < (int)vec.size()-1; i++)
	{
		for(int k = 0; k < i; k++)
			file << "      \t";
		for(int j = i+1; j < (int)vec.size(); j++)
			file << fixed << setprecision(4) 
			<< (1-CosineSimilarity(vec[i].get_fea_vec(), vec[j].get_fea_vec()))
			<< "\t";
		file << endl;
	}
	file.close();

	fstream list_file;	
	list_file.open("name_list.txt", ios::out);
	if(!list_file)
	{
		cout << "ERROR: matrix file can not be opened" << endl;
		exit(0);
	}	
	for(int i = 0; i < (int)vec.size(); i++)
			list_file << i+1 << " " << vec[i].get_name() << endl;
	list_file.close();	
}

void Data::getClass(vector<string> &vClass) const
{
	vector<string> vSiteTmp;

	for(int i = 0; i < (int)vec.size(); i++)
	{
		vSiteTmp.clear();
		vSiteTmp = vec[i].get_site_vec();
		for(int j = 0; j < (int)vSiteTmp.size(); j++)
		{
			if(find(vClass.begin(), vClass.end(), vSiteTmp[j]) == vClass.end())
				vClass.push_back(vSiteTmp[j]);
		}
	}
}

// bool check_with_label(const char *file_name)
// {
// 	fstream file;
// 	string str;
// 	int name_dis = 0;
// 
// 	file.open(file_name, ios_base::in);
// 	if(!file)
// 	{
// 		system("hostname");
// 		perror(file_name);
// 		exit(0);
// 	}
// 	getline(file, str);
// 	if (str.substr(0, 1) != ">")
// 	{
// 		cout << "[ERROR] input file = " << file_name << " is not fasta format" << endl;
// 		cout << str[0] << " != > " << endl;
// 		exit(1);
// 	}	
// 	
// 	while(!file.eof())
// 	{
// 		getline(file, str);
// 		if((int)str.length() > 0)
// 		{
// 			if(str.substr(0, 1) == ">")
// 			{
// 				break;	
// 			}
// 			else
// 			{
// 				name_dis++;
// 			}
// 		}
// 	}		
// 	file.close();
// 
// 	return (name_dis==2);
// }

// Data::Data(const char *file_name, const struct PSLDoc_PARAMETER *param)
// {
// 	fstream file;
// 	string str, name, seq;
// 	vector<string> site;		
// 	bool with_label = check_with_label(file_name);
// 
// 	file.open(file_name, ios_base::in);
// 	if(!file)
// 	{
// 		system("hostname");
// 		perror(file_name);
// 		exit(0);
// 	}
// 	while(!file.eof())
// 	{
// 		getline(file, str);
// 		if((int)str.length() > 0)
// 		{
// 			name = str.substr(1);
// 			if(with_label)
// 			{
// 				getline(file, str);
// 				site.clear();
// 				if(str.find(',') == string::npos)
// 					site.push_back(str);			
// 			}
// 			getline(file, seq);
// 			vec.push_back(Protein(name, seq, site));
// 		}
// 	}
// 	file.close();
// 	file.clear();
// 
// 	data_path = param->data_path;
// /*mk dir for self name
// 	data_path = param->data_path + (strrchr(file_name,'/')?strrchr(file_name,'/')+1:file_name) + "/";
// 	string str_make_dir = "mkdir -p " + data_path;
// 	if(system(str_make_dir.c_str()) != 0)
// 	{
// 		cout << "(ERROR) system cmd : " << str_make_dir << endl;
// 		exit(1);
// 	}
// */
// }