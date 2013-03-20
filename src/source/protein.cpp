#include <string>
#include <fstream>
#include <iostream>
#include <iomanip>

#include "config.h"
#include "term.h"
#include "pssm.h"
#include "PSLDoc.h"
#include "protein.h"
#include "util.h"

Protein::Protein()
{
}

Protein::Protein(string name_input, string seq_input, vector<string> site_input)
{
	name = name_input;
	seq = seq_input;
	for(int i = 0; i < (int)site_input.size(); i++)
		site_vec.push_back(site_input[i]);
}

Protein::~Protein()
{
}

void perform_blast(string tmp, string pssm_file_name, const struct PSLDoc_PARAMETER *param)
{
  	//fast,insensitive by book "BLAST" page 146, July 2003 version
	//blast version: -F "m S" -f 999 -M BLOSUM80 -G 9 -E 2 -e 1e-5
	string fast_par=" -matrix BLOSUM80 -evalue 1e-5 -gapopen 9 -gapextend 2 -threshold 999 -seg yes -soft_masking true -num_iterations " + param->PSIBLAST_iteration;
	string default_par=" -num_iterations "+param->PSIBLAST_iteration+" -evalue "+param->PSIBLAST_e_value;
        string blast_par = (getenv ("BLAST_MODE")=="fast"?fast_par:default_par);
 
//check whether blast is ok or not
	string blast_db = (getenv ("BLAST_DB")!= NULL?string(getenv ("BLAST_DB")): "uniprot");
	string str_db = (string)getenv("BLAST_PATH") + "/" + blast_db;
	string str_check_db = (string)getenv("BLAST_PATH") + "/" + blast_db + ".00.phr";

	if(check_pg_exist(PSIBLAST_CMD) == 0)
	{
		cout << "(ERROR) psi-blast cmd does not exist: " << PSIBLAST_CMD << endl;
		exit(1);
	}
	else if(!fexists(str_check_db))
	{
		cout << "(ERROR) psi-blast db " << str_db << " does not exist"  << endl;
		cout << "\tPlease pecify the path of blast database by the following cmd" << endl;
		cout << "\texport BLAST_PATH=/your_blast_path" << endl;
		exit(1);
	}

 	string cmd = (string)PSIBLAST_CMD + blast_par + " -db " + str_db + " -query \"" + tmp + "\" -out_ascii_pssm \"" + pssm_file_name + "\" &> /dev/null";

	cout << "RUN psiblast : " << cmd << endl;
	system(cmd.c_str());
	if(!fexists(pssm_file_name))
	{
	 cmd = (string)PSIBLAST_CMD + default_par + " -db " + str_db + " -query \"" + tmp + "\" -out_ascii_pssm \"" + pssm_file_name + "\" &> /dev/null";
	 cout << "RERUN psiblast : " << cmd << endl;
	 system(cmd.c_str());
	}
	
	cmd = "/bin/rm " + tmp;
	system(cmd.c_str());
}

void Protein::gen_PSSM(const char *dic_path, const struct PSLDoc_PARAMETER *param) const
{
	fstream file;
	string tmp = name + ".tmp";
	string pssm_file_name = dic_path + name + ".pssm";

	if(!fexists(pssm_file_name.c_str()))
	{
		file.clear(); file.open(tmp.c_str(), ios_base::out);
		file << ">" << name << endl << seq << endl;
		file.close();
		perform_blast(tmp, pssm_file_name, param);
	}
	else
		cout << "pssm_file exist: " << pssm_file_name << endl;
}

void Protein::gen_TFPSSM(const char *dic_path, Term term, const struct PSLDoc_PARAMETER *param) const
{
	int i;
	float termWei, weiSum=0;
	string pssm_file_name = dic_path + name + ".pssm";
	string tfpssm_file_name = dic_path + name + ".tfpssm";
	string tfpssm_csv_name = dic_path + name + ".tfpssm.csv";

	if(!fexists(tfpssm_file_name.c_str()) || (param->overwrite_tfpssm > 0))
	{
		cout << "process " << tfpssm_file_name << endl;
		float *saveVector = new float [(int)term.size()];
		
		PSSM pssm((int)seq.length(), pssm_file_name.c_str());
		//pssm.normal_smooth();
		//pssm.window_smooth();
		
		for(i = 0; i < (int)term.size(); i++)
		{
			termWei = pssm.gap_dip_fre(term.get_term(i));
			saveVector[i] = termWei;
			weiSum += termWei;
		}	
		//Normalization
		for(i = 0; i < (int)term.size(); i++)
				saveVector[i] = saveVector[i]/weiSum;
	
//output vector for binary
		fstream file;
		file.open(tfpssm_file_name.c_str(), ios::out|ios::binary);
		if(!file)
		{
			perror(tfpssm_file_name.c_str());
			exit(0);
		}
		file.write((char*)saveVector,sizeof(float)*(int)term.size());
		file.close();

//output vector for csv
		fstream csvFile;
		csvFile.open(tfpssm_csv_name.c_str(), ios::out);
		if(!csvFile)
		{
			perror(tfpssm_csv_name.c_str());
			exit(0);
		}
		csvFile << name << ",";
		for(i = 0; i < (int)site_vec.size(); i++)
			csvFile << site_vec[i] << ",";
		for(i = 0; i < (int)term.size(); i++)
			csvFile << saveVector[i] << ",";
		csvFile << endl;
		csvFile.close();

		delete [] saveVector;
	}
	else if(fexists(tfpssm_file_name))
		cout << "tfpssm_file exist: " << tfpssm_file_name << endl;
	
// 	for(int i = 0; i < 10; i++)
// 	{
// 	   if(fexists(tfpssm_file_name))
// 	     break;
// 	   sleep(10);
// 	}	
}

//From File
void Protein::construct_fea_vec(const char*tfpssm_dic)
{
	fstream file;	
	float *tmp_vec = new float [gapped_dipeptide_term.size()];
	string tfpssm_file = tfpssm_dic + name + ".tfpssm";
		
	file.open(tfpssm_file.c_str(), ios::in|ios::binary);
	if(!file)
	{
		perror(tfpssm_file.c_str());
		exit(0);
	}	
	file.read((char*)tmp_vec,sizeof(float)*gapped_dipeptide_term.size());
	file.close();

	for(int i = 0; i < gapped_dipeptide_term.size(); i++)
		fea_vec.push_back(tmp_vec[i]);
							
	delete [] tmp_vec;
}

//From PLSA reduced vector
void Protein::construct_fea_vec(const gsl_vector *gsl_vec)
{
	fea_vec.clear();
	
	for(int i = 0; i < (int)gsl_vec->size; i++)
		fea_vec.push_back(gsl_vector_get(gsl_vec,i));
}

void Protein::set_pred(Protein pro)
{
	PRE_ELEMENT tmp_pre;
	for(int i = 0; i < (int)pro.site_vec.size(); i++)
	{
		tmp_pre.pre_label = pro.site_vec[i];
		pre_vec.push_back(tmp_pre);
	}
}

void Protein::show_fea(fstream& out) const
{
	for(int i = 0; i < (int)fea_vec.size(); i++)
		out << i+1 << ":" << fixed << setprecision(6) << fea_vec[i] << " ";
}

void Protein::show_site(fstream& out) const
{
	for(int j =0; j < (int)site_vec.size(); j++)
		out << site_vec[j] << " ";
}

void Protein::show_pred(fstream& out) const
{
	for(int j =0; j < (int)pre_vec.size(); j++)
		out << pre_vec[j].pre_label << " ";
}

void Protein::show_probability(fstream& out) const
{
	for(int i = 0; i < (int)pre_vec.size(); i++)
		for(int j = 0; j < (int)pre_vec[i].confidence.size(); j++)
			out << pre_vec[i].confidence[j] << ",";
}

int Protein::pred_correct() const
{
	int i;
	map<string, int> site_map;
	int count;
	string pred_site;
	
	count = 0;
	if((int)site_vec.size() == 1)
	{
		for(i =0; i < (int)pre_vec.size(); i++)
		{
			pred_site = pre_vec[i].pre_label;
			if(site_map.find(pred_site) == site_map.end())
				site_map.insert(make_pair(pred_site,1));
			else
				site_map[pred_site]++;
		}
		for(map<string, int>::iterator iter = site_map.begin(); iter != site_map.end(); iter++)
		{
			if(iter->second > count)
			{
				count = iter->second;
				pred_site = iter->first;
			}
		}
		
		if(pred_site == site_vec.front())
			return 1;
		else
			return 0;
	}
	else
	{
			for(i =0; i < (int)site_vec.size(); i++)
			{
				if(site_vec[i] != pre_vec[i].pre_label)
					return 0;
			}
			return 1;
	}
}
