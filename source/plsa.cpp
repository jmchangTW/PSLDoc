#include <iostream>
#include <algorithm>
#include "PSLDoc.h"
#include "plsa.h"
#include "data.h"
#include "term.h"

#define RANDOM_MATRIX true
//Equal to matlabe function : normalize(m, 1)
//function [M, z] = normalise(A, dim)
//  z = sum(A);
//  s = z + (z==0);
//  M = A ./ repmat(s, size(A,1), 1);
Matrix::Matrix()
{
}

Matrix::Matrix(int size1, int size2)
{
	int i, j;
	matrix = gsl_matrix_alloc(size1, size2);
	gsl_matrix *norlMatrix = gsl_matrix_alloc(size1, size2);
	gsl_vector *normalV = gsl_vector_alloc(size2);
	float columSum;
	float rankNumber;
  	
	if(RANDOM_MATRIX)
	{
		srand(time(NULL));        	
		for(i=0; i < size2; i++)
		{
			columSum = 0;
			for(j=0; j < size1; j++)
			{
				rankNumber = rand()/(RAND_MAX+1.0);
				columSum += rankNumber;
				gsl_matrix_set(matrix, j, i, rankNumber);
			}
			if(columSum != 0)
				gsl_vector_set(normalV, i, columSum);
			else
				gsl_vector_set(normalV, i, 1);		
		}
	
		for(i=0; i < size1; i++)
			gsl_matrix_set_row(norlMatrix, i, normalV);
		gsl_matrix_div_elements(matrix, norlMatrix);
		gsl_vector_free(normalV);
		gsl_matrix_free(norlMatrix);	
	}
	else
	{
		for(i=0; i < size2; i++)
		{
			for(j=0; j < size1; j++)
				gsl_matrix_set(matrix, j, i, (float)1/size1);
		}
	}
}

Matrix::Matrix(int size1, int size2, const char *file_name)
{
	matrix = gsl_matrix_alloc(size1, size2);	
	FILE *file = fopen(file_name, "rb");
	gsl_matrix_fread(file, matrix);
	fclose(file);
}

Matrix::Matrix(int size1, int size2, const Data *data)
{
	matrix = gsl_matrix_alloc(size1, size2);
	
	for(int i = 0; i < (int)data->size(); i++)
	{
			for(int j = 0; j< (int)(*data)[i].fea_size(); j++)
				gsl_matrix_set(matrix, j, i, (*data)[i].get_fea(j));
	}
}

Matrix::~Matrix()
{
//	gsl_matrix_free(matrix); //add this code will cause "*** glibc detected *** double free or corruption"
}

void Matrix::show() const
{
	for(int i = 0; i < (int)matrix->size1; i++)
	{
		cout << i << "\t";
		for(int j = 0; j < (int)matrix->size2; j++)
			cout << gsl_matrix_get(matrix, i, j) << " ";
		cout << endl;			
	}
}


void Matrix::write_file2binary_format(const char *file_name)
{
	FILE *file = fopen(file_name, "wb");
	if(!file)
	{
//		system("hostname");
		cout << "ERROR: open file failed : " << file_name << endl;
		exit(0);
	}	
	gsl_matrix_fwrite(file, matrix);
	fclose(file);
}

void Matrix::write_file2csv_format(const char *file_name)
{
	int iDimOneSize = (int)matrix->size1;
	int iDimTwoSize = (int)matrix->size2;
	int i, j;
	float fElement;
	fstream file ;
	file.open(file_name, ios::out);

	if(!file)
	{
//		system("hostname");
		cout << "ERROR: open file failed : " << file_name << endl;
		exit(0);
	}	

	for(j = 0; j < iDimTwoSize; j++)
		file << ",Dim2-" << j << " ";
	file << endl;

	for(i = 0; i < iDimOneSize; i++)
	{
		file << "Dim1-" << i << " ";
		for(j = 0; j < iDimTwoSize; j++)
		{
			fElement = gsl_matrix_get(matrix, i, j);
			file << ","<< fElement << " ";
		}
		file << endl;
	}
	file << endl;
	file.close();
}

void Matrix::normalization()
{
	gsl_matrix *norlMatrix = gsl_matrix_alloc(matrix->size1, matrix->size2);	
	gsl_vector *normalV = gsl_vector_alloc(matrix->size2);
	float columSum;
        	
	for(int i=0; i < (int)matrix->size2; i++)
	{
		columSum = 0;
		for(int j=0; j < (int)matrix->size1; j++)
			columSum += gsl_matrix_get(matrix, j, i);

		if(columSum != 0)
			gsl_vector_set(normalV, i, columSum);
		else
			gsl_vector_set(normalV, i, 1);			
	}
	
	for(int i=0; i < (int)matrix->size1; i++)
		gsl_matrix_set_row(norlMatrix, i, normalV);

	gsl_matrix_div_elements(matrix, norlMatrix);
	gsl_vector_free(normalV);
	gsl_matrix_free(norlMatrix);
}

void Matrix::mergeByClass(Matrix *averageMatrix, const Data *data, const vector<string> &vClass)
{
	vector<string> vSite;
	int i, j, k, l;
	map<string, int> mapCount;
	map<string, float> mapSum;
   	map<string, int>::iterator iterCount;
	map<string, float>::iterator iterSum;
	float fElement;

	//initialization	
	for(i = 0; i < (int)vClass.size(); i++)
	{
		mapCount.insert(make_pair(vClass[i], 0));
		mapSum.insert(make_pair(vClass[i], 0));
	}
	
	for(i = 0; i < (int)matrix->size1; i++)
	{
		//initialization		
		for(iterCount = mapCount.begin(); iterCount != mapCount.end(); iterCount++)
			iterCount->second = 0;
		for(iterSum = mapSum.begin(); iterSum != mapSum.end(); iterSum++)
			iterSum->second = 0;
		vSite.clear();

		for(k = 0; k < (int)matrix->size2; k++)
		{
			vSite = (*data)[k].get_site_vec();
			fElement = gsl_matrix_get(matrix, i, k);
			for(l = 0; l < (int)vSite.size(); l++)
			{
				mapCount[vSite[l] ]++;				
				mapSum[vSite[l] ] += fElement;
			}
		}

		for(j = 0; j < (int)vClass.size(); j++)
			averageMatrix->set_element(i, j, (float)mapSum[vClass[j] ]/mapCount[vClass[j] ]);
	}
}

template <class T>
void Matrix::get_row_vector(int iRowIndex, vector<T> &vInput) const
{
	for(int i = 0; i < (int)matrix->size2; i++)
		vInput.push_back(gsl_matrix_get(matrix, iRowIndex, i));
};

bool cmp(const PAIR& a, const PAIR& b)
{
    return a.value > b.value;
}

void Matrix::select(vector<vector<PAIR> > &vvReturn)
{
	int i, j, putIndex;
	PAIR pTmp;
	int iDimOneSize = (int)matrix->size1;
	int iDimTwoSize = (int)matrix->size2;

	vvReturn.resize(iDimTwoSize);	
	for(i = 0; i < iDimOneSize; i++)
	{	
		pTmp.value = DBL_MIN;
		pTmp.index = i;
		putIndex = 0;
		for(j = 0; j < iDimTwoSize; j++)
		{
			if(gsl_matrix_get(matrix, i, j) > pTmp.value)
			{
				pTmp.value = gsl_matrix_get(matrix, i, j);
				putIndex = j;
			}
		}
		vvReturn[putIndex].push_back(pTmp);
	}
	for(i = 0; i < iDimTwoSize; i++)
		sort(vvReturn[i].begin(), vvReturn[i].end(), cmp);
};

void PLSA::EM(bool learn)
{
	gsl_matrix *wt_t = gsl_matrix_alloc(topic_size, word_size);
	gsl_matrix *td_t = gsl_matrix_alloc(doc_size, topic_size);
	gsl_matrix *wd_cal = gsl_matrix_alloc(word_size, doc_size);
	gsl_matrix *wd_div = gsl_matrix_alloc(word_size, doc_size);;
	gsl_matrix *wt_n = gsl_matrix_alloc(topic_size, doc_size);
	gsl_matrix *td_n = gsl_matrix_alloc(word_size, topic_size); 	
	if((wt_t == 0)||(td_t == 0)||(wd_cal == 0)||(wd_div == 0)||(wt_n == 0)||(td_n == 0))
	{
		cout << "ERROR in plsa_EM: gsl_matrix_alloc fail " << endl;
	}
		
//Original Matlab code		       
// while ~done; iter = iter+ 1;         
//	td = normalize( td .* ( wt' * ( counts ./ ( wt*td + eps ) ) ), 1);   
//   	if learn;
//    	   wt = normalize( wt.* (( counts ./ ( wt*td + eps ) ) *td' ), 1);
//   	end	         
//	end	

	cout << "(Start) EM-algorithm for "  << iteration << " iterations" << endl;
	for(int i = 0; i < iteration; i++)
	{
		cout << i+1 << ",";
		if((i+1)%30 == 0)
			cout << endl;
		else if((i+1)%10 == 0)
			cout << " ";

		gsl_matrix_transpose_memcpy(wt_t, wt.matrix);
		gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, wt.matrix, td.matrix, 0.0, wd_cal);
		gsl_matrix_add_constant(wd_cal, DBL_EPSILON); //DBL_EPSILON = eps, Machine epsilon, epsmch, is defined as the smallest positive number such that 1.0 + epsmch is not equal to 1.0.
		
		gsl_matrix_memcpy(wd_div, wd.matrix);
		gsl_matrix_div_elements(wd_div, wd_cal);
		
		gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, wt_t, wd_div, 0.0, wt_n);
		
		gsl_matrix_mul_elements(td.matrix, wt_n);
		td.normalization();
		
		if(learn)
		{
				gsl_matrix_transpose_memcpy(td_t, td.matrix);
				gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, wt.matrix, td.matrix, 0.0, wd_cal);		
				
				gsl_matrix_add_constant(wd_cal, DBL_EPSILON); //DBL_EPSILON = eps, Machine epsilon, epsmch, is defined as the smallest positive number such that 1.0 + epsmch is not equal to 1.0.				
				gsl_matrix_memcpy(wd_div, wd.matrix);
				gsl_matrix_div_elements(wd_div, wd_cal);
				
				gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, wd_div, td_t, 0.0, td_n);
				
				gsl_matrix_mul_elements(wt.matrix, td_n);
				wt.normalization();
		}
	}
	
	gsl_matrix_free(wt_t);
	gsl_matrix_free(td_t);
	gsl_matrix_free(wt_n);
	gsl_matrix_free(td_n);
	gsl_matrix_free(wd_cal);	
	gsl_matrix_free(wd_div);

	cout << endl;
	cout << "(Finish) EM-algorithm for "  << iteration << " iterations" << endl;	
}

PLSA::PLSA(const Data *data, const struct PSLDoc_PARAMETER *param, int word_size, int doc_size)
:td(param->feature_size, doc_size), wd(word_size, doc_size, data), wt(word_size, param->feature_size)
{
	iteration = param->iteration;
	this->word_size = word_size;
	this->topic_size = param->feature_size;
	this->doc_size = doc_size;
}

PLSA::PLSA(const Data *data, const char *wt_file_name, const struct PSLDoc_PARAMETER *param, int word_size, int doc_size):td(param->feature_size, doc_size), wd(word_size, doc_size, data), wt(word_size, param->feature_size, wt_file_name)
{
	iteration = param->iteration;
	topic_size = param->feature_size;
	this->word_size = word_size;
	this->doc_size = doc_size;
}

PLSA::PLSA(const Data *data, const char *wt_file_name, const char *td_file_name, const struct PSLDoc_PARAMETER *param, int word_size, int doc_size):td(param->feature_size, doc_size, td_file_name), wd(word_size, doc_size, data), wt(word_size, param->feature_size, wt_file_name)
{
	iteration = param->iteration;
	topic_size = param->feature_size;
	this->word_size = word_size;
	this->doc_size = doc_size;
}

PLSA::~PLSA()
{
}

void PLSA::write_matrix(string path, string format)
{
	string wt_path;
	string td_path;

	if (format == "binary")
	{
		wt_path = path + ".wt";
		td_path = path + ".td";
		cout << "Output wt matrix = " << wt_path << endl
			<< "Output td matrix = " << td_path << endl;
		wt.write_file2binary_format(wt_path.c_str());
		td.write_file2binary_format(td_path.c_str());
	}
	else if(format == "csv")
	{
		wt_path = path + "_wt.csv";
		td_path = path + "_td.csv";
		cout << "Output wt matrix = " << wt_path << endl
			<< "Output td matrix = " << td_path << endl;
		wt.write_file2csv_format(wt_path.c_str());
		td.write_file2csv_format(td_path.c_str());
	}
	else
	{
		cout << "[ERROR] unknown format: " << format << " for PLSA write_matrix function" << endl;
		exit(1);
	}
}

//         WT            TD
//       t               d
//     --------       --------
//    -              - 
//   w-             t- 
//    -              -
//
void PLSA::analyze(string sFilePath, const Data *data, const vector<string> &vClass)
{
	int SELECT_TOP_NUM = 5;
	int SELECT_WORD_NUM = 10; //# of signature = SELECT_TOP_NUM * SELECT_WORD_NUM
	int iTopicNumLoop,
	    iWordNumLoop,
 	    iSelectTopic,
            CLASS_NUM;
	gapped_dipeptide gdTmp;
	fstream file;

	cout << "STRAT function : PLSA::analyze" << endl;

	CLASS_NUM = (int)vClass.size();
	file.open(sFilePath.c_str(), ios::out);
	if(!file)	
	{
		perror(sFilePath.c_str());
		exit (0);
	}

	vector<vector<PAIR> > vvTW, vvCT;
	Matrix tc(td.get_size1(), CLASS_NUM);

	wt.select(vvTW);
	cout << "FINISH : wt select" << endl;

	td.mergeByClass(&tc, data, vClass);
	tc.write_file2csv_format("tc.csv");
	cout << "FINISH : tc merge" << endl;

	tc.select(vvCT);
	cout << "FINISH : tc select" << endl;

	//Output signature gapped-depeptide	
	//For each class
	for(int i = 0; i < CLASS_NUM; i++)
	{
		//For each topic
		file <<  "Class :" << vClass[i] << endl;
		if(SELECT_TOP_NUM > (int)vvCT[i].size())
			cout << "[NOTICE] # of preferred TOPIC < # set selected TOPIC: " << (int)vvCT[i].size() << " < " << SELECT_TOP_NUM << endl;
		iTopicNumLoop = (SELECT_TOP_NUM < (int)vvCT[i].size())?SELECT_TOP_NUM:(int)vvCT[i].size();
        	for(int j = 0; j < iTopicNumLoop; j++)
		{
			//For each word
			iSelectTopic = vvCT[i][j].index;
			if(SELECT_WORD_NUM > (int)vvTW[iSelectTopic].size())
				cout << "[NOTICE] # of preferred WORD < # set selected WORD: " << (int)vvTW[iSelectTopic].size() << " < " << SELECT_WORD_NUM << endl;
			iWordNumLoop = (SELECT_WORD_NUM < (int)vvTW[iSelectTopic].size())?SELECT_WORD_NUM:(int)vvTW[iSelectTopic].size();
			for(int k = 0; k < iWordNumLoop; k++)
			{
				gdTmp = gapped_dipeptide_term.get_term(vvTW[iSelectTopic][k].index);
				file << gdTmp.head << gdTmp.gapLen << gdTmp.tail << ", ";
			}
			file << endl;
		}
		file << endl;		
	}
	file.close();

	cout << "FINISH : output gapped-dipeptide" << endl;
}
