#ifndef _PLSA_H
#define _PLSA_H

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_sort_double.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_math.h>
#include <fstream>
#include <map>

using namespace std;

struct PAIR
{
	int index;
	float value;
};

class Data;

class Matrix
{
public:
	Matrix();
	Matrix(int size1, int size2);
	Matrix(int size1, int size2, const char *file_name);
	Matrix(int size1, int size2, const Data *data);
	~Matrix();
	void write_file2binary_format(const char *file);
	void write_file2csv_format(const char *file);
	void show_size() const { cout << "Size1 " << matrix->size1 << endl << "Size2 " << matrix->size2 << endl;	};
	int get_size1() const { return matrix->size1;	};
	int get_size2() const { return matrix->size2;	};
	void show() const;
	template <class T>
	void get_row_vector(int iRowIndex, vector<T> &vInput) const;
	void normalization();
	void set_element(int i, int j, float value){gsl_matrix_set(matrix, i, j, value);};
	void mergeByClass(Matrix *averageMatrix, const Data *data, const vector<string> &vClass);
	void select(vector<vector<PAIR> > &vvReturn);
	friend class PLSA;
	friend class Data;
private:
	gsl_matrix *matrix;
};


class PLSA
{
public:
	PLSA(const Data *data, const struct PSLDoc_PARAMETER *param, int w_num, int d_num);
	PLSA(const Data *data, const char *wt_name, const struct PSLDoc_PARAMETER *param, int w_num, int d_num);
	PLSA(const Data *data, const char *wt_name, const char *td_name, const struct PSLDoc_PARAMETER *param, int w_num, int d_num);
	~PLSA();
	void EM(bool learn);
	void write_matrix(string path, string format);
	void analyze(string sFilepPath, const Data *data, const vector<string> &vClass);
	Matrix td;
private:	
	Matrix wd;
	Matrix wt;
	int iteration;
	int word_size;
	int topic_size;
	int doc_size;
};

#endif /* _PLSA_H */
