#ifndef _DATA_H
#define _DATA_H

#include <string>
#include <vector>
#include "protein.h"

using namespace std;

class Matrix;

//===Data===
class Data
{
public:
	Data();
	Data(const char *fileName, const struct PSLDoc_PARAMETER *param);
	~Data();
	int size() const{return (int)vec.size();} ;
	Protein operator[](int i) const;
	string get_data_path(){return data_path;};
	void protein2feature(int i, const struct PSLDoc_PARAMETER *param, Term term) const;
	void gen_distance_matrix();
	void show_pro(int i){
		cout << vec[i].get_name() << endl;
//		if(i<(int)vec.size()) 
//			vec[i].show_fea();
	};
	void update_vec(Matrix td);
	void gen_svm_file(const char* filename) const;
	void getClass(vector<string> &vClass) const;
protected:
	string data_path;
	vector<Protein> vec;
};

#endif /* _DATA_H */
