#ifndef _PSLDoc_H
#define _PSLDoc_H

#include <string>
#include "term.h"

using namespace std;

struct PSLDoc_PARAMETER
{
	int gapped_distance;
	string data_path;
	string PSIBLAST_iteration; 	//Set PSIBLAST parameter
	string PSIBLAST_e_value;
	int feature_size; 	//PLSA parameter
	int iteration;
	int overwrite_tfpssm;
	int pred_type;
};

extern Term gapped_dipeptide_term;

const int PARALLEL_NODE = 0;

#endif /* _PSLDoc_H */
