#include <iostream>
#include <fstream>
#include "term.h"
#include "protein.h"
#include "train_data.h"

Train_data::Train_data(const char *file_name, bool bSelectAllData, const struct PSLDoc_PARAMETER *param,  int index, const int *nr_fold)
:Data(file_name, param)
{
	int i;
	vector<Protein>::iterator iter;
		
	for(iter = vec.begin(), i = 0; iter != vec.end(); iter++, i++) 
	{
		if(bSelectAllData||((i%*nr_fold) != index)) //bSelectAllData: select all
			iter->construct_fea_vec(data_path.c_str());
		else
			vec.erase(iter--);
	}
}

Train_data::Train_data(const char *file_name, const struct PSLDoc_PARAMETER *param)
:Data(file_name, param)
{
	int i;
	vector<Protein>::iterator iter;
		
	for(iter = vec.begin(), i = 0; iter != vec.end(); iter++, i++) 
		iter->construct_fea_vec(data_path.c_str());
}
