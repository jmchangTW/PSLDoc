#ifndef _TEST_DATA_H
#define _TEST_DATA_H

#include <string>
#include <vector>
#include "data.h"
#include "train_data.h"

using namespace std;

class Test_data : public Data
{
public:
	Test_data(const char *fileName, bool bSelectAllData, const struct PSLDoc_PARAMETER *param, int index, const int *nr_fold);
	Test_data(const char *fileName, const struct PSLDoc_PARAMETER *param);
	void pred_kNN(Train_data train_data, int kNN);
	void show_pred(const char *) const;
	void parse_svm_result(const char * filename);
//	void Evaluation();
private:
	int class_num;
};

#endif /* _TEST_DATA_H */
