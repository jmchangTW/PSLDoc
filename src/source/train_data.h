#ifndef _TRAIN_DATA_H
#define _TRAIN_DATA_H

#include <string>
#include <vector>
#include "data.h"

using namespace std;

//===Data===
class Train_data : public Data
{
public:
  Train_data(const char *fileName, bool bSelectAllData, const struct PSLDoc_PARAMETER *param, int index, const int *nr_fold);
  Train_data(const char *fileName, const struct PSLDoc_PARAMETER *param);
};

#endif /* _TRAIN_DATA_H */
