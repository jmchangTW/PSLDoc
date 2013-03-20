#ifndef _CONFIG_H
#define _CONFIG_H

#ifndef DATA_PATH
#define DATA_PATH "/projects/cjiaming/2010-01_classification/test/"
#endif

#ifndef PSIBLAST_CMD
#define PSIBLAST_CMD "/users/cn/jchang/program/ncbi-blast-2.2.22+/bin/psiblast"
#endif

#ifndef PSIBLAST_DATABASE_PATH
#define PSIBLAST_DATABASE_PATH "/projects/cjiaming/database/4blast/"
#endif

#ifndef PSIBLAST_DATABASE
#define PSIBLAST_DATABASE "uniprot"
#endif

#ifndef LIBSVM_EASYPY_CMD
#define LIBSVM_EASYPY_CMD "/users/cn/jchang/program/libsvm-2.89/tools/easy.py"
#endif

#ifndef SVM_SCALE_CMD
#define SVM_SCALE_CMD "svm-scale"
#endif

#ifndef SVM_PREDICT_CMD
#define SVM_PREDICT_CMD "svm-predict"
#endif

#endif /* _CONFIG_H */
