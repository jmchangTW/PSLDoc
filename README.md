PSLDoc uses gapped-dipeptides and probabilistic latent semantic analysis solve 
to prediction protein subcellular localization.

Table of Contents
=================

- Installation
- Data Format
- 'PSLDoc-prepare' Usage
- 'PSLDoc-train' Usage
- 'PSLDoc-test' Usage
- 'PSLDoc-analyze' Usage
- Examples
- Additional Information

Installation
============================

1.Install needed softwares

PSLDoc needs the following programs.
 - GNU gsl:	http://www.gnu.org/software/gsl/
 - psiblast & database
 - libsvm:	http://www.csie.ntu.edu.tw/~cjlin/libsvm/

2.Specify their paths

* makefile GNU gsl  
Please modify the following line in makefile according to your gsl path

		FLAGS		= -Wall -D NDEBUG -O2 -I/soft/general/gsl-1.12/include/ -L/soft/general/gsl-1.12/lib/

* PSLDoc.cfg  
Please specify the path which PSSM and TFPSSM files store.

		DATA_PATH = your_path

* psiblast  
Please make sure that blastpgp is available in shell command or you could modify the following line in the file.

		PSIBLAST_CMD = blastpgp
 		PSIBLAST_DATABASE = nr
 		PSIBLAST_DATABASE_PATH = your_path # if you have different setting with default blast setting, please sepecify and umark

* libsvm  
Please make sure that svm-scale and svm-predict are available in shell command or you could modify the following lines in the file.

		LIBSVM_EASYPY_CMD = /users/yourname/program/libsvm-2.89/tools/easy.py
		SVM_SCALE_CMD = svm-scale
		SVM_PREDICT_CMD = svm-predict

* grid.py  
You also have modify svmscale_exe, svmtrain_exe, svmpredict_exe and grid_py paths in two files, easy.py and grid.py, in libsvm "tools" directory.
For instance,

		svmscale_exe = "/users/your_name/program/libsvm-2.89/svm-scale"
		svmtrain_exe = "/users/your_name/program/libsvm-2.89/svm-train"
		svmpredict_exe = "/users/your_name/program/libsvm-2.89/svm-predict"
		grid_py = "/users/your_name/program/libsvm-2.89/tools/grid.py"

* easy.py  
Train SVM model for supporting probability estimates, please replace the following line

		cmd = '%s -c %s -g %s "%s" "%s"' % (svmtrain_exe,c,g,scaled_file,model_file)
by

		cmd = '%s -b 1 -c %s -g %s "%s" "%s"' % (svmtrain_exe,c,g,scaled_file,model_file)

3.On Unix systems,  
type `./configure.pl` to generate config.h file  
type `make` to build the PSLDoc-prepare, PSLDoc-train and PSLDoc-test programs.  
Run them without arguments to show the usages of them.

Data Format
============================

The format of training and testing data file is like FAST format:

For training file:

\>protein_name class=class_label;  
protein_seq

For example:

\>3122878 class=1;  
MPLDLYNTLTRRKERFEPMTPDRVGMYVCG...

For testing file with label or without label:

\>protein_name class=class_label;  
protein_seq

or

\>protein_name  
protein_seq

For example:

\>3914018 class=2;  
MSISMTTKLSYGFGAFGKDFAIGIVYMYLMY...

or 

\>3914018  
MSISMTTKLSYGFGAFGKDFAIGIVYMYLMY...

'PSLDoc-prepare' Usage
=================

	PSLDoc-prepare [options] data_set_file  
	options:  
        	-o overwrite_TFPSSM     overwrite TFPSSM data (default 0)  
                	0 -- No  
                	1 -- Yes	
        	-s signature    set signature file (NULL)  
        	-d distance     set gapped-dipeptide distance (default 13)  
        	-j loop         set the number of the Loop in PSIBLAST (default 3)  
        	-e value        set the e-value of PSIBLAST (default 0.01)  

For instance:
>PSLDoc-prepare -e 0.001 data/simple_train


'PSLDoc-train' Usage
===================

	PSLDoc-train [options] train_file
	options:
        -r reduction : feature reduction or not (default 1)
                0 -- Do not perform feature reduction
                1 -- Perform feature reduction by PLSA
        -s signature : set signature file (NULL)
        -d distance : set gapped-dipeptide distance (default 13)
        -f size : set the reduced feature size of PLSA (default 80)
        -i iteration : set the iteration of PLSA (default 300)
        -e exist_data : use existing data (default 1)
                0 -- No
                1 -- Yes

For instance:
>PSLDoc-train data/simple_train


'PSLDoc-test' Usage
=================

	PSLDoc-test [options] model_file test_data_file
	options:
        -r reduction : feature reduction or not (default 1)
                0 -- Do not perform feature reduction
                1 -- Perform feature reduction by PLSA
        -p predict_method : set prediction method (default 1)
                0 -- k Nearest Neighbor
                1 -- SVM
        -s signature : set signature file (NULL)
        -k # of neighgors : set the number of Nearest Neighbors (default 1)
        -d distance : set gapped-dipeptide distance (default 13)
        -f size : set the reduced feature size of PLSA (default 80)
        -i iteration : set the iteration of PLSA (default 300)
        -b probability : whether to output the probabiliy of SVM prediction (default 1)
                0 -- Do not output
                1 -- Output

For instance:
>PSLDoc-test data/simple_train data/simple_test
<pre>
Output:
        SVM input = data/simple_test.svm_input
        prediction result = data/simple_test.predict
        prediction result (csv format) = data/simple_test.csv

For instance:
>PSLDoc-test data/simple_train data/simple_test

NOTICE:
Predictions may be different with and without probability estimation.
http://www.csie.ntu.edu.tw/~cjlin/libsvm/faq.html#f501
Q: Why using svm-predict -b 0 and -b 1 gives different accuracy values?
Let's just consider two-class classification here. After probability information is obtained in training, we do not have
prob > = 0.5 if and only if decision value >= 0.
So predictions may be different with -b 0 and 1.


'PSLDoc-analyze' Usage
=================

	PSLDoc-analyze [options] data_set_file
	options:
        -s signature : set signature file (NULL)
        -d distance : set gapped-dipeptide distance (default 13)

For instance
>PSLDoc-analyze data/simple_train  
<pre>
Output:
       	Topic vs localization class = tc.csv
       	Analyze result = data/simple_train.analysis

See 'Examples' in this file for examples.

Examples
========
You have two files.
One is a training file, data/simple_train.
Other one is a testing file, data/simple_test.

Use the following four steps to perform PSLDoc prediction for the testing file based on the training file.

Generate PSSM and TFPSSM files for data/simple_train and data/simple_test
> PSLDoc-prepare data/simple_train
> PSLDoc-prepare data/simple_test

Perform PLSA and SVM training on data/simple_train
> PSLDoc-train data/simple_train

Perform PLSA fold in and SVM testing on data/simple_test based on the previous training result, data/simple_train
> PSLDoc-test data/simple_train data/simple_test

Perform PLSA analysis for selecting preferred topics for each class based on previous trained wt and td matrix, that is, you have to make sure two binary files, data/simple_train.wt and data/simple_train.td, exist.
> PSLDoc-analyze data/simple_train

Visulize wt, td and tc matrix.
You have to run 'PSLDoc-train' and 'PSLDoc-analysis' to get three wt, td and tc matrix in csv format.

After the following two steps.
> PSLDoc-train data/simple_train
> PSLDoc-analyze data/simple_train

You could use 'tools/visMatrix.R' and then get three jpg files, 'wt.jpg', 'td.jpg' and 'tc.jpg'.
visMatrix('data/simple_train_wt.csv', 'data/simple_train_td.csv', 'tc.csv')

Additional Information
======================

If you find PSLDoc helpful, please cite it as  
Chang, J.M., Su, Emily C.Y., Lo, A., Chiu, H.S., Sung, T.Y. and Hsu, W.L. (2008) Protein Subcellular Localization Prediction based on based on gapped-dipeptides and probabilistic latent semantic analysis. PROTEINS: Structure, Function, and Bioinformatics,72, 693-710.

Creative Commons License 
Attribution-Noncommercial-Share Alike 2.5 Taiwan License. 
http://creativecommons.org/licenses/by-nc-sa/2.5/tw/
