.KEEP_STAT:

all: PSLDoc-prepare PSLDoc-train PSLDoc-test PSLDoc-analyze

Compiler	= g++
FLAGS		= -Wall -fopenmp -D NDEBUG -O2 -I/soft/general/gsl-1.12/include/ -L/soft/general/gsl-1.12/lib/
# FLAGS		= -Wall -fopenmp -g -O2 -I/soft/general/gsl-1.12/include/ -L/soft/general/gsl-1.12/lib/
LIB		= -lgsl -lgslcblas -lm

OBJECT4PREPARE		= PSLDoc-prepare.o data.o protein.o term.o pssm.o test_data.o util.o

OBJECT4TRAIN		= PSLDoc-train.o data.o protein.o term.o pssm.o test_data.o train_data.o plsa.o util.o

OBJECT4TEST		= PSLDoc-test.o data.o protein.o term.o pssm.o test_data.o train_data.o plsa.o util.o

OBJECT4ANALYZE		= PSLDoc-analyze.o data.o protein.o term.o pssm.o test_data.o train_data.o plsa.o util.o

PSLDoc-prepare:	$(OBJECT4PREPARE)
			$(Compiler) $(FLAGS) $(OBJECT4PREPARE) -o PSLDoc-prepare $(LIB)
			mv -f PSLDoc-prepare ../.

PSLDoc-train: $(OBJECT4TRAIN)
			$(Compiler) $(FLAGS) $(OBJECT4TRAIN) -o PSLDoc-train $(LIB)
			mv -f PSLDoc-train ../.

PSLDoc-test: $(OBJECT4TEST)
			$(Compiler) $(FLAGS) $(OBJECT4TEST) -o PSLDoc-test $(LIB)
			mv -f PSLDoc-test ../.

PSLDoc-analyze: $(OBJECT4ANALYZE)
			$(Compiler) $(FLAGS) $(OBJECT4ANALYZE) -o PSLDoc-analyze $(LIB)
			mv -f PSLDoc-analyze ../.

%.o:		%.cpp
			$(Compiler) $(FLAGS) -c $<

clean:
			rm -f *.o

install:
			mv -f ../PSLDoc-prepare ~/bin/.
			mv -f ../PSLDoc-train ~/bin/.
			mv -f ../PSLDoc-test ~/bin/.
			mv -f ../PSLDoc-analyze ~/bin/.
