CPPFLAGS=-Wall -W -g -O0
CXX ?= g++

kaggleTitanic: string.o svm.o list.o bayesianNetwork.o laa.o BayesianInformationCriterion.o
clean:
	rm *.o preProcessing.exe