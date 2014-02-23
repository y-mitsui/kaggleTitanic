CPPFLAGS=-Wall -W -g -O2
CXX ?= g++

preProcessing: string.o svm.o list.o bayesianNetwork.o laa.o BayesianInformationCriterion.o
clean:
	rm *.o preProcessing.exe