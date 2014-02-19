CPPFLAGS=-Wall -W -g
CXX ?= g++

preProcessing: string.o svm.o list.o naiveBeise.o bayesianNetwork.o laa.o BayesianInformationCriterion.o
string:
