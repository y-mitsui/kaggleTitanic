#include <stdio.h>
#include <stdlib.h>

void fatalError(char *message,char *filename,int line){
	fprintf(stderr,"%s %s:%d\n",message,filename,line);
	exit(EXIT_FAILURE);
}
