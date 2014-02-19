#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "kaggleTitanic.h"

int explode(char ***tokens,char *str,const char demeliter){
	int open=0;
	int n=0,i;
	int tokensNum=0;
	int len=strlen(str);
	char clumn[1024];
	char **token;
	int limit=100;

	token=Malloc(char*,limit);
	for(i=0;i<len;i++){
		if(str[i]=='\"'){
			open=!open;
		}else if(str[i]==demeliter && !open){
			clumn[n]=0;
			token[tokensNum++]=strdup(clumn);
			n=0;
		}else
			clumn[n++]=str[i];
	}
	clumn[n++]=0;
	token[tokensNum++]=strdup(clumn);

	*tokens=token;
	return tokensNum;
}
