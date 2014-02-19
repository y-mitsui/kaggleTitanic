#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <limits.h>
#include "kaggleTitanic.h"
#include "error.h"
#include "list.h"
#include "laa.h"

double bic(double logLike,int numSample,int numModel);

typedef struct _baysianNode{
	const char *name;
	int value;
	int numVariablePattern;
	double *probability;
}baysianNode;

typedef struct {
	int *edge;
	double score;
}bayesianNetScore;
typedef struct{
	bayesianNetScore **models;
	int numModels;
	Associate ***CPT;
}bayesianNetwork;
bayesianNetScore* makeBNScore(int *edge,int numEdge,double score){
	bayesianNetScore *r=Malloc(bayesianNetScore,1);
	r->edge=Malloc(int,numEdge);
	memcpy(r->edge,edge,sizeof(int)*numEdge);
	r->score=score;
	return r;
}
void printfBayesianNetwork(baysianNode *node,int numNode,int *edge){
	int i;
	puts("---------");
	for(i=0;i<numNode*numNode;i++){
		if(edge[i]==1){
			printf("%s -> %s\n",node[i%numNode].name,node[(int)i/numNode].name);
		}
	}
	puts("");
	
}
void printSquereMatrix(int *matrix,int numMatrix){
	int width=sqrt(numMatrix);
	int i,j;
	puts("--------------");
	printf("  ");
	for(i=0;i<width;i++){
		printf("%c ",'A'+i);
	}
	puts("");
	for(i=0;i<width;i++){
		printf("%c ",'A'+i);
		for(j=0;j<width;j++){
			printf("%d ",matrix[i*width+j]);
		}
		puts("");
	}
	puts("--------------");
}

int getNumParents(int *edge,int currentNode,int numNode){
	int i,num=0;
	for(i=0;i<numNode;i++){
		if(edge[currentNode*numNode+i]==1){
			num++;
			num+=getNumParents(edge,i,numNode);
		}
	}
	return num;
}
double likelihood(int numNode,int *edge,double **likelies){
	double like;
	int i,j,sum;
	unsigned char parentPattern;

	for(like=0.0,i=0;i<numNode;i++){
		for(parentPattern=0,j=0;j<numNode;j++){
			if(edge[i*numNode+j]==1) parentPattern=parentPattern|(1<<j);
		}
		like+=likelies[i][parentPattern];
	}
	
	/* 親の総数*/
	for(sum=0,i=0;i<numNode;i++){
		for(j=0;j<numNode;j++){
			if(edge[i*numNode+j]==1) sum++;
		}
	}
	return -2*like+sum;
}
int isExistChildrenChain(int *edge,int numNode,int child,int key){
	int i;
	if(edge[key*numNode+child]==1) return 0;
	for(i=0;i<numNode;i++){
		if(edge[i*numNode+child]==1){
			if(isExistChildrenChain(edge,numNode,i,key)==0) return 0;
		}
	}
	return 1;
}
int isEnableArc(int *edge,int numNode,int idx){
	int parent=idx%numNode;
	int child=idx/numNode;

	if(parent==child) return 0;
	if(edge[idx]==1) return 0;
	return isExistChildrenChain(edge,numNode,child,parent);
	
}

int __buildBayesianNet(baysianNode *node,int numNode,int *edge,int numLine,bayesianNetScore **likely,double **likelies,int start,int rank){
	int numLikely=0;
	int i;
	double tmp;
	if(rank==numLine){
		tmp=likelihood(numNode,edge,likelies);
		likely[numLikely++]=makeBNScore(edge,numNode*numNode,tmp);
		return numLikely;
	}else{
		for(i=start;i<numNode*numNode;i++){
			if(isEnableArc(edge,numNode,i)){
				edge[i]=1;
				numLikely+=__buildBayesianNet(node,numNode,edge,numLine,&likely[numLikely],likelies,i+1,rank+1);
				edge[i]=0;
			}
		}
	}
	return numLikely;
}
int combination(int n,int m){
	double sum=1.0,sum2=1.0;
	int i;
	for(i=0;i<m;i++){
		sum*=n-i;
	}
	for(i=0;i<m;i++){
		sum2*=(m-i);
	}
	return (int)(sum/sum2);
}
unsigned long long getNumberOfBayesianNode(int limit){
	int i;
	unsigned long long sum;
	if(limit==0) return 1;
	for(sum=0,i=1;i<=limit;i++){
		unsigned long long nn=pow(-1,i+1)*combination(limit,i)*pow(2,i*(limit-i));
		if(limit-i > 1) nn*=getNumberOfBayesianNode(limit-i);
		sum+=nn;
	}
	return sum;
}
int *int2Number(int num){
	int *r=Malloc(int,1);
	*r=num;
	return r;
}
int sortCmp(const void *x,const void *y){
	bayesianNetScore **a=(bayesianNetScore **)x;
	bayesianNetScore **b=(bayesianNetScore **)y;
	if((*a)->score > (*b)->score) return 1;
	else if((*a)->score < (*b)->score) return -1;
	return 0;
}
void localScore(baysianNode *node,int numNode,int currentNode,int numLine,Associate *cHist,unsigned char parentPattern,double *likelies,Associate **CPT,int rank,int start){
	int *pattern;
	int *nums;
	Associate *tmp;
	double like;
	int i,j,patternCt,total;

	if(numLine==rank){
		/* 子ノードcurrentNodeにおける親パターンedgeについて、条件付き確率を計算する */
		tmp=makeLAA(cHist->numKeys,numLine);	//親値パターンと子ノード値ごとの総数
		pattern=Calloc(int,numNode);
		for(i=0;i<cHist->numKeys;i++){	
			for(patternCt=0,j=0;j<numNode;j++){
				if((parentPattern&((unsigned char)1<<j))!=0){
					pattern[patternCt++]=cHist->keys[i*cHist->patternSize+j];
				}
			}
			if(!(nums=(int*)getLAA(tmp,pattern))) nums=Calloc(int,node->numVariablePattern);
			for(j=0;j<node->numVariablePattern;j++){
				nums[j]+=((int*)cHist->array[i])[j];
			}
			setLAA(tmp,pattern,nums);
		}
		/* 条件付き確立を計算*/
		CPT[parentPattern]=makeLAA(tmp->numKeys,numLine);
		
		for(i=0;i<tmp->numKeys;i++){
			double *prop=Malloc(double,node->numVariablePattern);
			for(total=0,j=0;j<node->numVariablePattern;j++){
				total+=((int*)tmp->array[i])[j];
			}
			for(j=0;j<node->numVariablePattern;j++){
				prop[j]=(double)((int*)tmp->array[i])[j]/total;
				if(prop[j]==0.0) prop[j]=0.00000001;
			}			
			setLAA(CPT[parentPattern],&tmp->keys[i*tmp->patternSize],prop);
		}
		
		/* 尤度を計算 */
		like=0.0;
		total=0;
		for(i=0;i<tmp->numKeys;i++){
			double* prop=(double*)getLAA(CPT[parentPattern],&tmp->keys[i*tmp->patternSize]);
			for(j=0;j<node->numVariablePattern;j++){
				like+=log(prop[j])*((int*)tmp->array[i])[j];
				total+=((int*)tmp->array[i])[j];
			}
			free(tmp->array[i]);
		}
		likelies[parentPattern]=like;
		free(pattern);
		freeLAA(tmp);
		free(tmp);
	}else{
		for(i=start;i<numNode;i++){
			if(i!=currentNode){
				localScore(node,numNode,currentNode,numLine,cHist,parentPattern | (1<<i),likelies,CPT,rank+1,i+1);
			}
		}
	}
}
void __freeCPT(Associate **CPT,int numNode,int currentNode,unsigned char parentPattern,int numLine,int rank,int start){
	int i;
	if(numLine==rank){
		for(i=0;i<CPT[parentPattern]->numKeys;i++){
			free(CPT[parentPattern]->array[i]);
		}
		freeLAA(CPT[parentPattern]);
		free(CPT[parentPattern]);
	}
	for(i=start;i<numNode;i++){
		if(i!=currentNode){
			__freeCPT(CPT,numNode,currentNode,parentPattern | (1<<i),numLine,rank+1,i+1);
		}
	}
}
void freeCPT(Associate ***CPT,int numVariable){
	int i,j;

	for(i=0;i<numVariable;i++){
		for(j=0;j<numVariable;j++){
			__freeCPT(CPT[i],numVariable,i,0,j,0,0);
		}
		free(CPT[i]);
	}
}
bayesianNetwork* bayesianNetworkTrain(list_t *passengerList){
	Associate *frequency;
	int i,j,k;
	int *pattern;
	int *edge;
	int numLikely;
	cell_t *cur;
	bayesianNetwork *model;

	/* 分割表を作成 */
	frequency=makeLAA(passengerList->size,NUM_VARIABLE);
	pattern=Malloc(int,NUM_VARIABLE);
	for(cur=passengerList->first;cur;cur=cur->next){
		passenger *human=(passenger *)cur->data;
		pattern[0]=(human->survived==-1) ? 0 : 1;
		pattern[1]=human->sex;
		pattern[2]=human->rank-1;
		pattern[3]=(human->age < 10) ? 0 : (human->age < 20) ? 1 :  (human->age < 30) ? 2  : (human->age < 40) ?  3  : (human->age < 50) ?  4 : (human->age < 60) ?  5 : 6;
		pattern[4]=human->prop1;

		int* cur=(int*)getLAA(frequency,pattern);
		if(cur==NULL){
			setLAA(frequency,pattern,int2Number(1));
		}else{
			(*cur)++;
		}
	}
	/* 条件付き頻度表を作成*/	
	Associate **cHist=Malloc(Associate*,NUM_VARIABLE);
	for(i=0;i<NUM_VARIABLE;i++){
		cHist[i]=makeLAA(passengerList->size,NUM_VARIABLE);
		for(j=0;j<frequency->numKeys;j++){
			for(k=0;k<frequency->patternSize;k++){
				pattern[k]=(i==k) ? 0 : frequency->keys[j*frequency->patternSize+k];
			}
			int val=frequency->keys[j*frequency->patternSize+i];
			int *cur=(int*)getLAA(cHist[i],pattern);
			if(cur) cur[val]+=*((int*)frequency->array[j]);
			else{
				cur=Calloc(int,node[i].numVariablePattern);
				cur[val]=*((int*)frequency->array[j]);
			}
			setLAA(cHist[i],pattern,cur);
		}
	}

	/* ローカルスコアを計算*/
	Associate ***CPT=Malloc(Associate **,NUM_VARIABLE);
	double **likelies=Malloc(double *,NUM_VARIABLE);
	for(i=0;i<NUM_VARIABLE;i++){
		likelies[i]=Malloc(double,pow(2,NUM_VARIABLE));
		CPT[i]=Malloc(Associate *,pow(2,NUM_VARIABLE));
		for(j=0;j<NUM_VARIABLE;j++){
			localScore(&node[i],NUM_VARIABLE,i,j,cHist[i],0,likelies[i],CPT[i],0,0);
		}
	}

	int  numCombinationPattern=getNumberOfBayesianNode(NUM_VARIABLE);	//組み合わせ数を計算
	bayesianNetScore **arr=Malloc(bayesianNetScore*,numCombinationPattern);
	edge=Calloc(int,NUM_VARIABLE*NUM_VARIABLE);
	for(numLikely=0,i=0;i<=10;i++){
		numLikely+=__buildBayesianNet(node,NUM_VARIABLE,edge,i,&arr[numLikely],likelies,0,0);	
	}

	free(edge);
	for(i=0;i<NUM_VARIABLE;i++){
		free(likelies[i]);
		for(j=0;j<cHist[i]->numKeys;j++){
			free(cHist[i]->array[j]);
		}
		freeLAA(cHist[i]);
		free(cHist[i]);
	}
	free(likelies);
	free(cHist);
	for(i=0;i<frequency->numKeys;i++)
		free(frequency->array[i]);
	freeLAA(frequency);
	free(frequency);
	free(pattern);


	qsort(arr,numLikely,sizeof(bayesianNetScore*),sortCmp);

	model=Malloc(bayesianNetwork,1);
	model->models=arr;
	model->numModels=numLikely;
	model->CPT=CPT;
	return model;
}
void *bayesianNetworkPredict(bayesianNetwork *model,passenger *human){
	pattern[0]=human->sex;
	pattern[1]=human->rank-1;
	pattern[2]=(human->age < 10) ? 0 : (human->age < 20) ? 1 :  (human->age < 30) ? 2  : (human->age < 40) ?  3  : (human->age < 50) ?  4 : (human->age < 60) ?  5 : 6;
	pattern[3]=human->prop1;
	double *prop=(double*)getLAA(model->CPT[0][2|4|8|16],pattern);
	if(prop==NULL){
		fprintf(stderr,"prop is NULL");
		exit(1);
	}
	int realServived=(human->survived==-1) ? 0 : 1;
	like+=log(prop[realServived]);
	int serv=(prop[0] >= prop[1]) ? 0 : 1;
}
void freeBayesianNetwork(bayesianNetwork *model){
	for(i=0;i<model->numModels;i++){
		free(model->models[i]->edge);
		free(model->models[i]);
	}
	free(model->models);
	freeCPT(model->CPT,NUM_VARIABLE);
	free(model->CPT);
	
}
int bayesianNetworkIsEqual(int *x,passenger *human){
	return *x==human->servived;
}
double clossValidation(list_t *data, void* (*train)(list_t *),void* (*predict)(void*),int isEqual(void*,void*),void (*freeModel)(void*)){
	void *model;
	void *val;
	int hit;
	for(hit=0,cur=data->first;cur;cur=cur->next){
		removeList(data,cur);
		model=train(list);
		val=predict(model,cur->data);
		freeModel(model);
		if(isEqual(val,cur->data)) hit++;
		revivalList(list,cur);
	}
	return (double)hit/data->size;
}
#define NUM_VARIABLE 5
void bayesianNetwork(list_t *passengerList,list_t *testPassengerList){
	baysianNode node[NUM_VARIABLE];
	clock_t time1;

	memset(&node,0,sizeof(node[0])*NUM_VARIABLE);
	node[0].name="SURVIVED";
	node[0].numVariablePattern=2;
	node[1].name="SEX";
	node[1].numVariablePattern=2;
	node[2].name="RANK";
	node[2].numVariablePattern=3;
	node[3].name="AGE";
	node[3].numVariablePattern=7;
	node[4].name="FRIEND";
	node[4].numVariablePattern=3;

	double cv=clossValidation(passengerList,bayesianNetworkTrain,bayesianNetworkPredict,bayesianNetworkIsEqual,freeBayesianNetwork);
	printf("CV:%lf\n",cv);
	time1=clock();
	
}
