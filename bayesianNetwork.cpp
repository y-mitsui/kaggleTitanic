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

#define NUM_LIKELIES 50

double bic(double logLike,int numSample,int numModel);
void setVariable(int *pattern,passenger *human);

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
typedef struct{
	baysianNode *node;
	int numVariable;
}bayesianNetworkOption;
bayesianNetScore* makeBNScore(int *edge,int numEdge,double score){
	bayesianNetScore *r=Malloc(bayesianNetScore,1);
	r->edge=Malloc(int,numEdge);
	memcpy(r->edge,edge,sizeof(int)*numEdge);
	r->score=score;
	return r;
}

clock_t times[3];
clock_t first_time;

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

double likelihood(int numNode,int *edge,double **likelies){
	double like;
	int i,j,sum;
	unsigned short parentPattern;

	for(sum=0,like=0.0,i=0;i<numNode;i++){
		for(parentPattern=0,j=0;j<numNode;j++){
			if(edge[i*numNode+j]==1){
				parentPattern=parentPattern|(1<<j);
				sum++; // 親の総数
			}
		}
		like+=likelies[i][parentPattern];
	}
	
	
	return -2*like+sum;
}
int isExistChildrenChain(int *edge,int numNode,int child,int key){
	int i;
	if(edge[key*numNode+child]==1) return 0;
	register int *pEdge=edge+child;
	for(i=0;i<numNode;i++,pEdge+=numNode){
		if(*pEdge==1){
		//if(edge[i*numNode+child]==1){
			if(isExistChildrenChain(edge,numNode,i,key)==0) return 0;
		}
	}
	return 1;
}
int isEnableArc(int *edge,int numNode,int idx){

	first_time=clock();

	int parent=idx%numNode;
	int child=idx/numNode;

	if(parent==child) return 0;
	if(edge[idx]==1) return 0;
	int tmp=isExistChildrenChain(edge,numNode,child,parent);

	times[2]+=clock()-first_time;

	return tmp;	
}

int __buildBayesianNet(int numNode,int *edge,int numLine,bayesianNetScore **likely,double **likelies,int limitLikely,int start,int rank){
	int numLikely=0;
	int i,j;
	double tmp;
	if(rank==numLine){

		first_time=clock();

		tmp=likelihood(numNode,edge,likelies);
		
		times[0]+=clock()-first_time;
		first_time=clock();

		if(likely[limitLikely-1]->score > tmp){
			for(i=limitLikely-1;i>0;i--){
				if(likely[i]->score < tmp) break;
			}
			for(j=i+1;j<limitLikely;j++){
				likely[j]=likely[j-1];
			}
			likely[i]=makeBNScore(edge,numNode*numNode,tmp);
		}

		times[1]+=clock()-first_time;
		

		return 1;
	}else{
		for(i=start;i<numNode*numNode;i++){
			if(isEnableArc(edge,numNode,i)){

				edge[i]=1;
				numLikely+=__buildBayesianNet(numNode,edge,numLine,likely,likelies,limitLikely,i+1,rank+1);
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
long long getNumberOfBayesianNode(int limit){
	int i;
	long long sum;
	if(limit==0) return 1;
	for(sum=0,i=1;i<=limit;i++){
		long long nn=pow(-1,i+1)*combination(limit,i)*pow(2,i*(limit-i));
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
void localScore(baysianNode *node,int numNode,int currentNode,int numLine,Associate *cHist,unsigned short parentPattern,double *likelies,Associate **CPT,int rank,int start){
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
				if((parentPattern&((unsigned short)1<<j))!=0){
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
void __freeCPT(Associate **CPT,int numNode,int currentNode,unsigned short parentPattern,int numLine,int rank,int start){
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
bayesianNetwork* bayesianNetworkTrain(list_t *passengerList,bayesianNetworkOption *opt){
	Associate *frequency;
	int i,j,k;
	int *pattern;
	int *edge;
	cell_t *cur;
	bayesianNetwork *model;

	/* 分割表を作成 */
	frequency=makeLAA(passengerList->size,opt->numVariable);
	pattern=Malloc(int,opt->numVariable);
	for(cur=passengerList->first;cur;cur=cur->next){
		passenger *human=(passenger *)cur->data;
		pattern[0]=human->survived;
		setVariable(&pattern[1],human);

		int* cur=(int*)getLAA(frequency,pattern);
		if(cur==NULL){
			setLAA(frequency,pattern,int2Number(1));
		}else{
			(*cur)++;
		}
	}
	/* 条件付き頻度表を作成*/	
	Associate **cHist=Malloc(Associate*,opt->numVariable);
	for(i=0;i<opt->numVariable;i++){
		cHist[i]=makeLAA(passengerList->size,opt->numVariable);
		for(j=0;j<frequency->numKeys;j++){
			for(k=0;k<frequency->patternSize;k++){
				pattern[k]=(i==k) ? 0 : frequency->keys[j*frequency->patternSize+k];
			}
			int val=frequency->keys[j*frequency->patternSize+i];
			int *cur=(int*)getLAA(cHist[i],pattern);
			if(cur) cur[val]+=*((int*)frequency->array[j]);
			else{
				cur=Calloc(int,opt->node[i].numVariablePattern);
				cur[val]=*((int*)frequency->array[j]);
			}
			setLAA(cHist[i],pattern,cur);
		}
	}	
	/* ローカルスコアを計算*/
	Associate ***CPT=Malloc(Associate **,opt->numVariable);
	double **likelies=Malloc(double *,opt->numVariable);
	for(i=0;i<opt->numVariable;i++){
		likelies[i]=Malloc(double,pow(2,opt->numVariable));
		CPT[i]=Malloc(Associate *,pow(2,opt->numVariable));
		for(j=0;j<opt->numVariable;j++){
			localScore(&opt->node[i],opt->numVariable,i,j,cHist[i],0,likelies[i],CPT[i],0,0);
		}
	}

	model=Malloc(bayesianNetwork,1);
	model->CPT=CPT;
	
#ifdef BUILD_BAYESIAN_NET
	long long numLikely=0;
	long long  numCombinationPattern=getNumberOfBayesianNode(opt->numVariable);	//組み合わせ数を計算
	edge=Calloc(int,opt->numVariable*opt->numVariable);
	bayesianNetScore **arr=Malloc(bayesianNetScore*,NUM_LIKELIES);
	for(i=0;i<NUM_LIKELIES;i++){
		arr[i]=makeBNScore(edge,opt->numVariable*opt->numVariable,99999999);
	}
	for(i=0;numCombinationPattern!=numLikely;i++){
		numLikely+=__buildBayesianNet(opt->numVariable,edge,i,arr,likelies,NUM_LIKELIES,0,0);
#ifndef MAIN_PREDICTION
		printf("complite %015lld  / %lld \n",numLikely,numCombinationPattern);
#endif
	}
	
	model->models=arr;
	model->numModels=NUM_LIKELIES;

	free(edge);
#else
	model->numModels=0;
#endif
	
	for(i=0;i<opt->numVariable;i++){
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
	
	return model;
}
void freeBayesianNetwork(bayesianNetwork *model,bayesianNetworkOption *opt){
	int i;
	for(i=0;i<model->numModels;i++){
		free(model->models[i]->edge);
		free(model->models[i]);
	}
	free(model->models);
	freeCPT(model->CPT,opt->numVariable);
	free(model->CPT);
	
}
int bayesianNetworkIsEqual(int *x,passenger *human){
	return *x==human->survived;
}
double clossValidation(list_t *data,void* arg, void* (*train)(list_t *,void*),void* (*predict)(void*,void*,void*,double*),int isEqual(void*,void*),void (*freeModel)(void*,void*)){
	void *model;
	void *val;
	int hit;
	cell_t *cur;
	int num,i;

	num=(data->size < 100) ? data->size : 100;
	for(i=0,hit=0,cur=data->first;cur && i<num;cur=cur->next,i++){
		removeList(data,cur);
		model=train(data,arg);
		val=predict(model,arg,cur->data,NULL);
		freeModel(model,arg);
		if(isEqual(val,cur->data)) hit++;
		revivalList(data,cur);
		free(val);
	}
	return (double)hit/(double)num;
}

//相互情報量をベースとした情報量基準を計算
inline double mutualInformationMeasure(double *probabilityTarget,double *probabilityPredictor,double *jointProbability,int numA,int numB,int *base){
	double sum=0.0;
	int i,j;
	
	for(i=0;i<numA;i++){
		for(j=0;j<numB;j++){
			if((probabilityTarget[i]!=0 && probabilityPredictor[j])!=0 && jointProbability[i*numA+j]!=0) 
				sum+=jointProbability[i*numA+j]*log2(jointProbability[i*numA+j]/(probabilityTarget[i]*probabilityPredictor[j]));
		}
	}
	double mean=0.0;
	mean+=(double)((base[0])*(base[0]));
	for(i=0;i<numB-1;i++){
		mean+=(double)((base[i+1]-base[i])*(base[i+1]-base[i]));
	}
	mean+=(double)((80-base[i])*(80-base[i]));
	mean/=(double)numB;
	return 2*sum-numB-32*log(mean);
}
inline double func(int *passengers,int *sourcies,int numPassengers,int *base,int numSplits){
	int numTarget=2,i,j;
	double *jointProbability,*probabilityTarget,*probabilityPredictor;
	int *sumJoint,*sumTarget,*sumPredictor;

	sumJoint=Calloc(int,numTarget*(numSplits+1));
	sumTarget=Calloc(int,numTarget);
	sumPredictor=Calloc(int,(numSplits+1));

	jointProbability=Malloc(double,numTarget*(numSplits+1));
	probabilityTarget=Malloc(double,numTarget);
	probabilityPredictor=Malloc(double,(numSplits+1));
	for(j=0;j<numPassengers;j++){
		register int source=sourcies[j];
		for(i=0;i<numSplits;i++) if(source < base[i]) break;
		source=i;
		register int target=passengers[j];
		sumJoint[target*numTarget+source]++;
		sumTarget[target]++;
		sumPredictor[source]++;
	}
	for(i=0;i<numSplits;i++){
		probabilityPredictor[i]=(double)sumPredictor[i]/(double)numPassengers;
	}
	for(i=0;i<numTarget;i++){
		probabilityTarget[i]=(double)sumTarget[i]/(double)numPassengers;
		for(j=0;j<numSplits;j++){
			jointProbability[i*numTarget+j]=(double)sumJoint[i*numTarget+j]/(double)numPassengers;
		}
	}
	double r=mutualInformationMeasure(probabilityTarget,probabilityPredictor,jointProbability,numTarget,numSplits,base);
	free(jointProbability);
	free(probabilityTarget);
	free(probabilityPredictor);
	free(sumJoint);
	free(sumTarget);
	free(sumPredictor);
	return r;

}
void setBase(int *targets,int *sourcies,int numPassengers,int *base,int num,int max,int *maxBase,double *maxScore,int *maxNum,int rank,int start){
	int i;
	double score;
	if(rank==num){
		score=func(targets,sourcies,numPassengers,base,num);
		if(num > 1 && score > *maxScore){
			*maxScore=score;
			memcpy(maxBase,base,sizeof(int)*num);
			*maxNum=num;
		}
		/*printf("Score:%lf  ",score);
		for(i=0;i<num;i++){
			printf("%d ",base[i]);
		}puts("");*/
		//printf("num:%d %lf\n",num,score);
	}else{
		for(i=start;i<max;i++){
			base[rank]=i;
			setBase(targets,sourcies,numPassengers,base,num,max,maxBase,maxScore,maxNum,rank+1,i+1);
		}
	}
}
void *bayesianNetworkPredict(bayesianNetwork *model,bayesianNetworkOption *opt,passenger *human,double *likeli){
	int *pattern,i;
	unsigned short parentPattern=0;

	pattern=Malloc(int,opt->numVariable);
	setVariable(pattern,human);

	for(i=1;i<opt->numVariable;i++) parentPattern=parentPattern | (1<<i);
	double *prop=(double*)getLAA(model->CPT[0][parentPattern],pattern);
	if(prop==NULL){
		fprintf(stderr,"prop is NULL\n");
		return int2Number(0);
	}
	if(likeli){
		*likeli+=log(prop[human->survived]);
	}
	int serv=(prop[0] >= prop[1]) ? 0 : 1;
	free(pattern);
	return int2Number(serv);
}

void preProccess(list_t *passengerList,int *ageBorderLine,int numAgeBorderLine){
	int i;
	cell_t *cur;
	for(cur=passengerList->first;cur;cur=cur->next){
		passenger *human=(passenger *)cur->data;
		
		for(i=0;i<numAgeBorderLine;i++) if(human->age < ageBorderLine[i]) break;
		human->age=i;
		human->extention=(strcmp(human->name->honorific,"Don")==0 || strcmp(human->name->honorific,"Dr")==0 || strcmp(human->name->honorific,"Major")==0
			 || strcmp(human->name->honorific,"Jonkheer")==0 || strcmp(human->name->honorific,"Col")==0 || strcmp(human->name->honorific,"Capt")==0
			 || strcmp(human->name->honorific,"Countess")==0 || strcmp(human->name->honorific,"Rev")==0) ? 0 : 1;
		human->extention2=(*human->cabin==0) ? 2 : (human->cabin[strlen(human->cabin)-1]-'0')%2;
		//human->floor=(*human->cabin==0) ? 0 : 'G'-*human->cabin+1;
	}
}
void setVariable(int *pattern,passenger *human){
	
	pattern[0]=human->sex;
	pattern[1]=human->rank-1;
	pattern[2]=human->age;
	pattern[3]=human->prop1;
	pattern[4]=human->extention2;
	pattern[5]=human->alone;
	pattern[6]=human->floor;
}

#define NUM_VARIABLE 8

void ML_bayesianNetwork(list_t *passengerList,list_t *testPassengerList){
	baysianNode node[NUM_VARIABLE];
	bayesianNetworkOption opt;
	cell_t *cur;
	int i;
	
	
//量的変数の最適分割ラインを計算する
#ifdef CALCLATION_SPLITS
	#define NUM_SPLITS 5
	int maxBase[NUM_SPLITS];
	int base[NUM_SPLITS];
	int max,maxNum;
	max=0;
	double maxScore=-9999999.0;
	int *targets=Malloc(int,passengerList->size);
	int *sourcies=Malloc(int,passengerList->size);
	for(i=0,cur=passengerList->first;cur;cur=cur->next,i++){
		passenger *human=(passenger *)cur->data;
		targets[i]=human->survived;
		sourcies[i]=human->age;
		if(max < human->age) max=human->age;
	}
	for(i=0;i<NUM_SPLITS;i++){
		setBase(targets,sourcies,passengerList->size,base,i+1,max,maxBase,&maxScore,&maxNum,0,0);
	}
	
	printf("maxScore:%lf\n",maxScore);
	for(i=0;i<maxNum;i++){
		printf("%d ",maxBase[i]);
	}
	puts("");
	exit(0);
#else
	#define NUM_SPLITS 5
	int maxBase[NUM_SPLITS]={12,24,38,52,66};
	int maxNum=NUM_SPLITS;
#endif
	
	
	memset(&node,0,sizeof(node[0])*NUM_VARIABLE);
	node[0].name="SURVIVED";
	node[0].numVariablePattern=2;
	node[1].name="SEX";
	node[1].numVariablePattern=2;
	node[2].name="RANK";
	node[2].numVariablePattern=3;
	node[3].name="AGE";
	node[3].numVariablePattern=maxNum+1;
	node[4].name="FRIEND";
	node[4].numVariablePattern=3;
	node[5].name="POSITION";
	node[5].numVariablePattern=3;
	node[6].name="ALONE";
	node[6].numVariablePattern=2;
	node[7].name="FOOR";
	node[7].numVariablePattern=8;

	opt.numVariable=NUM_VARIABLE;
	opt.node=node;
	
	preProccess(passengerList,maxBase,maxNum);

//交差確認法でテストをする
#ifdef CLOSS_VALIDATION
	double cv=clossValidation(passengerList,&opt,(void* (*)(list_t*, void*))bayesianNetworkTrain,
				(void* (*)(void*,void*,void*,double*))bayesianNetworkPredict,(int (*)(void*, void*))bayesianNetworkIsEqual,(void (*)(void*, void*))freeBayesianNetwork);
	printf("CV:%lf\n",cv);
#endif

//トレーニングデータを使ってテストをする
#ifdef TEST
	
	bayesianNetwork *model=bayesianNetworkTrain(passengerList,&opt);
	for(i=0;i<model->numModels;i++){
		printf("SCORE:%lf\n",model->models[i]->score);
		printfBayesianNetwork(node,NUM_VARIABLE,model->models[i]->edge);
	}
	double likeli=0.0;
	int hit=0;
	for(cur=passengerList->first;cur;cur=cur->next){
		passenger* human=(passenger*)cur->data;
		
		int *val=(int*)bayesianNetworkPredict(model,&opt,(passenger*)human,&likeli);
		if(*val==human->survived) hit++;
	}
	printf("BIC:%lf  Accuary:%lf\n",bic(likeli,passengerList->size,NUM_VARIABLE-1),(double)hit/(double)passengerList->size);
	//freeBayesianNetwork(model,&opt);
#endif 

//本番予測
#ifdef MAIN_PREDICTION
	preProccess(testPassengerList,maxBase,maxNum);
	
	bayesianNetwork *modelMain=bayesianNetworkTrain(passengerList,&opt);
	puts("PassengerId,Survived");
	for(cur=testPassengerList->first;cur;cur=cur->next){
		passenger* human=(passenger*)cur->data;
		int *val=(int*)bayesianNetworkPredict(modelMain,&opt,(passenger*)human,NULL);
		printf("%d,%d\n",human->passengerId,*val);
	}
	freeBayesianNetwork(modelMain,&opt);
#else
	for(i=0;i<3;i++)
		printf("times[%d]:%lf\n",i,(double)times[i] / (double)CLOCKS_PER_SEC);
#endif


}
