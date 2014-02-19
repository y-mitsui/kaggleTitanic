#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>

#include "kaggleTitanic.h"
#include "list.h"

int isTest=1;

enum FIELD_NAME{
	SEX,
	TICKET,
	RANK,
	
};
enum RANGE_FIELD{
	AGE=0,
	FARE,
};
typedef struct{
	double *source;
	int n;
	double var;
	double mean;
	double bandWidth;
}kdeModel;

double kernel(double x){
	double var=0.1;
	double mean=0.0;
	return 1.0/sqrt(2.0*M_PI*var)*exp(-(x-mean)*(x-mean)/(2.0*var));
}
double variance(double *data,int num){
 int i;
 double mean=0,sum=0;
  
 for(i=0;i<num;i++){
  mean+=data[i];
 }
 mean/=num;
 
 for(i=0;i<num;i++){
  sum+=pow(mean-data[i],2);
 }
 
 return sum/num;
}
//プラグイン法によるバンド幅
//data:データ配列
//num:配列の要素数
//戻り値：バンド幅
double calcBandWidth(double *data,int num){
 return 1.06*sqrt(variance(data,num))*pow(num,-0.2);
}

kdeModel *makeKernelDentisyEstimationModel(double *source,int n,double bandWidth){
	kdeModel *model;
	model=Malloc(kdeModel,1);
	double mean,var,sum;
	double *data;
	int i;

	data=Malloc(double,n);
	memcpy(data,source,sizeof(double)*n);
	sum=0.0;
	for(i=0;i<n;i++){
		sum+=data[i];
	}
	mean=sum/(double)n;
	sum=0.0;
	for(i=0;i<n;i++){
		sum+=(data[i]-mean)*(data[i]-mean);
	}
	var=sum/n;
	var=sqrt(var);
	for(i=0;i<n;i++){
		data[i]=(data[i]-mean)/var;
	}
	
	model->source=data;
	model->var=var;
	model->mean=mean;
	model->bandWidth=(bandWidth >= 0.0) ?  bandWidth : calcBandWidth(data,n);
	model->n=n;
	return model;
}
double kernelDentisyEstimation(double *source,int n,double target,double bandWidth){
	double sum=0.0;
	int i;

	for(i=0;i<n;i++){
		sum+=(1/(n*bandWidth))*kernel((target-source[i])/bandWidth);
	}
	return sum;
}
double getRangeProp(kdeModel *model,double from,double to){
	double t,sum2;
	from=(from-model->mean)/model->var;
	to=(to-model->mean)/model->var;

	for(t=from,sum2=0.0;t<to;t+=0.001){
		double diff=kernelDentisyEstimation(model->source,model->n,t,model->bandWidth)
				+kernelDentisyEstimation(model->source,model->n,t+0.001,model->bandWidth);
		diff/=2.0;
		diff=0.001*diff;
		sum2+=diff;
	}
	return sum2;
}
void message(const char* format, ...){
	va_list args;
	va_start( args, format );
	if(!isTest) vprintf(format,args);
}
void debug(const char* format, ...){
	va_list args;
	va_start( args, format );
	if(isTest) vprintf(format,args);	
}

double dmax(double x,double y){
	return (x > y) ? x : y;
}

void convertGrouping(double *data,int n,int numGroup){
	int i;
	double max,division;
	
	for(max=0.0,i=0;i<n;i++){
		if(max < data[i]) max=data[i];
	}
	division=(max+1)/numGroup;
	for(i=0;i<n;i++){
		data[i]=(double)((int)(data[i]/division));
	}
}
void total(double *data,int n,int *hinst){
	int i;
	for(i=0;i<n;i++){
		hinst[(int)data[i]]++;
	}
}
double sendRangeData(int n,void *data,void *arg){
	list_t *passengerList=(list_t *)arg;
	cell_t *cur;
	int allDead=0,hasFriend=0;
	passenger *human=(passenger *)data;
	switch(n){
	case AGE:  return (double)human->age; break;
	case FARE: return (double)human->fare; break;
	default: exit(1); break;
	}
}
double sendData(int n,void *data,void *arg){
	list_t *passengerList=(list_t *)arg;
	cell_t *cur;
	int allDead=0,hasFriend=0;
	passenger *human=(passenger *)data;
	switch(n){
	case SEX: return (double)human->sex; break;
	case TICKET: 
		for(cur=passengerList->first;cur;cur=cur->next){
			passenger *human2=(passenger *)cur->data;
			if(human!=human2 && strcmp(human->ticketNo,human2->ticketNo)==0){
				hasFriend++;
				if(human2->survived==-1)
					allDead++;
			}
		}
		return (hasFriend >2 && hasFriend==allDead) ? 1 : 0;
		break;	
	default: exit(1); break;
	}
}
double** groupsPropability(list_t *motherList,int numDataTypes,int *numGrouping,double (*getData)(int,void*,void*),void *getDataArg,double **data){
	cell_t *cur;
	double *groupsData;
	int *hinst;
	double **propabirities;
	int i,j;

	groupsData=Malloc(double,motherList->size*numDataTypes);
	for(i=0,cur=motherList->first;cur;cur=cur->next,i++){
		for(j=0;j<numDataTypes;j++){
			groupsData[j*motherList->size+i]=getData(j,cur->data,getDataArg);
		}
	}
	propabirities=Malloc(double*,numDataTypes);
	for(i=0;i<numDataTypes;i++){
		convertGrouping(&groupsData[i*motherList->size],motherList->size,numGrouping[i]);
		hinst=Calloc(int,numGrouping[i]);
		total(&groupsData[i*motherList->size],motherList->size,hinst);
		propabirities[i]=Malloc(double,numGrouping[i]);
		for(j=0;j<numGrouping[i];j++)
			propabirities[i][j]=(double)hinst[j]/(double)motherList->size;
		free(hinst);
	}
	if(data) *data=groupsData;
	else free(groupsData);

	return propabirities;
}

kdeModel **rangePropability(list_t *motherList,int numDataTypes,double *bandWidth,double (*getData)(int,void*,void*),void *getDataArg,double **data){
	double *rangeData;
	kdeModel **models;
	int i,j;
	cell_t *cur;

	rangeData=Malloc(double,motherList->size*numDataTypes);
	for(i=0,cur=motherList->first;cur;cur=cur->next,i++){
		for(j=0;j<numDataTypes;j++){
			rangeData[j*motherList->size+i]=getData(j,cur->data,getDataArg);
		}
	}
	models=Malloc(kdeModel*,numDataTypes);
	for(i=0;i<numDataTypes;i++){
		models[i]=makeKernelDentisyEstimationModel(&rangeData[i*motherList->size],motherList->size,bandWidth[i]);
	}
	if(data) *data=rangeData;
	else free(rangeData);
	return models;
}
#define NUM_DATA_TYPE 2
#define NUM_RANGE_DATA_TYPE 2
void naiveBayesClassifier2(list_t *passengerList,list_t *testPassengerList){
	
}
void naiveBayesClassifier(list_t *passengerList,list_t *testPassengerList){
	double **allPropabirity,**deadPropabirity,**alivePropabirity,*data;
	double deadUser,aliveUser,deadProp;
	list_t deadPassenger,alivePassenger;
	cell_t *cur;
	int i,j,hit,serv;
	int groupSize[NUM_DATA_TYPE];
	kdeModel **deadRangePropabirity,**aliveRangePropabirity,**allRangePropabirity;
	double rangeProp,deadRangeProp,aliveRangeProp,*rangeData;
	/* 範囲確率設定*/
	double bandWidth[NUM_RANGE_DATA_TYPE]={-1,-1};
	double range[NUM_RANGE_DATA_TYPE]={1.0,10.0};

	groupSize[SEX]=2;	
	groupSize[TICKET]=2;

	memset(&deadPassenger,0,sizeof(list_t));
	memset(&alivePassenger,0,sizeof(list_t));
	
	for(cur=passengerList->first;cur;cur=cur->next){
		passenger *human=(passenger *)cur->data;
		if(human->survived==-1)
			addList(&deadPassenger,human);
		else
			addList(&alivePassenger,human);
	}

	allPropabirity=groupsPropability(passengerList,NUM_DATA_TYPE,groupSize,sendData,passengerList,&data);
	deadPropabirity=groupsPropability(&deadPassenger,NUM_DATA_TYPE,groupSize,sendData,&passengerList,NULL);
	alivePropabirity=groupsPropability(&alivePassenger,NUM_DATA_TYPE,groupSize,sendData,&passengerList,NULL);
	
	deadProp=(double)deadPassenger.size/(double)passengerList->size;

	allRangePropabirity=rangePropability(passengerList,NUM_RANGE_DATA_TYPE,bandWidth,sendRangeData,&passengerList,&rangeData);
	deadRangePropabirity=rangePropability(&deadPassenger,NUM_RANGE_DATA_TYPE,bandWidth,sendRangeData,&passengerList,NULL);
	aliveRangePropabirity=rangePropability(&alivePassenger,NUM_RANGE_DATA_TYPE,bandWidth,sendRangeData,&passengerList,NULL);




	list_t *accuary=(isTest) ? passengerList : testPassengerList;
	message("PassengerId,Survived\n");
	for(hit=0,i=0,cur=accuary->first;cur;cur=cur->next,i++){
		deadUser=deadProp;
		aliveUser=1.0-deadProp;
		passenger *human=(passenger *)cur->data;
		for(j=0;j<NUM_DATA_TYPE;j++){
			int d=(int)data[j*passengerList->size+i];
			debug("[%d] parameter:%d befor:%lf after:%lf\n",j,d,deadUser,deadUser*deadPropabirity[j][d]);
			debug("[%d] parameter:%d befor:%lf after:%lf\n",j,d,aliveUser,aliveUser*alivePropabirity[j][d]);
			deadUser*=deadPropabirity[j][d];
			aliveUser*=alivePropabirity[j][d];
		}

		for(j=0;j<NUM_RANGE_DATA_TYPE;j++){
			int d=(int)rangeData[j*passengerList->size+i];

			deadRangeProp=dmax(getRangeProp(deadRangePropabirity[j],(double)d-range[d],(double)d+range[d]),0.001);
			aliveRangeProp=dmax(getRangeProp(aliveRangePropabirity[j],(double)d-range[d],(double)d+range[d]),0.001);
			debug("[%d] parameter:%d befor:%lf after:%lf\n",j,d,deadUser,deadUser*deadRangeProp);
			debug("[%d] parameter:%d befor:%lf after:%lf\n",j,d,aliveUser,aliveUser*aliveRangeProp);
			deadUser*=deadRangeProp;
			aliveUser*=aliveRangeProp;
		}
		debug("--------------\n");

		serv=(deadUser>aliveUser) ? ((isTest) ? -1 : 0) : 1;
		message("%d,%d\n",human->passengerId,serv);
		if(serv==human->survived) hit++;
	}
	debug("Accuary=%lf\n",(double)hit/(double)passengerList->size);

	for(i=0;i<NUM_DATA_TYPE;i++){
		free(allPropabirity[i]);
		free(deadPropabirity[i]);
	}
	free(allPropabirity);
	free(deadPropabirity);
	free(data);
	free(rangeData);
	allDeleteList(deadPassenger.first);
	allDeleteList(alivePassenger.first);
}


/* 女性率*/
double womanPropabirity(list_t *passengerList){
	int sum=0;
	cell_t *cur;
	for(cur=passengerList->first;cur;cur=cur->next){
		passenger *human=(passenger *)cur->data;
		if(human->sex==1) sum++;
	}
	return (double)sum/(double)passengerList->size;
}
void naiveBeise(list_t *passengerList,list_t *testPassengerList){
	return naiveBayesClassifier(passengerList,testPassengerList);
	list_t deadPassenger;
	double deadProp,womanProp,womanDeadProp,isDead;
	double *fare,*deadFare;
	double *age,*deadAge,sum2;
	double isAllDeadProp,isAllServProp,deadIsAllDeadProp,deadIsAllServProp;
	int sum,serv,hit,i;
	double rankProp[3],deadRankProp[3];
	int rankCt[3];
	int hasFriendCt,isAllDeadCt,isAllServCt;
	cell_t *cur,*cur2;
	double bandWidth=0.97;
	int parFriendDead,parHasFriend;
	double parFriendDeadProp,deadParFriendDeadProp;

	hasFriendCt=isAllDeadCt=isAllServCt=0;
	//parDead=parHasFriend=0;
	for(cur=passengerList->first;cur;cur=cur->next){
		passenger *human=(passenger *)cur->data;
		int isDead=0;
		int isServ=0;
		for(cur2=passengerList->first;cur2;cur2=cur2->next){
			passenger *human2=(passenger *)cur2->data;
			if(human!=human2 && strcmp(human->ticketNo,human2->ticketNo)==0){
				human->ufriend.hasFriend++;
				if(human2->survived==-1) isDead=1;
				else isServ=1;
			}
		}
		if(human->ufriend.hasFriend>0){
			if(isDead) parFriendDead++;
			parHasFriend++;
		}
		if(human->ufriend.hasFriend>1){
			if(isServ==0){ human->ufriend.isAllDead=1; isAllDeadCt++; }
			if(isDead==0){ human->ufriend.isAllServ=1; isAllServCt++; }
			hasFriendCt++;
		}
	}
	isAllDeadProp=(double)isAllDeadCt/(double)hasFriendCt;
	isAllServProp=(double)isAllServCt/(double)hasFriendCt;
	//printf("%d %d %d %lf %lf\n",isAllDeadCt,isAllServCt,hasFriendCt,isAllDeadProp,isAllServProp);
	parFriendDeadProp=(double)parFriendDead/deadParFriendDeadProp;

	memset(&deadPassenger,0,sizeof(list_t));
	age=Malloc(double,passengerList->size);
	fare=Malloc(double,passengerList->size);
	rankCt[0]=rankCt[1]=rankCt[2]=0;
	for(i=0,cur=passengerList->first;cur;cur=cur->next,i++){
		passenger *human=(passenger *)cur->data;
		if(human->survived==-1){
			addList(&deadPassenger,human);
		}
		age[i]=(double)human->age;
		fare[i]=(double)human->fare;
		rankCt[human->rank-1]++;
	}
	kdeModel *ageModel=makeKernelDentisyEstimationModel(age,i,bandWidth);
	kdeModel *fareModel=makeKernelDentisyEstimationModel(fare,i,0.45);

	womanProp=womanPropabirity(passengerList);
	
	deadProp=(double)deadPassenger.size/(double)passengerList->size;
	rankProp[0]=(double)rankCt[0]/(double)passengerList->size;
	rankProp[1]=(double)rankCt[1]/(double)passengerList->size;
	rankProp[2]=(double)rankCt[2]/(double)passengerList->size;

	deadAge=Malloc(double,deadPassenger.size);
	deadFare=Malloc(double,deadPassenger.size);
	hasFriendCt=isAllDeadCt=isAllServCt=0;
	rankCt[0]=rankCt[1]=rankCt[2]=0;
	for(i=0,sum=0,cur=deadPassenger.first;cur;cur=cur->next,i++){
		passenger *human=(passenger *)cur->data;
		if(human->sex==1) sum++;
		deadAge[i]=human->age;
		deadFare[i]=human->fare;
		if(human->ufriend.hasFriend>1){
			if(human->ufriend.isAllDead==1){ isAllDeadCt++; }
			else if(human->ufriend.isAllServ==1){ isAllServCt++; }
			hasFriendCt++;
		}
		rankCt[human->rank-1]++;
	}
	deadRankProp[0]=(double)rankCt[0]/(double)deadPassenger.size;
	deadRankProp[1]=(double)rankCt[1]/(double)deadPassenger.size;
	deadRankProp[2]=(double)rankCt[2]/(double)deadPassenger.size;
	deadIsAllDeadProp=(double)isAllDeadCt/(double)hasFriendCt;
	deadIsAllServProp=(double)isAllServCt/(double)hasFriendCt;

	kdeModel *deadAgeModel=makeKernelDentisyEstimationModel(deadAge,i,bandWidth);
	kdeModel *deadFareModel=makeKernelDentisyEstimationModel(deadFare,i,0.45);

	sum2=0.0;
	/*for(i=0;i<512;i++){
		//sum2+=getAgeProp(deadFareModel,(double)i-0.5,(double)i+0.5);
		printf("%d %lf  %lf\n",i,getAgeProp(deadFareModel,(double)i-0.5,(double)i+0.5),getAgeProp(fareModel,(double)i-0.5,(double)i+0.5));	
	}*/
	//printf("sum2:%lf\n",sum2);

	/* 死者のうち女性である確立*/
	womanDeadProp=(double)sum/(double)deadPassenger.size;

	hit=0;
	//puts("PassengerId,Survived");	
	for(i=0,cur=passengerList->first;cur;cur=cur->next,i++){
		isDead=deadProp;
		passenger *human=(passenger *)cur->data;
		if(human->sex==1) isDead=(womanDeadProp*isDead)/womanProp;
		else isDead=((1.0-womanDeadProp)*isDead)/(1.0-womanProp);


		isDead=isDead*(deadRankProp[human->rank-1]/rankProp[human->rank-1]);

		/*if(human->ufriend.hasFriend>1){
			if(human->ufriend.isAllDead==1){
				//puts("friends all dead");
				isDead=(deadIsAllDeadProp*isDead)/isAllDeadProp;
			}else if(human->ufriend.isAllServ==1){
				//puts("friends all servied");
				isDead=(deadIsAllServProp*isDead)/isAllServProp;
			}
		}
		double ageProp=dmax(getAgeProp(ageModel,(double)human->age-0.5,(double)human->age+0.5),0.001);
		double deadAgeProp=dmax(getAgeProp(deadAgeModel,(double)human->age-0.5,(double)human->age+0.5),0.001);
		isDead=(deadAgeProp*isDead)/ageProp;*/

		/*double fareProp=dmax(getAgeProp(fareModel,(double)human->fare-0.5,(double)human->fare+0.5),0.0001);
		double deadFareProp=dmax(getAgeProp(deadFareModel,(double)human->fare-0.5,(double)human->fare+0.5),0.0001);
		isDead=(deadFareProp*isDead)/fareProp;*/

		

		serv=(isDead>0.5) ? -1 : 1;
		//printf("%d,%d\n",human->passengerId,serv);
		if(serv==human->survived) hit++;
	}
	printf("Accuary2=%lf\n",(double)hit/(double)passengerList->size);
}
