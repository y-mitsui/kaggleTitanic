#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <dirent.h>
#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>
#include "kaggleTitanic.h"
#include "error.h"
#include "list.h"
#include "mystring.h"
#include "svm.h"

void naiveBeise(list_t *passengerList,list_t *testPassengerList);
void ML_bayesianNetwork(list_t *passengerList,list_t *testPassengerList);

int famiryNo=1;

nameData *parseName(char *name){
	char *names[3];
	char tmp[1024];
	int tmpNum=0,namesNum=0,i;
	nameData *result=Malloc(nameData,1);

	int len=strlen(name);
	names[0]=names[1]=names[2]=NULL;
	
	for(i=0;i<len;i++){
		switch(name[i]){
		case ' ': break;
		case ',':
		case '.':
			tmp[tmpNum]='\0';
			tmpNum=0;
			names[namesNum++]=strdup(tmp);
			break;
		case '(': goto LOOP_OUT;
		default:
			tmp[tmpNum++]=name[i];
		}
	}
LOOP_OUT:
	tmp[tmpNum]='\0';
	names[namesNum++]=strdup(tmp);
	result->first=names[2];
	result->honorific=names[1];
	result->second=names[1];
	return result;
}
void setFamiryNo(list_t *passengerList){
	cell_t *cur,*cur2;
	for(cur=passengerList->first;cur;cur=cur->next){
		passenger *src=(passenger *)cur->data;
		for(cur2=passengerList->first;cur2;cur2=cur2->next){
			passenger *target=(passenger *)cur2->data;
			if(!strcmp(src->name->second,target->name->second)
			&& strcmp(src->name->honorific,target->name->honorific)){
				src->famiryNo=target->famiryNo;
			}
			if(src->famiryNo==0) src->famiryNo=famiryNo++;
		}
	}
}

double sameTicketScale(list_t *passengerList){
	int dead,deadMe=0,mother=0;
	cell_t *cur,*cur2;
	int flag;

	for(cur=passengerList->first;cur;cur=cur->next){
		passenger *source=(passenger *)cur->data;
		dead=0;
		flag=0;
		for(cur2=passengerList->first;cur2;cur2=cur2->next){
			passenger *target=(passenger *)cur2->data;
			
			if(target->passengerId!=source->passengerId
			&& !strcmp(target->ticketNo,source->ticketNo)){
				flag++;
				if(target->survived==1){
					dead++;
					
				}
			}
		}
		if(flag && flag==dead) mother++;
		if(flag && flag==dead) source->prop1=2;	//‘Sˆõ•‚©‚Á‚½ê‡
		else if(flag && dead==0) source->prop1=0;	//‘SˆõŽ€‚ñ‚¾ê‡
		else source->prop1=1;
		/*
		if(flag && flag==dead) mother++;
		if(flag && flag==dead && source->survived==1){
			deadMe++;
			source->prop1=1000;
		}else if(flag && dead==0)
			source->prop1=0;
		else
			source->prop1=500;*/
		/*
		if(flag>0 && dead>0){
			if(source->survived==1){
				deadMe++;
			}
		}*/
	}
	
	return (double)deadMe/(double)mother;
}
double sameFamiryScale(list_t *passengerList){
	int dead,deadMe=0,mother=0;
	cell_t *cur,*cur2;
	int flag;

	for(cur=passengerList->first;cur;cur=cur->next){
		passenger *source=(passenger *)cur->data;
		dead=0;
		flag=0;
		for(cur2=passengerList->first;cur2;cur2=cur2->next){
			passenger *target=(passenger *)cur2->data;
			
			if(target->passengerId!=source->passengerId
			&& target->famiryNo==source->famiryNo){
				flag++;
				if(target->survived==-1){
					dead++;
				}
			}
		}
		if(flag) mother++;
		if(flag && flag==dead){			
			if(source->survived==-1){
				deadMe++;
			}
		}
	}
	
	return (double)deadMe/(double)mother;
}
void scale(struct svm_node **x,int n,struct svm_node **x2,int n2,int w){
	int i,j;
	double *min,*max;

	min=Malloc(double,w);
	max=Malloc(double,w);
	
	for(i=0;i<w;i++){
		min[i]=-100000000;
		max[i]=100000000;
	}
	for(i=0;i<n;i++){
		for(j=0;j<w;j++){
			if(min[j] > x[i][j].value) min[j]=x[i][j].value;
			if(max[j] < x[i][j].value) max[j]=x[i][j].value;
		}
	}
	for(i=0;i<n2;i++){
		for(j=0;j<w;j++){
			if(min[j] > x2[i][j].value) min[j]=x2[i][j].value;
			if(max[j] < x2[i][j].value) max[j]=x2[i][j].value;
		}
	}
	for(i=0;i<n;i++){
		for(j=0;j<w;j++){
			x[i][j].value=(x[i][j].value-min[j])/(max[j]-min[j])*2.0-1.0;
		}
	}
	for(i=0;i<n2;i++){
		for(j=0;j<w;j++){
			x2[i][j].value=(x2[i][j].value-min[j])/(max[j]-min[j])*2.0-1.0;
		}
	}
}
double cv(struct svm_parameter *param,struct svm_problem *prob){
	struct svm_problem prob2;
	struct svm_model *model;
	int i,j;
	double predict_label,target_label;
	double sum,sum2,mean;
	prob2.l=prob->l;
	prob2.y = Malloc(double,prob->l);
	prob2.x = Malloc(struct svm_node *,prob->l);
	double error = 0;
	double sump = 0, sumt = 0, sumpp = 0, sumtt = 0, sumpt = 0;
	int total = 0;
	prob2.l--;

	int hit=0;

	for(sum=0.0,i=0;i<prob->l;i++){
		sum+=prob->y[i];
	}
	mean=sum/prob->l;
	printf("mean:%lf\n",mean);
	for(sum=0.0,sum2=0.0,i=0;i<prob->l;i++){
		for(j=0;j<i;j++){
			prob2.x[j]=prob->x[j];
			prob2.y[j]=prob->y[j];
		}
		for(j=i;j<prob2.l;j++){
			prob2.x[j]=prob->x[j+1];
			prob2.y[j]=prob->y[j+1];
		}
		model = svm_train(&prob2,param);		
		predict_label = svm_predict(model,prob->x[i]);
		//predict_label=mean;
		svm_free_and_destroy_model(&model);
		sum+=(predict_label-prob->y[i])*(predict_label-prob->y[i]);
		sum2+=(prob->y[i]-mean)*(prob->y[i]-mean);

		if(predict_label==prob->y[i]) hit++;
		printf("predict_label:%lf prob->y[i]:%lf\n",predict_label,prob->y[i]);
		target_label=prob->y[i];
		error += (predict_label-target_label)*(predict_label-target_label);
		sump += predict_label;
		sumt += target_label;
		sumpp += predict_label*predict_label;
		sumtt += target_label*target_label;
		sumpt += predict_label*target_label;
		++total;
	}
	return error/(double)total;
	return ((total*sumpt-sump*sumt)*(total*sumpt-sump*sumt))/ ((total*sumpp-sump*sump)*(total*sumtt-sumt*sumt));	
	return 1.0-sum/sum2;
	return (double)hit/(double)prob->l;
	
}
double getNameNo(char *name){
	if(!strcmp("Mrs",name)) return 3.0;
	if(!strcmp("Mr",name)) return 2.0;
	if(!strcmp("Miss",name)) return 1.0;
	if(!strcmp("Master",name)) return 0.0;
	return 1.5;
}
void fillAge(list_t *passengerList,list_t *testPassengerList){
	list_t inAge,notAge;
	struct svm_parameter param;		// set by parse_command_line
	struct svm_problem prob,prob2;		// set by read_problem
	struct svm_model *model;
	struct svm_node *x_space;
	int max_index=5;
	int i;
	cell_t *cur;
	passenger *human;

	// default values
	param.svm_type = NU_SVR;
	param.kernel_type = LINEAR;
	param.degree = 3;
	param.gamma = 0;	// 1/num_features
	param.coef0 = 0;
	param.nu = 0.5;
	param.cache_size = 0;
	param.C = 1;
	param.eps = 1e-2;
	param.p = 0.1;
	param.shrinking = 1;
	param.probability = 0;
	param.nr_weight = 0;
	param.weight_label = NULL;
	param.weight = NULL;

	memset(&inAge,0,sizeof(list_t));
	memset(&notAge,0,sizeof(list_t));
	for(cur=passengerList->first;cur;cur=cur->next){
		human=(passenger *)cur->data;
		if(human->age==-1) addList(&notAge,human);
		else addList(&inAge,human);
	}
	for(cur=testPassengerList->first;cur;cur=cur->next){
		human=(passenger *)cur->data;
		if(human->age==-1) addList(&notAge,human);
		else addList(&inAge,human);
	}
	
	prob.l=inAge.size;
	prob.y = Malloc(double,prob.l);
	prob.x = Malloc(struct svm_node *,prob.l);
	for(i=0,cur=inAge.first;cur;cur=cur->next,i++){
		x_space = Malloc(struct svm_node,6);
		passenger *human=(passenger *)cur->data;
		x_space[0].index=1;
		x_space[0].value=(double)human->sex;
		x_space[1].index=2;
		x_space[1].value=(double)human->rank;
		x_space[2].index=3;
		x_space[2].value=(double)human->fare;
		x_space[3].index=4;
		x_space[3].value=(double)human->prop1;
		x_space[4].index=5;
		x_space[4].value=getNameNo(human->name->honorific);
		
		x_space[max_index].index = -1;
		prob.x[i]=x_space;
		prob.y[i]=(double)human->age;
	}
	prob2.l=notAge.size;
	prob2.x = Malloc(struct svm_node *,prob2.l);
	passenger **tmp= Malloc(passenger *,prob2.l);
	for(i=0,cur=notAge.first;cur;cur=cur->next,i++){
		x_space = Malloc(struct svm_node,6);
		passenger *human=(passenger *)cur->data;
		x_space[0].index=1;
		x_space[0].value=(double)human->sex;
		x_space[1].index=2;
		x_space[1].value=(double)human->rank;
		x_space[2].index=3;
		x_space[2].value=(double)human->fare;
		x_space[3].index=4;
		x_space[3].value=(double)human->prop1;
		x_space[4].index=5;
		x_space[4].value=getNameNo(human->name->honorific);
		/*x_space[5].index=6;
		x_space[5].value=(human->survived==-1) ? 1.0:0.0;*/
		x_space[max_index].index = -1;
		prob2.x[i]=x_space;
		tmp[i]=human;
	}

	param.gamma=1.0/(double)max_index;
	scale(prob.x,prob.l,prob2.x,prob2.l,max_index);
	//printf("%lf\n",cv(&param,&prob));
	model = svm_train(&prob,&param);
	
	for(i=0;i<prob2.l;i++){
		tmp[i]->age=(int)svm_predict(model,prob2.x[i]);
	}
	svm_free_model_content(model);
	svm_free_and_destroy_model(&model);
	svm_destroy_param(&param);

	
	
	for(i=0,cur=inAge.first;cur;cur=cur->next,i++){
		free(prob.x[i]);
	}
	for(i=0,cur=notAge.first;cur;cur=cur->next,i++){
		free(prob2.x[i]);
	}
	free(prob.x);
	free(prob.y);
	free(prob2.x);

}
void print_null(const char *s) { s++; }
int main(void){
	struct svm_parameter param;		// set by parse_command_line
	struct svm_problem prob,prob2;		// set by read_problem
	struct svm_model *model;
	struct svm_node *x_space;
	double predict_label;
	list_t passengerList,testPassengerList;
	FILE *fp;
	cell_t *cur;
	int i,j,hit;

	// default values
	param.svm_type = C_SVC;
	param.kernel_type = RBF;
	param.degree = 3;
	param.gamma = 0;	// 1/num_features
	param.coef0 = 0;
	param.nu = 0.5;
	param.cache_size = 100;
	param.C = 1;
	param.eps = 1e-3;
	param.p = 0.1;
	param.shrinking = 1;
	param.probability = 0;
	param.nr_weight = 0;
	param.weight_label = NULL;
	param.weight = NULL;

	svm_set_print_string_function(print_null);

	memset(&passengerList,0,sizeof(passengerList));
	memset(&testPassengerList,0,sizeof(passengerList));

	char buf[1024];
	char **token;

	dfopen(fp,"data/train.csv","r",exit(EXIT_FAILURE));
	fgets(buf,sizeof(buf),fp);

	while(fgets(buf,sizeof(buf),fp)){		
		int n=explode(&token,buf,',');
		passenger *human=Calloc(passenger,1);
		human->passengerId=atoi(token[0]);
		human->survived=atoi(token[1]);
		human->sex=(strcmp("male",token[4])==0) ? 0 : 1;
		human->age=(*token[5]=='\0') ? -1 : atoi(token[5]);
		human->rank=atoi(token[2]);
		human->fare=atoi(token[9]);
		human->name=parseName(token[3]);
		human->ticketNo=strdup(token[8]);
		human->cabin=strdup(token[10]);
		addList(&passengerList,human);
		for(i=0;i<n;i++){
			free(token[i]);
		}
		free(token);
	}
	fclose(fp);
	//setFamiryNo(&passengerList);
	sameTicketScale(&passengerList);
	

	dfopen(fp,"data/test.csv","r",exit(EXIT_FAILURE));
	fgets(buf,sizeof(buf),fp);
	while(fgets(buf,sizeof(buf),fp)){
		explode(&token,strdup(buf),',');
		passenger *human=Calloc(passenger,1);
		human->passengerId=atoi(token[0]);
		human->sex=(strcmp("male",token[3])==0) ? 0 : 1;
		human->age=(*token[4]=='\0') ? -1 : 100-atoi(token[5]);
		human->rank=atoi(token[1]);
		human->fare=atoi(token[8]);
		human->name=parseName(token[2]);
		human->ticketNo=strdup(token[7]);
		human->cabin=strdup(token[9]);
		addList(&testPassengerList,human);
	}
	fclose(fp);

	//setFamiryNo(&testPassengerList);
	sameTicketScale(&testPassengerList);

	fillAge(&passengerList,&testPassengerList);

	ML_bayesianNetwork(&passengerList,&testPassengerList);
	exit(0);

#define MAX_INDEX 5
	prob.l=passengerList.size;
	prob.x = Malloc(struct svm_node *,prob.l);
	prob.y = Malloc(double,prob.l);
	for(i=0,cur=passengerList.first;cur;cur=cur->next,i++){
		x_space = Malloc(struct svm_node,6);
		passenger *human=(passenger *)cur->data;
		x_space[0].index=1;
		x_space[0].value=(double)human->age;
		x_space[1].index=2;
		x_space[1].value=(double)human->sex;
		x_space[2].index=3;
		x_space[2].value=(double)human->rank;
		x_space[3].index=4;
		x_space[3].value=(double)human->fare;
		x_space[4].index=5;
		x_space[4].value=(double)human->prop1;
		x_space[MAX_INDEX].index = -1;
		prob.x[i]=x_space;
		prob.y[i]=human->survived;
	}
	param.gamma=1.0/(double)MAX_INDEX;

	prob2.l=testPassengerList.size;
	passenger **tmp = Malloc(passenger *,prob2.l);
	prob2.x = Malloc(struct svm_node *,prob2.l);
	
	for(i=0,cur=testPassengerList.first;cur;cur=cur->next,i++){
		x_space = Malloc(struct svm_node,6);
		passenger *human=(passenger *)cur->data;
		x_space[0].index=1;
		x_space[0].value=(double)human->age;
		x_space[1].index=2;
		x_space[1].value=(double)human->sex;
		x_space[2].index=3;
		x_space[2].value=(double)human->rank;
		x_space[3].index=4;
		x_space[3].value=(double)human->fare;
		x_space[4].index=5;
		x_space[4].value=(double)human->prop1;
		x_space[MAX_INDEX].index = -1;
		prob2.x[i]=x_space;
		tmp[i]=human;
	}
	
	naiveBeise(&passengerList,&testPassengerList);
	exit(0);

	scale(prob.x,prob.l,prob2.x,prob2.l,MAX_INDEX);

	

	prob2.l=prob.l;
	prob2.y = Malloc(double,prob.l);
	prob2.x = Malloc(struct svm_node *,prob.l);

	prob2.l--;

	hit=0;
	for(i=0;i<prob.l;i++){
		for(j=0;j<i;j++){
			prob2.x[j]=prob.x[j];
			prob2.y[j]=prob.y[j];
		}
		for(j=i;j<prob2.l;j++){
			prob2.x[j]=prob.x[j+1];
			prob2.y[j]=prob.y[j+1];
		}
		model = svm_train(&prob2,&param);		
		predict_label = svm_predict(model,prob.x[i]);
		svm_free_and_destroy_model(&model);
		
		if(predict_label==prob.y[i]) hit++;
		//printf("predict_label:%lf\n",predict_label);
	}
	printf("%lf\n",(double)hit/(double)prob.l);
	exit(1);



	model = svm_train(&prob,&param);
	puts("PassengerId,Survived");	
	for(i=0;i<prob2.l;i++){
		predict_label = svm_predict(model,prob2.x[i]);
		//printf("%d\n",(int)predict_label);
		printf("%d,%d\n",tmp[i]->passengerId,(predict_label==-1.0) ? (int)0 : (int)1);
	}
	
	return 0;
}

