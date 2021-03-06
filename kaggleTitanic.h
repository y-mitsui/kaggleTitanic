
typedef struct {
	int hasFriend;
	int isAllDead;
	int isAllServ;
}friendInfo;
typedef struct {
	char *first;
	char *second;
	char *honorific;
}nameData;

typedef struct {
	int survived;
	int passengerId;
	nameData *name;
	int sex;
	int age;
	int rank;
	int fare;
	char* cabin;
	int famiryNo;
	char *ticketNo;
	char *proof;
	int prop1;
	int extention;
	int extention2;
	int alone;
	int floor;
	friendInfo ufriend;
}passenger;

#define Malloc(type,n) (type *)malloc((n)*sizeof(type))
#define Calloc(type,n) (type *)calloc(1,(n)*sizeof(type))
#define LAMBDA(rettype, ARG_LIST, BODY)          \
({                                               \
   rettype __lambda_funcion__ ARG_LIST { BODY; } \
   __lambda_funcion__;                           \
})

#define MAX(dataType,da,db)				\
({							\
	dataType __max__ (dataType a,dataType b) {	\
		return a > b ? a : b;			\
	}						\
	__max__;					\
})(da,db)						\

#define dfopen(fp,fileName,mode,errorProcess)	{	\
		if((fp=fopen(fileName,mode))==NULL){	\
			perror(fileName);		\
			errorProcess;			\
		}					\
}

#define dmalloc(p,size)	{				\
		if((p=malloc(size))==NULL){		\
			fatalError("out of memory");	\
		}					\
}
#define drealloc(p,p2,size)	{			\
		if((p=realloc(p2,size))==NULL){		\
			fatalError("out of memory");	\
		}					\
}

#define dcalloc(p,size)	{				\
		if((p=calloc(1,size))==NULL){		\
			fatalError("out of memory");	\
		}					\
}

#define dstrdup(p,str)	{				\
		if((p=strdup(str))==NULL){		\
			fatalError("out of memory");	\
		}					\
}

//交差確認法をする
//#define CLOSS_VALIDATION 1

//トレーニングデータを使ってテストをする
#define TEST 1

//本番予測
//#define MAIN_PREDICTION 1

//量的変数の最適分割ラインを計算する
//#define CALCLATION_SPLITS 1

//ベイジアンネットの構造学習をする
//#define BUILD_BAYESIAN_NET 1