typedef struct {
	void* (*train)(list_t *);
	void* (*predict)(void*);
	int isEqual(void*,void*);
	void (*freeModel)(void*);
}machineLearning;

