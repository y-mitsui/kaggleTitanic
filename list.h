#include <stdio.h>

typedef struct _cell_t{
	struct _cell_t *next;
	struct _cell_t *prev;
	void *data;
}cell_t;

typedef struct {
	cell_t *first;
	cell_t *current;
	size_t size;
}list_t;

void* listSearch(list_t *list,void *arg,int (*func)(void *,void *));
void addList(list_t *list,void *data);
void allDeleteList(cell_t *p);
void removeList(list_t *list,cell_t *cur);
void revivalList(list_t *list,cell_t *cur);
