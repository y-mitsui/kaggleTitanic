#include <stdio.h>
#include <stdlib.h>
#include "list.h"
#include "memory.h"
#include "kaggleTitanic.h"
#include "error.h"

void* listSearch(list_t *list,void *arg,int (*func)(void *,void *)){
	cell_t *cur;
	for(cur=list->first;cur;cur=cur->next){
		if(func(arg,cur->data)) return cur->data;
	}
	return NULL;
}
void addList(list_t *list,void *data){
	cell_t *tmp;

	tmp=(cell_t*)calloc(1,sizeof(cell_t));
	tmp->data=data;
	tmp->prev=list->current;
	if(list->current){
		list->current->next=tmp;
	}else
		list->first=tmp;
	list->current=tmp;
	list->size++;
}
void removeList(list_t *list,cell_t *cur){
	if(cur->prev) cur->prev->next=cur->next;
	else list->first=cur->next;
	if(cur->next) cur->next->prev=cur->prev;
	else list->current=cur->prev;
}
void revivalList(list_t *list,cell_t *cur){
	if(cur->prev) cur->prev->next=cur;
	else list->first=cur;
	if(cur->next) cur->next->prev=cur;
	else list->current=cur;
	
}
/*!
 ******************************************************************************
 * \fn int setList(cell_t* p, void* data)
 * \param *p �}���������Z���̃|�C���^�i���̒���ɃZ����}������j
 * \param data �f�[�^
 * \
 *
 ******************************************************************************
 */
cell_t*
setList(cell_t* p, void* data)
{
  cell_t *tmp;
	tmp=(cell_t*)malloc(sizeof(cell_t));

  tmp->next = p->next;
  tmp->data = data;
  p->next = tmp;
  return tmp;
}

/*!
 ******************************************************************************
 * \fn int deleteList(cell_t* p)
 * \brief �����Ŏ󂯎�����Z���̎��̃|�C���^�̃f�[�^���폜
 * \param *p �폜�������Z�����w���Ă���|�C���^
 * \return -1:�G���[ 0:����I��
 *
 ******************************************************************************
 */
int
deleteList(cell_t* p)
{
  cell_t *tmp;
  if(p->next == NULL)
  {
    return(-1);
  }
  tmp = p->next;
  p->next = tmp->next;
  free(tmp);
  return(0);
}


/*!
 ******************************************************************************
 * \fn void allDeleteList(cell_t *p)
 * \brief �f�[�^�̑S�폜
 * \param *p �擪�̃Z���̃|�C���^
 * \return None
 *
 ******************************************************************************
 */
void
allDeleteList(cell_t *p)
{
  cell_t *tmp;
  while (p->next != NULL)
  {
    tmp = p->next;
    p->next = tmp->next;
    free(tmp);
  }
}
