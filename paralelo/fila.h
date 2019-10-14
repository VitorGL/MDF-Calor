#include <stdio.h>
#include <stdlib.h>

typedef struct list
{
	struct list *prev;
	struct list *next;
	void *chave;
}List;

void alocList(List **fila)
{
	*fila = (List *) malloc(sizeof(List));
}

void exib(List **fila)
{
	if(*fila == NULL)
	{
		printf("Nulo\n");
		return;	
	}
	printf(	"-----\n"
			"| %d |\n"
			"-----\n", (int *)((*fila)->chave));
	if((*fila)->next != NULL)
		exib(&(*fila)->next);

	return;
}

int emptyList(List **no)
{
	if(*no == NULL)
	{
		return 1;
	}
	return 0;
}

void insert(List **no, void *chave)
{
	
	if(emptyList(no))
	{
		alocList(no);

		(*no)->next = NULL;
		(*no)->prev = NULL;
		(*no)->chave = chave;

		return;
	}

	if(emptyList(&(*no)->next))
	{
		alocList(&(*no)->next);

		(*no)->next->next = NULL;
		(*no)->next->prev = *no;
		(*no)->next->chave = chave;

		return;
	}

	insert(&(*no)->next, chave);
	return;
}

void invert(List **entrada, List **saida)
{
	if(*entrada == NULL)
		return;

	if((*entrada)->next != NULL)
		invert(&(*entrada)->next, saida);

	insert(saida, (*entrada)->chave);

	(*entrada)->next = NULL;
	(*entrada)->prev = NULL;
	free(*entrada);
	*entrada = NULL;
}

void *remover(List **item)
{
	if((*item)->next == NULL)
	{
		void *retorno = (void *)((*item)->chave);

		// (*item)->prev->next = NULL;

		*item = NULL;
		return retorno;
	}
	remover(&(*item)->next);
}

void *pop(List **fila)
{
	List *saida = NULL;

	invert(fila, &saida);

	void *retorno = remover(&saida);

	invert(&saida, fila);

	return retorno;
}