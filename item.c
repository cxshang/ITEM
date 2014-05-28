
//  gcc -o item item.c ./docitem.so  -I/usr/include/python2.7 -lpython2.7


#include <Python.h>
#include <stdio.h>
#include <string.h>
#include "docitem.h"

int main(int argc, char *argv[]) 
	{
	char tips[]="     ITEM processes graph in undirected, unweighted form. Edgelist must be two values separated with space, and nodes index must start from 0. Indexs of nodes and edges must be sequential without gap. Item can output the resulting communities to a file. In output file, each line is a community which contains some nodes.\nUsages: ./item fn st th\n     In above line, fn is the file name of input network. The st value indicates whether the output indexs starting from 0 or 1. st can only be 0 or 1. The th value is the cut threshold when selecting a seed, if the ratio of overlapped nodes between the seed and an already selected seed is greater than th, the seed will be deleted. The th value is among {0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0}.\nExamples:\nitem filename 0 0.5 > outputfilename\n";
		
	if (argc != 4)
		{fprintf(stderr, tips);  exit(1);}


	if (!(atoi(argv[2]) ==0 || atoi(argv[2]) ==1))
		{fprintf(stderr, tips);  exit(1);}
		

	if (atof(argv[3])<0.1 || atof(argv[3])>1.0) 
		{fprintf(stderr, tips);  exit(1);}

		
	fprintf(stderr, "ITEM starting\n");
	Py_Initialize();
	initdocitem();
	item(argv[1],atoi(argv[2]),atof(argv[3]));
	Py_Finalize();
	return 0;
	}
	
