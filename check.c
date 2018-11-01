#include "hmm.h"
#include <string.h>
#include <math.h>
#include <assert.h>
#define ITER 2

int main(int argc, char** argv)
{
	assert(argc == 4);
	char* file1 = argv[1], *file2 = argv[2];
	char name1[50], name2[50];

	FILE* fp = open_or_die( file1, "r" );
	FILE* fp2 = open_or_die( file2, "r" );
	int same = 0, total = 0;
	int error[5][5]={}, i, j;
	while( fscanf(fp, "%s", name1) > 0){
		fscanf(fp2, "%s", name2);
		if(strcmp(name1, name2)==0)
			same++;
		else{
			// model_01.txt
			i = name1[7]-'1';
			j = name2[7]-'1';
			// if(i==4 || j==4)
			// 	continue;
			error[i][j]++;
		}
		total++;
	}
	fclose(fp);
	fclose(fp2);
	printf("result: %d same out of %d. accuracy: %.6f\n", same, total, (double)same/total);
	for(i = 0; i < 5; i++){
		for(j = 0; j < 5; j++){
			printf("%3d ", error[i][j]);
		}
		puts("");
	}
	return 0;
}