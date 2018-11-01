#include "hmm.h"
#include <string.h>
#include <math.h>
#include <assert.h>

void seq_to_state();
void make_Delta();

int main(int argc, char** argv)
{

	HMM hmms[5];
	char* listname;
	char* dataname;
	char* resultname;
	char accu = 0;
	// parse arguments
	if( argc < 4 ){
		fprintf(stderr, "usage: ./test modellist.txt testing_data.txt result.txt\n");
		return 1;
	}
	else{
		listname = argv[1];
		dataname = argv[2];
		resultname = argv[3];
		if (argc>4)
			accu = argv[4][0]=='1';
	}

	int model_count = load_models( listname, hmms, 5);
	int N = 6, T = 50;
	char seq[MAX_SEQ]={};
	FILE* fp = open_or_die( dataname, "r" );
	FILE* fp2 = open_or_die( resultname, "w" );
	while( fscanf(fp, "%s", seq) > 0){
		T = strlen(seq);
		seq_to_state(seq, T);
		int max_index = -1;
		double max_delta_t = -1;
		for(int mod = 0; mod < model_count; mod++){
			double Delta[MAX_SEQ][MAX_STATE];			
			HMM* hmm = &(hmms[mod]);
			make_Delta(Delta, hmm, seq, N, T);
			for(int j = 0; j < N; j++){
				if(Delta[T-1][j]>max_delta_t){
					max_delta_t = Delta[T-1][j];
					max_index = mod;
				}
			}
		}

		if(accu)
			fprintf(fp2, "%s %e\n", hmms[max_index].model_name, max_delta_t);
		else
			fprintf(fp2, "%s\n", hmms[max_index].model_name);
	}
	fclose(fp);
	fclose(fp2);
	return 0;
}








void seq_to_state(char* seq, int len){
	for(int i = 0; i < len; i++){
		seq[i] = seq[i] - 'A';
	}
	return;
}
void make_Delta(double Delta[MAX_SEQ][MAX_STATE], HMM* hmm, char* Seq, int N, int T){
	for(int i = 0; i < N; i++){
		Delta[0][i] = hmm->initial[i]*hmm->observation[Seq[0]][i];
	}
	for(int t = 1; t < T; t++){
		for(int j = 0; j < N; j++){
			double max_del_a = -1;
			for(int i = 0; i < N; i++){
				double cur = Delta[t-1][i]*hmm->transition[i][j];
				if (cur > max_del_a)
					max_del_a = cur;
			}
			Delta[t][j] = hmm->observation[Seq[t]][j]*max_del_a;
		}
	}
}