#include "hmm.h"
#include <string.h>
#include <math.h>
#include <assert.h>

void seq_to_state();
void make_Alpha_Beta();
void make_Gamma();
void make_Epsilon();
void clean_sum();
void update_sum();
void update_model();

int main(int argc, char** argv)
{

	HMM hmm_body;
	char* init_name;
	char* seq_name;
	char* model_name;
	int ITER = 1;
	
	// parse arguments
	if( argc != 5 ){
		fprintf(stderr, "usage: ./train iteration model_init.txt seq_model_01.txt model_01.txt\n");
		return 1;
	}
	else{
		sscanf(argv[1], "%d", &ITER);
		init_name = argv[2];
		seq_name = argv[3];
		model_name = argv[4];
		// printf("%d %s %s %s\n", ITER, init_name, seq_name, model_name);
	}
	
	HMM* hmm = &(hmm_body);
	loadHMM( hmm, init_name );

	while(ITER--){
		char seq[MAX_SEQ]={};
		int samples = 0, N = hmm->state_num, T = 50;

		double new_pi[MAX_STATE];
		double sum_epsilon[MAX_STATE][MAX_STATE];
		double sum_gamma[MAX_STATE];
		double sum_gamma_match[MAX_OBSERV][MAX_STATE];

		clean_sum(new_pi, sum_epsilon, sum_gamma, sum_gamma_match);

		FILE* fp = open_or_die( seq_name, "r" );
		while( fscanf(fp, "%s", seq) > 0){
			samples++;
			double Alpha[MAX_SEQ][MAX_STATE]={};
			double Beta[MAX_SEQ][MAX_STATE]={};
			double Gamma[MAX_SEQ][MAX_STATE]={};
			double Epsilon[MAX_SEQ][MAX_STATE][MAX_STATE]={};

			T = strlen(seq);
			seq_to_state(seq, T);

			make_Alpha_Beta(Alpha, Beta, hmm, seq, N, T);
			// make_Beta(Beta, hmm, seq, N, T);
			make_Gamma(Gamma, Alpha, Beta, N, T);
			make_Epsilon(Epsilon, Alpha, Beta, hmm, seq, N, T);

			update_sum(new_pi, sum_epsilon, sum_gamma, sum_gamma_match, Gamma, Epsilon, seq, N, T);

		}
		fclose(fp);
		
		update_model(hmm, new_pi, sum_epsilon, sum_gamma, sum_gamma_match, N, (double)samples);

		// dumpHMM( stderr, hmm );
	}

	FILE* fp = open_or_die( model_name, "w" );
	dumpHMM( fp, hmm );
	fclose(fp);	


	return 0;
}








void seq_to_state(char* seq, int len){
	for(int i = 0; i < len; i++){
		seq[i] = seq[i] - 'A';
	}
}
void make_Alpha_Beta(double Alpha[MAX_SEQ][MAX_STATE], double Beta[MAX_SEQ][MAX_STATE], HMM* hmm, char* Seq, int N, int T){
	double sum[MAX_SEQ] = {0};
	
	// make Alpha
	for(int i = 0; i < N; i++){
		Alpha[0][i] = hmm->initial[i]*hmm->observation[Seq[0]][i];
		// sum[0] += Alpha[0][i];
	}	
	for(int t = 0; t < T-1; t++){
		for(int j = 0; j < N; j++){
			double sum_alpha_a = 0.0;
			for(int i = 0; i < N; i++){
				sum_alpha_a += Alpha[t][i]*hmm->transition[i][j];
			}
			Alpha[t+1][j] = sum_alpha_a*hmm->observation[Seq[t+1]][j];
			// sum[t+1] += Alpha[t+1][j];
		}
	}


	// make Beta 
	for(int i = 0; i < N; i++){
		Beta[T-1][i] = 1.0;
	}
	for(int t = T-2; t >= 0; t--){
		for(int i = 0; i < N; i++){
			Beta[t][i] = 0;
			for(int j = 0; j < N; j++){
				Beta[t][i] += hmm->transition[i][j]*hmm->observation[Seq[t+1]][j]*Beta[t+1][j];
			}
		}
	}

	// Normalization
	// for(int t = 0; t < T; t++){
	// 	if(sum[t]>0){
	// 		for(int i = 0; i < N; i++){
	// 			Beta[t][i] /= sum[t];
	// 			Alpha[t][i] /= sum[t];
	// 		}
	// 	}
	// }
}
void make_Gamma(double Gamma[MAX_SEQ][MAX_STATE], double Alpha[MAX_SEQ][MAX_STATE], double Beta[MAX_SEQ][MAX_STATE], int N, int T){
	for(int t = 0; t < T; t++){
		double sum = 0.0;
		for(int i = 0; i < N; i++){
			Gamma[t][i] = Alpha[t][i]*Beta[t][i];
			sum += Gamma[t][i];
		}

		if(sum>0){
			for(int i = 0; i < N; i++){
				Gamma[t][i] /= sum;
			}
		}
	}
}
void make_Epsilon(double Epsilon[MAX_SEQ][MAX_STATE][MAX_STATE], double Alpha[MAX_SEQ][MAX_STATE], double Beta[MAX_SEQ][MAX_STATE], HMM* hmm, char* Seq, int N, int T){
	for(int t = 0; t < T-1; t++){
		double sum = 0.0;
		for(int i = 0; i < N; i++){
			for(int j = 0; j < N; j++){
				Epsilon[t][i][j] = Alpha[t][i]*hmm->transition[i][j]*hmm->observation[Seq[t+1]][j]*Beta[t+1][j];
				sum += Epsilon[t][i][j];
			}
		}
				// *sigh*
		if(sum>0){
			for(int i = 0; i < N; i++){
				for(int j = 0; j < N; j++){
					Epsilon[t][i][j] /= sum;
				}
			}
		}
	}
}
void clean_sum(double new_pi[MAX_STATE],
	double sum_epsilon[MAX_STATE][MAX_STATE],
	double sum_gamma[MAX_STATE],
	double sum_gamma_match[MAX_OBSERV][MAX_STATE]){

	for(int i = 0; i < MAX_STATE; i++){
		new_pi[i] = 0;
		sum_gamma[i] = 0;
		for(int j = 0; j < MAX_STATE; j++){
			sum_epsilon[i][j] = 0;
		}
		for(int j = 0; j < MAX_OBSERV; j++){
			sum_gamma_match[j][i] = 0;
		}
	}
}
void update_sum(double new_pi[MAX_STATE],
	double sum_epsilon[MAX_STATE][MAX_STATE],
	double sum_gamma[MAX_STATE],
	double sum_gamma_match[MAX_OBSERV][MAX_STATE],
	double Gamma[MAX_SEQ][MAX_STATE],
	double Epsilon[MAX_SEQ][MAX_STATE][MAX_STATE], char* Seq, int N, int T){
	
	for( int i = 0; i < N; i++){	
		new_pi[i] += Gamma[0][i];
		for(int j = 0; j < N; j++){
			for(int t = 0; t < T-1; t++){
				sum_epsilon[i][j] += Epsilon[t][i][j];
			}
		}

		for(int t = 0; t < T - 1; t++){
			sum_gamma[i] += Gamma[t][i];
					// for(int k = 0; k < N; k++){
					// 	if(k==seq[t])
					// 		+= Gamma[t][i]
					// }

			sum_gamma_match[Seq[t]][i] += Gamma[t][i];
		}
	}
}
void update_model(HMM* hmm, double new_pi[MAX_STATE],
	double sum_epsilon[MAX_STATE][MAX_STATE],
	double sum_gamma[MAX_STATE],
	double sum_gamma_match[MAX_OBSERV][MAX_STATE], int N, double samples){
	// Update HMM
		// initial
	for( int i = 0; i < N; i++){
		hmm->initial[i] = new_pi[i]/samples;
	}
		// transition
	for( int i = 0; i < N; i++){
		for(int j = 0; j < N; j++){
			double upper = sum_epsilon[i][j];
			double lower = sum_gamma[i];
			hmm->transition[i][j] = upper>0? (upper / lower) : 0;
		}
	}
		// emission
	for( int j = 0; j < N; j++){
		for( int k = 0; k < hmm->observ_num; k++){
			double upper = sum_gamma_match[k][j];
			double lower = sum_gamma[j];
			hmm->observation[k][j] = upper>0? (upper / lower) : 0;
		}
	}
}