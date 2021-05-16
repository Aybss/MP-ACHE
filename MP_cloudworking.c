#include <string>
#include <iostream>
#include <algorithm>
#include <utility>
#include <tfhe/tfhe.h>
#include <tfhe/tfhe_io.h>
#include <stdio.h>
#include <time.h>
#include <vector>
#include <cassert>
#include<sys/time.h>

#define T_FILE "averagestandard.txt"

void addmp(LweSample *sum, const LweSample *x, const LweSample *y, const LweSample *c, const int32_t nb_bits, const TFheGateBootstrappingCloudKeySet *keyset)
{
        const LweParams *in_out_params = keyset->params->in_out_params;

        LweSample *carry = new_LweSample_array(1, in_out_params);
        LweSample *axc = new_LweSample_array(1, in_out_params);
        LweSample *bxc = new_LweSample_array(1, in_out_params);

        bootsCOPY(carry, c, keyset);

        for(int32_t  i = 0; i < nb_bits; i++)
        {
                #pragma omp parallel sections num_threads(2)
                {
                        #pragma omp section
                        bootsXOR(axc, x + i, carry, keyset);
                        #pragma omp section
                        bootsXOR(bxc, y + i, carry, keyset);
                }
                #pragma omp parallel sections num_threads(2)
                {
                        #pragma omp section
                        bootsXOR(sum + i, x + i, bxc, keyset);
                        #pragma omp section
                        bootsAND(axc, axc, bxc, keyset);
                }
                bootsXOR(carry, carry, axc, keyset);
        }
}

void mulmp(LweSample *result, LweSample *a, LweSample *b, const LweSample *c, const int32_t nb_bits, TFheGateBootstrappingCloudKeySet *keyset)
{
	const LweParams *in_out_params = keyset->params->in_out_params;

	LweSample **p = new LweSample*[32];

	#pragma omp parallel for num_threads(32)
	for (int32_t i = 0; i < nb_bits; i++)
	{
		p[i] = new_LweSample_array(32, in_out_params);
		for (int32_t k = 0; k < nb_bits; k++)
		{
			// Initialize the arrays to 0
			bootsCONSTANT(p[i] + k, 0, keyset);
		}
	}

	#pragma omp parallel for num_threads(32)
	for (int32_t i = 0; i < nb_bits; i++)
	{
		for (int32_t k = 0; k < nb_bits - i; k++)
		{
			bootsAND(&p[i][i+k], &a[k], &b[i], keyset);
		}
	}

	#pragma omp parallel sections num_threads(16)
	{
		#pragma omp section
		addmp(p[0], p[0], p[1], c, 32, keyset);
		#pragma omp section
                addmp(p[2], p[2], p[3], c, 32, keyset);
                #pragma omp section
                addmp(p[4], p[4], p[5], c, 32, keyset);
                #pragma omp section
                addmp(p[6], p[6], p[7], c, 32, keyset);

                #pragma omp section
                addmp(p[8], p[8], p[9], c, 32, keyset);
                #pragma omp section
                addmp(p[10], p[10], p[11], c, 32, keyset);
                #pragma omp section
                addmp(p[12], p[12], p[13], c, 32, keyset);
                #pragma omp section
                addmp(p[14], p[14], p[15], c, 32, keyset);

                #pragma omp section
                addmp(p[16], p[16], p[17], c, 32, keyset);
                #pragma omp section
                addmp(p[18], p[18], p[19], c, 32, keyset);
                #pragma omp section
                addmp(p[20], p[20], p[21], c, 32, keyset);
                #pragma omp section
                addmp(p[22], p[22], p[23], c, 32, keyset);

                #pragma omp section
                addmp(p[24], p[24], p[25], c, 32, keyset);
                #pragma omp section
                addmp(p[26], p[26], p[27], c, 32, keyset);
                #pragma omp section
                addmp(p[28], p[28], p[29], c, 32, keyset);
                #pragma omp section
                addmp(p[30], p[30], p[31], c, 32, keyset);
	}

	#pragma omp parallel sections num_threads(8)
	{
		#pragma omp section
                addmp(p[0], p[0], p[2], c, 32, keyset);
                #pragma omp section
                addmp(p[4], p[4], p[6], c, 32, keyset);
                #pragma omp section
                addmp(p[8], p[8], p[10], c, 32, keyset);
                #pragma omp section
                addmp(p[12], p[12], p[14], c, 32, keyset);

                #pragma omp section
                addmp(p[16], p[16], p[18], c, 32, keyset);
                #pragma omp section
                addmp(p[20], p[20], p[22], c, 32, keyset);
                #pragma omp section
                addmp(p[24], p[24], p[26], c, 32, keyset);
                #pragma omp section
                addmp(p[28], p[28], p[30], c, 32, keyset);
	}

	#pragma omp parallel sections num_threads(4)
	{
                #pragma omp section
                addmp(p[0], p[0], p[4], c, 32, keyset);
                #pragma omp section
                addmp(p[8], p[8], p[12], c, 32, keyset);

                #pragma omp section
                addmp(p[16], p[16], p[20], c, 32, keyset);
                #pragma omp section
                addmp(p[24], p[24], p[28], c, 32, keyset);
	}

	#pragma omp parallel sections num_threads(2)
	{
                #pragma omp section
                addmp(p[0], p[0], p[8], c, 32, keyset);
                #pragma omp section
                addmp(p[16], p[16], p[24], c, 32, keyset);
	}
	addmp(result, p[0], p[16], c, 32, keyset);

        for(int i = 0; i < nb_bits; i++)
                delete_LweSample_array(32, p[i]);
        delete[] p;
}

int main() {

	printf("reading the key...\n");

	// reads the cloud key from file
	FILE* cloud_key = fopen("cloud.key", "rb");
	TFheGateBootstrappingCloudKeySet* bk = new_tfheGateBootstrappingCloudKeySet_fromFile(cloud_key);
	fclose(cloud_key);

	// if necessary, the params are inside the key
	const TFheGateBootstrappingParameterSet* params = bk->params;

	printf("reading the input...\n");

	// read the 2x32 ciphertexts
	LweSample* ciphertext1 = new_gate_bootstrapping_ciphertext_array(32, params);
	LweSample* ciphertext2 = new_gate_bootstrapping_ciphertext_array(32, params);
	LweSample* ciphertext3 = new_gate_bootstrapping_ciphertext_array(32, params);

	// reads the 2x16 ciphertexts from the cloud file
	FILE* cloud_data = fopen("cloud.data", "rb");
	for (int i = 0; i<32; i++)
		import_gate_bootstrapping_ciphertext_fromFile(cloud_data, &ciphertext1[i], params);

	for (int i = 0; i<32; i++)
		import_gate_bootstrapping_ciphertext_fromFile(cloud_data, &ciphertext2[i], params);

	for (int i = 0; i<32; i++)
		import_gate_bootstrapping_ciphertext_fromFile(cloud_data, &ciphertext3[i], params);

	fclose(cloud_data);

	printf("doing the homomorphic computation...\n");

	// do some operations on the ciphertexts: here, we will compute the
	// product of the two
	LweSample* result = new_gate_bootstrapping_ciphertext_array(32, params);

	struct timeval start, end;
	double get_time;
	gettimeofday(&start, NULL);
	mulmp(result, ciphertext1, ciphertext2, ciphertext3, 32, bk);
	gettimeofday(&end, NULL);
        get_time = (end.tv_sec - start.tv_sec) + (end.tv_usec - start.tv_usec) * 1.0E-6;
        printf("Computation Time: %lf[sec]\n", get_time);

	// write computation time to file
	FILE *t_file;
	t_file = fopen(T_FILE, "a");
	fprintf(t_file, "%lf\n", get_time);
	fclose(t_file);

	printf("writing the answer to file...\n");

	// export the 32 ciphertexts to a file (for the cloud)
	FILE* answer_data = fopen("answer.data", "wb");
	for (int i = 0; i<32; i++)
		export_gate_bootstrapping_ciphertext_toFile(answer_data, &result[i], params);
	fclose(answer_data);

	// clean up all pointers
	delete_gate_bootstrapping_ciphertext_array(32, result);
	delete_gate_bootstrapping_ciphertext_array(32, ciphertext1);
	delete_gate_bootstrapping_ciphertext_array(32, ciphertext2);
	delete_gate_bootstrapping_ciphertext_array(32, ciphertext3);
	delete_gate_bootstrapping_cloud_keyset(bk);
}
