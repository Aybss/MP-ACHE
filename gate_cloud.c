#include <tfhe/tfhe.h>
#include <tfhe/tfhe_io.h>
#include <stdio.h>
#include <time.h>
#include <sys/time.h>

void add(LweSample *sum, const LweSample *x, const LweSample *y, const LweSample *c, const int32_t nb_bits, const TFheGateBootstrappingCloudKeySet *keyset)
{
        const LweParams *in_out_params = keyset->params->in_out_params;

        //LweSample *carry = new_LweSample_array(1, in_out_params);
        //LweSample *a = new_LweSample_array(32, in_out_params);
        //LweSample *b = new_LweSample_array(32, in_out_params);

	//bootsCOPY(sum, x, keyset);

	for(int32_t i = 0; i < 3; i++)
	{
		bootsCOPY(sum + i, x + i, keyset);
	}

        //delete_LweSample_array(1, carry);
        //delete_LweSample_array(32, a);
        //delete_LweSample_array(32, b);
}

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

	delete_LweSample_array(1, carry);
	delete_LweSample_array(1, axc);
	delete_LweSample_array(1, bxc);
}

int main() {

    printf("reading the key...\n");

    // reads the cloud key from file
    FILE* cloud_key = fopen("cloud.key","rb");
    TFheGateBootstrappingCloudKeySet* bk = new_tfheGateBootstrappingCloudKeySet_fromFile(cloud_key);
    fclose(cloud_key);

    // if necessary, the params are inside the key
    const TFheGateBootstrappingParameterSet* params = bk->params;

    printf("reading the input...\n");

    // read the 2x32 ciphertexts
    LweSample* ciphertext1 = new_gate_bootstrapping_ciphertext_array(32, params);
    LweSample* ciphertext2 = new_gate_bootstrapping_ciphertext_array(32, params);
    LweSample* ciphertext3 = new_gate_bootstrapping_ciphertext_array(32, params);

    // reads the 2x32 ciphertexts from the cloud file
    FILE* cloud_data = fopen("cloud.data","rb");
    for (int i=0; i<32; i++)
        import_gate_bootstrapping_ciphertext_fromFile(cloud_data, &ciphertext1[i], params);

    for (int i=0; i<32; i++)
        import_gate_bootstrapping_ciphertext_fromFile(cloud_data, &ciphertext2[i], params);

    for (int i = 0; i<32; i++)
        import_gate_bootstrapping_ciphertext_fromFile(cloud_data, &ciphertext3[i], params);

    fclose(cloud_data);

    printf("doing the homomorphic computation...\n");

    // do some operations on the ciphertexts: here, we will compute the
    // addition of the two
    LweSample* result = new_gate_bootstrapping_ciphertext_array(32, params);

    struct timeval start, end;
    double get_time;
    gettimeofday(&start, NULL);
    add(result, ciphertext1, ciphertext2, ciphertext3, 32, bk);
    gettimeofday(&end, NULL);
    get_time = (end.tv_sec - start.tv_sec) + (end.tv_usec - start.tv_usec) * 1.0E-6;
    printf("Computation Time: %lf[sec]\n", get_time);

    printf("writing the answer to file...\n");

    // export the 32 ciphertexts to a file (for the cloud)
    FILE* answer_data = fopen("answer.data","wb");
    for (int i=0; i<3; i++)
        export_gate_bootstrapping_ciphertext_toFile(answer_data, &result[i], params);
    fclose(answer_data);

    // clean up all pointers
    delete_gate_bootstrapping_ciphertext_array(32, result);
    delete_gate_bootstrapping_ciphertext_array(32, ciphertext2);
    delete_gate_bootstrapping_ciphertext_array(32, ciphertext1);
    delete_gate_bootstrapping_ciphertext_array(32, ciphertext3);
    delete_gate_bootstrapping_cloud_keyset(bk);
}
