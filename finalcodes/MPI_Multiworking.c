#include <mpi.h>
#include <iostream>
#include <algorithm>
#include <utility>
#include <tfhe/tfhe.h>
#include <tfhe/tfhe_io.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <vector>
#include <cassert>
#include <string>
#include <sys/time.h>

#define T_FILE "averageMPI.txt"

void add_full(LweSample *sum, const LweSample *x, const LweSample *y, const LweSample *c, const int32_t nb_bits, const TFheGateBootstrappingCloudKeySet *keyset)
{
	const LweParams *in_out_params = keyset->params->in_out_params;
	if (nb_bits == 0)
	{
		return;
	}
	LweSample *carry = new_LweSample_array(1, in_out_params);
	LweSample *axc = new_LweSample_array(1, in_out_params);
	LweSample *bxc = new_LweSample_array(1, in_out_params);
	bootsCOPY(carry, c, keyset);
	for(int32_t i = 0; i < nb_bits; ++i)
	{
		bootsCOPY(axc, x + i, keyset);
		bootsCOPY(bxc, y + i, keyset);
		bootsXOR(axc, axc, carry, keyset);
		bootsXOR(bxc, bxc, carry, keyset);
		bootsXOR(sum + i, x + i, bxc, keyset);
		bootsAND(axc, axc, bxc, keyset);
		bootsXOR(carry, carry, axc, keyset);
	}

	delete_LweSample_array(1, carry);
	delete_LweSample_array(1, axc);
	delete_LweSample_array(1, bxc);
}

void MPImul(LweSample *result, LweSample *x, LweSample *y, const LweSample *c, const int32_t nb_bits, int rank, const TFheGateBootstrappingCloudKeySet *keyset) 
{
	const LweParams *in_out_params = keyset->params->in_out_params;
	LweSample *tmp = new_LweSample_array(32, in_out_params);

	for(int i = 0; i < nb_bits; i++)
	{
		bootsCONSTANT(tmp + i, 0, keyset);
		bootsCONSTANT(result + i, 0, keyset);
	}

	for(int k = 0; k < nb_bits; k++)
	{
		bootsAND(result + k,  x + k, y + rank, keyset);
	}

	delete_LweSample_array(32, tmp);
}

void mul(LweSample *result, LweSample *a, LweSample *b, const LweSample *c, const int32_t nb_bits, TFheGateBootstrappingCloudKeySet *keyset) 
{
	const LweParams *in_out_params = keyset->params->in_out_params;
	LweSample *sum1 = new_LweSample_array(32, in_out_params); // Current Sum
	LweSample *sum2 = new_LweSample_array(32, in_out_params);
	LweSample *tmp = new_LweSample_array(32, in_out_params); // A AND B

	for (int32_t i = 0; i < nb_bits; ++i)
	{
		bootsCONSTANT(sum1 + i, 0, keyset);
		bootsCONSTANT(result + i, 0, keyset);
		bootsCONSTANT(tmp + i, 0, keyset);
	}

	for (int32_t i = 0; i < nb_bits; ++i)
	{
		for (int32_t k = 0; k < nb_bits; ++k)
		{
			bootsAND(tmp + k, b + k, a + i, keyset);
		}

		add_full(sum1 + i, sum1 + i, tmp, c, nb_bits - i, keyset);
	}

	for (int32_t i = 0; i < nb_bits; ++i)
	{
		bootsCOPY(result + i, sum1 + i, keyset);
	}

	delete_LweSample_array(32, sum1);
	delete_LweSample_array(32, tmp);
}

void dieDrama(LweSample* result, const LweSample* a, const TFheGateBootstrappingCloudKeySet* keyset, const size_t size) 
{
		for(int i = 0; i < size; i++)
        {
                bootsCOPY(result + i, a + i, keyset);
        }
}

void leftShift(LweSample* result, const LweSample* a, const TFheGateBootstrappingCloudKeySet* keyset, const size_t size, int amnt) 
{
	const LweParams * in_out_params = keyset->params->in_out_params;
	LweSample *tmp = new_LweSample_array(32, in_out_params);
        int shifted = size - amnt;

	for(int32_t i = 0; i < size; i++)
	{
		bootsCOPY(tmp + i, result + i, keyset);
	}
        for(int32_t i = 0; i < size; i++)
        {
                bootsCONSTANT(result + i, 0, keyset);
        }
        dieDrama(result + amnt, tmp, keyset, shifted);
}

void cleanup(LweSample* result, const int32_t nb_bits, TFheGateBootstrappingCloudKeySet* keyset) 
{
	for(int i = 0; i < nb_bits; i++)
	{
		bootsCONSTANT(result + i, 0, keyset);
	}
}

int main(int argc, char **argv) {

	MPI_Init(&argc,&argv);
	int world_rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
	int world_size;
	MPI_Comm_size(MPI_COMM_WORLD, &world_size);

	//MPI test
	printf("Hello World from process %d of %d\n", world_rank, world_size);
	MPI_Barrier(MPI_COMM_WORLD);

	printf("reading the key...\n");

	//reads the cloud key from file
	FILE* cloud_key = fopen("cloud.key", "rb");
	TFheGateBootstrappingCloudKeySet* bk = new_tfheGateBootstrappingCloudKeySet_fromFile(cloud_key);
        fclose(cloud_key);

	//if necessary, the params are inside the key
	const TFheGateBootstrappingParameterSet* params = bk->params;

	printf("reading the input...\n");

	//Variables for ciphertexts
	LweSample* ciphertext1 = new_gate_bootstrapping_ciphertext_array(32, params);
	LweSample* ciphertext2 = new_gate_bootstrapping_ciphertext_array(32, params);
	LweSample* ciphertext3 = new_gate_bootstrapping_ciphertext_array(32, params);
	LweSample* anstext1 = new_gate_bootstrapping_ciphertext_array(32, params);
    	LweSample* anstext2 = new_gate_bootstrapping_ciphertext_array(32, params);
    	LweSample* result1 = new_gate_bootstrapping_ciphertext_array(32, params);


	//reads the 2x16 ciphertexts from the cloud file
	FILE* cloud_data = fopen("cloud.data", "rb");
	for (int i = 0; i<32; i++)
		import_gate_bootstrapping_ciphertext_fromFile(cloud_data, &ciphertext1[i], params);

	for (int i = 0; i<32; i++)
		import_gate_bootstrapping_ciphertext_fromFile(cloud_data, &ciphertext2[i], params);

	for (int i = 0; i<32; i++)
		import_gate_bootstrapping_ciphertext_fromFile(cloud_data, &ciphertext3[i], params);

	fclose(cloud_data);

	printf("doing the homomorphic computation...on %d\n", world_rank);

	LweSample* result = new_gate_bootstrapping_ciphertext_array(32, params);

	struct timeval start, end;
	double get_time;
	gettimeofday(&start, NULL);

	MPImul(result, ciphertext1, ciphertext2, ciphertext3, 32, world_rank, bk);
	leftShift(result, result, bk, 32, world_rank);

	printf("writing the answer to file...for %d, Calculating\n", world_rank);

	//Dynamically create files to write output from AND
	std::string r = "r";
	std::string s = std::to_string(world_rank);
	s = r + s;
	std::string name = "file";
	char const *fileparam = s.c_str();
	FILE* answer_data = fopen(fileparam, "wb");

	for(int i = 0; i < 32; i++)
		export_gate_bootstrapping_ciphertext_toFile(answer_data, &result[i], params);
	fclose(answer_data);

	cleanup(result, 32, bk);

	MPI_Barrier(MPI_COMM_WORLD);

	if(world_rank % 2 == 0)
	{
		std::string r_source = "r";
		std::string s1 = std::to_string(world_rank);
		std::string s2 = std::to_string(world_rank + 1);
		s1 = r + s1;
		char const *fileparam1 = s1.c_str();

		FILE* ans1 = fopen(fileparam1, "rb");
		for (int i = 0; i < 32; i++)
			import_gate_bootstrapping_ciphertext_fromFile(ans1, &anstext1[i], params);
		fclose(ans1);

		s2 = r + s2;
		char const *fileparam2 = s2.c_str();

		FILE* ans2 = fopen(fileparam2, "rb");
		for (int i = 0; i < 32; i++)
			import_gate_bootstrapping_ciphertext_fromFile(ans2, &anstext2[i], params);
		fclose(ans2);

		add_full(result1, anstext1, anstext2, ciphertext3, 32, bk);

		std::string ans = "loop1answer.data";
		std::string loc = std::to_string(world_rank);
		ans = loc + ans;
		char const *writefileparam = ans.c_str();

		FILE* answer_data = fopen(writefileparam,"wb");
        	for (int i=0; i<32; i++)
            		export_gate_bootstrapping_ciphertext_toFile(answer_data, &result1[i], params);
        	fclose(answer_data);

		cleanup(anstext1, 32, bk);
		cleanup(anstext2, 32, bk);
		cleanup(result1, 32, bk);
	}

	MPI_Barrier(MPI_COMM_WORLD);

	if(world_rank % 4 == 0)
	{
		std::string r_source = "loop1answer.data";
		std::string s1 = std::to_string(world_rank);
		std::string s2 = std::to_string(world_rank + 2);
		s1 = s1 + r_source;
		char const *fileparam1 = s1.c_str();

		FILE* ans1 = fopen(fileparam1, "rb");
		for (int i = 0; i < 32; i++)
		{
			import_gate_bootstrapping_ciphertext_fromFile(ans1, &anstext1[i], params);
		}
		fclose(ans1);

		s2 = s2 + r_source;
		char const *fileparam2 = s2.c_str();

		FILE* ans2 = fopen(fileparam2, "rb");
		for (int i = 0; i < 32; i++)
		{
			import_gate_bootstrapping_ciphertext_fromFile(ans2, &anstext2[i], params);
		}
		fclose(ans2);

		add_full(result1, anstext1, anstext2, ciphertext3, 32, bk);

		std::string ans = "loop2answer.data";
		std::string loc = std::to_string(world_rank);
		ans = loc + ans;
		char const *writefileparam = ans.c_str();

		FILE* answer_data = fopen(writefileparam, "wb");
		for (int i = 0; i < 32; i++)
			export_gate_bootstrapping_ciphertext_toFile(answer_data, &result1[i], params);
		fclose(answer_data);

	        cleanup(anstext1, 32, bk);
	       	cleanup(anstext2, 32, bk);
	        cleanup(result1, 32, bk);
	}

	MPI_Barrier(MPI_COMM_WORLD);

	if(world_rank % 8 == 0)
	{
		std::string r_source = "loop2answer.data";
		std::string s1 = std::to_string(world_rank);
		std::string s2 = std::to_string(world_rank + 4);
		s1 = s1 + r_source;
		char const *fileparam1 = s1.c_str();

		FILE* ans1 = fopen(fileparam1, "rb");
		for (int i = 0; i < 32; i++)
		{
			import_gate_bootstrapping_ciphertext_fromFile(ans1, &anstext1[i], params);
		}
		fclose(ans1);

		s2 = s2 + r_source;
		char const *fileparam2 = s2.c_str();

		FILE* ans2 = fopen(fileparam2, "rb");
		for (int i = 0; i < 32; i++)
		{
			import_gate_bootstrapping_ciphertext_fromFile(ans2, &anstext2[i], params);
		}
		fclose(ans2);

		add_full(result1, anstext1, anstext2, ciphertext3, 32, bk);

		std::string ans = "loop3answer.data";
		std::string loc = std::to_string(world_rank);
		ans = loc + ans;
		char const *writefileparam = ans.c_str();

		FILE* answer_data = fopen(writefileparam, "wb");
		for (int i =0; i < 32; i++)
			export_gate_bootstrapping_ciphertext_toFile(answer_data, &result1[i], params);
		fclose(answer_data);

	        cleanup(anstext1, 32, bk);
	        cleanup(anstext2, 32, bk);
	        cleanup(result1, 32, bk);
	}

	MPI_Barrier(MPI_COMM_WORLD);

	if(world_rank % 16 == 0)
	{
		std::string r_source = "loop3answer.data";
	        std::string s1 = std::to_string(world_rank);
	        std::string s2 = std::to_string(world_rank + 8);
	        s1 = s1 + r_source;
	        char const *fileparam1 = s1.c_str();

		FILE* ans1 = fopen(fileparam1, "rb");
	        for (int i = 0; i < 32; i++)
        		import_gate_bootstrapping_ciphertext_fromFile(ans1, &anstext1[i], params);
		fclose(ans1);

	        s2 = s2 + r_source;
	        char const *fileparam2 = s2.c_str();

		FILE* ans2 = fopen(fileparam2, "rb");
	        for (int i = 0; i < 32; i++)
			import_gate_bootstrapping_ciphertext_fromFile(ans2, &anstext2[i], params);
		fclose(ans2);

        	add_full(result1, anstext1, anstext2, ciphertext3, 32, bk);

		std::string ans = "loop4answer.data";
		std::string loc = std::to_string(world_rank);
		ans = loc + ans;
		char const *writefileparam = ans.c_str();

        	FILE* answer_data = fopen(writefileparam, "wb");
	       	for (int i =0; i < 32; i++)
			export_gate_bootstrapping_ciphertext_toFile(answer_data, &result1[i], params);
        	fclose(answer_data);

        	cleanup(anstext1, 32, bk);
        	cleanup(anstext2, 32, bk);
        	cleanup(result1, 32, bk);
	}

	if(world_rank % 32 == 0)
	{
                std::string r_source = "loop4answer.data";
                std::string s1 = std::to_string(world_rank);
                std::string s2 = std::to_string(world_rank + 16);
                s1 = s1 + r_source;
                char const *fileparam1 = s1.c_str();

                FILE* ans1 = fopen(fileparam1, "rb");
                for (int i = 0; i < 32; i++)
                        import_gate_bootstrapping_ciphertext_fromFile(ans1, &anstext1[i], params);
                fclose(ans1);

                s2 = s2 + r_source;
                char const *fileparam2 = s2.c_str();

                FILE* ans2 = fopen(fileparam2, "rb");
                for (int i = 0; i < 32; i++)
                        import_gate_bootstrapping_ciphertext_fromFile(ans2, &anstext2[i], params);
                fclose(ans2);

                add_full(result1, anstext1, anstext2, ciphertext3, 32, bk);

	        FILE* answer_data = fopen("answer.data", "wb");
       		for (int i =0; i < 32; i++)
            		export_gate_bootstrapping_ciphertext_toFile(answer_data, &result1[i], params);
        	fclose(answer_data);

        	cleanup(anstext1, 32, bk);
        	cleanup(anstext2, 32, bk);
        	cleanup(result1, 32, bk);

	}

	MPI_Barrier(MPI_COMM_WORLD);

	gettimeofday(&end, NULL);

	get_time = (end.tv_sec - start.tv_sec) + (end.tv_usec - start.tv_usec) * 1.0E-6;

	printf("Computation Time: %lf[sec]\n", get_time);
	MPI_Barrier(MPI_COMM_WORLD);

	double time_taken;
	MPI_Reduce(&get_time, &time_taken, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

	FILE *t_file;
	if(world_rank == 0)
	{
		printf("avg = %lf\n", time_taken / 32);
		t_file = fopen(T_FILE, "a");
		fprintf(t_file, "%lf\n", time_taken / 32);
		fclose(t_file);
	}

	//clean up all pointers
	delete_gate_bootstrapping_ciphertext_array(32, result);
	delete_gate_bootstrapping_ciphertext_array(32, ciphertext2);
	delete_gate_bootstrapping_ciphertext_array(32, ciphertext1);
	delete_gate_bootstrapping_ciphertext_array(32, ciphertext3);
	delete_gate_bootstrapping_ciphertext_array(32, anstext1);
	delete_gate_bootstrapping_ciphertext_array(32, anstext2);
	delete_gate_bootstrapping_ciphertext_array(32, result1);
	delete_gate_bootstrapping_cloud_keyset(bk);

	MPI_Finalize();
}
