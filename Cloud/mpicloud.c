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
#include <sys/time.h>
#include <omp.h>
#include <fstream>
#include <mpi.h>
#include <string>
using namespace std;
ifstream read;

#define SIZE 1024
#define T_FILE "averageMPI.txt"

void add_addition(LweSample *sum, LweSample *carryover, const LweSample *x, const LweSample *y, const LweSample *c, const int32_t nb_bits, const TFheGateBootstrappingCloudKeySet *keyset)
{
    const LweParams *in_out_params = keyset->params->in_out_params;
    LweSample *carry = new_LweSample_array(1, in_out_params);
    LweSample *axc = new_LweSample_array(1, in_out_params);
    LweSample *bxc = new_LweSample_array(1, in_out_params);
    bootsCOPY(carry, c, keyset);
    for (int32_t i = 0; i < nb_bits; i++)
    {
        bootsXOR(axc, x + i, carry, keyset);
        bootsXOR(bxc, y + i, carry, keyset);
        bootsXOR(sum + i, x + i, bxc, keyset);
        bootsAND(axc, axc, bxc, keyset);
        bootsXOR(carry, carry, axc, keyset);
    }

    bootsCOPY(carryover, carry, keyset);

    delete_LweSample_array(1, carry);
    delete_LweSample_array(1, axc);
    delete_LweSample_array(1, bxc);
}

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
    for (int32_t i = 0; i < nb_bits; ++i)
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

void cleanup(LweSample *result, const int32_t nb_bits, TFheGateBootstrappingCloudKeySet *keyset)
{
    for (int i = 0; i < nb_bits; i++)
    {
        bootsCONSTANT(result + i, 0, keyset);
    }
}

void MPImul(LweSample *result, LweSample *x, LweSample *y, const LweSample *c, const int32_t nb_bits, int rank, const TFheGateBootstrappingCloudKeySet *keyset)
{
    const LweParams *in_out_params = keyset->params->in_out_params;
    LweSample *tmp = new_LweSample_array(32, in_out_params);

    for (int i = 0; i < nb_bits; i++)
    {
        bootsCONSTANT(tmp + i, 0, keyset);
        bootsCONSTANT(result + i, 0, keyset);
    }

    for (int k = 0; k < nb_bits; k++)
    {
        bootsAND(result + k, x + k, y + rank, keyset);
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

void dieDrama(LweSample *result, const LweSample *a, const TFheGateBootstrappingCloudKeySet *keyset, const size_t size)
{
    for (int i = 0; i < size; i++)
    {
        bootsCOPY(result + i, a + i, keyset);
    }
}

void leftShift(LweSample *result, const LweSample *a, const TFheGateBootstrappingCloudKeySet *keyset, const size_t size, int amnt)
{
    const LweParams *in_out_params = keyset->params->in_out_params;
    LweSample *tmp = new_LweSample_array(32, in_out_params);
    int shifted = size - amnt;

    for (int32_t i = 0; i < size; i++)
    {
        bootsCOPY(tmp + i, result + i, keyset);
    }
    for (int32_t i = 0; i < size; i++)
    {
        bootsCONSTANT(result + i, 0, keyset);
    }
    dieDrama(result + amnt, tmp, keyset, shifted);
}

void zero(LweSample *result, const TFheGateBootstrappingCloudKeySet *keyset, const size_t size)
{
    for (int i = 0; i < size; i++)
    {
        bootsCONSTANT(result + i, 0, keyset);
    }
}

void NOT(LweSample *result, const LweSample *x, const TFheGateBootstrappingCloudKeySet *keyset, const size_t size)
{
    for (int i = 0; i < size; i++)
    {
        bootsNOT(result + i, x + i, keyset);
    }
}

int main(int argc, char **argv)
{
    char test[SIZE];
    FILE *f = fopen("operator.txt", "r");
    fgets(test, SIZE, f);
    int x = atoi(test);
    printf("%d", x);

    FILE *cloud_key = fopen("cloud.key", "rb");
    TFheGateBootstrappingCloudKeySet *bk = new_tfheGateBootstrappingCloudKeySet_fromFile(cloud_key);
    fclose(cloud_key);

    FILE *nbit_key = fopen("nbit.key", "rb");
    TFheGateBootstrappingSecretKeySet *nbitkey = new_tfheGateBootstrappingSecretKeySet_fromFile(nbit_key);
    fclose(nbit_key);

    const TFheGateBootstrappingParameterSet *params = bk->params;

    const TFheGateBootstrappingParameterSet *nbitparams = nbitkey->params;

    LweSample *ciphertextbit = new_gate_bootstrapping_ciphertext_array(32, nbitparams);

    LweSample *ciphertextnegative1 = new_gate_bootstrapping_ciphertext_array(32, nbitparams);
    LweSample *ciphertextbit1 = new_gate_bootstrapping_ciphertext_array(32, nbitparams);

    LweSample *ciphertextnegative2 = new_gate_bootstrapping_ciphertext_array(32, nbitparams);
    LweSample *ciphertextbit2 = new_gate_bootstrapping_ciphertext_array(32, nbitparams);

    LweSample *ciphertext1 = new_gate_bootstrapping_ciphertext_array(32, params);
    LweSample *ciphertext2 = new_gate_bootstrapping_ciphertext_array(32, params);
    LweSample *ciphertext3 = new_gate_bootstrapping_ciphertext_array(32, params);
    LweSample *ciphertext4 = new_gate_bootstrapping_ciphertext_array(32, params);
    LweSample *ciphertext5 = new_gate_bootstrapping_ciphertext_array(32, params);
    LweSample *ciphertext6 = new_gate_bootstrapping_ciphertext_array(32, params);
    LweSample *ciphertext7 = new_gate_bootstrapping_ciphertext_array(32, params);
    LweSample *ciphertext8 = new_gate_bootstrapping_ciphertext_array(32, params);
    LweSample *ciphertext9 = new_gate_bootstrapping_ciphertext_array(32, params);
    LweSample *ciphertext10 = new_gate_bootstrapping_ciphertext_array(32, params);
    LweSample *ciphertext11 = new_gate_bootstrapping_ciphertext_array(32, params);
    LweSample *ciphertext12 = new_gate_bootstrapping_ciphertext_array(32, params);
    LweSample *ciphertext13 = new_gate_bootstrapping_ciphertext_array(32, params);
    LweSample *ciphertext14 = new_gate_bootstrapping_ciphertext_array(32, params);
    LweSample *ciphertext15 = new_gate_bootstrapping_ciphertext_array(32, params);
    LweSample *ciphertext16 = new_gate_bootstrapping_ciphertext_array(32, params);

    LweSample *ciphertextcarry1 = new_gate_bootstrapping_ciphertext_array(32, params);
    LweSample *ciphertextcarry2 = new_gate_bootstrapping_ciphertext_array(32, params);

    printf("Reading input 1...\n");

    FILE *cloud_data = fopen("cloud.data", "rb");
    for (int i = 0; i < 32; i++) // 
        import_gate_bootstrapping_ciphertext_fromFile(cloud_data, &ciphertextnegative1[i], nbitparams);
    for (int i = 0; i < 32; i++) // 
        import_gate_bootstrapping_ciphertext_fromFile(cloud_data, &ciphertextbit1[i], nbitparams);

    int32_t int_bit1 = 0;
    for (int i = 0; i < 32; i++)
    {
        int ai = bootsSymDecrypt(&ciphertextbit1[i], nbitkey) > 0;
        int_bit1 |= (ai << i);
    }

    for (int i = 0; i < 32; i++)
        import_gate_bootstrapping_ciphertext_fromFile(cloud_data, &ciphertext1[i], params);
    for (int i = 0; i < 32; i++)
        import_gate_bootstrapping_ciphertext_fromFile(cloud_data, &ciphertext2[i], params);
    for (int i = 0; i < 32; i++)
        import_gate_bootstrapping_ciphertext_fromFile(cloud_data, &ciphertext3[i], params);
    for (int i = 0; i < 32; i++)
        import_gate_bootstrapping_ciphertext_fromFile(cloud_data, &ciphertext4[i], params);
    for (int i = 0; i < 32; i++)
        import_gate_bootstrapping_ciphertext_fromFile(cloud_data, &ciphertext5[i], params);
    for (int i = 0; i < 32; i++)
        import_gate_bootstrapping_ciphertext_fromFile(cloud_data, &ciphertext6[i], params);
    for (int i = 0; i < 32; i++)
        import_gate_bootstrapping_ciphertext_fromFile(cloud_data, &ciphertext7[i], params);
    for (int i = 0; i < 32; i++)
        import_gate_bootstrapping_ciphertext_fromFile(cloud_data, &ciphertext8[i], params);

    for (int i = 0; i < 32; i++) // 
        import_gate_bootstrapping_ciphertext_fromFile(cloud_data, &ciphertextcarry1[i], params);

    printf("Reading input 2...\n");

    for (int i = 0; i < 32; i++) // 
        import_gate_bootstrapping_ciphertext_fromFile(cloud_data, &ciphertextnegative2[i], nbitparams);
    for (int i = 0; i < 32; i++) // 
        import_gate_bootstrapping_ciphertext_fromFile(cloud_data, &ciphertextbit2[i], nbitparams);

    int32_t int_bit2 = 0;
    for (int i = 0; i < 32; i++)
    {
        int ai = bootsSymDecrypt(&ciphertextbit2[i], nbitkey) > 0;
        int_bit2 |= (ai << i);
    }

    for (int i = 0; i < 32; i++)
        import_gate_bootstrapping_ciphertext_fromFile(cloud_data, &ciphertext9[i], params);
    for (int i = 0; i < 32; i++)
        import_gate_bootstrapping_ciphertext_fromFile(cloud_data, &ciphertext10[i], params);
    for (int i = 0; i < 32; i++)
        import_gate_bootstrapping_ciphertext_fromFile(cloud_data, &ciphertext11[i], params);
    for (int i = 0; i < 32; i++)
        import_gate_bootstrapping_ciphertext_fromFile(cloud_data, &ciphertext12[i], params);
    for (int i = 0; i < 32; i++)
        import_gate_bootstrapping_ciphertext_fromFile(cloud_data, &ciphertext13[i], params);
    for (int i = 0; i < 32; i++)
        import_gate_bootstrapping_ciphertext_fromFile(cloud_data, &ciphertext14[i], params);
    for (int i = 0; i < 32; i++)
        import_gate_bootstrapping_ciphertext_fromFile(cloud_data, &ciphertext15[i], params);
    for (int i = 0; i < 32; i++)
        import_gate_bootstrapping_ciphertext_fromFile(cloud_data, &ciphertext16[i], params);

    for (int i = 0; i < 32; i++) // line21
        import_gate_bootstrapping_ciphertext_fromFile(cloud_data, &ciphertextcarry2[i], params);

    printf("Reading operation code...\n");

    LweSample *ciphertextnegative = new_gate_bootstrapping_ciphertext_array(32, nbitparams);
    LweSample *carry1 = new_gate_bootstrapping_ciphertext_array(32, params);

    int32_t int_negative1 = 0;
    for (int i = 0; i < 32; i++)
    {
        int ai = bootsSymDecrypt(&ciphertextnegative1[i], nbitkey) > 0;
        int_negative1 |= (ai << i);
    }
    std::cout << int_negative1 << " => negative1"
              << "\n";

    if (int_negative1 == 2)
    {
        int_negative1 = 1;
    }

    // Decrypts Negative2
    int32_t int_negative2 = 0;
    for (int i = 0; i < 32; i++)
    {
        int ai = bootsSymDecrypt(&ciphertextnegative2[i], nbitkey) > 0;
        int_negative2 |= (ai << i);
    }
    std::cout << int_negative2 << " => negative2"
              << "\n";
    int32_t int_negative;
    int_negative = (int_negative1 + int_negative2);

    FILE *answer_data = fopen("answer.data", "wb");

    int32_t ciphernegative = 0;
    if (int_negative == 1)
    {
        ciphernegative = 1;
    }
    if (int_negative == 2)
    {
        ciphernegative = 2;
    }
    if (int_negative == 3)
    {
        ciphernegative = 4;
    }
    for (int i = 0; i < 32; i++)
    {
        bootsSymEncrypt(&ciphertextnegative[i], (ciphernegative >> i) & 1, nbitkey);
    }
    for (int i = 0; i < 32; i++)
        export_gate_bootstrapping_ciphertext_toFile(answer_data, &ciphertextnegative[i], nbitparams);
    std::cout << ciphernegative << " => total negatives"
              << "\n";

    delete_gate_bootstrapping_ciphertext_array(32, ciphertextnegative);

    int32_t int_op;
    read.open("operator.txt");
    read >> int_op;

    int32_t int_bit = 0;
    if (int_op == 4)
    {
        if (int_bit1 >= int_bit2)
        {
            int_bit = (int_bit1 * 2);
        }
        else
        {
            int_bit = (int_bit2 * 2);
        }
        for (int i = 0; i < 32; i++)
        {
            bootsSymEncrypt(&ciphertextbit[i], (int_bit >> i) & 1, nbitkey);
        }
        for (int i = 0; i < 32; i++)
            export_gate_bootstrapping_ciphertext_toFile(answer_data, &ciphertextbit[i], nbitparams);
        std::cout << int_bit << " written to answer.data"
                  << "\n";
        if (int_bit1 >= int_bit2)
        {
            int_bit = int_bit1;
        }
        else
        {
            int_bit = int_bit2;
        }
    }
    else if (int_bit1 >= int_bit2)
    {
        int_bit = int_bit1;
        for (int i = 0; i < 32; i++)
            export_gate_bootstrapping_ciphertext_toFile(answer_data, &ciphertextbit1[i], nbitparams);
        std::cout << int_bit << " written to answer.data"
                  << "\n";
    }
    else
    {
        int_bit = int_bit2;
        for (int i = 0; i < 32; i++)
            export_gate_bootstrapping_ciphertext_toFile(answer_data, &ciphertextbit2[i], nbitparams);
        std::cout << int_bit << " written to answer.data"
                  << "\n";
    }

    fclose(cloud_data);

    printf("MPICloud.c is being computed [mpicloud.c]\n");
    if ((int_op == 1 && (int_negative != 1 && int_negative != 2)) || (int_op == 2 && (int_negative == 1 || int_negative == 2)))
    {
        if (int_op == 1)
        {
            std::cout << int_bit << " bit Addition computation"
                      << "\n";
        }
        else
        {
            std::cout << int_bit << " bit Subtraction computation"
                      << "\n";
        }

        std::cout << int_bit << " is the int_bit";
        printf("\n");
        if (int_bit == 32)
        {
            printf("MPI Addition is being computed 32 bit[mpicloud.c]");
            MPI_Init(&argc, &argv);
            int world_rank;
            MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
            int world_size;
            MPI_Comm_size(MPI_COMM_WORLD, &world_size);


            printf("Hello World from process %d of %d\n", world_rank, world_size);

            LweSample *result = new_gate_bootstrapping_ciphertext_array(32, params);
            LweSample *carry1 = new_gate_bootstrapping_ciphertext_array(32, params);

            struct timeval start, end;
            double get_time;
            gettimeofday(&start, NULL);

            printf("Doing the homomorphic computation...\n");

            add_addition(result, carry1, ciphertext1, ciphertext9, ciphertextcarry1, 32, bk);

            gettimeofday(&end, NULL);
            get_time = (end.tv_sec - start.tv_sec) + (end.tv_usec - start.tv_usec) * 1.0E-6;
            printf("Computation Time: %lf[sec]\n", get_time);

            for (int i = 0; i < 32; i++) // 
                export_gate_bootstrapping_ciphertext_toFile(answer_data, &result[i], params);
            for (int i = 0; i < 32; i++) // 2
                export_gate_bootstrapping_ciphertext_toFile(answer_data, &ciphertextcarry1[i], params);
            for (int i = 0; i < 32; i++) // 3
                export_gate_bootstrapping_ciphertext_toFile(answer_data, &ciphertextcarry1[i], params);
            for (int i = 0; i < 32; i++) // 4
                export_gate_bootstrapping_ciphertext_toFile(answer_data, &ciphertextcarry1[i], params);
            for (int i = 0; i < 32; i++) // 5
                export_gate_bootstrapping_ciphertext_toFile(answer_data, &ciphertextcarry1[i], params);
            for (int i = 0; i < 32; i++) // 6
                export_gate_bootstrapping_ciphertext_toFile(answer_data, &ciphertextcarry1[i], params);
            for (int i = 0; i < 32; i++) // 7
                export_gate_bootstrapping_ciphertext_toFile(answer_data, &ciphertextcarry1[i], params);
            for (int i = 0; i < 32; i++) // 8
                export_gate_bootstrapping_ciphertext_toFile(answer_data, &ciphertextcarry1[i], params);
            for (int i = 0; i < 32; i++) // carry
                export_gate_bootstrapping_ciphertext_toFile(answer_data, &ciphertextcarry1[i], params);
            fclose(answer_data);

            printf("writing the answer to file...\n");


            printf("writing the answer to file...for %d, Calculating\n", world_rank);

            get_time = (end.tv_sec - start.tv_sec) + (end.tv_usec - start.tv_usec) * 1.0E-6;

            printf("Computation Time: %lf[sec]\n", get_time);
            MPI_Barrier(MPI_COMM_WORLD);

            double time_taken;
            MPI_Reduce(&get_time, &time_taken, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

            FILE *t_file;
            if (world_rank == 0)
            {
                printf("avg = %lf\n", time_taken / 32);
                t_file = fopen(T_FILE, "a");
                fprintf(t_file, "%lf\n", time_taken / 32);
                fclose(t_file);
            }

            delete_gate_bootstrapping_ciphertext_array(32, result);
            delete_gate_bootstrapping_ciphertext_array(32, carry1);
            delete_gate_bootstrapping_ciphertext_array(32, ciphertextbit1);
            delete_gate_bootstrapping_ciphertext_array(32, ciphertextnegative1);
            delete_gate_bootstrapping_ciphertext_array(32, ciphertextbit2);
            delete_gate_bootstrapping_ciphertext_array(32, ciphertextnegative2);
            delete_gate_bootstrapping_ciphertext_array(32, ciphertext1);
            delete_gate_bootstrapping_ciphertext_array(32, ciphertext9);
            delete_gate_bootstrapping_ciphertext_array(32, ciphertextcarry1);

            delete_gate_bootstrapping_cloud_keyset(bk);
            delete_gate_bootstrapping_secret_keyset(nbitkey);
            MPI_Finalize();
        }
        else if (int_bit == 64)
        {
            printf("MPI Addition is being computed 64 bit [mpicloud.c]\n");
            MPI_Init(&argc, &argv);
            int world_rank;
            MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
            int world_size;
            MPI_Comm_size(MPI_COMM_WORLD, &world_size);


            printf("Hello World from process %d of %d\n", world_rank, world_size);
            MPI_Barrier(MPI_COMM_WORLD);

            LweSample *result = new_gate_bootstrapping_ciphertext_array(32, params);
            LweSample *result2 = new_gate_bootstrapping_ciphertext_array(32, params);
            LweSample *result3 = new_gate_bootstrapping_ciphertext_array(32, params);
            LweSample *carry1 = new_gate_bootstrapping_ciphertext_array(32, params);
            LweSample *carry2 = new_gate_bootstrapping_ciphertext_array(32, params);
            LweSample *carry3 = new_gate_bootstrapping_ciphertext_array(32, params);
            LweSample *temp1 = new_gate_bootstrapping_ciphertext_array(32, params);
            LweSample *temp2 = new_gate_bootstrapping_ciphertext_array(32, params);
            LweSample *temp3 = new_gate_bootstrapping_ciphertext_array(32, params);
            LweSample *temp4 = new_gate_bootstrapping_ciphertext_array(32, params);

            LweSample *anstext1 = new_gate_bootstrapping_ciphertext_array(32, params);
            LweSample *anstext2 = new_gate_bootstrapping_ciphertext_array(32, params);

            struct timeval start, end;
            double get_time;
            gettimeofday(&start, NULL);

            printf("Doing the homomorphic computation (64bit Addition)...\n");

            zero(temp1, bk, 32);

            bootsCONSTANT(temp2, 1, bk);

            add_addition(result, carry1, ciphertext1, ciphertext9, ciphertextcarry1, 32, bk);
            MPI_Barrier(MPI_COMM_WORLD);

            //carry = 0
            if (world_rank == 0)
            {
                printf("Rank (i=0) %d\n", world_rank);
                time_t t;
                time(&t);
                printf("current time is : %s", ctime(&t));

                add_addition(result2, carry2, ciphertext2, ciphertext10, temp1, 32, bk);

                std::string name = "result2.data";
                char const *fileparam = name.c_str();
                FILE *answer1 = fopen(fileparam, "wb");
                for (int i = 0; i < 32; i++)
                    export_gate_bootstrapping_ciphertext_toFile(answer1, &result2[i], params);
                fclose(answer1);

                cleanup(result2, 32, bk);
            }

            //carry == 1
            if (world_rank == 1)
            {
                printf("Rank (i=1) %d\n", world_rank);
                time_t t;
                time(&t);
                printf("current time is : %s", ctime(&t));

                add_addition(result3, carry3, ciphertext2, ciphertext10, temp2, 32, bk);

                std::string name1 = "result3.data";
                char const *fileparam1 = name1.c_str();
                FILE *answer2 = fopen(fileparam1, "wb");
                for (int i = 0; i < 32; i++)
                    export_gate_bootstrapping_ciphertext_toFile(answer2, &result3[i], params);
                fclose(answer2);

                cleanup(result3, 32, bk);
            }

            MPI_Barrier(MPI_COMM_WORLD);

            std::string name2 = "result2.data";
            char const *fileparam3 = name2.c_str();
            FILE *ans1 = fopen(fileparam3, "rb");
            for (int i = 0; i < 32; i++)
                import_gate_bootstrapping_ciphertext_fromFile(ans1, &anstext1[i], params);
            fclose(ans1);

            std::string name3 = "result3.data";
            char const *fileparam4 = name3.c_str();
            FILE *ans2 = fopen(fileparam4, "rb");
            for (int i = 0; i < 32; i++)
                import_gate_bootstrapping_ciphertext_fromFile(ans2, &anstext2[i], params);
            fclose(ans2);

            // Timings
            gettimeofday(&end, NULL);
            get_time = (end.tv_sec - start.tv_sec) + (end.tv_usec - start.tv_usec) * 1.0E-6;
            printf("Computation Time: %lf[sec]\n", get_time);

            // export the result ciphertexts to a file (for the cloud)
            for (int i = 0; i < 32; i++) // result1 behind
                export_gate_bootstrapping_ciphertext_toFile(answer_data, &result[i], params);
            for (int i = 0; i < 32; i++) // result1 infront
                export_gate_bootstrapping_ciphertext_toFile(answer_data, &carry1[i], params);
            for (int i = 0; i < 32; i++) // result2 binary2
                export_gate_bootstrapping_ciphertext_toFile(answer_data, &anstext1[i], params);
            for (int i = 0; i < 32; i++) // result3 binary3
                export_gate_bootstrapping_ciphertext_toFile(answer_data, &anstext2[i], params);
            /* for (int i = 0; i < 32; i++) // 5
                export_gate_bootstrapping_ciphertext_toFile(answer_data, &ciphertextcarry1[i], params);
            for (int i = 0; i < 32; i++) // 6
                export_gate_bootstrapping_ciphertext_toFile(answer_data, &ciphertextcarry1[i], params);*/
            for (int i = 0; i < 32; i++) // 5
                export_gate_bootstrapping_ciphertext_toFile(answer_data, &ciphertextcarry1[i], params);
            for (int i = 0; i < 32; i++) // 6
                export_gate_bootstrapping_ciphertext_toFile(answer_data, &ciphertextcarry1[i], params);
            for (int i = 0; i < 32; i++) // 7
                export_gate_bootstrapping_ciphertext_toFile(answer_data, &ciphertextcarry1[i], params);
            for (int i = 0; i < 32; i++) // 8
                export_gate_bootstrapping_ciphertext_toFile(answer_data, &ciphertextcarry1[i], params);
            for (int i = 0; i < 32; i++) // carry
                export_gate_bootstrapping_ciphertext_toFile(answer_data, &ciphertextcarry1[i], params);
            fclose(answer_data);

            MPI_Barrier(MPI_COMM_WORLD);

            printf("writing the answer to file...\n");

            //Clean up

            printf("writing the answer to file...for %d, Calculating\n", world_rank);

            get_time = (end.tv_sec - start.tv_sec) + (end.tv_usec - start.tv_usec) * 1.0E-6;

            printf("Computation Time: %lf[sec]\n", get_time);
            MPI_Barrier(MPI_COMM_WORLD);

            double time_taken;
            MPI_Reduce(&get_time, &time_taken, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

            FILE *t_file;
            if (world_rank == 0)
            {
                printf("avg = %lf\n", time_taken / 32);
                t_file = fopen(T_FILE, "a");
                fprintf(t_file, "%lf\n", time_taken / 32);
                fclose(t_file);
            }

            delete_gate_bootstrapping_ciphertext_array(32, result);
            delete_gate_bootstrapping_ciphertext_array(32, result2);
            delete_gate_bootstrapping_ciphertext_array(32, result3);
            delete_gate_bootstrapping_ciphertext_array(32, carry1);
            delete_gate_bootstrapping_ciphertext_array(32, carry2);
            delete_gate_bootstrapping_ciphertext_array(32, carry3);
            delete_gate_bootstrapping_ciphertext_array(32, ciphertextbit1);
            delete_gate_bootstrapping_ciphertext_array(32, ciphertextnegative1);
            delete_gate_bootstrapping_ciphertext_array(32, ciphertextbit2);
            delete_gate_bootstrapping_ciphertext_array(32, ciphertextnegative2);
            delete_gate_bootstrapping_ciphertext_array(32, ciphertext1);
            delete_gate_bootstrapping_ciphertext_array(32, ciphertext2);
            delete_gate_bootstrapping_ciphertext_array(32, ciphertext9);
            delete_gate_bootstrapping_ciphertext_array(32, ciphertext10);
            delete_gate_bootstrapping_ciphertext_array(32, temp1);
            delete_gate_bootstrapping_ciphertext_array(32, temp2);
            delete_gate_bootstrapping_ciphertext_array(32, temp3);
            delete_gate_bootstrapping_ciphertext_array(32, temp4);
            delete_gate_bootstrapping_ciphertext_array(32, anstext1);
            delete_gate_bootstrapping_ciphertext_array(32, anstext2);
            delete_gate_bootstrapping_cloud_keyset(bk);
            delete_gate_bootstrapping_secret_keyset(nbitkey);

            MPI_Finalize();
        }

        else if (int_bit == 128)
        {
            printf("MPI Addition is being computed 128 bit [mpicloud.c]\n");
            MPI_Init(&argc, &argv);
            int world_rank;
            MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
            int world_size;
            MPI_Comm_size(MPI_COMM_WORLD, &world_size);

            //MPI test

            printf("Hello World from process %d of %d\n", world_rank, world_size);
            MPI_Barrier(MPI_COMM_WORLD);

            // Ciphertext to hold the result and carry
            LweSample *result = new_gate_bootstrapping_ciphertext_array(32, params);
            LweSample *result2 = new_gate_bootstrapping_ciphertext_array(32, params);
            LweSample *result3 = new_gate_bootstrapping_ciphertext_array(32, params);
            LweSample *result4 = new_gate_bootstrapping_ciphertext_array(32, params);
            LweSample *result5 = new_gate_bootstrapping_ciphertext_array(32, params);
            LweSample *result6 = new_gate_bootstrapping_ciphertext_array(32, params);
            LweSample *result7 = new_gate_bootstrapping_ciphertext_array(32, params);

            LweSample *carry1 = new_gate_bootstrapping_ciphertext_array(32, params);
            LweSample *carry2 = new_gate_bootstrapping_ciphertext_array(32, params);
            LweSample *carry3 = new_gate_bootstrapping_ciphertext_array(32, params);
            LweSample *carry4 = new_gate_bootstrapping_ciphertext_array(32, params);
            LweSample *carry5 = new_gate_bootstrapping_ciphertext_array(32, params);
            LweSample *carry6 = new_gate_bootstrapping_ciphertext_array(32, params);
            LweSample *carry7 = new_gate_bootstrapping_ciphertext_array(32, params);

            LweSample *temp1 = new_gate_bootstrapping_ciphertext_array(32, params);
            LweSample *temp2 = new_gate_bootstrapping_ciphertext_array(32, params);
            LweSample *temp3 = new_gate_bootstrapping_ciphertext_array(32, params);
            LweSample *temp4 = new_gate_bootstrapping_ciphertext_array(32, params);
            LweSample *temp5 = new_gate_bootstrapping_ciphertext_array(32, params);
            LweSample *temp6 = new_gate_bootstrapping_ciphertext_array(32, params);

            LweSample *anstext1 = new_gate_bootstrapping_ciphertext_array(32, params);
            LweSample *anstext2 = new_gate_bootstrapping_ciphertext_array(32, params);
            LweSample *anstext3 = new_gate_bootstrapping_ciphertext_array(32, params);
            LweSample *anstext4 = new_gate_bootstrapping_ciphertext_array(32, params);
            LweSample *anstext5 = new_gate_bootstrapping_ciphertext_array(32, params);
            LweSample *anstext6 = new_gate_bootstrapping_ciphertext_array(32, params);

            LweSample *anstextImportCarry2 = new_gate_bootstrapping_ciphertext_array(32, params);
            LweSample *anstextImportCarry3 = new_gate_bootstrapping_ciphertext_array(32, params);
            LweSample *anstextImportCarry4 = new_gate_bootstrapping_ciphertext_array(32, params);
            LweSample *anstextImportCarry5 = new_gate_bootstrapping_ciphertext_array(32, params);
            LweSample *anstextImportCarry6 = new_gate_bootstrapping_ciphertext_array(32, params);
            LweSample *anstextImportCarry7 = new_gate_bootstrapping_ciphertext_array(32, params);

            printf("Doing the homomorphic computation (64bit Addition)...\n");

            //making carry = 0
            zero(temp1, bk, 32);
            zero(temp3, bk, 32);
            zero(temp5, bk, 32);

            //making carry = 1
            bootsCONSTANT(temp2, 1, bk);
            bootsCONSTANT(temp4, 1, bk);
            bootsCONSTANT(temp6, 1, bk);

            struct timeval start, end;
            double get_time;
            gettimeofday(&start, NULL);

            //first 32 bit addition
            add_addition(result, carry1, ciphertext1, ciphertext9, ciphertextcarry1, 32, bk);

            //add_addition(result2, carry2, ciphertext2, ciphertext10, carry1, 32, bk);
            /* add_addition(result2, carry2, ciphertext2, ciphertext10, carry1, 32, bk);
            add_addition(result3, carry3, ciphertext3, ciphertext11, carry2, 32, bk);
            add_addition(result4, carry4, ciphertext4, ciphertext12, carry3, 32, bk);*/

            MPI_Barrier(MPI_COMM_WORLD);

            //carry = 0 result2
            if (world_rank == 0)
            {
                printf("Rank (i=0) %d\n", world_rank);
                time_t t;
                time(&t);
                printf("current time is : %s", ctime(&t));

                add_addition(result2, carry2, ciphertext2, ciphertext10, temp1, 32, bk);

                std::string name = "result2.data";
                char const *fileparam = name.c_str();
                FILE *answer1 = fopen(fileparam, "wb");
                for (int i = 0; i < 32; i++)
                    export_gate_bootstrapping_ciphertext_toFile(answer1, &result2[i], params);
                fclose(answer1);

                std::string nameCarry1 = "carry2.data";
                char const *fileparamCarry1 = nameCarry1.c_str();
                FILE *answerCarry1 = fopen(fileparamCarry1, "wb");
                for (int i = 0; i < 32; i++)
                    export_gate_bootstrapping_ciphertext_toFile(answerCarry1, &carry2[i], params);
                fclose(answerCarry1);

                cleanup(result2, 32, bk);
                cleanup(carry2, 32, bk);
            }
            //carry == 1 result3
            if (world_rank == 1)
            {
                printf("Rank (i=1) %d\n", world_rank);
                time_t t;
                time(&t);
                printf("current time is : %s", ctime(&t));

                add_addition(result3, carry3, ciphertext2, ciphertext10, temp2, 32, bk);

                std::string name1 = "result3.data";
                char const *fileparam1 = name1.c_str();
                FILE *answer2 = fopen(fileparam1, "wb");
                for (int i = 0; i < 32; i++)
                    export_gate_bootstrapping_ciphertext_toFile(answer2, &result3[i], params);
                fclose(answer2);

                std::string nameCarry2 = "carry3.data";
                char const *fileparamCarry2 = nameCarry2.c_str();
                FILE *answerCarry2 = fopen(fileparamCarry2, "wb");
                for (int i = 0; i < 32; i++)
                    export_gate_bootstrapping_ciphertext_toFile(answerCarry2, &carry3[i], params);
                fclose(answerCarry2);

                cleanup(result3, 32, bk);
                cleanup(carry3, 32, bk);
            }

            //carry = 0 result4
            if (world_rank == 2)
            {
                printf("Rank (i=2) %d\n", world_rank);
                time_t t;
                time(&t);
                printf("current time is : %s", ctime(&t));

                add_addition(result4, carry4, ciphertext3, ciphertext11, temp3, 32, bk);

                std::string name2 = "result4.data";
                char const *fileparam2 = name2.c_str();
                FILE *answer3 = fopen(fileparam2, "wb");
                for (int i = 0; i < 32; i++)
                    export_gate_bootstrapping_ciphertext_toFile(answer3, &result4[i], params);
                fclose(answer3);

                std::string nameCarry3 = "carry4.data";
                char const *fileparamCarry3 = nameCarry3.c_str();
                FILE *answerCarry3 = fopen(fileparamCarry3, "wb");
                for (int i = 0; i < 32; i++)
                    export_gate_bootstrapping_ciphertext_toFile(answerCarry3, &carry4[i], params);
                fclose(answerCarry3);

                cleanup(result4, 32, bk);
                cleanup(carry4, 32, bk);
            }
            if (world_rank == 3)
            {
                printf("Rank (i=3) %d\n", world_rank);
                time_t t;
                time(&t);
                printf("current time is : %s", ctime(&t));

                add_addition(result5, carry5, ciphertext3, ciphertext11, temp4, 32, bk);

                std::string name3 = "result5.data";
                char const *fileparam3 = name3.c_str();
                FILE *answer4 = fopen(fileparam3, "wb");
                for (int i = 0; i < 32; i++)
                    export_gate_bootstrapping_ciphertext_toFile(answer4, &result5[i], params);
                fclose(answer4);

                std::string nameCarry4 = "carry5.data";
                char const *fileparamCarry4 = nameCarry4.c_str();
                FILE *answerCarry4 = fopen(fileparamCarry4, "wb");
                for (int i = 0; i < 32; i++)
                    export_gate_bootstrapping_ciphertext_toFile(answerCarry4, &carry5[i], params);
                fclose(answerCarry4);

                cleanup(result5, 32, bk);
                cleanup(carry5, 32, bk);
            }

            if (world_rank == 4)
            {
                printf("Rank (i=4) %d\n", world_rank);
                time_t t;
                time(&t);
                printf("current time is : %s", ctime(&t));

                add_addition(result6, carry6, ciphertext4, ciphertext12, temp5, 32, bk);

                std::string name4 = "result6.data";
                char const *fileparam4 = name4.c_str();
                FILE *answer5 = fopen(fileparam4, "wb");
                for (int i = 0; i < 32; i++)
                    export_gate_bootstrapping_ciphertext_toFile(answer5, &result6[i], params);
                fclose(answer5);

                std::string nameCarry5 = "carry6.data";
                char const *fileparamCarry5 = nameCarry5.c_str();
                FILE *answerCarry5 = fopen(fileparamCarry5, "wb");
                for (int i = 0; i < 32; i++)
                    export_gate_bootstrapping_ciphertext_toFile(answerCarry5, &carry6[i], params);
                fclose(answerCarry5);

                cleanup(result5, 32, bk);
                cleanup(carry6, 32, bk);
            }
            //carry = 1 result7
            if (world_rank == 5)
            {
                printf("Rank (i=5) %d\n", world_rank);
                time_t t;
                time(&t);
                printf("current time is : %s", ctime(&t));

                add_addition(result7, carry7, ciphertext4, ciphertext12, temp6, 32, bk);

                std::string name5 = "result7.data";
                char const *fileparam5 = name5.c_str();
                FILE *answer6 = fopen(fileparam5, "wb");
                for (int i = 0; i < 32; i++)
                    export_gate_bootstrapping_ciphertext_toFile(answer6, &result7[i], params);
                fclose(answer6);

                std::string nameCarry6 = "carry7.data";
                char const *fileparamCarry6 = nameCarry6.c_str();
                FILE *answerCarry6 = fopen(fileparamCarry6, "wb");
                for (int i = 0; i < 32; i++)
                    export_gate_bootstrapping_ciphertext_toFile(answerCarry6, &carry7[i], params);
                fclose(answerCarry6);

                cleanup(result5, 32, bk);
                cleanup(carry7, 32, bk);
            }

            MPI_Barrier(MPI_COMM_WORLD);

            // Timings
            gettimeofday(&end, NULL);
            get_time = (end.tv_sec - start.tv_sec) + (end.tv_usec - start.tv_usec) * 1.0E-6;
            printf("Computation Time [MPI]: %lf[sec]\n", get_time);

            std::string name2 = "result2.data";
            char const *fileparam3 = name2.c_str();
            FILE *ans1 = fopen(fileparam3, "rb");
            for (int i = 0; i < 32; i++)
                import_gate_bootstrapping_ciphertext_fromFile(ans1, &anstext1[i], params);
            fclose(ans1);

            std::string name3 = "result3.data";
            char const *fileparam4 = name3.c_str();
            FILE *ans2 = fopen(fileparam4, "rb");
            for (int i = 0; i < 32; i++)
                import_gate_bootstrapping_ciphertext_fromFile(ans2, &anstext2[i], params);
            fclose(ans2);

            std::string name4 = "result4.data";
            char const *fileparam5 = name4.c_str();
            FILE *ans3 = fopen(fileparam5, "rb");
            for (int i = 0; i < 32; i++)
                import_gate_bootstrapping_ciphertext_fromFile(ans3, &anstext3[i], params);
            fclose(ans3);

            std::string name5 = "result5.data";
            char const *fileparam6 = name5.c_str();
            FILE *ans4 = fopen(fileparam6, "rb");
            for (int i = 0; i < 32; i++)
                import_gate_bootstrapping_ciphertext_fromFile(ans4, &anstext4[i], params);
            fclose(ans4);

            std::string name6 = "result6.data";
            char const *fileparam7 = name6.c_str();
            FILE *ans5 = fopen(fileparam7, "rb");
            for (int i = 0; i < 32; i++)
                import_gate_bootstrapping_ciphertext_fromFile(ans5, &anstext5[i], params);
            fclose(ans5);

            std::string name7 = "result7.data";
            char const *fileparam8 = name7.c_str();
            FILE *ans6 = fopen(fileparam8, "rb");
            for (int i = 0; i < 32; i++)
                import_gate_bootstrapping_ciphertext_fromFile(ans6, &anstext6[i], params);
            fclose(ans6);

            ////////////////////////////////////////////////
            //Carry////////////////////////////////////////
            //////////////////////////////////////////////

            std::string nameImportCarry2 = "carry2.data";
            char const *fileparamImportCarry2 = nameImportCarry2.c_str();
            FILE *ansImportCarry1 = fopen(fileparamImportCarry2, "rb");
            for (int i = 0; i < 32; i++)
                import_gate_bootstrapping_ciphertext_fromFile(ansImportCarry1, &anstextImportCarry2[i], params);
            fclose(ansImportCarry1);

            std::string nameImportCarry3 = "carry3.data";
            char const *fileparamImportCarry3 = nameImportCarry3.c_str();
            FILE *ansImportCarry3 = fopen(fileparamImportCarry3, "rb");
            for (int i = 0; i < 32; i++)
                import_gate_bootstrapping_ciphertext_fromFile(ansImportCarry3, &anstextImportCarry3[i], params);
            fclose(ansImportCarry3);

            std::string nameImportCarry4 = "carry4.data";
            char const *fileparamImportCarry4 = nameImportCarry4.c_str();
            FILE *ansImportCarry4 = fopen(fileparamImportCarry4, "rb");
            for (int i = 0; i < 32; i++)
                import_gate_bootstrapping_ciphertext_fromFile(ansImportCarry4, &anstextImportCarry4[i], params);
            fclose(ansImportCarry4);

            std::string nameImportCarry5 = "carry5.data";
            char const *fileparamImportCarry5 = nameImportCarry5.c_str();
            FILE *ansImportCarry5 = fopen(fileparamImportCarry5, "rb");
            for (int i = 0; i < 32; i++)
                import_gate_bootstrapping_ciphertext_fromFile(ansImportCarry5, &anstextImportCarry5[i], params);
            fclose(ansImportCarry5);

            std::string nameImportCarry6 = "carry6.data";
            char const *fileparamImportCarry6 = nameImportCarry6.c_str();
            FILE *ansImportCarry6 = fopen(fileparamImportCarry6, "rb");
            for (int i = 0; i < 32; i++)
                import_gate_bootstrapping_ciphertext_fromFile(ansImportCarry6, &anstextImportCarry6[i], params);
            fclose(ansImportCarry6);

            std::string nameImportCarry7 = "carry7.data";
            char const *fileparamImportCarry7 = nameImportCarry7.c_str();
            FILE *ansImportCarry7 = fopen(fileparamImportCarry7, "rb");
            for (int i = 0; i < 32; i++)
                import_gate_bootstrapping_ciphertext_fromFile(ansImportCarry7, &anstextImportCarry7[i], params);
            fclose(ansImportCarry7);

            for (int i = 0; i < 32; i++) //  behind
                export_gate_bootstrapping_ciphertext_toFile(answer_data, &result[i], params);
            for (int i = 0; i < 32; i++) //  behind
                export_gate_bootstrapping_ciphertext_toFile(answer_data, &carry1[i], params);
            for (int i = 0; i < 32; i++) //  infront
                export_gate_bootstrapping_ciphertext_toFile(answer_data, &anstextImportCarry2[i], params);
            for (int i = 0; i < 32; i++) // result1 infront
                export_gate_bootstrapping_ciphertext_toFile(answer_data, &anstextImportCarry3[i], params);
            for (int i = 0; i < 32; i++) //  infront
                export_gate_bootstrapping_ciphertext_toFile(answer_data, &anstextImportCarry4[i], params);
            for (int i = 0; i < 32; i++) //  infron
                export_gate_bootstrapping_ciphertext_toFile(answer_data, &anstextImportCarry5[i], params);
            for (int i = 0; i < 32; i++) // result1 infront
                export_gate_bootstrapping_ciphertext_toFile(answer_data, &anstextImportCarry6[i], params);
            for (int i = 0; i < 32; i++) //  infront
                export_gate_bootstrapping_ciphertext_toFile(answer_data, &anstextImportCarry7[i], params);
            for (int i = 0; i < 32; i++) //  binary2
                export_gate_bootstrapping_ciphertext_toFile(answer_data, &anstext1[i], params);
            for (int i = 0; i < 32; i++) //  binary3
                export_gate_bootstrapping_ciphertext_toFile(answer_data, &anstext2[i], params);
            for (int i = 0; i < 32; i++) // 5
                export_gate_bootstrapping_ciphertext_toFile(answer_data, &anstext3[i], params);
            for (int i = 0; i < 32; i++) // 6
                export_gate_bootstrapping_ciphertext_toFile(answer_data, &anstext4[i], params);
            for (int i = 0; i < 32; i++) // 7
                export_gate_bootstrapping_ciphertext_toFile(answer_data, &anstext5[i], params);
            for (int i = 0; i < 32; i++) // 8
                export_gate_bootstrapping_ciphertext_toFile(answer_data, &anstext6[i], params);
            for (int i = 0; i < 32; i++) // carry
                export_gate_bootstrapping_ciphertext_toFile(answer_data, &ciphertextcarry1[i], params);
            fclose(answer_data);

            MPI_Barrier(MPI_COMM_WORLD);

            printf("writing the answer to file...\n");

            //Clean up
            printf("writing the answer to file...for %d, Calculating\n", world_rank);
            MPI_Barrier(MPI_COMM_WORLD);

            double time_taken;
            MPI_Reduce(&get_time, &time_taken, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

            FILE *t_file;
            if (world_rank == 0)
            {
                printf("avg = %lf\n", time_taken / 32);
                t_file = fopen(T_FILE, "a");
                fprintf(t_file, "%lf\n", time_taken / 32);
                fclose(t_file);
            }

            delete_gate_bootstrapping_ciphertext_array(32, result);
            delete_gate_bootstrapping_ciphertext_array(32, result2);
            delete_gate_bootstrapping_ciphertext_array(32, result3);
            delete_gate_bootstrapping_ciphertext_array(32, result4);
            delete_gate_bootstrapping_ciphertext_array(32, result5);
            delete_gate_bootstrapping_ciphertext_array(32, result6);
            delete_gate_bootstrapping_ciphertext_array(32, result7);
            delete_gate_bootstrapping_ciphertext_array(32, carry1);
            delete_gate_bootstrapping_ciphertext_array(32, carry2);
            delete_gate_bootstrapping_ciphertext_array(32, carry3);
            delete_gate_bootstrapping_ciphertext_array(32, carry4);
            delete_gate_bootstrapping_ciphertext_array(32, carry5);
            delete_gate_bootstrapping_ciphertext_array(32, carry6);
            delete_gate_bootstrapping_ciphertext_array(32, carry7);
            delete_gate_bootstrapping_ciphertext_array(32, ciphertextbit1);
            delete_gate_bootstrapping_ciphertext_array(32, ciphertextnegative1);
            delete_gate_bootstrapping_ciphertext_array(32, ciphertextbit2);
            delete_gate_bootstrapping_ciphertext_array(32, ciphertextnegative2);
            delete_gate_bootstrapping_ciphertext_array(32, ciphertext1);
            delete_gate_bootstrapping_ciphertext_array(32, ciphertext2);
            delete_gate_bootstrapping_ciphertext_array(32, ciphertext3);
            delete_gate_bootstrapping_ciphertext_array(32, ciphertext4);
            delete_gate_bootstrapping_ciphertext_array(32, ciphertext9);
            delete_gate_bootstrapping_ciphertext_array(32, ciphertext10);
            delete_gate_bootstrapping_ciphertext_array(32, ciphertext11);
            delete_gate_bootstrapping_ciphertext_array(32, ciphertext12);
            delete_gate_bootstrapping_ciphertext_array(32, temp1);
            delete_gate_bootstrapping_ciphertext_array(32, temp2);
            delete_gate_bootstrapping_ciphertext_array(32, temp3);
            delete_gate_bootstrapping_ciphertext_array(32, temp4);
            delete_gate_bootstrapping_ciphertext_array(32, anstext1);
            delete_gate_bootstrapping_ciphertext_array(32, anstext2);
            delete_gate_bootstrapping_ciphertext_array(32, anstext3);
            delete_gate_bootstrapping_ciphertext_array(32, anstext4);
            delete_gate_bootstrapping_ciphertext_array(32, anstext5);
            delete_gate_bootstrapping_ciphertext_array(32, anstext6);
            delete_gate_bootstrapping_ciphertext_array(32, anstextImportCarry2);
            delete_gate_bootstrapping_ciphertext_array(32, anstextImportCarry3);
            delete_gate_bootstrapping_ciphertext_array(32, anstextImportCarry4);
            delete_gate_bootstrapping_ciphertext_array(32, anstextImportCarry5);
            delete_gate_bootstrapping_ciphertext_array(32, anstextImportCarry6);
            delete_gate_bootstrapping_ciphertext_array(32, anstextImportCarry7);

            delete_gate_bootstrapping_cloud_keyset(bk);
            delete_gate_bootstrapping_secret_keyset(nbitkey);

            MPI_Finalize();
        }
    }

    //MPI Sub
    else if (int_op == 2 || (int_op == 1 && (int_negative == 1 || int_negative == 2)))
    {
        if ((int_op == 2 && int_negative == 0) || (int_op == 1 && int_negative == 2))
        {
            if (int_op == 2)
            {
                std::cout << int_bit << " bit Subtraction computation"
                          << "\n";
            }
            else
            {
                std::cout << int_bit << " bit Addition computation with 2nd value negative"
                          << "\n";
            }

            if (int_bit == 32)
            {
                MPI_Init(&argc, &argv);
                int world_rank;
                MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
                int world_size;
                MPI_Comm_size(MPI_COMM_WORLD, &world_size);

                printf("Doing the homomorphic computation...\n");

                LweSample *temp = new_gate_bootstrapping_ciphertext_array(32, params);
                LweSample *inverse1 = new_gate_bootstrapping_ciphertext_array(32, params);

                LweSample *tempcarry1 = new_gate_bootstrapping_ciphertext_array(32, params);

                LweSample *twosresult1 = new_gate_bootstrapping_ciphertext_array(32, params);

                LweSample *twoscarry1 = new_gate_bootstrapping_ciphertext_array(32, params);

                struct timeval start, end;
                double get_time;

                gettimeofday(&start, NULL);


                NOT(inverse1, ciphertext9, bk, 32);

                zero(temp, bk, 32);
                zero(tempcarry1, bk, 32);

                bootsCONSTANT(temp, 1, bk);

                add_addition(twosresult1, twoscarry1, inverse1, temp, tempcarry1, 32, bk);

                LweSample *result1 = new_gate_bootstrapping_ciphertext_array(32, params);

                LweSample *carry1 = new_gate_bootstrapping_ciphertext_array(32, params);

                add_addition(result1, carry1, ciphertext1, twosresult1, ciphertextcarry1, 32, bk);

                gettimeofday(&end, NULL);

                get_time = (end.tv_sec - start.tv_sec) + (end.tv_usec - start.tv_usec) * 1.0E-6;
                printf("Computation Time: %lf[sec]\n", get_time);

                printf("writing the answer to file...\n");

                for (int i = 0; i < 32; i++) //result1
                    export_gate_bootstrapping_ciphertext_toFile(answer_data, &result1[i], params);
                for (int i = 0; i < 32; i++) // 2
                    export_gate_bootstrapping_ciphertext_toFile(answer_data, &ciphertextcarry1[i], params);
                for (int i = 0; i < 32; i++) // 3
                    export_gate_bootstrapping_ciphertext_toFile(answer_data, &ciphertextcarry1[i], params);
                for (int i = 0; i < 32; i++) // 4
                    export_gate_bootstrapping_ciphertext_toFile(answer_data, &ciphertextcarry1[i], params);
                for (int i = 0; i < 32; i++) // 5
                    export_gate_bootstrapping_ciphertext_toFile(answer_data, &ciphertextcarry1[i], params);
                for (int i = 0; i < 32; i++) // 6
                    export_gate_bootstrapping_ciphertext_toFile(answer_data, &ciphertextcarry1[i], params);
                for (int i = 0; i < 32; i++) // 7
                    export_gate_bootstrapping_ciphertext_toFile(answer_data, &ciphertextcarry1[i], params);
                for (int i = 0; i < 32; i++) // 8
                    export_gate_bootstrapping_ciphertext_toFile(answer_data, &ciphertextcarry1[i], params);
                for (int i = 0; i < 32; i++) // carry
                    export_gate_bootstrapping_ciphertext_toFile(answer_data, &ciphertextcarry1[i], params);
                fclose(answer_data);

                delete_gate_bootstrapping_ciphertext_array(32, temp);
                delete_gate_bootstrapping_ciphertext_array(32, inverse1);

                delete_gate_bootstrapping_ciphertext_array(32, tempcarry1);

                delete_gate_bootstrapping_ciphertext_array(32, twosresult1);

                delete_gate_bootstrapping_ciphertext_array(32, twoscarry1);

                delete_gate_bootstrapping_ciphertext_array(32, carry1);

                delete_gate_bootstrapping_ciphertext_array(32, result1);

                delete_gate_bootstrapping_ciphertext_array(32, ciphertextbit1);
                delete_gate_bootstrapping_ciphertext_array(32, ciphertextnegative1);
                delete_gate_bootstrapping_ciphertext_array(32, ciphertextbit2);
                delete_gate_bootstrapping_ciphertext_array(32, ciphertextnegative2);
                delete_gate_bootstrapping_ciphertext_array(32, ciphertext1);
                delete_gate_bootstrapping_ciphertext_array(32, ciphertext9);
                delete_gate_bootstrapping_ciphertext_array(32, ciphertextcarry1);

                delete_gate_bootstrapping_cloud_keyset(bk);
                delete_gate_bootstrapping_secret_keyset(nbitkey);

                MPI_Finalize();
            }
            else if (int_bit == 64)
            {
                MPI_Init(&argc, &argv);
                int world_rank;
                MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
                int world_size;
                MPI_Comm_size(MPI_COMM_WORLD, &world_size);

                printf("Doing the homomorphic computation (64bit Subtraction)...\n");

                LweSample *temp = new_gate_bootstrapping_ciphertext_array(32, params);

                LweSample *inverse1 = new_gate_bootstrapping_ciphertext_array(32, params);
                LweSample *inverse2 = new_gate_bootstrapping_ciphertext_array(32, params);

                LweSample *tempcarry1 = new_gate_bootstrapping_ciphertext_array(32, params);
                LweSample *tempcarry2 = new_gate_bootstrapping_ciphertext_array(32, params);

                LweSample *twosresult1 = new_gate_bootstrapping_ciphertext_array(32, params);
                LweSample *twosresult2 = new_gate_bootstrapping_ciphertext_array(32, params);

                LweSample *twoscarry1 = new_gate_bootstrapping_ciphertext_array(32, params);
                LweSample *twoscarry2 = new_gate_bootstrapping_ciphertext_array(32, params);

                struct timeval start, end;
                double get_time;
                gettimeofday(&start, NULL);

   
                NOT(inverse1, ciphertext9, bk, 32);
                NOT(inverse2, ciphertext10, bk, 32);

                zero(temp, bk, 32);

                zero(tempcarry1, bk, 32);
                zero(tempcarry2, bk, 32);

                bootsCONSTANT(temp, 1, bk);

                add_addition(twosresult1, twoscarry1, inverse1, temp, tempcarry1, 32, bk);

                //10
                add_addition(twosresult2, twoscarry2, inverse2, tempcarry2, twoscarry1, 32, bk);

                ///
                LweSample *result1 = new_gate_bootstrapping_ciphertext_array(32, params);
                LweSample *result2 = new_gate_bootstrapping_ciphertext_array(32, params);
                LweSample *result3 = new_gate_bootstrapping_ciphertext_array(32, params);

                LweSample *carry1 = new_gate_bootstrapping_ciphertext_array(32, params);
                LweSample *carry2 = new_gate_bootstrapping_ciphertext_array(32, params);

  
                add_addition(result1, carry1, ciphertext1, twosresult1, ciphertextcarry1, 32, bk);

                LweSample *makeitzero = new_gate_bootstrapping_ciphertext_array(32, params);
                LweSample *makeitone = new_gate_bootstrapping_ciphertext_array(32, params);

                LweSample *anstext1 = new_gate_bootstrapping_ciphertext_array(32, params);
                LweSample *anstext2 = new_gate_bootstrapping_ciphertext_array(32, params);

                zero(makeitzero, bk, 32);
                bootsCONSTANT(makeitone, 1, bk);

                MPI_Barrier(MPI_COMM_WORLD);

                if (world_rank == 0)
                {
                    printf("Rank (i=0) %d\n", world_rank);
                    time_t t;
                    time(&t);
                    printf("current time is : %s", ctime(&t));

                    add_addition(result2, carry2, ciphertext2, twosresult2, makeitzero, 32, bk);

                    std::string name = "result2sub.data";
                    char const *fileparam = name.c_str();
                    FILE *answer1 = fopen(fileparam, "wb");
                    for (int i = 0; i < 32; i++)
                        export_gate_bootstrapping_ciphertext_toFile(answer1, &result2[i], params);
                    fclose(answer1);

                    cleanup(result2, 32, bk);
                }

                //carry = 1
                if (world_rank == 1)
                {
                    printf("Rank (i=0) %d\n", world_rank);
                    time_t t;
                    time(&t);
                    printf("current time is : %s", ctime(&t));

                    add_addition(result3, carry2, ciphertext2, twosresult2, makeitone, 32, bk);

                    std::string name1 = "result3sub.data";
                    char const *fileparam1 = name1.c_str();
                    FILE *answer2 = fopen(fileparam1, "wb");
                    for (int i = 0; i < 32; i++)
                        export_gate_bootstrapping_ciphertext_toFile(answer2, &result3[i], params);
                    fclose(answer2);

                    cleanup(result3, 32, bk);
                }

                MPI_Barrier(MPI_COMM_WORLD);
                gettimeofday(&end, NULL);
                get_time = (end.tv_sec - start.tv_sec) + (end.tv_usec - start.tv_usec) * 1.0E-6;
                printf("Computation Time [MPI]: %lf[sec]\n", get_time);

                std::string name2 = "result2sub.data";
                char const *fileparam3 = name2.c_str();
                FILE *ans1 = fopen(fileparam3, "rb");
                for (int i = 0; i < 32; i++)
                    import_gate_bootstrapping_ciphertext_fromFile(ans1, &anstext1[i], params);
                fclose(ans1);

                std::string name3 = "result3sub.data";
                char const *fileparam4 = name3.c_str();
                FILE *ans2 = fopen(fileparam4, "rb");
                for (int i = 0; i < 32; i++)
                    import_gate_bootstrapping_ciphertext_fromFile(ans2, &anstext2[i], params);
                fclose(ans2);

                printf("Exporting 64bit Sub... \n");

                printf("writing the answer to file...\n");

                //export the 32 ciphertexts to a file (for the cloud)
                for (int i = 0; i < 32; i++) //result1
                    export_gate_bootstrapping_ciphertext_toFile(answer_data, &result1[i], params);
                for (int i = 0; i < 32; i++) // 2
                    export_gate_bootstrapping_ciphertext_toFile(answer_data, &carry1[i], params);
                for (int i = 0; i < 32; i++) // 3
                    export_gate_bootstrapping_ciphertext_toFile(answer_data, &anstext1[i], params);
                for (int i = 0; i < 32; i++) // 4
                    export_gate_bootstrapping_ciphertext_toFile(answer_data, &anstext2[i], params);
                for (int i = 0; i < 32; i++) // 5
                    export_gate_bootstrapping_ciphertext_toFile(answer_data, &ciphertextcarry1[i], params);
                for (int i = 0; i < 32; i++) // 6
                    export_gate_bootstrapping_ciphertext_toFile(answer_data, &ciphertextcarry1[i], params);
                for (int i = 0; i < 32; i++) // 7
                    export_gate_bootstrapping_ciphertext_toFile(answer_data, &ciphertextcarry1[i], params);
                for (int i = 0; i < 32; i++) // 8
                    export_gate_bootstrapping_ciphertext_toFile(answer_data, &ciphertextcarry1[i], params);
                for (int i = 0; i < 32; i++) // carry
                    export_gate_bootstrapping_ciphertext_toFile(answer_data, &ciphertextcarry1[i], params);
                fclose(answer_data);

                delete_gate_bootstrapping_ciphertext_array(32, temp);
                delete_gate_bootstrapping_ciphertext_array(32, inverse1);
                delete_gate_bootstrapping_ciphertext_array(32, inverse2);
                delete_gate_bootstrapping_ciphertext_array(32, tempcarry1);
                delete_gate_bootstrapping_ciphertext_array(32, tempcarry2);
                delete_gate_bootstrapping_ciphertext_array(32, twosresult1);
                delete_gate_bootstrapping_ciphertext_array(32, twosresult2);
                delete_gate_bootstrapping_ciphertext_array(32, twoscarry1);
                delete_gate_bootstrapping_ciphertext_array(32, twoscarry2);
                delete_gate_bootstrapping_ciphertext_array(32, carry1);
                delete_gate_bootstrapping_ciphertext_array(32, carry2);
                delete_gate_bootstrapping_ciphertext_array(32, result1);
                delete_gate_bootstrapping_ciphertext_array(32, result2);
                delete_gate_bootstrapping_ciphertext_array(32, result3);
                delete_gate_bootstrapping_ciphertext_array(32, anstext1);
                delete_gate_bootstrapping_ciphertext_array(32, anstext2);

                delete_gate_bootstrapping_ciphertext_array(32, ciphertextbit1);
                delete_gate_bootstrapping_ciphertext_array(32, ciphertextnegative1);
                delete_gate_bootstrapping_ciphertext_array(32, ciphertextbit2);
                delete_gate_bootstrapping_ciphertext_array(32, ciphertextnegative2);
                delete_gate_bootstrapping_ciphertext_array(32, ciphertext1);
                delete_gate_bootstrapping_ciphertext_array(32, ciphertext2);
                delete_gate_bootstrapping_ciphertext_array(32, ciphertext9);
                delete_gate_bootstrapping_ciphertext_array(32, ciphertext10);
                delete_gate_bootstrapping_ciphertext_array(32, ciphertextcarry1);

                delete_gate_bootstrapping_cloud_keyset(bk);
                delete_gate_bootstrapping_secret_keyset(nbitkey);

                MPI_Finalize();
            }
            //sub
            else if (int_bit == 128)
            {
                printf("MPI Subtraction is being computed 129 bit [mpicloud.c]\n");
                MPI_Init(&argc, &argv);
                int world_rank;
                MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
                int world_size;
                MPI_Comm_size(MPI_COMM_WORLD, &world_size);

                //MPI test

                printf("Hello World from process %d of %d\n", world_rank, world_size);
                MPI_Barrier(MPI_COMM_WORLD);

                // Ciphertext to hold the result and carry
                LweSample *result = new_gate_bootstrapping_ciphertext_array(32, params);
                LweSample *result2 = new_gate_bootstrapping_ciphertext_array(32, params);
                LweSample *result3 = new_gate_bootstrapping_ciphertext_array(32, params);
                LweSample *result4 = new_gate_bootstrapping_ciphertext_array(32, params);
                LweSample *result5 = new_gate_bootstrapping_ciphertext_array(32, params);
                LweSample *result6 = new_gate_bootstrapping_ciphertext_array(32, params);
                LweSample *result7 = new_gate_bootstrapping_ciphertext_array(32, params);

                LweSample *carry1 = new_gate_bootstrapping_ciphertext_array(32, params);
                LweSample *carry2 = new_gate_bootstrapping_ciphertext_array(32, params);
                LweSample *carry3 = new_gate_bootstrapping_ciphertext_array(32, params);
                LweSample *carry4 = new_gate_bootstrapping_ciphertext_array(32, params);
                LweSample *carry5 = new_gate_bootstrapping_ciphertext_array(32, params);
                LweSample *carry6 = new_gate_bootstrapping_ciphertext_array(32, params);
                LweSample *carry7 = new_gate_bootstrapping_ciphertext_array(32, params);

                LweSample *temp1 = new_gate_bootstrapping_ciphertext_array(32, params);
                LweSample *temp2 = new_gate_bootstrapping_ciphertext_array(32, params);
                LweSample *temp3 = new_gate_bootstrapping_ciphertext_array(32, params);
                LweSample *temp4 = new_gate_bootstrapping_ciphertext_array(32, params);
                LweSample *temp5 = new_gate_bootstrapping_ciphertext_array(32, params);
                LweSample *temp6 = new_gate_bootstrapping_ciphertext_array(32, params);

                LweSample *anstext1 = new_gate_bootstrapping_ciphertext_array(32, params);
                LweSample *anstext2 = new_gate_bootstrapping_ciphertext_array(32, params);
                LweSample *anstext3 = new_gate_bootstrapping_ciphertext_array(32, params);
                LweSample *anstext4 = new_gate_bootstrapping_ciphertext_array(32, params);
                LweSample *anstext5 = new_gate_bootstrapping_ciphertext_array(32, params);
                LweSample *anstext6 = new_gate_bootstrapping_ciphertext_array(32, params);

                LweSample *anstextImportCarry2 = new_gate_bootstrapping_ciphertext_array(32, params);
                LweSample *anstextImportCarry3 = new_gate_bootstrapping_ciphertext_array(32, params);
                LweSample *anstextImportCarry4 = new_gate_bootstrapping_ciphertext_array(32, params);
                LweSample *anstextImportCarry5 = new_gate_bootstrapping_ciphertext_array(32, params);
                LweSample *anstextImportCarry6 = new_gate_bootstrapping_ciphertext_array(32, params);
                LweSample *anstextImportCarry7 = new_gate_bootstrapping_ciphertext_array(32, params);

                LweSample *tempStart = new_gate_bootstrapping_ciphertext_array(32, params);

                LweSample *inverse1 = new_gate_bootstrapping_ciphertext_array(32, params);
                LweSample *inverse2 = new_gate_bootstrapping_ciphertext_array(32, params);
                LweSample *inverse3 = new_gate_bootstrapping_ciphertext_array(32, params);
                LweSample *inverse4 = new_gate_bootstrapping_ciphertext_array(32, params);

                LweSample *twosresult1 = new_gate_bootstrapping_ciphertext_array(32, params);
                LweSample *twosresult2 = new_gate_bootstrapping_ciphertext_array(32, params);
                LweSample *twosresult3 = new_gate_bootstrapping_ciphertext_array(32, params);
                LweSample *twosresult4 = new_gate_bootstrapping_ciphertext_array(32, params);

                LweSample *tempcarry1 = new_gate_bootstrapping_ciphertext_array(32, params);
                LweSample *tempcarry2 = new_gate_bootstrapping_ciphertext_array(32, params);
                LweSample *tempcarry3 = new_gate_bootstrapping_ciphertext_array(32, params);
                LweSample *tempcarry4 = new_gate_bootstrapping_ciphertext_array(32, params);

                LweSample *twoscarry1 = new_gate_bootstrapping_ciphertext_array(32, params);
                LweSample *twoscarry2 = new_gate_bootstrapping_ciphertext_array(32, params);
                LweSample *twoscarry3 = new_gate_bootstrapping_ciphertext_array(32, params);
                LweSample *twoscarry4 = new_gate_bootstrapping_ciphertext_array(32, params);

                struct timeval start, end;
                double get_time;
                gettimeofday(&start, NULL);

                printf("Doing the homomorphic computation (128bit Subtraction)...\n");

                NOT(inverse1, ciphertext9, bk, 32);
                NOT(inverse2, ciphertext10, bk, 32);
                NOT(inverse3, ciphertext11, bk, 32);
                NOT(inverse4, ciphertext12, bk, 32);

                zero(tempStart, bk, 32);
                zero(tempcarry1, bk, 32);
                zero(tempcarry2, bk, 32);
                zero(tempcarry3, bk, 32);
                zero(tempcarry4, bk, 32);

                bootsCONSTANT(tempStart, 1, bk);

                add_addition(twosresult1, twoscarry1, inverse1, tempStart, tempcarry1, 32, bk);

                add_addition(twosresult2, twoscarry2, inverse2, tempcarry2, twoscarry1, 32, bk);
                add_addition(twosresult3, twoscarry3, inverse3, tempcarry3, twoscarry2, 32, bk);
                add_addition(twosresult4, twoscarry4, inverse4, tempcarry4, twoscarry3, 32, bk);

                ////////////////////////////////////

                zero(temp1, bk, 32);
                zero(temp3, bk, 32);
                zero(temp5, bk, 32);

                //making carry = 1
                bootsCONSTANT(temp2, 1, bk);
                bootsCONSTANT(temp4, 1, bk);
                bootsCONSTANT(temp6, 1, bk);

                //first 32 bit addition
                add_addition(result, carry1, ciphertext1, twosresult1, ciphertextcarry1, 32, bk);

                MPI_Barrier(MPI_COMM_WORLD);



                if (world_rank == 0)
                {
                    printf("Rank (i=0) %d\n", world_rank);
                    time_t t;
                    time(&t);
                    printf("current time is : %s", ctime(&t));

                    add_addition(result2, carry2, ciphertext2, twosresult2, temp1, 32, bk);

                    std::string name = "result2sub.data";
                    char const *fileparam = name.c_str();
                    FILE *answer1 = fopen(fileparam, "wb");
                    for (int i = 0; i < 32; i++)
                        export_gate_bootstrapping_ciphertext_toFile(answer1, &result2[i], params);
                    fclose(answer1);

                    std::string nameCarry1 = "carry2sub.data";
                    char const *fileparamCarry1 = nameCarry1.c_str();
                    FILE *answerCarry1 = fopen(fileparamCarry1, "wb");
                    for (int i = 0; i < 32; i++)
                        export_gate_bootstrapping_ciphertext_toFile(answerCarry1, &carry2[i], params);
                    fclose(answerCarry1);

                    cleanup(result2, 32, bk);
                    cleanup(carry2, 32, bk);
                }
                if (world_rank == 1)
                {
                    printf("Rank (i=1) %d\n", world_rank);
                    time_t t;
                    time(&t);
                    printf("current time is : %s", ctime(&t));

                    add_addition(result3, carry3, ciphertext2, twosresult2, temp2, 32, bk);

                    std::string name1 = "result3sub.data";
                    char const *fileparam1 = name1.c_str();
                    FILE *answer2 = fopen(fileparam1, "wb");
                    for (int i = 0; i < 32; i++)
                        export_gate_bootstrapping_ciphertext_toFile(answer2, &result3[i], params);
                    fclose(answer2);

                    std::string nameCarry2 = "carry3sub.data";
                    char const *fileparamCarry2 = nameCarry2.c_str();
                    FILE *answerCarry2 = fopen(fileparamCarry2, "wb");
                    for (int i = 0; i < 32; i++)
                        export_gate_bootstrapping_ciphertext_toFile(answerCarry2, &carry3[i], params);
                    fclose(answerCarry2);

                    cleanup(result3, 32, bk);
                    cleanup(carry3, 32, bk);
                }

                if (world_rank == 2)
                {
                    printf("Rank (i=2) %d\n", world_rank);
                    time_t t;
                    time(&t);
                    printf("current time is : %s", ctime(&t));

                    add_addition(result4, carry4, ciphertext3, twosresult3, temp3, 32, bk);

                    std::string name2 = "result4sub.data";
                    char const *fileparam2 = name2.c_str();
                    FILE *answer3 = fopen(fileparam2, "wb");
                    for (int i = 0; i < 32; i++)
                        export_gate_bootstrapping_ciphertext_toFile(answer3, &result4[i], params);
                    fclose(answer3);

                    std::string nameCarry3 = "carry4sub.data";
                    char const *fileparamCarry3 = nameCarry3.c_str();
                    FILE *answerCarry3 = fopen(fileparamCarry3, "wb");
                    for (int i = 0; i < 32; i++)
                        export_gate_bootstrapping_ciphertext_toFile(answerCarry3, &carry4[i], params);
                    fclose(answerCarry3);

                    cleanup(result4, 32, bk);
                    cleanup(carry4, 32, bk);
                }
                //carry = 1 result5
                if (world_rank == 3)
                {
                    printf("Rank (i=3) %d\n", world_rank);
                    time_t t;
                    time(&t);
                    printf("current time is : %s", ctime(&t));

                    add_addition(result5, carry5, ciphertext3, twosresult3, temp4, 32, bk);

                    std::string name3 = "result5sub.data";
                    char const *fileparam3 = name3.c_str();
                    FILE *answer4 = fopen(fileparam3, "wb");
                    for (int i = 0; i < 32; i++)
                        export_gate_bootstrapping_ciphertext_toFile(answer4, &result5[i], params);
                    fclose(answer4);

                    std::string nameCarry4 = "carry5sub.data";
                    char const *fileparamCarry4 = nameCarry4.c_str();
                    FILE *answerCarry4 = fopen(fileparamCarry4, "wb");
                    for (int i = 0; i < 32; i++)
                        export_gate_bootstrapping_ciphertext_toFile(answerCarry4, &carry5[i], params);
                    fclose(answerCarry4);

                    cleanup(result5, 32, bk);
                    cleanup(carry5, 32, bk);
                }

                if (world_rank == 4)
                {
                    printf("Rank (i=4) %d\n", world_rank);
                    time_t t;
                    time(&t);
                    printf("current time is : %s", ctime(&t));

                    add_addition(result6, carry6, ciphertext4, twosresult4, temp5, 32, bk);

                    std::string name4 = "result6sub.data";
                    char const *fileparam4 = name4.c_str();
                    FILE *answer5 = fopen(fileparam4, "wb");
                    for (int i = 0; i < 32; i++)
                        export_gate_bootstrapping_ciphertext_toFile(answer5, &result6[i], params);
                    fclose(answer5);

                    std::string nameCarry5 = "carry6sub.data";
                    char const *fileparamCarry5 = nameCarry5.c_str();
                    FILE *answerCarry5 = fopen(fileparamCarry5, "wb");
                    for (int i = 0; i < 32; i++)
                        export_gate_bootstrapping_ciphertext_toFile(answerCarry5, &carry6[i], params);
                    fclose(answerCarry5);

                    cleanup(result5, 32, bk);
                    cleanup(carry6, 32, bk);
                }
                //carry = 1 result7
                if (world_rank == 5)
                {
                    printf("Rank (i=5) %d\n", world_rank);
                    time_t t;
                    time(&t);
                    printf("current time is : %s", ctime(&t));

                    add_addition(result7, carry7, ciphertext4, twosresult4, temp6, 32, bk);

                    std::string name5 = "result7sub.data";
                    char const *fileparam5 = name5.c_str();
                    FILE *answer6 = fopen(fileparam5, "wb");
                    for (int i = 0; i < 32; i++)
                        export_gate_bootstrapping_ciphertext_toFile(answer6, &result7[i], params);
                    fclose(answer6);

                    std::string nameCarry6 = "carry7sub.data";
                    char const *fileparamCarry6 = nameCarry6.c_str();
                    FILE *answerCarry6 = fopen(fileparamCarry6, "wb");
                    for (int i = 0; i < 32; i++)
                        export_gate_bootstrapping_ciphertext_toFile(answerCarry6, &carry7[i], params);
                    fclose(answerCarry6);

                    cleanup(result5, 32, bk);
                    cleanup(carry7, 32, bk);
                }

                MPI_Barrier(MPI_COMM_WORLD);
                // Timings
                gettimeofday(&end, NULL);
                get_time = (end.tv_sec - start.tv_sec) + (end.tv_usec - start.tv_usec) * 1.0E-6;
                printf("Computation Time [MPI]: %lf[sec]\n", get_time);

                std::string name2 = "result2sub.data";
                char const *fileparam3 = name2.c_str();
                FILE *ans1 = fopen(fileparam3, "rb");
                for (int i = 0; i < 32; i++)
                    import_gate_bootstrapping_ciphertext_fromFile(ans1, &anstext1[i], params);
                fclose(ans1);

                std::string name3 = "result3sub.data";
                char const *fileparam4 = name3.c_str();
                FILE *ans2 = fopen(fileparam4, "rb");
                for (int i = 0; i < 32; i++)
                    import_gate_bootstrapping_ciphertext_fromFile(ans2, &anstext2[i], params);
                fclose(ans2);

                std::string name4 = "result4sub.data";
                char const *fileparam5 = name4.c_str();
                FILE *ans3 = fopen(fileparam5, "rb");
                for (int i = 0; i < 32; i++)
                    import_gate_bootstrapping_ciphertext_fromFile(ans3, &anstext3[i], params);
                fclose(ans3);

                std::string name5 = "result5sub.data";
                char const *fileparam6 = name5.c_str();
                FILE *ans4 = fopen(fileparam6, "rb");
                for (int i = 0; i < 32; i++)
                    import_gate_bootstrapping_ciphertext_fromFile(ans4, &anstext4[i], params);
                fclose(ans4);

                std::string name6 = "result6sub.data";
                char const *fileparam7 = name6.c_str();
                FILE *ans5 = fopen(fileparam7, "rb");
                for (int i = 0; i < 32; i++)
                    import_gate_bootstrapping_ciphertext_fromFile(ans5, &anstext5[i], params);
                fclose(ans5);

                std::string name7 = "result7sub.data";
                char const *fileparam8 = name7.c_str();
                FILE *ans6 = fopen(fileparam8, "rb");
                for (int i = 0; i < 32; i++)
                    import_gate_bootstrapping_ciphertext_fromFile(ans6, &anstext6[i], params);
                fclose(ans6);

     

                std::string nameImportCarry2 = "carry2sub.data";
                char const *fileparamImportCarry2 = nameImportCarry2.c_str();
                FILE *ansImportCarry1 = fopen(fileparamImportCarry2, "rb");
                for (int i = 0; i < 32; i++)
                    import_gate_bootstrapping_ciphertext_fromFile(ansImportCarry1, &anstextImportCarry2[i], params);
                fclose(ansImportCarry1);

                std::string nameImportCarry3 = "carry3sub.data";
                char const *fileparamImportCarry3 = nameImportCarry3.c_str();
                FILE *ansImportCarry3 = fopen(fileparamImportCarry3, "rb");
                for (int i = 0; i < 32; i++)
                    import_gate_bootstrapping_ciphertext_fromFile(ansImportCarry3, &anstextImportCarry3[i], params);
                fclose(ansImportCarry3);

                std::string nameImportCarry4 = "carry4sub.data";
                char const *fileparamImportCarry4 = nameImportCarry4.c_str();
                FILE *ansImportCarry4 = fopen(fileparamImportCarry4, "rb");
                for (int i = 0; i < 32; i++)
                    import_gate_bootstrapping_ciphertext_fromFile(ansImportCarry4, &anstextImportCarry4[i], params);
                fclose(ansImportCarry4);

                std::string nameImportCarry5 = "carry5sub.data";
                char const *fileparamImportCarry5 = nameImportCarry5.c_str();
                FILE *ansImportCarry5 = fopen(fileparamImportCarry5, "rb");
                for (int i = 0; i < 32; i++)
                    import_gate_bootstrapping_ciphertext_fromFile(ansImportCarry5, &anstextImportCarry5[i], params);
                fclose(ansImportCarry5);

                std::string nameImportCarry6 = "carry6sub.data";
                char const *fileparamImportCarry6 = nameImportCarry6.c_str();
                FILE *ansImportCarry6 = fopen(fileparamImportCarry6, "rb");
                for (int i = 0; i < 32; i++)
                    import_gate_bootstrapping_ciphertext_fromFile(ansImportCarry6, &anstextImportCarry6[i], params);
                fclose(ansImportCarry6);

                std::string nameImportCarry7 = "carry7sub.data";
                char const *fileparamImportCarry7 = nameImportCarry7.c_str();
                FILE *ansImportCarry7 = fopen(fileparamImportCarry7, "rb");
                for (int i = 0; i < 32; i++)
                    import_gate_bootstrapping_ciphertext_fromFile(ansImportCarry7, &anstextImportCarry7[i], params);
                fclose(ansImportCarry7);

                // export the result ciphertexts to a file (for the cloud)
                for (int i = 0; i < 32; i++) // result1 behind
                    export_gate_bootstrapping_ciphertext_toFile(answer_data, &result[i], params);
                for (int i = 0; i < 32; i++) // result1 behind
                    export_gate_bootstrapping_ciphertext_toFile(answer_data, &carry1[i], params);
                for (int i = 0; i < 32; i++) // result1 infront
                    export_gate_bootstrapping_ciphertext_toFile(answer_data, &anstextImportCarry2[i], params);
                for (int i = 0; i < 32; i++) // result1 infront
                    export_gate_bootstrapping_ciphertext_toFile(answer_data, &anstextImportCarry3[i], params);
                for (int i = 0; i < 32; i++) // result1 infront
                    export_gate_bootstrapping_ciphertext_toFile(answer_data, &anstextImportCarry4[i], params);
                for (int i = 0; i < 32; i++) // result1 infron
                    export_gate_bootstrapping_ciphertext_toFile(answer_data, &anstextImportCarry5[i], params);
                for (int i = 0; i < 32; i++) // result1 infront
                    export_gate_bootstrapping_ciphertext_toFile(answer_data, &anstextImportCarry6[i], params);
                for (int i = 0; i < 32; i++) // result1 infront
                    export_gate_bootstrapping_ciphertext_toFile(answer_data, &anstextImportCarry7[i], params);
                for (int i = 0; i < 32; i++) // result2 binary2
                    export_gate_bootstrapping_ciphertext_toFile(answer_data, &anstext1[i], params);
                for (int i = 0; i < 32; i++) // result3 binary3
                    export_gate_bootstrapping_ciphertext_toFile(answer_data, &anstext2[i], params);
                for (int i = 0; i < 32; i++) // 5
                    export_gate_bootstrapping_ciphertext_toFile(answer_data, &anstext3[i], params);
                for (int i = 0; i < 32; i++) // 6
                    export_gate_bootstrapping_ciphertext_toFile(answer_data, &anstext4[i], params);
                for (int i = 0; i < 32; i++) // 7
                    export_gate_bootstrapping_ciphertext_toFile(answer_data, &anstext5[i], params);
                for (int i = 0; i < 32; i++) // 8
                    export_gate_bootstrapping_ciphertext_toFile(answer_data, &anstext6[i], params);
                for (int i = 0; i < 32; i++) // carry
                    export_gate_bootstrapping_ciphertext_toFile(answer_data, &ciphertextcarry1[i], params);
                fclose(answer_data);

                MPI_Barrier(MPI_COMM_WORLD);

                printf("writing the answer to file...\n");


                printf("writing the answer to file...for %d, Calculating\n", world_rank);
                MPI_Barrier(MPI_COMM_WORLD);

                double time_taken;
                MPI_Reduce(&get_time, &time_taken, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

                FILE *t_file;
                if (world_rank == 0)
                {
                    printf("avg = %lf\n", time_taken / 32);
                    t_file = fopen(T_FILE, "a");
                    fprintf(t_file, "%lf\n", time_taken / 32);
                    fclose(t_file);
                }

                delete_gate_bootstrapping_ciphertext_array(32, result);
                delete_gate_bootstrapping_ciphertext_array(32, result2);
                delete_gate_bootstrapping_ciphertext_array(32, result3);
                delete_gate_bootstrapping_ciphertext_array(32, result4);
                delete_gate_bootstrapping_ciphertext_array(32, result5);
                delete_gate_bootstrapping_ciphertext_array(32, result6);
                delete_gate_bootstrapping_ciphertext_array(32, result7);
                delete_gate_bootstrapping_ciphertext_array(32, carry1);
                delete_gate_bootstrapping_ciphertext_array(32, carry2);
                delete_gate_bootstrapping_ciphertext_array(32, carry3);
                delete_gate_bootstrapping_ciphertext_array(32, carry4);
                delete_gate_bootstrapping_ciphertext_array(32, carry5);
                delete_gate_bootstrapping_ciphertext_array(32, carry6);
                delete_gate_bootstrapping_ciphertext_array(32, carry7);
                delete_gate_bootstrapping_ciphertext_array(32, ciphertextbit1);
                delete_gate_bootstrapping_ciphertext_array(32, ciphertextnegative1);
                delete_gate_bootstrapping_ciphertext_array(32, ciphertextbit2);
                delete_gate_bootstrapping_ciphertext_array(32, ciphertextnegative2);
                delete_gate_bootstrapping_ciphertext_array(32, ciphertext1);
                delete_gate_bootstrapping_ciphertext_array(32, ciphertext2);
                delete_gate_bootstrapping_ciphertext_array(32, ciphertext3);
                delete_gate_bootstrapping_ciphertext_array(32, ciphertext4);
                delete_gate_bootstrapping_ciphertext_array(32, ciphertext9);
                delete_gate_bootstrapping_ciphertext_array(32, ciphertext10);
                delete_gate_bootstrapping_ciphertext_array(32, ciphertext11);
                delete_gate_bootstrapping_ciphertext_array(32, ciphertext12);
                delete_gate_bootstrapping_ciphertext_array(32, temp1);
                delete_gate_bootstrapping_ciphertext_array(32, temp2);
                delete_gate_bootstrapping_ciphertext_array(32, temp3);
                delete_gate_bootstrapping_ciphertext_array(32, temp4);
                delete_gate_bootstrapping_ciphertext_array(32, anstext1);
                delete_gate_bootstrapping_ciphertext_array(32, anstext2);
                delete_gate_bootstrapping_ciphertext_array(32, anstext3);
                delete_gate_bootstrapping_ciphertext_array(32, anstext4);
                delete_gate_bootstrapping_ciphertext_array(32, anstext5);
                delete_gate_bootstrapping_ciphertext_array(32, anstext6);
                delete_gate_bootstrapping_ciphertext_array(32, anstextImportCarry2);
                delete_gate_bootstrapping_ciphertext_array(32, anstextImportCarry3);
                delete_gate_bootstrapping_ciphertext_array(32, anstextImportCarry4);
                delete_gate_bootstrapping_ciphertext_array(32, anstextImportCarry5);
                delete_gate_bootstrapping_ciphertext_array(32, anstextImportCarry6);
                delete_gate_bootstrapping_ciphertext_array(32, anstextImportCarry7);

                delete_gate_bootstrapping_ciphertext_array(32, tempStart);
                delete_gate_bootstrapping_ciphertext_array(32, inverse1);
                delete_gate_bootstrapping_ciphertext_array(32, inverse2);
                delete_gate_bootstrapping_ciphertext_array(32, inverse3);
                delete_gate_bootstrapping_ciphertext_array(32, inverse4);

                delete_gate_bootstrapping_ciphertext_array(32, tempcarry1);
                delete_gate_bootstrapping_ciphertext_array(32, tempcarry2);
                delete_gate_bootstrapping_ciphertext_array(32, tempcarry3);
                delete_gate_bootstrapping_ciphertext_array(32, tempcarry4);

                delete_gate_bootstrapping_ciphertext_array(32, twosresult1);
                delete_gate_bootstrapping_ciphertext_array(32, twosresult2);
                delete_gate_bootstrapping_ciphertext_array(32, twosresult3);
                delete_gate_bootstrapping_ciphertext_array(32, twosresult4);

                delete_gate_bootstrapping_ciphertext_array(32, twoscarry1);
                delete_gate_bootstrapping_ciphertext_array(32, twoscarry2);
                delete_gate_bootstrapping_ciphertext_array(32, twoscarry3);
                delete_gate_bootstrapping_ciphertext_array(32, twoscarry4);

                delete_gate_bootstrapping_cloud_keyset(bk);
                delete_gate_bootstrapping_secret_keyset(nbitkey);

                MPI_Finalize();
            }
        }

        else
        {
            if (int_op == 2)
            {
                std::cout << int_bit << " bit Subtraction computation"
                          << "\n";
            }
            else
            {
                std::cout << int_bit << " bit Addition computation with 1st value negative"
                          << "\n";
            }

            if (int_bit == 32)
            {
                MPI_Init(&argc, &argv);
                int world_rank;
                MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
                int world_size;
                MPI_Comm_size(MPI_COMM_WORLD, &world_size);

                printf("Doing the homomorphic computation...\n");

                LweSample *temp = new_gate_bootstrapping_ciphertext_array(32, params);
                LweSample *inverse1 = new_gate_bootstrapping_ciphertext_array(32, params);

                LweSample *tempcarry1 = new_gate_bootstrapping_ciphertext_array(32, params);

                LweSample *twosresult1 = new_gate_bootstrapping_ciphertext_array(32, params);

                LweSample *twoscarry1 = new_gate_bootstrapping_ciphertext_array(32, params);

                struct timeval start, end;
                double get_time;

                gettimeofday(&start, NULL);

                NOT(inverse1, ciphertext1, bk, 32);


                zero(temp, bk, 32);
                zero(tempcarry1, bk, 32);

                bootsCONSTANT(temp, 1, bk);

                add_addition(twosresult1, twoscarry1, inverse1, temp, tempcarry1, 32, bk);

                LweSample *result1 = new_gate_bootstrapping_ciphertext_array(32, params);

                LweSample *carry1 = new_gate_bootstrapping_ciphertext_array(32, params);

                add_addition(result1, carry1, ciphertext9, twosresult1, ciphertextcarry1, 32, bk);

                gettimeofday(&end, NULL);

                get_time = (end.tv_sec - start.tv_sec) + (end.tv_usec - start.tv_usec) * 1.0E-6;
                printf("Computation Time: %lf[sec]\n", get_time);

                printf("writing the answer to file...\n");

                for (int i = 0; i < 32; i++) //result1
                    export_gate_bootstrapping_ciphertext_toFile(answer_data, &result1[i], params);
                for (int i = 0; i < 32; i++) // 2
                    export_gate_bootstrapping_ciphertext_toFile(answer_data, &ciphertextcarry1[i], params);
                for (int i = 0; i < 32; i++) // 3
                    export_gate_bootstrapping_ciphertext_toFile(answer_data, &ciphertextcarry1[i], params);
                for (int i = 0; i < 32; i++) // 4
                    export_gate_bootstrapping_ciphertext_toFile(answer_data, &ciphertextcarry1[i], params);
                for (int i = 0; i < 32; i++) // 5
                    export_gate_bootstrapping_ciphertext_toFile(answer_data, &ciphertextcarry1[i], params);
                for (int i = 0; i < 32; i++) // 6
                    export_gate_bootstrapping_ciphertext_toFile(answer_data, &ciphertextcarry1[i], params);
                for (int i = 0; i < 32; i++) // 7
                    export_gate_bootstrapping_ciphertext_toFile(answer_data, &ciphertextcarry1[i], params);
                for (int i = 0; i < 32; i++) // 8
                    export_gate_bootstrapping_ciphertext_toFile(answer_data, &ciphertextcarry1[i], params);
                for (int i = 0; i < 32; i++) // carry
                    export_gate_bootstrapping_ciphertext_toFile(answer_data, &ciphertextcarry1[i], params);
                fclose(answer_data);

                delete_gate_bootstrapping_ciphertext_array(32, temp);
                delete_gate_bootstrapping_ciphertext_array(32, inverse1);

                delete_gate_bootstrapping_ciphertext_array(32, tempcarry1);

                delete_gate_bootstrapping_ciphertext_array(32, twosresult1);

                delete_gate_bootstrapping_ciphertext_array(32, twoscarry1);

                delete_gate_bootstrapping_ciphertext_array(32, carry1);

                delete_gate_bootstrapping_ciphertext_array(32, result1);

                delete_gate_bootstrapping_ciphertext_array(32, ciphertextbit1);
                delete_gate_bootstrapping_ciphertext_array(32, ciphertextnegative1);
                delete_gate_bootstrapping_ciphertext_array(32, ciphertextbit2);
                delete_gate_bootstrapping_ciphertext_array(32, ciphertextnegative2);
                delete_gate_bootstrapping_ciphertext_array(32, ciphertext1);
                delete_gate_bootstrapping_ciphertext_array(32, ciphertext9);
                delete_gate_bootstrapping_ciphertext_array(32, ciphertextcarry1);

                delete_gate_bootstrapping_cloud_keyset(bk);
                delete_gate_bootstrapping_secret_keyset(nbitkey);

                MPI_Finalize();
            }
            else if (int_bit == 64)
            {
                MPI_Init(&argc, &argv);
                int world_rank;
                MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
                int world_size;
                MPI_Comm_size(MPI_COMM_WORLD, &world_size);

                printf("Doing the homomorphic computation (64bit Subtraction) [(-A)+B] ...\n");

                LweSample *temp = new_gate_bootstrapping_ciphertext_array(32, params);

                LweSample *inverse1 = new_gate_bootstrapping_ciphertext_array(32, params);
                LweSample *inverse2 = new_gate_bootstrapping_ciphertext_array(32, params);

                LweSample *tempcarry1 = new_gate_bootstrapping_ciphertext_array(32, params);
                LweSample *tempcarry2 = new_gate_bootstrapping_ciphertext_array(32, params);

                LweSample *twosresult1 = new_gate_bootstrapping_ciphertext_array(32, params);
                LweSample *twosresult2 = new_gate_bootstrapping_ciphertext_array(32, params);

                LweSample *twoscarry1 = new_gate_bootstrapping_ciphertext_array(32, params);
                LweSample *twoscarry2 = new_gate_bootstrapping_ciphertext_array(32, params);

                struct timeval start, end;
                double get_time;
                gettimeofday(&start, NULL);


                NOT(inverse1, ciphertext1, bk, 32);
                NOT(inverse2, ciphertext2, bk, 32);

                zero(temp, bk, 32);

                zero(tempcarry1, bk, 32);
                zero(tempcarry2, bk, 32);

                bootsCONSTANT(temp, 1, bk);

   
                add_addition(twosresult1, twoscarry1, inverse1, temp, tempcarry1, 32, bk);

                add_addition(twosresult2, twoscarry2, inverse2, tempcarry2, twoscarry1, 32, bk);

                LweSample *result1 = new_gate_bootstrapping_ciphertext_array(32, params);
                LweSample *result2 = new_gate_bootstrapping_ciphertext_array(32, params);
                LweSample *result3 = new_gate_bootstrapping_ciphertext_array(32, params);

                LweSample *carry1 = new_gate_bootstrapping_ciphertext_array(32, params);
                LweSample *carry2 = new_gate_bootstrapping_ciphertext_array(32, params);


                add_addition(result1, carry1, ciphertext9, twosresult1, ciphertextcarry1, 32, bk);

                LweSample *makeitzero = new_gate_bootstrapping_ciphertext_array(32, params);
                LweSample *makeitone = new_gate_bootstrapping_ciphertext_array(32, params);

                LweSample *anstext1 = new_gate_bootstrapping_ciphertext_array(32, params);
                LweSample *anstext2 = new_gate_bootstrapping_ciphertext_array(32, params);

                zero(makeitzero, bk, 32);
                bootsCONSTANT(makeitone, 1, bk);

                MPI_Barrier(MPI_COMM_WORLD);

                if (world_rank == 0)
                {
                    printf("Rank (i=0) %d\n", world_rank);
                    time_t t;
                    time(&t);
                    printf("current time is : %s", ctime(&t));

                    add_addition(result2, carry2, ciphertext10, twosresult2, makeitzero, 32, bk);

                    std::string name = "result2sub.data";
                    char const *fileparam = name.c_str();
                    FILE *answer1 = fopen(fileparam, "wb");
                    for (int i = 0; i < 32; i++)
                        export_gate_bootstrapping_ciphertext_toFile(answer1, &result2[i], params);
                    fclose(answer1);

                    cleanup(result2, 32, bk);
                }

                //carry = 1
                if (world_rank == 1)
                {
                    printf("Rank (i=0) %d\n", world_rank);
                    time_t t;
                    time(&t);
                    printf("current time is : %s", ctime(&t));

                    add_addition(result3, carry2, ciphertext10, twosresult2, makeitone, 32, bk);

                    std::string name1 = "result3sub.data";
                    char const *fileparam1 = name1.c_str();
                    FILE *answer2 = fopen(fileparam1, "wb");
                    for (int i = 0; i < 32; i++)
                        export_gate_bootstrapping_ciphertext_toFile(answer2, &result3[i], params);
                    fclose(answer2);

                    cleanup(result3, 32, bk);
                }

                MPI_Barrier(MPI_COMM_WORLD);

                gettimeofday(&end, NULL);
                get_time = (end.tv_sec - start.tv_sec) + (end.tv_usec - start.tv_usec) * 1.0E-6;
                printf("Computation Time [MPI]: %lf[sec]\n", get_time);

                std::string name2 = "result2sub.data";
                char const *fileparam3 = name2.c_str();
                FILE *ans1 = fopen(fileparam3, "rb");
                for (int i = 0; i < 32; i++)
                    import_gate_bootstrapping_ciphertext_fromFile(ans1, &anstext1[i], params);
                fclose(ans1);

                std::string name3 = "result3sub.data";
                char const *fileparam4 = name3.c_str();
                FILE *ans2 = fopen(fileparam4, "rb");
                for (int i = 0; i < 32; i++)
                    import_gate_bootstrapping_ciphertext_fromFile(ans2, &anstext2[i], params);
                fclose(ans2);

                printf("Exporting 64bit Sub... \n");

                printf("writing the answer to file...\n");

                //export the 32 ciphertexts to a file (for the cloud)
                for (int i = 0; i < 32; i++) //result1
                    export_gate_bootstrapping_ciphertext_toFile(answer_data, &result1[i], params);
                for (int i = 0; i < 32; i++) // 2
                    export_gate_bootstrapping_ciphertext_toFile(answer_data, &carry1[i], params);
                for (int i = 0; i < 32; i++) // 3
                    export_gate_bootstrapping_ciphertext_toFile(answer_data, &anstext1[i], params);
                for (int i = 0; i < 32; i++) // 4
                    export_gate_bootstrapping_ciphertext_toFile(answer_data, &anstext2[i], params);
                for (int i = 0; i < 32; i++) // 5
                    export_gate_bootstrapping_ciphertext_toFile(answer_data, &ciphertextcarry1[i], params);
                for (int i = 0; i < 32; i++) // 6
                    export_gate_bootstrapping_ciphertext_toFile(answer_data, &ciphertextcarry1[i], params);
                for (int i = 0; i < 32; i++) // 7
                    export_gate_bootstrapping_ciphertext_toFile(answer_data, &ciphertextcarry1[i], params);
                for (int i = 0; i < 32; i++) // 8
                    export_gate_bootstrapping_ciphertext_toFile(answer_data, &ciphertextcarry1[i], params);
                for (int i = 0; i < 32; i++) // carry
                    export_gate_bootstrapping_ciphertext_toFile(answer_data, &ciphertextcarry1[i], params);
                fclose(answer_data);

                delete_gate_bootstrapping_ciphertext_array(32, temp);
                delete_gate_bootstrapping_ciphertext_array(32, inverse1);
                delete_gate_bootstrapping_ciphertext_array(32, inverse2);
                delete_gate_bootstrapping_ciphertext_array(32, tempcarry1);
                delete_gate_bootstrapping_ciphertext_array(32, tempcarry2);
                delete_gate_bootstrapping_ciphertext_array(32, twosresult1);
                delete_gate_bootstrapping_ciphertext_array(32, twosresult2);
                delete_gate_bootstrapping_ciphertext_array(32, twoscarry1);
                delete_gate_bootstrapping_ciphertext_array(32, twoscarry2);
                delete_gate_bootstrapping_ciphertext_array(32, carry1);
                delete_gate_bootstrapping_ciphertext_array(32, carry2);
                delete_gate_bootstrapping_ciphertext_array(32, result1);
                delete_gate_bootstrapping_ciphertext_array(32, result2);
                delete_gate_bootstrapping_ciphertext_array(32, result3);
                delete_gate_bootstrapping_ciphertext_array(32, anstext1);
                delete_gate_bootstrapping_ciphertext_array(32, anstext2);

                delete_gate_bootstrapping_ciphertext_array(32, ciphertextbit1);
                delete_gate_bootstrapping_ciphertext_array(32, ciphertextnegative1);
                delete_gate_bootstrapping_ciphertext_array(32, ciphertextbit2);
                delete_gate_bootstrapping_ciphertext_array(32, ciphertextnegative2);
                delete_gate_bootstrapping_ciphertext_array(32, ciphertext1);
                delete_gate_bootstrapping_ciphertext_array(32, ciphertext2);
                delete_gate_bootstrapping_ciphertext_array(32, ciphertext9);
                delete_gate_bootstrapping_ciphertext_array(32, ciphertext10);
                delete_gate_bootstrapping_ciphertext_array(32, ciphertextcarry1);

                delete_gate_bootstrapping_cloud_keyset(bk);
                delete_gate_bootstrapping_secret_keyset(nbitkey);

                MPI_Finalize();
            }

            else if (int_bit == 128)
            {
                printf("MPI Subtraction is being computed 128 bit [mpicloud.c] [(-A)+B] \n");
                MPI_Init(&argc, &argv);
                int world_rank;
                MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
                int world_size;
                MPI_Comm_size(MPI_COMM_WORLD, &world_size);


                printf("Hello World from process %d of %d\n", world_rank, world_size);
                MPI_Barrier(MPI_COMM_WORLD);

                LweSample *result = new_gate_bootstrapping_ciphertext_array(32, params);
                LweSample *result2 = new_gate_bootstrapping_ciphertext_array(32, params);
                LweSample *result3 = new_gate_bootstrapping_ciphertext_array(32, params);
                LweSample *result4 = new_gate_bootstrapping_ciphertext_array(32, params);
                LweSample *result5 = new_gate_bootstrapping_ciphertext_array(32, params);
                LweSample *result6 = new_gate_bootstrapping_ciphertext_array(32, params);
                LweSample *result7 = new_gate_bootstrapping_ciphertext_array(32, params);

                LweSample *carry1 = new_gate_bootstrapping_ciphertext_array(32, params);
                LweSample *carry2 = new_gate_bootstrapping_ciphertext_array(32, params);
                LweSample *carry3 = new_gate_bootstrapping_ciphertext_array(32, params);
                LweSample *carry4 = new_gate_bootstrapping_ciphertext_array(32, params);
                LweSample *carry5 = new_gate_bootstrapping_ciphertext_array(32, params);
                LweSample *carry6 = new_gate_bootstrapping_ciphertext_array(32, params);
                LweSample *carry7 = new_gate_bootstrapping_ciphertext_array(32, params);

                LweSample *temp1 = new_gate_bootstrapping_ciphertext_array(32, params);
                LweSample *temp2 = new_gate_bootstrapping_ciphertext_array(32, params);
                LweSample *temp3 = new_gate_bootstrapping_ciphertext_array(32, params);
                LweSample *temp4 = new_gate_bootstrapping_ciphertext_array(32, params);
                LweSample *temp5 = new_gate_bootstrapping_ciphertext_array(32, params);
                LweSample *temp6 = new_gate_bootstrapping_ciphertext_array(32, params);

                LweSample *anstext1 = new_gate_bootstrapping_ciphertext_array(32, params);
                LweSample *anstext2 = new_gate_bootstrapping_ciphertext_array(32, params);
                LweSample *anstext3 = new_gate_bootstrapping_ciphertext_array(32, params);
                LweSample *anstext4 = new_gate_bootstrapping_ciphertext_array(32, params);
                LweSample *anstext5 = new_gate_bootstrapping_ciphertext_array(32, params);
                LweSample *anstext6 = new_gate_bootstrapping_ciphertext_array(32, params);

                LweSample *anstextImportCarry2 = new_gate_bootstrapping_ciphertext_array(32, params);
                LweSample *anstextImportCarry3 = new_gate_bootstrapping_ciphertext_array(32, params);
                LweSample *anstextImportCarry4 = new_gate_bootstrapping_ciphertext_array(32, params);
                LweSample *anstextImportCarry5 = new_gate_bootstrapping_ciphertext_array(32, params);
                LweSample *anstextImportCarry6 = new_gate_bootstrapping_ciphertext_array(32, params);
                LweSample *anstextImportCarry7 = new_gate_bootstrapping_ciphertext_array(32, params);

                LweSample *tempStart = new_gate_bootstrapping_ciphertext_array(32, params);

                LweSample *inverse1 = new_gate_bootstrapping_ciphertext_array(32, params);
                LweSample *inverse2 = new_gate_bootstrapping_ciphertext_array(32, params);
                LweSample *inverse3 = new_gate_bootstrapping_ciphertext_array(32, params);
                LweSample *inverse4 = new_gate_bootstrapping_ciphertext_array(32, params);

                LweSample *twosresult1 = new_gate_bootstrapping_ciphertext_array(32, params);
                LweSample *twosresult2 = new_gate_bootstrapping_ciphertext_array(32, params);
                LweSample *twosresult3 = new_gate_bootstrapping_ciphertext_array(32, params);
                LweSample *twosresult4 = new_gate_bootstrapping_ciphertext_array(32, params);

                LweSample *tempcarry1 = new_gate_bootstrapping_ciphertext_array(32, params);
                LweSample *tempcarry2 = new_gate_bootstrapping_ciphertext_array(32, params);
                LweSample *tempcarry3 = new_gate_bootstrapping_ciphertext_array(32, params);
                LweSample *tempcarry4 = new_gate_bootstrapping_ciphertext_array(32, params);

                LweSample *twoscarry1 = new_gate_bootstrapping_ciphertext_array(32, params);
                LweSample *twoscarry2 = new_gate_bootstrapping_ciphertext_array(32, params);
                LweSample *twoscarry3 = new_gate_bootstrapping_ciphertext_array(32, params);
                LweSample *twoscarry4 = new_gate_bootstrapping_ciphertext_array(32, params);

                struct timeval start, end;
                double get_time;
                gettimeofday(&start, NULL);

                printf("Doing the homomorphic computation (128bit Subtraction)...\n");

                NOT(inverse1, ciphertext1, bk, 32);
                NOT(inverse2, ciphertext2, bk, 32);
                NOT(inverse3, ciphertext3, bk, 32);
                NOT(inverse4, ciphertext4, bk, 32);

                //iniailize tempcarry and temp carry to 0
                zero(tempStart, bk, 32);
                zero(tempcarry1, bk, 32);
                zero(tempcarry2, bk, 32);
                zero(tempcarry3, bk, 32);
                zero(tempcarry4, bk, 32);

                bootsCONSTANT(tempStart, 1, bk);

                add_addition(twosresult1, twoscarry1, inverse1, tempStart, tempcarry1, 32, bk);

                add_addition(twosresult2, twoscarry2, inverse2, tempcarry2, twoscarry1, 32, bk);
                add_addition(twosresult3, twoscarry3, inverse3, tempcarry3, twoscarry2, 32, bk);
                add_addition(twosresult4, twoscarry4, inverse4, tempcarry4, twoscarry3, 32, bk);


                zero(temp1, bk, 32);
                zero(temp3, bk, 32);
                zero(temp5, bk, 32);

                bootsCONSTANT(temp2, 1, bk);
                bootsCONSTANT(temp4, 1, bk);
                bootsCONSTANT(temp6, 1, bk);

                add_addition(result, carry1, ciphertext9, twosresult1, ciphertextcarry1, 32, bk);

                MPI_Barrier(MPI_COMM_WORLD);

        
                //carry = 0 result2
                if (world_rank == 0)
                {
                    printf("Rank (i=0) %d\n", world_rank);
                    time_t t;
                    time(&t);
                    printf("current time is : %s", ctime(&t));

                    add_addition(result2, carry2, ciphertext10, twosresult2, temp1, 32, bk);

                    std::string name = "result2sub1.data";
                    char const *fileparam = name.c_str();
                    FILE *answer1 = fopen(fileparam, "wb");
                    for (int i = 0; i < 32; i++)
                        export_gate_bootstrapping_ciphertext_toFile(answer1, &result2[i], params);
                    fclose(answer1);

                    std::string nameCarry1 = "carry2sub1.data";
                    char const *fileparamCarry1 = nameCarry1.c_str();
                    FILE *answerCarry1 = fopen(fileparamCarry1, "wb");
                    for (int i = 0; i < 32; i++)
                        export_gate_bootstrapping_ciphertext_toFile(answerCarry1, &carry2[i], params);
                    fclose(answerCarry1);

                    cleanup(result2, 32, bk);
                    cleanup(carry2, 32, bk);
                }
                //carry == 1 result3
                if (world_rank == 1)
                {
                    printf("Rank (i=1) %d\n", world_rank);
                    time_t t;
                    time(&t);
                    printf("current time is : %s", ctime(&t));

                    add_addition(result3, carry3, ciphertext10, twosresult2, temp2, 32, bk);

                    std::string name1 = "result3sub1.data";
                    char const *fileparam1 = name1.c_str();
                    FILE *answer2 = fopen(fileparam1, "wb");
                    for (int i = 0; i < 32; i++)
                        export_gate_bootstrapping_ciphertext_toFile(answer2, &result3[i], params);
                    fclose(answer2);

                    std::string nameCarry2 = "carry3sub1.data";
                    char const *fileparamCarry2 = nameCarry2.c_str();
                    FILE *answerCarry2 = fopen(fileparamCarry2, "wb");
                    for (int i = 0; i < 32; i++)
                        export_gate_bootstrapping_ciphertext_toFile(answerCarry2, &carry3[i], params);
                    fclose(answerCarry2);

                    cleanup(result3, 32, bk);
                    cleanup(carry3, 32, bk);
                }

                //carry = 0 result4
                if (world_rank == 2)
                {
                    printf("Rank (i=2) %d\n", world_rank);
                    time_t t;
                    time(&t);
                    printf("current time is : %s", ctime(&t));

                    add_addition(result4, carry4, ciphertext11, twosresult3, temp3, 32, bk);

                    std::string name2 = "result4sub1.data";
                    char const *fileparam2 = name2.c_str();
                    FILE *answer3 = fopen(fileparam2, "wb");
                    for (int i = 0; i < 32; i++)
                        export_gate_bootstrapping_ciphertext_toFile(answer3, &result4[i], params);
                    fclose(answer3);

                    std::string nameCarry3 = "carry4sub1.data";
                    char const *fileparamCarry3 = nameCarry3.c_str();
                    FILE *answerCarry3 = fopen(fileparamCarry3, "wb");
                    for (int i = 0; i < 32; i++)
                        export_gate_bootstrapping_ciphertext_toFile(answerCarry3, &carry4[i], params);
                    fclose(answerCarry3);

                    cleanup(result4, 32, bk);
                    cleanup(carry4, 32, bk);
                }
                //carry = 1 result5
                if (world_rank == 3)
                {
                    printf("Rank (i=3) %d\n", world_rank);
                    time_t t;
                    time(&t);
                    printf("current time is : %s", ctime(&t));

                    add_addition(result5, carry5, ciphertext11, twosresult3, temp4, 32, bk);

                    std::string name3 = "result5sub1.data";
                    char const *fileparam3 = name3.c_str();
                    FILE *answer4 = fopen(fileparam3, "wb");
                    for (int i = 0; i < 32; i++)
                        export_gate_bootstrapping_ciphertext_toFile(answer4, &result5[i], params);
                    fclose(answer4);

                    std::string nameCarry4 = "carry5sub1.data";
                    char const *fileparamCarry4 = nameCarry4.c_str();
                    FILE *answerCarry4 = fopen(fileparamCarry4, "wb");
                    for (int i = 0; i < 32; i++)
                        export_gate_bootstrapping_ciphertext_toFile(answerCarry4, &carry5[i], params);
                    fclose(answerCarry4);

                    cleanup(result5, 32, bk);
                    cleanup(carry5, 32, bk);
                }

                //carry = 0 result6
                if (world_rank == 4)
                {
                    printf("Rank (i=4) %d\n", world_rank);
                    time_t t;
                    time(&t);
                    printf("current time is : %s", ctime(&t));

                    add_addition(result6, carry6, ciphertext12, twosresult4, temp5, 32, bk);

                    std::string name4 = "result6sub1.data";
                    char const *fileparam4 = name4.c_str();
                    FILE *answer5 = fopen(fileparam4, "wb");
                    for (int i = 0; i < 32; i++)
                        export_gate_bootstrapping_ciphertext_toFile(answer5, &result6[i], params);
                    fclose(answer5);

                    std::string nameCarry5 = "carry6sub1.data";
                    char const *fileparamCarry5 = nameCarry5.c_str();
                    FILE *answerCarry5 = fopen(fileparamCarry5, "wb");
                    for (int i = 0; i < 32; i++)
                        export_gate_bootstrapping_ciphertext_toFile(answerCarry5, &carry6[i], params);
                    fclose(answerCarry5);

                    cleanup(result5, 32, bk);
                    cleanup(carry6, 32, bk);
                }
                //carry = 1 result7
                if (world_rank == 5)
                {
                    printf("Rank (i=5) %d\n", world_rank);
                    time_t t;
                    time(&t);
                    printf("current time is : %s", ctime(&t));

                    add_addition(result7, carry7, ciphertext12, twosresult4, temp6, 32, bk);

                    std::string name5 = "result7sub1.data";
                    char const *fileparam5 = name5.c_str();
                    FILE *answer6 = fopen(fileparam5, "wb");
                    for (int i = 0; i < 32; i++)
                        export_gate_bootstrapping_ciphertext_toFile(answer6, &result7[i], params);
                    fclose(answer6);

                    std::string nameCarry6 = "carry7sub1.data";
                    char const *fileparamCarry6 = nameCarry6.c_str();
                    FILE *answerCarry6 = fopen(fileparamCarry6, "wb");
                    for (int i = 0; i < 32; i++)
                        export_gate_bootstrapping_ciphertext_toFile(answerCarry6, &carry7[i], params);
                    fclose(answerCarry6);

                    cleanup(result5, 32, bk);
                    cleanup(carry7, 32, bk);
                }

                MPI_Barrier(MPI_COMM_WORLD);
                // Timings
                gettimeofday(&end, NULL);
                get_time = (end.tv_sec - start.tv_sec) + (end.tv_usec - start.tv_usec) * 1.0E-6;
                printf("Computation Time [MPI]: %lf[sec]\n", get_time);

                std::string name2 = "result2sub1.data";
                char const *fileparam3 = name2.c_str();
                FILE *ans1 = fopen(fileparam3, "rb");
                for (int i = 0; i < 32; i++)
                    import_gate_bootstrapping_ciphertext_fromFile(ans1, &anstext1[i], params);
                fclose(ans1);

                std::string name3 = "result3sub1.data";
                char const *fileparam4 = name3.c_str();
                FILE *ans2 = fopen(fileparam4, "rb");
                for (int i = 0; i < 32; i++)
                    import_gate_bootstrapping_ciphertext_fromFile(ans2, &anstext2[i], params);
                fclose(ans2);

                std::string name4 = "result4sub1.data";
                char const *fileparam5 = name4.c_str();
                FILE *ans3 = fopen(fileparam5, "rb");
                for (int i = 0; i < 32; i++)
                    import_gate_bootstrapping_ciphertext_fromFile(ans3, &anstext3[i], params);
                fclose(ans3);

                std::string name5 = "result5sub1.data";
                char const *fileparam6 = name5.c_str();
                FILE *ans4 = fopen(fileparam6, "rb");
                for (int i = 0; i < 32; i++)
                    import_gate_bootstrapping_ciphertext_fromFile(ans4, &anstext4[i], params);
                fclose(ans4);

                std::string name6 = "result6sub1.data";
                char const *fileparam7 = name6.c_str();
                FILE *ans5 = fopen(fileparam7, "rb");
                for (int i = 0; i < 32; i++)
                    import_gate_bootstrapping_ciphertext_fromFile(ans5, &anstext5[i], params);
                fclose(ans5);

                std::string name7 = "result7sub1.data";
                char const *fileparam8 = name7.c_str();
                FILE *ans6 = fopen(fileparam8, "rb");
                for (int i = 0; i < 32; i++)
                    import_gate_bootstrapping_ciphertext_fromFile(ans6, &anstext6[i], params);
                fclose(ans6);

                ////////////////////////////////////////////////
                //Carry////////////////////////////////////////
                //////////////////////////////////////////////

                std::string nameImportCarry2 = "carry2sub1.data";
                char const *fileparamImportCarry2 = nameImportCarry2.c_str();
                FILE *ansImportCarry1 = fopen(fileparamImportCarry2, "rb");
                for (int i = 0; i < 32; i++)
                    import_gate_bootstrapping_ciphertext_fromFile(ansImportCarry1, &anstextImportCarry2[i], params);
                fclose(ansImportCarry1);

                std::string nameImportCarry3 = "carry3sub1.data";
                char const *fileparamImportCarry3 = nameImportCarry3.c_str();
                FILE *ansImportCarry3 = fopen(fileparamImportCarry3, "rb");
                for (int i = 0; i < 32; i++)
                    import_gate_bootstrapping_ciphertext_fromFile(ansImportCarry3, &anstextImportCarry3[i], params);
                fclose(ansImportCarry3);

                std::string nameImportCarry4 = "carry4sub1.data";
                char const *fileparamImportCarry4 = nameImportCarry4.c_str();
                FILE *ansImportCarry4 = fopen(fileparamImportCarry4, "rb");
                for (int i = 0; i < 32; i++)
                    import_gate_bootstrapping_ciphertext_fromFile(ansImportCarry4, &anstextImportCarry4[i], params);
                fclose(ansImportCarry4);

                std::string nameImportCarry5 = "carry5sub1.data";
                char const *fileparamImportCarry5 = nameImportCarry5.c_str();
                FILE *ansImportCarry5 = fopen(fileparamImportCarry5, "rb");
                for (int i = 0; i < 32; i++)
                    import_gate_bootstrapping_ciphertext_fromFile(ansImportCarry5, &anstextImportCarry5[i], params);
                fclose(ansImportCarry5);

                std::string nameImportCarry6 = "carry6sub1.data";
                char const *fileparamImportCarry6 = nameImportCarry6.c_str();
                FILE *ansImportCarry6 = fopen(fileparamImportCarry6, "rb");
                for (int i = 0; i < 32; i++)
                    import_gate_bootstrapping_ciphertext_fromFile(ansImportCarry6, &anstextImportCarry6[i], params);
                fclose(ansImportCarry6);

                std::string nameImportCarry7 = "carry7sub1.data";
                char const *fileparamImportCarry7 = nameImportCarry7.c_str();
                FILE *ansImportCarry7 = fopen(fileparamImportCarry7, "rb");
                for (int i = 0; i < 32; i++)
                    import_gate_bootstrapping_ciphertext_fromFile(ansImportCarry7, &anstextImportCarry7[i], params);
                fclose(ansImportCarry7);

                // export the result ciphertexts to a file (for the cloud)
                for (int i = 0; i < 32; i++) // result1 behind
                    export_gate_bootstrapping_ciphertext_toFile(answer_data, &result[i], params);
                for (int i = 0; i < 32; i++) // result1 behind
                    export_gate_bootstrapping_ciphertext_toFile(answer_data, &carry1[i], params);
                for (int i = 0; i < 32; i++) // result1 infront
                    export_gate_bootstrapping_ciphertext_toFile(answer_data, &anstextImportCarry2[i], params);
                for (int i = 0; i < 32; i++) // result1 infront
                    export_gate_bootstrapping_ciphertext_toFile(answer_data, &anstextImportCarry3[i], params);
                for (int i = 0; i < 32; i++) // result1 infront
                    export_gate_bootstrapping_ciphertext_toFile(answer_data, &anstextImportCarry4[i], params);
                for (int i = 0; i < 32; i++) // result1 infron
                    export_gate_bootstrapping_ciphertext_toFile(answer_data, &anstextImportCarry5[i], params);
                for (int i = 0; i < 32; i++) // result1 infront
                    export_gate_bootstrapping_ciphertext_toFile(answer_data, &anstextImportCarry6[i], params);
                for (int i = 0; i < 32; i++) // result1 infront
                    export_gate_bootstrapping_ciphertext_toFile(answer_data, &anstextImportCarry7[i], params);
                for (int i = 0; i < 32; i++) // result2 binary2
                    export_gate_bootstrapping_ciphertext_toFile(answer_data, &anstext1[i], params);
                for (int i = 0; i < 32; i++) // result3 binary3
                    export_gate_bootstrapping_ciphertext_toFile(answer_data, &anstext2[i], params);
                for (int i = 0; i < 32; i++) // 5
                    export_gate_bootstrapping_ciphertext_toFile(answer_data, &anstext3[i], params);
                for (int i = 0; i < 32; i++) // 6
                    export_gate_bootstrapping_ciphertext_toFile(answer_data, &anstext4[i], params);
                for (int i = 0; i < 32; i++) // 7
                    export_gate_bootstrapping_ciphertext_toFile(answer_data, &anstext5[i], params);
                for (int i = 0; i < 32; i++) // 8
                    export_gate_bootstrapping_ciphertext_toFile(answer_data, &anstext6[i], params);
                for (int i = 0; i < 32; i++) // carry
                    export_gate_bootstrapping_ciphertext_toFile(answer_data, &ciphertextcarry1[i], params);
                fclose(answer_data);

                MPI_Barrier(MPI_COMM_WORLD);

                printf("writing the answer to file...\n");

                //Clean up

                printf("writing the answer to file...for %d, Calculating\n", world_rank);
                MPI_Barrier(MPI_COMM_WORLD);

                double time_taken;
                MPI_Reduce(&get_time, &time_taken, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

                FILE *t_file;
                if (world_rank == 0)
                {
                    printf("avg = %lf\n", time_taken / 32);
                    t_file = fopen(T_FILE, "a");
                    fprintf(t_file, "%lf\n", time_taken / 32);
                    fclose(t_file);
                }

                delete_gate_bootstrapping_ciphertext_array(32, result);
                delete_gate_bootstrapping_ciphertext_array(32, result2);
                delete_gate_bootstrapping_ciphertext_array(32, result3);
                delete_gate_bootstrapping_ciphertext_array(32, result4);
                delete_gate_bootstrapping_ciphertext_array(32, result5);
                delete_gate_bootstrapping_ciphertext_array(32, result6);
                delete_gate_bootstrapping_ciphertext_array(32, result7);
                delete_gate_bootstrapping_ciphertext_array(32, carry1);
                delete_gate_bootstrapping_ciphertext_array(32, carry2);
                delete_gate_bootstrapping_ciphertext_array(32, carry3);
                delete_gate_bootstrapping_ciphertext_array(32, carry4);
                delete_gate_bootstrapping_ciphertext_array(32, carry5);
                delete_gate_bootstrapping_ciphertext_array(32, carry6);
                delete_gate_bootstrapping_ciphertext_array(32, carry7);
                delete_gate_bootstrapping_ciphertext_array(32, ciphertextbit1);
                delete_gate_bootstrapping_ciphertext_array(32, ciphertextnegative1);
                delete_gate_bootstrapping_ciphertext_array(32, ciphertextbit2);
                delete_gate_bootstrapping_ciphertext_array(32, ciphertextnegative2);
                delete_gate_bootstrapping_ciphertext_array(32, ciphertext1);
                delete_gate_bootstrapping_ciphertext_array(32, ciphertext2);
                delete_gate_bootstrapping_ciphertext_array(32, ciphertext3);
                delete_gate_bootstrapping_ciphertext_array(32, ciphertext4);
                delete_gate_bootstrapping_ciphertext_array(32, ciphertext9);
                delete_gate_bootstrapping_ciphertext_array(32, ciphertext10);
                delete_gate_bootstrapping_ciphertext_array(32, ciphertext11);
                delete_gate_bootstrapping_ciphertext_array(32, ciphertext12);
                delete_gate_bootstrapping_ciphertext_array(32, temp1);
                delete_gate_bootstrapping_ciphertext_array(32, temp2);
                delete_gate_bootstrapping_ciphertext_array(32, temp3);
                delete_gate_bootstrapping_ciphertext_array(32, temp4);
                delete_gate_bootstrapping_ciphertext_array(32, anstext1);
                delete_gate_bootstrapping_ciphertext_array(32, anstext2);
                delete_gate_bootstrapping_ciphertext_array(32, anstext3);
                delete_gate_bootstrapping_ciphertext_array(32, anstext4);
                delete_gate_bootstrapping_ciphertext_array(32, anstext5);
                delete_gate_bootstrapping_ciphertext_array(32, anstext6);
                delete_gate_bootstrapping_ciphertext_array(32, anstextImportCarry2);
                delete_gate_bootstrapping_ciphertext_array(32, anstextImportCarry3);
                delete_gate_bootstrapping_ciphertext_array(32, anstextImportCarry4);
                delete_gate_bootstrapping_ciphertext_array(32, anstextImportCarry5);
                delete_gate_bootstrapping_ciphertext_array(32, anstextImportCarry6);
                delete_gate_bootstrapping_ciphertext_array(32, anstextImportCarry7);

                delete_gate_bootstrapping_ciphertext_array(32, tempStart);
                delete_gate_bootstrapping_ciphertext_array(32, inverse1);
                delete_gate_bootstrapping_ciphertext_array(32, inverse2);
                delete_gate_bootstrapping_ciphertext_array(32, inverse3);
                delete_gate_bootstrapping_ciphertext_array(32, inverse4);

                delete_gate_bootstrapping_ciphertext_array(32, tempcarry1);
                delete_gate_bootstrapping_ciphertext_array(32, tempcarry2);
                delete_gate_bootstrapping_ciphertext_array(32, tempcarry3);
                delete_gate_bootstrapping_ciphertext_array(32, tempcarry4);

                delete_gate_bootstrapping_ciphertext_array(32, twosresult1);
                delete_gate_bootstrapping_ciphertext_array(32, twosresult2);
                delete_gate_bootstrapping_ciphertext_array(32, twosresult3);
                delete_gate_bootstrapping_ciphertext_array(32, twosresult4);

                delete_gate_bootstrapping_ciphertext_array(32, twoscarry1);
                delete_gate_bootstrapping_ciphertext_array(32, twoscarry2);
                delete_gate_bootstrapping_ciphertext_array(32, twoscarry3);
                delete_gate_bootstrapping_ciphertext_array(32, twoscarry4);

                delete_gate_bootstrapping_cloud_keyset(bk);
                delete_gate_bootstrapping_secret_keyset(nbitkey);

                MPI_Finalize();
            }
        }
    }

    else if (x == 4)
    {
        MPI_Init(&argc, &argv);
        int world_rank;
        MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
        int world_size;
        MPI_Comm_size(MPI_COMM_WORLD, &world_size);

        printf("Hello World from process %d of %d\n", world_rank, world_size);
        MPI_Barrier(MPI_COMM_WORLD);

        printf("doing the homomorphic computation...on %d\n", world_rank);

        LweSample *result = new_gate_bootstrapping_ciphertext_array(32, params);
        LweSample *result1 = new_gate_bootstrapping_ciphertext_array(32, params);
        LweSample *anstext1 = new_gate_bootstrapping_ciphertext_array(32, params);
        LweSample *anstext2 = new_gate_bootstrapping_ciphertext_array(32, params);
        LweSample *carry1 = new_gate_bootstrapping_ciphertext_array(32, params);

        struct timeval start, end;
        double get_time;
        gettimeofday(&start, NULL);

        MPImul(result, ciphertext1, ciphertext9, ciphertextcarry1, 32, world_rank, bk);
        leftShift(result, result, bk, 32, world_rank);

        printf("writing the answer to file...for %d, Calculating\n", world_rank);

        std::string r = "r";
        std::string s = std::to_string(world_rank);
        s = r + s;
        std::string name = "file";
        char const *fileparam = s.c_str();
        FILE *answer_data1 = fopen(fileparam, "wb");

        for (int i = 0; i < 32; i++)
            export_gate_bootstrapping_ciphertext_toFile(answer_data1, &result[i], params);
        fclose(answer_data1);

        cleanup(result, 32, bk);

        MPI_Barrier(MPI_COMM_WORLD);

        if (world_rank % 2 == 0)
        {
            printf("Rank2\n");
            std::string r_source = "r";
            std::string s1 = std::to_string(world_rank);
            std::string s2 = std::to_string(world_rank + 1);
            s1 = r + s1;
            char const *fileparam1 = s1.c_str();

            FILE *ans1 = fopen(fileparam1, "rb");
            for (int i = 0; i < 32; i++)
                import_gate_bootstrapping_ciphertext_fromFile(ans1, &anstext1[i], params);
            fclose(ans1);

            s2 = r + s2;
            char const *fileparam2 = s2.c_str();

            FILE *ans2 = fopen(fileparam2, "rb");
            for (int i = 0; i < 32; i++)
                import_gate_bootstrapping_ciphertext_fromFile(ans2, &anstext2[i], params);
            fclose(ans2);

            add_addition(result1, carry1, anstext1, anstext2, ciphertextcarry1, 32, bk);

            std::string ans = "loop1answer.data";
            std::string loc = std::to_string(world_rank);
            ans = loc + ans;
            char const *writefileparam = ans.c_str();

            FILE *answer_data1 = fopen(writefileparam, "wb");
            for (int i = 0; i < 32; i++)
                export_gate_bootstrapping_ciphertext_toFile(answer_data1, &result1[i], params);
            fclose(answer_data1);

            cleanup(anstext1, 32, bk);
            cleanup(anstext2, 32, bk);
            cleanup(result1, 32, bk);
        }

        MPI_Barrier(MPI_COMM_WORLD);

        if (world_rank % 4 == 0)
        {
            printf("Rank4\n");
            std::string r_source = "loop1answer.data";
            std::string s1 = std::to_string(world_rank);
            std::string s2 = std::to_string(world_rank + 2);
            s1 = s1 + r_source;
            char const *fileparam1 = s1.c_str();

            FILE *ans1 = fopen(fileparam1, "rb");
            for (int i = 0; i < 32; i++)
            {
                import_gate_bootstrapping_ciphertext_fromFile(ans1, &anstext1[i], params);
            }
            fclose(ans1);

            s2 = s2 + r_source;
            char const *fileparam2 = s2.c_str();

            FILE *ans2 = fopen(fileparam2, "rb");
            for (int i = 0; i < 32; i++)
            {
                import_gate_bootstrapping_ciphertext_fromFile(ans2, &anstext2[i], params);
            }
            fclose(ans2);

            add_addition(result1, carry1, anstext1, anstext2, ciphertextcarry1, 32, bk);

            std::string ans = "loop2answer.data";
            std::string loc = std::to_string(world_rank);
            ans = loc + ans;
            char const *writefileparam = ans.c_str();

            FILE *answer_data1 = fopen(writefileparam, "wb");
            for (int i = 0; i < 32; i++)
                export_gate_bootstrapping_ciphertext_toFile(answer_data1, &result1[i], params);
            fclose(answer_data1);

            cleanup(anstext1, 32, bk);
            cleanup(anstext2, 32, bk);
            cleanup(result1, 32, bk);
        }

        MPI_Barrier(MPI_COMM_WORLD);

        if (world_rank % 8 == 0)
        {
            printf("Rank8\n");
            std::string r_source = "loop2answer.data";
            std::string s1 = std::to_string(world_rank);
            std::string s2 = std::to_string(world_rank + 4);
            s1 = s1 + r_source;
            char const *fileparam1 = s1.c_str();

            FILE *ans1 = fopen(fileparam1, "rb");
            for (int i = 0; i < 32; i++)
            {
                import_gate_bootstrapping_ciphertext_fromFile(ans1, &anstext1[i], params);
            }
            fclose(ans1);

            s2 = s2 + r_source;
            char const *fileparam2 = s2.c_str();

            FILE *ans2 = fopen(fileparam2, "rb");
            for (int i = 0; i < 32; i++)
            {
                import_gate_bootstrapping_ciphertext_fromFile(ans2, &anstext2[i], params);
            }
            fclose(ans2);

            add_addition(result1, carry1, anstext1, anstext2, ciphertextcarry1, 32, bk);

            std::string ans = "loop3answer.data";
            std::string loc = std::to_string(world_rank);
            ans = loc + ans;
            char const *writefileparam = ans.c_str();

            FILE *answer_data1 = fopen(writefileparam, "wb");
            for (int i = 0; i < 32; i++)
                export_gate_bootstrapping_ciphertext_toFile(answer_data1, &result1[i], params);
            fclose(answer_data1);

            cleanup(anstext1, 32, bk);
            cleanup(anstext2, 32, bk);
            cleanup(result1, 32, bk);
        }

        MPI_Barrier(MPI_COMM_WORLD);

        if (world_rank % 16 == 0)
        {
            printf("Rank16\n");
            std::string r_source = "loop3answer.data";
            std::string s1 = std::to_string(world_rank);
            std::string s2 = std::to_string(world_rank + 8);
            s1 = s1 + r_source;
            char const *fileparam1 = s1.c_str();

            FILE *ans1 = fopen(fileparam1, "rb");
            for (int i = 0; i < 32; i++)
                import_gate_bootstrapping_ciphertext_fromFile(ans1, &anstext1[i], params);
            fclose(ans1);

            s2 = s2 + r_source;
            char const *fileparam2 = s2.c_str();

            FILE *ans2 = fopen(fileparam2, "rb");
            for (int i = 0; i < 32; i++)
                import_gate_bootstrapping_ciphertext_fromFile(ans2, &anstext2[i], params);
            fclose(ans2);

            add_addition(result1, carry1, anstext1, anstext2, ciphertextcarry1, 32, bk);

            std::string ans = "loop4answer.data";
            std::string loc = std::to_string(world_rank);
            ans = loc + ans;
            char const *writefileparam = ans.c_str();

            FILE *answer_data1 = fopen(writefileparam, "wb");
            for (int i = 0; i < 32; i++)
                export_gate_bootstrapping_ciphertext_toFile(answer_data1, &result1[i], params);
            fclose(answer_data1);

            cleanup(anstext1, 32, bk);
            cleanup(anstext2, 32, bk);
            cleanup(result1, 32, bk);
        }

        if (world_rank % 32 == 0)
        {

            printf("Rank32\n");
            std::string r_source = "loop4answer.data";
            std::string s1 = std::to_string(world_rank);
            std::string s2 = std::to_string(world_rank + 16);
            s1 = s1 + r_source;
            char const *fileparam1 = s1.c_str();

            FILE *ans1 = fopen(fileparam1, "rb");
            for (int i = 0; i < 32; i++)
                import_gate_bootstrapping_ciphertext_fromFile(ans1, &anstext1[i], params);
            fclose(ans1);

            s2 = s2 + r_source;
            char const *fileparam2 = s2.c_str();

            FILE *ans2 = fopen(fileparam2, "rb");
            for (int i = 0; i < 32; i++)
                import_gate_bootstrapping_ciphertext_fromFile(ans2, &anstext2[i], params);
            fclose(ans2);

            add_addition(result1, carry1, anstext1, anstext2, ciphertextcarry1, 32, bk);

            printf("Exporting\n");
            //export the 32 ciphertexts to a file (for the cloud)
            for (int i = 0; i < 32; i++) // result1
                export_gate_bootstrapping_ciphertext_toFile(answer_data, &result1[i], params);
            for (int i = 0; i < 32; i++) // 2
                export_gate_bootstrapping_ciphertext_toFile(answer_data, &ciphertextcarry1[i], params);
            for (int i = 0; i < 32; i++) // 3
                export_gate_bootstrapping_ciphertext_toFile(answer_data, &ciphertextcarry1[i], params);
            for (int i = 0; i < 32; i++) // 4
                export_gate_bootstrapping_ciphertext_toFile(answer_data, &ciphertextcarry1[i], params);
            for (int i = 0; i < 32; i++) // 5
                export_gate_bootstrapping_ciphertext_toFile(answer_data, &ciphertextcarry1[i], params);
            for (int i = 0; i < 32; i++) // 6
                export_gate_bootstrapping_ciphertext_toFile(answer_data, &ciphertextcarry1[i], params);
            for (int i = 0; i < 32; i++) // 7
                export_gate_bootstrapping_ciphertext_toFile(answer_data, &ciphertextcarry1[i], params);
            for (int i = 0; i < 32; i++) // 8
                export_gate_bootstrapping_ciphertext_toFile(answer_data, &ciphertextcarry1[i], params);
            for (int i = 0; i < 32; i++) // carry
                export_gate_bootstrapping_ciphertext_toFile(answer_data, &ciphertextcarry1[i], params);
            fclose(answer_data);

            cleanup(anstext1, 32, bk);
            cleanup(anstext2, 32, bk);
        }

        MPI_Barrier(MPI_COMM_WORLD);

        gettimeofday(&end, NULL);
        get_time = (end.tv_sec - start.tv_sec) + (end.tv_usec - start.tv_usec) * 1.0E-6;
        printf("Computation Time: %lf[sec]\n", get_time);
        
        MPI_Barrier(MPI_COMM_WORLD);

        double time_taken;
        MPI_Reduce(&get_time, &time_taken, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

        FILE *t_file;
        if (world_rank == 0)
        {
            printf("avg = %lf\n", time_taken / 32);
            t_file = fopen(T_FILE, "a");
            fprintf(t_file, "%lf\n", time_taken / 32);
            fclose(t_file);
        }

        delete_gate_bootstrapping_ciphertext_array(32, result);
        delete_gate_bootstrapping_ciphertext_array(32, ciphertext9);
        delete_gate_bootstrapping_ciphertext_array(32, ciphertext1);
        delete_gate_bootstrapping_ciphertext_array(32, anstext1);
        delete_gate_bootstrapping_ciphertext_array(32, anstext2);
        delete_gate_bootstrapping_ciphertext_array(32, result1);

        delete_gate_bootstrapping_ciphertext_array(32, ciphertextbit1);
        delete_gate_bootstrapping_ciphertext_array(32, ciphertextnegative1);
        delete_gate_bootstrapping_ciphertext_array(32, ciphertextbit2);
        delete_gate_bootstrapping_ciphertext_array(32, ciphertextnegative2);
        delete_gate_bootstrapping_secret_keyset(nbitkey);
        delete_gate_bootstrapping_ciphertext_array(32, ciphertextcarry1);
        delete_gate_bootstrapping_cloud_keyset(bk);

        MPI_Finalize();
    }

    fclose(f);
    return 0;
}
