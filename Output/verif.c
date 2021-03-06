#include tfhetfhe.h
#include tfhetfhe_io.h
#include stdio.h
#include stdlib.h
#include bitset
#include iostream
#include iomanip
#include boostmultiprecisioncpp_int.hpp
#include boostdynamic_bitset.hpp
#include boostlexical_cast.hpp
#include string
#include systime.h
#include time.h
using namespace std;
using namespace boostmultiprecision;
using boostmultiprecisionint256_t;

#define SIZE 1024

int main()
{

    char mpioropenmp[SIZE];
    FILE f = fopen(mpioropenmp.txt, rb);
    fgets(mpioropenmp, SIZE, f);
    printf(This is workingn);
    printf(%sn, mpioropenmp);
    printf(%dn, strcmp(mpioropenmp, OpenMP));
    printf(%dn, strcmp(mpioropenmp, MPI));

    if (strcmp(mpioropenmp, OpenMP) == 0)
    {
        printf(OpenMP is being computed);
         Addition
         reads the secret key from file
        FILE secret_key = fopen(secret.key, rb);
        TFheGateBootstrappingSecretKeySet key = new_tfheGateBootstrappingSecretKeySet_fromFile(secret_key);
        fclose(secret_key);

         reads the nbit key from file
        FILE nbit_key = fopen(nbit.key, rb);
        TFheGateBootstrappingSecretKeySet nbitkey = new_tfheGateBootstrappingSecretKeySet_fromFile(nbit_key);
        fclose(nbit_key);

         if necessary, the params are inside the key
        const TFheGateBootstrappingParameterSet params = key-params;

         if necessary, the params are inside the key
        const TFheGateBootstrappingParameterSet nbitparams = nbitkey-params;

        struct timeval start, end;
        double get_time;
        gettimeofday(&start, NULL);

        read the 16 ciphertexts of the result
        LweSample negative = new_gate_bootstrapping_ciphertext_array(32, nbitparams);
        LweSample bit = new_gate_bootstrapping_ciphertext_array(32, nbitparams);

        FILE answer_data = fopen(answer.data, rb);
        for (int i = 0; i  32; i++)
            import_gate_bootstrapping_ciphertext_fromFile(answer_data, &negative[i], nbitparams);

        for (int i = 0; i  32; i++)
            import_gate_bootstrapping_ciphertext_fromFile(answer_data, &bit[i], nbitparams);

        decrypt and rebuild the answer

         negativity code
        int32_t int_negative = 0;
        for (int i = 0; i  32; i++)
        {
            int ai = bootsSymDecrypt(&negative[i], nbitkey)  0;
            int_negative = (ai  i);
        }
        stdcout  Negative   int_negative  n;

         TODO Obtain opcode from file
        int32_t int_op;
        FILE fptr;
        fptr = fopen(operator.txt, r);
        fscanf(fptr, %d, &int_op);
        stdcout  Opcode   int_op  n
                   n;

         Bit count
        int32_t int_bit = 0;
        for (int i = 0; i  32; i++)
        {
            int ai = bootsSymDecrypt(&bit[i], nbitkey)  0;
            int_bit = (ai  i);
        }

        if (int_op == 1)
        {
            stdcout  Result for   int_bit   bit Addition computation
                       n
                       n;

            if (int_bit == 32)
            {
                LweSample result = new_gate_bootstrapping_ciphertext_array(32, params);

                for (int i = 0; i  32; i++)
                    import_gate_bootstrapping_ciphertext_fromFile(answer_data, &result[i], params);
                fclose(answer_data);

                 decrypt and rebuild the answer
                int32_t int_answer1 = 0;
                for (int i = 0; i  32; i++)
                {
                    int ai = bootsSymDecrypt(&result[i], key)  0;
                    int_answer1 = (ai  i);
                }

                stdstring binary1 = stdbitset32(int_answer1).to_string();

                stdstring binary_combined = binary1;

                if (int_negative != 4)
                {
                    stdcout  The result in binary form is
                               n;
                    stdcout  binary_combined  n;
                }

                Length is the number of bits
                int length = binary_combined.length();
                int lastchar = length - 1;
                char char_array[length + 1];
                strcpy(char_array, binary_combined.c_str());

                int256_t total = 0;

                Positive or negative checking
                char charvalue = char_array[0];
                stdstring stringvalue(1, charvalue);
                int intvalue = stdstoi(stringvalue);

                if the result is either double positive or negative
                if (int_negative == 0  int_negative == 4)
                {
                    for (int i = 0; i  length; ++i)
                    {
                        char charvalue = char_array[i];
                        stdstring stringvalue(1, charvalue);
                        int intvalue = stdstoi(stringvalue);
                        int256_t calc1 = total  2;
                        total = calc1 + intvalue;
                    }
                    if either one of the values are negative, check for negative
                }
                else if (int_negative != 4)
                {
                    if (length == 32 && intvalue == 1)
                    {
                        int256_t lastcharvalue = 2;
                        for (int i = 0; i  lastchar - 1; ++i)
                        {
                            lastcharvalue = lastcharvalue  2;
                        }

                        for (int i = 1; i  length; ++i)
                        {
                            char charvalue = char_array[i];
                            stdstring stringvalue(1, charvalue);
                            int intvalue = stdstoi(stringvalue);
                            int256_t calc1 = total  2;
                            total = calc1 + intvalue;
                        }
                        total = total - lastcharvalue;
                    }
                    else
                    {
                        for (int i = 0; i  length; ++i)
                        {
                            char charvalue = char_array[i];
                            stdstring stringvalue(1, charvalue);
                            int intvalue = stdstoi(stringvalue);
                            int256_t calc1 = total  2;
                            total = calc1 + intvalue;
                        }
                    }
                }

                if its a double negative then inverse it
                if (int_negative == 4)
                {
                    int256_t invertedtotal = total  -1;
                    stdcout  The result in decimal form is
                               n;
                    stdcout  invertedtotal;
                }
                else
                {
                    stdcout  n
                               n;
                    stdcout  The result in decimal form is
                               n;
                    stdcout  total;
                }

                printf(n);
                printf(n);
                gettimeofday(&end, NULL);
                get_time = (end.tv_sec - start.tv_sec) + (end.tv_usec - start.tv_usec)  1.0E-6;
                printf(Computation Time %lf[sec]n, get_time);
                printf(n);

                printf(I hope you remembered what calculation you performed!n);

                 clean up all pointers
                delete_gate_bootstrapping_ciphertext_array(32, result);
                delete_gate_bootstrapping_secret_keyset(key);
                delete_gate_bootstrapping_secret_keyset(nbitkey);
            }

            else if (int_bit == 64)
            {
                LweSample result = new_gate_bootstrapping_ciphertext_array(32, params);
                LweSample result2 = new_gate_bootstrapping_ciphertext_array(32, params);

                 export the 64 ciphertexts to a file (for the cloud)

                for (int i = 0; i  32; i++)
                    import_gate_bootstrapping_ciphertext_fromFile(answer_data, &result[i], params);

                for (int i = 0; i  32; i++)
                    import_gate_bootstrapping_ciphertext_fromFile(answer_data, &result2[i], params);

                fclose(answer_data);

                 decrypt and rebuild the answer
                int32_t int_answer1 = 0;
                for (int i = 0; i  32; i++)
                {
                    int ai = bootsSymDecrypt(&result[i], key)  0;
                    int_answer1 = (ai  i);
                }

                int32_t int_answer2 = 0;
                for (int i = 0; i  32; i++)
                {
                    int ai = bootsSymDecrypt(&result2[i], key)  0;
                    int_answer2 = (ai  i);
                }

                stdstring binary1 = stdbitset32(int_answer1).to_string();
                stdstring binary2 = stdbitset32(int_answer2).to_string();
                stdstring binary_combined = binary2 + binary1;

                if (int_negative != 4)
                {
                    stdcout  The result in binary form is
                               n;
                    stdcout  binary_combined  n;
                }

                Length is the number of bits
                int length = binary_combined.length();
                int lastchar = length - 1;
                char char_array[length + 1];
                strcpy(char_array, binary_combined.c_str());

                int256_t total = 0;

                Positive or negative checking
                char charvalue = char_array[0];
                stdstring stringvalue(1, charvalue);
                int intvalue = stdstoi(stringvalue);
                if the result is either double positive or negative
                if (int_negative == 0  int_negative == 4)
                {
                    for (int i = 0; i  length; ++i)
                    {
                        char charvalue = char_array[i];
                        stdstring stringvalue(1, charvalue);
                        int intvalue = stdstoi(stringvalue);
                        int256_t calc1 = total  2;
                        total = calc1 + intvalue;
                    }
                    if either one of the values are negative, check for negative
                }
                else if (int_negative != 4)
                {
                    if (length == 64 && intvalue == 1)
                    {
                        int256_t lastcharvalue = 2;
                        for (int i = 0; i  lastchar - 1; ++i)
                        {
                            lastcharvalue = lastcharvalue  2;
                        }

                        for (int i = 1; i  length; ++i)
                        {
                            char charvalue = char_array[i];
                            stdstring stringvalue(1, charvalue);
                            int intvalue = stdstoi(stringvalue);
                            int256_t calc1 = total  2;
                            total = calc1 + intvalue;
                        }

                        total = total - lastcharvalue;
                    }
                    else
                    {
                        for (int i = 0; i  length; ++i)
                        {
                            char charvalue = char_array[i];
                            stdstring stringvalue(1, charvalue);
                            int intvalue = stdstoi(stringvalue);
                            int256_t calc1 = total  2;
                            total = calc1 + intvalue;
                        }
                    }
                }
                if (int_negative == 4)
                {
                    int256_t invertedtotal = total  -1;
                    stdcout  The result in decimal form is
                               n;
                    stdcout  invertedtotal;
                }
                else
                {
                    stdcout  n
                               n;
                    stdcout  The result in decimal form is
                               n;
                    stdcout  total;
                }

                printf(n);
                printf(n);
                gettimeofday(&end, NULL);
                get_time = (end.tv_sec - start.tv_sec) + (end.tv_usec - start.tv_usec)  1.0E-6;
                printf(Computation Time %lf[sec]n, get_time);
                printf(n);

                printf(I hope you remembered what calculation you performed!n);

                 clean up all pointers
                delete_gate_bootstrapping_ciphertext_array(32, result);
                delete_gate_bootstrapping_ciphertext_array(32, result2);
                delete_gate_bootstrapping_secret_keyset(key);
                delete_gate_bootstrapping_secret_keyset(nbitkey);
            }

            else if (int_bit == 128)
            {
                LweSample result = new_gate_bootstrapping_ciphertext_array(32, params);
                LweSample result2 = new_gate_bootstrapping_ciphertext_array(32, params);
                LweSample result3 = new_gate_bootstrapping_ciphertext_array(32, params);
                LweSample result4 = new_gate_bootstrapping_ciphertext_array(32, params);

                 export the 64 ciphertexts to a file (for the cloud)

                for (int i = 0; i  32; i++)
                    import_gate_bootstrapping_ciphertext_fromFile(answer_data, &result[i], params);

                for (int i = 0; i  32; i++)
                    import_gate_bootstrapping_ciphertext_fromFile(answer_data, &result2[i], params);

                for (int i = 0; i  32; i++)
                    import_gate_bootstrapping_ciphertext_fromFile(answer_data, &result3[i], params);

                for (int i = 0; i  32; i++)
                    import_gate_bootstrapping_ciphertext_fromFile(answer_data, &result4[i], params);

                fclose(answer_data);

                 decrypt and rebuild the answer
                int32_t int_answer1 = 0;
                for (int i = 0; i  32; i++)
                {
                    int ai = bootsSymDecrypt(&result[i], key)  0;
                    int_answer1 = (ai  i);
                }

                int32_t int_answer2 = 0;
                for (int i = 0; i  32; i++)
                {
                    int ai = bootsSymDecrypt(&result2[i], key)  0;
                    int_answer2 = (ai  i);
                }

                int32_t int_answer3 = 0;
                for (int i = 0; i  32; i++)
                {
                    int ai = bootsSymDecrypt(&result3[i], key)  0;
                    int_answer3 = (ai  i);
                }

                int32_t int_answer4 = 0;
                for (int i = 0; i  32; i++)
                {
                    int ai = bootsSymDecrypt(&result4[i], key)  0;
                    int_answer4 = (ai  i);
                }

                stdstring binary1 = stdbitset32(int_answer1).to_string();
                stdstring binary2 = stdbitset32(int_answer2).to_string();
                stdstring binary3 = stdbitset32(int_answer3).to_string();
                stdstring binary4 = stdbitset32(int_answer4).to_string();
                stdstring binary_combined = binary4 + binary3 + binary2 + binary1;

                if (int_negative != 4)
                {
                    stdcout  The result in binary form is
                               n;
                    stdcout  binary_combined  n;
                }
                Length is the number of bits
                int length = binary_combined.length();
                int lastchar = length - 1;
                char char_array[length + 1];
                strcpy(char_array, binary_combined.c_str());

                int256_t total = 0;

                Positive or negative checking
                char charvalue = char_array[0];
                stdstring stringvalue(1, charvalue);
                int intvalue = stdstoi(stringvalue);

                if the result is either double positive or negative
                if (int_negative == 0  int_negative == 4)
                {
                    for (int i = 0; i  length; ++i)
                    {
                        char charvalue = char_array[i];
                        stdstring stringvalue(1, charvalue);
                        int intvalue = stdstoi(stringvalue);
                        int256_t calc1 = total  2;
                        total = calc1 + intvalue;
                    }
                }
                if either one of the values are negative, check for negative
                else if (int_negative != 4)
                {
                    if (length == 128 && intvalue == 1)
                    {
                        int256_t lastcharvalue = 2;

                        for (int i = 0; i  lastchar - 1; ++i)
                        {
                            lastcharvalue = lastcharvalue  2;
                        }

                        for (int i = 1; i  length; ++i)
                        {
                            char charvalue = char_array[i];
                            stdstring stringvalue(1, charvalue);
                            int intvalue = stdstoi(stringvalue);
                            int256_t calc1 = total  2;
                            total = calc1 + intvalue;
                        }
                        total = total - lastcharvalue;
                    }
                    else
                    {
                        for (int i = 0; i  length; ++i)
                        {
                            char charvalue = char_array[i];
                            stdstring stringvalue(1, charvalue);
                            int intvalue = stdstoi(stringvalue);
                            int256_t calc1 = total  2;
                            total = calc1 + intvalue;
                        }
                    }
                }
                if (int_negative == 4)
                {
                    int256_t invertedtotal = total  -1;
                    stdcout  The result in decimal form is
                               n;
                    stdcout  invertedtotal;
                }
                else
                {
                    stdcout  n
                               n;
                    stdcout  The result in decimal form is
                               n;
                    stdcout  total;
                }
                printf(n);
                printf(n);
                gettimeofday(&end, NULL);
                get_time = (end.tv_sec - start.tv_sec) + (end.tv_usec - start.tv_usec)  1.0E-6;
                printf(Computation Time %lf[sec]n, get_time);
                printf(n);

                printf(I hope you remembered what calculation you performed!n);

                 clean up all pointers
                delete_gate_bootstrapping_ciphertext_array(32, result);
                delete_gate_bootstrapping_ciphertext_array(32, result2);
                delete_gate_bootstrapping_ciphertext_array(32, result3);
                delete_gate_bootstrapping_ciphertext_array(32, result4);
                delete_gate_bootstrapping_secret_keyset(key);
                delete_gate_bootstrapping_secret_keyset(nbitkey);
            }

            else if (int_bit == 256)
            {
                LweSample result = new_gate_bootstrapping_ciphertext_array(32, params);
                LweSample result2 = new_gate_bootstrapping_ciphertext_array(32, params);
                LweSample result3 = new_gate_bootstrapping_ciphertext_array(32, params);
                LweSample result4 = new_gate_bootstrapping_ciphertext_array(32, params);
                LweSample result5 = new_gate_bootstrapping_ciphertext_array(32, params);
                LweSample result6 = new_gate_bootstrapping_ciphertext_array(32, params);
                LweSample result7 = new_gate_bootstrapping_ciphertext_array(32, params);
                LweSample result8 = new_gate_bootstrapping_ciphertext_array(32, params);

                 export the 64 ciphertexts to a file (for the cloud)

                for (int i = 0; i  32; i++)
                    import_gate_bootstrapping_ciphertext_fromFile(answer_data, &result[i], params);

                for (int i = 0; i  32; i++)
                    import_gate_bootstrapping_ciphertext_fromFile(answer_data, &result2[i], params);

                for (int i = 0; i  32; i++)
                    import_gate_bootstrapping_ciphertext_fromFile(answer_data, &result3[i], params);

                for (int i = 0; i  32; i++)
                    import_gate_bootstrapping_ciphertext_fromFile(answer_data, &result4[i], params);

                for (int i = 0; i  32; i++)
                    import_gate_bootstrapping_ciphertext_fromFile(answer_data, &result5[i], params);

                for (int i = 0; i  32; i++)
                    import_gate_bootstrapping_ciphertext_fromFile(answer_data, &result6[i], params);

                for (int i = 0; i  32; i++)
                    import_gate_bootstrapping_ciphertext_fromFile(answer_data, &result7[i], params);

                for (int i = 0; i  32; i++)
                    import_gate_bootstrapping_ciphertext_fromFile(answer_data, &result8[i], params);
                fclose(answer_data);

                 decrypt and rebuild the answer
                int32_t int_answer1 = 0;
                for (int i = 0; i  32; i++)
                {
                    int ai = bootsSymDecrypt(&result[i], key)  0;
                    int_answer1 = (ai  i);
                }

                int32_t int_answer2 = 0;
                for (int i = 0; i  32; i++)
                {
                    int ai = bootsSymDecrypt(&result2[i], key)  0;
                    int_answer2 = (ai  i);
                }

                int32_t int_answer3 = 0;
                for (int i = 0; i  32; i++)
                {
                    int ai = bootsSymDecrypt(&result3[i], key)  0;
                    int_answer3 = (ai  i);
                }

                int32_t int_answer4 = 0;
                for (int i = 0; i  32; i++)
                {
                    int ai = bootsSymDecrypt(&result4[i], key)  0;
                    int_answer4 = (ai  i);
                }

                int32_t int_answer5 = 0;
                for (int i = 0; i  32; i++)
                {
                    int ai = bootsSymDecrypt(&result5[i], key)  0;
                    int_answer5 = (ai  i);
                }

                int32_t int_answer6 = 0;
                for (int i = 0; i  32; i++)
                {
                    int ai = bootsSymDecrypt(&result6[i], key)  0;
                    int_answer6 = (ai  i);
                }

                int32_t int_answer7 = 0;
                for (int i = 0; i  32; i++)
                {
                    int ai = bootsSymDecrypt(&result7[i], key)  0;
                    int_answer7 = (ai  i);
                }

                int32_t int_answer8 = 0;
                for (int i = 0; i  32; i++)
                {
                    int ai = bootsSymDecrypt(&result8[i], key)  0;
                    int_answer8 = (ai  i);
                }

                stdstring binary1 = stdbitset32(int_answer1).to_string();
                stdstring binary2 = stdbitset32(int_answer2).to_string();
                stdstring binary3 = stdbitset32(int_answer3).to_string();
                stdstring binary4 = stdbitset32(int_answer4).to_string();
                stdstring binary5 = stdbitset32(int_answer5).to_string();
                stdstring binary6 = stdbitset32(int_answer6).to_string();
                stdstring binary7 = stdbitset32(int_answer7).to_string();
                stdstring binary8 = stdbitset32(int_answer8).to_string();
                stdstring binary_combined = binary8 + binary7 + binary6 + binary5 + binary4 + binary3 + binary2 + binary1;

                if (int_negative != 4)
                {
                    stdcout  The result in binary form is
                               n;
                    stdcout  binary_combined  n;
                }

                Length is the number of bits
                int length = binary_combined.length();
                int lastchar = length - 1;
                char char_array[length + 1];
                strcpy(char_array, binary_combined.c_str());

                int256_t total = 0;

                Positive or negative checking
                char charvalue = char_array[0];
                stdstring stringvalue(1, charvalue);
                int intvalue = stdstoi(stringvalue);

                if the result is either double positive or negative
                if (int_negative == 0  int_negative == 4)
                {
                    for (int i = 0; i  length; ++i)
                    {
                        char charvalue = char_array[i];
                        stdstring stringvalue(1, charvalue);
                        int intvalue = stdstoi(stringvalue);
                        int256_t calc1 = total  2;
                        total = calc1 + intvalue;
                    }
                }
                if either one of the values are negative, check for negative
                else if (int_negative == 1  int_negative == 2)
                {
                    if (length == 256 && intvalue == 1)
                    {
                        int256_t lastcharvalue = 2;

                        for (int i = 0; i  lastchar - 1; ++i)
                        {
                            lastcharvalue = lastcharvalue  2;
                        }

                        for (int i = 1; i  length; ++i)
                        {
                            char charvalue = char_array[i];
                            stdstring stringvalue(1, charvalue);
                            int intvalue = stdstoi(stringvalue);
                            int256_t calc1 = total  2;
                            total = calc1 + intvalue;
                        }
                        total = total - lastcharvalue;
                    }
                    else
                    {
                        for (int i = 0; i  length; ++i)
                        {
                            char charvalue = char_array[i];
                            stdstring stringvalue(1, charvalue);
                            int intvalue = stdstoi(stringvalue);
                            int256_t calc1 = total  2;
                            total = calc1 + intvalue;
                        }
                    }
                }

                if (int_negative == 4)
                {
                    int256_t invertedtotal = total  -1;
                    stdcout  The result in decimal form is
                               n;
                    stdcout  invertedtotal;
                }
                else
                {
                    stdcout  n
                               n;
                    stdcout  The result in decimal form is
                               n;
                    stdcout  total;
                }

                printf(n);
                printf(n);
                gettimeofday(&end, NULL);
                get_time = (end.tv_sec - start.tv_sec) + (end.tv_usec - start.tv_usec)  1.0E-6;
                printf(Computation Time %lf[sec]n, get_time);
                printf(n);

                printf(I hope you remembered what calculation you performed!n);

                 clean up all pointers
                delete_gate_bootstrapping_ciphertext_array(32, result);
                delete_gate_bootstrapping_ciphertext_array(32, result2);
                delete_gate_bootstrapping_ciphertext_array(32, result3);
                delete_gate_bootstrapping_ciphertext_array(32, result4);
                delete_gate_bootstrapping_ciphertext_array(32, result5);
                delete_gate_bootstrapping_ciphertext_array(32, result6);
                delete_gate_bootstrapping_ciphertext_array(32, result7);
                delete_gate_bootstrapping_ciphertext_array(32, result8);
                delete_gate_bootstrapping_secret_keyset(key);
                delete_gate_bootstrapping_secret_keyset(nbitkey);
            }
        }

         Subtraction
        else if (int_op == 2)
        {
            stdcout  Result for   int_bit   bit Subtraction computation
                       n
                       n;

            if (int_bit == 32)
            {
                LweSample result = new_gate_bootstrapping_ciphertext_array(32, params);
                for (int i = 0; i  32; i++)
                    import_gate_bootstrapping_ciphertext_fromFile(answer_data, &result[i], params);
                fclose(answer_data);

                 decrypt and rebuild the answer
                int32_t int_answer1 = 0;
                for (int i = 0; i  32; i++)
                {
                    int ai = bootsSymDecrypt(&result[i], key)  0;
                    int_answer1 = (ai  i);
                }

                stdstring binary1 = stdbitset32(int_answer1).to_string();
                stdstring binary_combined = binary1;

                if (int_negative != 1)
                {
                    stdcout  The result in binary form is
                               n;
                    stdcout  binary_combined;
                }

                Length is the number of bits
                int length = binary_combined.length();
                int lastchar = length - 1;
                char char_array[length + 1];
                strcpy(char_array, binary_combined.c_str());
                int256_t total = 0;

                Positive or negative checking
                char charvalue = char_array[0];
                stdstring stringvalue(1, charvalue);
                int intvalue = stdstoi(stringvalue);

                if (int_negative == 2)
                {

                    for (int i = 0; i  length; ++i)
                    {
                        char charvalue = char_array[i];
                        stdstring stringvalue(1, charvalue);
                        int intvalue = stdstoi(stringvalue);
                        int256_t calc1 = total  2;
                        total = calc1 + intvalue;
                    }
                }
                else
                {
                     int_negative is not equals to 2
                    if (length == 32 && intvalue == 1)
                    {
                        int256_t lastcharvalue = 2;

                        for (int i = 0; i  lastchar - 1; ++i)
                        {
                            lastcharvalue = lastcharvalue  2;
                        }

                        for (int i = 1; i  length; ++i)
                        {
                            char charvalue = char_array[i];
                            stdstring stringvalue(1, charvalue);
                            int intvalue = stdstoi(stringvalue);
                            int256_t calc1 = total  2;
                            total = calc1 + intvalue;
                        }
                        total = total - lastcharvalue;
                    }

                    else
                    {
                        for (int i = 0; i  length; ++i)
                        {
                            char charvalue = char_array[i];
                            stdstring stringvalue(1, charvalue);
                            int intvalue = stdstoi(stringvalue);
                            int256_t calc1 = total  2;
                            total = calc1 + intvalue;
                        }
                    }
                }

                if (int_negative == 1)
                {
                    int256_t invertedtotal = total  -1;
                    stdcout  The result in decimal form is
                               n;
                    stdcout  invertedtotal;
                }
                else
                {
                    stdcout  n
                               n;
                    stdcout  The result in decimal form is
                               n;
                    stdcout  total;
                }
                printf(n);
                printf(n);
                gettimeofday(&end, NULL);
                get_time = (end.tv_sec - start.tv_sec) + (end.tv_usec - start.tv_usec)  1.0E-6;
                printf(Computation Time %lf[sec]n, get_time);
                printf(n);

                printf(I hope you remembered what calculation you performed!n);

                 clean up all pointers
                delete_gate_bootstrapping_ciphertext_array(32, result);
                delete_gate_bootstrapping_secret_keyset(key);
                delete_gate_bootstrapping_secret_keyset(nbitkey);
            }

            else if (int_bit == 64)
            {
                LweSample result = new_gate_bootstrapping_ciphertext_array(32, params);
                LweSample result2 = new_gate_bootstrapping_ciphertext_array(32, params);

                 export the 64 ciphertexts to a file (for the cloud)

                for (int i = 0; i  32; i++)
                    import_gate_bootstrapping_ciphertext_fromFile(answer_data, &result[i], params);

                for (int i = 0; i  32; i++)
                    import_gate_bootstrapping_ciphertext_fromFile(answer_data, &result2[i], params);

                fclose(answer_data);

                 decrypt and rebuild the answer
                int32_t int_answer1 = 0;
                for (int i = 0; i  32; i++)
                {
                    int ai = bootsSymDecrypt(&result[i], key)  0;
                    int_answer1 = (ai  i);
                }

                int32_t int_answer2 = 0;
                for (int i = 0; i  32; i++)
                {
                    int ai = bootsSymDecrypt(&result2[i], key)  0;
                    int_answer2 = (ai  i);
                }

                stdstring binary1 = stdbitset32(int_answer1).to_string();
                stdstring binary2 = stdbitset32(int_answer2).to_string();

                stdstring binary_combined = binary2 + binary1;
                if (int_negative != 1)
                {
                    stdcout  The result in binary form is
                               n;
                    stdcout  binary_combined;
                }

                Length is the number of bits
                int length = binary_combined.length();
                int lastchar = length - 1;
                char char_array[length + 1];
                strcpy(char_array, binary_combined.c_str());

                int256_t total = 0;

                Positive or negative checking
                char charvalue = char_array[0];
                stdstring stringvalue(1, charvalue);
                int intvalue = stdstoi(stringvalue);

                if (int_negative == 2)
                {
                    for (int i = 0; i  length; ++i)
                    {
                        char charvalue = char_array[i];
                        stdstring stringvalue(1, charvalue);
                        int intvalue = stdstoi(stringvalue);
                        int256_t calc1 = total  2;
                        total = calc1 + intvalue;
                    }
                }
                else
                {
                    if (length == 64 && intvalue == 1)
                    {
                        int256_t lastcharvalue = 2;

                        for (int i = 0; i  lastchar - 1; ++i)
                        {
                            lastcharvalue = lastcharvalue  2;
                        }

                        for (int i = 1; i  length; ++i)
                        {
                            char charvalue = char_array[i];
                            stdstring stringvalue(1, charvalue);
                            int intvalue = stdstoi(stringvalue);
                            int256_t calc1 = total  2;
                            total = calc1 + intvalue;
                        }
                        total = total - lastcharvalue;
                    }
                    else
                    {
                        for (int i = 0; i  length; ++i)
                        {
                            char charvalue = char_array[i];
                            stdstring stringvalue(1, charvalue);
                            int intvalue = stdstoi(stringvalue);
                            int256_t calc1 = total  2;
                            total = calc1 + intvalue;
                        }
                    }
                }
                if (int_negative == 1)
                {
                    int256_t invertedtotal = total  -1;
                    stdcout  The result in decimal form is
                               n;
                    stdcout  invertedtotal;
                }
                else
                {
                    stdcout  n
                               n;
                    stdcout  The result in decimal form is
                               n;
                    stdcout  total;
                }

                printf(n);
                printf(n);
                gettimeofday(&end, NULL);
                get_time = (end.tv_sec - start.tv_sec) + (end.tv_usec - start.tv_usec)  1.0E-6;
                printf(Computation Time %lf[sec]n, get_time);
                printf(n);

                printf(I hope you remembered what calculation you performed!n);

                 clean up all pointers
                delete_gate_bootstrapping_ciphertext_array(32, result);
                delete_gate_bootstrapping_ciphertext_array(32, result2);
                delete_gate_bootstrapping_secret_keyset(key);
                delete_gate_bootstrapping_secret_keyset(nbitkey);
            }

            else if (int_bit == 128)
            {
                LweSample result = new_gate_bootstrapping_ciphertext_array(32, params);
                LweSample result2 = new_gate_bootstrapping_ciphertext_array(32, params);
                LweSample result3 = new_gate_bootstrapping_ciphertext_array(32, params);
                LweSample result4 = new_gate_bootstrapping_ciphertext_array(32, params);

                 export the 64 ciphertexts to a file (for the cloud)

                for (int i = 0; i  32; i++)
                    import_gate_bootstrapping_ciphertext_fromFile(answer_data, &result[i], params);

                for (int i = 0; i  32; i++)
                    import_gate_bootstrapping_ciphertext_fromFile(answer_data, &result2[i], params);

                for (int i = 0; i  32; i++)
                    import_gate_bootstrapping_ciphertext_fromFile(answer_data, &result3[i], params);

                for (int i = 0; i  32; i++)
                    import_gate_bootstrapping_ciphertext_fromFile(answer_data, &result4[i], params);

                fclose(answer_data);

                 decrypt and rebuild the answer
                int32_t int_answer1 = 0;
                for (int i = 0; i  32; i++)
                {
                    int ai = bootsSymDecrypt(&result[i], key)  0;
                    int_answer1 = (ai  i);
                }

                int32_t int_answer2 = 0;
                for (int i = 0; i  32; i++)
                {
                    int ai = bootsSymDecrypt(&result2[i], key)  0;
                    int_answer2 = (ai  i);
                }

                int32_t int_answer3 = 0;
                for (int i = 0; i  32; i++)
                {
                    int ai = bootsSymDecrypt(&result3[i], key)  0;
                    int_answer3 = (ai  i);
                }

                int32_t int_answer4 = 0;
                for (int i = 0; i  32; i++)
                {
                    int ai = bootsSymDecrypt(&result4[i], key)  0;
                    int_answer4 = (ai  i);
                }

                stdstring binary1 = stdbitset32(int_answer1).to_string();
                stdstring binary2 = stdbitset32(int_answer2).to_string();
                stdstring binary3 = stdbitset32(int_answer3).to_string();
                stdstring binary4 = stdbitset32(int_answer4).to_string();
                stdstring binary_combined = binary4 + binary3 + binary2 + binary1;

                if (int_negative != 1)
                {
                    stdcout  The result in binary form is
                               n;
                    stdcout  binary_combined;
                }

                Length is the number of bits
                int length = binary_combined.length();
                int lastchar = length - 1;
                char char_array[length + 1];
                strcpy(char_array, binary_combined.c_str());

                int256_t total = 0;

                Positive or negative checking
                char charvalue = char_array[0];
                stdstring stringvalue(1, charvalue);
                int intvalue = stdstoi(stringvalue);

                if (int_negative == 2)
                {
                    for (int i = 0; i  length; ++i)
                    {
                        char charvalue = char_array[i];
                        stdstring stringvalue(1, charvalue);
                        int intvalue = stdstoi(stringvalue);
                        int256_t calc1 = total  2;
                        total = calc1 + intvalue;
                    }
                }
                else
                {
                    if (length == 128 && intvalue == 1)
                    {
                        int256_t lastcharvalue = 2;

                        for (int i = 0; i  lastchar - 1; ++i)
                        {
                            lastcharvalue = lastcharvalue  2;
                        }

                        for (int i = 1; i  length; ++i)
                        {
                            char charvalue = char_array[i];
                            stdstring stringvalue(1, charvalue);
                            int intvalue = stdstoi(stringvalue);
                            int256_t calc1 = total  2;
                            total = calc1 + intvalue;
                        }

                        total = total - lastcharvalue;
                    }
                    else
                    {
                        for (int i = 0; i  length; ++i)
                        {
                            char charvalue = char_array[i];
                            stdstring stringvalue(1, charvalue);
                            int intvalue = stdstoi(stringvalue);
                            int256_t calc1 = total  2;
                            total = calc1 + intvalue;
                        }
                    }
                }
                if (int_negative == 1)
                {
                    int256_t invertedtotal = total  -1;
                    stdcout  The result in decimal form is
                               n;
                    stdcout  invertedtotal;
                }
                else
                {
                    stdcout  n
                               n;
                    stdcout  The result in decimal form is
                               n;
                    stdcout  total;
                }

                printf(n);
                printf(n);
                gettimeofday(&end, NULL);
                get_time = (end.tv_sec - start.tv_sec) + (end.tv_usec - start.tv_usec)  1.0E-6;
                printf(Computation Time %lf[sec]n, get_time);
                printf(n);

                printf(I hope you remembered what calculation you performed!n);

                 clean up all pointers
                delete_gate_bootstrapping_ciphertext_array(32, result);
                delete_gate_bootstrapping_ciphertext_array(32, result2);
                delete_gate_bootstrapping_ciphertext_array(32, result3);
                delete_gate_bootstrapping_ciphertext_array(32, result4);
                delete_gate_bootstrapping_secret_keyset(key);
                delete_gate_bootstrapping_secret_keyset(nbitkey);
            }

            else if (int_bit == 256)
            {
                LweSample result = new_gate_bootstrapping_ciphertext_array(32, params);
                LweSample result2 = new_gate_bootstrapping_ciphertext_array(32, params);
                LweSample result3 = new_gate_bootstrapping_ciphertext_array(32, params);
                LweSample result4 = new_gate_bootstrapping_ciphertext_array(32, params);
                LweSample result5 = new_gate_bootstrapping_ciphertext_array(32, params);
                LweSample result6 = new_gate_bootstrapping_ciphertext_array(32, params);
                LweSample result7 = new_gate_bootstrapping_ciphertext_array(32, params);
                LweSample result8 = new_gate_bootstrapping_ciphertext_array(32, params);

                 export the 64 ciphertexts to a file (for the cloud)

                for (int i = 0; i  32; i++)
                    import_gate_bootstrapping_ciphertext_fromFile(answer_data, &result[i], params);

                for (int i = 0; i  32; i++)
                    import_gate_bootstrapping_ciphertext_fromFile(answer_data, &result2[i], params);

                for (int i = 0; i  32; i++)
                    import_gate_bootstrapping_ciphertext_fromFile(answer_data, &result3[i], params);

                for (int i = 0; i  32; i++)
                    import_gate_bootstrapping_ciphertext_fromFile(answer_data, &result4[i], params);

                for (int i = 0; i  32; i++)
                    import_gate_bootstrapping_ciphertext_fromFile(answer_data, &result5[i], params);

                for (int i = 0; i  32; i++)
                    import_gate_bootstrapping_ciphertext_fromFile(answer_data, &result6[i], params);

                for (int i = 0; i  32; i++)
                    import_gate_bootstrapping_ciphertext_fromFile(answer_data, &result7[i], params);

                for (int i = 0; i  32; i++)
                    import_gate_bootstrapping_ciphertext_fromFile(answer_data, &result8[i], params);
                fclose(answer_data);

                 decrypt and rebuild the answer
                int32_t int_answer1 = 0;
                for (int i = 0; i  32; i++)
                {
                    int ai = bootsSymDecrypt(&result[i], key)  0;
                    int_answer1 = (ai  i);
                }

                int32_t int_answer2 = 0;
                for (int i = 0; i  32; i++)
                {
                    int ai = bootsSymDecrypt(&result2[i], key)  0;
                    int_answer2 = (ai  i);
                }

                int32_t int_answer3 = 0;
                for (int i = 0; i  32; i++)
                {
                    int ai = bootsSymDecrypt(&result3[i], key)  0;
                    int_answer3 = (ai  i);
                }

                int32_t int_answer4 = 0;
                for (int i = 0; i  32; i++)
                {
                    int ai = bootsSymDecrypt(&result4[i], key)  0;
                    int_answer4 = (ai  i);
                }

                int32_t int_answer5 = 0;
                for (int i = 0; i  32; i++)
                {
                    int ai = bootsSymDecrypt(&result5[i], key)  0;
                    int_answer5 = (ai  i);
                }

                int32_t int_answer6 = 0;
                for (int i = 0; i  32; i++)
                {
                    int ai = bootsSymDecrypt(&result6[i], key)  0;
                    int_answer6 = (ai  i);
                }

                int32_t int_answer7 = 0;
                for (int i = 0; i  32; i++)
                {
                    int ai = bootsSymDecrypt(&result7[i], key)  0;
                    int_answer7 = (ai  i);
                }

                int32_t int_answer8 = 0;
                for (int i = 0; i  32; i++)
                {
                    int ai = bootsSymDecrypt(&result8[i], key)  0;
                    int_answer8 = (ai  i);
                }

                stdstring binary1 = stdbitset32(int_answer1).to_string();
                stdstring binary2 = stdbitset32(int_answer2).to_string();
                stdstring binary3 = stdbitset32(int_answer3).to_string();
                stdstring binary4 = stdbitset32(int_answer4).to_string();
                stdstring binary5 = stdbitset32(int_answer5).to_string();
                stdstring binary6 = stdbitset32(int_answer6).to_string();
                stdstring binary7 = stdbitset32(int_answer7).to_string();
                stdstring binary8 = stdbitset32(int_answer8).to_string();
                stdstring binary_combined = binary8 + binary7 + binary6 + binary5 + binary4 + binary3 + binary2 + binary1;

                if (int_negative != 1)
                {
                    stdcout  The result in binary form is
                               n;
                    stdcout  binary_combined;
                }

                Length is the number of bits
                int length = binary_combined.length();
                int lastchar = length - 1;
                char char_array[length + 1];
                strcpy(char_array, binary_combined.c_str());

                int256_t total = 0;

                Positive or negative checking
                char charvalue = char_array[0];
                stdstring stringvalue(1, charvalue);
                int intvalue = stdstoi(stringvalue);

                if (int_negative == 2)
                {
                    for (int i = 0; i  length; ++i)
                    {
                        char charvalue = char_array[i];
                        stdstring stringvalue(1, charvalue);
                        int intvalue = stdstoi(stringvalue);
                        int256_t calc1 = total  2;
                        total = calc1 + intvalue;
                    }
                }
                else
                {
                    if (length == 256 && intvalue == 1)
                    {
                        int256_t lastcharvalue = 2;

                        for (int i = 0; i  lastchar - 1; ++i)
                        {
                            lastcharvalue = lastcharvalue  2;
                        }

                        for (int i = 1; i  length; ++i)
                        {
                            char charvalue = char_array[i];
                            stdstring stringvalue(1, charvalue);
                            int intvalue = stdstoi(stringvalue);
                            int256_t calc1 = total  2;
                            total = calc1 + intvalue;
                        }
                        total = total - lastcharvalue;
                    }
                    else
                    {
                        for (int i = 0; i  length; ++i)
                        {
                            char charvalue = char_array[i];
                            stdstring stringvalue(1, charvalue);
                            int intvalue = stdstoi(stringvalue);
                            int256_t calc1 = total  2;
                            total = calc1 + intvalue;
                        }
                    }
                }

                if (int_negative == 1)
                {
                    int256_t invertedtotal = total  -1;
                    stdcout  The result in decimal form is
                               n;
                    stdcout  invertedtotal;
                }
                else
                {
                    stdcout  n
                               n;
                    stdcout  The result in decimal form is
                               n;
                    stdcout  total;
                }
                printf(n);
                printf(n);
                gettimeofday(&end, NULL);
                get_time = (end.tv_sec - start.tv_sec) + (end.tv_usec - start.tv_usec)  1.0E-6;
                printf(Computation Time %lf[sec]n, get_time);
                printf(n);

                printf(I hope you remembered what calculation you performed!n);

                 clean up all pointers
                delete_gate_bootstrapping_ciphertext_array(32, result);
                delete_gate_bootstrapping_ciphertext_array(32, result2);
                delete_gate_bootstrapping_ciphertext_array(32, result3);
                delete_gate_bootstrapping_ciphertext_array(32, result4);
                delete_gate_bootstrapping_ciphertext_array(32, result5);
                delete_gate_bootstrapping_ciphertext_array(32, result6);
                delete_gate_bootstrapping_ciphertext_array(32, result7);
                delete_gate_bootstrapping_ciphertext_array(32, result8);
                delete_gate_bootstrapping_secret_keyset(key);
                delete_gate_bootstrapping_secret_keyset(nbitkey);
            }
        }

         Multiplication
        else if (int_op == 4)
        {
            stdcout  Result for   int_bit   bit Multiplication computation
                       n
                       n;
            if (int_bit == 256)
            {
                LweSample finalresult1 = new_gate_bootstrapping_ciphertext_array(32, params);
                LweSample finalresult2 = new_gate_bootstrapping_ciphertext_array(32, params);
                LweSample finalresult3 = new_gate_bootstrapping_ciphertext_array(32, params);
                LweSample finalresult4 = new_gate_bootstrapping_ciphertext_array(32, params);
                LweSample finalresult5 = new_gate_bootstrapping_ciphertext_array(32, params);
                LweSample finalresult6 = new_gate_bootstrapping_ciphertext_array(32, params);
                LweSample finalresult7 = new_gate_bootstrapping_ciphertext_array(32, params);
                LweSample finalresult8 = new_gate_bootstrapping_ciphertext_array(32, params);

                export the 32 ciphertexts to a file (for the cloud)

                for (int i = 0; i  32; i++)
                    import_gate_bootstrapping_ciphertext_fromFile(answer_data, &finalresult1[i], params);

                for (int i = 0; i  32; i++)
                    import_gate_bootstrapping_ciphertext_fromFile(answer_data, &finalresult2[i], params);

                for (int i = 0; i  32; i++)
                    import_gate_bootstrapping_ciphertext_fromFile(answer_data, &finalresult3[i], params);

                for (int i = 0; i  32; i++)
                    import_gate_bootstrapping_ciphertext_fromFile(answer_data, &finalresult4[i], params);

                for (int i = 0; i  32; i++)
                    import_gate_bootstrapping_ciphertext_fromFile(answer_data, &finalresult5[i], params);

                for (int i = 0; i  32; i++)
                    import_gate_bootstrapping_ciphertext_fromFile(answer_data, &finalresult6[i], params);

                for (int i = 0; i  32; i++)
                    import_gate_bootstrapping_ciphertext_fromFile(answer_data, &finalresult7[i], params);

                for (int i = 0; i  32; i++)
                    import_gate_bootstrapping_ciphertext_fromFile(answer_data, &finalresult8[i], params);
                fclose(answer_data);

                decrypt and rebuild the answer
                int32_t int_answer = 0;
                for (int i = 0; i  32; i++)
                {
                    int ai = bootsSymDecrypt(&finalresult1[i], key)  0;
                    int_answer = (ai  i);
                }

                int32_t int_answer2 = 0;
                for (int i = 0; i  32; i++)
                {
                    int ai = bootsSymDecrypt(&finalresult2[i], key)  0;
                    int_answer2 = (ai  i);
                }

                int32_t int_answer3 = 0;
                for (int i = 0; i  32; i++)
                {
                    int ai = bootsSymDecrypt(&finalresult3[i], key)  0;
                    int_answer3 = (ai  i);
                }

                int32_t int_answer4 = 0;
                for (int i = 0; i  32; i++)
                {
                    int ai = bootsSymDecrypt(&finalresult4[i], key)  0;
                    int_answer4 = (ai  i);
                }

                int32_t int_answer5 = 0;
                for (int i = 0; i  32; i++)
                {
                    int ai = bootsSymDecrypt(&finalresult5[i], key)  0;
                    int_answer5 = (ai  i);
                }

                int32_t int_answer6 = 0;
                for (int i = 0; i  32; i++)
                {
                    int ai = bootsSymDecrypt(&finalresult6[i], key)  0;
                    int_answer6 = (ai  i);
                }

                int32_t int_answer7 = 0;
                for (int i = 0; i  32; i++)
                {
                    int ai = bootsSymDecrypt(&finalresult7[i], key)  0;
                    int_answer7 = (ai  i);
                }

                int32_t int_answer8 = 0;
                for (int i = 0; i  32; i++)
                {
                    int ai = bootsSymDecrypt(&finalresult8[i], key)  0;
                    int_answer8 = (ai  i);
                }

                stdstring binary1 = stdbitset32(int_answer).to_string();
                stdstring binary2 = stdbitset32(int_answer2).to_string();
                stdstring binary3 = stdbitset32(int_answer3).to_string();
                stdstring binary4 = stdbitset32(int_answer4).to_string();
                stdstring binary5 = stdbitset32(int_answer5).to_string();
                stdstring binary6 = stdbitset32(int_answer6).to_string();
                stdstring binary7 = stdbitset32(int_answer7).to_string();
                stdstring binary8 = stdbitset32(int_answer8).to_string();
                stdstring binary_combined = binary8 + binary7 + binary6 + binary5 + binary4 + binary3 + binary2 + binary1;

                int length = binary_combined.length();
                char char_array[length + 1];
                strcpy(char_array, binary_combined.c_str());
                int256_t total = 0;

                if (int_negative == 0  int_negative == 4)
                {
                    stdcout  The result in binary form is
                               n;
                    stdcout  binary_combined;
                }

                for (int i = 0; i  length; ++i)
                {
                    char charvalue = char_array[i];
                    stdstring stringvalue(1, charvalue);
                    int intvalue = stdstoi(stringvalue);
                    int256_t calc1 = total  2;
                    total = calc1 + intvalue;
                }

                if (int_negative == 1  int_negative == 2)
                {
                    int256_t invertedtotal = total  -1;
                    stdcout  The result in decimal form is
                               n;
                    stdcout  invertedtotal;
                }
                else
                {
                    stdcout  n
                               n;
                    stdcout  The result in decimal form is
                               n;
                    stdcout  total;
                }

                printf(n);
                printf(n);
                gettimeofday(&end, NULL);
                get_time = (end.tv_sec - start.tv_sec) + (end.tv_usec - start.tv_usec)  1.0E-6;
                printf(Computation Time %lf[sec]n, get_time);
                printf(n);

                printf(I hope you remembered what calculation you performed!n);

                clean up all pointers
                delete_gate_bootstrapping_ciphertext_array(32, finalresult1);
                delete_gate_bootstrapping_ciphertext_array(32, finalresult2);
                delete_gate_bootstrapping_ciphertext_array(32, finalresult3);
                delete_gate_bootstrapping_ciphertext_array(32, finalresult4);
                delete_gate_bootstrapping_ciphertext_array(32, finalresult5);
                delete_gate_bootstrapping_ciphertext_array(32, finalresult6);
                delete_gate_bootstrapping_ciphertext_array(32, finalresult7);
                delete_gate_bootstrapping_ciphertext_array(32, finalresult8);
                delete_gate_bootstrapping_secret_keyset(key);
                delete_gate_bootstrapping_secret_keyset(nbitkey);
            }

            else if (int_bit == 128)
            {
                LweSample finalresult = new_gate_bootstrapping_ciphertext_array(32, params);
                LweSample finalresult2 = new_gate_bootstrapping_ciphertext_array(32, params);
                LweSample finalresult3 = new_gate_bootstrapping_ciphertext_array(32, params);
                LweSample finalresult4 = new_gate_bootstrapping_ciphertext_array(32, params);

                export the 32 ciphertexts to a file (for the cloud)

                for (int i = 0; i  32; i++)
                    import_gate_bootstrapping_ciphertext_fromFile(answer_data, &finalresult[i], params);

                for (int i = 0; i  32; i++)
                    import_gate_bootstrapping_ciphertext_fromFile(answer_data, &finalresult2[i], params);

                for (int i = 0; i  32; i++)
                    import_gate_bootstrapping_ciphertext_fromFile(answer_data, &finalresult3[i], params);

                for (int i = 0; i  32; i++)
                    import_gate_bootstrapping_ciphertext_fromFile(answer_data, &finalresult4[i], params);

                fclose(answer_data);

                decrypt and rebuild the answer
                int32_t int_answer = 0;
                for (int i = 0; i  32; i++)
                {
                    int ai = bootsSymDecrypt(&finalresult[i], key)  0;
                    int_answer = (ai  i);
                }

                int32_t int_answer2 = 0;
                for (int i = 0; i  32; i++)
                {
                    int ai = bootsSymDecrypt(&finalresult2[i], key)  0;
                    int_answer2 = (ai  i);
                }

                int32_t int_answer3 = 0;
                for (int i = 0; i  32; i++)
                {
                    int ai = bootsSymDecrypt(&finalresult3[i], key)  0;
                    int_answer3 = (ai  i);
                }

                int32_t int_answer4 = 0;
                for (int i = 0; i  32; i++)
                {
                    int ai = bootsSymDecrypt(&finalresult4[i], key)  0;
                    int_answer4 = (ai  i);
                }

                stdstring binary1 = stdbitset32(int_answer).to_string();
                stdstring binary2 = stdbitset32(int_answer2).to_string();
                stdstring binary3 = stdbitset32(int_answer3).to_string();
                stdstring binary4 = stdbitset32(int_answer4).to_string();
                stdstring binary_combined = binary4 + binary3 + binary2 + binary1;

                int length = binary_combined.length();
                char char_array[length + 1];
                strcpy(char_array, binary_combined.c_str());
                int256_t total = 0;

                if (int_negative == 0  int_negative == 4)
                {
                    stdcout  The result in binary form is
                               n;
                    stdcout  binary_combined;
                }

                for (int i = 0; i  length; ++i)
                {
                    char charvalue = char_array[i];
                    stdstring stringvalue(1, charvalue);
                    int intvalue = stdstoi(stringvalue);
                    int256_t calc1 = total  2;
                    total = calc1 + intvalue;
                }

                if (int_negative == 1  int_negative == 2)
                {
                    int256_t invertedtotal = total  -1;
                    stdcout  The result in decimal form is
                               n;
                    stdcout  invertedtotal;
                }
                else
                {
                    stdcout  n
                               n;
                    stdcout  The result in decimal form is
                               n;
                    stdcout  total;
                }

                printf(n);
                printf(n);
                gettimeofday(&end, NULL);
                get_time = (end.tv_sec - start.tv_sec) + (end.tv_usec - start.tv_usec)  1.0E-6;
                printf(Computation Time %lf[sec]n, get_time);
                printf(n);
                printf(I hope you remembered what calculation you performed!n);

                clean up all pointers
                delete_gate_bootstrapping_ciphertext_array(32, finalresult);
                delete_gate_bootstrapping_ciphertext_array(32, finalresult2);
                delete_gate_bootstrapping_ciphertext_array(32, finalresult3);
                delete_gate_bootstrapping_ciphertext_array(32, finalresult4);
                delete_gate_bootstrapping_secret_keyset(key);
                delete_gate_bootstrapping_secret_keyset(nbitkey);
            }

            else if (int_bit == 64)
            {
                LweSample finalresult = new_gate_bootstrapping_ciphertext_array(32, params);
                LweSample finalresult2 = new_gate_bootstrapping_ciphertext_array(32, params);

                export the 32 ciphertexts to a file (for the cloud)

                for (int i = 0; i  32; i++)
                    import_gate_bootstrapping_ciphertext_fromFile(answer_data, &finalresult[i], params);

                for (int i = 0; i  32; i++)
                    import_gate_bootstrapping_ciphertext_fromFile(answer_data, &finalresult2[i], params);

                fclose(answer_data);

                decrypt and rebuild the answer
                int32_t int_answer = 0;
                for (int i = 0; i  32; i++)
                {
                    int ai = bootsSymDecrypt(&finalresult[i], key)  0;
                    int_answer = (ai  i);
                }

                int32_t int_answer2 = 0;
                for (int i = 0; i  32; i++)
                {
                    int ai = bootsSymDecrypt(&finalresult2[i], key)  0;
                    int_answer2 = (ai  i);
                }

                stdstring binary1 = stdbitset32(int_answer).to_string();
                stdstring binary2 = stdbitset32(int_answer2).to_string();
                stdstring binary_combined = binary2 + binary1;

                int length = binary_combined.length();
                char char_array[length + 1];
                strcpy(char_array, binary_combined.c_str());
                int256_t total = 0;
                if (int_negative == 0  int_negative == 4)
                {
                    stdcout  The result in binary form is
                               n;
                    stdcout  binary_combined;
                }

                for (int i = 0; i  length; ++i)
                {
                    char charvalue = char_array[i];
                    stdstring stringvalue(1, charvalue);
                    int intvalue = stdstoi(stringvalue);
                    int256_t calc1 = total  2;
                    total = calc1 + intvalue;
                }

                if (int_negative == 1  int_negative == 2)
                {
                    int256_t invertedtotal = total  -1;
                    stdcout  The result in decimal form is
                               n;
                    stdcout  invertedtotal;
                }
                else
                {
                    stdcout  n
                               n;
                    stdcout  The result in decimal form is
                               n;
                    stdcout  total;
                }

                printf(n);
                printf(n);
                gettimeofday(&end, NULL);
                get_time = (end.tv_sec - start.tv_sec) + (end.tv_usec - start.tv_usec)  1.0E-6;
                printf(Computation Time %lf[sec]n, get_time);
                printf(n);

                printf(I hope you remembered what calculation you performed!n);

                clean up all pointers
                delete_gate_bootstrapping_ciphertext_array(32, finalresult);
                delete_gate_bootstrapping_ciphertext_array(32, finalresult2);
                delete_gate_bootstrapping_secret_keyset(key);
                delete_gate_bootstrapping_secret_keyset(nbitkey);
            }
        }
    }

    else if (strcmp(mpioropenmp, MPI) == 0)
    {
        printf(MPI is being computedn);
         Addition
         reads the secret key from file
        FILE secret_key = fopen(secret.key, rb);
        TFheGateBootstrappingSecretKeySet key = new_tfheGateBootstrappingSecretKeySet_fromFile(secret_key);
        fclose(secret_key);

         reads the nbit key from file
        FILE nbit_key = fopen(nbit.key, rb);
        TFheGateBootstrappingSecretKeySet nbitkey = new_tfheGateBootstrappingSecretKeySet_fromFile(nbit_key);
        fclose(nbit_key);

         if necessary, the params are inside the key
        const TFheGateBootstrappingParameterSet params = key-params;

         if necessary, the params are inside the key
        const TFheGateBootstrappingParameterSet nbitparams = nbitkey-params;

        struct timeval start, end;
        double get_time;
        gettimeofday(&start, NULL);

        read the 16 ciphertexts of the result
        LweSample negative = new_gate_bootstrapping_ciphertext_array(32, nbitparams);
        LweSample bit = new_gate_bootstrapping_ciphertext_array(32, nbitparams);

        FILE answer_data = fopen(answer.data, rb);
        for (int i = 0; i  32; i++)
            import_gate_bootstrapping_ciphertext_fromFile(answer_data, &negative[i], nbitparams);

        for (int i = 0; i  32; i++)
            import_gate_bootstrapping_ciphertext_fromFile(answer_data, &bit[i], nbitparams);

        decrypt and rebuild the answer

         negativity code
        int32_t int_negative = 0;
        for (int i = 0; i  32; i++)
        {
            int ai = bootsSymDecrypt(&negative[i], nbitkey)  0;
            int_negative = (ai  i);
        }
        stdcout  Negative   int_negative  n;

         TODO Obtain opcode from file
        int32_t int_op;
        FILE fptr;
        fptr = fopen(operator.txt, r);
        fscanf(fptr, %d, &int_op);
        stdcout  Opcode   int_op  n
                   n;

         Bit count
        int32_t int_bit = 0;
        for (int i = 0; i  32; i++)
        {
            int ai = bootsSymDecrypt(&bit[i], nbitkey)  0;
            int_bit = (ai  i);
        }

        if (int_op == 1)
        {
            stdcout  Result for   int_bit   bit Addition computation
                       n
                       n;

            if (int_bit == 32)
            {
                LweSample result = new_gate_bootstrapping_ciphertext_array(32, params);

                for (int i = 0; i  32; i++)
                    import_gate_bootstrapping_ciphertext_fromFile(answer_data, &result[i], params);
                fclose(answer_data);

                 decrypt and rebuild the answer
                int32_t int_answer1 = 0;
                for (int i = 0; i  32; i++)
                {
                    int ai = bootsSymDecrypt(&result[i], key)  0;
                    int_answer1 = (ai  i);
                }

                stdstring binary1 = stdbitset32(int_answer1).to_string();

                stdstring binary_combined = binary1;

                if (int_negative != 4)
                {
                    stdcout  The result in binary form is
                               n;
                    stdcout  binary_combined  n;
                }

                Length is the number of bits
                int length = binary_combined.length();
                int lastchar = length - 1;
                char char_array[length + 1];
                strcpy(char_array, binary_combined.c_str());

                int256_t total = 0;

                Positive or negative checking
                char charvalue = char_array[0];
                stdstring stringvalue(1, charvalue);
                int intvalue = stdstoi(stringvalue);

                if the result is either double positive or negative
                if (int_negative == 0  int_negative == 4)
                {
                    for (int i = 0; i  length; ++i)
                    {
                        char charvalue = char_array[i];
                        stdstring stringvalue(1, charvalue);
                        int intvalue = stdstoi(stringvalue);
                        int256_t calc1 = total  2;
                        total = calc1 + intvalue;
                    }
                    if either one of the values are negative, check for negative
                }
                else if (int_negative != 4)
                {
                    if (length == 32 && intvalue == 1)
                    {
                        int256_t lastcharvalue = 2;
                        for (int i = 0; i  lastchar - 1; ++i)
                        {
                            lastcharvalue = lastcharvalue  2;
                        }

                        for (int i = 1; i  length; ++i)
                        {
                            char charvalue = char_array[i];
                            stdstring stringvalue(1, charvalue);
                            int intvalue = stdstoi(stringvalue);
                            int256_t calc1 = total  2;
                            total = calc1 + intvalue;
                        }
                        total = total - lastcharvalue;
                    }
                    else
                    {
                        for (int i = 0; i  length; ++i)
                        {
                            char charvalue = char_array[i];
                            stdstring stringvalue(1, charvalue);
                            int intvalue = stdstoi(stringvalue);
                            int256_t calc1 = total  2;
                            total = calc1 + intvalue;
                        }
                    }
                }

                if its a double negative then inverse it
                if (int_negative == 4)
                {
                    int256_t invertedtotal = total  -1;
                    stdcout  The result in decimal form is
                               n;
                    stdcout  invertedtotal;
                }
                else
                {
                    stdcout  n
                               n;
                    stdcout  The result in decimal form is
                               n;
                    stdcout  total;
                }

                printf(n);
                printf(n);
                gettimeofday(&end, NULL);
                get_time = (end.tv_sec - start.tv_sec) + (end.tv_usec - start.tv_usec)  1.0E-6;
                printf(Computation Time for MPI %lf[sec]n, get_time);
                printf(n);

                printf(I hope you remembered what calculation you performed!n);

                 clean up all pointers
                delete_gate_bootstrapping_ciphertext_array(32, result);
                delete_gate_bootstrapping_secret_keyset(key);
                delete_gate_bootstrapping_secret_keyset(nbitkey);
            }
            else if (int_bit == 64)
            {
                printf(Verifying 64 bit addition on MPI n);
                 read the 32 ciphertexts of the result
                LweSample result = new_gate_bootstrapping_ciphertext_array(32, params);
                LweSample result2 = new_gate_bootstrapping_ciphertext_array(32, params);
                LweSample result3 = new_gate_bootstrapping_ciphertext_array(32, params);
                LweSample carry1 = new_gate_bootstrapping_ciphertext_array(32, params);

                 export the 64 ciphertexts to a file (for the cloud)
                for (int i = 0; i  32; i++)
                    import_gate_bootstrapping_ciphertext_fromFile(answer_data, &result[i], params);
                for (int i = 0; i  32; i++)
                    import_gate_bootstrapping_ciphertext_fromFile(answer_data, &carry1[i], params);
                for (int i = 0; i  32; i++)
                    import_gate_bootstrapping_ciphertext_fromFile(answer_data, &result2[i], params);
                for (int i = 0; i  32; i++)
                    import_gate_bootstrapping_ciphertext_fromFile(answer_data, &result3[i], params);

                fclose(answer_data);
                printf(Exportingn);

                 decrypt and rebuild the answer
                int32_t firstResult = 0;
                for (int i = 0; i  32; i++)
                {
                    int ai = bootsSymDecrypt(&result[i], key)  0;
                    firstResult = (ai  i);
                }

                int32_t VerifyCarry = 0;
                for (int i = 0; i  32; i++)
                {
                    int ai = bootsSymDecrypt(&carry1[i], key)  0;
                    VerifyCarry = (ai  i);
                }

                int32_t ifCarryisZero = 0;
                for (int i = 0; i  32; i++)
                {
                    int ai = bootsSymDecrypt(&result2[i], key)  0;
                    ifCarryisZero = (ai  i);
                }

                int32_t ifCarryisOne = 0;
                for (int i = 0; i  32; i++)
                {
                    int ai = bootsSymDecrypt(&result3[i], key)  0;
                    ifCarryisOne = (ai  i);
                }

                last 32
                when testing == ciphertext1
                stdstring binary1 = stdbitset32(firstResult).to_string();
                carry1 in binary
                stdstring binary4 = stdbitset32(VerifyCarry).to_string();
                32 if carry1 = 0
                stdstring binary2 = stdbitset32(ifCarryisZero).to_string();
                32 if carry1 = 1
                stdstring binary3 = stdbitset32(ifCarryisOne).to_string();

                stdcout  First 32 (Binary1) ;
                stdcout  binary1;
                stdcout  n
                           n;
                stdcout  Carry1 (Carry1 in binary) ;
                stdcout  binary4;
                stdcout  n
                           n;
                stdcout  If carry is zero (Binary2) ;
                stdcout  binary2;
                stdcout  n
                           n;
                stdcout  If carry is one (Binary3) ;
                stdcout  binary3;
                stdcout  n
                           n;

                printf(Carry1 is %dn, VerifyCarry);

                if (VerifyCarry == 0)
                {
                    printf(Binary2 + Binary1 %dn, VerifyCarry);
                    stdstring binary_combined = binary2 + binary1;

                    if (int_negative != 4)
                    {
                        stdcout  The result in binary form is
                                   n;
                        stdcout  binary_combined  n;
                    }

                    int length = binary_combined.length();
                    int lastchar = length - 1;
                    char char_array[length + 1];
                    strcpy(char_array, binary_combined.c_str());

                    int256_t total = 0;

                    Positive or negative checking
                    char charvalue = char_array[0];
                    stdstring stringvalue(1, charvalue);
                    int intvalue = stdstoi(stringvalue);
                    if the result is either double positive or negative
                    if (int_negative == 0  int_negative == 4)
                    {
                        for (int i = 0; i  length; ++i)
                        {
                            char charvalue = char_array[i];
                            stdstring stringvalue(1, charvalue);
                            int intvalue = stdstoi(stringvalue);
                            int256_t calc1 = total  2;
                            total = calc1 + intvalue;
                        }
                        if either one of the values are negative, check for negative
                    }
                    else if (int_negative != 4)
                    {
                        if (length == 64 && intvalue == 1)
                        {
                            int256_t lastcharvalue = 2;
                            for (int i = 0; i  lastchar - 1; ++i)
                            {
                                lastcharvalue = lastcharvalue  2;
                            }

                            for (int i = 1; i  length; ++i)
                            {
                                char charvalue = char_array[i];
                                stdstring stringvalue(1, charvalue);
                                int intvalue = stdstoi(stringvalue);
                                int256_t calc1 = total  2;
                                total = calc1 + intvalue;
                            }

                            total = total - lastcharvalue;
                        }
                        else
                        {
                            for (int i = 0; i  length; ++i)
                            {
                                char charvalue = char_array[i];
                                stdstring stringvalue(1, charvalue);
                                int intvalue = stdstoi(stringvalue);
                                int256_t calc1 = total  2;
                                total = calc1 + intvalue;
                            }
                        }
                    }
                    if (int_negative == 4)
                    {
                        int256_t invertedtotal = total  -1;
                        stdcout  The result in decimal form is
                                   n;
                        stdcout  invertedtotal;
                    }
                    else
                    {
                        stdcout  n
                                   n;
                        stdcout  The result in decimal form is
                                   n;
                        stdcout  total;
                    }
                }

                if carry = 1
                else if (VerifyCarry == 1)
                {
                    printf(Binary3 + Binary1 %dn, VerifyCarry);
                    stdstring binary_combined2 = binary3 + binary1;

                    if (int_negative != 4)
                    {
                        stdcout  The result in binary form is
                                   n;
                        stdcout  binary_combined2  n;
                    }

                    int length = binary_combined2.length();
                    int lastchar = length - 1;
                    char char_array[length + 1];
                    strcpy(char_array, binary_combined2.c_str());

                    int256_t total = 0;

                    Positive or negative checking
                    char charvalue = char_array[0];
                    stdstring stringvalue(1, charvalue);
                    int intvalue = stdstoi(stringvalue);
                    if the result is either double positive or negative
                    if (int_negative == 0  int_negative == 4)
                    {
                        for (int i = 0; i  length; ++i)
                        {
                            char charvalue = char_array[i];
                            stdstring stringvalue(1, charvalue);
                            int intvalue = stdstoi(stringvalue);
                            int256_t calc1 = total  2;
                            total = calc1 + intvalue;
                        }
                        if either one of the values are negative, check for negative
                    }
                    else if (int_negative != 4)
                    {
                        if (length == 64 && intvalue == 1)
                        {
                            int256_t lastcharvalue = 2;
                            for (int i = 0; i  lastchar - 1; ++i)
                            {
                                lastcharvalue = lastcharvalue  2;
                            }

                            for (int i = 1; i  length; ++i)
                            {
                                char charvalue = char_array[i];
                                stdstring stringvalue(1, charvalue);
                                int intvalue = stdstoi(stringvalue);
                                int256_t calc1 = total  2;
                                total = calc1 + intvalue;
                            }

                            total = total - lastcharvalue;
                        }
                        else
                        {
                            for (int i = 0; i  length; ++i)
                            {
                                char charvalue = char_array[i];
                                stdstring stringvalue(1, charvalue);
                                int intvalue = stdstoi(stringvalue);
                                int256_t calc1 = total  2;
                                total = calc1 + intvalue;
                            }
                        }
                    }
                    if (int_negative == 4)
                    {
                        int256_t invertedtotal = total  -1;
                        stdcout  The result in decimal form is
                                   n;
                        stdcout  invertedtotal;
                    }
                    else
                    {
                        stdcout  n
                                   n;
                        stdcout  The result in decimal form is
                                   n;
                        stdcout  total;
                    }
                }

                else
                {
                    printf(Failn);
                }

                printf(n);
                printf(n);
                printf(n);
                printf(I hope you remembered what calculation you performed!n);

                delete_gate_bootstrapping_ciphertext_array(32, result);
                delete_gate_bootstrapping_ciphertext_array(32, result2);
                delete_gate_bootstrapping_ciphertext_array(32, result3);
                delete_gate_bootstrapping_ciphertext_array(32, carry1);
                delete_gate_bootstrapping_secret_keyset(key);
                delete_gate_bootstrapping_secret_keyset(nbitkey);
            }

            else if (int_bit == 128)
            {
                printf(Verifying 128 bit addition on MPI n);
                 read the 32 ciphertexts of the result
                LweSample result = new_gate_bootstrapping_ciphertext_array(32, params);
                LweSample result2 = new_gate_bootstrapping_ciphertext_array(32, params);
                LweSample result3 = new_gate_bootstrapping_ciphertext_array(32, params);
                LweSample result4 = new_gate_bootstrapping_ciphertext_array(32, params);
                LweSample result5 = new_gate_bootstrapping_ciphertext_array(32, params);
                LweSample result6 = new_gate_bootstrapping_ciphertext_array(32, params);
                LweSample result7 = new_gate_bootstrapping_ciphertext_array(32, params);
                LweSample carry1 = new_gate_bootstrapping_ciphertext_array(32, params);
                LweSample carry2 = new_gate_bootstrapping_ciphertext_array(32, params);
                LweSample carry3 = new_gate_bootstrapping_ciphertext_array(32, params);
                LweSample carry4 = new_gate_bootstrapping_ciphertext_array(32, params);
                LweSample carry5 = new_gate_bootstrapping_ciphertext_array(32, params);
                LweSample carry6 = new_gate_bootstrapping_ciphertext_array(32, params);
                LweSample carry7 = new_gate_bootstrapping_ciphertext_array(32, params);

                 export the 64 ciphertexts to a file (for the cloud)
                for (int i = 0; i  32; i++)
                    import_gate_bootstrapping_ciphertext_fromFile(answer_data, &result[i], params);
                for (int i = 0; i  32; i++)
                    import_gate_bootstrapping_ciphertext_fromFile(answer_data, &carry1[i], params);
                for (int i = 0; i  32; i++)
                    import_gate_bootstrapping_ciphertext_fromFile(answer_data, &carry2[i], params);
                for (int i = 0; i  32; i++)
                    import_gate_bootstrapping_ciphertext_fromFile(answer_data, &carry3[i], params);
                for (int i = 0; i  32; i++)
                    import_gate_bootstrapping_ciphertext_fromFile(answer_data, &carry4[i], params);
                for (int i = 0; i  32; i++)
                    import_gate_bootstrapping_ciphertext_fromFile(answer_data, &carry5[i], params);
                for (int i = 0; i  32; i++)
                    import_gate_bootstrapping_ciphertext_fromFile(answer_data, &carry6[i], params);
                for (int i = 0; i  32; i++)
                    import_gate_bootstrapping_ciphertext_fromFile(answer_data, &carry7[i], params);
                for (int i = 0; i  32; i++)
                    import_gate_bootstrapping_ciphertext_fromFile(answer_data, &result2[i], params);
                for (int i = 0; i  32; i++)
                    import_gate_bootstrapping_ciphertext_fromFile(answer_data, &result3[i], params);
                for (int i = 0; i  32; i++)
                    import_gate_bootstrapping_ciphertext_fromFile(answer_data, &result4[i], params);
                for (int i = 0; i  32; i++)
                    import_gate_bootstrapping_ciphertext_fromFile(answer_data, &result5[i], params);
                for (int i = 0; i  32; i++)
                    import_gate_bootstrapping_ciphertext_fromFile(answer_data, &result6[i], params);
                for (int i = 0; i  32; i++)
                    import_gate_bootstrapping_ciphertext_fromFile(answer_data, &result7[i], params);
                fclose(answer_data);
                printf(Exportingn);

                 decrypt and rebuild the answer
                int32_t firstResult = 0;
                for (int i = 0; i  32; i++)
                {
                    int ai = bootsSymDecrypt(&result[i], key)  0;
                    firstResult = (ai  i);
                }
                int32_t firstcarry = 0;
                for (int i = 0; i  32; i++)
                {
                    int ai = bootsSymDecrypt(&carry1[i], key)  0;
                    firstcarry = (ai  i);
                }
                int32_t secondCarryZero = 0;
                for (int i = 0; i  32; i++)
                {
                    int ai = bootsSymDecrypt(&carry2[i], key)  0;
                    secondCarryZero = (ai  i);
                }
                int32_t secondCarryOne = 0;
                for (int i = 0; i  32; i++)
                {
                    int ai = bootsSymDecrypt(&carry3[i], key)  0;
                    secondCarryOne = (ai  i);
                }

                int32_t thirdCarryZero = 0;
                for (int i = 0; i  32; i++)
                {
                    int ai = bootsSymDecrypt(&carry4[i], key)  0;
                    thirdCarryZero = (ai  i);
                }
                int32_t thirdCarryOne = 0;
                for (int i = 0; i  32; i++)
                {
                    int ai = bootsSymDecrypt(&carry5[i], key)  0;
                    thirdCarryOne = (ai  i);
                }
                int32_t fourthCarryZero = 0;
                for (int i = 0; i  32; i++)
                {
                    int ai = bootsSymDecrypt(&carry6[i], key)  0;
                    fourthCarryZero = (ai  i);
                }
                int32_t fourthCarryOne = 0;
                for (int i = 0; i  32; i++)
                {
                    int ai = bootsSymDecrypt(&carry7[i], key)  0;
                    fourthCarryOne = (ai  i);
                }
                
                int32_t secondResultZero = 0;
                for (int i = 0; i  32; i++)
                {
                    int ai = bootsSymDecrypt(&result2[i], key)  0;
                    secondResultZero = (ai  i);
                }
                int32_t secondResultOne = 0;
                for (int i = 0; i  32; i++)
                {
                    int ai = bootsSymDecrypt(&result3[i], key)  0;
                    secondResultOne = (ai  i);
                }
                int32_t thirdResultZero = 0;
                for (int i = 0; i  32; i++)
                {
                    int ai = bootsSymDecrypt(&result4[i], key)  0;
                    thirdResultZero = (ai  i);
                }
                int32_t thirdResultOne = 0;
                for (int i = 0; i  32; i++)
                {
                    int ai = bootsSymDecrypt(&result5[i], key)  0;
                    thirdResultOne = (ai  i);
                }
                int32_t fourthResultZero = 0;
                for (int i = 0; i  32; i++)
                {
                    int ai = bootsSymDecrypt(&result6[i], key)  0;
                    fourthResultZero = (ai  i);
                }
                int32_t fourthResultOne = 0;
                for (int i = 0; i  32; i++)
                {
                    int ai = bootsSymDecrypt(&result7[i], key)  0;
                    fourthResultOne = (ai  i);
                }

                printf(carry1  %dn, firstcarry);
                printf(carry2 secondCarryZero  %dn, secondCarryZero);
                printf(carry3 secondCarryOne  %dn, secondCarryOne);
                printf(carry4 thirdCarryZero  %dn, thirdCarryZero);
                printf(carry5 thirdCarryOne  %dn, thirdCarryOne);
                printf(carry6 fourthCarryZero  %dn, fourthCarryZero);
                printf(carry7 fourthCarryOne  %dn, fourthCarryOne);

                last 32 bits
                stdstring binary1 = stdbitset32(firstResult).to_string();

                carry of ciphertext1 and 9
                stdstring binary2 = stdbitset32(firstcarry).to_string();

                carry=0
                carry of ciphertext2 and 10
                stdstring binary3 = stdbitset32(secondCarryZero).to_string();

                carry=1
                carry of ciphertext2 and 10
                stdstring binary4 = stdbitset32(secondCarryOne).to_string();

                carry=0
                carry of ciphertext3 and 11
                stdstring binary5 = stdbitset32(thirdCarryZero).to_string();

                carry=1
                carry of ciphertext3 and 11
                stdstring binary6 = stdbitset32(thirdCarryOne).to_string();

                carry=0
                carry of ciphertext4 and 12
                stdstring binary7 = stdbitset32(fourthCarryZero).to_string();

                carry=1
                carry of ciphertext4 and 12
                stdstring binary8 = stdbitset32(fourthCarryOne).to_string();

                result of next 32 if carry is 0        2
                stdstring binary9 = stdbitset32(secondResultZero).to_string();
                result of next 32 if carry is 1        3
                stdstring binary10 = stdbitset32(secondResultOne).to_string();
                64 - 96 if carry is 0                  4
                stdstring binary11 = stdbitset32(thirdResultZero).to_string();
                64 - 96 if carry is 1                  5
                stdstring binary12 = stdbitset32(thirdResultOne).to_string();
                96 -128 if carry is 0                  6
                stdstring binary13 = stdbitset32(fourthResultZero).to_string();
                96 -128 carry is 1                     7
                stdstring binary14 = stdbitset32(fourthResultOne).to_string();


                printf(Verifying.............................. n);

                first carry =0 - 2nd carry = 0 - third carry = 0 (13 11 9) check
                if (firstcarry == 0 && secondCarryZero == 0 && thirdCarryZero == 0)
                {
                    printf(13 + 11 + 9 + 1 n);
                    stdstring binary_combined1 = binary13 + binary11 + binary9 + binary1;
                    if (int_negative != 4)
                    {
                        stdcout  The result in binary form is
                                   n;
                        stdcout  binary_combined1  n;
                    }
                    Length is the number of bits
                    int length = binary_combined1.length();
                    int lastchar = length - 1;
                    char char_array[length + 1];
                    strcpy(char_array, binary_combined1.c_str());

                    int256_t total = 0;

                    Positive or negative checking
                    char charvalue = char_array[0];
                    stdstring stringvalue(1, charvalue);
                    int intvalue = stdstoi(stringvalue);

                    if the result is either double positive or negative
                    if (int_negative == 0  int_negative == 4)
                    {
                        for (int i = 0; i  length; ++i)
                        {
                            char charvalue = char_array[i];
                            stdstring stringvalue(1, charvalue);
                            int intvalue = stdstoi(stringvalue);
                            int256_t calc1 = total  2;
                            total = calc1 + intvalue;
                        }
                    }
                    if either one of the values are negative, check for negative
                    else if (int_negative != 4)
                    {
                        if (length == 128 && intvalue == 1)
                        {
                            int256_t lastcharvalue = 2;

                            for (int i = 0; i  lastchar - 1; ++i)
                            {
                                lastcharvalue = lastcharvalue  2;
                            }

                            for (int i = 1; i  length; ++i)
                            {
                                char charvalue = char_array[i];
                                stdstring stringvalue(1, charvalue);
                                int intvalue = stdstoi(stringvalue);
                                int256_t calc1 = total  2;
                                total = calc1 + intvalue;
                            }
                            total = total - lastcharvalue;
                        }
                        else
                        {
                            for (int i = 0; i  length; ++i)
                            {
                                char charvalue = char_array[i];
                                stdstring stringvalue(1, charvalue);
                                int intvalue = stdstoi(stringvalue);
                                int256_t calc1 = total  2;
                                total = calc1 + intvalue;
                            }
                        }
                    }
                    if (int_negative == 4)
                    {
                        int256_t invertedtotal = total  -1;
                        stdcout  The result in decimal form is
                                   n;
                        stdcout  invertedtotal;
                    }
                    else
                    {
                        stdcout  n
                                   n;
                        stdcout  The result in decimal form is
                                   n;
                        stdcout  total;
                    }

                    printf(n);
                    printf(n);
                    gettimeofday(&end, NULL);
                    get_time = (end.tv_sec - start.tv_sec) + (end.tv_usec - start.tv_usec)  1.0E-6;
                    printf(Computation Time %lf[sec]n, get_time);
                    printf(n);

                    printf(I hope you remembered what calculation you performed!n);
                }

                first carry =0 - 2nd carry = 0 - third carry = 1 (14 11 9) check
                else if (firstcarry == 0 && secondCarryZero == 0 && thirdCarryZero == 1)
                {
                    printf(14 + 11 + 9 + 1 n);
                    stdstring binary_combined2 = binary14 + binary11 + binary9 + binary1;
                    if (int_negative != 4)
                    {
                        stdcout  The result in binary form is
                                   n;
                        stdcout  binary_combined2  n;
                    }
                    Length is the number of bits
                    int length = binary_combined2.length();
                    int lastchar = length - 1;
                    char char_array[length + 1];
                    strcpy(char_array, binary_combined2.c_str());

                    int256_t total = 0;

                    Positive or negative checking
                    char charvalue = char_array[0];
                    stdstring stringvalue(1, charvalue);
                    int intvalue = stdstoi(stringvalue);

                    if the result is either double positive or negative
                    if (int_negative == 0  int_negative == 4)
                    {
                        for (int i = 0; i  length; ++i)
                        {
                            char charvalue = char_array[i];
                            stdstring stringvalue(1, charvalue);
                            int intvalue = stdstoi(stringvalue);
                            int256_t calc1 = total  2;
                            total = calc1 + intvalue;
                        }
                    }
                    if either one of the values are negative, check for negative
                    else if (int_negative != 4)
                    {
                        if (length == 128 && intvalue == 1)
                        {
                            int256_t lastcharvalue = 2;

                            for (int i = 0; i  lastchar - 1; ++i)
                            {
                                lastcharvalue = lastcharvalue  2;
                            }

                            for (int i = 1; i  length; ++i)
                            {
                                char charvalue = char_array[i];
                                stdstring stringvalue(1, charvalue);
                                int intvalue = stdstoi(stringvalue);
                                int256_t calc1 = total  2;
                                total = calc1 + intvalue;
                            }
                            total = total - lastcharvalue;
                        }
                        else
                        {
                            for (int i = 0; i  length; ++i)
                            {
                                char charvalue = char_array[i];
                                stdstring stringvalue(1, charvalue);
                                int intvalue = stdstoi(stringvalue);
                                int256_t calc1 = total  2;
                                total = calc1 + intvalue;
                            }
                        }
                    }
                    if (int_negative == 4)
                    {
                        int256_t invertedtotal = total  -1;
                        stdcout  The result in decimal form is
                                   n;
                        stdcout  invertedtotal;
                    }
                    else
                    {
                        stdcout  n
                                   n;
                        stdcout  The result in decimal form is
                                   n;
                        stdcout  total;
                    }

                    printf(n);
                    printf(n);
                    gettimeofday(&end, NULL);
                    get_time = (end.tv_sec - start.tv_sec) + (end.tv_usec - start.tv_usec)  1.0E-6;
                    printf(Computation Time %lf[sec]n, get_time);
                    printf(n);

                    printf(I hope you remembered what calculation you performed!n);
                }

                first carry =0 - 2nd carry = 1 - third carry = 1 (14 12 9) check
                else if (firstcarry == 0 && secondCarryZero == 1 && thirdCarryOne == 1)
                {
                    printf(14 + 12 + 9 + 1 n);
                    stdstring binary_combined3 = binary14 + binary12 + binary9 + binary1;
                    if (int_negative != 4)
                    {
                        stdcout  The result in binary form is
                                   n;
                        stdcout  binary_combined3  n;
                    }
                    Length is the number of bits
                    int length = binary_combined3.length();
                    int lastchar = length - 1;
                    char char_array[length + 1];
                    strcpy(char_array, binary_combined3.c_str());

                    int256_t total = 0;

                    Positive or negative checking
                    char charvalue = char_array[0];
                    stdstring stringvalue(1, charvalue);
                    int intvalue = stdstoi(stringvalue);

                    if the result is either double positive or negative
                    if (int_negative == 0  int_negative == 4)
                    {
                        for (int i = 0; i  length; ++i)
                        {
                            char charvalue = char_array[i];
                            stdstring stringvalue(1, charvalue);
                            int intvalue = stdstoi(stringvalue);
                            int256_t calc1 = total  2;
                            total = calc1 + intvalue;
                        }
                    }
                    if either one of the values are negative, check for negative
                    else if (int_negative != 4)
                    {
                        if (length == 128 && intvalue == 1)
                        {
                            int256_t lastcharvalue = 2;

                            for (int i = 0; i  lastchar - 1; ++i)
                            {
                                lastcharvalue = lastcharvalue  2;
                            }

                            for (int i = 1; i  length; ++i)
                            {
                                char charvalue = char_array[i];
                                stdstring stringvalue(1, charvalue);
                                int intvalue = stdstoi(stringvalue);
                                int256_t calc1 = total  2;
                                total = calc1 + intvalue;
                            }
                            total = total - lastcharvalue;
                        }
                        else
                        {
                            for (int i = 0; i  length; ++i)
                            {
                                char charvalue = char_array[i];
                                stdstring stringvalue(1, charvalue);
                                int intvalue = stdstoi(stringvalue);
                                int256_t calc1 = total  2;
                                total = calc1 + intvalue;
                            }
                        }
                    }
                    if (int_negative == 4)
                    {
                        int256_t invertedtotal = total  -1;
                        stdcout  The result in decimal form is
                                   n;
                        stdcout  invertedtotal;
                    }
                    else
                    {
                        stdcout  n
                                   n;
                        stdcout  The result in decimal form is
                                   n;
                        stdcout  total;
                    }

                    printf(n);
                    printf(n);
                    gettimeofday(&end, NULL);
                    get_time = (end.tv_sec - start.tv_sec) + (end.tv_usec - start.tv_usec)  1.0E-6;
                    printf(Computation Time %lf[sec]n, get_time);
                    printf(n);

                    printf(I hope you remembered what calculation you performed!n);
                }

                first carry =1 - 2nd carry = 1 - third carry = 1 (14 12 10) check
                else if (firstcarry == 1 && secondCarryOne == 1 && thirdCarryOne == 1)
                {
                    printf(14 + 12 + 10 + 1 n);
                    stdstring binary_combined4 = binary14 + binary12 + binary10 + binary1;
                    if (int_negative != 4)
                    {
                        stdcout  The result in binary form is
                                   n;
                        stdcout  binary_combined4  n;
                    }
                    Length is the number of bits
                    int length = binary_combined4.length();
                    int lastchar = length - 1;
                    char char_array[length + 1];
                    strcpy(char_array, binary_combined4.c_str());

                    int256_t total = 0;

                    Positive or negative checking
                    char charvalue = char_array[0];
                    stdstring stringvalue(1, charvalue);
                    int intvalue = stdstoi(stringvalue);

                    if the result is either double positive or negative
                    if (int_negative == 0  int_negative == 4)
                    {
                        for (int i = 0; i  length; ++i)
                        {
                            char charvalue = char_array[i];
                            stdstring stringvalue(1, charvalue);
                            int intvalue = stdstoi(stringvalue);
                            int256_t calc1 = total  2;
                            total = calc1 + intvalue;
                        }
                    }
                    if either one of the values are negative, check for negative
                    else if (int_negative != 4)
                    {
                        if (length == 128 && intvalue == 1)
                        {
                            int256_t lastcharvalue = 2;

                            for (int i = 0; i  lastchar - 1; ++i)
                            {
                                lastcharvalue = lastcharvalue  2;
                            }

                            for (int i = 1; i  length; ++i)
                            {
                                char charvalue = char_array[i];
                                stdstring stringvalue(1, charvalue);
                                int intvalue = stdstoi(stringvalue);
                                int256_t calc1 = total  2;
                                total = calc1 + intvalue;
                            }
                            total = total - lastcharvalue;
                        }
                        else
                        {
                            for (int i = 0; i  length; ++i)
                            {
                                char charvalue = char_array[i];
                                stdstring stringvalue(1, charvalue);
                                int intvalue = stdstoi(stringvalue);
                                int256_t calc1 = total  2;
                                total = calc1 + intvalue;
                            }
                        }
                    }
                    if (int_negative == 4)
                    {
                        int256_t invertedtotal = total  -1;
                        stdcout  The result in decimal form is
                                   n;
                        stdcout  invertedtotal;
                    }
                    else
                    {
                        stdcout  n
                                   n;
                        stdcout  The result in decimal form is
                                   n;
                        stdcout  total;
                    }

                    printf(n);
                    printf(n);
                    gettimeofday(&end, NULL);
                    get_time = (end.tv_sec - start.tv_sec) + (end.tv_usec - start.tv_usec)  1.0E-6;
                    printf(Computation Time %lf[sec]n, get_time);
                    printf(n);

                    printf(I hope you remembered what calculation you performed!n);
                }

                first carry =1 - 2nd carry = 0 - third carry = 0 (13 11 10) check
                else if (firstcarry == 1 && secondCarryOne == 0 && thirdCarryZero == 0)
                {
                    printf(13 + 11 + 10 + 1 n);
                    stdstring binary_combined5 = binary13 + binary11 + binary10 + binary1;
                    if (int_negative != 4)
                    {
                        stdcout  The result in binary form is
                                   n;
                        stdcout  binary_combined5  n;
                    }
                    Length is the number of bits
                    int length = binary_combined5.length();
                    int lastchar = length - 1;
                    char char_array[length + 1];
                    strcpy(char_array, binary_combined5.c_str());

                    int256_t total = 0;

                    Positive or negative checking
                    char charvalue = char_array[0];
                    stdstring stringvalue(1, charvalue);
                    int intvalue = stdstoi(stringvalue);

                    if the result is either double positive or negative
                    if (int_negative == 0  int_negative == 4)
                    {
                        for (int i = 0; i  length; ++i)
                        {
                            char charvalue = char_array[i];
                            stdstring stringvalue(1, charvalue);
                            int intvalue = stdstoi(stringvalue);
                            int256_t calc1 = total  2;
                            total = calc1 + intvalue;
                        }
                    }
                    if either one of the values are negative, check for negative
                    else if (int_negative != 4)
                    {
                        if (length == 128 && intvalue == 1)
                        {
                            int256_t lastcharvalue = 2;

                            for (int i = 0; i  lastchar - 1; ++i)
                            {
                                lastcharvalue = lastcharvalue  2;
                            }

                            for (int i = 1; i  length; ++i)
                            {
                                char charvalue = char_array[i];
                                stdstring stringvalue(1, charvalue);
                                int intvalue = stdstoi(stringvalue);
                                int256_t calc1 = total  2;
                                total = calc1 + intvalue;
                            }
                            total = total - lastcharvalue;
                        }
                        else
                        {
                            for (int i = 0; i  length; ++i)
                            {
                                char charvalue = char_array[i];
                                stdstring stringvalue(1, charvalue);
                                int intvalue = stdstoi(stringvalue);
                                int256_t calc1 = total  2;
                                total = calc1 + intvalue;
                            }
                        }
                    }
                    if (int_negative == 4)
                    {
                        int256_t invertedtotal = total  -1;
                        stdcout  The result in decimal form is
                                   n;
                        stdcout  invertedtotal;
                    }
                    else
                    {
                        stdcout  n
                                   n;
                        stdcout  The result in decimal form is
                                   n;
                        stdcout  total;
                    }

                    printf(n);
                    printf(n);
                    gettimeofday(&end, NULL);
                    get_time = (end.tv_sec - start.tv_sec) + (end.tv_usec - start.tv_usec)  1.0E-6;
                    printf(Computation Time %lf[sec]n, get_time);
                    printf(n);

                    printf(I hope you remembered what calculation you performed!n);
                }

                first carry =1 - 2nd carry = 0 - third carry = 1 (14 11 10) check
                else if (firstcarry == 1 && secondCarryOne == 0 && thirdCarryZero == 1)
                {
                    printf(14 + 11 + 10 + 1 n);
                    stdstring binary_combined6 = binary14 + binary11 + binary10 + binary1;
                    if (int_negative != 4)
                    {
                        stdcout  The result in binary form is
                                   n;
                        stdcout  binary_combined6  n;
                    }
                    Length is the number of bits
                    int length = binary_combined6.length();
                    int lastchar = length - 1;
                    char char_array[length + 1];
                    strcpy(char_array, binary_combined6.c_str());

                    int256_t total = 0;

                    Positive or negative checking
                    char charvalue = char_array[0];
                    stdstring stringvalue(1, charvalue);
                    int intvalue = stdstoi(stringvalue);

                    if the result is either double positive or negative
                    if (int_negative == 0  int_negative == 4)
                    {
                        for (int i = 0; i  length; ++i)
                        {
                            char charvalue = char_array[i];
                            stdstring stringvalue(1, charvalue);
                            int intvalue = stdstoi(stringvalue);
                            int256_t calc1 = total  2;
                            total = calc1 + intvalue;
                        }
                    }
                    if either one of the values are negative, check for negative
                    else if (int_negative != 4)
                    {
                        if (length == 128 && intvalue == 1)
                        {
                            int256_t lastcharvalue = 2;

                            for (int i = 0; i  lastchar - 1; ++i)
                            {
                                lastcharvalue = lastcharvalue  2;
                            }

                            for (int i = 1; i  length; ++i)
                            {
                                char charvalue = char_array[i];
                                stdstring stringvalue(1, charvalue);
                                int intvalue = stdstoi(stringvalue);
                                int256_t calc1 = total  2;
                                total = calc1 + intvalue;
                            }
                            total = total - lastcharvalue;
                        }
                        else
                        {
                            for (int i = 0; i  length; ++i)
                            {
                                char charvalue = char_array[i];
                                stdstring stringvalue(1, charvalue);
                                int intvalue = stdstoi(stringvalue);
                                int256_t calc1 = total  2;
                                total = calc1 + intvalue;
                            }
                        }
                    }
                    if (int_negative == 4)
                    {
                        int256_t invertedtotal = total  -1;
                        stdcout  The result in decimal form is
                                   n;
                        stdcout  invertedtotal;
                    }
                    else
                    {
                        stdcout  n
                                   n;
                        stdcout  The result in decimal form is
                                   n;
                        stdcout  total;
                    }

                    printf(n);
                    printf(n);
                    gettimeofday(&end, NULL);
                    get_time = (end.tv_sec - start.tv_sec) + (end.tv_usec - start.tv_usec)  1.0E-6;
                    printf(Computation Time %lf[sec]n, get_time);
                    printf(n);

                    printf(I hope you remembered what calculation you performed!n);
                }

                first carry =1 - 2nd carry = 1 - third carry = 0 (13 12 10) check
                else if (firstcarry == 1 && secondCarryOne == 1 && thirdCarryOne == 0)
                {
                    printf(13 + 12 + 10 + 1 n);
                    stdstring binary_combined7 = binary13 + binary12 + binary10 + binary1;
                    if (int_negative != 4)
                    {
                        stdcout  The result in binary form is
                                   n;
                        stdcout  binary_combined7  n;
                    }
                    Length is the number of bits
                    int length = binary_combined7.length();
                    int lastchar = length - 1;
                    char char_array[length + 1];
                    strcpy(char_array, binary_combined7.c_str());

                    int256_t total = 0;

                    Positive or negative checking
                    char charvalue = char_array[0];
                    stdstring stringvalue(1, charvalue);
                    int intvalue = stdstoi(stringvalue);

                    if the result is either double positive or negative
                    if (int_negative == 0  int_negative == 4)
                    {
                        for (int i = 0; i  length; ++i)
                        {
                            char charvalue = char_array[i];
                            stdstring stringvalue(1, charvalue);
                            int intvalue = stdstoi(stringvalue);
                            int256_t calc1 = total  2;
                            total = calc1 + intvalue;
                        }
                    }
                    if either one of the values are negative, check for negative
                    else if (int_negative != 4)
                    {
                        if (length == 128 && intvalue == 1)
                        {
                            int256_t lastcharvalue = 2;

                            for (int i = 0; i  lastchar - 1; ++i)
                            {
                                lastcharvalue = lastcharvalue  2;
                            }

                            for (int i = 1; i  length; ++i)
                            {
                                char charvalue = char_array[i];
                                stdstring stringvalue(1, charvalue);
                                int intvalue = stdstoi(stringvalue);
                                int256_t calc1 = total  2;
                                total = calc1 + intvalue;
                            }
                            total = total - lastcharvalue;
                        }
                        else
                        {
                            for (int i = 0; i  length; ++i)
                            {
                                char charvalue = char_array[i];
                                stdstring stringvalue(1, charvalue);
                                int intvalue = stdstoi(stringvalue);
                                int256_t calc1 = total  2;
                                total = calc1 + intvalue;
                            }
                        }
                    }
                    if (int_negative == 4)
                    {
                        int256_t invertedtotal = total  -1;
                        stdcout  The result in decimal form is
                                   n;
                        stdcout  invertedtotal;
                    }
                    else
                    {
                        stdcout  n
                                   n;
                        stdcout  The result in decimal form is
                                   n;
                        stdcout  total;
                    }

                    printf(n);
                    printf(n);
                    gettimeofday(&end, NULL);
                    get_time = (end.tv_sec - start.tv_sec) + (end.tv_usec - start.tv_usec)  1.0E-6;
                    printf(Computation Time %lf[sec]n, get_time);
                    printf(n);

                    printf(I hope you remembered what calculation you performed!n);
                }

                first carry =0 - 2nd carry = 1 - third carry = 0 (13 12 9) check
                else if (firstcarry == 0 && secondCarryZero == 1 && thirdCarryOne == 0)
                {
                    printf(13 + 12 + 9 + 1 n);
                    stdstring binary_combined8 = binary13 + binary12 + binary9 + binary1;
                    if (int_negative != 4)
                    {
                        stdcout  The result in binary form is
                                   n;
                        stdcout  binary_combined8  n;
                    }
                    Length is the number of bits
                    int length = binary_combined8.length();
                    int lastchar = length - 1;
                    char char_array[length + 1];
                    strcpy(char_array, binary_combined8.c_str());

                    int256_t total = 0;

                    Positive or negative checking
                    char charvalue = char_array[0];
                    stdstring stringvalue(1, charvalue);
                    int intvalue = stdstoi(stringvalue);

                    if the result is either double positive or negative
                    if (int_negative == 0  int_negative == 4)
                    {
                        for (int i = 0; i  length; ++i)
                        {
                            char charvalue = char_array[i];
                            stdstring stringvalue(1, charvalue);
                            int intvalue = stdstoi(stringvalue);
                            int256_t calc1 = total  2;
                            total = calc1 + intvalue;
                        }
                    }
                    if either one of the values are negative, check for negative
                    else if (int_negative != 4)
                    {
                        if (length == 128 && intvalue == 1)
                        {
                            int256_t lastcharvalue = 2;

                            for (int i = 0; i  lastchar - 1; ++i)
                            {
                                lastcharvalue = lastcharvalue  2;
                            }

                            for (int i = 1; i  length; ++i)
                            {
                                char charvalue = char_array[i];
                                stdstring stringvalue(1, charvalue);
                                int intvalue = stdstoi(stringvalue);
                                int256_t calc1 = total  2;
                                total = calc1 + intvalue;
                            }
                            total = total - lastcharvalue;
                        }
                        else
                        {
                            for (int i = 0; i  length; ++i)
                            {
                                char charvalue = char_array[i];
                                stdstring stringvalue(1, charvalue);
                                int intvalue = stdstoi(stringvalue);
                                int256_t calc1 = total  2;
                                total = calc1 + intvalue;
                            }
                        }
                    }
                    if (int_negative == 4)
                    {
                        int256_t invertedtotal = total  -1;
                        stdcout  The result in decimal form is
                                   n;
                        stdcout  invertedtotal;
                    }
                    else
                    {
                        stdcout  n
                                   n;
                        stdcout  The result in decimal form is
                                   n;
                        stdcout  total;
                    }

                    printf(n);
                    printf(n);
                    gettimeofday(&end, NULL);
                    get_time = (end.tv_sec - start.tv_sec) + (end.tv_usec - start.tv_usec)  1.0E-6;
                    printf(Computation Time %lf[sec]n, get_time);
                    printf(n);

                    printf(I hope you remembered what calculation you performed!n);
                }

                else
                {
                    printf(Something went wrong... n);
                }
            }
        }

        else if (int_op == 2)
        {
            stdcout  Result for   int_bit   bit Subtraction computation
                       n
                       n;

            if (int_bit == 32)
            {
                LweSample result = new_gate_bootstrapping_ciphertext_array(32, params);
                for (int i = 0; i  32; i++)
                    import_gate_bootstrapping_ciphertext_fromFile(answer_data, &result[i], params);
                fclose(answer_data);

                 decrypt and rebuild the answer
                int32_t int_answer1 = 0;
                for (int i = 0; i  32; i++)
                {
                    int ai = bootsSymDecrypt(&result[i], key)  0;
                    int_answer1 = (ai  i);
                }

                stdstring binary1 = stdbitset32(int_answer1).to_string();
                stdstring binary_combined = binary1;

                if (int_negative != 1)
                {
                    stdcout  The result in binary form is
                               n;
                    stdcout  binary_combined;
                }

                Length is the number of bits
                int length = binary_combined.length();
                int lastchar = length - 1;
                char char_array[length + 1];
                strcpy(char_array, binary_combined.c_str());
                int256_t total = 0;

                Positive or negative checking
                char charvalue = char_array[0];
                stdstring stringvalue(1, charvalue);
                int intvalue = stdstoi(stringvalue);

                if (int_negative == 2)
                {

                    for (int i = 0; i  length; ++i)
                    {
                        char charvalue = char_array[i];
                        stdstring stringvalue(1, charvalue);
                        int intvalue = stdstoi(stringvalue);
                        int256_t calc1 = total  2;
                        total = calc1 + intvalue;
                    }
                }
                else
                {
                     int_negative is not equals to 2
                    if (length == 32 && intvalue == 1)
                    {
                        int256_t lastcharvalue = 2;

                        for (int i = 0; i  lastchar - 1; ++i)
                        {
                            lastcharvalue = lastcharvalue  2;
                        }

                        for (int i = 1; i  length; ++i)
                        {
                            char charvalue = char_array[i];
                            stdstring stringvalue(1, charvalue);
                            int intvalue = stdstoi(stringvalue);
                            int256_t calc1 = total  2;
                            total = calc1 + intvalue;
                        }
                        total = total - lastcharvalue;
                    }

                    else
                    {
                        for (int i = 0; i  length; ++i)
                        {
                            char charvalue = char_array[i];
                            stdstring stringvalue(1, charvalue);
                            int intvalue = stdstoi(stringvalue);
                            int256_t calc1 = total  2;
                            total = calc1 + intvalue;
                        }
                    }
                }

                if (int_negative == 1)
                {
                    int256_t invertedtotal = total  -1;
                    stdcout  The result in decimal form is
                               n;
                    stdcout  invertedtotal;
                }
                else
                {
                    stdcout  n
                               n;
                    stdcout  The result in decimal form is
                               n;
                    stdcout  total;
                }
                printf(n);
                printf(n);
                gettimeofday(&end, NULL);
                get_time = (end.tv_sec - start.tv_sec) + (end.tv_usec - start.tv_usec)  1.0E-6;
                printf(Computation Time %lf[sec]n, get_time);
                printf(n);

                printf(I hope you remembered what calculation you performed!n);

                 clean up all pointers
                delete_gate_bootstrapping_ciphertext_array(32, result);
                delete_gate_bootstrapping_secret_keyset(key);
                delete_gate_bootstrapping_secret_keyset(nbitkey);
            }
            else if (int_bit == 64)
            {
                printf(Verifying 64 bit Subtraction on MPI n);
                 read the 32 ciphertexts of the result
                LweSample result = new_gate_bootstrapping_ciphertext_array(32, params);
                LweSample result2 = new_gate_bootstrapping_ciphertext_array(32, params);
                LweSample result3 = new_gate_bootstrapping_ciphertext_array(32, params);
                LweSample carry1 = new_gate_bootstrapping_ciphertext_array(32, params);

                 export the 64 ciphertexts to a file (for the cloud)
                for (int i = 0; i  32; i++)
                    import_gate_bootstrapping_ciphertext_fromFile(answer_data, &result[i], params);
                for (int i = 0; i  32; i++)
                    import_gate_bootstrapping_ciphertext_fromFile(answer_data, &carry1[i], params);
                for (int i = 0; i  32; i++)
                    import_gate_bootstrapping_ciphertext_fromFile(answer_data, &result2[i], params);
                for (int i = 0; i  32; i++)
                    import_gate_bootstrapping_ciphertext_fromFile(answer_data, &result3[i], params);
                fclose(answer_data);
                printf(Exportingn);

                 decrypt and rebuild the answer
                int32_t firstResult = 0;
                for (int i = 0; i  32; i++)
                {
                    int ai = bootsSymDecrypt(&result[i], key)  0;
                    firstResult = (ai  i);
                }

                int32_t VerifyCarry = 0;
                for (int i = 0; i  32; i++)
                {
                    int ai = bootsSymDecrypt(&carry1[i], key)  0;
                    VerifyCarry = (ai  i);
                }

                int32_t ifCarryisZero = 0;
                for (int i = 0; i  32; i++)
                {
                    int ai = bootsSymDecrypt(&result2[i], key)  0;
                    ifCarryisZero = (ai  i);
                }

                int32_t ifCarryisOne = 0;
                for (int i = 0; i  32; i++)
                {
                    int ai = bootsSymDecrypt(&result3[i], key)  0;
                    ifCarryisOne = (ai  i);
                }

                last 32
                when testing == ciphertext1
                stdstring binary1 = stdbitset32(firstResult).to_string();
                carry1 in binary
                stdstring binary4 = stdbitset32(VerifyCarry).to_string();
                32 if carry1 = 0
                stdstring binary2 = stdbitset32(ifCarryisZero).to_string();
                32 if carry1 = 1
                stdstring binary3 = stdbitset32(ifCarryisOne).to_string();

                stdcout  First 32 (Binary1) ;
                stdcout  binary1;
                stdcout  n
                           n;
                stdcout  Carry1 (Carry1 in binary) ;
                stdcout  binary4;
                stdcout  n
                           n;
                stdcout  If carry is zero (Binary2) ;
                stdcout  binary2;
                stdcout  n
                           n;
                stdcout  If carry is one (Binary3) ;
                stdcout  binary3;
                stdcout  n
                           n;

                printf(Carry1 is %dn, VerifyCarry);

                if (VerifyCarry == 0)
                {
                    printf(Binary2 + Binary1 %dn, VerifyCarry);
                    stdstring binary_combined = binary2 + binary1;

                    if (int_negative != 1)
                    {
                        stdcout  The result in binary form is
                                   n;
                        stdcout  binary_combined;
                    }

                    int length = binary_combined.length();
                    int lastchar = length - 1;
                    char char_array[length + 1];
                    strcpy(char_array, binary_combined.c_str());

                    int256_t total = 0;

                    Positive or negative checking
                    char charvalue = char_array[0];
                    stdstring stringvalue(1, charvalue);
                    int intvalue = stdstoi(stringvalue);

                    if (int_negative == 2)
                    {
                        for (int i = 0; i  length; ++i)
                        {
                            char charvalue = char_array[i];
                            stdstring stringvalue(1, charvalue);
                            int intvalue = stdstoi(stringvalue);
                            int256_t calc1 = total  2;
                            total = calc1 + intvalue;
                        }
                    }
                    else
                    {
                        if (length == 64 && intvalue == 1)
                        {
                            int256_t lastcharvalue = 2;

                            for (int i = 0; i  lastchar - 1; ++i)
                            {
                                lastcharvalue = lastcharvalue  2;
                            }

                            for (int i = 1; i  length; ++i)
                            {
                                char charvalue = char_array[i];
                                stdstring stringvalue(1, charvalue);
                                int intvalue = stdstoi(stringvalue);
                                int256_t calc1 = total  2;
                                total = calc1 + intvalue;
                            }
                            total = total - lastcharvalue;
                        }
                        else
                        {
                            for (int i = 0; i  length; ++i)
                            {
                                char charvalue = char_array[i];
                                stdstring stringvalue(1, charvalue);
                                int intvalue = stdstoi(stringvalue);
                                int256_t calc1 = total  2;
                                total = calc1 + intvalue;
                            }
                        }
                    }
                    if (int_negative == 1)
                    {
                        int256_t invertedtotal = total  -1;
                        stdcout  The result in decimal form is
                                   n;
                        stdcout  invertedtotal;
                    }
                    else
                    {
                        stdcout  n
                                   n;
                        stdcout  The result in decimal form is
                                   n;
                        stdcout  total;
                    }
                }
                if carry = 1
                else if (VerifyCarry == 1)
                {
                    printf(Binary3 + Binary1 %dn, VerifyCarry);
                    stdstring binary_combined2 = binary3 + binary1;

                    if (int_negative != 1)
                    {
                        stdcout  The result in binary form is
                                   n;
                        stdcout  binary_combined2;
                    }

                    int length = binary_combined2.length();
                    int lastchar = length - 1;
                    char char_array[length + 1];
                    strcpy(char_array, binary_combined2.c_str());

                    int256_t total = 0;

                    Positive or negative checking
                    char charvalue = char_array[0];
                    stdstring stringvalue(1, charvalue);
                    int intvalue = stdstoi(stringvalue);

                    if (int_negative == 2)
                    {
                        for (int i = 0; i  length; ++i)
                        {
                            char charvalue = char_array[i];
                            stdstring stringvalue(1, charvalue);
                            int intvalue = stdstoi(stringvalue);
                            int256_t calc1 = total  2;
                            total = calc1 + intvalue;
                        }
                    }
                    else
                    {
                        if (length == 64 && intvalue == 1)
                        {
                            int256_t lastcharvalue = 2;

                            for (int i = 0; i  lastchar - 1; ++i)
                            {
                                lastcharvalue = lastcharvalue  2;
                            }

                            for (int i = 1; i  length; ++i)
                            {
                                char charvalue = char_array[i];
                                stdstring stringvalue(1, charvalue);
                                int intvalue = stdstoi(stringvalue);
                                int256_t calc1 = total  2;
                                total = calc1 + intvalue;
                            }
                            total = total - lastcharvalue;
                        }
                        else
                        {
                            for (int i = 0; i  length; ++i)
                            {
                                char charvalue = char_array[i];
                                stdstring stringvalue(1, charvalue);
                                int intvalue = stdstoi(stringvalue);
                                int256_t calc1 = total  2;
                                total = calc1 + intvalue;
                            }
                        }
                    }
                    if (int_negative == 1)
                    {
                        int256_t invertedtotal = total  -1;
                        stdcout  The result in decimal form is
                                   n;
                        stdcout  invertedtotal;
                    }
                    else
                    {
                        stdcout  n
                                   n;
                        stdcout  The result in decimal form is
                                   n;
                        stdcout  total;
                    }
                }

                else
                {
                    printf(Failn);
                }

                printf(n);
                printf(n);
                printf(n);
                printf(I hope you remembered what calculation you performed!n);

                delete_gate_bootstrapping_ciphertext_array(32, result);
                delete_gate_bootstrapping_ciphertext_array(32, result2);
                delete_gate_bootstrapping_ciphertext_array(32, result3);
                delete_gate_bootstrapping_ciphertext_array(32, carry1);
                delete_gate_bootstrapping_secret_keyset(key);
                delete_gate_bootstrapping_secret_keyset(nbitkey);
            }
            else if (int_bit == 128)
            {
                printf(Verifying 128 bit subtraction on MPI n);
                 read the 32 ciphertexts of the result
                LweSample result = new_gate_bootstrapping_ciphertext_array(32, params);
                LweSample result2 = new_gate_bootstrapping_ciphertext_array(32, params);
                LweSample result3 = new_gate_bootstrapping_ciphertext_array(32, params);
                LweSample result4 = new_gate_bootstrapping_ciphertext_array(32, params);
                LweSample result5 = new_gate_bootstrapping_ciphertext_array(32, params);
                LweSample result6 = new_gate_bootstrapping_ciphertext_array(32, params);
                LweSample result7 = new_gate_bootstrapping_ciphertext_array(32, params);
                LweSample carry1 = new_gate_bootstrapping_ciphertext_array(32, params);
                LweSample carry2 = new_gate_bootstrapping_ciphertext_array(32, params);
                LweSample carry3 = new_gate_bootstrapping_ciphertext_array(32, params);
                LweSample carry4 = new_gate_bootstrapping_ciphertext_array(32, params);
                LweSample carry5 = new_gate_bootstrapping_ciphertext_array(32, params);
                LweSample carry6 = new_gate_bootstrapping_ciphertext_array(32, params);
                LweSample carry7 = new_gate_bootstrapping_ciphertext_array(32, params);

                 export the 64 ciphertexts to a file (for the cloud)
                for (int i = 0; i  32; i++)
                    import_gate_bootstrapping_ciphertext_fromFile(answer_data, &result[i], params);
                for (int i = 0; i  32; i++)
                    import_gate_bootstrapping_ciphertext_fromFile(answer_data, &carry1[i], params);
                for (int i = 0; i  32; i++)
                    import_gate_bootstrapping_ciphertext_fromFile(answer_data, &carry2[i], params);
                for (int i = 0; i  32; i++)
                    import_gate_bootstrapping_ciphertext_fromFile(answer_data, &carry3[i], params);
                for (int i = 0; i  32; i++)
                    import_gate_bootstrapping_ciphertext_fromFile(answer_data, &carry4[i], params);
                for (int i = 0; i  32; i++)
                    import_gate_bootstrapping_ciphertext_fromFile(answer_data, &carry5[i], params);
                for (int i = 0; i  32; i++)
                    import_gate_bootstrapping_ciphertext_fromFile(answer_data, &carry6[i], params);
                for (int i = 0; i  32; i++)
                    import_gate_bootstrapping_ciphertext_fromFile(answer_data, &carry7[i], params);
                for (int i = 0; i  32; i++)
                    import_gate_bootstrapping_ciphertext_fromFile(answer_data, &result2[i], params);
                for (int i = 0; i  32; i++)
                    import_gate_bootstrapping_ciphertext_fromFile(answer_data, &result3[i], params);
                for (int i = 0; i  32; i++)
                    import_gate_bootstrapping_ciphertext_fromFile(answer_data, &result4[i], params);
                for (int i = 0; i  32; i++)
                    import_gate_bootstrapping_ciphertext_fromFile(answer_data, &result5[i], params);
                for (int i = 0; i  32; i++)
                    import_gate_bootstrapping_ciphertext_fromFile(answer_data, &result6[i], params);
                for (int i = 0; i  32; i++)
                    import_gate_bootstrapping_ciphertext_fromFile(answer_data, &result7[i], params);
                fclose(answer_data);
                printf(Exportingn);

                 decrypt and rebuild the answer
                int32_t firstResult = 0;
                for (int i = 0; i  32; i++)
                {
                    int ai = bootsSymDecrypt(&result[i], key)  0;
                    firstResult = (ai  i);
                }
                int32_t firstcarry = 0;
                for (int i = 0; i  32; i++)
                {
                    int ai = bootsSymDecrypt(&carry1[i], key)  0;
                    firstcarry = (ai  i);
                }
                int32_t secondCarryZero = 0;
                for (int i = 0; i  32; i++)
                {
                    int ai = bootsSymDecrypt(&carry2[i], key)  0;
                    secondCarryZero = (ai  i);
                }
                int32_t secondCarryOne = 0;
                for (int i = 0; i  32; i++)
                {
                    int ai = bootsSymDecrypt(&carry3[i], key)  0;
                    secondCarryOne = (ai  i);
                }

                int32_t thirdCarryZero = 0;
                for (int i = 0; i  32; i++)
                {
                    int ai = bootsSymDecrypt(&carry4[i], key)  0;
                    thirdCarryZero = (ai  i);
                }
                int32_t thirdCarryOne = 0;
                for (int i = 0; i  32; i++)
                {
                    int ai = bootsSymDecrypt(&carry5[i], key)  0;
                    thirdCarryOne = (ai  i);
                }
                int32_t fourthCarryZero = 0;
                for (int i = 0; i  32; i++)
                {
                    int ai = bootsSymDecrypt(&carry6[i], key)  0;
                    fourthCarryZero = (ai  i);
                }
                int32_t fourthCarryOne = 0;
                for (int i = 0; i  32; i++)
                {
                    int ai = bootsSymDecrypt(&carry7[i], key)  0;
                    fourthCarryOne = (ai  i);
                }
                
                int32_t secondResultZero = 0;
                for (int i = 0; i  32; i++)
                {
                    int ai = bootsSymDecrypt(&result2[i], key)  0;
                    secondResultZero = (ai  i);
                }
                int32_t secondResultOne = 0;
                for (int i = 0; i  32; i++)
                {
                    int ai = bootsSymDecrypt(&result3[i], key)  0;
                    secondResultOne = (ai  i);
                }
                int32_t thirdResultZero = 0;
                for (int i = 0; i  32; i++)
                {
                    int ai = bootsSymDecrypt(&result4[i], key)  0;
                    thirdResultZero = (ai  i);
                }
                int32_t thirdResultOne = 0;
                for (int i = 0; i  32; i++)
                {
                    int ai = bootsSymDecrypt(&result5[i], key)  0;
                    thirdResultOne = (ai  i);
                }
                int32_t fourthResultZero = 0;
                for (int i = 0; i  32; i++)
                {
                    int ai = bootsSymDecrypt(&result6[i], key)  0;
                    fourthResultZero = (ai  i);
                }
                int32_t fourthResultOne = 0;
                for (int i = 0; i  32; i++)
                {
                    int ai = bootsSymDecrypt(&result7[i], key)  0;
                    fourthResultOne = (ai  i);
                }

                printf(carry1  %dn, firstcarry);
                printf(carry2 secondCarryZero  %dn, secondCarryZero);
                printf(carry3 secondCarryOne  %dn, secondCarryOne);
                printf(carry4 thirdCarryZero  %dn, thirdCarryZero);
                printf(carry5 thirdCarryOne  %dn, thirdCarryOne);
                printf(carry6 fourthCarryZero  %dn, fourthCarryZero);
                printf(carry7 fourthCarryOne  %dn, fourthCarryOne);

                last 32 bits
                stdstring binary1 = stdbitset32(firstResult).to_string();

                carry of ciphertext1 and 9
                stdstring binary2 = stdbitset32(firstcarry).to_string();

                carry=0
                carry of ciphertext2 and 10
                stdstring binary3 = stdbitset32(secondCarryZero).to_string();

                carry=1
                carry of ciphertext2 and 10
                stdstring binary4 = stdbitset32(secondCarryOne).to_string();

                carry=0
                carry of ciphertext3 and 11
                stdstring binary5 = stdbitset32(thirdCarryZero).to_string();

                carry=1
                carry of ciphertext3 and 11
                stdstring binary6 = stdbitset32(thirdCarryOne).to_string();

                carry=0
                carry of ciphertext4 and 12
                stdstring binary7 = stdbitset32(fourthCarryZero).to_string();

                carry=1
                carry of ciphertext4 and 12
                stdstring binary8 = stdbitset32(fourthCarryOne).to_string();

                result of next 32 if carry is 0        2
                stdstring binary9 = stdbitset32(secondResultZero).to_string();
                result of next 32 if carry is 1        3
                stdstring binary10 = stdbitset32(secondResultOne).to_string();
                64 - 96 if carry is 0                  4
                stdstring binary11 = stdbitset32(thirdResultZero).to_string();
                64 - 96 if carry is 1                  5
                stdstring binary12 = stdbitset32(thirdResultOne).to_string();
                96 -128 if carry is 0                  6
                stdstring binary13 = stdbitset32(fourthResultZero).to_string();
                96 -128 carry is 1                     7
                stdstring binary14 = stdbitset32(fourthResultOne).to_string();


                printf(Verifying.............................. n);

                first carry =0 - 2nd carry = 0 - third carry = 0 (13 11 9) check
                if (firstcarry == 0 && secondCarryZero == 0 && thirdCarryZero == 0)
                {
                    printf(13 + 11 + 9 + 1 n);
                    stdstring binary_combined1 = binary13 + binary11 + binary9 + binary1;
                    if (int_negative != 1)
                    {
                        stdcout  The result in binary form is
                                   n;
                        stdcout  binary_combined1;
                    }

                    Length is the number of bits
                    int length = binary_combined1.length();
                    int lastchar = length - 1;
                    char char_array[length + 1];
                    strcpy(char_array, binary_combined1.c_str());

                    int256_t total = 0;

                    Positive or negative checking
                    char charvalue = char_array[0];
                    stdstring stringvalue(1, charvalue);
                    int intvalue = stdstoi(stringvalue);

                    if (int_negative == 2)
                    {
                        for (int i = 0; i  length; ++i)
                        {
                            char charvalue = char_array[i];
                            stdstring stringvalue(1, charvalue);
                            int intvalue = stdstoi(stringvalue);
                            int256_t calc1 = total  2;
                            total = calc1 + intvalue;
                        }
                    }
                    else
                    {
                        if (length == 128 && intvalue == 1)
                        {
                            int256_t lastcharvalue = 2;

                            for (int i = 0; i  lastchar - 1; ++i)
                            {
                                lastcharvalue = lastcharvalue  2;
                            }

                            for (int i = 1; i  length; ++i)
                            {
                                char charvalue = char_array[i];
                                stdstring stringvalue(1, charvalue);
                                int intvalue = stdstoi(stringvalue);
                                int256_t calc1 = total  2;
                                total = calc1 + intvalue;
                            }

                            total = total - lastcharvalue;
                        }
                        else
                        {
                            for (int i = 0; i  length; ++i)
                            {
                                char charvalue = char_array[i];
                                stdstring stringvalue(1, charvalue);
                                int intvalue = stdstoi(stringvalue);
                                int256_t calc1 = total  2;
                                total = calc1 + intvalue;
                            }
                        }
                    }
                    if (int_negative == 1)
                    {
                        int256_t invertedtotal = total  -1;
                        stdcout  The result in decimal form is
                                   n;
                        stdcout  invertedtotal;
                    }
                    else
                    {
                        stdcout  n
                                   n;
                        stdcout  The result in decimal form is
                                   n;
                        stdcout  total;
                    }

                    printf(n);
                    printf(n);
                    gettimeofday(&end, NULL);
                    get_time = (end.tv_sec - start.tv_sec) + (end.tv_usec - start.tv_usec)  1.0E-6;
                    printf(Computation Time %lf[sec]n, get_time);
                    printf(n);

                    printf(I hope you remembered what calculation you performed!n);
                }

                first carry =0 - 2nd carry = 0 - third carry = 1 (14 11 9) check
                else if (firstcarry == 0 && secondCarryZero == 0 && thirdCarryZero == 1)
                {
                    printf(14 + 11 + 9 + 1 n);
                    stdstring binary_combined2 = binary14 + binary11 + binary9 + binary1;
                    if (int_negative != 1)
                    {
                        stdcout  The result in binary form is
                                   n;
                        stdcout  binary_combined2;
                    }

                    Length is the number of bits
                    int length = binary_combined2.length();
                    int lastchar = length - 1;
                    char char_array[length + 1];
                    strcpy(char_array, binary_combined2.c_str());

                    int256_t total = 0;

                    Positive or negative checking
                    char charvalue = char_array[0];
                    stdstring stringvalue(1, charvalue);
                    int intvalue = stdstoi(stringvalue);

                    if (int_negative == 2)
                    {
                        for (int i = 0; i  length; ++i)
                        {
                            char charvalue = char_array[i];
                            stdstring stringvalue(1, charvalue);
                            int intvalue = stdstoi(stringvalue);
                            int256_t calc1 = total  2;
                            total = calc1 + intvalue;
                        }
                    }
                    else
                    {
                        if (length == 128 && intvalue == 1)
                        {
                            int256_t lastcharvalue = 2;

                            for (int i = 0; i  lastchar - 1; ++i)
                            {
                                lastcharvalue = lastcharvalue  2;
                            }

                            for (int i = 1; i  length; ++i)
                            {
                                char charvalue = char_array[i];
                                stdstring stringvalue(1, charvalue);
                                int intvalue = stdstoi(stringvalue);
                                int256_t calc1 = total  2;
                                total = calc1 + intvalue;
                            }

                            total = total - lastcharvalue;
                        }
                        else
                        {
                            for (int i = 0; i  length; ++i)
                            {
                                char charvalue = char_array[i];
                                stdstring stringvalue(1, charvalue);
                                int intvalue = stdstoi(stringvalue);
                                int256_t calc1 = total  2;
                                total = calc1 + intvalue;
                            }
                        }
                    }
                    if (int_negative == 1)
                    {
                        int256_t invertedtotal = total  -1;
                        stdcout  The result in decimal form is
                                   n;
                        stdcout  invertedtotal;
                    }
                    else
                    {
                        stdcout  n
                                   n;
                        stdcout  The result in decimal form is
                                   n;
                        stdcout  total;
                    }

                    printf(n);
                    printf(n);
                    gettimeofday(&end, NULL);
                    get_time = (end.tv_sec - start.tv_sec) + (end.tv_usec - start.tv_usec)  1.0E-6;
                    printf(Computation Time %lf[sec]n, get_time);
                    printf(n);

                    printf(I hope you remembered what calculation you performed!n);
                }

                first carry =0 - 2nd carry = 1 - third carry = 1 (14 12 9) check
                else if (firstcarry == 0 && secondCarryZero == 1 && thirdCarryOne == 1)
                {
                    printf(14 + 12 + 9 + 1 n);
                    stdstring binary_combined3 = binary14 + binary12 + binary9 + binary1;
                    if (int_negative != 1)
                    {
                        stdcout  The result in binary form is
                                   n;
                        stdcout  binary_combined3;
                    }

                    Length is the number of bits
                    int length = binary_combined3.length();
                    int lastchar = length - 1;
                    char char_array[length + 1];
                    strcpy(char_array, binary_combined3.c_str());

                    int256_t total = 0;

                    Positive or negative checking
                    char charvalue = char_array[0];
                    stdstring stringvalue(1, charvalue);
                    int intvalue = stdstoi(stringvalue);

                    if (int_negative == 2)
                    {
                        for (int i = 0; i  length; ++i)
                        {
                            char charvalue = char_array[i];
                            stdstring stringvalue(1, charvalue);
                            int intvalue = stdstoi(stringvalue);
                            int256_t calc1 = total  2;
                            total = calc1 + intvalue;
                        }
                    }
                    else
                    {
                        if (length == 128 && intvalue == 1)
                        {
                            int256_t lastcharvalue = 2;

                            for (int i = 0; i  lastchar - 1; ++i)
                            {
                                lastcharvalue = lastcharvalue  2;
                            }

                            for (int i = 1; i  length; ++i)
                            {
                                char charvalue = char_array[i];
                                stdstring stringvalue(1, charvalue);
                                int intvalue = stdstoi(stringvalue);
                                int256_t calc1 = total  2;
                                total = calc1 + intvalue;
                            }

                            total = total - lastcharvalue;
                        }
                        else
                        {
                            for (int i = 0; i  length; ++i)
                            {
                                char charvalue = char_array[i];
                                stdstring stringvalue(1, charvalue);
                                int intvalue = stdstoi(stringvalue);
                                int256_t calc1 = total  2;
                                total = calc1 + intvalue;
                            }
                        }
                    }
                    if (int_negative == 1)
                    {
                        int256_t invertedtotal = total  -1;
                        stdcout  The result in decimal form is
                                   n;
                        stdcout  invertedtotal;
                    }
                    else
                    {
                        stdcout  n
                                   n;
                        stdcout  The result in decimal form is
                                   n;
                        stdcout  total;
                    }

                    printf(n);
                    printf(n);
                    gettimeofday(&end, NULL);
                    get_time = (end.tv_sec - start.tv_sec) + (end.tv_usec - start.tv_usec)  1.0E-6;
                    printf(Computation Time %lf[sec]n, get_time);
                    printf(n);

                    printf(I hope you remembered what calculation you performed!n);
                }

                first carry =1 - 2nd carry = 1 - third carry = 1 (14 12 10) check
                else if (firstcarry == 1 && secondCarryOne == 1 && thirdCarryOne == 1)
                {
                    printf(14 + 12 + 10 + 1 n);
                    stdstring binary_combined4 = binary14 + binary12 + binary10 + binary1;
                    if (int_negative != 1)
                    {
                        stdcout  The result in binary form is
                                   n;
                        stdcout  binary_combined4;
                    }

                    Length is the number of bits
                    int length = binary_combined4.length();
                    int lastchar = length - 1;
                    char char_array[length + 1];
                    strcpy(char_array, binary_combined4.c_str());

                    int256_t total = 0;

                    Positive or negative checking
                    char charvalue = char_array[0];
                    stdstring stringvalue(1, charvalue);
                    int intvalue = stdstoi(stringvalue);

                    if (int_negative == 2)
                    {
                        for (int i = 0; i  length; ++i)
                        {
                            char charvalue = char_array[i];
                            stdstring stringvalue(1, charvalue);
                            int intvalue = stdstoi(stringvalue);
                            int256_t calc1 = total  2;
                            total = calc1 + intvalue;
                        }
                    }
                    else
                    {
                        if (length == 128 && intvalue == 1)
                        {
                            int256_t lastcharvalue = 2;

                            for (int i = 0; i  lastchar - 1; ++i)
                            {
                                lastcharvalue = lastcharvalue  2;
                            }

                            for (int i = 1; i  length; ++i)
                            {
                                char charvalue = char_array[i];
                                stdstring stringvalue(1, charvalue);
                                int intvalue = stdstoi(stringvalue);
                                int256_t calc1 = total  2;
                                total = calc1 + intvalue;
                            }

                            total = total - lastcharvalue;
                        }
                        else
                        {
                            for (int i = 0; i  length; ++i)
                            {
                                char charvalue = char_array[i];
                                stdstring stringvalue(1, charvalue);
                                int intvalue = stdstoi(stringvalue);
                                int256_t calc1 = total  2;
                                total = calc1 + intvalue;
                            }
                        }
                    }
                    if (int_negative == 1)
                    {
                        int256_t invertedtotal = total  -1;
                        stdcout  The result in decimal form is
                                   n;
                        stdcout  invertedtotal;
                    }
                    else
                    {
                        stdcout  n
                                   n;
                        stdcout  The result in decimal form is
                                   n;
                        stdcout  total;
                    }

                    printf(n);
                    printf(n);
                    gettimeofday(&end, NULL);
                    get_time = (end.tv_sec - start.tv_sec) + (end.tv_usec - start.tv_usec)  1.0E-6;
                    printf(Computation Time %lf[sec]n, get_time);
                    printf(n);

                    printf(I hope you remembered what calculation you performed!n);
                }

                first carry =1 - 2nd carry = 0 - third carry = 0 (13 11 10) check
                else if (firstcarry == 1 && secondCarryOne == 0 && thirdCarryZero == 0)
                {
                    printf(13 + 11 + 10 + 1 n);
                    stdstring binary_combined5 = binary13 + binary11 + binary10 + binary1;
                    if (int_negative != 1)
                    {
                        stdcout  The result in binary form is
                                   n;
                        stdcout  binary_combined5;
                    }

                    Length is the number of bits
                    int length = binary_combined5.length();
                    int lastchar = length - 1;
                    char char_array[length + 1];
                    strcpy(char_array, binary_combined5.c_str());

                    int256_t total = 0;

                    Positive or negative checking
                    char charvalue = char_array[0];
                    stdstring stringvalue(1, charvalue);
                    int intvalue = stdstoi(stringvalue);

                    if (int_negative == 2)
                    {
                        for (int i = 0; i  length; ++i)
                        {
                            char charvalue = char_array[i];
                            stdstring stringvalue(1, charvalue);
                            int intvalue = stdstoi(stringvalue);
                            int256_t calc1 = total  2;
                            total = calc1 + intvalue;
                        }
                    }
                    else
                    {
                        if (length == 128 && intvalue == 1)
                        {
                            int256_t lastcharvalue = 2;

                            for (int i = 0; i  lastchar - 1; ++i)
                            {
                                lastcharvalue = lastcharvalue  2;
                            }

                            for (int i = 1; i  length; ++i)
                            {
                                char charvalue = char_array[i];
                                stdstring stringvalue(1, charvalue);
                                int intvalue = stdstoi(stringvalue);
                                int256_t calc1 = total  2;
                                total = calc1 + intvalue;
                            }

                            total = total - lastcharvalue;
                        }
                        else
                        {
                            for (int i = 0; i  length; ++i)
                            {
                                char charvalue = char_array[i];
                                stdstring stringvalue(1, charvalue);
                                int intvalue = stdstoi(stringvalue);
                                int256_t calc1 = total  2;
                                total = calc1 + intvalue;
                            }
                        }
                    }
                    if (int_negative == 1)
                    {
                        int256_t invertedtotal = total  -1;
                        stdcout  The result in decimal form is
                                   n;
                        stdcout  invertedtotal;
                    }
                    else
                    {
                        stdcout  n
                                   n;
                        stdcout  The result in decimal form is
                                   n;
                        stdcout  total;
                    }

                    printf(n);
                    printf(n);
                    gettimeofday(&end, NULL);
                    get_time = (end.tv_sec - start.tv_sec) + (end.tv_usec - start.tv_usec)  1.0E-6;
                    printf(Computation Time %lf[sec]n, get_time);
                    printf(n);

                    printf(I hope you remembered what calculation you performed!n);
                }

                first carry =1 - 2nd carry = 0 - third carry = 1 (14 11 10) check
                else if (firstcarry == 1 && secondCarryOne == 0 && thirdCarryZero == 1)
                {
                    printf(14 + 11 + 10 + 1 n);
                    stdstring binary_combined6 = binary14 + binary11 + binary10 + binary1;
                    if (int_negative != 1)
                    {
                        stdcout  The result in binary form is
                                   n;
                        stdcout  binary_combined6;
                    }

                    Length is the number of bits
                    int length = binary_combined6.length();
                    int lastchar = length - 1;
                    char char_array[length + 1];
                    strcpy(char_array, binary_combined6.c_str());

                    int256_t total = 0;

                    Positive or negative checking
                    char charvalue = char_array[0];
                    stdstring stringvalue(1, charvalue);
                    int intvalue = stdstoi(stringvalue);

                    if (int_negative == 2)
                    {
                        for (int i = 0; i  length; ++i)
                        {
                            char charvalue = char_array[i];
                            stdstring stringvalue(1, charvalue);
                            int intvalue = stdstoi(stringvalue);
                            int256_t calc1 = total  2;
                            total = calc1 + intvalue;
                        }
                    }
                    else
                    {
                        if (length == 128 && intvalue == 1)
                        {
                            int256_t lastcharvalue = 2;

                            for (int i = 0; i  lastchar - 1; ++i)
                            {
                                lastcharvalue = lastcharvalue  2;
                            }

                            for (int i = 1; i  length; ++i)
                            {
                                char charvalue = char_array[i];
                                stdstring stringvalue(1, charvalue);
                                int intvalue = stdstoi(stringvalue);
                                int256_t calc1 = total  2;
                                total = calc1 + intvalue;
                            }

                            total = total - lastcharvalue;
                        }
                        else
                        {
                            for (int i = 0; i  length; ++i)
                            {
                                char charvalue = char_array[i];
                                stdstring stringvalue(1, charvalue);
                                int intvalue = stdstoi(stringvalue);
                                int256_t calc1 = total  2;
                                total = calc1 + intvalue;
                            }
                        }
                    }
                    if (int_negative == 1)
                    {
                        int256_t invertedtotal = total  -1;
                        stdcout  The result in decimal form is
                                   n;
                        stdcout  invertedtotal;
                    }
                    else
                    {
                        stdcout  n
                                   n;
                        stdcout  The result in decimal form is
                                   n;
                        stdcout  total;
                    }

                    printf(n);
                    printf(n);
                    gettimeofday(&end, NULL);
                    get_time = (end.tv_sec - start.tv_sec) + (end.tv_usec - start.tv_usec)  1.0E-6;
                    printf(Computation Time %lf[sec]n, get_time);
                    printf(n);

                    printf(I hope you remembered what calculation you performed!n);
                }

                first carry =1 - 2nd carry = 1 - third carry = 0 (13 12 10) check
                else if (firstcarry == 1 && secondCarryOne == 1 && thirdCarryOne == 0)
                {
                    printf(13 + 12 + 10 + 1 n);
                    stdstring binary_combined7 = binary13 + binary12 + binary10 + binary1;
                    if (int_negative != 1)
                    {
                        stdcout  The result in binary form is
                                   n;
                        stdcout  binary_combined7;
                    }

                    Length is the number of bits
                    int length = binary_combined7.length();
                    int lastchar = length - 1;
                    char char_array[length + 1];
                    strcpy(char_array, binary_combined7.c_str());

                    int256_t total = 0;

                    Positive or negative checking
                    char charvalue = char_array[0];
                    stdstring stringvalue(1, charvalue);
                    int intvalue = stdstoi(stringvalue);

                    if (int_negative == 2)
                    {
                        for (int i = 0; i  length; ++i)
                        {
                            char charvalue = char_array[i];
                            stdstring stringvalue(1, charvalue);
                            int intvalue = stdstoi(stringvalue);
                            int256_t calc1 = total  2;
                            total = calc1 + intvalue;
                        }
                    }
                    else
                    {
                        if (length == 128 && intvalue == 1)
                        {
                            int256_t lastcharvalue = 2;

                            for (int i = 0; i  lastchar - 1; ++i)
                            {
                                lastcharvalue = lastcharvalue  2;
                            }

                            for (int i = 1; i  length; ++i)
                            {
                                char charvalue = char_array[i];
                                stdstring stringvalue(1, charvalue);
                                int intvalue = stdstoi(stringvalue);
                                int256_t calc1 = total  2;
                                total = calc1 + intvalue;
                            }

                            total = total - lastcharvalue;
                        }
                        else
                        {
                            for (int i = 0; i  length; ++i)
                            {
                                char charvalue = char_array[i];
                                stdstring stringvalue(1, charvalue);
                                int intvalue = stdstoi(stringvalue);
                                int256_t calc1 = total  2;
                                total = calc1 + intvalue;
                            }
                        }
                    }
                    if (int_negative == 1)
                    {
                        int256_t invertedtotal = total  -1;
                        stdcout  The result in decimal form is
                                   n;
                        stdcout  invertedtotal;
                    }
                    else
                    {
                        stdcout  n
                                   n;
                        stdcout  The result in decimal form is
                                   n;
                        stdcout  total;
                    }

                    printf(n);
                    printf(n);
                    gettimeofday(&end, NULL);
                    get_time = (end.tv_sec - start.tv_sec) + (end.tv_usec - start.tv_usec)  1.0E-6;
                    printf(Computation Time %lf[sec]n, get_time);
                    printf(n);

                    printf(I hope you remembered what calculation you performed!n);
                }

                first carry =0 - 2nd carry = 1 - third carry = 0 (13 12 9) check
                else if (firstcarry == 0 && secondCarryZero == 1 && thirdCarryOne == 0)
                {
                    printf(13 + 12 + 9 + 1 n);
                    stdstring binary_combined8 = binary13 + binary12 + binary9 + binary1;
                    if (int_negative != 1)
                    {
                        stdcout  The result in binary form is
                                   n;
                        stdcout  binary_combined8;
                    }

                    Length is the number of bits
                    int length = binary_combined8.length();
                    int lastchar = length - 1;
                    char char_array[length + 1];
                    strcpy(char_array, binary_combined8.c_str());

                    int256_t total = 0;

                    Positive or negative checking
                    char charvalue = char_array[0];
                    stdstring stringvalue(1, charvalue);
                    int intvalue = stdstoi(stringvalue);

                    if (int_negative == 2)
                    {
                        for (int i = 0; i  length; ++i)
                        {
                            char charvalue = char_array[i];
                            stdstring stringvalue(1, charvalue);
                            int intvalue = stdstoi(stringvalue);
                            int256_t calc1 = total  2;
                            total = calc1 + intvalue;
                        }
                    }
                    else
                    {
                        if (length == 128 && intvalue == 1)
                        {
                            int256_t lastcharvalue = 2;

                            for (int i = 0; i  lastchar - 1; ++i)
                            {
                                lastcharvalue = lastcharvalue  2;
                            }

                            for (int i = 1; i  length; ++i)
                            {
                                char charvalue = char_array[i];
                                stdstring stringvalue(1, charvalue);
                                int intvalue = stdstoi(stringvalue);
                                int256_t calc1 = total  2;
                                total = calc1 + intvalue;
                            }

                            total = total - lastcharvalue;
                        }
                        else
                        {
                            for (int i = 0; i  length; ++i)
                            {
                                char charvalue = char_array[i];
                                stdstring stringvalue(1, charvalue);
                                int intvalue = stdstoi(stringvalue);
                                int256_t calc1 = total  2;
                                total = calc1 + intvalue;
                            }
                        }
                    }
                    if (int_negative == 1)
                    {
                        int256_t invertedtotal = total  -1;
                        stdcout  The result in decimal form is
                                   n;
                        stdcout  invertedtotal;
                    }
                    else
                    {
                        stdcout  n
                                   n;
                        stdcout  The result in decimal form is
                                   n;
                        stdcout  total;
                    }
                    printf(n);
                    printf(n);
                    gettimeofday(&end, NULL);
                    get_time = (end.tv_sec - start.tv_sec) + (end.tv_usec - start.tv_usec)  1.0E-6;
                    printf(Computation Time %lf[sec]n, get_time);
                    printf(n);

                    printf(I hope you remembered what calculation you performed!n);
                }

                else
                {
                    printf(Something went wrong... n);
                }
            }

            else if (int_bit == 256)
            {
                LweSample result = new_gate_bootstrapping_ciphertext_array(32, params);
                LweSample result2 = new_gate_bootstrapping_ciphertext_array(32, params);
                LweSample result3 = new_gate_bootstrapping_ciphertext_array(32, params);
                LweSample result4 = new_gate_bootstrapping_ciphertext_array(32, params);
                LweSample result5 = new_gate_bootstrapping_ciphertext_array(32, params);
                LweSample result6 = new_gate_bootstrapping_ciphertext_array(32, params);
                LweSample result7 = new_gate_bootstrapping_ciphertext_array(32, params);
                LweSample result8 = new_gate_bootstrapping_ciphertext_array(32, params);

                 export the 64 ciphertexts to a file (for the cloud)

                for (int i = 0; i  32; i++)
                    import_gate_bootstrapping_ciphertext_fromFile(answer_data, &result[i], params);

                for (int i = 0; i  32; i++)
                    import_gate_bootstrapping_ciphertext_fromFile(answer_data, &result2[i], params);

                for (int i = 0; i  32; i++)
                    import_gate_bootstrapping_ciphertext_fromFile(answer_data, &result3[i], params);

                for (int i = 0; i  32; i++)
                    import_gate_bootstrapping_ciphertext_fromFile(answer_data, &result4[i], params);

                for (int i = 0; i  32; i++)
                    import_gate_bootstrapping_ciphertext_fromFile(answer_data, &result5[i], params);

                for (int i = 0; i  32; i++)
                    import_gate_bootstrapping_ciphertext_fromFile(answer_data, &result6[i], params);

                for (int i = 0; i  32; i++)
                    import_gate_bootstrapping_ciphertext_fromFile(answer_data, &result7[i], params);

                for (int i = 0; i  32; i++)
                    import_gate_bootstrapping_ciphertext_fromFile(answer_data, &result8[i], params);
                fclose(answer_data);

                 decrypt and rebuild the answer
                int32_t int_answer1 = 0;
                for (int i = 0; i  32; i++)
                {
                    int ai = bootsSymDecrypt(&result[i], key)  0;
                    int_answer1 = (ai  i);
                }

                int32_t int_answer2 = 0;
                for (int i = 0; i  32; i++)
                {
                    int ai = bootsSymDecrypt(&result2[i], key)  0;
                    int_answer2 = (ai  i);
                }

                int32_t int_answer3 = 0;
                for (int i = 0; i  32; i++)
                {
                    int ai = bootsSymDecrypt(&result3[i], key)  0;
                    int_answer3 = (ai  i);
                }

                int32_t int_answer4 = 0;
                for (int i = 0; i  32; i++)
                {
                    int ai = bootsSymDecrypt(&result4[i], key)  0;
                    int_answer4 = (ai  i);
                }

                int32_t int_answer5 = 0;
                for (int i = 0; i  32; i++)
                {
                    int ai = bootsSymDecrypt(&result5[i], key)  0;
                    int_answer5 = (ai  i);
                }

                int32_t int_answer6 = 0;
                for (int i = 0; i  32; i++)
                {
                    int ai = bootsSymDecrypt(&result6[i], key)  0;
                    int_answer6 = (ai  i);
                }

                int32_t int_answer7 = 0;
                for (int i = 0; i  32; i++)
                {
                    int ai = bootsSymDecrypt(&result7[i], key)  0;
                    int_answer7 = (ai  i);
                }

                int32_t int_answer8 = 0;
                for (int i = 0; i  32; i++)
                {
                    int ai = bootsSymDecrypt(&result8[i], key)  0;
                    int_answer8 = (ai  i);
                }

                stdstring binary1 = stdbitset32(int_answer1).to_string();
                stdstring binary2 = stdbitset32(int_answer2).to_string();
                stdstring binary3 = stdbitset32(int_answer3).to_string();
                stdstring binary4 = stdbitset32(int_answer4).to_string();
                stdstring binary5 = stdbitset32(int_answer5).to_string();
                stdstring binary6 = stdbitset32(int_answer6).to_string();
                stdstring binary7 = stdbitset32(int_answer7).to_string();
                stdstring binary8 = stdbitset32(int_answer8).to_string();
                stdstring binary_combined = binary8 + binary7 + binary6 + binary5 + binary4 + binary3 + binary2 + binary1;

                if (int_negative != 1)
                {
                    stdcout  The result in binary form is
                               n;
                    stdcout  binary_combined;
                }

                Length is the number of bits
                int length = binary_combined.length();
                int lastchar = length - 1;
                char char_array[length + 1];
                strcpy(char_array, binary_combined.c_str());

                int256_t total = 0;

                Positive or negative checking
                char charvalue = char_array[0];
                stdstring stringvalue(1, charvalue);
                int intvalue = stdstoi(stringvalue);

                if (int_negative == 2)
                {
                    for (int i = 0; i  length; ++i)
                    {
                        char charvalue = char_array[i];
                        stdstring stringvalue(1, charvalue);
                        int intvalue = stdstoi(stringvalue);
                        int256_t calc1 = total  2;
                        total = calc1 + intvalue;
                    }
                }
                else
                {
                    if (length == 256 && intvalue == 1)
                    {
                        int256_t lastcharvalue = 2;

                        for (int i = 0; i  lastchar - 1; ++i)
                        {
                            lastcharvalue = lastcharvalue  2;
                        }

                        for (int i = 1; i  length; ++i)
                        {
                            char charvalue = char_array[i];
                            stdstring stringvalue(1, charvalue);
                            int intvalue = stdstoi(stringvalue);
                            int256_t calc1 = total  2;
                            total = calc1 + intvalue;
                        }
                        total = total - lastcharvalue;
                    }
                    else
                    {
                        for (int i = 0; i  length; ++i)
                        {
                            char charvalue = char_array[i];
                            stdstring stringvalue(1, charvalue);
                            int intvalue = stdstoi(stringvalue);
                            int256_t calc1 = total  2;
                            total = calc1 + intvalue;
                        }
                    }
                }

                if (int_negative == 1)
                {
                    int256_t invertedtotal = total  -1;
                    stdcout  The result in decimal form is
                               n;
                    stdcout  invertedtotal;
                }
                else
                {
                    stdcout  n
                               n;
                    stdcout  The result in decimal form is
                               n;
                    stdcout  total;
                }
                printf(n);
                printf(n);
                gettimeofday(&end, NULL);
                get_time = (end.tv_sec - start.tv_sec) + (end.tv_usec - start.tv_usec)  1.0E-6;
                printf(Computation Time %lf[sec]n, get_time);
                printf(n);

                printf(I hope you remembered what calculation you performed!n);

                 clean up all pointers
                delete_gate_bootstrapping_ciphertext_array(32, result);
                delete_gate_bootstrapping_ciphertext_array(32, result2);
                delete_gate_bootstrapping_ciphertext_array(32, result3);
                delete_gate_bootstrapping_ciphertext_array(32, result4);
                delete_gate_bootstrapping_ciphertext_array(32, result5);
                delete_gate_bootstrapping_ciphertext_array(32, result6);
                delete_gate_bootstrapping_ciphertext_array(32, result7);
                delete_gate_bootstrapping_ciphertext_array(32, result8);
                delete_gate_bootstrapping_secret_keyset(key);
                delete_gate_bootstrapping_secret_keyset(nbitkey);
            }
        }

        else if (int_op == 4)
        {
            stdcout  Result for   int_bit   bit Multiplication computation
                       n
                       n;

            if (int_bit == 256)
            {
                LweSample finalresult1 = new_gate_bootstrapping_ciphertext_array(32, params);
                LweSample finalresult2 = new_gate_bootstrapping_ciphertext_array(32, params);
                LweSample finalresult3 = new_gate_bootstrapping_ciphertext_array(32, params);
                LweSample finalresult4 = new_gate_bootstrapping_ciphertext_array(32, params);
                LweSample finalresult5 = new_gate_bootstrapping_ciphertext_array(32, params);
                LweSample finalresult6 = new_gate_bootstrapping_ciphertext_array(32, params);
                LweSample finalresult7 = new_gate_bootstrapping_ciphertext_array(32, params);
                LweSample finalresult8 = new_gate_bootstrapping_ciphertext_array(32, params);

                export the 32 ciphertexts to a file (for the cloud)

                for (int i = 0; i  32; i++)
                    import_gate_bootstrapping_ciphertext_fromFile(answer_data, &finalresult1[i], params);

                for (int i = 0; i  32; i++)
                    import_gate_bootstrapping_ciphertext_fromFile(answer_data, &finalresult2[i], params);

                for (int i = 0; i  32; i++)
                    import_gate_bootstrapping_ciphertext_fromFile(answer_data, &finalresult3[i], params);

                for (int i = 0; i  32; i++)
                    import_gate_bootstrapping_ciphertext_fromFile(answer_data, &finalresult4[i], params);

                for (int i = 0; i  32; i++)
                    import_gate_bootstrapping_ciphertext_fromFile(answer_data, &finalresult5[i], params);

                for (int i = 0; i  32; i++)
                    import_gate_bootstrapping_ciphertext_fromFile(answer_data, &finalresult6[i], params);

                for (int i = 0; i  32; i++)
                    import_gate_bootstrapping_ciphertext_fromFile(answer_data, &finalresult7[i], params);

                for (int i = 0; i  32; i++)
                    import_gate_bootstrapping_ciphertext_fromFile(answer_data, &finalresult8[i], params);
                fclose(answer_data);

                decrypt and rebuild the answer
                int32_t int_answer = 0;
                for (int i = 0; i  32; i++)
                {
                    int ai = bootsSymDecrypt(&finalresult1[i], key)  0;
                    int_answer = (ai  i);
                }

                int32_t int_answer2 = 0;
                for (int i = 0; i  32; i++)
                {
                    int ai = bootsSymDecrypt(&finalresult2[i], key)  0;
                    int_answer2 = (ai  i);
                }

                int32_t int_answer3 = 0;
                for (int i = 0; i  32; i++)
                {
                    int ai = bootsSymDecrypt(&finalresult3[i], key)  0;
                    int_answer3 = (ai  i);
                }

                int32_t int_answer4 = 0;
                for (int i = 0; i  32; i++)
                {
                    int ai = bootsSymDecrypt(&finalresult4[i], key)  0;
                    int_answer4 = (ai  i);
                }

                int32_t int_answer5 = 0;
                for (int i = 0; i  32; i++)
                {
                    int ai = bootsSymDecrypt(&finalresult5[i], key)  0;
                    int_answer5 = (ai  i);
                }

                int32_t int_answer6 = 0;
                for (int i = 0; i  32; i++)
                {
                    int ai = bootsSymDecrypt(&finalresult6[i], key)  0;
                    int_answer6 = (ai  i);
                }

                int32_t int_answer7 = 0;
                for (int i = 0; i  32; i++)
                {
                    int ai = bootsSymDecrypt(&finalresult7[i], key)  0;
                    int_answer7 = (ai  i);
                }

                int32_t int_answer8 = 0;
                for (int i = 0; i  32; i++)
                {
                    int ai = bootsSymDecrypt(&finalresult8[i], key)  0;
                    int_answer8 = (ai  i);
                }

                stdstring binary1 = stdbitset32(int_answer).to_string();
                stdstring binary2 = stdbitset32(int_answer2).to_string();
                stdstring binary3 = stdbitset32(int_answer3).to_string();
                stdstring binary4 = stdbitset32(int_answer4).to_string();
                stdstring binary5 = stdbitset32(int_answer5).to_string();
                stdstring binary6 = stdbitset32(int_answer6).to_string();
                stdstring binary7 = stdbitset32(int_answer7).to_string();
                stdstring binary8 = stdbitset32(int_answer8).to_string();
                stdstring binary_combined = binary8 + binary7 + binary6 + binary5 + binary4 + binary3 + binary2 + binary1;

                int length = binary_combined.length();
                char char_array[length + 1];
                strcpy(char_array, binary_combined.c_str());
                int256_t total = 0;

                if (int_negative == 0  int_negative == 4)
                {
                    stdcout  The result in binary form is
                               n;
                    stdcout  binary_combined;
                }

                for (int i = 0; i  length; ++i)
                {
                    char charvalue = char_array[i];
                    stdstring stringvalue(1, charvalue);
                    int intvalue = stdstoi(stringvalue);
                    int256_t calc1 = total  2;
                    total = calc1 + intvalue;
                }

                if (int_negative == 1  int_negative == 2)
                {
                    int256_t invertedtotal = total  -1;
                    stdcout  The result in decimal form is
                               n;
                    stdcout  invertedtotal;
                }
                else
                {
                    stdcout  n
                               n;
                    stdcout  The result in decimal form is
                               n;
                    stdcout  total;
                }

                printf(n);
                printf(n);
                gettimeofday(&end, NULL);
                get_time = (end.tv_sec - start.tv_sec) + (end.tv_usec - start.tv_usec)  1.0E-6;
                printf(Computation Time %lf[sec]n, get_time);
                printf(n);

                printf(I hope you remembered what calculation you performed!n);

                clean up all pointers
                delete_gate_bootstrapping_ciphertext_array(32, finalresult1);
                delete_gate_bootstrapping_ciphertext_array(32, finalresult2);
                delete_gate_bootstrapping_ciphertext_array(32, finalresult3);
                delete_gate_bootstrapping_ciphertext_array(32, finalresult4);
                delete_gate_bootstrapping_ciphertext_array(32, finalresult5);
                delete_gate_bootstrapping_ciphertext_array(32, finalresult6);
                delete_gate_bootstrapping_ciphertext_array(32, finalresult7);
                delete_gate_bootstrapping_ciphertext_array(32, finalresult8);
                delete_gate_bootstrapping_secret_keyset(key);
                delete_gate_bootstrapping_secret_keyset(nbitkey);
            }

            else if (int_bit == 128)
            {
                LweSample finalresult = new_gate_bootstrapping_ciphertext_array(32, params);
                LweSample finalresult2 = new_gate_bootstrapping_ciphertext_array(32, params);
                LweSample finalresult3 = new_gate_bootstrapping_ciphertext_array(32, params);
                LweSample finalresult4 = new_gate_bootstrapping_ciphertext_array(32, params);

                export the 32 ciphertexts to a file (for the cloud)

                for (int i = 0; i  32; i++)
                    import_gate_bootstrapping_ciphertext_fromFile(answer_data, &finalresult[i], params);

                for (int i = 0; i  32; i++)
                    import_gate_bootstrapping_ciphertext_fromFile(answer_data, &finalresult2[i], params);

                for (int i = 0; i  32; i++)
                    import_gate_bootstrapping_ciphertext_fromFile(answer_data, &finalresult3[i], params);

                for (int i = 0; i  32; i++)
                    import_gate_bootstrapping_ciphertext_fromFile(answer_data, &finalresult4[i], params);

                fclose(answer_data);

                decrypt and rebuild the answer
                int32_t int_answer = 0;
                for (int i = 0; i  32; i++)
                {
                    int ai = bootsSymDecrypt(&finalresult[i], key)  0;
                    int_answer = (ai  i);
                }

                int32_t int_answer2 = 0;
                for (int i = 0; i  32; i++)
                {
                    int ai = bootsSymDecrypt(&finalresult2[i], key)  0;
                    int_answer2 = (ai  i);
                }

                int32_t int_answer3 = 0;
                for (int i = 0; i  32; i++)
                {
                    int ai = bootsSymDecrypt(&finalresult3[i], key)  0;
                    int_answer3 = (ai  i);
                }

                int32_t int_answer4 = 0;
                for (int i = 0; i  32; i++)
                {
                    int ai = bootsSymDecrypt(&finalresult4[i], key)  0;
                    int_answer4 = (ai  i);
                }

                stdstring binary1 = stdbitset32(int_answer).to_string();
                stdstring binary2 = stdbitset32(int_answer2).to_string();
                stdstring binary3 = stdbitset32(int_answer3).to_string();
                stdstring binary4 = stdbitset32(int_answer4).to_string();
                stdstring binary_combined = binary4 + binary3 + binary2 + binary1;

                int length = binary_combined.length();
                char char_array[length + 1];
                strcpy(char_array, binary_combined.c_str());
                int256_t total = 0;

                if (int_negative == 0  int_negative == 4)
                {
                    stdcout  The result in binary form is
                               n;
                    stdcout  binary_combined;
                }

                for (int i = 0; i  length; ++i)
                {
                    char charvalue = char_array[i];
                    stdstring stringvalue(1, charvalue);
                    int intvalue = stdstoi(stringvalue);
                    int256_t calc1 = total  2;
                    total = calc1 + intvalue;
                }

                if (int_negative == 1  int_negative == 2)
                {
                    int256_t invertedtotal = total  -1;
                    stdcout  The result in decimal form is
                               n;
                    stdcout  invertedtotal;
                }
                else
                {
                    stdcout  n
                               n;
                    stdcout  The result in decimal form is
                               n;
                    stdcout  total;
                }

                printf(n);
                printf(n);
                gettimeofday(&end, NULL);
                get_time = (end.tv_sec - start.tv_sec) + (end.tv_usec - start.tv_usec)  1.0E-6;
                printf(Computation Time %lf[sec]n, get_time);
                printf(n);
                printf(I hope you remembered what calculation you performed!n);

                clean up all pointers
                delete_gate_bootstrapping_ciphertext_array(32, finalresult);
                delete_gate_bootstrapping_ciphertext_array(32, finalresult2);
                delete_gate_bootstrapping_ciphertext_array(32, finalresult3);
                delete_gate_bootstrapping_ciphertext_array(32, finalresult4);
                delete_gate_bootstrapping_secret_keyset(key);
                delete_gate_bootstrapping_secret_keyset(nbitkey);
            }

            else if (int_bit == 64)
            {
                printf(Importing ciphertext...n);

                LweSample finalresult = new_gate_bootstrapping_ciphertext_array(32, params);
                LweSample finalresult2 = new_gate_bootstrapping_ciphertext_array(32, params);

                LweSample result = new_gate_bootstrapping_ciphertext_array(32, params);

                export the 32 ciphertexts to a file (for the cloud)

                for (int i = 0; i  32; i++)
                    import_gate_bootstrapping_ciphertext_fromFile(answer_data, &finalresult[i], params);

                for (int i = 0; i  32; i++)
                    import_gate_bootstrapping_ciphertext_fromFile(answer_data, &finalresult2[i], params);

                fclose(answer_data);

                decrypt and rebuild the answer
                int32_t int_answer = 0;
                for (int i = 0; i  32; i++)
                {
                    int ai = bootsSymDecrypt(&finalresult[i], key)  0;
                    int_answer = (ai  i);
                }

                int32_t int_answer2 = 0;
                for (int i = 0; i  32; i++)
                {
                    int ai = bootsSymDecrypt(&finalresult2[i], key)  0;
                    int_answer2 = (ai  i);
                }

                stdstring binary1 = stdbitset32(int_answer).to_string();
                stdstring binary2 = stdbitset32(int_answer2).to_string();
                stdstring binary_combined = binary2 + binary1;

                int length = binary_combined.length();
                char char_array[length + 1];
                strcpy(char_array, binary_combined.c_str());
                int256_t total = 0;

                stdcout  The result in binary form is
                           n;
                stdcout  binary_combined;

                for (int i = 0; i  length; ++i)
                {
                    char charvalue = char_array[i];
                    stdstring stringvalue(1, charvalue);
                    int intvalue = stdstoi(stringvalue);
                    int256_t calc1 = total  2;
                    total = calc1 + intvalue;
                }

                stdcout  n
                           n;
                stdcout  The result in decimal form is
                           n;
                stdcout  total;

                printf(n);
                printf(n);
                gettimeofday(&end, NULL);
                get_time = (end.tv_sec - start.tv_sec) + (end.tv_usec - start.tv_usec)  1.0E-6;
                printf(Computation Time %lf[sec]n, get_time);
                printf(n);

                printf(I hope you remembered what calculation you performed!n);

                printf(And the result is %dn, int_answer);
                printf(Plaintext1 after decryption is %dn, int_ciphertext1);
                printf(Plaintext2 after decryption is %dn, int_ciphertext2);
                printf(I hope you remembered what calculation you performed!n);

                 clean up all pointers
                delete_gate_bootstrapping_ciphertext_array(32, int_answer);
                delete_gate_bootstrapping_secret_keyset(key);
            }
        }
        else
        {
            printf(FAILn);
        }
    }

    fclose(f);
    return 0;
}