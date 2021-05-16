#include <tfhe/tfhe.h>
#include <tfhe/tfhe_io.h>
#include <stdio.h>


int main() {

	// reads the cloud key from file
	FILE* secret_key = fopen("secret.key","rb");
	TFheGateBootstrappingSecretKeySet* key = new_tfheGateBootstrappingSecretKeySet_fromFile(secret_key);
    	fclose(secret_key);

	// if necessary, the params are inside the key
    	const TFheGateBootstrappingParameterSet* params = key->params;

	// read the 32 ciphertexts of the result
    	LweSample* answer = new_gate_bootstrapping_ciphertext_array(32, params);
    	LweSample* checkans = new_gate_bootstrapping_ciphertext_array(32, params);
    	LweSample* checkdenominator = new_gate_bootstrapping_ciphertext_array(32, params);

	// export the 32 ciphertexts to a file (for the cloud)
    	FILE* answer_data = fopen("answer.data","rb");
    	for (int i=0; i<32; i++)
        	import_gate_bootstrapping_ciphertext_fromFile(answer_data, &answer[i], params);

    	for (int i=0; i<32; i++)
        	import_gate_bootstrapping_ciphertext_fromFile(answer_data, &checkans[i], params);

    	for (int i=0; i<32; i++)
        	import_gate_bootstrapping_ciphertext_fromFile(answer_data, &checkdenominator[i], params);

    	fclose(answer_data);

	//decrypt and rebuild the answer
    	uint32_t int_answer = 0;
    	uint32_t check_answer = 0;
	uint32_t check_denominator = 0;
    	for (int i=0; i<32; i++) {
        	int ai = bootsSymDecrypt(&answer[i], key)>0;
        	int_answer |= (ai<<i);
    	}

    	for (int i=0; i<32; i++) {
        	int ai = bootsSymDecrypt(&checkans[i], key)>0;
        	check_answer |= (ai<<i);
    	}

    	for (int i=0; i<32; i++) {
        	int ai = bootsSymDecrypt(&checkdenominator[i], key)>0;
        	check_denominator |= (ai<<i);
	}

    	printf("The Quotient is: %d\n",int_answer);
    	printf("The Remainder is: %d\n",check_answer);
    	printf("The Denominator is: %d\n",check_denominator);

	float n, sum;
	n  = (float)check_answer/(float)check_denominator;
	sum = int_answer + (float)n;
  	printf("%d / %d = %.3f \n",check_answer,check_denominator,(float)n);
  	printf("%d + (%d / %d) = %.3f \n", int_answer, check_answer, check_denominator, (float)sum);
    	printf("I hope you remembered what calculation you performed!\n");

	//clean up all pointers
	delete_gate_bootstrapping_ciphertext_array(32, answer);
	delete_gate_bootstrapping_secret_keyset(key);

}
