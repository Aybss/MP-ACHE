#include <tfhe/tfhe.h>
#include <tfhe/tfhe_io.h>
#include <stdio.h>

int main() {

	// generate a keyset
	const int minimum_lambda = 110;
	TFheGateBootstrappingParameterSet* params = new_default_gate_bootstrapping_parameters(minimum_lambda);

	// generate a random key
	uint32_t seed[] = { 314, 1592, 657 };
	tfhe_random_generator_setSeed(seed,3);
	TFheGateBootstrappingSecretKeySet* key = new_random_gate_bootstrapping_secret_keyset(params);

	//generate encrypt the 32 bits
	uint32_t plaintext1 = 20;
	LweSample* ciphertext1 = new_gate_bootstrapping_ciphertext_array(32, params);
	for (int i=0; i<32; i++)
	{
		bootsSymEncrypt(&ciphertext1[i], (plaintext1>>i)&1, key);
	}

	// generate encrypt the 32 bits
	uint32_t plaintext2 = 5;
	LweSample* ciphertext2 = new_gate_bootstrapping_ciphertext_array(32, params);
	for (int i=0; i<32; i++)
	{
		bootsSymEncrypt(&ciphertext2[i], (plaintext2>>i)&1, key);
	}

	//carry
	uint16_t plaintext3 = 0;
	LweSample* ciphertext3 = new_gate_bootstrapping_ciphertext_array(16, params);
	for (int i = 0; i<16; i++) {
		bootsSymEncrypt(&ciphertext3[i], (plaintext3 >> i) & 1, key);
	}

	printf("Hi there! Today, I will ask the cloud the calculation results of the two data you input.\n");


	//export the secret key to file for later use
	FILE* secret_key = fopen("secret.key","wb");
	export_tfheGateBootstrappingSecretKeySet_toFile(secret_key, key);
	fclose(secret_key);

	//export the cloud key to a file (for the cloud)
	FILE* cloud_key = fopen("cloud.key","wb");
	export_tfheGateBootstrappingCloudKeySet_toFile(cloud_key, &key->cloud);
	fclose(cloud_key);

	//export the 32 ciphertexts to a file (for the cloud)
	FILE* cloud_data = fopen("cloud.data","wb");
	for (int i=0; i<32; i++)
		export_gate_bootstrapping_ciphertext_toFile(cloud_data, &ciphertext1[i], params);
	for (int i=0; i<32; i++)
 		export_gate_bootstrapping_ciphertext_toFile(cloud_data, &ciphertext2[i], params);
	for (int i = 0; i<16; i++)
		export_gate_bootstrapping_ciphertext_toFile(cloud_data, &ciphertext3[i], params);
	fclose(cloud_data);

	//clean up all pointers
	delete_gate_bootstrapping_ciphertext_array(32, ciphertext2);
	delete_gate_bootstrapping_ciphertext_array(32, ciphertext1);
	delete_gate_bootstrapping_ciphertext_array(16, ciphertext3);
	delete_gate_bootstrapping_secret_keyset(key);
	delete_gate_bootstrapping_parameters(params);

}
