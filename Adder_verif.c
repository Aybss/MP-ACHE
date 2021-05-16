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

    // read the 64 ciphertexts of the result
    LweSample* result = new_gate_bootstrapping_ciphertext_array(64, params);
    LweSample* ciphertext1 = new_gate_bootstrapping_ciphertext_array(64, params);
    LweSample* ciphertext2 = new_gate_bootstrapping_ciphertext_array(64, params);

    // export the 64 ciphertexts to a file (for the cloud)
    FILE* answer_data = fopen("answer.data","rb");
    for (int i=0; i<64; i++)
        import_gate_bootstrapping_ciphertext_fromFile(answer_data, &result[i], params);
    fclose(answer_data);

    FILE* ciphertext1_data = fopen("ciphertext1.data","rb");
    for (int i=0; i<64; i++)
        import_gate_bootstrapping_ciphertext_fromFile(ciphertext1_data, &ciphertext1[i], params);
    fclose(ciphertext1_data);

    FILE* ciphertext2_data = fopen("ciphertext2.data","rb");
    for (int i=0; i<64; i++)
        import_gate_bootstrapping_ciphertext_fromFile(ciphertext2_data, &ciphertext2[i], params);
    fclose(ciphertext2_data);

    // decrypt and rebuild the answer
    int64_t int_answer = 0;
    for (int i=0; i<64; i++) {
        int ai = bootsSymDecrypt(&result[i], key)>0;
        int_answer |= (ai<<i);
    }

    int64_t int_ciphertext1 = 0;
    for (int i=0; i<64; i++) {
        int bi = bootsSymDecrypt(&ciphertext1[i], key)>0;
        int_ciphertext1 |= (bi<<i);
    }

    int64_t int_ciphertext2 = 0;
    for (int i=0; i<64; i++) {
        int ci = bootsSymDecrypt(&ciphertext2[i], key)>0;
        int_ciphertext2 |= (ci<<i);
    }

    printf("And the result is: %lu\n",int_answer);
    printf("Plaintext1 after decryption is: %lu\n",int_ciphertext1);
    printf("Plaintext2 after decryption is: %lu\n",int_ciphertext2);
    printf("I hope you remembered what calculation you performed!\n");

    // clean up all pointers
    delete_gate_bootstrapping_ciphertext_array(64, result);
    delete_gate_bootstrapping_secret_keyset(key);
}
