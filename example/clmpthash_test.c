#include "../include/clmpthash.h"

int main(int argc, char** argv) {
    clmpthash_config cfg;
    LVA* lvas;
    PhysicalAddr* pas;
    LVA* querys;
    uint64_t num_lva, num_querys;
    parse_configuration("lmpthash_config", &cfg, &lvas, &pas, &num_lva, &querys, &num_querys);

    void* index = build_index(lvas, pas, num_lva, &cfg);

    PhysicalAddr pa;
    for(int i = 0; i < num_querys; ++i) {
        int ret = get_pa(querys[i], index, &pa);
        
        // check if the result is correct: first 8 bytes should be equal to the original lva
        if(ret==0){
            LVA _lva=0;
            for(int j=0;j<8;j++){
                _lva=(_lva<<8)+pa.data[j];
            }
            if(_lva!=querys[i]){
                printf("wrong result: should be 0x%lx, but got 0x%lx\n", querys[i], _lva);
                return -1;
            }
        }else{
            printf("error getting pa\n");
            return -1;
        }
    }
    printf("test passed\n");
    return 0;

}