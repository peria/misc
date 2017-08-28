#include"test.h"

void tn_rand_init(unsigned long seed, int type){
	switch(type){
		case 0 :
			srand((int)seed);
			break;
		case 1 :
			init_genrand((unsigned long)seed);
			break;
	}
};


long  tn_rand_gen(type){
	switch(type){
		case 0 :
			return rand();
		case 1 :
			return genrand_int32();
	}
};
