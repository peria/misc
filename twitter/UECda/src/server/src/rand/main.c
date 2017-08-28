#include <stdio.h>
#include "test.h"

int main(void)
{
    int i;
    tn_rand_init(0xffffffUL,0);
    for (i=0; i<100; i++) {
      printf("%10lu ", tn_rand_gen(0));
      if (i%5==4) printf("\n");
    }
    return 0;
}
