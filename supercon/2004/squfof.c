#include <stdio.h>
#include <ctype.h>
#include <string.h>

typedef unsigned long long uint64;

typedef struct {
  uint64 h;  /* high */
  uint64 l;  /* low */
} Uint128;

/* c = a * b */
void Uint128Mul(Uint128* a, int b, Uint128* c) {

}

void Parse(char* str, Uint128* input) {
  input->h = 0;
  input->l = 0;

  if (strncmp(str, "0x", 2) == 0 || strncmp(str, "0X", 2) == 0) {
    /* hex mode */
    char* begin = str + 2;
    char* end = begin;
    for (; isxdigit(*end); ++end) {}
    char* lbegin = end - 16;
    if (lbegin < begin)
      lbegin = begin;
    sscanf(lbegin, "%llx", &input->l);

    if (begin == lbegin) {
      input->h = 0;
    }
    *lbegin = '\0';
    sscanf(begin, "%llx", &input->h);

    return;
  }

  /* decimal */

}

int main(int argc, char* argv[]) {
  if (argc < 2) {
    fprintf(stderr, "Usage: %s number\n", argv[0]);
    return 0;
  }

  Uint128 input;
  Parse(argv[1], &input);

  return 0;
}
