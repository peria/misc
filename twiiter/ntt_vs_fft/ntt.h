#include "fmt.h"

class NTT {
 public:
  static void SetGlobals(void* work, void* table, void*);
  static bool Validate(uint64* data);

  static void InitTable(int log2n);
  static void Forward(int log2n, int n, uint64* data);
  static void Backward(int log2n, int n, uint64* data);

 private:
  static void Radix2(const int width, const int height,
                     uint64* ptr, uint64* x, uint64* y);
};
