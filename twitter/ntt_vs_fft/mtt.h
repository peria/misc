#include "fmt.h"

// MTT stands for Montgomery NTT
class MTT {
 public:
  static void SetGlobals(void* work, void* table, void*);
  static bool Validate(uint64* data);
  static void Square(uint64* data, int n);

  static void InitTable(int log2n);
  static void Forward(int log2n, int n, uint64* data);
  static void Backward(int log2n, int n, uint64* data);

 private:
  static void Core(int log2n, int n, uint64* table,
                   uint64* work, uint64* data);
  static void Radix2(const int width, const int height,
                     uint64* ptr, uint64* x, uint64* y);
};
