#include <iostream>
using namespace std;

using int64 = long long;

void Solve(int64 b) {
  int64 low = 1, up = 1, diff = -b;
  while (diff) {
    while (diff < 0)
      diff += up++;
    while (diff > 0)
      diff -= low++;
  }
  cout << low << " " << (up - low) << "\n";
}

int main() {
  for (int64 b; cin >> b && b;) {
    Solve(b);
  }
  return 0;
}
