#include <algorithm>
#include <iostream>
#include <numeric>
#include <vector>
using namespace std;

using int64 = long long;

int main() {
  for (int64 n; cin >> n && n;) {
    vector<int64> a(n);
    for (int64 i = 0; i < n; ++i)
      cin >> a[i];
    int64 sum = accumulate(a.begin(), a.end(), 0LL);
    int64 avg = sum / n;
    sort(a.begin(), a.end());
    auto itr = upper_bound(a.begin(), a.end(), avg);
    cout << (itr - a.begin()) << "\n";
  }
  return 0;
}
