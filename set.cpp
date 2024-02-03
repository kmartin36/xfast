#include <iostream>
#include <utility>
#include <random>
#include <set>
#include <chrono>
#include <ctime>
#include <cmath>
#include <sys/resource.h>

// Compile with: g++ -O3 -std=c++11 -o set set.cpp
// Run with: ./set 1000000

int main(int argc, char **argv) {
  std::set<unsigned long long> b;
//   std::random_device rd;
  std::mt19937_64 gen(55);
  std::chrono::time_point<std::chrono::high_resolution_clock> start, end;
  std::chrono::duration<double> elapsed_b(0.0);
  std::cout << "N rb_insert(us) rb_lower_bound(us) rusage(kB)" << std::endl;
  
  long tcount = 0;
  while (*++argv) {
    long count = atol(*argv) - tcount;
    tcount += count;
    for (long i=0; i<count; i++) {
      unsigned long long k = gen();
      start = std::chrono::high_resolution_clock::now();
      b.emplace(k);
      end = std::chrono::high_resolution_clock::now();
      elapsed_b += end - start;
    }
    std::cout << tcount << ' ' << elapsed_b.count()*1000000/count << ' ';
    elapsed_b = std::chrono::duration<double>::zero();
    unsigned long long good = 0;
    for (int i=0; i<1000000; i++) {
      unsigned long long x = gen(), k0, k1;
      end = std::chrono::high_resolution_clock::now();
      auto it = b.lower_bound(x);
      start = std::chrono::high_resolution_clock::now();
      good += *it;
      elapsed_b += start - end;
    }
    struct rusage usage;
    getrusage(RUSAGE_SELF, &usage);
    std::cout << elapsed_b.count() << ' ' << usage.ru_maxrss << ' ' << (good != 0) << std::endl;
    elapsed_b = std::chrono::duration<double>::zero();
  }
  
  return 0;
}
