#include <iostream>
#include <utility>
#include <random>
#include <vector>
#include <limits>
#include <set>
#include <chrono>
#include <ctime>
#include <cmath>
#include <sys/resource.h>

// Compile with: g++ -O3 -std=c++11 -o ktrie ktrie.cpp
// Run with: ./ktrie 1000000

using std::cout;
using std::endl;

template <class key>
class ktrie {
public:
  ktrie(std::size_t expected = 1) {
    expected = std::max(expected, (std::size_t)1);
    shift = __builtin_clzll(expected-1);
    mask = ~((1ULL << shift) - 1);
    count = 0;
    tbl.resize(a(mask) + 1);
  }
  bool contains(key x) {return tbl[a(x)].count() > 0;};
  key lower_bound(key x) {
    auto i = tbl[a(x)].lower_bound(b(x));
    if (i != tbl[a(x)].end())
      return *i | (mask & x);
    for (auto i=tbl.begin() + a(x)+1; i!=tbl.end(); i++)
      if (!i->empty())
        return *i->begin() | ((i-tbl.begin())<<shift);
    return std::numeric_limits<key>::max();
  }
  key upper_bound(key x) {
    auto i = tbl[a(x)].upper_bound(b(x));
    if (i != tbl[a(x)].end())
      return *i | (mask & x);
    for (auto i=tbl.begin() + a(x)+1; i!=tbl.end(); i++)
      if (!i->empty())
        return *i->begin() | ((i-tbl.begin())<<shift);
    return std::numeric_limits<key>::max();
  }
  bool insert(key x) {
    if (!tbl[a(x)].emplace(b(x)).second)
      return false;
    if (count == (mask >> shift)) {
      // resize the array
      shift--;
      mask |= 1ULL << shift;
      std::size_t s = tbl.size();
      tbl.resize(a(mask)+1);
      for (std::size_t i=s; --i;) {
        for (auto j : tbl[i]) {
          key c = (i << (shift+1)) | j;
          tbl[a(c)].emplace_hint(tbl[a(c)].end(), b(c));
        }
        tbl[i].clear();
      }
      auto i = tbl[0].lower_bound(((key)1) << shift);
      tbl[1].insert(i, tbl[0].end());
      tbl[0].erase(i, tbl[0].end());
    }
    count++;
    return true;
  }
  bool remove(key x) {
    std::size_t s = tbl[a(x)].erase(b(x));
    if (s == 0)
      return false;
    if (mask && count == (mask ^ (mask >> 1))) {
      // resize the array
      mask <<= 1;
      shift++;
      s = tbl.size();
      for (std::size_t i=1; i<s; i++) {
        for (auto j : tbl[i]) {
          key c = (i << (shift-1)) | j;
          tbl[i>>1].emplace_hint(tbl[i>>1].end(), b(c));
        }
        tbl[i].clear();
      }
      tbl.resize(a(mask)+1);
    }
    count--;
    return true;
  }
  std::size_t getDepth() {
    std::size_t s=0;
    for (auto const &i : tbl)
      s = std::max(s, i.size());
    return s;
  }

private:
  std::vector<std::set<key> > tbl;
  key mask;
  std::size_t count;
  uint8_t shift;
  key a(key x) {return (x & mask) >> shift;};
  key b(key x) {return x & (~mask);};
};

int main(int argc, char **argv) {
  // std::set<unsigned long long> b;
  ktrie<unsigned long long> c;
//   std::random_device rd;
  std::mt19937_64 gen(55);
  std::chrono::time_point<std::chrono::high_resolution_clock> start, end;
  std::chrono::duration<double> elapsed_b(0.0), elapsed_c(0.0);
  std::cout << "N rb_insert(us) kt_insert(us) rb_lower_bound(us) kt_lower_bound(us) rusage(kB) tree_depth" << std::endl;
  
  long tcount = 0;
  while (*++argv) {
    long count = atol(*argv) - tcount;
    tcount += count;
    for (long i=0; i<count; i++) {
      unsigned long long k = gen();
      // start = std::chrono::high_resolution_clock::now();
      // b.emplace(k);
      end = std::chrono::high_resolution_clock::now();
      // elapsed_b += end - start;
      c.insert(k);
      start = std::chrono::high_resolution_clock::now();
      elapsed_c += start - end;
    }
    std::cout << tcount << ' ' << elapsed_b.count()*1000000/count << ' ' << elapsed_c.count()*1000000/count << ' ';
    elapsed_b = elapsed_c = std::chrono::duration<double>::zero();
    unsigned long long good = 0;
    for (int i=0; i<1000000; i++) {
      unsigned long long x = gen(), k0, k1;
      start = std::chrono::high_resolution_clock::now();
      k0 = c.lower_bound(x);
      end = std::chrono::high_resolution_clock::now();
      good += k0;
      elapsed_c += end - start;
      // auto it = b.lower_bound(x);
      // start = std::chrono::high_resolution_clock::now();
      // elapsed_b += start - end;
      // if (it == b.end())
      //   k1 = std::numeric_limits<unsigned long long>::max();
      // else
      //   k1 = *it;
      // if (k0 == k1)
      //   good++;
      // else
      //   printf("%016llx %016llx %016llx\n", k0, k1, x);
    }
    struct rusage usage;
    getrusage(RUSAGE_SELF, &usage);
    std::cout << elapsed_b.count() << ' ' << elapsed_c.count() << ' ' << usage.ru_maxrss << ' ' << c.getDepth() << ' ' << (good != 0) << std::endl;
    elapsed_c = std::chrono::duration<double>::zero();
  }
  
  return 0;
}
