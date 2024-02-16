#include <iostream>
#include <utility>
#include <random>
#include <vector>
#include <limits>
#include <set>
#include <chrono>
#include <ctime>
#include <cmath>
#include <iterator>
#include <sys/resource.h>

// Compile with: g++ -O3 -std=c++17 -o ktrie ktrie.cpp
// Run with: ./ktrie 1000000

using std::cout;
using std::endl;

template <class T>
class extVec {
  class extIter {
public:
    using iterator_category = std::random_access_iterator_tag;
    using difference_type   = std::ptrdiff_t;
    using value_type        = T;
    using pointer           = T*;
    using reference         = T&;
    extIter() : i(0), a(nullptr) {}
    extIter(std::size_t Idx, extVec<T> *Vec) : i(Idx), a(Vec) {}
    extIter(const extIter &it) : i(it.i), a(it,a) {}
    extIter operator=(const extIter &it) {i=it.i; a=it.a;}
    reference operator*() const { return (*a)[i]; }
    pointer operator->() { return &(*a)[i]; }
    reference operator[](difference_type n) const { return (*a)[i+n]; }
    extIter& operator++() { i++; return *this; }
    extIter operator++(int) { extIter tmp = *this; i++; return tmp; }
    extIter& operator--() { i--; return *this; }
    extIter operator--(int) { extIter tmp = *this; i--; return tmp; }
    extIter& operator+=(difference_type n) { i+=n; return *this; }
    extIter& operator-=(difference_type n) { i-=n; return *this; }
    friend bool operator==(const extIter& a, const extIter& b) { return a.a == b.a && a.i == b.i; }
    friend bool operator!=(const extIter& a, const extIter& b) { return a.a != b.a || a.i != b.i; }
    friend bool operator<(const extIter& a, const extIter& b) { return a.i < b.i; }
    friend bool operator>(const extIter& a, const extIter& b) { return a.i > b.i; }
    friend bool operator<=(const extIter& a, const extIter& b) { return a.i <= b.i; }
    friend bool operator>=(const extIter& a, const extIter& b) { return a.i >= b.i; }
    friend extIter operator+(const extIter& it, difference_type n) { extIter tmp = it; tmp+=n; return tmp; }
    friend extIter operator+(difference_type n, const extIter& it) { return it + n; }
    friend extIter operator-(const extIter& it, difference_type n) { extIter tmp = it; tmp-=n; return tmp; }
    friend extIter operator-(difference_type n, const extIter& it) { return it - n; }
    friend difference_type operator-(const extIter& lhs, const extIter& rhs) { return lhs.i - rhs.i; }
private:
    std::size_t i;
    extVec<T> *a;
  };
public:
  extVec(std::size_t count) {tbl.reserve(8*sizeof(std::size_t) + 1);resize(count);}
  ~extVec() {tbl.clear();}
  std::size_t size() {return tbl.empty() ? 0 : a();}
  void resize(std::size_t count) {
    if (!count) {
      tbl.clear();
      return;
    }
    if (tbl.empty())
      tbl.push_back(new T [1]);
    for (std::size_t sz = a(); sz < count; sz <<= 1) {
      tbl.push_back(new T [sz]);
    }
    for (std::size_t sz = a(); sz > count; sz >>= 1) {
      delete[] tbl.back();
      tbl.pop_back();
    }
  }
  T& operator[](std::size_t pos) {
    if (pos == 0)
      return tbl[0][0];
    uint8_t shift = __builtin_clzll(pos);
    return tbl[64-shift][pos & ((~0ULL >> shift) >> 1)];
  }
  extIter begin() { return extIter((std::size_t)0, this); }
  extIter end() { return extIter(size(), this); }
private:
  std::vector<T*> tbl;
  constexpr std::size_t a() {return ((std::size_t)1)<<(tbl.size()-1);};
};

template <class key>
class ktrie {
public:
  ktrie() : shift(8 * sizeof(key)), mask(0), count(0), tbl(1) {}
  constexpr bool contains(key x) {return tbl[a(x)].count() > 0;}
  constexpr key lower_bound(key x) {
    auto &t = tbl[a(x)];
    auto i = t.lower_bound(b(x));
    return inner_bound(mask & x, t, i);
  }
  constexpr key upper_bound(key x) {
    auto &t = tbl[a(x)];
    auto i = t.upper_bound(b(x));
    return inner_bound(mask & x, t, i);
  }
  bool insert(key x) {
    if (!tbl[a(x)].emplace(b(x)).second)
      return false;
    if (count == (mask >> (shift-1))+1) {
      // resize the array
      shift--;
      mask |= ((key)1) << shift;
      tbl.resize(a(mask)+1);
      for (std::size_t i = ((mask>>shift)>>1)+1; --i;) {
        for (auto j : tbl[i]) {
          key c = (i << (shift+1)) | j;
          tbl[a(c)].emplace_hint(tbl[a(c)].end(), b(c));
        }
        tbl[i].clear();
      }
      auto i = tbl[0].lower_bound(-mask);
      for (auto j=i; j!=tbl[0].end(); j++)
        tbl[1].emplace_hint(tbl[1].end(), b(*j));
      tbl[0].erase(i, tbl[0].end());
    }
    count++;
    return true;
  }
  bool remove(key x) {
    std::size_t s = tbl[a(x)].erase(b(x));
    if (s == 0)
      return false;
    if (mask && count == -(mask >> 1)) {
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
    auto f = [] (const std::set<key> &i) {return i.size();};
    return std::transform_reduce(tbl.begin(), tbl.end(), 0, std::max<std::size_t>, f);
    return 1;
  }

private:
  extVec<std::set<key> > tbl;
  // std::vector<std::set<key> > tbl;
  key mask;
  std::size_t count, shift;
  constexpr key a(key x) {return (x & mask) >> shift;};
  constexpr key b(key x) {return x & (~mask);};
  typedef typename std::set<key>::const_iterator set_iter;
  constexpr key inner_bound(key x, const std::set<key> &t, set_iter &i) {
    if (i != t.end())
      return *i | x;
    for (x-=mask; x; x-=mask)
      if (!tbl[x>>shift].empty())
        return *tbl[x>>shift].begin() | x;
    return std::numeric_limits<key>::max();
  }
};

#undef doTest
#define doTest

int main(int argc, char **argv) {
  using namespace std::chrono_literals;
  using std::chrono::high_resolution_clock;
  typedef unsigned long long key;
  #ifdef doTest
  std::set<key> b;
  #endif
  ktrie<key> c;
  // std::random_device rd;
  std::mt19937_64 gen(55);
  high_resolution_clock::time_point start = high_resolution_clock::now(), end;
  std::chrono::duration<double> elapsed_b(0.0), elapsed_c(0.0);
  cout << "N rb_insert(us) kt_insert(us) rb_lower_bound(us) kt_lower_bound(us) "
          "rusage(kB) tree_depth all_good" << endl;
  
  long tcount = 0;
  while (*++argv) {
    long count = atol(*argv) - tcount;
    tcount += count;
    for (long i=0; i<count; i++) {
      key k = gen();
      #ifdef doTest
      start = high_resolution_clock::now();
      b.insert(k);
      #endif
      end = high_resolution_clock::now();
      elapsed_b += end - start;
      c.insert(k);
      start = high_resolution_clock::now();
      elapsed_c += start - end;
    }
    double scale = 1000000.0 / count;
    cout << tcount << ' ' << elapsed_b.count()*scale << ' ' << elapsed_c.count()*scale << ' ';
    elapsed_b = elapsed_c = 0s;
    bool good = true;
    int toPrint = 10;
    for (int i=0; i<1000000; i++) {
      key x = gen(), k0, k1;
      start = high_resolution_clock::now();
      k0 = c.lower_bound(x);
      end = high_resolution_clock::now();
      elapsed_c += end - start;
      #ifndef doTest
      good &= k0 != 0;
      #else
      auto it = b.lower_bound(x);
      start = high_resolution_clock::now();
      elapsed_b += start - end;
      k1 = it == b.end() ? std::numeric_limits<key>::max() : *it;
      good &= k0 == k1;
      if (k0 != k1 && toPrint--)
        printf("%016llx %016llx %016llx\n", k0, k1, x);
      #endif
    }
    struct rusage usage;
    getrusage(RUSAGE_SELF, &usage);
    cout << elapsed_b.count() << ' ' << elapsed_c.count() << ' ' << usage.ru_maxrss << ' '
         << c.getDepth() << ' ' << good << endl;
    elapsed_b = elapsed_c = 0s;
  }
  
  return 0;
}
