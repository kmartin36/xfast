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
public:
  class iterator {
public:
    using iterator_category = std::random_access_iterator_tag;
    using difference_type   = std::ptrdiff_t;
    using value_type        = T;
    using pointer           = T*;
    using reference         = T&;
    iterator() : i(0), a(nullptr) {}
    iterator(std::size_t Idx, extVec<T> *Vec) : i(Idx), a(Vec) {}
    iterator(const iterator &it) : i(it.i), a(it.a) {}
    iterator operator=(const iterator &it) {i=it.i; a=it.a; return *this; }
    reference operator*() const { return (*a)[i]; }
    pointer operator->() { return &(*a)[i]; }
    reference operator[](difference_type n) const { return (*a)[i+n]; }
    iterator& operator++() { i++; return *this; }
    iterator operator++(int) { iterator tmp = *this; i++; return tmp; }
    iterator& operator--() { i--; return *this; }
    iterator operator--(int) { iterator tmp = *this; i--; return tmp; }
    iterator& operator+=(difference_type n) { i+=n; return *this; }
    iterator& operator-=(difference_type n) { i-=n; return *this; }
    friend bool operator==(const iterator& a, const iterator& b) { return a.a == b.a && a.i == b.i; }
    friend bool operator!=(const iterator& a, const iterator& b) { return a.a != b.a || a.i != b.i; }
    friend bool operator<(const iterator& a, const iterator& b) { return a.i < b.i; }
    friend bool operator>(const iterator& a, const iterator& b) { return a.i > b.i; }
    friend bool operator<=(const iterator& a, const iterator& b) { return a.i <= b.i; }
    friend bool operator>=(const iterator& a, const iterator& b) { return a.i >= b.i; }
    friend iterator operator+(const iterator& it, difference_type n) { iterator tmp = it; tmp+=n; return tmp; }
    friend iterator operator+(difference_type n, const iterator& it) { return it + n; }
    friend iterator operator-(const iterator& it, difference_type n) { iterator tmp = it; tmp-=n; return tmp; }
    friend iterator operator-(difference_type n, const iterator& it) { return it - n; }
    friend difference_type operator-(const iterator& lhs, const iterator& rhs) { return lhs.i - rhs.i; }
private:
    std::size_t i;
    extVec<T> *a;
  };
  extVec(std::size_t count) {tbl.reserve(8*sizeof(std::size_t) - 12);resize(count);}
  ~extVec() {clear();}
  std::size_t size() {return tbl.empty() ? t0.size() : (4096ULL<<tbl.size());}
  void clear() {resize(0);}
  void resize(std::size_t count) {
    if (count <= 4096)
      t0.resize(count);
    for (std::size_t sz = size(); sz < count; sz <<= 1) {
      tbl.push_back(new T [sz]);
    }
    while (size() > count) {
      delete[] tbl.back();
      tbl.pop_back();
    }
  }
  T& operator[](std::size_t pos) {
    if (pos < 4096)
      return t0[pos];
    uint8_t shift = __builtin_clzll(pos >> 12);
    return tbl.at(63-shift)[pos & (~0ULL >> (shift-11))];
  }
  iterator begin() { return iterator((std::size_t)0, this); }
  iterator end() { return iterator(size(), this); }
  iterator rbegin() {return std::reverse_iterator(end());}
  iterator rend() {return std::reverse_iterator(begin());}
  void swap(extVec &other) {tbl.swap(other.tbl);}
private:
  std::vector<T*> tbl;
  std::vector<T> t0;
};

template <class key>
class ktrie {
private:
  typedef extVec<std::set<key> > vecT;
  // typedef std::vector<std::set<key> > vecT;

public:
  class iterator {
public:
    typedef typename vecT::iterator T0;
    typedef typename std::set<key>::iterator T1;
    using iterator_category = std::bidirectional_iterator_tag;
    using difference_type   = std::ptrdiff_t;
    using value_type        = key;
    using pointer           = key*;
    using reference         = key&;
    iterator() {}
    iterator(T0 iTop, T1 iBot, ktrie &ikt) : i0(iTop), i1(iBot), kt(ikt) {normalize();}
    iterator(const iterator &it) : i0(it.i0), i1(it.i1), kt(it.kt) {normalize();}
    iterator operator=(const iterator &it) {i0=it.i0; i1=it.i1; kt = i0.kt; normalize(); return *this;}
    reference operator*() {k = i0-kt.tbl.end(); k = (k << kt.shift) | *i1; return k;}
    pointer operator->() {return &*this;}
    iterator& operator++() {
      if (i0 == kt.tbl.end())
        return *this;
      if (i1 != i0->end())
        ++i1;
      while (i1==i0->end() && ++i0 != kt.tbl.end())
        i1 = i0->begin();
      return *this;
    }
    iterator operator++(int) {iterator tmp = *this; ++*this; return tmp;}
    iterator& operator--() {
      while (i1 == i0->begin() && i0 != kt.tbl.begin())
        i1 = (--i0)->end();
      if (i1 != i0->begin())
        --i1;
      return *this;
    }
    iterator operator--(int) {iterator tmp = *this; --*this; return tmp;}
    friend bool operator==(const iterator& a, const iterator& b) { return a.i0 == b.i0 && (a.i1 == b.i1 || a.i0 == a.kt.tbl.end()); }
    friend bool operator!=(const iterator& a, const iterator& b) { return !(a == b); }
    T0 i0;
    T1 i1;
    key k;
    ktrie &kt;
    void normalize() {
      while (i0 != kt.tbl.end() && i1 == i0->end() && ++i0 != kt.tbl.end())
        i1 = i0->begin();
    }
  };
  ktrie() : shift(8 * sizeof(key)), mask(0), sz(0), tbl(1) {}
  template<class InputIt>
  ktrie(InputIt first, InputIt last) : shift(8 * sizeof(key)), mask(0), sz(0), tbl(1) {
    using category = typename std::iterator_traits<InputIt>::iterator_category;
    static_assert(std::is_base_of_v<std::input_iterator_tag, category>);
    if constexpr (std::is_base_of_v<std::random_access_iterator_tag, category>) {
      resize(last - first);
      for (auto i = first; i != last; i++)
        tbl[a(*i)].emplace(b(*i));
    } else {
      insert(first, last);
    }
  }
  ktrie(std::initializer_list<key> init) : tbl(1) {
    resize(init.size());
    for (auto &i : init)
      tbl[a(i)].emplace(b(i));
  }
  ~ktrie() {clear();}
  ktrie& operator=(const ktrie& other) {tbl=other.tbl; mask=other.mask; sz=other.sz; shift=other.shift; return *this;}
  ktrie& operator=(const ktrie&& other) {tbl=std::move(other.tbl); mask=other.mask; sz=other.sz; shift=other.shift; return *this;}
  ktrie& operator=(std::initializer_list<key> init) {
    resize(init.size());
    for (auto &i : init)
      tbl[a(i)].emplace(b(i));
    return *this;
  }
  iterator begin() {return iterator(tbl.begin(), tbl[0].begin(), *this);}
  iterator end() {return iterator(tbl.end(), tbl[mask >> shift].end(), *this);}
  iterator rbegin() {return std::reverse_iterator(end());}
  iterator rend() {return std::reverse_iterator(begin());}
  bool empty() const {return sz == 0;}
  std::size_t size() const {return sz;}
  void clear() {tbl.clear(); shift=8 * sizeof(key); mask=sz=0;}
  std::pair<iterator, bool> insert(key x) {
    auto i0 = tbl.begin() + a(x);
    auto i1 = i0->emplace(b(x));
    if (!i1.second)
      return std::make_pair(iterator(i0, i1.first, *this), false);
    if (sz == (mask >> (shift-1))+1) {
      // resize the array
      shift--;
      mask |= ((key)1) << shift;
      tbl.resize(a(mask)+1);
      for (std::size_t i = ((mask>>shift)>>1)+1; --i;) {
        for (auto &j : tbl[i]) {
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
    sz++;
    i0 = tbl.begin() + a(x);
    return std::make_pair(iterator(i0, i0->find(b(x)), *this), true);
  }
  template<class InputIt>
  void insert(InputIt first, InputIt last) {for (InputIt i=first; i!=last; i++) {insert(*i);}}
  void insert(std::initializer_list<key> ilist) {for (auto i : ilist) {insert(i);}}
  iterator erase(iterator pos) {
    key x = *pos;
    pos.i0->erase(pos.i1);
    checkShrink();
    sz--;
    return upper_bound(x);
  }
  std::size_t erase(const key& x) {
    std::size_t s = tbl[a(x)].erase(b(x));
    if (s)
      checkShrink();
    sz -= s;
    return s;
  }
  void swap(ktrie &other) {
    tbl.swap(other.tbl);
    std::swap(mask, other.mask);
    std::swap(sz, other.sz);
    std::swap(shift, other.shift);
  }
  std::size_t count(const key& x) const {return contains(x);}
  iterator find(const key& x) {auto i = tbl.begin() + a(x); return iterator(i, i->find(x), *this);}
  bool contains(const key& x) const {return tbl[a(x)].count(x) > 0;}
  std::pair<iterator, iterator> equal_range(const key& x) {
    auto i=lower_bound(x);
    return std::make_pair(i, *i != x ? i : i+1);
  }
  iterator lower_bound(const key& x) {
    auto i = tbl.begin() + a(x);
    return iterator(i, i->lower_bound(b(x)), *this);
  }
  iterator upper_bound(const key& x) {
    auto i = tbl.begin() + a(x);
    return iterator(i, i->upper_bound(b(x)), *this);
  }
  std::size_t getDepth() {
    auto f = [] (const std::set<key> &i) {return i.size();};
    return std::transform_reduce(tbl.begin(), tbl.end(), 0, std::max<std::size_t>, f);
    return 1;
  }

private:
  vecT tbl;
  key mask;
  std::size_t sz, shift;
  void resize(std::size_t count) {
    sz = count;
    shift = __builtin_clzll(std::max(sz, (std::size_t)1)-1);
    mask = ~0ULL << shift;
    tbl.resize((mask >> shift) + 1);
  }
  void checkShrink() {
    if (mask && sz == -(mask >> 1)) {
      // resize the array
      mask <<= 1;
      shift++;
      for (std::size_t i=1; i<tbl.size(); i++) {
        for (auto j : tbl[i]) {
          key c = (i << (shift-1)) | j;
          tbl[i>>1].emplace_hint(tbl[i>>1].end(), b(c));
        }
        tbl[i].clear();
      }
      tbl.resize(a(mask)+1);
    }
  }
  constexpr key a(key x) {return (x & mask) >> shift;};
  constexpr key b(key x) {return x & (~mask);};
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
      auto it0 = c.lower_bound(x);
      end = high_resolution_clock::now();
      elapsed_c += end - start;
      k0 = it0 == c.end() ? std::numeric_limits<key>::max() : *it0;
      #ifndef doTest
      good &= k0 != 0;
      #else
      auto it1 = b.lower_bound(x);
      start = high_resolution_clock::now();
      elapsed_b += start - end;
      k1 = it1 == b.end() ? std::numeric_limits<key>::max() : *it1;
      good &= k0 == k1;
      if (k0 != k1 && toPrint) {
        printf("%016llx %016llx %016llx\n", k0, k1, x);
        --toPrint;
      }
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
