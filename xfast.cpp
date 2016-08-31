#include <unordered_map>
#include <iostream>
#include <utility>
#include <random>
#include <vector>
#include <limits>
#include <map>
#include <chrono>
#include <ctime>
#include <cmath>
#include <sparsehash/sparse_hash_map>
#include <sys/resource.h>

using std::cout;
using std::endl;

template <class key, class T>
class xfast {
public:
  xfast() {accesses = 0;}
  bool contains(key x) {accesses++; return t0.count(x) > 0;};
  key nearest(key x);
  key lower_bound(key x);
  key upper_bound(key x);
  bool insert(key x, T y);
  bool remove(key x);
  T operator[](key x) {accesses++; return t0[x].second;};
  unsigned long long getOperationsCount() {return accesses;};
  std::size_t getNumElements() {return t0.size();};
  std::size_t getStorageUsed();

private:
  static const int w = 8 * sizeof(key);
  std::unordered_map<key, std::pair<key, key> > tables[w];
  std::unordered_map<key, std::pair<std::pair<key, key>, T> > t0;
  std::pair<std::pair<key, key>, int> inner(key x);
  key closer(key x, std::pair<key, key> u);
  key near(key x);
  unsigned long long accesses;
};

template <class key, class T>
std::pair<std::pair<key, key>, int> xfast<key, T>::inner(key x) {
  if (t0.empty())
    throw;
  int l = 0, h = w-1;
  while (h-l > 1) {
    int i = (l+h) / 2;
    key q = (x >> (w-1-i)) >> 1;
    accesses++;
    if (tables[i].count(q))
      l = i;
    else
      h = i;
  }
  accesses++;
  return std::make_pair(tables[l][(x >> (w-1-l)) >> 1], l);
}

template <class key, class T>
key xfast<key, T>::closer(key x, std::pair<key, key> u) {
  if (x < u.first)
    return u.first;
  return u.second;
}

template <class key, class T>
key xfast<key, T>::near(key x) {
  if (t0.empty())
    throw;
  return closer(x, inner(x).first);
}

template <class key, class T>
key xfast<key, T>::nearest(key x) {
  if (getNumElements() == 0)
    return std::numeric_limits<key>::max();
  key k = near(x);
  if (!contains(k))
    return std::numeric_limits<key>::max();
  accesses++;
  std::pair<key, key> u = t0[k].first;
  if (k < x)
    if (x - k > u.second - x)
      return u.second;
    else
      return k;
  if (k - x > x - u.first)
    return u.first;
  return k;
}

template <class key, class T>
key xfast<key, T>::lower_bound(key x) {
  if (getNumElements() == 0)
    return std::numeric_limits<key>::max();
  key k = near(x);
  if (k < x) {
    accesses++;
    k = t0[k].first.second;
  }
  return k;
}

template <class key, class T>
key xfast<key, T>::upper_bound(key x) {
  if (getNumElements() == 0)
    return std::numeric_limits<key>::max();
  key k = near(x);
  if (k <= x) {
    accesses++;
    k = t0[k].first.second;
  }
  return k;
}

template <class key, class T>
bool xfast<key, T>::insert(key x, T y) {
  if (t0.empty()) {
    accesses++;
    t0.insert({x, std::make_pair(std::make_pair(x, x), y)});
    for (int i=0; i<w; i++) {
      accesses++;
      tables[i].insert({(x >> (w-1-i)) >> 1, std::make_pair(x, x)});
    }
    return true;
  }
  if (contains(x))
    return false;
  
  key kl = near(x);
  std::pair<key, key> u = t0[kl].first;
  key kh = (kl < x) ? u.second : u.first;
  if (kl > x)
    std::swap(kl, kh);
  t0[kl].first.second = x;
  accesses++;
  t0[kh].first.first = x;
  accesses++;
  t0.insert({x, std::make_pair(std::make_pair(kl, kh), y)});
  
  for (int i=w; i--; ) {
    accesses++;
    if (tables[i].count((x >> (w-1-i)) >> 1)) {
      if (x < tables[i][(x >> (w-1-i)) >> 1].first)
        tables[i][(x >> (w-1-i)) >> 1].first = x;
      else if (x > tables[i][(x >> (w-1-i)) >> 1].second)
        tables[i][(x >> (w-1-i)) >> 1].second = x;
      else
        break;
    } else {
      tables[i].insert({(x >> (w-1-i)) >> 1, std::make_pair(x, x)});
    }
  }
  
  return true;
}

template <class key, class T>
bool xfast<key, T>::remove(key x) {
  if (t0.empty() || !contains(x))
    return false;
  
  std::pair<key, key> u = t0[x].first;
  accesses++;
  t0[u.first].first.second = u.second;
  accesses++;
  t0[u.second].first.first = u.first;
  t0.erase(x);
  
  key x1 = x;
  for (int i=w-1; i>=0; i--) {
    accesses++;
    x1 >>= 1;
    std::pair<key, key> v = tables[i][x1];
    if (v.first == v.second)
      tables[i].erase(x1);
    else if (v.first == x)
      tables[i][x1].first = u.second;
    else if (v.second == x)
      tables[i][x1].second = u.first;
  }
  
  return true;
}

template <class key, class T>
std::size_t xfast<key, T>::getStorageUsed() {
  std::size_t s = getNumElements();
  for (auto t : tables)
    s += t.size();
  return s;
}

template <class key, class T>
class kfast {
public:
  kfast() {accesses = 0;}
  bool contains(key x) {accesses++; return t0.count(x) > 0;};
  key nearest(key x);
  key lower_bound(key x);
  key upper_bound(key x);
  bool insert(key x, T y);
  bool remove(key x);
  T operator[](key x) {accesses++; return t0[x].second;};
  unsigned long long getOperationsCount() {return accesses;};
  std::size_t getNumElements() {return t0.size();};
  std::size_t getStorageUsed();
  std::size_t getDepth() {return tables.size();};
  int clz(key x);

private:
  std::vector<google::sparse_hash_map<key, std::pair<key, key> > > tables;
  google::sparse_hash_map<key, std::pair<std::pair<key, key>, T> > t0;
  std::pair<std::pair<key, key>, int> inner(key x);
  key closer(key x, std::pair<key, key> u);
  key near(key x);
  unsigned long long accesses;
};

template <class key, class T>
std::pair<std::pair<key, key>, int> kfast<key, T>::inner(key x) {
  if (t0.empty())
    throw;
  int l = 0, h = tables.size()-1;
  const int w = 8 * sizeof(key);
  while (h-l > 1) {
    int i = (l+h) / 2;
    key q = (x >> (w-1-i)) >> 1;
    accesses++;
    if (tables[i].count(q))
      l = i;
    else
      h = i;
  }
  key q = (x >> (w-1-h)) >> 1;
  accesses++;
  if (tables[h].count(q))
      l = h;
  accesses++;
  return std::make_pair(tables[l][(x >> (w-1-l)) >> 1], l);
}

template <class key, class T>
key kfast<key, T>::closer(key x, std::pair<key, key> u) {
  int i = clz(x ^ u.first);
  int j = clz(x ^ u.second);
  if (i > j)
    return u.first;
  if (i < j)
    return u.second;
  if (x < u.first)
    return u.first;
  return u.second;
}

template <class key, class T>
key kfast<key, T>::near(key x) {
  if (t0.empty())
    throw;
  return closer(x, inner(x).first);
}

template <class key, class T>
int kfast<key, T>::clz(key x) {
  if (x == 0)
    return 8*sizeof(key);
  int i = 4*sizeof(key);
  int j = 8*sizeof(key)-1;
  while (i) {
    if (x >> i) {
      x >>= i;
      j -= i;
    }
    i >>= 1;
  }
  return j;
}

template <class key, class T>
key kfast<key, T>::nearest(key x) {
  if (getNumElements() == 0)
    return std::numeric_limits<key>::max();
  key k = near(x);
  accesses++;
  std::pair<key, key> u = t0[k].first;
  if (k < x)
    if (x - k > u.second - x)
      return u.second;
    else
      return k;
  if (k - x > x - u.first)
    return u.first;
  return k;
}

template <class key, class T>
key kfast<key, T>::lower_bound(key x) {
  if (getNumElements() == 0)
    return std::numeric_limits<key>::max();
  key k = near(x);
  if (k < x) {
    accesses++;
    k = t0[k].first.second;
  }
  return k;
}

template <class key, class T>
key kfast<key, T>::upper_bound(key x) {
  if (getNumElements() == 0)
    return std::numeric_limits<key>::max();
  key k = near(x);
  if (k <= x) {
    accesses++;
    k = t0[k].first.second;
  }
  return k;
}

template <class key, class T>
bool kfast<key, T>::insert(key x, T y) {
  if (t0.empty()) {
    accesses++;
    t0.insert({x, std::make_pair(std::make_pair(x, x), y)});
    tables.resize(tables.size() + 1);
    tables.back()[0] = std::make_pair(x, x);
    return true;
  }
  if (contains(x))
    return false;
  
  auto in = inner(x);
  key k = closer(x, in.first);
  int i = in.second;
  int j = clz(x ^ k);
  i = std::min((int)tables.size()-1, j);
  
  std::pair<key, key> u = t0[k].first;
  key kh = (k < x) ? u.second : u.first;
  key kl = k;
  if (kl > x)
    std::swap(kl, kh);
  t0[kl].first.second = x;
  accesses++;
  t0[kh].first.first = x;
  accesses++;
  t0.insert({x, std::make_pair(std::make_pair(kl, kh), y)});
  
  const int w = 8 * sizeof(key);
  key x1 = x >> (w-1-i);
  std::pair<key, key> v(std::min(x, k), std::max(x, k));
  for (; i>=0; i--) {
    accesses++;
    x1 >>= 1;
    if (tables[i].count(x1)) {
      if (x < tables[i][x1].first)
        tables[i][x1].first = x;
      else if (x > tables[i][x1].second)
        tables[i][x1].second = x;
      else
        break;
    } else {
      tables[i].insert({x1, v});
    }
  }
  for (i=tables.size(); i<=j; i++) {
    tables.resize(tables.size() + 1);
    tables.back()[(x >> (w-1-i)) >> 1] = v;
  }
  
  return true;
}

template <class key, class T>
bool kfast<key, T>::remove(key x) {
  if (t0.empty() || !contains(x))
    return false;
  
  if (t0.size() == 1) {
    accesses++;
    t0.clear();
    tables.clear();
    return true;
  }
  
  std::pair<key, key> u = t0[x].first;
  accesses++;
  t0[u.first].first.second = u.second;
  accesses++;
  t0[u.second].first.first = u.first;
  t0.erase(x);
  
  key k = closer(x, u);
  int i = clz(x ^ k);
  const int w = 8 * sizeof(key);
  key x1 = x >> (w-1-i);
  for (; i>=0; i--) {
    accesses++;
    x1 >>= 1;
    std::pair<key, key> v = tables[i][x1];
    if ((v.first == x && v.second == u.second) || (v.second == x && v.first == u.first))
      tables[i].erase(x1);
    else if (v.first == x)
      tables[i][x1].first = u.second;
    else if (v.second == x)
      tables[i][x1].second = u.first;
    else
      break;
  }
  
  return true;
}

template <class key, class T>
std::size_t kfast<key, T>::getStorageUsed() {
  std::size_t s = getNumElements();
  for (auto t : tables)
    s += t.size();
  return s;
}

int main(int argc, char **argv) {
  std::map<unsigned long long, unsigned long long> b;
  kfast<unsigned long long, unsigned long long> c;
  std::random_device rd;
  std::mt19937_64 gen(rd());
//   unsigned long long k = gen();
  std::chrono::time_point<std::chrono::high_resolution_clock> start, end;
  std::chrono::duration<double> elapsed_b(0.0), elapsed_c(0.0);
  
  int count = atoi(argv[1]);
  for (int i=0; i<count; i++) {
//     unsigned long long k = gen(), v = i;
    unsigned long long k = gen();
//     k = k * 2862933555777941757UL + 3037000493UL;
    start = std::chrono::high_resolution_clock::now();
    b.emplace(k, k);
    end = std::chrono::high_resolution_clock::now();
    elapsed_b += end - start;
    c.insert(k, k);
    start = std::chrono::high_resolution_clock::now();
    elapsed_c += start - end;
  }
//   std::cout << c.getNumElements() << ' ' << c.getStorageUsed() << ' ' << c.getOperationsCount() << ' ' << elapsed_b.count() << ' ' << elapsed_c.count() << std::endl;
  std::cout << count << ' ' << elapsed_b.count()*1000000/count << ' ' << elapsed_c.count()*1000000/count << ' ';
  
  elapsed_b = std::chrono::duration<double>::zero();
  elapsed_c = std::chrono::duration<double>::zero();
  unsigned long long ops = c.getOperationsCount();
  unsigned long long good = 0;
  for (int i=0; i<1000000; i++) {
    unsigned long long x = gen(), k0, k1;
//     unsigned long long x = k, k0, k1;
//     k = k * 2862933555777941757UL + 3037000493UL;
    start = std::chrono::high_resolution_clock::now();
    k0 = c.lower_bound(x);
    end = std::chrono::high_resolution_clock::now();
    elapsed_c += end - start;
    auto it = b.lower_bound(x);
    start = std::chrono::high_resolution_clock::now();
    elapsed_b += start - end;
    if (it == b.end())
      it = b.begin();
    k1 = it->first;
    if (k0 == k1)
      good++;
    else
      printf("%016llx %016llx %016llx\n", k0, k1, x);
  }
//   std::cout << good << ' ' << c.getNumElements() << ' ' << c.getOperationsCount()-ops << ' ' << c.getDepth() << ' ' << elapsed_b.count() << ' ' << elapsed_c.count() << std::endl;
  struct rusage usage;
  getrusage(RUSAGE_SELF, &usage);
  std::cout << elapsed_b.count() << ' ' << elapsed_c.count() << ' ' << log(usage.ru_maxrss*1000.0/128)/10.0 << std::endl;
  
//   elapsed_b = std::chrono::duration<double>::zero();
//   elapsed_c = std::chrono::duration<double>::zero();
//   ops = c.getOperationsCount();
//   for (int i=0; i<count; i++) {
//     unsigned long long x = gen();
//     auto it = b.lower_bound(x);
//     if (it == b.end())
//       it = b.begin();
//     start = std::chrono::high_resolution_clock::now();
//     c.remove(it->first);
//     end = std::chrono::high_resolution_clock::now();
//     elapsed_c += end - start;
//     b.erase(it->first);
//     start = std::chrono::high_resolution_clock::now();
//     elapsed_b += start - end;
//   }
//   std::cout << c.getNumElements() << ' ' << c.getStorageUsed() << ' ' << c.getOperationsCount()-ops << ' ' << elapsed_b.count() << ' ' << elapsed_c.count() << std::endl;
//   std::cout << elapsed_b.count()/count << ' ' << elapsed_c.count()/count << std::endl;
  
  return 0;
}
