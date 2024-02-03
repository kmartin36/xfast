#include <iostream>
#include <utility>
#include <random>
#include <vector>
#include <limits>
#include <set>
#include <chrono>
#include <ctime>
#include <cmath>
#include <unordered_map>
#include <sparsehash/sparse_hash_map>
#include <sys/resource.h>

// 100000000 3.05715 2.84022 2.79974 2.39814 31930504
// 110000000 3.85684 4.94752 2.8591 2.34912 35699908
// 120000000 3.86588 2.97757 2.85305 2.33472 38602060
// 130000000 3.90349 3.006 2.87935 2.35273 41459860
// 140000000 3.95209 3.0478 2.92428 2.39232 44278060
// 150000000 3.97941 3.07911 2.94388 2.40543 47059036
// 160000000 4.01087 3.10131 3.04504 2.47508 49805692
// 170000000 4.07793 3.15478 3.00763 2.45849 52520132

// Compile with: g++ -O3 -std=c++11 -I. -o vebt vebt.cpp
// Run with: ./vebt 1000000

template <typename key>
struct vebtNext {
  typedef uint8_t halfType;
};

template <>
struct vebtNext<uint16_t> {
  typedef uint8_t halfType;
};

template <>
struct vebtNext<uint32_t> {
  typedef uint16_t halfType;
};

template <>
struct vebtNext<uint64_t> {
  typedef uint32_t halfType;
};

template <typename key>
using vebtNextType = typename vebtNext<key>::halfType;

template <typename key>
class vebt {
public:
  vebt() {
    min = max = std::numeric_limits<key>::max();
    count = 0;
    summary = NULL;
    // printf("Init size %lx\n", (unsigned long)mask+1);
  }
  bool contains(key x) {
    if (!count)
      return false;
    if (x == min || x == max)
      return true;
    if (!subtrees.count(a(x)))
      return false;
    return subtrees[a(x)].contains(b(x));
  };
  key front() {return min;};
  key back() {return max;};
  std::size_t size() {return count;};
  key lower_bound(key x) {
    if (x <= min)
      return min;
    if (x == max)
      return max;
    if (x > max)
      return std::numeric_limits<key>::max();
    if (!summary || summary->back() < a(x))
      return max;
    if (summary->back() == a(x) && subtrees[a(x)].back() < b(x))
      return max;
    if (subtrees.count(a(x)) && subtrees[a(x)].back() >= b(x))
      return (x & ~mask) | subtrees[a(x)].lower_bound(b(x));
    x = (key)summary->lower_bound(a(x)+1) << shift;
    return x | subtrees[a(x)].front();
  }
  key upper_bound(key x) {
    if (x < min)
      return min;
    if (x >= max)
      return std::numeric_limits<key>::max();
    if (!summary || summary->back() < a(x))
      return max;
    if (summary->back() == a(x) && subtrees[a(x)].back() <= b(x))
      return max;
    if (subtrees.count(a(x)) && subtrees[a(x)].back() > b(x))
      return (x & ~mask) | subtrees[a(x)].upper_bound(b(x));
    x = (key)summary->upper_bound(a(x)) << shift;
    return x | subtrees[a(x)].front();
  }
  bool insert(key x) {
    if (count == 0) {
      min = max = x;
      count++;
      // printf("Inserted %016lx into empty\n", (unsigned long)x);
      return true;
    } else if (count == 1) {
      if (min == x)
        return false;
      min = std::min(x, min);
      max = std::max(x, max);
      count++;
      // printf("Inserted %016lx into size(1)\n", (unsigned long)x);
      return true;
    } 
    if (contains(x))
      return false;
    // printf("Inserted %016lx into size %ld\n", (unsigned long)x, count);
    if (x < min)
      std::swap(x, min);
    else if (x > max)
      std::swap(x, max);
    if (summary == NULL)
      summary = new vebt<vebtNextType<key>>();
    summary->insert(a(x));
    subtrees[a(x)].insert(b(x));
    count++;
    return true;
  }
  bool remove(key x) {
    if (count == 0 || !contains(x))
      return false;
    count--;
    if (count == 0)
      min = max = std::numeric_limits<key>::max();
    else if (count == 1 && x == min)
      min = max;
    else if (count == 1)
      max = min;
    if (count > 1) {
      if (x == min)
        x = min = ((key)summary->front() << shift) | subtrees[summary->front()].front();
      else if (x == max)
        x = max = ((key)summary->back() << shift) | subtrees[summary->back()].back();
      subtrees[a(x)].remove(b(x));
      if (subtrees[a(x)].size() == 0) {
        subtrees.erase(a(x));
        summary->remove(a(x));
      }
      if (subtrees.empty()) {
        delete summary;
        summary = NULL;
      }
    }
    return true;
  }

private:
  vebt<vebtNextType<key>> *summary;
  std::unordered_map<vebtNextType<key>, vebt<vebtNextType<key>> > subtrees;
  key min, max;
  std::size_t count;
  static constexpr uint8_t shift = sizeof(key) * 4;
  static constexpr key mask = ((key)1 << shift) - 1;
  vebtNextType<key> a(key x) {return x >> shift;};
  vebtNextType<key> b(key x) {return x & mask;};
};

template <>
class vebt<uint8_t> {
public:
  vebt() {
    min = max = std::numeric_limits<uint8_t>::max();
    count = 0;
    summary = NULL;
    // printf("Init size %lx\n", (unsigned long)mask+1);
  }
  bool contains(uint8_t x) {
    if (!count)
      return false;
    if (x == min || x == max)
      return true;
    if (!subtrees.count(a(x)))
      return false;
    return subtrees[a(x)].count(b(x)) > 0;
  };
  uint8_t front() {return min;};
  uint8_t back() {return max;};
  std::size_t size() {return count;};
  uint8_t lower_bound(uint8_t x) {
    if (x <= min)
      return min;
    if (x == max)
      return max;
    if (x > max)
      return mask;
    if (!summary || *summary->rbegin() < a(x))
      return max;
    if (*summary->rbegin() == a(x) && *subtrees[a(x)].rbegin() < b(x))
      return max;
    if (subtrees.count(a(x)) && *subtrees[a(x)].rbegin() >= b(x))
      return (x & ~mask) | *subtrees[a(x)].lower_bound(b(x));
    x = *summary->lower_bound(a(x)+1) << shift;
    return x | *subtrees[a(x)].begin();
  }
  uint8_t upper_bound(uint8_t x) {
    if (x < min)
      return min;
    if (x >= max)
      return std::numeric_limits<uint8_t>::max();
    if (!summary || *summary->rbegin() < a(x))
      return max;
    if (*summary->rbegin() == a(x) && *subtrees[a(x)].rbegin() <= b(x))
      return max;
    if (subtrees.count(a(x)) && *subtrees[a(x)].rbegin() > b(x))
      return (x & ~mask) | *subtrees[a(x)].upper_bound(b(x));
    x = *summary->upper_bound(a(x)) << shift;
    return x | *subtrees[a(x)].begin();
  }
  bool insert(uint8_t x) {
    if (count == 0) {
      min = max = x;
      count++;
      // printf("Inserted %016lx into empty\n", (unsigned long)x);
      return true;
    } else if (count == 1) {
      if (min == x)
        return false;
      min = std::min(x, min);
      max = std::max(x, max);
      count++;
      // printf("Inserted %016lx into size(1)\n", (unsigned long)x);
      return true;
    } 
    if (contains(x))
      return false;
    // printf("Inserted %016lx into size %ld\n", (unsigned long)x, count);
    if (x < min)
      std::swap(x, min);
    else if (x > max)
      std::swap(x, max);
    if (summary == NULL)
      summary = new std::set<uint8_t>();
    summary->insert(a(x));
    subtrees[a(x)].insert(b(x));
    count++;
    return true;
  }
  bool remove(uint8_t x) {
    if (count == 0 || !contains(x))
      return false;
    count--;
    if (count == 0)
      min = max = std::numeric_limits<uint8_t>::max();
    else if (count == 1 && x == min)
      min = max;
    else if (count == 1)
      max = min;
    if (count > 1) {
      if (x == min)
        x = min = (*summary->begin() << shift) | *subtrees[*summary->begin()].begin();
      else if (x == max)
        x = max = (*summary->rbegin() << shift) | *subtrees[*summary->rbegin()].rbegin();
      subtrees[a(x)].erase(b(x));
      if (subtrees[a(x)].size() == 0) {
        subtrees.erase(a(x));
        summary->erase(a(x));
      }
      if (subtrees.empty()) {
        delete summary;
        summary = NULL;
      }
    }
    return true;
  }
  
private:
  std::set<uint8_t> *summary;
  std::unordered_map<uint8_t, std::set<uint8_t> > subtrees;
  uint8_t min, max;
  std::size_t count;
  static constexpr uint8_t shift = sizeof(uint8_t) * 4;
  static constexpr uint8_t mask = ((uint8_t)1 << shift) - 1;
  uint8_t a(uint8_t x) {return x >> shift;};
  uint8_t b(uint8_t x) {return x & mask;};
};

int main(int argc, char **argv) {
  // std::set<uint64_t> b;
  vebt<uint64_t> c;
//   std::random_device rd;
  std::mt19937_64 gen(0x55aa55aa);
  std::chrono::time_point<std::chrono::high_resolution_clock> start, end;
  std::chrono::duration<double> elapsed_b(0.0), elapsed_c(0.0);
  std::cout << "N rb_insert(us) vebt_insert(us) rb_lower_bound(us) vebt_lower_bound(us) rusage(kB)" << std::endl;
  
  long tcount = 0;
  while (*++argv) {
    long count = atol(*argv) - tcount;
    tcount += count;
    for (long i=0; i<count; i++) {
      uint64_t k = gen();
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
    std::size_t good = 0, bad = 0;
    for (int i=0; i<1000000; i++) {
      uint64_t x = gen(), k0, k1;
      start = std::chrono::high_resolution_clock::now();
      k0 = c.lower_bound(x);
      end = std::chrono::high_resolution_clock::now();
      elapsed_c += end - start;
      good += k0;
      // auto it = b.lower_bound(x);
      // start = std::chrono::high_resolution_clock::now();
      // elapsed_b += start - end;
      // if (it == b.end())
      //   k1 = std::numeric_limits<uint64_t>::max();
      // else
      //   k1 = *it;
      // if (k0 == k1)
      //   good++;
      // else if (++bad < 10)
      //   printf("%016lx %016lx %016lx\n", k0, k1, x);
    }
    struct rusage usage;
    getrusage(RUSAGE_SELF, &usage);
    std::cout << elapsed_b.count() << ' ' << elapsed_c.count() << ' ' << usage.ru_maxrss << ' ' << (good != 0) << std::endl;
    elapsed_c = std::chrono::duration<double>::zero();
  }
  
  return 0;
}
