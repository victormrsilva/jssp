#ifndef HASH_H
#define HASH_H
#include "Instance.hpp"
#include "lp.hpp"
#include <unordered_set>
#include <vector>
extern "C"{
#include "cgraph/cgraph.h"
#include "cgraph/clique_separation.h"
}

    struct S {
        int i;
        int j;
        int t;
        int var;

        S(){

        }

        S(int i1, int j1, int t1, int var1){
          i = i1;
          j = j1;
          t = t1;
          var = var1;
        }

        S& operator =(const S& a)
        {
            i = a.i;
            j = a.j;
            t = a.t;
            var = a.var;
            return *this;
        }

        inline bool operator==(const S a) const {
            return a.var==this->var;
        }

        bool operator<(const S& rhs) const { return var < rhs.var; }
    };

namespace std {
    template <> struct hash<std::vector<int>> {
      size_t operator()(const vector<int>& v) const {
        std::hash<int> hasher;
        std::size_t seed = 0;
        for (int i : v) {
            seed ^= hasher(i) + 0x9e3779b9 + (seed<<6) + (seed>>2);
        }
        return seed;
      }
    };
    template <> struct hash<S> {
      size_t operator()(S v) const {
        return v.var + (v.i<<6);
      }
    };
    template <>
    struct hash<std::vector<S>> {
      size_t operator()(const vector<S>& v) const {
        std::size_t seed = 0;
        for (S i : v) {
          seed ^= i.var + 0x9e3779b9 + (seed<<6) + (seed>>2);
        }
        return seed;
      }
    };
};
#endif