/*
 *  SPDX-License-Identifier: MIT
 */

#include "../parameters.h"
#include "../macros.h"

#include <algorithm>
#include <array>
#include <cctype>
#include <cstdint>
#include <cstdlib>
#include <filesystem>
#include <functional>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <numeric>
#include <optional>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

using namespace std;

namespace {
  class big_int {
    vector<uint64_t> words_;

    void normalize() {
      while (!words_.empty() && words_.back() == 0) {
        words_.pop_back();
      }
    }

  public:
    big_int() {}

    big_int(uint64_t value) {
      if (value) {
        words_.push_back(value);
      }
    }

    bool is_zero() const {
      return words_.empty();
    }

    uint64_t word_at(size_t i) const {
      return i < words_.size() ? words_[i] : 0;
    }

    bool test(unsigned int i) const {
      return (word_at(i / 64) >> (i % 64)) & 1;
    }

    void set(unsigned int i) {
      words_.resize(max(words_.size(), size_t(i / 64 + 1)), 0);
      words_[i / 64] |= uint64_t{1} << (i % 64);
    }

    void flip(unsigned int i) {
      words_.resize(max(words_.size(), size_t(i / 64 + 1)), 0);
      words_[i / 64] ^= uint64_t{1} << (i % 64);
      normalize();
    }

    unsigned int lsb() const {
      if (is_zero()) {
        throw invalid_argument("lsb of zero");
      }

      size_t i = 0;
      while (!words_[i]) {
        ++i;
      }

#if __has_builtin(__builtin_ctzll)
      auto b = __builtin_ctzll(words_[i]);
#else
      auto w         = words_[i];
      unsigned int b = 0;
      while (!(w & 1)) {
        w >>= 1;
        ++b;
      }
#endif

      return 64 * i + b;
    }

    unsigned int msb() const {
      if (is_zero()) {
        throw invalid_argument("msb of zero");
      }

#if __has_builtin(__builtin_clzll)
      auto b = 63 - __builtin_clzll(words_.back());
#else
      auto w         = words_.back();
      unsigned int b = 0;
      while (w >>= 1) {
        ++b;
      }
#endif

      return 64 * (words_.size() - 1) + b;
    }

    big_int& operator^=(const big_int& other) {
      words_.resize(max(words_.size(), other.words_.size()), 0);
      for (size_t i = 0; i < other.words_.size(); ++i) {
        words_[i] ^= other.words_[i];
      }
      normalize();
      return *this;
    }

    big_int& operator|=(const big_int& other) {
      words_.resize(max(words_.size(), other.words_.size()), 0);
      for (size_t i = 0; i < other.words_.size(); ++i) {
        words_[i] |= other.words_[i];
      }
      return *this;
    }

    big_int& operator&=(const big_int& other) {
      words_.resize(min(words_.size(), other.words_.size()));
      for (size_t i = 0; i < words_.size(); ++i) {
        words_[i] &= other.words_[i];
      }
      normalize();
      return *this;
    }

    friend big_int operator<<(const big_int& a, unsigned int shift) {
      if (a.is_zero()) {
        return 0;
      }
      const auto offset = shift / 64;
      const auto bits   = shift % 64;
      big_int result;
      result.words_.resize(a.words_.size() + offset + (bits != 0), 0);
      for (size_t i = 0; i < a.words_.size(); ++i) {
        result.words_[i + offset] |= a.words_[i] << bits;
        if (bits) {
          result.words_[i + offset + 1] |= a.words_[i] >> (64 - bits);
        }
      }
      result.normalize();
      return result;
    }

    big_int& operator<<=(unsigned int shift) {
      return *this = *this << shift;
    }

    // Polynomial multiplication only needs a single-bit right shift.
    void shift_right_one() {
      for (size_t i = 0; i < words_.size(); ++i) {
        words_[i] = (words_[i] >> 1) | (word_at(i + 1) << 63);
      }
      normalize();
    }

    big_int& operator++() {
      for (auto& w : words_) {
        if (++w) {
          return *this;
        }
      }
      words_.push_back(1);
      return *this;
    }

    friend bool operator==(const big_int& a, const big_int& b) {
      return a.words_ == b.words_;
    }

    friend bool operator<(const big_int& a, const big_int& b) {
      if (a.words_.size() != b.words_.size()) {
        return a.words_.size() < b.words_.size();
      }
      return lexicographical_compare(a.words_.rbegin(), a.words_.rend(), b.words_.rbegin(),
                                     b.words_.rend());
    }
  };

  big_int operator^(big_int a, const big_int& b) {
    return a ^= b;
  }

  big_int operator|(big_int a, const big_int& b) {
    return a |= b;
  }

  big_int operator&(big_int a, const big_int& b) {
    return a &= b;
  }

  bool operator!=(const big_int& a, const big_int& b) {
    return !(a == b);
  }

  [[noreturn]] void fail(const string& msg) {
    throw std::runtime_error(msg);
  }

  void check(bool cond, const string& msg) {
    if (!cond) {
      fail(msg);
    }
  }

  template <typename R, typename... Args>
  class cached_eval {
    map<tuple<Args...>, R> cache_;
    function<R(Args...)> f_;

  public:
    cached_eval(function<R(Args...)> f) : f_{f} {}

    const R& operator()(Args... args) {
      const auto key = make_tuple(args...);
      auto it        = cache_.find(key);
      if (it != cache_.end()) {
        return it->second;
      }

      auto result = f_(args...);
      cache_[key] = result;
      return cache_[key];
    }
  };

  big_int bit(unsigned int i) {
    big_int ret{0};
    ret.set(i);
    return ret;
  }

  template <typename Fn>
  void for_each_bit(const big_int& x, Fn fn) {
    if (x.is_zero()) {
      return;
    }

    const auto msb = x.msb();
    for (unsigned int i = x.lsb(); i <= msb; ++i) {
      if (x.test(i)) {
        fn(i);
      }
    }
  }

  optional<unsigned int> pdeg(const big_int& a) {
    if (a.is_zero()) {
      return nullopt;
    }
    return a.msb();
  }

  big_int pmul(big_int a, big_int b) {
    big_int r = 0;
    while (!b.is_zero()) {
      if (b.test(0)) {
        r ^= a;
      }
      b.shift_right_one();
      a <<= 1;
    }
    return r;
  }

  big_int pmod(big_int a, const big_int& m) {
    const auto d = pdeg(m);
    check(d.has_value(), "polynomial modulus must be nonzero");
    while (!a.is_zero() && *pdeg(a) >= *d) {
      a ^= m << (*pdeg(a) - *d);
    }
    return a;
  }

  big_int pmulmod(const big_int& a, const big_int& b, const big_int& m) {
    return pmod(pmul(a, b), m);
  }

  big_int pgcd(big_int a, big_int b) {
    while (!b.is_zero()) {
      const big_int r = pmod(a, b);
      a               = b;
      b               = r;
    }
    return a;
  }

  pair<big_int, big_int> pdivmod(big_int a, const big_int& b) {
    big_int q     = 0;
    const auto db = pdeg(b);
    check(db.has_value(), "polynomial divisor must be nonzero");
    while (!a.is_zero() && *pdeg(a) >= *db) {
      const int s = *pdeg(a) - *db;
      q.set(s);
      a ^= b << s;
    }
    return {q, a};
  }

  big_int pinvmod(const big_int& a, const big_int& m) {
    big_int r0{m};
    big_int r1{pmod(a, m)};
    big_int s0{0};
    big_int s1{1};
    while (!r1.is_zero()) {
      const auto [q, r] = pdivmod(r0, r1);
      r0                = r1;
      r1                = r;
      const big_int ns  = s0 ^ pmul(q, s1);
      s0                = s1;
      s1                = ns;
    }
    check(r0 == 1, "pinvmod: polynomials are not coprime");
    return pmod(s0, m);
  }

  big_int poly_pow(const big_int& p, unsigned int e) {
    big_int r = 1;
    for (unsigned int i = 0; i < e; ++i) {
      r = pmul(r, p);
    }
    return r;
  }

  big_int mtree_poly(const vector<big_int>& tree_moduli) {
    return accumulate(tree_moduli.begin(), tree_moduli.end(), big_int{1},
                      [](auto left, auto right) { return pmul(left, right); });
  }

  vector<big_int> crt_lift_cols(const vector<big_int>& tree_moduli) {
    big_int mtree = mtree_poly(tree_moduli);

    vector<big_int> cols;
    for (const big_int& m : tree_moduli) {
      const auto [q, rem] = pdivmod(mtree, m);
      check(rem.is_zero(), "tree modulus does not divide M_tree");

      const big_int e = pmulmod(q, pinvmod(q, m), mtree);
      const auto deg  = pdeg(m).value_or(0);
      for (unsigned int b = 0; b < deg; ++b) {
        cols.emplace_back(pmulmod(e, bit(b), mtree));
      }
    }
    return cols;
  }

  big_int x_pow_2k(const big_int& f, unsigned int k) {
    big_int r = pmod(2, f);
    for (unsigned int i = 0; i < k; ++i) {
      r = pmulmod(r, r, f);
    }
    return r;
  }

  set<unsigned int> prime_factors(unsigned int n) {
    set<unsigned int> fs;
    for (unsigned int d = 2; d * d <= n; ++d) {
      while (!(n % d)) {
        fs.insert(d);
        n /= d;
      }
    }
    if (n > 1) {
      fs.insert(n);
    }
    return fs;
  }

  bool is_irreducible(const big_int& f) {
    const auto n = pdeg(f).value_or(0);
    if (n <= 0) {
      return false;
    }
    if (n == 1) {
      return true;
    }
    if (!f.test(0)) {
      return false;
    }
    if (x_pow_2k(f, n) != pmod(2, f)) {
      return false;
    }

    for (auto q : prime_factors(n)) {
      const big_int h = pmod(x_pow_2k(f, n / q) ^ 2, f);
      if (pgcd(f, h) != 1) {
        return false;
      }
    }
    return true;
  }

  vector<big_int> irreducibles_of_degree_impl(unsigned int d) {
    vector<big_int> out;
    for (big_int f = bit(d) | 1; !f.test(d + 1); ++f) {
      if (is_irreducible(f)) {
        out.push_back(f);
      }
    }
    return out;
  }

  cached_eval<vector<big_int>, unsigned int> irreducibles_of_degree{irreducibles_of_degree_impl};

  big_int faest_modulus(unsigned int lambda) {
    big_int modulus{1 | (1 << 2)};
    switch (lambda) {
    case 128: {
      modulus.set(128);
      modulus.set(7);
      modulus.set(1);
      break;
    }
    case 192: {
      modulus.set(192);
      modulus.set(7);
      modulus.set(1);
      break;
    }
    case 256: {
      modulus.set(256);
      modulus.set(10);
      modulus.set(5);
      break;
    }
    default:
      fail("unsupported FAEST modulus");
    }
    return modulus;
  }

  big_int combine(const big_int& mask, const vector<big_int>& rows) {
    big_int out = 0;
    for_each_bit(mask, [&](unsigned long i) {
      check(i < rows.size(), "row-combination mask exceeds row count");
      out ^= rows[i];
    });
    return out;
  }

  unsigned int parity(const big_int& x) {
    bool p = false;
    for_each_bit(x, [&](unsigned long) { p = !p; });
    return p ? 1 : 0;
  }

  vector<big_int> reduction_rows(const big_int& q, unsigned int ncols) {
    const auto d = pdeg(q);
    vector<big_int> rows(d.value_or(0), 0);
    big_int cur = 1;
    for (unsigned int j = 0; j < ncols; ++j) {
      for_each_bit(cur, [&rows, j](unsigned long i) {
        check(i < rows.size(), "reduction row index out of range");
        rows[i].set(j);
      });
      cur = pmod(cur << 1, q);
    }
    return rows;
  }

  vector<big_int> pascal_rows(unsigned int e) {
    vector<big_int> rows;
    rows.reserve(e);
    for (unsigned int i = 0; i < e; ++i) {
      big_int row = 0;
      for (unsigned int j = i; j < e; ++j) {
        if ((i & j) == i) {
          row.set(j);
        }
      }
      rows.emplace_back(row);
    }
    return rows;
  }

  vector<big_int> left_inverse(const vector<big_int>& rows, unsigned int ncols) {
    const unsigned int m = rows.size();
    vector<pair<big_int, big_int>> work;
    work.reserve(rows.size());
    for (unsigned int i = 0; i < m; ++i) {
      work.push_back({rows[i], bit(i)});
    }

    vector<bool> used(m, false);
    vector<unsigned int> piv_of_col(ncols, -1);
    for (unsigned int c = 0; c < ncols; ++c) {
      optional<unsigned int> piv = nullopt;
      for (unsigned int i = 0; i < m; ++i) {
        if (!used[i] && work[i].first.test(c)) {
          piv = i;
          break;
        }
      }
      check(piv.has_value(), "evaluation map is not full rank");
      used[*piv]     = true;
      piv_of_col[c]  = *piv;
      const auto& pr = work[*piv].first;
      const auto& pt = work[*piv].second;
      for (unsigned int i = 0; i < m; ++i) {
        if (i != piv && work[i].first.test(c)) {
          work[i].first ^= pr;
          work[i].second ^= pt;
        }
      }
    }

    vector<big_int> out;
    out.reserve(ncols);
    for (unsigned int c = 0; c < ncols; ++c) {
      out.push_back(work[piv_of_col[c]].second);
    }
    return out;
  }

  constexpr unsigned int F4M[4][4] = {
      {0, 0, 0, 0},
      {0, 1, 2, 3},
      {0, 2, 3, 1},
      {0, 3, 1, 2},
  };

  constexpr unsigned int F4INV[4] = {0, 1, 3, 2};

  constexpr unsigned int f4_mul(unsigned int a, unsigned int b) {
    return F4M[a][b];
  }

  constexpr unsigned int f4_pow(unsigned int t, unsigned int j) {
    unsigned int r = 1;
    for (unsigned int i = 0; i < j; ++i) {
      r = f4_mul(r, t);
    }
    return r;
  }

  vector<unsigned int> f4_pmul(const vector<unsigned int>& a, const vector<unsigned int>& b) {
    vector<unsigned int> out(a.size() + b.size() - 1, 0);
    for (size_t i = 0; i < a.size(); ++i) {
      if (!a[i]) {
        continue;
      }
      for (size_t j = 0; j < b.size(); ++j) {
        out[i + j] ^= f4_mul(a[i], b[j]);
      }
    }
    return out;
  }

  vector<unsigned int> f4_pmod(const vector<unsigned int>& a, const vector<unsigned int>& q) {
    check(q.size(), "q needs to be non-zero");
    const unsigned int dq     = q.size() - 1;
    vector<unsigned int> work = a;
    for (int i = work.size() - 1; i >= (int)dq; --i) {
      const auto c = work[i];
      if (!c) {
        continue;
      }
      for (unsigned int j = 0; j <= dq; ++j) {
        work[i - dq + j] ^= f4_mul(c, q[j]);
      }
    }
    work.resize(dq, 0);
    return work;
  }

  unsigned int f4_peval(const vector<unsigned int>& p, unsigned int t) {
    return accumulate(p.rbegin(), p.rend(), 0u,
                      [&t](unsigned int l, unsigned int r) { return f4_mul(l, t) ^ r; });
  }

  template <typename Fn>
  void f4_tuples_rec(unsigned int len, vector<unsigned int>& cur, Fn fn) {
    if (cur.size() == len) {
      fn(cur);
    } else {
      for (unsigned int v = 0; v < 4; ++v) {
        cur.push_back(v);
        f4_tuples_rec(len, cur, fn);
        cur.pop_back();
      }
    }
  }

  bool f4_poly_has_roots(const vector<unsigned int>& p) {
    for (unsigned int u = 0; u != 4; ++u) {
      if (!f4_peval(p, u)) {
        return true;
      }
    }
    return false;
  }

  vector<vector<unsigned int>> f4_irreducibles_impl(unsigned int m) {
    vector<vector<unsigned int>> quads;
    if (m == 4) {
      vector<unsigned int> cur;
      f4_tuples_rec(2, cur, [&](const auto& t) {
        vector<unsigned int> p = t;
        p.push_back(1);
        if (!f4_poly_has_roots(p)) {
          quads.push_back(p);
        }
      });
    }

    vector<vector<unsigned int>> out;
    vector<unsigned int> cur;
    f4_tuples_rec(m, cur, [&](const auto& t) {
      auto p = t;
      p.push_back(1);
      if (f4_poly_has_roots(p)) {
        return;
      }
      if (m == 4) {
        for (const auto& q2 : quads) {
          const auto rem = f4_pmod(p, q2);
          if (all_of(rem.begin(), rem.end(), logical_not{})) {
            return;
          }
        }
      }
      out.push_back(p);
    });
    return out;
  }

  cached_eval<vector<vector<unsigned int>>, unsigned int> f4_irreducibles{f4_irreducibles_impl};

  vector<vector<unsigned int>> f4_xpow_cols(const vector<unsigned int>& q, unsigned int count) {
    const unsigned int dq = q.size() - 1;
    vector<unsigned int> cur(dq, 0);
    cur[0] = 1;
    vector<vector<unsigned int>> cols;
    cols.reserve(count);
    for (unsigned int i = 0; i < count; ++i) {
      cols.push_back(cur);
      vector<unsigned int> shifted;
      shifted.reserve(cur.size() + 1);
      shifted.push_back(0);
      shifted.insert(shifted.end(), cur.begin(), cur.end());
      cur = f4_pmod(shifted, q);
    }
    return cols;
  }

  vector<vector<unsigned int>> f4_matinv(const vector<vector<unsigned int>>& m) {
    const unsigned int n = m.size();
    vector<vector<unsigned int>> a;
    a.reserve(m.size());
    for (unsigned int i = 0; i < n; ++i) {
      auto row = m[i];
      row.resize(2 * n, 0);
      row[n + i] = 1;
      a.push_back(std::move(row));
    }

    for (unsigned int c = 0; c < n; ++c) {
      optional<unsigned int> piv = nullopt;
      for (unsigned int i = c; i < n; ++i) {
        if (a[i][c]) {
          piv = i;
          break;
        }
      }
      check(piv.has_value(), "F4 matrix is singular");
      swap(a[c], a[*piv]);
      const auto inv = F4INV[a[c][c]];
      for (unsigned int j = 0; j < 2 * n; ++j) {
        a[c][j] = f4_mul(inv, a[c][j]);
      }
      for (unsigned int i = 0; i < n; ++i) {
        if (i == c || !a[i][c]) {
          continue;
        }
        const auto factor = a[i][c];
        for (unsigned int j = 0; j < 2 * n; ++j) {
          a[i][j] ^= f4_mul(factor, a[c][j]);
        }
      }
    }

    vector<vector<unsigned int>> inv;
    inv.reserve(m.size());
    for (unsigned int i = 0; i < n; ++i) {
      inv.emplace_back(a[i].begin() + n, a[i].end());
    }
    return inv;
  }

  struct Alg {
    unsigned int na = 0;
    unsigned int nb = 0;
    vector<pair<big_int, big_int>> gates;
    vector<big_int> w;
  };

  Alg precompose(const Alg& alg, const vector<big_int>& ar, const vector<big_int>& br,
                 unsigned int na, unsigned int nb) {
    Alg out;
    out.na = na;
    out.nb = nb;
    out.gates.reserve(alg.gates.size());
    for (const auto& [fa, fb] : alg.gates) {
      out.gates.push_back({combine(fa, ar), combine(fb, br)});
    }
    out.w = alg.w;
    return out;
  }

  struct Place {
    enum class Kind { Inf, Lin, Irr };

    Kind kind      = Kind::Inf;
    big_int poly   = 0;
    unsigned int e = 1;
  };

  struct Plan {
    enum class Kind { Base, Mont5, Split, Crt, Full, ShortSplit };

    Kind kind      = Kind::Base;
    unsigned int s = 0;
    unsigned int m = 0;
    vector<Place> picks;
  };

  struct CostPlan {
    unsigned int gates = 0;
    Plan plan;
  };

  pair<Alg, vector<big_int>> place_alg(const Place& pk, unsigned int na, unsigned int nb);
  Alg full_alg_impl(unsigned int n);
  Alg short_alg_impl(unsigned int n);
  cached_eval<Alg, unsigned int> full_alg{full_alg_impl};
  cached_eval<Alg, unsigned int> short_alg{short_alg_impl};

  constexpr unsigned int TOWER_FIELD_6 = 15;
  constexpr unsigned int TOWER_FIELD_8 = 24;

  unsigned int tower_field_cost(unsigned int d) {
    if (d == 6) {
      return TOWER_FIELD_6;
    }
    if (d == 8) {
      return TOWER_FIELD_8;
    }
    return 0;
  }

  pair<vector<vector<unsigned int>>, vector<vector<unsigned int>>>
  f4_full_product_alg(unsigned int m) {
    vector<vector<unsigned int>> gates4;
    vector<vector<unsigned int>> ev;
    vector<vector<pair<unsigned int, unsigned int>>> r4;

    const unsigned int npts = m == 2 ? 2 : 4;
    for (unsigned int t = 0; t < npts; ++t) {
      r4.push_back({{gates4.size(), 1}});
      vector<unsigned int> gate;
      vector<unsigned int> row;
      for (unsigned int j = 0; j < m; ++j) {
        gate.push_back(f4_pow(t, j));
      }
      for (unsigned int j = 0; j < 2 * m - 1; ++j) {
        row.push_back(f4_pow(t, j));
      }
      gates4.push_back(std::move(gate));
      ev.push_back(std::move(row));
    }

    r4.push_back({{gates4.size(), 1}});
    vector<unsigned int> inf_gate(m, 0);
    inf_gate[m - 1] = 1;
    gates4.push_back(inf_gate);
    vector<unsigned int> inf_row(2 * m - 1, 0);
    inf_row[2 * m - 2] = 1;
    ev.push_back(inf_row);

    if (m == 4) {
      const auto q2      = f4_irreducibles(2)[0];
      const auto beta    = q2[0];
      const auto alpha   = q2[1];
      const auto cols_in = f4_xpow_cols(q2, m);
      vector<unsigned int> r0;
      vector<unsigned int> r1;
      vector<unsigned int> rx;
      for (const auto& col : cols_in) {
        r0.push_back(col[0]);
        r1.push_back(col[1]);
        rx.push_back(col[0] ^ col[1]);
      }
      const unsigned int ga = gates4.size();
      gates4.push_back(r0);
      gates4.push_back(r1);
      gates4.push_back(std::move(rx));

      const auto cols_out = f4_xpow_cols(q2, 2 * m - 1);
      vector<unsigned int> ev0;
      vector<unsigned int> ev1;
      for (const auto& col : cols_out) {
        ev0.push_back(col[0]);
        ev1.push_back(col[1]);
      }
      ev.push_back(std::move(ev0));
      ev.push_back(std::move(ev1));
      r4.push_back({{ga, 1}, {ga + 1, beta}});
      r4.push_back({{ga, 1}, {ga + 1, 1 ^ alpha}, {ga + 2, 1}});
    }

    const auto inv        = f4_matinv(ev);
    const unsigned int ng = gates4.size();
    vector<vector<unsigned int>> w4(2 * m - 1, vector<unsigned int>(ng, 0));
    for (unsigned int k = 0; k < 2 * m - 1; ++k) {
      for (unsigned int j = 0; j < r4.size(); ++j) {
        const auto c = inv[k][j];
        if (!c) {
          continue;
        }
        for (const auto& [g, coeff] : r4[j]) {
          w4[k][g] ^= f4_mul(c, coeff);
        }
      }
    }
    return {gates4, w4};
  }

  big_int t_to_bits(const vector<unsigned int>& t) {
    big_int out{0};
    for (unsigned int j = 0; j < t.size(); ++j) {
      out |= big_int(t[j]) << (2 * j);
    }
    return out;
  }

  pair<Alg, vector<big_int>> tower_place_alg(const big_int& q, unsigned int na, unsigned int nb) {
    const unsigned int d     = pdeg(q).value_or(0);
    const unsigned int m     = d / 2;
    const auto tower_modulus = f4_irreducibles(m)[0];
    const auto [gates4, w4]  = f4_full_product_alg(m);

    const auto rq          = f4_xpow_cols(tower_modulus, 2 * m - 1);
    const unsigned int ng4 = gates4.size();
    vector<vector<unsigned int>> w4m(m, vector<unsigned int>(ng4, 0));
    for (unsigned int k = 0; k < m; ++k) {
      for (unsigned int j = 0; j < 2 * m - 1; ++j) {
        const auto c = rq[j][k];
        for (unsigned int g = 0; g < ng4; ++g) {
          w4m[k][g] ^= f4_mul(c, w4[j][g]);
        }
      }
    }

    optional<vector<unsigned int>> rho;
    vector<unsigned int> cand;
    f4_tuples_rec(m, cand, [&](const auto& c) {
      if (rho.has_value()) {
        return;
      }
      vector<unsigned int> acc(m, 0);
      for (int i = d; i >= 0; --i) {
        acc = f4_pmod(f4_pmul(acc, c), tower_modulus);
        if (q.test(i)) {
          acc[0] ^= 1;
        }
      }
      if (all_of(acc.begin(), acc.end(), [](auto v) { return !v; })) {
        rho = c;
      }
    });
    check(rho.has_value(), "no root of q in F4 tower field");

    vector<big_int> cols;
    vector<unsigned int> power(1, 1);
    for (unsigned int i = 0; i < d; ++i) {
      cols.push_back(t_to_bits(power));
      power = f4_pmod(f4_pmul(power, *rho), tower_modulus);
    }

    vector<big_int> mrows(d, 0);
    for (unsigned int r = 0; r < d; ++r) {
      big_int row = 0;
      for (unsigned int i = 0; i < d; ++i) {
        if (cols[i].test(r)) {
          row.set(i);
        }
      }
      mrows[r] = row;
    }
    const auto minv = left_inverse(mrows, d);

    const auto red_a = reduction_rows(q, na);
    const auto red_b = reduction_rows(q, nb);
    vector<big_int> trow_a(d);
    vector<big_int> trow_b(d);
    for (unsigned int r = 0; r < d; ++r) {
      trow_a[r] = combine(mrows[r], red_a);
      trow_b[r] = combine(mrows[r], red_b);
    }

    auto form_bits = [&](const auto& l, const auto& trow) {
      big_int p = 0;
      big_int r = 0;
      for (unsigned int j = 0; j < l.size(); ++j) {
        const auto c = l[j];
        if (!c) {
          continue;
        }

        const auto& pj = trow[2 * j];
        const auto& rj = trow[2 * j + 1];
        if (c & 1) {
          p ^= pj;
          r ^= rj;
        }
        if (c & 2) {
          p ^= rj;
          r ^= pj ^ rj;
        }
      }
      return make_pair(p, r);
    };

    Alg alg;
    alg.na = na;
    alg.nb = nb;
    vector<pair<big_int, big_int>> comp;
    for (const auto& l : gates4) {
      const auto [pa, ra]  = form_bits(l, trow_a);
      const auto [pb, rb]  = form_bits(l, trow_b);
      const unsigned int b = alg.gates.size();
      alg.gates.push_back({pa, pb});
      alg.gates.push_back({ra, rb});
      alg.gates.push_back({pa ^ ra, pb ^ rb});
      comp.push_back({bit(b) | bit(b + 1), bit(b) | bit(b + 2)});
    }

    vector<big_int> towout(d, 0);
    for (unsigned int k = 0; k < m; ++k) {
      big_int o0 = 0;
      big_int o1 = 0;
      for (unsigned int g = 0; g < ng4; ++g) {
        const auto c = w4m[k][g];
        if (!c) {
          continue;
        }

        const auto [g0, g1] = comp[g];
        if ((c & 1)) {
          o0 ^= g0;
          o1 ^= g1;
        }
        if ((c & 2)) {
          o0 ^= g1;
          o1 ^= g0 ^ g1;
        }
      }
      towout[2 * k]     = o0;
      towout[2 * k + 1] = o1;
    }

    alg.w.reserve(d);
    for (unsigned int r = 0; r < d; ++r) {
      alg.w.push_back(combine(minv[r], towout));
    }
    return {alg, reduction_rows(q, na + nb - 1)};
  }

  const array<big_int, 13> MONT5_GATES = {
      big_int(0b11111), big_int(0b11101), big_int(0b10111), big_int(0b11011), big_int(0b01101),
      big_int(0b10110), big_int(0b11000), big_int(0b00011), big_int(0b10001), big_int(0b10000),
      big_int(0b01000), big_int(0b00010), big_int(0b00001),
  };

  const array<vector<unsigned int>, 9> MONT5_W = {
      vector<unsigned int>{12},
      vector<unsigned int>{7, 11, 12},
      vector<unsigned int>{2, 5, 7, 8, 9, 12},
      vector<unsigned int>{0, 1, 3, 6, 8, 9},
      vector<unsigned int>{0, 4, 5, 6, 7, 9, 10, 11, 12},
      vector<unsigned int>{0, 2, 3, 7, 8, 12},
      vector<unsigned int>{1, 4, 6, 8, 9, 12},
      vector<unsigned int>{6, 9, 10},
      vector<unsigned int>{9},
  };

  Alg mont5_alg() {
    Alg alg;
    alg.na = 5;
    alg.nb = 5;
    for (const big_int& mask : MONT5_GATES) {
      alg.gates.push_back({mask, mask});
    }
    for (const auto& row : MONT5_W) {
      big_int w = 0;
      for (auto g : row) {
        w.set(g);
      }
      alg.w.push_back(w);
    }
    return alg;
  }

  CostPlan short_cost_impl(unsigned int n);
  cached_eval<CostPlan, unsigned int> short_cost{short_cost_impl};

  unsigned int place_gate_cost(const Place& pk);

  struct CtrSearchOption {
    unsigned int degree = 0;
    unsigned int cost   = 0;
    vector<Place> picks;
  };

  optional<pair<unsigned int, vector<Place>>>
  crt_search(unsigned int target, unsigned int maxe_lin, unsigned int cap_irr,
             const map<unsigned int, vector<big_int>>& irr_pools) {
    vector<vector<CtrSearchOption>> groups;

    for (const auto& spec : {pair<Place::Kind, big_int>{Place::Kind::Inf, 0},
                             pair<Place::Kind, big_int>{Place::Kind::Lin, 0b10},
                             pair<Place::Kind, big_int>{Place::Kind::Lin, 0b11}}) {
      vector<CtrSearchOption> opts;
      for (unsigned int e = 1; e <= maxe_lin; ++e) {
        Place pk{spec.first, spec.second, e};
        opts.push_back({e, short_cost(e).gates, {pk}});
      }
      if (!opts.empty()) {
        groups.push_back(std::move(opts));
      }
    }

    for (const auto& [d, pool] : irr_pools) {
      if (pool.empty()) {
        continue;
      }
      if (d == 2) {
        const big_int& p = pool[0];
        vector<CtrSearchOption> opts;
        for (unsigned int e = 1; d * e <= cap_irr; ++e) {
          Place pk{Place::Kind::Irr, p, e};
          opts.push_back({d * e, place_gate_cost(pk), {pk}});
        }
        if (!opts.empty()) {
          groups.push_back(std::move(opts));
        }
      } else {
        Place first{Place::Kind::Irr, pool[0], 1};
        const auto per  = place_gate_cost(first);
        const auto maxt = min<size_t>(pool.size(), target / d + 1);
        vector<CtrSearchOption> opts;
        for (unsigned int t = 1; t <= maxt; ++t) {
          vector<Place> picks;
          for (unsigned int i = 0; i < t; ++i) {
            picks.push_back({Place::Kind::Irr, pool[i], 1});
          }
          opts.push_back({t * d, t * per, std::move(picks)});
        }
        if (!opts.empty()) {
          groups.push_back(std::move(opts));
        }
      }
    }

    map<unsigned int, pair<unsigned int, vector<Place>>> dp;
    dp[0] = {0, {}};
    for (const auto& opts : groups) {
      auto ndp = dp;
      for (const auto& [deg0, state] : dp) {
        const auto& [cost0, picks0] = state;
        for (const auto& opt : opts) {
          const auto nd = min(deg0 + opt.degree, target);
          const auto nc = cost0 + opt.cost;
          auto it       = ndp.find(nd);
          if (it == ndp.end() || nc < it->second.first) {
            vector<Place> picks = picks0;
            picks.insert(picks.end(), opt.picks.begin(), opt.picks.end());
            ndp[nd] = {nc, std::move(picks)};
          }
        }
      }
      dp = std::move(ndp);
    }

    auto it = dp.find(target);
    if (it == dp.end()) {
      return nullopt;
    }
    return it->second;
  }

  CostPlan full_cost_impl(unsigned int n);
  cached_eval<CostPlan, unsigned int> full_cost{full_cost_impl};

  CostPlan full_cost_impl(unsigned int n) {
    optional<CostPlan> best;
    if (n == 1) {
      best = CostPlan{1, Plan{Plan::Kind::Base, 0, 0, {}}};
    } else {
      if (n == 5) {
        best = CostPlan{13, Plan{Plan::Kind::Mont5, 0, 0, {}}};
      }
      for (unsigned int s = 2; s < n; ++s) {
        const unsigned int m = (n + s - 1) / s;
        const unsigned int c = full_cost(s).gates * full_cost(m).gates;
        if (!best.has_value() || c < best->gates) {
          best = CostPlan{c, Plan{Plan::Kind::Split, s, m, {}}};
        }
      }
      map<unsigned int, vector<big_int>> pools;
      for (unsigned int d = 2; d < n; ++d) {
        pools[d] = irreducibles_of_degree(d);
      }
      auto res = crt_search(2 * n - 1, n - 1, n - 1, pools);
      if (res.has_value()) {
        const auto& [cost, picks] = *res;
        if (!best.has_value() || cost < best->gates) {
          best = CostPlan{cost, Plan{Plan::Kind::Crt, 0, 0, picks}};
        }
      }
    }
    check(best.has_value(), "full_cost failed");
    return *best;
  }

  CostPlan short_cost_impl(unsigned int n) {
    CostPlan best;
    if (n == 1) {
      best = CostPlan{1, Plan{Plan::Kind::Base, 0, 0, {}}};
    } else {
      best = CostPlan{full_cost(n).gates, Plan{Plan::Kind::Full, 0, 0, {}}};
      for (unsigned int m = (n + 1) / 2; m != n; ++m) {
        const unsigned int c = full_cost(m).gates + 2 * short_cost(n - m).gates;
        if (c < best.gates) {
          best = CostPlan{c, Plan{Plan::Kind::ShortSplit, 0, m, {}}};
        }
      }
    }
    return best;
  }

  unsigned int place_gate_cost(const Place& pk) {
    if (pk.kind == Place::Kind::Inf || pk.kind == Place::Kind::Lin) {
      return short_cost(pk.e).gates;
    }
    const unsigned int d          = pdeg(pk.poly).value_or(0);
    const unsigned int tower_cost = tower_field_cost(d);
    if (pk.e == 1 && tower_cost) {
      return min(tower_cost, full_cost(d).gates);
    }
    return full_cost(d * pk.e).gates;
  }

  Alg full_alg_impl(unsigned int n) {
    const Plan plan = full_cost(n).plan;
    switch (plan.kind) {
    case Plan::Kind::Base: {
      return Alg{1, 1, {{1, 1}}, {1}};
    }
    case Plan::Kind::Mont5: {
      return mont5_alg();
    }
    case Plan::Kind::Split: {
      const auto s           = plan.s;
      const auto m           = plan.m;
      const auto oa          = full_alg(s);
      const auto ia          = full_alg(m);
      const unsigned int ngi = ia.gates.size();

      vector<pair<big_int, big_int>> gates;
      for (const auto& [fo, go] : oa.gates) {
        vector<big_int> ar(m, 0);
        vector<big_int> br(m, 0);
        for (unsigned int i = 0; i < m; ++i) {
          big_int& fa = ar[i];
          big_int& fb = br[i];
          for_each_bit(fo, [&](unsigned long t) {
            if (t * m + i < n) {
              fa.set(t * m + i);
            }
          });
          for_each_bit(go, [&](unsigned long t) {
            if (t * m + i < n) {
              fb.set(t * m + i);
            }
          });
        }
        for (const auto& [fi, gi] : ia.gates) {
          gates.push_back({combine(fi, ar), combine(gi, br)});
        }
      }

      vector<big_int> w;
      for (unsigned int q = 0; q < 2 * n - 1; ++q) {
        big_int row = 0;
        for (unsigned int r = 0; r < 2 * s - 1; ++r) {
          if (r * m > q) {
            continue;
          }
          const unsigned int i = q - r * m;
          if (i <= 2 * m - 2) {
            for_each_bit(oa.w[r], [&](unsigned long o) { row ^= ia.w[i] << (o * ngi); });
          }
        }
        w.push_back(row);
      }
      return Alg{n, n, std::move(gates), std::move(w)};
    }
    case Plan::Kind::Crt: {
      vector<pair<big_int, big_int>> gates;
      vector<big_int> resmasks;
      vector<big_int> evrows;
      for (const Place& pk : plan.picks) {
        const auto [palg, ev]  = place_alg(pk, n, n);
        const unsigned int off = gates.size();
        gates.insert(gates.end(), palg.gates.begin(), palg.gates.end());
        for (const big_int& row : palg.w) {
          resmasks.push_back(row << off);
        }
        evrows.insert(evrows.end(), ev.begin(), ev.end());
      }
      const auto x = left_inverse(evrows, 2 * n - 1);
      vector<big_int> w;
      for (unsigned int q = 0; q < 2 * n - 1; ++q) {
        w.push_back(combine(x[q], resmasks));
      }
      return Alg{n, n, std::move(gates), std::move(w)};
    }
    default:
      fail("invalid plan");
    }
  }

  Alg short_alg_impl(unsigned int n) {
    const Plan plan = short_cost(n).plan;
    switch (plan.kind) {
    case Plan::Kind::Base: {
      return Alg{1, 1, {{1, 1}}, {1}};
    }
    case Plan::Kind::Full: {
      return full_alg(n);
    }
    case Plan::Kind::ShortSplit: {
      const auto m = plan.m;
      vector<big_int> lo;
      for (unsigned int i = 0; i < m; ++i) {
        lo.push_back(bit(i));
      }
      const Alg g1 = precompose(full_alg(m), lo, lo, n, n);
      const Alg s  = short_alg(n - m);
      vector<big_int> lo2;
      vector<big_int> hi2;
      for (unsigned int i = 0; i < n - m; ++i) {
        lo2.push_back(bit(i));
        hi2.push_back(bit(m + i));
      }
      const Alg g2          = precompose(s, lo2, hi2, n, n);
      const Alg g3          = precompose(s, hi2, lo2, n, n);
      const unsigned int o2 = g1.gates.size();
      const unsigned int o3 = o2 + g2.gates.size();

      vector<pair<big_int, big_int>> gates = g1.gates;
      gates.reserve(gates.size() + g2.gates.size() + g3.gates.size());
      gates.insert(gates.end(), g2.gates.begin(), g2.gates.end());
      gates.insert(gates.end(), g3.gates.begin(), g3.gates.end());

      vector<big_int> w;
      for (unsigned int j = 0; j < n; ++j) {
        big_int row{0};
        if (j <= 2 * m - 2) {
          row ^= g1.w[j];
        }
        if (j >= m) {
          const auto t = j - m;
          row ^= (g2.w[t] << o2) ^ (g3.w[t] << o3);
        }
        w.push_back(row);
      }

      return Alg{n, n, std::move(gates), std::move(w)};
    }

    default:
      fail("invalid plan");
    }
  }

  pair<Alg, vector<big_int>> place_alg(const Place& pk, unsigned int na, unsigned int nb) {
    const auto n2 = na + nb - 1;
    if (pk.kind == Place::Kind::Inf) {
      check(pk.e <= min(na, nb), "infinity multiplicity exceeds input width");
      vector<big_int> rows_a;
      vector<big_int> rows_b;
      vector<big_int> ev;
      for (unsigned int i = 0; i < pk.e; ++i) {
        rows_a.push_back(bit(na - 1 - i));
        rows_b.push_back(bit(nb - 1 - i));
        ev.push_back(bit(n2 - 1 - i));
      }
      return {precompose(short_alg(pk.e), rows_a, rows_b, na, nb), ev};
    }

    if (pk.kind == Place::Kind::Lin) {
      const auto q   = poly_pow(pk.poly, pk.e);
      const auto ra  = reduction_rows(q, na);
      const auto rb  = reduction_rows(q, nb);
      const auto rev = reduction_rows(q, n2);
      if (pk.poly == 0b10) {
        return {precompose(short_alg(pk.e), ra, rb, na, nb), rev};
      }

      const auto p = pascal_rows(pk.e);
      vector<big_int> ina;
      vector<big_int> inb;
      vector<big_int> ev;
      for (const auto& mask : p) {
        ina.push_back(combine(mask, ra));
        inb.push_back(combine(mask, rb));
        ev.push_back(combine(mask, rev));
      }
      return {precompose(short_alg(pk.e), ina, inb, na, nb), ev};
    }

    const auto tower_cost = tower_field_cost(pdeg(pk.poly).value_or(0));
    if (pk.e == 1 && tower_cost && tower_cost < full_cost(pdeg(pk.poly).value_or(0)).gates) {
      return tower_place_alg(pk.poly, na, nb);
    }

    const big_int q = poly_pow(pk.poly, pk.e);
    const auto d    = pdeg(q);
    const Alg full  = full_alg(*d);
    Alg alg         = precompose(full, reduction_rows(q, na), reduction_rows(q, nb), na, nb);
    const auto rout = reduction_rows(q, 2 * *d - 1);
    vector<big_int> w;
    for (unsigned int i = 0; (int)i < d; ++i) {
      w.push_back(combine(rout[i], alg.w));
    }
    alg.w = std::move(w);
    return {alg, reduction_rows(q, n2)};
  }

  struct Tables {
    vector<big_int> f;
    vector<big_int> g;
    vector<big_int> w_tree;
    vector<big_int> w_gate;
    unsigned int n_tree = 0;
    vector<pair<Place, unsigned int>> report;
  };

  Tables build_tables(unsigned int lambda, unsigned int wgrind, const big_int& p,
                      const vector<big_int>& tree_moduli, const vector<Place>& picks) {
    const unsigned int na = lambda;
    const unsigned int nb = lambda - wgrind;
    unsigned int ntree    = 0;
    for (const big_int& m : tree_moduli) {
      ntree += pdeg(m).value_or(0);
    }
    check(ntree == nb, "tree degrees must sum to lambda - wgrind");

    const auto n2 = na + nb - 1;
    vector<big_int> evrows;
    for (const big_int& m : tree_moduli) {
      const auto rows = reduction_rows(m, n2);
      evrows.insert(evrows.end(), rows.begin(), rows.end());
    }

    vector<pair<big_int, big_int>> gates;
    vector<big_int> resmasks;
    vector<pair<Place, unsigned int>> report;
    for (const Place& pk : picks) {
      const auto [palg, ev]  = place_alg(pk, na, nb);
      const unsigned int off = gates.size();
      gates.insert(gates.end(), palg.gates.begin(), palg.gates.end());
      for (const big_int& row : palg.w) {
        resmasks.emplace_back(row << off);
      }
      evrows.insert(evrows.end(), ev.begin(), ev.end());
      report.push_back({pk, palg.gates.size()});
    }

    const auto x   = left_inverse(evrows, n2);
    const auto red = reduction_rows(p, n2);
    vector<big_int> wt(lambda, 0);
    vector<big_int> wg(lambda, 0);
    for (unsigned int r = 0; r < lambda; ++r) {
      const big_int sel = combine(red[r], x);
      for_each_bit(sel, [&](unsigned long idx) {
        if (idx < ntree) {
          wt[r].set(idx);
        } else {
          wg[r] ^= resmasks[idx - ntree];
        }
      });
    }

    vector<big_int> f;
    vector<big_int> g;
    f.reserve(gates.size());
    g.reserve(gates.size());
    const auto cols = crt_lift_cols(tree_moduli);
    for (const auto& [fa, fb] : gates) {
      f.push_back(fa);
      big_int row;
      for (unsigned int j = 0; j < nb; ++j) {
        if (parity(fb & cols[j])) {
          row.set(j);
        }
      }
      g.emplace_back(row);
    }
    return {std::move(f), std::move(g), std::move(wt), std::move(wg), ntree, std::move(report)};
  }

  vector<big_int> wcrt_rows(const vector<big_int>& tree_moduli, unsigned int lambda) {
    const auto cols = crt_lift_cols(tree_moduli);
    vector<big_int> rows;
    rows.reserve(lambda);
    for (unsigned int r = 0; r < lambda; ++r) {
      big_int row;
      for (unsigned int j = 0; j < cols.size(); ++j) {
        if (cols[j].test(r)) {
          row.set(j);
        }
      }
      rows.emplace_back(row);
    }
    return rows;
  }

  void prune(vector<big_int>& f, vector<big_int>& g, vector<big_int>& wg) {
    map<pair<big_int, big_int>, int> canon;
    for (unsigned int i = 0; i < f.size(); ++i) {
      if (!f[i].is_zero() && !g[i].is_zero()) {
        const auto key = make_pair(f[i], g[i]);
        if (canon.find(key) == canon.end()) {
          canon[key] = i;
        }
      }
    }

    vector<big_int> wg1;
    wg1.reserve(wg.size());
    for (const big_int& row : wg) {
      big_int nr;
      for_each_bit(row, [&](unsigned long i) {
        if (!f[i].is_zero() && !g[i].is_zero()) {
          nr.flip(canon[make_pair(f[i], g[i])]);
        }
      });
      wg1.emplace_back(nr);
    }

    big_int used = 0;
    for (const big_int& row : wg1) {
      used |= row;
    }

    map<unsigned int, unsigned int> newidx;
    vector<big_int> f2;
    vector<big_int> g2;
    for (unsigned int i = 0; i < f.size(); ++i) {
      if (used.test(i)) {
        newidx[i] = f2.size();
        f2.push_back(f[i]);
        g2.push_back(g[i]);
      }
    }

    vector<big_int> wg2;
    wg2.reserve(wg1.size());
    for (const big_int& row : wg1) {
      big_int nr = 0;
      for_each_bit(row, [&](unsigned long i) { nr.set(newidx[i]); });
      wg2.emplace_back(nr);
    }

    f  = std::move(f2);
    g  = std::move(g2);
    wg = std::move(wg2);
  }

  uint64_t word_at(const big_int& v, unsigned int word) {
    return v.word_at(word);
  }

  constexpr unsigned int words_of(unsigned int width) {
    return (width + 63) / 64;
  }

  string upper_identifier(string s) {
    for (char& c : s) {
      if (c == '-') {
        c = '_';
      } else {
        c = static_cast<char>(toupper(static_cast<unsigned char>(c)));
      }
    }
    return s;
  }

  string prefix_for_name(const string& name) {
    return "FAEST_" + upper_identifier(name);
  }

  struct uint64_printer {
    uint64_t value;
  };

  ostream& operator<<(ostream& ofs, const uint64_printer& p) {
    const auto flags = ofs.flags();
    const auto fill  = ofs.fill();
    ofs << "UINT64_C(0x" << hex << nouppercase << setfill('0') << setw(16) << p.value << ")";
    ofs.flags(flags);
    ofs.fill(fill);
    return ofs;
  }

  struct uint16_printer {
    uint16_t value;
  };

  ostream& operator<<(ostream& ofs, const uint16_printer& p) {
    const auto flags = ofs.flags();
    const auto fill  = ofs.fill();
    ofs << "UINT16_C(0x" << hex << nouppercase << setfill('0') << setw(4) << p.value << ")";
    ofs.flags(flags);
    ofs.fill(fill);
    return ofs;
  }

  string row_literal(const big_int& r, unsigned int words) {
    ostringstream ss;
    ss << "{ ";
    for (unsigned int w = 0; w < words; ++w) {
      ss << uint64_printer{word_at(r, w)} << ", ";
    }
    ss << " }";
    return ss.str();
  }

  struct EmitTable {
    string suffix;
    const vector<big_int>& rows;
    unsigned int words;
    string rows_macro;
    string words_macro;
  };

  string emit_c(const string& path, const string& name, const vector<big_int>& f,
                const vector<big_int>& g, const vector<big_int>& wt, const vector<big_int>& wg,
                const vector<big_int>& tree_moduli, unsigned int lambda, unsigned int wgrind,
                unsigned int ntree) {
    const vector<big_int> wcrt = wcrt_rows(tree_moduli, lambda);
    const big_int m_tree       = mtree_poly(tree_moduli);
    const unsigned int ng      = f.size();
    const unsigned int tau     = tree_moduli.size();
    for (const big_int& m : tree_moduli) {
      check(pdeg(m).has_value() && pdeg(m).value() < 16,
            "tree modulus degree >= 16 does not fit one uint16_t");
    }

    const auto w_f  = words_of(lambda);
    const auto w_g  = words_of(lambda - wgrind);
    const auto w_wt = words_of(ntree);
    const auto w_wg = words_of(ng);
    const auto w_wc = w_g;
    const auto w_mt = words_of(*pdeg(m_tree) + 1);

    const string pfx = prefix_for_name(name);
    string src_path;
    if (path.size() >= 2 && path.substr(path.size() - 2) == ".h") {
      src_path = path.substr(0, path.size() - 2) + ".c";
    } else {
      src_path = path + ".c";
    }
    const string header_name = filesystem::path(path).filename().string();

    const array<EmitTable, 5> tables{{
        {"F", f, w_f, pfx + "_NGATES", pfx + "_F_WORDS"},
        {"G", g, w_g, pfx + "_NGATES", pfx + "_G_WORDS"},
        {"W_TREE", wt, w_wt, pfx + "_LAMBDA", pfx + "_W_TREE_WORDS"},
        {"W_GATE", wg, w_wg, pfx + "_LAMBDA", pfx + "_W_GATE_WORDS"},
        {"W_CRT", wcrt, w_wc, pfx + "_LAMBDA", pfx + "_W_CRT_WORDS"},
    }};

    ofstream out(path);
    check(out.good(), "failed to open header for writing: " + path);
    out << "/* generated by vole_mult_tables.cpp: " << name << "\n"
        << " *\n"
        << " * Bilinear VOLE-multiplication tables.\n"
        << " *\n"
        << " * Encoding: each table row is an F_2 bit-vector packed\n"
        << " * little-endian into uint64 words -- bit i sits in word i/64\n"
        << " * at bit position i%64, and word 0 holds the low bits.\n"
        << " */\n";
    out << "#ifndef " << pfx << "_TABLES_H\n#define " << pfx << "_TABLES_H\n\n";
    out << "#include <stdint.h>\n\n";
    out << "#define " << pfx << "_NGATES " << ng << "\n";
    out << "#define " << pfx << "_WGRIND " << wgrind << "\n";
    out << "#define " << pfx << "_F_WORDS " << w_f << "\n";
    out << "#define " << pfx << "_G_WORDS " << w_g << "\n";
    out << "#define " << pfx << "_W_TREE_WORDS " << w_wt << "\n";
    out << "#define " << pfx << "_W_GATE_WORDS " << w_wg << "\n";
    out << "#define " << pfx << "_W_CRT_WORDS " << w_wc << "\n";
    out << "#define " << pfx << "_M_TREE_WORDS " << w_mt << "\n\n";

    for (const auto& tbl : tables) {
      out << "static const uint64_t " << pfx << "_" << tbl.suffix << "[" << tbl.rows_macro << "]["
          << tbl.words_macro << "] = {\n";
      for (const big_int& row : tbl.rows) {
        out << "  " << row_literal(row, tbl.words) << ",\n";
      }
      out << "};\n\n";
    }

    out << "static const uint16_t " << pfx << "_TREE_MODULI[" << pfx << "_TAU] = { ";
    for (unsigned int i = 0; i < tau; ++i) {
      out << uint16_printer{static_cast<uint16_t>(word_at(tree_moduli[i], 0))} << ", ";
    }
    out << " };\n\n";

    out << "static const uint64_t " << pfx << "_M_TREE[" << pfx << "_M_TREE_WORDS] = { ";
    for (unsigned int w = 0; w < w_mt; ++w) {
      out << uint64_printer{word_at(m_tree, w)} << ", ";
    }
    out << " };\n\n";
    out << "#endif\n";

    return src_path;
  }

  struct Preset {
    unsigned int lambda;
    unsigned int tau;
    unsigned int wgrind;
  };

  const map<string, Preset> PRESETS = {
      {"128s", {FAEST_128_LAMBDA, FAEST_128S_TAU, FAEST_128S_W_GRIND}},
      {"128f", {FAEST_128_LAMBDA, FAEST_128F_TAU, FAEST_128F_W_GRIND}},
      {"192s", {FAEST_192_LAMBDA, FAEST_192S_TAU, FAEST_192S_W_GRIND}},
      {"192f", {FAEST_192_LAMBDA, FAEST_192F_TAU, FAEST_192F_W_GRIND}},
      {"256s", {FAEST_256_LAMBDA, FAEST_256S_TAU, FAEST_256S_W_GRIND}},
      {"256f", {FAEST_256_LAMBDA, FAEST_256F_TAU, FAEST_256F_W_GRIND}},
      {"em_192s", {FAEST_192_LAMBDA, FAEST_EM_192S_TAU, FAEST_EM_192F_W_GRIND}},
      {"em_192f", {FAEST_192_LAMBDA, FAEST_EM_192F_TAU, FAEST_EM_192S_W_GRIND}},
  };

  vector<pair<unsigned int, unsigned int>> faest_tree_spec(unsigned int lambda, unsigned int tau,
                                                           unsigned int wgrind) {
    const unsigned int n    = lambda - wgrind;
    const unsigned int d1   = n / tau + 1;
    const unsigned int tau1 = n % tau;
    vector<pair<unsigned int, unsigned int>> spec;
    if (tau1) {
      spec.push_back({d1, tau1});
    }
    spec.push_back({d1 - 1, tau - tau1});
    return spec;
  }

  struct Args {
    string preset;
    string emit_c_path;
  };

  void print_help(const char* prog) {
    cout << "Usage:\n"
         << "  " << prog << " --preset NAME --emit-c FILE.h\n\n"
         << "Search + frozen C table generation for CRT-based F_{2^lambda} VOLE multiplication.\n\n"
         << "Options:\n"
         << "  --preset NAME         one of 128s, 128f, 192s, 192f, 256s, 256f, em_192s, em_192f\n"
         << "  --emit-c FILE.h       emit C header and companion source\n";
  }

  Args parse_args(int argc, char** argv) {
    Args args;
    for (int i = 1; i < argc; ++i) {
      const string a  = argv[i];
      auto need_value = [&](const string& opt) {
        if (i + 1 >= argc) {
          fail(opt + " requires a value");
        }
        return argv[++i];
      };

      if (a == "--help" || a == "-h") {
        print_help(argv[0]);
        exit(0);
      } else if (a == "--preset") {
        args.preset = need_value(a);
      } else if (a == "--emit-c") {
        args.emit_c_path = need_value(a);
      } else {
        fail("unknown argument: " + a);
      }
    }
    return args;
  }

  constexpr unsigned int MAX_E = 10;
  constexpr unsigned int MAX_D = 10;

  void run_set(const string& name, unsigned int lambda, unsigned int wgrind,
               const vector<pair<unsigned int, unsigned int>>& tree_spec, const string& c_path) {
    unsigned int tree_sum = 0;
    for (const auto& [d, c] : tree_spec) {
      tree_sum += d * c;
    }
    check(tree_sum == lambda - wgrind,
          "tree degrees sum to " + to_string(tree_sum) +
              ", expected lambda - wgrind = " + to_string(lambda - wgrind));

    const unsigned int n2      = lambda + (lambda - wgrind) - 1;
    const unsigned int deficit = n2 - tree_sum;

    vector<big_int> tree_moduli;
    for (const auto& [d, c] : tree_spec) {
      const auto& pool = irreducibles_of_degree(d);
      if (c > pool.size()) {
        fail("only " + to_string(pool.size()) + " irreducibles of degree " + to_string(d) +
             " over F_2; cannot pick " + to_string(c));
      }
      tree_moduli.insert(tree_moduli.end(), pool.begin(), pool.begin() + c);
    }
    set<big_int> treeset(tree_moduli.begin(), tree_moduli.end());
    map<unsigned int, vector<big_int>> pools;
    for (unsigned int d = 2; d <= MAX_D; ++d) {
      for (const big_int& q : irreducibles_of_degree(d)) {
        if (treeset.find(q) == treeset.end()) {
          pools[d].push_back(q);
        }
      }
    }

    const big_int p = faest_modulus(lambda);
    auto res        = crt_search(deficit, MAX_E, MAX_D, pools);
    check(res.has_value(), "portfolio search failed; raise MAX_E, MAX_D");
    const vector<Place>& picks = res->second;

    Tables tables = build_tables(lambda, wgrind, p, tree_moduli, picks);
    prune(tables.f, tables.g, tables.w_gate);

    map<pair<string, unsigned int>, array<unsigned int, 3>> agg;
    for (const auto& [pk, gates] : tables.report) {
      const auto d =
          pk.kind == Place::Kind::Inf || pk.kind == Place::Kind::Lin ? pk.e : *pdeg(pk.poly) * pk.e;
      string type;
      if (pk.kind == Place::Kind::Irr) {
        type = "irr deg " + to_string(*pdeg(pk.poly));
      } else if (pk.kind == Place::Kind::Inf) {
        type = "inf";
      } else {
        type = "lin";
      }
      auto& row = agg[{type, pk.e}];
      row[0] += 1;
      row[1] += d;
      row[2] += gates;
    }

    const auto wcrt = wcrt_rows(tree_moduli, lambda);
    emit_c(c_path, name, tables.f, tables.g, tables.w_tree, tables.w_gate, tree_moduli, lambda,
           wgrind, tables.n_tree);
  }
} // namespace

int main(int argc, char** argv) {
  try {
    Args args = parse_args(argc, argv);
    if (args.preset.empty() || args.emit_c_path.empty()) {
      print_help(argv[0]);
      return 0;
    }

    auto it = PRESETS.find(args.preset);
    if (it == PRESETS.end()) {
      fail("unknown preset: " + args.preset);
    }
    const auto tree_spec = faest_tree_spec(it->second.lambda, it->second.tau, it->second.wgrind);
    run_set(args.preset, it->second.lambda, it->second.wgrind, tree_spec, args.emit_c_path);
    return 0;
  } catch (const std::exception& e) {
    std::cerr << "vole_mult_tables: " << e.what() << "\n";
    return 1;
  }
}
