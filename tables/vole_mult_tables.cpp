/*
 *  SPDX-License-Identifier: MIT
 */

// TODO: This is a pretty direct translation from the Python script and should be cleaned up at some
// point.

#include <boost/multiprecision/cpp_int.hpp>
#include <boost/multiprecision/integer.hpp>

#include <algorithm>
#include <array>
#include <cstdint>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <optional>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

using Int = boost::multiprecision::cpp_int;

namespace {
  [[noreturn]] void fail(const std::string& msg) {
    throw std::runtime_error(msg);
  }

  void check(bool cond, const std::string& msg) {
    if (!cond) {
      fail(msg);
    }
  }

  Int bit(unsigned int i) {
    return Int(1) << i;
  }

  bool test_bit(const Int& x, unsigned int i) {
    return ((x >> i) & 1) != 0;
  }

  template <typename Fn>
  void for_each_bit(Int x, Fn fn) {
    while (x != 0) {
      const auto b = boost::multiprecision::lsb(x);
      fn(b);
      x ^= bit(b);
    }
  }

  std::optional<unsigned int> pdeg(const Int& a) {
    if (a == 0) {
      return std::nullopt;
    }
    return boost::multiprecision::msb(a);
  }

  Int pmul(Int a, Int b) {
    Int r = 0;
    while (b != 0) {
      if ((b & 1) != 0) {
        r ^= a;
      }
      b >>= 1;
      a <<= 1;
    }
    return r;
  }

  Int pmod(Int a, const Int& m) {
    const auto d = pdeg(m);
    check(d.has_value(), "polynomial modulus must be nonzero");
    while (a != 0 && *pdeg(a) >= *d) {
      a ^= m << (*pdeg(a) - *d);
    }
    return a;
  }

  Int pmulmod(const Int& a, const Int& b, const Int& m) {
    return pmod(pmul(a, b), m);
  }

  Int pgcd(Int a, Int b) {
    while (b != 0) {
      const Int r = pmod(a, b);
      a           = b;
      b           = r;
    }
    return a;
  }

  std::pair<Int, Int> pdivmod(Int a, const Int& b) {
    Int q         = 0;
    const auto db = pdeg(b);
    check(db.has_value(), "polynomial divisor must be nonzero");
    while (a != 0 && *pdeg(a) >= *db) {
      const int s = *pdeg(a) - *db;
      q |= bit(s);
      a ^= b << s;
    }
    return {q, a};
  }

  Int pinvmod(const Int& a, const Int& m) {
    Int r0 = m;
    Int r1 = pmod(a, m);
    Int s0 = 0;
    Int s1 = 1;
    while (r1 != 0) {
      const auto [q, r] = pdivmod(r0, r1);
      r0                = r1;
      r1                = r;
      const Int ns      = s0 ^ pmul(q, s1);
      s0                = s1;
      s1                = ns;
    }
    check(r0 == 1, "pinvmod: polynomials are not coprime");
    return pmod(s0, m);
  }

  Int poly_pow(const Int& p, unsigned int e) {
    Int r = 1;
    for (unsigned int i = 0; i < e; ++i) {
      r = pmul(r, p);
    }
    return r;
  }

  Int mtree_poly(const std::vector<Int>& tree_moduli) {
    Int m = 1;
    for (const Int& q : tree_moduli) {
      m = pmul(m, q);
    }
    return m;
  }

  std::vector<Int> crt_lift_cols(const std::vector<Int>& tree_moduli) {
    Int mtree = 1;
    for (const Int& m : tree_moduli) {
      mtree = pmul(mtree, m);
    }

    std::vector<Int> cols;
    for (const Int& m : tree_moduli) {
      const auto [q, rem] = pdivmod(mtree, m);
      check(rem == 0, "tree modulus does not divide M_tree");
      const Int e = pmulmod(q, pinvmod(q, m), mtree);
      for (unsigned int b = 0; b < pdeg(m).value_or(0); ++b) {
        cols.push_back(pmulmod(e, bit(b), mtree));
      }
    }
    return cols;
  }

  Int x_pow_2k(const Int& f, unsigned int k) {
    Int r = pmod(2, f);
    for (unsigned int i = 0; i < k; ++i) {
      r = pmulmod(r, r, f);
    }
    return r;
  }

  std::set<unsigned int> prime_factors(unsigned int n) {
    std::set<unsigned int> fs;
    for (unsigned int d = 2; d * d <= n; ++d) {
      while (n % d == 0) {
        fs.insert(d);
        n /= d;
      }
    }
    if (n > 1) {
      fs.insert(n);
    }
    return fs;
  }

  bool is_irreducible(const Int& f) {
    const auto n = pdeg(f);
    if (n.value_or(0) <= 0) {
      return false;
    }
    if (*n == 1) {
      return true;
    }
    if ((f & 1) == 0) {
      return false;
    }
    if (x_pow_2k(f, *n) != pmod(2, f)) {
      return false;
    }
    for (auto q : prime_factors(*n)) {
      const Int h = pmod(x_pow_2k(f, *n / q) ^ 2, f);
      if (pgcd(f, h) != 1) {
        return false;
      }
    }
    return true;
  }

  const std::vector<Int>& irreducibles_of_degree(unsigned int d) {
    static std::map<unsigned int, std::vector<Int>> cache;
    auto it = cache.find(d);
    if (it != cache.end()) {
      return it->second;
    }

    std::vector<Int> out;
    for (Int f = bit(d); f < bit(d + 1); ++f) {
      if (is_irreducible(f)) {
        out.push_back(f);
      }
    }
    auto inserted = cache.emplace(d, std::move(out));
    return inserted.first->second;
  }

  Int faest_modulus(unsigned int lam) {
    switch (lam) {
    case 128:
      return bit(128) | bit(7) | bit(2) | bit(1) | 1;
    case 192:
      return bit(192) | bit(7) | bit(2) | bit(1) | 1;
    case 256:
      return bit(256) | bit(10) | bit(5) | bit(2) | 1;
    default:
      fail("unsupported FAEST modulus");
    }
  }

  Int combine(const Int& mask, const std::vector<Int>& rows) {
    Int out = 0;
    for_each_bit(mask, [&](unsigned long i) {
      check(i < rows.size(), "row-combination mask exceeds row count");
      out ^= rows[i];
    });
    return out;
  }

  unsigned int parity(Int x) {
    bool p = false;
    for_each_bit(x, [&](unsigned long) { p = !p; });
    return p ? 1 : 0;
  }

  std::vector<Int> reduction_rows(const Int& q, unsigned int ncols) {
    const auto d = pdeg(q);
    std::vector<Int> rows(d.value_or(0), 0);
    Int cur = 1;
    for (unsigned int j = 0; j < ncols; ++j) {
      for_each_bit(cur, [&](unsigned long i) {
        check((int)i < d, "reduction row index out of range");
        rows[i] |= bit(j);
      });
      cur = pmod(cur << 1, q);
    }
    return rows;
  }

  std::vector<Int> pascal_rows(unsigned int e) {
    std::vector<Int> rows;
    rows.reserve(e);
    for (unsigned int i = 0; i < e; ++i) {
      Int row = 0;
      for (unsigned int j = 0; j < e; ++j) {
        if ((i & j) == i) {
          row |= bit(j);
        }
      }
      rows.push_back(row);
    }
    return rows;
  }

  std::vector<Int> left_inverse(const std::vector<Int>& rows, unsigned int ncols) {
    const unsigned int m = rows.size();
    std::vector<std::pair<Int, Int>> work;
    work.reserve(rows.size());
    for (unsigned int i = 0; i < m; ++i) {
      work.push_back({rows[i], bit(i)});
    }

    std::vector<bool> used(m, false);
    std::vector<unsigned int> piv_of_col(ncols, -1);
    for (unsigned int c = 0; c < ncols; ++c) {
      std::optional<unsigned int> piv = std::nullopt;
      for (unsigned int i = 0; i < m; ++i) {
        if (!used[i] && test_bit(work[i].first, c)) {
          piv = i;
          break;
        }
      }
      if (!piv.has_value()) {
        fail("evaluation map is not full rank");
      }
      used[*piv]    = true;
      piv_of_col[c] = *piv;
      const Int pr  = work[*piv].first;
      const Int pt  = work[*piv].second;
      for (unsigned int i = 0; i < m; ++i) {
        if (i != piv && test_bit(work[i].first, c)) {
          work[i].first ^= pr;
          work[i].second ^= pt;
        }
      }
    }

    std::vector<Int> out;
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

  std::vector<unsigned int> f4_pmul(const std::vector<unsigned int>& a,
                                    const std::vector<unsigned int>& b) {
    std::vector<unsigned int> out(a.size() + b.size() - 1, 0);
    for (std::size_t i = 0; i < a.size(); ++i) {
      if (a[i] == 0) {
        continue;
      }
      for (std::size_t j = 0; j < b.size(); ++j) {
        out[i + j] ^= f4_mul(a[i], b[j]);
      }
    }
    return out;
  }

  std::vector<unsigned int> f4_pmod(const std::vector<unsigned int>& a,
                                    const std::vector<unsigned int>& q) {
    const int dq                   = q.size() - 1;
    std::vector<unsigned int> work = a;
    for (int i = work.size() - 1; i >= dq; --i) {
      const auto c = work[i];
      if (c == 0) {
        continue;
      }
      for (unsigned int j = 0; (int)j <= dq; ++j) {
        work[i - dq + j] ^= f4_mul(c, q[j]);
      }
    }
    work.resize(dq, 0);
    return work;
  }

  unsigned int f4_peval(const std::vector<unsigned int>& p, unsigned int t) {
    unsigned int r = 0;
    for (auto it = p.rbegin(); it != p.rend(); ++it) {
      r = f4_mul(r, t) ^ *it;
    }
    return r;
  }

  template <typename Fn>
  void f4_tuples_rec(unsigned int len, std::vector<unsigned int>& cur, Fn fn) {
    if (cur.size() == len) {
      fn(cur);
      return;
    }
    for (unsigned int v = 0; v < 4; ++v) {
      cur.push_back(v);
      f4_tuples_rec(len, cur, fn);
      cur.pop_back();
    }
  }

  std::vector<std::vector<unsigned int>> f4_irreducibles(unsigned int m) {
    std::vector<std::vector<unsigned int>> quads;
    if (m == 4) {
      std::vector<unsigned int> cur;
      f4_tuples_rec(2, cur, [&](const auto& t) {
        std::vector<unsigned int> p = t;
        p.push_back(1);
        bool no_roots = true;
        for (unsigned int u = 0; u < 4; ++u) {
          if (f4_peval(p, u) == 0) {
            no_roots = false;
            break;
          }
        }
        if (no_roots) {
          quads.push_back(p);
        }
      });
    }

    std::vector<std::vector<unsigned int>> out;
    std::vector<unsigned int> cur;
    f4_tuples_rec(m, cur, [&](const auto& t) {
      std::vector<unsigned int> p = t;
      p.push_back(1);
      for (unsigned int u = 0; u < 4; ++u) {
        if (f4_peval(p, u) == 0) {
          return;
        }
      }
      if (m == 4) {
        for (const auto& q2 : quads) {
          const auto rem = f4_pmod(p, q2);
          if (std::all_of(rem.begin(), rem.end(), [](auto v) { return v == 0; })) {
            return;
          }
        }
      }
      out.push_back(p);
    });
    return out;
  }

  std::vector<std::vector<unsigned int>> f4_xpow_cols(const std::vector<unsigned int>& q,
                                                      unsigned int count) {
    const unsigned int dq = q.size() - 1;
    std::vector<unsigned int> cur(dq, 0);
    cur[0] = 1;
    std::vector<std::vector<unsigned int>> cols;
    cols.reserve(count);
    for (unsigned int i = 0; i < count; ++i) {
      cols.push_back(cur);
      std::vector<unsigned int> shifted;
      shifted.reserve(cur.size() + 1);
      shifted.push_back(0);
      shifted.insert(shifted.end(), cur.begin(), cur.end());
      cur = f4_pmod(shifted, q);
    }
    return cols;
  }

  std::vector<std::vector<unsigned int>>
  f4_matinv(const std::vector<std::vector<unsigned int>>& m) {
    const unsigned int n            = m.size();
    constexpr unsigned int F4INV[4] = {0, 1, 3, 2};
    std::vector<std::vector<unsigned int>> a;
    a.reserve(m.size());
    for (unsigned int i = 0; i < n; ++i) {
      auto row = m[i];
      row.resize(2 * n, 0);
      row[n + i] = 1;
      a.push_back(std::move(row));
    }

    for (unsigned int c = 0; c < n; ++c) {
      std::optional<unsigned int> piv = std::nullopt;
      for (unsigned int i = c; i < n; ++i) {
        if (a[i][c] != 0) {
          piv = i;
          break;
        }
      }
      check(piv.has_value(), "F4 matrix is singular");
      std::swap(a[c], a[*piv]);
      const auto inv = F4INV[a[c][c]];
      for (unsigned int j = 0; j < 2 * n; ++j) {
        a[c][j] = f4_mul(inv, a[c][j]);
      }
      for (unsigned int i = 0; i < n; ++i) {
        if (i == c || a[i][c] == 0) {
          continue;
        }
        const auto factor = a[i][c];
        for (unsigned int j = 0; j < 2 * n; ++j) {
          a[i][j] ^= f4_mul(factor, a[c][j]);
        }
      }
    }

    std::vector<std::vector<unsigned int>> inv;
    inv.reserve(m.size());
    for (unsigned int i = 0; i < n; ++i) {
      inv.emplace_back(a[i].begin() + n, a[i].end());
    }
    return inv;
  }

  struct Alg {
    unsigned int na = 0;
    unsigned int nb = 0;
    std::vector<std::pair<Int, Int>> gates;
    std::vector<Int> w;
  };

  Alg precompose(const Alg& alg, const std::vector<Int>& ar, const std::vector<Int>& br,
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
    Int poly       = 0;
    unsigned int e = 1;
  };

  struct Plan {
    enum class Kind { Base, Mont5, Split, Crt, Full, ShortSplit };

    Kind kind      = Kind::Base;
    unsigned int s = 0;
    unsigned int m = 0;
    std::vector<Place> picks;
  };

  struct CostPlan {
    unsigned int gates = 0;
    Plan plan;
  };

  std::pair<Alg, std::vector<Int>> place_alg(const Place& pk, unsigned int na, unsigned int nb);
  Alg full_alg(unsigned int n);
  Alg short_alg(unsigned int n);

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

  std::pair<std::vector<std::vector<unsigned int>>, std::vector<std::vector<unsigned int>>>
  f4_full_product_alg(unsigned int m) {
    std::vector<std::vector<unsigned int>> gates4;
    std::vector<std::vector<unsigned int>> ev;
    std::vector<std::vector<std::pair<unsigned int, unsigned int>>> r4;

    const unsigned int npts = m == 2 ? 2 : 4;
    for (unsigned int t = 0; t < npts; ++t) {
      r4.push_back({{gates4.size(), 1}});
      std::vector<unsigned int> gate;
      std::vector<unsigned int> row;
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
    std::vector<unsigned int> inf_gate(m, 0);
    inf_gate[m - 1] = 1;
    gates4.push_back(inf_gate);
    std::vector<unsigned int> inf_row(2 * m - 1, 0);
    inf_row[2 * m - 2] = 1;
    ev.push_back(inf_row);

    if (m == 4) {
      const auto q2      = f4_irreducibles(2)[0];
      const auto beta    = q2[0];
      const auto alpha   = q2[1];
      const auto cols_in = f4_xpow_cols(q2, m);
      std::vector<unsigned int> r0;
      std::vector<unsigned int> r1;
      for (const auto& col : cols_in) {
        r0.push_back(col[0]);
        r1.push_back(col[1]);
      }
      const unsigned int ga = gates4.size();
      gates4.push_back(r0);
      gates4.push_back(r1);
      std::vector<unsigned int> rx;
      for (unsigned int i = 0; i < m; ++i) {
        rx.push_back(r0[i] ^ r1[i]);
      }
      gates4.push_back(std::move(rx));

      const auto cols_out = f4_xpow_cols(q2, 2 * m - 1);
      std::vector<unsigned int> ev0;
      std::vector<unsigned int> ev1;
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
    std::vector<std::vector<unsigned int>> w4(2 * m - 1, std::vector<unsigned int>(ng, 0));
    for (unsigned int k = 0; k < 2 * m - 1; ++k) {
      for (unsigned int j = 0; j < r4.size(); ++j) {
        const auto c = inv[k][j];
        if (c == 0) {
          continue;
        }
        for (const auto& [g, coeff] : r4[j]) {
          w4[k][g] ^= f4_mul(c, coeff);
        }
      }
    }
    return {gates4, w4};
  }

  Int t_to_bits(const std::vector<unsigned int>& t) {
    Int out = 0;
    for (unsigned int j = 0; j < t.size(); ++j) {
      out |= Int(t[j]) << (2 * j);
    }
    return out;
  }

  std::pair<Alg, std::vector<Int>> tower_place_alg(const Int& q, unsigned int na, unsigned int nb) {
    const unsigned int d     = pdeg(q).value_or(0);
    const unsigned int m     = d / 2;
    const auto tower_modulus = f4_irreducibles(m)[0];
    const auto [gates4, w4]  = f4_full_product_alg(m);

    const auto rq          = f4_xpow_cols(tower_modulus, 2 * m - 1);
    const unsigned int ng4 = gates4.size();
    std::vector<std::vector<unsigned int>> w4m(m, std::vector<unsigned int>(ng4, 0));
    for (unsigned int k = 0; k < m; ++k) {
      for (unsigned int j = 0; j < 2 * m - 1; ++j) {
        const auto c = rq[j][k];
        if (c == 0) {
          continue;
        }
        for (unsigned int g = 0; g < ng4; ++g) {
          w4m[k][g] ^= f4_mul(c, w4[j][g]);
        }
      }
    }

    const std::vector<unsigned int> one = [&] {
      std::vector<unsigned int> v(m, 0);
      v[0] = 1;
      return v;
    }();

    auto tmul = [&](const auto& a, const auto& b) { return f4_pmod(f4_pmul(a, b), tower_modulus); };

    std::optional<std::vector<unsigned int>> rho;
    std::vector<unsigned int> cand;
    f4_tuples_rec(m, cand, [&](const auto& c) {
      if (rho.has_value()) {
        return;
      }
      std::vector<unsigned int> acc(m, 0);
      for (int i = d; i >= 0; --i) {
        acc = tmul(acc, c);
        if (test_bit(q, i)) {
          for (unsigned int j = 0; j < m; ++j) {
            acc[j] ^= one[j];
          }
        }
      }
      if (std::all_of(acc.begin(), acc.end(), [](auto v) { return v == 0; })) {
        rho = c;
      }
    });
    check(rho.has_value(), "no root of q in F4 tower field");

    std::vector<Int> cols;
    auto power = one;
    for (unsigned int i = 0; i < d; ++i) {
      cols.push_back(t_to_bits(power));
      power = tmul(power, *rho);
    }

    std::vector<Int> mrows(d, 0);
    for (unsigned int r = 0; r < d; ++r) {
      Int row = 0;
      for (unsigned int i = 0; i < d; ++i) {
        if (test_bit(cols[i], r)) {
          row |= bit(i);
        }
      }
      mrows[r] = row;
    }
    const auto minv = left_inverse(mrows, d);

    const auto red_a = reduction_rows(q, na);
    const auto red_b = reduction_rows(q, nb);
    std::vector<Int> trow_a(d);
    std::vector<Int> trow_b(d);
    for (unsigned int r = 0; r < d; ++r) {
      trow_a[r] = combine(mrows[r], red_a);
      trow_b[r] = combine(mrows[r], red_b);
    }

    auto form_bits = [&](const auto& l, const auto& trow) {
      Int p = 0;
      Int r = 0;
      for (unsigned int j = 0; j < l.size(); ++j) {
        const auto c = l[j];
        if (c == 0) {
          continue;
        }
        const auto& pj = trow[2 * j];
        const auto& rj = trow[2 * j + 1];
        if ((c & 1) != 0) {
          p ^= pj;
          r ^= rj;
        }
        if ((c & 2) != 0) {
          p ^= rj;
          r ^= pj ^ rj;
        }
      }
      return std::pair<Int, Int>{p, r};
    };

    Alg alg;
    alg.na = na;
    alg.nb = nb;
    std::vector<std::pair<Int, Int>> comp;
    for (const auto& l : gates4) {
      const auto [pa, ra]  = form_bits(l, trow_a);
      const auto [pb, rb]  = form_bits(l, trow_b);
      const unsigned int b = alg.gates.size();
      alg.gates.push_back({pa, pb});
      alg.gates.push_back({ra, rb});
      alg.gates.push_back({pa ^ ra, pb ^ rb});
      comp.push_back({bit(b) | bit(b + 1), bit(b) | bit(b + 2)});
    }

    std::vector<Int> towout(d, 0);
    for (unsigned int k = 0; k < m; ++k) {
      Int o0 = 0;
      Int o1 = 0;
      for (unsigned int g = 0; g < ng4; ++g) {
        const auto c = w4m[k][g];
        if (c == 0) {
          continue;
        }
        const auto [g0, g1] = comp[g];
        if ((c & 1) != 0) {
          o0 ^= g0;
          o1 ^= g1;
        }
        if ((c & 2) != 0) {
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

  const std::array<Int, 13> MONT5_GATES = {
      Int(0b11111), Int(0b11101), Int(0b10111), Int(0b11011), Int(0b01101),
      Int(0b10110), Int(0b11000), Int(0b00011), Int(0b10001), Int(0b10000),
      Int(0b01000), Int(0b00010), Int(0b00001),
  };

  const std::array<std::vector<unsigned int>, 9> MONT5_W = {
      std::vector<unsigned int>{12},
      std::vector<unsigned int>{7, 11, 12},
      std::vector<unsigned int>{2, 5, 7, 8, 9, 12},
      std::vector<unsigned int>{0, 1, 3, 6, 8, 9},
      std::vector<unsigned int>{0, 4, 5, 6, 7, 9, 10, 11, 12},
      std::vector<unsigned int>{0, 2, 3, 7, 8, 12},
      std::vector<unsigned int>{1, 4, 6, 8, 9, 12},
      std::vector<unsigned int>{6, 9, 10},
      std::vector<unsigned int>{9},
  };

  Alg mont5_alg() {
    Alg alg;
    alg.na = 5;
    alg.nb = 5;
    for (const Int& mask : MONT5_GATES) {
      alg.gates.push_back({mask, mask});
    }
    for (const auto& row : MONT5_W) {
      Int w = 0;
      for (auto g : row) {
        w |= bit(g);
      }
      alg.w.push_back(w);
    }
    return alg;
  }

  unsigned int place_gate_cost(const Place& pk);
  CostPlan short_cost(unsigned int n);

  struct CtrSearchOption {
    unsigned int degree = 0;
    unsigned int cost   = 0;
    std::vector<Place> picks;
  };

  std::optional<std::pair<unsigned int, std::vector<Place>>>
  crt_search(unsigned int target, unsigned int maxe_lin, unsigned int cap_irr,
             const std::map<unsigned int, std::vector<Int>>& irr_pools) {
    std::vector<std::vector<CtrSearchOption>> groups;

    for (const auto& spec : {std::pair<Place::Kind, Int>{Place::Kind::Inf, 0},
                             std::pair<Place::Kind, Int>{Place::Kind::Lin, 0b10},
                             std::pair<Place::Kind, Int>{Place::Kind::Lin, 0b11}}) {
      std::vector<CtrSearchOption> opts;
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
        const Int p = pool[0];
        std::vector<CtrSearchOption> opts;
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
        const auto maxt = std::min<std::size_t>(pool.size(), target / d + 1);
        std::vector<CtrSearchOption> opts;
        for (unsigned int t = 1; t <= maxt; ++t) {
          std::vector<Place> picks;
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

    std::map<unsigned int, std::pair<unsigned int, std::vector<Place>>> dp;
    dp[0] = {0, {}};
    for (const auto& opts : groups) {
      auto ndp = dp;
      for (const auto& [deg0, state] : dp) {
        const auto& [cost0, picks0] = state;
        for (const auto& opt : opts) {
          const auto nd = std::min(deg0 + opt.degree, target);
          const auto nc = cost0 + opt.cost;
          auto it       = ndp.find(nd);
          if (it == ndp.end() || nc < it->second.first) {
            std::vector<Place> picks = picks0;
            picks.insert(picks.end(), opt.picks.begin(), opt.picks.end());
            ndp[nd] = {nc, std::move(picks)};
          }
        }
      }
      dp = std::move(ndp);
    }

    auto it = dp.find(target);
    if (it == dp.end()) {
      return std::nullopt;
    }
    return it->second;
  }

  CostPlan full_cost(unsigned int n) {
    static std::map<unsigned int, CostPlan> cache;
    auto it = cache.find(n);
    if (it != cache.end()) {
      return it->second;
    }

    std::optional<CostPlan> best;
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
      std::map<unsigned int, std::vector<Int>> pools;
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
    cache[n] = *best;
    return *best;
  }

  CostPlan short_cost(unsigned int n) {
    static std::map<unsigned int, CostPlan> cache;
    auto it = cache.find(n);
    if (it != cache.end()) {
      return it->second;
    }

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
    cache[n] = best;
    return best;
  }

  unsigned int place_gate_cost(const Place& pk) {
    if (pk.kind == Place::Kind::Inf || pk.kind == Place::Kind::Lin) {
      return short_cost(pk.e).gates;
    }
    const unsigned int d          = pdeg(pk.poly).value_or(0);
    const unsigned int tower_cost = tower_field_cost(d);
    if (pk.e == 1 && tower_cost != 0) {
      return std::min(tower_cost, full_cost(d).gates);
    }
    return full_cost(d * pk.e).gates;
  }

  Alg full_alg(unsigned int n) {
    static std::map<unsigned int, Alg> cache;
    auto it = cache.find(n);
    if (it != cache.end()) {
      return it->second;
    }

    const Plan plan = full_cost(n).plan;
    Alg alg;
    if (plan.kind == Plan::Kind::Base) {
      alg = Alg{1, 1, {{1, 1}}, {1}};
    } else if (plan.kind == Plan::Kind::Mont5) {
      alg = mont5_alg();
    } else if (plan.kind == Plan::Kind::Split) {
      const auto s           = plan.s;
      const auto m           = plan.m;
      const auto oa          = full_alg(s);
      const auto ia          = full_alg(m);
      const unsigned int ngi = ia.gates.size();

      std::vector<std::pair<Int, Int>> gates;
      for (const auto& [fo, go] : oa.gates) {
        std::vector<Int> ar(m, 0);
        std::vector<Int> br(m, 0);
        for (unsigned int i = 0; i < m; ++i) {
          Int fa = 0;
          Int fb = 0;
          for_each_bit(fo, [&](unsigned long t) {
            if (t * m + i < n) {
              fa |= bit(t * m + i);
            }
          });
          for_each_bit(go, [&](unsigned long t) {
            if (t * m + i < n) {
              fb |= bit(t * m + i);
            }
          });
          ar[i] = fa;
          br[i] = fb;
        }
        for (const auto& [fi, gi] : ia.gates) {
          gates.push_back({combine(fi, ar), combine(gi, br)});
        }
      }

      std::vector<Int> w;
      for (unsigned int q = 0; q < 2 * n - 1; ++q) {
        Int row = 0;
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
      alg = Alg{n, n, std::move(gates), std::move(w)};
    } else if (plan.kind == Plan::Kind::Crt) {
      std::vector<std::pair<Int, Int>> gates;
      std::vector<Int> resmasks;
      std::vector<Int> evrows;
      for (const Place& pk : plan.picks) {
        const auto [palg, ev]  = place_alg(pk, n, n);
        const unsigned int off = gates.size();
        gates.insert(gates.end(), palg.gates.begin(), palg.gates.end());
        for (const Int& row : palg.w) {
          resmasks.push_back(row << off);
        }
        evrows.insert(evrows.end(), ev.begin(), ev.end());
      }
      const auto x = left_inverse(evrows, 2 * n - 1);
      std::vector<Int> w;
      for (unsigned int q = 0; q < 2 * n - 1; ++q) {
        w.push_back(combine(x[q], resmasks));
      }
      alg = Alg{n, n, std::move(gates), std::move(w)};
    } else {
      fail("invalid full_alg plan");
    }

    cache[n] = alg;
    return alg;
  }

  Alg short_alg(unsigned int n) {
    static std::map<unsigned int, Alg> cache;
    auto it = cache.find(n);
    if (it != cache.end()) {
      return it->second;
    }

    const Plan plan = short_cost(n).plan;
    Alg alg;
    if (plan.kind == Plan::Kind::Base) {
      alg = Alg{1, 1, {{1, 1}}, {1}};
    } else if (plan.kind == Plan::Kind::Full) {
      const Alg a = full_alg(n);
      alg.na      = n;
      alg.nb      = n;
      alg.gates   = a.gates;
      alg.w.assign(a.w.begin(), a.w.begin() + n);
    } else if (plan.kind == Plan::Kind::ShortSplit) {
      const auto m = plan.m;
      std::vector<Int> lo;
      for (unsigned int i = 0; i < m; ++i) {
        lo.push_back(bit(i));
      }
      const Alg g1 = precompose(full_alg(m), lo, lo, n, n);
      const Alg s  = short_alg(n - m);
      std::vector<Int> lo2;
      std::vector<Int> hi2;
      for (unsigned int i = 0; i < n - m; ++i) {
        lo2.push_back(bit(i));
        hi2.push_back(bit(m + i));
      }
      const Alg g2          = precompose(s, lo2, hi2, n, n);
      const Alg g3          = precompose(s, hi2, lo2, n, n);
      const unsigned int o2 = g1.gates.size();
      const unsigned int o3 = o2 + g2.gates.size();

      std::vector<std::pair<Int, Int>> gates = g1.gates;
      gates.insert(gates.end(), g2.gates.begin(), g2.gates.end());
      gates.insert(gates.end(), g3.gates.begin(), g3.gates.end());
      std::vector<Int> w;
      for (unsigned int j = 0; j < n; ++j) {
        Int row = 0;
        if (j <= 2 * m - 2) {
          row ^= g1.w[j];
        }
        if (j >= m) {
          const auto t = j - m;
          row ^= (g2.w[t] << o2) ^ (g3.w[t] << o3);
        }
        w.push_back(row);
      }
      alg = Alg{n, n, std::move(gates), std::move(w)};
    } else {
      fail("invalid short_alg plan");
    }

    cache[n] = alg;
    return alg;
  }

  std::pair<Alg, std::vector<Int>> place_alg(const Place& pk, unsigned int na, unsigned int nb) {
    const auto n2 = na + nb - 1;
    if (pk.kind == Place::Kind::Inf) {
      check(pk.e <= std::min(na, nb), "infinity multiplicity exceeds input width");
      std::vector<Int> rows_a;
      std::vector<Int> rows_b;
      std::vector<Int> ev;
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
      std::vector<Int> ina;
      std::vector<Int> inb;
      std::vector<Int> ev;
      for (unsigned int i = 0; i < pk.e; ++i) {
        ina.push_back(combine(p[i], ra));
        inb.push_back(combine(p[i], rb));
        ev.push_back(combine(p[i], rev));
      }
      return {precompose(short_alg(pk.e), ina, inb, na, nb), ev};
    }

    const auto tower_cost = tower_field_cost(pdeg(pk.poly).value_or(0));
    if (pk.e == 1 && tower_cost != 0 && tower_cost < full_cost(pdeg(pk.poly).value_or(0)).gates) {
      return tower_place_alg(pk.poly, na, nb);
    }

    const Int q     = poly_pow(pk.poly, pk.e);
    const auto d    = pdeg(q);
    const Alg full  = full_alg(*d);
    Alg alg         = precompose(full, reduction_rows(q, na), reduction_rows(q, nb), na, nb);
    const auto rout = reduction_rows(q, 2 * *d - 1);
    std::vector<Int> w;
    for (unsigned int i = 0; (int)i < d; ++i) {
      w.push_back(combine(rout[i], alg.w));
    }
    alg.w = std::move(w);
    return {alg, reduction_rows(q, n2)};
  }

  struct Tables {
    std::vector<Int> f;
    std::vector<Int> g;
    std::vector<Int> w_tree;
    std::vector<Int> w_gate;
    unsigned int n_tree = 0;
    std::vector<std::pair<Place, unsigned int>> report;
  };

  Tables build_tables(unsigned int lam, unsigned int wgrind, const Int& p,
                      const std::vector<Int>& tree_moduli, const std::vector<Place>& picks) {
    const unsigned int na = lam;
    const unsigned int nb = lam - wgrind;
    unsigned int ntree    = 0;
    for (const Int& m : tree_moduli) {
      ntree += pdeg(m).value_or(0);
    }
    check(ntree == nb, "tree degrees must sum to lambda - wgrind");

    const auto n2 = na + nb - 1;
    std::vector<Int> evrows;
    for (const Int& m : tree_moduli) {
      const auto rows = reduction_rows(m, n2);
      evrows.insert(evrows.end(), rows.begin(), rows.end());
    }

    std::vector<std::pair<Int, Int>> gates;
    std::vector<Int> resmasks;
    std::vector<std::pair<Place, unsigned int>> report;
    for (const Place& pk : picks) {
      const auto [palg, ev]  = place_alg(pk, na, nb);
      const unsigned int off = gates.size();
      gates.insert(gates.end(), palg.gates.begin(), palg.gates.end());
      for (const Int& row : palg.w) {
        resmasks.push_back(row << off);
      }
      evrows.insert(evrows.end(), ev.begin(), ev.end());
      report.push_back({pk, palg.gates.size()});
    }

    const auto x   = left_inverse(evrows, n2);
    const auto red = reduction_rows(p, n2);
    std::vector<Int> wt(lam, 0);
    std::vector<Int> wg(lam, 0);
    for (unsigned int r = 0; r < lam; ++r) {
      const Int sel = combine(red[r], x);
      for_each_bit(sel, [&](unsigned long idx) {
        if (idx < ntree) {
          wt[r] |= bit(idx);
        } else {
          wg[r] ^= resmasks[idx - ntree];
        }
      });
    }

    std::vector<Int> f;
    std::vector<Int> g;
    f.reserve(gates.size());
    g.reserve(gates.size());
    const auto cols = crt_lift_cols(tree_moduli);
    for (const auto& [fa, fb] : gates) {
      f.push_back(fa);
      Int row = 0;
      for (unsigned int j = 0; j < nb; ++j) {
        if (parity(fb & cols[j]) != 0) {
          row |= bit(j);
        }
      }
      g.push_back(row);
    }
    return {std::move(f), std::move(g), std::move(wt), std::move(wg), ntree, std::move(report)};
  }

  std::vector<Int> wcrt_rows(const std::vector<Int>& tree_moduli, unsigned int lam) {
    const auto cols = crt_lift_cols(tree_moduli);
    std::vector<Int> rows;
    rows.reserve(lam);
    for (unsigned int r = 0; r < lam; ++r) {
      Int row = 0;
      for (unsigned int j = 0; j < cols.size(); ++j) {
        if (test_bit(cols[j], r)) {
          row |= bit(j);
        }
      }
      rows.push_back(row);
    }
    return rows;
  }

  void prune(std::vector<Int>& f, std::vector<Int>& g, std::vector<Int>& wg) {
    std::map<std::pair<Int, Int>, int> canon;
    for (unsigned int i = 0; i < f.size(); ++i) {
      if (f[i] != 0 && g[i] != 0) {
        const auto key = std::make_pair(f[i], g[i]);
        if (canon.find(key) == canon.end()) {
          canon[key] = i;
        }
      }
    }

    std::vector<Int> wg1;
    wg1.reserve(wg.size());
    for (const Int& row : wg) {
      Int nr = 0;
      for_each_bit(row, [&](unsigned long i) {
        if (f[i] != 0 && g[i] != 0) {
          nr ^= bit(canon[std::make_pair(f[i], g[i])]);
        }
      });
      wg1.push_back(nr);
    }

    Int used = 0;
    for (const Int& row : wg1) {
      used |= row;
    }

    std::map<unsigned int, unsigned int> newidx;
    std::vector<Int> f2;
    std::vector<Int> g2;
    for (unsigned int i = 0; i < f.size(); ++i) {
      if (test_bit(used, i)) {
        newidx[i] = f2.size();
        f2.push_back(f[i]);
        g2.push_back(g[i]);
      }
    }

    std::vector<Int> wg2;
    wg2.reserve(wg1.size());
    for (const Int& row : wg1) {
      Int nr = 0;
      for_each_bit(row, [&](int i) { nr |= bit(newidx[i]); });
      wg2.push_back(nr);
    }

    f  = std::move(f2);
    g  = std::move(g2);
    wg = std::move(wg2);
  }

  std::string to_hex(Int v) {
    if (v == 0) {
      return "0";
    }
    static constexpr char digits[] = "0123456789abcdef";
    std::string out;
    while (v != 0) {
      out.push_back(digits[(v & 0xf).convert_to<unsigned>()]);
      v >>= 4;
    }
    std::reverse(out.begin(), out.end());
    return out;
  }

  std::string word_hex(uint64_t v) {
    std::ostringstream ss;
    ss << "0x" << std::hex << std::nouppercase << std::setfill('0') << std::setw(16) << v << "ULL";
    return ss.str();
  }

  uint64_t word_at(const Int& v, unsigned int word) {
    return ((v >> (64 * word)) & ((Int(1) << 64) - 1)).convert_to<uint64_t>();
  }

  constexpr unsigned int words_of(unsigned int width) {
    return (width + 63) / 64;
  }

  std::string upper_identifier(std::string s) {
    for (char& c : s) {
      if (c == '-') {
        c = '_';
      } else {
        c = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
      }
    }
    return s;
  }

  std::string prefix_for_name(const std::string& name) {
    return "FAEST_" + upper_identifier(name);
  }

  std::string row_literal(const Int& r, unsigned int words) {
    std::ostringstream ss;
    ss << "{ ";
    for (unsigned int w = 0; w < words; ++w) {
      if (w != 0) {
        ss << ", ";
      }
      ss << word_hex(word_at(r, w));
    }
    ss << " }";
    return ss.str();
  }

  struct EmitTable {
    std::string suffix;
    const std::vector<Int>& rows;
    int words = 0;
    std::string rows_macro;
    std::string words_macro;
  };

  std::string emit_c(const std::string& path, const std::string& name, const std::vector<Int>& f,
                     const std::vector<Int>& g, const std::vector<Int>& wt,
                     const std::vector<Int>& wg, const std::vector<Int>& tree_moduli,
                     unsigned int lam, unsigned int wgrind, unsigned int ntree) {
    const std::vector<Int> wcrt = wcrt_rows(tree_moduli, lam);
    const Int m_tree            = mtree_poly(tree_moduli);
    const unsigned int ng       = f.size();
    const unsigned int nb       = lam - wgrind;
    const unsigned int tau      = tree_moduli.size();
    for (const Int& m : tree_moduli) {
      check(pdeg(m).has_value() && pdeg(m).value() < 64,
            "tree modulus degree >= 64 does not fit one uint64 word");
    }

    const int w_f  = words_of(lam);
    const int w_g  = words_of(nb);
    const int w_wt = words_of(ntree);
    const int w_wg = words_of(ng);
    const int w_wc = words_of(nb);
    const int w_mt = words_of(*pdeg(m_tree) + 1);

    const std::string pfx   = prefix_for_name(name);
    const std::string guard = pfx + "_TABLES_H";
    std::string src_path;
    if (path.size() >= 2 && path.substr(path.size() - 2) == ".h") {
      src_path = path.substr(0, path.size() - 2) + ".c";
    } else {
      src_path = path + ".c";
    }
    const std::string header_name = std::filesystem::path(path).filename().string();

    const std::array<EmitTable, 5> tables = {{
        {"F", f, w_f, pfx + "_NGATES", pfx + "_F_WORDS"},
        {"G", g, w_g, pfx + "_NGATES", pfx + "_G_WORDS"},
        {"W_TREE", wt, w_wt, pfx + "_LAMBDA", pfx + "_W_TREE_WORDS"},
        {"W_GATE", wg, w_wg, pfx + "_LAMBDA", pfx + "_W_GATE_WORDS"},
        {"W_CRT", wcrt, w_wc, pfx + "_LAMBDA", pfx + "_W_CRT_WORDS"},
    }};

    std::ofstream out(path);
    check(out.good(), "failed to open header for writing: " + path);
    out << "/* generated by vole_mult_tables.cpp: " << name << "\n"
        << " *\n"
        << " * Bilinear VOLE-multiplication tables.\n"
        << " *\n"
        << " * Encoding: each table row is an F_2 bit-vector packed\n"
        << " * little-endian into uint64 words -- bit i sits in word i/64\n"
        << " * at bit position i%64, and word 0 holds the low bits.\n"
        << " */\n";
    out << "#ifndef " << guard << "\n#define " << guard << "\n\n";
    out << "#include <stdint.h>\n\n";
    out << "#define " << pfx << "_NGATES " << ng << "\n";
    out << "#define " << pfx << "_WGRIND " << wgrind << "\n";
    out << "#define " << pfx << "_NDELTA_BITS " << nb << "\n";
    out << "#define " << pfx << "_NTREE_BITS " << ntree << "\n";
    out << "#define " << pfx << "_F_WORDS " << w_f << "\n";
    out << "#define " << pfx << "_G_WORDS " << w_g << "\n";
    out << "#define " << pfx << "_W_TREE_WORDS " << w_wt << "\n";
    out << "#define " << pfx << "_W_GATE_WORDS " << w_wg << "\n";
    out << "#define " << pfx << "_W_CRT_WORDS " << w_wc << "\n";
    out << "#define " << pfx << "_M_TREE_WORDS " << w_mt << "\n\n";

    for (const auto& tbl : tables) {
      out << "static const uint64_t " << pfx << "_" << tbl.suffix << "[" << tbl.rows_macro << "]["
          << tbl.words_macro << "] = {\n";
      for (const Int& row : tbl.rows) {
        out << "  " << row_literal(row, tbl.words) << ",\n";
      }
      out << "};\n\n";
    }

    out << "static const uint64_t " << pfx << "_TREE_MODULI[" << pfx << "_TAU] = { ";
    for (unsigned int i = 0; i < tau; ++i) {
      if (i != 0) {
        out << ", ";
      }
      out << "0x" << to_hex(tree_moduli[i]) << "ULL";
    }
    out << " };\n\n";

    out << "static const uint64_t " << pfx << "_M_TREE[" << pfx << "_M_TREE_WORDS] = { ";
    for (int w = 0; w < w_mt; ++w) {
      if (w != 0) {
        out << ", ";
      }
      out << word_hex(word_at(m_tree, w));
    }
    out << " };\n\n";
    out << "#endif\n";

    return src_path;
  }

  struct Preset {
    unsigned int lam    = 0;
    unsigned int tau    = 0;
    unsigned int wgrind = 0;
  };

  const std::map<std::string, Preset> PRESETS = {
      {"128s", {128, 11, 7}},    {"128f", {128, 16, 8}}, {"192s", {192, 16, 12}},
      {"192f", {192, 24, 8}},    {"256s", {256, 22, 6}}, {"256f", {256, 32, 8}},
      {"em_192s", {192, 16, 8}},
  };

  std::vector<std::pair<unsigned int, unsigned int>>
  faest_tree_spec(unsigned int lam, unsigned int tau, unsigned int wgrind) {
    const unsigned int n    = lam - wgrind;
    const unsigned int d1   = n / tau + 1;
    const unsigned int tau1 = n % tau;
    std::vector<std::pair<unsigned int, unsigned int>> spec;
    if (tau1 != 0) {
      spec.push_back({d1, tau1});
    }
    spec.push_back({d1 - 1, tau - tau1});
    return spec;
  }

  struct Args {
    int maxd = 10;
    int maxe = 10;
    std::string preset;
    std::string emit_c_path;
  };

  void print_help(const char* prog) {
    std::cout << "Usage:\n"
              << "  " << prog << " --preset NAME --emit-c FILE.h\n\n"
              << "Search + frozen C table generation for CRT-based F_{2^lambda} VOLE "
                 "multiplication.\n\n"
              << "Options:\n"
              << "  --preset NAME         one of 128s, 128f, 192s, 192f, 256s, 256f, 192s_em\n"
              << "  --emit-c FILE.h       emit C header and companion source\n";
  }

  Args parse_args(int argc, char** argv) {
    Args args;
    for (int i = 1; i < argc; ++i) {
      const std::string a = argv[i];
      auto need_value     = [&](const std::string& opt) {
        if (i + 1 >= argc) {
          fail(opt + " requires a value");
        }
        return std::string(argv[++i]);
      };

      if (a == "--help" || a == "-h") {
        print_help(argv[0]);
        std::exit(0);
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

  void run_set(const std::string& name, unsigned int lam, unsigned int wgrind,
               const std::vector<std::pair<unsigned int, unsigned int>>& tree_spec,
               const Args& args, const std::string& c_path) {
    unsigned int tree_sum = 0;
    for (const auto& [d, c] : tree_spec) {
      tree_sum += d * c;
    }
    check(tree_sum == lam - wgrind,
          "tree degrees sum to " + std::to_string(tree_sum) +
              ", expected lambda - wgrind = " + std::to_string(lam - wgrind));

    const unsigned int n2      = lam + (lam - wgrind) - 1;
    const unsigned int deficit = n2 - tree_sum;

    std::vector<Int> tree_moduli;
    for (const auto& [d, c] : tree_spec) {
      const auto& pool = irreducibles_of_degree(d);
      if (c > pool.size()) {
        fail("only " + std::to_string(pool.size()) + " irreducibles of degree " +
             std::to_string(d) + " over F_2; cannot pick " + std::to_string(c));
      }
      tree_moduli.insert(tree_moduli.end(), pool.begin(), pool.begin() + c);
    }
    std::set<Int> treeset(tree_moduli.begin(), tree_moduli.end());
    std::map<unsigned int, std::vector<Int>> pools;
    for (int d = 2; d <= args.maxd; ++d) {
      for (const Int& q : irreducibles_of_degree(d)) {
        if (treeset.find(q) == treeset.end()) {
          pools[d].push_back(q);
        }
      }
    }

    const Int p = faest_modulus(lam);
    auto res    = crt_search(deficit, args.maxe, args.maxd, pools);
    check(res.has_value(), "portfolio search failed; raise --maxd/--maxe");
    const std::vector<Place>& picks = res->second;

    Tables tables = build_tables(lam, wgrind, p, tree_moduli, picks);
    prune(tables.f, tables.g, tables.w_gate);

    std::map<std::pair<std::string, unsigned int>, std::array<unsigned int, 3>> agg;
    for (const auto& [pk, gates] : tables.report) {
      const auto d =
          pk.kind == Place::Kind::Inf || pk.kind == Place::Kind::Lin ? pk.e : *pdeg(pk.poly) * pk.e;
      std::string type;
      if (pk.kind == Place::Kind::Irr) {
        type = "irr deg " + std::to_string(*pdeg(pk.poly));
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

    const std::vector<Int> wcrt = wcrt_rows(tree_moduli, lam);

    emit_c(c_path, name, tables.f, tables.g, tables.w_tree, tables.w_gate, tree_moduli, lam, wgrind,
           tables.n_tree);
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
    const auto tree_spec = faest_tree_spec(it->second.lam, it->second.tau, it->second.wgrind);
    run_set(args.preset, it->second.lam, it->second.wgrind, tree_spec, args, args.emit_c_path);
    return 0;
  } catch (const std::exception& e) {
    std::cerr << "vole_mult_tables: " << e.what() << "\n";
    return 1;
  }
}
