#include "e_lib.h"
#define assert(x)

#include "modp_data.h"

struct gmp_div_inverse
{
  /* Normalization shift count. */
  unsigned shift;
  /* Normalized divisor (d0 unused for mpn_div_qr_1) */
  mp_limb_t d1, d0;
  /* Inverse, for 2/1 or 3/2. */
  mp_limb_t di;
};

static modp_indata_t inbuf;
static volatile unsigned wait_flag;

#define Q_LEN 7
mp_limb_t primorial[Q_LEN] = { 3990313926u, 1301064672u, 203676164u, 3309914335u, 2220064684u, 2929319840u, 153406374u }; 

#define CHAR_BIT 8
#define GMP_LIMB_BITS (sizeof(mp_limb_t) * CHAR_BIT)

#define GMP_LIMB_MAX (~ (mp_limb_t) 0)
#define GMP_LIMB_HIGHBIT ((mp_limb_t) 1 << (GMP_LIMB_BITS - 1))

#define GMP_HLIMB_BIT ((mp_limb_t) 1 << (GMP_LIMB_BITS / 2))
#define GMP_LLIMB_MASK (GMP_HLIMB_BIT - 1)

#define gmp_clz(count, xx) do {                                         \
    unsigned config = 0x1;                                              \
    unsigned _x = (xx);                                                 \
    unsigned half = 0x3f000000;                                         \
    unsigned _res;                                                      \
    __asm__ __volatile__ (                                              \
      "movts config, %[config]\n\t"                                     \
      "nop\n\t"                                                         \
      "float %[x], %[x]\n\t"                                            \
      "lsl %[config], %[config], #19\n\t"                               \
      "fadd %[x], %[x], %[half]\n\t"                                    \
      "mov %[res], #158\n\t"                                            \
      "lsr %[x], %[x], #23\n\t"                                         \
      "sub %[res], %[res], %[x]\n\t"                                    \
      "movts config, %[config]" :                                       \
      [res] "=r" (_res),                                                \
      [config] "+r" (config), [x] "+r" (_x) :                           \
      [half] "r" (half));                                               \
    (count) = _res;                                                     \
  } while (0)

#define gmp_add_ssaaaa(sh, sl, ah, al, bh, bl) \
  do {                                                                  \
    mp_limb_t __x;                                                      \
    __x = (al) + (bl);                                                  \
    (sh) = (ah) + (bh) + (__x < (al));                                  \
    (sl) = __x;                                                         \
  } while (0)

#define gmp_umul_ppmm(w1, w0, u, v)                                     \
  do {                                                                  \
    mp_limb_t __x0, __x1, __x2, __x3;                                   \
    unsigned __ul, __vl, __uh, __vh;                                    \
    mp_limb_t __u = (u), __v = (v);                                     \
                                                                        \
    __ul = __u & GMP_LLIMB_MASK;                                        \
    __uh = __u >> (GMP_LIMB_BITS / 2);                                  \
    __vl = __v & GMP_LLIMB_MASK;                                        \
    __vh = __v >> (GMP_LIMB_BITS / 2);                                  \
                                                                        \
    __x0 = (mp_limb_t) __ul * __vl;                                     \
    __x1 = (mp_limb_t) __ul * __vh;                                     \
    __x2 = (mp_limb_t) __uh * __vl;                                     \
    __x3 = (mp_limb_t) __uh * __vh;                                     \
                                                                        \
    __x1 += __x0 >> (GMP_LIMB_BITS / 2);/* this can't give carry */     \
    __x1 += __x2;               /* but this indeed can */               \
    if (__x1 < __x2)            /* did we get it? */                    \
      __x3 += GMP_HLIMB_BIT;    /* yes, add it in the proper pos. */    \
                                                                        \
    (w1) = __x3 + (__x1 >> (GMP_LIMB_BITS / 2));                        \
    (w0) = (__x1 << (GMP_LIMB_BITS / 2)) + (__x0 & GMP_LLIMB_MASK);     \
  } while (0)

#define gmp_udiv_rnnd_preinv(r, nh, nl, d, di)                      \
  do {                                                                  \
    mp_limb_t _qh, _ql, _r, _mask;                                      \
    gmp_umul_ppmm (_qh, _ql, (nh), (di));                               \
    gmp_add_ssaaaa (_qh, _ql, _qh, _ql, (nh) + 1, (nl));                \
    _r = (nl) - _qh * (d);                                              \
    _mask = -(mp_limb_t) (_r > _ql); /* both > and >= are OK */         \
    _qh += _mask;                                                       \
    _r += _mask & (d);                                                  \
    if (_r >= (d))                                                      \
      {                                                                 \
        _r -= (d);                                                      \
        _qh++;                                                          \
      }                                                                 \
                                                                        \
    (r) = _r;                                                           \
  } while (0)

#define mpn_invert_limb(x) mpn_invert_3by2 ((x), 0)

// Divide algorithm modified version of the builtin udivsi3
// original source: libgcc/config/epiphany/udivsi3-float.c
unsigned int
udiv (unsigned int a, unsigned int b)
{
  unsigned int d, t, s0, s1, r0, r1;

  if (b > a)
    return 0;
  if (b & 0x80000000)
    return 1;

  /* Assuming B is nonzero, compute S0 such that 0 <= S0,
     (B << S0+1) does not overflow,
     A < 4.01 * (B << S0), with S0 chosen as small as possible
     without taking to much time calculating.  */
#if 1
  {
    unsigned config = 0x1;
    unsigned half = 0x3f000000;
    unsigned allbits = 32;
    unsigned tmp;
    asm volatile (
      "movts config, %[config]\n\t"
      "nop\n\t"
      "float %[s0], %[a]\n\t"
      "float %[s1], %[b]\n\t"
      "lsl %[config], %[config], #19\n\t"
      "lsr %[s0], %[s0], #23\n\t"
      "lsr %[s1], %[s1], #23\n\t"
      "movts config, %[config]\n\t"
      "sub %[s0], %[s0], #126\n\t"
      "  isub %[config], %[config], %[config]\n\t"
      "sub %[s1], %[s1], #125\n\t"
      "sub %[tmp], %[a], #0\n\t"
      "movlt %[s0], %[allbits]\n\t"
      "sub %[s0], %[s0], %[s1]\n\t"
      "movlt %[s0], %[config]\n\t" :
      [config] "+r" (config), [s0] "=r" (s0), [s1] "=r" (s1), [tmp] "=r" (tmp) :
      [half] "r" (half), [allbits] "r" (allbits), [a] "r" (a), [b] "r" (b));
  }
#else
  gmp_clz(s0, a);
  gmp_clz(s1, b);
  s0 = 32 - s0;
  s1 = (32 - s1) + 1;
  if (s0 < s1) s0 = 0;
  else s0 -= s1;
#endif

  b <<= s0;
  r1 = 0;

  r0 = 1 << s0;
  a = ((t=a) - b);
  if (a <= t)
    {
      r1 += r0;
      a = ((t=a) - b);
      if (a <= t)
      do {
        r1 += r0;
        a = ((t=a) - b);
      } while (a <= t);
    }
  a += b;
  d = b - 1;

#define STEP(n) case n: a += a; t = a - d; if (t <= a) a = t;
  switch (s0)
    {
    STEP (31)
    STEP (30)
    STEP (29)
    STEP (28)
    STEP (27)
    STEP (26)
    STEP (25)
    STEP (24)
    STEP (23)
    STEP (22)
    STEP (21)
    STEP (20)
    STEP (19)
    STEP (18)
    STEP (17)
    STEP (16)
    STEP (15)
    STEP (14)
    STEP (13)
    STEP (12)
    STEP (11)
    STEP (10)
    STEP (9)
    STEP (8)
    STEP (7)
    STEP (6)
    STEP (5)
    STEP (4)
    STEP (3)
    STEP (2)
    STEP (1)
    case 0: ;
    }
#undef STEP
  r0 = r1 | ((r0-1) & a);
  return r0;
}

static mp_limb_t
mpn_invert_3by2 (mp_limb_t u1, mp_limb_t u0)
{
  mp_limb_t r, p, m;
  unsigned ul, uh;
  unsigned ql, qh;

  /* First, do a 2/1 inverse. */
  /* The inverse m is defined as floor( (B^2 - 1 - u1)/u1 ), so that 0 <
   * B^2 - (B + m) u1 <= u1 */
  assert (u1 >= GMP_LIMB_HIGHBIT);

  ul = u1 & GMP_LLIMB_MASK;
  uh = u1 >> (GMP_LIMB_BITS / 2);

  qh = udiv(~u1, uh);
  r = ((~u1 - (mp_limb_t) qh * uh) << (GMP_LIMB_BITS / 2)) | GMP_LLIMB_MASK;

  p = (mp_limb_t) qh * ul;
  /* Adjustment steps taken from udiv_qrnnd_c */
  if (r < p)
    {
      qh--;
      r += u1;
      if (r >= u1) /* i.e. we didn't get carry when adding to r */
        if (r < p)
          {
            qh--;
            r += u1;
          }
    }
  r -= p;

  /* Do a 3/2 division (with half limb size) */
  p = (r >> (GMP_LIMB_BITS / 2)) * qh + r;
  ql = (p >> (GMP_LIMB_BITS / 2)) + 1;

  /* By the 3/2 method, we don't need the high half limb. */
  r = (r << (GMP_LIMB_BITS / 2)) + GMP_LLIMB_MASK - ql * u1;

  if (r >= (p << (GMP_LIMB_BITS / 2)))
    {
      ql--;
      r += u1;
    }
  m = ((mp_limb_t) qh << (GMP_LIMB_BITS / 2)) + ql;
  if (r >= u1)
    {
      m++;
      r -= u1;
    }

  if (u0 > 0)
    {
      mp_limb_t th, tl;
      r = ~r;
      r += u0;
      if (r < u0)
        {
          m--;
          if (r >= u1)
            {
              m--;
              r -= u1;
            }
          r -= u1;
        }
      gmp_umul_ppmm (th, tl, u0, m);
      r += th;
      if (r < th)
        {
          m--;
          m -= ((r > u1) | ((r == u1) & (tl > u0)));
        }
    }

  return m;
}

static mp_limb_t
mpn_lshift (mp_ptr rp, mp_srcptr up, mp_size_t n, unsigned int cnt)
{
  mp_limb_t high_limb, low_limb;
  unsigned int tnc, cntmul;
  mp_size_t i;
  mp_limb_t retval;

  assert (n >= 1);
  assert (cnt >= 1);
  assert (cnt < GMP_LIMB_BITS);

  up += n;
  rp += n;

  tnc = GMP_LIMB_BITS - cnt;
  cntmul = 1<<cnt;
  low_limb = *--up;
  retval = low_limb >> tnc;
  high_limb = (low_limb * cntmul);

  for (i = n; --i != 0;)
    {
      low_limb = *--up;
      *--rp = high_limb | (low_limb >> tnc);
      high_limb = (low_limb * cntmul);
    }
  *--rp = high_limb;

  return retval;
}


static mp_limb_t
mpn_div_r_1_preinv_ns(mp_srcptr np, mp_size_t nn,
                      const struct gmp_div_inverse *inv)
{
  mp_limb_t d, di;
  mp_limb_t r = 0;

  d = inv->d1;
  di = inv->di;
  while (nn-- > 0)
    {
      gmp_udiv_rnnd_preinv (r, r, np[nn], d, di);
    }

  return r >> inv->shift;
}

static void
mpn_div_r_1_preinv_ns_2(mp_limb_t* rp1, mp_limb_t* rp2,
                        mp_srcptr np, mp_size_t nn,
                        const struct gmp_div_inverse *inv1,
                        const struct gmp_div_inverse *inv2)
{
  mp_limb_t d1, di1, d2, di2;
  mp_limb_t r1 = 0, r2 = 0;

  assert(inv1->shift == inv2->shift);

  d1 = inv1->d1;
  di1 = inv1->di;
  d2 = inv2->d1;
  di2 = inv2->di;
  while (nn-- > 0)
    {
      gmp_udiv_rnnd_preinv (r1, r1, np[nn], d1, di1);
      gmp_udiv_rnnd_preinv (r2, r2, np[nn], d2, di2);
    }

  *rp1 = r1 >> inv1->shift;
  *rp2 = r2 >> inv2->shift;
}

	  
static void
mpn_div_qr_1_invert (struct gmp_div_inverse *inv, mp_limb_t d)
{
  unsigned shift;

  assert (d > 0);
  gmp_clz (shift, d);
  inv->shift = shift;
  inv->d1 = d << shift;
  inv->di = mpn_invert_limb (inv->d1);
}


// return t such that at = 1 mod b
// a, b < 2^31.
// Algorithm from http://www.ucl.ac.uk/~ucahcjm/combopt/ext_gcd_python_programs.pdf
// with modification to keep s, t in bounds.
static unsigned int inverse(unsigned int a, unsigned int b)
{
  int alpha, beta;
  int u, v, s, t;
  u = 1; v = 0; s = 0; t = 1;
  alpha = a; beta = b;

  if (a == 0)
    return 0;

  // Keep a = u * alpha + v * beta
  while ((a&1) == 0)
  {
    a >>= 1;
    if ((u|v) & 1)
    {
      u = (u + beta) >> 1;
      v = (v - alpha) >> 1;
    }
    else
    {
      u >>= 1;
      v >>= 1;
    }
  }
  while (a!=b)
  {
    if ((b&1)==0)
    {
      b >>= 1;
      if ((s|t) & 1)
      {
        s = (s + beta) >> 1;
        t = (t - alpha) >> 1;
      }
      else
      {
        s >>= 1;
        t >>= 1;
      }
    }
    else if (b < a)
    {
      int tmp;
      tmp = a;
      a = b;
      b = tmp;
      tmp = u;
      u = s;
      s = tmp;
      tmp = v;
      v = t;
      t = tmp;
    }
    else
    {
      b = b - a;
      s = s - u;
      t = t - v;

      // Fix up to avoid out of range s,t
      if (u > 0)
      {
        if (s <= -beta)
        {
          s += beta;
          t -= alpha;
        }
      }
      else
      {
        if (s >= beta)
        {
          s -= beta;
          t += alpha;
        }
      }
    }
  }
  if (s < 0) s += beta;
  return s;
}

static unsigned mulmod(mp_limb_t a, mp_limb_t b, struct gmp_div_inverse* inv)
{
  mp_limb_t th, tl, r;

  // th:tl = a*b
  gmp_umul_ppmm(th, tl, a, b);

  // Shift to inverse normalized form
  // Can't overflow if a, b < p.
  th = (th << inv->shift) | (tl >> (GMP_LIMB_BITS - inv->shift));
  tl <<= inv->shift;

  // Divide
  {
    mp_limb_t d1, di;
    d1 = inv->d1;
    di = inv->di;
    gmp_udiv_rnnd_preinv(r, th, tl, d1, di);
  }

  return r >> inv->shift;
}

void null_isr(int);

static void doshift(mp_ptr nshifted, mp_size_t* nshiftedn, 
                    mp_ptr qshifted, mp_size_t* qshiftedn, unsigned shift)
{
              mp_limb_t h;
              h = mpn_lshift(nshifted, inbuf.n, inbuf.nn, shift);
              if (h) 
              {
                nshifted[inbuf.nn] = h;
                *nshiftedn = inbuf.nn + 1;
              }
              else
                *nshiftedn = inbuf.nn;

              h = mpn_lshift(qshifted, primorial, Q_LEN, shift);
              if (h) 
              {
                qshifted[Q_LEN] = h;
                *qshiftedn = Q_LEN + 1;
              }
              else
                *qshiftedn = Q_LEN;
}

int main()
{
  e_coreid_t coreid;
  coreid = e_get_coreid();
  unsigned int row, col, core;
  e_coords_from_coreid(coreid, &row, &col);
  core = row * 4 + col;

  modp_indata_t* in = (modp_indata_t*)(0x8f000000+(SHARED_MEM_PER_CORE*core));
  modp_outdata_t* outbufs[2];
  outbufs[0] = (modp_outdata_t*)((char*)in + sizeof(modp_indata_t));
  outbufs[1] = (modp_outdata_t*)((char*)outbufs[0] + sizeof(modp_outdata_t));
  int buffer;

  // Must set up a null interrupt routine to allow host to wake us from idle
  e_irq_attach(E_SYNC, null_isr);
  e_irq_mask(E_SYNC, E_FALSE);
  e_irq_global_mask(E_FALSE);
  
  while (1)
  {
    buffer = 0;
    e_dma_copy(&inbuf, in, sizeof(modp_indata_t));

    // Ensure we don't run this block again.
    in->pbase = 0;

    unsigned num_results;
    unsigned i = 0;
    mp_fixed_len_num nshifted, qshifted;
    mp_size_t nshiftedn = 0, qshiftedn = 0;
    unsigned lastshift = 0xffffffff;

    do
    {
      modp_outdata_t* out = outbufs[buffer];
      buffer ^= 1;

      // Wait for buffer ready, without spamming reads to the
      // off core memory.
      while(out->results_status != 0)
      {
        wait_flag = 1;
        while (wait_flag == 1);
      }

      // Writing
      out->results_status = 1;
      num_results = 0;
      mp_limb_t p[2];
      unsigned pidx = 0;
      for (; i < MODP_E_SIEVE_SIZE && num_results + pidx < MODP_RESULTS_PER_PAGE; ++i)
      {
        if ((inbuf.sieve[i>>5] & (1<<(i&0x1f))) == 0)
        {
          p[pidx++] = inbuf.pbase + (i<<1);
          if (pidx == 2)
          {
            struct gmp_div_inverse inv1, inv2;
            mpn_div_qr_1_invert (&inv1, p[0]);
            mpn_div_qr_1_invert (&inv2, p[1]);

            if (inv1.shift != lastshift)
            {
              doshift(nshifted, &nshiftedn, qshifted, &qshiftedn, inv1.shift);
              lastshift = inv1.shift;
            }

            if (inv1.shift == inv2.shift)
            {
              modp_result_t result1, result2;
              mpn_div_r_1_preinv_ns_2(&result1.r, &result2.r, nshifted, nshiftedn, &inv1, &inv2);
              mp_limb_t q1, q2;
              mpn_div_r_1_preinv_ns_2(&q1, &q2, qshifted, qshiftedn, &inv1, &inv2);
#ifdef MODP_RESULT_DEBUG
              result1.p = p[0];
              result2.p = p[1];
              result1.q = q1;
              result2.q = q2;
              result1.x = result1.r;
              result2.x = result2.r;
#endif
              q1 = inverse(q1, p[0]);
              q2 = inverse(q2, p[1]);
              result1.r = mulmod(result1.r, q1, &inv1);
              result2.r = mulmod(result2.r, q2, &inv2);
              q1 <<= 1;
              if (q1 >= p[0]) q1 -= p[0];
              q2 <<= 1;
              if (q2 >= p[1]) q2 -= p[1];
              result1.twoqinv = q1;
              result2.twoqinv = q2;
              out->result[num_results++] = result1;
              out->result[num_results++] = result2;

              pidx = 0;
            }
            else
            {
              modp_result_t result;
              result.r = mpn_div_r_1_preinv_ns(nshifted, nshiftedn, &inv1);
              mp_limb_t q = mpn_div_r_1_preinv_ns(qshifted, qshiftedn, &inv1);
#ifdef MODP_RESULT_DEBUG
              result.p = p[0];
              result.q = q;
              result.x = result.r;
#endif
              q = inverse(q, p[0]);
              result.r = mulmod(result.r, q, &inv1);
              q <<= 1;
              if (q >= p[0]) q -= p[0];
              result.twoqinv = q;
              out->result[num_results++] = result;

              p[0] = p[1];
              pidx = 1;
            }
          }
        }
      }

      if (pidx == 1)
      {
        struct gmp_div_inverse inv;
        mpn_div_qr_1_invert (&inv, p[0]);
        if (lastshift != inv.shift)
        {
          doshift(nshifted, &nshiftedn, qshifted, &qshiftedn, inv.shift);
          lastshift = inv.shift;
        }
        modp_result_t result;
        result.r = mpn_div_r_1_preinv_ns(nshifted, nshiftedn, &inv);
        mp_limb_t q = mpn_div_r_1_preinv_ns(qshifted, qshiftedn, &inv);
#ifdef MODP_RESULT_DEBUG
        result.p = p[0];
        result.q = q;
        result.x = result.r;
#endif
        q = inverse(q, p[0]);
        result.r = mulmod(result.r, q, &inv);
        q <<= 1;
        if (q >= p[0]) q -= p[0];
        result.twoqinv = q;
        out->result[num_results++] = result;
      }

      out->num_results = num_results;
      out->results_status = 2;
    } while(num_results == MODP_RESULTS_PER_PAGE);

    while (in->pbase == 0) __asm__ __volatile__ ("idle");
  }

  return 0;
}

void __attribute__((interrupt)) null_isr(int x) 
{ 
  wait_flag = 0;
  return;
}
