#include "e_lib.h"
#define assert(x)

#include "ptest_data.h"

struct gmp_div_inverse
{
  /* Normalization shift count. */
  unsigned shift;
  /* Normalized divisor (d0 unused for mpn_div_qr_1) */
  mp_limb_t d1, d0;
  /* Inverse, for 2/1 or 3/2. */
  mp_limb_t di;
};

static ptest_indata_t inbuf;
static volatile unsigned wait_flag;

#define Q_LEN 7
mp_limb_t primorial[Q_LEN] = { 3990313926u, 1301064672u, 203676164u, 3309914335u, 2220064684u, 2929319840u, 153406374u };

#define CHAR_BIT 8
#define GMP_LIMB_BITS (sizeof(mp_limb_t) * CHAR_BIT)

#define GMP_LIMB_MAX (~ (mp_limb_t) 0)
#define GMP_LIMB_HIGHBIT ((mp_limb_t) 1 << (GMP_LIMB_BITS - 1))

#define GMP_HLIMB_BIT ((mp_limb_t) 1 << (GMP_LIMB_BITS / 2))
#define GMP_LLIMB_MASK (GMP_HLIMB_BIT - 1)

#undef MPN_REDC_1
#define MPN_REDC_1(rp, up, mp, n, invm)                                 \
  do {                                                                  \
    mp_limb_t cy;                                                       \
    cy = mpn_redc_1 (rp, up, mp, n, invm);                              \
    if (cy != 0)                                                        \
      mpn_sub_n (rp, rp, mshifted, n);                                  \
  } while (0)

const unsigned char  binvert_limb_table[128] = {
  0x01, 0xAB, 0xCD, 0xB7, 0x39, 0xA3, 0xC5, 0xEF,
  0xF1, 0x1B, 0x3D, 0xA7, 0x29, 0x13, 0x35, 0xDF,
  0xE1, 0x8B, 0xAD, 0x97, 0x19, 0x83, 0xA5, 0xCF,
  0xD1, 0xFB, 0x1D, 0x87, 0x09, 0xF3, 0x15, 0xBF,
  0xC1, 0x6B, 0x8D, 0x77, 0xF9, 0x63, 0x85, 0xAF,
  0xB1, 0xDB, 0xFD, 0x67, 0xE9, 0xD3, 0xF5, 0x9F,
  0xA1, 0x4B, 0x6D, 0x57, 0xD9, 0x43, 0x65, 0x8F,
  0x91, 0xBB, 0xDD, 0x47, 0xC9, 0xB3, 0xD5, 0x7F,
  0x81, 0x2B, 0x4D, 0x37, 0xB9, 0x23, 0x45, 0x6F,
  0x71, 0x9B, 0xBD, 0x27, 0xA9, 0x93, 0xB5, 0x5F,
  0x61, 0x0B, 0x2D, 0x17, 0x99, 0x03, 0x25, 0x4F,
  0x51, 0x7B, 0x9D, 0x07, 0x89, 0x73, 0x95, 0x3F,
  0x41, 0xEB, 0x0D, 0xF7, 0x79, 0xE3, 0x05, 0x2F,
  0x31, 0x5B, 0x7D, 0xE7, 0x69, 0x53, 0x75, 0x1F,
  0x21, 0xCB, 0xED, 0xD7, 0x59, 0xC3, 0xE5, 0x0F,
  0x11, 0x3B, 0x5D, 0xC7, 0x49, 0x33, 0x55, 0xFF
};

#define binvert_limb(inv,n)                                             \
  do {                                                                  \
    mp_limb_t  __n = (n);                                               \
    mp_limb_t  __inv;                                                   \
    assert ((__n & 1) == 1);                                            \
                                                                        \
    __inv = binvert_limb_table[(__n/2) & 0x7F]; /*  8 */                \
    __inv = 2 * __inv - __inv * __inv * __n;                            \
    __inv = 2 * __inv - __inv * __inv * __n;                            \
                                                                        \
    assert ((__inv * __n) == 1);                                        \
    (inv) = __inv;                                                      \
  } while (0)

#define gmp_clz(count, x) do {                                          \
    unsigned __clz_x = (x);                                             \
    union {unsigned asInt; float asFloat;} _u;                          \
    __clz_x &= ~(__clz_x >> 1);                                         \
    _u.asFloat = (float)__clz_x + 0.5f;                                 \
    (count) = 158 - (_u.asInt >> 23);                                   \
  } while (0)

#define gmp_add_ssaaaa(sh, sl, ah, al, bh, bl) \
  do {                                                                  \
    mp_limb_t __x;                                                      \
    __x = (al) + (bl);                                                  \
    (sh) = (ah) + (bh) + (__x < (al));                                  \
    (sl) = __x;                                                         \
  } while (0)

#define gmp_sub_ddmmss(sh, sl, ah, al, bh, bl) \
  do {                                                                  \
    mp_limb_t __x;                                                      \
    __x = (al) - (bl);                                                  \
    (sh) = (ah) - (bh) - ((al) < (bl));                                 \
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

#define gmp_udiv_qr_3by2(q, r1, r0, n2, n1, n0, d1, d0, dinv)           \
  do {                                                                  \
    mp_limb_t _q0, _t1, _t0, _mask;                                     \
    gmp_umul_ppmm ((q), _q0, (n2), (dinv));                             \
    gmp_add_ssaaaa ((q), _q0, (q), _q0, (n2), (n1));                    \
                                                                        \
    /* Compute the two most significant limbs of n - q'd */             \
    (r1) = (n1) - (d1) * (q);                                           \
    gmp_sub_ddmmss ((r1), (r0), (r1), (n0), (d1), (d0));                \
    gmp_umul_ppmm (_t1, _t0, (d0), (q));                                \
    gmp_sub_ddmmss ((r1), (r0), (r1), (r0), _t1, _t0);                  \
    (q)++;                                                              \
                                                                        \
    /* Conditionally adjust q and the remainders */                     \
    _mask = - (mp_limb_t) ((r1) >= _q0);                                \
    (q) += _mask;                                                       \
    gmp_add_ssaaaa ((r1), (r0), (r1), (r0), _mask & (d1), _mask & (d0)); \
    if ((r1) >= (d1))                                                   \
      {                                                                 \
        if ((r1) > (d1) || (r0) >= (d0))                                \
          {                                                             \
            (q)++;                                                      \
            gmp_sub_ddmmss ((r1), (r0), (r1), (r0), (d1), (d0));        \
          }                                                             \
      }                                                                 \
  } while (0)

#define mpn_invert_limb(x) mpn_invert_3by2 ((x), 0)

#if 0
static mp_size_t
mpn_normalized_size (mp_srcptr xp, mp_size_t n)
{
  for (; n > 0 && xp[n-1] == 0; n--)
    ;
  return n;
}
#endif

static mp_limb_t
mpn_sub_1 (mp_ptr rp, mp_srcptr ap, mp_size_t n, mp_limb_t b)
{
  mp_size_t i;

  assert (n > 0);

  i = 0;
  do
    {
      mp_limb_t a = ap[i];
      /* Carry out */
      mp_limb_t cy = a < b;;
      rp[i] = a - b;
      b = cy;
    }
  while (++i < n);

  return b;
}

static mp_limb_t
mpn_sub_n (mp_ptr rp, mp_srcptr ap, mp_srcptr bp, mp_size_t n)
{
  mp_size_t i;
  mp_limb_t cy;

  for (i = 0, cy = 0; i < n; i++)
    {
      mp_limb_t a, b;
      a = ap[i]; b = bp[i];
      b += cy;
      cy = (b < cy);
      cy += (a < b);
      rp[i] = a - b;
    }
  return cy;
}

static mp_limb_t
mpn_add_1_inplace(mp_ptr np, mp_size_t n, mp_limb_t a)
{
  for (mp_size_t i = 0; i < n && a; ++i)
  {
    mp_limb_t b = np[i];
    mp_limb_t r = b + a;
    a = r < b;
    np[i] = r;
  }
  return a;
}

static mp_limb_t
mpn_add_n (mp_ptr rp, mp_srcptr ap, mp_srcptr bp, mp_size_t n)
{
  mp_size_t i;
  mp_limb_t cy;

  for (i = 0, cy = 0; i < n; i++)
    {
      mp_limb_t a, b, r;
      a = ap[i]; b = bp[i];
      r = a + cy;
      cy = (r < cy);
      r += b;
      cy += (r < b);
      rp[i] = r;
    }
  return cy;
}

static mp_limb_t
mpn_add_n2 (mp_ptr rp, mp_srcptr ap, mp_size_t an, mp_srcptr bp, mp_size_t bn)
{
  mp_size_t i;
  mp_limb_t cy;
  assert(an >= bn);

  for (i = 0, cy = 0; i < bn; i++)
    {
      mp_limb_t a, b, r;
      a = ap[i]; b = bp[i];
      r = a + cy;
      cy = (r < cy);
      r += b;
      cy += (r < b);
      rp[i] = r;
    }

  for (; i < an; i++)
    {
      mp_limb_t a, r;
      a = ap[i];
      r = a + cy;
      cy = (r < cy);
      rp[i] = r;
    }

  return cy;
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

  qh = ~u1 / uh;
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

static void
mpn_div_qr_invert (struct gmp_div_inverse *inv,
                   mp_srcptr dp, mp_size_t dn)
{
  assert (dn > 2);

    {
      unsigned shift;
      mp_limb_t d1, d0;

      d1 = dp[dn-1];
      d0 = dp[dn-2];
      assert (d1 > 0);
      gmp_clz (shift, d1);
      inv->shift = shift;
      if (shift > 0)
        {
          d1 = (d1 << shift) | (d0 >> (GMP_LIMB_BITS - shift));
          d0 = (d0 << shift) | (dp[dn-3] >> (GMP_LIMB_BITS - shift));
        }
      inv->d1 = d1;
      inv->d0 = d0;
      inv->di = mpn_invert_3by2 (d1, d0);
    }
}

static mp_limb_t
mpn_mul_1 (mp_ptr rp, mp_srcptr up, mp_size_t n, mp_limb_t vl)
{
  mp_limb_t ul, cl, hpl, lpl;

  assert (n >= 1);

  cl = 0;
  do
    {
      ul = *up++;
      gmp_umul_ppmm (hpl, lpl, ul, vl);

      lpl += cl;
      cl = (lpl < cl) + hpl;

      *rp++ = lpl;
    }
  while (--n != 0);

  return cl;
}

static mp_limb_t
mpn_addmul_1 (mp_ptr rp, mp_srcptr up, mp_size_t n, mp_limb_t vl)
{
  mp_limb_t ul, cl, hpl, lpl, rl;
  mp_limb_t ul2, cl2, hpl2, lpl2, rl2;
  mp_ptr up2, rp2;
  mp_size_t halfn = n>>1;

  assert (n >= 1);

  cl = 0;
  if (n & 1)
    {
      ul = *up++;
      gmp_umul_ppmm (hpl, lpl, ul, vl);

      lpl += cl;
      cl = (lpl < cl) + hpl;

      rl = *rp;
      lpl = rl + lpl;
      cl += lpl < rl;
      *rp++ = lpl;

      n--;
    }

  up2 = up + (n>>1);
  rp2 = rp + (n>>1);
  cl2 = 0;

  while (n != 0)
    {
      ul = *up++;
      gmp_umul_ppmm (hpl, lpl, ul, vl);
        ul2 = *up2++;
        gmp_umul_ppmm (hpl2, lpl2, ul2, vl);

      lpl += cl;
      cl = (lpl < cl) + hpl;
        lpl2 += cl2;
        cl2 = (lpl2 < cl2) + hpl2;

      rl = *rp;
      lpl = rl + lpl;
        rl2 = *rp2;
        lpl2 = rl2 + lpl2;
      cl += lpl < rl;
      *rp++ = lpl;
        cl2 += lpl2 < rl2;
        *rp2++ = lpl2;

      n -= 2;
    }

  cl2 += mpn_add_1_inplace(rp, halfn, cl);

  return cl2;
}

static mp_limb_t
mpn_submul_1 (mp_ptr rp, mp_srcptr up, mp_size_t n, mp_limb_t vl)
{
  mp_limb_t ul, cl, hpl, lpl, rl;

  assert (n >= 1);

  cl = 0;
  do
    {
      ul = *up++;
      gmp_umul_ppmm (hpl, lpl, ul, vl);

      lpl += cl;
      cl = (lpl < cl) + hpl;

      rl = *rp;
      lpl = rl - lpl;
      cl += lpl > rl;
      *rp++ = lpl;
    }
  while (--n != 0);

  return cl;
}

#if 0
static mp_limb_t
mpn_mul (mp_ptr rp, mp_srcptr up, mp_size_t un, mp_srcptr vp, mp_size_t vn)
{
  assert (un >= vn);
  assert (vn >= 1);

  /* We first multiply by the low order limb. This result can be
     stored, not added, to rp. We also avoid a loop for zeroing this
     way. */

  rp[un] = mpn_mul_1 (rp, up, un, vp[0]);
  rp += 1, vp += 1, vn -= 1;

  /* Now accumulate the product of up[] and the next higher limb from
     vp[]. */

  while (vn >= 1)
    {
      rp[un] = mpn_addmul_1 (rp, up, un, vp[0]);
      rp += 1, vp += 1, vn -= 1;
    }
  return rp[un - 1];
}
#endif

static mp_limb_t
mpn_redc_1 (mp_ptr rp, mp_ptr up, mp_srcptr mp, mp_size_t n, mp_limb_t invm)
{
  mp_size_t j;
  mp_limb_t cy;

  assert (n > 0);

  for (j = n - 1; j >= 0; j--)
    {
      cy = mpn_addmul_1 (up, mp, n, up[0] * invm);
      assert (up[0] == 0);
      up[0] = cy;
      up++;
    }

  cy = mpn_add_n (rp, up, up - n, n);
  return cy;
}

static void
mpn_div_r_pi1 (mp_ptr np, mp_size_t nn, mp_limb_t n1,
               mp_srcptr dp, mp_size_t dn,
               mp_limb_t dinv)
{
  mp_size_t i;

  mp_limb_t d1, d0;
  mp_limb_t cy, cy1;
  mp_limb_t q;

  assert (dn > 2);
  assert (nn >= dn);

  d1 = dp[dn - 1];
  d0 = dp[dn - 2];

  assert ((d1 & GMP_LIMB_HIGHBIT) != 0);
  /* Iteration variable is the index of the q limb.
   *
   * We divide <n1, np[dn-1+i], np[dn-2+i], np[dn-3+i],..., np[i]>
   * by            <d1,          d0,        dp[dn-3],  ..., dp[0] >
   */

  i = nn - dn;
  do
    {
      mp_limb_t n0 = np[dn-1+i];

      if (n1 == d1 && n0 == d0)
        {
          q = GMP_LIMB_MAX;
          mpn_submul_1 (np+i, dp, dn, q);
          n1 = np[dn-1+i];      /* update n1, last loop's value will now be invalid */
        }
      else
        {
          gmp_udiv_qr_3by2 (q, n1, n0, n1, n0, np[dn-2+i], d1, d0, dinv);

          cy = mpn_submul_1 (np + i, dp, dn-2, q);

          cy1 = n0 < cy;
          n0 = n0 - cy;
          cy = n1 < cy1;
          n1 = n1 - cy1;
          np[dn-2+i] = n0;

          if (cy != 0)
            {
              n1 += d1 + mpn_add_n (np + i, np + i, dp, dn - 1);
              q--;
            }
        }
    }
  while (--i >= 0);

  np[dn - 1] = n1;
}


static void
mpn_div_r_preinv_ns (mp_ptr np, mp_size_t nn,
                     mp_srcptr dp, mp_size_t dn,
                     const struct gmp_div_inverse *inv)
{
  assert (dn > 2);
  assert (nn >= dn);

    {
      mp_limb_t nh;

      assert (inv->d1 == dp[dn-1]);
      assert (inv->d0 == dp[dn-2]);
      assert ((inv->d1 & GMP_LIMB_HIGHBIT) != 0);

      nh = 0;

      mpn_div_r_pi1 (np, nn, nh, dp, dn, inv->di);
    }
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
mpn_rshift (mp_ptr rp, mp_srcptr up, mp_size_t n, unsigned int cnt)
{
  mp_limb_t high_limb, low_limb;
  unsigned int tnc;
  mp_size_t i;
  mp_limb_t retval;

  assert (n >= 1);
  assert (cnt >= 1);
  assert (cnt < GMP_LIMB_BITS);

  tnc = GMP_LIMB_BITS - cnt;
  high_limb = *up++;
  retval = (high_limb << tnc);
  low_limb = high_limb >> cnt;

  for (i = n; --i != 0;)
    {
      high_limb = *up++;
      *rp++ = low_limb | (high_limb << tnc);
      low_limb = high_limb >> cnt;
    }
  *rp = low_limb;

  return retval;
}

static void
mpn_sqr (mp_ptr rp, mp_srcptr up, mp_size_t n)
{
  mp_size_t i;
  {
    mp_limb_t ul, lpl;
    ul = up[0];
    gmp_umul_ppmm (rp[1], lpl, ul, ul);
    rp[0] = lpl;
  }
  if (n > 1)
    {
      mp_double_fixed_len_num tp;
      mp_limb_t cy;

      cy = mpn_mul_1 (tp, up + 1, n - 1, up[0]);
      tp[n - 1] = cy;
      for (i = 2; i < n; i++)
        {
          mp_limb_t cy;
          cy = mpn_addmul_1 (&tp[2 * i - 2], up + i, n - i, up[i - 1]);
          tp[n + i - 2] = cy;
        }

      for (i = 0; i < n; i++)
      {
        mp_limb_t ul, lpl;
        ul = up[i];
        gmp_umul_ppmm (rp[2 * i + 1], lpl, ul, ul);
        rp[2 * i] = lpl;
      }
      cy = mpn_lshift (tp, tp, 2 * n - 2, 1);
      cy += mpn_add_n (rp + 1, rp + 1, tp, 2 * n - 2);
      rp[2 * n - 1] += cy;
    }
}

#if 0
static int
mpn_cmp (mp_srcptr ap, mp_srcptr bp, mp_size_t n)
{
  while (--n >= 0)
    {
      if (ap[n] != bp[n])
        return ap[n] > bp[n] ? 1 : -1;
    }
  return 0;
}
#endif

static int
my_fermat_test (const mp_srcptr msp, mp_size_t mn)
{
  mp_fixed_len_num r, ep, mshifted;
  mp_double_fixed_len_num t;
  mp_ptr mp, tp, rp;
  mp_size_t en;
  mp_limb_t startbit;
  struct gmp_div_inverse minv;
  unsigned i;

  // REDCify: r = B^n * 2 % M
  mp = msp;
  mpn_div_qr_invert (&minv, mp, mn);

  if (minv.shift > 0)
  {
    mpn_lshift(mshifted, mp, mn, minv.shift);
    mp = mshifted;
  }
  else
  {
    for (i = 0; i < mn; ++i) mshifted[i] = mp[i];
  }

  for (i = 0; i < mn; ++i) r[i] = 0;
  r[mn] = 1 << minv.shift;
  mpn_div_r_preinv_ns(r, mn+1, mp, mn, &minv);

  if (minv.shift > 0)
  {
    mpn_rshift(r, r, mn, minv.shift);
    mp = msp;
  }

  en = mn;
  mpn_sub_1(ep, mp, mn, 1);
  gmp_clz(startbit, ep[en-1]);
  startbit = GMP_LIMB_HIGHBIT >> startbit;
  mp_limb_t mi;
  binvert_limb(mi, mp[0]);
  mi = -mi;

  tp = t;
  rp = r;

  while (en-- > 0)
    {
      mp_limb_t w = ep[en];
      mp_limb_t bit;

      bit = startbit;
      startbit = GMP_LIMB_HIGHBIT;
      do
        {
          mpn_sqr (tp, rp, mn);
          MPN_REDC_1(rp, tp, mp, mn, mi);
          if (w & bit)
          {
            mp_limb_t carry = mpn_lshift (rp, rp, mn, 1);
            while (carry)
            {
              carry -= mpn_sub_n(rp, rp, mshifted, mn);
            }
          }
          bit >>= 1;
        }
      while (bit > 0);
    }

  // DeREDCify - necessary as rp can have a large
  //             multiple of m in it (although I'm not 100% sure
  //             why it can't after this redc!)
  for (i = 0; i < mn; ++i) t[i] = r[i];
  for (; i < mn*2; ++i) t[i] = 0;
  MPN_REDC_1(rp, tp, mp, mn, mi);

  if (rp[mn-1] != 0)
  {
    // Compare to m+1
    mpn_sub_1(tp, rp, mn, 1);
    for (i = 0; i < mn; ++i) if (tp[i] != mp[i]) return 0;
    return 1;
  }

  // COmpare to 1
  for (i = 1; i < mn-1; ++i) if (rp[i] != 0) return 0;
  return rp[0] == 1;
}

void null_isr(int);

int main()
{
  e_coreid_t coreid;
  coreid = e_get_coreid();
  unsigned int row, col, core;
  e_coords_from_coreid(coreid, &row, &col);
  core = row * 4 + col;

  ptest_indata_t* in = (ptest_indata_t*)(0x8f000000+(0x00020000*core));
  ptest_outdata_t* out = (ptest_outdata_t*)((char*)in + sizeof(ptest_indata_t));

  // Must set up a null interrupt routine to allow host to wake us from idle
  e_irq_attach(E_SYNC, null_isr);
  e_irq_mask(E_SYNC, E_FALSE);
  e_irq_global_mask(E_FALSE);

  while (1)
  {
    e_dma_copy(&inbuf, in, sizeof(ptest_indata_t));

    // Ensure we don't run this block again.
    in->num_candidates = 0;

    mp_fixed_len_num c;
    mp_size_t cn;
    unsigned num_results = 0;
    for (unsigned i = 0; i < inbuf.num_candidates; ++i)
    {
      unsigned primes = 0;

      c[Q_LEN] = mpn_mul_1(c, primorial, Q_LEN, inbuf.k[i]);
      if (c[Q_LEN]) cn = Q_LEN+1;
      else          cn = Q_LEN;
      mpn_add_n2(c, inbuf.n, inbuf.nn, c, cn);
      cn = inbuf.nn;

      if (my_fermat_test(c, cn)) primes++;
      if (primes < 1) continue;

      mpn_add_1_inplace(c, cn, 4);
      if (my_fermat_test(c, cn)) primes++;

      mpn_add_1_inplace(c, cn, 2);
      if (my_fermat_test(c, cn)) primes++;

      mpn_add_1_inplace(c, cn, 4);
      if (my_fermat_test(c, cn)) primes++;
      if (primes < 2) continue;

      mpn_add_1_inplace(c, cn, 2);
      if (my_fermat_test(c, cn)) primes++;

      if (primes >= 3)
      {
        mpn_add_1_inplace(c, cn, 4);
        if (my_fermat_test(c, cn)) primes++;
      }

      ptest_result_t result;
      result.k = inbuf.k[i];
      result.primes = primes;
      out->result[num_results++] = result;
    }
    out->num_results = num_results;

    while (in->num_candidates == 0) __asm__ __volatile__ ("idle");
  }

  return 0;
}

void __attribute__((interrupt)) null_isr(int x)
{
  wait_flag = 0;
  return;
}
