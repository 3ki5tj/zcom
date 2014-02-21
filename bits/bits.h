#include "util.h" /* to define uint32_t etc. */
#ifndef BITS_H__
#define BITS_H__



/* these masks can be generated from the helper python script mkmask.py */
uint32_t bitsmask32_[33] = {0,
 0x00000001, 0x00000003, 0x00000007, 0x0000000f,
 0x0000001f, 0x0000003f, 0x0000007f, 0x000000ff,
 0x000001ff, 0x000003ff, 0x000007ff, 0x00000fff,
 0x00001fff, 0x00003fff, 0x00007fff, 0x0000ffff,
 0x0001ffff, 0x0003ffff, 0x0007ffff, 0x000fffff,
 0x001fffff, 0x003fffff, 0x007fffff, 0x00ffffff,
 0x01ffffff, 0x03ffffff, 0x07ffffff, 0x0fffffff,
 0x1fffffff, 0x3fffffff, 0x7fffffff, 0xffffffff};



uint32_t invbitsmask32_[33] = {0xffffffff,
 0xfffffffe, 0xfffffffc, 0xfffffff8, 0xfffffff0,
 0xffffffe0, 0xffffffc0, 0xffffff80, 0xffffff00,
 0xfffffe00, 0xfffffc00, 0xfffff800, 0xfffff000,
 0xffffe000, 0xffffc000, 0xffff8000, 0xffff0000,
 0xfffe0000, 0xfffc0000, 0xfff80000, 0xfff00000,
 0xffe00000, 0xffc00000, 0xff800000, 0xff000000,
 0xfe000000, 0xfc000000, 0xf8000000, 0xf0000000,
 0xe0000000, 0xc0000000, 0x80000000,        0x0};



/* make the ith bit */
#define MKBIT32(n) (((uint32_t) 1u) << (uint32_t) (n))

/* make a mask with the lowest n bits being 1s, other bits being 0s */
#define MKBITSMASK32(n) bitsmask32_[n]

/* make a mask with the lowest n bits being 0s, other bits being 1s */
#define MKINVBITSMASK32(n) invbitsmask32_[n]



/* if the lowest two bits of x are zeros, and x has n nonzero bits,
 * then x, x+1, x+2, x+3 have respectivly n, n+1, n+1, n+2 nonzero bits */
#define BITB2_(n)        n,         n+1,         n+1,         n+2
#define BITB4_(n) BITB2_(n), BITB2_(n+1), BITB2_(n+1), BITB2_(n+2)
#define BITB6_(n) BITB4_(n), BITB4_(n+1), BITB4_(n+1), BITB4_(n+2)

/* how many bits in 0..255 */
const unsigned char bits256_[256] = {
  BITB6_(0), BITB6_(1), BITB6_(1), BITB6_(2) };



/* count the number of 1 bits in x
 * http://graphics.stanford.edu/~seander/bithacks.html#CountBitsSetTable */
#define BITCOUNT8(x) bits256_[ ((unsigned char *)&(x))[0] ]

#define BITCOUNT16(x) ( BITCOUNT8(x) \
  + bits256_[ ((unsigned char *)&(x))[1] ] )

#define BITCOUNT24(x) ( BITCOUNT16(x) \
  + bits256_[ ((unsigned char *)&(x))[2] ] )

#define BITCOUNT32(x) ( BITCOUNT24(x) \
  + bits256_[ ((unsigned char *)&(x))[3] ] )

#define BITCOUNT40(x) ( BITCOUNT32(x) \
  + bits256_[ ((unsigned char *)&(x))[4] ] )

#define BITCOUNT48(x) ( BITCOUNT40(x) \
  + bits256_[ ((unsigned char *)&(x))[5] ] )

#define BITCOUNT56(x) ( BITCOUNT48(x) \
  + bits256_[ ((unsigned char *)&(x))[6] ] )

#define BITCOUNT64(x) ( BITCOUNT56(x) \
  + bits256_[ ((unsigned char *)&(x))[7] ] )

INLINE int bitcount32(uint32_t x) { return BITCOUNT32(x); }


/* invert the lowest n bits */
#define BITNINV32(x, n) ((x) ^ MKBITSMASK32(n))

INLINE uint32_t bitninv32(uint32_t x, int n)
{ return BITNINV32(x, n); }



#define BITR2_(n)        n,         n + 2*64,         n + 1*64,         n + 3*64
#define BITR4_(n) BITR2_(n), BITR2_(n + 2*16), BITR2_(n + 1*16), BITR2_(n + 3*16)
#define BITR6_(n) BITR4_(n), BITR4_(n + 2*4 ), BITR4_(n + 1*4 ), BITR4_(n + 3*4 )

const unsigned char bitrev256_[256] = {
  BITR6_(0), BITR6_(2), BITR6_(1), BITR6_(3) };

/* reverse bits a = rev(b) */
#define BITREV32(a, b) { \
  ((unsigned char *)&(a))[0] = bitrev256_[ ((unsigned char *)&(b))[3] ]; \
  ((unsigned char *)&(a))[1] = bitrev256_[ ((unsigned char *)&(b))[2] ]; \
  ((unsigned char *)&(a))[2] = bitrev256_[ ((unsigned char *)&(b))[1] ]; \
  ((unsigned char *)&(a))[3] = bitrev256_[ ((unsigned char *)&(b))[0] ]; }

INLINE uint32_t bitrev32(uint32_t x)
{ uint32_t c; BITREV32(c, x) return c; }




const int bruijn32_[32] =
  { 0, 1, 28, 2, 29, 14, 24, 3, 30, 22, 20, 15, 25, 17, 4, 8,
    31, 27, 13, 23, 21, 19, 16, 7, 26, 12, 18, 6, 11, 5, 10, 9};

#define BIT2ID32(b) bruijn32_[((b) * 0x077CB531) >> 27]



/* macro version of bitfirstlow(); */
#define BITFIRSTLOW32(id, x, b) { \
  (b) = (x) & (-(x)); \
  (id) = BIT2ID32(b); }

/* find the index of the lowest 1 bit
 * on return *b is the nonzero bit
 * http://graphics.stanford.edu/~seander/bithacks.html#ZerosOnRightMultLookup
 * ``Using de Bruijn Sequences to Index 1 in a Computer Word''
 * by Charles E. Leiserson, Harald Prokof, and Keith H. Randall.
 * http://supertech.csail.mit.edu/papers/debruijn.pdf */
INLINE int bitfirstlow32(uint32_t x, uint32_t *b)
{
  (*b) = x & (-x); /* such that only the lowest 1-bit survives */
  return BIT2ID32(*b);
}



/* find the least significant nonzero bit
 * http://en.wikipedia.org/wiki/Find_first_set#Tool_and_library_support
 * http://chessprogramming.wikispaces.com/BitScan */
#ifdef __INTEL_COMPILER

/* directly map to the intrinsic function
 * use it only when you are sure x != 0
 * the intrinsic for the most significant nonzero bit
 * is _bit_scan_reverse(x) */
#define BITFIRSTNZ32(x)  _bit_scan_forward(x)
#define BITFIRST32(id, x) id = (x) ? BITFIRSTNZ32(x) : 0;
#define bitfirst32(x) ((x) ? BITFIRSTNZ32(x) : 0)

#elif defined(__GNUC__) && !defined(__INTEL_COMPILER)

/* the intrinsic for the least significant nonzero bit
 * is __builtin_clz */
#define BITFIRSTNZ32(x)  __builtin_ctz(x)
#define BITFIRST32(id, x) id = (x) ? BITFIRSTNZ32(x) : 0;
#define bitfirst32(x)  ((x) ? BITFIRSTNZ32(x) : 0)

#else /* general */

#define BITFIRSTNZ32(x)   bitfirst32(x)
#define BITFIRST32(id, x) { uint32_t b32_; BITFIRSTLOW32(id, x, b32_) }

/* index of nonzero bit */
INLINE int bitfirst32(uint32_t x)
{
  uint32_t b;
  return bitfirstlow32(x, &b);
}

#endif /* defined(__INTEL_COMPILER/__GNUC__) */



uint64_t bitsmask64_[65] = {0,
 CU64(0x0000000000000001), CU64(0x0000000000000003), CU64(0x0000000000000007), CU64(0x000000000000000f),
 CU64(0x000000000000001f), CU64(0x000000000000003f), CU64(0x000000000000007f), CU64(0x00000000000000ff),
 CU64(0x00000000000001ff), CU64(0x00000000000003ff), CU64(0x00000000000007ff), CU64(0x0000000000000fff),
 CU64(0x0000000000001fff), CU64(0x0000000000003fff), CU64(0x0000000000007fff), CU64(0x000000000000ffff),
 CU64(0x000000000001ffff), CU64(0x000000000003ffff), CU64(0x000000000007ffff), CU64(0x00000000000fffff),
 CU64(0x00000000001fffff), CU64(0x00000000003fffff), CU64(0x00000000007fffff), CU64(0x0000000000ffffff),
 CU64(0x0000000001ffffff), CU64(0x0000000003ffffff), CU64(0x0000000007ffffff), CU64(0x000000000fffffff),
 CU64(0x000000001fffffff), CU64(0x000000003fffffff), CU64(0x000000007fffffff), CU64(0x00000000ffffffff),
 CU64(0x00000001ffffffff), CU64(0x00000003ffffffff), CU64(0x00000007ffffffff), CU64(0x0000000fffffffff),
 CU64(0x0000001fffffffff), CU64(0x0000003fffffffff), CU64(0x0000007fffffffff), CU64(0x000000ffffffffff),
 CU64(0x000001ffffffffff), CU64(0x000003ffffffffff), CU64(0x000007ffffffffff), CU64(0x00000fffffffffff),
 CU64(0x00001fffffffffff), CU64(0x00003fffffffffff), CU64(0x00007fffffffffff), CU64(0x0000ffffffffffff),
 CU64(0x0001ffffffffffff), CU64(0x0003ffffffffffff), CU64(0x0007ffffffffffff), CU64(0x000fffffffffffff),
 CU64(0x001fffffffffffff), CU64(0x003fffffffffffff), CU64(0x007fffffffffffff), CU64(0x00ffffffffffffff),
 CU64(0x01ffffffffffffff), CU64(0x03ffffffffffffff), CU64(0x07ffffffffffffff), CU64(0x0fffffffffffffff),
 CU64(0x1fffffffffffffff), CU64(0x3fffffffffffffff), CU64(0x7fffffffffffffff), CU64(0xffffffffffffffff)};

uint64_t invbitsmask64_[65] = {CU64(0xffffffffffffffff),
 CU64(0xfffffffffffffffe), CU64(0xfffffffffffffffc), CU64(0xfffffffffffffff8), CU64(0xfffffffffffffff0),
 CU64(0xffffffffffffffe0), CU64(0xffffffffffffffc0), CU64(0xffffffffffffff80), CU64(0xffffffffffffff00),
 CU64(0xfffffffffffffe00), CU64(0xfffffffffffffc00), CU64(0xfffffffffffff800), CU64(0xfffffffffffff000),
 CU64(0xffffffffffffe000), CU64(0xffffffffffffc000), CU64(0xffffffffffff8000), CU64(0xffffffffffff0000),
 CU64(0xfffffffffffe0000), CU64(0xfffffffffffc0000), CU64(0xfffffffffff80000), CU64(0xfffffffffff00000),
 CU64(0xffffffffffe00000), CU64(0xffffffffffc00000), CU64(0xffffffffff800000), CU64(0xffffffffff000000),
 CU64(0xfffffffffe000000), CU64(0xfffffffffc000000), CU64(0xfffffffff8000000), CU64(0xfffffffff0000000),
 CU64(0xffffffffe0000000), CU64(0xffffffffc0000000), CU64(0xffffffff80000000), CU64(0xffffffff00000000),
 CU64(0xfffffffe00000000), CU64(0xfffffffc00000000), CU64(0xfffffff800000000), CU64(0xfffffff000000000),
 CU64(0xffffffe000000000), CU64(0xffffffc000000000), CU64(0xffffff8000000000), CU64(0xffffff0000000000),
 CU64(0xfffffe0000000000), CU64(0xfffffc0000000000), CU64(0xfffff80000000000), CU64(0xfffff00000000000),
 CU64(0xffffe00000000000), CU64(0xffffc00000000000), CU64(0xffff800000000000), CU64(0xffff000000000000),
 CU64(0xfffe000000000000), CU64(0xfffc000000000000), CU64(0xfff8000000000000), CU64(0xfff0000000000000),
 CU64(0xffe0000000000000), CU64(0xffc0000000000000), CU64(0xff80000000000000), CU64(0xff00000000000000),
 CU64(0xfe00000000000000), CU64(0xfc00000000000000), CU64(0xf800000000000000), CU64(0xf000000000000000),
 CU64(0xe000000000000000), CU64(0xc000000000000000), CU64(0x8000000000000000), CU64(0x0000000000000000)};



/* make the ith bit */
#define MKBIT64(n) (((uint64_t) 1u) << (uint64_t) (n))

/* make a mask with the lowest n bits being 1s, other bits being 0s */
#define MKBITSMASK64(n) bitsmask64_[n]

/* make a mask with the lowest n bits being 0s, other bits being 1s */
#define MKINVBITSMASK64(n) invbitsmask64_[n]



INLINE int bitcount64(uint64_t x) { return BITCOUNT64(x); }



/* invert the lowest k bits */
#define BITNINV64(x, n) ((x) ^ MKBITSMASK64(n))

INLINE uint64_t bitninv64(uint64_t x, int n)
{ return BITNINV64(x, n); }



/* reverse bits a = rev(b) */
#define BITREV64(a, b) { \
  ((unsigned char *)&(a))[0] = bitrev256_[ ((unsigned char *)&(b))[7] ]; \
  ((unsigned char *)&(a))[1] = bitrev256_[ ((unsigned char *)&(b))[6] ]; \
  ((unsigned char *)&(a))[2] = bitrev256_[ ((unsigned char *)&(b))[5] ]; \
  ((unsigned char *)&(a))[3] = bitrev256_[ ((unsigned char *)&(b))[4] ]; \
  ((unsigned char *)&(a))[4] = bitrev256_[ ((unsigned char *)&(b))[3] ]; \
  ((unsigned char *)&(a))[5] = bitrev256_[ ((unsigned char *)&(b))[2] ]; \
  ((unsigned char *)&(a))[6] = bitrev256_[ ((unsigned char *)&(b))[1] ]; \
  ((unsigned char *)&(a))[7] = bitrev256_[ ((unsigned char *)&(b))[0] ]; }

INLINE uint64_t bitrev64(uint64_t x)
{ uint64_t c; BITREV64(c, x) return c; }



const int bruijn64_[64] =
  {  0,  1,  2,  7,  3, 13,  8, 19,  4, 25, 14, 28,  9, 34, 20, 40,
     5, 17, 26, 38, 15, 46, 29, 48, 10, 31, 35, 54, 21, 50, 41, 57,
    63,  6, 12, 18, 24, 27, 33, 39, 16, 37, 45, 47, 30, 53, 49, 56,
    62, 11, 23, 32, 36, 44, 52, 55, 61, 22, 43, 51, 60, 42, 59, 58};

#define BIT2ID64(b) bruijn64_[((b) * CU64(0x218A392CD3D5DBF)) >> 58]




/* macro version of bitfirstlow(); */
#define BITFIRSTLOW64(id, x, b) { \
  (b) = (x) & (-(x)); \
  (id) = BIT2ID64(b); }

/* find the index of the lowest 1 bit
 * on return *b is the nonzero bit
 * http://graphics.stanford.edu/~seander/bithacks.html#ZerosOnRightMultLookup
 * ``Using de Bruijn Sequences to Index 1 in a Computer Word''
 * by Charles E. Leiserson, Harald Prokof, and Keith H. Randall.
 * http://supertech.csail.mit.edu/papers/debruijn.pdf */
INLINE int bitfirstlow64(uint64_t x, uint64_t *b)
{
  (*b) = x & (-x); /* such that only the lowest 1-bit survives */
  return BIT2ID64(*b);
}



#define BITFIRSTNZ64(x) bitfirst64(x)

#define BITFIRST64(id, x) { uint64_t b64_; BITFIRSTLOW64(id, x, b64_) }

/* index of nonzero bit */
INLINE int bitfirst64(uint64_t x)
{
  uint64_t b;
  return bitfirstlow64(x, &b);
}


#ifndef WORDBITS_
#define WORDBITS_ 32
#endif /* defined(WORDBITS_) */


/* one word operations */
#if WORDBITS_ == 64

  #ifndef BITS_WORD_T
  #define BITS_WORD_T uint64_t
  #endif
  #define MKBIT(n)                MKBIT64(n)
  #define MKBITSMASK(n)           MKBITSMASK64(n)
  #define MKINVBITSMASK(n)        MKINVBITSMASK64(n)
  #define BITCOUNT(x)             BITCOUNT64(x)
  #define bitcount(x)             bitcount64(x)
  #define BITNINV(x, k)           BITNINV64(x, k)
  #define bitninv(x, k)           bitninv64(x, k)
  #define BITREV(x)               BITREV64(x)
  #define bitrev(x)               bitrev64(x)
  #define BIT2ID(b)               BIT2ID64(b)
  #define BITFIRSTLOW(id, x, b)   BITFIRSTLOW64(id, x, b)
  #define bitfirstlow(x, b)       bitfirstlow64(x, b)
  #define BITFIRSTNZ(x)           BITFIRSTNZ64(x)
  #define BITFIRST(id, x)         BITFIRST64(id, x)
  #define bitfirst(x)             bitfirst64(x)

#elif WORDBITS_ == 32

  #ifndef BITS_WORD_T
  #define BITS_WORD_T uint32_t
  #endif
  #define MKBIT(n)                MKBIT32(n)
  #define MKBITSMASK(n)           MKBITSMASK32(n)
  #define MKINVBITSMASK(n)        MKINVBITSMASK32(n)
  #define BITCOUNT(x)             BITCOUNT32(x)
  #define bitcount(x)             bitcount32(x)
  #define BITNINV(x, k)           BITNINV32(x, k)
  #define bitninv(x, k)           bitninv32(x, k)
  #define BITREV(x)               BITREV32(x)
  #define bitrev(x)               bitrev32(x)
  #define BIT2ID(b)               BIT2ID32(b)
  #define BITFIRSTLOW(id, x, b)   BITFIRSTLOW32(id, x, b)
  #define bitfirstlow(x, b)       bitfirstlow32(x, b)
  #define BITFIRSTNZ(x)           BITFIRSTNZ32(x)
  #define BITFIRST(id, x)         BITFIRST32(id, x)
  #define bitfirst(x)             bitfirst32(x)

#else

  #error bad WORDBITS_ definition

#endif /* WORDBITS_ == 32 */



/* multiple-word bit-set operations */

#if WORDBITS_ == 32
  #define BITS_GETNW(n)      (((n) + 31) >> 5)
  #define BITS_GETNQ(n)      ((n) >> 5)            /* n / 32 */
  #define BITS_GETNR(n)      ((n) & 0x1f)          /* n % 32 */
#elif WORDBITS_ == 64
  #define BITS_GETNW(n)      (((n) + 63) >> 6)
  #define BITS_GETNQ(n)      ((n) >> 6)            /* n / 64 */
  #define BITS_GETNR(n)      ((n) & 0x3f)          /* n % 64 */
#else
  /* get the number of words to represent n elements */
  #define BITS_GETNW(n)      (((n) + WORDBITS_ - 1) / WORDBITS_)
  #define BITS_GETNQ(n)      ((n) / WORDBITS_)
  #define BITS_GETNR(n)      ((n) % WORDBITS_)
#endif
#define BITS_GETNB(n)      MKBIT( BITS_GETNR(n) )


#define BITS_CPY(a, b, nw) { int iw_; \
  for (iw_ = 0; iw_ < nw; iw_++) \
    a[iw_] = b[iw_]; }

#define BITS_CLEAR(bs, nw) { int iw_; \
  for (iw_ = 0; iw_ < nw; iw_++) (bs)[iw_] = 0; }

#define BITS_MKBIT(i, bi, iq) { iq = BITS_GETNQ(i), bi = MKBIT(DG_IR(i)); }
/* construct the set `bs' with a single vertex */
#define BITS_ONEBIT(bs, nw, bi, iq) { BITS_CLEAR(bs, nw); bs[iq] = bi; }

/* make a mask for the lowest k elements, clear higher bits up to nw */
#define BITS_MKBITSMASK(bs, k, nw) { int iw_; \
  int kq_ = BITS_GETNQ(k), kr_ = BITS_GETNR(k); \
  for (iw_ = 0; iw_ < kq_; iw_++) \
    (bs)[iw_] = MKINVBITSMASK(0); \
  if (kr_ != 0) (bs)[kq_++] = MKBITSMASK(kr_); \
  for (iw_ = kq_; iw_ < nw; iw_++) \
    (bs)[iw_] = 0; }

INLINE BITS_WORD_T *bits_mkbitsmask(BITS_WORD_T *bs, int k, int nw)
{ BITS_MKBITSMASK(bs, k, nw) return bs; }

/* a similar set of macros BITS_MKINVSET() are defined later
 * after we have defined BITS_REMOVE() */



#define BITS_NONZERO(nz, bs, nw) { int iw_; \
  for (nz = 0, iw_ = 0; iw_ < nw; iw_++) \
    if ( (bs)[iw_] ) nz = 1; }

INLINE int bits_nonzero(BITS_WORD_T *bs, int nw)
{  int nz; BITS_NONZERO(nz, bs, nw) return nz; }

#define bits_iszero(bs, nw) ( !bits_nonzero(bs, nw) )

#define BITS_NOTEQUAL(ne, a, b, nw) { int iw_; \
  for (ne = 0, iw_ = 0; iw_ < nw; iw_++) \
    if ( (a)[iw_] != (b)[iw_] ) ne = 1; }

INLINE int bits_notequal(BITS_WORD_T *a, BITS_WORD_T *b, int nw)
{ int ne; BITS_NOTEQUAL(ne, a, b, nw) return ne; }

#define bits_equal(a, b, nw) ( !bits_notequal(a, b, nw) )

#define BITS_COUNT(cnt, bs, nw) { int iw_; \
  for (cnt = 0, iw_ = 0; iw_ < nw; iw_++) \
    cnt += bitcount( (bs)[iw_] ); }

INLINE int bits_count(BITS_WORD_T *bs, int nw) \
{ int cnt; BITS_COUNT(cnt, bs, nw) return cnt; }

/* invert the lowest n bits of b, keep the rest */
#define BITS_NINV(a, b, n, nw) { int iw_; \
  for (iw_ = 0; iw_ < nw - 1; iw_++) (a)[iw_] = ~(b)[iw_]; \
  (a)[iw_] = BITNINV((b)[iw_], n - iw_ * WORDBITS_); }



/* save the bit of `bs' corresponding to the first element to `b'
 * `bs' is a set, `b' is a word, `bq' is the word offset */
#define BITS_FIRSTBIT(bs, nw, b, bq) { \
  BITS_WORD_T w_; \
  (b) = 0; \
  for (bq = 0; bq < nw; bq++) { \
    if ((w_ = (bs)[bq]) != 0) { \
      (b) = w_ & (-w_); \
      break; \
    } } }

/* save the index of the first element to `id'
 * and the bit corresponding to `id' is `b'
 * `bs' is a set, `b' is a word, `bq' is the word offset */
#define BITS_FIRSTLOW(id, bs, nw, b, bq) { \
  BITS_WORD_T w_; \
  b = 0; id = 0; \
  for (bq = 0; bq < nw; bq++) { \
    if ((w_ = (bs)[bq]) != 0) { \
      BITFIRSTLOW(id, w_, b); \
      id += bq * WORDBITS_; \
      break; \
    } } }

/* return the index of the first element in the set */
INLINE int bits_first(BITS_WORD_T *bs, int nw)
{
  BITS_WORD_T b;
  int id, bq;
  BITS_FIRSTLOW(id, bs, nw, b, bq);
  return id;
}

/* assuming `bs' contains a single bit, return the index */
INLINE int bits_bit2id(BITS_WORD_T *bs, int nw)
{
  BITS_WORD_T b;
  int bq;

  for (bq = 0; bq < nw; bq++)
    if ((b = bs[bq]) != 0)
      return bq * WORDBITS_ * BIT2ID(b);
  return 0;
}



/* bitwise/set operations */
/* a ^= b */
#define BITS_XOR(a, b, nw) { int iw_; \
  for (iw_ = 0; iw_ < nw; iw_++) (a)[iw_] ^= (b)[iw_]; }

#define BITS_XOR1(a, b, iq)     (a)[iq] ^= b;

/* a = b ^ c */
#define BITS_XOR2(a, b, c, nw) { int iw_; \
  for (iw_ = 0; iw_ < nw; iw_++) \
    (a)[iw_] = (b)[iw_] ^ (c)[iw_]; }

INLINE BITS_WORD_T *bits_xor2(BITS_WORD_T *a, BITS_WORD_T *b, BITS_WORD_T *c, int nw)
{ BITS_XOR2(a, b, c, nw); return a; }

/* a |= b */
#define BITS_OR(a, b, nw) { int iw_; \
  for (iw_ = 0; iw_ < nw; iw_++) (a)[iw_] |= (b)[iw_]; }

#define BITS_OR1(a, b, iq)      (a)[iq] |= b;

/* a = b | c */
#define BITS_OR2(a, b, c, nw) { int iw_; \
  for (iw_ = 0; iw_ < nw; iw_++) \
    (a)[iw_] = (b)[iw_] | (c)[iw_]; }

INLINE BITS_WORD_T *bits_or2(BITS_WORD_T *a, BITS_WORD_T *b, BITS_WORD_T *c, int nw)
{ BITS_OR2(a, b, c, nw); return a; }

/* a &= b */
#define BITS_AND(a, b, nw) { int iw_; \
  for (iw_ = 0; iw_ < nw; iw_++) \
    (a)[iw_] &= (b)[iw_]; }

#define BITS_AND1(a, b, iq)     (a)[iq] &= b;

/* a = b & c */
#define BITS_AND2(a, b, c, nw) { int iw_; \
  for (iw_ = 0; iw_ < nw; iw_++) \
    (a)[iw_] = (b)[iw_] & (c)[iw_]; }

INLINE BITS_WORD_T *bits_and2(BITS_WORD_T *a, BITS_WORD_T *b, BITS_WORD_T *c, int nw)
{ BITS_AND2(a, b, c, nw); return a; }

/* a &= ~b or a -= b */
#define BITS_MINUS(a, b, nw) { int iw_; \
  for (iw_ = 0; iw_ < nw; iw_++) \
    (a)[iw_] &= ~(b)[iw_]; }

#define BITS_MINUS1(a, b, iq)  (a)[iq] &= ~(b);

/* a = b & ~c or a = b - c */
#define BITS_MINUS2(a, b, c, nw) { int iw_; \
  for (iw_ = 0; iw_ < nw; iw_++) \
    (a)[iw_] = (b)[iw_] & ~(c)[iw_]; }

/* a = b & ~c or a = b - c */
INLINE BITS_WORD_T *bits_minus2(BITS_WORD_T *a, BITS_WORD_T *b, BITS_WORD_T *c, int nw)
{ BITS_MINUS2(a, b, c, nw) return a; }


/* add element `i' into the set `bs' */
#define BITS_ADD(bs, i)         BITS_OR1(bs, BITS_GETNB(i), BITS_GETNQ(i))

/* remove element `i' from the set `bs' */
#define BITS_REMOVE(bs, i)      BITS_MINUS1(bs, BITS_GETNB(i), BITS_GETNQ(i))

/* flip the state of element `i' in the set `bs' */
#define BITS_FLIP(bs, i)        BITS_XOR1(bs, BITS_GETNB(i), BITS_GETNQ(i))

/* test if the set `bs' has element `i' */
#define BITS_HAS(bs, i)         (int) (((bs)[BITS_GETNQ(i)] >> BITS_GETNR(i)) & 0x1)

/* test if the set `bs' has bit `bi' at offset `iq' */
#define BITS_HASBIT(bs, bi, iq) (((bs)[iq] & (bi)) != 0)

/* right shift one bit */
#define BITS_RSHIFT1(a, b, nw) { int iw_; \
  for (iw_ = 0; iw_ < nw - 1; iw_++) { \
    a[iw_] = (b[iw_] >> 1) | (b[iw_ + 1] << (WORDBITS_ - 1)); \
  } a[iw_] = b[iw_] >> 1; }

/* left shift one bit */
#define BITS_LSHIFT1(a, b, nw) { int iw_; \
  for (iw_ = nw - 1; iw_ > 0; iw_--) { \
    a[iw_] = (b[iw_] << 1) | ( b[iw_ - 1] >> (WORDBITS_ - 1) ); \
  } a[0] = b[0] << 1; }



/* make a set of lowest k - 1 bits, without the ith element */
#define BITS_MKINVSET(bs, k, nw, i)  \
  BITS_MKBITSMASK(bs, k, nw) BITS_REMOVE(bs, i)

INLINE BITS_WORD_T *bits_mkinvset(BITS_WORD_T *bs, int n, int nw, int i)
{ BITS_MKINVSET(bs, n, nw, i) return bs; }

/* make a set of the lowest k - 2 bits, without the ith and jth elements */
#define BITS_MKINVSET2(bs, k, nw, i, j)  \
  BITS_MKBITSMASK(bs, k, nw) BITS_REMOVE(bs, i) BITS_REMOVE(bs, j)

INLINE BITS_WORD_T *bits_mkinvset2(BITS_WORD_T *bs, int n, int nw, int i, int j)
{ BITS_MKINVSET2(bs, n, nw, i, j) return bs; }



#define bits_print(bs, n, name) bits_fprint(bs, n, stdout, name)

/* print a set as star or blank pattern */
INLINE void bits_fprint(BITS_WORD_T *bs, int n, FILE *fp, const char *name)
{
  int i;

  if (name != NULL) fprintf(fp, "%-8s: [", name);
  for (i = 0; i < n; i++)
    fprintf(fp, "%c", BITS_HAS(bs, i) ? '*' : ' ');
  fprintf(fp, "]\n");
}



#define bits_printn(bs, n, name) bits_fprintn(bs, n, stdout, name)

/* print a set as numbers */
INLINE void bits_fprintn(BITS_WORD_T *bs, int n, FILE *fp, const char *name)
{
  int i;

  if (name != NULL) fprintf(fp, "%-8s:", name);
  for (i = 0; i < n; i++)
    if ( BITS_HAS(bs, i) )
      fprintf(fp, " %d", i);
  fprintf(fp, "\n");
}


#endif /* BITS_H__ */

