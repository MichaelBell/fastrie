#include"global.h"
#include <chrono>
#include <assert.h>
#include "tsqueue.hpp"
#include <math.h>
#include <primesieve.hpp>

#include <immintrin.h>

union xmmreg
{
  uint32_t v[4];
  __m128i m128;
};

#define zeroesBeforeHashInPrime	8

#define DEBUG 0

#define DEBUG_TIMING 0

#if DEBUG
#define DPRINTF(fmt, args...) do { printf("line %d: " fmt, __LINE__, ##args); fflush(stdout); } while(0)
#else
#define DPRINTF(fmt, ...) do { } while(0)
#endif

static constexpr int NONCE_REGEN_SECONDS = 195;
static constexpr uint32_t riecoin_sieveBits = 23; /* normally 22. 8 million, or 1MB, tuned for Haswell L3 */
static constexpr uint32_t riecoin_sieveSize = (1UL<<riecoin_sieveBits); /* 1MB, tuned for L3 of Haswell */
static constexpr uint32_t riecoin_sieveWords64 = riecoin_sieveSize/64;
static constexpr uint32_t riecoin_sieveWords128 = riecoin_sieveSize/128;

bool search_for_shares;
uint32_t riecoin_primeTestLimit;
uint32_t num_entries_per_segment = 0;
int N_THREADS = 4;
int N_SIEVE_WORKERS = 2;

uint32 riecoin_primorialNumber = 40; /* 15 is the 64 bit limit */
/* 39 gives us 2^36 increments */
/* Based on the primorial (40 is 226 bits), we only have about 2^29
 * increments before overflowing the 256 bit nonce field in Riecoin.
 * Each loop goes through riecoin_sieveSize increments, which means
 * that our max loop iter count is 2^29/riecoin_sieveSize.  With
 * the current settings of 8M for sieveSize, this means 64
 * iterations.
 */
static constexpr uint64_t max_increments = (1ULL<<29); /* 36 = 39.  29 = 40 */
static constexpr uint64_t maxiter = (max_increments/riecoin_sieveSize);

//static const uint32_t primorial_offset = 97;
static constexpr uint32_t primorial_offset = 16057; /* For > 26 or so */

//static uint32 primeTupleBias[6] = {0,4,6,10,12,16};
static uint32 primeTupleOffset[6] = {0, 4, 2, 4, 2, 4};

static constexpr uint32_t riecoin_denseLimit = 16384; /* A few cachelines */
uint32_t* riecoin_primeTestTable;
uint32_t riecoin_primeTestSize;
uint32_t riecoin_primeTestStoreOffsetsSize;
uint32_t *inverts;
mpz_t  z_primorial;
uint32_t startingPrimeIndex;

static constexpr int WORK_INDEXES = 64;
uint32_t *segment_counts;

enum JobType { TYPE_CHECK, TYPE_MOD, TYPE_SIEVE };

struct riecoinPrimeTestWork {
  /* prime = z_target+z_remainderPrimorial + ((loop * riecoin_sieveSize) + i) * z_primorial */
  /* Shared r/o access to z_target, z_remainderPrimorial, sievesize, primorial */
  JobType type;
  union {
    struct {
      uint32_t loop; 
      uint32_t n_indexes;
      uint32_t indexes[WORK_INDEXES];
    } testWork;
    struct {
      uint32_t start;
      uint32_t end;
    } modWork;
    struct {
      uint32_t start;
      uint32_t end;
      uint32_t sieveId;
    } sieveWork;
  };
};

ts_queue<riecoinPrimeTestWork, 1024> verifyWorkQueue;
ts_queue<int, 3096> workerDoneQueue;
ts_queue<int, 3096> testDoneQueue;
thread_local bool i_am_master = false;
bool there_is_a_master = false;
CRITICAL_SECTION master_lock;
CRITICAL_SECTION bucket_lock; /* careful */

/* These are globals that are written only by the 'master' thread,
 * but that are read by all of the verifiers */

mpz_t z_verify_target, z_verify_remainderPrimorial;
minerRiecoinBlock_t* verify_block;

void riecoin_init(uint64_t sieveMax, int numThreads, bool solo)
{
        N_THREADS = numThreads;
	N_SIEVE_WORKERS = std::max(1, numThreads/4);
	N_SIEVE_WORKERS = std::min(N_SIEVE_WORKERS, 8);
        InitializeCriticalSection(&master_lock);
	InitializeCriticalSection(&bucket_lock);
        mpz_init(z_verify_target);
	mpz_init(z_verify_remainderPrimorial);

	search_for_shares = !solo;
        riecoin_primeTestLimit = sieveMax;
	printf("Generating table of small primes for Riecoin...\n");
	// generate prime table
	{
		std::vector<uint32_t> pvec;
		primesieve::generate_primes(riecoin_primeTestLimit, &pvec);
		riecoin_primeTestTable = (uint32_t*)malloc(sizeof(uint32)*pvec.size());
		for (size_t i = 0; i < pvec.size(); i++) {
			riecoin_primeTestTable[i] = pvec[i];
		}
		riecoin_primeTestSize = pvec.size();
	}

	DPRINTF("Table with %d entries generated\n", riecoin_primeTestSize);

	// generate primorial
	mpz_init_set_ui(z_primorial, riecoin_primeTestTable[0]);
	for(uint32_t i=1; i<riecoin_primorialNumber; i++)
	{
		mpz_mul_ui(z_primorial, z_primorial, riecoin_primeTestTable[i]);
	}
#if DEBUG
	gmp_printf("z_primorial: %Zd\n", z_primorial);
#endif
	inverts = (uint32_t *)calloc(sizeof(uint32_t), riecoin_primeTestSize);
	if (inverts == NULL) {
	  perror("could not malloc inverts");
	  exit(-1);
	}

	mpz_t z_tmp, z_p;
	mpz_init(z_tmp);
	mpz_init(z_p);
	for (uint32_t i = 5; i < riecoin_primeTestSize; i++) {
	  mpz_set_ui(z_p, riecoin_primeTestTable[i]);
	  mpz_invert(z_tmp, z_primorial, z_p);
	  inverts[i] = mpz_get_ui(z_tmp);
	}
	mpz_clear(z_p);
	mpz_clear(z_tmp);

	uint64_t high_segment_entries = 0;
	double high_floats = 0.0;
	riecoin_primeTestStoreOffsetsSize = 0;
	for (uint32_t i = 5; i < riecoin_primeTestSize; i++) {
	  uint32_t p = riecoin_primeTestTable[i];
	  if (p < max_increments) {
	    riecoin_primeTestStoreOffsetsSize++;
	  } else {
	    high_floats += ((6.0f * max_increments) / (double)p);
	  }
	}
	high_segment_entries = ceil(high_floats);
	if (high_segment_entries == 0) {
	  num_entries_per_segment = 1;
	} else {
	  num_entries_per_segment = high_segment_entries/maxiter + 4; /* rounding up a bit */
	  num_entries_per_segment = (num_entries_per_segment + (num_entries_per_segment>>3));
	}
	segment_counts = (uint32_t *)malloc(sizeof(uint32_t) * maxiter);
	if (segment_counts == NULL) {
	  perror("could not malloc segment_counts");
	  exit(-1);
	}
}

typedef uint32_t sixoff[6];

thread_local uint8_t* riecoin_sieve = NULL;
sixoff *offsets = NULL;
uint32_t *segment_hits[maxiter];

uint8_t **riecoin_sieves;


uint32 _getHexDigitValue(uint8 c)
{
	if( c >= '0' && c <= '9' )
		return c-'0';
	else if( c >= 'a' && c <= 'f' )
		return c-'a'+10;
	else if( c >= 'A' && c <= 'F' )
		return c-'A'+10;
	return 0;
}

/*
 * Parses a hex string
 * Length should be a multiple of 2
 */
void debug_parseHexString(char* hexString, uint32 length, uint8* output)
{
	uint32 lengthBytes = length / 2;
	for(uint32 i=0; i<lengthBytes; i++)
	{
		// high digit
		uint32 d1 = _getHexDigitValue(hexString[i*2+0]);
		// low digit
		uint32 d2 = _getHexDigitValue(hexString[i*2+1]);
		// build byte
		output[i] = (uint8)((d1<<4)|(d2));	
	}
}

void debug_parseHexStringLE(char* hexString, uint32 length, uint8* output)
{
	uint32 lengthBytes = length / 2;
	for(uint32 i=0; i<lengthBytes; i++)
	{
		// high digit
		uint32 d1 = _getHexDigitValue(hexString[i*2+0]);
		// low digit
		uint32 d2 = _getHexDigitValue(hexString[i*2+1]);
		// build byte
		output[lengthBytes-i-1] = (uint8)((d1<<4)|(d2));	
	}
}

inline void silly_sort_indexes(uint32_t indexes[6]) {
  for (int i = 0; i < 5; i++) {
    for (int j = i+1; j < 6; j++) {
      if (indexes[j] < indexes[i]) {
	std::swap(indexes[i], indexes[j]);
      }
    }
  }
}

#define PENDING_SIZE 32

inline void init_pending(uint32_t pending[PENDING_SIZE]) {
  for (int i = 0; i < PENDING_SIZE; i++) {
    pending[i] = 0;
  }
}

inline void add_to_pending(uint8_t *sieve, uint32_t pending[PENDING_SIZE], uint32_t &pos, uint32_t ent) {
  __builtin_prefetch(&(sieve[ent>>3]));
  uint32_t old = pending[pos];
  if (old != 0) {
    assert(old < riecoin_sieveSize);
    sieve[old>>3] |= (1<<(old&7));
  }
  pending[pos] = ent;
  pos++;
  pos &= PENDING_SIZE-1;
}

inline void add6_to_pending(uint8_t *sieve, uint32_t pending[PENDING_SIZE], uint32_t &pos, xmmreg ent1, xmmreg ent2) {
  xmmreg addr;
#define ADD_REG(reg, cpts) \
  addr.m128 = _mm_srli_epi32(reg.m128, 3); \
  for (size_t i = 0; i < cpts; ++i) \
  { \
    if (reg.v[i] < riecoin_sieveSize) { \
      uint32_t old = pending[pos]; \
      if (old != 0) { \
        assert(old < riecoin_sieveSize); \
        sieve[old>>3] |= (1<<(old&7)); \
      } \
      __builtin_prefetch(&(sieve[addr.v[i]])); \
      pending[pos] = reg.v[i]; \
      pos++; \
      pos &= PENDING_SIZE-1; \
    } \
  }
  ADD_REG(ent1, 2)
  ADD_REG(ent2, 4)
#undef ADD_REG
}

void put_offsets_in_segments(uint32_t *offsets, int n_offsets) {
  EnterCriticalSection(&bucket_lock);
  for (int i = 0; i < n_offsets; i++) {
    uint32_t index = offsets[i];
    uint32_t segment = index>>riecoin_sieveBits;
    uint32_t sc = segment_counts[segment]++;
    if (sc >= num_entries_per_segment) { 
      printf("EEEEK segment %u  %u with index %u is > %u\n", segment, sc, index, num_entries_per_segment); exit(-1);
    }
    segment_hits[segment][sc] = index & (riecoin_sieveSize-1);
  }
  LeaveCriticalSection(&bucket_lock);
}

thread_local uint32_t *offset_stack = NULL;

void update_remainders(uint32_t start_i, uint32_t end_i) {
  mpz_t tar;
  mpz_init_set(tar, z_verify_target);
  mpz_add(tar, tar, z_verify_remainderPrimorial);
  int n_offsets = 0;
  static const int OFFSET_STACK_SIZE = 16384;
  if (offset_stack == NULL) {
    offset_stack = (uint32_t *)malloc(sizeof(uint32_t) * OFFSET_STACK_SIZE);
  }

  for (auto i = start_i; i < end_i; i++) {
    uint32_t p = riecoin_primeTestTable[i];
    uint32_t remainder = mpz_tdiv_ui(tar, p);
    bool is_once_only = false;

    /* Also update the offsets unless once only */
    if (p >= max_increments) {
      is_once_only = true;
    }
     
    uint64_t inverted = inverts[i];
    uint64_t pa = p - remainder;
    uint64_t index = pa*inverted;
    index %= p;
    uint64_t inverted2 = inverted << 1;
    if (inverted2 > p) inverted2 -= p;
    uint64_t inverted4 = inverted2 << 1;
    if (inverted4 > p) inverted4 -= p;
    if (!is_once_only) {
      offsets[i][0] = index;
      if (index < inverted4) index += p;
      index -= inverted4;
      offsets[i][1] = index;
      for (uint32_t f = 2; f < 6; f += 2)
      {
        if (index < inverted2) index += p;
        index -= inverted2;
        offsets[i][f] = index;
        if (index < inverted4) index += p;
        index -= inverted4;
        offsets[i][f+1] = index;
      }
    }
    else {
      if (n_offsets + 6 >= OFFSET_STACK_SIZE)
      {
        put_offsets_in_segments(offset_stack, n_offsets);
	n_offsets = 0;
      }
      if (index < max_increments) offset_stack[n_offsets++] = index;
      if (index < inverted4) index += p;
      index -= inverted4;
      if (index < max_increments) offset_stack[n_offsets++] = index;
      for (uint32_t f = 0; f < 2; ++f)
      {
        if (index < inverted2) index += p;
        index -= inverted2;
        if (index < max_increments) offset_stack[n_offsets++] = index;
        if (index < inverted4) index += p;
        index -= inverted4;
        if (index < max_increments) offset_stack[n_offsets++] = index;
      }
    }
  }
  if (n_offsets > 0) {
    put_offsets_in_segments(offset_stack, n_offsets);
    n_offsets = 0;
  }
  mpz_clear(tar);
}


void process_sieve(uint8_t *sieve, uint32_t start_i, uint32_t end_i) {
  uint32_t pending[PENDING_SIZE];
  uint32_t pending_pos = 0;
  init_pending(pending);
  
  xmmreg opnomax;
  opnomax.m128 = _mm_set1_epi32(riecoin_sieveSize);

  for (unsigned int i = start_i; i < end_i; i++) {
    uint32_t pno = i+startingPrimeIndex;
    xmmreg p;
    xmmreg opno1, opno2, nextopno1, nextopno2;
    xmmreg cmpres1, cmpres2;
    p.m128 = _mm_set1_epi32(riecoin_primeTestTable[pno]);
    opno1.m128 = _mm_loadu_si128((__m128i const *)&offsets[pno][0]);
    opno2.m128 = _mm_loadu_si128((__m128i const *)&offsets[pno][2]);

    while (true) {
      cmpres1.m128 = _mm_cmpgt_epi32(opnomax.m128, opno1.m128);
      cmpres2.m128 = _mm_cmpgt_epi32(opnomax.m128, opno2.m128);
      if ((_mm_movemask_epi8(cmpres1.m128) == 0) && (_mm_movemask_epi8(cmpres2.m128) == 0)) break;
      add6_to_pending(sieve, pending, pending_pos, opno1, opno2);
      nextopno1.m128 = _mm_and_si128(cmpres1.m128, p.m128);
      nextopno2.m128 = _mm_and_si128(cmpres2.m128, p.m128);
      opno1.m128 = _mm_add_epi32(opno1.m128, nextopno1.m128);
      opno2.m128 = _mm_add_epi32(opno2.m128, nextopno2.m128);
    }
    opno1.m128 = _mm_sub_epi32(opno1.m128, opnomax.m128);
    opno2.m128 = _mm_sub_epi32(opno2.m128, opnomax.m128);
    offsets[pno][0] = opno1.v[0];
    offsets[pno][1] = opno1.v[1];
    _mm_storeu_si128((__m128i*)&offsets[pno][2], opno2.m128);
  }

  for (unsigned int i = 0; i < PENDING_SIZE; i++) {
    const uint32_t old = pending[i];
    if (old != 0) {
      assert(old < riecoin_sieveSize);
      sieve[old>>3] |= (1<<(old&7));
    }
  }
}

#if DEBUG_TIMING
volatile float modTime = 0.f;
volatile float sieveTime = 0.f;
volatile float checkTime = 0.f;
#endif

void verify_thread() {
  /* Check for a prime cluster.  A "share" on ypool is any
   * four or more of the elements prime, but for speed,
   * check further only if the first passes the primality
   * test.  The first test is the bottleneck for the
   * miner.
   *
   * Uses the fermat test - jh's code noted that it is slightly faster.
   * Could do an MR test as a follow-up, but the server can do this too
   * for the one-in-a-whatever case that Fermat is wrong.
   */

  mpz_t z_ft_r, z_ft_n, z_temp, z_temp2, z_ploop;
  mpz_inits(z_ft_r, z_ft_n, z_temp, z_temp2, z_ploop, NULL);
  mpz_t z_ft_b;
  mpz_init_set_ui(z_ft_b, 2);

  while (1) {
    auto job = verifyWorkQueue.pop_front();

#if DEBUG_TIMING
    auto start = std::chrono::system_clock::now();
#endif

    if (job.type == TYPE_MOD) {
      update_remainders(job.modWork.start, job.modWork.end);
      DPRINTF("wdq-push-mod %d-%d start\n", job.modWork.start, job.modWork.end);

      workerDoneQueue.push_back(1);
      DPRINTF("wdq-push-mod %d-%d done\n", job.modWork.start, job.modWork.end);

#if DEBUG_TIMING
      auto end = std::chrono::system_clock::now();
      auto dur = std::chrono::duration_cast<std::chrono::microseconds>(end-start).count();
      modTime += dur * 0.000001f;
#endif

      continue;
    }

    if (job.type == TYPE_SIEVE) {
      process_sieve(riecoin_sieves[job.sieveWork.sieveId], job.sieveWork.start, job.sieveWork.end);
      DPRINTF("wdq-push-sieve %d-%d start\n", job.modWork.start, job.modWork.end);
      workerDoneQueue.push_back(1);
      DPRINTF("wdq-push-sieve %d-%d done\n", job.modWork.start, job.modWork.end);

#if DEBUG_TIMING
      auto end = std::chrono::system_clock::now();
      auto dur = std::chrono::duration_cast<std::chrono::microseconds>(end-start).count();
      sieveTime += dur * 0.000001f;
#endif

      continue;
    }
    /* fallthrough:  job.type == TYPE_CHECK */
    if (job.type == TYPE_CHECK) {
      mpz_mul_ui(z_ploop, z_primorial, job.testWork.loop * riecoin_sieveSize);
      mpz_add(z_ploop, z_ploop, z_verify_remainderPrimorial);
      mpz_add(z_ploop, z_ploop, z_verify_target);
      
      for (unsigned int idx = 0; idx < job.testWork.n_indexes; idx++) {
	int nPrimes = 0;
	mpz_mul_ui(z_temp, z_primorial, job.testWork.indexes[idx]);
	mpz_add(z_temp, z_temp, z_ploop);
	
	mpz_sub_ui(z_ft_n, z_temp, 1);
	mpz_powm(z_ft_r, z_ft_b, z_ft_n, z_temp);
	if (mpz_cmp_ui(z_ft_r, 1) != 0) {
	  continue;
	}
	
	nPrimes++;
	totalChainCount[nPrimes]++;
	
	/* New definition of shares:  Any 4+ valid primes.  Search method 
	 * is for 1st + any 3 to avoid doing too much primality testing.
	 */
	
	/* Note start at 1 - we've already tested bias 0 */
	for (int i = 1; i < 6; i++) {
	  mpz_add_ui(z_temp, z_temp, primeTupleOffset[i]);
	  mpz_sub_ui(z_ft_n, z_temp, 1);
	  mpz_powm(z_ft_r, z_ft_b, z_ft_n, z_temp);
	  if (mpz_cmp_ui(z_ft_r, 1) == 0) {
	    nPrimes++;
	    totalChainCount[nPrimes]++;
	  }
	  else if (!search_for_shares) break;
	  int candidatesRemaining = 5-i;
	  if ((nPrimes + candidatesRemaining) < 4) { continue; }
	}
	
	/* The statistics are a little confusing because of the interaction
	 * with early-exit above.  They overcount relative to finding consecutive
	 * primes, but undercount relative to counting all primes.  But they're
	 * still useful for benchmarking within a variant of the program with
	 * all else held equal. */
	
	if (nPrimes < 4) continue;
	
	mpz_set(z_temp, z_primorial);
	mpz_mul_ui(z_temp, z_temp, job.testWork.loop);
	mpz_mul_ui(z_temp, z_temp, riecoin_sieveSize);
	mpz_set(z_temp2, z_primorial);
	mpz_mul_ui(z_temp2, z_temp2, job.testWork.indexes[idx]);
	mpz_add(z_temp, z_temp, z_temp2);
	mpz_add(z_temp, z_temp, z_verify_remainderPrimorial);
	mpz_add(z_temp, z_temp, z_verify_target);
	
	mpz_sub(z_temp2, z_temp, z_verify_target); // offset = tested - target
	// submit share
	uint8 nOffset[32];
	memset(nOffset, 0x00, 32);
#if defined _WIN64 || __X86_64__
	for(uint32 d=0; d<std::min(32/8, z_temp2->_mp_size); d++)
	  {
	    *(uint64*)(nOffset+d*8) = z_temp2->_mp_d[d];
	  }
#elif defined _WIN32 
	for(uint32 d=0; d<std::min(32/4, z_temp2->_mp_size); d++)
	  {
	    *(uint32*)(nOffset+d*4) = z_temp2->_mp_d[d];
	  }
#elif defined __GNUC__
#ifdef	__x86_64__
	for(uint32 d=0; d<std::min(32/8, z_temp2->_mp_size); d++)
	  {
	    *(uint64*)(nOffset+d*8) = z_temp2->_mp_d[d];
	  }
#else  
	for(uint32 d=0; d<std::min(32/4, z_temp2->_mp_size); d++)
	  {
	    *(uint32*)(nOffset+d*4) = z_temp2->_mp_d[d];
	  }
#endif
#endif
	totalShareCount++;
	xptMiner_submitShare(verify_block, nOffset);
      }
      DPRINTF("tdq-push start\n");
      testDoneQueue.push_back(1);
      DPRINTF("tdq-push done\n");

#if DEBUG_TIMING
      auto end = std::chrono::system_clock::now();
      auto dur = std::chrono::duration_cast<std::chrono::microseconds>(end-start).count();
      checkTime += dur * 0.000001f;
#endif
    }
  }
}


void riecoin_process(minerRiecoinBlock_t* block)
{
#if DEDUB_TIMING
  modTime = sieveTime = checkTime = 0;
#endif

	uint32 searchBits = block->targetCompact;

	if (!there_is_a_master) {
	  EnterCriticalSection(&master_lock);
	  if (!there_is_a_master) {
	    there_is_a_master = true;
	    i_am_master = true;
	    DPRINTF("I have become the master thread\n");
	  }
	  LeaveCriticalSection(&master_lock);
	}

	if (!i_am_master) {
	  DPRINTF("I am a verifier thread\n");
	  verify_thread(); /* Runs forever */
	  printf("If you see this, life is VERY bad\n");fflush(stdout);
	  return;
	}

	if( !riecoin_sieve ) {
	  riecoin_sieve = (uint8*)malloc(riecoin_sieveSize/8);
	  riecoin_sieves = (uint8_t **)malloc(sizeof(uint8_t *) * N_SIEVE_WORKERS);
	  DPRINTF("allocating sieves for %d workers\n", N_SIEVE_WORKERS);
	  for (int i = 0; i < N_SIEVE_WORKERS; i++) {
	    riecoin_sieves[i] = (uint8_t*)malloc(riecoin_sieveSize/8);
	    if (riecoin_sieves[i] == NULL) {
	      perror("malloc error on sieve allocation");
	      printf("Eek:  Could not allocate sieve %d!\n", i);
	      exit(-1);
	    }
	  }
	  size_t offsize = sizeof(sixoff) * (riecoin_primeTestStoreOffsetsSize+1024);
	  offsets = (sixoff *)malloc(offsize);
	  if (offsets == NULL) {
	    perror("malloc error on offsets allocation");
	    printf("could not allocate %lu bytes for offsets\n", offsize);
	    exit(-1);
	  }
	  memset(offsets, 0, offsize);
	  for (unsigned int i = 0; i < maxiter; i++) {
	    segment_hits[i] = (uint32_t *)malloc(sizeof(uint32_t) * num_entries_per_segment);
	    if (segment_hits[i] == NULL) {
	      perror("malloc error on segment_hits allocation");
	      printf("could not allocate %lu bytes for segment hits %d\n", sizeof(uint32_t) * num_entries_per_segment, i);
	      exit(-1);
	    }
	  }
	}
	uint8* sieve = riecoin_sieve;

	time_t start_time = time(NULL);
#if DEBUG
	auto start = std::chrono::system_clock::now();
#endif

	// test data
	// getblock 16ee31c116b75d0299dc03cab2b6cbcb885aa29adf292b2697625bc9d28b2b64
	//debug_parseHexStringLE("c59ba5357285de73b878fed43039a37f85887c8960e66bcb6e86bdad565924bd", 64, block->merkleRoot);
	//block->version = 2;
	//debug_parseHexStringLE("c64673c670fb327c2e009b3b626d2def01d51ad4131a7a1040e9cef7bfa34838", 64, block->prevBlockHash);
	//block->nTime = 1392151955;
	//block->nBits = 0x02013000;
	//debug_parseHexStringLE("0000000000000000000000000000000000000000000000000000000070b67515", 64, block->nOffset);
	// generate PoW hash (version to nBits)
	uint8 powHash[32];
	sha256_ctx ctx;
	sha256_init(&ctx);
	sha256_update(&ctx, (uint8*)block, 80);
	sha256_final(&ctx, powHash);
	sha256_init(&ctx);
	sha256_update(&ctx, powHash, 32);
	sha256_final(&ctx, powHash);
	// generatePrimeBase
	uint32* powHashU32 = (uint32*)powHash;
	mpz_t z_target;
	mpz_t z_temp;
	mpz_init(z_temp);
	mpz_t z_remainderPrimorial;
	mpz_init(z_remainderPrimorial);

	mpz_init_set_ui(z_target, 1);
	mpz_mul_2exp(z_target, z_target, zeroesBeforeHashInPrime);
	for(uint32 i=0; i<256; i++)
	{
		mpz_mul_2exp(z_target, z_target, 1);
		if( (powHashU32[i/32]>>(i%32))&1 )
			z_target->_mp_d[0]++;
	}
	unsigned int trailingZeros = searchBits - 1 - zeroesBeforeHashInPrime - 256;
	mpz_mul_2exp(z_target, z_target, trailingZeros);
	// find first offset where target%primorial = primorial_offset

	mpz_tdiv_r(z_remainderPrimorial, z_target, z_primorial);
	mpz_abs(z_remainderPrimorial, z_remainderPrimorial);
	mpz_sub(z_remainderPrimorial, z_primorial, z_remainderPrimorial);
	mpz_tdiv_r(z_remainderPrimorial, z_remainderPrimorial, z_primorial);
	mpz_abs(z_remainderPrimorial, z_remainderPrimorial);
	mpz_add_ui(z_remainderPrimorial, z_remainderPrimorial, primorial_offset);

	mpz_add(z_temp, z_target, z_remainderPrimorial);

	uint32_t primeIndex = riecoin_primorialNumber;	
	
	{
	  startingPrimeIndex = primeIndex;
	  mpz_set(z_verify_target, z_target);
	  mpz_set(z_verify_remainderPrimorial, z_remainderPrimorial);
	  verify_block = block;
	}

	//	  printf("Boss is working1\n");fflush(stdout);

	uint32_t n_dense = 0;
	uint64_t n_sparse = 0;

	/* Get some stuff done while the workers are working */
	for (unsigned int i = 0; i < maxiter; i++) {
	  segment_counts[i] = 0;
	}


	uint32_t incr = riecoin_primeTestSize/128;
	riecoinPrimeTestWork wi;
	//	  printf("Boss is working2\n");fflush(stdout);
	wi.type = TYPE_MOD;
	int n_workers = 0;
	for (auto base = primeIndex; base < riecoin_primeTestSize; base += incr) {
	  uint32_t lim = std::min(riecoin_primeTestSize, base+incr);
	  wi.modWork.start = base;
	  wi.modWork.end = lim;
	  DPRINTF("master vwq-push start\n");
	  verifyWorkQueue.push_back(wi);
	  DPRINTF("master vwq-push done\n");
	  n_workers++;
	}


	DPRINTF("master wdq-wait complete start\n");

	for (int i = 0; i < n_workers; i++) {
	  workerDoneQueue.pop_front();
	}
	DPRINTF("master wdq-wait complete done\n");

	//printf("primeIndex: %d  is %d\n", primeIndex, riecoin_primeTestTable[primeIndex]);
	/* XXX - move this all into init - it's only done once */
	for( ; primeIndex < riecoin_primeTestSize; primeIndex++)
	{
	  uint32 p = riecoin_primeTestTable[primeIndex];
	  if (p < riecoin_denseLimit) {
	    n_dense++;
	  } else if (p < max_increments) {
	    n_sparse++;
	  }
	}
	uint32_t round = 8 - (n_dense & 0x7);
	n_dense += round;
	n_sparse -= round;

	//printf("\n");
#if DEBUG

	auto end = std::chrono::system_clock::now();
	auto dur = std::chrono::duration_cast<std::chrono::milliseconds>(end-start).count();
	printf("Loop start offset compute time:  %ld ms\n", dur);
#endif

	/* Main processing loop:
	 * 1)  Sieve "dense" primes;
	 * 2)  Sieve "sparse" primes;
	 * 3)  Sieve "so sparse they happen at most once" primes;
	 * 4)  Scan sieve for candidates, test, report
	 */

	uint32 countCandidates = 0;
	uint32_t outstandingTests = 0;

	for (unsigned int loop = 0; loop < maxiter; loop++) {
	    __sync_synchronize(); /* gcc specific - memory barrier for checking height */
	    if( block->height != monitorCurrentBlockHeight ) {
	      break;
	    }
	    time_t cur_time = time(NULL);
	    if ((cur_time - start_time) > NONCE_REGEN_SECONDS) {
	      break;
	    }
	    DPRINTF("Loop %d after %d seconds\n", loop, (cur_time-start_time));

	    for (int i = 0; i < N_SIEVE_WORKERS; i++) {
	      memset(riecoin_sieves[i], 0, riecoin_sieveSize/8);
	    }
	    
	    DPRINTF("done zeroing sieves\n");

	    wi.type = TYPE_SIEVE;
	    n_workers = 0;
	    incr = ((n_sparse)/N_SIEVE_WORKERS);
            round = 8 - (incr & 0x7);
            incr += round;
	    int which_sieve = 0;
	    //printf("n_dense: %u  sparse: %u\n", n_dense, n_sparse);
	    for (auto base = n_dense; base < (n_dense+n_sparse); base += incr) {
	      uint32_t lim = std::min((uint32_t)(n_dense+n_sparse), (uint32_t)base+incr);
	      if ((lim + 1000) > (n_dense+n_sparse)) {
		lim = (n_dense+n_sparse);
	      }
	      DPRINTF("Sieving from %d to %d\n", base, lim);
	      wi.sieveWork.start = base;
	      wi.sieveWork.end = lim;
	      wi.sieveWork.sieveId = which_sieve;
	      /* Need to do something for thread to sieve affinity */
	      DPRINTF("master vwq-push start to sieve %d\n", which_sieve);
	      verifyWorkQueue.push_front(wi);
	      DPRINTF("master vwq-push done\n");
	      which_sieve++;
	      which_sieve %= N_SIEVE_WORKERS;
	      n_workers++;
	      if ((lim + 1000) > (n_dense+n_sparse)) {
		break;
	      }
	    }

	    DPRINTF("zeroing master sieve\n");

	    memset(sieve, 0, riecoin_sieveSize/8);
	    for (unsigned int i = 0; i < n_dense; i++) {
	      uint32_t pno = i+startingPrimeIndex;
	      silly_sort_indexes(offsets[pno]);
	      uint32_t p = riecoin_primeTestTable[pno];
	      while (offsets[pno][5] < riecoin_sieveSize) {
#define MARK(x) sieve[offsets[pno][x]>>3] |= (1<<(offsets[pno][x]&7)); offsets[pno][x] += p;
		      MARK(0);
		      MARK(1);
		      MARK(2);
		      MARK(3);
		      MARK(4);
		      MARK(5);
	      }
	      for (uint32 f = 0; f < 6; f++) {
		auto opnof = offsets[pno][f];
		if (opnof < riecoin_sieveSize) {
		  sieve[opnof>>3] |= (1<<(opnof&7));
		  opnof += p;
		}
		offsets[pno][f] = (opnof - riecoin_sieveSize);
	      }
	    }

	    DPRINTF("master tdq clear start\n");

	    outstandingTests -= testDoneQueue.clear();

	    DPRINTF("master wdq-verify pop start\n");
	    for (int i = 0; i < n_workers; i++) {
	      workerDoneQueue.pop_front();
	    }
	    DPRINTF("master wdq-verify pop done\n");

	    if (workerDoneQueue.size() != 0) {
	      printf("ERROR:  workerDoneQueue has grown too large!  Report this error\n");fflush(stdout);
	      exit(-1);
	    }
	    __m128i *sp128 = (__m128i *)sieve;

	    for (unsigned int i = 0; i < riecoin_sieveWords128; i++) {
              __m128i s128;
              s128 = _mm_loadu_si128(&sp128[i]);
	      for (int j = 0; j < N_SIEVE_WORKERS; j++) {
		__m128i ws128;
                ws128 = _mm_loadu_si128(&((__m128i *)(riecoin_sieves[j]))[i]);
                s128 = _mm_or_si128(s128, ws128);
	      }
              _mm_storeu_si128(&sp128[i], s128);
	    }
	    DPRINTF("master - done merging sieves\n");
		  
	    //process_sieve(sieve, n_dense, (n_dense+n_sparse));

	    uint32_t pending[PENDING_SIZE];
	    init_pending(pending);
	    uint32_t pending_pos = 0;
	    for (uint32_t i = 0; i < segment_counts[loop]; i++) {
	      add_to_pending(sieve, pending, pending_pos, segment_hits[loop][i]);
	    }

	    DPRINTF("master - done adding segment counts to workers\n");

	    for (unsigned int i = 0; i < PENDING_SIZE; i++) {
	      uint32_t old = pending[i];
	      if (old != 0) {
		assert(old < riecoin_sieveSize);
		sieve[old>>3] |= (1<<(old&7));
	      }
	    }
	  

	    // scan for candidates
	    riecoinPrimeTestWork w;
	    w.testWork.n_indexes = 0;
	    w.testWork.loop = loop;
	    w.type = TYPE_CHECK;

	    bool do_reset = false;
	    uint64_t *sieve64 = (uint64_t *)sieve;
	    for(uint32 b = 0; !do_reset && b < riecoin_sieveWords64; b++) {
	      uint64_t sb = ~sieve64[b];
	      
	      int sb_process_count = 0;
	      while (sb != 0) {
		sb_process_count++;
		if (sb_process_count > 65) {
		  printf("Impossible:  process count too high.  Bug bug\n");
		  fflush(stdout);
		  exit(-1);
		}
		uint32_t highsb = __builtin_clzll(sb);
		uint32_t i = (b*64) + (63-highsb);
		sb &= ~(1ULL<<(63-highsb));
		
		countCandidates++;
		
		w.testWork.indexes[w.testWork.n_indexes] = i;
		w.testWork.n_indexes++;
		outstandingTests -= testDoneQueue.clear();

		if (w.testWork.n_indexes == WORK_INDEXES) {
		  DPRINTF("master vwq push new work start sb 0x%lx\n", sb);
		  verifyWorkQueue.push_back(w);
		  DPRINTF("master vwq push new work done\n");
		  w.testWork.n_indexes = 0;
		  outstandingTests++;
		}
		outstandingTests -= testDoneQueue.clear();
		
		/* Low overhead but still often enough */
		if( block->height != monitorCurrentBlockHeight ) {
		  outstandingTests -= verifyWorkQueue.clear();
		  do_reset = true;
		  break;
		}
	      }
	      DPRINTF("sieveWord done\n");

	    }
	    

	    if (w.testWork.n_indexes > 0) {
	      DPRINTF("master vwq push remaining new work start\n");
	      verifyWorkQueue.push_back(w);
	      outstandingTests++;
	      w.testWork.n_indexes = 0;
	      DPRINTF("master vwq push remaining new work done\n");
	    }

	}
	DPRINTF("Loop end:  vwQ size: %d\n", verifyWorkQueue.size());
	outstandingTests -= testDoneQueue.clear();
	DPRINTF("master waiting on all work completion\n");

	while (outstandingTests > 0) {
	  testDoneQueue.pop_front();
	  outstandingTests--;
	  if( block->height != monitorCurrentBlockHeight ) {
	    outstandingTests -= verifyWorkQueue.clear();
	  }
	}

	printf("Total candidates evaluated: %d\n", countCandidates);
#if DEBUG_TIMING
	float sumTime = (modTime + sieveTime + checkTime) * .01f;
	printf("Thread timing: %.3f %.3f %.3f\n", modTime / sumTime, sieveTime / sumTime, checkTime / sumTime);
#endif
	mpz_clears(z_target, z_temp, z_remainderPrimorial, NULL);
}
