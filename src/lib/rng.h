#ifndef RNG_H
#define RNG_H

// __________________________________________________________________________
// rng.h - a Random Number Generator Class
// rng.C - contains the non-inline class methods

// __________________________________________________________________________
// CAUTIONS:

// 1. Some of this code might not work correctly on 64-bit machines.  I
// have hacked the 32 bit version try and make it work, but the 64-bit
// version is not extensively tested.
//
// 2. This generator should NOT be used as in the following line.
// for (int32_t i = 0; i < 100; ++i) { RNG x; cout << x.uniform() << endl; }
// The problem is that each time through the loop, a new RNG 'x' is
// created, and that RNG is used to generate exactly one random number.
// While the results may be satisfactory, the class is designed to
// produce quality random numbers by having a single (or a few) RNGs
// called repeatedly.
// The better way to do the above loop is:
// RNG x; for (int32_t i = 0; i < 100; ++i) { cout << x.uniform() << endl; }

// __________________________________________________________________________
// This C++ code uses the simple, fast "KISS" (Keep It Simple Stupid)
// random number generator suggested by George Marsaglia in a Usenet
// posting from 1999.  He describes it as "one of my favorite
// generators".  It generates high-quality random numbers that
// apparently pass all commonly used tests for randomness.  In fact, it
// generates random numbers by combining the results of three simple
// random number generators that have different periods and are
// constructed from completely different algorithms.  It does not have
// the ultra-long period of some other generators - a "problem" that can
// be fixed fairly easily - but that seems to be its only potential
// problem.  The period is about 2^123.

// The KISS algorithm is only used directly in the function rand_int32.
// rand_int32 is then used (directly or indirectly) by every other
// member function of the class that generates random numbers.  For
// faster random numbers, one can redefine rand_int32 to return either
// WMC(), CONG(), or SHR3().  The speed will be two to three times
// faster, and the quality of the random numbers should be  sufficient
// for many purposes.  The three alternatives are comparable in terms of
// both speed and quality.

// The ziggurat method of Marsaglia is used to generate exponential and
// normal variates.  The method as well as source code can be found in
// the article "The Ziggurat Method for Generating Random Variables" by
// Marsaglia and Tsang, Journal of Statistical Software 5, 2000.

// The method for generating gamma variables appears in "A Simple Method
// for Generating Gamma Variables" by Marsaglia and Tsang, ACM
// Transactions on Mathematical Software, Vol. 26, No 3, Sep 2000, pages
// 363-372.
// __________________________________________________________________________

#include <cstdlib>
#include <cstdint>
#include <ctime>
#include <cmath>
#include <climits>
#include <vector>

using std::vector;

static const double PI = 3.1415926535897932;
static const double AD_l = 0.6931471805599453;
static const double AD_a = 5.7133631526454228;
static const double AD_b = 3.4142135623730950;
static const double AD_c = -1.6734053240284925;
static const double AD_p = 0.9802581434685472;
static const double AD_A = 5.6005707569738080;
static const double AD_B = 3.3468106480569850;
static const double AD_H = 0.0026106723602095;
static const double AD_D = 0.0857864376269050;

class RNG
{
 private:
  static uint32_t tm; // Used to ensure different RNGs have different seeds.
  uint32_t z, w, jsr, jcong; // Seeds

  uint32_t kn[128], ke[256];
  double wn[128], fn[128], we[256], fe[256];

 public:
  RNG()
  {
      init();
      zigset();
  }

  RNG(uint32_t x_) :
      z(x_), w(x_), jsr(x_), jcong(x_)
  {
      zigset();
  }

  RNG(uint32_t z_, uint32_t w_, uint32_t jsr_, uint32_t jcong_) :
      z(z_), w(w_), jsr(jsr_), jcong(jcong_)
  {
      zigset();
  }

  ~RNG()
  {
  }

  uint32_t znew()
  {
      return (z = 36969 * (z & 0xfffful) + (z >> 16));
  }

  uint32_t wnew()
  {
      return (w = 18000 * (w & 0xfffful) + (w >> 16));
  }

  uint32_t MWC()
  {
      return ((znew() << 16) + wnew());
  }

  uint32_t SHR3()
  {
      jsr ^= (jsr << 17);
      jsr ^= (jsr >> 13);
      return (jsr ^= (jsr << 5));
  }

  uint32_t CONG()
  {
      return (jcong = 69069 * jcong + 1234567);
  }

  uint32_t rand_int32()       // [0,2^32-1]
  {
      return ((MWC() ^ CONG()) + SHR3());
  }

  uint32_t rand_int()         // [0,2^32-1]
  {
      return ((MWC() ^ CONG()) + SHR3());
  }

  double RNOR()
  {
      int32_t h = int32_t(rand_int32()), i = h & 127;
      return (((uint32_t) std::abs(h) < kn[i]) ? h * wn[i] : nfix(h, i));
  }

  double REXP()
  {
      uint32_t j = rand_int32(), i = j & 255;
      return ((j < ke[i]) ? j * we[i] : efix(j, i));
  }

  double nfix(int32_t h, uint32_t i);

  double efix(uint32_t j, uint32_t i);

  void zigset();

  void init()
  {
      z = w = jsr = jcong = uint32_t(time(0)) + tm;
      tm += 123457;
  }

  void init(uint32_t z_, uint32_t w_, uint32_t jsr_, uint32_t jcong_)
  {
      z = z_;
      w = w_;
      jsr = jsr_;
      jcong = jcong_;
  }

  // For a faster but lower quality RNG, uncomment the following
  // line, and comment out the original definition of rand_int above.
  // In practice, the faster RNG will be fine for simulations
  // that do not simulate more than a few billion random numbers.
  // uint32_t rand_int() { return SHR3(); }

  int32_t rand_int31()          // [0,2^31-1]
  {
      return (int32_t) (rand_int32() >> 1);
  }

  double rand_closed01()     // [0,1]
  {
      return ((double) rand_int() / double(ULONG_MAX));
  }

  double rand_open01()       // (0,1)
  {
      return (((double) rand_int() + 1.0) / (ULONG_MAX + 2.0));
  }

  double rand_halfclosed01() // [0,1)
  {
      return ((double) rand_int() / (ULONG_MAX + 1.0));
  }

  double rand_halfopen01()   // (0,1]
  {
      return (((double) rand_int() + 1.0) / (ULONG_MAX + 1.0));
  }

  // Continuous Distributions
  double uniform(double x = 0.0, double y = 1.0)
  {
      return rand_closed01() * (y - x) + x;
  }

  double normal(double mu = 0.0, double sd = 1.0)
  {
      return RNOR() * sd + mu;
  }

  double exponential(double lambda = 1)
  {
      return REXP() / lambda;
  }

  double gamma(double shape = 1, double scale = 1);

  double chi_square(double df)
  {
      return gamma(df / 2.0, 0.5);
  }

  double beta(double a1, double a2)
  {
      const double x1 = gamma(a1, 1);
      return (x1 / (x1 + gamma(a2, 1)));
  }

  void uniform(vector<double> &res, double x = 0.0, double y = 1.0)
  {
      for (vector<double>::iterator i = res.begin(); i != res.end(); ++i)
      {
          *i = uniform(x, y);
      }
  }

  void normal(vector<double> &res, double mu = 0.0, double sd = 1.0)
  {
      for (vector<double>::iterator i = res.begin(); i != res.end(); ++i)
      {
          *i = normal(mu, sd);
      }
  }

  void exponential(vector<double> &res, double lambda = 1)
  {
      for (vector<double>::iterator i = res.begin(); i != res.end(); ++i)
      {
          *i = exponential(lambda);
      }
  }

  void gamma(vector<double> &res, double shape = 1, double scale = 1)
  {
      for (vector<double>::iterator i = res.begin(); i != res.end(); ++i)
      {
          *i = gamma(shape, scale);
      }
  }

  void chi_square(vector<double> &res, double df)
  {
      for (vector<double>::iterator i = res.begin(); i != res.end(); ++i)
      {
          *i = chi_square(df);
      }
  }

  void beta(vector<double> &res, double a1, double a2)
  {
      for (vector<double>::iterator i = res.begin(); i != res.end(); ++i)
      {
          *i = beta(a1, a2);
      }
  }

  // Discrete Distributions
  int32_t poisson(double mu);

  int32_t binomial(double p, int32_t n);

  void multinom(uint32_t n, const vector<double> &probs, vector<uint32_t> &samp);

  void multinom(uint32_t n, const double *prob, uint32_t K, uint32_t *samp);

  void poisson(vector<int32_t> &res, double lambda)
  {
      for (vector<int32_t>::iterator i = res.begin(); i != res.end(); ++i)
      {
          *i = poisson(lambda);
      }
  }

  void binomial(vector<int32_t> &res, double p, int32_t n)
  {
      for (vector<int32_t>::iterator i = res.begin(); i != res.end(); ++i)
      {
          *i = binomial(p, n);
      }
  }

}; // class RNG

#endif // RNG_H
