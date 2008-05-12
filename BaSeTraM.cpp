#include <boost/program_options.hpp>
#include <fstream>
#include <iostream>
#include <string>
#include <list>
#include <boost/iostreams/device/file.hpp>
#include <boost/filesystem.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/regex.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/math/special_functions/beta.hpp>
#include <boost/math/special_functions/binomial.hpp>
#include <boost/circular_buffer.hpp>
#include <boost/lambda/lambda.hpp>
#include <boost/lambda/bind.hpp>
#include <cmath>
#include <assert.h>
#include "../parsegenbank/GenbankParser.hpp"
#include <sys/time.h>
#include <time.h>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/date_time/microsec_time_clock.hpp>

#ifndef __GNUC__
#define USUALLY(x) x
#define RARELY(x) x
#else
#define USUALLY(x) __builtin_expect(!!(x), 1)
#define RARELY(x) __builtin_expect(!!(x), 0)
#endif

namespace po = boost::program_options;
namespace io = boost::iostreams;
namespace fs = boost::filesystem;
namespace ublas = boost::numeric::ublas;
namespace ll = boost::lambda;
namespace pt = boost::posix_time;

enum Base
{
  BASE_A = 0,
  BASE_C,
  BASE_G,
  BASE_T,
  BASE_N /* for undefined */
};

pt::ptime gStartupTime;

class Motif
{
public:
  Motif(const std::string& aAccession, ublas::matrix<double> aRawFreqD)
    : mAccession(aAccession)
  {
    uint32_t max_n = 0;

    mLength = aRawFreqD.size1();
    mBaseProb = new double[mLength * 4];

    for (uint32_t i = 0; i < mLength; i++)
    {
      uint32_t n_i = 0;
      uint32_t rawFreq[4];

      for (uint32_t b = 0; b < 4; b++)
      {
        rawFreq[b] = static_cast<uint32_t>(::round(aRawFreqD(i, b)));
        n_i += rawFreq[b];
      }

      max_n = std::max(n_i, max_n);

      for (uint32_t b = 0; b < 4; b++)
      {
        uint32_t f_ib(rawFreq[b]);
        mBaseProb[4 * i + b] =
          boost::math::beta<double>(f_ib + 2.0, n_i - f_ib + 1.0)
          /
          boost::math::beta<double>(f_ib + 1.0, n_i - f_ib + 1.0);
      }
    }

    // XXX should we make the prior configurable or something? For now, just
    // assume that of 4 billion possible positions, we have found 1/100th of
    // them. This is unavoidably a total guess.
    mH1 = max_n * (100.0 / 4.0E9);
    mH0 = 1.0 - mH1;
  }

  ~Motif()
  {
    delete [] mBaseProb;
  }

  double
  getPosterior(
               double pD_given_H0, double pD_given_H1
              )
  {

    return mH1 * pD_given_H1 / (mH1 * pD_given_H1 + mH0 * pD_given_H0);
  }

  double getAlternativeProbability(uint32_t position, Base b)
  {
    return mBaseProb[(position << 2) + b];
  }

  uint32_t
  getLength() const
  {
    return mLength;
  }

  const std::string&
  getAccession()
  {
    return mAccession;
  }

private:
  std::string mAccession;
  uint32_t mLength;
  double* mBaseProb;
  double mH1, mH0;
};

class BayesianSearcher
  : public GenBankSink
{
public:
  BayesianSearcher(std::istream& aMatrices)
    : mGBP(NewGenBankParser()), mBaseQueue(1 + kMinorLocality * 2 + kMajorLocality * 2)
  {
    static const boost::regex AcPat("^AC[ \\t]+(.*)$");
    static const boost::regex BaPat("^BA[ \\t]+[ \\t]+([0-9]+)[ \\t]+.*$");
    static const boost::regex FrPat("^([0-9]+)[ \\t]+([0-9\\.]+)[ \\t]+([0-9\\.]+)"
                                    "[ \\t]+([0-9\\.]+)[ \\t]+([0-9\\.]+)[ \\t]+.*$");

    mGBP->SetSink(this);

    std::string l;

    std::string acc;
    uint32_t basis = 0;
    ublas::matrix<double> matrix;
    uint32_t max = 0;

    mMaxLength = 0;

    while (aMatrices.good())
    {
      std::getline(aMatrices, l);

      boost::smatch res;
      if (l == "//")
      {
        // Sum the first row to see if it looks like it is normalised to a percentage...
        if (matrix.size1() > 0)
        {
          double sum1 =
            ublas::sum(ublas::matrix_row<ublas::matrix<double> >(matrix, 0));
          if (sum1 >= 96 && sum1 <= 104)
          {
            // It is probably a percentage matrix. Can we normalise it?
            if (basis != 0)
            {
              uint32_t i, j;
              for (i = 0; i < matrix.size1(); i++)
              {
                sum1 =
                  ublas::sum(ublas::matrix_row<ublas::matrix<double> >(matrix, i));
                ublas::matrix_row<ublas::matrix<double> >(matrix, i) *=
                  (basis / sum1);
              }
            }
          }

          // We now have a matrix of doubles. Put it on the list...
          motifs.push_back(new Motif(acc, matrix));
          mMaxLength = std::max(mMaxLength, matrix.size1());
        }

        basis = 0;
        acc = "";
        max = 0;
        matrix.resize(0, 4);
      }
      else if (boost::regex_match(l, res, FrPat))
      {
        uint32_t i = strtoul(res[1].str().c_str(), NULL, 10);
        double count_a = strtod(res[2].str().c_str(), NULL),
          count_c = strtod(res[3].str().c_str(), NULL),
          count_g = strtod(res[4].str().c_str(), NULL),
          count_t = strtod(res[5].str().c_str(), NULL);

        if (i == 0)
          continue;
        if (i > max)
        {
          matrix.resize(i, 4);
          max = i;
        }

        matrix(i - 1, 0) = count_a;
        matrix(i - 1, 1) = count_c;
        matrix(i - 1, 2) = count_g;
        matrix(i - 1, 3) = count_t;
      }
      else if (boost::regex_match(l, res, BaPat))
        basis = strtoul(res[1].str().c_str(), NULL, 10);
      else if (boost::regex_match(l, res, AcPat))
        acc = res[1];
    }

    mAltProbs = new double[motifs.size()];
    std::sort(motifs.begin(), motifs.end(),
              ll::bind(&Motif::getLength, ll::_1) <
              ll::bind(&Motif::getLength, ll::_2));

    std::vector<Motif*>::iterator i(motifs.begin());
    for (uint32_t x = 0; x < mMaxLength; x++)
    {
      while (i != motifs.end() && (*i)->getLength() - 1 < x)
        i++;
      mMotifLengthEnds.push_back(i);
    }
    mMotifLengthEnds.push_back(motifs.end());
  }

  ~BayesianSearcher()
  {
    std::vector<Motif*>::iterator i;
    for (i = motifs.begin(); i != motifs.end(); i++)
      delete (*i);

    delete mGBP;

    delete [] mAltProbs;
  }

  void
  search(const std::string& aFile)
  {
    TextSource* ts = NewBufferedFileSource(aFile.c_str());
    mGBP->SetSource(ts);
    try
    {
      mGBP->Parse();
    }
    catch (ParserException& pe)
    {
      std::cout << "Parse error: " << pe.what() << std::endl;
    }

    mGBP->SetSource(NULL);
    delete ts;
  }

  void
  OpenKeyword(const char* name, const char* value)
  {
    if (!::strcmp(name, "LOCUS"))
    {
      if (mBaseQueue.size() != 0)
        flushBaseQueue();

      mBaseQueue.clear();

      memset(mMajorCounts, 0, sizeof(mMajorCounts));
      memset(mMinorCounts, 0, sizeof(mMinorCounts));
      mOffsetInto = 0;
      std::cout << "LOCUS       " << value << std::endl
                << "FEATURES             Location/Qualifiers" << std::endl;
    }
  }

  void
  CloseKeyword()
  {
  }

  void
  OpenFeature(const char* name, const char* location)
  {
  }

  void
  CloseFeature()
  {
  }

  void
  Qualifier(const char* name, const char* value)
  {
  }

  void
  CodingData(const char* data)
  {
    char base;

    // std::cout << "Batch of data:" << std::endl;
    while ((base = *data++))
    {
      Base b;
      switch (base)
      {
      case 'a': case 'A':
        b = BASE_A;
        break;
      case 'c': case 'C':
        b = BASE_C;
        break;
      case 'g': case 'G':
        b = BASE_G;
        break;
      case 't': case 'T':
        b = BASE_T;
        break;
      default:
        continue;
      }

      pushBase(b);
    }
    // std::cout << std::endl;
  }
private:
  void
  pushBase(Base aBase)
  {
    uint32_t * major, * minor;

    if (mBaseQueue.size() >= kMajorLocality * 2 + kMinorLocality * 2 + 1)
    {
      mMajorCounts[BASE_N*4 + mBaseQueue.front()]--;
      mMajorCounts[mBaseQueue.front() * 4 + *(++mBaseQueue.begin())]--;
    }

    mMajorCounts[BASE_N*4 + aBase]++;
    if (mBaseQueue.size() > 0)
      mMajorCounts[mBaseQueue.back() * 4 + aBase]++;

    if (mBaseQueue.size() >= kMajorLocality + kMinorLocality * 2 + 1)
    {
      Base outgoing
        (mBaseQueue[mBaseQueue.size() - (kMajorLocality +
                                         kMinorLocality * 2 + 1)]);
      Base outgoing1
        (mBaseQueue[mBaseQueue.size() - (kMajorLocality +
                                         kMinorLocality * 2 + 1) + 1]);
      mMinorCounts[BASE_N*4 + outgoing1]--;
      mMinorCounts[outgoing * 4 + outgoing1]--;
    }

    mBaseQueue.push_back(aBase);

    if (mBaseQueue.size() == kMajorLocality + kMinorLocality + 1)
    {
      for (uint32_t i = 0; i < kMinorLocality + 1; i++)
      {
        Base incoming(mBaseQueue[i]);
        Base incoming1(mBaseQueue[i + 1]);
        mMinorCounts[BASE_N*4 + incoming1]++;
        mMinorCounts[incoming*4 + incoming1]++;
      }
    }
    else if (mBaseQueue.size() > kMajorLocality + kMinorLocality + 1)
    {
      Base incoming(mBaseQueue[mBaseQueue.size() - kMajorLocality - 1]);
      Base incoming1(mBaseQueue[mBaseQueue.size() - kMajorLocality]);

      mMinorCounts[BASE_N*4 + incoming1]++;
      mMinorCounts[incoming * 4 + incoming1]++;
    }

    if (mBaseQueue.size() >= kMajorLocality + kMinorLocality + 1)
    {
      processFrameAt(mBaseQueue.size() - (kMajorLocality + kMinorLocality + 1));
    }
  }

  void
  flushBaseQueue()
  {
    if (mBaseQueue.size() < kMajorLocality + kMinorLocality)
      return;

    uint32_t index = mBaseQueue.size() - kMajorLocality - kMinorLocality;
    uint32_t *major, *minor;

    for (; index < mBaseQueue.size(); index++)
    {
      if (index >= kMajorLocality + kMinorLocality + 2)
      {
        Base outgoing(mBaseQueue[index - kMajorLocality - kMinorLocality - 2]);
        Base outgoing1(mBaseQueue[index - kMajorLocality - kMinorLocality - 1]);
        mMajorCounts[BASE_N*4 + outgoing1]--;
        mMajorCounts[outgoing*4 + outgoing1]--;
      }

      if (index >= kMinorLocality + 2)
      {
        Base outgoing(mBaseQueue[index - kMinorLocality - 2]);
        Base outgoing1(mBaseQueue[index - kMinorLocality - 1]);
        mMinorCounts[BASE_N*4 + outgoing1]--;
        mMinorCounts[outgoing*4 + outgoing1]--;
      }

      processFrameAt(index);

      index++;
    }
  }

  void
  processFrameAt(uint32_t index)
  {
    // Firstly, we convert our 5*4 tables into frequencies...
    for (uint32_t row = BASE_A; USUALLY(row <= BASE_N); row++)
    {
      uint32_t totMaj = 0, totMin = 0;
      for (uint32_t col = BASE_A; USUALLY(col < BASE_N); col++)
      {
        totMaj += mMajorCounts[row * 4 + col];
        totMin += mMinorCounts[row * 4 + col];
      }

      for (uint32_t col = BASE_A; USUALLY(col < BASE_N); col++)
        mFreqTab[row * 4 + col] =
          (
           static_cast<double>(mMajorCounts[row * 4 + col]) /
           totMaj
           +
           static_cast<double>(mMinorCounts[row * 4 + col]) /
           totMin
          ) / 2.0;
    }

    mOffsetInto++;

#ifdef TIME_PROFILING
    if (mOffsetInto > 1000000)
    {
      pt::ptime stopTime(pt::microsec_clock::universal_time());
      pt::time_duration runTime(stopTime - gStartupTime);
      std::cout << "1 megabase in "
                << runTime.hours() << " hours, "
                << runTime.minutes() << " minutes, "
                << runTime.seconds() << " seconds, "
                << runTime.fractional_seconds() << " microseconds."
                << std::endl;
      
      exit(0);
    }
#endif

    // Next we go through each base and build the probability of the data
    // given the null hypothesis, and for each motif, the probability of the
    // alternative hypothesis.
    boost::circular_buffer<Base>::iterator bi(mBaseQueue.begin() + index);
    std::vector<std::vector<Motif*>::iterator>::iterator
      mlei(mMotifLengthEnds.begin());

    Base prev = BASE_N;
    double pBackground = 1.0;
    const double* e = mAltProbs + motifs.size();
    for (double* p = mAltProbs; USUALLY(p < e); p++)
      *p = 1.0;

    if (RARELY(mMaxLength >= mBaseQueue.size()))
      return;

    for (uint32_t l = 0;
         USUALLY(l < mMaxLength);
         l++, mlei++)
    {
      Base cur = *bi++;
      pBackground *= mFreqTab[prev * 4 + cur];
      prev = cur;

      double* ap = mAltProbs + (*mlei - motifs.begin());
      
      for (std::vector<Motif*>::iterator i = *mlei, e = motifs.end();
           USUALLY(i != e);
           i++)
        *ap++ *= (*i)->getAlternativeProbability(l, cur);

      for (std::vector<Motif*>::iterator i = *mlei, e = *(mlei + 1);
           USUALLY(i != e); i++)
      {
        double pp = (*i)->getPosterior(pBackground, mAltProbs[i - motifs.begin()]);
        if (pp > mPosteriorCutoff)
        {
          const char* t = "ACGT";
          std::cout << "     TFBS            " << mOffsetInto << ".."
                    << (mOffsetInto + (*i)->getLength()) << std::endl
                    << "                     /probability=\""
                    << pp << "\"" << std::endl
                    << "                     /db_xref=\"TRANSFAC:"
                    << (*i)->getAccession() << "\"" << std::endl;
#if 0
          std::cout << "Context: "
                    << t[mBaseQueue[index]]
                    << t[mBaseQueue[index + 1]]
                    << t[mBaseQueue[index + 2]]
                    << t[mBaseQueue[index + 3]]
                    << t[mBaseQueue[index + 4]]
                    << t[mBaseQueue[index + 5]]
                    << t[mBaseQueue[index + 6]]
                    << t[mBaseQueue[index + 7]]
                    << t[mBaseQueue[index + 8]]
                    << t[mBaseQueue[index + 9]]
                    << t[mBaseQueue[index + 10]]
                    << t[mBaseQueue[index + 11]]
                    << t[mBaseQueue[index + 12]]
                    << t[mBaseQueue[index + 13]]
                    << t[mBaseQueue[index + 14]]
                    << t[mBaseQueue[index + 15]]
                    << t[mBaseQueue[index + 16]]
                    << t[mBaseQueue[index + 17]]
                    << t[mBaseQueue[index + 18]]
                    << t[mBaseQueue[index + 19]]
                    << std::endl;
#endif
        }
      }
    }
  }

  // Look-ahead an look-behind contexts for background frequency table.
  static const size_t kMajorLocality = 1000;
  static const size_t kMinorLocality = 250;
  static const double mPosteriorCutoff = 0.5;
  uint32_t mMajorCounts[5 * 4], mMinorCounts[5 * 4], mOffsetInto;
  double mFreqTab[5 * 4];
  boost::circular_buffer<Base> mBaseQueue;
  std::vector<std::vector<Motif*>::iterator > mMotifLengthEnds;
  size_t mMaxLength;
  double *mAltProbs;

  std::vector<Motif*> motifs;
  GenBankParser* mGBP;
};

int
main(int argc, char** argv)
{
  std::string matrices, genbank;

  po::options_description desc("Allowed options");
  desc.add_options()
    ("matrices", po::value<std::string>(&matrices), "TRANSFAC format matrices file")
    ("genbank", po::value<std::string>(&genbank), "Directory containing GenBank files")
    ;

  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc), vm);
  po::notify(vm);

  if (!fs::is_directory(genbank))
  {
    std::cerr << "Supplied GenBank 'directory' is not a valid directory."
              << std::endl;
    return 1;
  }

  io::filtering_istream fsmatrices;
  fsmatrices.push(io::file_source(matrices));

  BayesianSearcher searcher(fsmatrices);

  gStartupTime = pt::microsec_clock::universal_time();

  for (fs::directory_iterator it(genbank); it != fs::directory_iterator(); it++)
  {
    if (fs::extension(it->path()) != ".gbk")
      continue;
    searcher.search(it->path().string());
  }
}
