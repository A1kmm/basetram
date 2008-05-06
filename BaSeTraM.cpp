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
#include <cmath>
#include <assert.h>
#include "../parsegenbank/GenbankParser.hpp"

namespace po = boost::program_options;
namespace io = boost::iostreams;
namespace fs = boost::filesystem;
namespace ublas = boost::numeric::ublas;

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
    // assume that of 4 billion possible positions, we have found 1/10th of
    // them. This is unavoidably a total guess.
    mH1 = max_n * (10.0 / 4.0E9);
    mH0 = 1.0 - mH1;
  }

  ~Motif()
  {
    delete [] mBaseProb;
  }

  double
  getPosterior(
               uint32_t aIndex,
               boost::circular_buffer<char>& aBaseQueue,
               double aPA, double aPG, double aPC, double aPT
              )
  {
    // Impossible if we don't have enough data...
    if (aIndex + mLength >= aBaseQueue.size())
      return 0.0;

    double pD_given_H1 = 1.0, pD_given_H0 = 1.0;

    double* bb = mBaseProb;

    boost::circular_buffer<char>::iterator i = aBaseQueue.begin() + aIndex,
      e = i + mLength;

    for (; i != e; i++, bb += 4)
    {
      char base = *i;
      switch (base)
      {
      case 'A': case 'a':
        pD_given_H0 *= aPA;
        pD_given_H1 *= bb[0];
        continue;
      case 'C': case 'c':
        pD_given_H0 *= aPC;
        pD_given_H1 *= bb[1];
        continue;
      case 'G': case 'g':
        pD_given_H0 *= aPG;
        pD_given_H1 *= bb[2];
        continue;
      case 'T': case 't':
        pD_given_H0 *= aPT;
        pD_given_H1 *= bb[3];
        continue;
      }
    }

    return mH1 * pD_given_H1 / (mH1 * pD_given_H1 + mH0 * pD_given_H0);
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
    
  }

  ~BayesianSearcher()
  {
    std::list<Motif*>::iterator i;
    for (i = motifs.begin(); i != motifs.end(); i++)
      delete (*i);

    delete mGBP;
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
      mAMajor = mAMinor = mGMajor = mGMinor = mCMajor = mCMinor = mTMajor =
        mTMinor = mOffsetInto = 0;
      std::cout << "Locus = " << value << std::endl;
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
      if (base != 'a' && base != 'A' && base != 'g' && base != 'G' &&
          base != 'c' && base != 'C' && base != 't' && base != 'T')
        continue;

      // std::cout << base;

      pushBase(base);
    }
    // std::cout << std::endl;
  }
private:
  void
  getMajorCounts(char aBase, uint32_t** aMajor)
  {
    switch (aBase)
    {
    case 'a': case 'A':
      *aMajor = &mAMajor;
      return;
    case 'g': case 'G':
      *aMajor = &mGMajor;
      return;
    case 'c': case 'C':
      *aMajor = &mCMajor;
      return;
    case 't': case 'T':
      *aMajor = &mTMajor;
      return;
    }
    assert(0); // aBase precondition violated.
  }

  void
  getMinorCounts(char aBase, uint32_t** aMinor)
  {
    switch (aBase)
    {
    case 'a': case 'A':
      *aMinor = &mAMinor;
      return;
    case 'g': case 'G':
      *aMinor = &mGMinor;
      return;
    case 'c': case 'C':
      *aMinor = &mCMinor;
      return;
    case 't': case 'T':
      *aMinor = &mTMinor;
      return;
    }
    assert(0); // aBase precondition violated.
  }

  void
  pushBase(char aBase)
  {
    uint32_t * major, * minor;

    if (mBaseQueue.size() >= kMajorLocality * 2 + kMinorLocality * 2 + 1)
    {
      getMajorCounts(mBaseQueue[0], &major);
      (*major)--;
    }

    getMajorCounts(aBase, &major);
    (*major)++;

    if (mBaseQueue.size() >= kMajorLocality + kMinorLocality * 2 + 1)
    {
      char outgoing
        (mBaseQueue[mBaseQueue.size() - (kMajorLocality +
                                         kMinorLocality * 2 + 1)]);
      getMinorCounts(outgoing, &minor);
      (*minor)--;
    }

    mBaseQueue.push_back(aBase);

    if (mBaseQueue.size() == kMajorLocality + kMinorLocality + 1)
    {
      for (uint32_t i = 0; i < kMinorLocality + 1; i++)
      {
        getMinorCounts(mBaseQueue[i], &minor);
        (*minor)++;
      }
    }
    else if (mBaseQueue.size() > kMajorLocality + kMinorLocality + 1)
    {
      getMinorCounts(mBaseQueue[mBaseQueue.size() - kMajorLocality - 1], &minor);
      (*minor)++;
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
        getMajorCounts(mBaseQueue[index - kMajorLocality - kMinorLocality - 2], &major);
        (*major)--;
      }

      if (index >= kMinorLocality + 2)
      {
        getMinorCounts(mBaseQueue[index - kMinorLocality - 2], &minor);
        (*minor)--;
      }

      processFrameAt(index);

      index++;
    }
  }

  void
  processFrameAt(uint32_t index)
  {
    // Now compute the local base frequencies as the average of the major and
    // minor moving averages...
    uint32_t majorCount(mBaseQueue.size());
    uint32_t minorCount(std::min(mBaseQueue.size(), kMinorLocality * 2 + 1));
    double imajorCount2(1.0 / (2.0 * majorCount));
    double iminorCount2(1.0 / (2.0 * minorCount));
    
    double freqA(imajorCount2 * mAMajor + iminorCount2 * mAMinor);
    double freqG(imajorCount2 * mGMajor + iminorCount2 * mGMinor);
    double freqC(imajorCount2 * mCMajor + iminorCount2 * mCMinor);
    double freqT(imajorCount2 * mTMajor + iminorCount2 * mTMinor);

    mOffsetInto++;

    // Now we go through all the motifs and look for matches...
    std::list<Motif*>::iterator i;
    for (i = motifs.begin(); i != motifs.end(); i++)
    {
      double pp;
      pp = (*i)->getPosterior(index, mBaseQueue, freqA, freqG, freqC, freqT);
      if (pp > mPosteriorCutoff)
      {
        std::cout << "At offset " << mOffsetInto << ": Match "
                  << (*i)->getAccession() << " with probability "
                  << pp << std::endl;
        std::cout << "Context: "
                  << mBaseQueue[index]
                  << mBaseQueue[index + 1]
                  << mBaseQueue[index + 2]
                  << mBaseQueue[index + 3]
                  << mBaseQueue[index + 4]
                  << mBaseQueue[index + 5]
                  << mBaseQueue[index + 6]
                  << mBaseQueue[index + 7]
                  << mBaseQueue[index + 8]
                  << mBaseQueue[index + 9]
                  << mBaseQueue[index + 10]
                  << mBaseQueue[index + 11]
                  << mBaseQueue[index + 12]
                  << mBaseQueue[index + 13]
                  << mBaseQueue[index + 14]
                  << mBaseQueue[index + 15]
                  << mBaseQueue[index + 16]
                  << mBaseQueue[index + 17]
                  << mBaseQueue[index + 18]
                  << mBaseQueue[index + 19]
                  << std::endl;

      }
    }

    // std::cout << freqA << "," << freqG << "," << freqC << "," << freqT
    // << std::endl;
  }

  // Look-ahead an look-behind contexts for background frequency table.
  static const unsigned int kMajorLocality = 500;
  static const unsigned int kMinorLocality = 50;
  static const double mPosteriorCutoff = 0.5;
  uint32_t mAMajor, mAMinor, mGMajor, mGMinor, mCMajor, mCMinor,
           mTMajor, mTMinor, mOffsetInto;
  boost::circular_buffer<char> mBaseQueue;

  std::list<Motif*> motifs;
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

  for (fs::directory_iterator it(genbank); it != fs::directory_iterator(); it++)
  {
    if (fs::extension(it->path()) != ".gbk")
      continue;
    searcher.search(it->path().string());
  }
}
