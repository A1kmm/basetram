#define NON_MPI

/*
    Bayesian Search for Transcription factor Motifs.
    Copyright (C) 2008-2009  Andrew Miller

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU Affero General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Affero General Public License for more details.

    You should have received a copy of the GNU Affero General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/
#define BOOST_MATH_NO_LONG_DOUBLE_MATH_FUNCTIONS
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
#ifndef NON_MPI
#include <boost/mpi/environment.hpp>
#include <boost/mpi/communicator.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/utility.hpp>
#endif
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
#ifndef NON_MPI
namespace mpi = boost::mpi;
#endif

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
  Motif(const std::string& aAccession, ublas::matrix<float> aRawFreqD)
    : mAccession(aAccession), mReverseComplement(false), mNotChunked(true)
  {
    uint32_t max_n = 0;

    mLength = aRawFreqD.size1();
    mBaseProb = new float[mLength * 4];

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
          boost::math::beta<float>(f_ib + 2.0, n_i - f_ib + 1.0)
          /
          boost::math::beta<float>(f_ib + 1.0, n_i - f_ib + 1.0);
      }
    }

    // XXX should we make the prior configurable or something? For now, just
    // assume that of 4 billion possible positions, we have found 1/100th of
    // them. This is unavoidably a total guess.
    mH1 = max_n * (100.0 / 4.0E9);
    mH0 = 1.0 - mH1;
  }

  Motif(const Motif& aM, bool aComplement)
    : mAccession(aM.getAccession()), mLength(aM.getLength()),
      mH1(aM.getH1()), mH0(aM.getH0()),
      mReverseComplement(aComplement),
      mNotChunked(true)
  {
    mBaseProb = new float[mLength * 4];
    if (!aComplement)
      memcpy(mBaseProb, aM.getBaseProb(), mLength * 4 * sizeof(float));
    else
    {
      const float* oldProbs = aM.getBaseProb();

      const Base complements[] =
        {
          /* BASE_A */BASE_T,
          /* BASE_C */BASE_G,
          /* BASE_G */BASE_C,
          /* BASE_T */BASE_A
        };
      
      for (uint32_t i = 0; i < mLength; i++)
        for (uint32_t j = 0; j < 4; j++)
        {
          mBaseProb[(i << 2) + complements[j]] =
            oldProbs[((mLength - i - 1) << 2) + j];
        }
    }
  }

  ~Motif()
  {
    if (mNotChunked)
      delete [] mBaseProb;
  }

  float
  getPosterior(
               float pD_given_H0, float pD_given_H1
              ) const
  {
    return mH1 * pD_given_H1 / (mH1 * pD_given_H1 + mH0 * pD_given_H0);
  }

  float getAlternativeProbability(uint32_t position, Base b) const
  {
    return mBaseProb[(position << 2) + b];
  }

  uint32_t
  getLength() const
  {
    return mLength;
  }

  const std::string&
  getAccession() const
  {
    return mAccession;
  }

  const float* getBaseProb() const
  {
    return mBaseProb;
  }

  float getH1() const
  {
    return mH1;
  }

  float getH0() const
  {
    return mH0;
  }

  bool isReverseComplement() const
  {
    return mReverseComplement;
  }

  void useConsolidatedMatrixMemory(char** aSpace)
  {
    mNotChunked = false;
    float* oldData = mBaseProb;

    mBaseProb = reinterpret_cast<float*>(*aSpace);
    *aSpace += mLength * 4 * sizeof(float);

    memcpy(mBaseProb, oldData, sizeof(float) * 4 * mLength);
    
    delete [] oldData;
  }

private:
  std::string mAccession;
  uint32_t mLength;
  float* mBaseProb;
  float mH1, mH0;
  bool mReverseComplement, mNotChunked;
};

#ifndef NON_MPI
class WorkDoneException
{
};
#endif

class BayesianSearcher
  : public GenBankSink
{
public:
  BayesianSearcher(std::istream& aMatrices, const fs::path& aOutputDir)
    : mGBP(NewGenBankParser()), mBaseQueue(1 + kMinorLocality * 2 + kMajorLocality * 2),
      mOutputDirectory(aOutputDir)
  {
    static const boost::regex AcPat("^AC[ \\t]+(.*)$");
    static const boost::regex BaPat("^BA[ \\t]+[ \\t]+([0-9]+)[ \\t]+.*$");
    static const boost::regex FrPat("^([0-9]+)[ \\t]+([0-9\\.]+)[ \\t]+([0-9\\.]+)"
                                    "[ \\t]+([0-9\\.]+)[ \\t]+([0-9\\.]+)[ \\t]+.*$");
    static const boost::regex MxPat("^M([ACGT])[ \\t]+([0-9\\. \\t]+)$");

    mGBP->SetSink(this);

    std::string l;

    std::string acc;
    uint32_t basis = 0;
    ublas::matrix<float> matrix;
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
          float sum1 =
            ublas::sum(ublas::matrix_row<ublas::matrix<float> >(matrix, 0));
          if (sum1 >= 96 && sum1 <= 104)
          {
            // It is probably a percentage matrix. Can we normalise it?
            if (basis != 0)
            {
              uint32_t i, j;
              for (i = 0; i < matrix.size1(); i++)
              {
                sum1 =
                  ublas::sum(ublas::matrix_row<ublas::matrix<float> >(matrix, i));
                ublas::matrix_row<ublas::matrix<float> >(matrix, i) *=
                  (basis / sum1);
              }
            }
          }

          // We now have a matrix of floats. Put it on the list...
          Motif *mForward = new Motif(acc, matrix);
          motifs.push_back(mForward);
          Motif *mReverse = new Motif(*mForward, true);
          motifs.push_back(mReverse);
          mMaxLength = std::max(mMaxLength, matrix.size1());
        }

        basis = 0;
        acc = "";
        max = 0;
        matrix.resize(0, 4);
      }
      else if (boost::regex_match(l, res, MxPat))
      {
        uint32_t base;
        switch (res[1].str().c_str()[0])
        {
        case 'A':
          base = 0;
          break;
        case 'C':
          base = 1;
          break;
        case 'G':
          base = 2;
          break;
        default:
          base = 3;
        }

        typedef boost::tokenizer<boost::char_separator<char> > tok;
        boost::char_separator<char> sep(" \t", "");
        tok mytok(res[2].str(), sep);
        uint32_t j = 0;
        for (tok::iterator i = mytok.begin(); i != mytok.end(); i++, j++)
        {
          if (j > max)
          {
            matrix.resize(j, 4);
            max = j;
          }
          matrix(j, base) = strtod((*i).c_str(), NULL);
        }
      }
      else if (boost::regex_match(l, res, FrPat))
      {
        uint32_t i = strtoul(res[1].str().c_str(), NULL, 10);
        float count_a = strtod(res[2].str().c_str(), NULL),
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

    // Now we go through all motifs and try to get the matrices on the same
    // page so as to avoid cache misses later.
    size_t needCapacity = sizeof(float) * motifs.size() +
      sizeof(float*) * motifs.size();

    for (i = motifs.begin(); i != motifs.end(); i++)
      needCapacity += (*i)->getLength() * 4 * sizeof(float);

    posix_memalign(reinterpret_cast<void**>(&mMatrixSpace),
                   65536, needCapacity);

    char *mm = mMatrixSpace;
    mAltProbs = reinterpret_cast<float*>(mm);
    mm += sizeof(float) * motifs.size();

    mMatrixShortcuts = reinterpret_cast<const float**>(mm);
    reinterpret_cast<char*>(mm) += sizeof(float*) * motifs.size();

    const float** p = mMatrixShortcuts;

    for (i = motifs.begin(); i != motifs.end(); i++)
    {
      (*i)->useConsolidatedMatrixMemory(&mm);
      *p++ = (*i)->getBaseProb();
    }
  }

  ~BayesianSearcher()
  {
    std::vector<Motif*>::iterator i;
    for (i = motifs.begin(); i != motifs.end(); i++)
      delete (*i);

    delete mGBP;

    if (mSegmentOutput.is_open())
      mSegmentOutput.close();

    free(mMatrixSpace);
  }

  void
  search(const std::string& aFile, uint32_t aSeekTo = 0)
  {
    if (mSegmentOutput.is_open())
      mSegmentOutput.close();

    mChromosomeOutput = mOutputDirectory;
    mChromosomeOutput /= fs::basename(aFile);
    fs::create_directories(mChromosomeOutput);

#ifndef NON_MPI
    mSeenLocus = false;
#endif

    TextSource* ts = NewBufferedFileSource(aFile.c_str(), aSeekTo);
    mGBP->SetSource(ts);
    try
    {
      mGBP->Parse();
    }
    catch (ParserException& pe)
    {
      std::cout << "Parse error: " << pe.what() << std::endl;
    }
#ifndef NON_MPI
    catch (WorkDoneException& wde)
    {
    }
#endif

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

      if (mSegmentOutput.is_open())
        mSegmentOutput.close();

#ifndef NON_MPI
      if (mSeenLocus)
        throw WorkDoneException();
      mSeenLocus = true;
#endif

      fs::path mSegmentFile = mChromosomeOutput;
      std::string locus(value);
      size_t pos = locus.find(" ");
      mSegmentFile /= locus.substr(0, pos);
      mSegmentOutput.open(mSegmentFile.string().c_str());

      mOffsetInto = 0;
      mSegmentOutput << "LOCUS       " << value << std::endl
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
           static_cast<float>(mMajorCounts[row * 4 + col]) /
           totMaj
           +
           static_cast<float>(mMinorCounts[row * 4 + col]) /
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
    float pBackground = 1.0;
    const float* e = mAltProbs + motifs.size();
    for (float* p = mAltProbs; USUALLY(p < e); p++)
      *p = 1.0;

    if (RARELY(mMaxLength >= (mBaseQueue.size() - index)))
      return;

    for (uint32_t l = 0;
         USUALLY(l < mMaxLength);
         l++, mlei++)
    {
      Base cur = *bi++;
      pBackground *= mFreqTab[prev * 4 + cur];
      prev = cur;

      float* ap = mAltProbs + (*mlei - motifs.begin());
      const float** ms = mMatrixShortcuts + (*mlei - motifs.begin()),
           ** mse = mMatrixShortcuts + motifs.size();

      // This is the most time consuming part of the program, hence why it has
      // so many ugly optimisations to help ensure everything can be done with
      // the registers + the L1 cache.
      for (; USUALLY(ms != mse); ms++)
        *ap++ *= (*ms)[(l << 2) + cur];

      for (std::vector<Motif*>::iterator i = *mlei, e = *(mlei + 1);
           USUALLY(i != e); i++)
      {
        float pp = (*i)->getPosterior(pBackground, mAltProbs[i - motifs.begin()]);
        if (RARELY(pp > mPosteriorCutoff))
        {
          const char* t = "ACGT";

          if (!(*i)->isReverseComplement())
            mSegmentOutput << "     TFBS            " << mOffsetInto << ".."
                           << (mOffsetInto + (*i)->getLength()) << std::endl;
          else
            mSegmentOutput << "     TFBS            complement(" << mOffsetInto << ".."
                           << (mOffsetInto + (*i)->getLength()) << ")" << std::endl;

          mSegmentOutput << "                     /probability=\""
                         << pp << "\"" << std::endl
                         << "                     /db_xref=\"TRANSFAC:"
                         << (*i)->getAccession() << "\"" << std::endl;
        }
      }
    }
  }

  // Look-ahead an look-behind contexts for background frequency table.
  static const size_t kMajorLocality = 1000;
  static const size_t kMinorLocality = 250;
  static const float mPosteriorCutoff = 0.5;
  uint32_t mMajorCounts[5 * 4], mMinorCounts[5 * 4], mOffsetInto;
  float mFreqTab[5 * 4];
  boost::circular_buffer<Base> mBaseQueue;
  std::vector<std::vector<Motif*>::iterator > mMotifLengthEnds;
  size_t mMaxLength;
  float *mAltProbs;

  const fs::path& mOutputDirectory;
  fs::path mChromosomeOutput;
  std::ofstream mSegmentOutput;
  char* mMatrixSpace;
  const float** mMatrixShortcuts;

  std::vector<Motif*> motifs;
  GenBankParser* mGBP;

#ifndef NON_MPI
  bool mSeenLocus;
#endif
};

class LocusOffsetFinder
  : public GenBankSink
{
public:
  LocusOffsetFinder()
    : mGBP(NewGenBankParser())
  {
    mGBP->SetSink(this);
  }

  std::vector<uint32_t>&
  search(const std::string& aFile)
  {
    mOffsets.clear();
    mTS = NewBufferedFileSource(aFile.c_str());
    mGBP->SetSource(mTS);
    try
    {
      mGBP->Parse();
    }
    catch (ParserException& pe)
    {
      std::cout << "Parse error: " << pe.what() << std::endl;
    }

    mGBP->SetSource(NULL);
    delete mTS;

    return mOffsets;
  }

  void
  OpenKeyword(const char* name, const char* value)
  {
    if (!::strcmp(name, "LOCUS"))
    {
      size_t offset
        ((mTS->getOffset() - mGBP->CountBufferedBytes())
         - 14 - strlen(value));

      mOffsets.push_back(offset);
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
  }

private:
  GenBankParser* mGBP;
  TextSource* mTS;
  std::vector<uint32_t> mOffsets;
};

int
main(int argc, char** argv)
{
#ifndef NON_MPI
  mpi::environment env(argc, argv);
  mpi::communicator world;
#endif

  std::string matrices, genbank, outdir;

  po::options_description desc("Allowed options");
  desc.add_options()
    ("matrices", po::value<std::string>(&matrices), "TRANSFAC format matrices file")
    ("genbank", po::value<std::string>(&genbank), "Directory containing GenBank files")
    ("outdir", po::value<std::string>(&outdir), "Top-level directory for results")
    ("help", "produce help message")
    ;

  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc), vm);
  po::notify(vm);

  std::string wrong;
  if (!vm.count("help"))
  {
    if (!vm.count("matrices"))
      wrong = "matrices";
    else if (!vm.count("genbank"))
      wrong = "genbank";
    else if (!vm.count("outdir"))
      wrong = "outdir";
  }

  if (wrong != "")
    std::cerr << "Missing option: " << wrong << std::endl;
  if (vm.count("help") || wrong != "")
  {
    std::cout << desc << std::endl;
    return 1;
  }

  if (!fs::is_directory(genbank))
  {
    std::cerr << "Supplied GenBank 'directory' is not a valid directory."
              << std::endl;
    return 1;
  }

  if (!fs::is_directory(outdir))
  {
    std::cerr << "Supplied output 'directory' is not a valid directory."
              << std::endl;
    return 1;
  }

#ifdef NON_MPI
  io::filtering_istream fsmatrices;
  fsmatrices.push(io::file_source(matrices));

  fs::path outdirp(outdir);

  BayesianSearcher searcher(fsmatrices, outdirp);

  gStartupTime = pt::microsec_clock::universal_time();

  for (fs::directory_iterator it(genbank); it != fs::directory_iterator(); it++)
  {
    if (fs::extension(it->path()) != ".gbk")
      continue;
    searcher.search(it->path().string());
  }
#else
  if (world.size() < 2)
  {
    std::cerr << "Sorry, you must have at least two MPI processes - one for "
              << "control and one or more worker processes."
              << std::endl;
    return 1;
  }

  if (world.rank() == 0)
  {
    uint32_t nWorkers = world.size() - 1;

    // Make a list of chromosomes and assign workers to them...
    uint32_t index = 0;
    
    uint32_t expectedResponses = 0;

    for (fs::directory_iterator it(genbank); it != fs::directory_iterator();
         it++)
    {
      if (fs::extension(it->path()) != ".gbk")
        continue;
      uint32_t worker = (index++ % nWorkers) + 1;

      expectedResponses++;

      std::cout << "Starting worker " << worker << " on scan of "
                << it->path().string() << std::endl;
      world.send(worker, 0, true);
      world.send(worker, 0, it->path().string());
    }

    // Let all workers know the scan phase is over and they should prepare for
    // messages to start searching contigs...
    for (uint32_t worker = 1; worker <= nWorkers; worker++)
      world.send(worker, 0, false);

    std::list<std::pair<std::string, uint32_t> > contigs;
    while (expectedResponses)
    {
      std::pair<std::string, std::vector<uint32_t> > data;

      while (!world.iprobe(mpi::any_source, 0))
        sleep(1);

      mpi::status s(world.recv(mpi::any_source, 0, data));

      std::cout << "Controller got response from " << s.source() << std::endl;
      expectedResponses--;

      for (std::vector<uint32_t>::iterator di = data.second.begin();
           di != data.second.end();
           di++)
        contigs.push_back(std::pair<std::string, uint32_t>
                          (data.first, *di));
    }

    std::cout << "Preliminary scan of contigs complete."
              << std::endl;

    std::list<uint32_t> availableWorkers;
    for (uint32_t i = 1; i <= nWorkers; i++)
      availableWorkers.push_back(i);

    while (!contigs.empty() || (availableWorkers.size() != nWorkers))
    {
      while (!availableWorkers.empty() &&
             !contigs.empty())
      {
        uint32_t worker = availableWorkers.back();
        availableWorkers.pop_back();
        world.send(worker, 1, true);
        std::cout << "Starting worker " << worker
                  << " on search of " << contigs.front().first
                  << " contig at " << contigs.front().second
                  << std::endl;
        world.send(worker, 1, contigs.front());
        contigs.pop_front();
      }

      while (!world.iprobe(mpi::any_source, 1))
        sleep(1);

      mpi::status s(world.recv(mpi::any_source, 1));
      uint32_t worker(s.source());

      std::cout << "Worker " << worker << " is now free."
                << std::endl;
      availableWorkers.push_back(worker);
    }

    for (uint32_t worker = 1; worker <= nWorkers; worker++)
      world.send(worker, 1, false);

    return 0;
  }

  // If we get here, we are a worker process. Carry out any scan commands we
  // get...
  {
    LocusOffsetFinder lof;
    while (true)
    {
      bool stillScanning;
      world.recv(0, 0, stillScanning);
      
      if (!stillScanning)
        break;

      std::string chromosome;
      world.recv(0, 0, chromosome);

      std::cout << "[" << world.rank() << "]: Received instruction to scan "
                << chromosome << std::endl;

      std::vector<uint32_t>& r(lof.search(chromosome));

      world.send(0, 0, std::pair<std::string, std::vector<uint32_t> >
                 (chromosome, r));
    }
  }

  {
    io::filtering_istream fsmatrices;
    fsmatrices.push(io::file_source(matrices));
    
    fs::path outdirp(outdir);
    
    BayesianSearcher searcher(fsmatrices, outdirp);

    while (true)
    {
      bool stillSearching;
      world.recv(0, 1, stillSearching);
      
      if (!stillSearching)
        break;

      std::pair<std::string, uint32_t> contig;
      world.recv(0, 1, contig);

      std::cout << "[" << world.rank() << "]: Received instruction to search "
                << contig.first << " offset "
		<< contig.second << std::endl;

      searcher.search(contig.first, contig.second);

      // and tell the controller we finished...
      world.send(0, 1);
    }
  }
#endif
}
