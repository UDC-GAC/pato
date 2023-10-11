// ==========================================================================
//                                triplexator
// ==========================================================================
// Copyright (c) 2011,2012, Fabian Buske, UQ
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of Fabian Buske or the University of Queensland nor
//       the names of its contributors may be used to endorse or promote
//       products derived from this software without specific prior written
//       permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL KNUT REINERT OR THE FU BERLIN BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
// OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.
//
// ==========================================================================
// Author: Fabian Buske <fbuske@uq.edu.au>
// ==========================================================================

#ifndef TRIPLEX_PATTERN_HPP
#define TRIPLEX_PATTERN_HPP

#include <sstream>

#include <seqan/basic.h>
#include <seqan/modifier.h>
#include <seqan/sequence.h>

#include "triplex_functors.hpp"

namespace seqan {

template <typename THost_, typename TString> class ModStringTriplex {
public:
  // clang-format off
  typedef typename Value<TString>::Type TAlphabet;
  typedef typename Host<ModStringTriplex>::Type THost;
  typedef typename Infix<THost>::Type TSegment;
  typedef unsigned TId;
  typedef Pair<TId, typename Position<ModStringTriplex<THost_, TString>>::Type> TDuplicate;
  typedef ::std::vector<TDuplicate> TDuplicates;

  typedef ModifiedString<TSegment, ModView<FunctorTCMotif>> TtcMotif;
  typedef ModifiedString<TSegment, ModView<FunctorGTMotif>> TgtMotif;
  typedef ModifiedString<ModifiedString<TSegment, ModView<FunctorGTMotif>>, ModReverse> TgtMotifRev;
  typedef ModifiedString<ModifiedString<TSegment, ModView<FunctorGAMotif>>, ModReverse> TgaMotif;
  typedef ModifiedString<TSegment, ModView<FunctorTTSMotif>> TttsMotif;
  typedef ModifiedString<ModifiedString<TSegment, ModView<FunctorTTSMotifCompl>>, ModReverse> TttsMotifRevComp;
  // clang-format on

  TString mask_string; // a string which masks the host segment
  bool parallel;       // whether the motif string is parallel or antiparallel
                       // oriented to the host
  TSegment segment;    // the host segment masked by the string
  TId seqNo;           // the sequence id
  int copies;          // number of copies
  bool isTFO; // indicates whether this is a TFO (=true) or a TTS (=false)
  char motif; // the motif of the triplex
  unsigned int score;
  TDuplicates duplicates;

  void _updateMaskString() {
    if (!isTFO) {
      if (motif == '+')
        mask_string = TttsMotif(segment);
      else
        mask_string = TttsMotifRevComp(segment);
    } else {
      if (parallel) {
        if (motif == 'M')
          mask_string = TgtMotif(segment);
        else if (motif == 'Y')
          mask_string = TtcMotif(segment);
      } else {
        if (motif == 'M')
          mask_string = TgtMotifRev(segment);
        else if (motif == 'R')
          mask_string = TgaMotif(segment);
      }
    }
  }

public:
  ModStringTriplex(typename Parameter_<THost>::Type _host,
                   bool _parallel_orientation, unsigned _seqNo, bool _isTFO,
                   char _motif)
      : parallel(_parallel_orientation), segment(_host), seqNo(_seqNo),
        copies(-1), isTFO(_isTFO), motif(_motif) {
    _updateMaskString();
  }

  ModStringTriplex(TSegment _segment, bool _parallel_orientation,
                   unsigned _seqNo, bool _isTFO, char _motif)
      : parallel(_parallel_orientation), segment(_segment), seqNo(_seqNo),
        copies(-1), isTFO(_isTFO), motif(_motif) {
    _updateMaskString();
  }

  ModStringTriplex(typename Parameter_<THost>::Type _host,
                   typename Position<THost>::Type _begin_index,
                   typename Position<THost>::Type _end_index,
                   bool _parallel_orientation, unsigned _seqNo, bool _isTFO,
                   char _motif)
      : parallel(_parallel_orientation),
        segment(_host, _begin_index, _end_index), seqNo(_seqNo), copies(-1),
        isTFO(_isTFO), motif(_motif) {
    _updateMaskString();
  }

  ModStringTriplex(typename Parameter_<THost>::Type _host,
                   typename Position<THost>::Type _begin_index,
                   typename Position<THost>::Type _end_index,
                   bool _parallel_orientation, unsigned _seqNo, bool _isTFO,
                   char _motif, unsigned _copies)
      : parallel(_parallel_orientation),
        segment(_host, _begin_index, _end_index), seqNo(_seqNo),
        copies(_copies), isTFO(_isTFO), motif(_motif) {
    _updateMaskString();
  }

  ModStringTriplex(typename Parameter_<THost>::Type _host,
                   typename Iterator<THost, Standard>::Type _begin,
                   typename Iterator<THost, Standard>::Type _end,
                   bool _parallel_orientation, unsigned _seqNo, bool _isTFO,
                   char _motif, int _copies)
      : parallel(_parallel_orientation), segment(_host, _begin, _end),
        seqNo(_seqNo), copies(_copies), isTFO(_isTFO), motif(_motif) {
    _updateMaskString();
  }

  template <typename TSource>
  inline ModStringTriplex &operator=(TSource const &source) {
    assign(*this, source);
    return *this;
  }

  template <typename TPos>
  inline typename Reference<ModStringTriplex>::Type operator[](TPos pos) {
    return value(*this, pos);
  }

  template <typename TPos>
  inline typename Reference<ModStringTriplex const>::Type
  operator[](TPos pos) const {
    return value(*this, pos);
  }
};

template <typename THost, typename TString, typename TPos>
inline typename Reference<ModStringTriplex<THost, TString>>::Type
value(ModStringTriplex<THost, TString> &me, TPos pos) {
  return *(begin(me.mask_string, Standard()) + pos);
}

template <typename THost, typename TString, typename TPos>
inline typename Reference<ModStringTriplex<THost, TString> const>::Type
value(ModStringTriplex<THost, TString> const &me, TPos pos) {
  return *(begin(me.mask_string, Standard()) + pos);
}

template <typename THost_, typename TString>
inline typename Parameter_<THost_>::Type
host(ModStringTriplex<THost_, TString> &me) {
  return host(me.segment);
}

template <typename THost_, typename TString>
inline typename Parameter_<THost_>::Type
host(ModStringTriplex<THost_, TString> const &me) {
  return host(me.segment);
}

template <typename TStream, typename THost, typename TString>
inline TStream &operator<<(TStream &target,
                           ModStringTriplex<THost, TString> const &source) {
  write(target, ttsString(source));
  return target;
}

template <typename THost_, typename TString>
inline typename Iterator<ModStringTriplex<THost_, TString>, Standard>::Type
begin(ModStringTriplex<THost_, TString> &me, Standard) {
  return begin(me.mask_string);
}

template <typename THost_, typename TString>
inline
    typename Iterator<ModStringTriplex<THost_, TString> const, Standard>::Type
    begin(ModStringTriplex<THost_, TString> const &me, Standard) {
  return begin(me.mask_string);
}

template <typename THost_, typename TString>
bool merge(ModStringTriplex<THost_, TString> &m1,
           ModStringTriplex<THost_, TString> &m2) {
  if (m1.seqNo != m2.seqNo || m1.motif != m2.motif ||
      m1.parallel != m2.parallel)
    return false;
  else {
    setBeginPosition(m1, std::min(beginPosition(m1), beginPosition(m2)));
    setEndPosition(m1, std::max(endPosition(m1), endPosition(m2)));
    m1._updateMaskString();
    return true;
  }
}

template <typename THost_, typename TString>
inline typename Position<ModStringTriplex<THost_, TString>>::Type
beginPosition(ModStringTriplex<THost_, TString> &me) {
  return beginPosition(me.segment);
}
template <typename THost_, typename TString>
inline typename Position<ModStringTriplex<THost_, TString> const>::Type
beginPosition(ModStringTriplex<THost_, TString> const &me) {
  return beginPosition(me.segment);
}

template <typename THost_, typename TString>
inline bool isParallel(ModStringTriplex<THost_, TString> &me) {
  return me.parallel;
}

template <typename THost_, typename TString>
inline bool isParallel(ModStringTriplex<THost_, TString> const &me) {
  return me.parallel;
}

template <typename THost_, typename TString, typename TId, typename TPos>
inline void addDuplicate(ModStringTriplex<THost_, TString> &me,
                         TId const &seqnr, TPos const &pos) {
  typename ModStringTriplex<THost_, TString>::TDuplicate d(seqnr, pos);
  appendValue(me.duplicates, d);
}

template <typename THost_, typename TString>
inline typename ModStringTriplex<THost_, TString>::TDuplicates
getDuplicates(ModStringTriplex<THost_, TString> &me) {
  return me.duplicates;
}

template <typename THost_, typename TString>
inline typename ModStringTriplex<THost_, TString>::TDuplicates
getDuplicates(ModStringTriplex<THost_, TString> const &me) {
  return me.duplicates;
}

template <typename THost_, typename TString, typename TIndex>
inline typename ModStringTriplex<THost_, TString>::TDuplicate
getDuplicateAt(ModStringTriplex<THost_, TString> &me, TIndex const &pos) {
  return value(me.duplicates, pos);
}

template <typename THost_, typename TString, typename TIndex>
inline typename ModStringTriplex<THost_, TString>::TDuplicate
getDuplicateAt(ModStringTriplex<THost_, TString> const &me, TIndex const &pos) {
  return value(me.duplicates, pos);
}

template <typename THost_, typename TString>
inline bool isTFO(ModStringTriplex<THost_, TString> &me) {
  return me.isTFO;
}

template <typename THost_, typename TString>
inline bool isTFO(ModStringTriplex<THost_, TString> const &me) {
  return me.isTFO;
}

template <typename THost_, typename TString>
inline char getMotif(ModStringTriplex<THost_, TString> &me) {
  return me.motif;
}

template <typename THost_, typename TString>
inline char getMotif(ModStringTriplex<THost_, TString> const &me) {
  return me.motif;
}

template <typename THost_, typename TString>
inline typename Id<ModStringTriplex<THost_, TString>>::Type
getSequenceNo(ModStringTriplex<THost_, TString> &me) {
  return me.seqNo;
}

template <typename THost_, typename TString>
inline typename Id<ModStringTriplex<THost_, TString>>::Type
getSequenceNo(ModStringTriplex<THost_, TString> const &me) {
  return me.seqNo;
}

template <typename THost_, typename TString>
inline void
setSequenceNo(ModStringTriplex<THost_, TString> &me,
              typename Id<ModStringTriplex<THost_, TString>>::Type new_seqno) {
  me.seqNo = new_seqno;
}

template <typename THost_, typename TString, typename TIterator>
inline void setBegin(ModStringTriplex<THost_, TString> &me,
                     TIterator new_begin) {
  me.data_begin_position = new_begin - begin(host(me)); //, Standard());
}

template <typename THost_, typename TString_, typename TPosition_>
inline void setBeginPosition(ModStringTriplex<THost_, TString_> &me,
                             TPosition_ new_begin) {
  setBeginPosition(me.segment, new_begin);
}

template <typename THost_, typename TString_>
inline typename Iterator<ModStringTriplex<THost_, TString_>, Standard>::Type
end(ModStringTriplex<THost_, TString_> &me, Standard) {
  return end(me.mask_string);
}

template <typename THost_, typename TString_>
inline
    typename Iterator<ModStringTriplex<THost_, TString_> const, Standard>::Type
    end(ModStringTriplex<THost_, TString_> const &me, Standard) {
  return end(me.mask_string);
}

template <typename THost_, typename TString>
inline typename Position<ModStringTriplex<THost_, TString>>::Type
endPosition(ModStringTriplex<THost_, TString> &me) {
  return endPosition(me.segment);
}

template <typename THost_, typename TString>
inline typename Position<ModStringTriplex<THost_, TString>>::Type
endPosition(ModStringTriplex<THost_, TString> const &me) {
  return endPosition(me.segment);
}

template <typename THost_, typename TString>
inline unsigned int score(ModStringTriplex<THost_, TString> &me) {
  return me.score;
}

template <typename THost_, typename TString>
inline unsigned int score(ModStringTriplex<THost_, TString> const &me) {
  return me.score;
}

template <typename THost_, typename TString_, typename TScore>
inline void setScore(ModStringTriplex<THost_, TString_> &me, TScore score) {
  me.score = score;
}

template <typename THost, typename TString>
inline typename Infix<THost>::Type
getSegment(ModStringTriplex<THost, TString> &me) {
  return me.segment;
}

template <typename THost, typename TString>
inline typename Infix<THost>::Type const
getSegment(ModStringTriplex<THost, TString> const &me) {
  return me.segment;
}

template <typename THost, typename TString>
inline double guanineRate(ModStringTriplex<THost, TString> &me) {
  double guanines = 0.;
  for (unsigned i = 0; i < length(me); ++i) {
    if (value(me.mask_string, i) == 'G' || value(me.mask_string, i) == 'g')
      ++guanines;
  }
  return guanines / length(me);
}

template <typename THost, typename TString>
inline double guanineRate(ModStringTriplex<THost, TString> const &me) {
  double guanines = 0.;
  for (unsigned i = 0; i < length(me); ++i) {
    if (value(me.mask_string, i) == 'G' || value(me.mask_string, i) == 'g')
      ++guanines;
  }
  return guanines / length(me);
}

template <typename THost, typename TString>
inline CharString errorString(ModStringTriplex<THost, TString> &me) {
  std::ostringstream errors;

  CharString ps = prettyString(me);
  if (me.isTFO) {
    for (unsigned i = 0; i < length(ps); ++i) {
      if (!isupper(value(ps, i)))
        errors << 'o' << i;
    }
  } else {
    if (me.motif == '+') {
      for (unsigned i = 0; i < length(ps); ++i) {
        if (!isupper(value(ps, i)))
          errors << 'd' << i;
      }
    } else {
      for (unsigned i = length(ps); i > 0; --i) {
        if (!isupper(value(ps, i - 1)))
          errors << 'd' << (length(ps) - i);
      }
    }
  }
  return errors.str();
}

template <typename THost, typename TString>
inline CharString const
errorString(ModStringTriplex<THost, TString> const &me) {
  std::ostringstream errors;
  CharString ps = prettyString(me);
  if (me.isTFO) {
    for (unsigned i = 0; i < length(ps); ++i) {
      if (!isupper(value(ps, i)))
        errors << 'o' << i;
    }
  } else {
    if (me.motif == '+') {
      for (unsigned i = 0; i < length(ps); ++i) {
        if (!isupper(value(ps, i)))
          errors << 'd' << i;
      }
    } else {
      for (unsigned i = length(ps); i > 0; --i) {
        if (!isupper(value(ps, i - 1)))
          errors << 'd' << (length(ps) - i);
      }
    }
  }
  return errors.str();
}

template <typename THost, typename TString>
inline TString ttsString(ModStringTriplex<THost, TString> &me) {
  return me.mask_string;
}

template <typename THost, typename TString>
inline TString const ttsString(ModStringTriplex<THost, TString> const &me) {
  return me.mask_string;
}

template <typename THost, typename TString>
inline THost tfoString(ModStringTriplex<THost, TString> &me) {
  THost tempSeq2(me.segment);
  return tempSeq2;
}

template <typename THost, typename TString>
inline THost const tfoString(ModStringTriplex<THost, TString> const &me) {
  THost tempSeq2(me.segment);
  return tempSeq2;
}

template <typename THost, typename TString>
inline CharString outputString(ModStringTriplex<THost, TString> &me) {
  typedef typename Infix<THost>::Type TSegment;
  typedef ModifiedString<TSegment, ModView<FunctorTTSMotifOutput>>
      TttsMotifOutput;
  typedef ModifiedString<
      ModifiedString<TSegment, ModView<FunctorTTSMotifComplOutput>>, ModReverse>
      TttsMotifRevCompOutput;
  typedef ModifiedString<TSegment, ModView<FunctorTCMotifOutput>>
      TtcMotifOutput;
  typedef ModifiedString<TSegment, ModView<FunctorGTMotifOutput>>
      TgtMotifOutput;
  typedef ModifiedString<TSegment, ModView<FunctorGAMotifOutput>>
      TgaMotifOutput;

  if (!me.isTFO) {
    if (me.motif == '+') {
      return TttsMotifOutput(me.segment);
    } else {
      return TttsMotifRevCompOutput(me.segment);
    }
  } else {
    if (me.motif == 'M')
      return TgtMotifOutput(me.segment);
    else if (me.motif == 'Y')
      return TtcMotifOutput(me.segment);
    else if (me.motif == 'R')
      return TgaMotifOutput(me.segment);
    else
      return "";
  }
}

template <typename THost, typename TString>
inline CharString const
outputString(ModStringTriplex<THost, TString> const &me) {
  typedef typename Infix<THost>::Type TSegment;
  typedef ModifiedString<TSegment, ModView<FunctorTTSMotifOutput>>
      TttsMotifOutput;
  typedef ModifiedString<
      ModifiedString<TSegment, ModView<FunctorTTSMotifComplOutput>>, ModReverse>
      TttsMotifRevCompOutput;
  typedef ModifiedString<TSegment, ModView<FunctorTCMotifOutput>>
      TtcMotifOutput;
  typedef ModifiedString<TSegment, ModView<FunctorGTMotifOutput>>
      TgtMotifOutput;
  typedef ModifiedString<TSegment, ModView<FunctorGAMotifOutput>>
      TgaMotifOutput;

  if (!me.isTFO) {
    if (me.motif == '+') {
      return TttsMotifOutput(me.segment);
    } else {
      return TttsMotifRevCompOutput(me.segment);
    }
  } else {
    if (me.motif == 'M')
      return TgtMotifOutput(me.segment);
    else if (me.motif == 'Y')
      return TtcMotifOutput(me.segment);
    else if (me.motif == 'R')
      return TgaMotifOutput(me.segment);
    else
      return "";
  }
}

template <typename THost, typename TString>
inline CharString prettyString(ModStringTriplex<THost, TString> &me) {
  typedef typename Infix<THost>::Type TSegment;
  typedef ModifiedString<TSegment, ModView<FunctorTTSMotifPretty>>
      TttsMotifOutput;
  typedef ModifiedString<
      ModifiedString<TSegment, ModView<FunctorTTSMotifComplPretty>>, ModReverse>
      TttsMotifRevCompOutput;
  typedef ModifiedString<TSegment, ModView<FunctorTCMotifPretty>>
      TtcMotifOutput;
  typedef ModifiedString<TSegment, ModView<FunctorGTMotifPretty>>
      TgtMotifOutput;
  typedef ModifiedString<TSegment, ModView<FunctorGAMotifPretty>>
      TgaMotifOutput;

  if (!me.isTFO) {
    if (me.motif == '+') {
      return TttsMotifOutput(me.segment);
    } else {
      return TttsMotifRevCompOutput(me.segment);
    }
  } else {
    if (me.motif == 'M')
      return TgtMotifOutput(me.segment);
    else if (me.motif == 'Y')
      return TtcMotifOutput(me.segment);
    else if (me.motif == 'R')
      return TgaMotifOutput(me.segment);
    else
      return "";
  }
}

template <typename THost, typename TString>
inline CharString const
prettyString(ModStringTriplex<THost, TString> const &me) {
  typedef typename Infix<THost>::Type TSegment;
  typedef ModifiedString<TSegment, ModView<FunctorTTSMotifPretty>>
      TttsMotifOutput;
  typedef ModifiedString<
      ModifiedString<TSegment, ModView<FunctorTTSMotifComplPretty>>, ModReverse>
      TttsMotifRevCompOutput;
  typedef ModifiedString<TSegment, ModView<FunctorTCMotifPretty>>
      TtcMotifOutput;
  typedef ModifiedString<TSegment, ModView<FunctorGTMotifPretty>>
      TgtMotifOutput;
  typedef ModifiedString<TSegment, ModView<FunctorGAMotifPretty>>
      TgaMotifOutput;

  if (!me.isTFO) {
    if (me.motif == '+') {
      return TttsMotifOutput(me.segment);
    } else {
      return TttsMotifRevCompOutput(me.segment);
    }
  } else {
    if (me.motif == 'M')
      return TgtMotifOutput(me.segment);
    else if (me.motif == 'Y')
      return TtcMotifOutput(me.segment);
    else if (me.motif == 'R')
      return TgaMotifOutput(me.segment);
    else
      return "";
  }
}

template <typename THost_, typename TString, typename TIterator>
inline void setEnd(ModStringTriplex<THost_, TString> &me, TIterator new_end) {
  me.data_end_position = new_end - begin(host(me)); //, Standard());
}

template <typename THost_, typename TString_, typename TPosition_>
inline void setEndPosition(ModStringTriplex<THost_, TString_> &me,
                           TPosition_ new_end) {
  setEndPosition(me.segment, new_end);
}

template <typename THost_, typename TString>
inline void _setLength(ModStringTriplex<THost_, TString> &me,
                       typename Size<THost_>::Type new_length) {
  me.data_end_position = me.data_begin_position + new_length;
}

template <typename THost_, typename TString>
inline void setHost(ModStringTriplex<THost_, TString> &me,
                    typename Parameter_<THost_>::Type _host) {
  me.data_host = _toPointer(_host);
}

template <typename THost, typename TString>
inline typename Size<ModStringTriplex<THost, TString> const>::Type
length(ModStringTriplex<THost, TString> const &me) {
  return endPosition(me) - beginPosition(me);
}

template <typename THost, typename TString>
inline TString substr(ModStringTriplex<THost, TString> &me, int beginPos,
                      int endPos) {
  typedef typename Iterator<TString, Rooted>::Type TIter;

  ::std::string tmp_sub;
  tmp_sub.reserve(endPos - beginPos);

  TIter it = begin(me.mask_string, Rooted());
  TIter itEnd = begin(me.mask_string, Rooted());
  goFurther(it, beginPos);
  goFurther(itEnd, endPos);
  for (; it != itEnd; ++it) {
    tmp_sub += *(it);
  }
  TString tmp_tstring(tmp_sub);
  return tmp_tstring;
}

template <typename THost_, typename TString>
struct Value<ModStringTriplex<THost_, TString>> {
  typedef typename Value<TString>::Type Type;
};

template <typename THost_, typename TString>
struct Value<ModStringTriplex<THost_, TString> const>
    : public Value<ModStringTriplex<THost_, TString>> {};

template <typename THost_, typename TString>
struct Host<ModStringTriplex<THost_, TString>> {
  typedef THost_ Type;
};

template <typename THost_, typename TString>
struct Host<ModStringTriplex<THost_, TString> const>
    : public Host<ModStringTriplex<THost_, TString>> {};

template <typename THost_, typename TString>
struct IsSequence<ModStringTriplex<THost_, TString>> {
  typedef True Type;
  enum { VALUE = true };
};

template <typename THost_, typename TString>
inline char orientationString(ModStringTriplex<THost_, TString> &me) {
  if (me.parallel)
    return 'P';
  else
    return 'A';
}

template <typename THost_, typename TString>
inline char orientationString(ModStringTriplex<THost_, TString> const &me) {
  if (me.parallel)
    return 'P';
  else
    return 'A';
}

template <typename THost_, typename TString>
inline int duplicates(ModStringTriplex<THost_, TString> &me) {
  return me.copies;
}

template <typename THost_, typename TString>
inline int duplicates(ModStringTriplex<THost_, TString> const &me) {
  return me.copies;
}

template <typename THost_, typename TString, typename TPos>
inline void duplicates(ModStringTriplex<THost_, TString> &me, TPos &copies) {
  me.copies = copies;
}

template <typename THost_, typename TString>
struct GetValue<ModStringTriplex<THost_, TString>> {
  typedef typename GetValue<TString>::Type Type;
};

template <typename THost_, typename TString>
struct GetValue<ModStringTriplex<THost_, TString> const> {
  typedef typename GetValue<TString const>::Type Type;
};

template <typename THost_, typename TString>
struct Iterator<ModStringTriplex<THost_, TString>, Rooted> {
  typedef ModStringTriplex<THost_, TString> TSequence;
  typedef typename Iterator<TString, Standard>::Type TIterator;
  typedef Iter<TSequence, AdaptorIterator<TIterator>> Type;
};

template <typename THost_, typename TString>
struct Iterator<ModStringTriplex<THost_, TString> const, Rooted> {
  typedef ModStringTriplex<THost_, TString> const TSequence;
  typedef typename Iterator<TString const, Standard>::Type TIterator;
  typedef Iter<TSequence, AdaptorIterator<TIterator>> Type;
};

template <typename THost_, typename TString>
struct Iterator<ModStringTriplex<THost_, TString>, Standard>
    : Iterator<TString, Standard> {};

template <typename THost_, typename TString>
struct Iterator<ModStringTriplex<THost_, TString> const, Standard>
    : Iterator<TString, Standard> {};

template <typename THost_, typename TString>
struct Size<ModStringTriplex<THost_, TString>> {
  typedef typename Size<THost_>::Type Type;
};

template <typename THost_, typename TString>
struct Size<ModStringTriplex<THost_, TString> const> {
  typedef typename Size<THost_>::Type Type;
};

template <typename THost_, typename TString>
struct Position<ModStringTriplex<THost_, TString>> {
  typedef typename Position<THost_>::Type Type;
};

template <typename THost_, typename TString>
struct Position<ModStringTriplex<THost_, TString> const> {
  typedef typename Position<THost_>::Type Type;
};

template <typename THost_, typename TString>
struct Id<ModStringTriplex<THost_, TString>> {
  typedef typename ModStringTriplex<THost_, TString>::TId Type;
};
template <typename THost_, typename TString>
struct Id<ModStringTriplex<THost_, TString> const> {
  typedef typename ModStringTriplex<THost_, TString>::TId Type;
};

} // namespace seqan

#endif
