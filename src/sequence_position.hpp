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

#ifndef SEQUENCE_POSITION_HPP
#define SEQUENCE_POSITION_HPP

namespace seqan {

template <typename TSeq, typename TPos> struct SeqPos {
  TSeq seqnr;
  TPos position;

  bool operator==(const SeqPos<TSeq, TPos> &b) const;
  bool operator!=(const SeqPos<TSeq, TPos> &b) const;
  bool operator<(const SeqPos<TSeq, TPos> &b) const;
  bool operator>(const SeqPos<TSeq, TPos> &b) const;

  SeqPos(TSeq _seqnr, TPos _pos) : seqnr(_seqnr), position(_pos) {}

  SeqPos() {}

  template <typename TSource> inline SeqPos &operator=(TSource const &source) {
    assign(*this, source);
    return *this;
  }
};

template <typename TSeq, typename TPos>
bool SeqPos<TSeq, TPos>::operator==(const SeqPos<TSeq, TPos> &b) const {
  if (seqnr != b.seqnr)
    return false;
  if (position != b.position)
    return false;
  return true;
}

template <typename TSeq, typename TPos>
bool SeqPos<TSeq, TPos>::operator!=(const SeqPos<TSeq, TPos> &b) const {
  return !(*this == b);
}

template <typename TSeq, typename TPos>
bool SeqPos<TSeq, TPos>::operator<(const SeqPos<TSeq, TPos> &b) const {
  if (seqnr < b.seqnr)
    return true;
  if (seqnr > b.seqnr)
    return false;
  if (position < b.position)
    return true;
  return false;
}

template <typename TSeq, typename TPos>
bool SeqPos<TSeq, TPos>::operator>(const SeqPos<TSeq, TPos> &b) const {
  return b < *this;
}

template <typename TSeq, typename TPos>
inline TSeq getSequenceNo(SeqPos<TSeq, TPos> &me) {
  return me.seqnr;
}

template <typename TSeq, typename TPos>
inline TSeq getSequenceNo(SeqPos<TSeq, TPos> const &me) {
  return me.seqnr;
}

template <typename TSeq, typename TPos>
inline TPos getPosition(SeqPos<TSeq, TPos> &me) {
  return me.position;
}

template <typename TSeq, typename TPos>
inline TPos getPosition(SeqPos<TSeq, TPos> const &me) {
  return me.position;
}

} // namespace seqan

#endif
