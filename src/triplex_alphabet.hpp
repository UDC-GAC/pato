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
//       the names of its contributors may be used to endorse or promote products 
//       derived from this software without specific prior written permission.
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

#ifndef TRIPLEX_ALPHABET_HPP
#define TRIPLEX_ALPHABET_HPP

#include <functional>

#include <seqan/basic.h>
#include <seqan/modifier.h>
#include <seqan/sequence.h>

namespace seqan
{

typedef char Ascii;
typedef wchar_t Unicode;

template <typename T = void>
struct TranslateTableTriplex2Ascii_
{
	static char const VALUE[9];
};
template <typename T>
char const TranslateTableTriplex2Ascii_<T>::VALUE[9] = {'A', 'C', 'G', 'T', 'R', 'Y', 'K', 'M', 'N'};

template <typename T = void>
struct TranslateTableDna2Triplex_
{
	static char const VALUE[4];
};
template <typename T>
char const TranslateTableDna2Triplex_<T>::VALUE[4] =
{
	'0', //'A'
	'1', //'C'
	'2', //'G'
	'3'  //'T'
};

template <typename T = void>
struct TranslateTableDna52Triplex_
{
	static char const VALUE[5];
};
template <typename T>
char const TranslateTableDna52Triplex_<T>::VALUE[5] =
{
	'0', //'A'
	'1', //'C'
	'2', //'G'
	'3', //'T'
	'8'  //'N'
};

template <typename T = void>
struct TranslateTableIupac2Triplex_
{
	static char const VALUE[16];
};
template <typename T>
char const TranslateTableIupac2Triplex_<T>::VALUE[16] =
{
	3, //'U'
	3, //'T'
	0, //'A'
	8, //'W' = TA
	1, //'C'
	5, //'Y' = TC
	7, //'M' = AC
	8, //'H' = not-G
	2, //'G'
	6, //'K' = TG
	4, //'R' = AG
	8, //'D' = not-C
	8, //'S' = CG
	8, //'B' = non-A
	8, //'V' = non-T
	8  //'N' = any
};

template <typename T = void>
struct TranslateTableAscii2Triplex_
{
	static char const VALUE[256];
};
template <typename T>
char const TranslateTableAscii2Triplex_<T>::VALUE[256] =
{
	8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8, //0
	8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8, //1
	8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8, //2
	8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8, //3

	8,   0,   8,   1,   8,   8,   8,   2,   8,   8,   8,   6,   8,   7,   8,   8, //4
//	 ,   A,   B,   C,   D,   E,   D,   G,   H,   I,   J,   K,   L,   M,   N,   O,

	8,   8,   4,   8,   3,   3,   8,   8,   8,   5,   8,   8,   8,   8,   8,   8, //5
//	P,   Q,   R,   S,   T,   U,   V,   W,   X,   Y,   Z,    ,    ,    ,    ,

	8,   0,   8,   1,   8,   8,   8,   2,   8,   8,   8,   6,   8,   7,   8,   8, //6
//   ,   a,   b,   c,   d,   e,   f,   g,   h,   i,   j,   k,   l,   m,   n,   o,

	8,   8,   4,   8,   3,   3,   8,   8,   8,   5,   8,   8,   8,   8,   8,   8, //7
//  p,   q,   r,   s,   t,   u,   v,   w,   x,   y,   z,    ,    ,    ,    ,

	8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8, //8
	8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8, //9
	8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8, //10
	8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8, //11
	8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8, //12
	8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8, //13
	8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8, //14
	8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8  //15
};

template <typename T = void>
struct TranslateTableByte2Triplex_
{
	static char const VALUE[256];
};
template <typename T>
	char const TranslateTableByte2Triplex_<T>::VALUE[256] = {
	0,   1,   8,   3,   4,   5,   6,   7,   8,   8,   8,   8,   8,   8,   8,   8, //0
	8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8, //1
	8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8, //2
	8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8, //3
	8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8, //4
	8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8, //5
	8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8, //6
	8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8, //7
	8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8, //8
	8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8, //9
	8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8, //10
	8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8, //11
	8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8, //12
	8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8, //13
	8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8, //14
	8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8  //15
};

struct Triplex_ {};
typedef SimpleType<unsigned char, Triplex_> Triplex;

template <> struct ValueSize< Triplex > { 
	typedef uint8_t Type;
	static const Type VALUE = 9;
};
template <> struct BitsPerValue< Triplex > { 
	typedef uint8_t Type;
    static const Type VALUE = 4;
};

inline void assign(Ascii & c_target,
				   Triplex const & source)
{
	c_target = TranslateTableTriplex2Ascii_<>::VALUE[source.value];
}

template <>
struct CompareType<Triplex, uint8_t> { typedef Triplex Type; };
inline void assign(Triplex & target, uint8_t c_source)
{
	target.value = TranslateTableByte2Triplex_<>::VALUE[c_source];
}

template <>
struct CompareType<Triplex, Ascii> { typedef Triplex Type; };
inline void assign(Triplex & target, Ascii c_source)
{
	target.value = TranslateTableAscii2Triplex_<>::VALUE[(unsigned char)c_source];
}

template <>
struct CompareType<Triplex, Unicode> { typedef Triplex Type; };
inline void assign(Triplex & target, Unicode c_source)
{
	target.value = TranslateTableAscii2Triplex_<>::VALUE[(unsigned char) c_source];
}

template <>
struct CompareType<Triplex, Dna> { typedef Triplex Type; };
inline void assign(Triplex & target, Dna const & c_source)
{
	target.value = TranslateTableDna2Triplex_<>::VALUE[(unsigned char) c_source];
}

template <>
struct CompareType<Triplex, Dna5> { typedef Triplex Type; };
inline void assign(Triplex & target, Dna5 const & c_source)
{
	target.value = TranslateTableDna52Triplex_<>::VALUE[(unsigned char) c_source];
}

template <>
struct CompareType<Triplex, Iupac> { typedef Triplex Type; };
inline void assign(Triplex & target, Iupac const & source)
{
	target.value = TranslateTableIupac2Triplex_<>::VALUE[source.value];
}

template <typename T = void>
struct TranslateTableTriplex2TriplexComplement_
{
	static char const VALUE[9];
};
template <typename T>
	char const TranslateTableTriplex2TriplexComplement_<T>::VALUE[9] = {'T', 'G', 'C', 'A', 'Y', 'R', 'M', 'K'};

inline bool isMatch(Triplex val1, Triplex val2)
{
	if (val1.value == 8 || val2.value == 8) // 'N'
		return false;
	if (val1.value == val2.value)
		return true;
	else 
		return false;
}

template <typename TAlphabet>
inline bool isMatch(TAlphabet val1, TAlphabet val2)
{
	if (val1 == val2)
		return true;
	else 
		return false;
}

inline bool _repeatMaskValue(Triplex const &val) 
{
	return val == unknownValue<Triplex>(); // 'N'
}

typedef String<Triplex, Alloc<void> > TriplexString;

typedef ModView< FunctorComplement<Triplex> >	ModComplementTriplex;

typedef ModifiedString<TriplexString, ModView< FunctorComplement<Triplex> > >		TriplexStringComplement;

typedef ModifiedString<
			ModifiedString<	TriplexString, ModView< FunctorComplement<Triplex> > >,
			ModReverse
		>	TriplexStringReverseComplement;

template <>
struct FunctorComplement<Triplex> : std::function<Triplex(Triplex)>
{
	inline Triplex operator()(Triplex x) const {
		return TranslateTableTriplex2TriplexComplement_<>::VALUE[x.value];
	}
};

}

#endif
