/*
 * MIT License
 *
 * Copyright (c) 2022 Iñaki Amatria-Barral
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

#ifndef TRIPLEX_DEFINITIONS_HPP
#define TRIPLEX_DEFINITIONS_HPP

#include <unordered_map>
#include <vector>

#include <seqan/basic.h>
#include <seqan/graph_types.h>
#include <seqan/index.h>
#include <seqan/map.h>
#include <seqan/modifier.h>
#include <seqan/sequence.h>

#include "sequence_position.hpp"
#include "triplex_alphabet.hpp"
#include "triplex_functors.hpp"
#include "triplex_match.hpp"
#include "triplex_pattern.hpp"

// clang-format off
struct pair_hash_t {
  template <class T1, class T2>
  std::size_t operator()(const std::pair<T1, T2> &pair) const {
    return std::hash<T1>()(pair.first) ^ std::hash<T2>()(pair.second);
  }
};

struct _tts;
struct _tfo;
struct _mixed_motif;
struct _purine_motif;
struct _pyrimidine_motif;

typedef seqan::Tag<_tts> tts_t;
typedef seqan::Tag<_tfo> tfo_t;
typedef seqan::Tag<_mixed_motif> mixed_motif_t;
typedef seqan::Tag<_purine_motif> purine_motif_t;
typedef seqan::Tag<_pyrimidine_motif> pyrimidine_motif_t;

typedef seqan::TriplexString triplex_t;
typedef seqan::ModStringTriplex<triplex_t, triplex_t> motif_t;
typedef seqan::Repeat<unsigned int, unsigned int> repeat_t;
typedef typename seqan::Infix<triplex_t>::Type segment_t;
typedef seqan::ModifiedString<motif_t, seqan::ModView<seqan::FunctorRYFilter>> filter_t;
typedef seqan::TriplexMatch<seqan::Difference<seqan::TriplexString>::Type, unsigned int, unsigned int> match_t;
typedef seqan::TriplexPotential<std::pair<unsigned int, unsigned int>> potential_t;
typedef seqan::TriplexPotential<unsigned int> motif_potential_t;

typedef std::vector<triplex_t> triplex_set_t;
typedef std::vector<seqan::CharString> name_set_t;
typedef std::vector<motif_t> motif_set_t;
typedef std::vector<repeat_t> repeat_set_t;
typedef std::vector<segment_t> segment_set_t;
typedef std::vector<char> char_set_t;
typedef std::vector<char_set_t> char_set_set_t;
typedef std::vector<match_t> match_set_t;
typedef std::vector<unsigned int> index_set_t;
typedef std::vector<motif_potential_t> motif_potential_set_t;
typedef std::unordered_map<std::pair<unsigned int, unsigned int>, potential_t, pair_hash_t> potential_set_t;

#if defined(_OPENMP)
typedef std::vector<match_set_t> match_set_set_t;
#endif

typedef seqan::Graph<seqan::Automaton<seqan::Triplex, seqan::Triplex>> graph_t;
typedef seqan::VertexDescriptor<graph_t>::Type vertex_descriptor_t;

typedef seqan::StringSet<triplex_t> triplex_string_set_t;
typedef seqan::Index<triplex_string_set_t> index_t;
typedef seqan::Finder<index_t> finder_t;
typedef typename seqan::Position<motif_t>::Type duplicate_pos_t;
typedef typename seqan::Id<motif_t>::Type duplicate_id_t;
typedef seqan::SeqPos<duplicate_id_t, duplicate_pos_t> duplicate_t;
typedef seqan::Map<duplicate_t, seqan::Skiplist<>> duplicate_map_t;
typedef typename seqan::Host<motif_t>::Type host_t;
typedef typename seqan::Infix<host_t>::Type host_infix_t;
typedef typename seqan::Value<motif_t>::Type host_value_t;
typedef seqan::ModifiedString<seqan::ModifiedString<host_infix_t, seqan::ModView<seqan::FunctorComplement<host_value_t>>>, seqan::ModReverse> mod_rev_t;
// clang-format on

#endif
