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

#include "segment_parser.hpp"

#include <seqan/graph_types.h>
#include <seqan/sequence.h>

void make_parser(graph_t &parser, triplex_t &valid_chars,
                 triplex_t &invalid_chars, unsigned int max_interrupts) {
  vertex_descriptor_t root = seqan::addVertex(parser);
  seqan::assignRoot(parser, root);

  vertex_descriptor_t valid_child = seqan::addVertex(parser);
  for (auto valid_char : valid_chars) {
    seqan::addEdge(parser, root, valid_child, valid_char);
    seqan::addEdge(parser, valid_child, valid_child, valid_char);
  }

  vertex_descriptor_t last_child = valid_child;
  for (unsigned int i = 0; i < max_interrupts; i++) {
    vertex_descriptor_t invalid_child = seqan::addVertex(parser);

    for (auto invalid_char : invalid_chars) {
      seqan::addEdge(parser, last_child, invalid_child, invalid_char);
    }
    for (auto valid_char : valid_chars) {
      seqan::addEdge(parser, invalid_child, valid_child, valid_char);
    }

    last_child = invalid_child;
  }
}

void parse_segments(graph_t &parser, segment_set_t &segments,
                    triplex_t &sequence, unsigned int max_interrupts,
                    int min_length) {
  typedef typename seqan::Iterator<triplex_t>::Type iterator_t;

  iterator_t it = seqan::begin(sequence);
  iterator_t run_it = seqan::begin(sequence);
  iterator_t end_it = seqan::end(sequence);

  vertex_descriptor_t root = seqan::getRoot(parser);
  seqan::parseString(parser, root, run_it, end_it);

  while (run_it != end_it) {
    unsigned int shift =
        std::min(max_interrupts, static_cast<unsigned int>(run_it - it));
    run_it -= shift;

    unsigned int size = std::max(run_it - run_it, run_it - it);
    if (static_cast<int>(size) >= min_length) {
      segments.push_back(seqan::infix(sequence, it, run_it));
    }

    run_it += shift + 1;
    it = run_it;

    if (run_it != end_it) {
      seqan::parseString(parser, root, run_it, end_it);
    }
  }

  unsigned int size = std::max(end_it - end_it, end_it - it);
  if (static_cast<int>(size) >= min_length) {
    segments.push_back(seqan::infix(sequence, it, end_it));
  }
}
