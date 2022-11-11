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
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

#ifndef OUTPUT_WRITER_HPP
#define OUTPUT_WRITER_HPP

#include <cstdio>

#include "options.hpp"
#include "triplex_definitions.hpp"

struct output_writer_state_t
{
    bool summary_only;

    std::FILE* output_file;
    std::FILE* summary_file;
};

bool create_output_state(output_writer_state_t& state, const options& opts);
void destroy_output_state(output_writer_state_t& state);

void print_motifs(motif_set_t& motifs,
                  name_set_t& names,
                  output_writer_state_t& state,
                  const options& opts);
void print_summary(motif_potential_set_t& potentials,
                   name_set_t& names,
                   output_writer_state_t& state,
                   const options& opts);

#if !defined(_OPENMP)
void print_triplex_pairs(match_set_t& matches,
                         motif_set_t& tfo_motifs,
                         name_set_t& tfo_names,
                         motif_set_t& tts_motifs,
                         name_set_t& tts_names,
                         output_writer_state_t& state,
                         const options& opts);
#else
void print_triplex_pairs(match_set_set_t& matches,
                         motif_set_t& tfo_motifs,
                         name_set_t& tfo_names,
                         motif_set_t& tts_motifs,
                         name_set_t& tts_names,
                         output_writer_state_t& state,
                         const options& opts);
#endif
void print_triplex_summary(potential_set_t& potentials,
                           name_set_t& tfo_names,
                           name_set_t& tts_names,
                           output_writer_state_t& state);

#endif
