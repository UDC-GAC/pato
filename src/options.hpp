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

#ifndef OPTIONS_HPP
#define OPTIONS_HPP

#include <seqan/sequence.h>
#include <seqan/stream.h>

#include "triplex_enums.hpp"

struct options {
  seqan::CharString tfo_file;
  seqan::CharString tts_file;
  seqan::CharString output_file;

  double error_rate;
  double min_guanine_rate;
  double max_guanine_rate;
  double mixed_parallel_max_guanine;
  double mixed_antiparallel_min_guanine;

  int min_length;
  int max_length;
  int maximal_error;

  unsigned int chunk_size;
  unsigned int min_block_run;
  unsigned int min_repeat_length;
  unsigned int max_repeat_period;
  unsigned int max_interruptions;

  run_mode_t run_mode;
  output_format_t output_format;
  error_reference_t error_reference;

  bool tc_motif;
  bool ga_motif;
  bool gt_p_motif;
  bool gt_a_motif;
  bool all_matches;
  bool pretty_output;
  bool filter_repeats;
  bool merge_features;
};

void print_options(const options &opts);

#endif
