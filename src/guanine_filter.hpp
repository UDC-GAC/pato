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

#ifndef GUANINE_FILTER_HPP
#define GUANINE_FILTER_HPP

#include <algorithm>
#include <cmath>

#include "motif_set_reducer.hpp"
#include "options.hpp"
#include "triplex_definitions.hpp"
#include "triplex_enums.hpp"
#include "triplex_pattern.hpp"

struct filter_arguments {
  motif_set_t &motifs;

  char_set_set_t &block_runs;
  char_set_set_t &encoded_seq;

  orientation_t ornt;

  bool reduce_set;

  char filter_char;
  char interrupt_char;

  filter_arguments(motif_set_t &_motifs, char_set_set_t &_block_runs,
                   char_set_set_t &_encoded_seq)
      : motifs(_motifs), block_runs(_block_runs), encoded_seq(_encoded_seq) {}
};

template <typename string_t>
void encode_sequence(string_t &motif, char filter_char, char interrupt_char,
                     char_set_set_t &block_runs, char_set_set_t &encoded_seq,
                     unsigned int min_block_run) {
  if (encoded_seq.empty() || seqan::length(motif) > encoded_seq[0].size()) {
    encoded_seq.clear();
    encoded_seq.resize(3, char_set_t(seqan::length(motif) * 2, false));
  } else {
    encoded_seq[0].assign(seqan::length(motif), false);
    encoded_seq[1].assign(seqan::length(motif), false);
    encoded_seq[2].assign(seqan::length(motif), false);
  }

  unsigned int counter = 0;
  unsigned int run_counter = 0;

  for (auto m : motif) {
    if (m == filter_char) {
      encoded_seq[0][counter] = true;
    } else if (m == interrupt_char) {
      encoded_seq[1][counter] = true;
      encoded_seq[2][counter] = true;

      if (counter - run_counter >= min_block_run) {
        for (unsigned int i = 0; i <= counter - min_block_run; i++) {
          for (unsigned int j = std::max(i, run_counter) + min_block_run;
               j <= seqan::length(motif); j++) {
            block_runs[i][j] = true;
          }
        }
      }

      run_counter = counter + 1;
    } else {
      encoded_seq[2][counter] = true;
    }
    counter++;
  }

  if (counter - run_counter >= min_block_run) {
    for (unsigned int i = 0; i <= counter - min_block_run; i++) {
      for (unsigned int j = std::max(i, run_counter) + min_block_run;
           j <= seqan::length(motif); j++) {
        block_runs[i][j] = true;
      }
    }
  }
}

void increase_right(char_set_set_t &encoded_seq, unsigned int &pos,
                    unsigned int &filter_chars, unsigned int &interrupt_chars,
                    unsigned int &non_filter_chars);

void increase_left(char_set_set_t &encoded_seq, unsigned int &pos,
                   unsigned int &filter_chars, unsigned int &interrupt_chars,
                   unsigned int &non_filter_chars);

bool is_interrupt_char(char_set_set_t &encoded_seq, unsigned int pos);

template <typename tag_t>
bool motif_specific_constraint(__attribute__((unused)) double filter_rate,
                               __attribute__((unused)) orientation_t ornt,
                               __attribute__((unused)) const tag_t &tag,
                               __attribute__((unused)) const options &opts) {
  return true;
}

bool motif_specific_constraint(double filter_rate, orientation_t ornt,
                               const mixed_motif_t &tag, const options &opts);

void add_match(motif_set_t &motifs, motif_t &motif, unsigned int start,
               unsigned int end, unsigned int errors, const tfo_t &tag);

void add_match(motif_set_t &motifs, motif_t &motif, unsigned int start,
               unsigned int end, unsigned int errors, const mixed_motif_t &tag);

void add_match(motif_set_t &motifs, motif_t &motif, unsigned int start,
               unsigned int end, unsigned int errors,
               const purine_motif_t &tag);

void add_match(motif_set_t &motifs, motif_t &motif, unsigned int start,
               unsigned int end, unsigned int errors,
               const pyrimidine_motif_t &tag);

void add_match(motif_set_t &motifs, motif_t &motif, unsigned int start,
               unsigned int end, unsigned int errors, const tts_t &tag);

template <typename tag_t>
unsigned int filter_guanine_error_rate(motif_t &motif, filter_arguments &args,
                                       const tag_t &tag, const options &opts) {
  motif_set_t tmp_set;
  motif_set_t &motifs_ptr = args.reduce_set ? tmp_set : args.motifs;

  auto motif_length = seqan::length(motif);
  if (args.block_runs.empty() ||
      motif_length - opts.min_block_run + 1 > args.block_runs.size()) {
    args.block_runs.clear();
    args.block_runs.resize((motif_length - opts.min_block_run + 1) * 2,
                           char_set_t((motif_length + 1) * 2, false));
  } else {
    for (unsigned int i = 0; i < motif_length - opts.min_block_run + 1; i++) {
      args.block_runs[i].assign(motif_length + 1, false);
    }
  }

  char filter_char = args.filter_char;
  char interrupt_char = args.interrupt_char;
  if (opts.min_guanine_rate <= 0.0) {
    filter_t filtered_sequence(motif);
    filter_char = filter_char == 'G' ? 'R' : 'Y';
    encode_sequence(filtered_sequence, filter_char, interrupt_char,
                    args.block_runs, args.encoded_seq, opts.min_block_run);
  } else {
    encode_sequence(motif, filter_char, interrupt_char, args.block_runs,
                    args.encoded_seq, opts.min_block_run);
  }

  double max_error = std::floor(motif_length * opts.error_rate);
  double max_tolerated =
      std::floor(motif_length * (1.0 - opts.min_guanine_rate));
  if (opts.maximal_error >= 0) {
    max_error = std::min(max_error, static_cast<double>(opts.maximal_error));
  }

  unsigned int filter_chars = 0;
  unsigned int interrupt_chars = 0;
  unsigned int non_filter_chars = 0;

  unsigned int max_length = motif_length;
  if (opts.max_length >= opts.min_length) {
    max_length = static_cast<unsigned int>(opts.max_length);
  }

  unsigned int tmp_start = 0;
  unsigned int tmp_end = 0;
  unsigned int tmp_errors = 0;
  unsigned int covered_end = 0;

  unsigned int left = 0;
  unsigned int right = 0;

  unsigned int matches = 0;
  while (args.block_runs[left][motif_length] &&
         left + opts.min_length <= motif_length) {
    while (static_cast<int>(right - left) < opts.min_length &&
           right < motif_length) {
      while (static_cast<int>(right - left) < opts.min_length &&
             right < motif_length) {
        increase_right(args.encoded_seq, right, filter_chars, interrupt_chars,
                       non_filter_chars);
      }
      while (interrupt_chars > max_error) {
        increase_left(args.encoded_seq, left, filter_chars, interrupt_chars,
                      non_filter_chars);
      }
      while (non_filter_chars > max_tolerated) {
        increase_left(args.encoded_seq, left, filter_chars, interrupt_chars,
                      non_filter_chars);
      }
      while (left < motif_length && is_interrupt_char(args.encoded_seq, left)) {
        increase_left(args.encoded_seq, left, filter_chars, interrupt_chars,
                      non_filter_chars);
      }

      if (right < left) {
        right = left;

        filter_chars = 0;
        interrupt_chars = 0;
        non_filter_chars = 0;
      }
    }

    if (static_cast<int>(right - left) < opts.min_length) {
      break;
    }

    bool is_match = false;

    while (interrupt_chars <= max_error && non_filter_chars <= max_tolerated &&
           right - left <= max_length) {
      double filter_chars_rate =
          static_cast<double>(filter_chars) / (right - left);
      double interrupt_chars_rate =
          static_cast<double>(interrupt_chars) / (right - left);

      if (args.block_runs[left][right] &&
          !is_interrupt_char(args.encoded_seq, right - 1) &&
          interrupt_chars_rate <= opts.error_rate &&
          opts.min_guanine_rate <= filter_chars_rate &&
          filter_chars_rate <= opts.max_guanine_rate &&
          motif_specific_constraint(filter_chars_rate, args.ornt, tag, opts)) {
        is_match = true;
        matches++;

        tmp_start = left;
        tmp_end = right;
        tmp_errors = interrupt_chars;

        if (opts.all_matches) {
          covered_end = tmp_end;
          is_match = false;
          add_match(motifs_ptr, motif, tmp_start, tmp_end, tmp_errors, tag);
        }
      }

      if (right < motif_length) {
        increase_right(args.encoded_seq, right, filter_chars, interrupt_chars,
                       non_filter_chars);
      } else {
        break;
      }
    }

    if (is_match && tmp_end > covered_end) {
      covered_end = tmp_end;
      add_match(motifs_ptr, motif, tmp_start, tmp_end, tmp_errors, tag);
    }

    left++;
    while (left < motif_length && is_interrupt_char(args.encoded_seq, left)) {
      left++;
    }

    right = left;

    filter_chars = 0;
    interrupt_chars = 0;
    non_filter_chars = 0;
  }

  if (args.reduce_set) {
    reduce_motif_set(args.motifs, motifs_ptr);
  }

  return matches;
}

#endif
