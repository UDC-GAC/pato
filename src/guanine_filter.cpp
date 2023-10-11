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

#include "guanine_filter.hpp"

void increase_right(char_set_set_t &encoded_seq, unsigned int &pos,
                    unsigned int &filter_chars, unsigned int &interrupt_chars,
                    unsigned int &non_filter_chars) {
  filter_chars += encoded_seq[0][pos];
  interrupt_chars += encoded_seq[1][pos];
  non_filter_chars += encoded_seq[2][pos];

  pos++;
}

void increase_left(char_set_set_t &encoded_seq, unsigned int &pos,
                   unsigned int &filter_chars, unsigned int &interrupt_chars,
                   unsigned int &non_filter_chars) {
  filter_chars -= encoded_seq[0][pos];
  interrupt_chars -= encoded_seq[1][pos];
  non_filter_chars -= encoded_seq[2][pos];

  pos++;
}

bool is_interrupt_char(char_set_set_t &encoded_seq, unsigned int pos) {
  return encoded_seq[1][pos];
}

bool motif_specific_constraint(double filter_rate, orientation_t ornt,
                               __attribute__((unused)) const mixed_motif_t &tag,
                               const options &opts) {
  if (ornt != orientation_t::parallel &&
      filter_rate >= opts.mixed_antiparallel_min_guanine) {
    return true;
  } else if (ornt != orientation_t::antiparallel &&
             filter_rate <= opts.mixed_parallel_max_guanine) {
    return true;
  }
  return false;
}

void add_match(motif_set_t &motifs, motif_t &motif, unsigned int start,
               unsigned int end, unsigned int errors,
               __attribute__((unused)) const tfo_t &tag) {
  auto motif_length = seqan::length(motif);
  auto st = seqan::isParallel(motif) ? start : motif_length - end;
  auto nd = seqan::isParallel(motif) ? end : motif_length - start;

  motif_t tmp_motif(seqan::host(motif), seqan::beginPosition(motif) + st,
                    seqan::beginPosition(motif) + nd, seqan::isParallel(motif),
                    seqan::getSequenceNo(motif), seqan::isTFO(motif),
                    seqan::getMotif(motif));
  seqan::setScore(tmp_motif, end - start - errors);

  motifs.push_back(std::move(tmp_motif));
}

void add_match(motif_set_t &motifs, motif_t &motif, unsigned int start,
               unsigned int end, unsigned int errors,
               __attribute__((unused)) const mixed_motif_t &tag) {
  add_match(motifs, motif, start, end, errors, tfo_t());
}

void add_match(motif_set_t &motifs, motif_t &motif, unsigned int start,
               unsigned int end, unsigned int errors,
               __attribute__((unused)) const purine_motif_t &tag) {
  add_match(motifs, motif, start, end, errors, tfo_t());
}

void add_match(motif_set_t &motifs, motif_t &motif, unsigned int start,
               unsigned int end, unsigned int errors,
               __attribute__((unused)) const pyrimidine_motif_t &tag) {
  add_match(motifs, motif, start, end, errors, tfo_t());
}

void add_match(motif_set_t &motifs, motif_t &motif, unsigned int start,
               unsigned int end, unsigned int errors,
               __attribute__((unused)) const tts_t &tag) {
  auto motif_length = seqan::length(motif);
  auto st = seqan::getMotif(motif) == '+' ? start : motif_length - end;
  auto nd = seqan::getMotif(motif) == '+' ? end : motif_length - start;

  motif_t tmp_motif(seqan::host(motif), seqan::beginPosition(motif) + st,
                    seqan::beginPosition(motif) + nd, seqan::isParallel(motif),
                    seqan::getSequenceNo(motif), seqan::isTFO(motif),
                    seqan::getMotif(motif));
  seqan::setScore(tmp_motif, end - start - errors);

  motifs.push_back(std::move(tmp_motif));
}
