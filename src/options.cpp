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

#include "options.hpp"

#include <iostream>

#include <seqan/sequence.h>

void print_options(const options &opts) {
  // clang-format off
  std::cout << "\033[1mPATO: high PerformAnce TriplexatOr\033[0m -- ";
  switch (opts.run_mode) {
  case run_mode_t::tfo_search:
    std::cout << "TFO search\n\n";
    break;
  case run_mode_t::tts_search:
    std::cout << "TTS search\n\n";
    break;
  case run_mode_t::tpx_search:
    std::cout << "TPX search\n\n";
    break;
  }

  if (opts.run_mode != run_mode_t::tts_search) {
    std::cout << "  (-ss) \033[1mTFO file:\033[0m " << opts.tfo_file << "\n";
  }
  if (opts.run_mode != run_mode_t::tfo_search) {
    std::cout << "  (-ds) \033[1mTTS file:\033[0m " << opts.tts_file << "\n";
  }
  std::cout << "\n";

  std::cout << "   (-l) \033[1mMinimum triplex length:\033[0m " << opts.min_length << "\n";
  std::cout << "   (-L) \033[1mMaximum triplex length:\033[0m " << opts.max_length << "\n";
  std::cout << "   (-e) \033[1mError rate:\033[0m " << opts.error_rate * 100.0 << "%\n";
  std::cout << "   (-E) \033[1mMaximal error rate:\033[0m " << opts.maximal_error << "\n";
  std::cout << "   (-c) \033[1mMaximum consecutive errors:\033[0m " << opts.max_interruptions << "\n";
  std::cout << "   (-g) \033[1mMinimal guanine proportion:\033[0m " << opts.min_guanine_rate * 100.0 << "%\n";
  std::cout << "   (-G) \033[1mMaximal guanine proportion:\033[0m " << opts.max_guanine_rate * 100.0 << "%\n";
  if (opts.run_mode != run_mode_t::tts_search) {
    std::cout << "   (-m) \033[1mTriplex motifs allowed\033[0m \n";
    std::cout << "        - \033[1mTC:\033[0m " << (opts.tc_motif ? "on" : "off") << "\n";
    std::cout << "        - \033[1mGA:\033[0m " << (opts.ga_motif ? "on" : "off") << "\n";
    std::cout << "        - \033[1mParallel GT:\033[0m " << (opts.gt_p_motif ? "on" : "off") << "\n";
    std::cout << "        - \033[1mAntiprallel GT:\033[0m " << (opts.gt_a_motif ? "on" : "off") << "\n";
  }
  std::cout << "(-mpmg) \033[1mMaximum guanine content in a parallel mixed-motif:\033[0m " << opts.mixed_parallel_max_guanine * 100.0 << "%\n";
  std::cout << "(-mamg) \033[1mMinimum guanine content in an anti-parallel mixed-motif:\033[0m " << opts.mixed_antiparallel_min_guanine * 100.0 << "%\n";
  std::cout << "   (-b) \033[1mRequired number of consecutive matches:\033[0m " << opts.min_block_run << "\n";
  std::cout << "   (-a) \033[1mProcess all sub-matches:\033[0m " << (opts.all_matches ? "on" : "off") << "\n\n";

  std::cout << "  (-fr) \033[1mDisregard low-complex regions:\033[0m " << (opts.filter_repeats ? "on" : "off") << "\n";
  if (opts.filter_repeats) {
    std::cout << " (-mrl) \033[1mMinimum length to filter a low-complex region:\033[0m " << opts.min_repeat_length << "\n";
    std::cout << " (-mrp) \033[1mMaximum repeat period to filter a low-complex region:\033[0m " << opts.max_repeat_period << "\n";
  }
  if (opts.run_mode != run_mode_t::tpx_search) {
    std::cout << "  (-mf) \033[1mMerge overlapping features:\033[0m " << (opts.merge_features ? "on" : "off") << "\n";
  }
  std::cout << "\n";

  std::cout << "   (-o) \033[1mOutput file name:\033[0m " << opts.output_file << "\n";
  std::cout << "  (-of) \033[1mOutput format:\033[0m ";
  if (opts.output_format == output_format_t::summary) {
    std::cout << "summary only\n";
  } else if (opts.output_format == output_format_t::bed) {
    std::cout << "BED\n";
  } else if (opts.run_mode == run_mode_t::tpx_search) {
    std::cout << "Triplexator\n";
  } else {
    std::cout << "FASTA\n";
  }
  std::cout << "  (-po) \033[1mPretty output:\033[0m " << (opts.pretty_output ? "on" : "off") << "\n";
  std::cout << "  (-er) \033[1mError reference:\033[0m ";
  switch (opts.error_reference) {
  case error_reference_t::watson_strand:
    std::cout << "watson strand\n\n";
    break;
  case error_reference_t::purine_strand:
    std::cout << "purine strand\n\n";
    break;
  case error_reference_t::last:
  case error_reference_t::third_strand:
    std::cout << "third strand\n\n";
    break;
  }

  if (opts.run_mode != run_mode_t::tfo_search) {
    std::cout << "  (-cs) \033[1mTTS window size\033[0m: " << opts.chunk_size << "\n\n";
  }
  // clang-format on
}
