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

#include "tfo_finder.hpp"

#include <algorithm>
#include <iostream>
#include <numeric>

#include <seqan/parallel/parallel_macros.h>

#include "guanine_filter.hpp"
#include "output_writer.hpp"
#include "repeat_filter.hpp"
#include "segment_parser.hpp"
#include "sequence_loader.hpp"
#include "triplex_enums.hpp"
#include "triplex_match.hpp"

struct tfo_arguments {
  graph_t tc_parser;
  graph_t ga_parser;
  graph_t gt_parser;

  repeat_set_t repeats;

#if !defined(_OPENMP)
  motif_set_t &motifs;
  motif_potential_set_t &potentials;
#else
  motif_set_t motifs;
  motif_potential_set_t potentials;
#endif

  char_set_set_t block_runs;
  char_set_set_t encoded_seq;

  segment_set_t segments;

  filter_arguments filter_args;

#if !defined(_OPENMP)
  explicit tfo_arguments(motif_set_t &_motifs,
                         motif_potential_set_t &_potentials)
      : motifs(_motifs), potentials(_potentials),
        filter_args(motifs, block_runs, encoded_seq) {
    filter_args.filter_char = 'G';
    filter_args.interrupt_char = 'N';
    filter_args.reduce_set = true;
  }
#else
  tfo_arguments() : filter_args(motifs, block_runs, encoded_seq) {
    filter_args.filter_char = 'G';
    filter_args.interrupt_char = 'N';
    filter_args.reduce_set = true;
  }
#endif
};

void make_tfo_parsers(tfo_arguments &args, unsigned int max_interrupts) {
  triplex_t valid_chars, invalid_chars;

  valid_chars = "TCY";
  invalid_chars = "GARN";
  make_parser(args.tc_parser, valid_chars, invalid_chars, max_interrupts);

  valid_chars = "GAR";
  invalid_chars = "TCYN";
  make_parser(args.ga_parser, valid_chars, invalid_chars, max_interrupts);

  valid_chars = "GTK";
  invalid_chars = "CAMN";
  make_parser(args.gt_parser, valid_chars, invalid_chars, max_interrupts);
}

void find_tfo_motifs(triplex_t &sequence, unsigned int id,
                     tfo_arguments &tfo_args, const options &opts) {
  unsigned int matches_y = 0;
  unsigned int matches_r = 0;
  unsigned int matches_m = 0;

  // TC motif
  if (opts.tc_motif) {
    parse_segments(tfo_args.tc_parser, tfo_args.segments, sequence,
                   opts.max_interruptions, opts.min_length);

    tfo_args.filter_args.ornt = orientation_t::parallel;

    for (auto &segment : tfo_args.segments) {
      motif_t motif(segment, true, id, true, 'Y');
      matches_y += filter_guanine_error_rate(motif, tfo_args.filter_args,
                                             pyrimidine_motif_t(), opts);
    }
    tfo_args.segments.clear();
  }

  // GA motif
  if (opts.ga_motif) {
    parse_segments(tfo_args.ga_parser, tfo_args.segments, sequence,
                   opts.max_interruptions, opts.min_length);

    tfo_args.filter_args.ornt = orientation_t::antiparallel;

    for (auto &segment : tfo_args.segments) {
      motif_t motif(segment, false, id, true, 'R');
      matches_r += filter_guanine_error_rate(motif, tfo_args.filter_args,
                                             purine_motif_t(), opts);
    }
    tfo_args.segments.clear();
  }

  // GT motif
  if (opts.gt_a_motif || opts.gt_p_motif) {
    parse_segments(tfo_args.gt_parser, tfo_args.segments, sequence,
                   opts.max_interruptions, opts.min_length);

    if (opts.gt_a_motif && opts.gt_p_motif &&
        opts.run_mode != run_mode_t::tfo_search) {
      tfo_args.filter_args.ornt = orientation_t::both;
    } else if (opts.gt_p_motif || opts.run_mode == run_mode_t::tfo_search) {
      tfo_args.filter_args.ornt = orientation_t::parallel;
    } else {
      tfo_args.filter_args.ornt = orientation_t::antiparallel;
    }

    for (auto &segment : tfo_args.segments) {
      if ((tfo_args.filter_args.ornt == orientation_t::both ||
           tfo_args.filter_args.ornt == orientation_t::parallel) &&
          opts.mixed_parallel_max_guanine >= opts.min_guanine_rate) {
        motif_t motif(segment, true, id, true, 'M');
        matches_m += filter_guanine_error_rate(motif, tfo_args.filter_args,
                                               mixed_motif_t(), opts);
      }
      if ((tfo_args.filter_args.ornt == orientation_t::both ||
           tfo_args.filter_args.ornt == orientation_t::antiparallel) &&
          opts.mixed_antiparallel_min_guanine <= opts.max_guanine_rate) {
        motif_t motif(segment, false, id, true, 'M');
        matches_m += filter_guanine_error_rate(motif, tfo_args.filter_args,
                                               mixed_motif_t(), opts);
      }
    }
    tfo_args.segments.clear();
  }

  if (opts.run_mode == run_mode_t::tfo_search) {
    motif_potential_t potential(id);

    seqan::addCount(potential, matches_y, 'Y');
    seqan::addCount(potential, matches_r, 'R');
    seqan::addCount(potential, matches_m, 'M');

    seqan::setNorm(potential, seqan::length(sequence), opts);

    tfo_args.potentials.push_back(potential);
  }
}

void find_tfo_motifs(motif_set_t &motifs, motif_potential_set_t &potentials,
                     triplex_set_t &sequences, const options &opts) {
  // clang-format off
  index_set_t indices(sequences.size(), 0);
  std::iota(indices.begin(), indices.end(), 0);
  std::sort(indices.begin(), indices.end(),
            [&](unsigned int i, unsigned int j) -> bool {
              return seqan::length(sequences[i]) > seqan::length(sequences[j]);
            });

SEQAN_OMP_PRAGMA(parallel)
{
#if !defined(_OPENMP)
  tfo_arguments tfo_args(motifs, potentials);
#else
  tfo_arguments tfo_args;
#endif
  if (opts.run_mode == run_mode_t::tfo_search) {
    tfo_args.filter_args.reduce_set = opts.merge_features;
  }

  make_tfo_parsers(tfo_args, opts.max_interruptions);

SEQAN_OMP_PRAGMA(for schedule(dynamic) nowait)
  for (unsigned int i = 0; i < sequences.size(); i++) {
    if (opts.filter_repeats) {
      filter_repeats(tfo_args.repeats, sequences[indices[i]],
                     opts.min_repeat_length, opts.max_repeat_period);
    }
    find_tfo_motifs(sequences[indices[i]], indices[i], tfo_args, opts);
  }

#if defined(_OPENMP)
#pragma omp critical(motifs)
{
  motifs.reserve(motifs.size() + tfo_args.motifs.size());
  std::move(tfo_args.motifs.begin(), tfo_args.motifs.end(),
            std::back_inserter(motifs));
} // #pragma omp critical
  if (opts.run_mode == run_mode_t::tfo_search) {
#pragma omp critical(potentials)
{
    potentials.reserve(potentials.size() + tfo_args.potentials.size());
    std::move(tfo_args.potentials.begin(), tfo_args.potentials.end(),
              std::back_inserter(potentials));
} // #pragma omp critical
  }
#endif
} // #pragma omp parallel
  // clang-format on
}

void find_tfo_motifs(const options &opts) {
  // clang-format off
  if (!file_exists(seqan::toCString(opts.tfo_file))) {
    std::cerr << "PATO: error opening input file '" << opts.tfo_file << "'\n";
    return;
  }

  output_writer_state_t tfo_output_file_state;
  if (!create_output_state(tfo_output_file_state, opts)) {
    return;
  }

  name_set_t tfo_names;
  motif_set_t tfo_motifs;
  triplex_set_t tfo_sequences;
  motif_potential_set_t tfo_potentials;

  if (!load_sequences(tfo_sequences, tfo_names,
                      seqan::toCString(opts.tfo_file))) {
    return;
  }
  find_tfo_motifs(tfo_motifs, tfo_potentials, tfo_sequences, opts);

SEQAN_OMP_PRAGMA(parallel sections num_threads(2))
{
SEQAN_OMP_PRAGMA(section)
  print_motifs(tfo_motifs, tfo_names, tfo_output_file_state, opts);
SEQAN_OMP_PRAGMA(section)
  print_summary(tfo_potentials, tfo_names, tfo_output_file_state, opts);
} // #pragma omp parallel sections num_threads(2)

  destroy_output_state(tfo_output_file_state);

  std::cout << "\033[1mTFO search:\033[0m done\n";
  // clang-format on
}
