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

#include "output_writer.hpp"

#include <cctype>
#include <iomanip>
#include <sstream>
#include <iostream>

#include <seqan/modifier.h>
#include <seqan/sequence.h>

#include "triplex_enums.hpp"
#include "triplex_match.hpp"
#include "triplex_pattern.hpp"

seqan::CharString triplex_alignment_string(match_t& match,
                                           motif_set_t& tfo_motifs,
                                           motif_set_t& tts_motifs)
{
    std::ostringstream alignment;

    motif_t tfo_motif(seqan::host(tfo_motifs[match.tfoNo]),
                      match.oBegin,
                      match.oEnd,
                      match.parallel,
                      tfo_motifs[match.tfoNo].seqNo,
                      true,
                      match.motif);
    motif_t tts_motif(seqan::host(tts_motifs[match.ttsNo]),
                      match.dBegin,
                      match.dEnd,
                      match.parallel,
                      match.ttsSeqNo,
                      false,
                      match.strand);

    seqan::CharString tfo_ps = seqan::prettyString(tfo_motif);
    seqan::CharString tts_ps = seqan::prettyString(tts_motif);

    seqan::CharString opp(tts_ps);
    seqan::complement(opp);

    alignment << "\n";
    if (match.strand == '-') {
        seqan::reverse(opp);
        seqan::reverse(tts_ps);
        alignment << "     5'- " << opp << " -3'" << "\n";
        alignment << "TTS: 3'- " << tts_ps << " -5'" << "\n";
        alignment << "         ";

        auto tts_it = seqan::begin(tts_motif);
        auto tts_end = seqan::end(tts_motif);
        auto tfo_it = seqan::begin(tfo_motif);
        auto tfo_end = seqan::end(tfo_motif);
        while (tts_end != tts_it && tfo_end != tfo_it) {
            tts_end--;
            tfo_end--;
            if (*tts_end == *tfo_end) {
                alignment << "|";
            } else {
                alignment << "*";
            }
        }
        alignment << "\n";

        if (!tfo_motifs[match.tfoNo].parallel) {
            alignment << "TFO: 5'- " << tfo_ps << " -3'" << "\n";
        } else {
            seqan::reverse(tfo_ps);
            alignment << "TFO: 3'- " << tfo_ps << " -5'" << "\n";
        }
    } else {
        if (!tfo_motifs[match.tfoNo].parallel) {
            seqan::reverse(tfo_ps);
            alignment << "TFO: 3'- " << tfo_ps << " -5'" << "\n";
        } else {
            alignment << "TFO: 5'- " << tfo_ps << " -3'" << "\n";
        }
        alignment << "         ";

        auto tts_it = seqan::begin(tts_motif);
        auto tts_end = seqan::end(tts_motif);
        auto tfo_it = seqan::begin(tfo_motif);
        auto tfo_end = seqan::end(tfo_motif);
        while (tts_it != tts_end && tfo_it != tfo_end) {
            if (*tts_it == *tfo_it) {
                alignment << "|";
            } else {
                alignment << "*";
            }
            tts_it++;
            tfo_it++;
        }
        alignment << "\n";

        alignment << "TTS: 5'- " << tts_ps << " -3'" << "\n";
        alignment << "     3'- " << opp << " -5'" << "\n";
    }

    return alignment.str();
}

seqan::CharString triplex_error_string(match_t& match,
                                       motif_set_t& tfo_motifs,
                                       motif_set_t& tts_motifs,
                                       const options& opts)
{
    std::ostringstream errors;

    motif_t tfo_motif(seqan::host(tfo_motifs[match.tfoNo]),
                      match.oBegin,
                      match.oEnd,
                      match.parallel,
                      tfo_motifs[match.tfoNo].seqNo,
                      true,
                      match.motif);
    motif_t tts_motif(seqan::host(tts_motifs[match.ttsNo]),
                      match.dBegin,
                      match.dEnd,
                      match.parallel,
                      match.ttsSeqNo,
                      false,
                      match.strand);

    seqan::CharString tfo_ps = seqan::prettyString(tfo_motif);
    seqan::CharString tts_ps = seqan::prettyString(tts_motif);

    if (opts.error_reference == error_reference_t::purine_strand) {
        if (!tfo_motifs[match.tfoNo].parallel) {
            seqan::reverse(tfo_ps);
        }
    } else if (opts.error_reference == error_reference_t::third_strand) {
        if (!tfo_motifs[match.tfoNo].parallel) {
            seqan::reverse(tts_ps);
        }
    } else {
        if (match.strand == '-') {
            seqan::reverse(tts_ps);
            if (tfo_motifs[match.tfoNo].parallel) {
                seqan::reverse(tfo_ps);
            }
        } else {
            if (!tfo_motifs[match.tfoNo].parallel) {
                seqan::reverse(tfo_ps);
            }
        }
    }

    auto tts_it = seqan::begin(tts_motif);
    auto tts_end = seqan::end(tts_motif);
    auto tfo_it = seqan::begin(tfo_motif);
    auto tfo_end = seqan::end(tfo_motif);
    unsigned int i = 0;

    if (opts.error_reference == error_reference_t::purine_strand
        || (opts.error_reference == error_reference_t::watson_strand && match.strand == '+')
        || (opts.error_reference == error_reference_t::third_strand && tfo_motifs[match.tfoNo].parallel)) {
        while (tts_it != tts_end && tfo_it != tfo_end) {
            if (*tts_it != *tfo_it) {
                if (!isupper(seqan::value(tts_ps, i))
                    && !isupper(seqan::value(tfo_ps, i))) {
                    errors << "b" << i;
                } else if (!isupper(seqan::value(tts_ps, i))) {
                    errors << "d" << i;
                } else if (!isupper(seqan::value(tfo_ps, i))) {
                    errors << "o" << i;
                } else {
                    errors << "t" << i;
                }
            }
            tts_it++;
            tfo_it++;
            i++;
        }
    } else {
        while (tts_end != tts_it && tfo_end != tfo_it) {
            tts_end--;
            tfo_end--;
            if (*tts_end != *tfo_end) {
                if (!isupper(seqan::value(tts_ps, i))
                    && !isupper(seqan::value(tfo_ps, i))) {
                    errors << "b" << i;
                } else if (!isupper(seqan::value(tts_ps, i))) {
                    errors << "d" << i;
                } else if (!isupper(seqan::value(tfo_ps, i))) {
                    errors << "o" << i;
                } else {
                    errors << "t" << i;
                }
            }
            i++;
        }
    }

    return errors.str();
}

bool create_output_state(output_writer_state_t& state, const options& opts)
{
    {
        seqan::CharString output_file_name;
        seqan::append(output_file_name, opts.output_file);
        seqan::append(output_file_name, ".summary");

        state.summary_file = std::fopen(seqan::toCString(output_file_name), "w");
        if (!state.summary_file) {
            std::cerr << "PATO: error opening output file '"
                      << output_file_name << "'\n";
            return false;
        }

        if (opts.run_mode == run_mode_t::tfo_search) {
            std::fprintf(state.summary_file,
                         "# Sequence-ID\tTFOs (abs)\tTFOs (rel)\tGA (abs)\tGA ("
                         "rel)\tTC (abs)\tTC (rel)\tGT (abs)\tGT (rel)\n");
        } else if (opts.run_mode == run_mode_t::tts_search) {
            std::fprintf(state.summary_file,
                         "# Duplex-ID\tTTSs (abs)\tTTSs (rel)\n");
        } else {
            std::fprintf(state.summary_file,
                         "# Duplex-ID\tSequence-ID\tTotal (abs)\tTotal (rel)\tG"
                         "A (abs)\tGA (rel)\tTC (abs)\tTC (rel)\tGT (abs)\tGT ("
                         "rel)\n");
        }
    }

    state.summary_only = opts.output_format == output_format_t::summary;
    if (opts.output_format == output_format_t::summary) {
        return true;
    }

    {
        seqan::CharString output_file_name;
        seqan::append(output_file_name, opts.output_file);
        seqan::append(output_file_name, ".out");

        state.output_file = std::fopen(seqan::toCString(output_file_name), "w");
        if (!state.output_file) {
            std::cerr << "PATO: error opening output file '"
                      << output_file_name << "'\n";
            return false;
        }

        if (opts.run_mode == run_mode_t::tpx_search) {
            std::fprintf(state.output_file,
                         "# Sequence-ID\tTFO start\tTFO end\tDuplex-ID\tTTS sta"
                         "rt\tTTS end\tScore\tError-rate\tErrors\tMotif\tStrand"
                         "\tOrientation\tGuanine-rate\n");
        } else {
            if (opts.output_format == output_format_t::bed) {
                if (opts.run_mode == run_mode_t::tfo_search) {
                    std::fprintf(state.output_file,
                                 "# Sequence-ID\tStart\tEnd\tScore\tMotif\tErro"
                                 "r-rate\tErrors\tGuanine-rate\tDuplicates\tTFO"
                                 "\tDuplicate locations\n");
                } else {
                    std::fprintf(state.output_file,
                                 "# Duplex-ID\tStart\tEnd\tScore\tStrand\tError"
                                 "-rate\tErrors\tGuanine-rate\tDuplicates\tTTS"
                                 "\tDuplicate locations\n");
                }
            }
        }
    }

    return true;
}

void destroy_output_state(output_writer_state_t& state)
{
    std::fclose(state.summary_file);
    if (!state.summary_only) {
        std::fclose(state.output_file);
    }
}

void print_motifs(motif_set_t& motifs,
                  name_set_t& names,
                  output_writer_state_t& state,
                  const options& opts)
{
    if (opts.output_format == output_format_t::summary || motifs.empty()) {
        return;
    }

    unsigned int counter = 1;
    unsigned int last_sequence_id = seqan::getSequenceNo(motifs[0]);

    for (auto& m : motifs) {
        if (opts.output_format == output_format_t::bed) {
            std::fprintf(state.output_file,
                         "%s\t%lu\t%lu\t%u\t%c\t%.2g\t%s\t%.2g\t%d\t%s\t-\n",
                         seqan::toCString(names[seqan::getSequenceNo(m)]),
                         seqan::beginPosition(m),
                         seqan::endPosition(m),
                         seqan::score(m),
                         seqan::getMotif(m),
                         1.0 - static_cast<double>(seqan::score(m)) / (seqan::endPosition(m) - seqan::beginPosition(m)),
                         seqan::toCString(seqan::errorString(m)),
                         seqan::guanineRate(m),
                         seqan::duplicates(m),
                         seqan::toCString(opts.pretty_output ? seqan::prettyString(m) : seqan::outputString(m)));
        } else {
            if (last_sequence_id != seqan::getSequenceNo(m)) {
                counter = 1;
                last_sequence_id = seqan::getSequenceNo(m);
            }

            std::fprintf(state.output_file,
                         ">%s_%u\t%lu-%lu %c\t%u\t%s\t%g\t%d\t-\n%s\n",
                         seqan::toCString(names[seqan::getSequenceNo(m)]),
                         counter++,
                         seqan::beginPosition(m),
                         seqan::endPosition(m),
                         seqan::getMotif(m),
                         seqan::score(m),
                         seqan::toCString(seqan::errorString(m)),
                         seqan::guanineRate(m),
                         seqan::duplicates(m),
                         seqan::toCString(opts.pretty_output ? seqan::prettyString(m) : seqan::outputString(m)));
        }
    }
}

void print_summary(motif_potential_set_t& potentials,
                   name_set_t& names,
                   output_writer_state_t& state,
                   const options& opts)
{
    for (auto& potential : potentials) {
        if (seqan::hasCount(potential)) {
            std::fprintf(state.summary_file,
                         "%s\t%u\t%.3g",
                         seqan::toCString(names[seqan::getKey(potential)]),
                         seqan::getCounts(potential),
                         seqan::getCounts(potential) / seqan::getNorm(potential));
            if (opts.run_mode == run_mode_t::tfo_search) {
                std::fprintf(state.summary_file,
                             "\t%u\t%.3g\t%u\t%.3g\t%u\t%.3g\t",
                             seqan::getCount(potential, 'R'),
                             seqan::getCount(potential, 'R') / seqan::getNorm(potential),
                             seqan::getCount(potential, 'Y'),
                             seqan::getCount(potential, 'Y') / seqan::getNorm(potential),
                             seqan::getCount(potential, 'M'),
                             seqan::getCount(potential, 'M') / seqan::getNorm(potential));
            }
            std::fprintf(state.summary_file, "\n");
        }
    }
}

#if !defined(_OPENMP)
void print_triplex_pairs(match_set_t& matches,
                         motif_set_t& tfo_motifs,
                         name_set_t& tfo_names,
                         motif_set_t& tts_motifs,
                         name_set_t& tts_names,
                         output_writer_state_t& state,
                         const options& opts)
#else
void print_triplex_pairs(match_set_set_t& matches,
                         motif_set_t& tfo_motifs,
                         name_set_t& tfo_names,
                         motif_set_t& tts_motifs,
                         name_set_t& tts_names,
                         output_writer_state_t& state,
                         const options& opts)
#endif
{
    if (opts.output_format == output_format_t::summary) {
        return;
    }

#if !defined(_OPENMP)
    for (auto& match: matches) {
#else
    for (auto& local_matches : matches) {
        for (auto& match : local_matches) {
#endif
            auto tfo_seq_id = tfo_motifs[match.tfoNo].seqNo;
            auto tts_seq_id = match.ttsSeqNo;

            std::fprintf(state.output_file,
                         "%s\t%lu\t%lu\t%s\t%lu\t%lu\t%u\t%.2g\t%s\t%c\t%c\t%c\t%.2g",
                         seqan::toCString(tfo_names[tfo_seq_id]),
                         match.oBegin,
                         match.oEnd,
                         seqan::toCString(tts_names[tts_seq_id]),
                         match.dBegin,
                         match.dEnd,
                         match.mScore,
                         1.0 - static_cast<double>(match.mScore) / (match.dEnd - match.dBegin),
                         seqan::toCString(triplex_error_string(match, tfo_motifs, tts_motifs, opts)),
                         match.motif,
                         match.strand,
                         match.parallel ? 'P' : 'A',
                         static_cast<double>(match.guanines) / (match.dEnd - match.dBegin));
            if (opts.output_format == output_format_t::triplex) {
                std::fprintf(state.output_file,
                             "%s",
                             seqan::toCString(triplex_alignment_string(match, tfo_motifs, tts_motifs)));
            }
            std::fprintf(state.output_file, "\n");
#if defined(_OPENMP)
        }
#endif
    }
}

void print_triplex_summary(potential_set_t& potentials,
                           name_set_t& tfo_names,
                           name_set_t& tts_names,
                           output_writer_state_t& state)
{
    for (auto& potential_entry : potentials) {
        auto& potential = potential_entry.second;
        if (seqan::hasCount(potential)) {
            std::fprintf(state.summary_file,
                         "%s\t%s\t%u\t%.3g\t%u\t%.3g\t%u\t%.3g\t%u\t%.3g\t\n",
                         seqan::toCString(tts_names[seqan::getKey(potential).second]),
                         seqan::toCString(tfo_names[seqan::getKey(potential).first]),
                         seqan::getCounts(potential),
                         seqan::getCounts(potential) / seqan::getNorm(potential),
                         seqan::getCount(potential, 'R'),
                         seqan::getCount(potential, 'R') / seqan::getNorm(potential),
                         seqan::getCount(potential, 'Y'),
                         seqan::getCount(potential, 'Y') / seqan::getNorm(potential),
                         seqan::getCount(potential, 'M'),
                         seqan::getCount(potential, 'M') / seqan::getNorm(potential));
        }
    }
}
