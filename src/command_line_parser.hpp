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

#ifndef COMMAND_LINE_PARSER_HPP
#define COMMAND_LINE_PARSER_HPP

#include <cmath>
#include <string>
#include <iostream>
#include <algorithm>

#include <seqan/sequence.h>
#include <seqan/arg_parse.h>

#include "triplex_enums.hpp"

void parse_motifs(options& opts, const std::string& motifs)
{
    opts.tc_motif = false;
    opts.ga_motif = false;
    opts.gt_p_motif = false;
    opts.gt_a_motif = false;

    if (motifs.find('A') != std::string::npos) {
        opts.ga_motif = true;
        opts.gt_a_motif = true;
    }
    if (motifs.find('M') != std::string::npos) {
        opts.gt_p_motif = true;
        opts.gt_a_motif = true;
    }
    if (motifs.find('P') != std::string::npos) {
        opts.tc_motif = true;
        opts.gt_p_motif = true;
    }
    if (motifs.find('R') != std::string::npos) {
        opts.ga_motif = true;
    }
    if (motifs.find('Y') != std::string::npos) {
        opts.tc_motif = true;
    }
}

bool parse_command_line(options& opts, int argc, char *argv[])
{
    seqan::ArgumentParser parser("PATO");

    seqan::setShortDescription(parser, "high PerformAnce TriplexatOr");
    seqan::addUsageLine(parser, "[options] {-ss tfo_file | -ds tts_file | -ss tfo_file -ds tts_file}");
    seqan::addDescription(parser, "PATO is a high performance tool for the fast and efficient detection of triple helices and triplex features in nucleotide sequences. PATO is a more accurate and modern altenative to Triplexator and TDF and functions as a drop in replacement to accelerate the triplex analyses in multicore computers.");

    seqan::setDate(parser, "October 2022");
    seqan::setVersion(parser, "v0.0.0");
    seqan::setUrl(parser, "https://github.com/amatria/pato");
    seqan::setShortCopyright(parser, "2022 Iñaki Amatria-Barral.");
    seqan::setCompilationOpts(parser, "COMPILATIONOPTS_PLACEHOLDER");

    seqan::addSection(parser, "Input options");
    seqan::addOption(parser, seqan::ArgParseOption("ss", "single-strand-file", "File in FASTA format that is searched for TFOs (e.g. RNA or DNA).", seqan::ArgParseOption::STRING));
    seqan::addOption(parser, seqan::ArgParseOption("ds", "duplex-file", "File in FASTA format that is searched for TTSs (e.g. DNA).", seqan::ArgParseOption::STRING));
    seqan::addSection(parser, "Main options");
    seqan::addOption(parser, seqan::ArgParseOption("l", "lower-length-bound", "Minimum triplex feature length required.", seqan::ArgParseOption::INTEGER));
    seqan::addOption(parser, seqan::ArgParseOption("L", "upper-length-bound", "Maximum triplex feature length permitted (disable with -1).", seqan::ArgParseOption::INTEGER));
    seqan::addOption(parser, seqan::ArgParseOption("e", "error-rate", "Set the maximal error rate tolerated in %.", seqan::ArgParseOption::DOUBLE));
    seqan::addOption(parser, seqan::ArgParseOption("E", "maximal-error", "Set the maximal overall error tolerated (disable with -1).", seqan::ArgParseOption::INTEGER));
    seqan::addOption(parser, seqan::ArgParseOption("c", "consecutive-errors", "Maximum number of consecutive errors.", seqan::ArgParseOption::INTEGER));
    seqan::addOption(parser, seqan::ArgParseOption("g", "min-guanine", "Set the minimal guanine proportion required in %.", seqan::ArgParseOption::DOUBLE));
    seqan::addOption(parser, seqan::ArgParseOption("G", "max-guanine", "Set the maximal guanine proportion allowed in %.", seqan::ArgParseOption::DOUBLE));
    seqan::addOption(parser, seqan::ArgParseOption("m", "triplex-motifs", "Triplex motifs allowed [R,Y,M,P,A].", seqan::ArgParseOption::STRING));
    seqan::addOption(parser, seqan::ArgParseOption("mpmg", "mixed-parallel-max-guanine", "Maximum guanine content to consider parallel binding in a mixed-motif in %.", seqan::ArgParseOption::DOUBLE));
    seqan::addOption(parser, seqan::ArgParseOption("mamg", "mixed-antiparallel-min-guanine", "Minimum guanine content to consider anti-parallel binding in a mixed-motif in %.", seqan::ArgParseOption::DOUBLE));
    seqan::addOption(parser, seqan::ArgParseOption("b", "minimum-block-run", "Required number of consecutive matches.", seqan::ArgParseOption::INTEGER));
    seqan::addOption(parser, seqan::ArgParseOption("a", "all-matches", "Process and report all sub-matches in addition to the longest match.", seqan::ArgParseOption::BOOL));
    seqan::addSection(parser, "Filtering options");
    seqan::addOption(parser, seqan::ArgParseOption("fr", "filter-repeats", "Disregards repeated and low-complex regions if enabled.", seqan::ArgParseOption::BOOL));
    seqan::addOption(parser, seqan::ArgParseOption("mrl", "minimum-repeat-length", "Minimum length requirement for low-complex regions to be filtered.", seqan::ArgParseOption::INTEGER));
    seqan::addOption(parser, seqan::ArgParseOption("mrp", "maximum-repeat-period", "Maximum repeat period for low-complex regions to be filtered.", seqan::ArgParseOption::INTEGER));
    seqan::addOption(parser, seqan::ArgParseOption("mf", "merge-features", "Merge overlapping features into a cluster and report the spanning region.", seqan::ArgParseOption::BOOL));
    seqan::addSection(parser, "Output options");
    seqan::addOption(parser, seqan::ArgParseOption("o", "output", "Output file name", seqan::ArgParseOption::STRING));
    seqan::addOption(parser, seqan::ArgParseOption("of", "output-format", "Set the output format [0,1,2].", seqan::ArgParseOption::INTEGER));
    seqan::addOption(parser, seqan::ArgParseOption("po", "pretty-output", "Indicate matching/mismatching characters with upper/lower case.", seqan::ArgParseOption::BOOL));
    seqan::addOption(parser, seqan::ArgParseOption("er", "error-reference", "Reference to which the error should correspond [0,1,2].", seqan::ArgParseOption::INTEGER));
    seqan::addSection(parser, "Performance options");
    seqan::addOption(parser, seqan::ArgParseOption("cs", "chunk-size", "Set the TTS window processing size (lower size equals less memory usage but implies less parallelism).", seqan::ArgParseOption::INTEGER));

    // input options
    seqan::setDefaultValue(parser, "ss", "(null)");
    seqan::setDefaultValue(parser, "ds", "(null)");

    // main options
    seqan::setDefaultValue(parser, "l", 16);
    seqan::setDefaultValue(parser, "L", 30);
    seqan::setDefaultValue(parser, "e", 5.0);
    seqan::setDefaultValue(parser, "E", -1);
    seqan::setDefaultValue(parser, "c", 1);
    seqan::setDefaultValue(parser, "g", 10.0);
    seqan::setDefaultValue(parser, "G", 100.0);
    seqan::setDefaultValue(parser, "m", "R,Y,M,P,A");
    seqan::setDefaultValue(parser, "mpmg", 100.0);
    seqan::setDefaultValue(parser, "mamg", 0.0);
    seqan::setDefaultValue(parser, "b", 1);
    seqan::setDefaultValue(parser, "a", false);

    // filtering options
    seqan::setDefaultValue(parser, "fr", true);
    seqan::setDefaultValue(parser, "mrl", 10);
    seqan::setDefaultValue(parser, "mrp", 4);
    seqan::setDefaultValue(parser, "mf", false);

    // output options
    seqan::setDefaultValue(parser, "o", "pato");
    seqan::setDefaultValue(parser, "of", static_cast<unsigned int>(output_format_t::bed));
    seqan::setDefaultValue(parser, "po", false);
    seqan::setDefaultValue(parser, "er", static_cast<unsigned int>(error_reference_t::watson_strand));

    // performance options
    seqan::setDefaultValue(parser, "cs", 128);

    if (seqan::parse(parser, argc, argv) != seqan::ArgumentParser::PARSE_OK) {
        return false;
    }

    // input options
    seqan::getOptionValue(opts.tfo_file, parser, "ss");
    seqan::getOptionValue(opts.tts_file, parser, "ds");

    // main options
    seqan::getOptionValue(opts.min_length, parser, "l");
    seqan::getOptionValue(opts.max_length, parser, "L");
    seqan::getOptionValue(opts.error_rate, parser, "e");
    seqan::getOptionValue(opts.maximal_error, parser, "E");
    seqan::getOptionValue(opts.max_interruptions, parser, "c");
    seqan::getOptionValue(opts.min_guanine_rate, parser, "g");
    seqan::getOptionValue(opts.max_guanine_rate, parser, "G");
    seqan::getOptionValue(opts.mixed_parallel_max_guanine, parser, "mpmg");
    seqan::getOptionValue(opts.mixed_antiparallel_min_guanine, parser, "mamg");
    seqan::getOptionValue(opts.min_block_run, parser, "b");
    seqan::getOptionValue(opts.all_matches, parser, "a");

    std::string motifs;
    seqan::getOptionValue(motifs, parser, "m");
    parse_motifs(opts, motifs);

    // filtering options
    seqan::getOptionValue(opts.filter_repeats, parser, "fr");
    seqan::getOptionValue(opts.min_repeat_length, parser, "mrl");
    seqan::getOptionValue(opts.max_repeat_period, parser, "mrp");
    seqan::getOptionValue(opts.merge_features, parser, "mf");

    // output options
    seqan::getOptionValue(opts.output_file, parser, "o");
    seqan::getOptionValue(opts.pretty_output, parser, "po");

    unsigned int tmp;
    seqan::getOptionValue(tmp, parser, "of");
    opts.output_format = output_format_t(tmp);
    seqan::getOptionValue(tmp, parser, "er");
    opts.error_reference = error_reference_t(tmp);

    // performance options
    seqan::getOptionValue(opts.chunk_size, parser, "cs");

    // check options
    if (opts.tfo_file == "(null)" && opts.tts_file == "(null)") {
        std::cerr << "PATO: at least one type of input file has to be supplied\n";
        return false;
    }
    if (opts.min_length < 5) {
        std::cerr << "PATO: the minimum triplex length must be greater or equal to 5\n";
        return false;
    }
    if (opts.max_length > 1000) {
        std::cerr << "PATO: the maximum triplex length must be smaller or equal to 1000\n";
        return false;
    }
    if (opts.error_rate < 0.0 || opts.error_rate > 20.0) {
        std::cerr << "PATO: the error rate must be a value between 0 and 20\n";
        return false;
    }
    if (opts.max_interruptions > 3) {
        std::cerr << "PATO: the maximum consecutive interruptions must be smaller or equal to 3\n";
        return false;
    }
    if (opts.min_guanine_rate < 0.0 || opts.min_guanine_rate > 100.0) {
        std::cerr << "PATO: the minimum guanine proportion in the triplex target site must be a value between 0 and 100\n";
        return false;
    }
    if (opts.max_guanine_rate < 0.0 || opts.max_guanine_rate > 100.0) {
        std::cerr << "PATO: the maximum guanine proportion in the triplex target site must be a value between 0 and 100\n";
        return false;
    }
    if (opts.min_guanine_rate > opts.max_guanine_rate) {
        std::cerr << "PATO: the maximum guanine proportion cannot be smaller than the minimum guanine proportion\n";
        return false;
    }
    if (opts.mixed_parallel_max_guanine < 0.0 || opts.mixed_parallel_max_guanine > 100.0) {
        std::cerr << "PATO: the maximum guanine proportion in a parallel mixed motif must be a value between 0 and 100\n";
        return false;
    }
    if (opts.mixed_antiparallel_min_guanine < 0.0 || opts.mixed_antiparallel_min_guanine > 100.0) {
        std::cerr << "PATO: the minimum guanine proportion in an anti-parallel mixed motif must be a value between 0 and 100\n";
        return false;
    }
    if (opts.error_reference >= error_reference_t::last) {
        std::cerr << "PATO: error reference not known\n";
        return false;
    }
    if (opts.output_format >= output_format_t::last) {
        std::cerr << "PATO: output format not known\n";
        return false;
    }
    if (opts.chunk_size == 0) {
        std::cerr << "PATO: chunk size must be a positive number\n";
        return false;
    }

    unsigned int tolerated_error = static_cast<unsigned int>(std::floor(opts.error_rate * opts.min_length));
    if (opts.min_block_run > opts.min_length - 2 * tolerated_error) {
        std::cerr << "PATO: block match too large given minimum length constraint and error rate\n";
        return false;
    }

    // prepare options
    if (opts.tts_file == "(null)") {
        opts.run_mode = run_mode_t::tfo_search;
    } else if (opts.tfo_file == "(null)") {
        opts.run_mode = run_mode_t::tts_search;
    } else {
        opts.run_mode = run_mode_t::tpx_search;
    }

    opts.error_rate /= 100.0;
    opts.min_guanine_rate /= 100.0;
    opts.max_guanine_rate /= 100.0;
    opts.mixed_parallel_max_guanine /= 100.0;
    opts.mixed_antiparallel_min_guanine /= 100.0;

    if (opts.error_rate == 0.0 || opts.maximal_error == 0) {
        opts.max_interruptions = 0;
    }

    if (opts.max_length >= opts.min_length) {
        if (opts.maximal_error < 0) {
            opts.maximal_error = static_cast<int>(opts.error_rate * opts.max_length);
        } else {
            opts.maximal_error = std::min(opts.maximal_error,
                                          static_cast<int>(opts.error_rate * opts.max_length));
        }
    }

    return true;
}

#endif
