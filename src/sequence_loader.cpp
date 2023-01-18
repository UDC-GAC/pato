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

#include "sequence_loader.hpp"

#include <string>
#include <iostream>
#include <algorithm>

#include <seqan/stream.h>
#include <seqan/sequence.h>

void crop_sequence_name(seqan::CharString& name)
{
    std::string tmp_name(name.data_begin, seqan::length(name));
    std::size_t num_chars = std::min(tmp_name.find_first_of(' '), tmp_name.size());
    name = tmp_name.substr(0, num_chars);
}

bool file_exists(const char *file_name)
{
    seqan::SeqFileIn fasta_file;
    return seqan::open(fasta_file, file_name);
}

bool create_loader_state(sequence_loader_state_t& state, const char *file_name)
{
    if (!seqan::open(state.fasta_file, file_name)) {
        std::cerr << "PATO: error opening input file '" << file_name << "'\n";
        return false;
    }
    return true;
}

void destroy_loader_state(sequence_loader_state_t& state)
{
    seqan::close(state.fasta_file);
}

bool load_sequences(triplex_set_t& sequences,
                    name_set_t& names,
                    const char *file_name)
{
    seqan::SeqFileIn fasta_file;
    if (!seqan::open(fasta_file, file_name)) {
        std::cerr << "PATO: error opening input file '" << file_name << "'\n";
        return false;
    }

    seqan::readRecords(names, sequences, fasta_file);
    for (auto& name : names) {
        crop_sequence_name(name);
    }

    return sequences.size() > 0;
}

bool load_sequences(triplex_set_t& sequences,
                    name_set_t& names,
                    sequence_loader_state_t& state,
                    const options& opts)
{
    seqan::readRecords(names, sequences, state.fasta_file, opts.chunk_size);
    for (auto& name : names) {
        crop_sequence_name(name);
    }

    return sequences.size() > 0;
}
