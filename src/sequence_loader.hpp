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

#ifndef SEQUENCE_LOADER_HPP
#define SEQUENCE_LOADER_HPP

#include <seqan/seq_io.h>

#include "options.hpp"
#include "triplex_definitions.hpp"

struct sequence_loader_state_t {
  seqan::SeqFileIn fasta_file;
};

bool file_exists(const char *file_name);

bool create_loader_state(sequence_loader_state_t &state, const char *file_name);
void destroy_loader_state(sequence_loader_state_t &state);

bool load_sequences(triplex_set_t &sequences, name_set_t &names,
                    const char *file_name);
bool load_sequences(triplex_set_t &sequences, name_set_t &names,
                    sequence_loader_state_t &state, const options &opts);

#endif
