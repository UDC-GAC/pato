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

#include "command_line_parser.hpp"
#include "options.hpp"
#include "tfo_finder.hpp"
#include "triplex_enums.hpp"
#include "triplex_finder.hpp"
#include "tts_finder.hpp"

int main(int argc, char *argv[]) {
  options opts;
  if (!parse_command_line(opts, argc, argv)) {
    return 0;
  }

  print_options(opts);
  switch (opts.run_mode) {
  case run_mode_t::tfo_search:
    find_tfo_motifs(opts);
    break;
  case run_mode_t::tts_search:
    find_tts_motifs(opts);
    break;
  case run_mode_t::tpx_search:
    find_triplexes(opts);
    break;
  }

  return 0;
}
