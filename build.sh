#!/usr/bin/env bash
set -euo pipefail

include_flags=(-I. -Iinclude)
[[ -d third_party/bit_vector/include ]] &&
    include_flags+=(-Ithird_party/bit_vector/include)
[[ -d third_party/GSL/include ]] &&
    include_flags+=(-Ithird_party/GSL/include)
[[ -d "$HOME/.local/include" ]] && include_flags+=(-I"$HOME/.local/include")

library_flags=()
[[ -d "$HOME/.local/lib" ]] && library_flags+=(-L"$HOME/.local/lib")

find_header() {
    local header=$1
    shift
    for directory in "$@"; do
        [[ -f "$directory/$header" ]] && return 0
    done
    return 1
}

find_header pasta/bit_vector/bit_vector.hpp \
    . include third_party/bit_vector/include "$HOME/.local/include" /usr/include || {
    echo "Missing PASTA header: pasta/bit_vector/bit_vector.hpp" >&2
    exit 1
}
find_header gsl/span \
    . include third_party/GSL/include "$HOME/.local/include" /usr/include || {
    echo "Missing Microsoft GSL header: gsl/span" >&2
    exit 1
}
find_header sdsl/bit_vectors.hpp \
    "$HOME/.local/include" /usr/local/include /usr/include || {
    echo "Missing SDSL headers (for example sdsl/bit_vectors.hpp)" >&2
    exit 1
}

g++ -O3 -DNDEBUG -std=c++2a -ffunction-sections -fdata-sections \
    "${include_flags[@]}" \
    RLZ_Rank_sdsl.cpp Rank.cpp TextParser.cpp \
    WordBasedText.cpp main.cpp libsais.c libsais64.c \
    -Wl,--gc-sections \
    "${library_flags[@]}" \
    -lsdsl -ldivsufsort -ldivsufsort64 -lz -ldl \
    -o rlz_benchmark
