#pragma once

#include <cxxabi.h>
#include <stdint.h>
#include <array>
#include <cassert>
#include <chrono>
#include <exception>
#include <fstream>
#include <iostream>
#include <numeric>
#include <regex>
#include <sstream>
#include <string>
#include <vector>

#include <sdsl/bit_vectors.hpp>

namespace sketch_search {

static constexpr int MAX_BITS = 8;
static constexpr int MAX_DIM = 64;

enum class node_reps : int { HYBRID = 1, DHT = 2, LIST = 3 };

inline std::string get_rep_name(node_reps rep) {
    switch (rep) {
        case node_reps::HYBRID:
            return "HYBRID";
        case node_reps::DHT:
            return "DHT";
        case node_reps::LIST:
            return "LIST";
    }
    return "????????";
}

struct config_t {
    int dim;
    int bits;
    int blocks;
    float suf_thr;  // for super sparse layer
    node_reps rep_type;
};

struct score_t {
    uint32_t id;
    int errs;
};

template <class It>
inline void print_ints(std::ostream& os, It beg, It end, const char* title) {
    if (title) {
        os << title << ": ";
    }
    for (auto it = beg; it != end; ++it) {
        std::cout << static_cast<int32_t>(*it) << " ";
    }
    std::cout << std::endl;
}

template <class Vec>
inline void print_ints(std::ostream& os, const Vec& vec, const char* title) {
    print_ints(os, vec.begin(), vec.end(), title);
}

template <class T>
inline T get_max_value(uint32_t width) {
    assert(width < sizeof(T) * 8);
    return (T(1) << width) - T(1);
}

template <class T, size_t Num>
constexpr size_t array_size(const T (&array)[Num]) {
    return Num;
}

struct stat_t {
    size_t num_cands = 0;
    size_t num_actnodes = 0;
};

struct entry_t {
    const uint8_t* key;
    std::vector<uint32_t> ids;
};

inline std::vector<entry_t> make_entries(const std::vector<const uint8_t*>& keys, int dim) {
    std::vector<uint32_t> perms(keys.size());
    std::iota(perms.begin(), perms.end(), 0);
    std::sort(perms.begin(), perms.end(), [&](uint32_t i1, uint32_t i2) {
        int cmp = std::memcmp(keys[i1], keys[i2], dim);
        return cmp < 0;
    });

    std::vector<entry_t> entries;
    entries.reserve(keys.size());

    auto push_entries = [&](uint32_t beg, uint32_t end) {
        entry_t entry;
        entry.key = keys[perms[beg]];
        entry.ids.resize(end - beg);
        for (uint32_t i = beg; i < end; ++i) {
            entry.ids[i - beg] = perms[i];
        }
        entries.emplace_back(std::move(entry));
    };

    uint32_t beg = 0;
    for (uint32_t i = 1; i < keys.size(); ++i) {
        uint32_t p1 = perms[i - 1];
        uint32_t p2 = perms[i];
        if (std::memcmp(keys[p1], keys[p2], dim) == 0) {
            continue;
        }
        push_entries(beg, i);
        beg = i;
    }
    push_entries(beg, keys.size());

    entries.shrink_to_fit();
    return entries;
}

inline std::vector<std::vector<uint32_t>> parse_trie(const std::vector<entry_t>& entries, int dim) {
    std::vector<std::vector<uint32_t>> node_begs(dim + 1);
    node_begs[0] = std::vector<uint32_t>{0, static_cast<uint32_t>(entries.size())};

    for (int h = 0; h < dim; ++h) {
        node_begs[h + 1].push_back(0);
        for (uint32_t i = 1; i < node_begs[h].size(); ++i) {
            uint32_t e_beg = node_begs[h][i - 1];
            uint32_t e_end = node_begs[h][i];
            uint8_t prev_c = entries[e_beg].key[h];
            for (uint32_t j = e_beg + 1; j < e_end; ++j) {
                uint8_t cur_c = entries[j].key[h];
                assert(prev_c <= cur_c);
                if (prev_c != cur_c) {
                    node_begs[h + 1].push_back(j);
                    prev_c = cur_c;
                }
            }
            node_begs[h + 1].push_back(e_end);
        }
    }

    return node_begs;
}

template <class LhsIt, class RhsIt>
inline int get_hamdist(LhsIt lhs, RhsIt rhs, int dim, int max_errs = std::numeric_limits<int>::max()) {
    int errs = 0;
    for (int i = 0; i < dim; ++i) {
        if (lhs[i] != rhs[i]) {
            ++errs;
            if (errs > max_errs) {
                break;
            }
        }
    }
    return errs;
}

template <class LhsIt, class RhsIt>
inline int get_hamdist_v(LhsIt lhs, RhsIt rhs, int bits, int max_errs = std::numeric_limits<int>::max()) {
    int errs = 0;
    uint64_t cumdiff = 0;
    for (int j = 0; j < bits; ++j) {
        uint64_t diff = lhs[j] ^ rhs[j];
        cumdiff |= diff;
        errs = int(sdsl::bits::cnt(cumdiff));
        if (errs > max_errs) {
            return errs;
        }
    }
    return errs;
}

inline void to_vertical_code(const uint8_t* code, int bits, int dim, uint64_t* vcode) {
    for (int j = 0; j < bits; ++j) {
        uint64_t vc = 0;
        for (int i = 0; i < dim; ++i) {
            uint64_t b = (code[i] >> j) & 1ULL;
            vc |= (b << i);
        }
        vcode[j] = vc;
    }
}

// in bvecs format
inline std::vector<uint8_t> load_sketches(const std::string& fn, const config_t& conf) {
    std::ios::sync_with_stdio(false);

    std::ifstream ifs(fn);
    if (!ifs) {
        std::cerr << "open error: " << fn << '\n';
        exit(1);
    }

    uint8_t buf[MAX_DIM];
    std::vector<uint8_t> sketches;

    while (true) {
        uint32_t dim = 0;
        ifs.read(reinterpret_cast<char*>(&dim), sizeof(dim));

        if (ifs.eof()) {
            break;
        }

        if (int(dim) < conf.dim) {
            std::cerr << "error: dim < conf.dim => " << dim << " < " << conf.dim << std::endl;
            exit(1);
        }
        if (MAX_DIM < dim) {
            std::cerr << "error: MAX_DIM < dim => " << MAX_DIM << " < " << dim << std::endl;
            exit(1);
        }

        ifs.read(reinterpret_cast<char*>(buf), dim);
        std::copy(buf, buf + conf.dim, std::back_inserter(sketches));
    }
    sketches.shrink_to_fit();

    int mask = (1 << conf.bits) - 1;
    for (size_t i = 0; i < sketches.size(); ++i) {
        sketches[i] = static_cast<uint8_t>(sketches[i] & mask);
    }

    return sketches;
}

inline std::vector<const uint8_t*> extract_ptrs(const std::vector<uint8_t>& sketches, const config_t& conf) {
    if (sketches.size() % conf.dim != 0) {
        std::cerr << "error: sketches.size() % conf.dim != 0" << std::endl;
        exit(1);
    }

    std::vector<const uint8_t*> ptrs(sketches.size() / conf.dim);
    for (size_t i = 0; i < ptrs.size(); ++i) {
        ptrs[i] = sketches.data() + (i * size_t(conf.dim));
    }
    return ptrs;
}

inline std::string get_ext(const std::string& fn) {
    return fn.substr(fn.find_last_of(".") + 1);
}
inline bool is_file_exist(const std::string& fn) {
    if (fn.empty()) {
        return false;
    }
    std::ifstream ifs(fn);
    return ifs.good();
}

template <typename T>
inline std::string realname() {
    int status;
    char* name = abi::__cxa_demangle(typeid(T).name(), nullptr, nullptr, &status);

    std::string ret;
    if (name) {
        if (status == 0) {
            ret = std::string(name);
        }
        free(name);
    }
    return ret;
}

template <typename T>
inline std::string short_realname() {
    auto name = realname<T>();
    name = std::regex_replace(name, std::regex{R"( |sketch_search::)"}, "");
    name = std::regex_replace(name, std::regex{R"((\d+)ul{0,2})"}, "$1");
    return name;
}

}  // namespace sketch_search