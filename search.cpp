#include <chrono>
#include <iostream>

#include "hash_table.hpp"
#include "multi_index.hpp"
#include "sketch_trie.hpp"

#include "cmdline.h"

using namespace sketch_search;

constexpr double ABORT_BORDER_IN_MS = 1000.0;

class timer {
  public:
    using hrc = std::chrono::high_resolution_clock;

    timer() = default;

    template <class Duration>
    double get() const {
        return std::chrono::duration_cast<Duration>(hrc::now() - tp_).count();
    }

  private:
    hrc::time_point tp_ = hrc::now();
};

std::vector<std::string> string_split(const std::string& s, char delim) {
    std::vector<std::string> elems;
    std::string item;
    for (char ch : s) {
        if (ch == delim) {
            if (!item.empty()) elems.push_back(item);
            item.clear();
        } else {
            item += ch;
        }
    }
    if (!item.empty()) elems.push_back(item);
    return elems;
}

std::tuple<int, int, int> parse_range(const std::string& range_str) {
    auto elems = string_split(range_str, ':');
    if (elems.size() == 1) {
        int max = std::stoi(elems[0]);
        return {0, max, 1};
    }
    if (elems.size() == 2) {
        int min = std::stoi(elems[0]);
        int max = std::stoi(elems[1]);
        return {min, max, 1};
    }
    if (elems.size() == 3) {
        int min = std::stoi(elems[0]);
        int max = std::stoi(elems[1]);
        int stp = std::stoi(elems[2]);
        return {min, max, stp};
    }

    std::cerr << "error: invalid format of range string " << range_str << std::endl;
    exit(1);
}

template <class Index>
int bench_index(const cmdline::parser& p) {
    auto name = p.get<std::string>("name");
    auto index_fn = p.get<std::string>("index_fn");
    auto base_fn = p.get<std::string>("base_fn");
    auto query_fn = p.get<std::string>("query_fn");
    auto dim = p.get<int>("dim");
    auto bits = p.get<int>("bits");
    auto blocks = p.get<int>("blocks");
    auto errs_range = p.get<std::string>("errs_range");
    auto validation = p.get<bool>("validation");
    auto suf_thr = p.get<float>("suf_thr");

    if (dim == 0 or MAX_DIM < dim) {
        std::cerr << "error: dim == 0 or MAX_DIM < dim" << std::endl;
        return 1;
    }
    if (bits == 0 or MAX_BITS < bits) {
        std::cerr << "error: bits == 0 or MAX_BITS < bits" << std::endl;
        return 1;
    }

    std::cout << "### " << short_realname<Index>() << " ###" << std::endl;

    Index index;
    std::vector<uint8_t> keys_buf;
    std::vector<const uint8_t*> keys;
    std::vector<uint8_t> queries_buf;
    std::vector<const uint8_t*> queries;

    config_t conf;
    conf.dim = dim;
    conf.bits = bits;
    conf.blocks = blocks;
    conf.suf_thr = suf_thr;
    conf.rep_type = node_reps::HYBRID;

    if (is_file_exist(base_fn)) {
        std::cout << "Now loading keys..." << std::endl;
        timer t;
        keys_buf = load_sketches(base_fn, conf);
        keys = extract_ptrs(keys_buf, conf);
        double elapsed = t.get<std::chrono::seconds>();
        std::cout << "--> " << keys.size() << " keys" << std::endl;
        std::cout << "--> " << elapsed << " sec" << std::endl;
    }

    if (!index_fn.empty()) {
        std::ostringstream oss;
        oss << index_fn << "." << dim << "m" << bits << "b" << blocks << "B." << name;
        index_fn = oss.str();
    }

    if (is_file_exist(index_fn)) {
        std::cout << "Now loading index" << std::endl;
        sdsl::load_from_file(index, index_fn);
    } else {
        if (keys.empty()) {
            std::cerr << "error: keys is empty" << std::endl;
            return 1;
        }
        std::cout << "Now constructing index" << std::endl;
        timer t;
        index.build(keys, conf);
        double elapsed = t.get<std::chrono::seconds>();
        std::cout << "--> " << elapsed << " sec" << std::endl;
        if (!index_fn.empty()) {
            std::cout << "Now writing " << index_fn << std::endl;
            sdsl::store_to_file(index, index_fn);
        }
    }

    {
        size_t bytes = sdsl::size_in_bytes(index);
        std::cout << "--> " << bytes << " bytes; " << bytes / (1024.0 * 1024.0) << " MiB" << std::endl;
    }

    index.show_stats(std::cout);

    std::cout << "Now loading queries..." << std::endl;
    queries_buf = load_sketches(query_fn, conf);
    queries = extract_ptrs(queries_buf, conf);
    std::cout << "--> " << queries.size() << " queries" << std::endl;

    int min_errs, max_errs, err_step;
    std::tie(min_errs, max_errs, err_step) = parse_range(errs_range);

    auto searcher = index.make_searcher();

    if (validation) {
        if (keys.empty()) {
            std::cerr << "error: keys is empty" << std::endl;
            return 1;
        }

        std::cout << "Now validating with " << (min_errs + max_errs) / 2 << " errs..." << std::endl;

        stat_t stat;
        for (size_t j = 0; j < queries.size(); ++j) {
            auto& ret = searcher(queries[j], min_errs, stat);
            std::vector<score_t> searched_ans(ret.begin(), ret.end());

            std::vector<score_t> true_ans;
            for (uint32_t i = 0; i < keys.size(); ++i) {
                int hamdist = get_hamdist(keys[i], queries[j], dim, min_errs);
                if (hamdist <= min_errs) {
                    true_ans.push_back({i, hamdist});
                }
            }

            if (searched_ans.size() != true_ans.size()) {
                std::cerr << "validation error: searched_ans.size() != true_ans.size() -> " << searched_ans.size()
                          << " != " << true_ans.size() << std::endl;
                std::cerr << "  at " << j << "-th query: ";
                print_ints(std::cerr, queries[j], queries[j] + dim, nullptr);
                return 1;
            }

            std::sort(searched_ans.begin(), searched_ans.end(),
                      [](const score_t& lhs, const score_t& rhs) { return lhs.id < rhs.id; });

            for (uint32_t i = 0; i < searched_ans.size(); ++i) {
                if (searched_ans[i].id != true_ans[i].id or searched_ans[i].errs != true_ans[i].errs) {
                    std::cerr << "validation error: searched_ans[i].id != true_ans[i].id for i = " << i << std::endl;
                    std::cerr << "  at " << j << "-th query: ";
                    print_ints(std::cerr, queries[j], queries[j] + dim, nullptr);
                    return 1;
                }
            }
        }
        std::cout << "--> No problem!!" << std::endl;
        return 0;
    }

    {
        std::cout << "Now simlarity searching..." << std::endl;

        for (int errs = min_errs; errs <= max_errs; errs += err_step) {
            size_t num_ans = 0;
            stat_t stat;

            timer t;
            for (uint32_t i = 0; i < queries.size(); ++i) {
                auto ret = searcher(queries[i], errs, stat);
                num_ans += ret.size();
            }
            double elapsed = t.get<std::chrono::milliseconds>();

            std::cout << "--> " << errs << " errs; " << double(num_ans) / queries.size() << " ans; ";
            std::cout << double(stat.num_cands) / queries.size() << " cands; ";
            // std::cout << double(stat.num_actnodes) / queries.size() << " actnodes; ";
            std::cout << elapsed / queries.size() << " ms" << std::endl;

            if (ABORT_BORDER_IN_MS * queries.size() < elapsed) {
                std::cout << "**** forced termination due to ABORT_BORDER_IN_MS!! ****" << std::endl;
                break;
            }
        }
    }

    return 0;
}

int main(int argc, char* argv[]) {
    cmdline::parser p;
    p.add<std::string>("name", 'n', "index name (hash | trie)", true);
    p.add<std::string>("index_fn", 'i', "input/output file name of index", true);
    p.add<std::string>("base_fn", 'd', "input file name of database sketches", true);
    p.add<std::string>("query_fn", 'q', "input file name of query sketches", true);
    p.add<int>("dim", 'm', "dimension (<= 64)", false, 32);
    p.add<int>("bits", 'b', "#bits of alphabet (<= 8)", false, 2);
    p.add<int>("blocks", 'B', "#blocks (B=1 means to use single index)", false, 1);
    p.add<std::string>("errs_range", 'e', "range of errs (min:max:step)", false, "1:5:1");
    p.add<bool>("validation", 'v', "validation", false, false);
    p.add<float>("suf_thr", 's', "suf_thr", false, 2.0);
    p.parse_check(argc, argv);

    auto name = p.get<std::string>("name");
    auto blocks = p.get<int>("blocks");

    if (blocks == 1) {
        if (name == "hash") {
            return bench_index<hash_table>(p);
        }
        if (name == "trie") {
            return bench_index<sketch_trie>(p);
        }
    } else {
        if (name == "hash") {
            return bench_index<multi_index<hash_table>>(p);
        }
        if (name == "trie") {
            return bench_index<multi_index<sketch_trie>>(p);
        }
    }

    return 1;
}
