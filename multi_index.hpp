#pragma once

#include <functional>
#include <numeric>

#include "misc.hpp"

namespace sketch_search {

template <class Index>
class multi_index {
  public:
    using this_type = multi_index<Index>;
    using index_type = Index;
    using size_type = uint64_t;

  public:
    multi_index() = default;
    ~multi_index() = default;

    void build(const std::vector<const uint8_t*>& keys, const config_t& conf) {
        m_conf = conf;

        if (m_conf.blocks < 2) {
            std::cerr << "error: blocks < 2" << std::endl;
            exit(1);
        }

        m_dims.resize(m_conf.blocks);
        m_indexes.resize(m_conf.blocks);

        std::vector<const uint8_t*> sub_keys(keys.size());
        int dim_beg = 0;

        config_t conf_b = m_conf;
        for (int b = 0; b < m_conf.blocks; ++b) {
            m_dims[b] = (int(m_conf.dim) + b) / m_conf.blocks;

            for (size_t i = 0; i < keys.size(); ++i) {
                sub_keys[i] = keys[i] + dim_beg;
            }
            conf_b.dim = m_dims[b];
            m_indexes[b].build(sub_keys, conf_b);
            dim_beg += m_dims[b];
        }

        uint64_t vcode[MAX_BITS];
        m_vert_codes = sdsl::int_vector<>(keys.size() * uint64_t(conf.bits), 0, conf.dim);
        for (size_t i = 0; i < keys.size(); ++i) {
            to_vertical_code(keys[i], conf.bits, conf.dim, vcode);
            std::copy(vcode, vcode + conf.bits, m_vert_codes.begin() + i * size_t(conf.bits));
        }
    }

    class searcher {
      public:
        using index_searcher_type = typename index_type::searcher;

        searcher() = default;

        const std::vector<score_t>& operator()(const uint8_t* q, int max_errs, stat_t& stat) {
            m_score.clear();
            reset_dupflags_();

            static uint64_t vq[MAX_BITS];
            to_vertical_code(q, m_obj->m_conf.bits, m_obj->m_conf.dim, vq);

            int blocks = m_obj->num_blocks();
            {
                float gph_errs = max_errs - blocks + 1;
                for (int b = 0; b < blocks; ++b) {
                    sub_errs_[b] = std::floor((gph_errs + b) / blocks);
                }
                assert(std::accumulate(sub_errs_.begin(), sub_errs_.end(), 0) == int(gph_errs));
            }

            for (int b = 0; b < blocks; ++b) {
                const uint8_t* sub_q = q + dim_begs_[b];

                // cand0, err0, cand1, err1, cand2, err2, ...
                const std::vector<score_t>& cands = index_searchers_[b](sub_q, sub_errs_[b], stat);

                for (size_t i = 0; i < cands.size(); ++i) {
                    uint32_t cand = cands[i].id;

                    if (get_dupflag_(cand)) {
                        continue;
                    }

                    ++stat.num_cands;

                    uint64_t offset = cand * uint64_t(m_obj->m_conf.bits);
                    const auto vcode = m_obj->m_vert_codes.begin() + offset;
                    int hamdist = get_hamdist_v(vq, vcode, m_obj->m_conf.bits, max_errs);

                    if (hamdist <= max_errs) {
                        m_score.push_back({cand, hamdist});
                    }
                    set_dupflag_(cand);
                }
            }

            return m_score;
        }

      private:
        const this_type* m_obj = nullptr;
        std::vector<score_t> m_score;
        std::vector<uint64_t> dupflags_;
        std::vector<int> sub_errs_;
        std::vector<int> dim_begs_;
        std::vector<index_searcher_type> index_searchers_;

        explicit searcher(const this_type* obj) : m_obj(obj) {
            int blocks = m_obj->num_blocks();

            m_score.reserve(1U << 10);
            dupflags_.resize(m_obj->num_keys() / 64 + 1);
            sub_errs_.resize(blocks);
            dim_begs_.resize(blocks + 1);

            int dim_beg = 0;
            for (int b = 0; b < blocks; ++b) {
                dim_begs_[b] = dim_beg;
                dim_beg += m_obj->m_dims[b];
                index_searchers_.emplace_back(m_obj->m_indexes[b].make_searcher());
            }
            dim_begs_[blocks] = dim_beg;
        }

        bool get_dupflag_(uint64_t i) const {
            assert(i / 64 < dupflags_.size());
            return (dupflags_[i / 64] & (1ULL << (i % 64))) != 0;
        }
        void set_dupflag_(uint64_t i) {
            assert(i / 64 < dupflags_.size());
            dupflags_[i / 64] |= (1ULL << (i % 64));
        }
        void reset_dupflags_() {
            for (uint64_t i = 0; i < dupflags_.size(); i++) {
                dupflags_[i] = 0;
            }
        }

        friend class multi_index;
    };  // searcher

    searcher make_searcher() const {
        return searcher(this);
    }

    void debug_dump(std::ostream& os) const {}

    uint64_t num_keys() const {
        return m_indexes[0].num_keys();
    }
    int num_blocks() const {
        return m_conf.blocks;
    }
    config_t get_config() const {
        return m_conf;
    }

    void show_stats(std::ostream& os) const {
        for (int b = 0; b < m_conf.blocks; ++b) {
            m_indexes[b].show_stats(os);
        }
    }

    size_type serialize(std::ostream& out, sdsl::structure_tree_node* v = nullptr, std::string name = "") const {
        auto child = sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(*this));
        size_type written_bytes = 0;
        written_bytes += sdsl::serialize(m_conf, out, child, "m_conf");
        written_bytes += sdsl::serialize(m_dims, out, child, "m_dims");
        written_bytes += sdsl::serialize(m_indexes, out, child, "m_indexes");
        written_bytes += sdsl::serialize(m_vert_codes, out, child, "m_vert_codes");
        sdsl::structure_tree::add_size(child, written_bytes);
        return written_bytes;
    }

    void load(std::istream& in) {
        sdsl::load(m_conf, in);
        sdsl::load(m_dims, in);
        sdsl::load(m_indexes, in);
        sdsl::load(m_vert_codes, in);
    }

    multi_index(const multi_index&) = delete;
    multi_index& operator=(const multi_index&) = delete;

    multi_index(multi_index&& rhs) noexcept : multi_index() {
        *this = std::move(rhs);
    }
    multi_index& operator=(multi_index&& rhs) noexcept {
        if (this != &rhs) {
            m_conf = std::move(rhs.m_conf);
            m_dims = std::move(rhs.m_dims);
            m_indexes = std::move(rhs.m_indexes);
            m_vert_codes = std::move(rhs.m_vert_codes);
        }
        return *this;
    }

  private:
    config_t m_conf;
    std::vector<int> m_dims;
    std::vector<index_type> m_indexes;
    sdsl::int_vector<> m_vert_codes;
};

}  // namespace sketch_search