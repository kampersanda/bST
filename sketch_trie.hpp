#pragma once

#include "bit_vector.hpp"
#include "misc.hpp"

namespace sketch_search {

// #define UNDEFINE_DENSE_LAYER

class sketch_trie {
  public:
    using size_type = uint64_t;  // for sdsl

    sketch_trie() = default;
    ~sketch_trie() = default;

    void build(std::vector<const uint8_t*>& keys, const config_t& conf) {
        m_conf = conf;
        build_trie(keys);
    }

    class searcher {
      public:
        searcher() = default;

        const std::vector<score_t>& operator()(const uint8_t* q, int max_errs, stat_t& stat) {
            m_score.clear();
            if (max_errs < 0) {
                return m_score;
            }

            m_q = q;
            m_max_errs = max_errs;

            if (m_obj->m_suf_dim != 0) {
                to_vertical_code(m_q + m_trie_height, m_obj->m_conf.bits, m_obj->m_suf_dim, m_q_vert_suf);
            }
            ph_traverse_(0, 0, 0);

            return m_score;
        }

      private:
        const sketch_trie* m_obj = nullptr;
        const uint8_t* m_q = nullptr;
        uint64_t m_q_vert_suf[MAX_BITS];
        const int m_sigma = 0;
        const int m_trie_height = 0;
        int m_max_errs = 0;
        std::vector<score_t> m_score;

        searcher(const sketch_trie* obj)
            : m_obj(obj), m_sigma(1 << obj->m_conf.bits), m_trie_height(obj->m_conf.dim - obj->m_suf_dim) {
            m_score.reserve(1U << 10);
        }

        void ph_traverse_(int h, int errs, uint64_t rank) {
            if (h == m_obj->m_perf_height) {
                traverse_(h, errs, rank);
                return;
            }

            const int c = int(m_q[h]);
            rank *= m_sigma;

            if (errs == m_max_errs) {
                ph_traverse_(h + 1, errs, rank + c);
            } else {
                for (int i = 0; i < m_sigma; ++i) {
                    if (i == c) {
                        ph_traverse_(h + 1, errs, rank + i);
                    } else {
                        ph_traverse_(h + 1, errs + 1, rank + i);
                    }
                }
            }
        }

        void traverse_(int h, int errs, uint64_t rank) {
            assert(0 <= errs and errs <= m_max_errs);

            if (h == m_trie_height) {
                if (m_obj->m_suf_dim != 0) {
                    uint64_t suf_beg = m_obj->m_suf_begs.select(rank);
                    uint64_t suf_end = suf_beg;

                    assert(suf_beg + 1 < m_obj->m_suf_begs.size());

                    do {
                        const auto vert_suf = m_obj->m_vert_sufs.begin() + suf_end * m_obj->m_conf.bits;
                        int hamdist = get_hamdist_v(vert_suf, m_q_vert_suf, m_obj->m_conf.bits, m_max_errs - errs);

                        if (errs + hamdist <= m_max_errs) {
                            uint64_t id_beg = m_obj->m_id_begs.select(suf_end);
                            uint64_t id_end = id_beg;

                            assert(id_beg + 1 < m_obj->m_id_begs.size());

                            int e = errs + hamdist;
                            do {
                                m_score.push_back({static_cast<uint32_t>(m_obj->m_ids[id_end]), e});
                            } while (!m_obj->m_id_begs[++id_end]);
                        }
                    } while (!m_obj->m_suf_begs[++suf_end]);
                } else {
                    uint64_t id_beg = m_obj->m_id_begs.select(rank);
                    uint64_t id_end = id_beg;

                    assert(id_beg + 1 < m_obj->m_id_begs.size());

                    do {
                        m_score.push_back({static_cast<uint32_t>(m_obj->m_ids[id_end]), errs});
                    } while (!m_obj->m_id_begs[++id_end]);
                }
                return;
            }

            const medium_aux_t& med_aux = m_obj->m_medium_auxes[h - m_obj->m_perf_height];
            uint64_t c = m_q[h];

            if (med_aux.nd_type == DHT) {  // DHT
                uint64_t pos_beg = med_aux.begin + (rank << m_obj->m_conf.bits);
                assert(pos_beg + m_sigma <= m_obj->m_dhts.size());

                if (errs == m_max_errs) {
                    uint64_t pos = pos_beg + c;
                    if (!m_obj->m_dhts[pos]) {
                        return;
                    }
                    uint64_t next_rank = m_obj->m_dhts.rank(pos) - med_aux.prefix_sum;
                    traverse_(h + 1, errs, next_rank);
                    return;
                }

                uint64_t next_rank = m_obj->m_dhts.rank(pos_beg) - med_aux.prefix_sum;

                for (uint64_t i = 0; i < uint64_t(m_sigma); ++i) {
                    uint64_t pos = pos_beg + i;
                    if (!m_obj->m_dhts[pos]) {
                        continue;
                    }
                    traverse_(h + 1, i == c ? errs : errs + 1, next_rank++);
                }
            } else {  // List
                uint64_t pos = m_obj->m_list_bits.select(rank + med_aux.prefix_sum);
                if (errs == m_max_errs) {
                    do {
                        if (m_obj->m_list_chars[pos] == c) {
                            traverse_(h + 1, errs, pos - med_aux.begin);
                        }
                    } while (!m_obj->m_list_bits[++pos]);
                    return;
                }
                do {
                    if (m_obj->m_list_chars[pos] == c) {
                        traverse_(h + 1, errs, pos - med_aux.begin);
                    } else {
                        traverse_(h + 1, errs + 1, pos - med_aux.begin);
                    }
                } while (!m_obj->m_list_bits[++pos]);
            }
        }

        friend class sketch_trie;
    };  // searcher

    searcher make_searcher() const {
        return searcher(this);
    }

    uint64_t num_keys() const {
        return m_ids.size();
    }
    config_t get_config() const {
        return m_conf;
    }

    uint64_t get_trie_memory() const {
        return sdsl::size_in_bytes(*this) - (sdsl::size_in_bytes(m_ids) + sdsl::size_in_bytes(m_id_begs));
    }

    void show_stats(std::ostream& os) const {
        os << "Statistics of sketch_trie\n";
        os << "--> perf_height: " << m_perf_height << '\n';
        os << "--> suff_dim: " << m_suf_dim << '\n';
        os << "--> suf_thr: " << m_conf.suf_thr << '\n';
        os << "--> rep_type: " << get_rep_name(m_conf.rep_type) << std::endl;
    }

    size_type serialize(std::ostream& out, sdsl::structure_tree_node* v = nullptr, std::string name = "") const {
        auto child = sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(*this));
        size_type written_bytes = 0;
        written_bytes += sdsl::serialize(m_conf, out, child, "m_conf");
        written_bytes += sdsl::serialize(m_perf_height, out, child, "m_perf_height");
        written_bytes += sdsl::serialize(m_medium_auxes, out, child, "m_medium_auxes");
        written_bytes += sdsl::serialize(m_dhts, out, child, "m_dhts");
        written_bytes += sdsl::serialize(m_list_bits, out, child, "m_list_bits");
        written_bytes += sdsl::serialize(m_list_chars, out, child, "m_list_chars");
        written_bytes += sdsl::serialize(m_suf_dim, out, child, "m_suf_dim");
        written_bytes += sdsl::serialize(m_vert_sufs, out, child, "m_vert_sufs");
        written_bytes += sdsl::serialize(m_suf_begs, out, child, "m_suf_begs");
        written_bytes += sdsl::serialize(m_ids, out, child, "m_ids");
        written_bytes += sdsl::serialize(m_id_begs, out, child, "m_id_begs");
        sdsl::structure_tree::add_size(child, written_bytes);
        return written_bytes;
    }

    void load(std::istream& in) {
        sdsl::load(m_conf, in);
        sdsl::load(m_perf_height, in);
        sdsl::load(m_medium_auxes, in);
        sdsl::load(m_dhts, in);
        sdsl::load(m_list_bits, in);
        sdsl::load(m_list_chars, in);
        sdsl::load(m_suf_dim, in);
        sdsl::load(m_vert_sufs, in);
        sdsl::load(m_suf_begs, in);
        sdsl::load(m_ids, in);
        sdsl::load(m_id_begs, in);
    }

    sketch_trie(const sketch_trie&) = delete;
    sketch_trie& operator=(const sketch_trie&) = delete;

    sketch_trie(sketch_trie&& rhs) noexcept : sketch_trie() {
        *this = std::move(rhs);
    }
    sketch_trie& operator=(sketch_trie&& rhs) noexcept {
        if (this != &rhs) {
            m_conf = std::move(rhs.m_conf);
            m_perf_height = std::move(rhs.m_perf_height);
            m_medium_auxes = std::move(rhs.m_medium_auxes);
            m_dhts = std::move(rhs.m_dhts);
            m_list_bits = std::move(rhs.m_list_bits);
            m_list_chars = std::move(rhs.m_list_chars);
            m_suf_dim = std::move(rhs.m_suf_dim);
            m_vert_sufs = std::move(rhs.m_vert_sufs);
            m_suf_begs = std::move(rhs.m_suf_begs);
            m_ids = std::move(rhs.m_ids);
            m_id_begs = std::move(rhs.m_id_begs);
        }
        return *this;
    }

  private:
    enum ds_types : uint8_t { DHT, LIST };

    struct medium_aux_t {
        ds_types nd_type;
        uint64_t begin;
        uint64_t prefix_sum;
    };

    config_t m_conf;

    // Super dense layer
    int m_perf_height = 0;

    // Medium layer
    std::vector<medium_aux_t> m_medium_auxes;
    bit_vector m_dhts;
    bit_vector m_list_bits;
    sdsl::int_vector<> m_list_chars;

    // Super sparse layer
    int m_suf_dim = 0;
    sdsl::int_vector<> m_vert_sufs;  // in vcodes
    bit_vector m_suf_begs;  // suffix to ids

    // ID Lists
    sdsl::int_vector<> m_ids;
    bit_vector m_id_begs;  // suffix to ids

    void build_trie(std::vector<const uint8_t*>& keys) {
        auto entries = make_entries(keys, m_conf.dim);
        auto node_begs = parse_trie(entries, m_conf.dim);

        auto num_leaves = [&](int h) -> uint64_t { return node_begs[h].size() - 1; };

        // 1. Super dense layer
        int h = 0;
#ifdef UNDEFINE_DENSE_LAYER
        std::cerr << "!! UNDEFINE_DENSE_LAYER !!" << std::endl;
#else
        for (; h < m_conf.dim; ++h) {
            if ((num_leaves(h) << m_conf.bits) != num_leaves(h + 1)) {
                break;
            }
        }
#endif
        m_perf_height = h;

        // 2. Medium dense layer
        {
            std::vector<bool> dhts;
            std::vector<bool> list_bits;
            std::vector<uint8_t> list_chars;

            medium_aux_t dht_aux = {DHT, 0, 0};
            medium_aux_t list_aux = {LIST, 0, 0};

            float ds_thr = 0.0;  // dense or sparse

            switch (m_conf.rep_type) {
                case node_reps::HYBRID:
                    ds_thr = float(1 << m_conf.bits) / (m_conf.bits + 1);
                    break;
                case node_reps::DHT:
                    ds_thr = 0.0;
                    break;
                case node_reps::LIST:
                    ds_thr = float(1 << m_conf.bits) + 1;
                    break;
                default:
                    std::cerr << "invalid rep type" << std::endl;
                    exit(1);
            }

            for (; h < m_conf.dim; ++h) {
                if (num_leaves(h + 1) * m_conf.suf_thr > entries.size()) {
                    break;
                }

                const float ave_degree = float(num_leaves(h + 1)) / num_leaves(h);
                const ds_types ds_type = (ave_degree >= ds_thr) ? DHT : LIST;

                uint64_t dht_beg = dhts.size();
                if (ds_type == DHT) {
                    dhts.resize(dhts.size() + (num_leaves(h) << m_conf.bits));
                }

                const auto& prev_begs = node_begs[h];

                for (uint32_t i = 1; i < prev_begs.size(); ++i) {
                    uint32_t e_beg = prev_begs[i - 1];
                    uint32_t e_end = prev_begs[i];
                    uint8_t prev_c = entries[e_beg].key[h];

                    uint64_t _dht_beg = dht_beg + ((i - 1) << m_conf.bits);

                    if (ds_type == LIST) {
                        list_bits.push_back(true);
                    }

                    for (uint32_t j = e_beg + 1; j < e_end; ++j) {
                        uint8_t cur_c = entries[j].key[h];

                        if (prev_c != cur_c) {
                            if (cur_c < prev_c) {
                                std::cerr << "error: non lex" << std::endl;
                                exit(1);
                            }

                            if (ds_type == DHT) {
                                assert(_dht_beg + prev_c < dhts.size());
                                assert(!dhts[_dht_beg + prev_c]);
                                dhts[_dht_beg + prev_c] = true;
                            } else {  // ds_type == LIST
                                list_bits.push_back(false);
                                list_chars.push_back(prev_c);
                            }

                            prev_c = cur_c;
                        }
                    }

                    if (ds_type == DHT) {
                        assert(_dht_beg + prev_c < dhts.size());
                        assert(!dhts[_dht_beg + prev_c]);
                        dhts[_dht_beg + prev_c] = true;
                    } else {  // ds_type == LIST
                        list_chars.push_back(prev_c);
                    }
                }

                if (ds_type == DHT) {  // DHT
                    m_medium_auxes.push_back(dht_aux);
                    dht_aux.begin = dhts.size();
                    dht_aux.prefix_sum += num_leaves(h + 1);
                } else {  // List
                    m_medium_auxes.push_back(list_aux);
                    list_aux.begin = list_bits.size();
                    list_aux.prefix_sum += num_leaves(h);
                }
            }

            list_bits.push_back(true);
            list_chars.push_back('\0');

            m_medium_auxes.shrink_to_fit();
            m_dhts.build(dhts, true);
            m_list_bits.build(list_bits, false, true);
            m_list_chars = sdsl::int_vector<>(list_chars.size(), 0, m_conf.bits);
            std::copy(list_chars.begin(), list_chars.end(), m_list_chars.begin());
        }

        // 3. suffixes
        m_suf_dim = m_conf.dim - h;

        sdsl::bit_vector suf_begs;
        sdsl::bit_vector id_begs;

        if (m_suf_dim != 0) {
            m_vert_sufs = sdsl::int_vector<>(entries.size() * m_conf.bits, 0, m_suf_dim);
            suf_begs = sdsl::bit_vector(entries.size() + 1);
        }

        m_ids = sdsl::int_vector<>(keys.size(), 0, sdsl::bits::hi(keys.size()) + 1);
        id_begs = sdsl::bit_vector(keys.size() + 1);

        size_t ids_size = 0;
        size_t sufs_size = 0;

        const auto& prev_begs = node_begs[h];

        uint64_t vsuf[MAX_BITS];
        for (uint32_t i = 1; i < prev_begs.size(); ++i) {
            uint32_t e_beg = prev_begs[i - 1];
            uint32_t e_end = prev_begs[i];

            if (m_suf_dim != 0) {
                suf_begs[sufs_size] = 1;
            }

            for (uint32_t j = e_beg; j < e_end; ++j) {
                const auto& e = entries[j];

                if (m_suf_dim != 0) {
                    to_vertical_code(e.key + h, m_conf.bits, m_suf_dim, vsuf);
                    std::copy(vsuf, vsuf + m_conf.bits, m_vert_sufs.begin() + sufs_size++ * m_conf.bits);
                }

                id_begs[ids_size] = 1;
                for (uint32_t id : e.ids) {
                    m_ids[ids_size++] = id;
                }
            }
        }

        assert(ids_size == m_ids.size());

        if (m_suf_dim != 0) {
            suf_begs[sufs_size++] = 1;
            m_suf_begs.build(std::move(suf_begs), false, true);
        }

        id_begs[ids_size++] = 1;
        m_id_begs.build(std::move(id_begs), false, true);
    }
};

}  // namespace sketch_search