#pragma once

#include "misc.hpp"
#include "sig_generator.hpp"
#include "sig_size.hpp"

namespace sketch_search {

static constexpr size_t SIG_LIMIT = 100000000;  // 100M

class hash_table {
  public:
    using size_type = uint64_t;
    static constexpr double load_factor = 1.5;

    hash_table() = default;
    ~hash_table() = default;

    void build(std::vector<const uint8_t*>& keys, const config_t& conf) {
        m_conf = conf;
        build_(keys);
    }

    class searcher {
      public:
        searcher() = default;

        const std::vector<score_t>& operator()(const uint8_t* q, int max_errs, stat_t& stat) {
            m_score.clear();
            if (max_errs < 0) {
                return m_score;
            }

            // If # of signatures exceeds # of keys, then we stop the query process
            // because the running time will be slower than that of plain linear scan
            if (get_sigsize(m_obj->m_conf.bits, m_obj->m_conf.dim, max_errs) >= SIG_LIMIT) {
                std::cerr << "**** forced termination due to massive signatures!! ****" << std::endl;
                exit(1);
            }

            for (int errs = 0; errs <= max_errs; ++errs) {
                m_gen.set(q, m_obj->m_conf.dim, m_obj->m_conf.dim, m_obj->m_conf.bits, errs);
                while (m_gen.has_next()) {
                    m_q = m_gen.next();
                    find_(errs);
                }
            }
            return m_score;
        }

      private:
        const hash_table* m_obj = nullptr;
        const uint8_t* m_q = nullptr;
        sig_generator m_gen;
        std::vector<score_t> m_score;

        searcher(const hash_table* obj) : m_obj(obj) {
            m_score.reserve(1U << 10);
        }

        void find_(int errs) {
            size_t pos = fnv1a_hash_(m_q, m_obj->m_conf.dim) % m_obj->m_table.size();

            // Linear probing
            while (true) {
                if (m_obj->m_table[pos].key_pos == UINT32_MAX) {
                    break;
                }

                auto key = m_obj->m_keys.begin() + (m_obj->m_table[pos].key_pos * m_obj->m_conf.dim);
                if (std::equal(m_q, m_q + m_obj->m_conf.dim, key)) {
                    for (uint32_t i = m_obj->m_table[pos].id_beg; i < m_obj->m_table[pos].id_end; ++i) {
                        m_score.push_back({static_cast<uint32_t>(m_obj->m_ids[i]), errs});
                    }
                    return;
                }
                ++pos;
                if (pos == m_obj->m_table.size()) {
                    pos = 0;
                }
            }
        }

        friend class hash_table;
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

    void show_stats(std::ostream& os) const {}

    size_type serialize(std::ostream& out, sdsl::structure_tree_node* v = nullptr, std::string name = "") const {
        auto child = sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(*this));
        size_type written_bytes = 0;
        written_bytes += sdsl::serialize(m_conf, out, child, "m_conf");
        written_bytes += sdsl::serialize(m_table, out, child, "m_table");
        written_bytes += sdsl::serialize(m_keys, out, child, "m_keys");
        written_bytes += sdsl::serialize(m_ids, out, child, "m_ids");
        sdsl::structure_tree::add_size(child, written_bytes);
        return written_bytes;
    }

    void load(std::istream& in) {
        sdsl::load(m_conf, in);
        sdsl::load(m_table, in);
        sdsl::load(m_keys, in);
        sdsl::load(m_ids, in);
    }

    hash_table(const hash_table&) = delete;
    hash_table& operator=(const hash_table&) = delete;

    hash_table(hash_table&& rhs) noexcept : hash_table() {
        *this = std::move(rhs);
    }
    hash_table& operator=(hash_table&& rhs) noexcept {
        if (this != &rhs) {
            m_conf = std::move(rhs.m_conf);
            m_table = std::move(rhs.m_table);
            m_keys = std::move(rhs.m_keys);
            m_ids = std::move(rhs.m_ids);
        }
        return *this;
    }

  private:
    struct element_t {
        uint32_t key_pos;
        uint32_t id_beg;
        uint32_t id_end;
    };
    config_t m_conf;
    std::vector<element_t> m_table;
    sdsl::int_vector<> m_keys;
    sdsl::int_vector<> m_ids;

    void build_(std::vector<const uint8_t*>& keys) {
        const auto entries = make_entries(keys, m_conf.dim);
        const size_t num_elems = entries.size() * load_factor;

        m_table.resize(num_elems, element_t{UINT32_MAX, 0, 0});
        m_keys = sdsl::int_vector<>(entries.size() * m_conf.dim, 0, m_conf.bits);
        m_ids = sdsl::int_vector<>(keys.size(), 0, sdsl::bits::hi(keys.size()) + 1);

        size_t id_beg = 0;

        for (size_t i = 0; i < entries.size(); ++i) {
            const auto& e = entries[i];
            size_t pos = fnv1a_hash_(e.key, m_conf.dim) % num_elems;

            // Linear probing
            while (m_table[pos].key_pos != UINT32_MAX) {
                ++pos;
                if (pos == num_elems) {
                    pos = 0;
                }
            }

            m_table[pos].key_pos = static_cast<uint32_t>(i);
            std::copy(e.key, e.key + m_conf.dim, m_keys.begin() + (i * m_conf.dim));

            m_table[pos].id_beg = id_beg;
            for (uint32_t id : e.ids) {
                m_ids[id_beg++] = id;
            }
            m_table[pos].id_end = id_beg;
        }
    }

    static size_t fnv1a_hash_(const uint8_t* key, size_t length) {
        static const size_t init = size_t((sizeof(size_t) == 8) ? 0xcbf29ce484222325 : 0x811c9dc5);
        static const size_t multiplier = size_t((sizeof(size_t) == 8) ? 0x100000001b3 : 0x1000193);

        size_t hash = init;
        for (size_t i = 0; i < length; ++i) {
            hash ^= key[i];
            hash *= multiplier;
        }
        return hash;
    }
};

}  // namespace sketch_search