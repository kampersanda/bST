#pragma once

#include <bitset>
#include <iostream>

#include "misc.hpp"

namespace sketch_search {

class sig_generator {
  public:
    sig_generator() = default;

    void set(const uint8_t* base, int pfx_dim, int dim, int bits, int errs) {
        if (64 < dim) {
            std::cerr << "error: 64 < dim" << std::endl;
            exit(1);
        }

        m_base = base;
        m_pfx_dim = pfx_dim;
        m_dim = dim;
        m_mask = (1 << bits) - 1;
        m_errs = errs;

        for (int i = 0; i < errs; ++i) {
            m_power[i] = i;
        }
        m_power[errs] = m_pfx_dim + 1;
        m_bit = errs - 1;
        m_bitstr = 0;
    }

    bool has_next() const {
        return m_gen_ints or (m_bit != m_errs);
    }
    const uint8_t* next() {
        assert(has_next());

        if (m_gen_ints) {
            return next_ints();
        }

        while (m_bit != -1) {
            if (m_power[m_bit] == m_bit) {
                m_bitstr ^= 1ULL << m_power[m_bit];
            } else {
                m_bitstr ^= 3ULL << (m_power[m_bit] - 1);
            }
            ++m_power[m_bit];
            --m_bit;
        }

        uint64_t tmpstr = m_bitstr;
        while ((++m_bit < m_errs) and (m_power[m_bit] == m_power[m_bit + 1] - 1)) {
            assert(m_power[m_bit] > 0);
            m_bitstr ^= 1ULL << (m_power[m_bit] - 1);
            m_power[m_bit] = m_bit;
        }

        m_gen_ints = true;

        int r = 0;
        for (uint8_t i = 0; r < m_errs; ++i) {
            if (((tmpstr >> i) & 1ULL) == 1ULL) {
                m_combs[r] = i;
                m_chars[r] = m_base[i];
                m_cntrs[r] = 1;
                ++r;
            }
        }

        return next_ints();
    }

  private:
    const uint8_t* m_base = nullptr;
    int m_dim = 0;
    int m_pfx_dim = 0;
    int m_mask = 0;
    int m_errs = 0;

    uint8_t m_sig[64];

    // For combination
    uint64_t m_bitstr = 0;
    int m_bit = 0;
    int m_power[64];

    // For int codes
    bool m_gen_ints = false;
    uint8_t m_combs[64];
    uint8_t m_chars[64];
    uint8_t m_cntrs[64];

    const uint8_t* next_ints() {
        std::memcpy(m_sig, m_base, m_dim);

        int r = 0;
        for (r = 0; r < m_errs; ++r) {
            m_sig[m_combs[r]] = uint8_t((m_chars[r] + m_cntrs[r]) & m_mask);
        }

        for (r = 0; r < m_errs; ++r) {
            if (m_cntrs[r] < m_mask) {
                ++m_cntrs[r];
                break;
            }
            m_cntrs[r] = 1;
        }
        if (r == m_errs) {
            m_gen_ints = false;
        }
        return m_sig;
    }
};

}  // namespace sketch_search