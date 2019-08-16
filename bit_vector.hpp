#pragma once

#include "misc.hpp"

namespace sketch_search {

class bit_vector {
  public:
    using this_type = bit_vector;
    using size_type = uint64_t;

    bit_vector() = default;
    ~bit_vector() = default;

    void build(sdsl::bit_vector&& bits, bool use_rank = false, bool use_select = false) {
        m_bits = std::move(bits);
        if (use_rank) {
            sdsl::util::init_support(m_bits_r1, &m_bits);
        }
        if (use_select) {
            sdsl::util::init_support(m_bits_s1, &m_bits);
        }
    }

    void build(const std::vector<bool>& bits, bool use_rank = false, bool use_select = false) {
        m_bits = sdsl::bit_vector(bits.size());
        for (size_type i = 0; i < bits.size(); ++i) {
            m_bits[i] = bits[i];
        }
        if (use_rank) {
            sdsl::util::init_support(m_bits_r1, &m_bits);
        }
        if (use_select) {
            sdsl::util::init_support(m_bits_s1, &m_bits);
        }
    }

    bool operator[](size_type i) const {
        return get_bit(i);
    }
    bool get_bit(size_type i) const {
        return m_bits[i];
    }

    size_type rank(size_type i) const {
        return m_bits_r1(i);
    }
    size_type rank0(size_type i) const {
        return i - m_bits_r1(i);
    }
    size_type select(size_type i) const {
        return m_bits_s1(i + 1);
    }

    size_type size() const {
        return m_bits.size();
    }

    size_type serialize(std::ostream& out, sdsl::structure_tree_node* v = nullptr, std::string name = "") const {
        auto child = sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(*this));
        size_type written_bytes = 0;
        written_bytes += sdsl::serialize(m_bits, out, child, "m_bits");
        written_bytes += sdsl::serialize(m_bits_r1, out, child, "m_bits_r1");
        written_bytes += sdsl::serialize(m_bits_s1, out, child, "m_bits_s1");
        sdsl::structure_tree::add_size(child, written_bytes);
        return written_bytes;
    }

    void load(std::istream& in) {
        sdsl::load(m_bits, in);
        sdsl::load(m_bits_r1, in);
        m_bits_r1.set_vector(&m_bits);
        sdsl::load(m_bits_s1, in);
        m_bits_s1.set_vector(&m_bits);
    }

    bit_vector(const bit_vector&) = delete;
    bit_vector& operator=(const bit_vector&) = delete;

    bit_vector(bit_vector&& rhs) noexcept : bit_vector() {
        *this = std::move(rhs);
    }
    bit_vector& operator=(bit_vector&& rhs) noexcept {
        if (this != &rhs) {
            m_bits = std::move(rhs.m_bits);
            m_bits_r1 = std::move(rhs.m_bits_r1);
            m_bits_r1.set_vector(&m_bits);
            m_bits_s1 = std::move(rhs.m_bits_s1);
            m_bits_s1.set_vector(&m_bits);
        }
        return *this;
    }

  private:
    sdsl::bit_vector m_bits;
    sdsl::bit_vector::rank_1_type m_bits_r1;
    sdsl::bit_vector::select_1_type m_bits_s1;
};

}  // namespace sketch_search