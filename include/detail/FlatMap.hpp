#pragma once

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <utility>
#include <vector>

namespace klc::detail {

template <typename Key, typename Value, typename Hasher>
class FlatMap {
public:
    explicit FlatMap(size_t initial_cap = 4096) {
        capacity_ = 1;
        while (capacity_ < initial_cap) {
            if (capacity_ > std::numeric_limits<size_t>::max() / 2)
                throw std::length_error("Cache capacity overflow");
            capacity_ <<= 1;
        }
        keys_.resize(capacity_);
        values_.resize(capacity_);
        occupied_.assign(capacity_, 0);
    }

    void clear() {
        std::fill(occupied_.begin(), occupied_.end(), 0);
        size_ = 0;
    }

    size_t size() const { return size_; }

    // Returned pointers are invalidated by insertion that grows the table.

    Value* find(const Key& key) {
        if (size_ == 0) return nullptr;
        size_t mask = capacity_ - 1;
        size_t idx = hasher_(key) & mask;
        while (occupied_[idx]) {
            if (keys_[idx] == key) {
                return &values_[idx];
            }
            idx = (idx + 1) & mask;
        }
        return nullptr;
    }

    void insert_or_assign(const Key& key, const Value& val) {
        if (size_ * 10 >= capacity_ * 7) {
            const Value stable_value = val; // val may refer to an entry in this map.
            grow();
            insert_or_assign(key, stable_value);
            return;
        }
        size_t mask = capacity_ - 1;
        size_t idx = hasher_(key) & mask;
        while (occupied_[idx]) {
            if (keys_[idx] == key) {
                values_[idx] = val;
                return;
            }
            idx = (idx + 1) & mask;
        }
        keys_[idx] = key;
        values_[idx] = val;
        occupied_[idx] = 1;
        size_++;
    }

private:
    void grow() {
        if (capacity_ > std::numeric_limits<size_t>::max() / 2)
            throw std::length_error("Cache capacity overflow");
        // Allocate and rehash before committing; allocation failure leaves this map usable.
        FlatMap replacement(capacity_ * 2);
        replacement.hasher_ = hasher_;
        for (size_t i = 0; i < capacity_; ++i)
            if (occupied_[i]) replacement.insert_or_assign(keys_[i], values_[i]);
        keys_.swap(replacement.keys_);
        values_.swap(replacement.values_);
        occupied_.swap(replacement.occupied_);
        std::swap(capacity_, replacement.capacity_);
        std::swap(size_, replacement.size_);
    }

    size_t capacity_ = 0;
    size_t size_ = 0;
    std::vector<Key> keys_;
    std::vector<Value> values_;
    std::vector<uint8_t> occupied_;
    Hasher hasher_{};
};


} // namespace klc::detail
