#pragma once

#include <vector>
#include <functional>
#include <cstdint>
#include <utility>
#include <limits>
#include <stdexcept>

template <typename Key, typename Value, typename Hasher = std::hash<Key>>
class FastMap {
public:
    FastMap(size_t initial_capacity = 1024) {
        resize(initial_capacity);
    }

    void clear() {
        std::fill(occupied.begin(), occupied.end(), 0); 
        _size = 0;
    }

    Value& operator[](const Key& key) {
        if (_size * 2 > _capacity) {
            grow();
        }

        size_t idx = hasher(key) & _mask;
        while (occupied[idx]) {
            if (keys[idx] == key) {
                return values[idx];
            }
            idx = (idx + 1) & _mask;
        }

        keys[idx] = key;
        occupied[idx] = 1;
        _size++;
        return values[idx];
    }

    bool count(const Key& key) const {
        if (_capacity == 0) return false;
        size_t idx = hasher(key) & _mask;
        while (occupied[idx]) {
            if (keys[idx] == key) {
                return true;
            }
            idx = (idx + 1) & _mask;
        }
        return false;
    }

    const Value& at(const Key& key) const {
        if (_capacity > 0) {
            size_t idx = hasher(key) & _mask;
            while (occupied[idx]) {
                if (keys[idx] == key) {
                    return values[idx];
                }
                idx = (idx + 1) & _mask;
            }
        }
        throw std::runtime_error("Key not found");
    }

private:
    std::vector<Key> keys;
    std::vector<Value> values;
    std::vector<uint8_t> occupied;
    
    size_t _capacity = 0;
    size_t _size = 0;
    size_t _mask = 0;
    Hasher hasher;

    void resize(size_t new_cap) {
        size_t cap = 1;
        while (cap < new_cap) cap <<= 1;
        
        std::vector<Key> new_keys(cap);
        std::vector<Value> new_values(cap);
        std::vector<uint8_t> new_occupied(cap, 0);
        
        _capacity = cap;
        _mask = cap - 1;

        for (size_t i = 0; i < occupied.size(); ++i) {
            if (occupied[i]) {
                size_t idx = hasher(keys[i]) & _mask;
                while (new_occupied[idx]) {
                    idx = (idx + 1) & _mask;
                }
                new_keys[idx] = std::move(keys[i]);
                new_values[idx] = std::move(values[i]);
                new_occupied[idx] = 1;
            }
        }
        keys = std::move(new_keys);
        values = std::move(new_values);
        occupied = std::move(new_occupied);
    }

    void grow() {
        resize(_capacity == 0 ? 1024 : _capacity * 2);
    }
};
