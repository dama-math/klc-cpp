#pragma once

#include "Polynomial.hpp"
#include <vector>
#include <unordered_map>
#include <cstdint>
#include <cassert>

/**
 * @brief Manages a pool of unique Polynomial objects to reduce memory usage.
 * 
 * Implements the Flyweight pattern by storing unique polynomials in a single
 * pool and returning an identifier (ID) to reference them, rather than passing
 * around copies of identical polynomials.
 */
class PolynomialPool {
public:
    using IdType = uint32_t;

    /**
     * @brief Default constructor for PolynomialPool.
     */
    PolynomialPool() {
        // ID 0 can be implicitly used as the starting ID.
    }

    /**
     * @brief Gets the unique ID for a given polynomial.
     * 
     * If the polynomial does not exist in the pool, it is added, and a new
     * unique ID is generated and returned.
     * 
     * @param p The polynomial to look up or store.
     * @return The unique ID for the given polynomial.
     */
    IdType get_id(const Polynomial& p) {
        if (auto it = map.find(p); it != map.end()) {
            return it->second;
        }
        IdType id = static_cast<IdType>(pool.size());
        pool.push_back(p);
        map.emplace(p, id);
        return id;
    }

    /**
     * @brief Retrieves a polynomial from the pool by its unique ID.
     * 
     * @param id The ID of the polynomial to retrieve.
     * @return A constant reference to the requested polynomial.
     */
    const Polynomial& get_poly(IdType id) const {
        assert(id < pool.size() && "PolynomialPool: invalid polynomial ID");
        return pool[id];
    }
    
    /**
     * @brief Clears the pool and internal mappings.
     */
    void clear() {
        pool.clear();
        map.clear();
    }
    
    /**
     * @brief Returns the total number of unique polynomials in the pool.
     * 
     * @return The current size of the pool.
     */
    size_t size() const { return pool.size(); }

private:
    std::vector<Polynomial> pool;                                 ///< Stores the unique polynomials.
    std::unordered_map<Polynomial, IdType, PolynomialHash> map;   ///< Maps polynomials to their respective IDs.
};
