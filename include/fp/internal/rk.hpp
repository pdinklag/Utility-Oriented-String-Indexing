/**
 * rk.hpp
 * part of pdinklag/fp
 * 
 * MIT License
 * 
 * Copyright (c) Patrick Dinklage, Jonas Ellert and Alexander Herlez
 * 
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 * 
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

#ifndef _FP_RK_HPP
#define _FP_RK_HPP

#include <bit>
#include <concepts>
#include <cstddef>
#include <iostream>
#include <limits>

namespace fp::internal {

template<std::unsigned_integral _Fingerprint, typename Product, int _exponent>
class RabinKarp {
public:
    /**
     * \brief Alias for the fingerprint type
     */
    using Fingerprint = _Fingerprint;

private:
    static_assert(sizeof(Product) >= 2 * sizeof(_Fingerprint));

    static constexpr Fingerprint EXPONENT = _exponent;
    static constexpr Fingerprint PRIME = (Fingerprint(1ULL) << _exponent) - 1;

    static constexpr Fingerprint mod(Product const value) {
        Fingerprint const i = (value & PRIME) + (value >> EXPONENT);
        return i - PRIME * (i >= PRIME);
    }

    static constexpr Fingerprint power(Fingerprint base, Fingerprint exponent) {
        Fingerprint result = 1;
        while (exponent > 0) {
            if (exponent & 1ULL) {
                result = mod(Product(base) * result);
            }
            base = mod(Product(base) * base);
            exponent >>= 1;
        }
        return result;
    }

    Fingerprint base_;
    Fingerprint table_[256][256] = {};

public:
    RabinKarp() : base_(0) {
    }

    /**
     * \brief Constructs a fingerprinter for the given base and window size
     * 
     * \param base the base for fingerprint computation
     * \param window the window size for rolling fingerprints (pass zero if not needed)
     */
    RabinKarp(Fingerprint const base, Fingerprint const window = 0) : base_(mod(base)) {
        if(std::has_single_bit(base)) {
            std::cerr << "[pdinklag/fp] Your Rabin-Karp fingerprint base (" << base << ") is a power of two, which is strongly recommended against. Please consult the readme for more details." << std::endl;
        }

        if(window > 0) {
            Fingerprint const max_exponent = power(base_, window);
            for (size_t i = 0; i < 256; ++i) {
                table_[i][0] = mod(PRIME - mod(i * (Product)max_exponent));
                for (size_t j = 1; j < 256; ++j) {
                    table_[i][j] = mod(table_[i][j - 1] + 1);
                }
            }
        }
    };

    RabinKarp(RabinKarp&&) = default;
    RabinKarp& operator=(RabinKarp&&) = default;

    RabinKarp(RabinKarp const&) = default;
    RabinKarp& operator=(RabinKarp const&) = default;

    /**
     * \brief Extends a fingerprint
     * 
     * \param fp the fingerprint to extend
     * \param push_right the character to append on the right
     * \return the new fingerprint
     */
    __attribute__((always_inline)) 
    inline Fingerprint push(Fingerprint const fp, uint8_t const push_right) const {
        return mod(((Product)base_) * fp + push_right);
    }

    /**
     * \brief Rolls a fingerprint
     * 
     * \param fp the fingerprint to roll
     * \param pop_left the character to remove from the beginning of the string
     * \param push_right the character to append to the end of the string
     * \return the new fingerprint
     */
    __attribute__((always_inline)) 
    inline Fingerprint roll(Fingerprint const fp, uint8_t const pop_left, uint8_t const push_right) const {
        auto const lookup = table_[pop_left][push_right];
        return mod(((Product)base_) * fp + lookup);
    }
};

}

 #endif
