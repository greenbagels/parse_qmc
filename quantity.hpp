/*! @file quantity.hpp
 *  @brief Interface file for the math quantity template class
 *  @author Sameed Pervaiz (pervaiz.8@osu.edu)
 *  @date 2021-06-03
 *  @copyright GPLv3
 */

#ifndef QUANTITY_HPP
#define QUANTITY_HPP

#include <array>
#include <cmath>
#include <ostream>
#include <stdexcept>
#include <vector>

namespace qmc
{
    /*! Wrapper class for scientific measurements */
    template <typename T>
    class quantity
    {
        public:
            using value_type = T;
            using reference = value_type&;
            using const_reference = const value_type&;

            quantity(T v = 0, T uc = 0) : _value(v), _uncertainty(uc)
            {
            }

            quantity(std::initializer_list<T> l)
            {
                if (l.size() != 2)
                {
                    throw std::invalid_argument("Initializer list must have size 2!");
                }

                _value = l.begin()[0];
                _uncertainty = l.begin()[1];
            }

            /*! Returns the quantity's actual value */
            reference value()
            {
                return _value;
            }

            /*! Returns the quantity's actual value */
            const_reference value() const
            {
                return _value;
            }

            /*! Return the quantity's uncertainty (typically 1 stdev) */
            reference uncertainty()
            {
                return _uncertainty;
            }

            /*! Return the quantity's uncertainty (typically 1 stdev) */
            const_reference uncertainty() const
            {
                return _uncertainty;
            }

            /*! Add other quantity to self with propagated error */
            quantity& operator+=(const quantity& rhs)
            {
                _value += rhs.value();
                _uncertainty = std::hypot(_uncertainty, rhs.uncertainty());
                return *this;
            }

            /*! Reuse compound assignment to propagate sum errors */
            friend quantity operator+(quantity lhs, const quantity &rhs)
            {
                lhs += rhs;
                return lhs;
            }

            /*! Subtract another quantity to self with propagated error */
            quantity& operator-=(const quantity& rhs)
            {
                _value -= rhs.value();
                _uncertainty = std::hypot(_uncertainty, rhs.uncertainty());
                return *this;
            }

            /*! Reuse compound assignment to propagate difference errors */
            friend quantity operator-(quantity lhs, const quantity &rhs)
            {
                lhs -= rhs;
                return lhs;
            }

            /*! Negate a quantity */
            friend quantity operator-(const quantity& rhs)
            {
                return quantity(-rhs.value(), rhs.uncertainty());
            }

            quantity& operator*=(const quantity &rhs)
            {
                // First, grab the dX / X terms and square
                auto term1 = _uncertainty / _value;
                auto term2 = rhs.uncertainty() / rhs.value();
                // Now, update the value and uncertainty
                _value *= rhs.value();
                _uncertainty = _value * std::hypot(term1, term2);
                return *this;
            }

            friend quantity operator*(quantity lhs, const quantity &rhs)
            {
                lhs *= rhs;
                return lhs;
            }

            friend quantity operator*(const T &lhs, quantity rhs)
            {
                rhs.value() *= lhs;
                rhs.uncertainty() *= std::abs(lhs);
                return rhs;
            }

            friend quantity operator*(quantity lhs, const T &rhs)
            {
                return rhs * lhs;
            }

            quantity& operator/=(const quantity &rhs)
            {
                // First, grab the dX / X terms and square
                auto term1 = _uncertainty / _value;
                auto term2 = rhs.uncertainty() / rhs.value();
                // Now, update the value and uncertainty
                _value /= rhs.value();
                _uncertainty = _value * std::hypot(term1, term2);
                return *this;
            }

            friend quantity operator/(quantity lhs, const quantity &rhs)
            {
                lhs /= rhs;
                return lhs;
            }

            friend std::ostream& operator<<(std::ostream& out, const quantity& qty)
            {
                out << qty.value() << " " << qty.uncertainty();
                return out;
            }

        private:
            double _value, _uncertainty;
    };

    // T requirements: has .value(), .quantity() of numeric type
    template <typename T>
    T average(std::vector<T> list)
    {
        typename T::value_type weighted_sum(0.), sum_of_weights(0.);
        for (auto &qty : list)
        {
            if (qty.uncertainty() == 0.)
                return qty;
            auto weight = 1. / (qty.uncertainty() * qty.uncertainty());
            weighted_sum += weight * qty.value();
            sum_of_weights += weight;
        }
        return T(weighted_sum / sum_of_weights, 1. / std::sqrt(sum_of_weights));
    }

}

#endif
