/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

#include "detray/definitions/detail/accessor.hpp"
#include "detray/definitions/qualifiers.hpp"

namespace detray {

template <template <typename...> class tuple_t, typename id_type,
          typename... Ts>
class base_container {

    public:
    template <typename... Args>
    using tuple_type = tuple_t<Args...>;
    using container_type = vtuple::tuple<Ts...>;

    DETRAY_HOST_DEVICE
    base_container(container_type &&container)
        : m_container(std::move(container)) {}

    template <id_type ID>
    DETRAY_HOST_DEVICE size_t size() const {
        return detail::get<ID>(m_container).size();
    }

    DETRAY_HOST_DEVICE constexpr std::size_t size() const {
        return sizeof...(Ts);
    }

    template <id_type ID>
    DETRAY_HOST_DEVICE bool empty() const {
        return detail::get<ID>(m_container).empty();
    }

    template <id_type ID>
    DETRAY_HOST_DEVICE constexpr auto &group() {
        return detail::get<ID>(m_container);
    }

    template <id_type ID>
    DETRAY_HOST_DEVICE constexpr const auto &group() const {
        return detail::get<ID>(m_container);
    }

    template <std::size_t ref_idx = 0>
    DETRAY_HOST_DEVICE static constexpr id_type to_id(const std::size_t index) {
        if (ref_idx == index) {
            // Produce a more helpful error than the usual tuple index error
            static_assert(
                ref_idx < sizeof...(Ts),
                "Index out of range: Please make sure that indices and type "
                "enums match the number of types in container.");
            return static_cast<id_type>(index);
        }
        if constexpr (ref_idx < sizeof...(Ts) - 1) {
            return to_id<ref_idx + 1>(index);
        }
        // This produces a compiler error when used in type unrolling code
        return static_cast<id_type>(sizeof...(Ts));
    }

    protected:
    container_type m_container;
};

}  // namespace detray
