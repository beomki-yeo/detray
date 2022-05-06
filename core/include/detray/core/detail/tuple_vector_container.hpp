/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Detray include(s)
#include "detray/core/detail/base_container.hpp"
#include "detray/definitions/indexing.hpp"
#include "detray/definitions/qualifiers.hpp"

// Vecmem include(s)
#include <vecmem/containers/device_vector.hpp>
#include <vecmem/memory/memory_resource.hpp>

// System include(s)
#include <type_traits>
#include <utility>

namespace detray {

template <template <typename...> class tuple_t,
          template <typename> class vector_t, typename id_type, typename... Ts>
class tuple_vector_container
    : public base_container<tuple_t, id_type, vector_t<Ts>...> {

    public:
    using base_type = base_container<tuple_t, id_type, vector_t<Ts>...>;
    using base_type::base_type;

    template <typename... Args>
    using tuple_type = typename base_type::template tuple_type<Args...>;

    template <typename T>
    using vector_type = vector_t<T>;

    using container_type = typename base_type::container_type;
    using container_data_type = tuple_type<vecmem::data::vector_view<Ts>...>;

    DETRAY_HOST
    tuple_vector_container(vecmem::memory_resource &resource)
        : base_type(container_type(vector_type<Ts>{&resource}...)) {}

    template <
        typename container_data_t,
        std::enable_if_t<
            !std::is_base_of_v<vecmem::memory_resource, container_data_t> &&
                !std::is_same_v<tuple_vector_container, container_data_t>,
            bool> = true>
    DETRAY_DEVICE tuple_vector_container(container_data_t &container_data)
        : base_type(initialize_device_vectors(
              container_data, std::make_index_sequence<sizeof...(Ts)>{})) {}

    template <typename container_data_t, std::size_t... ints>
    DETRAY_DEVICE container_type
    initialize_device_vectors(container_data_t &container_data,
                              std::index_sequence<ints...> /*seq*/) {
        return vtuple::make_tuple(
            vector_type<Ts>(detail::get<ints>(container_data.m_data))...);
    }

    template <std::size_t... ints>
    DETRAY_HOST container_data_type
    initialize_data_vectors(std::index_sequence<ints...> /*seq*/) {
        return detail::make_tuple<tuple_type>(vecmem::data::vector_view<Ts>(
            vecmem::get_data(detail::get<ints>(this->m_container)))...);
    }

    template <id_type ID, typename... Args>
    DETRAY_HOST auto &add_element(Args &&... args) noexcept(false) {

        auto &gr = detail::get<ID>(this->m_container);

        return gr.emplace_back(std::forward<Args>(args)...);
    }

    template <std::size_t current_id = 0, typename T>
    DETRAY_HOST inline void add_vector(vector_type<T> &vector) noexcept(false) {

        auto &gr = detail::get<current_id>(this->m_container);

        if constexpr (std::is_same_v<decltype(vector), decltype(gr)>) {

            gr.reserve(gr.size() + vector.size());
            gr.insert(gr.end(), vector.begin(), vector.end());
        }

        if constexpr (current_id < sizeof...(Ts) - 1) {
            return add_vector<current_id + 1>(vector);
        }
    }

    template <std::size_t current_id = 0, typename T>
    DETRAY_HOST inline void add_vector(vector_type<T> &&vector) noexcept(
        false) {
        auto &gr = detail::get<current_id>(this->m_container);

        if constexpr (std::is_same_v<decltype(vector), decltype(gr)>) {
            gr.reserve(gr.size() + vector.size());
            gr.insert(gr.end(), std::make_move_iterator(vector.begin()),
                      std::make_move_iterator(vector.end()));
        }

        if constexpr (current_id < sizeof...(Ts) - 1) {
            return add_vector<current_id + 1>(vector);
        }
    }

    template <std::size_t current_id = 0>
    DETRAY_HOST inline void append_container(tuple_vector_container &&other) {
        auto &gr = detail::get<current_id>(other);
        add_vector(gr);

        if constexpr (current_id < sizeof...(Ts) - 1) {
            return append_container<current_id + 1>(other);
        }
    }
};

template <typename tuple_vector_container_t>
struct tuple_vector_container_data {

    using container_data_type =
        typename tuple_vector_container_t::container_data_type;

    template <std::size_t... ints>
    DETRAY_HOST tuple_vector_container_data(tuple_vector_container_t &container,
                                            std::index_sequence<ints...> seq)
        : m_data(container.initialize_data_vectors(seq)) {}

    container_data_type m_data;
};

template <template <typename...> class tuple_t,
          template <typename...> class vector_t, typename id_type,
          typename... Ts>
inline tuple_vector_container_data<
    tuple_vector_container<tuple_t, vector_t, id_type, Ts...>>
get_data(tuple_vector_container<tuple_t, vector_t, id_type, Ts...> &container) {
    return tuple_vector_container_data<
        tuple_vector_container<tuple_t, vector_t, id_type, Ts...>>(
        container, std::make_index_sequence<sizeof...(Ts)>{});
}

}  // namespace detray