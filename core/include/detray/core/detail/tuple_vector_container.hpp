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
#include <utility>

namespace detray {

template <template <typename...> class tuple_t,
          template <typename...> class vector_t, typename id_type,
          typename... Ts>
struct tuple_vector_container
    : public base_container<tuple_t, id_type, vector_t<Ts>...> {

    public:
    template <typename T>
    using vector_type = vector_t<T>;

    template <typename... Args>
    using tuple_type = tuple_t<Args...>;

    using base_type = base_container<tuple_type, id_type, vector_type<Ts>...>;

    using container_type = typename base_type::container_type;

    using container_data_type = tuple_type<vecmem::data::vector_view<Ts>...>;

    tuple_vector_container() = delete;

    DETRAY_HOST
    tuple_vector_container(vecmem::memory_resource &resource)
        : _container(vector_type<Ts>{&resource}...) {}

    template <
        typename container_data_t,
        std::enable_if_t<
            !std::is_base_of_v<vecmem::memory_resource, container_data_t> &&
                !std::is_same_v<tuple_vector_container, container_data_t>,
            bool> = true>
    DETRAY_DEVICE tuple_vector_container(container_data_t &container_data)
        : _container(initialize_device_vectors(
              container_data, std::make_index_sequence<sizeof...(Ts)>{})) {}

    template <typename container_data_t, std::size_t... ints>
    DETRAY_DEVICE container_type
    initialize_device_vectors(container_data_t &container_data,
                              std::index_sequence<ints...> /*seq*/) {
        return vtuple::make_tuple(
            vector_type<Ts>(detail::get<ints>(container_data._data))...);
    }

    template <std::size_t... ints>
    DETRAY_HOST container_data_type
    initialize_data_vectors(std::index_sequence<ints...> /*seq*/) {
        return detail::make_tuple<tuple_type>(vecmem::data::vector_view<Ts>(
            vecmem::get_data(detail::get<ints>(_container)))...);
    }

    template <id_type ID, typename... e_types>
    DETRAY_HOST auto &add_elements(e_types &&... elements) noexcept(false) {

        auto &group = detail::get<ID>(_container);

        return group.emplace_back(std::forward<e_types>(elements)...);
    }

    template <std::size_t current_id = 0, typename T>
    DETRAY_HOST inline void add_vector(vector_type<T> &vector) noexcept(false) {

        auto &group = detail::get<current_id>(_container);

        if constexpr (std::is_same_v<decltype(vector), decltype(group)>) {

            group.reserve(group.size() + vector.size());
            group.insert(group.end(), vector.begin(), vector.end());
        }

        // Next mask type
        if constexpr (current_id < sizeof...(Ts) - 1) {
            return add_vector<current_id + 1>(vector);
        }
    }

    template <std::size_t current_id = 0, typename T>
    DETRAY_HOST inline void add_vector(vector_type<T> &&vector) noexcept(
        false) {
        // Get the mask group that will be updated
        auto &group = detail::get<current_id>(_container);

        if constexpr (std::is_same_v<decltype(vector), decltype(group)>) {
            // Reserve memory and copy new masks
            group.reserve(group.size() + vector.size());
            group.insert(group.end(), std::make_move_iterator(vector.begin()),
                         std::make_move_iterator(vector.end()));
        }

        // Next mask type
        if constexpr (current_id < sizeof...(Ts) - 1) {
            return add_vector<current_id + 1>(vector);
        }
    }

    template <std::size_t current_id = 0>
    DETRAY_HOST inline void append_container(tuple_vector_container &&other) {
        // Add masks to current group
        auto &group = detail::get<current_id>(other);
        add_vector(group);

        // Next mask type
        if constexpr (current_id < sizeof...(Ts) - 1) {
            return append_container<current_id + 1>(other);
        }
    }

    private:
    container_type _container;
};

}  // namespace detray