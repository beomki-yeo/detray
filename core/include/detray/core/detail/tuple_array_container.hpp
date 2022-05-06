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
          template <typename, std::size_t> class array_t, typename ID,
          typename... Ts>
class tuple_array_container
    : public base_container<tuple_t, ID,
                            array_t<typename Ts::object_type, Ts::N>...> {

    public:
    using base_type =
        base_container<tuple_t, ID,
                       array_t<typename Ts::object_type, Ts::N>...>;
    using base_type::base_type;

    template <typename... Args>
    using tuple_type = typename base_type::template tuple_type<Args...>;

    template <typename T, std::size_t I>
    using array_type = array_t<T, I>;

    using container_type = typename base_type::container_type;

    /// Host container constructor
    DETRAY_HOST
    tuple_array_container(vecmem::memory_resource& resource)
        : base_type(
              container_type(initialize_host_arrays<typename Ts::object_type>(
                  resource, std::make_index_sequence<Ts::N>{})...)) {}

    template <typename T, std::size_t... ints>
    DETRAY_HOST array_type<T, sizeof...(ints)> initialize_host_arrays(
        vecmem::memory_resource& resource,
        std::index_sequence<ints...> /*seq*/) {

        array_type<vecmem::memory_resource*, sizeof...(ints)> resources;

        std::fill(resources.begin(), resources.end(), &resource);

        return array_type<T, sizeof...(ints)>{T(resources[ints])...};
    }

    /// Device container constructor
    template <typename container_data_t,
              std::enable_if_t<
                  !std::is_base_of_v<vecmem::memory_resource, container_data_t>,
                  bool> = true>
    DETRAY_HOST_DEVICE tuple_array_container(
        const container_data_t& container_data)
        : base_type(initialize_device_arrays<Ts::object_type>(
              container_data, std::make_index_sequence<Ts::N>{})...) {}

    template <typename T, std::size_t... ints, typename container_data_t>
    DETRAY_HOST_DEVICE array_type<T, sizeof...(ints)> initialize_device_arrays(
        const container_data_t& container_data,
        std::index_sequence<ints...> /*seq*/) {
        return array_type<T, sizeof...(ints)>{container_data.data[ints]...};
    }
};

}  // namespace detray