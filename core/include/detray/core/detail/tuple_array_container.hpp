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

    static constexpr std::size_t m_tuple_size = sizeof...(Ts);

    template <typename T, std::size_t I>
    using array_type = array_t<T, I>;

    using container_type = typename base_type::container_type;
    using container_data_type =
        tuple_type<array_type<typename Ts::data_type, Ts::N>...>;
    using container_view_type =
        tuple_type<array_type<typename Ts::view_type, Ts::N>...>;

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
        : base_type(initialize_device_arrays(
              container_data, std::make_index_sequence<sizeof...(Ts)>{})) {}

    template <typename container_data_t, std::size_t... ints>
    DETRAY_HOST_DEVICE container_type
    initialize_device_arrays(const container_data_t& container_data,
                             std::index_sequence<ints...> /*seq*/) {

        return container_type(to_device<Ts::object_type>(
            vtuple::get<ints>(container_data.m_data))...);
    }

    template <typename object_t, typename data_t, std::size_t N>
    DETRAY_HOST_DEVICE array_type<object_t, N> to_device(
        const array_type<data_t, N>& A) {
        to_device<object_t>(A, std::make_index_sequence<N>{});
    }

    template <typename object_t, typename data_t, std::size_t N,
              std::size_t... ints>
    DETRAY_HOST_DEVICE array_type<object_t, N> to_device(
        const array_type<data_t, N>& A, std::index_sequence<ints...> /*seq*/) {
        return array_type<object_t, N>{object_t(A[ints])...};
    }

    template <std::size_t... ints>
    DETRAY_HOST container_data_type
    initialize_data_arrays(vecmem::memory_resource& resource,
                           std::index_sequence<ints...> /*seq*/) {

        return container_data_type{array_type<typename Ts::data_type, Ts::N>(
            to_data<typename Ts::data_type>(
                detail::get<ints>(this->m_container), resource))...};
    }

    template <typename data_t, typename object_t, std::size_t N>
    DETRAY_HOST array_type<data_t, N> to_data(
        array_type<object_t, N>& A, vecmem::memory_resource& resource) {
        return to_data<data_t>(A, resource, std::make_index_sequence<N>{});
    }

    template <typename data_t, typename object_t, std::size_t N,
              std::size_t... ints>
    DETRAY_HOST array_type<data_t, N> to_data(
        array_type<object_t, N>& A, vecmem::memory_resource& resource,
        std::index_sequence<ints...> /*seq*/) {

        (void)resource;
        return array_type<data_t, N>{vecmem::get_data(A[ints])...};
    }
};

template <typename container_t>
struct tuple_array_container_data {

    template <typename... Args>
    using tuple_type = typename container_t::template tuple_type<Args...>;

    using container_view_type = typename container_t::container_view_type;

    tuple_array_container_data(container_t& container,
                               vecmem::memory_resource& resource)
        : m_data(container.initialize_data_arrays(
              resource, std::make_index_sequence<container.size()>{})) {}

    template <std::size_t... ints>
    container_view_type initialize_view_arrays(
        std::index_sequence<ints...> /*seq*/) {
        return container_view_type(
            detail::make_tuple<tuple_type>(detail::get<ints>(m_data)...));
    }

    typename container_t::container_data_type m_data;
};

template <typename container_t>
struct tuple_array_container_view {

    template <typename... Args>
    using tuple_type = typename container_t::template tuple_type<Args...>;

    tuple_array_container_view(
        tuple_array_container_data<container_t>& container_data)
        : m_data(container_data.initialize_view_arrays(
              std::make_index_sequence<container_t::m_tuple_size>{})) {}

    typename container_t::container_view_type m_data;
};

template <template <typename...> class tuple_t,
          template <typename, std::size_t> class array_t, typename ID,
          typename... Ts>
inline tuple_array_container_data<
    tuple_array_container<tuple_t, array_t, ID, Ts...>>
get_data(tuple_array_container<tuple_t, array_t, ID, Ts...>& container,
         vecmem::memory_resource& resource) {
    return {container, resource};
}

}  // namespace detray