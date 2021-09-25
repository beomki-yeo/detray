/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2020 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

#include <vecmem/memory/memory_resource.hpp>

#include "definitions/invalid_values.hpp"

namespace detray {

/** A two-dimensional grid for object storage
 *
 * @tparam populator_type  is a prescription what to do when a bin gets
 *pupulated, it broadcasts also the value type
 * @tparam tparam axis_p0_type the type of the first axis
 * @tparam tparam axis_p1_type the type of the second axis
 * @tparam serialzier_type  type of the serializer to the storage represenations
 *
 **/
template <typename populator_type, typename axis_p0_type, typename axis_p1_type,
          typename serializer_type,
          template <typename, unsigned int> class array_type = darray,
          template <typename...> class tuple_type = dtuple,
          template <typename> class vector_type = dvector>
class grid2 {

    public:
    using populator_t = populator_type;
    using axis_p0_t = axis_p0_type;
    using axis_p1_t = axis_p1_type;

    using bare_value = typename populator_type::bare_value;
    using serialized_storage = typename populator_t::serialized_storage;
    using point2 = __plugin::point2;

    template <typename neighbor_t>
    using neighborhood = array_type<array_type<neighbor_t, 2>, 2>;

    static constexpr array_type<dindex, 2> hermit1 = {0u, 0u};
    static constexpr neighborhood<dindex> hermit2 = {hermit1, hermit1};

    /** Constructor from axes - copy semantics
     *
     * @param axis_p0 is the axis in the first coordinate
     * @param axis_p1 is the axis in the second coordinate
     *
     **/
    grid2(const axis_p0_type &axis_p0, const axis_p1_type &axis_p1,
          vecmem::memory_resource &mr,
          const bare_value kInvalid = invalid_value<bare_value>())
        : _axis_p0(axis_p0),
          _axis_p1(axis_p1),
          _data_serialized(&mr),
          _populator(kInvalid) {
        _data_serialized = serialized_storage(_axis_p0.bins() * _axis_p1.bins(),
                                              _populator.init());
    }

    /** Constructor from axes - move semantics
     *
     * @param axis_p0 is the axis in the first coordinate
     * @param axis_p1 is the axis in the second coordinate
     *
     **/
    grid2(axis_p0_type &&axis_p0, axis_p1_type &&axis_p1,
          vecmem::memory_resource &mr,
          const bare_value kInvalid = invalid_value<bare_value>())
        : _axis_p0(std::move(axis_p0)),
          _axis_p1(std::move(axis_p1)),
          _populator(kInvalid),
          _data_serialized(&mr) {
        _data_serialized = serialized_storage(_axis_p0.bins() * _axis_p1.bins(),
                                              _populator.init());
    }

    /** Allow for grid shift, when using a centralized store and indices
     *
     * @param offset is the applied offset shift
     *
     **/
    void shift(const typename populator_type::bare_value &offset) {
        std::for_each(_data_serialized.begin(), _data_serialized.end(),
                      [&](auto &ds) { _populator.shift(ds, offset); });
    }

    /** Fill/populate operation
     *
     * @tparam point2_type the 2D local point type
     *
     * @param p2 the point in p2 local frame
     * @param fvalue is a single fill value to be filled
     *
     **/
    template <typename point2_type>
    void populate(const point2_type &p2,
                  typename populator_type::bare_value &&fvalue) {
        auto sbin = _serializer.template serialize<axis_p0_type, axis_p1_type>(
            _axis_p0, _axis_p1, _axis_p0.bin(p2[0]), _axis_p1.bin(p2[1]));
        _populator(_data_serialized[sbin], std::move(fvalue));
    }

    /** Fill/populate operation - with bin entry
     *
     * @param bin The two-dimensional bin2
     * @param fvalue is a single fill value to be filled
     *
     **/
    void populate(dindex bin0, dindex bin1,
                  typename populator_type::bare_value &&fvalue) {
        auto sbin = _serializer.template serialize<axis_p0_type, axis_p1_type>(
            _axis_p0, _axis_p1, bin0, bin1);
        _populator(_data_serialized[sbin], std::move(fvalue));
    }

    /** Return the value of a single bin - with direct bin acess
     *
     * @param bin0 the index of bin 0
     * @param bin1 the index of bin 1
     *
     * @return the const reference to the value in this bin
     **/
    const auto &bin(dindex bin0, dindex bin1) const {
        return _data_serialized[_serializer.template serialize<
            axis_p0_type, axis_p1_type>(_axis_p0, _axis_p1, bin0, bin1)];
    }

    /** Return the value of a single bin
     *
     * @param p2 is point in the local frame
     *
     * @return the const reference to the value in this bin
     **/
    const auto &bin(const point2 &p2) const {
        return _data_serialized[_serializer.template serialize<axis_p0_type,
                                                               axis_p1_type>(
            _axis_p0, _axis_p1, _axis_p0.bin(p2[0]), _axis_p1.bin(p2[1]))];
    }

    /** Return the value of a single bin - non-const access
     *
     * @param p2 is point in the local frame
     *
     * @return the const reference to the value in this bin
     **/
    auto &bin(const point2 &p2) {
        return _data_serialized[_serializer.template serialize<axis_p0_type,
                                                               axis_p1_type>(
            _axis_p0, _axis_p1, _axis_p0.bin(p2[0]), _axis_p1.bin(p2[1]))];
    }

    /** Return a zone around a single bin, either with binned or scalar
     *neighborhood
     *
     * The zone is done with a neighborhood around the bin which is defined by
     *p2
     *
     * @param p2 is point in the local frame
     * @param nhood is the binned/scalar neighborhood
     * @param sort is a directive whether to sort or not
     *
     * @return the sequence of values
     **/
    template <typename neighbor_t>
    vector_type<typename populator_type::bare_value> zone_t(
        const point2 &p2, const neighborhood<neighbor_t> &nhood,
        bool sort) const {
        auto zone0 = _axis_p0.zone(p2[0], nhood[0]);
        auto zone1 = _axis_p1.zone(p2[1], nhood[1]);

        vector_type<typename populator_type::bare_value> zone;

        // Specialization for bare value equal to store value
        if constexpr (std::is_same_v<typename populator_type::bare_value,
                                     typename populator_type::store_value>) {
            unsigned int iz = 0;
            zone = vector_type<typename populator_type::bare_value>(
                zone0.size() * zone1.size(), {});
            for (const auto z1 : zone1) {
                for (const auto z0 : zone0) {
                    auto sbin =
                        _serializer
                            .template serialize<axis_p0_type, axis_p1_type>(
                                _axis_p0, _axis_p1, z0, z1);
                    zone[iz++] = _data_serialized[sbin];
                }
            }
        } else {
            zone.reserve(10);
            for (const auto z1 : zone1) {
                for (const auto z0 : zone0) {
                    auto sbin =
                        _serializer
                            .template serialize<axis_p0_type, axis_p1_type>(
                                _axis_p0, _axis_p1, z0, z1);
                    auto bin_data = _data_serialized[sbin];
                    auto bin_content = _populator.sequence(bin_data);
                    zone.insert(zone.end(), bin_content.begin(),
                                bin_content.end());
                }
            }
        }

        if (sort) {
            std::sort(zone.begin(), zone.end());
        }
        return zone;
    }

    /** Return a zone around a single bin, either with binned neighborhood
     *
     * The zone is done with a neighborhood around the bin which is defined by
     *p2
     *
     * @param p2 is point in the local frame
     * @param nhood is the binned neighborhood
     * @param sort is a directive whether to sort or not
     *
     * @return the sequence of values
     **/
    vector_type<typename populator_type::bare_value> zone(
        const point2 &p2, const neighborhood<dindex> &nhood = hermit2,
        bool sort = false) const {
        return zone_t<dindex>(p2, nhood, sort);
    }

    /** Return a zone around a single bin, either with scalar neighborhood
     *
     * The zone is done with a neighborhood around the bin which is defined by
     *p2
     *
     * @param p2 is point in the local frame
     * @param nhood is the binned neighborhood
     * @param sort is a directive whether to sort or not
     *
     * @return the sequence of values
     **/
    vector_type<typename populator_type::bare_value> zone(
        const point2 &p2, const neighborhood<scalar> &nhood,
        bool sort = false) const {
        return zone_t<scalar>(p2, nhood, sort);
    }

    /** Const access to axis p0  */
    const axis_p0_type &axis_p0() const { return _axis_p0; }

    /** Const access to axis p1 */
    const axis_p1_type &axis_p1() const { return _axis_p1; }

    /* Copy of axes in a tuple */
    tuple_type<axis_p0_type, axis_p1_type> axes() const {
        return std::tie(_axis_p0, _axis_p1);
    }

    /** Const acess to the serializer */
    const serializer_type &serializer() const { return _serializer; }

    /** Const acess to the polulator */
    const populator_type &populator() const { return _populator; }

    serialized_storage &data() { return _data_serialized; }

    private:
    serialized_storage _data_serialized;
    axis_p0_type _axis_p0;
    axis_p1_type _axis_p1;
    serializer_type _serializer;
    populator_type _populator;
};

/** A two-dimensional grid data for gpu device usage
 *
 * @tparam grid_t type of grid
 **/
template <typename grid_t>
struct grid2_data {

    using populator_t = typename grid_t::populator_t;
    using vector_view_t = typename populator_t::vector_view_t;
    using axis_p0_t = typename grid_t::axis_p0_t;
    using axis_p1_t = typename grid_t::axis_p1_t;

    grid2_data(grid_t &grid, vecmem::memory_resource *resource = nullptr)
        : _axis_p0(grid.axis_p0()),
          _axis_p1(grid.axis_p1()),
          _data_serialized(populator_t::get_data(grid.data(), resource)) {}
    
    vector_view_t _data_serialized;
    const axis_p0_t _axis_p0;
    const axis_p1_t _axis_p1;
};

/** A two-dimensional grid buffer for gpu device usage
 *
 * @tparam grid_t type of grid
 **/
template <typename grid_t>
struct grid2_buffer {

    using axis_p0_t = typename grid_t::axis_p0_t;
    using axis_p1_t = typename grid_t::axis_p1_t;
    using populator_t = typename grid_t::populator_t;
    using vector_buffer_t = typename populator_t::vector_buffer_t;
    using buffer_size_t = typename populator_t::buffer_size_t;    
    
    grid2_buffer(const axis_p0_t& axis_p0,
		 const axis_p0_t& axis_p1,
		 buffer_size_t& sizes,
		 buffer_size_t& capacities,
		 vecmem::memory_resource& resource)
        : _axis_p0(axis_p0),
          _axis_p1(axis_p1),
          _buffer(sizes, capacities, resource) {}
    
    vector_buffer_t _buffer;
    const axis_p0_t _axis_p0;
    const axis_p1_t _axis_p1;
};
        
}  // namespace detray
