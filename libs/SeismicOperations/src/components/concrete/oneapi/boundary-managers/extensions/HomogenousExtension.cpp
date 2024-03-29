/**
 * Copyright (C) 2021 by Brightskies inc
 *
 * This file is part of SeismicToolbox.
 *
 * SeismicToolbox is free software: you can redistribute it and/or modify it
 * under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * SeismicToolbox is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with GEDLIB. If not, see <http://www.gnu.org/licenses/>.
 */

#include <vector>

#include <operations/components/independents/concrete/boundary-managers/extensions/HomogenousExtension.hpp>
#include <bs/base/backend/Backend.hpp>

using namespace std;
using namespace cl::sycl;
using namespace bs::base::backend;
using namespace operations::components;
using namespace operations::components::addons;
using namespace operations::dataunits;

HomogenousExtension::HomogenousExtension(bool use_top_layer) {
    this->mUseTop = use_top_layer;
}

void HomogenousExtension::VelocityExtensionHelper(float *property_array,
                                                  int start_x, int start_y, int start_z,
                                                  int end_x, int end_y, int end_z,
                                                  int nx, int ny, int nz,
                                                  uint boundary_length) {
    /*!
     * change the values of velocities at boundaries (HALF_LENGTH excluded) to
     * zeros the start for x , y and z is at HALF_LENGTH and the end is at (nx -
     * HALF_LENGTH) or (ny - HALF_LENGTH) or (nz- HALF_LENGTH)
     */
    int nz_nx = nx * nz;

    // In case of 2D
    if (ny == 1) {
        end_y = 1;
        start_y = 0;
    } else {
        // general case for 3D
        /*!putting the nearest property_array adjacent to the boundary as the value
         * for all velocities at the boundaries for y and with all x and z */
        vector<uint> gridDims{(uint) end_x - start_x, boundary_length, (uint) end_z - start_z};
        vector<uint> blockDims{1, 1, 1};
        auto configs = Backend::GetInstance()->CreateKernelConfiguration(gridDims, blockDims);
        Backend::GetInstance()->GetDeviceQueue()->submit([&](handler &cgh) {
            auto global_nd_range = nd_range<3>(configs.mGridDimensions, configs.mBlockDimensions);

            cgh.parallel_for(global_nd_range, [=](nd_item<3> it) {
                int column = it.get_global_id(0) + start_x;
                int depth = it.get_global_id(1);
                int row = it.get_global_id(2) + start_z;

                /*!for values from y = HALF_LENGTH TO y= HALF_LENGTH +BOUND_LENGTH*/
                int p_idx = (depth + start_y) * nz_nx + row * nx + column;
                int p2_idx =
                        (boundary_length + start_y) * nz_nx + row * nx + column;
                property_array[p_idx] = property_array[p2_idx];

                /*!for values from y = ny-HALF_LENGTH TO y =
                 * ny-HALF_LENGTH-BOUND_LENGTH*/
                p_idx = (end_y - 1 - depth) * nz_nx + row * nx + column;
                p2_idx = (end_y - 1 - boundary_length) * nz_nx + row * nx + column;
                property_array[p_idx] = property_array[p2_idx];
            });
        });
        Backend::GetInstance()->GetDeviceQueue()->wait();
    }

    /*!putting the nearest property_array adjacent to the boundary as the value
     * for all velocities at the boundaries for x and with all z and y */
    vector<uint> gridDims{boundary_length, (uint) end_y - start_y, (uint) end_z - start_z};
    vector<uint> blockDims{1, 1, 1};
    auto configs = Backend::GetInstance()->CreateKernelConfiguration(gridDims, blockDims);
    Backend::GetInstance()->GetDeviceQueue()->submit([&](handler &cgh) {
        auto global_nd_range = nd_range<3>(configs.mGridDimensions, configs.mBlockDimensions);

        cgh.parallel_for(global_nd_range, [=](nd_item<3> it) {
            int column = it.get_global_id(0);
            int depth = it.get_global_id(1) + start_y;
            int row = it.get_global_id(2) + start_z;

            /*!for values from x = HALF_LENGTH TO x= HALF_LENGTH +BOUND_LENGTH*/
            int p_idx = depth * nz_nx + row * nx + column + start_x;
            int p2_idx = depth * nz_nx + row * nx + boundary_length + start_x;
            property_array[p_idx] = property_array[p2_idx];

            /*!for values from x = nx-HALF_LENGTH TO x =
             * nx-HALF_LENGTH-BOUND_LENGTH*/
            p_idx = depth * nz_nx + row * nx + (end_x - 1 - column);
            p2_idx = depth * nz_nx + row * nx + (end_x - 1 - boundary_length);
            property_array[p_idx] = property_array[p2_idx];
        });
    });
    Backend::GetInstance()->GetDeviceQueue()->wait();

    bool extend_top = this->mUseTop;
    /*!putting the nearest property_array adjacent to the boundary as the value
     * for all velocities at the boundaries for z and with all x and y */
    gridDims = {(uint) end_x - start_x, (uint) end_y - start_y, boundary_length};
    configs = Backend::GetInstance()->CreateKernelConfiguration(gridDims, blockDims);
    Backend::GetInstance()->GetDeviceQueue()->submit([&](handler &cgh) {
        auto global_nd_range = nd_range<3>(configs.mGridDimensions, configs.mBlockDimensions);

        cgh.parallel_for(global_nd_range, [=](nd_item<3> it) {
            int column = it.get_global_id(0) + start_x;
            int depth = it.get_global_id(1) + start_y;
            int row = it.get_global_id(2);

            /*!for values from z = HALF_LENGTH TO z = HALF_LENGTH +BOUND_LENGTH */
            int p_idx = depth * nz_nx + (start_z + row) * nx + column;
            int p2_idx =
                    depth * nz_nx + (start_z + boundary_length) * nx + column;
            if (extend_top) {
                property_array[p_idx] = property_array[p2_idx];
            }

            /*!for values from z = nz-HALF_LENGTH TO z =
             * nz-HALF_LENGTH-BOUND_LENGTH*/
            p_idx = depth * nz_nx + (end_z - 1 - row) * nx + column;
            p2_idx = depth * nz_nx + (end_z - 1 - boundary_length) * nx + column;
            property_array[p_idx] = property_array[p2_idx];
        });
    });
    Backend::GetInstance()->GetDeviceQueue()->wait();
}

void HomogenousExtension::TopLayerExtensionHelper(float *property_array,
                                                  int start_x, int start_y, int start_z,
                                                  int end_x, int end_y, int end_z,
                                                  int nx, int ny, int nz, uint boundary_length) {
    if (this->mUseTop) {
        int nz_nx = nx * nz;
        if (ny == 1) {
            start_y = 0;
            end_y = 1;
        }
        /*!putting the nearest property_array adjacent to the boundary as the value
         * for all velocities at the boundaries for z and with all x and y */
        vector<uint> gridDims{(uint) end_x - start_x, (uint) end_y - start_y, boundary_length};
        vector<uint> blockDims{1, 1, 1};
        auto configs = Backend::GetInstance()->CreateKernelConfiguration(gridDims, blockDims);
        Backend::GetInstance()->GetDeviceQueue()->submit([&](handler &cgh) {
            auto global_nd_range = nd_range<3>(configs.mGridDimensions, configs.mBlockDimensions);

            cgh.parallel_for(global_nd_range, [=](nd_item<3> it) {
                int column = it.get_global_id(0) + start_x;
                int depth = it.get_global_id(1) + start_y;
                int row = it.get_global_id(2);

                /*!for values from z = HALF_LENGTH TO z = HALF_LENGTH +BOUND_LENGTH */
                int p_idx = depth * nz_nx + (start_z + row) * nx + column;
                int p2_idx =
                        depth * nz_nx + (start_z + boundary_length) * nx + column;
                property_array[p_idx] = property_array[p2_idx];
            });
        });
        Backend::GetInstance()->GetDeviceQueue()->wait();
    }
}

void HomogenousExtension::TopLayerRemoverHelper(float *property_array,
                                                int start_x, int start_y, int start_z,
                                                int end_x, int end_y, int end_z,
                                                int nx, int ny, int nz, uint boundary_length) {
    if (this->mUseTop) {
        if (ny == 1) {
            start_y = 0;
            end_y = 1;
        }
        int nz_nx = nx * nz;
        /*!putting the nearest property_array adjacent to the boundary as the value
         * for all velocities at the boundaries for z and with all x and y */
        vector<uint> gridDims{(uint) end_x - start_x, (uint) end_y - start_y, boundary_length};
        vector<uint> blockDims{1, 1, 1};
        auto configs = Backend::GetInstance()->CreateKernelConfiguration(gridDims, blockDims);
        Backend::GetInstance()->GetDeviceQueue()->submit([&](handler &cgh) {
            auto global_nd_range = nd_range<3>(configs.mGridDimensions, configs.mBlockDimensions);

            cgh.parallel_for(global_nd_range, [=](nd_item<3> it) {
                int column = it.get_global_id(0) + start_x;
                int depth = it.get_global_id(1) + start_y;
                int row = it.get_global_id(2);

                /*!for values from z = HALF_LENGTH TO z = HALF_LENGTH +BOUND_LENGTH */
                int p_idx = depth * nz_nx + (start_z + row) * nx + column;
                property_array[p_idx] = 0;
            });
        });
        Backend::GetInstance()->GetDeviceQueue()->wait();
    }
}
