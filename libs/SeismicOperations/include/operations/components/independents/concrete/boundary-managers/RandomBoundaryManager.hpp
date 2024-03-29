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

#ifndef OPERATIONS_LIB_COMPONENTS_BOUNDARY_MANAGERS_RANDOM_BOUNDARY_MANAGER_HPP
#define OPERATIONS_LIB_COMPONENTS_BOUNDARY_MANAGERS_RANDOM_BOUNDARY_MANAGER_HPP

#include <vector>

#include <operations/components/independents/concrete/boundary-managers/extensions/Extension.hpp>
#include <operations/components/independents/primitive/BoundaryManager.hpp>
#include <operations/components/dependency/concrete/HasNoDependents.hpp>


namespace operations {
    namespace components {

        /**
         * @note
         * Random Boundary Manager is based on the following paper:
         * https://library.seg.org/doi/abs/10.1190/1.3255432.
         */
        class RandomBoundaryManager : public BoundaryManager,
                                      public dependency::HasNoDependents {
        public:
            explicit RandomBoundaryManager(bs::base::configurations::ConfigurationMap *apConfigurationMap);

            ~RandomBoundaryManager() override;

            void ApplyBoundary(uint kernel_id) override;

            void ExtendModel() override;

            void ReExtendModel() override;

            void SetComputationParameters(common::ComputationParameters *apParameters) override;

            void SetGridBox(dataunits::GridBox *apGridBox) override;

            void AdjustModelForBackward() override;

            void AcquireConfiguration() override;

        private:
            void InitializeExtensions();

        private:
            common::ComputationParameters *mpParameters = nullptr;

            dataunits::GridBox *mpGridBox = nullptr;

            std::vector<addons::Extension *> mvExtensions;

            int mGrainSideLength;
        };

    }//namespace components
}//namespace operations

#endif // OPERATIONS_LIB_COMPONENTS_BOUNDARY_MANAGERS_RANDOM_BOUNDARY_MANAGER_HPP
