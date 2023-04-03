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

#ifndef OPERATIONS_LIB_COMPONENTS_FORWARD_COLLECTORS_TWO_PROPAGATION_HPP
#define OPERATIONS_LIB_COMPONENTS_FORWARD_COLLECTORS_TWO_PROPAGATION_HPP

#include <cstdlib>
#include <cstring>
#include <fstream>
#include <sstream>
#include <unistd.h>

#include <bs/base/memory/MemoryManager.hpp>

#include <operations/components/independents/concrete/forward-collectors/file-handler/file_handler.h>
#include <operations/components/dependents/concrete/memory-handlers/WaveFieldsMemoryHandler.hpp>
#include <operations/components/independents/primitive/ForwardCollector.hpp>
#include <operations/components/dependency/concrete/HasDependents.hpp>

// #ifdef BUILD_FOR_NVIDIA
#include <fstream>
#include "nvcomp/lz4.hpp"
#include "nvcomp.hpp"
#include "nvcomp/nvcompManagerFactory.hpp"
#include "cusz.h"
#include "cuszapi.hh"
#include "veloc.hpp"
#include "veloc.h"

#define checkCuda(ans) { checkCudaFunc((ans), __FILE__, __LINE__); }
inline void checkCudaFunc(cudaError_t code, const char *file, int line, bool abort=true) {
   if (code != cudaSuccess) {
      fprintf(stderr,"========= GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
      if (abort) exit(code);
   }
}

// #endif

namespace operations {
    namespace components {

        class TwoPropagation : public ForwardCollector,
                               public dependency::HasDependents {
        public:
            explicit TwoPropagation(bs::base::configurations::ConfigurationMap *apConfigurationMap, const std::string &velocClient);

            ~TwoPropagation() override;

            void SetComputationParameters(common::ComputationParameters *apParameters) override;

            void SetGridBox(dataunits::GridBox *apGridBox) override;

            void SetDependentComponents(
                    operations::helpers::ComponentsMap<DependentComponent> *apDependentComponentsMap) override;

            void FetchForward() override;

            void SaveForward() override;

            void ResetGrid(bool aIsForwardRun) override;

            dataunits::GridBox *GetForwardGrid() override;

            void AcquireConfiguration() override;


        private:
            common::ComputationParameters *mpParameters = nullptr;

            dataunits::GridBox *mpMainGridBox = nullptr;

            dataunits::GridBox *mpInternalGridBox = nullptr;

            WaveFieldsMemoryHandler *mpWaveFieldsMemoryHandler = nullptr;

            dataunits::FrameBuffer<float> *mpForwardPressure = nullptr;

            float *mpForwardPressureHostMemory = nullptr;

            float *mpTempPrev = nullptr;

            float *mpTempCurr = nullptr;

            float *mpTempNext = nullptr;

            bool mIsMemoryFit;

            uint mTimeCounter;

            unsigned long long mMaxNT;

            unsigned long long mMaxDeviceNT;

            unsigned int mpMaxNTRatio;

            std::string mWritePath;

            bool mIsCompression;

            /* ZFP Properties. */

            int mZFP_Parallel;

            bool mZFP_IsRelative;

            float mZFP_Tolerance;

            // vector<uint64_t> cmpr_times;
            std::vector<size_t> iter_counts;
            std::vector<uint64_t> func_start_times;
            std::vector<uint64_t> func_end_times;
            std::vector<uint64_t> nvcomp_times;
            std::vector<uint64_t> nvcomp_sizes;
            std::vector<uint64_t> data_sizes;

            // std::vector<uint64_t> ckpt_times;
            // std::vector<uint64_t> func_times;
            // std::vector<uint64_t> d2d_times;
            // std::vector<uint64_t> d2h_times;
            // std::vector<uint64_t> h2f_times;
            // std::vector<uint64_t> data_sizes;
            uint64_t prev_ckpt_time = 0;
        };
    }//namespace components
}//namespace operations

#endif //OPERATIONS_LIB_COMPONENTS_FORWARD_COLLECTORS_TWO_PROPAGATION_HPP
