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

#include <operations/components/independents/concrete/trace-writers/SeismicTraceWriter.hpp>


using namespace std;
using namespace operations::components;
using namespace operations::common;
using namespace operations::dataunits;

void SeismicTraceWriter::RecordTrace(uint time_step) {
    int trace_size = mTraceNumber;

    int wnx = this->mpGridBox->GetWindowAxis()->GetXAxis().GetActualAxisSize();
    int wnz_wnx = this->mpGridBox->GetWindowAxis()->GetZAxis().GetActualAxisSize() * wnx;

    int std_offset = (mpParameters->GetBoundaryLength() + mpParameters->GetHalfLength()) * wnx;
    float current_time = (time_step - 1) * this->mpGridBox->GetDT();
    uint trace_step = uint(current_time / mTraceSampling);
    if (time_step > 1) {
        float previous_time = (time_step - 2) * this->mpGridBox->GetDT();
        uint previous_trace_step = uint(previous_time / mTraceSampling);
        if (previous_trace_step == trace_step) {
            return;
        }
    }
    trace_step = min(trace_step, mSampleNumber - 1);
    auto positions_y = this->mpDPositionsY.GetNativePointer();
    auto positions_x = this->mpDPositionsX.GetNativePointer();
    auto values = this->mpDTraces.GetNativePointer();
    float *pressure = this->mpGridBox->Get(WAVE | GB_PRSS | CURR | DIR_Z)->GetNativePointer();
    for (int i = 0; i < trace_size; i++) {
        int offset = positions_y[i] * wnz_wnx + std_offset + positions_x[i];
        values[(trace_step) * trace_size + i] = pressure[offset];
    }
}
