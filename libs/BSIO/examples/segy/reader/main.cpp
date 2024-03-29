/**
 * Copyright (C) 2021 by Brightskies inc
 *
 * This file is part of BS I/O.
 *
 * BS I/O is free software: you can redistribute it and/or modify it
 * under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * BS I/O is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with GEDLIB. If not, see <http://www.gnu.org/licenses/>.
 */

#include <iostream>

#include <prerequisites/libraries/nlohmann/json.hpp>

#include <bs/base/configurations/concrete/JSONConfigurationMap.hpp>

#include <bs/io/streams/concrete/readers/SegyReader.hpp>
#include <bs/io/utils/displayers/Displayer.hpp>
#include <bs/io/utils/timer/ExecutionTimer.hpp>
#include <bs/io/configurations/MapKeys.h>


using namespace std;
using namespace bs::base::configurations;
using namespace bs::io::streams;
using namespace bs::io::dataunits;
using namespace bs::io::utils::timer;
using namespace bs::io::utils::displayers;


int main(int argc, char *argv[]) {
    nlohmann::json configuration_map;
    configuration_map[IO_K_PROPERTIES][IO_K_TEXT_HEADERS_ONLY] = false;
    configuration_map[IO_K_PROPERTIES][IO_K_TEXT_HEADERS_STORE] = false;

    std::vector<TraceHeaderKey> gather_keys = {TraceHeaderKey::FLDR};
    std::vector<std::pair<TraceHeaderKey, Gather::SortDirection>> sorting_keys;
    std::vector<std::string> paths = {DATA_PATH "/shots0601_0800.segy",
                                      DATA_PATH "/vel_z6.25m_x12.5m_exact.segy"};

    SegyReader r(new JSONConfigurationMap(configuration_map));
    r.AcquireConfiguration();
    r.Initialize(gather_keys, sorting_keys, paths);

    ExecutionTimer::Evaluate([&]() {
        r.ReadAll();
    }, true);

    Displayer::PrintTextHeader(r.GetTextHeader());
    if (r.HasExtendedTextHeader()) {
        Displayer::PrintTextHeader(r.GetExtendedTextHeader());
    }
}
