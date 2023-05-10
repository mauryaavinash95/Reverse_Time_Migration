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

///
/// @brief This should contain the main function that drives the
/// Seismic Engine execution.
///

#include <bs/base/logger/concrete/LoggerSystem.hpp>
#include <bs/base/logger/concrete/FileLogger.hpp>
#include <bs/base/logger/concrete/ConsoleLogger.hpp>

#include <stbx/parsers/Parser.hpp>
#include <stbx/parsers/ArgumentsParser.hpp>
#include <stbx/generators/Generator.hpp>

using namespace std;
using namespace stbx::parsers;
using namespace stbx::generators;
using namespace stbx::writers;
using namespace operations::dataunits;
using namespace operations::engines;
using namespace bs::base::logger;
using namespace bs::timer::configurations;


int main(int argc, char *argv[]) {
    string parameter_file = WORKLOAD_PATH "/computation_parameters.json";
    string configuration_file = WORKLOAD_PATH "/engine_configuration.json";
    string callback_file = WORKLOAD_PATH "/callback_configuration.json";
    string system_file = WORKLOAD_PATH "/system_configuration.json";
    string veloc_config = WORKLOAD_PATH "/veloc_config.cfg";
    string write_path = WRITE_PATH;

    auto logger = LoggerSystem::GetInstance();

    /* Registering the needed loggers. */
    logger->RegisterLogger(new FileLogger(write_path + "/log.txt"));
    logger->RegisterLogger(new ConsoleLogger());

    /* Configuring the logger. */
    logger->ConfigureLoggers("Main logger", FILE_CONSOLE, DATE_TIME);
    logger->Info() << "Starting Seismic Engine..." << '\n';

    ArgumentsParser::Parse(parameter_file,
                           configuration_file,
                           callback_file,
                           system_file,
                           write_path,
                           argc, argv);

    auto parser = Parser::GetInstance();
    parser->RegisterFile(parameter_file);
    parser->RegisterFile(configuration_file);
    parser->RegisterFile(callback_file);
    parser->RegisterFile(system_file);
    parser->BuildMap();

    auto generator = new Generator(parser->GetMap());
    auto engine = generator->GenerateEngine(write_path);

    TimerManager::GetInstance()->Configure(generator->GenerateTimerConfiguration());

    auto agent = generator->GenerateAgent();
    agent->AssignEngine(engine);
    agent->AssignArgs(argc, argv);
    auto md = agent->Execute(veloc_config);

    delete engine;

    auto writer = generator->GenerateWriter();
    writer->AssignMigrationData(md);
    writer->Write(write_path);

    TimerManager::GetInstance()->Terminate(true);
    TimerManager::Kill();

    delete parser;
    delete generator;

    return EXIT_SUCCESS;
}
