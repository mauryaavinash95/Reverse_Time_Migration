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

#include <sys/stat.h>

#include <bs/base/api/cpp/BSBase.hpp>
#include <bs/timer/api/cpp/BSTimer.hpp>

#include <operations/components/independents/concrete/forward-collectors/TwoPropagation.hpp>
#include <operations/configurations/MapKeys.h>
#include <operations/utils/compressor/Compressor.hpp>
#include <chrono>


using namespace std;
using namespace bs::timer;
using namespace bs::base::logger;
using namespace bs::base::memory;
using namespace operations::helpers;
using namespace operations::components;
using namespace operations::components::helpers;
using namespace operations::common;
using namespace operations::dataunits;
using namespace operations::utils::compressors;
using namespace nvcomp;

static float *initial_internalGridbox_curr = nullptr;

// TwoPropagation::TwoPropagation(bs::base::configurations::ConfigurationMap *apConfigurationMap, const std::string &velocConfig) {
TwoPropagation::TwoPropagation(bs::base::configurations::ConfigurationMap *apConfigurationMap) {
    this->mpConfigurationMap = apConfigurationMap;
    this->mpInternalGridBox = new GridBox();
    this->mpForwardPressure = nullptr;
    this->mIsMemoryFit = false;
    this->mTimeCounter = 0;
    this->mIsCompression = false;
    this->mZFP_Tolerance = 0.01f;
    this->mZFP_Parallel = true;
    this->mZFP_IsRelative = false;
    this->mMaxNT = 0;
    // #ifdef BUILD_FOR_NVIDIA
    // checkCuda(cudaStreamCreate(&stream));
    // #endif
}

TwoPropagation::~TwoPropagation() {
    LoggerSystem *Logger = LoggerSystem::GetInstance();
    // for(int i=0; i<data_gen_rate.size(); i++) {
    //     Logger->Info() << i << " = " << data_gen_rate[i] << "\n";
    // }

    std::ofstream logfile;
    uint64_t curr_time = std::chrono::duration_cast<std::chrono::seconds>(std::chrono::system_clock::now().time_since_epoch()).count();
    std::string filename = std::string(std::to_string(curr_time) + "-logs.csv");
    logfile.open (filename);
    int n = data_sizes.size();
    logfile << "D2H transfer intval, " << this->mMaxDeviceNT <<  "H2F trf intval, " << this->mMaxNT << "\n";
    logfile << "Ckpt No., Curr ckpt start, Prev ckpt end, Time diff, Datasize, NVComp time, Cmpr size, Cmpr ratio, Ckpt gen rate\n";
    for(int i=0; i<n; i++) {
        uint64_t prev_func_end_times = 0;
        if (i>0)
            prev_func_end_times = func_end_times[i-1];
        logfile << iter_counts[i] << ", "
            << func_start_times[i] << ", " 
            << func_end_times[i-1] << ", " 
            << (func_start_times[i]-prev_func_end_times) << ", " 
            << data_sizes[i] << ", "
            << nvcomp_times[i] << ", " 
            << nvcomp_sizes[i] << ", " 
            << (double)data_sizes[i]/(double)nvcomp_sizes[i] << ", "
            << (double)data_sizes[i]/(double)(func_start_times[i]-prev_func_end_times) << "\n";
    }

    logfile.close();

    if (this->mpForwardPressureHostMemory != nullptr) {
        mem_free(this->mpForwardPressureHostMemory);
    }
    delete this->mpForwardPressure;
    this->mpInternalGridBox->Set(WAVE | GB_PRSS | CURR | DIR_Z, initial_internalGridbox_curr);
    this->mpWaveFieldsMemoryHandler->FreeWaveFields(this->mpInternalGridBox);
    delete this->mpInternalGridBox;
}

void TwoPropagation::AcquireConfiguration() {
    if (this->mpConfigurationMap->Contains(OP_K_PROPRIETIES, OP_K_COMPRESSION)) {
        this->mIsCompression = this->mpConfigurationMap->GetValue(OP_K_PROPRIETIES, OP_K_COMPRESSION,
                                                                  this->mIsCompression);
    }
    if (this->mpConfigurationMap->Contains(OP_K_PROPRIETIES, OP_K_WRITE_PATH)) {
        std::string write_path = this->mpConfigurationMap->GetValue(OP_K_PROPRIETIES, OP_K_WRITE_PATH,
                                                                    this->mWritePath);
        mkdir(write_path.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
        this->mWritePath = write_path + "/two_prop";
        mkdir(this->mWritePath.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    }
    if (this->mpConfigurationMap->Contains(OP_K_PROPRIETIES, OP_K_ZFP_TOLERANCE)) {
        this->mZFP_Tolerance = this->mpConfigurationMap->GetValue(OP_K_PROPRIETIES, OP_K_ZFP_TOLERANCE,
                                                                  this->mZFP_Tolerance);
    }
    if (this->mpConfigurationMap->Contains(OP_K_PROPRIETIES, OP_K_ZFP_PARALLEL)) {
        this->mZFP_Parallel = this->mpConfigurationMap->GetValue(OP_K_PROPRIETIES, OP_K_ZFP_PARALLEL,
                                                                 this->mZFP_Parallel);
    }
    if (this->mpConfigurationMap->Contains(OP_K_PROPRIETIES, OP_K_ZFP_RELATIVE)) {
        this->mZFP_IsRelative = this->mpConfigurationMap->GetValue(OP_K_PROPRIETIES, OP_K_ZFP_RELATIVE,
                                                                   this->mZFP_IsRelative);
    }
}

// void TwoPropagation::FetchForward() {
void TwoPropagation::FetchForward(std::string &ckpt_name) {

    uint wnx = this->mpMainGridBox->GetWindowAxis()->GetXAxis().GetActualAxisSize();
    uint wny = this->mpMainGridBox->GetWindowAxis()->GetYAxis().GetActualAxisSize();
    uint wnz = this->mpMainGridBox->GetWindowAxis()->GetZAxis().GetActualAxisSize();


    uint const window_size = wnx * wny * wnz;
    // Retrieve data from files to host buffer
    if ((this->mTimeCounter + 1) % this->mMaxNT == 0) {
        if (this->mIsCompression) {
            string str = this->mWritePath + "/temp_" + to_string(this->mTimeCounter / this->mMaxNT);
            {
                ScopeTimer t("ForwardCollector::Decompression");
                Compressor::Decompress(this->mpForwardPressureHostMemory, wnx, wny, wnz,
                                       this->mMaxNT,
                                       (double) this->mZFP_Tolerance,
                                       this->mZFP_Parallel,
                                       str.c_str(),
                                       this->mZFP_IsRelative);
            }
        } else {
            string str = this->mWritePath + "/temp_" + to_string(this->mTimeCounter / this->mMaxNT);
            {
                ScopeTimer t("IO::ReadForward");
                bin_file_load(str.c_str(), this->mpForwardPressureHostMemory, this->mMaxNT * window_size);
            }
        }
    }
    // Retrieve data from host buffer
    if ((this->mTimeCounter + 1) % this->mMaxDeviceNT == 0) {

        int host_index = (this->mTimeCounter + 1) / this->mMaxDeviceNT - 1;

        Device::MemCpy(this->mpForwardPressure->GetNativePointer(),
                       this->mpForwardPressureHostMemory +
                       (host_index % this->mpMaxNTRatio) * (this->mMaxDeviceNT * window_size),
                       this->mMaxDeviceNT * window_size * sizeof(float),
                       Device::COPY_HOST_TO_DEVICE);
    }
    this->mpInternalGridBox->Set(WAVE | GB_PRSS | CURR | DIR_Z,
                                 this->mpForwardPressure->GetNativePointer() +
                                 ((this->mTimeCounter) % this->mMaxDeviceNT) * window_size);
    this->mTimeCounter--;
}

void TwoPropagation::ResetGrid(bool aIsForwardRun) {


    uint wnx = this->mpMainGridBox->GetWindowAxis()->GetXAxis().GetActualAxisSize();
    uint wny = this->mpMainGridBox->GetWindowAxis()->GetYAxis().GetActualAxisSize();
    uint wnz = this->mpMainGridBox->GetWindowAxis()->GetZAxis().GetActualAxisSize();

    uint const window_size = wnx * wny * wnz;

    if (aIsForwardRun) {
        this->mpMainGridBox->CloneMetaData(this->mpInternalGridBox);
        this->mpMainGridBox->CloneParameters(this->mpInternalGridBox);

        if (this->mpInternalGridBox->GetWaveFields().empty()) {
            this->mpWaveFieldsMemoryHandler->CloneWaveFields(this->mpMainGridBox,
                                                             this->mpInternalGridBox);

            // save the pressure pointer for deletion afterwards
            initial_internalGridbox_curr = mpInternalGridBox->Get(WAVE | GB_PRSS | CURR | DIR_Z)->GetNativePointer();

        } else {
            this->mpWaveFieldsMemoryHandler->CopyWaveFields(this->mpMainGridBox,
                                                            this->mpInternalGridBox);
        }


        this->mTimeCounter = 0;
        if (this->mpForwardPressureHostMemory == nullptr) {
            /// Add one for empty timeframe at the start of the simulation
            /// (The first previous) since SaveForward is called before each step.
            this->mMaxNT = 10000;  // this->mpMainGridBox->GetNT() + 1;


            this->mMaxDeviceNT = 300; // save 100 frames in the Device memory, then reflect to host memory

            this->mpForwardPressureHostMemory = (float *) mem_allocate(
                    (sizeof(float)), this->mMaxNT * window_size, "forward_pressure");

            this->mpForwardPressure = new FrameBuffer<float>();
            this->mpForwardPressure->Allocate(window_size * this->mMaxDeviceNT);

            if (this->mpForwardPressureHostMemory != nullptr) {
                this->mIsMemoryFit = true;
            } else {
                this->mIsMemoryFit = false;
                while (this->mpForwardPressureHostMemory == nullptr) {
                    this->mMaxNT = this->mMaxNT / 2;
                    this->mpForwardPressureHostMemory = (float *) mem_allocate(
                            (sizeof(float)), this->mMaxNT * window_size, "forward_pressure");
                }

                mem_free(this->mpForwardPressureHostMemory);

                // another iteration as a safety measure
                this->mMaxNT = this->mMaxNT / 2;
                this->mpForwardPressureHostMemory = (float *) mem_allocate(
                        (sizeof(float)), this->mMaxNT * window_size, "forward_pressure");
            }

            mpMaxNTRatio = mMaxNT / this->mMaxDeviceNT;

        }

        this->mpTempCurr = this->mpMainGridBox->Get(WAVE | GB_PRSS | CURR | DIR_Z)->GetNativePointer();
        this->mpTempNext = this->mpMainGridBox->Get(WAVE | GB_PRSS | NEXT | DIR_Z)->GetNativePointer();

        Device::MemSet(this->mpForwardPressure->GetNativePointer(), 0.0f,
                       window_size * sizeof(float));
        Device::MemSet(this->mpForwardPressure->GetNativePointer() + window_size, 0.0f,
                       window_size * sizeof(float));


        if (this->mpParameters->GetEquationOrder() == SECOND) {
            this->mpTempPrev = this->mpMainGridBox->Get(WAVE | GB_PRSS | PREV | DIR_Z)->GetNativePointer();

            Device::MemSet(this->mpForwardPressure->GetNativePointer() + 2 * window_size, 0.0f,
                           window_size * sizeof(float));
        }


        if (this->mpParameters->GetEquationOrder() == SECOND) {
            this->mpMainGridBox->Set(WAVE | GB_PRSS | PREV | DIR_Z, this->mpForwardPressure->GetNativePointer());
            this->mpMainGridBox->Set(WAVE | GB_PRSS | CURR | DIR_Z,
                                     this->mpForwardPressure->GetNativePointer() + window_size);
            // Save forward is called before the kernel in the engine.
            // When called will advance pressure next to the right point.
            this->mpMainGridBox->Set(WAVE | GB_PRSS | NEXT | DIR_Z,
                                     this->mpForwardPressure->GetNativePointer() + window_size);
        } else {
            this->mpMainGridBox->Set(WAVE | GB_PRSS | CURR | DIR_Z, this->mpForwardPressure->GetNativePointer());
            this->mpMainGridBox->Set(WAVE | GB_PRSS | NEXT | DIR_Z,
                                     this->mpForwardPressure->GetNativePointer() + window_size);
        }
    } else {
        Device::MemSet(this->mpTempCurr, 0.0f, window_size * sizeof(float));
        if (this->mpParameters->GetEquationOrder() == SECOND) {
            Device::MemSet(this->mpTempPrev, 0.0f, window_size * sizeof(float));
        }

        for (auto const &wave_field : this->mpMainGridBox->GetWaveFields()) {
            if (GridBox::Includes(wave_field.first, GB_PRTC)) {
                Device::MemSet(wave_field.second->GetNativePointer(), 0.0f, window_size * sizeof(float));
            }
        }

        if (!this->mIsMemoryFit) {
            this->mTimeCounter++;
            this->mpInternalGridBox->Set(WAVE | GB_PRSS | CURR | DIR_Z,
                                         this->mpMainGridBox->Get(WAVE | GB_PRSS | CURR | DIR_Z)->GetNativePointer());
        } else {
            // Pressure size will be minimized in FetchForward call at first step.
            this->mpInternalGridBox->Set(WAVE | GB_PRSS | CURR | DIR_Z,
                                         this->mpMainGridBox->Get(WAVE | GB_PRSS | CURR | DIR_Z)->GetNativePointer() +
                                         window_size);
        }

        if (this->mpParameters->GetEquationOrder() == SECOND) {
            this->mpMainGridBox->Set(WAVE | GB_PRSS | PREV | DIR_Z, this->mpTempPrev);
        }
        this->mpMainGridBox->Set(WAVE | GB_PRSS | CURR | DIR_Z, this->mpTempCurr);
        this->mpMainGridBox->Set(WAVE | GB_PRSS | NEXT | DIR_Z, this->mpTempNext);
    }
}

// void TwoPropagation::SaveForward() {
void TwoPropagation::SaveForward(std::string &ckpt_name) {

    uint64_t func_start = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::system_clock::now().time_since_epoch()).count();

    uint wnx = this->mpMainGridBox->GetWindowAxis()->GetXAxis().GetActualAxisSize();
    uint wny = this->mpMainGridBox->GetWindowAxis()->GetYAxis().GetActualAxisSize();
    uint wnz = this->mpMainGridBox->GetWindowAxis()->GetZAxis().GetActualAxisSize();

    uint const window_size = wnx * wny * wnz;

    this->mTimeCounter++;

    // Transfer from Device memory to host memory
    LoggerSystem *Logger = LoggerSystem::GetInstance();
    // veloc_client->mem_protect(0, ptr, window_size, sizeof(float), DEFAULT);
    if (this->mTimeCounter <= 1) {
        Logger->Info() <<"Total number of checkpoints: " << this->mpMainGridBox->GetNT() << " \n";
        Logger->Info() << "Size of memory region: " << window_size*sizeof(float) << ", block size: " << this->mMaxDeviceNT << ", total: " << (uint64_t)this->mpMainGridBox->GetNT()*(uint64_t)window_size*(uint64_t)sizeof(float) << "\n";
        for(int i=(this->mpMainGridBox->GetNT()/this->mMaxDeviceNT)-1; i>=0; i--) {
            VELOC_Prefetch_enqueue(ckpt_name.c_str(), i, 0);
        }
        
        VELOC_Mem_protect(0, this->mpForwardPressure->GetNativePointer(), window_size*this->mMaxDeviceNT, sizeof(float), should_compress);
        ckpt_id = 0;
    }
    
    if ((this->mTimeCounter + 1) % this->mMaxDeviceNT == 0) {
        VELOC_Checkpoint(ckpt_name.c_str(), ckpt_id);
	    ckpt_id++;
    }
    uint64_t d2h_end = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::system_clock::now().time_since_epoch()).count();

    // Save host memory to file
    uint64_t h2f_start = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::system_clock::now().time_since_epoch()).count();
    // Logger->Info() << "In save forward function from host " << this->mTimeCounter << " maxdeviceNT " << this->mMaxNT << " compr " << this->mIsCompression << '\n';
    if ((this->mTimeCounter + 1) % this->mMaxNT == 0) {
        // if (false) { // Always write uncompressed data to file
        if (this->mIsCompression) {
            string str = this->mWritePath + "/temp_" + to_string(this->mTimeCounter / this->mMaxNT);
            {
                ScopeTimer t("ForwardCollector::Compression");;
                Compressor::Compress(this->mpForwardPressureHostMemory, wnx, wny, wnz,
                                     this->mMaxNT,
                                     (double) this->mZFP_Tolerance,
                                     this->mZFP_Parallel,
                                     str.c_str(),
                                     this->mZFP_IsRelative);
            }
        } else {
            string str =
                    this->mWritePath + "/temp_" + to_string(this->mTimeCounter / this->mMaxNT);
            Logger->Info() << "Saving file at " << this->mWritePath << "\n";
            {
                ScopeTimer t("IO::WriteForward");
                bin_file_save(str.c_str(), this->mpForwardPressureHostMemory, this->mMaxNT * window_size);
            }
        }
    }
    uint64_t h2f_end = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::system_clock::now().time_since_epoch()).count();

    uint64_t d2d_start = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::system_clock::now().time_since_epoch()).count();
    this->mpMainGridBox->Set(WAVE | GB_PRSS | CURR | DIR_Z,
                             this->mpForwardPressure->GetNativePointer() +
                             ((this->mTimeCounter) % this->mMaxDeviceNT) * window_size);
    this->mpMainGridBox->Set(WAVE | GB_PRSS | NEXT | DIR_Z,
                             this->mpForwardPressure->GetNativePointer() +
                             ((this->mTimeCounter + 1) % this->mMaxDeviceNT) * window_size);
    if (this->mpParameters->GetEquationOrder() == SECOND) {
        this->mpMainGridBox->Set(WAVE | GB_PRSS | PREV | DIR_Z,
                                 this->mpForwardPressure->GetNativePointer() +
                                 ((this->mTimeCounter - 1) % this->mMaxDeviceNT) * window_size);
    }
    uint64_t d2d_end = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::system_clock::now().time_since_epoch()).count();

    uint64_t func_end = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::system_clock::now().time_since_epoch()).count();

    float nvcomp_ms = 0;
    checkCuda(cudaEventSynchronize(nvcomp_end));
    checkCuda(cudaEventElapsedTime(&nvcomp_ms, nvcomp_start, nvcomp_end));
    uint64_t nvcomp_time = nvcomp_ms*1e6;
    if ( ((this->mTimeCounter + 1) % this->mMaxDeviceNT == 0) || ((this->mTimeCounter + 1) % this->mMaxNT == 0) ) {  
        // Logger->Info() << "Ckpt: " << this->mTimeCounter+1 
        //     << " Fstart: " << func_start 
        //     << " Fend: " << func_end
        //     << " D2D: " << d2d_end-d2d_start 
        //     << " D2H: " << d2h_end-d2h_start
        //     << " H2F: " << h2f_end-h2f_start
        //     << " Func: " << func_end-func_start << "\n";
        iter_counts.push_back(this->mTimeCounter+1);
        func_start_times.push_back(func_start);
        func_end_times.push_back(func_end);
        // nvcomp_times.push_back(nvcomp_end-nvcomp_start);
        nvcomp_times.push_back(nvcomp_time);
        nvcomp_sizes.push_back(cmpr_size);
        data_sizes.push_back(data_size);
    }
}

void TwoPropagation::SetComputationParameters(ComputationParameters *apParameters) {
    LoggerSystem *Logger = LoggerSystem::GetInstance();
    this->mpParameters = (ComputationParameters *) apParameters;
    if (this->mpParameters == nullptr) {
        Logger->Error() << "No computation parameters provided... Terminating..." << '\n';
        exit(EXIT_FAILURE);
    }
}

void TwoPropagation::SetGridBox(GridBox *apGridBox) {
    LoggerSystem *Logger = LoggerSystem::GetInstance();
    this->mpMainGridBox = apGridBox;
    if (this->mpMainGridBox == nullptr) {
        Logger->Error() << "Not a compatible GridBox... Terminating..." << '\n';
        exit(EXIT_FAILURE);
    }

    /*
     * in case of two propagation the next buffer and previous buffer get separated
     * so we need to register a new wave field with actual allocated memory
     */
    auto framebuffer = new FrameBuffer<float>();
    this->mpMainGridBox->RegisterWaveField(WAVE | GB_PRSS | NEXT | DIR_Z, framebuffer);


    framebuffer->Allocate(
            this->mpMainGridBox->GetWindowAxis()->GetXAxis().GetActualAxisSize() *
            this->mpMainGridBox->GetWindowAxis()->GetYAxis().GetActualAxisSize() *
            this->mpMainGridBox->GetWindowAxis()->GetZAxis().GetActualAxisSize(),
            mpParameters->GetHalfLength(),
            "next pressure");
}

void TwoPropagation::SetDependentComponents(
        ComponentsMap<DependentComponent> *apDependentComponentsMap) {
    LoggerSystem *Logger = LoggerSystem::GetInstance();
    HasDependents::SetDependentComponents(apDependentComponentsMap);

    this->mpWaveFieldsMemoryHandler =
            (WaveFieldsMemoryHandler *)
                    this->GetDependentComponentsMap()->Get(MEMORY_HANDLER);
    if (this->mpWaveFieldsMemoryHandler == nullptr) {
        Logger->Error() << "No Wave Fields Memory Handler provided... " << "Terminating..." << '\n';
        exit(EXIT_FAILURE);
    }
}

GridBox *TwoPropagation::GetForwardGrid() {
    return this->mpInternalGridBox;
}