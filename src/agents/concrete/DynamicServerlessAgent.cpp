//
// Created by mennatallah on 6/9/20.
//

#if defined(USING_MPI)

#include <stbx/agents/concrete/DynamicServerlessAgent.hpp>

#include <mpi.h>

using namespace std;
using namespace stbx::agents;
using namespace operations::dataunits;

DynamicServerlessAgent::~DynamicServerlessAgent() {
    for (auto it : this->mStacks) {
        delete[] it;
    }
}

GridBox *DynamicServerlessAgent::Initialize() {
    this->mpGridBox = mpEngine->Initialize();

    int provided;
    MPI_Init_thread(&this->argc, &this->argv, MPI_THREAD_FUNNELED, &provided);

    if (provided != MPI_THREAD_FUNNELED) {
        std::cerr << "Warning MPI did not provide MPI_THREAD_FUNNELED..." << std::endl;
    }

    this->mCommunication = MPI_COMM_WORLD;
    MPI_Comm_rank(this->mCommunication, &this->self);
    MPI_Comm_size(this->mCommunication, &this->mProcessCount);

    return this->mpGridBox;
}

void DynamicServerlessAgent::BeforeMigration() {
    if (this->self == 0) {
        this->mPossibleShots = mpEngine->GetValidShots();
        this->mShotsSize = this->mPossibleShots.size();

        // The server is assigned the first shot ID to work on.
        this->mMasterCurrentShotID = this->mPossibleShots[this->mShotTracker];

        /*Incrementing the shot ID tracker.*/
        this->mShotTracker++;

        // Incrementing the count that indicates the first shots iterator
        this->mShotsCount++;

        // The server(process with rank 0) allocates the mig_flag to each client process:
        // 0 for a process if (process rank < # of shots) -->it will be working on a shot.
        // 1 if (process rank >=# of shots) -->there is no enough shots for this process.
        for (int i = 1; i < this->mProcessCount; i++) {
            // process_rank >= number of shots
            if (this->mShotsCount == this->mShotsSize) {
                this->flag[4] = 1;
            } else {
                //  process rank < # of shots
                this->mShotsCount++;
            }
            // Server sends a mig_flag to each client process
            MPI_Send(&this->flag[4], 1, MPI_INT, i, 20, this->mCommunication);
        }

        // Server sends one shot ID to each client process that will be working on a shot.
        for (int i = 1; i < min<int>(this->mProcessCount, this->mShotsSize); ++i) {
            flag[0] = this->mPossibleShots[i];
            MPI_Send(&this->flag[0], 1, MPI_INT, i, 2, this->mCommunication);
            // Server updates the tracking of the coming shot ID to be processed
            this->mShotTracker++;
        }
    } else {
        // Each client process receives the mig_flag form the client
        MPI_Recv(&this->flag[4], 1, MPI_INT, 0, 20, mCommunication, &mStatus);

        //Getting the valid shots.
        vector<uint> temporary_shots;
        temporary_shots = mpEngine->GetValidShots();
    }
}

void DynamicServerlessAgent::AfterMigration() {
    if (this->self == 0) {
        while (this->mShotTracker < this->mShotsSize) {
            // Server is assigned the coming shot ID to work on
            this->mMasterCurrentShotID = this->mPossibleShots[mShotTracker];

            // Server updates the tracking of the coming shot ID to be processed
            this->mShotTracker++;

            if (this->mShotTracker == this->mShotsSize) {
                this->mHasShotDynamic = true;
                break;
            } else {
                for (int m = 0; m < this->mProcessCount - 1; ++m) {
                    /*The server is sending to all client processes a shot_flag of 0 which means there are still shots available*/
                    MPI_Send(&this->flag[1], 1, MPI_INT, m + 1, 3, this->mCommunication);
                    /*The server receives the av_flag and rank from the finished client*/
                    MPI_Recv(&this->flag[2], 2, MPI_INT, m + 1, 14, this->mCommunication, &this->mStatus);

                    /*If the av_flag equals 1*/
                    if (this->flag[2] == 1) {

                        this->flag[0] = this->mPossibleShots[mShotTracker];
                        int finished = this->flag[3];
                        // Server sends the coming shot ID to be processed
                        // to the client who send the av_flag
                        MPI_Send(&this->flag[0], 1, MPI_INT, finished, 2, this->mCommunication);
                        /*The server updates the tracking of the coming shot ID to be processed*/
                        this->mShotTracker++;
                        /*If there are no shots available*/
                        if (this->mShotTracker == this->mShotsSize) {
                            break;
                        }
                    }
                }
            }
            return;
        }
        // As there are no shots available the shot_flag will be 1
        this->flag[1] = 1;
        // The server sends to all client processes a shot_flag of 1
        for (int i = 1; i < min<int>(this->mProcessCount, this->mShotsSize); ++i) {
            MPI_Send(&this->flag[1], 1, MPI_INT, i, 3, this->mCommunication);
        }
    } else {
        // Each client process receives the shot_flag from the server.
        MPI_Recv(&this->flag[1], 1, MPI_INT, 0, 3, this->mCommunication, &this->mStatus);

        if (this->flag[1] == 0) {
            // Indicates that it is now available
            this->flag[2] = 1;
            // Current process rank
            this->flag[3] = self;

            // The client who finishes, sends av_flag (flag[2])
            // to the server to state that it is now available.
            // It sends also its rank to the server as an identification.
            MPI_Send(&this->flag[2], 2, MPI_INT, 0, 14, this->mCommunication);
        }
    }
}

void DynamicServerlessAgent::BeforeFinalize() {}

MigrationData *DynamicServerlessAgent::AfterFinalize(MigrationData *apMigrationData) {
    MigrationData *md = apMigrationData;

    uint size = md->GetGridSize(X_AXIS) *
                md->GetGridSize(Y_AXIS) *
                md->GetGridSize(Z_AXIS) *
                md->GetGatherDimension();

    for (int i = 0; i < md->GetResults().size(); ++i) {
        mStacks.push_back(new float[size]);
        MPI_Reduce(md->GetResultAt(i)->GetData(), this->mStacks[i],
                   size, MPI_FLOAT, MPI_SUM, 0, this->mCommunication);
    }
    if (this->self == 0) {
        for (int i = 0; i < md->GetResults().size(); ++i) {
            md->SetResults(i, new Result(this->mStacks[i]));
        }
    }
    MPI_Finalize();
    if (this->self != 0) {
        exit(0);
    }
    return md;
}

bool DynamicServerlessAgent::HasNextShot() {
    if (self == 0) {
        if (this->mHasShotDynamic && this->flag[1] == 1) {
            this->mCount++;
            return this->mCount > 1 ? false : true;
        } else if (this->flag[1] == 1) {
            return false;
        } else {
            return true;
        }
    } else {
        if (this->flag[4] == 1) {
            return false;
        } else {
            if (this->flag[1] == 0) {
                // Each client process receives the shot ID to work on from the server.
                MPI_Recv(&this->flag[0], 1, MPI_INT, 0, 2, this->mCommunication, &this->mStatus);
                return true;
            } else {
                return false;
            }
        }
    }
}

vector<uint> DynamicServerlessAgent::GetNextShot() {
    if (this->self == 0) {
        this->mProcessShots.clear();
        this->mProcessShots.push_back(this->mMasterCurrentShotID);
        return this->mProcessShots;
    } else {
        this->mPossibleShots.clear();
        this->mPossibleShots.push_back(this->flag[0]);
        return this->mPossibleShots;
    }
}

#endif
