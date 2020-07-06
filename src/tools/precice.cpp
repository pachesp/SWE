#include "tools/precice.hh"

void write_preCICE(SolverInterface &interface, SWE_Block &wavePropagationBlock, PreciceData *data, int size, int columNr){

    for(int i = 0; i < size ; i++){ //Data stored column-wise
      data->snd_height_db[i] = (double)wavePropagationBlock.getWaterHeight()[columNr][i];
      data->snd_hu_db[i] = (double)wavePropagationBlock.getDischarge_hu()[columNr][i];
      data->snd_hv_db[i] = (double)wavePropagationBlock.getDischarge_hv()[columNr][i];
    }

    interface.writeBlockScalarData(data->snd_heightId, size, data->vertexIDs, data->snd_height_db);
    interface.writeBlockScalarData(data->snd_huId, size, data->vertexIDs, data->snd_hu_db);
    interface.writeBlockScalarData(data->snd_hvId, size, data->vertexIDs, data->snd_hv_db);
}

void read_preCICE(SolverInterface &interface, SWE_Block &wavePropagationBlock,
                  SWE_Block1D* ghoshtBlock, PreciceData *data, int size, int columNr){

    interface.readBlockScalarData(data->recv_heightId, size, data->vertexIDs, data->recv_height_db);
    interface.readBlockScalarData(data->recv_huId, size, data->vertexIDs, data->recv_hu_db);
    interface.readBlockScalarData(data->recv_hvId, size, data->vertexIDs, data->recv_hv_db);

    float* h = doublePointer2floatPointer(data->recv_height_db, size);
    float* hu = doublePointer2floatPointer(data->recv_hu_db, size);
    float* hv = doublePointer2floatPointer(data->recv_hv_db, size);

    // Data sent by left neighbour from precice
    SWE_Block1D* newBlock = new SWE_Block1D{h, hu, hv, 1};
    ghoshtBlock->copyFrom(newBlock, size );

    if(newBlock) delete newBlock;
    if(h) delete h;
    if(hu) delete hu;
    if(hv) delete hv;
}

void writeCheckpoint(PreciceData *data, SWE_Block &wavePropagationBlock, float time, float &time_CP, int size, int columNr){

    if(data->CP_height_f2d) delete data->CP_height_f2d;
    if(data->CP_hu_f2d) delete data->CP_hu_f2d;
    if(data->CP_hv_f2d) delete data->CP_hv_f2d;


    time_CP = time;

    data->CP_height_f2d = new Float2D(wavePropagationBlock.getWaterHeight(), false);
    data->CP_hu_f2d = new Float2D(wavePropagationBlock.getDischarge_hu(), false);
    data->CP_hv_f2d = new Float2D(wavePropagationBlock.getDischarge_hv(), false);

}

void restoreCheckpoint(PreciceData *data, SWE_Block &wavePropagationBlock, float &time, float time_CP, int size, int columNr){
    time = time_CP;

    for(int j = 0; j < size ; j++){
        for(int i = 0; i < size ; i++){

            wavePropagationBlock.getWaterHeight()[i][j]  =  (*(data->CP_height_f2d))[i][j];
            wavePropagationBlock.getDischarge_hu()[i][j] =  (*(data->CP_hu_f2d))[i][j];
            wavePropagationBlock.getDischarge_hv()[i][j] =  (*(data->CP_hv_f2d))[i][j];
        }
    }
}

// void storeData_preCICE(SolverInterface &interface, SWE_Block &wavePropagationBlock, PreciceData *data, int columNr, int size){
//     for(int i = 0; i < size ; i++){
//       data->snd_height_db[i] = wavePropagationBlock.getWaterHeight().float2D2doublePointer()[columNr * size + i];
//       data->snd_hu_db[i] = wavePropagationBlock.getDischarge_hu().float2D2doublePointer()[columNr * size + i];
//       data->snd_hv_db[i] = wavePropagationBlock.getDischarge_hv().float2D2doublePointer()[columNr * size + i];
//     }
// }
//
// void writeData_preCICE(SolverInterface &interface, SWE_Block &wavePropagationBlock, PreciceData *data, int columNr, int size){
//     interface.writeBlockScalarData(data->snd_heightId, size, data->vertexIDs, data->snd_height_db);
//     interface.writeBlockScalarData(data->snd_huId, size, data->vertexIDs, data->snd_hu_db);
//     interface.writeBlockScalarData(data->snd_hvId, size, data->vertexIDs, data->snd_hv_db);
// }
