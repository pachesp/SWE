#ifndef TOOLS_PRECICE_H
#define TOOLS_PRECICE_H

#include "precice/SolverInterface.hpp"
using namespace precice;
using namespace precice::constants;

#include "blocks/SWE_WavePropagationBlock.hh"
#include "tools/help.hh"

struct PreciceData{

  int  snd_heightId;
  int  snd_huId;
  int  snd_hvId;

  double* snd_height_db;
  double* snd_hu_db;
  double* snd_hv_db;

  int recv_heightId;
  int recv_huId;
  int recv_hvId;

  double* recv_height_db;
  double* recv_hu_db;
  double* recv_hv_db;

  int* vertexIDs;

  Float2D* CP_height_f2d;
  Float2D* CP_hu_f2d ;
  Float2D* CP_hv_f2d;

};

void write_preCICE(SolverInterface &interface, SWE_Block &wavePropagationBlock, PreciceData *data,
                  int size, int columNr);

void read_preCICE(SolverInterface &interface, SWE_Block &wavePropagationBlock, SWE_Block1D* ghoshtBlock,
                  PreciceData *data, int size, int columNr = 0);

void writeCheckpoint(PreciceData *data, SWE_Block &wavePropagationBlock, float time, float &time_CP, int size, int columNr);

void restoreCheckpoint(PreciceData *data, SWE_Block &wavePropagationBlock, float &time, float time_CP, int size, int columNr);

// void storeData_preCICE(SolverInterface &interface, SWE_Block &wavePropagationBlock, PreciceData *data,
//                      int columNr, int size);
//
// void writeData_preCICE(SolverInterface &interface, SWE_Block &wavePropagationBlock, PreciceData *data,
//                 int columNr, int size);
#endif // TOOLS_PRECICE_H
