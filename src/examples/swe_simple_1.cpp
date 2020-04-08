/**
 * @file
 * This file is part of SWE.
 *
 * @author Alexander Breuer (breuera AT in.tum.de, http://www5.in.tum.de/wiki/index.php/Dipl.-Math._Alexander_Breuer)
 *         Michael Bader (bader AT in.tum.de, http://www5.in.tum.de/wiki/index.php/Univ.-Prof._Dr._Michael_Bader)
 *
 * @section LICENSE
 *
 * SWE is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * SWE is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with SWE.  If not, see <http://www.gnu.org/licenses/>.
 *
 *
 * @section DESCRIPTION
 *
 * Basic setting of SWE, which uses a wave propagation solver and an artificial or ASAGI scenario on a single block.
 */
#include <cassert>
#include <cstdlib>
#include <string>
#include <iostream>


#include "blocks/SWE_WavePropagationBlock.hh"
#include "writer/VtkWriter.hh"
#include "scenarios/SWE_simple_scenarios.hh"


#include "tools/args.hh"
#include "tools/help.hh"
#include "tools/Logger.hh"
#include "tools/ProgressBar.hh"

//precice
#include "precice/SolverInterface.hpp"
using namespace precice;
using namespace precice::constants;


/**
 * Main program for the simulation on a single SWE_WavePropagationBlock.
 */
int main( int argc, char** argv ) {
  /**
   * Initialization.
   */
  // Parse command line parameters
  tools::Args args;

  args.addOption("grid-size-x", 'x', "Number of cells in x direction");
  args.addOption("grid-size-y", 'y', "Number of cells in y direction");
  args.addOption("output-basepath", 'o', "Output base file name");

  tools::Args::Result ret = args.parse(argc, argv);

  switch (ret)
  {
  case tools::Args::Error:
	  return 1;
  case tools::Args::Help:
	  return 0;
  }

  //! number of grid cells in x- and y-direction.
  int l_nX, l_nY;

  //! l_baseName of the plots.
  std::string l_baseName;

  // read command line parameters
  l_nX = args.getArgument<int>("grid-size-x");
  l_nY = args.getArgument<int>("grid-size-y");
  l_baseName = args.getArgument<std::string>("output-basepath");

  // create a simple artificial scenario
  SWE_FranciscoScenario l_scenario;

  //! number of checkpoints for visualization (at each checkpoint in time, an output file is written).
  // int l_numberOfCheckPoints = 50;
  int l_numberOfCheckPoints = l_scenario.setNumberCheckpoints();

  //! size of a single cell in x- and y-direction
  float l_dX, l_dY;

  // compute the size of a single cell
  l_dX = (l_scenario.getBoundaryPos(BND_RIGHT) - l_scenario.getBoundaryPos(BND_LEFT) )/l_nX;
  l_dY = (l_scenario.getBoundaryPos(BND_TOP) - l_scenario.getBoundaryPos(BND_BOTTOM) )/l_nY;

  // create a single wave propagation block
  SWE_WavePropagationBlock l_wavePropgationBlock(l_nX,l_nY,l_dX,l_dY);

  //! origin of the simulation domain in x- and y-direction
  float l_originX, l_originY;

  // get the origin from the scenario
  l_originX = l_scenario.getBoundaryPos(BND_LEFT);
  l_originY = l_scenario.getBoundaryPos(BND_BOTTOM);

  //***************preCICE**************************
  // *
  // *
  std::string configFileName("precice-config.xml");
  std::string solverName = "Solver1";
  SolverInterface interface(solverName, configFileName, 0, 1);
  int dimensions = interface.getDimensions();
  int meshID = interface.getMeshID("Solver1_Nodes");
  int heightS1Id = interface.getDataID("heightS1", meshID);
  int huS1Id = interface.getDataID("huS1", meshID);
  int hvS1Id = interface.getDataID("hvS1", meshID);
  int heightS2Id = interface.getDataID("heightS2", meshID);
  int huS2Id = interface.getDataID("huS2", meshID);
  int hvS2Id = interface.getDataID("hvS2", meshID);
  int* vertexIDs;
  vertexIDs = new int[(l_nY + 2)];
  double* grid;
  grid = new double[dimensions * (l_nY + 2)];
  int count=0;
  for (int j=0; j <= l_nY + 1 ; j++){
      grid[count++] = 0;
      grid[count++] = j;
    }
  interface.setMeshVertices(meshID, (l_nY+2) , grid, vertexIDs);
  cout << "Initialize preCICE..." << endl;
  float precice_dt = interface.initialize();

  double *heightS1_db = new double[(l_nX+2)];
  double *huS1_db = new double[(l_nX+2)];
  double *hvS1_db = new double[(l_nX+2)];
  //
  double *heightS2_db = new double[(l_nX+2)];
  double *huS2_db = new double[(l_nX+2)];
  double *hvS2_db = new double[(l_nX+2)];
  // *
  // *
  //***************preCICE**************************


  // initialize the wave propagation block
  l_wavePropgationBlock.initScenario(l_originX, l_originY, l_scenario);

  //! time when the simulation ends.
  float l_endSimulation = l_scenario.endSimulation();
  // float l_endSimulation = 200.0;

  //! checkpoints when output files are written.
  float* l_checkPoints = new float[l_numberOfCheckPoints+1];

  // compute the checkpoints in time
  for(int cp = 0; cp <= l_numberOfCheckPoints; cp++) {
     l_checkPoints[cp] = cp*(l_endSimulation/l_numberOfCheckPoints);
  }



  //***************preCICE SEND begin**************************
    // for(int i = 0; i < l_nX +2 ; i++){
    //   heightS1_db[i] = l_wavePropgationBlock.getWaterHeight().float2D2doublePointer()[l_nY * (l_nX+2)+(i)];
    //   huS1_db[i] = l_wavePropgationBlock.getDischarge_hu().float2D2doublePointer()[l_nY * (l_nX+2)+(i)];
    //   hvS1_db[i] = l_wavePropgationBlock.getDischarge_hv().float2D2doublePointer()[l_nY * (l_nX+2)+(i)];
    // }
    //
    // interface.writeBlockScalarData(heightS1Id, (l_nX+2), vertexIDs, heightS1_db);
    // interface.writeBlockScalarData(huS1Id, (l_nX+2), vertexIDs, huS1_db);
    // interface.writeBlockScalarData(hvS1Id, (l_nX+2), vertexIDs, hvS1_db);
  //***************preCICE SEND end**************************


  // // Debbugging
  //  std::cout << "output1" << '\n';
  //  for(int j = 0; j < l_nX +2 ; j++){
  //     for(int i = 0; i < l_nY + 2; i++){
  //       std::cout << l_wavePropgationBlock.getWaterHeight().float2D2doublePointer()[i*(l_nX+2)+(j)] << "\t";
  //   }
  //   std::cout <<"\n";
  // }
  //
  //  std::cout << "output 2" << '\n';
  //  for(int i = 0; i < l_nY +2; i++){
  //    std::cout << heightS1_db[i] << '\n';
  //  }

  //
  // //***************preCICE RCV begin**************************
  //   interface.readBlockScalarData(heightS2Id, (l_nX + 2), vertexIDs, heightS2_db);
  //   interface.readBlockScalarData(huS2Id, (l_nX + 2), vertexIDs, huS2_db);
  //   interface.readBlockScalarData(hvS2Id, (l_nX + 2), vertexIDs, hvS2_db);
  //
  //   // Data sent by left neighbour from precice
  //   SWE_Block1D* rightNeighbourDataPreCICE = new SWE_Block1D{ doublePointer2floatPointer(heightS2_db, l_nX + 2),
  //                                                             doublePointer2floatPointer(huS2_db, l_nX + 2),
  //                                                             doublePointer2floatPointer(hvS2_db, l_nX + 2), 1 };
  //
  //   deepSWE_Block1DCopy(rightNeighbourDataPreCICE, l_rightInflow, l_nX+2);
  // //***************preCICE RCV end**************************


  // Init fancy progressbar
  tools::ProgressBar progressBar(l_endSimulation);

  // write the output at time zero
  tools::Logger::logger.printOutputTime((float) 0.);
  progressBar.update(0.);

  std::string l_fileName = generateBaseFileName(l_baseName,0,0);
  //boundary size of the ghost layers
  io::BoundarySize l_boundarySize = {{1, 1, 1, 1}};

  // consturct a VtkWriter
  io::VtkWriter l_writer( l_fileName,
		  l_wavePropgationBlock.getBathymetry(),
		  l_boundarySize,
		  l_nX, l_nY,
		  l_dX, l_dY );

  // Write zero time step
  l_writer.writeTimeStep( l_wavePropgationBlock.getWaterHeight(),
                          l_wavePropgationBlock.getDischarge_hu(),
                          l_wavePropgationBlock.getDischarge_hv(),
                          // l_wavePropgationBlock.getDummy(),
                          (float) 0.);

  /**
   * Simulation.
   */
  // print the start message and reset the wall clock time
  progressBar.clear();
  tools::Logger::logger.printStartMessage();
  tools::Logger::logger.initWallClockTime(time(NULL));

  //! simulation time.
  float l_t = 0.0;
  progressBar.update(l_t);

  SWE_Block1D* l_rightInflow  = l_wavePropgationBlock.grabGhostLayer(BND_RIGHT);

  unsigned int l_iterations = 0;
  int c=1;
  while(interface.isCouplingOngoing()){

    // loop over checkpoints
    if(l_t < l_checkPoints[c]){

    //***************preCICE SEND begin**************************
      for(int i = 0; i < l_nX +2 ; i++){
        heightS1_db[i] = l_wavePropgationBlock.getWaterHeight().float2D2doublePointer()[l_nY * (l_nX+2)+(i)];
        huS1_db[i] = l_wavePropgationBlock.getDischarge_hu().float2D2doublePointer()[l_nY * (l_nX+2)+(i)];
        hvS1_db[i] = l_wavePropgationBlock.getDischarge_hv().float2D2doublePointer()[l_nY * (l_nX+2)+(i)];
      }

      interface.writeBlockScalarData(heightS1Id, (l_nX+2), vertexIDs, heightS1_db);
      interface.writeBlockScalarData(huS1Id, (l_nX+2), vertexIDs, huS1_db);
      interface.writeBlockScalarData(hvS1Id, (l_nX+2), vertexIDs, hvS1_db);
    //***************preCICE SEND end**************************

    //***************preCICE RCV begin**************************
      interface.readBlockScalarData(heightS2Id, (l_nY + 2), vertexIDs, heightS2_db);
      interface.readBlockScalarData(huS2Id, (l_nY + 2), vertexIDs, huS2_db);
      interface.readBlockScalarData(hvS2Id, (l_nY + 2), vertexIDs, hvS2_db);

      // Data sent by left neighbour from precice
      SWE_Block1D* rightNeighbourDataPreCICE = new SWE_Block1D{ doublePointer2floatPointer(heightS2_db, l_nY + 2),
                                                              doublePointer2floatPointer(huS2_db, l_nY + 2),
                                                              doublePointer2floatPointer(hvS2_db, l_nY + 2), l_nY + 2 };

      // deepSWE_Block1DCopy(rightNeighbourDataPreCICE, l_rightInflow, l_nY+2);
      l_rightInflow->copyFrom(rightNeighbourDataPreCICE, l_nY+2 );
    //***************preCICE RCV end**************************

      // set values in ghost cells:
      l_wavePropgationBlock.setGhostLayer();

      // reset the cpu clock
      tools::Logger::logger.resetClockToCurrentTime("Cpu");

      // approximate the maximum time step
      // TODO: This calculation should be replaced by the usage of the wave speeds occuring during the flux computation
      // Remark: The code is executed on the CPU, therefore a "valid result" depends on the CPU-GPU-synchronization.
      l_wavePropgationBlock.computeMaxTimestep();

      // compute numerical flux on each edge
      l_wavePropgationBlock.computeNumericalFluxes();

      //! maximum allowed time step width.
      float l_maxTimeStepWidth = l_wavePropgationBlock.getMaxTimestep();

      // update the cell values
      l_wavePropgationBlock.updateUnknowns(l_maxTimeStepWidth);

      l_maxTimeStepWidth = std::min(l_maxTimeStepWidth, precice_dt );

      precice_dt = interface.advance(l_maxTimeStepWidth);

      // update the cpu time in the logger
      tools::Logger::logger.updateTime("Cpu");

      // update simulation time with time step width.
      l_t += l_maxTimeStepWidth;
      l_iterations++;

      // print the current simulation time
      progressBar.clear();
      tools::Logger::logger.printSimulationTime(l_t);
      progressBar.update(l_t);
    }
    else{

      // print current simulation time of the output
      progressBar.clear();
      tools::Logger::logger.printOutputTime(l_t);
      progressBar.update(l_t);


      // write output
      l_writer.writeTimeStep( l_wavePropgationBlock.getWaterHeight(),
                              l_wavePropgationBlock.getDischarge_hu(),
                              l_wavePropgationBlock.getDischarge_hv(),
                              // dummyFloat2D_n,
                              l_t);
      c++;
    }
  }

  interface.finalize();

  /**
   * Finalize.
   */
  // write the statistics message
  progressBar.clear();
  tools::Logger::logger.printStatisticsMessage();

  // print the cpu time
  tools::Logger::logger.printTime("Cpu", "CPU time");

  // print the wall clock time (includes plotting)
  tools::Logger::logger.printWallClockTime(time(NULL));

  // printer iteration counter
  tools::Logger::logger.printIterationsDone(l_iterations);

  return 0;
}

// // Debbugging
//  std::cout << "output1" << '\n';
//  for(int j = 0; j < l_nX +2 ; j++){
//     for(int i = 0; i < l_nY + 2; i++){
//       std::cout << l_wavePropgationBlock.getWaterHeight().float2D2doublePointer()[i*(l_nX+2)+(j)] << "\t";
//   }
//   std::cout <<"\n";
// }
//
//  std::cout << "output 2" << '\n';
//  for(int i = 0; i < l_nY +2; i++){
//    std::cout << heightS1_db[i] << '\n';
//  }
