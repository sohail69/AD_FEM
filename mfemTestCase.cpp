#include <fstream>
#include <iostream>
#include <vector>
#include <iomanip>
#include <cmath>
#include "mfem.hpp"
#include "include/nlOperator/tADNonLinearForm.hpp"
#include "include/UtilityObjects/Visualisation.hpp"

// Mathematical objects
#include "include/templatedMathObjs/dualNumber.hpp"
#include "include/templatedMathObjs/tVector.hpp"
#include "include/templatedMathObjs/tMultiVarVector.hpp"

// Mathematical Functions
#include "include/templatedMaths/tCmath.hpp"


int main(){
  // 1. Initialize MPI and HYPRE.
  //    and Parse command-line options
  Mpi::Init();
  const int myid = Mpi::WorldRank();
  int ref_levels=-1, order=2;
  const char *mesh_file = "data/star.mesh";
  const char *device_config = "cpu";
  bool use_dev=false;
  mfem::Device device(device_config);
  mfem::MemoryType mt = device.GetMemoryType();

  // 2. Read the mesh from the given mesh file, and refine once uniformly.
  Mesh mesh(mesh_file);
  int dim = mesh.Dimension();
  if (ref_levels == -1) ref_levels = (int)floor(log(10000./mesh.GetNE())/log(2.)/dim);
  for (int l = 0; l < ref_levels; l++) mesh.UniformRefinement();
  ParMesh pmesh(MPI_COMM_WORLD, mesh);

  // 3. Define a finite element space on the mesh. Here we use H1 continuous
  //    high-order Lagrange finite elements of the given order.
  H1_FECollection fec(order, dim);
  ParFiniteElementSpace fespace(&pmesh, &fec);
  std::vector<mfem::ParGridFunction*> gFuncs;
  std::vector<std::string>            FieldNames;
  gFuncs.push_back(new mfem::ParGridFunction(&fespace)); FieldNames.push_back("field_1");
  gFuncs.push_back(new mfem::ParGridFunction(&fespace)); FieldNames.push_back("field_2");

  int NEQs=0;
  for(int I=0; I<gFuncs.size(); I++) NEQs += gFuncs[I]->ParFESpace()->GetTrueVSize();
  mfem::Vector x(NEQs,mt), y(NEQs,mt); 

  //Test case for the non-linear
  //form (make sure it compiles 
  //and runs)
  tADNLForm<mfem::real_t> nlProb(gFuncs, device, mt, use_dev);
  nlProb.buildJacobian(x);
  nlProb.Mult(x,y);

  // 5. Output the vector data
  //    into paraview 
  ParaViewVisualise("testNLProblem",gFuncs,FieldNames,order,&pmesh,0.00);

  // Delete the objects
  // an clean-up
  for(int I=0; I<gFuncs.size(); I++) delete gFuncs[I];
  gFuncs.clear();
  return 0;
};
