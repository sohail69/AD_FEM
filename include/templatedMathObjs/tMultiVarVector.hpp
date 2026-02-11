#pragma once
#include "../UtilityObjects/macros.hpp"
#include "../UtilityObjects/lowLevelMFEM.hpp"
#include "mfem.hpp"


/**
 This code is used for traversing densely populated multi-variate 
 local data with variable-tensor ranks. This is primarily for dealing
 with element level sampling to get Vars used for non-linear
 coefficients
**/


/** The residual calculation properties for different assembly types
|======================================================================================|
| CASE           |     Problem Characteristic    | Description                         |
|                |            Size               |                                     |
|================|===============================|=====================================|
|FULL    Assembly| nEQ                           | Store r-global directly             |
|----------------|-------------------------------|-------------------------------------|
|ELEMENT Assembly| nDOF*nELMs                    | Store r-Element directly            |
|----------------|-------------------------------|-------------------------------------|
|PARTIAL Assembly| Sdim*Edim*nCoeff*nQuads*nELMs | Store r-coefficients at the         |
|                |                               | integration-points without          !
|                |                               | multiplying by the basis function   |
|----------------|-------------------------------|-------------------------------------|

For non-linear problems
D = f(B.x);
**/



/**
 This part of the code build an extendable
 iterator struct for traversing mulit-tensor
 unstrtuctured environments, however not
 necessarily suitable for GPU, but useful for
 IO
**/

//The variable iterator
//data for problem size
template<typename uint>
struct Var{
  // TRanks   : Tensor rank of Variables
  // sizes    : Sizes of the tensor ranks (contiguous vector)
  uint ParentTrueVar, TRank;
  std::vector<uint> sizes;
};


//The variable iterator
//data for problem size
template<typename uint>
struct VarIterData{
  // TRanks   : Tensor rank of Variables [nVars]
  // sizes    : Sizes of the tensor ranks (contiguous vector) [nVars]
  // Soffsets : Offsets for the size vector of a Vars Tensor dimensions  [sum_{nVars}(Ranks)]
  // Voffsets : Offsets for the Vars starting points [nVars]
  // sDim     : Spatial dimension of the problem
  uint Tsize, sDim;
  std::vector<uint> TRanks, sizes, Soffsets, Voffsets;
};


//Clears the iterator of all its
//data for alternative usage
template<typename uint>
void clearIterator(VarIterData<uint> & data)
{
  data.Tsize=0;
  data.sDim=0;
  data.TRanks.clear();
  data.sizes.clear();
  data.Soffsets.clear();
  data.Voffsets.clear();
};


//Add data to the multi-variate iterator
//for multi-tensor objects
template<typename uint>
void AddVarIteratorDat(VarIterData<uint> & data, const std::vector<uint> TRank, const std::vector<uint> sizes)
{
  for(uint I=0; I<data.TRank.size(); I++) data.TRanks.push_back(TRank[I]);
  for(uint I=0; I<data.TRank.size(); I++){
    for(uint J=0; J<data.TRank[I]; J++){
      data.sizes.push_back(sizes[I]);
    }
  }
};

/**
 This part of the code passes the Iterator data
 collected from the IO class and passes it to a
 more GPU friendly class built for MFEM
**/

//The variable iterator data 
//for problem size using MFEM
//data layout
template<typename uint>
struct MFEMVarIterData{
  // TRanks   : Tensor rank of Variables [nVars]
  // sizes    : Sizes of the tensor ranks (contiguous vector) [nVars]
  // Soffsets : Offsets for the size vector of a Vars Tensor dimensions  [sum_{nVars}(Ranks)]
  // Voffsets : Offsets for the Vars starting points [nVars]
  // sDim     : Spatial dimension of the problem
  uint Tsize, sDim;
  mfem::Memory<uint> TRanks, sizes, Soffsets, Voffsets;
};

// Takes the base IO-iterator and passes
// it to a mfem-suitable data structure
template<typename uint>
void MakeMultiVarMFEMIter(const VarIterData<uint> & IterBase, MFEMVarIterData<uint> & IterMFEM)
{
  IterMFEM.Tsize = IterBase.Tsize;
  IterMFEM.sDim  = IterBase.sDim;
};


//Forward iterator for the
//multi-dimensional var-data
// [Almost exclusively used]
template<typename uint>
FORCE_INLINE uint MultiVarFwdIterator(const MFEMVarIterData<uint> data, const uint VarID, uint Iters[])
{
  uint Iter1D=0, Iter_tmp=0;
  #pragma unroll
  for(uint I=0; I < data.TRank[VarID]; I++){
    Iter_tmp = Iters[I];
    #pragma unroll
    for(uint J=I+1; J < data.TRank[VarID]; I++){
      Iter_tmp *= data.sizes[data.Soffsets[VarID] + J];
    }
    Iter1D += Iter_tmp;
  }
  return Iter1D + data.Voffsets[VarID];
};

//Inverse iterator for the
//multi-dimensional var-data
// [Unlikely to be ever used]
template<typename uint>
FORCE_INLINE void MultiVarInvIterator(const MFEMVarIterData<uint> data
                                    , const uint VarID
                                    , const uint vec_Iter
                                    , uint Iters[])
{
  uint Iter_tmp, Iter1D = vec_Iter - data.Voffsets[VarID];
  #pragma unroll
  for(uint I=0; I < data.TRank[VarID]; I++){
    Iter_tmp = Iter1D;
    #pragma unroll
    for(uint J=I+1; J < data.TRank[VarID]; I++){
      Iter_tmp /= data.sizes[data.Soffsets[VarID] + J];
    }
    Iters[I] = Iter_tmp;
  }
};
