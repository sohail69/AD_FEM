#pragma once
#include "../macros.hpp"
#include "../../../../MFEM_STUFF/mfem-4.7/build/include/mfem/general/array.hpp"
#include "../../../../MFEM_STUFF/mfem-4.7/build/include/mfem/general/globals.hpp"
#include "../../../../MFEM_STUFF/mfem-4.7/build/include/mfem/general/mem_manager.hpp"
#include "../../../../MFEM_STUFF/mfem-4.7/build/include/mfem/general/device.hpp"

#include "mfem.hpp"

//
// templated vector numbers 
//
template<typename Numeric>
struct tVector{
  //Definition of templated vector
  mfem::Memory<Numeric> val;
  UINT64 size, Iter;

  tVector(){};
  tVector(int size_, mfem::MemoryType mt): val(size_, mt), size(size_){};
  ~tVector(){val.Delete();};


  /***************************************\
  !
  !  Vector Memory and allocation operations
  !
  \***************************************/
  /// Shortcut for mfem::Read(vec.GetMemory(), vec.Size(), on_dev).
  const Numeric *Read(bool on_dev = true) const
  { return mfem::Read(val, size, on_dev); }

  /// Shortcut for mfem::Read(vec.GetMemory(), vec.Size(), false).
  const Numeric *HostRead() const
  { return mfem::Read(val, size, false); }

  /// Shortcut for mfem::Write(vec.GetMemory(), vec.Size(), on_dev).
  Numeric *Write(bool on_dev = true)
  { return mfem::Write(val, size, on_dev); }

  /// Shortcut for mfem::Write(vec.GetMemory(), vec.Size(), false).
  Numeric *HostWrite()
  { return mfem::Write(val, size, false); }

  /// Shortcut for mfem::ReadWrite(vec.GetMemory(), vec.Size(), on_dev).
  Numeric *ReadWrite(bool on_dev = true)
  { return mfem::ReadWrite(val, size, on_dev); }

  /// Shortcut for mfem::ReadWrite(vec.GetMemory(), vec.Size(), false).
  Numeric *HostReadWrite()
  { return mfem::ReadWrite(val, size, false); }

   inline Numeric &operator[](int i) { return (*this)(i); };

   inline Numeric &operator()(int i) { return val[i]; };
};
