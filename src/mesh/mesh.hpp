// Copyright (c) 2010-2022, Lawrence Livermore National Security, LLC. Produced
// at the Lawrence Livermore National Laboratory. All Rights reserved. See files
// LICENSE and NOTICE for details. LLNL-CODE-806117.
//
// This file is part of the MFEM library. For more information and source code
// availability visit https://mfem.org.
//
// MFEM is free software; you can redistribute it and/or modify it under the
// terms of the BSD-3 license. We welcome feedback and contributions, see file
// CONTRIBUTING.md for details.

#ifndef PPM_MESH
#define PPM_MESH

#include "mfem/config/config.hpp"
#include "mfem/general/stable3d.hpp"
#include "mfem/general/globals.hpp"
#include "mfem/mesh/triangle.hpp"
#include "mfem/mesh/tetrahedron.hpp"
#include "mfem/mesh/vertex.hpp"
#include "mfem/mesh/vtk.hpp"
#include "mfem/mesh/ncmesh.hpp"
#include "mfem/fem/eltrans.hpp"
#include "mfem/fem/coefficient.hpp"
#include "mfem/general/zstr.hpp"
#ifdef MFEM_USE_ADIOS2
#include "mfem/general/adios2stream.hpp"
#endif
#include <iostream>

namespace ppm
{

class Mesh : public mfem::Mesh
{
public:
   void PrintGmsh(std::ostream &os);
   void PrintGmsh(std::string &this_mesh);
};

#endif
