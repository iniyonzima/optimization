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

#ifndef PPM_DATACOLLECTION
#define PPM_DATACOLLECTION

#include "mfem/config/config.hpp"
#include "mfem/fem/gridfunc.hpp"
#include "mfem/fem/qfunction.hpp"
#include "mfem/fem/datacollection.hpp"
#include <string>
#include <map>
#include <fstream>

namespace ppm
{
// Begining of code added by inno
//===============================
// Helper class for Gmsh visualization data
class GmshDataCollection : public mfem::DataCollection
{
private:
   std::fstream gmsh_stream;

protected:
   void WriteGmshHeader(std::ostream &out);
   void WriteGmshFooter(std::ostream &out, const std::string &vtu_prefix);
   void SaveDataGmsh(std::ostream &out, int ref);
   void SaveGFieldGmsh(std::ostream& out, int ref_, const FieldMapIterator& it);

   std::string GenerateCollectionPath();
   //Mesh *ppm_mesh;

public:
   GmshDataCollection(const std::string& collection_name,
                          mfem::Mesh *mesh_ = NULL);
   virtual void Save() override;
   virtual void SaveMesh() override;
   void AddSolution(int time_step, double time_value, std::string this_name, std::string &this_string, int flag_data_type, int ncomp);
   inline void AddSolution(std::string this_name, std::string &this_string, int flag_data_type, int ncomp)
   {
      AddSolution(GetCycle(), GetTime(), this_name, this_string, flag_data_type, ncomp);
   }
   void Save(const std::string &field_name, std::string &solution_string);

   void AddSolution(int time_step, double time_value, std::string this_name, 
      std::ostringstream &this_ostringstream, int flag_data_type, int ncomp, int Flag_AMR, int Flag_Map_Per_TimeStep);
   inline void AddSolution(std::string this_name, std::ostringstream &this_stream, int flag_data_type, int ncomp, int flag_AMR = 0, int Flag_Map_Per_TimeStep = 0)
   {
      AddSolution(GetCycle(), GetTime(), this_name, this_stream, flag_data_type, ncomp, flag_AMR, Flag_Map_Per_TimeStep);
   }
   void Save(const std::string &field_name, std::ostringstream &solution_stream);

   void Create_Directory();

   void PrintGmsh(std::string &this_mesh);
   void PrintGmsh(std::ostream &os);

};
// End of code added by inno
//==========================

}
#endif
