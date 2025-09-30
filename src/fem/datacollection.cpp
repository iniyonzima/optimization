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



// #include "mfem/fem/fem.hpp"
// #include "../mesh/nurbs.hpp"
// #include "../mesh/vtk.hpp"
// #include "../general/binaryio.hpp"
// #include "../general/text.hpp"
// #include "picojson.h"

#include <cerrno>      // errno
#include <sstream>
#include <iostream>
#include <regex>

#include <string>
#include <iostream>

#ifndef _WIN32
#include <sys/stat.h>  // mkdir
#else
#include <direct.h>    // _mkdir
#define mkdir(dir, mode) _mkdir(dir)
#endif

#include<unistd.h>  
//#include<datacollection.h>  

#include "mfem.hpp"
// #include "mfem/fem/fem.hpp"
// #include "mfem/linalg/vector.hpp"
// #include "mfem/fem/datacollection.hpp"

// #include "mfem/mesh/nurbs.hpp"
// #include "mfem/mesh/vtk.hpp"
// #include "mfem/general/binaryio.hpp"
// #include "mfem/general/text.hpp"
// //#include "../../mfem-4.8/fem/picojson.h"
// //#include "../mesh/mesh.hpp"

#include "datacollection.hpp"

using namespace mfem;
using namespace std;

namespace ppm
{
   GmshDataCollection::GmshDataCollection(const std::string& collection_name, mfem::Mesh *mesh_)
      : mfem::DataCollection(collection_name, mesh_) {}

   std::string GmshDataCollection::GenerateCollectionPath()
   {
      return prefix_path + DataCollection::GetCollectionName();
   }

   void GmshDataCollection::Save()
   {
      for (FieldMapIterator it = field_map.begin(); it != field_map.end(); ++it)
      {
         std::string col_path = GenerateCollectionPath();
         int stat;
         stat = mkdir(col_path.c_str(), 0777);
         if (stat)
         {
            error = WRITE_ERROR;
            MFEM_WARNING("Error creating directory: " << col_path << " with error " << stat);
            return; // do not even try to write the mesh
         }
         std::string solution_name = (it->first);
         std::ofstream os(col_path + solution_name + ".msh");
         mfem::Mesh * mesh = GetMesh();
         PrintGmsh(os);      

         std::string solution_name_new = '"' + solution_name + '"';
         const char * scalar_view_name = solution_name_new.c_str();

         double time_value = GetTime();
         int time_step = GetCycle();
         int num_component = 1;
         int num_vertices = mesh->GetNV();
         os << "$NodeData\n" 
         << 1 <<'\n' 
         << scalar_view_name << '\n' 
         << 1 << '\n' 
         << time_value << '\n' 
         << 3 << '\n' 
         << time_step << '\n' 
         << num_component << '\n' 
         << num_vertices << '\n';

         GridFunction *this_gf = new GridFunction(*(GetField(it->first)));
         Vector *this_Vector = new Vector();
         
         int vdim = 1;
         this_gf->GetNodalValues(*(this_Vector), vdim);
         for (int i = 0; i < num_vertices; i++) 
         {
            int num_vertex = (i+1);
            double nodal_value = this_Vector->Elem(i);
            os << num_vertex << ' ' << nodal_value << '\n';
         }
         os << "$EndNodeData\n";
      }
      return; 
   }

   void GmshDataCollection::Create_Directory() 
   {
      std::string col_path = GenerateCollectionPath();
      std::cout << "col_path = " << col_path << std::endl;
      int stat;
      stat = mkdir(col_path.c_str(), 0777);
      if (stat)
      {
         error = WRITE_ERROR;
         MFEM_WARNING("Error creating directory: " << col_path << " with error " << stat);
         return; // do not even try to write the mesh
      }
   }

   void GmshDataCollection::Save(const std::string &this_field_name, std::string &solution_string) 
   {

      //MFEM_WARNING("Start saving the solution : " << this_field_name << "...");

      std::string col_path = GenerateCollectionPath();
      //std::cout << "col_path = " << col_path << std::endl;
      int status;

      struct stat info;
      if( stat(col_path.c_str(), &info) != 0 ){   
         MFEM_WARNING("Directory " << col_path << " does not exist. It is going to be created.");
         status = mkdir(col_path.c_str(), 0777);
      }
      if (status){
         error = WRITE_ERROR;
         MFEM_WARNING("Error creating directory: " << col_path << " with error " << status);
         return; // do not even try to write the mesh
      }
      std::string solution_name = (std::string) this_field_name;
      std::ofstream os((col_path + solution_name + ".msh"));
      os << solution_string;
      os.close();
      //MFEM_WARNING("End saving the solution : " << this_field_name << "...");
      return; 
   }

   void GmshDataCollection::AddSolution(int time_step, double time_value, std::string this_name, 
      std::string &this_string, int Flag_Data_Type, int ncomp)
   {
      // //std::cout << "Start adding the solution ("  << this_name << ") for the time step " << time_step << std::endl;
      // if(this_string.empty()){
      //    mfem::Mesh * mesh = GetMesh();
      //    PrintGmsh(this_string);
      //    //std::cout << "Writing the mesh for the field ("  << this_name << ")." << std::endl;
      // }
      int num_vertices = mesh->GetNV();
      int num_elements = mesh->GetNE();
      std::string solution_name_new = '"' + this_name + '"';
      const char * scalar_view_name = solution_name_new.c_str();

      std::setprecision(16);

      if(Flag_Data_Type == 1){
         int num_component = 1;
         this_string = this_string + "$NodeData\n" + "1\n" + scalar_view_name + "\n" + "1\n" + std::to_string(time_value) + "\n" + std::to_string(ncomp) + "\n";
         this_string = this_string + std::to_string(time_step) + "\n" + std::to_string(num_component) + "\n" + std::to_string(num_vertices) + "\n";

         GridFunction *this_gf = GetField(this_name.c_str());
         Vector *this_Vector = new Vector();
         
         int vdim = 1;
         this_gf->GetNodalValues(*(this_Vector), vdim);
         for (int i = 0; i < num_vertices; i++) {
            int num_vertex = (i+1);
            double nodal_value = this_Vector->Elem(i);
            this_string = this_string + std::to_string(num_vertex) + ' ' + std::to_string(nodal_value) + '\n';
         }
         this_string  = this_string + "$EndNodeData\n";
      }
      else if (Flag_Data_Type == 2){
         int num_component = 1;
         this_string = this_string + "$ElementNodeData\n" + "1\n" + scalar_view_name + "\n" + "1\n" + std::to_string(time_value) + "\n" + "3\n";
         this_string = this_string + std::to_string(time_step) + "\n" + std::to_string(num_component) + "\n" + std::to_string(num_elements) + "\n";

         GridFunction *this_gf = GetField(this_name.c_str());
         Vector *this_Vector = new Vector();
         
         int vdim = 1;
         this_gf->GetNodalValues(*(this_Vector), vdim);
         for (int i = 0; i < num_elements; i++) 
         {
            Element *this_element = mesh->GetElement(i);
            int *this_vertices_array = this_element->GetVertices();
            this_string = this_string + to_string(i+1) + ' ' + to_string(this_element->GetNVertices());
            for (int j = 0; j < 3 ; j++)
            {
               int num_vertex = this_vertices_array[j];
               long double nodal_value = this_Vector->Elem(num_vertex);            
               this_string = this_string + ' ' + std::to_string(nodal_value);
            }
            this_string = this_string + '\n';
         }
         this_string  = this_string + "$EndElementNodeData\n";      
      }
      else if (Flag_Data_Type == 3){
         int num_component = 3;
         this_string = this_string + "$ElementNodeData\n" + "1\n" + scalar_view_name + "\n" + "1\n" + std::to_string(time_value) + "\n" + "3\n";
         this_string = this_string + std::to_string(time_step) + "\n" + std::to_string(num_component) + "\n" + std::to_string(num_elements) + "\n";

         GridFunction *this_gf = GetField(this_name.c_str());
         Vector *this_Vector = new Vector();
         
         int vdim = 1;
         this_gf->GetNodalValues(*(this_Vector), vdim);
         for (int i = 0; i < num_elements; i++) 
         {
            Element *this_element = mesh->GetElement(i);
            int *this_vertices_array = this_element->GetVertices();
            this_string = this_string + to_string(i+1) + ' ' + to_string(this_element->GetNVertices());
            for (int j = 0; j < 3 ; j++)
            {
               int num_vertex = this_vertices_array[j];
               long double nodal_value = this_Vector->Elem(num_vertex);            
               this_string = this_string + ' ' + to_string(0.0) + ' ' + to_string(0.0) + ' ' + std::to_string(nodal_value);
            }
            this_string = this_string + '\n';
         }
         this_string  = this_string + "$EndElementNodeData\n";      
      }

      std::cout << "End adding the solution ("  << this_name << ") for the time step " << time_step << std::endl;

      return;
   }

   void GmshDataCollection::Save(const std::string &this_field_name, std::ostringstream &solution_stream) 
   {
      // std::cout << "Start saving the solution : " << this_field_name << "..." << std::endl;
      std::string col_path = GenerateCollectionPath();
      //std::cout << "col_path = " << col_path << std::endl;
      int status;

      struct stat info;
      
      if( stat(col_path.c_str(), &info) != 0 ) {   
         MFEM_WARNING("Directory " << col_path << " does not exist. It is created now...");
         status = mkdir(col_path.c_str(), 0777);
      }
      if (status) {
         error = WRITE_ERROR;
      }
      std::string solution_name = (std::string) this_field_name;
      std::ofstream os((col_path + solution_name + ".msh"));
      std::cout.setf( std::ios::fixed, std:: ios::scientific );
      os << solution_stream.str();
      os.close();
      // std::cout << "End saving the solution : " << this_field_name << "..." << std::endl;
      return; 
   }

   void GmshDataCollection::AddSolution(int time_step, double time_value, std::string this_name, 
      std::ostringstream &this_stream, int Flag_Data_Type, int ncomp, int Flag_AMR, int Flag_Map_Per_TimeStep)
   {
      // std::cout << "Start adding the solution " << this_name << " for the time step " << time_step << std::endl;
      {
         if(this_stream.tellp() == 0){
            // mfem::Mesh * mesh = GetMesh();
            //std::cout << "Start writing the unrefined mesh for the field (" << this_name << ")..." << std::endl;
            PrintGmsh(this_stream);
            //std::cout << "End writing the unrefined mesh for the field (" << this_name << ")..." << std::endl;
         }
         else {
            if(Flag_Map_Per_TimeStep != 0) 
            {
               this_stream.str() = "";
               // mfem::Mesh * mesh = GetMesh();
               //std::cout << "Start writing the refined mesh for the field (" << this_name << ")..." << std::endl;
               PrintGmsh(this_stream);
               //std::cout << "End writing the refined mesh for the field (" << this_name << ")..." << std::endl;
            }
         }
      }
      int num_vertices = mesh->GetNV();
      int num_elements = mesh->GetNE();
      std::string solution_name_new = '"' + this_name + '"';
      const char * scalar_view_name = solution_name_new.c_str();

      std::setprecision(16);

      if(Flag_Data_Type == 1)
      {
         int num_component = 1;
         this_stream << "$NodeData\n" << "1\n" << scalar_view_name << "\n" << "1\n" << time_value << "\n";
         this_stream << "3\n" << time_step << "\n" << num_component << "\n" << num_vertices << "\n";

         GridFunction *this_gf = GetField(this_name.c_str());
         Vector *this_Vector = new Vector();
         
         int vdim = 1;
         this_gf->GetNodalValues(*(this_Vector), vdim);
         for (int i = 0; i < num_vertices; i++) {
            int num_vertex = (i+1);
            double nodal_value = this_Vector->Elem(i);
            this_stream << num_vertex << ' ' << nodal_value << '\n';
         }
         this_stream << "$EndNodeData\n";

      }
      else if (Flag_Data_Type == 2){
         int num_component = 1;
         this_stream << "$ElementNodeData\n" << "1\n" << scalar_view_name << "\n" << "1\n" << time_value << "\n";
         this_stream << "3\n" << time_step << "\n" << num_component << "\n" << num_elements << "\n";

         GridFunction *this_gf = GetField(this_name.c_str());
         Vector *this_Vector = new Vector();
         
         int vdim = 1;
         this_gf->GetNodalValues(*(this_Vector), vdim);
         for (int i = 0; i < num_elements; i++) {
            Element *this_element = mesh->GetElement(i);
            int *this_vertices_array = this_element->GetVertices();
            this_stream << (i+1) << ' ' << this_element->GetNVertices();
            for (int j = 0; j < 3 ; j++){
               int num_vertex = this_vertices_array[j];
               long double nodal_value = this_Vector->Elem(num_vertex);            
               this_stream << ' ' << nodal_value;
            }
            this_stream << '\n';
         }
         this_stream << "$EndElementNodeData\n";      
      }
      else if (Flag_Data_Type == 3) {
         int num_component = 3;
         this_stream << "$ElementNodeData\n" << "1\n" << scalar_view_name << "\n" << "1\n" << time_value << "\n";
         this_stream  << "3\n" << time_step << "\n" << num_component << "\n" << num_elements << "\n";

         GridFunction *this_gf = GetField(this_name.c_str());
         Vector *this_Vector = new Vector();
         
         int vdim = 1;
         this_gf->GetNodalValues(*(this_Vector), vdim);
         for (int i = 0; i < num_elements; i++) {
            Element *this_element = mesh->GetElement(i);
            int *this_vertices_array = this_element->GetVertices();
            this_stream << (i+1) << ' ' << this_element->GetNVertices();
            for (int j = 0; j < 3 ; j++){
               int num_vertex = this_vertices_array[j];
               long double nodal_value = this_Vector->Elem(num_vertex);            
               this_stream << ' ' << 0.0 << ' ' << 0.0 << ' ' << nodal_value;
            }
            this_stream << '\n';
         }
         this_stream << "$EndElementNodeData\n";      
      }
      // std::cout << "End adding the solution "  << this_name << " for the time step " << time_step << std::endl;
      return;
   }

   void GmshDataCollection::SaveMesh()
   {
      std::cout << "function member GmshDataCollection::SaveMesh() is under construction" << std::endl;
   }

   void GmshDataCollection::WriteGmshHeader(std::ostream &os)
   {
      std::cout << "function member GmshDataCollection::WriteGmshHeader() is under construction" << std::endl;
   }

   void GmshDataCollection::WriteGmshFooter(std::ostream &os, const std::string &vtu_prefix)
   {
      std::cout << "function member GmshDataCollection::WriteGmshFooter() is under construction" << std::endl;
   }

   void GmshDataCollection::SaveDataGmsh(std::ostream &os, int ref)
   {
      std::cout << "function member GmshDataCollection::SaveDataGmsh() is under construction" << std::endl;
   }

   void GmshDataCollection::SaveGFieldGmsh(std::ostream &os, int ref_, const FieldMapIterator &it)
   {
      std::cout << "function member GmshDataCollection::SaveGFieldGmsh() is under construction" << std::endl;
   }

   void GmshDataCollection::PrintGmsh(std::string &this_string)
   {
      mfem::Mesh *mesh = GetMesh();
      this_string = "";
      this_string = this_string + "$MeshFormat\n" + "2.2 0 8\n" + "$EndMeshFormat\n" + "$Nodes\n";
      this_string = this_string + to_string(mesh->GetNV()) + "\n";   

      for (int i = 0; i < mesh->GetNV(); i++){
         const double *v = mesh->GetVertex(i);
         this_string = this_string + to_string(i+1) + " " + to_string(v[0]) + " " + to_string(v[1]) + " " + to_string(v[2]) + "\n";
      }

      this_string = this_string + "$EndNodes\n" + "$Elements\n" + to_string(mesh->GetNE()) + "\n";
      
      for (int i = 0; i < mesh->GetNE(); i++){
         int int_quad = 1;
         const mfem::Element *ele = mesh->GetElement(i);
         int ele_attribute = ele->GetAttribute();
         mfem::Element::Type ele_type = ele->GetType();
         int element_type;
         if(ele_type == mfem::Element::SEGMENT){
            element_type = 1;
         }
         else if(ele_type == mfem::Element::TRIANGLE){
            element_type = 2;
         }
         else if(ele_type == mfem::Element::QUADRILATERAL){
            element_type = 3;
         }
         else if(ele_type == mfem::Element::TETRAHEDRON){
            element_type = 4;
         }
         else if(ele_type == mfem::Element::HEXAHEDRON){
            element_type = 5;
         }
         else if(ele_type == mfem::Element::WEDGE){
            element_type = 6;         
         }
         else if(ele_type == mfem::Element::PYRAMID){
            element_type = 7;         
         }

         const int *v = ele->GetVertices();
         if (element_type == 2){
            this_string = this_string + to_string(i+1) + " " + to_string(element_type) + " 2 " + to_string(ele_attribute) + " ";
            this_string = this_string + to_string(ele_attribute) + " " + to_string(v[0]+1) + " " + to_string(v[1]+1) + " " + to_string(v[2]+1) + "\n";
         }
         else if (element_type == 3){
            if(int_quad < 100)
            {
               std::cout << "v[0]+1 = " << v[0]+1 << ", v[1]+1 = " << v[1]+1 << ", v[2]+1 = " << v[2]+1 << ", v[3]+1 = " << v[3]+1 << std::endl; 
               int_quad = int_quad + 1;
            }
            this_string = this_string + to_string(i+1) + " " + to_string(element_type) + " 2 " + to_string(ele_attribute) + " ";
            this_string = this_string + to_string(ele_attribute) + " " + to_string(v[0]+1) + " " + to_string(v[1]+1) + " " + to_string(v[2]+1) + " " + to_string(v[3]+1) + "\n";
         }
      }
      this_string = this_string + "$EndElements\n";
      return;
   }
   void GmshDataCollection::PrintGmsh(std::ostream &os)
   {
      //std::cout << "Start writing the mesh ..." << std::endl;
      if(os.tellp() == 0) {
         os << "$MeshFormat\n" 
         "2.2 0 8\n"
         "$EndMeshFormat\n";
      }
      os << "$Nodes\n";
      
      std::cout << std::setprecision(16);
      mfem::Mesh *mesh = GetMesh();
      os << mesh->GetNV() << "\n";   
      for (int i = 0; i < mesh->GetNV(); i++){
         const double *v = mesh->GetVertex(i);
         if(mesh->Dimension() == 1) { os << (i+1) << " " << v[0] << " " << 0 << " " << 0 << "\n"; }
         else if(mesh->Dimension() == 2) { os << (i+1) << " " << v[0] << " " << v[1] << " " << 0 << "\n"; }
         else { os << (i+1) << " " << v[0] << " " << v[1] << " " << v[2] << "\n"; }
      }
      os << "$EndNodes\n";
      os << "$Elements\n";
      // int num_elts = mesh->GetNE();
      os << mesh->GetNE() << "\n";
      
      int ele_attribute;
      mfem::Element::Type ele_type;
      int element_type;
      for (int i = 0; i < mesh->GetNE(); i++){
         const mfem::Element *ele = mesh->GetElement(i);
         const int *v = ele->GetVertices();

         ele_attribute = ele->GetAttribute();
         ele_type = ele->GetType();
         if(ele_type == Element::SEGMENT){
            element_type = 1;
         }
         else if(ele_type == Element::TRIANGLE){
            element_type = 2;
         }
         else if(ele_type == Element::QUADRILATERAL){
            element_type = 3;
         }
         else if(ele_type == Element::TETRAHEDRON){
            element_type = 4;
         }
         else if(ele_type == Element::HEXAHEDRON){
            element_type = 5;
         }
         else if(ele_type == Element::WEDGE){
            element_type = 6;         
         }
         else if(ele_type == Element::PYRAMID){
            element_type = 7;         
         }
         if (element_type == 2){
            os << (i+1) << " " << element_type << " 2 " << ele_attribute << " " << ele_attribute << " " << v[0]+1 << " " << v[1]+1 << " " << v[2]+1 << "\n";
         }
         else if (element_type == 3){
            os << (i+1) << " " << element_type << " 2 " << ele_attribute << " " << ele_attribute << " " << v[0]+1 << " " << v[1]+1 << " " << v[2]+1 << " " << v[3]+1 << "\n";
         }
      }
      os << "$EndElements\n";
      //std::cout << "End writing the mesh ..." << std::endl;
      return;
   }
}  // end namespace ppm
