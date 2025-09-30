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

// Implementation of data type mesh

#include <iostream>
#include <sstream>
#include <fstream>
#include <limits>
#include <cmath>
#include <cstring>
#include <ctime>
#include <functional>
#include <map>
#include <set>

#include "mfem/mesh/mesh_headers.hpp"
#include "mfem/fem/fem.hpp"
#include "mfem/general/sort_pairs.hpp"
#include "mfem/general/binaryio.hpp"
#include "mfem/general/text.hpp"
#include "mfem/general/device.hpp"
#include "mfem/general/tic_toc.hpp"
#include "mfem/general/gecko.hpp"
#include "mfem/fem/quadinterpolator.hpp"

#include "mesh.hpp"

namespace ppm
{

void Mesh::PrintGmsh(std::string &this_string)
{
   // Written by Innocent Niyonzima
   //==============================
   this_string = "";
   //std::string mesh_str;
   this_string = this_string + "$MeshFormat\n" + "2.2 0 8\n" + "$EndMeshFormat\n" + "$Nodes\n";
   this_string = this_string + to_string(GetNV()) + "\n";   
   for (int i = 0; i < GetNV(); i++)
   {
      const double *v = GetVertex(i);
      this_string = this_string + to_string(i+1) + " " + to_string(v[0]) + " " + to_string(v[1]) + " " + to_string(v[2]) + "\n";
   }

   this_string = this_string + "$EndNodes\n" + "$Elements\n" + to_string(GetNE()) + "\n";
   
   //int ele_attribute;
   //Element::Type ele_type;
   //int element_type;
   for (int i = 0; i < GetNE(); i++)
   {

      int int_quad = 1;

      const Element *ele = GetElement(i);
      int ele_attribute = ele->GetAttribute();
      Element::Type ele_type = ele->GetType();
      int element_type;
      
      if(ele_type == Element::SEGMENT)
      {
         element_type = 1;
      }
      else if(ele_type == Element::TRIANGLE)
      {
         element_type = 2;
      }
      else if(ele_type == Element::QUADRILATERAL)
      {
         element_type = 3;
      }
      else if(ele_type == Element::TETRAHEDRON)
      {
         element_type = 4;
      }
      else if(ele_type == Element::HEXAHEDRON)
      {
         element_type = 5;
      }
      else if(ele_type == Element::WEDGE)
      {
         element_type = 6;         
      }
      else if(ele_type == Element::PYRAMID)
      {
         element_type = 7;         
      }

      const int *v = ele->GetVertices();
      if (element_type == 2)
      {
         this_string = this_string + to_string(i+1) + " " + to_string(element_type) + " 2 " + to_string(ele_attribute) + " ";
         this_string = this_string + to_string(ele_attribute) + " " + to_string(v[0]+1) + " " + to_string(v[1]+1) + " " + to_string(v[2]+1) + "\n";
      }
      else if (element_type == 3)
      {
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

void Mesh::PrintGmsh(std::ostream &os)
{
   std::cout << "Start writing the mesh ..." << std::endl;

   // Written by Inno
   //================
   if(os.tellp() == 0) 
   {
      os << "$MeshFormat\n" 
      "2.2 0 8\n"
      "$EndMeshFormat\n";
   }
   os << "$Nodes\n";
   
   cout << setprecision(16);
   
   os << GetNV() << "\n";   
   for (int i = 0; i < GetNV(); i++)
   {
      const double *v = GetVertex(i);
      if(Dimension() == 1) { os << (i+1) << " " << v[0] << " " << 0 << " " << 0 << "\n"; }
      else if(Dimension() == 2) { os << (i+1) << " " << v[0] << " " << v[1] << " " << 0 << "\n"; }
      else { os << (i+1) << " " << v[0] << " " << v[1] << " " << v[2] << "\n"; }
   }
   os << "$EndNodes\n";
   os << "$Elements\n";
   int num_elts = GetNE();
   os << GetNE() << "\n";
   
   int ele_attribute;
   Element::Type ele_type;
   int element_type;
   for (int i = 0; i < GetNE(); i++)
   {
      const Element *ele = GetElement(i);
      const int *v = ele->GetVertices();

      ele_attribute = ele->GetAttribute();
      ele_type = ele->GetType();
      if(ele_type == Element::SEGMENT)
      {
         element_type = 1;
      }
      else if(ele_type == Element::TRIANGLE)
      {
         element_type = 2;
      }
      else if(ele_type == Element::QUADRILATERAL)
      {
         element_type = 3;
      }
      else if(ele_type == Element::TETRAHEDRON)
      {
         element_type = 4;
      }
      else if(ele_type == Element::HEXAHEDRON)
      {
         element_type = 5;
      }
      else if(ele_type == Element::WEDGE)
      {
         element_type = 6;         
      }
      else if(ele_type == Element::PYRAMID)
      {
         element_type = 7;         
      }
      //os << (i+1) << " " << element_type << " 2 " << ele_attribute << " 2 " << v[0]+1 << " " << v[1]+1 << " " << v[2]+1 << "\n";
      if (element_type == 2)
      {
         os << (i+1) << " " << element_type << " 2 " << ele_attribute << " " << ele_attribute << " " << v[0]+1 << " " << v[1]+1 << " " << v[2]+1 << "\n";
      }
      else if (element_type == 3)
      {
         os << (i+1) << " " << element_type << " 2 " << ele_attribute << " " << ele_attribute << " " << v[0]+1 << " " << v[1]+1 << " " << v[2]+1 << " " << v[3]+1 << "\n";
      }
   }
   os << "$EndElements\n";
   std::cout << "End writing the mesh ..." << std::endl;
   return;
}

}
