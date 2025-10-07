//                                MFEM Example 0
//
// Compile with: make my_capteur_magsta
//
//
// Description: This example code demonstrates the most basic usage of MFEM to
//              define a simple finite element discretization of the Laplace
//              problem -Delta u = 1 with zero Dirichlet boundary conditions.
//              General 2D/3D mesh files and finite element polynomial degrees
//              can be specified by command line options.

#include "mfem.hpp"
#include <fstream>
#include <iostream>
#include <cmath>
#include <stdio.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <filesystem>
#include <string>
#include <math.h>
#include <algorithm>    // std::min

#ifdef OPTIMIZATION_USE_MPI
   #include <mpi.h>
#endif

#include "../../src/fem/datacollection.hpp"

const char *gd_path = "/home/innocent/softwares/src/optimization/build/bin/Gmsh_Streamers/";
const char *mesh_file = "/home/innocent/softwares/src/ppm/data/mesh/capteur/capteur.msh";

using namespace std;
using namespace mfem;
using namespace ppm;

const double PI = 3.141592653589793238463;

class MagneticFlux
{
private:
   ParGridFunction az_;
   PWConstCoefficient j0_;
public:
   MagneticFlux( ParGridFunction az, PWConstCoefficient j0)
      : az_(az), j0_(j0) {}
   double ComputeFlux() 
   {
      ParFiniteElementSpace *FESpace_tmp = dynamic_cast<ParFiniteElementSpace *>(const_cast<FiniteElementSpace *>(az_.FESpace()));
      if (!FESpace_tmp) {
         std::cerr << "FESpace is not a ParFiniteElementSpace!" << std::endl;
         std::abort(); // or throw
      }
      ParLinearForm j0_lf(FESpace_tmp);
      j0_lf.AddDomainIntegrator(new DomainLFIntegrator(j0_));
      j0_lf.Assemble();
      return 0.5 * j0_lf(az_);
   }
};

class MagSta_Problem
{
private:
   ParFiniteElementSpace &fes;
   Coefficient *lhs_coeff;
   Coefficient *rhs_coeff;
   Coefficient *design_char_coeff;
   ParBilinearForm *K;
   ParLinearForm *F;
public:
   MagSta_Problem(ParFiniteElementSpace &fes_, Coefficient *lhs_coeff_, Coefficient *rhs_coeff_) : fes(fes_), 
   lhs_coeff(lhs_coeff_), rhs_coeff(rhs_coeff_), design_char_coeff(NULL)
   {
      K = new ParBilinearForm(&fes);
      K->AddDomainIntegrator(new DiffusionIntegrator(*lhs_coeff));
      F = new ParLinearForm(&fes);
      F->AddDomainIntegrator(new DomainLFIntegrator(*rhs_coeff));
   }
   MagSta_Problem(ParFiniteElementSpace &fes_, Coefficient *lhs_coeff_, Coefficient *rhs_coeff_, Coefficient *design_char_coeff_) : fes(fes_), 
   lhs_coeff(lhs_coeff_), rhs_coeff(rhs_coeff_), design_char_coeff(design_char_coeff_)
   {
      K = new ParBilinearForm(&fes);
      F = new ParLinearForm(&fes);
   }

   Coefficient *get_LHS_Coeff()
   {
      return lhs_coeff;
   }
   void set_LHS_Coeff(Coefficient *This_Coeff)
   {
      lhs_coeff = This_Coeff;
   }
   void set_Char_Coeff(Coefficient *This_Char_Coeff)
   {
      design_char_coeff = This_Char_Coeff;
   }

   void UpdateCoefficient(Coefficient &h_gf_coeff, int pow)
   {
      PowerCoefficient *h_pow = new PowerCoefficient(h_gf_coeff, pow);
      double mu_0 = 4 * PI * 1.0e-7, nu_0 = 1.0/mu_0, mu_r = 1000, nu_r = 1.0/mu_r, nu_fe = nu_r * nu_0, nu_rel = (nu_fe - nu_0);
      ConstantCoefficient *nu_0_c = new ConstantCoefficient(nu_0);
      ConstantCoefficient *nu_rel_c = new ConstantCoefficient(nu_rel);
      ProductCoefficient *nu_hat_fe_c = new ProductCoefficient(*nu_rel_c, *h_pow);
      SumCoefficient *nu_hat = new SumCoefficient(*nu_0_c, *nu_hat_fe_c);
      set_LHS_Coeff(nu_hat);
      // return nu_hat;
      return;
   }

   void ComputeBilinearForm(int Flag_Update_Matrix)
   {
      BilinearFormIntegrator *integ = new DiffusionIntegrator(*lhs_coeff);
      ParBilinearForm *K_Temp = new ParBilinearForm(&fes);
      if(Flag_Update_Matrix)
      {
         delete K;
         K = new ParBilinearForm(&fes, K_Temp);
         K->AddDomainIntegrator(integ);
      }
      else 
      {
         K->AddDomainIntegrator(integ);
      }
      K->Assemble();
      delete integ; 
      delete K_Temp;
      return;
   }

   void ComputeLinearForm(int Flag_Update_RHS) 
   {
      if(Flag_Update_RHS)
      {
         if(F)
         { 
            delete F; 
         }
         F = new ParLinearForm(&fes);
         F->AddDomainIntegrator(new DomainLFIntegrator(*rhs_coeff));
      }
      else 
      {
         F->AddDomainIntegrator(new DomainLFIntegrator(*rhs_coeff));
      }
      F->Assemble();
      return;
   }

   void Solve( ParGridFunction  &x, Array<int> ess_tdof_list) 
   {
      // if(0)
      // {
      //    SparseMatrix A;
      //    Vector B, X;
      //    K->FormLinearSystem(ess_tdof_list, x, *F, A, X, B);

      //    #ifndef MFEM_USE_SUITESPARSE
      //       DSmoother M(A);         
      //       PCG(A, M, B, X, 0, 5000, 1e-12, 0.0);
      //    #else
      //       UMFPackSolver umf_solver;
      //       umf_solver.Control[UMFPACK_ORDERING] = UMFPACK_ORDERING_METIS;
      //       umf_solver.SetOperator(A);
      //       umf_solver.Mult(B, X);
      //    #endif
      //    K->RecoverFEMSolution(X, *F, x);
      // }
      // else
      // {
         HypreParMatrix A;
         Vector B, X;
         K->FormLinearSystem(ess_tdof_list, x, *F, A, X, B);

         // 11. Solve the system using PCG with hypre's BoomerAMG preconditioner.
         HypreBoomerAMG P(A);
         // BiCGSTABSolver cg(MPI_COMM_WORLD);
         CGSolver cg(MPI_COMM_WORLD);
         cg.SetRelTol(1e-12);
         cg.SetMaxIter(2000);
         cg.SetPrintLevel(0);
         cg.SetPreconditioner(P);
         cg.SetOperator(A);
         cg.Mult(B, X);

         // Reconstruct total solution
         //===========================
         K->RecoverFEMSolution(X, *F, x);

      // }
      return;
   }

   double ComputeFlux( ParGridFunction  &x) 
   {
      ParLinearForm j0_lf(&fes);
      j0_lf.AddDomainIntegrator(new DomainLFIntegrator(*rhs_coeff) );
      j0_lf.Assemble();
      return j0_lf(x);
   }

   double ComputeFlux( ParGridFunction  &x, PWConstCoefficient &j0) 
   {
      MagneticFlux *MagFlux = new MagneticFlux(x, j0);
      return MagFlux->ComputeFlux();
   }

   double ComputeVolume(Coefficient &h_GF_Coeff) 
   {
      ParLinearForm j0_lf(&fes);
      ConstantCoefficient lambda_coeff(1);
      j0_lf.AddDomainIntegrator(new DomainLFIntegrator(lambda_coeff) );
      ParGridFunction *vol_GF = new ParGridFunction(&fes);
      vol_GF->ProjectCoefficient(h_GF_Coeff);
      j0_lf.Assemble();
      double vol_val = j0_lf(*vol_GF);
      delete vol_GF;
      return vol_val;
   }

   GridFunctionCoefficient ConvertToGridFunctionCoefficient(Coefficient &this_coeff)
   {
      ParGridFunction *gf_tmp = new ParGridFunction(&fes);
      gf_tmp->ProjectCoefficient(this_coeff);
      GridFunctionCoefficient this_gf_coeff(gf_tmp);
      return this_gf_coeff;
   }

   double Evaluate_Objective_Function( ParGridFunction  &x, double lambda, Coefficient &h_GF_Coeff) 
   {
      double flux, volume;
      flux = ComputeFlux(x);
      volume = ComputeVolume(h_GF_Coeff);
      std::cout << "Done in the function 'Evaluate_Objective_Function'. The flux = " << flux << ". The volume = " << volume << std::endl; 
      double obj_func = flux + lambda * volume;
      return obj_func;
   }

   void PrintGridFunction_UsingGmsh(GmshDataCollection *_gd_magsta, int _flag_data_type, int _ncomp)
   {
      std::ostringstream solution_name_x;
      std::ostringstream solution_name_x_adj;

      //int ncomp = 1;
      int time_step_gmsh = 0;
      double time_gmsh = 0.0;
      {
         cout << " GmshDataCollection step = " << time_step_gmsh << ", t = " << time_gmsh << endl;
         _gd_magsta->SetCycle (time_step_gmsh);
         _gd_magsta->SetTime(time_gmsh);
         int flag_data_type = 1;
         int ncomp = 1;

         std::string name_field_direct = "az_dir";
         std::string name_field_adjoint = "az_adj";

         _gd_magsta->AddSolution(name_field_direct, solution_name_x, flag_data_type, ncomp);
         std::cout << "Saving the solution " << "az...\n" << std::endl;
         _gd_magsta->Save(name_field_direct, solution_name_x);
         std::cout << "End saving the solution " << "az...\n" << std::endl;

         _gd_magsta->AddSolution(name_field_adjoint, solution_name_x_adj, flag_data_type, ncomp);
         std::cout << "Saving the solution " << "az_adj...\n" << std::endl;
         _gd_magsta->Save(name_field_adjoint, solution_name_x_adj);
         std::cout << "End saving the solution " << "az_adj...\n" << std::endl;
      }
      return;
   }

   void PrintGridFunction_UsingGmsh(GmshDataCollection *_gd_magsta, ParGridFunction &this_GF, string old_name, string new_name) 
   {
      _gd_magsta->DeregisterField(old_name);
      _gd_magsta->RegisterField(new_name, &this_GF);

      std::ostringstream solution_name_x;
      std::ostringstream solution_name_x_adj;

      int time_step_gmsh = 0;
      double time_gmsh = 0.0;
      {
         cout << " GmshDataCollection step = " << time_step_gmsh << ", t = " << time_gmsh << endl;
         _gd_magsta->SetCycle (time_step_gmsh);
         _gd_magsta->SetTime(time_gmsh);
         int flag_data_type = 1;
         int ncomp = 1;

         _gd_magsta->AddSolution(new_name, solution_name_x, flag_data_type, ncomp);
         std::cout << "Saving the solution " << new_name << "\n" << std::endl;
         _gd_magsta->Save(new_name, solution_name_x);
         std::cout << "End saving the solution " << new_name << "\n" << std::endl;
      }
      return;
   }

   void PrintGridFunctionCoefficient_UsingGmsh(GmshDataCollection *_gd_magsta_gfc, 
      GridFunctionCoefficient &gf_coeff, string name_field, int _flag_data_type, int _ncomp)
   {
      std::ostringstream solution_name_x;
      std::ostringstream solution_name_x_adj;

      int time_step_gmsh = 0;
      double time_gmsh = 0.0;
      ParFiniteElementSpace *FESpace_tmp = dynamic_cast<ParFiniteElementSpace *>(const_cast<FiniteElementSpace *>(gf_coeff.GetGridFunction()->FESpace()));
      ParGridFunction *gf_tmp = dynamic_cast<ParGridFunction *>(const_cast<GridFunction *>(gf_coeff.GetGridFunction()));
      if ((!FESpace_tmp) || (!gf_tmp)) {
         // Handle error: the cast failed
         std::cerr << "FESpace is not a ParFiniteElementSpace! of gf_tmp is not ParGridFunction! " << std::endl;
         exit(1); // or throw
      }

      // ParGridFunction *x_tmp = new ParGridFunction(FESpace_tmp, *(gf_coeff.GetGridFunction()));
      ParGridFunction *x_tmp = new ParGridFunction(FESpace_tmp, *(gf_tmp));
      _gd_magsta_gfc->RegisterField(name_field, x_tmp);
      {
         cout << " GmshDataCollection step = " << time_step_gmsh << ", t = " << time_gmsh << endl;
         _gd_magsta_gfc->SetCycle (time_step_gmsh);
         _gd_magsta_gfc->SetTime(time_gmsh);
         int flag_data_type = 1;
         int ncomp = 1;
         _gd_magsta_gfc->AddSolution(name_field, solution_name_x, flag_data_type, ncomp);
         std::cout << "Saving the solution ... \n" << std::endl;
         _gd_magsta_gfc->Save(name_field, solution_name_x);
         std::cout << "End saving the solution ... \n" << std::endl;
      }
      _gd_magsta_gfc->DeregisterField(name_field);
      delete x_tmp;
      return;
   }

   void PrintGridFunctionCoefficient_UsingGmsh(GmshDataCollection *_gd_magsta_gfc, 
      Coefficient &this_coeff, string name_field, int _flag_data_type, int _ncomp)
   {
      ParGridFunction *gf_tmp = new ParGridFunction(&fes);
      gf_tmp->ProjectCoefficient(this_coeff);
      GridFunctionCoefficient this_gf_coeff(gf_tmp);
      PrintGridFunctionCoefficient_UsingGmsh(_gd_magsta_gfc, this_gf_coeff, name_field, _flag_data_type, _ncomp);
      return;
   }

   void PrintGridFunctionCoefficient_UsingParaview(ParaViewDataCollection *_pd_magsta_gfc, Coefficient &h_gf_coeff, string name_field)
   {
      std::ostringstream solution_name_x;
      std::ostringstream solution_name_x_adj;

      int time_step_pv = 0;
      double time_pv = 0.0;

      ParGridFunction *x_tmp = new  ParGridFunction (&fes);
      x_tmp->ProjectCoefficient(h_gf_coeff);
      _pd_magsta_gfc->RegisterField(name_field, x_tmp);
      {
         cout << " GmshDataCollection step = " << time_step_pv << ", t = " << time_pv << endl;
         _pd_magsta_gfc->SetCycle (time_step_pv);
         _pd_magsta_gfc->SetTime(time_pv);
         _pd_magsta_gfc->Save();
      }
     _pd_magsta_gfc->DeregisterField(name_field);
     delete x_tmp;
     return;
   }

   void ComputeSensitivity( ParGridFunction  &dx, ParGridFunction &dx_adj, Coefficient &nu_hat, ParGridFunction &sensitivity_gridfunction)
   {
      VectorGridFunctionCoefficient grad_x(&dx);
      VectorGridFunctionCoefficient grad_x_adj(&dx_adj);
      InnerProductCoefficient *sensitivity_unsigned = new InnerProductCoefficient(grad_x, grad_x_adj);
      ProductCoefficient *sensitivity_unsigned_2 = new ProductCoefficient(*sensitivity_unsigned, nu_hat);
      ConstantCoefficient minus_one(-1.0);

      ProductCoefficient sensitivity_pc_minus(minus_one, *sensitivity_unsigned_2);
      ProductCoefficient sensitivity_pc(sensitivity_pc_minus, *design_char_coeff);

      sensitivity_gridfunction.ProjectCoefficient(sensitivity_pc);
      delete sensitivity_unsigned;
      delete sensitivity_unsigned_2;
      return;
   }

   ~MagSta_Problem() 
   {
      delete K;
      delete F;
   }
};


int main(int argc, char *argv[])
{
   // 1. Initialize MPI and HYPRE
   //============================
   Mpi::Init(argc, argv);
   int num_procs = Mpi::WorldSize();
   int myid = Mpi::WorldRank();
   Hypre::Init();

   // 1. Parse command line options
   //==============================
   int order = 1;
   OptionsParser args(argc, argv);
   args.AddOption(&mesh_file, "-m", "--mesh", "Mesh file to use.");
   args.AddOption(&order, "-o", "--order", "Finite element polynomial degree");

   // args.ParseCheck();

   args.Parse();
   if (!args.Good()) {
      if (myid == 0) {
         args.PrintUsage(cout);
      }
      MPI_Finalize();
      return 1;
   }
   if (myid == 0) {
      mfem::out << num_procs << " number of process created.\n";
      args.PrintOptions(cout);
   }

   // 2. Read the mesh from the given mesh file, and refine once uniformly.
   //======================================================================
   Mesh *mesh = new Mesh(mesh_file);
   // mesh->UniformRefinement(1);   
   cout << "There are " << mesh->attributes.Size() << " physical region in the mesh" << endl;
   mesh->Finalize();

   // 5. Define a parallel mesh by partitioning the serial mesh.  Once the
   //    parallel mesh is defined, the serial mesh can be deleted.
   ParMesh *pmesh = new ParMesh(MPI_COMM_WORLD, *mesh);
   delete mesh;


   // 3. Define a finite element space on the mesh. Here we use H1 continuous
   //    high-order Lagrange finite elements of the given order.
   //========================================================================
   H1_FECollection fec(order, pmesh->Dimension());
   ParFiniteElementSpace fespace(pmesh, &fec);

   int order_vol = 0;
   L2_FECollection fec_vol(order_vol, pmesh->Dimension());
   ParFiniteElementSpace fespace_vol(pmesh, &fec_vol);

   H1_FECollection fec_adj(order, pmesh->Dimension());
   ParFiniteElementSpace fespace_adj(pmesh, &fec_adj);

   ND_FECollection fec_nd(order, pmesh->Dimension());
   ParFiniteElementSpace fespace_nd(pmesh, &fec_nd);

   ND_FECollection fec_nd_adj(order, pmesh->Dimension());
   ParFiniteElementSpace fespace_nd_adj(pmesh, &fec_nd_adj);

   // 4. Extract the list of all the boundary DOFs. These will be marked as
   //    Dirichlet in order to enforce zero boundary conditions.  
   //====================================================================== 
   Array<int> ess_tdof_list;
   Array<int> ess_bdr;
   if (pmesh->bdr_attributes.Size())
   {
      ess_bdr.SetSize(pmesh->bdr_attributes.Max());
      ess_bdr = 0;
      ess_bdr[999] = 1;
      fespace.GetEssentialTrueDofs(ess_bdr, ess_tdof_list);
   }
   Vector x_Diri(pmesh->bdr_attributes.Max() );
   x_Diri = 0.0;
   x_Diri(999) = 0.0;
   PWConstCoefficient x_Diri_c(x_Diri);

   Array<int> ess_tdof_list_adj;
   Array<int> ess_bdr_adj;
   if (pmesh->bdr_attributes.Size())
   {
      ess_bdr_adj.SetSize(pmesh->bdr_attributes.Max());
      ess_bdr_adj = 0;
      ess_bdr_adj[999] = 1;
      fespace_adj.GetEssentialTrueDofs(ess_bdr_adj, ess_tdof_list_adj);
   }
   Vector x_Diri_adj(pmesh->bdr_attributes.Max() );
   x_Diri_adj = 0.0;
   x_Diri_adj(999) = 0.0;
   PWConstCoefficient x_Diri_c_adj(x_Diri_adj);

   // 6. Define the coefficient of the RHS used in the LForm of the direct problem.
   //==============================================================================
   Vector js(pmesh->attributes.Max());
   js = 0.0;
   js(1999) = M_PI * 1.0e8;
   PWConstCoefficient js_c(js);

   // 7. Define the coefficient of the RHS used in the LForm of the adjoint problem.
   //===============================================================================
   double surf_indu;
   surf_indu = 8e-6;
   Vector j0_vec(pmesh->attributes.Max());
   j0_vec = 0.0;
   j0_vec(2000) = -1.0/surf_indu;
   j0_vec(2001) = 1.0/surf_indu;
   j0_vec(2002) = 1.0/surf_indu;
   j0_vec(2003) = -1.0/surf_indu;
   PWConstCoefficient j0_coeff(j0_vec);
   ParGridFunction j0_gf(&fespace);
   j0_gf.ProjectCoefficient(j0_coeff);

   // 8. Define the solution x as a finite element grid function in fespace. Set
   //    the initial guess to zero, which also sets the boundary conditions.
   //===========================================================================
   ParGridFunction x(&fespace);
   x = 0.0;
   x.ProjectBdrCoefficient(x_Diri_c, ess_bdr);

   ParGridFunction dx(&fespace_nd);
   dx = 0.0;
   DiscreteLinearOperator grad(&fespace, &fespace_nd);
   grad.AddDomainInterpolator(new GradientInterpolator);
   
   ParGridFunction x_adj(&fespace_adj);
   x_adj = 0.0;
   x_adj.ProjectBdrCoefficient(x_Diri_c_adj, ess_bdr_adj);

   ParGridFunction dx_adj(&fespace_nd_adj);
   dx_adj = 0.0;
   DiscreteLinearOperator grad_adj(&fespace_adj, &fespace_nd_adj);
   grad_adj.AddDomainInterpolator(new GradientInterpolator);
   
   // 9. Define Datacollections (Paraview and Gmsh) used for the postoperation.
   //==========================================================================
   ParaViewDataCollection *pd_magsta = NULL;
   pd_magsta = new ParaViewDataCollection("Capteur_MagSta", pmesh);
   pd_magsta->SetPrefixPath("ParaView_Capteur");
   pd_magsta->RegisterField("solution", &x);
   pd_magsta->RegisterField("adjoint", &x_adj);
   pd_magsta->RegisterField("der_solution", &dx);
   pd_magsta->RegisterField("der_adjoint", &dx_adj);
   pd_magsta->RegisterField("j0_gf", &j0_gf);
   pd_magsta->SetLevelsOfDetail(order);
   pd_magsta->SetDataFormat(VTKFormat::BINARY);
   pd_magsta->SetHighOrderOutput(true);

   ParaViewDataCollection *pd_magsta_gfc = NULL;
   pd_magsta_gfc = new ParaViewDataCollection("Capteur_MagSta", pmesh);
   pd_magsta_gfc->SetPrefixPath("ParaView_Capteur");
   pd_magsta_gfc->RegisterField("solution", &x);
   pd_magsta_gfc->SetLevelsOfDetail(order);
   pd_magsta_gfc->SetDataFormat(VTKFormat::BINARY);
   pd_magsta_gfc->SetHighOrderOutput(true);

   GmshDataCollection *gd_magsta = NULL;
   gd_magsta = new GmshDataCollection("", pmesh);
   gd_magsta->SetPrefixPath(gd_path);
   gd_magsta->RegisterField("az_dir", &x);
   gd_magsta->RegisterField("az_dir_0", &x);
   gd_magsta->RegisterField("az_adj", &x_adj);

   GmshDataCollection *gd_magsta_gfc = NULL;
   gd_magsta_gfc = new GmshDataCollection("", pmesh);
   gd_magsta_gfc->SetPrefixPath(gd_path);

   std::cout << "... " << std::endl;
   std::cout << "Start initializing the optimization problem... " << std::endl;
   std::cout << "... " << std::endl;
   
   double epsilon = 1e-10, tolerance = 1e-6;
   double step = 1.0;
   double min_step = 1.0e-3, max_step = 1.0e3;
   double lambda = 1.0e-2;
   int pow = 3; // 7;
   int num_max_ti_optim = 100;
   int it_linear_search_max = 9;

   double objective_value, new_objective_value;

   std::string old_name = "az_dir_0";
   std::string new_name = "az_dir_0";

   ParGridFunction x_old(&fespace);
   x_old.ProjectGridFunction(x);

   // 1. Initialization
   //==================

   // 1.1. h <-- h_init
   //==================
   Vector h_init(pmesh->attributes.Max());
   h_init = 0.0;
   h_init(2004) = 0.5;
   h_init(2005) = 0.5;
   PWConstCoefficient h_init_c(h_init);

   Vector nu(pmesh->attributes.Max());
   double mu0 = 4*M_PI*1.0e-7, mu_r = 1.0e0;
   nu = 1.0/mu0;
   nu(2004) = 1.0/(mu_r*mu0);
   nu(2005) = 1.0/(mu_r*mu0);
   PWConstCoefficient nu_c(nu);


   ParGridFunction h_gf_old(&fespace_vol); // Attention changer l'espace fonctionnel
   h_gf_old.ProjectCoefficient(h_init_c);
   GridFunctionCoefficient h_gf_coeff_old(&h_gf_old);

   ParGridFunction h_gf_new(&fespace_vol); // Attention changer l'espace fonctionnel
   h_gf_new.ProjectCoefficient(h_init_c);
   GridFunctionCoefficient h_gf_coeff_new(&h_gf_new);

   ParGridFunction sensitivity_gf(&fespace_vol); // Attention changer l'espace fonctionnel 
   sensitivity_gf.ProjectCoefficient(h_init_c);
   GridFunctionCoefficient sensitivity(&sensitivity_gf);

   // 1.2. Solve direct problem 
   //==========================
   int Flag_Update_Matrix_Direct = 1;
   int Flag_Update_RHS_Direct = 1;
   std::cout << "... " << std::endl;
   std::cout << "Start solving the direct problem... " << std::endl;
   std::cout << "... " << std::endl;
   MagSta_Problem *direct_Problem = new MagSta_Problem(fespace, &h_init_c, &js_c, &h_init_c);
   // Coefficient *nu_hat = direct_Problem->UpdateCoefficient(h_init_c, pow);
   direct_Problem->UpdateCoefficient(h_init_c, pow);
   direct_Problem->ComputeBilinearForm(Flag_Update_Matrix_Direct);
   direct_Problem->ComputeLinearForm(Flag_Update_RHS_Direct);
   direct_Problem->Solve(x, ess_tdof_list);
   std::cout << "... " << std::endl;
   std::cout << "End solving the direct problem... " << std::endl;
   std::cout << "... " << std::endl;

   // 1.3. Compute objective function
   //================================
   new_objective_value = direct_Problem->Evaluate_Objective_Function(x, lambda, h_gf_coeff_new); // ATTENTION : compute volume à corriger
   objective_value = new_objective_value;
   std::cout << "The initial value of the objective function objective_value = " << objective_value << std::endl;

   std::cout << "Before. Done Printing the fields in Gmsh = " << std::endl;
   int _flag_data_type = 1;
   int _ncomp = 1;

   // Print the solution
   //===================
   old_name = new_name;
   
   std::string string_num_it_optim;
   string_num_it_optim = to_string(000);
   new_name = "az_dir_" + std::to_string(0) + "_" + string_num_it_optim;
   direct_Problem->PrintGridFunction_UsingGmsh(gd_magsta, _flag_data_type, _ncomp);
   direct_Problem->PrintGridFunction_UsingGmsh(gd_magsta, x, old_name, new_name);

   // Print the initial design
   //=========================
   std::string name_field = "h_design";
   direct_Problem->PrintGridFunctionCoefficient_UsingGmsh(gd_magsta_gfc, h_gf_coeff_new, name_field,_flag_data_type, _ncomp);
   std::cout << "After. Done Printing the fields in Gmsh = " << std::endl;

   // 2. Steepest descent optimization. Optimization loop
   //====================================================
   for (int it_optim = 0; it_optim <= num_max_ti_optim; it_optim++)
   {
      std::cout << "Optimization iteration number " << it_optim << "\n" << std::endl;

      // u_old = u_new
      //==============
      x_old.ProjectGridFunction(x);

      // Project h_gf_coeff_new on both h_gf_new and h_gf_old
      //===================================================== 
      h_gf_new.ProjectCoefficient(h_gf_coeff_new);// Mise à jour de h_gf_coef_new !!!!
      h_gf_old.ProjectCoefficient(h_gf_coeff_new); 

      // Compute onbjective function
      //============================
      new_objective_value = direct_Problem->Evaluate_Objective_Function(x, lambda, h_gf_coeff_new);

      // Compute adjoint solution
      //========================
      std::cout << "Start solving the adjoint problem... " << "\n" << std::endl;
      int Flag_Update_Matrix_Adjoint = 1; 
      int Flag_Update_RHS_Adjoint = 1; 

      MagSta_Problem *adjoint_Problem = new MagSta_Problem(fespace_adj, &h_init_c, &j0_coeff);
      // Coefficient *nu_hat = adjoint_Problem->UpdateCoefficient(h_gf_coeff_new, pow);
      adjoint_Problem->UpdateCoefficient(h_gf_coeff_new, pow);
      adjoint_Problem->ComputeBilinearForm(Flag_Update_Matrix_Adjoint);
      adjoint_Problem->ComputeLinearForm(Flag_Update_RHS_Adjoint);
      adjoint_Problem->Solve(x_adj, ess_tdof_list_adj);
      double flux_value = adjoint_Problem->ComputeFlux(x);

      std::cout << "End solving the adjoint problem... " << "\n" << std::endl;
      std::cout << "The magnetic flux Phi = " << flux_value << " Wb. \n" << std::endl;

      grad.Assemble();
      grad.Finalize();
      grad.Mult(x, dx);

      grad_adj.Assemble();
      grad_adj.Finalize();
      grad_adj.Mult(x_adj, dx_adj);

      int ti = 1;
      double t = 0;
      cout << " ParaViewDataCollection step = " << ti << ", t = " << t << endl;
      pd_magsta->SetCycle(ti);
      pd_magsta->SetTime(t);
      pd_magsta->Save();

      // Compute and save sensitivity
      //====================
      // direct_Problem->ComputeSensitivity(dx, dx_adj, *nu_hat, sensitivity_gf);
      direct_Problem->ComputeSensitivity(dx, dx_adj, *(adjoint_Problem->get_LHS_Coeff()), sensitivity_gf);

      int _flag_data_type = 1;
      int _ncomp = 1;
      direct_Problem->PrintGridFunction_UsingGmsh(gd_magsta, _flag_data_type, _ncomp);
      std::string name_field = "sensitivity_";

      // Prepare the char used to write the coefficients
      //================================================
      std::string string_num_it_optim;
      if((it_optim >= 0) && (it_optim < 10))
      {
         string_num_it_optim = "00" + to_string(it_optim);
      }
      else if(((it_optim >= 10) && (it_optim < 100)))
      {
         string_num_it_optim = "0" + to_string(it_optim);
      }
      else
      {
         string_num_it_optim = "" + to_string(it_optim);            
      }

      direct_Problem->PrintGridFunctionCoefficient_UsingGmsh(gd_magsta_gfc, sensitivity, name_field + string_num_it_optim,_flag_data_type, _ncomp);

      // Line search 
      //============
      int flag_line_search = 0;
      if(flag_line_search)
      {      

      }
      else
      {
         // for (int it_linear_search = 0; it_linear_search <= it_linear_search_max; it_linear_search++) 
         for (int it_linear_search = 0; it_linear_search <= 1; it_linear_search++) 
         {

            // Compute h_gf_new_coeff
            //=======================

            // Compute tau
            //============
            double tau;
            if(1)
            {
               // tau = 1.0e4;
               tau = 1.0e02;
            }
            else
            {
               // Compute tau and h = h_old_gd_coeff - tau * sensitivity
               //=======================================================
               ParGridFunction this_GF = ParGridFunction(&fespace);
               this_GF.ProjectCoefficient(sensitivity);
               Vector vec_sensitivity;
               this_GF.SetFromTrueDofs(vec_sensitivity);
               // tau = 1.0e-4 * step/(epsilon + vec_sensitivity.Normlinf() );
               tau = 1.0e-3 * step/(epsilon + vec_sensitivity.Normlinf() );
            }

            // Compute h_gf_new_coeff
            //=======================
            SumCoefficient *h_gf_coeff_new_tmp_tmp = new SumCoefficient(h_gf_coeff_old, sensitivity, 1, -tau);
            h_gf_new.ProjectCoefficient(*h_gf_coeff_new_tmp_tmp);
            Vector h_vector;
            h_gf_new.GetTrueDofs(h_vector);
            for (int i  = 0; i < h_vector.Size(); i++)
            {
               if(h_vector.Elem(i) > 1)
               {
                  h_vector.Elem(i) = 1;
               }
               else if(h_vector.Elem(i) < 0)
               {
                  h_vector.Elem(i) = 0;
               }
            }
            h_gf_new.SetFromTrueDofs(h_vector);
            GridFunctionCoefficient *h_gf_coeff_new_tmp = new GridFunctionCoefficient(&h_gf_new);

            // Prepare the character used for writing coefficients
            //====================================================
            std::string string_num_it_optim = "";
            if((it_optim >= 0) && (it_optim < 10)) {
               string_num_it_optim = "00" + to_string(it_optim) + to_string(it_linear_search);
            }
            else if(((it_optim >= 10) && (it_optim < 100))) {
               string_num_it_optim = "0" + to_string(it_optim) + to_string(it_linear_search);
            }
            else {
               string_num_it_optim = "" + to_string(it_optim) + to_string(it_linear_search);            
            }

            direct_Problem->PrintGridFunctionCoefficient_UsingGmsh(gd_magsta_gfc, *h_gf_coeff_new_tmp, "h_design_" + string_num_it_optim,_flag_data_type, _ncomp);

            Vector tmp_x_new_vec;
            x.GetTrueDofs(tmp_x_new_vec);
            x_old.SetFromTrueDofs(tmp_x_new_vec);

            // Solve the direct problem and save the solution
            //===============================================
            std::cout << "Start solving the direct problem ..." << std::endl;
            direct_Problem->UpdateCoefficient(*h_gf_coeff_new_tmp, pow);
            direct_Problem->ComputeBilinearForm(Flag_Update_Matrix_Direct);
            direct_Problem->ComputeLinearForm(Flag_Update_RHS_Direct);
            direct_Problem->Solve(x, ess_tdof_list);
            std::cout << "End solving the direct problem ..." << std::endl;

            old_name = new_name;
            new_name = "az_dir_" + string_num_it_optim;
            direct_Problem->PrintGridFunction_UsingGmsh(gd_magsta, x, old_name, new_name);

            // compute the value of the objective function
            //============================================
            new_objective_value = direct_Problem->Evaluate_Objective_Function(x, lambda, *h_gf_coeff_new_tmp);

            // // Adjust new step
            // //================
            // if( (new_objective_value < (objective_value + tolerance * new_objective_value) ) || (step < min_step) ) {
            //    step = std::min(1.1 * step, max_step);
            //    std::cout << "[objective_new < (objective_old + tol * objective_new) ] || [step < min_step] " << std::endl;
            //    break;
            // }
            // else {
            //    step = std::max(0.6 * step, min_step);
            //    Vector tmp_x_old_vec;
            //    x_old.GetTrueDofs(tmp_x_old_vec);
            //    x.SetFromTrueDofs(tmp_x_old_vec);
            //    h_gf_coeff_new = h_gf_coeff_old;
            //    new_objective_value = objective_value;
            //    std::cout << "NOT ([objective_new < (objective_old + tol * objective_new) ] || [step < min_step] ) " << std::endl;
            // }

            // sleep(1.0);

            // Delete all objects created
            //===========================
            delete h_gf_coeff_new_tmp;
            std::cout << "End computing for the optimization iteration number " << it_optim << "\n" << std::endl;
         }
      }
   } 
   
   std::ostringstream solution_name_x;
   std::ostringstream solution_name_x_adj;
   int ti = 1;
   double t = 0;
   {
       cout << " ParaViewDataCollection step = " << ti << ", t = " << t << endl;
      {
         pd_magsta->SetCycle(ti);
         pd_magsta->SetTime(t);
         pd_magsta->Save();
      }
   }

   int time_step_gmsh = 0;
   double time_gmsh = 0.0;
   {
      cout << " GmshDataCollection step = " << ti << ", t = " << t << endl;
      gd_magsta->SetCycle (time_step_gmsh);
      gd_magsta->SetTime(time_gmsh);
      int flag_data_type = 1;
      int ncomp = 1;
      gd_magsta->AddSolution("az", solution_name_x, flag_data_type, ncomp);
      std::cout << "Saving the solution " << "az...\n" << std::endl;
      gd_magsta->Save("az", solution_name_x);
      std::cout << "End saving the solution " << "az...\n" << std::endl;

      gd_magsta->AddSolution("az_adj", solution_name_x_adj, flag_data_type, ncomp);
      std::cout << "Saving the solution " << "az_adj...\n" << std::endl;
      gd_magsta->Save("az_adj", solution_name_x_adj);
      std::cout << "End saving the solution " << "az_adj...\n" << std::endl;
   }
   return 0;
}