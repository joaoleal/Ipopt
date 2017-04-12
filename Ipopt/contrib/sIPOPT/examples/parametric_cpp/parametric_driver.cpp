// Copyright 2010, 2011 Hans Pirnay
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// Date   : 2010-01-05

#include "parametricTNLP.hpp"

#include "IpIpoptApplication.hpp"
#include "SensApplication.hpp"
#include "IpPDSearchDirCalc.hpp"
#include "IpIpoptAlg.hpp"
#include "SensRegOp.hpp"

int main(int argv, char**argc)
{
  using namespace Ipopt;

  SmartPtr<IpoptApplication> app_ipopt = new IpoptApplication();

  // Register sIPOPT options
  RegisterOptions_sIPOPT(app_ipopt->RegOptions());
  app_ipopt->Options()->SetRegisteredOptions(app_ipopt->RegOptions());
  //app_ipopt->Options()->SetIntegerValue("print_level", 11);

  // Call Initialize the first time to create a journalist, but ignore
  // any options file
  ApplicationReturnStatus retval;
  retval = app_ipopt->Initialize("");
  if (retval != Solve_Succeeded) {
    //printf("ampl_ipopt.cpp: Error in first Initialize!!!!\n");
    exit(-100);
  }
  app_ipopt->Initialize();

  // create AmplSensTNLP from argc. This is an nlp because we are using our own TNLP Adapter
  SmartPtr<TNLP> sens_tnlp = new ParametricTNLP();

  retval = app_ipopt->OptimizeTNLP(sens_tnlp);

  /**
   * give pointers to Ipopt algorithm objects to Sens Application
   */
  SmartPtr<SensApplication> app_sens = new SensApplication(app_ipopt->Jnlst(),
                                                           app_ipopt->Options(),
                                                           app_ipopt->RegOptions());
  app_sens->Initialize();

  app_sens->SetIpoptAlgorithmObjects(app_ipopt, retval);

  /* In this function, the indices for the parametric computations are set.
   * To keep track of the parameters, each parameter gets an index from 1 to n_parameters.
   * In this case, [1] eta_1, [2] eta_2.
   * The following metadata vectors are important:
   */

  Index n = 5;
  /* In this list of Numbers (=doubles), the perturbed
   *    values for the parameters are set.
   */
  std::vector<Number> sens_state_value(n,0);
  sens_state_value[3] = 4.5;
  sens_state_value[4] = 1.0;

  app_sens->SetUpdatedParameters(sens_state_value);

  printf("\n");
  printf("#-------------------------------------------\n");
  printf("# Sensitivity without bound checking\n");
  printf("#-------------------------------------------\n");

  app_sens->Run();
  
  printf("\n");
  printf("#-------------------------------------------\n");
  printf("# Sensitivity with bound checking\n");
  printf("#-------------------------------------------\n");
  app_ipopt->Options()->SetStringValue("sens_boundcheck", "yes");
  app_sens->Run();
}
