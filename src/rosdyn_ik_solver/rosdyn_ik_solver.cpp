/*
Copyright (c) 2022, JRL-CARI CNR-STIIMA/UNIBS
Manuel Beschi manuel.beschi@unibs.it
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
    * Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.
    * Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in the
      documentation and/or other materials provided with the distribution.
    * Neither the name of the <organization> nor the
      names of its contributors may be used to endorse or promote products
      derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL <COPYRIGHT HOLDER> BE LIABLE FOR ANY
DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#include <rosdyn_ik_solver/rosdyn_ik_solver.h>
#include <pluginlib/class_list_macros.h>

PLUGINLIB_EXPORT_CLASS(ik_solver::RosdynIkSolver, ik_solver::IkSolver)

namespace ik_solver
{

inline bool RosdynIkSolver::customConfig()
{
  Eigen::Vector3d gravity;
  gravity << 0,0,-9.806;
  chain_ = rosdyn::createChain(model_,base_frame_,flange_frame_,gravity);
  return chain_->setInputJointsName(joint_names_);

}


std::vector<Eigen::VectorXd> RosdynIkSolver::getIk(const Eigen::Affine3d& T_base_flange,
                                                   const std::vector<Eigen::VectorXd> & seeds,
                                                   const int& desired_solutions,
                                                   const int& max_stall_iterations)
{
  std::vector<Eigen::VectorXd > solutions;

  unsigned int n_seed = seeds.size();
  bool found = false;



  Eigen::VectorXd center=(ub_+lb_)*0.5;
  Eigen::VectorXd width=(ub_-lb_)*0.5;

  int stall=0;
  int iter=0;
  for (; iter<max_iter_; iter++)
  {
    if (stall>max_stall_iterations)
    {
      if (solutions.size()==0)
        ROS_DEBUG("reach stall generations without any solution");
      else
        ROS_DEBUG("reach stall generation (iteration %u)",iter);
      break;
    }
    if (solutions.size()>desired_solutions)
    {
      ROS_DEBUG("reach desired solutions number");
      break;
    }
    Eigen::VectorXd js;
    Eigen::VectorXd start;
    Eigen::VectorXd r(lb_.size());
    r.setRandom(); // between [-1,1]

    if (iter<n_seed)
    {
      start=seeds.at(iter);
    }
    else
    {
      start=center+width.cwiseProduct(r);
    }

    if (outOfBound(start))
      continue;

    stall++;

    if (chain_->computeLocalIk(js,T_base_flange,start,1e-6,ros::Duration(0.005)))
    {

      if (outOfBound(js))
        continue;
      bool is_diff = true;

      for (const Eigen::VectorXd& sol: solutions)
      {
        if ((sol-js).norm()<TOLERANCE)
        {
          is_diff = false;
          break;
        }
      }
      if (not is_diff)
        continue;

      stall=0;
      std::vector<Eigen::VectorXd> multiturn = chain_->getMultiplicity(js);
      for (const Eigen::VectorXd& tmp: multiturn)
      {
        if (outOfBound(tmp))
          continue;

        is_diff = true;
        for (const Eigen::VectorXd& sol: solutions)
        {
          if ((sol-tmp).norm()<TOLERANCE)
          {
            is_diff = false;
            break;
          }
        }
        if (not is_diff)
          continue;

        solutions.push_back(tmp);
      }
      found = true;


    }
  }

  return solutions;
}


}   // end namespace ik_solver
