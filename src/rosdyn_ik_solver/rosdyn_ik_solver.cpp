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
#if ROS_VERSION == 1
  #include <pluginlib/class_list_macros.h>
  #define DEBUG(s) ROS_DEBUG_STREAM(s)
  #define WARN(s) ROS_WARN_STREAM(s)
#elif ROS_VERSION == 2
  #include <pluginlib/class_list_macros.hpp>
  #define DEBUG(s) RCLCPP_DEBUG_STREAM(rclcpp::get_logger("RosdynIkSolver"), s)
  #define WARN(s) RCLCPP_WARN_STREAM(rclcpp::get_logger("RosdynIkSolver"), s)
#endif

PLUGINLIB_EXPORT_CLASS(ik_solver::RosdynIkSolver, ik_solver::IkSolver)

namespace ik_solver
{
inline bool RosdynIkSolver::config(const std::string& params_ns)
{
  if(!IkSolver::config(params_ns))
  {
    WARN("Ik solver initial config FAILED");
    return false;
  }
  Eigen::Vector3d gravity;
  gravity << 0,0,-9.806;
  chain_ = rdyn::createChain(*model_,base_frame_,flange_frame_,gravity);
  std::string what;
  return chain_->setInputJointsName(joint_names_, what);
}


Solutions RosdynIkSolver::getIk(const Eigen::Affine3d& T_base_flange,
                                const Configurations & seeds,
                                const int& desired_solutions,
                                const int& min_stall_iterations,
                                const int& max_stall_iterations)
{
  Solutions solutions;
  solutions.clear();

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
      if (solutions.configurations().size()==0)
        DEBUG("reach stall generations without any solution");
      else
        DEBUG("reach stall generation (iteration " << iter << ")");
      break;
    }
    if (solutions.configurations().size()>desired_solutions)
    {
      DEBUG("reach desired solutions number");
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
    std::vector<int> out_of_bound = outOfBound(start, ub_, lb_);
    if (std::find(out_of_bound.begin(), out_of_bound.end(), 0) != out_of_bound.end())
      continue;
    stall++;

    if (chain_->computeLocalIk(js,T_base_flange,start,1e-6, 0.005))
    {
      out_of_bound = outOfBound(js, ub_, lb_);
      if (std::find(out_of_bound.begin(), out_of_bound.end(), 0) != out_of_bound.end())
        continue;
      bool is_diff = true;

      for (const Eigen::VectorXd& sol: solutions.configurations())
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
        out_of_bound = outOfBound(tmp, ub_, lb_);
        if (std::find(out_of_bound.begin(), out_of_bound.end(), 0) != out_of_bound.end())
          continue;

        is_diff = true;
        for (const Eigen::VectorXd& sol: solutions.configurations())
        {
          if ((sol-tmp).norm()<TOLERANCE)
          {
            is_diff = false;
            break;
          }
        }
        if (not is_diff)
          continue;

        solutions.configurations().push_back(tmp);
        Eigen::Affine3d Tfk = chain_->getTransformation(tmp);
        solutions.translation_residuals().push_back((Tfk.translation() - T_base_flange.translation()).norm());
        solutions.rotation_residuals().push_back(Eigen::AngleAxisd(T_base_flange.linear().inverse() * Tfk.linear()).angle());
      }
      found = true;

    }
  }

  solutions.message() = "Solution found: " + std::to_string(found? "True":"False") + ". Iterations: " + std::to_string(iter+1) + " (" + std::to_string(min_stall_iter_) + " - " + std::to_string(max_stall_iter_) + "). Seeds: " + std::to_string(n_seed);

  return solutions;
}

Eigen::Affine3d RosdynIkSolver::getFK(const Configuration& s)
{
  return chain_->getTransformation(s);
}


}   // end namespace ik_solver
