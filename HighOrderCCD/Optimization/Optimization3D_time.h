#ifndef OPTIMIZATION3D_TIME_H
#define OPTIMIZATION3D_TIME_H

#include "HighOrderCCD/Config/Config.h"

#include "HighOrderCCD/Distance.h"
#include "HighOrderCCD/BVH/BVH.h"
#include "HighOrderCCD/CCD/CCD.h"

#include "HighOrderCCD/Subdivide.h"
#include "HighOrderCCD/Energy.h"
#include "HighOrderCCD/Gradient.h"
#include "HighOrderCCD/Step.h"

#include <vector>
#include <ctime>
#include <Eigen/SparseCholesky>
#include <unsupported/Eigen/AutoDiff>
#include <unsupported/Eigen/KroneckerProduct>
//#include "gurobi_c++.h"
PRJ_BEGIN

class Optimization3D_time 
{
public:

  typedef Eigen::MatrixXd Data;
 //Eigen::Dynamic
  typedef Eigen::Matrix<double,Eigen::Dynamic,1> inner_derivative_t;//3*(order_num+1)
  typedef Eigen::AutoDiffScalar<inner_derivative_t> inner_scalar_t;
  typedef Eigen::Matrix<inner_scalar_t,Eigen::Dynamic,1> derivative_t;
  typedef Eigen::AutoDiffScalar<derivative_t> scalar_t;
  typedef Eigen::Matrix<scalar_t,Eigen::Dynamic,1> Vec12;
  typedef Eigen::Matrix<scalar_t,1,3> Vec3;

  typedef Eigen::SparseMatrix<double> SpMat; // declares a column-major sparse matrix type of double

  typedef std::vector< std::tuple< int, std::pair<double,double>, Eigen::MatrixXd > > Tree;
  typedef std::tuple< int, std::pair<double,double>, Eigen::MatrixXd >  Node;


  static void optimization(Data& spline, double& piece_time,
                           const Eigen::MatrixXd & V,const Eigen::MatrixXi& F,
                           BVH& bvh)
  {
      Data direction;
      double t_direction;
      
      if(max_step<1e-10)
      {
        Subdivide::update_tree(spline,V,F,bvh);
      }
      if(adaptive_change)
      {
        Subdivide::update_tree(spline,V,F,bvh);
        Subdivide::update_vel(spline,piece_time,bvh);
        Subdivide::update_acc(spline,piece_time,bvh);
      }
      
    
      bvh.InitTrajectory(spline);
      
      std::cout<<"tree_size: "<<subdivide_tree.size()<<std::endl;
      if(adaptive_change)
      {
        std::cout<<"vel_size: "<<vel_tree.size()<<std::endl;
        std::cout<<"acc_size: "<<acc_tree.size()<<std::endl;
      }
      /*
      for(unsigned int i=0;i<subdivide_tree.size();i++)
      {
        std::cout<<std::get<1>(subdivide_tree[i]).first<<" "<<std::get<1>(subdivide_tree[i]).second<<std::endl;
      }
      */
     
      //if(iter>0 && iter%100==0 && lambda <= 1/margin2)
      //lambda*=2;
      Energy::barrier_energy(spline,V,F,bvh);
      
      clock_t time0 = clock();
      descent_direction(spline, direction, piece_time,t_direction, V,F,bvh);
      clock_t time1 = clock();     

      
      line_search(spline, direction,piece_time, t_direction, V,F, bvh, false);//vtxs
      
      //update_time(spline, bvh);

      std::cout<<"piece_time:"<<piece_time<<std::endl;

      
      clock_t time2 = clock();
      
      
      std::cout<<std::endl<<"time10:"<<(time1-time0)/(CLOCKS_PER_SEC/1000)<<std::endl<<std::endl;
      std::cout<<"time21:"<<(time2-time1)/(CLOCKS_PER_SEC/1000)<<std::endl<<std::endl;

  }


  static void descent_direction(const Data& spline, Data& direction, const double& piece_time, double& t_direction,
                                const Eigen::MatrixXd & V,const Eigen::MatrixXi& F,
                                BVH& bvh)
  {
    int t_n=trajectory_num;

    Eigen::VectorXd grad,g1,g2,g3, g4, partgrad;
    Eigen::MatrixXd hessian,h1,h2,h3;
    double g_t,h_t,g_t1,g_t2,h_t1,h_t2;
    //clock_t time0 = clock();
    Gradient::dynamic_gradient(spline,piece_time,g1,h1,g_t1,h_t1);

    clock_t time1 = clock();
    Gradient::fast_barrier_gradient(spline, g2,h2, V,F,bvh);//fast_
    clock_t time2 = clock();
    
    Gradient::bound_gradient(spline,piece_time, g3,h3, g_t2, h_t2, partgrad);
    //Gradient::auto_bound_gradient(spline,piece_time, g3,h3);
    clock_t time3 = clock();

    std::cout<<"time_g2:"<<(time2-time1)/(CLOCKS_PER_SEC/1000)<<std::endl<<std::endl;
    std::cout<<"time_g3:"<<(time3-time2)/(CLOCKS_PER_SEC/1000)<<std::endl<<std::endl;

    //double g_t,h_t;
    //double dy_energy=Energy::dynamic_energy(spline,piece_time);
    //Gradient::time_gradient(spline,piece_time, dy_energy, g_t, h_t, partgrad);

    //std::cout<<g1<<"\n\n";
    //std::cout<<g_t<<"\n\n";
    //std::cout<<h_t<<"\n\n";
    grad=ks*g1 + lambda*g2 + lambda*g3;
    hessian =ks*h1 + lambda*h2 + lambda*h3;

    g_t=g_t1+lambda*g_t2;
    h_t=h_t1+lambda*h_t2;

    int fixed_size=fixed_points.size();
    for(int i=0;i<fixed_size;i++)
    {
      int id=fixed_points[i]*(order_num-2)+1;
      for(int j=0;j<3;j++)
      {
        grad.row(3*(id-1)+j)-=0.5*grad.row(3*id+j);
        grad.row(3*(id+1)+j)-=0.5*grad.row(3*id+j);

        g1.row(3*(id-1)+j)-=0.5*g1.row(3*id+j);
        g1.row(3*(id+1)+j)-=0.5*g1.row(3*id+j);

        partgrad.row(3*(id-1)+j)-=0.5*partgrad.row(3*id+j);
        partgrad.row(3*(id+1)+j)-=0.5*partgrad.row(3*id+j);

        hessian.row(3*(id-1)+j)-=0.5*hessian.row(3*id+j);
        hessian.row(3*(id+1)+j)-=0.5*hessian.row(3*id+j);

        hessian.col(3*(id-1)+j)-=0.5*hessian.col(3*id+j);
        hessian.col(3*(id+1)+j)-=0.5*hessian.col(3*id+j);
      }
    }
    std::vector<int> id_list;
    id_list.push_back(1);
    for(int i=0;i<fixed_size;i++)
    {
      int id=fixed_points[i]*(order_num-2)+1;
      id_list.push_back(id);
    }
    id_list.push_back(t_n-2);

    Eigen::VectorXd ng((t_n-4-fixed_size)*3), g1_((t_n-4-fixed_size)*3), partgrad_((t_n-4-fixed_size)*3);
    Eigen::MatrixXd h((t_n-4-fixed_size)*3,(t_n-4-fixed_size)*3);

    int init0=0;
    for(int i=0;i<(int)id_list.size()-1;i++)
    {
      int id0=id_list[i];
      int id1=id_list[i+1];
      ng.segment(init0*3,(id1-id0-1)*3)=-grad.segment((id0+1)*3,(id1-id0-1)*3);
      g1_.segment(init0*3,(id1-id0-1)*3)=g1.segment((id0+1)*3,(id1-id0-1)*3);
      partgrad_.segment(init0*3,(id1-id0-1)*3)=partgrad.segment((id0+1)*3,(id1-id0-1)*3);
      
      int init1=0;
      for(int j=0;j<(int)id_list.size()-1;j++)
      {
        int id_0=id_list[j];
        int id_1=id_list[j+1];

        h.block(init0*3,init1*3,(id1-id0-1)*3,(id_1-id_0-1)*3) = hessian.block((id0+1)*3,(id_0+1)*3,(id1-id0-1)*3,(id_1-id_0-1)*3);
        init1+=id_1-id_0-1;
      }
      init0+=id1-id0-1;
    }
    
    Eigen::VectorXd ng0((t_n-4-fixed_size)*3+1);
    Eigen::MatrixXd h0((t_n-4-fixed_size)*3+1,(t_n-4-fixed_size)*3+1);
    ng0.head((t_n-4-fixed_size)*3)=ng;
    ng0((t_n-4-fixed_size)*3)=-g_t;

    h0.block(0,0,(t_n-4-fixed_size)*3,(t_n-4-fixed_size)*3)=h;
    h0((t_n-4-fixed_size)*3,(t_n-4-fixed_size)*3)=h_t;//h_t;

    g4=-ks*(2*der_num-1)*g1_/piece_time+lambda*partgrad_;

    h0.block(0,(t_n-4-fixed_size)*3,(t_n-4-fixed_size)*3,1)=g4;
    h0.block((t_n-4-fixed_size)*3,0,1,(t_n-4-fixed_size)*3)=g4.transpose();
  
    /*
    Eigen::VectorXd ng = -grad.segment(6,(t_n-4)*3);
    Eigen::MatrixXd h = hessian.block(6,6,(t_n-4)*3,(t_n-4)*3);
    
    Eigen::VectorXd ng0((t_n-4)*3+1);
    Eigen::MatrixXd h0((t_n-4)*3+1,(t_n-4)*3+1);
    ng0.head((t_n-4)*3)=ng;
    ng0((t_n-4)*3)=-g_t;

    h0.block(0,0,(t_n-4)*3,(t_n-4)*3)=h;
    h0((t_n-4)*3,(t_n-4)*3)=h_t;//h_t;

    g4=-ks*(2*der_num-1)*g1.segment(6,(t_n-4)*3)/std::pow(piece_time,2*der_num)+lambda*partgrad.segment(6,(t_n-4)*3);

    h0.block(0,(t_n-4)*3,(t_n-4)*3,1)=g4;
    h0.block((t_n-4)*3,0,1,(t_n-4)*3)=g4.transpose();
    */
    Eigen::VectorXd x, x0;
    
    Eigen::MatrixXd I=h0; I.setIdentity();
    Eigen::LLT<Eigen::MatrixXd> solver; 
    
    solver.compute(h0);
    
    if(solver.info() == Eigen::NumericalIssue)//50
    {
      Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigensolver(h0);
      Eigen::MatrixXd eigenvalue=eigensolver.eigenvalues();
      if(eigenvalue(0)<0)
      {
        //std::cout<<"eigenvalue:"<<eigenvalue(0)<<std::endl;
        h0=h0-eigenvalue(0)*I+0.01*I;
      }
      solver.compute(h0);    

    }
    
    
    x0 = solver.solve(ng0);
    //x0=ng0;
    wolfe=x0.dot(ng0);

    x=x0.head((t_n-4-fixed_size)*3);
    
    t_direction=x0((t_n-4-fixed_size)*3);
    /*
    SpMat H=h.sparseView();
    Eigen::SimplicialLLT<SpMat> solver;  // performs a Cholesky factorization of A
    solver.compute(H);

    SpMat I=H; I.setIdentity();
    while(solver.info()!=Eigen::Success)
    {
      H=H+I;
      solver.compute(H);
    }
    x = solver.solve(ng);

    wolfe=x.dot(ng);
    */
    
    //Eigen::MatrixXd d(Eigen::Map<Eigen::MatrixXd>(x.data(), 3,t_n-2));
    //Eigen::MatrixXd ngrad(Eigen::Map<Eigen::MatrixXd>(ng.data(), 3,t_n-2));
    //Eigen::MatrixXd d(Eigen::Map<Eigen::MatrixXd>(x.data(), 3,t_n-4));

    Eigen::MatrixXd d(Eigen::Map<Eigen::MatrixXd>(x.data(), 3,t_n-4-fixed_size));
    //Eigen::MatrixXd ngrad(Eigen::Map<Eigen::MatrixXd>(ng.data(), 3,t_n-4));
    Eigen::MatrixXd d_=d.transpose();
    
    direction.resize(t_n,3);
    
    //direction.block(1,0,t_n-2,3)=d.transpose();
    //direction.block(2,0,t_n-4,3)=d.transpose();

    init0=0;
    for(int i=0;i<(int)id_list.size()-1;i++)
    {
      int id0=id_list[i];
      int id1=id_list[i+1];
      direction.block(id0+1,0,id1-id0-1,3)=d_.block(init0,0,id1-id0-1,3);
      
      init0+=id1-id0-1;
    }

    for(int i=0;i<fixed_size;i++)
    {
      int id=fixed_points[i]*(order_num-2)+1;
      direction.row(id)=-0.5*(direction.row(id-1)+direction.row(id+1));
    }

    
    direction.row(0).setZero();
    direction.row(t_n-1).setZero();
    direction.row(1).setZero();
    direction.row(t_n-2).setZero();

    //direction.col(2).setZero();
    
    std::cout<<"gn:"<<ng0.norm()<<std::endl;
    std::cout<<"dn:"<<x0.norm()<<std::endl;
    std::cout<<"t_direction:"<<t_direction<<std::endl;
    
    gnorm=ng0.norm();
    //std::cout<<"gnorm:"<<gnorm<<std::endl;
  }

  static void line_search(Data& spline, const Data& direction, double& piece_time, const double& t_direction,
                          const Eigen::MatrixXd & V,const Eigen::MatrixXi& F,
                          BVH& bvh, bool exact)
  {
    
    double step=Step::position_step(spline, direction,V,F, bvh);
    //double time_step=Step::time_step(spline, direction, bvh);
    //if(time_step<step)
      //step=time_step;
    
    std::cout<<"highcdd:"<<step<<std::endl<<std::endl;
    if(piece_time+step*t_direction<=0)
    {
      step=-0.95*piece_time/t_direction;
    }

    std::cout.precision(10);
   
    std::cout<<"wolfe:"<<wolfe<<std::endl;


    // backtracking
    //time0 = clock();
    double e=Energy::fast_whole_energy(spline,piece_time,V,F,bvh);//fast_
    double init_time=piece_time;
   
    piece_time=init_time+step*t_direction;
    while(e-1e-4*wolfe*step<Energy::fast_whole_energy(spline+step*direction,piece_time,V,F,bvh))//fast_
    {
      step*=0.8;
      piece_time=init_time+step*t_direction;
    }
    //time1 = clock();
    //std::cout<<std::endl<<"searchtime:"<<(time1-time0)/(CLOCKS_PER_SEC/1000)<<std::endl;

    max_step=step;
    
    std::cout<<"step:"<<step<<std::endl;
    std::cout<<"result:"<<Energy::dynamic_energy(spline+step*direction,piece_time)<<
               " "<<lambda*Energy::barrier_energy(spline+step*direction,V,F,bvh)<<
               " "<<lambda*Energy::bound_energy(spline+step*direction,piece_time)<<
               " "<<kt*whole_weight*piece_time<<std::endl<<std::endl;
    //<<" "<<limit_energy(spline+step*direction)
    spline=spline+step*direction;

  }
};


PRJ_END

#endif
