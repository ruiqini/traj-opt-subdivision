#ifndef GRADIENT_H
#define GRADIENT_H

#include "Config/Config.h"

#include "HighOrderCCD/BVH/BVH.h"
#include "HighOrderCCD/CCD/CCD.h"
#include "HighOrderCCD/Distance.h"
#include "HighOrderCCD/Distance_der.h"

#include <vector>

#include <Eigen/SparseCholesky>
#include <unsupported/Eigen/AutoDiff>
#include <unsupported/Eigen/KroneckerProduct>

PRJ_BEGIN

class Gradient
{
  public:
    typedef Eigen::MatrixXd Data;
    typedef std::vector< std::tuple< int, std::pair<double,double>, Eigen::MatrixXd > > Tree;
    typedef std::tuple< int, std::pair<double,double>, Eigen::MatrixXd >  Node;

    typedef Eigen::Matrix<double,Eigen::Dynamic,1> inner_derivative_t;//3*(order_num+1)
    typedef Eigen::AutoDiffScalar<inner_derivative_t> inner_scalar_t;
    typedef Eigen::Matrix<inner_scalar_t,Eigen::Dynamic,1> derivative_t;
    typedef Eigen::AutoDiffScalar<derivative_t> scalar_t;
    typedef Eigen::Matrix<scalar_t,Eigen::Dynamic,1> Vec12;
    typedef Eigen::Matrix<scalar_t,1,3> Vec3;

    static void dynamic_gradient(const Data& spline, const double& piece_time, Eigen::VectorXd& grad, Eigen::MatrixXd& hessian, double& g_t, double& h_t)
    {
    
        int n=3*trajectory_num;
        
        grad.resize(n);
        grad.setZero();

        hessian.resize(n,n);
        hessian.setZero();
        
        Eigen::Matrix3d I;
        I.setIdentity();
        Eigen::MatrixXd B,M; 

        double dynamic_energy=0;

        //Dirichlet energy
        for(int sp_id=0;sp_id<piece_num;sp_id++)
        {
            Eigen::MatrixXd bz;
            int init=sp_id*(order_num-2);
            bz=spline.block<order_num+1,3>(init,0);
            Eigen::MatrixXd x=bz;
         
            M=convert_list[sp_id].transpose()*M_dynamic*convert_list[sp_id];
            M=convert_list[sp_id].transpose()*M_dynamic*convert_list[sp_id];

            for(int j=0;j<3;j++)
            {
                Eigen::VectorXd x0=bz.col(j);
                dynamic_energy+=1/std::pow(time_weight[sp_id]*piece_time,2*der_num-1)*0.5*x0.transpose()*M*x0;
            }

            x=M*x;
            B = Eigen::kroneckerProduct(M,I);
            x.transposeInPlace();
            Eigen::VectorXd v1(Eigen::Map<Eigen::VectorXd>(x.data(), 3*(order_num+1)));
            grad.segment<3*(order_num+1)>(3*init)+=1/std::pow(time_weight[sp_id]*piece_time,2*der_num-1)*v1;
            
            hessian.block<3*(order_num+1),3*(order_num+1)>(3*init,3*init)+=1/std::pow(time_weight[sp_id]*piece_time,2*der_num-1)*B;
        }   

        g_t=kt*whole_weight;
        g_t+=-ks*(2*der_num-1)*dynamic_energy/piece_time;

        h_t=ks*(2*der_num-1)*(2*der_num)*dynamic_energy/(piece_time*piece_time);   
    }

    static void fast_barrier_gradient(const Data& spline, 
                                    Eigen::VectorXd& grad, Eigen::MatrixXd& hessian,
                                    const Eigen::MatrixXd & V,const Eigen::MatrixXi& F,
                                    BVH& bvh)
    {
        
        int num=3*trajectory_num;
        
        grad.resize(num);
        grad.setZero();

        hessian.resize(num,num);
        hessian.setZero();

        //Eigen::VectorXd auto_grad=grad;
        //Eigen::MatrixXd auto_hessian=hessian;

        bvh.InitTrajectory(spline);

        std::vector<std::pair<unsigned int, unsigned int>> collision_pair, temp_pair;

        clock_t time1 = clock();
        bvh.CheckCollision(collision_pair,offset+margin);
        clock_t time2 = clock();
        std::cout<<"time_bvh:"<<(time2-time1)/(CLOCKS_PER_SEC/1000)<<std::endl<<std::endl;

        std::cout<<"bvh_size:"<<collision_pair.size()<<std::endl;
        double dmin=1.0;

        std::vector<std::vector<std::pair< int, int>>> segment_lists;

        std::vector<std::vector<int>> segment_ob_lists;

        segment_lists.resize(subdivide_tree.size());

        segment_ob_lists.resize(subdivide_tree.size());


        for(unsigned int i=0;i<collision_pair.size();i++)
        {            
            int tr_id=collision_pair[i].first;
            int ob_id=collision_pair[i].second;

            int sp_id=std::get<0>(subdivide_tree[tr_id]);
            double weight=std::get<1>(subdivide_tree[tr_id]).second-std::get<1>(subdivide_tree[tr_id]).first;
            Eigen::MatrixXd basis=std::get<2>(subdivide_tree[tr_id]);

            //std::cout<<"sp_id:"<<sp_id<<std::endl;
            Eigen::MatrixXd bz;
            int init=sp_id*(order_num-2);
            bz=spline.block<order_num+1,3>(init,0);

            std::vector<Eigen::RowVector3d> P(order_num+1);
            for(int j=0;j<=order_num;j++)
            {   
                P[j].setZero();
                for(int j0=0;j0<=order_num;j0++)
                {
                    P[j]+=basis(j,j0)*bz.row(j0);
                }
            }

            std::vector<Eigen::RowVector3d> S(3); 
            for(int j=0;j<3;j++)
            {
                S[j]=V.row(F(ob_id,j));
            }
            
            Eigen::RowVector3d C,C0,C1;
            double d;
        
            //point face
            
            for(int j=0;j<=order_num;j++)
            {
    
                Eigen::RowVector3d d_p;
                Eigen::Matrix3d h_p;
                Distance_der::point_face(P[j],S[0],S[1],S[2], d, d_p, h_p);
                
                d-=offset;
                if(d<margin)
                {
                    if(d<dmin)
                        dmin=d;
                    //-weight*(d-margin)*(d-margin)*log(d/margin)
                    //double e1=-weight*(2*(d-margin)*log(d/margin)+(d-margin)*(d-margin)/d);
                    double e1=-weight*((d*d-margin*margin)/(d*d)*log(d/margin)+(d-margin)*(d-margin)/(d*d));

                    Eigen::Matrix3d I; I.setIdentity();
                    Eigen::MatrixXd A=Eigen::kroneckerProduct(basis.row(j),I);
                    //std::cout<<A<<"\n";
                    Eigen::MatrixXd d_x=d_p*A;
                    
                    grad.segment(3*init,3*(order_num+1)) += e1*d_x.transpose();
                    
                    //double e2=-weight*(2*log(d/margin)+4*(d-margin)/d-(d-margin)*(d-margin)/(d*d));
                    double e2=-weight*(1/d+margin*margin/(d*d*d)*(2*log(d/margin)-1)+2*margin*(d-margin)/(d*d*d));

                    
                    hessian.block<3*(order_num+1),3*(order_num+1)>(3*init,3*init)+=e2*d_x.transpose()*d_x+e1*A.transpose()*h_p*A;
                }
            }

            //segment segment
            //for(int j=0;j<=order_num;j++)
            {
                int j=order_num;
                for(int k=0;k<3;k++)
                {
                    int k1=(k+1)%3;
                    
                    Distance<double,3>::segment_segment(P[j],P[(j+1)%(order_num+1)],
                                                        S[k],S[k1],d,C0,C1);
                    //std::cout<<"ss:"<<d<<std::endl;
                    d-=offset;
                    if(d<margin)
                    {
                        if(d<dmin)
                            dmin=d;
                        //segment_lists[sp_id].push_back(std::make_pair(j,k));
                        //segment_ob_lists[sp_id].push_back(ob_id);
                        //segment_tr_lists[sp_id].push_back(tr_id);
                        segment_lists[tr_id].push_back(std::make_pair(j,k));
                        segment_ob_lists[tr_id].push_back(ob_id);
                        
                    }
                }
            }
        }
        
        for(unsigned int tr_id=0;tr_id<subdivide_tree.size();tr_id++)
        {

            if(segment_lists[tr_id].size()==0)
            {
                continue;
            }

            int sp_id=std::get<0>(subdivide_tree[tr_id]);
            double weight=std::get<1>(subdivide_tree[tr_id]).second-std::get<1>(subdivide_tree[tr_id]).first;
            Eigen::MatrixXd basis=std::get<2>(subdivide_tree[tr_id]);

            //std::cout<<"sp_id:"<<sp_id<<std::endl;
            int init=sp_id*(order_num-2);
            
            Eigen::MatrixXd x=spline.block<order_num+1,3>(init,0);
            Vec12 X;
            X.resize(3*(order_num+1));
            for(int r=0;r<=order_num;r++)
            {
                X(3*r).value() = x(r,0);
                X(3*r+1).value() = x(r,1);
                X(3*r+2).value() = x(r,2);
            }
            //repeat partial derivatives for the inner AutoDiffScalar
            for(int id=0;id<3*(order_num+1);id++)
            {
                X(id).derivatives().resize(3*(order_num+1));
                X(id).derivatives().setZero();
                X(id).derivatives()(id)= 1;
                X(id).value().derivatives() = inner_derivative_t::Unit(3*(order_num+1),id);
            }
            //set the hessian matrix to zero
            for(int idx=0; idx<3*(order_num+1); idx++) {
                for(int id=0;id<3*(order_num+1);id++)
                {
                    X(id).derivatives()(idx).derivatives()  = inner_derivative_t::Zero(3*(order_num+1));
                }       
            }

            std::vector<Vec3> P(order_num+1); 

            for(int j=0;j<=order_num;j++)
            {
                P[j].setZero();
                for(int j0=0;j0<=order_num;j0++)
                {
                    P[j](0)+=basis(j,j0)*X(3*j0);
                    P[j](1)+=basis(j,j0)*X(3*j0+1);
                    P[j](2)+=basis(j,j0)*X(3*j0+2);
                }
            }
            
            Vec3 C,C0,C1;
            scalar_t d;

            std::vector<std::pair< int, int>> segment_list=segment_lists[tr_id];        

            //segment segment
        
            for(unsigned int j=0;j<segment_list.size();j++)
            {
                int ob_id=segment_ob_lists[tr_id][j];
                std::vector<Eigen::RowVector3d> S(3); 
                for(int i=0;i<3;i++)
                {
                    S[i]=V.row(F(ob_id,i));
                }
                int j0=segment_list[j].first;
                int k=segment_list[j].second;
                Distance<scalar_t,3>::segment_segment(P[j0],P[(j0+1)%(order_num+1)],
                                                        S[k],S[(k+1)%3],d,C0,C1);
                d=d-offset;
                //scalar_t e=-weight*(d-margin)*(d-margin)*log(d/margin);
                
                scalar_t e=-weight*(d-margin)*(d-margin)/d*log(d/margin);///d

                //scalar_t e=weight*(1-d/margin*d/margin)*(1-d/margin*d/margin);
                grad.segment(3*init,3*(order_num+1)) += e.value().derivatives();
                Eigen::Matrix<double,3*(order_num+1),3*(order_num+1)> B;

                for(int r=0;r<3*(order_num+1);r++)
                {
                    B.row(r)=e.derivatives()(r).derivatives().transpose();
                }
                hessian.block<3*(order_num+1),3*(order_num+1)>(3*init,3*init)+=B;
            }
        }
        
        std::cout<<"dmin:"<<dmin<<std::endl;
    }

    static void barrier_gradient(const Data& spline, 
                                 Eigen::VectorXd& grad, Eigen::MatrixXd& hessian,
                                 const Eigen::MatrixXd & V,const Eigen::MatrixXi& F,
                                 BVH& bvh)
    {
        
        int num=3*trajectory_num;
        
        grad.resize(num);
        grad.setZero();

        hessian.resize(num,num);
        hessian.setZero();

        //Eigen::VectorXd auto_grad=grad;
        //Eigen::MatrixXd auto_hessian=hessian;

        bvh.InitTrajectory(spline);

        std::vector<std::pair<unsigned int, unsigned int>> collision_pair, temp_pair;

        clock_t time1 = clock();
        bvh.CheckCollision(collision_pair,offset+margin);
        clock_t time2 = clock();
        std::cout<<"time_bvh:"<<(time2-time1)/(CLOCKS_PER_SEC/1000)<<std::endl<<std::endl;

        std::cout<<"bvh_size:"<<collision_pair.size()<<std::endl;
        double dmin=1.0;

        std::vector<std::vector<std::vector<int>>> segment_lists;

        std::vector<std::vector<std::vector<int>>> face_lists;

        std::vector<std::vector<int>> segment_ob_lists;

        std::vector<std::vector<int>> face_ob_lists;

        segment_lists.resize(subdivide_tree.size());

        face_lists.resize(subdivide_tree.size());

        segment_ob_lists.resize(subdivide_tree.size());

        face_ob_lists.resize(subdivide_tree.size());


        for(unsigned int i=0;i<collision_pair.size();i++)
        {            
            int tr_id=collision_pair[i].first;
            int ob_id=collision_pair[i].second;

            int sp_id=std::get<0>(subdivide_tree[tr_id]);
            double weight=std::get<1>(subdivide_tree[tr_id]).second-std::get<1>(subdivide_tree[tr_id]).first;
            Eigen::MatrixXd basis=std::get<2>(subdivide_tree[tr_id]);

            //std::cout<<"sp_id:"<<sp_id<<std::endl;
            Eigen::MatrixXd bz;
            int init=sp_id*(order_num-2);
            bz=spline.block<order_num+1,3>(init,0);

            std::vector<Eigen::RowVector3d> P(order_num+1);
            for(int j=0;j<=order_num;j++)
            {   
                P[j].setZero();
                for(int j0=0;j0<=order_num;j0++)
                {
                    P[j]+=basis(j,j0)*bz.row(j0);
                }
            }

            std::vector<Eigen::RowVector3d> S(3); 
            for(int j=0;j<3;j++)
            {
                S[j]=V.row(F(ob_id,j));
            }
            
            Eigen::RowVector3d C,C0,C1;
            double d;
        
            //point face
            
            for(int j=0;j<=order_num;j++)
            {
    
                Eigen::RowVector3d d_p;
                Eigen::Matrix3d h_p;
                Distance_der::point_face(P[j],S[0],S[1],S[2], d, d_p, h_p);
                
                d-=offset;
                if(d<margin)
                {
                    if(d<dmin)
                        dmin=d;
                    //-weight*(d-margin)*(d-margin)*log(d/margin)
                    //double e1=-weight*(2*(d-margin)*log(d/margin)+(d-margin)*(d-margin)/d);
                    double e1=-weight*((d*d-margin*margin)/(d*d)*log(d/margin)+(d-margin)*(d-margin)/(d*d));

                    Eigen::Matrix3d I; I.setIdentity();
                    Eigen::MatrixXd A=Eigen::kroneckerProduct(basis.row(j),I);
                    //std::cout<<A<<"\n";
                    Eigen::MatrixXd d_x=d_p*A;
                    
                    grad.segment(3*init,3*(order_num+1)) += e1*d_x.transpose();
                    
                    //double e2=-weight*(2*log(d/margin)+4*(d-margin)/d-(d-margin)*(d-margin)/(d*d));
                    double e2=-weight*(1/d+margin*margin/(d*d*d)*(2*log(d/margin)-1)+2*margin*(d-margin)/(d*d*d));

                    
                    hessian.block<3*(order_num+1),3*(order_num+1)>(3*init,3*init)+=e2*d_x.transpose()*d_x+e1*A.transpose()*h_p*A;
                }
            }

            //face point
            for(int j0=0;j0<order_num-1;j0++)
            {
                for(int j1=j0+1;j1<order_num;j1++)
                {
                    for(int j2=j1+1;j2<=order_num;j2++)
                    {
                        for(int k=0;k<3;k++)
                        {
                            Distance<double,3>::face_point(P[j0],P[j1],P[j2],
                                                            S[k], d, C);
                            d-=offset;
                            
                            if(d<margin)
                            { 
                               if(d<dmin)
                                  dmin=d;
                                std::vector<int> face_list; face_list.clear();
                                face_list.push_back(j0);
                                face_list.push_back(j1);
                                face_list.push_back(j2);
                                face_list.push_back(k);
                                
                                face_lists[tr_id].push_back(face_list);

                                face_ob_lists[tr_id].push_back(ob_id);  
                            }
                        }
                    }
                }
            }

            //segment segment

            for(int j0=0;j0<order_num;j0++)
            {
                for(int j1=j0+1;j1<=order_num;j1++)
                {
                    for(int k=0;k<3;k++)
                    {
                        int k1=(k+1)%3;
                        Distance<double,3>::segment_segment(P[j0],P[j1],
                                                            S[k],S[k1],d,C0,C1);
                        d-=offset;
                            
                        if(d<margin)
                        { 
                            if(d<dmin)
                               dmin=d;
                            std::vector<int> segment_list; segment_list.clear();
                            segment_list.push_back(j0);
                            segment_list.push_back(j1);
                            segment_list.push_back(k);

                            segment_lists[tr_id].push_back(segment_list);
                            segment_ob_lists[tr_id].push_back(ob_id);  
                        }
                    }
                }
            }
            
        }
        
        for(unsigned int tr_id=0;tr_id<subdivide_tree.size();tr_id++)
        {

            if(segment_lists[tr_id].size()==0)
            {
                continue;
            }

            int sp_id=std::get<0>(subdivide_tree[tr_id]);
            double weight=std::get<1>(subdivide_tree[tr_id]).second-std::get<1>(subdivide_tree[tr_id]).first;
            Eigen::MatrixXd basis=std::get<2>(subdivide_tree[tr_id]);

            //std::cout<<"sp_id:"<<sp_id<<std::endl;
            int init=sp_id*(order_num-2);
            
            Eigen::MatrixXd x=spline.block<order_num+1,3>(init,0);
            Vec12 X;
            X.resize(3*(order_num+1));
            for(int r=0;r<=order_num;r++)
            {
                X(3*r).value() = x(r,0);
                X(3*r+1).value() = x(r,1);
                X(3*r+2).value() = x(r,2);
            }
            //repeat partial derivatives for the inner AutoDiffScalar
            for(int id=0;id<3*(order_num+1);id++)
            {
                X(id).derivatives().resize(3*(order_num+1));
                X(id).derivatives().setZero();
                X(id).derivatives()(id)= 1;
                X(id).value().derivatives() = inner_derivative_t::Unit(3*(order_num+1),id);
            }
            //set the hessian matrix to zero
            for(int idx=0; idx<3*(order_num+1); idx++) {
                for(int id=0;id<3*(order_num+1);id++)
                {
                    X(id).derivatives()(idx).derivatives()  = inner_derivative_t::Zero(3*(order_num+1));
                }       
            }

            std::vector<Vec3> P(order_num+1); 

            for(int j=0;j<=order_num;j++)
            {
                P[j].setZero();
                for(int j0=0;j0<=order_num;j0++)
                {
                    P[j](0)+=basis(j,j0)*X(3*j0);
                    P[j](1)+=basis(j,j0)*X(3*j0+1);
                    P[j](2)+=basis(j,j0)*X(3*j0+2);
                }
            }
            
            Vec3 C,C0,C1;
            scalar_t d;

                

            //face point
            
            std::vector<std::vector<int>> face_list=face_lists[tr_id];    
        
            for(unsigned int j=0;j<face_list.size();j++)
            {
                int ob_id=face_ob_lists[tr_id][j];
                std::vector<Eigen::RowVector3d> S(3); 
                for(int i=0;i<3;i++)
                {
                    S[i]=V.row(F(ob_id,i));
                }
                int j0=face_list[j][0];
                int j1=face_list[j][1];
                int j2=face_list[j][2];
                int k=face_list[j][3];
          
                Distance<scalar_t,3>::face_point(P[j0],P[j1],P[j2],
                                                            S[k], d, C);
                d=d-offset;
                //scalar_t e=-weight*(d-margin)*(d-margin)*log(d/margin);
                
                scalar_t e=-weight*(d-margin)*(d-margin)/d*log(d/margin);///d

                //scalar_t e=weight*(1-d/margin*d/margin)*(1-d/margin*d/margin);
                grad.segment(3*init,3*(order_num+1)) += e.value().derivatives();
                Eigen::Matrix<double,3*(order_num+1),3*(order_num+1)> B;

                for(int r=0;r<3*(order_num+1);r++)
                {
                    B.row(r)=e.derivatives()(r).derivatives().transpose();
                }
                hessian.block<3*(order_num+1),3*(order_num+1)>(3*init,3*init)+=B;
            }
            
               

            //segment segment

            std::vector<std::vector<int>> segment_list=segment_lists[tr_id];     
        
            for(unsigned int j=0;j<segment_list.size();j++)
            {
                int ob_id=segment_ob_lists[tr_id][j];
                std::vector<Eigen::RowVector3d> S(3); 
                for(int i=0;i<3;i++)
                {
                    S[i]=V.row(F(ob_id,i));
                }
                int j0=segment_list[j][0];
                int j1=segment_list[j][1];
                int k=segment_list[j][2];
                Distance<scalar_t,3>::segment_segment(P[j0],P[j1],
                                                        S[k],S[(k+1)%3],d,C0,C1);
                d=d-offset;
                //scalar_t e=-weight*(d-margin)*(d-margin)*log(d/margin);
                
                scalar_t e=-weight*(d-margin)*(d-margin)/d*log(d/margin);///d

                //scalar_t e=weight*(1-d/margin*d/margin)*(1-d/margin*d/margin);
                grad.segment(3*init,3*(order_num+1)) += e.value().derivatives();
                Eigen::Matrix<double,3*(order_num+1),3*(order_num+1)> B;

                for(int r=0;r<3*(order_num+1);r++)
                {
                    B.row(r)=e.derivatives()(r).derivatives().transpose();
                }
                hessian.block<3*(order_num+1),3*(order_num+1)>(3*init,3*init)+=B;
            }
        }
        
        std::cout<<"dmin:"<<dmin<<std::endl;
    }

    static void fast_point_barrier_gradient(const Data& spline, 
                                       Eigen::VectorXd& grad, Eigen::MatrixXd& hessian,
                                       const Eigen::MatrixXd & V,
                                       BVH& bvh)
    {
        
        int num=3*trajectory_num;
        
        grad.resize(num);
        grad.setZero();

        hessian.resize(num,num);
        hessian.setZero();

        bvh.InitTrajectory(spline);

        std::vector<std::pair<unsigned int, unsigned int>> collision_pair, temp_pair;

        clock_t time1 = clock();
        bvh.pcCheckCollision(collision_pair,offset+margin);
        clock_t time2 = clock();
        std::cout<<"time_bvh:"<<(time2-time1)/(CLOCKS_PER_SEC/1000)<<std::endl<<std::endl;

        std::cout<<"bvh_size:"<<collision_pair.size()<<std::endl;
        double dmin=1.0;

        std::vector<std::vector<int>> segment_lists;

        std::vector<std::vector<int>> segment_ob_lists;

        segment_lists.resize(subdivide_tree.size());

        segment_ob_lists.resize(subdivide_tree.size());


        for(unsigned int i=0;i<collision_pair.size();i++)
        {            
            int tr_id=collision_pair[i].first;
            int ob_id=collision_pair[i].second;

            int sp_id=std::get<0>(subdivide_tree[tr_id]);
            double weight=std::get<1>(subdivide_tree[tr_id]).second-std::get<1>(subdivide_tree[tr_id]).first;
            Eigen::MatrixXd basis=std::get<2>(subdivide_tree[tr_id]);

            //std::cout<<"sp_id:"<<sp_id<<std::endl;
            Eigen::MatrixXd bz;
            int init=sp_id*(order_num-2);
            bz=spline.block<order_num+1,3>(init,0);

            std::vector<Eigen::RowVector3d> P(order_num+1);
            for(int j=0;j<=order_num;j++)
            {   
                P[j].setZero();
                for(int j0=0;j0<=order_num;j0++)
                {
                    P[j]+=basis(j,j0)*bz.row(j0);
                }
            }

            Eigen::RowVector3d S=V.row(ob_id); 
            
            
            Eigen::RowVector3d C,C0,C1;
            double d;
        
            //point point
            
            for(int j=0;j<=order_num;j++)
            {
    
                Eigen::RowVector3d d_p;
                Eigen::Matrix3d h_p;
                Distance_der::point_point(P[j],S, d, d_p, h_p);
                
                d-=offset;
                if(d<margin)
                {
                    if(d<dmin)
                        dmin=d;
                    //-weight*(d-margin)*(d-margin)*log(d/margin)
                    //double e1=-weight*(2*(d-margin)*log(d/margin)+(d-margin)*(d-margin)/d);
                    double e1=-weight*((d*d-margin*margin)/(d*d)*log(d/margin)+(d-margin)*(d-margin)/(d*d));

                    Eigen::Matrix3d I; I.setIdentity();
                    Eigen::MatrixXd A=Eigen::kroneckerProduct(basis.row(j),I);
                    //std::cout<<A<<"\n";
                    Eigen::MatrixXd d_x=d_p*A;
                    
                    grad.segment(3*init,3*(order_num+1)) += e1*d_x.transpose();
                    
                    //double e2=-weight*(2*log(d/margin)+4*(d-margin)/d-(d-margin)*(d-margin)/(d*d));
                    double e2=-weight*(1/d+margin*margin/(d*d*d)*(2*log(d/margin)-1)+2*margin*(d-margin)/(d*d*d));
                    
                    hessian.block<3*(order_num+1),3*(order_num+1)>(3*init,3*init)+=e2*d_x.transpose()*d_x+e1*A.transpose()*h_p*A;
                }
            }
            
            //segment segment
            
            //for(int j=0;j<=order_num;j++)
            {
                int j=order_num;
                
                Distance<double,3>::segment_constpoint(P[j],P[(j+1)%(order_num+1)],
                                                        S,d,C);
                //std::cout<<"ss:"<<d<<std::endl;
                d-=offset;
                if(d<margin)
                {
                    segment_lists[tr_id].push_back(j);
                    segment_ob_lists[tr_id].push_back(ob_id);
                    if(d<dmin)
                        dmin=d;
                }
                
            }

        }
        
        for(unsigned int tr_id=0;tr_id<subdivide_tree.size();tr_id++)
        {

            if(segment_lists[tr_id].size()==0)
            {
                continue;
            }

            int sp_id=std::get<0>(subdivide_tree[tr_id]);
            double weight=std::get<1>(subdivide_tree[tr_id]).second-std::get<1>(subdivide_tree[tr_id]).first;
            Eigen::MatrixXd basis=std::get<2>(subdivide_tree[tr_id]);

            //std::cout<<"sp_id:"<<sp_id<<std::endl;
            int init=sp_id*(order_num-2);
            
            Eigen::MatrixXd x=spline.block<order_num+1,3>(init,0);
            Vec12 X;
            X.resize(3*(order_num+1));
            for(int r=0;r<=order_num;r++)
            {
                X(3*r).value() = x(r,0);
                X(3*r+1).value() = x(r,1);
                X(3*r+2).value() = x(r,2);
            }
            //repeat partial derivatives for the inner AutoDiffScalar
            for(int id=0;id<3*(order_num+1);id++)
            {
                X(id).derivatives().resize(3*(order_num+1));
                X(id).derivatives().setZero();
                X(id).derivatives()(id)= 1;
                X(id).value().derivatives() = inner_derivative_t::Unit(3*(order_num+1),id);
            }
            //set the hessian matrix to zero
            for(int idx=0; idx<3*(order_num+1); idx++) {
                for(int id=0;id<3*(order_num+1);id++)
                {
                    X(id).derivatives()(idx).derivatives()  = inner_derivative_t::Zero(3*(order_num+1));
                }       
            }

            std::vector<Vec3> P(order_num+1); 

            for(int j=0;j<=order_num;j++)
            {
                P[j].setZero();
                for(int j0=0;j0<=order_num;j0++)
                {
                    P[j](0)+=basis(j,j0)*X(3*j0);
                    P[j](1)+=basis(j,j0)*X(3*j0+1);
                    P[j](2)+=basis(j,j0)*X(3*j0+2);
                }
            }
            
            Vec3 C,C0,C1;
            scalar_t d;

            //segment segment
            
            std::vector<int> segment_list=segment_lists[tr_id];     
            for(unsigned int j=0;j<segment_list.size();j++)
            {
                int ob_id=segment_ob_lists[tr_id][j];
                Eigen::RowVector3d S=V.row(ob_id); 
                
                int j0=segment_list[j];
                Distance<scalar_t,3>::segment_constpoint(P[j0],P[(j0+1)%(order_num+1)],
                                                         S,d,C);
                d=d-offset;
                scalar_t e=-weight*(d-margin)*(d-margin)/d*log(d/margin);
                //scalar_t e=weight*(1-d/margin*d/margin)*(1-d/margin*d/margin);
                grad.segment(3*init,3*(order_num+1)) += e.value().derivatives();
                Eigen::Matrix<double,3*(order_num+1),3*(order_num+1)> B;

                for(int r=0;r<3*(order_num+1);r++)
                {
                    B.row(r)=e.derivatives()(r).derivatives().transpose();
                }
                hessian.block<3*(order_num+1),3*(order_num+1)>(3*init,3*init)+=B;
            }
        }
        
        std::cout<<"dmin:"<<dmin<<std::endl;
    }
    
    static void point_barrier_gradient(const Data& spline, 
                                       Eigen::VectorXd& grad, Eigen::MatrixXd& hessian,
                                       const Eigen::MatrixXd & V,
                                       BVH& bvh)
    {
        
        int num=3*trajectory_num;
        
        grad.resize(num);
        grad.setZero();

        hessian.resize(num,num);
        hessian.setZero();

        //Eigen::VectorXd auto_grad=grad;
        //Eigen::MatrixXd auto_hessian=hessian;

        bvh.InitTrajectory(spline);

        std::vector<std::pair<unsigned int, unsigned int>> collision_pair, temp_pair;

        clock_t time1 = clock();
        bvh.pcCheckCollision(collision_pair,offset+margin);
        clock_t time2 = clock();
        std::cout<<"time_bvh:"<<(time2-time1)/(CLOCKS_PER_SEC/1000)<<std::endl<<std::endl;

        std::cout<<"bvh_size:"<<collision_pair.size()<<std::endl;
        double dmin=1.0;


        std::vector<std::vector<std::vector<int>>> face_lists;

        std::vector<std::vector<int>> face_ob_lists;

        face_lists.resize(subdivide_tree.size());

        face_ob_lists.resize(subdivide_tree.size());


        for(unsigned int i=0;i<collision_pair.size();i++)
        {            
            int tr_id=collision_pair[i].first;
            int ob_id=collision_pair[i].second;

            int sp_id=std::get<0>(subdivide_tree[tr_id]);
            //double weight=std::get<1>(subdivide_tree[tr_id]).second-std::get<1>(subdivide_tree[tr_id]).first;
            Eigen::MatrixXd basis=std::get<2>(subdivide_tree[tr_id]);

            //std::cout<<"sp_id:"<<sp_id<<std::endl;
            Eigen::MatrixXd bz;
            int init=sp_id*(order_num-2);
            bz=spline.block<order_num+1,3>(init,0);

            std::vector<Eigen::RowVector3d> P(order_num+1);
            for(int j=0;j<=order_num;j++)
            {   
                P[j].setZero();
                for(int j0=0;j0<=order_num;j0++)
                {
                    P[j]+=basis(j,j0)*bz.row(j0);
                }
            }

            Eigen::RowVector3d S=V.row(ob_id); 
            
            
            Eigen::RowVector3d C,C0,C1;
            double d;

           // face point

           for(int j0=0;j0<order_num-1;j0++)
            {
                for(int j1=j0+1;j1<order_num;j1++)
                {
                    for(int j2=j1+1;j2<=order_num;j2++)
                    {
                        
                        Distance<double,3>::face_point(P[j0],P[j1],P[j2],
                                                        S, d, C);
                        d-=offset;
                        
                        if(d<margin)
                        { 
                            if(d<dmin)
                                dmin=d;
                            std::vector<int> face_list; face_list.clear();
                            face_list.push_back(j0);
                            face_list.push_back(j1);
                            face_list.push_back(j2);
                            
                            face_lists[tr_id].push_back(face_list);

                            face_ob_lists[tr_id].push_back(ob_id);  
                        }
                        
                    }
                }
            }
           
        }
        
        for(unsigned int tr_id=0;tr_id<subdivide_tree.size();tr_id++)
        {

            if(face_lists[tr_id].size()==0)
            {
                continue;
            }

            int sp_id=std::get<0>(subdivide_tree[tr_id]);
            double weight=std::get<1>(subdivide_tree[tr_id]).second-std::get<1>(subdivide_tree[tr_id]).first;
            Eigen::MatrixXd basis=std::get<2>(subdivide_tree[tr_id]);

            //std::cout<<"sp_id:"<<sp_id<<std::endl;
            int init=sp_id*(order_num-2);
            
            Eigen::MatrixXd x=spline.block<order_num+1,3>(init,0);
            Vec12 X;
            X.resize(3*(order_num+1));
            for(int r=0;r<=order_num;r++)
            {
                X(3*r).value() = x(r,0);
                X(3*r+1).value() = x(r,1);
                X(3*r+2).value() = x(r,2);
            }
            //repeat partial derivatives for the inner AutoDiffScalar
            for(int id=0;id<3*(order_num+1);id++)
            {
                X(id).derivatives().resize(3*(order_num+1));
                X(id).derivatives().setZero();
                X(id).derivatives()(id)= 1;
                X(id).value().derivatives() = inner_derivative_t::Unit(3*(order_num+1),id);
            }
            //set the hessian matrix to zero
            for(int idx=0; idx<3*(order_num+1); idx++) {
                for(int id=0;id<3*(order_num+1);id++)
                {
                    X(id).derivatives()(idx).derivatives()  = inner_derivative_t::Zero(3*(order_num+1));
                }       
            }

            std::vector<Vec3> P(order_num+1); 

            for(int j=0;j<=order_num;j++)
            {
                P[j].setZero();
                for(int j0=0;j0<=order_num;j0++)
                {
                    P[j](0)+=basis(j,j0)*X(3*j0);
                    P[j](1)+=basis(j,j0)*X(3*j0+1);
                    P[j](2)+=basis(j,j0)*X(3*j0+2);
                }
            }
            
            Vec3 C,C0,C1;
            scalar_t d;

            //face point
            std::vector<std::vector<int>> face_list=face_lists[tr_id];    
        
            for(unsigned int j=0;j<face_list.size();j++)
            {
                int ob_id=face_ob_lists[tr_id][j];
                Eigen::RowVector3d S=V.row(ob_id);; 
                
                int j0=face_list[j][0];
                int j1=face_list[j][1];
                int j2=face_list[j][2];
          
                Distance<scalar_t,3>::face_point(P[j0],P[j1],P[j2],
                                                 S, d, C);
                d=d-offset;
                //scalar_t e=-weight*(d-margin)*(d-margin)*log(d/margin);
                
                scalar_t e=-weight*(d-margin)*(d-margin)/d*log(d/margin);///d

                //scalar_t e=weight*(1-d/margin*d/margin)*(1-d/margin*d/margin);
                grad.segment(3*init,3*(order_num+1)) += e.value().derivatives();
                Eigen::Matrix<double,3*(order_num+1),3*(order_num+1)> B;

                for(int r=0;r<3*(order_num+1);r++)
                {
                    B.row(r)=e.derivatives()(r).derivatives().transpose();
                }
                hessian.block<3*(order_num+1),3*(order_num+1)>(3*init,3*init)+=B;
            }
        }
        
        std::cout<<"dmin:"<<dmin<<std::endl;
    }

    static void bound_gradient(const Data& spline, const double& piece_time,
                               Eigen::VectorXd& grad, Eigen::MatrixXd& hessian,
                               double& g_t, double& h_t, 
                               Eigen::VectorXd& partgrad)
    {
        int num=3*trajectory_num;
        
        grad.resize(num);
        grad.setZero();

        hessian.resize(num,num);
        hessian.setZero();

        g_t=0;
        h_t=0;

        partgrad.resize(num);
        partgrad.setZero();
    
        double max_vel=0;
        double max_acc=0;

        for(unsigned int tr_id=0;tr_id<vel_tree.size();tr_id++)
        {
        
            int sp_id=std::get<0>(vel_tree[tr_id]);
            double weight=std::get<1>(vel_tree[tr_id]).second-std::get<1>(vel_tree[tr_id]).first;
            Eigen::MatrixXd basis=std::get<2>(vel_tree[tr_id]);
            
            int init=sp_id*(order_num-2);
            
            Eigen::MatrixXd bz;
            bz=spline.block<order_num+1,3>(init,0);

            std::vector<Eigen::RowVector3d> P(order_num+1);
            for(int j=0;j<=order_num;j++)
            {
                P[j].setZero();
                for(int j0=0;j0<=order_num;j0++)
                {
                  P[j]+=basis(j,j0)*bz.row(j0);
                } 
            }
                
            double d;
           
            for(int j=0;j<order_num;j++)
            {
                Eigen::RowVector3d P_=P[j+1]-P[j];
                Eigen::RowVector3d vel=order_num*P_;
                double v=vel.norm()/(weight);
                //d=vel_limit*piece_time-vel.norm()/weight;
                d=vel_limit-v/piece_time;

                if(v/piece_time>max_vel)
                  max_vel=v/piece_time;
                
                if(d<margin)
                { 
                   //g_t, h_t
                   g_t+=weight*v/(piece_time*piece_time)*(-2*(d-margin)*log(d/margin)-(d-margin)*(d-margin)/d);

                   h_t+=weight*(-2*v/std::pow(piece_time,3)*(-2*(d-margin)*log(d/margin)-(d-margin)*(d-margin)/d)+
                            v*v/std::pow(piece_time,4)*(-2*log(d/margin)-4*(d-margin)/d+(d-margin)*(d-margin)/(d*d)));
                    
                   //grad, hessian
                   //-weight*(d-margin)*(d-margin)*log(d/margin)
                    Eigen::Matrix3d I; I.setIdentity();
                    Eigen::MatrixXd A=Eigen::kroneckerProduct(basis.row(j+1),I)-
                                      Eigen::kroneckerProduct(basis.row(j),I);
                    //std::cout<<A<<"\n";
                    Eigen::RowVector3d d_p;
                    Eigen::Matrix3d h_p;

                    double d_=P_.norm();
                    
                    d_p=-order_num/(weight*piece_time)*P_/d_;
        
                    h_p=-order_num/(weight*piece_time)*(I/d_-P_.transpose()*P_/std::pow(d_,3));

                    Eigen::MatrixXd d_x=d_p*A;

                    double e1=-weight*(2*(d-margin)*log(d/margin)+(d-margin)*(d-margin)/d);
                    
                    grad.segment(3*init,3*(order_num+1)) += e1*d_x.transpose();
                    
                    double e2=-weight*(2*log(d/margin)+4*(d-margin)/d-(d-margin)*(d-margin)/(d*d));
                    
                    hessian.block<3*(order_num+1),3*(order_num+1)>(3*init,3*init)+=e2*d_x.transpose()*d_x+e1*A.transpose()*h_p*A;

                    //partgrad
  
                    //weight/piece_time*(d-vel_limit)*(2*(d-margin)*log(d/margin)+(d-margin)*(d-margin)/d);

                    double e3=weight/piece_time*(2*(d-margin)*log(d/margin)+(d-margin)*(d-margin)/d)
                             +weight/piece_time*(d-vel_limit)*(2*log(d/margin)+4*(d-margin)/d-(d-margin)*(d-margin)/(d*d));

                    partgrad.segment(3*init,3*(order_num+1)) += e3*d_x.transpose();

                }
            }
        }

        
        for(unsigned int tr_id=0;tr_id<acc_tree.size();tr_id++)
        {
        
            int sp_id=std::get<0>(acc_tree[tr_id]);
            double weight=std::get<1>(acc_tree[tr_id]).second-std::get<1>(acc_tree[tr_id]).first;
            Eigen::MatrixXd basis=std::get<2>(acc_tree[tr_id]);
            
            int init=sp_id*(order_num-2);

            Eigen::MatrixXd bz;
            bz=spline.block<order_num+1,3>(init,0);

            std::vector<Eigen::RowVector3d> P(order_num+1);
            for(int j=0;j<=order_num;j++)
            {
                P[j].setZero();
                for(int j0=0;j0<=order_num;j0++)
                {
                  P[j]+=basis(j,j0)*bz.row(j0);
                } 
            }
                
            double d;
           
            for(int j=0;j<order_num-1;j++)
            {
                Eigen::RowVector3d P_=P[j+2]-2*P[j+1]+P[j];
                Eigen::RowVector3d acc=order_num*(order_num-1)*P_;
                double a=acc.norm()/(weight*weight);
                
                //d=acc_limit*piece_time*piece_time-acc.norm()/(weight*weight);
                d=acc_limit-a/(piece_time*piece_time);

                if(a/(piece_time*piece_time)>max_acc)
                  max_acc=a/(piece_time*piece_time);
                
                if(d<margin)
                { 
                   //g_t, h_t
                   g_t+=weight*2*a/std::pow(piece_time,3)*(-2*(d-margin)*log(d/margin)-(d-margin)*(d-margin)/d);

                   h_t+=weight*(-6*a/std::pow(piece_time,4)*(-2*(d-margin)*log(d/margin)-(d-margin)*(d-margin)/d)
                                +4*a*a/std::pow(piece_time,6)*(-2*log(d/margin)-4*(d-margin)/d+(d-margin)*(d-margin)/(d*d)));

                   //grad, hessian
                   //-weight*(d-margin)*(d-margin)*log(d/margin)
                    Eigen::Matrix3d I; I.setIdentity();
                    Eigen::MatrixXd A=Eigen::kroneckerProduct(basis.row(j+2),I)-
                                      2*Eigen::kroneckerProduct(basis.row(j+1),I)+
                                      Eigen::kroneckerProduct(basis.row(j),I);
                    //std::cout<<A<<"\n";
                    Eigen::RowVector3d d_p;
                    Eigen::Matrix3d h_p;

                    double d_=P_.norm();
                    
                    d_p=-order_num*(order_num-1)/std::pow(weight*piece_time,2)*P_/d_;
        
                    h_p=-order_num*(order_num-1)/std::pow(weight*piece_time,2)*(I/d_-P_.transpose()*P_/std::pow(d_,3));

                    Eigen::MatrixXd d_x=d_p*A;

                    double e1=-weight*(2*(d-margin)*log(d/margin)+(d-margin)*(d-margin)/d);
                    
                    grad.segment(3*init,3*(order_num+1)) += e1*d_x.transpose();
                    
                    double e2=-weight*(2*log(d/margin)+4*(d-margin)/d-(d-margin)*(d-margin)/(d*d));
                    
                    hessian.block<3*(order_num+1),3*(order_num+1)>(3*init,3*init)+=e2*d_x.transpose()*d_x+e1*A.transpose()*h_p*A;

                    //partgrad

                    //2*weight/piece_time*(d-acc_limit)*(2*(d-margin)*log(d/margin)+(d-margin)*(d-margin)/d);

                    double e3=2*weight/piece_time*(2*(d-margin)*log(d/margin)+(d-margin)*(d-margin)/d)
                             +2*weight/piece_time*(d-acc_limit)*(2*log(d/margin)+4*(d-margin)/d-(d-margin)*(d-margin)/(d*d));

                    partgrad.segment(3*init,3*(order_num+1)) += e3*d_x.transpose();

                }
            }
        }

        std::cout<<"max_vel:"<<max_vel<<std::endl;
        std::cout<<"max_acc:"<<max_acc<<std::endl;
    }
    
};

PRJ_END

#endif