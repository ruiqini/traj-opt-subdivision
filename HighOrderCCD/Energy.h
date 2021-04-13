#ifndef ENERGY_H
#define ENERGY_H

#include "Config/Config.h"

#include "HighOrderCCD/BVH/BVH.h"
#include "HighOrderCCD/CCD/CCD.h"
#include "HighOrderCCD/Distance.h"
#include "HighOrderCCD/Distance_der.h"

#include <vector>

PRJ_BEGIN

class Energy
{
  public:
    typedef Eigen::MatrixXd Data;
    typedef std::vector< std::tuple< int, std::pair<double,double>, Eigen::MatrixXd > > Tree;
    typedef std::tuple< int, std::pair<double,double>, Eigen::MatrixXd >  Node;
    
    static double whole_energy(const Data& spline, const double& piece_time,
                       const Eigen::MatrixXd & V,const Eigen::MatrixXi& F,
                       BVH& bvh)
    { 
        return ks*dynamic_energy(spline,piece_time)+
               lambda*barrier_energy(spline,V,F,bvh)+ lambda*bound_energy(spline,piece_time)+
               kt*whole_weight*piece_time;
    }

    static double fast_whole_energy(const Data& spline, const double& piece_time,
                       const Eigen::MatrixXd & V,const Eigen::MatrixXi& F,
                       BVH& bvh)
    { 
        return ks*dynamic_energy(spline,piece_time)+
               lambda*fast_barrier_energy(spline,V,F,bvh)+ lambda*bound_energy(spline,piece_time)+
               kt*whole_weight*piece_time;
    }

    static double fast_point_whole_energy(const Data& spline, const double& piece_time,
                                        const Eigen::MatrixXd & V,
                                        BVH& bvh)
    { 
        return ks*dynamic_energy(spline,piece_time)+
               lambda*fast_point_barrier_energy(spline,V,bvh)+ lambda*bound_energy(spline,piece_time)+
               kt*whole_weight*piece_time;
    }

    static double point_whole_energy(const Data& spline, const double& piece_time,
                                        const Eigen::MatrixXd & V,
                                        BVH& bvh)
    { 
        return ks*dynamic_energy(spline,piece_time)+
               lambda*point_barrier_energy(spline,V,bvh)+ lambda*bound_energy(spline,piece_time)+
               kt*whole_weight*piece_time;
    }

    static double dynamic_energy(const Data& spline, const double& piece_time)
    {
        double energy=0;

        //dynamic energy
        Eigen::MatrixXd M;
        for(int sp_id=0;sp_id<piece_num;sp_id++)
        {
            Eigen::MatrixXd bz;
            bz=spline.block<order_num+1,3>(sp_id*(order_num-2),0);
            
            M=convert_list[sp_id].transpose()*M_dynamic*convert_list[sp_id];
            for(int j=0;j<3;j++)
            {
                Eigen::VectorXd x=bz.col(j);
                energy+=1/std::pow(time_weight[sp_id]*piece_time,2*der_num-1)*0.5*x.transpose()*M*x;
            }
        }
        return energy;
    }

    static double fast_barrier_energy(const Data& spline, 
                                 const Eigen::MatrixXd & V,const Eigen::MatrixXi& F,
                                 BVH& bvh)
    {
        double energy=0;
        bvh.InitTrajectory(spline);

        //bvh.UpdateTrajectory(spline);

        std::vector<std::pair<unsigned int, unsigned int>> collision_pair;
        //bvh.CheckCollision(collision_pair,margin);
        bvh.CheckCollision(collision_pair,offset+margin);

        for(unsigned int i=0;i<collision_pair.size();i++)
        {
            int tr_id=collision_pair[i].first;
            int ob_id=collision_pair[i].second;
        
            int sp_id=std::get<0>(subdivide_tree[tr_id]);
            double weight=std::get<1>(subdivide_tree[tr_id]).second-std::get<1>(subdivide_tree[tr_id]).first;
            Eigen::MatrixXd basis=std::get<2>(subdivide_tree[tr_id]);
            
            Eigen::MatrixXd bz;
            bz=spline.block<order_num+1,3>(sp_id*(order_num-2),0);
        
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
                Distance<double,3>::point_face(P[j],S[0],S[1],S[2], d, C);
                d-=offset;
                if(d<=0)
                    return INFINITY;
                
                if(d<margin)
                { 
                   energy+=-weight*(d-margin)*(d-margin)/d*log(d/margin); ///d

                //energy+=weight*(1-d/margin*d/margin)*(1-d/margin*d/margin);   
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
                    d-=offset;
                    if(d<=0)
                        return INFINITY;
                    if(d<margin)
                    {
                        energy+=-weight*(d-margin)*(d-margin)/d*log(d/margin);///d

                        //energy+=weight*(1-d/margin*d/margin)*(1-d/margin*d/margin);
                    }
                }
            }
           
        }

        return energy;  
    }

    static double barrier_energy(const Data& spline, 
                                 const Eigen::MatrixXd & V,const Eigen::MatrixXi& F,
                                 BVH& bvh)
    {
        double energy=0;
        bvh.InitTrajectory(spline);

        //bvh.UpdateTrajectory(spline);

        std::vector<std::pair<unsigned int, unsigned int>> collision_pair;
        //bvh.CheckCollision(collision_pair,margin);
        bvh.CheckCollision(collision_pair,offset+margin);

        for(unsigned int i=0;i<collision_pair.size();i++)
        {
            int tr_id=collision_pair[i].first;
            int ob_id=collision_pair[i].second;
        
            int sp_id=std::get<0>(subdivide_tree[tr_id]);
            double weight=std::get<1>(subdivide_tree[tr_id]).second-std::get<1>(subdivide_tree[tr_id]).first;
            Eigen::MatrixXd basis=std::get<2>(subdivide_tree[tr_id]);
            
            Eigen::MatrixXd bz;
            bz=spline.block<order_num+1,3>(sp_id*(order_num-2),0);
        
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
                Distance<double,3>::point_face(P[j],S[0],S[1],S[2], d, C);
                d-=offset;
                if(d<=0)
                {
                    return INFINITY;
                }
                    
                
                if(d<margin)
                { 
                   energy+=-weight*(d-margin)*(d-margin)/d*log(d/margin); ///d

                //energy+=weight*(1-d/margin*d/margin)*(1-d/margin*d/margin);   
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
                            if(d<=0)
                                return INFINITY;
                            if(d<margin)
                            { 
                               energy+=-weight*(d-margin)*(d-margin)/d*log(d/margin); ///d

                               //energy+=weight*(1-d/margin*d/margin)*(1-d/margin*d/margin);   
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
                        if(d<=0)
                        {
                            return INFINITY;
                        }
                            
                        if(d<margin)
                        { 
                            energy+=-weight*(d-margin)*(d-margin)/d*log(d/margin); ///d

                            //energy+=weight*(1-d/margin*d/margin)*(1-d/margin*d/margin);   
                        }
                    }
                }
            }
           
        }

        return energy;  
    }

    static double fast_point_barrier_energy(const Data& spline, 
                                       const Eigen::MatrixXd & V,
                                       BVH& bvh)
    {
        double energy=0;
        bvh.InitTrajectory(spline);

        //bvh.UpdateTrajectory(spline);

        std::vector<std::pair<unsigned int, unsigned int>> collision_pair;
        //bvh.CheckCollision(collision_pair,margin);
        bvh.pcCheckCollision(collision_pair,offset+margin);

        for(unsigned int i=0;i<collision_pair.size();i++)
        {
            int tr_id=collision_pair[i].first;
            int ob_id=collision_pair[i].second;
        
            int sp_id=std::get<0>(subdivide_tree[tr_id]);
            double weight=std::get<1>(subdivide_tree[tr_id]).second-std::get<1>(subdivide_tree[tr_id]).first;
            Eigen::MatrixXd basis=std::get<2>(subdivide_tree[tr_id]);
            
            Eigen::MatrixXd bz;
            bz=spline.block<order_num+1,3>(sp_id*(order_num-2),0);
        
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
                Distance<double,3>::point_point(P[j],S, d, C);
                d-=offset;
                if(d<=0)
                    return INFINITY;
                
                if(d<margin)
                { 
                   energy+=-weight*(d-margin)*(d-margin)/d*log(d/margin); 
                //energy+=weight*(1-d/margin*d/margin)*(1-d/margin*d/margin);   
                }
            }
            
            //segment segment
            
            //for(int j=0;j<=order_num;j++)
            {
                int j=order_num;
                
                Distance<double,3>::segment_constpoint(P[j],P[(j+1)%(order_num+1)],
                                                        S,d,C);
                d-=offset;
                if(d<=0)
                    return INFINITY;
                if(d<margin)
                {
                    energy+=-weight*(d-margin)*(d-margin)/d*log(d/margin);
                    //energy+=weight*(1-d/margin*d/margin)*(1-d/margin*d/margin);
                }
            }
           
        }

        return energy;  
    }

    static double point_barrier_energy(const Data& spline, 
                                       const Eigen::MatrixXd & V,
                                       BVH& bvh)
    {
        double energy=0;
        bvh.InitTrajectory(spline);

        //bvh.UpdateTrajectory(spline);

        std::vector<std::pair<unsigned int, unsigned int>> collision_pair;
        //bvh.CheckCollision(collision_pair,margin);
        bvh.pcCheckCollision(collision_pair,offset+margin);

        for(unsigned int i=0;i<collision_pair.size();i++)
        {
            int tr_id=collision_pair[i].first;
            int ob_id=collision_pair[i].second;
        
            int sp_id=std::get<0>(subdivide_tree[tr_id]);
            double weight=std::get<1>(subdivide_tree[tr_id]).second-std::get<1>(subdivide_tree[tr_id]).first;
            Eigen::MatrixXd basis=std::get<2>(subdivide_tree[tr_id]);
            
            Eigen::MatrixXd bz;
            bz=spline.block<order_num+1,3>(sp_id*(order_num-2),0);
        
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
                        if(d<=0)
                            return INFINITY;
                        if(d<margin)
                        { 
                            energy+=-weight*(d-margin)*(d-margin)/d*log(d/margin); ///d

                            //energy+=weight*(1-d/margin*d/margin)*(1-d/margin*d/margin);   
                        }
                        
                    }
                }
            }
        }

        return energy;  
    }

    static double bound_energy(const Data& spline,const double& piece_time)
    {
        double energy=0;

        for(unsigned int tr_id=0;tr_id<vel_tree.size();tr_id++)
        {
        
            int sp_id=std::get<0>(vel_tree[tr_id]);
            double weight=std::get<1>(vel_tree[tr_id]).second-std::get<1>(vel_tree[tr_id]).first;
            Eigen::MatrixXd basis=std::get<2>(vel_tree[tr_id]);
            
            Eigen::MatrixXd bz;
            bz=spline.block<order_num+1,3>(sp_id*(order_num-2),0);

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
                Eigen::RowVector3d vel=order_num*(P[j+1]-P[j]);
                //d=vel_limit*piece_time-vel.norm()/weight;
                d=vel_limit-vel.norm()/(weight*time_weight[sp_id]*piece_time);
                if(d<=0)
                    return INFINITY;
                
                if(d<margin)
                { 
                   energy+=-weight*(d-margin)*(d-margin)*log(d/margin); 
                }
            }
        }


        

        for(unsigned int tr_id=0;tr_id<acc_tree.size();tr_id++)
        {
        
            int sp_id=std::get<0>(acc_tree[tr_id]);
            double weight=std::get<1>(acc_tree[tr_id]).second-std::get<1>(acc_tree[tr_id]).first;
            Eigen::MatrixXd basis=std::get<2>(acc_tree[tr_id]);
            
            Eigen::MatrixXd bz;
            bz=spline.block<order_num+1,3>(sp_id*(order_num-2),0);

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
                Eigen::RowVector3d acc=order_num*(order_num-1)*(P[j+2]-2*P[j+1]+P[j]);
                //d=acc_limit*piece_time*piece_time-acc.norm()/(weight*weight);
                d=acc_limit-acc.norm()/(weight*weight*time_weight[sp_id]*time_weight[sp_id]*piece_time*piece_time);
                if(d<=0)
                    return INFINITY;
                
                if(d<margin)
                { 
                    //std::cout<<"e-acc:"<<acc<<std::endl;
                   energy+=-weight*(d-margin)*(d-margin)*log(d/margin); 
                }
            }
        }

        return energy;  

    }

};

PRJ_END

#endif