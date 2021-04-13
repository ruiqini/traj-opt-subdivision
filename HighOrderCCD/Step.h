#ifndef STEP_H
#define STEP_H

#include "Config/Config.h"

#include "HighOrderCCD/BVH/BVH.h"
#include "HighOrderCCD/CCD/CCD.h"
#include <vector>

PRJ_BEGIN

class Step
{
  public:
    typedef Eigen::MatrixXd Data;
    
    static double position_step(const Data& spline, const Data& direction, 
                                const Eigen::MatrixXd & V,const Eigen::MatrixXi& F,
                                BVH& bvh)
    {
        double step;
        if(step_choose)
        {
            step=1.0;
        }
        else
        {
            step=std::min(1.5*max_step,1.0);
        }

        std::vector<std::pair<unsigned int, unsigned int>> collision_pair,temp_pair;
        
        bvh.InitTrajectory(spline);
        bvh.CheckCollision(collision_pair, offset);
        int num_collision_pair=collision_pair.size();
        
        //clock_t time0 = clock();
        bvh.ccdInitTrajectory(spline, step*direction);
        bvh.CheckCollision(collision_pair, offset);
        std::cout<<"bvh_init: "<<collision_pair.size()<<std::endl;
        double init_step=step;
        if((int)collision_pair.size()>std::max(1000,2*num_collision_pair))
        {
            while(true)
            {
                init_step*=0.2;
                temp_pair=collision_pair;
                bvh.ccdInitTrajectory(spline, init_step*direction);
                bvh.CheckCollision(collision_pair, offset);
                std::cout<<std::endl<<"bvhsize:"<<collision_pair.size()<<std::endl;
                if((int)collision_pair.size()<=std::max(1000,2*num_collision_pair)&&(int)collision_pair.size()>num_collision_pair)
                {
                    break;
                }
                else if((int)collision_pair.size()<=std::max(1000,2*num_collision_pair)&&(int)collision_pair.size()<=num_collision_pair)
                {
                    init_step*=5;
                    collision_pair=temp_pair;
                    break;
                }
            }
        }
        step=init_step;

        //clock_t time1 = clock();
        //std::cout<<std::endl<<"bvhtime:"<<(time1-time0)/(CLOCKS_PER_SEC/1000)<<std::endl;
        bool is_collided=true;
        std::cout<<"bvh: "<<collision_pair.size()<<std::endl;
        
        double temp_step=step;
        double step0=0.0;
        //std::cout<<"start:"<<collision_pair.size()<<std::endl;
        //time0 = clock();

        for(int iter0=0;iter0<2;iter0++)//
        {
            if(collision_pair.size()>0)
            {
                is_collided=true;
            }
            else
            {
                is_collided=false;
            }   
            ///std::vector<std::pair<unsigned int, unsigned int>> temp_pair = collision_pair;

            temp_step=step-step0;

            for(unsigned int it = 0; it < collision_pair.size();it++)
            {
                int tr_id=collision_pair[it].first;
                int ob_id=collision_pair[it].second;

                int sp_id=std::get<0>(subdivide_tree[tr_id]);

                Eigen::MatrixXd basis=std::get<2>(subdivide_tree[tr_id]);
                
                Eigen::MatrixXd bz, bz_d;
                bz=spline.block<order_num+1,3>(sp_id*(order_num-2),0);
                bz_d=direction.block<order_num+1,3>(sp_id*(order_num-2),0);
                
                Eigen::MatrixXd P(order_num+1,3),D(order_num+1,3);
                P=basis*bz;
                D=basis*bz_d;
                
                int f0=F(ob_id,0); int f1=F(ob_id,1); int f2=F(ob_id,2);
                
                Eigen::Matrix3d _position;
                _position<<V.row(f0),V.row(f1),V.row(f2);

                
                is_collided= CCD::KDOPCCD(P,D,_position,offset,step0,step0+temp_step);                                    
                //std::cout<<"is_collided:"<<is_collided<<" \n";
                
                //is_collided=CCD::ExactCCD(P,D,_position);
                while(is_collided)
                {  
                    //is_collided= CCD::KDOPCCD(P,D,_position,offset,step0,step0+temp_step); 
                    is_collided= CCD::GJKCCD(P,D,_position, offset,step0,step0+temp_step);  //cgal
                    if(is_collided)
                    {
                        temp_step*=0.8;
                    }     
                }
            }
            
            if(temp_step/(step0+temp_step)<1e-8)
            {
                step0=step0+temp_step;
                temp_step=step0;
                std::cout<<"temp_step:"<<step0<<std::endl;
                break;
            }
            
            step0=step0+temp_step;
            temp_step=step0;
            std::cout<<"temp_step:"<<step0<<std::endl;
        }
        
        step=temp_step;

        return step;

    }
    
    static double point_position_step(const Data& spline, const Data& direction, 
                                      const Eigen::MatrixXd & V,
                                      BVH& bvh)
    {
        double step;
        if(step_choose)
        {
            step=1.0;
        }
        else
        {
            step=std::min(1.5*max_step,1.0);
        }


        std::vector<std::pair<unsigned int, unsigned int>> collision_pair,temp_pair;
        
        bvh.InitTrajectory(spline);
        bvh.pcCheckCollision(collision_pair, offset);
        int num_collision_pair=collision_pair.size();
        
        //clock_t time0 = clock();
        bvh.ccdInitTrajectory(spline, step*direction);
        bvh.pcCheckCollision(collision_pair, offset);
        std::cout<<"bvh_init: "<<collision_pair.size()<<std::endl;
        double init_step=step;
        if((int)collision_pair.size()>std::max(1000,2*num_collision_pair))
        {
            while(true)
            {
                init_step*=0.2;
                temp_pair=collision_pair;
                bvh.ccdInitTrajectory(spline, init_step*direction);
                bvh.pcCheckCollision(collision_pair, offset);
                std::cout<<std::endl<<"bvhsize:"<<collision_pair.size()<<std::endl;
                if((int)collision_pair.size()<=std::max(1000,2*num_collision_pair)&&(int)collision_pair.size()>num_collision_pair)
                {
                    break;
                }
                else if((int)collision_pair.size()<=std::max(1000,2*num_collision_pair)&&(int)collision_pair.size()<=num_collision_pair)
                {
                    init_step*=5;
                    collision_pair=temp_pair;
                    break;
                }
            }
        }
        step=init_step;

        //clock_t time1 = clock();
        //std::cout<<std::endl<<"bvhtime:"<<(time1-time0)/(CLOCKS_PER_SEC/1000)<<std::endl;
        bool is_collided=true;
        std::cout<<"bvh: "<<collision_pair.size()<<std::endl;
        
        double temp_step=step;
        double step0=0.0;
        //std::cout<<"start:"<<collision_pair.size()<<std::endl;
        //time0 = clock();

        for(int iter0=0;iter0<2;iter0++)//
        {
            if(collision_pair.size()>0)
            {
                is_collided=true;
            }
            else
            {
                is_collided=false;
            }   
            //std::vector<std::pair<unsigned int, unsigned int>> temp_pair = collision_pair;

            temp_step=step-step0;

            for(unsigned int it = 0; it < collision_pair.size();it++)
            {
                int tr_id=collision_pair[it].first;
                int ob_id=collision_pair[it].second;

                int sp_id=std::get<0>(subdivide_tree[tr_id]);

                Eigen::MatrixXd basis=std::get<2>(subdivide_tree[tr_id]);
                
                Eigen::MatrixXd bz, bz_d;
                bz=spline.block<order_num+1,3>(sp_id*(order_num-2),0);
                bz_d=direction.block<order_num+1,3>(sp_id*(order_num-2),0);
                
                Eigen::MatrixXd P(order_num+1,3),D(order_num+1,3);
                P=basis*bz;
                D=basis*bz_d;
                
                Eigen::RowVector3d _position;
                _position=V.row(ob_id);

                
                is_collided= CCD::KDOPCCD(P,D,_position,offset,step0,step0+temp_step);                                    
                
                while(is_collided)
                {  
                    //is_collided= CCD::KDOPCCD(P,D,_position,offset,step0,step0+temp_step); 
                    is_collided= CCD::GJKCCD(P,D,_position, offset,step0,step0+temp_step);  //cgal
                    //std::cout<<"is_collided:"<<it<<" "<<is_collided<<" \n";
                    if(is_collided)
                    {
                        temp_step*=0.8;
                    }     
                }
            }
           
            if(temp_step/(step0+temp_step)<1e-8)
            {
                step0=step0+temp_step;
                temp_step=step0;
                std::cout<<"temp_step:"<<step0<<std::endl;
                break;
            }
            
            step0=step0+temp_step;
            temp_step=step0;
            std::cout<<"temp_step:"<<step0<<std::endl;
        }
        
        step=temp_step;

        return step;

    }
    
};

PRJ_END

#endif