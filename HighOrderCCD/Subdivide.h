#ifndef SUBDIVIDE_H
#define SUBDIVIDE_H

#include "Config/Config.h"

#include "HighOrderCCD/BVH/BVH.h"
#include "HighOrderCCD/CCD/CCD.h"
#include <vector>

PRJ_BEGIN

class Subdivide
{
  public:
    typedef Eigen::MatrixXd Data;
    typedef std::vector< std::tuple< int, std::pair<double,double>, Eigen::MatrixXd > > Tree;
    typedef std::tuple< int, std::pair<double,double>, Eigen::MatrixXd >  Node;

    static void update_tree(const Data& spline, 
                         const Eigen::MatrixXd & V,const Eigen::MatrixXi& F,
                         BVH& bvh)
    {
      bvh.InitTrajectory(spline);

      std::vector<std::pair<unsigned int, unsigned int>> collision_pair;

      bvh.CheckCollision(collision_pair, offset+margin);
      
      std::vector<std::pair<int,std::vector<double>>> list_range;
      list_range.reserve(100);

      for(auto it = collision_pair.begin(); it != collision_pair.end();++it)
      {
        int tr_id=(*it).first;
        int ob_id=(*it).second;
        
        int sp_id=std::get<0>(subdivide_tree[tr_id]);
      
        int f0=F(ob_id,0); int f1=F(ob_id,1); int f2=F(ob_id,2);
        
        Eigen::Matrix3d _position;

        _position<<V.row(f0),V.row(f1),V.row(f2);

        Eigen::MatrixXd bz;
        bz=spline.block<order_num+1,3>(sp_id*(order_num-2),0);
        Eigen::MatrixXd x;
        /*
        if(sp_id==0)
        {
          x=M_head*bz;
        }
        else if(sp_id==piece_num-1)
        {
          x=M_tail*bz;
        }
        else
        {
          x=M_convert*bz;
        }
        */
        x=convert_list[sp_id]*bz;
        
        std::pair<double, double> range=std::get<1>(subdivide_tree[tr_id]);
        std::vector<double> range_list; range_list.resize(2);

        range_list[0]=range.first;
        range_list[1]=range.second;

        std::stack<std::pair<double, double>> ss;
        ss.push(range);

        while(!ss.empty()) 
        {      
          range=ss.top();
          ss.pop();

          double head=range.first;
          double tail=range.second;
          Eigen::MatrixXd basis;
          Blossom<order_num>::coefficient(basis, head, tail);
          Eigen::MatrixXd temp=basis*x;
          double diam2=(temp.row(0)-temp.row(order_num)).squaredNorm();
          if(CCD::KDOPDCD(temp, _position, offset+margin))
          {
            if(diam2>=epsilon2)
            {
              double mid = 0.5*(head+tail);
             ss.push(std::make_pair(head, mid));
             ss.push(std::make_pair(mid, tail));

             std::vector<double>::iterator it0,it1; 
    
              it0 = std::find (range_list.begin(), range_list.end(), head); 
              it1 = std::find (range_list.begin(), range_list.end(), mid); 
              if (it0 != range_list.end() && it1 == range_list.end()) 
              {
                range_list.insert(it0+1, mid); 
              }
            }
             
          }
        }
        auto it2 = std::find_if( list_range.begin(), list_range.end(),
                      [&](const std::pair<int, std::vector<double>>& element){ return element.first == tr_id;} );
        if(it2==list_range.end())
        {
          list_range.push_back(std::make_pair(tr_id,range_list));
        }
        //std::cout<<tr_id<<"\n";
        /*
        for(int i=0;i<range_list.size()-1;i++)
        {
          std::cout<<range_list[i]<<" "<<range_list[i+1]<<"\n";
        }
        */
      }
      Tree temp_tree=subdivide_tree;

      
      std::cout<<"list_range:"<<list_range.size()<<std::endl;
      std::sort(list_range.begin(),list_range.end());
      std::reverse(list_range.begin(),list_range.end());

      for(unsigned int i=0;i<list_range.size();i++)
      {
        //std::cout<<list_range[i].first<<std::endl;
        std::vector<double> range_list=list_range[i].second;
        int sp_id=std::get<0>(subdivide_tree[list_range[i].first]);
        int tree_size=range_list.size()-1;
        Tree new_tree; new_tree.reserve(tree_size);
        for (int j = 0; j < tree_size; j++) 
        {
          double head=range_list[j];
          double tail=range_list[j+1];
          Eigen::MatrixXd basis;

          Blossom<order_num>::coefficient(basis, head, tail);

          std::pair<double,double> range(head,tail);
          Node node;
          /*
          if(sp_id==0)
          {
            node=std::make_tuple(sp_id,range,basis*M_head);
          }
          else if(sp_id==piece_num-1)
          {
            node=std::make_tuple(sp_id,range,basis*M_tail);
          }
          else
          {
            node=std::make_tuple(sp_id,range,basis*M_convert);
          }
          */
          node=std::make_tuple(sp_id,range,basis*convert_list[sp_id]);
          //Node node=std::make_tuple(sp_id,range,basis);
          new_tree.push_back(node);
        }
        temp_tree.insert(temp_tree.begin()+list_range[i].first+1,new_tree.begin(),new_tree.end());
        temp_tree.erase(temp_tree.begin()+list_range[i].first);
      }

      subdivide_tree=temp_tree;
    
  }

  static void update_tree_point(const Data& spline, 
                         const Eigen::MatrixXd & V,
                         BVH& bvh)
    {
      bvh.InitTrajectory(spline);

      std::vector<std::pair<unsigned int, unsigned int>> collision_pair;

      bvh.pcCheckCollision(collision_pair, offset+margin);
      
      std::vector<std::pair<int,std::vector<double>>> list_range;
      list_range.reserve(100);

      for(auto it = collision_pair.begin(); it != collision_pair.end();++it)
      {
        int tr_id=(*it).first;
        int ob_id=(*it).second;
        
        int sp_id=std::get<0>(subdivide_tree[tr_id]);
      
        
        Eigen::RowVector3d _position=V.row(ob_id);

        Eigen::MatrixXd bz;
        bz=spline.block<order_num+1,3>(sp_id*(order_num-2),0);
        Eigen::MatrixXd x;
        /*
        if(sp_id==0)
        {
          x=M_head*bz;
        }
        else if(sp_id==piece_num-1)
        {
          x=M_tail*bz;
        }
        else
        {
          x=M_convert*bz;
        }
        */
        x=convert_list[sp_id]*bz;
        std::pair<double, double> range=std::get<1>(subdivide_tree[tr_id]);
        std::vector<double> range_list; range_list.resize(2);

        range_list[0]=range.first;
        range_list[1]=range.second;

        std::stack<std::pair<double, double>> ss;
        ss.push(range);

        while(!ss.empty()) 
        {      
          range=ss.top();
          ss.pop();

          double head=range.first;
          double tail=range.second;
          Eigen::MatrixXd basis;
          Blossom<order_num>::coefficient(basis, head, tail);
          Eigen::MatrixXd temp=basis*x;
          double diam2=(temp.row(0)-temp.row(order_num)).squaredNorm();
          if(CCD::KDOPDCD_point(temp, _position, offset+margin))
          {
            if(diam2>=epsilon2)
            {
              double mid = 0.5*(head+tail);
             ss.push(std::make_pair(head, mid));
             ss.push(std::make_pair(mid, tail));

             std::vector<double>::iterator it0,it1; 
    
              it0 = std::find (range_list.begin(), range_list.end(), head); 
              it1 = std::find (range_list.begin(), range_list.end(), mid); 
              if (it0 != range_list.end() && it1 == range_list.end()) 
              {
                range_list.insert(it0+1, mid); 
              }
            }
             
          }
        }
        auto it2 = std::find_if( list_range.begin(), list_range.end(),
                      [&](const std::pair<int, std::vector<double>>& element){ return element.first == tr_id;} );
        if(it2==list_range.end())
        {
          list_range.push_back(std::make_pair(tr_id,range_list));
        }
        //std::cout<<tr_id<<"\n";
        /*
        for(int i=0;i<range_list.size()-1;i++)
        {
          std::cout<<range_list[i]<<" "<<range_list[i+1]<<"\n";
        }
        */
      }
      Tree temp_tree=subdivide_tree;

      
      std::cout<<"list_range:"<<list_range.size()<<std::endl;
      std::sort(list_range.begin(),list_range.end());
      std::reverse(list_range.begin(),list_range.end());

      for(unsigned int i=0;i<list_range.size();i++)
      {
        //std::cout<<list_range[i].first<<std::endl;
        std::vector<double> range_list=list_range[i].second;
        int sp_id=std::get<0>(subdivide_tree[list_range[i].first]);
        int tree_size=range_list.size()-1;
        Tree new_tree; new_tree.reserve(tree_size);
        for (int j = 0; j < tree_size; j++) 
        {
          double head=range_list[j];
          double tail=range_list[j+1];
          Eigen::MatrixXd basis;

          Blossom<order_num>::coefficient(basis, head, tail);

          std::pair<double,double> range(head,tail);
          Node node;
          /*
          if(sp_id==0)
          {
            node=std::make_tuple(sp_id,range,basis*M_head);
          }
          else if(sp_id==piece_num-1)
          {
            node=std::make_tuple(sp_id,range,basis*M_tail);
          }
          else
          {
            node=std::make_tuple(sp_id,range,basis*M_convert);
          }
          */
          node=std::make_tuple(sp_id,range,basis*convert_list[sp_id]);
          //Node node=std::make_tuple(sp_id,range,basis);
          new_tree.push_back(node);
        }
        temp_tree.insert(temp_tree.begin()+list_range[i].first+1,new_tree.begin(),new_tree.end());
        temp_tree.erase(temp_tree.begin()+list_range[i].first);
      }

      subdivide_tree=temp_tree;
    
  }

  static void update_vel(const Data& spline, double piece_time,
                         BVH& bvh)
  {
      std::vector<std::pair<int,std::vector<double>>> list_range;
      list_range.reserve(100);

      for(unsigned int tr_id=0;tr_id<vel_tree.size();tr_id++)
      {

          int sp_id=std::get<0>(vel_tree[tr_id]);
          double weight=std::get<1>(vel_tree[tr_id]).second-std::get<1>(vel_tree[tr_id]).first;
          Eigen::MatrixXd basis=std::get<2>(vel_tree[tr_id]);

          //std::cout<<"sp_id:"<<sp_id<<std::endl;
          Eigen::MatrixXd bz;
          int init=sp_id*(order_num-2);
          bz=spline.block<order_num+1,3>(init,0);

          std::vector<Eigen::RowVector3d> P_(order_num+1);
          for(int j=0;j<=order_num;j++)
          {   
              P_[j].setZero();
              for(int j0=0;j0<=order_num;j0++)
              {
                  P_[j]+=basis(j,j0)*bz.row(j0);
              }
          }
          
          double d_;
          bool reach_bound=false;        
          for(int j=0;j<order_num;j++)
          {
              Eigen::RowVector3d vel=order_num*(P_[j+1]-P_[j]);
              d_=vel_limit-vel.norm()/(weight*time_weight[sp_id]*piece_time);

              if(d_<margin)
              {
                reach_bound=true;
              }
          }

          if(reach_bound)
          {
            std::pair<double, double> range=std::get<1>(vel_tree[tr_id]);
            std::vector<double> range_list; range_list.resize(2);

            range_list[0]=range.first;
            range_list[1]=range.second;

            std::stack<std::pair<double, double>> ss;
            ss.push(range);

            while(!ss.empty()) 
            {      
              range=ss.top();
              ss.pop();

              double head=range.first;
              double tail=range.second;

              weight=tail-head;
              
              Eigen::MatrixXd basis;
              Blossom<order_num>::coefficient(basis, head, tail);
              Eigen::MatrixXd temp=basis*bz;
              double diam2=(temp.row(0)-temp.row(order_num)).squaredNorm();

              reach_bound=false;        
              for(int j=0;j<order_num;j++)
              {
                  Eigen::RowVector3d vel=order_num*(temp.row(j+1)-temp.row(j));
                  d_=vel_limit-vel.norm()/(weight*time_weight[sp_id]*piece_time);

                  if(d_<margin)
                  {
                    reach_bound=true;
                  }
              }

              if(reach_bound)
              {
                if(diam2>=epsilon2)
                {
                  double mid = 0.5*(head+tail);
                  ss.push(std::make_pair(head, mid));
                  ss.push(std::make_pair(mid, tail));

                  std::vector<double>::iterator it0,it1; 
        
                  it0 = std::find (range_list.begin(), range_list.end(), head); 
                  it1 = std::find (range_list.begin(), range_list.end(), mid); 
                  if (it0 != range_list.end() && it1 == range_list.end()) 
                  {
                    range_list.insert(it0+1, mid); 
                  }
                }
                
              }
            }
            auto it2 = std::find_if( list_range.begin(), list_range.end(),
                          [&](const std::pair<int, std::vector<double>>& element){ return element.first == (int)tr_id;} );
            if(it2==list_range.end())
            {
              list_range.push_back(std::make_pair(tr_id,range_list));
            }
          }
      }

      Tree temp_tree=vel_tree;

      
      std::cout<<"vel_range:"<<list_range.size()<<std::endl;
      std::sort(list_range.begin(),list_range.end());
      std::reverse(list_range.begin(),list_range.end());

      for(unsigned int i=0;i<list_range.size();i++)
      {
        //std::cout<<list_range[i].first<<std::endl;
        std::vector<double> range_list=list_range[i].second;
        int sp_id=std::get<0>(vel_tree[list_range[i].first]);
        int tree_size=range_list.size()-1;
        Tree new_tree; new_tree.reserve(tree_size);
        for (int j = 0; j < tree_size; j++) 
        {
          double head=range_list[j];
          double tail=range_list[j+1];

          Eigen::MatrixXd basis;

          Blossom<order_num>::coefficient(basis, head, tail);

          std::pair<double,double> range(head,tail);
          Node node;
          /*
          if(sp_id==0)
          {
            node=std::make_tuple(sp_id,range,basis*M_head);
          }
          else if(sp_id==piece_num-1)
          {
            node=std::make_tuple(sp_id,range,basis*M_tail);
          }
          else
          {
            node=std::make_tuple(sp_id,range,basis*M_convert);
          }
          */
          node=std::make_tuple(sp_id,range,basis*convert_list[sp_id]);
          //Node node=std::make_tuple(sp_id,range,basis);
          new_tree.push_back(node);
        }
        temp_tree.insert(temp_tree.begin()+list_range[i].first+1,new_tree.begin(),new_tree.end());
        temp_tree.erase(temp_tree.begin()+list_range[i].first);
      }

      vel_tree=temp_tree;
  }

  static void update_acc(const Data& spline, const double& piece_time,
                         BVH& bvh)
  {
      std::vector<std::pair<int,std::vector<double>>> list_range;
      list_range.reserve(100);

      for(unsigned int tr_id=0;tr_id<acc_tree.size();tr_id++)
      {

          int sp_id=std::get<0>(acc_tree[tr_id]);
          double weight=std::get<1>(acc_tree[tr_id]).second-std::get<1>(acc_tree[tr_id]).first;
          Eigen::MatrixXd basis=std::get<2>(acc_tree[tr_id]);

          //std::cout<<"sp_id:"<<sp_id<<std::endl;
          Eigen::MatrixXd bz;
          int init=sp_id*(order_num-2);
          bz=spline.block<order_num+1,3>(init,0);

          std::vector<Eigen::RowVector3d> P_(order_num+1);
          for(int j=0;j<=order_num;j++)
          {   
              P_[j].setZero();
              for(int j0=0;j0<=order_num;j0++)
              {
                  P_[j]+=basis(j,j0)*bz.row(j0);
              }
          }
          
          double d_;
          bool reach_bound=false;        
          for(int j=0;j<order_num-1;j++)
          {
              Eigen::RowVector3d acc=order_num*(order_num-1)*(P_[j+2]-2*P_[j+1]+P_[j]);
                
              d_=acc_limit-acc.norm()/(weight*weight*time_weight[sp_id]*time_weight[sp_id]*piece_time*piece_time);

              if(d_<margin)
              {
                reach_bound=true;
              }
          }

          if(reach_bound)
          {
            std::pair<double, double> range=std::get<1>(acc_tree[tr_id]);
            std::vector<double> range_list; range_list.resize(2);

            range_list[0]=range.first;
            range_list[1]=range.second;

            std::stack<std::pair<double, double>> ss;
            ss.push(range);

            while(!ss.empty()) 
            {      
              range=ss.top();
              ss.pop();

              double head=range.first;
              double tail=range.second;
              weight=tail-head;

              Eigen::MatrixXd basis;
              Blossom<order_num>::coefficient(basis, head, tail);
              Eigen::MatrixXd temp=basis*bz;
              double diam2=(temp.row(0)-temp.row(order_num)).squaredNorm();

              reach_bound=false;        
              for(int j=0;j<order_num-1;j++)
              {
                  Eigen::RowVector3d acc=order_num*(order_num-1)*(temp.row(j+2)-2*temp.row(j+1)+temp.row(j));
                
                  d_=acc_limit-acc.norm()/(weight*weight*time_weight[sp_id]*time_weight[sp_id]*piece_time*piece_time);

                  if(d_<margin)
                  {
                    reach_bound=true;
                  }
              }

              if(reach_bound)
              {
                if(diam2>=epsilon2)
                {
                  double mid = 0.5*(head+tail);
                  ss.push(std::make_pair(head, mid));
                  ss.push(std::make_pair(mid, tail));

                  std::vector<double>::iterator it0,it1; 
        
                  it0 = std::find (range_list.begin(), range_list.end(), head); 
                  it1 = std::find (range_list.begin(), range_list.end(), mid); 
                  if (it0 != range_list.end() && it1 == range_list.end()) 
                  {
                    range_list.insert(it0+1, mid); 
                  }
                }
                
              }
            }
            auto it2 = std::find_if( list_range.begin(), list_range.end(),
                          [&](const std::pair<int, std::vector<double>>& element){ return element.first == (int)tr_id;} );
            if(it2==list_range.end())
            {
              list_range.push_back(std::make_pair(tr_id,range_list));
            }
          }
      }

      Tree temp_tree=acc_tree;

      
      std::cout<<"acc_range:"<<list_range.size()<<std::endl;
      std::sort(list_range.begin(),list_range.end());
      std::reverse(list_range.begin(),list_range.end());

      for(unsigned int i=0;i<list_range.size();i++)
      {
        //std::cout<<list_range[i].first<<std::endl;
        std::vector<double> range_list=list_range[i].second;
        int sp_id=std::get<0>(acc_tree[list_range[i].first]);
        int tree_size=range_list.size()-1;
        Tree new_tree; new_tree.reserve(tree_size);
        for (int j = 0; j < tree_size; j++) 
        {
          double head=range_list[j];
          double tail=range_list[j+1];
          Eigen::MatrixXd basis;

          Blossom<order_num>::coefficient(basis, head, tail);

          std::pair<double,double> range(head,tail);
          Node node;
          /*
          if(sp_id==0)
          {
            node=std::make_tuple(sp_id,range,basis*M_head);
          }
          else if(sp_id==piece_num-1)
          {
            node=std::make_tuple(sp_id,range,basis*M_tail);
          }
          else
          {
            node=std::make_tuple(sp_id,range,basis*M_convert);
          }
          */
          node=std::make_tuple(sp_id,range,basis*convert_list[sp_id]);
          //Node node=std::make_tuple(sp_id,range,basis);
          new_tree.push_back(node);
        }
        temp_tree.insert(temp_tree.begin()+list_range[i].first+1,new_tree.begin(),new_tree.end());
        temp_tree.erase(temp_tree.begin()+list_range[i].first);
      }

      acc_tree=temp_tree;
  }

};

PRJ_END

#endif