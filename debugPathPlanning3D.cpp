
#include "HighOrderCCD/Config/Config.h"

#include "HighOrderCCD/Optimization/Optimization3D_time.h"

#include "HighOrderCCD/Optimization/Optimization3D_point.h"

#include <vector>
#include <ctime>
#include <igl/opengl/glfw/Viewer.h>

#include <igl/read_triangle_mesh.h>
#include <igl/write_triangle_mesh.h>

#include "nlohmann/json.hpp"    // https://github.com/nlohmann/json/tree/develop/single_include/nlohmann/json.hpp
#include <fstream>              // ifstream, ofstream
#include <unistd.h>

USE_PRJ_NAMESPACE

using json = nlohmann::json;

Eigen::Vector3d getPosFromBezier(const Eigen::MatrixXd & polyCoeff, double t_now, int seg_now , int _traj_order) 
{
    Eigen::Vector3d ret = Eigen::VectorXd::Zero(3);
    Eigen::VectorXd ctrl_now = polyCoeff.row(seg_now);
    int ctrl_num1D = polyCoeff.cols() / 3;

    for(int i = 0; i < 3; i++)
        for(int j = 0; j < ctrl_num1D; j++)
            ret(i) += combination[_traj_order][j] * ctrl_now(i * ctrl_num1D + j) * pow(t_now, j) * pow((1 - t_now), (_traj_order - j) ); 

    return ret;  
}

void log_data(std::string meshfile, Eigen::MatrixXd spline, double piece_time)
{
              Eigen::VectorXd my_time(piece_num);
              Eigen::MatrixXd coeff(piece_num, 3*(order_num+1));
              for(int i=0;i<piece_num;i++)
              {
                my_time(i)=time_weight[i]*piece_time;
                Eigen::MatrixXd bz;
                bz=spline.block<order_num+1,3>(i*(order_num-2),0);
                
                bz=convert_list[i]*bz;
                coeff.block<1,order_num+1>(i,0)=bz.col(0).transpose();
                coeff.block<1,order_num+1>(i,order_num+1)=bz.col(1).transpose();
                coeff.block<1,order_num+1>(i,2*(order_num+1))=bz.col(2).transpose();
              }
              std::ofstream file(meshfile+"_ccd_test.txt");
              if (file.is_open())
              {
                file << "coeff:\n" << coeff << '\n';
                file << "time:" << '\n' <<  my_time << '\n';
              }
              file.close();

              std::vector<Eigen::Vector3d> ccd_traj;
              /*
              int segment_num  = coeff.rows();
              for(int i = 0; i < segment_num; i++ ){
                      Eigen::Vector3d state0 = getPosFromBezier( coeff, 0, i ,order_num);
                      Eigen::Vector3d cur0;
                      cur0(0) =   state0(0);
                      cur0(1) =   state0(1);
                      cur0(2) =   state0(2);
                  for (double t = 0.0; t < 1.0; t += 0.05 / my_time(i)){
                      Eigen::Vector3d state = getPosFromBezier( coeff, t, i ,order_num);
                      Eigen::Vector3d cur;
                      cur(0) =   state(0);
                      cur(1) =   state(1);
                      cur(2) =   state(2);
                      ccd_traj.push_back(cur);
                  }
              }
              */
             for (double t = 0.0; t < piece_num; t += 0.05 / piece_time){
                      int i = floor(t);
                      double cur_t=t-i;
                      Eigen::Vector3d state = getPosFromBezier( coeff, cur_t, i ,order_num);
                      Eigen::Vector3d cur;
                      cur(0) =   state(0);
                      cur(1) =   state(1);
                      cur(2) =   state(2);
                      ccd_traj.push_back(cur);
                  }
              std::cout<<"ccd time:"<<my_time.sum()<<std::endl;
              double  len_ccd=0;
              for(int i=0;i<(int)ccd_traj.size()-1;i++)
              {
                  len_ccd+=(ccd_traj[i+1]-ccd_traj[i]).norm();;
              }
              std::cout<<"ccd len:"<<len_ccd<<std::endl;
              

              Eigen::MatrixXd v_ccd(2*ccd_traj.size()-1,3);
              Eigen::MatrixXi f_ccd(ccd_traj.size()-1,3);
              for(int i=0;i<(int)ccd_traj.size();i++)
              {
                  v_ccd.row(i)=ccd_traj[i].transpose();
              }
              for(int i=0;i<(int)ccd_traj.size()-1;i++)
              {
                  v_ccd.row(i+ccd_traj.size())=0.5*(ccd_traj[i].transpose()+ccd_traj[i+1].transpose());
              }
              for(int i=0;i<(int)ccd_traj.size()-1;i++)
              {
                  f_ccd(i,0)=i; f_ccd(i,1)=i+1; f_ccd(i,2)=i+ccd_traj.size();
              }

              igl::write_triangle_mesh("ccd_traj_"+meshfile+".obj",v_ccd,f_ccd);
}

int main(int argc, char *argv[])
{
  if (argc < 2)
	{
		std::cerr << "Syntax: " << argv[0] << " <config file>" << std::endl;
		return -1;
  }

  #if 0
  GUI drawer;
  drawer.init();
  drawer.invoke();
  #else

  typedef Eigen::MatrixXd Data;

  
  igl::opengl::glfw::Viewer viewer;
  
  viewer.core().background_color<< 1.0f, 1.0f, 1.0f, 1.0f;
  viewer.core().is_animating = true;
  viewer.core().camera_zoom = 2.0f;
  viewer.data().line_width = 1.0f;
  viewer.data().point_size = 3.0f;

  const std::string config_file = argv[1];
  std::ifstream fin("Config File/"+config_file);   
  json j = json::parse(fin);
	fin.close();


  lambda=j["lambda"].get<double>();

  epsilon=j["epsilon"].get<double>();
  epsilon2 = epsilon*epsilon;
  margin=j["margin"].get<double>();//d hat
  margin2 = margin*margin;
  automove = j["auto"].get<int>();
  step_choose = j["step_choose"].get<int>();

  int init_spline=j["init"].get<int>();
  
  offset = j["offset"].get<double>();

  res = j["res"].get<int>();

  int if_exit=j["exit"].get<int>();

  int if_init_ob=j["init_ob"].get<int>();
  int if_init_wp=j["init_wp"].get<int>();

  int if_init_pc=j["init_pc"].get<int>();

  //lambda=1.0/margin2;

  double whole_time=0;
  
  adaptive_change=false;
  optimize_time=false;

  int dim = axis.size();
  for(int k=0;k<dim;k++)
  {
    axis[k].normalize();
  }

  dim = aabb_axis.size();
  for(int k=0;k<dim;k++)
  {
    aabb_axis[k].normalize();
  }
  //lambda/=epsilon;

  Eigen::MatrixXd V, BV;
  Eigen::MatrixXi F, BF;
  
  const std::string mesh_file = argv[2];
  igl::read_triangle_mesh(mesh_file,V,F);//32770 cylinder
  
  vel_limit=j["vel_limit"].get<double>();;//2
  acc_limit=j["acc_limit"].get<double>();;//2
    
  //energy_file.open ("energy/" +mesh_file + "_energy_file.txt");
  result_file.open ("result/" +mesh_file + "_result_file.txt");
  //gn_file.open ("gn/" +mesh_file + "_gn_file.txt");
  //step_file.open ("step/" +mesh_file + "_step_file.txt");
  //eigen_file.open("eigen/"+mesh_file + "_eigen_file.txt");
  //plot_file.open("plot/"+mesh_file + "_plot_file.txt");


  Data spline;
  if(init_spline==0)
  {
    double r=4.0, h=1.0;
    piece_num=5;
    trajectory_num = (order_num+1)+(piece_num-1)*(order_num+1-3);
    spline.resize(trajectory_num,3);
    std::vector<Eigen::RowVector3d> way_points;
    
    for(int i=0;i<=piece_num;i++)
    {
      Eigen::Vector3d p0(r*cos(i*M_PI/3.0),r*sin(i*M_PI/3.0),-(piece_num-i)/double(piece_num)*h+i/double(piece_num)*h);
      way_points.push_back(p0);
    }
    
    time_weight.resize(piece_num);
    whole_weight=0;
    for(int i=0;i<piece_num;i++)
    {
      time_weight[i]=1;
      //time_weight[i]=(way_points[i+1]-way_points[i]).norm()/vel_limit;
      whole_weight+=time_weight[i];
    }

    spline.row(0)=way_points[0];
    for(int i=0;i<piece_num;i++)
    {
      for(int j=0;j<=order_num-2;j++)
      {
        spline.row(j+i*(order_num-2)+1)=double(order_num-2-j)/(order_num-2)*way_points[i]+(double)j/(order_num-2)*way_points[i+1];
      }
    }
    spline.row(trajectory_num-1)=way_points[piece_num];
  }
  else if(init_spline==1)
  {
    std::string line;
    std::ifstream myfile ("init/" +mesh_file + "_init_file.txt");
    //V*=0.5;
    Eigen::Vector3d p0;
    std::vector<Eigen::RowVector3d> way_points;
    if (myfile.is_open())
    {
      int i=0;
      while ( getline (myfile,line) )
      {
        std::istringstream iss(line);
        
        iss>>p0(0)>>p0(1)>>p0(2);
        way_points.push_back(p0);
        
        i++;
      }
      myfile.close();
      piece_num=i-1;
    }
    time_weight.resize(piece_num);
    whole_weight=0;
    for(int i=0;i<piece_num;i++)
    {
      time_weight[i]=1;
      //time_weight[i]=(way_points[i+1]-way_points[i]).norm()/vel_limit;
      whole_weight+=time_weight[i];
    }
    
    trajectory_num = (order_num+1)+(piece_num-1)*(order_num+1-3);
    spline.resize(trajectory_num,3);

    spline.row(0)=way_points[0];
    for(int i=0;i<piece_num;i++)
    { 
      
      Eigen::RowVector3d head=0.9*way_points[i]+0.1*way_points[i+1];
      Eigen::RowVector3d tail=0.9*way_points[i+1]+0.1*way_points[i];
      spline.row(i*(order_num-2)+1)=way_points[i];
      for(int j=1;j<order_num-2;j++)
      {
        spline.row(j+i*(order_num-2)+1)=double(order_num-3-j)/(order_num-4)*head+(double)(j-1)/(order_num-4)*tail;
      }
      spline.row((i+1)*(order_num-2)+1)=way_points[i+1];
      /*
      for(int j=0;j<=order_num-2;j++)
      {
        spline.row(j+i*(order_num-2)+1)=double(order_num-2-j)/(order_num-2)*way_points[i]+(double)(j)/(order_num-2)*way_points[i+1];
      }
      */
      
    }
    spline.row(trajectory_num-1)=way_points[piece_num];

  }
  else if(init_spline==2)
  {
    std::string line;
    std::ifstream myfile ("init/b.txt");
    trajectory_num = (order_num+1)+(piece_num-1)*(order_num+1-3);
    spline.resize(trajectory_num,3);
    if (myfile.is_open())
    {
      int i=0;
      while ( getline (myfile,line) )
      {
        std::istringstream iss(line);
        if(i<trajectory_num)
        {
          iss>>spline(i,0)>>spline(i,1)>>spline(i,2);
        }
        i++;
      }
      myfile.close();
    }
    
  }
  if(if_init_wp)
  {
    std::string line;
    std::ifstream myfile ("init/" +mesh_file + "_fix.txt");
    if (myfile.is_open())
    {
      int i=0;
      while ( getline (myfile,line) )
      {
        std::istringstream iss(line);
        int id;
        iss>>id;
        fixed_points.push_back(id);
        i++;
      }
      myfile.close();
    }
  }
  else
  {
    fixed_points.clear();
  }
  

  //
  //std::cout<<spline<<"\n";
  //init_file.open ("init/" +mesh_file + "_init_file.txt");
  //spline*=10;
  //V*=10;

  max_step=1.0;



  combination = Combination<40>::value();

  //Eigen::MatrixXd V,O;
  
  //std::vector<Eigen::MatrixXd> p0,p1,d0;
  //V*=0.000000001;
  int num=1000000, turns=0;

  gnorm=1;
  
  ks=1e-3;//1e-3
  kt=1;
  double piece_time=20;//20

  Conversion<order_num>::convert_matrix();
  
  std::cout<<convert_list[0]<<std::endl;
  //std::cout<<M_head<<std::endl;
  //std::cout<<M_tail<<std::endl;
  //M_convert/=double(Factorial<order_num>::value());

  
  spline.row(1)=spline.row(0);
  spline.row(trajectory_num-2)=spline.row(trajectory_num-1);

  M_dynamic=Dynamic<order_num, der_num>::dynamic_matrix();//Dynamic<order_num,1>::dynamic_matrix()+
  
  subdivide_tree.resize(piece_num*res);
  Eigen::MatrixXd basis;
  
  for(int k=0;k<res;k++)
  {
    double a=k/double(res),b=(k+1)/double(res);

    Blossom<order_num>::coefficient(basis, a, b);

    for(int i=0;i<piece_num;i++)
    { 
      std::pair<double,double> range(a,b);
      subdivide_tree[i*res+k]=std::make_tuple(i,range,basis*convert_list[i]);
      
    }
  }

  vel_tree=subdivide_tree;
  acc_tree=subdivide_tree;
 
  for(int i=0;i<(int)fixed_points.size();i++)
  {
        int id=fixed_points[i]*(order_num-2)+1;
        std::cout<<0.25*(spline.row(id-1)+spline.row(id+1))+0.5*spline.row(id)<<"\n";
  }
  
  std::cout<<"before bvh init\n";
  BVH bvh;
  clock_t time1 = clock();
  if(if_init_ob)
  {
    if(if_init_pc)
    {
      bvh.InitPointcloud(V);
    }
    else
    {
      bvh.InitObstacle(V,F);
    }
  }
    

  clock_t time2 = clock();

  std::cout<<"time_obstacle:"<<(time2-time1)/(CLOCKS_PER_SEC/1000)<<std::endl<<std::endl;
  //std::cout<<F_<<std::endl;
  
  Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic> IR(640,800);
  Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic> IG(640,800);
  Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic> IB(640,800);
  Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic> IA(640,800);

  //bool draw_init=true;

  bool is_write=true;

  const auto &key_down = [&](igl::opengl::glfw::Viewer &,unsigned char key,int mod)->bool
  {
    switch(key)
    {
      case 'R':
      {
        iter=0;
        turns=0;
        break;
      }
      case 'C':
      {
        break;
      }
      case ' ':
      {
        //viewer.core().draw_buffer(viewer.data(),false,IR,IG,IB,IA);
        //igl::png::writePNG(IR,IG,IB,IA, "../png/" + std::to_string(turns) + ".png");
        turns=iter+1;
        break;
      }
      case '.':
      {
        automove=false;
        break;
      }
      case ',':
      {
        automove=true;
        break;
      }
      default:
        return false;
    }

    return true;
  };
  
  viewer.data().set_mesh(V,F);
  
  const auto &pre_draw = [&](igl::opengl::glfw::Viewer & )->bool
  {  
     if(iter<num) 
     {
        viewer.data().clear_edges();
        viewer.data().line_width = 1.0f;
 
                
        //viewer.data().line_width = 5.0f;
        for(unsigned int k=0;k<subdivide_tree.size();k++)
        {
          int sp_id=std::get<0>(subdivide_tree[k]);
          Eigen::MatrixXd basis=std::get<2>(subdivide_tree[k]);
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
          if(k%2==1)
          {
            for(int j=0;j<=order_num;j++)
            {
              viewer.data().add_edges(P[j], P[(j+1)%(order_num+1)], Eigen::RowVector3d(0.2,0.8,0.8));
            }
          }
          else
          {
            for(int j=0;j<=order_num;j++)
            {
              viewer.data().add_edges(P[j], P[(j+1)%(order_num+1)], Eigen::RowVector3d(0.8,0.2,0.8));
            }
          }
        }
        
        if(gnorm<1e-3 && is_write)//1e-3
        {
          if(optimize_time)
          {
            if(adaptive_change)
            {
              is_write = false;
              result_file<<iter<<std::endl;
              result_file<<whole_time<<std::endl;
              result_file<<gnorm<<std::endl;
              result_file<<V.rows()<<" "<<F.rows()<<std::endl;
              result_file<<spline<<std::endl;
              result_file<<piece_time<<std::endl;

              if(if_exit)
                exit(0);
              else
                automove=false;
              
              log_data(mesh_file, spline, piece_time);
            
            }

            adaptive_change=true;

          }
          else
          {
            result_file<<iter<<std::endl;
            result_file<<whole_time<<std::endl;
            result_file<<gnorm<<std::endl;
            result_file<<V.rows()<<" "<<F.rows()<<std::endl;
            result_file<<spline<<std::endl;

          
            optimize_time=true;
          }
          
          
        }
        
        if(iter<turns||automove)
        {

          clock_t time0 = clock();
          std::cout<<iter<<std::endl;
              

          
          if(if_init_pc)
          {
            Optimization3D_point::optimization(spline, piece_time,V, bvh);
          }
          else
          {
            Optimization3D_time::optimization(spline,piece_time, V, F, bvh);
          }
          
          
          
          clock_t time1 = clock();
          whole_time+=(time1-time0)/(CLOCKS_PER_SEC/1000);
          std::cout<<"time:"<<(time1-time0)/(CLOCKS_PER_SEC/1000)<<std::endl<<std::endl;
          //std::cout<<p0.size()<<std::endl;
          iter++;
        }
     }  
     
      return false;
  };
  
  viewer.callback_pre_draw = pre_draw;
  viewer.callback_key_down = key_down;

  viewer.launch();

  #endif
  return 0;
}

