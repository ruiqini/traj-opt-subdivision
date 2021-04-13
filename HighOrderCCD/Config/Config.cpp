#include "Config.h"

PRJ_BEGIN
std::vector<std::vector<long int>> combination;
//std::vector<double> element_len;
//int element_num;
int trajectory_num, piece_num;
int res;
int iter;
double epsilon2 , epsilon , 
       distance2, distance, 
       margin2, margin, 
       max_step, lambda, wolfe, offset, 
       gnorm, gtime,
       mu;
std::vector< std::tuple< int, std::pair<double,double>, Eigen::MatrixXd > > subdivide_tree, vel_tree, acc_tree;
bool automove, step_choose, adaptive_change, optimize_time;
//std::ofstream energy_file, result_file, gn_file, step_file, init_file, eigen_file, plot_file;
std::ofstream result_file;
Eigen::MatrixXd M_convert, M_head, M_tail, M_dynamic;

std::vector<double> time_weight;
double whole_weight;
std::vector<Eigen::MatrixXd> convert_list;

double piece_time, kt, ks;
double vel_limit, acc_limit;

std::vector<int> fixed_points;
std::vector<Eigen::Vector3d> aabb_axis={
                                          Eigen::Vector3d( 1, 0, 0),
                                          Eigen::Vector3d( 0, 1, 0),
                                          Eigen::Vector3d( 0, 0, 1),

                                          Eigen::Vector3d( 1, 1, 1),
                                          Eigen::Vector3d( 1,-1, 1),
                                          Eigen::Vector3d( 1, 1,-1),
                                          Eigen::Vector3d( 1,-1,-1),//
};
std::vector<Eigen::Vector3d> axis={
                                          Eigen::Vector3d( 1, 0, 0),
                                          Eigen::Vector3d( 0, 1, 0),
                                          Eigen::Vector3d( 0, 0, 1),//x y z

                                          Eigen::Vector3d( 1, 1, 1),
                                          Eigen::Vector3d( 1,-1, 1),
                                          Eigen::Vector3d( 1, 1,-1),
                                          Eigen::Vector3d( 1,-1,-1),//

                                          Eigen::Vector3d( 0, 1, 1),
                                          Eigen::Vector3d( 0, 1,-1),
                                          Eigen::Vector3d( 1, 0, 1),
                                          Eigen::Vector3d( 1, 0,-1),
                                          Eigen::Vector3d( 1, 1, 0),
                                          Eigen::Vector3d( 1,-1, 0),//
                                          
                                          
                                          Eigen::Vector3d( 0, 2, 1),
                                          Eigen::Vector3d( 0, 2,-1),
                                          Eigen::Vector3d( 0, 1, 2),
                                          Eigen::Vector3d( 0, 1,-2),

                                          Eigen::Vector3d( 2, 0, 1),
                                          Eigen::Vector3d( 2, 0,-1),
                                          Eigen::Vector3d( 1, 0, 2),
                                          Eigen::Vector3d( 1, 0,-2),

                                          Eigen::Vector3d( 2, 1, 0),
                                          Eigen::Vector3d( 2,-1, 0),
                                          Eigen::Vector3d( 1, 2, 0),
                                          Eigen::Vector3d( 1,-2, 0),
                                          
                                          Eigen::Vector3d( 1, 2, 1),
                                          Eigen::Vector3d( 1, 2,-1),
                                          Eigen::Vector3d( 1,-2, 1),
                                          Eigen::Vector3d(-1, 2, 1),

                                          Eigen::Vector3d( 1, 1, 2),
                                          Eigen::Vector3d( 1, 1,-2),
                                          Eigen::Vector3d( 1,-1, 2),
                                          Eigen::Vector3d(-1, 1, 2),

                                          Eigen::Vector3d( 2, 1, 1),
                                          Eigen::Vector3d( 2, 1,-1),
                                          Eigen::Vector3d( 2,-1, 1),
                                          Eigen::Vector3d(-2, 1, 1),

                                          Eigen::Vector3d( 2, 2, 1),
                                          Eigen::Vector3d( 2, 2,-1),
                                          Eigen::Vector3d( 2,-2, 1),
                                          Eigen::Vector3d(-2, 2, 1),

                                          Eigen::Vector3d( 2, 1, 2),
                                          Eigen::Vector3d( 2, 1,-2),
                                          Eigen::Vector3d( 2,-1, 2),
                                          Eigen::Vector3d(-2, 1, 2),

                                          Eigen::Vector3d( 1, 2, 2),
                                          Eigen::Vector3d( 1, 2,-2),
                                          Eigen::Vector3d( 1,-2, 2),
                                          Eigen::Vector3d(-1, 2, 2),

                                          };

PRJ_END