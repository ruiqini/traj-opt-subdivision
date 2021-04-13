#ifndef BVH_H
#define BVH_H

#include "HighOrderCCD/Config/Config.h"
//#include "HighOrderCCD/Element.h"
//#include "HighOrderCCD/Subdivide.h"
//#include "HighOrderCCD/ElementCD.h"
//#include "HighOrderCCD/Distance.h"
#include "src/AABB.h"

PRJ_BEGIN

class BVH 
{
  public:
    typedef std::vector< std::tuple< int, std::pair<double,double>, Eigen::MatrixXd > > SubdivideTree;
    typedef Eigen::MatrixXd Data;
    typedef std::pair<unsigned int, unsigned int> id_pair;

    aabb::Tree tr_tree;
    aabb::Tree ob_tree;
    aabb::Tree pc_tree;

    BVH();
	  ~BVH();

    void InitObstacle(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F);

    void InitPointcloud(const Eigen::MatrixXd& V);

    void InitTrajectory(const Data& spline);
    
    void ccdInitTrajectory(const Data& spline, const Data& direction);

    void UpdateTrajectory(const Data& spline);
    
    void ccdUpdateTrajectory(const Data& spline, const Data& direction);

    void CheckCollision(std::vector<id_pair>& collision_pair, double margin);
    
    void pcCheckCollision(std::vector<id_pair>& collision_pair, double margin);

};

PRJ_END

#endif