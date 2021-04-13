
#include "BVH.h"

PRJ_BEGIN

BVH::BVH()
{
	
}
BVH::~BVH()
{
	
}
void BVH::InitObstacle(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F)
{
    unsigned int n_ob=F.rows();
    int dim = aabb_axis.size();
    /*
    for(int k=0;k<dim;k++)
    {
      aabb_axis[k].normalize();
    }
    */
    ob_tree=aabb::Tree(dim,0.0, n_ob,true);
    
    std::cout << "\nInserting ob particles into AABB tree ...\n";
    for (unsigned int i=0;i<n_ob;i++)
    {
        std::vector<double> lowerBound(dim);
        std::vector<double> upperBound(dim);

        for(int k=0;k<dim;k++)
        {
          upperBound[k]=-INFINITY;
          lowerBound[k]=INFINITY;
        }

        for(int j=0;j<3;j++)
        {
          for(int k=0;k<dim;k++)
          {
            double level = aabb_axis[k].dot(V.row(F(i,j)));
            if(level<lowerBound[k])
              lowerBound[k]=level;
            if(level>upperBound[k])
              upperBound[k]=level;
          }

        }


        ob_tree.insertParticle(i, lowerBound, upperBound);
    }

}

void BVH::InitPointcloud(const Eigen::MatrixXd& V)
{
    unsigned int n_pc=V.rows();
    int dim = aabb_axis.size();
    /*
    for(int k=0;k<dim;k++)
    {
      aabb_axis[k].normalize();
    }
    */
    pc_tree=aabb::Tree(dim,0.0, n_pc,true);
    
    std::cout << "\nInserting pc particles into AABB tree ...\n";
    for (unsigned int i=0;i<n_pc;i++)
    {
        std::vector<double> lowerBound(dim);
        std::vector<double> upperBound(dim);

        for(int k=0;k<dim;k++)
        {
          upperBound[k]=-INFINITY;
          lowerBound[k]=INFINITY;
        }

       
        for(int k=0;k<dim;k++)
        {
          double level = aabb_axis[k].dot(V.row(i));
          if(level<lowerBound[k])
            lowerBound[k]=level;
          if(level>upperBound[k])
            upperBound[k]=level;
        }

        


        pc_tree.insertParticle(i, lowerBound, upperBound);
    }

}

void BVH::InitTrajectory(const Data& spline)
{
    unsigned int n_tr=subdivide_tree.size();
    int dim = aabb_axis.size();
    tr_tree=aabb::Tree(dim, 0.0, n_tr,true);
 
    //std::cout << n_tr<<"\nInserting tr particles into AABB tree ...\n";
    for (unsigned int i=0;i<n_tr;i++)
    {
        std::vector<double> lowerBound(dim);
        std::vector<double> upperBound(dim);

        int sp_id=std::get<0>(subdivide_tree[i]);
        Eigen::MatrixXd basis=std::get<2>(subdivide_tree[i]);
        std::vector<Eigen::RowVector3d> P(order_num+1);
        
        Eigen::MatrixXd bz;
        bz=spline.block<order_num+1,3>(sp_id*(order_num-2),0);
        /*
        if(sp_id>0)
        {
          bz=M_convert*bz;
        }
        */
        for(int k=0;k<dim;k++)
        {
          upperBound[k]=-INFINITY;
          lowerBound[k]=INFINITY;
        }

        for(int j=0;j<=order_num;j++)
        {

            P[j].setZero();
            for(int j0=0;j0<=order_num;j0++)
            {
              P[j]+=basis(j,j0)*bz.row(j0);
            } 

            for(int k=0;k<dim;k++)
            {
              double level = aabb_axis[k].dot(P[j]);
              if(level<lowerBound[k])
                lowerBound[k]=level;
              if(level>upperBound[k])
                upperBound[k]=level;
            }
           
        }
        tr_tree.insertParticle(i, lowerBound, upperBound);
    }

}

void BVH::ccdInitTrajectory(const Data& spline, const Data& direction)
{
    unsigned int n_tr=subdivide_tree.size();
    int dim = aabb_axis.size();
    tr_tree=aabb::Tree(dim, 0.0, n_tr,true);
   
    //std::cout << "\nUpdate tr particles in AABB tree ...\n";
    for (unsigned int i=0;i<n_tr;i++)
    {
        std::vector<double> lowerBound(dim);
        std::vector<double> upperBound(dim);

        int sp_id=std::get<0>(subdivide_tree[i]);
        Eigen::MatrixXd basis=std::get<2>(subdivide_tree[i]);
        std::vector<Eigen::RowVector3d> P(order_num+1), D(order_num+1);
        
        Eigen::MatrixXd bz, bz_d;
        bz=spline.block<order_num+1,3>(sp_id*(order_num-2),0);
        bz_d=direction.block<order_num+1,3>(sp_id*(order_num-2),0);
        
        for(int k=0;k<dim;k++)
        {
          upperBound[k]=-INFINITY;
          lowerBound[k]=INFINITY;
        }
        for(int j=0;j<=order_num;j++)
        {   
            P[j].setZero();
            D[j].setZero();
            for(int j0=0;j0<=order_num;j0++)
            {
              P[j]+=basis(j,j0)*bz.row(j0);
              D[j]+=basis(j,j0)*bz_d.row(j0);
            } 
            for(int k=0;k<dim;k++)
            {
              double level = aabb_axis[k].dot(P[j]);
              if(level<lowerBound[k])
                lowerBound[k]=level;
              if(level>upperBound[k])
                upperBound[k]=level;

              level = aabb_axis[k].dot(P[j]+D[j]);
              if(level<lowerBound[k])
                lowerBound[k]=level;
              if(level>upperBound[k])
                upperBound[k]=level;
            }
        }

        tr_tree.insertParticle(i, lowerBound, upperBound);
    }

}

void BVH::UpdateTrajectory(const Data& spline)
{
    unsigned int n_tr=subdivide_tree.size();
    //tr_tree=aabb::Tree(3, n_tr);
    int dim = aabb_axis.size();
    //std::cout << "\nUpdate tr particles in AABB tree ...\n";
    for (unsigned int i=0;i<n_tr;i++)
    {
        std::vector<double> lowerBound(dim);
        std::vector<double> upperBound(dim);

        int id=std::get<0>(subdivide_tree[i]);
        Eigen::MatrixXd basis=std::get<2>(subdivide_tree[i]);
        std::vector<Eigen::RowVector3d> P(order_num+1);
         
        for(int k=0;k<dim;k++)
        {
          upperBound[k]=-INFINITY;
          lowerBound[k]=INFINITY;
        }
        for(int j=0;j<=order_num;j++)
        {
            
            P[j].setZero();
            for(int j0=0;j0<=order_num;j0++)
            {
              P[j]+=basis(j,j0)*spline.row(id+j0);
            } 

            for(int k=0;k<dim;k++)
            {
              double level = aabb_axis[k].dot(P[j]);
              if(level<lowerBound[k])
                lowerBound[k]=level;
              if(level>upperBound[k])
                upperBound[k]=level;
            }
            
        }


        tr_tree.updateParticle(i, lowerBound, upperBound);
    }

}

void BVH::ccdUpdateTrajectory(const Data& spline, const Data& direction)
{
    unsigned int n_tr=subdivide_tree.size();
    //tr_tree=aabb::Tree(3, n_tr);
    int dim = aabb_axis.size();

    //std::cout << "\nUpdate tr particles in AABB tree ...\n";
    for (unsigned int i=0;i<n_tr;i++)
    {
        std::vector<double> lowerBound(dim);
        std::vector<double> upperBound(dim);

        int id=std::get<0>(subdivide_tree[i]);
        Eigen::MatrixXd basis=std::get<2>(subdivide_tree[i]);
        std::vector<Eigen::RowVector3d> P(order_num+1), D(order_num+1);
        
        for(int k=0;k<dim;k++)
        {
          upperBound[k]=-INFINITY;
          lowerBound[k]=INFINITY;
        }
        for(int j=0;j<=order_num;j++)
        {

            P[j].setZero();
            D[j].setZero();
            for(int j0=0;j0<=order_num;j0++)
            {
              P[j]+=basis(j,j0)*spline.row(id+j0);
              D[j]+=basis(j,j0)*direction.row(id+j0);
            } 
            for(int k=0;k<dim;k++)
            {
              double level = aabb_axis[k].dot(P[j]);
              if(level<lowerBound[k])
                lowerBound[k]=level;
              if(level>upperBound[k])
                upperBound[k]=level;

              level = aabb_axis[k].dot(P[j]+D[j]);
              if(level<lowerBound[k])
                lowerBound[k]=level;
              if(level>upperBound[k])
                upperBound[k]=level;
            }

        }


        tr_tree.updateParticle(i, lowerBound, upperBound);
    }

}

void BVH::CheckCollision(std::vector<id_pair>& collision_pair, double margin)
{
    collision_pair=tr_tree.query(ob_tree, margin);
}

void BVH::pcCheckCollision(std::vector<id_pair>& collision_pair, double margin)
{
    collision_pair=tr_tree.query(pc_tree, margin);
}





PRJ_END

