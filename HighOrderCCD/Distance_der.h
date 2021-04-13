#ifndef DISTANCE_DER_H
#define DISTANCE_DER_H

#include "Config/Config.h"

PRJ_BEGIN

class Distance_der
{
  public:
  typedef Eigen::RowVector3d Point;

  static void point_point(const Point& P, const Point& S, double &d, Eigen::RowVector3d &d_p, Eigen::Matrix3d &h_p)
  {
      d=(P-S).norm();
        
      d_p=(P-S)/d;
        
      Eigen::Matrix3d I; I.setIdentity();
      h_p=I/d-(P-S).transpose()*(P-S)/std::pow(d,3);
  }

  static void point_segment(const Point& P, const Point& S0, const Point& S1, double &d, Eigen::RowVector3d &d_p, Eigen::Matrix3d &h_p)
  {
      Point S=S1-S0;

      double s2=S.squaredNorm();
      if(s2==0)
      {
          point_point(P,S0,d,d_p,h_p);
      }
      else
      {
          double t = (P-S0).dot(S)/s2;
          if(t < 0)
          {
              point_point(P,S0,d,d_p,h_p);
          }
          else if(t > 1)
          {
              point_point(P,S1,d,d_p,h_p);
          }
          else
          {
              Eigen::Matrix3d I; I.setIdentity();
              Eigen::Matrix3d B = I-S.transpose()*S/s2;
              Point D=(P-S0)*B;
              d=D.norm();
               
              Point DB=D*B;
              d_p=DB/d;

              h_p=B*B/d-DB.transpose()*DB/std::pow(d,3);
          }

      } 
  }

  static void point_face(const Point& P, const Point& S0, const Point& S1, const Point& S2, double &d, Eigen::RowVector3d &d_p, Eigen::Matrix3d &h_p)
  {
      Point S10=S1-S0, S20=S2-S0;
      double s1=S10.norm(), s2=S20.norm();

      double a00=s1*s1;
      double a01=S10.dot(S20);
      double a10=a01;
      double a11=s2*s2;

      double b0=(P-S0).dot(S10);
      double b1=(P-S0).dot(S20);

      double det=a00*a11-a01*a10;
      double t0,t1,t2;
      
      if(det==0)
      {
          double d0, d1, d2;
          Eigen::RowVector3d d_p0, d_p1, d_p2;
          Eigen::Matrix3d h_p0, h_p1, h_p2;
          point_segment(P, S0, S1, d0, d_p0, h_p0);
          point_segment(P, S1, S2, d1, d_p1, h_p1);
          point_segment(P, S2, S0, d2, d_p2, h_p2);
          if(d0<=d1 && d0<=d2)
          {
              d=d0; d_p=d_p0; h_p=h_p0;
          }
          else if(d1<=d0 && d1<=d2)
          {
              d=d1; d_p=d_p1; h_p=h_p1;
          }
          else if(d2<=d0 && d2<=d1)
          {
              d=d2; d_p=d_p2; h_p=h_p2;
          }
      }
      else
      {  
        Point N=S10.cross(S20);
        N.normalize();

        t0=( a11*b0-a01*b1)/det;
        t1=(-a10*b0+a00*b1)/det;
        t2=1-t0-t1;
        
        if(t0>=0&&t1>=0&&t2>=0)
        {
            Eigen::Matrix3d B = N.transpose()*N;
            Point D=(P-S0)*B;
            d=D.norm();

            Point DB=D*B;
            d_p=DB/d;

            h_p=B*B/d-DB.transpose()*DB/std::pow(d,3);
        }
        else
        {
            double d0, d1, d2;
            Eigen::RowVector3d d_p0, d_p1, d_p2;
            Eigen::Matrix3d h_p0, h_p1, h_p2;
            point_segment(P, S0, S1, d0, d_p0, h_p0);
            point_segment(P, S1, S2, d1, d_p1, h_p1);
            point_segment(P, S2, S0, d2, d_p2, h_p2);
            if(d0<=d1 && d0<=d2)
            {
                d=d0; d_p=d_p0; h_p=h_p0;
            }
            else if(d1<=d0 && d1<=d2)
            {
                d=d1; d_p=d_p1; h_p=h_p1;
            }
            else if(d2<=d0 && d2<=d1)
            {
                d=d2; d_p=d_p2; h_p=h_p2;
            }
        }  
      }
  }
};

PRJ_END

#endif
