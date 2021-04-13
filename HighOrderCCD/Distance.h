#ifndef DISTANCE_H
#define DISTANCE_H

#include "Config/Config.h"

PRJ_BEGIN

template <typename Scalar, int ADIM>
class Distance 
{
  public:
  typedef Eigen::Matrix<Scalar,1,ADIM> Point;
  typedef Eigen::Matrix<double,1,ADIM> ConstPoint;

  static void point_point(const Point& P, const Point& S, Scalar &d, Point& C)
  { 
    C = S;
    d=(P-C).norm();
  }

  static void point_segment(const Point& P, const Point& S0, const Point& S1, Scalar &d, Point& C)
  { 
    Point S=S1-S0;
    Scalar t = (P-S0).dot(S)/S.squaredNorm();//((P(0,0)-S0(0,0))*S(0,0)+(P(0,1)-S0(0,1))*S(0,1))/S.squaredNorm();
    if(S.squaredNorm()==0)
    {
      C = S0;
      d=(P-C).norm();
    }
    else
    {
      if(t < 0)
      {
          C = S0;
          d=(P-C).norm();
      }
      else if(t>1)
      {
        C = S1;
        d=(P-C).norm();
      }
      else
      {
        C = S0+S*t;
        d=(P-C).norm();
      }
    }
  }

  static void point_constsegment(const Point& P, const ConstPoint& S0, const ConstPoint& S1, Scalar &d, Point& C)
  { 
    ConstPoint S=S1-S0;
    if(S.squaredNorm()==0)
    {
      for(int i=0;i<ADIM;i++)
      {
        C(i)=S0(i);
      }
      d=(P-C).norm();
    }
    else
    {
      Scalar t(0);
      for(int i=0;i<ADIM;i++)
      {
        t+=(P(i)-S0(i))*S(i);
      }
      t/=S.squaredNorm();//((P(0,0)-S0(0,0))*S(0,0)+(P(0,1)-S0(0,1))*S(0,1))/S.squaredNorm();
      
      if(t < 0)
      {
        for(int i=0;i<ADIM;i++)
        {
          C(i)=S0(i);
        }
        d=(P-C).norm();
      }
      else if(t>1)
      {
        for(int i=0;i<ADIM;i++)
        {
          C(i)=S1(i);
        }
        d=(P-C).norm();
      }
      else
      {
        for(int i=0;i<ADIM;i++)
        {
          C(i)=S0(i)+t*S(i);
        }
        d=(P-C).norm();
      }
    }
    
    
  }

  static void segment_segment_test(const Point& P0, const Point& P1, const ConstPoint& S0, const ConstPoint& S1, Scalar &d,Point& C0,Point& C1)
  {
    /*
    ConstPoint S=S1-S0;
    
    double len0=S.squaredNorm();
    //Scalar t0 = (P0-S0).dot(S)/len0; //((P0(0,0)-S0(0,0))*S(0,0)+(P0(0,1)-S0(0,1))*S(0,1))/len;
   // Scalar t1 = (P1-S0).dot(S)/len0; //((P1(0,0)-S0(0,0))*S(0,0)+(P1(0,1)-S0(0,1))*S(0,1))/len;
    if(len0==0)
      std::cout<<"sserr-------------------"<<std::endl;
    Scalar t0(0);
    for(int i=0;i<ADIM;i++)
    {
      t0+=(P0(i)-S0(i))*S(i);
    }
    t0/=len0;

    Scalar t1(0);
    for(int i=0;i<ADIM;i++)
    {
      t1+=(P1(i)-S0(i))*S(i);
    }
    t1/=len0;

    Point PlaneP0;
    for(int i=0;i<ADIM;i++)
    {
      PlaneP0(i)=P0(i)-t0*S(i);
    }
    Point PlaneP1;
    for(int i=0;i<ADIM;i++)
    {
      PlaneP1(i)=P1(i)-t1*S(i);
    }
    Point PlaneP=PlaneP1-PlaneP0;
    Scalar len1=PlaneP.squaredNorm();

    Scalar t(0);
    if(len1>1e-15)
    {
      for(int i=0;i<ADIM;i++)
      {
        t+=(S0(i)-PlaneP0(i))*PlaneP(i);
      }
      t/=len1;
      //t=(S0-PlaneP0).dot(PlaneP)/len1; //((S0(0,0)-PlaneP0(0,0))*PlaneP(0,0)+(S0(0,1)-PlaneP0(0,1))*PlaneP(0,1))/PlaneP.squaredNorm();

      Point P_;
      if(t<0)
      {
        P_=P0;
      }
      else if(t>1)
      {
        P_=P1;
      }
      else
      {
        P_=P0+t*(P1-P0);
      }
      //Point C0,C1;
      point_constsegment(P_,S0,S1,d,C0);
      point_segment(C0,P0,P1,d,C1);
    }
    else
    {
      std::cout<<"ss_err\n";
      d=PlaneP0.norm();
    }
    */
    
    Point P=P1-P0;
    
    Scalar len0=P.squaredNorm();

    ConstPoint S=S1-S0;
    
    double len2=S.squaredNorm();
    Scalar dot(0);
    for(int i=0;i<ADIM;i++)
    {
      dot+=S(i)*P(i);
    }
    
    //Scalar t0 = (P0-S0).dot(S)/len0; //((P0(0,0)-S0(0,0))*S(0,0)+(P0(0,1)-S0(0,1))*S(0,1))/len;
   // Scalar t1 = (P1-S0).dot(S)/len0; //((P1(0,0)-S0(0,0))*S(0,0)+(P1(0,1)-S0(0,1))*S(0,1))/len;
    if(len0==0)
      std::cout<<"sserr-------------------"<<std::endl;
    Scalar t0(0);
    for(int i=0;i<ADIM;i++)
    {
      t0+=(S0(i)-P0(i))*P(i);
    }
    t0/=len0;

    Scalar t1(0);
    for(int i=0;i<ADIM;i++)
    {
      t1+=(S1(i)-P0(i))*P(i);
    }
    t1/=len0;

    Point PlaneS0;
    for(int i=0;i<ADIM;i++)
    {
      PlaneS0(i)=S0(i)-t0*P(i);
    }
    Point PlaneS1;
    for(int i=0;i<ADIM;i++)
    {
      PlaneS1(i)=S1(i)-t1*P(i);
    }
    Point PlaneS=PlaneS1-PlaneS0;
    Scalar len1=PlaneS.squaredNorm();

    Scalar t(0);
    if(1-(dot*dot)/(len0*len2)>1e-8)
    {
      for(int i=0;i<ADIM;i++)
      {
        t+=(P0(i)-PlaneS0(i))*PlaneS(i);
      }
      t/=len1;
      //t=(S0-PlaneP0).dot(PlaneP)/len1; //((S0(0,0)-PlaneP0(0,0))*PlaneP(0,0)+(S0(0,1)-PlaneP0(0,1))*PlaneP(0,1))/PlaneP.squaredNorm();

      Point S_;
      if(t<0)
      {
        for(int i=0;i<ADIM;i++)
        {
          S_(i)=S0(i);
        }
      }
      else if(t>1)
      {
        for(int i=0;i<ADIM;i++)
        {
          S_(i)=S1(i);
        }
      }
      else
      {
        for(int i=0;i<ADIM;i++)
        {
          S_(i)=S0(i)+t*(S1-S0)(i);
        }
      }
      //Point C0,C1;
      point_segment(S_,P0,P1,d,C0);
      point_constsegment(C0,S0,S1,d,C1);
    }
    else
    {
      std::cout<<"ss_err\n";
      d=PlaneS0.norm();
    }
    
  }

  static void segment_segment(const Point& P0, const Point& P1, const ConstPoint& S0, const ConstPoint& S1, Scalar &d,Point& C0,Point& C1)
  { 
    
    ConstPoint S10=S1-S0;
    double S1x=S10(0); double S1y=S10(1); double S1z=S10(2);
    Scalar P1x=P1(0)-P0(0); Scalar P1y=P1(1)-P0(1); Scalar P1z=P1(2)-P0(2);

    Scalar a00=P1x*P1x+P1y*P1y+P1z*P1z;
    Scalar a01=-P1x*S1x-P1y*S1y-P1z*S1z;
    Scalar a10=a01;
    Scalar a11=Scalar(S1x*S1x+S1y*S1y+S1z*S1z);
    
    Scalar P0x=P0(0)-Scalar(S0(0)); 
    Scalar P0y=P0(1)-Scalar(S0(1));
    Scalar P0z=P0(2)-Scalar(S0(2));

    Scalar b0=-P1x*P0x-P1y*P0y-P1z*P0z;
    Scalar b1=S1x*P0x+S1y*P0y+S1z*P0z;

    Scalar det=a00*a11-a01*a10;
    
    Scalar t0,t1;

    if(det!=0)
    {
      t0=( a11*b0-a01*b1)/det;
      t1=(-a10*b0+a00*b1)/det;
      
      if(t0>=0&&t0<=1 && t1>=0&&t1<=1)
      {
        C0=P0+t0*(P1-P0);
        C1(0)=S0(0)+t1*S10(0);
        C1(1)=S0(1)+t1*S10(1);
        C1(2)=S0(2)+t1*S10(2);
        d=(C0-C1).norm();
      }
      else if(t0<0)
      {
        C0=P0;
        point_constsegment(C0, S0, S1,d,C1);
        point_segment(C1,P0,P1,d,C0);
      }
      else if(t0>1)
      {
        C0=P1;
        point_constsegment(C0, S0, S1,d,C1);
        point_segment(C1,P0,P1,d,C0);
      }
      else if(t1<0)
      {
        C1(0)=S0(0);
        C1(1)=S0(1);
        C1(2)=S0(2);
        point_segment(C1,P0,P1,d,C0);
        point_constsegment(C0, S0, S1,d,C1);

      }
      else if(t1>1)
      {
        C1(0)=S1(0);
        C1(1)=S1(1);
        C1(2)=S1(2);
        point_segment(C1,P0,P1,d,C0);
        point_constsegment(C0, S0, S1,d,C1);

      }
    }
    else
    {
      Scalar d0,d1;
      point_constsegment(P0, S0, S1,d0,C1);
      point_constsegment(P1, S0, S1,d1,C1);
      if(d0<d1)
      {
        d=d0;
      }
      else
      {
        d=d1;
      }
      
    }

  }

  static void segment_constpoint(const Point& P0, const Point& P1, const ConstPoint& S, Scalar &d, Point& C)
  { 
    Point P=P1-P0;
    Scalar t(0);
    for(int i=0;i<ADIM;i++)
    {
      t+=(S(i)-P0(i))*P(i);
    }
    t/=P.squaredNorm();//((P(0,0)-S0(0,0))*S(0,0)+(P(0,1)-S0(0,1))*S(0,1))/S.squaredNorm();
    
    Point D;
    if(t < 0)
    {
      C=P0;
      for(int i=0;i<ADIM;i++)
      {
        D(i)=S(i)-C(i);
      }
      d=D.norm();
    }
    else if(t>1)
    {
      C=P1;
      for(int i=0;i<ADIM;i++)
      {
        D(i)=S(i)-C(i);
      }
      d=D.norm();
    }
    else
    {
      C=P0+t*P;
      for(int i=0;i<ADIM;i++)
      {
        D(i)=S(i)-C(i);
      }
      d=D.norm();
    }
  }

  static void point_face(const Point& P, const ConstPoint& S0, const ConstPoint& S1, const ConstPoint& S2, Scalar &d,Point& C)
  { 
    ConstPoint S10=S1-S0;
    ConstPoint S20=S2-S0;

    double S1x=S10(0); double S1y=S10(1); double S1z=S10(2);
    double S2x=S20(0); double S2y=S20(1); double S2z=S20(2);

    double a00=S1x*S1x+S1y*S1y+S1z*S1z;
    double a01=S1x*S2x+S1y*S2y+S1z*S2z;
    double a10=a01;
    double a11=S2x*S2x+S2y*S2y+S2z*S2z;
    
    Scalar P0x=P(0)-Scalar(S0(0)); 
    Scalar P0y=P(1)-Scalar(S0(1));
    Scalar P0z=P(2)-Scalar(S0(2));

    Scalar b0=S1x*P0x+S1y*P0y+S1z*P0z;
    Scalar b1=S2x*P0x+S2y*P0y+S2z*P0z;

    double det=a00*a11-a01*a10;

    Scalar t0,t1,t2;

    if(det!=0)
    {
      t0=( a11*b0-a01*b1)/det;
      t1=(-a10*b0+a00*b1)/det;
      t2=1-t0-t1;
      /*
      std::cout.precision(50);
      std::cout<<t0<<std::endl;
      std::cout<<t1<<std::endl;
      std::cout<<t2<<std::endl;
      */
      if(t0>=0&&t1>=0&&t2>=0)
      {
        C(0)=S0(0)+t0*S1x+t1*S2x;
        C(1)=S0(1)+t0*S1y+t1*S2y;
        C(2)=S0(2)+t0*S1z+t1*S2z;
        d=(P-C).norm();
      }
      else
      {
        double cos0=(S1-S0).dot(S2-S0);
        double cos1=(S0-S1).dot(S2-S1);
        double cos2=(S1-S2).dot(S0-S2);
        if(cos0>0&&cos1>0&&cos2>0)
        {
          if(t0<=0)
          {
            point_constsegment(P, S0, S2,d,C);
          }
          else if(t1<=0)
          {
            point_constsegment(P, S0, S1,d,C);
          }
          else if(t2<=0)
          {
            point_constsegment(P, S1, S2,d,C);
          }
        }
        else 
        {
          Scalar d0, d1, d2;
          Point C0, C1, C2;
          point_constsegment(P, S0, S2,d0,C0);
          point_constsegment(P, S0, S1,d1,C1);
          point_constsegment(P, S1, S2,d2,C2);
          if(d0<=d1)
          {
            if(d0<=d2)
            {
              d=d0;
              C=C0;
            }
            else
            {
              d=d2;
              C=C2;
            }
          }
          else
          {
            if(d1<=d2)
            {
              d=d1;
              C=C1;
            }
            else
            {
              d=d2;
              C=C2;
            }
          }
        }

      }
      /*
      else if(t0<=0)
      {
        point_constsegment(P, S0, S2,d,C);
      }
      else if(t1<=0)
      {
        point_constsegment(P, S0, S1,d,C);
      }
      else if(t2<=0)
      {
        point_constsegment(P, S1, S2,d,C);
      }
      */
    }
    else
    {
      //std::cout<<"nf\n";
      /*
      ConstPoint S21=S2-S1;

      double s10=S10.squaredNorm();
      double s20=S20.squaredNorm();
      double s21=S21.squaredNorm();
      if(s10>=s20&&s10>=s21)
      {
        point_constsegment(P, S0, S1,d,C);
      }
      else if(s20>=s10&&s20>=s21)
      {
        point_constsegment(P, S0, S2,d,C);
      }
      else if(s21>=s20&&s21>=s10)
      {
        point_constsegment(P, S1, S2,d,C);
      }
      */
      Scalar d0, d1, d2;
      Point C0, C1, C2;
      point_constsegment(P, S0, S2,d0,C0);
      point_constsegment(P, S0, S1,d1,C1);
      point_constsegment(P, S1, S2,d2,C2);
      if(d0<=d1)
      {
        if(d0<=d2)
        {
          d=d0;
          C=C0;
        }
        else
        {
          d=d2;
          C=C2;
        }
      }
      else
      {
        if(d1<=d2)
        {
          d=d1;
          C=C1;
        }
        else
        {
          d=d2;
          C=C2;
        }
      }
    }
    
  }

  static void face_point(const Point& P0,const Point& P1,const Point& P2, const ConstPoint& S, Scalar &d,Point& C)
  { 

    Scalar P1x=P1(0)-P0(0); Scalar P1y=P1(1)-P0(1); Scalar P1z=P1(2)-P0(2);
    Scalar P2x=P2(0)-P0(0); Scalar P2y=P2(1)-P0(1); Scalar P2z=P2(2)-P0(2);

    Scalar a00=P1x*P1x+P1y*P1y+P1z*P1z;
    Scalar a01=P1x*P2x+P1y*P2y+P1z*P2z;
    Scalar a10=a01;
    Scalar a11=P2x*P2x+P2y*P2y+P2z*P2z;

    Scalar S0x=S(0)-P0(0); 
    Scalar S0y=S(1)-P0(1);
    Scalar S0z=S(2)-P0(2);

    Scalar b0=P1x*S0x+P1y*S0y+P1z*S0z;
    Scalar b1=P2x*S0x+P2y*S0y+P2z*S0z;

    Scalar det=a00*a11-a01*a10;
    
    Scalar t0,t1,t2;

    if(det!=0)
    {
      t0=( a11*b0-a01*b1)/det;
      t1=(-a10*b0+a00*b1)/det;
      t2=1-t0-t1;
      if(t0>=0&&t1>=0&&t2>=0)
      {
        C(0)=P0(0)+t0*P1x+t1*P2x;
        C(1)=P0(1)+t0*P1y+t1*P2y;
        C(2)=P0(2)+t0*P1z+t1*P2z;
        Point D;
        D(0)=S(0)-C(0); D(1)=S(1)-C(1); D(2)=S(2)-C(2);
        d=D.norm();
      }
      else
      {
        Scalar cos0=(P1-P0).dot(P2-P0);
        Scalar cos1=(P0-P1).dot(P2-P1);
        Scalar cos2=(P1-P2).dot(P0-P2);
        if(cos0>0&&cos1>0&&cos2>0)
        {
          if(t0<=0)
          {
            segment_constpoint(P0, P2, S, d, C);
          }
          else if(t1<=0)
          {
            segment_constpoint(P0, P1, S, d, C);
          }
          else if(t2<=0)
          {
            segment_constpoint(P1, P2, S, d, C);
          }
        }
        else
        {
          Scalar d0, d1, d2;
          Point C0, C1, C2;
          segment_constpoint(P0, P2, S, d0, C0);
          segment_constpoint(P0, P1, S, d1, C1);
          segment_constpoint(P1, P2, S, d2, C2);
          if(d0<=d1)
          {
            if(d0<=d2)
            {
              d=d0;
              C=C0;
            }
            else
            {
              d=d2;
              C=C2;
            }
          }
          else
          {
            if(d1<=d2)
            {
              d=d1;
              C=C1;
            }
            else
            {
              d=d2;
              C=C2;
            }
          }
        }
      }
      /*
      else if(t0<=0)
      {
        segment_constpoint(P0, P2, S, d, C);
      }
      else if(t1<=0)
      {
        segment_constpoint(P0, P1, S, d, C);
      }
      else if(t2<=0)
      {
        segment_constpoint(P1, P2, S, d, C);
      }
      */
    }
    else
    {
      Scalar d0, d1, d2;
      Point C0, C1, C2;
      segment_constpoint(P0, P2, S, d0, C0);
      segment_constpoint(P0, P1, S, d1, C1);
      segment_constpoint(P1, P2, S, d2, C2);
      if(d0<=d1)
      {
        if(d0<=d2)
        {
          d=d0;
          C=C0;
        }
        else
        {
          d=d2;
          C=C2;
        }
      }
      else
      {
        if(d1<=d2)
        {
          d=d1;
          C=C1;
        }
        else
        {
          d=d2;
          C=C2;
        }
      }
    
    }
  }

  static bool segment_face(const Point& P0, const Point& P1, const Point& S0, const Point& S1, const Point& S2)
  {
    Eigen::Matrix<Scalar,3,3> A;
    A.col(0)=(P1-P0).transpose();
    A.col(1)=(S0-S2).transpose();
    A.col(2)=(S1-S2).transpose();
    Eigen::Matrix<Scalar,3,1> b,x;
    b=(P1-S2).transpose();

    Scalar t,t0,t1;
    if(A.determinant()!=0)
    {
      x=A.inverse()*b;
      t=x(0);
      t0=x(1);
      t1=x(2);
      if(t>0 && t<1 && t0>0&&t1>0 && 1-t0-t1>0)
        return true;
    }
    else
    {
      Eigen::Matrix<Scalar,2,2> B;
      Eigen::Matrix<Scalar,2,1> b1,x1;

      B(0,0)=P1(0)-P0(0); B(0,1)=S1(0)-S0(0);
      B(1,0)=P1(1)-P0(1); B(1,1)=S1(1)-S0(1);
      b1(0)=P1(0)-S1(0);
      b1(1)=P1(1)-S1(1);
      if(B.determinant()!=0)
      {
        x1=B.inverse()*b1;
        t=x1(0);
        t0=x1(1);
        if(t>0 && t<1 && t0>0&&t0<1)
          return true;
      }
      

      B(0,0)=P1(0)-P0(0); B(0,1)=S2(0)-S0(0);
      B(1,0)=P1(1)-P0(1); B(1,1)=S2(1)-S0(1);
      b1(0)=P1(0)-S2(0);
      b1(1)=P1(1)-S2(1);
      if(B.determinant()!=0)
      {
        x1=B.inverse()*b1;
        t=x1(0);
        t0=x1(1);
        if(t>0 && t<1 && t0>0&&t0<1)
          return true;
      }
      
      B(0,0)=P1(0)-P0(0); B(0,1)=S2(0)-S1(0);
      B(1,0)=P1(1)-P0(1); B(1,1)=S2(1)-S1(1);
      b1(0)=P1(0)-S2(0);
      b1(1)=P1(1)-S2(1);
      if(B.determinant()!=0)
      {
        x1=B.inverse()*b1;
        t=x1(0);
        t0=x1(1);
        if(t>0 && t<1 && t0>0&&t0<1)
          return true;
      }

    }
    


    return false;
  }

};

PRJ_END

#endif
