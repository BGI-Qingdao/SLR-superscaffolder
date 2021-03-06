#ifndef __ALGORITHM_LINEAR_FITTING_LINEAR_H__
#define __ALGORITHM_LINEAR_FITTING_LINEAR_H__
/**********************************************************
 *
 * @Brief :
 *   A simple linear fitting code.
 *
 * ********************************************************/
#include <vector>
#include <cmath>
namespace BGIQD {
    namespace LINEARFITTING {

        // A line : 0 = ax + by + c
        template<class X , class Y>
            struct Linear
            {
                // 0 = ax + by + c
                double a ;
                double b ;
                double c ;

                Y getY(const X & x ) const 
                {
                    return ( -a *x -c ) / b ;
                }
                X getX(const Y & y ) const 
                {
                    return ( -b * y - c ) / a ;
                }
            };
        // A point (x,y)
        template<class X , class Y>
            struct Item
            {
                typedef X XType;
                typedef Y YType;
                XType x ;
                YType y ;
            };

        // LinearFit by least square method
        template<class Item>
            Linear<typename Item::XType , typename Item::YType> 
                lineFit(const std::vector<Item> &points)
            {
                Linear<typename Item::XType , typename Item::YType> ret ;
                int size = points.size();
                if(size < 2)
                {
                    ret.a = 0;
                    ret.b = 0;
                    ret.c = 0;
                    return ret;
                }
                double x_mean = 0;
                double y_mean = 0;
                for(int i = 0; i < size; i++)
                {
                    x_mean += points[i].x;
                    y_mean += points[i].y;
                }
                x_mean /= size;
                y_mean /= size;

                double Dxx = 0, Dxy = 0, Dyy = 0;

                for(int i = 0; i < size; i++)
                {
                    Dxx += (points[i].x - x_mean) * (points[i].x - x_mean);
                    Dxy += (points[i].x - x_mean) * (points[i].y - y_mean);
                    Dyy += (points[i].y - y_mean) * (points[i].y - y_mean);
                }
                double lambda = ( (Dxx + Dyy) - sqrt( (Dxx - Dyy) * (Dxx - Dyy) + 4 * Dxy * Dxy) ) / 2.0;
                double den = sqrt( Dxy * Dxy + (lambda - Dxx) * (lambda - Dxx) );
                ret.a = Dxy / den;
                ret.b = (lambda - Dxx) / den;
                ret.c = - ret.a * x_mean - ret.b * y_mean;
                return ret;
            } 
    }
}

#endif //__ALGORITHM_LINEAR_FITTING_LINEAR_H__
