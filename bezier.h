//
//  bezier.h
//
//
//  Created by Enrico Eberhard on 01/06/2018.
//
//
//	A class for creating and using bezier curves
//
//

#ifndef bezier_h
#define bezier_h

#include "math.h"


/*-----------------------------------------------------------------*/
/*                          Point Class                            */
/*-----------------------------------------------------------------*/

class Point {
public:
    float x, y;
    Point(void) : x(0), y(0) {}
    Point(float x, float y) : x(x), y(y) {}
};


Point operator + (Point a, Point b) {
    return Point(a.x + b.x, a.y + b.y);
}

Point operator - (Point a, Point b) {
    return Point(a.x - b.x, a.y - b.y);
}

Point operator * (float s, Point a) {
    return Point(s * a.x, s * a.y);
}


/*-----------------------------------------------------------------*/
/*                          Bezier Class                           */
/*-----------------------------------------------------------------*/

#define BezierMaxOrder 20		//maximum allowable order
#define BezierDefaultOrder 5	//default order
#define BezierMaxIter 50		//maximum iterations of root solver
#define BezierSolverTol 1e-6	//tolerance of root solver

class Bezier {
    
public:
    typedef float (Bezier::*Bez_h_fun)(float, float);
    
    Bezier(): order(BezierDefaultOrder) {calculateBinomials();};
    Bezier(int n): order(n) {calculateBinomials();};
    Bezier(float *x, float *y): order(BezierDefaultOrder) {calculateBinomials(); setControls(x, y);};
    Bezier(int n, float *x, float *y);
    
    void setOrder(int n) {order = n; calculateBinomials();};
    void setControls(float *x, float *y);
    Point getControl(int k) {if((k<=order)&&(k>=0)) return Controls[k]; else return Point();};
    void printControls(void);
    void printOrder(void) {printf("Order %i\n", order);};
    
    void setCorners(Point p0, float dydx0, Point pEnd, float dydxEnd, float roundness);
    
    Point getPoint(float t);
    float getX(float t);
    float getY(float t);
    
    Point getDerivativePoint(float t);
    float getDerivativeX(float t);
    float getDerivativeY(float t);
    
    
    //root finding
    float getTFromX(float x, float t0) {return newton_raphson(x, t0, &Bezier::h_fun_X);};
    float getTFromY(float y, float t0) {return newton_raphson(y, t0, &Bezier::h_fun_Y);};
    float getTFromX(float x) {return getTFromX(x, 0.5);};
    float getTFromY(float y) {return getTFromY(y, 0.5);};
    
    float getYFromX(float x, float t0) {return getY(getTFromX(x, t0));};
    float getXFromY(float y, float t0) {return getX(getTFromY(y, t0));};
    float getYFromX(float x) {return getYFromX(x, 0.5);};
    float getXFromY(float y) {return getXFromY(y, 0.5);};
    
    Point getPointFromX(float x, float t0) {return getPoint(getTFromX(x, t0));};
    Point getPointFromY(float y, float t0) {return getPoint(getTFromY(y, t0));};
    Point getPointFromX(float x) {return getPointFromX(x, 0.5);};
    Point getPointFromY(float y) {return getPointFromY(y, 0.5);};
    
    
    
private:
    int order;
    Point Controls[BezierMaxOrder + 1];
    int nchoosek[BezierMaxOrder + 1][BezierMaxOrder + 1];
    int factorial(int x);
    void calculateBinomials(void);
    float bernsteinBasis(int i, int n, float t);
    float h_fun_X(float t, float x);
    float h_fun_Y(float t, float y);
    float newton_raphson(float r, float t0, Bez_h_fun h_fun);
    
};

/* Constructor */
Bezier::Bezier(int n, float *x, float *y) {
    
    if (n > BezierMaxOrder)
        mju_warning("Requested Bezier curve order exceeds default maximum!");
    
    order = n;
    calculateBinomials();
    setControls(x, y);
}


/* Utilities */
void Bezier::printControls(void) {
    int i;
    
    printf("Control points (x,y): \n");
    for(i=0; i<(order + 1); i++)
        printf("%3.2f, %3.2f\n", Controls[i].x, Controls[i].y);
}

void Bezier::setControls(float *x, float *y) {
    
    int i;
    for (i = 0; i < (order + 1); i++) {
        Controls[i].x = x[i];
        Controls[i].y = y[i];
    }
}


/* Generate curve from corner points */

void Bezier::setCorners(Point p0, float dydx0, Point pEnd, float dydxEnd, float roundness) {
    
    Controls[0] = p0;
    
    Controls[1] = p0 + Point(1, dydx0);
    
    Controls[order - 1] = pEnd - Point(1, dydxEnd);
    
    Controls[order] = pEnd;
    
    
    
    
}




/* General Bezier curve formulas */
Point Bezier::getPoint(float t) {
    Point tmp[order + 1];
    memcpy(tmp, Controls, (order + 1) * sizeof(Point));
    
    
    int i;
    for (i = order; i > 0; i--) {
        for (int k = 0; k < i; k++)
            tmp[k] = tmp[k] + t * ( tmp[k+1] - tmp[k] );
    }
    
    
    Point p = Point(tmp[0].x, tmp[0].y);
    
    //delete tmp;
    
    return p;
}

float Bezier::getX(float t) {
    
    Point p = getPoint(t);
    
    return p.x;
}

float Bezier::getY(float t) {
    
    Point p = getPoint(t);
    
    return p.y;
}






/* Derivative of Bezier curve */
Point Bezier::getDerivativePoint(float t) {
    
    Point p;
    float b;
    
    int i;
    for (i = 0; i < order; i++) {
        b = bernsteinBasis(i,order-1,t);
        p = p + b * ( Controls[i+1] - Controls[i]);
    }
    
    return order * p;
}

float Bezier::getDerivativeX(float t) {
    
    Point p = getDerivativePoint(t);
    
    return p.x;
}

float Bezier::getDerivativeY(float t) {
    
    Point p = getDerivativePoint(t);
    
    return p.y;
}


/* Math functions */

void Bezier::calculateBinomials(void) {
    
    int i, n;
    
    for (n=1; n<(order + 1); n++) {
        for (i=0; i<=n; i++) {
            nchoosek[n][i] = factorial(n) / (factorial(i) * factorial(n-i));
            
            //printf("%i\t", nchoosek[n][i]);
        }
        //printf("\n");
    }
}

int Bezier::factorial(int x)
{
    if (x == 0)
        return 1;
    else
        return(x * factorial(x-1));
}

float Bezier::bernsteinBasis(int i, int n, float t) {
    
    if ((i < 0) || (i > n))
        return 0;
    
    return nchoosek[n][i] * mju_pow(t,i) * pow(1 - t, n - i);
    
}



/* Step functions */
// (step h as f(t) / (df(t)/dt))
float Bezier::h_fun_X(float t, float x) {
    float fx =	getX(t) - x;
    float dfx = getDerivativeX(t);
    
    return fx/dfx;
}

float Bezier::h_fun_Y(float t, float y) {
    float fx =	getY(t) - y;
    float dfx = getDerivativeY(t);
    
    return fx/dfx;
}



/* Root finding using netwon_raphson method */
float Bezier::newton_raphson(float r, float t0, Bez_h_fun h_fun) {
    
    int i, imax = BezierMaxIter;
    float h, tol = BezierSolverTol;
    
    float stepsize = 1;
    
    //initial guess for t
    if (t0 < 0)
        t0 = 0;
    else if (t0 > 1)
        t0 = 1;
    
    for (i=1; i<=imax; i++)
    {
        //calculate step h as f(t) / (df(t)/dt)
        h = (this->*h_fun)(t0, r);
        
        if (mju_abs(h) < tol) {
            return t0 - h;
            printf("h: %1.8f\n", mju_abs(h));
        }
        
        //take step
        t0 -= h * stepsize;
        
        //if searched t exceeds bounds 0 to 1, wrap over and try other side
        if (t0 < 0) {
            t0 = 1; //printf("wrap under (t = %1.4f, h = %1.4f)\n",t0,h);
            stepsize /= 2;
        }
        else if (t0 > 1) {
            t0 = 0; //printf("wrap over (t = %1.4f, h = %1.4f)\n",t0,h);
            stepsize /= 2;
        }
        
    }
    
    return -1;
}




/*-----------------------------------------------------------------*/
/*                     Segmented Bezier Class                      */
/*-----------------------------------------------------------------*/

#define BezierSetMax 10			//maximum number of curves in set
#define BezierSetOrder 5		//default order
#define BezierSetCount 5		//default number of curves in set


class BezierSet {
    
public:
    BezierSet(): count(BezierSetCount), order(BezierSetOrder) {};
    BezierSet(int k): count(k), order(BezierSetOrder) {};
    BezierSet(int k, int n): count(k), order(n) {};
    
    
private:
    int count;
    int order;
    Bezier Curves[BezierSetMax];
    
    int curveIndex(Point p);
    int curveIndexFromX(float x);
    int curveIndexFromY(float y);
    
    float getYFromX(float x);
    float getXFromY(float y);
};


int BezierSet::curveIndex(Point p) {
    
    Point p0, pEnd;
    int k;
    
    for(k=0; k<count; k++) {
        p0 = Curves[k].getPoint(0);
        pEnd = Curves[k].getPoint(order);
        
        if ((p.x >= p0.x) && (p.x <= pEnd.x) &&
            (p.y >= p0.y) && (p.y <= pEnd.y))
            return k;
    }
    
    return -1;
}



int BezierSet::curveIndexFromX(float x) {
    
    Point p0, pEnd;
    int k;
    
    for(k=0; k<count; k++) {
        p0 = Curves[k].getPoint(0);
        pEnd = Curves[k].getPoint(order);
        
        if ((x >= p0.x) && (x <= pEnd.x))
            return k;
    }
    
    return -1;
}

int BezierSet::curveIndexFromY(float y) {
    
    Point p0, pEnd;
    int k;
    
    for(k=0; k<count; k++) {
        p0 = Curves[k].getPoint(0);
        pEnd = Curves[k].getPoint(order);
        
        if ((y >= p0.y) && (y <= pEnd.y))
            return k;
    }
    
    return -1;
}




float BezierSet::getYFromX(float x) {
    
    int k = curveIndexFromX(x);
    
    return Curves[k].getYFromX(x);
}

float BezierSet::getXFromY(float y) {
    
    int k = curveIndexFromY(y);
    
    return Curves[k].getXFromY(y);
}


#endif /* bezier_h */
