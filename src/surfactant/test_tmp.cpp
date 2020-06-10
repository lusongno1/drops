/*************************************************************************
	> File Name: test_tmp.cpp
	> Author: lusongno1
	> Mail: lusong@lsec.cc.ac.cn
	> Created Time: 2020年06月10日 星期三 23时28分05秒
 ************************************************************************/
#include<iostream>
#include<surfactant/sfpde.h>
using namespace std;
int main()
{
    double a[3] = {1,3,4};
    double b[3] = {5,2,0};
    double r[3];
    double r1;
    vecMinus(a,b,r);
    for(int i=0; i<3; i++)
    {
        cout<<r[i]<<endl;
    }
    crossMul(a,b,r);
    for(int i=0; i<3; i++)
    {
        cout<<r[i]<<endl;
    }
    r1 = dotP3(a,b);
    cout<<r1<<endl;
    double tet[4][3] = {{0,0,0},
        {1,0,0},{0,1,0},{0,0,1}
    };
    double i = 0;
    double x = 0;
    double y = 0;
    double z = 0.2;
    double r2;
    r2 = getBaryCoord(tet,i,x,y,z);
    cout << r2 <<endl;
    exit(0);
}
