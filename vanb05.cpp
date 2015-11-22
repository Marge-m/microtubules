#include "libminim.h"
#include <string>
#include <iostream>
#include <cstdlib>
#include <armadillo>
#include <time.h>
#include <random>
#include <sys/timeb.h>
#include <fstream>
#include <math.h> 
#include <stdexcept> 
#include <ginac/ginac.h>
#include <sstream>

using namespace std;
using namespace arma;
using namespace GiNaC;


const double k_plus = 2.0*pow(double(10), double(6)); //2
const double d_0 = 4.0;
const double g_long = -9.4;//9.4
const double g_lat = -3.2;//3.2
const double g_curl = 3.7; //3.7
const double k_hydr = 0.95;
const double E = double(17)/9*pow(double(10), double(8)); // 17/9
const double k_lat = E/1000000000*2.7*4/5.15;
const double k_long = E/1000000000 * 5.15 * 2.7/4;
const double conc = 0.00001;
const double pi = 3.14159265359;
const double kbt = 4.11 / pow(double(10), double(21));
const double k_curl = 2.0 * g_curl/(22*pi/180)/(22*pi/180);


void define_lat_energy(int r, int c, vec& right_nei_x, vec& right_nei_y, vec& left_nei_x, vec& left_nei_y, mat& b_lat, mat& s, mat& s_x1, mat& s_y1, mat& s_z1, mat& s_x2, mat& s_y2, mat& s_z2)
{
	b_lat(r, c) = 0;
	
	double current_right_x_1 = s_x1(r, c) + right_nei_x(c);
	
	double current_right_y_1 = s_y1(r, c) + right_nei_y(c);
	double current_right_z_1 = s_z1(r, c) + double(6)/13;
	double current_left_x_1 =  s_x1(r, c) + left_nei_x(c);
	double current_left_y_1 =  s_y1(r, c) + left_nei_y(c);
	double current_left_z_1 =  s_z1(r, c) - double(6)/13;
	double current_right_x_2 =  s_x2(r, c) + right_nei_x(c);
	double current_right_y_2 =  s_y2(r, c) + right_nei_y(c);
	double current_right_z_2 =  s_z2(r, c) + double(6)/13;
	double current_left_x_2 =  s_x2(r, c) + left_nei_x(c);
	double current_left_y_2 =  s_y2(r, c) + left_nei_y(c);
	double current_left_z_2 =  s_z2(r, c) - double(6)/13;
	double neighbor_x, neighbor_y, neighbor_z;
	if (c == 12)
		{			
			if (s(r + 1, 0) > 0) // верхний мономер ПФ № 0  ряд r+1
			{	
				neighbor_x = s_x2(r + 1, 0) + left_nei_x(0);
				neighbor_y = s_y2(r + 1, 0) + left_nei_y(0);
				neighbor_z = s_z2(r + 1, 0) - double(6)/13;
				double ss = sqrt(pow((current_right_x_1 - neighbor_x),double(2)) + pow((current_right_y_1 - neighbor_y), double(2)) + pow((current_right_z_1 - neighbor_z), double(2)));
				b_lat(r, c) += min(double(0), (g_lat * kbt + 0.5* k_lat * ss*ss / pow(double(10), double(18)))/2);
//				if (r > 9)
//				{
//					cout << "later "<< r <<" " << c << " " << ss << " " << b_lat(r, c) << "\n";
//				}
				
			}
			if (s(r + 2, 0) > 0) // нижний мономер ПФ № 0  ряд r+2
			{	
				neighbor_x = s_x1(r + 2, 0) + left_nei_x(0);
				neighbor_y = s_y1(r + 2, 0) + left_nei_y(0);
				neighbor_z = s_z1(r + 2, 0) - double(6)/13;
				double ss = sqrt(pow((current_right_x_2 - neighbor_x),double(2)) + pow((current_right_y_2 - neighbor_y), double(2)) + pow((current_right_z_2 - neighbor_z), double(2)));
				b_lat(r, c) += min(double(0), (g_lat * kbt + 0.5* k_lat * ss*ss / pow(double(10), double(18)))/2);
//				if (r > 9)
//								{
//									cout << "later "<< r <<" " << c << " " << ss << " " << b_lat(r, c) << "\n";
//								}
			}
			
		}
		if (c == 0) // смотрим соседей слева
		{			
			if (s(r - 2, 12) > 0) // верхний мономер ПФ № 12  ряд r-2 (нижний монососед димера ПФ номер 0)
			{	
				neighbor_x = s_x2(r - 2, 12) + right_nei_x(12);
				neighbor_y = s_y2(r - 2, 12) + right_nei_y(12);
				neighbor_z = s_z2(r - 2, 12) + double(6)/13;
				double ss = sqrt(pow((current_left_x_1 - neighbor_x),double(2)) + pow((current_left_y_1 - neighbor_y), double(2)) + pow((current_left_z_1 - neighbor_z), double(2)));
//				cout << r << " j " << r<<" " << c << b_lat(r, c) << "\n";
//				cout << "cur " <<current_left_x_1<<" " <<  neighbor_x<<" " << s_x2(r - 2, 12)<<" " <<right_nei_x(12) << " " << s_x1(r, c) <<" " << left_nei_x(c) << "\n";
//				cout << "cur " <<current_left_y_1<<" " <<  neighbor_y<<" " << s_y2(r - 2, 12)<<" " <<right_nei_y(12) << " " << s_y1(r, c) <<" " << left_nei_y(c) << "\n";
//				cout << "cur " <<current_left_z_1<<" " <<  neighbor_z<<" " << s_z2(r - 2, 12)<<" " <<double(6)/13 << " " << s_z1(r, c) <<" " << -double(6)/13 << "\n";

				b_lat(r, c) += min(double(0), (g_lat * kbt + 0.5* k_lat * ss*ss / pow(double(10), double(18)))/2);
//				if (r > 9)
//								{
//									cout << "later "<< r <<" " << c << " " << ss << " " << b_lat(r, c) << "\n";
//								}
//				if (c==0 or c==1)
//				{
//					cout << "kj1 "<<c<<"\n"; 
//				cout << "ss " << ss << " " << (g_lat * kbt + 0.5* k_lat * ss*ss / pow(double(10), double(18)))/2 << "\n";
//				cout << "s1 " <<b_lat(r, c)<< "\n";
//				}
			}
			if (s(r - 1, 12) > 0) // нижний мономер ПФ № 12  ряд r - 1 (верхний сосед)
			{	
				neighbor_x = s_x1(r - 1, 12) + right_nei_x(12);
				neighbor_y = s_y1(r - 1, 12) + right_nei_y(12);
				neighbor_z = s_z1(r - 1, 12) + double(6)/13;
				double ss = sqrt(pow((current_left_x_2 - neighbor_x),double(2)) + pow((current_left_y_2 - neighbor_y), double(2)) + pow((current_left_z_2 - neighbor_z), double(2)));
				b_lat(r, c) += min(double(0), (g_lat * kbt + 0.5* k_lat * ss*ss / pow(double(10), double(18)))/2);
//				if (r > 9)
//								{
//									cout << "later "<< r <<" " << c << " " << ss << " " << b_lat(r, c) << "\n";
//								}
//				cout << r << " j " << c << b_lat(r, c) << "\n";
//				if (c==0 or c==1)
//				{
//					cout << "kj2 "<<c<<"\n"; 
//				cout << "ss " << ss << " " << (g_lat * kbt + 0.5* k_lat * ss*ss / pow(double(10), double(18)))/2 << "\n";
//				cout << "s1 " <<b_lat(r, c)<< "\n";
//				}
			}
			
		}
		if (c >= 0 and c <= 11)//соседи справа 
		{
			if (s(r, c + 1 ) > 0) //исправила
			{	
				neighbor_x = s_x1(r, c + 1) + left_nei_x(c + 1);
				neighbor_y = s_y1(r, c + 1) + left_nei_y(c + 1);
				neighbor_z = s_z1(r, c + 1) - double(6)/13;
				double ss = sqrt(pow((current_right_x_1 - neighbor_x),double(2)) + pow((current_right_y_1 - neighbor_y), double(2)) + pow((current_right_z_1 - neighbor_z), double(2)));
				b_lat(r, c) += min(double(0), (g_lat * kbt + 0.5* k_lat * ss*ss / pow(double(10), double(18)))/2);
//				if (r > 9)
//								{
//									cout << "later "<< r <<" " << c << " " << ss << " " << b_lat(r, c) << "\n";
//								}
				neighbor_x = s_x2(r, c + 1) + left_nei_x(c + 1);
				neighbor_y = s_y2(r, c + 1) + left_nei_y(c + 1);
				neighbor_z = s_z2(r, c + 1) - double(6)/13;
				ss = sqrt(pow((current_right_x_2 - neighbor_x),double(2)) + pow((current_right_y_2 - neighbor_y), double(2)) + pow((current_right_z_2 - neighbor_z), double(2)));
				b_lat(r, c) += min(double(0), (g_lat * kbt + 0.5* k_lat * ss*ss / pow(double(10), double(18)))/2);
//				if (r > 9)
//								{
//									cout << "later "<< r <<" " << c << " " << ss << " " << b_lat(r, c) << "\n";
//								}
//				if (c==0 or c==1)
//				{
//					cout << "kj3 "<<c<<"\n"; 
//				cout << "ss " << ss << " " << (g_lat * kbt + 0.5* k_lat * ss*ss / pow(double(10), double(18)))/2 << "\n";
//				cout << "s1 " <<b_lat(r, c)<< "\n";
//				}
			}

		}
		if (c >=1 and c <= 12)//соседи слева
		{
			if (s(r, c - 1) > 0) 
			{	
				
				neighbor_x = s_x1(r, c - 1) + right_nei_x(c - 1);
				neighbor_y = s_y1(r, c - 1) + right_nei_y(c - 1);
				neighbor_z = s_z1(r, c - 1) + double(6)/13;
				double ss = sqrt(pow((current_left_x_1 - neighbor_x),double(2)) + pow((current_left_y_1 - neighbor_y), double(2)) + pow((current_left_z_1 - neighbor_z), double(2)));
				b_lat(r, c) += min(double(0), (g_lat * kbt + 0.5* k_lat * ss*ss / pow(double(10), double(18)))/2);
//				if (r > 9)
//								{
//									cout << "later "<< r <<" " << c << " " << ss << " " << b_lat(r, c) << "\n";
//								}
				neighbor_x = s_x2(r, c - 1) + right_nei_x(c - 1);
				neighbor_y = s_y2(r, c - 1) + right_nei_y(c - 1);
				neighbor_z = s_z2(r, c - 1) + double(6)/13;
				ss = sqrt(pow((current_left_x_2 - neighbor_x),double(2)) + pow((current_left_y_2 - neighbor_y), double(2)) + pow((current_left_z_2 - neighbor_z), double(2)));
				b_lat(r, c) += min(double(0), (g_lat * kbt + 0.5* k_lat * ss*ss / pow(double(10), double(18)))/2);
//				if (r > 9)
//								{
//									cout << "later "<< r <<" " << c << " " << ss << " " << b_lat(r, c) << "\n";
//								}
//				if (c==0 or c==1)
//				{
//					cout << "kj4 "<<c<<"\n"; 
//				cout << "ss " << ss << " " << (g_lat * kbt + 0.5* k_lat * ss*ss / pow(double(10), double(18)))/2 << "\n";
//				cout << "s1 " <<b_lat(r, c)<< "\n";
//				}
			}
		}
}

void define_all_energies(int r, int c, vec& right_nei_x, vec& right_nei_y, vec& left_nei_x, vec& left_nei_y, mat& b_lat, mat& b_curl, mat& b_long, mat& s_tetta, mat& s_fi, mat& s_d, mat& s_x1, mat& s_y1, mat& s_z1, mat& s_x2, mat& s_y2, mat& s_z2, mat& s)
{
	
	    b_long(r, c) = 0;
		b_long(r, c) += min(double(0), g_long*kbt + 0.5*k_long*s_d(r, c)*s_d(r, c) / pow(double(10), double(18)));
		double l_p, m_p, n_p;
		if (s(r, c) == 2)
		{
			l_p = 0;
			m_p = 0;
			n_p = 1;
			
		}
		if (s(r, c) == 1)
		{
			l_p = sin(22*pi/180);
			m_p = 0;
			n_p = cos(22*pi/180);
			
		}
//		if (s(r, c) == 1)
//			{
//				cout << "klo "<<r<< " " << c<< " " << l_p<< "\n";
//			}
		double l, m, n;
		l = sin(s_tetta(r, c))*cos(s_fi(r, c));
		m = sin(s_tetta(r, c))*sin(s_fi(r, c));
		n = cos(s_tetta(r, c));
		
		double cos_f = l*l_p + m*m_p + n*n_p;
//		if (s(r, c) == 1)
//				{
//					cout << "klo2 "<< l <<" " << c<< " " << l_p<<" " << cos_f<< "\n";
//				}
		b_curl(r, c) = double(1)/2*k_curl*kbt*acos(cos_f)*acos(cos_f);
//		if (s(r, c) == 1)
//				{
//					cout << "klo3 "<< acos(cos_f)*acos(cos_f) <<" " <<  1/2*k_curl*kbt*acos(cos_f)*acos(cos_f)<< " " << b_curl(r, c)<< "\n";
//				}
		double current_right_x_1 = s_x1(r, c) + right_nei_x(c);
		double current_right_y_1 = s_y1(r, c) + right_nei_y(c);
	    double current_right_z_1 = s_z1(r, c) + double(6)/13;
		double current_left_x_1 =  s_x1(r, c) + left_nei_x(c);
		double current_left_y_1 =  s_y1(r, c) + left_nei_y(c);
		double current_left_z_1 =  s_z1(r, c) - double(6)/13;
		double current_right_x_2 =  s_x2(r, c) + right_nei_x(c);
		double current_right_y_2 =  s_y2(r, c) + right_nei_y(c);
		double current_right_z_2 =  s_z2(r, c) + double(6)/13;
		double current_left_x_2 =  s_x2(r, c) + left_nei_x(c);
		double current_left_y_2 =  s_y2(r, c) + left_nei_y(c);
		double current_left_z_2 =  s_z2(r, c) - double(6)/13;
		b_lat(r, c) = 0;
		double neighbor_x_3, neighbor_y_3, neighbor_z_3, neighbor_x_1, neighbor_y_1, neighbor_z_1, neighbor_x_4, neighbor_y_4, neighbor_z_4, neighbor_x_2, neighbor_y_2, neighbor_z_2;
		if (c == 12)
		{			
			if (s(r + 1, 0) > 0) // верхний мономер ПФ № 0  ряд r+1
			{	
				
				neighbor_x_2 = s_x2(r + 1, 0) + left_nei_x(0);
				neighbor_y_2 = s_y2(r + 1, 0) + left_nei_y(0);
				neighbor_z_2 = s_z2(r + 1, 0) - double(6)/13;
				double ss = sqrt(pow((current_right_x_1 - neighbor_x_2),double(2)) + pow((current_right_y_1 - neighbor_y_2), double(2)) + pow((current_right_z_1 - neighbor_z_2), double(2)));
				
				b_lat(r, c) += min(double(0), (g_lat * kbt + 0.5* k_lat * ss*ss / pow(double(10), double(18)))/2);
//				if (r > 9)
//				{
//					cout << "define "<< r <<" " << c << " " << ss << " " << b_lat(r, c) << "\n";
//				}
//				cout << "ss " << ss << " " << (g_lat * kbt + 0.5* k_lat * ss*ss / pow(double(10), double(18)))/2 << "\n";
			}
			if (s(r + 2, 0) > 0) // нижний мономер ПФ № 0  ряд r+2
			{	
				neighbor_x_4 = s_x1(r + 2, 0) + left_nei_x(0);
				neighbor_y_4 = s_y1(r + 2, 0) + left_nei_y(0);
				neighbor_z_4 = s_z1(r + 2, 0) - double(6)/13;
				double ss = sqrt(pow((current_right_x_2 - neighbor_x_4), double(2)) + pow((current_right_y_2 - neighbor_y_4), double(2)) + pow((current_right_z_2 - neighbor_z_4), double(2)));
				b_lat(r, c) += min(double(0), (g_lat * kbt + 0.5* k_lat * ss*ss / pow(double(10), double(18)))/2);
//				if (r > 9)
//				{
//					cout << "define "<< r <<" " << c << " " << ss << " " << b_lat(r, c) << "\n";
//				}
//				if (s(r, c) == 1)
//								{
//									cout << "kla2 "<< current_right_x_2 <<" " <<  current_right_y_2 << " " << current_right_z_2 <<" " <<ss<< " " << (g_lat * kbt + 0.5* k_lat * ss*ss / pow(double(10), double(18)))/2 <<" " << b_lat(r, c) << "\n";
//								}
//				cout << "ss " << ss << " " << (g_lat * kbt + 0.5* k_lat * ss*ss / pow(double(10), double(18)))/2 << "\n";
			}
			
		}
		if (c == 0) // смотрим соседей слева
		{			
			if (s(r - 2, 12) > 0) // верхний мономер ПФ № 12  ряд r-2 (нижний монососед димера ПФ номер 0)
			{	
				neighbor_x_1 = s_x2(r - 2, 12) + right_nei_x(12);
				neighbor_y_1 = s_y2(r - 2, 12) + right_nei_y(12);
				neighbor_z_1 = s_z2(r - 2, 12) + double(6)/13;
				double ss = sqrt(pow((current_left_x_1 - neighbor_x_1), double(2)) + pow((current_left_y_1 - neighbor_y_1), double(2)) + pow((current_left_z_1 - neighbor_z_1), double(2)) );
				b_lat(r, c) += min(double(0), (g_lat * kbt + 0.5* k_lat * ss*ss / pow(double(10), double(18)))/2);
//				if (r > 9)
//								{
//									cout << "define "<< r <<" " << c << " " << ss << " " << b_lat(r, c) << "\n";
//								}
//				if (s(r, c) == 1)
//								{
//									cout << "kla3 "<< current_left_x_1 <<" " <<  current_left_y_1 << " " << current_left_z_1 <<" " <<ss<< " " << (g_lat * kbt + 0.5* k_lat * ss*ss / pow(double(10), double(18)))/2 << " " << b_lat(r, c) <<"\n";
//								}
//				if (c==0 or c==1)
//				{
//					cout << "kj5 "<<c<<"\n"; 
//				cout << "ss " << ss << " " << (g_lat * kbt + 0.5* k_lat * ss*ss / pow(double(10), double(18)))/2 << "\n";
//				cout << "s1 " <<b_lat(r, c)<< "\n";
//				}
			}
			if (s(r - 1, 12) > 0) // нижний мономер ПФ № 12  ряд r - 1 (верхний сосед)
			{	
				neighbor_x_3 = s_x1(r - 1, 12) + right_nei_x(12);
				neighbor_y_3 = s_y1(r - 1, 12) + right_nei_y(12);
				neighbor_z_3 = s_z1(r - 1, 12) + double(6)/13;
				double ss = sqrt(pow((current_left_x_2 - neighbor_x_3), double(2)) + pow((current_left_y_2 - neighbor_y_3), double(2)) + pow((current_left_z_2 - neighbor_z_3), double(2)) );
				b_lat(r, c) += min(double(0), (g_lat * kbt + 0.5* k_lat * ss*ss / pow(double(10), double(18)))/2);
//				if (r > 9)
//								{
//									cout << "define "<< r <<" " << c << " " << ss << " " << b_lat(r, c) << "\n";
//								}
//				if (s(r, c) == 1)
//								{
//									cout << "kla4 "<< current_left_x_2 <<" " <<  current_left_y_2 << " " << current_left_z_2 <<" " <<ss<< " " << (g_lat * kbt + 0.5* k_lat * ss*ss / pow(double(10), double(18)))/2 <<" " << b_lat(r, c) << "\n";
//								}
//				cout << "ss " << ss << " " << (g_lat * kbt + 0.5* k_lat * ss*ss / pow(double(10), double(18)))/2 << "\n";
//				cout << "s2 "<< r<< " " << c << "\n";
//				if (c==0 or c==1)
//				{
//					cout << "kj6 "<<c<<"\n"; 
//					cout << "ss " << ss << " " << (g_lat * kbt + 0.5* k_lat * ss*ss / pow(double(10), double(18)))/2 << "\n";
//					cout << "s1 " <<b_lat(r, c)<< "\n";
//				}
			}
			
		}
		if (c >= 0 and c <= 11)//соседи справа
		{
			if (s(r, c + 1 ) > 0) 
			{	
				neighbor_x_2 = s_x1(r, c + 1) + left_nei_x(c + 1);
				neighbor_y_2 = s_y1(r, c + 1) + left_nei_y(c + 1);
				neighbor_z_2 = s_z1(r, c + 1) - double(6)/13;
				double ss = sqrt(pow((current_right_x_1 - neighbor_x_2), double(2)) + pow((current_right_y_1 - neighbor_y_2), double(2)) + pow((current_right_z_1 - neighbor_z_2), double(2)) );
				b_lat(r, c) += min(double(0), (g_lat * kbt + 0.5* k_lat * ss*ss /pow(double(10), double(18)))/2);
//				if (r > 9)
//								{
//									cout << "define "<< r <<" " << c << " " << ss << " " << b_lat(r, c) << "\n";
//									cout << "kla5 "<< current_right_x_1 <<" " <<  current_right_y_1 << " " << current_right_z_1 <<" " <<ss<< " " << neighbor_x_2 << " " << neighbor_y_2<<" " << neighbor_z_2 << " " << b_lat(r, c) << "\n";
//								}
//				if (s(r, c) == 1)
//								{
//								}
				neighbor_x_4 = s_x2(r, c + 1) + left_nei_x(c + 1);
				neighbor_y_4 = s_y2(r, c + 1) + left_nei_y(c + 1);
				neighbor_z_4 = s_z2(r, c + 1) - double(6)/13;
				ss = sqrt(pow((current_right_x_2 - neighbor_x_4), double(2)) + pow((current_right_y_2 - neighbor_y_4), double(2)) + pow((current_right_z_2 - neighbor_z_4), double(2)) );
				b_lat(r, c) += min(double(0), (g_lat * kbt + 0.5* k_lat * ss*ss / pow(double(10), double(18)))/2);
//				if (r > 9)
//								{
//									cout << "define "<< r <<" " << c << " " << ss << " " << b_lat(r, c) << "\n";
//									cout << "kla6 "<< current_right_x_2 <<" " <<  current_right_y_2 << " " << current_right_z_2 <<" " <<ss<< " " << neighbor_x_4 << " " << neighbor_y_4<<" " << neighbor_z_4 << " " << b_lat(r, c) << "\n";
//								}
//				if (s(r, c) == 1)
//								{
//									cout << "kla6 "<< current_right_x_2 <<" " <<  current_right_y_2 << " " << current_right_z_2 <<" " <<ss<< " " << (g_lat * kbt + 0.5* k_lat * ss*ss / pow(double(10), double(18)))/2 << " " << b_lat(r, c) <<"\n";
//								}
//				cout << "ss " << ss << " " << (g_lat * kbt + 0.5* k_lat * ss*ss / pow(double(10), double(18)))/2 << "\n";
//				if (c==0 or c==1)
//				{
//					cout << "kj7 "<<c<<"\n"; 
//					cout << "ss " << ss << " " << (g_lat * kbt + 0.5* k_lat * ss*ss / pow(double(10), double(18)))/2 << "\n";
//					cout << "s1 " <<b_lat(r, c)<< "\n";
//				}
			}

		}
		if (c >=1 and c <= 12)//соседи слева
		{
			if (s(r, c - 1) > 0) 
			{	
				neighbor_x_1 = s_x1(r, c - 1) + right_nei_x(c - 1);
				neighbor_y_1 = s_y1(r, c - 1) + right_nei_y(c - 1);
				neighbor_z_1 = s_z1(r, c - 1) + double(6)/13;
				double ss = sqrt(pow((current_left_x_1 - neighbor_x_1), double(2)) + pow((current_left_y_1 - neighbor_y_1), double(2)) + pow((current_left_z_1 - neighbor_z_1), double(2)) );
				b_lat(r, c) += min(double(0), (g_lat * kbt + 0.5* k_lat * ss*ss / pow(double(10), double(18)))/2);
//				if (r > 9)
//								{
//									cout << "define "<< r <<" " << c << " " << ss << " " << b_lat(r, c) << "\n";
//									cout << "kla7 "<< current_right_x_1 <<" " <<  current_right_y_1 << " " << current_right_z_1 <<" " <<ss<< " " << neighbor_x_1 << " " << neighbor_y_1<<" " << neighbor_z_1 << " " << b_lat(r, c) << "\n";
//
//								}
//				if (s(r, c) == 1)
//								{
//								}
//				if (c==0 or c==1)
//				{
//					cout << "kj81 "<<c<<"\n"; 
//				cout << "ss " << ss << " " << (g_lat * kbt + 0.5* k_lat * ss*ss / pow(double(10), double(18)))/2 << "\n";
//				cout << "s1 " <<b_lat(r, c)<< "\n";
//				}
				neighbor_x_3 = s_x2(r, c - 1) + right_nei_x(c - 1);
				neighbor_y_3 = s_y2(r, c - 1) + right_nei_y(c - 1);
				neighbor_z_3 = s_z2(r, c - 1) + double(6)/13;
				ss = sqrt(pow((current_left_x_2 - neighbor_x_3), double(2)) + pow((current_left_y_2 - neighbor_y_3), double(2)) + pow((current_left_z_2 - neighbor_z_3), double(2)) );
				b_lat(r, c) += min(double(0), (g_lat * kbt + 0.5* k_lat * ss*ss / pow(double(10), double(18)))/2);
//				if (r > 9)
//								{
//									cout << "define "<< r <<" " << c << " " << ss << " " << b_lat(r, c) << "\n";
//									cout << "kla8 "<< current_right_x_2 <<" " <<  current_right_y_2 << " " << current_right_z_2 <<" " <<ss<< " " << neighbor_x_3 << " " << neighbor_y_3<<" " << neighbor_z_3 << " " << b_lat(r, c) << "\n";
//								}
//				if (s(r, c) == 1)
//								{
//								}
//				cout << "ss " << ss << " " << (g_lat * kbt + 0.5* k_lat * ss*ss / pow(double(10), double(18)))/2 << "\n";
//				if (c==0 or c==1)
//				{
//					cout << "kj8 "<<c<<"\n"; 
//				cout << "ss " << ss << " " << (g_lat * kbt + 0.5* k_lat * ss*ss / pow(double(10), double(18)))/2 << "\n";
//				cout << "s1 " <<b_lat(r, c)<< "\n";
//				}
			}
		}
		
}




void minimize_energy(int r, int c, vec& right_nei_x, vec& right_nei_y, vec& left_nei_x, vec& left_nei_y, mat& b_lat, mat& b_curl, mat& b_long, mat& s_tetta, mat& s_fi, mat& s_d, mat& s_x1, mat& s_y1, mat& s_z1, mat& s_x2, mat& s_y2, mat& s_z2, mat& s, vec& leng)
{
	symbol current_tetta("x(1)"), current_fi("x(2)"), new_d("x(3)");
	ex curl_energy;
	ex lat_energy;
	ex long_energy;
	ex l = sin(current_tetta)*cos(current_fi);
	ex m = sin(current_tetta)*sin(current_fi);
	ex n = cos(current_tetta);
	
	
	double l_p, m_p, n_p;
	if (s(r, c) == 2)
	{
		l_p = 0;
		m_p = 0;
		n_p = 1;
		
	}
	if (s(r, c) == 1)
	{
		l_p = sin(22*pi/180);
		m_p = 0;
		n_p = cos(22*pi/180);
//		cout << "uio\n";
		
	}
	ex cos_f = l*l_p + m*m_p + n*n_p;
	curl_energy = double(1)/2*k_curl*kbt*acos(cos_f)*acos(cos_f);
	std::ostringstream buf;
	buf.clear();
	buf << curl_energy << " + ";
//	cout << "buf " << curl_energy << "\n";
//	total_energy += curl_energy;
	long_energy = g_long*kbt + 0.5*k_long*new_d*new_d / pow(double(10), double(18)); //размерность
	buf << "min(0.0, " << long_energy << ")";
//	ex min_long = minimal_dim(double(0), long_energy);
//	total_energy += min_long;
	double total_tetta = 0;
	double total_fi = 0;
	for (int ro = 10; ro < r; ro++)
	{
		total_tetta += s_tetta(ro, c);
		total_fi += s_fi(ro, c);
	} //суммарный угол точки, к которой крепится рассматриваемый вектор
	total_fi -= 2*pi/13*c; // посчитали абсолютные углы
	ex current_x_1 = s_x2(r - 1, c) + (new_d + 4)*sin(total_tetta + current_tetta)*cos(total_fi + current_fi);
	ex current_y_1 = s_y2(r - 1, c) + (new_d + 4)*sin(total_tetta + current_tetta)*sin(total_fi + current_fi);	
	ex current_z_1 = s_z2(r - 1, c) + (new_d + 4)*cos(total_tetta + current_tetta);
//	double curr_x_1_a = s_x2(r - 1, c);
//	double curr_x_1_b = right_nei_x(c);
//	double curr_y_1_a = s_y2(r - 1, c);
//	double curr_y_1_b = right_nei_y(c);
//	double curr_z_1_a = s_z2(r - 1, c);
//	double curr_x_1_c = left_nei_x(c);
//	double curr_y_1_c = left_nei_y(c);
	ex current_right_x_1 = current_x_1 + right_nei_x(c);
	ex current_right_y_1 = current_y_1 + right_nei_y(c);
	ex current_right_z_1 = current_z_1 + double(6)/13;
	ex current_left_x_1 = current_x_1 + left_nei_x(c);
	ex current_left_y_1 = current_y_1 + left_nei_y(c);
	ex current_left_z_1 = current_z_1 - double(6)/13;
	ex current_x_2 = s_x2(r - 1, c) + (new_d + 8)*sin(total_tetta + current_tetta)*cos(total_fi + current_fi);
	ex current_y_2 = s_y2(r - 1, c) + (new_d + 8)*sin(total_tetta + current_tetta)*sin(total_fi + current_fi);
	ex current_z_2 = s_z2(r - 1, c) + (new_d + 8)*cos(total_tetta + current_tetta);
	ex current_right_x_2 = current_x_2 + right_nei_x(c);
	ex current_right_y_2 = current_y_2 + right_nei_y(c);
	ex current_right_z_2 = current_z_2 + double(6)/13;
	ex current_left_x_2 = current_x_2 + left_nei_x(c);
	ex current_left_y_2 = current_y_2 + left_nei_y(c);
	ex current_left_z_2 = current_z_2 - double(6)/13;
	ex neighbor_x_3, neighbor_y_3, neighbor_z_3, neighbor_x_1, neighbor_y_1, neighbor_z_1, neighbor_x_4, neighbor_y_4, neighbor_z_4, neighbor_x_2, neighbor_y_2, neighbor_z_2;
	neighbor_x_1 = 0;
	neighbor_x_2 = 0;
	neighbor_x_3 = 0;
	neighbor_x_4 = 0;
	neighbor_y_1 = 0;
	neighbor_y_2 = 0;
	neighbor_y_3 = 0;
	neighbor_y_4 = 0;
	neighbor_z_1 = 0;
	neighbor_z_2 = 0;
	neighbor_z_3 = 0;
	neighbor_z_4 = 0;
//	double lat_1, lat_2, lat_3, lat_4;
//	lat_1 = 0;
//	lat_2 = 0;
//	lat_3 = 0;
//	lat_4 = 0;
	if (c == 12)
	{			
		if (s(r + 1, 0) > 0) // верхний мономер ПФ № 0  ряд r+1
		{	
//			lat_2 = 1;
			neighbor_x_2 = s_x2(r + 1, 0) + left_nei_x(0);
			neighbor_y_2 = s_y2(r + 1, 0) + left_nei_y(0);
			neighbor_z_2 = s_z2(r + 1, 0) - double(6)/13;
			ex s = sqrt((current_right_x_1 - neighbor_x_2)*(current_right_x_1 - neighbor_x_2) + (current_right_y_1 - neighbor_y_2)*(current_right_y_1 - neighbor_y_2) + (current_right_z_1 - neighbor_z_2)*(current_right_z_1 - neighbor_z_2));
//			total_energy += minimal_dim(double(0), (g_lat * kbt + 0.5* k_lat * s*s / pow(double(10), double(18)))/2);
			lat_energy = (g_lat * kbt + 0.5* k_lat * s*s / pow(double(10), double(18)))/2;
			buf << " + min(0.0, " << lat_energy << ")";
//			if (s(r, c) == 1)
//											{
//												cout << "kla1 "<< r<< " " << c << " " << current_right_x_1 <<" " <<  current_right_y_1 << " " << current_right_z_1 <<" " <<lat_energy<< "\n";
//											}
//			cout << "nei " << neighbor_x_2 <<"\n";
//									cout << "snei "<< s << "\n";
//									cout << "lat " << lat_energy << "\n";
//									buf << " + min(0.0, " << lat_energy << ")";
		}
		if (s(r + 2, 0) > 0) // нижний мономер ПФ № 0  ряд r+2
		{	
//			lat_4 = 1;
			neighbor_x_4 = s_x1(r + 2, 0) + left_nei_x(0);
			neighbor_y_4 = s_y1(r + 2, 0) + left_nei_y(0);
			neighbor_z_4 = s_z1(r + 2, 0) - double(6)/13;
			ex s = sqrt((current_right_x_2 - neighbor_x_4)*(current_right_x_2 - neighbor_x_4) + (current_right_y_2 - neighbor_y_4)* (current_right_y_2 - neighbor_y_4) + (current_right_z_2 - neighbor_z_4)*(current_right_z_2 - neighbor_z_4));
//			total_energy += minimal_dim(double(0), (g_lat * kbt + 0.5* k_lat * s*s / pow(double(10), double(18)))/2);
			lat_energy = (g_lat * kbt + 0.5* k_lat * s*s / pow(double(10), double(18)))/2;
			buf << " + min(0.0, " << lat_energy << ")";
//			if (s(r, c) == 1)
//						{
//							cout << "kla2 "<< r<< " " << c << " " <<lat_energy<< "\n";
//						}
//			cout << "nei " << neighbor_x_2 <<"\n";
//									cout << "snei "<< s << "\n";
//									cout << "lat " << lat_energy << "\n";
//									buf << " + min(0.0, " << lat_energy << ")";
		}
		
	}
	if (c == 0) // смотрим соседей слева
	{			
		if (s(r - 2, 12) > 0) // верхний мономер ПФ № 12  ряд r-2 (нижний монососед димера ПФ номер 0)
		{	
//			lat_1 = 1;
			neighbor_x_1 = s_x2(r - 2, 12) + right_nei_x(12);
			neighbor_y_1 = s_y2(r - 2, 12) + right_nei_y(12);
			neighbor_z_1 = s_z2(r - 2, 12) + double(6)/13;
			ex s = sqrt((current_left_x_1 - neighbor_x_1)*(current_left_x_1 - neighbor_x_1) + (current_left_y_1 - neighbor_y_1)*(current_left_y_1 - neighbor_y_1) + (current_left_z_1 - neighbor_z_1)*(current_left_z_1 - neighbor_z_1));
//			total_energy += minimal_dim(double(0), (g_lat * kbt + 0.5* k_lat * s*s / pow(double(10), double(18)))/2);
			lat_energy = (g_lat * kbt + 0.5* k_lat * s*s / pow(double(10), double(18)))/2;
			buf << " + min(0.0, " << lat_energy << ")";
//			if (s(r, c) == 1)
//				{
//					cout << "kla3 "<< r<< " " << c << " " <<lat_energy<< "\n";
//				}
//			cout << "nei " << current_left_x_1 << " " << neighbor_x_1 <<" " << current_left_y_1 << " " <<neighbor_y_1 << " "<<current_left_z_1 <<neighbor_z_1 << "\n";
//									cout << "snei "<< s << "\n";
//									cout << "lat " << lat_energy << "\n";
//									buf << " + min(0.0, " << lat_energy << ")";
		}
		if (s(r - 1, 12) > 0) // нижний мономер ПФ № 12  ряд r - 1 (верхний сосед)
		{	
//			lat_3 = 1;
			neighbor_x_3 = s_x1(r - 1, 12) + right_nei_x(12);
			neighbor_y_3 = s_y1(r - 1, 12) + right_nei_y(12);
			neighbor_z_3 = s_z1(r - 1, 12) + double(6)/13;
			ex s = sqrt((current_left_x_2 - neighbor_x_3)*(current_left_x_2 - neighbor_x_3) + (current_left_y_2 - neighbor_y_3)*(current_left_y_2 - neighbor_y_3) + (current_left_z_2 - neighbor_z_3)*(current_left_z_2 - neighbor_z_3) );
//			total_energy += minimal_dim(double(0), (g_lat * kbt + 0.5* k_lat * s*s / pow(double(10), double(18)))/2);
			lat_energy = (g_lat * kbt + 0.5* k_lat * s*s / pow(double(10), double(18)))/2;
			buf << " + min(0.0, " << lat_energy << ")";
//			if (s(r, c) == 1)
//						{
//							cout << "kla4 "<< r<< " " << c << " " <<lat_energy<< "\n";
//						}
//			cout << "nei " << current_left_x_2 << " " << neighbor_x_3<<" " << current_left_y_2 << " " <<neighbor_y_3 << " "<<current_left_z_2 <<neighbor_z_3 << "\n";
//									cout << "snei "<< s << "\n";
//									cout << "lat " << lat_energy << "\n";
//									buf << " + min(0.0, " << lat_energy << ")";
		}
		
	}
	if (c >= 0 and c <= 11)//соседи справа
	{
		if (s(r, c + 1 ) > 0) 
		{	
//			lat_4 = 1;
//			lat_2 = 1;
			neighbor_x_2 = s_x1(r, c + 1) + left_nei_x(c + 1);
			neighbor_y_2 = s_y1(r, c + 1) + left_nei_y(c + 1);
			neighbor_z_2 = s_z1(r, c + 1) - double(6)/13;
			ex s = sqrt((current_right_x_1 - neighbor_x_2)*(current_right_x_1 - neighbor_x_2) + (current_right_y_1 - neighbor_y_2)*(current_right_y_1 - neighbor_y_2) + (current_right_z_1 - neighbor_z_2)*(current_right_z_1 - neighbor_z_2));
//			total_energy += minimal_dim(double(0), (g_lat * kbt + 0.5* k_lat * s*s / pow(double(10), double(18)))/2);
			lat_energy = (g_lat * kbt + 0.5* k_lat * s*s / pow(double(10), double(18)))/2;
			buf << " + min(0.0, " << lat_energy << ")";
//			if (s(r, c) == 1)
//									{
//										cout << "kla5 "<< r<< " " << c << " " <<lat_energy<< "\n";
//									}
//			cout << "nei " << neighbor_x_2 <<"\n";
//									cout << "snei "<< s << "\n";
//									cout << "lat " << lat_energy << "\n";
//									buf << " + min(0.0, " << lat_energy << ")";
			neighbor_x_4 = s_x2(r, c + 1) + left_nei_x(c + 1);
			neighbor_y_4 = s_y2(r, c + 1) + left_nei_y(c + 1);
			neighbor_z_4 = s_z2(r, c + 1) - double(6)/13;
			s = sqrt((current_right_x_2 - neighbor_x_4)*(current_right_x_2 - neighbor_x_4) + (current_right_y_2 - neighbor_y_4)*(current_right_y_2 - neighbor_y_4) + (current_right_z_2 - neighbor_z_4)*(current_right_z_2 - neighbor_z_4) );
//			total_energy += minimal_dim(double(0), (g_lat * kbt + 0.5* k_lat * s*s / pow(double(10), double(18)))/2);
			lat_energy = (g_lat * kbt + 0.5* k_lat * s*s / pow(double(10), double(18)))/2;
			buf << " + min(0.0, " << lat_energy << ")";
		}

	}
	if (c >=1 and c <= 12)//соседи слева
	{
		if (s(r, c - 1) > 0) 
		{	
//			lat_1 = 1;
//			lat_3 = 1;
			neighbor_x_1 = s_x1(r, c - 1) + right_nei_x(c - 1);
			neighbor_y_1 = s_y1(r, c - 1) + right_nei_y(c - 1);
			neighbor_z_1 = s_z1(r, c - 1) + double(6)/13;
			ex s = sqrt((current_left_x_1 - neighbor_x_1)*(current_left_x_1 - neighbor_x_1) + (current_left_y_1 - neighbor_y_1)*(current_left_y_1 - neighbor_y_1) + (current_left_z_1 - neighbor_z_1)*(current_left_z_1 - neighbor_z_1) );
//			total_energy += minimal_dim(double(0), (g_lat * kbt + 0.5* k_lat * s*s / pow(double(10), double(18)))/2);
			lat_energy = (g_lat * kbt + 0.5* k_lat * s*s / pow(double(10), double(18)))/2;
			buf << " + min(0.0, " << lat_energy << ")";	
//			if (s(r, c) == 1)
//									{
//										cout << "kla2 "<< r<< " " << c << " " <<lat_energy<< "\n";
//									}
//			cout << "nei " << neighbor_x_2 <<"\n";
//									cout << "snei "<< s << "\n";
//									cout << "lat " << lat_energy << "\n";
//									buf << " + min(0.0, " << lat_energy << ")";
			neighbor_x_3 = s_x2(r, c - 1) + right_nei_x(c - 1);
			neighbor_y_3 = s_y2(r, c - 1) + right_nei_y(c - 1);
			neighbor_z_3 = s_z2(r, c - 1) + double(6)/13;
			s = sqrt((current_left_x_2 - neighbor_x_3)*(current_left_x_2 - neighbor_x_3) + (current_left_y_2 - neighbor_y_3)*(current_left_y_2 - neighbor_y_3) + (current_left_z_2 - neighbor_z_3)*(current_left_z_2 - neighbor_z_3) );
//			total_energy += minimal_dim(double(0), (g_lat * kbt + 0.5* k_lat * s*s / pow(double(10), double(18)))/2);
			lat_energy = (g_lat * kbt + 0.5* k_lat * s*s / pow(double(10), double(18)))/2;
			buf << " + min(0.0, " << lat_energy << ")";
		}
	}
	double add_to_total_tetta = 0;
	double add_to_total_fi = 0;
	for (int i = r + 1; i < leng(c) + 1; i++)
	{
		if (s(i, c) > 0)
		{
			ex new_current_x_1 = current_x_2 + (s_d(i, c) + 4)*sin(total_tetta + add_to_total_tetta + current_tetta + s_tetta(i, c))*cos(total_fi + add_to_total_fi + current_fi + s_fi(i, c));
			ex new_current_y_1 = current_y_2 + (s_d(i, c) + 4)*sin(total_tetta + add_to_total_tetta + current_tetta + s_tetta(i, c))*sin(total_fi + add_to_total_fi + current_fi +  s_fi(i, c));	
			ex new_current_z_1 = current_z_2 + (s_d(i, c) + 4)*cos(total_tetta + add_to_total_tetta + current_tetta + s_tetta(i, c));
			ex new_current_x_2 = current_x_2 + (s_d(i, c) + 8)*sin(total_tetta + add_to_total_tetta + current_tetta + s_tetta(i, c))*cos(total_fi + add_to_total_fi + current_fi + s_fi(i, c));
			ex new_current_y_2 = current_y_2 + (s_d(i, c) + 8)*sin(total_tetta + add_to_total_tetta + current_tetta + s_tetta(i, c))*sin(total_fi + add_to_total_fi + current_fi +  s_fi(i, c));	
			ex new_current_z_2 = current_z_2 + (s_d(i, c) + 8)*cos(total_tetta + add_to_total_tetta + current_tetta + s_tetta(i, c));
			ex new_current_right_x_1 = new_current_x_1 + right_nei_x(c);
			ex new_current_right_y_1 = new_current_y_1 + right_nei_y(c);
			ex new_current_right_z_1 = new_current_z_1 + double(6)/13;
			ex new_current_left_x_1 = new_current_x_1 + left_nei_x(c);
			ex new_current_left_y_1 = new_current_y_1 + left_nei_y(c);
			ex new_current_left_z_1 = new_current_z_1 - double(6)/13;
			ex new_current_right_x_2 = new_current_x_2 + right_nei_x(c);
			ex new_current_right_y_2 = new_current_y_2 + right_nei_y(c);
			ex new_current_right_z_2 = new_current_z_2 + double(6)/13;
			ex new_current_left_x_2 = new_current_x_2 + left_nei_x(c);
			ex new_current_left_y_2 = new_current_y_2 + left_nei_y(c);
			ex new_current_left_z_2 = new_current_z_2 - double(6)/13;
			if (c == 12)
				{			
					if (s(i + 1, 0) > 0) // верхний мономер ПФ № 0  ряд r+1
					{	
			//			lat_2 = 1;
						neighbor_x_2 = s_x2(i + 1, 0) + left_nei_x(0);
						neighbor_y_2 = s_y2(i + 1, 0) + left_nei_y(0);
						neighbor_z_2 = s_z2(i + 1, 0) - double(6)/13;
//						cout << "nei " << neighbor_x_2 <<"\n";
						ex s = sqrt((new_current_right_x_1 - neighbor_x_2)*(new_current_right_x_1 - neighbor_x_2) + (new_current_right_y_1 - neighbor_y_2)*(new_current_right_y_1 - neighbor_y_2) + (new_current_right_z_1 - neighbor_z_2)*(new_current_right_z_1 - neighbor_z_2));
//						cout << "snei "<< s << "\n";
//						total_energy += minimal_dim(double(0), (g_lat * kbt + 0.5* k_lat * s*s / pow(double(10), double(18)))/2);
						lat_energy = (g_lat * kbt + 0.5* k_lat * s*s / pow(double(10), double(18)))/2;
//						cout << "lat " << lat_energy << "\n";
						buf << " + min(0.0, " << lat_energy << ")";
					}
					if (s(i + 2, 0) > 0) // нижний мономер ПФ № 0  ряд r+2
					{	
			//			lat_4 = 1;
						neighbor_x_4 = s_x1(i + 2, 0) + left_nei_x(0);
						neighbor_y_4 = s_y1(i + 2, 0) + left_nei_y(0);
						neighbor_z_4 = s_z1(i + 2, 0) - double(6)/13;
						ex s = sqrt((new_current_right_x_2 - neighbor_x_4)*(new_current_right_x_2 - neighbor_x_4) + (new_current_right_y_2 - neighbor_y_4)*(new_current_right_y_2 - neighbor_y_4) + (new_current_right_z_2 - neighbor_z_4)*(new_current_right_z_2 - neighbor_z_4));
			//			total_energy += minimal_dim(double(0), (g_lat * kbt + 0.5* k_lat * s*s / pow(double(10), double(18)))/2);
						lat_energy = (g_lat * kbt + 0.5* k_lat * s*s / pow(double(10), double(18)))/2;
//						cout << "nei " << neighbor_x_2 <<"\n";
//						cout << "snei "<< s << "\n";
//						cout << "lat " << lat_energy << "\n";
						buf << " + min(0.0, " << lat_energy << ")";
					}
					
				}
				if (c == 0) // смотрим соседей слева
				{			
					if (s(i - 2, 12) > 0) // верхний мономер ПФ № 12  ряд r-2 (нижний монососед димера ПФ номер 0)
					{	
			//			lat_1 = 1;
						neighbor_x_1 = s_x2(i - 2, 12) + right_nei_x(12);
						neighbor_y_1 = s_y2(i - 2, 12) + right_nei_y(12);
						neighbor_z_1 = s_z2(i - 2, 12) + double(6)/13;
						ex s = sqrt((new_current_left_x_1 - neighbor_x_1)*(new_current_left_x_1 - neighbor_x_1) + (new_current_left_y_1 - neighbor_y_1)*(new_current_left_y_1 - neighbor_y_1) + (new_current_left_z_1 - neighbor_z_1)*(new_current_left_z_1 - neighbor_z_1) );
			//			total_energy += minimal_dim(double(0), (g_lat * kbt + 0.5* k_lat * s*s / pow(double(10), double(18)))/2);
						lat_energy = (g_lat * kbt + 0.5* k_lat * s*s / pow(double(10), double(18)))/2;
//						cout << "nei " << neighbor_x_2 <<"\n";
//						cout << "snei "<< s << "\n";
//						cout << "lat " << lat_energy << "\n";
						buf << " + min(0.0, " << lat_energy << ")";
					}
					if (s(i - 1, 12) > 0) // нижний мономер ПФ № 12  ряд r - 1 (верхний сосед)
					{	
			//			lat_3 = 1;
						neighbor_x_3 = s_x1(i - 1, 12) + right_nei_x(12);
						neighbor_y_3 = s_y1(i - 1, 12) + right_nei_y(12);
						neighbor_z_3 = s_z1(i - 1, 12) + double(6)/13;
						ex s = sqrt((new_current_left_x_2 - neighbor_x_3)*(new_current_left_x_2 - neighbor_x_3) + (new_current_left_y_2 - neighbor_y_3)*(new_current_left_y_2 - neighbor_y_3) + (new_current_left_z_2 - neighbor_z_3)* (new_current_left_z_2 - neighbor_z_3));
			//			total_energy += minimal_dim(double(0), (g_lat * kbt + 0.5* k_lat * s*s / pow(double(10), double(18)))/2);
						lat_energy = (g_lat * kbt + 0.5* k_lat * s*s / pow(double(10), double(18)))/2;
//						cout << "nei " << neighbor_x_2 <<"\n";
//												cout << "snei "<< s << "\n";
//												cout << "lat " << lat_energy << "\n";
						buf << " + min(0.0, " << lat_energy << ")";
					}
					
				}
				if (c >= 0 and c <= 11)//соседи справа
				{
					if (s(i, c + 1 ) > 0) 
					{	
			//			lat_4 = 1;
			//			lat_2 = 1;
						neighbor_x_2 = s_x1(i, c + 1) + left_nei_x(c + 1);
						neighbor_y_2 = s_y1(i, c + 1) + left_nei_y(c + 1);
						neighbor_z_2 = s_z1(i, c + 1) - double(6)/13;
						ex s = sqrt((new_current_right_x_1 - neighbor_x_2)*(new_current_right_x_1 - neighbor_x_2) + (new_current_right_y_1 - neighbor_y_2)*(new_current_right_y_1 - neighbor_y_2) + (new_current_right_z_1 - neighbor_z_2)* (new_current_right_z_1 - neighbor_z_2));
			//			total_energy += minimal_dim(double(0), (g_lat * kbt + 0.5* k_lat * s*s / pow(double(10), double(18)))/2);
						lat_energy = (g_lat * kbt + 0.5* k_lat * s*s / pow(double(10), double(18)))/2;
//						cout << "nei " << neighbor_x_2 <<"\n";
//												cout << "snei "<< s << "\n";
//												cout << "lat " << lat_energy << "\n";
						buf << " + min(0.0, " << lat_energy << ")";
						neighbor_x_4 = s_x2(i, c + 1) + left_nei_x(c + 1);
						neighbor_y_4 = s_y2(i, c + 1) + left_nei_y(c + 1);
						neighbor_z_4 = s_z2(i, c + 1) - double(6)/13;
						s = sqrt((new_current_right_x_2 - neighbor_x_4)*(new_current_right_x_2 - neighbor_x_4) + (new_current_right_y_2 - neighbor_y_4)*(new_current_right_y_2 - neighbor_y_4) + (new_current_right_z_2 - neighbor_z_4)*(new_current_right_z_2 - neighbor_z_4) );
			//			total_energy += minimal_dim(double(0), (g_lat * kbt + 0.5* k_lat * s*s / pow(double(10), double(18)))/2);
						lat_energy = (g_lat * kbt + 0.5* k_lat * s*s / pow(double(10), double(18)))/2;
						buf << " + min(0.0, " << lat_energy << ")";
					}

				}
				if (c >=1 and c <= 12)//соседи слева
				{
					if (s(i, c - 1) > 0) 
					{	
			//			lat_1 = 1;
			//			lat_3 = 1;
						neighbor_x_1 = s_x1(i, c - 1) + right_nei_x(c - 1);
						neighbor_y_1 = s_y1(i, c - 1) + right_nei_y(c - 1);
						neighbor_z_1 = s_z1(i, c - 1) + double(6)/13;
						ex s = sqrt((new_current_left_x_1 - neighbor_x_1)*(new_current_left_x_1 - neighbor_x_1) + (new_current_left_y_1 - neighbor_y_1)*(new_current_left_y_1 - neighbor_y_1) + (new_current_left_z_1 - neighbor_z_1)*(new_current_left_z_1 - neighbor_z_1) );
			//			total_energy += minimal_dim(double(0), (g_lat * kbt + 0.5* k_lat * s*s / pow(double(10), double(18)))/2);
						lat_energy = (g_lat * kbt + 0.5* k_lat * s*s / pow(double(10), double(18)))/2;
						buf << " + min(0.0, " << lat_energy << ")";		
//						cout << "nei " << neighbor_x_2 <<"\n";
//												cout << "snei "<< s << "\n";
//												cout << "lat " << lat_energy << "\n";
						neighbor_x_3 = s_x2(i, c - 1) + right_nei_x(c - 1);
						neighbor_y_3 = s_y2(i, c - 1) + right_nei_y(c - 1);
						neighbor_z_3 = s_z2(i, c - 1) + double(6)/13;
						s = sqrt((new_current_left_x_2 - neighbor_x_3)*(new_current_left_x_2 - neighbor_x_3) + (new_current_left_y_2 - neighbor_y_3)*(new_current_left_y_2 - neighbor_y_3) + (new_current_left_z_2 - neighbor_z_3)*(new_current_left_z_2 - neighbor_z_3) );
			//			total_energy += minimal_dim(double(0), (g_lat * kbt + 0.5* k_lat * s*s / pow(double(10), double(18)))/2);
						lat_energy = (g_lat * kbt + 0.5* k_lat * s*s / pow(double(10), double(18)))/2;
						buf << " + min(0.0, " << lat_energy << ")";
					}
				}
//				cout << "nei " << neighbor_x_2 <<"\n";
//										cout << "snei "<< s << "\n";
//										cout << "lat " << lat_energy << "\n";
				
		current_x_2 = new_current_x_2;
		current_y_2 = new_current_y_2;
		current_z_2 = new_current_z_2;
		add_to_total_tetta += s_tetta(i, c);
		add_to_total_fi += s_fi(i, c);
		}
	}
	
	
//	mwArray aa(l_p);
//	mwArray ab(m_p);
//	mwArray ac(n_p);
//	mwArray ad(lat_1);
//	mwArray ae(lat_2);
//	mwArray af(lat_3);
//	mwArray ag(lat_4);
//	mwArray ah(total_tetta);
//	mwArray ai(total_fi);
//	mwArray aj(curr_x_1_a);
//	mwArray ak(curr_x_1_b);
//	mwArray al(curr_y_1_a);
//	mwArray am(curr_y_1_b);
//	mwArray an(curr_z_1_a);
//	mwArray ao(curr_x_1_c);
//	mwArray ap(curr_y_1_c);
//	mwArray aq(neighbor_x_1);
//	mwArray ar(neighbor_x_2);
//	mwArray as(neighbor_x_3);
//	mwArray at(neighbor_x_4);
//	mwArray au(neighbor_y_1);
//	mwArray av(neighbor_y_2);
//	mwArray aw(neighbor_y_3);
//	mwArray ax(neighbor_y_4);
//	mwArray ay(neighbor_z_1);
//	mwArray az(neighbor_z_2);
//	mwArray bb(neighbor_z_3);
//	mwArray bc(neighbor_z_4);
	
//	if (s(r, c) == 1)
//	{
//	cout << r<<" " <<c<<" " << aa<<", "<< ab<<", "<< ac<<", "<< ad<<", "<< ae<<", "<< af<<", "<< ag<<", "<< ah<<", "<< ai<<", "<< aj<<", "<< ak<<", "<< al<<", "<< am<<", "<< an<<", "<< ao<<", "<< ap<<", "<< aq<<", "<< ar<<", "<< as<<", "<< at<<", "<< au<<", "<< av<<", "<< aw<<", "<< ax<<", "<< ay<<", "<< az<<", "<< bb<<", "<< bc <<"\n";
//	}
	
	
	string final_string = "10^30*( " + buf.str() + " )";
	char* final_function = new char[final_string.size() + 1];
	copy(final_string.begin(), final_string.end(), final_function);
	final_function[final_string.size()] = '\0';
	mwArray sum(1, 3, mxDOUBLE_CLASS);
//	mwArray[] sum = new mwArray[3];
//	mwArray sum;
	mwArray fina(final_function);
//	if (r > 9)
//	{
//	cout << "hydr " << r << " " << c << " " << final_function << "\n";
//	}
	
//	cout << final_function <<"\n";
	
	
//	minimize(1, sum, aa, ab, ac, ad, ae, af, ag, ah, ai, aj, ak, al, am, an, ao, ap, aq, ar, as, at, au, av, aw, ax, ay, az, bb, bc);
	minimize(1, sum, fina);
//	double *result = new double[3];
//	sum.GetData(result, 3);
//	double res1 = result[0];
//	cout << res1 <<  "\n";
//	sum = NULL;
	
	double tetta = sum(1);
	double fi = sum(2);
	double d = sum(3);
//	if (s(r, c) > 9)
//	{
//	std::cout << "hydr " << r<<" " <<c << " " << tetta << " " << fi << " " << d << "\n";
//	}
//	double tetta = 1;
//	double fi = 1;
//	double d = 1;
	delete[] final_function;
	s_d(r, c) = d;
	s_x1(r, c) = s_x2(r - 1, c) + (d + 4)*sin(total_tetta + tetta)*cos(total_fi + fi);
	s_y1(r, c) = s_y2(r - 1, c) + (d + 4)*sin(total_tetta + tetta)*sin(total_fi + fi);	
	s_z1(r, c) = s_z2(r - 1, c) + (d + 4)*cos(total_tetta + tetta);
	s_x2(r, c) = s_x2(r - 1, c) + (d + 8)*sin(total_tetta + tetta)*cos(total_fi + fi);
	s_y2(r, c) = s_y2(r - 1, c) + (d + 8)*sin(total_tetta + tetta)*sin(total_fi + fi);	
	s_z2(r, c) = s_z2(r - 1, c) + (d + 8)*cos(total_tetta + tetta);
	total_tetta += tetta;
	total_fi += fi;
	s_tetta(r, c) = tetta;
//	if (s(r, c) == 1.0000)
//	{
//		cout << "po " << s_tetta(r, c) << "\n";
//	}
	s_fi(r, c) = fi;
	define_all_energies(r, c, right_nei_x, right_nei_y, left_nei_x, left_nei_y, b_lat, b_curl, b_long, s_tetta, s_fi, s_d, s_x1, s_y1, s_z1, s_x2, s_y2, s_z2, s);
	if (c == 12)
	{			
		if (s(r + 1, 0) > 0) // верхний мономер ПФ № 0  ряд r+1
		{	
			define_lat_energy(r + 1, 0, right_nei_x, right_nei_y, left_nei_x, left_nei_y, b_lat, s, s_x1, s_y1, s_z1, s_x2, s_y2, s_z2);
		}
		if (s(r + 2, 0) > 0) // нижний мономер ПФ № 0  ряд r+2
		{	
			define_lat_energy(r + 2, 0, right_nei_x, right_nei_y, left_nei_x, left_nei_y, b_lat, s, s_x1, s_y1, s_z1, s_x2, s_y2, s_z2);
		}
		
	}
	if (c == 0) // смотрим соседей слева
	{			
		if (s(r - 2, 12) > 0) // верхний мономер ПФ № 12  ряд r-2 (нижний монососед димера ПФ номер 0)
		{	
			define_lat_energy(r - 2, 12, right_nei_x, right_nei_y, left_nei_x, left_nei_y, b_lat, s, s_x1, s_y1, s_z1, s_x2, s_y2, s_z2);
		}
		if (s(r - 1, 12) > 0) // нижний мономер ПФ № 12  ряд r - 1 (верхний сосед)
		{	
			define_lat_energy(r - 1, 12, right_nei_x, right_nei_y, left_nei_x, left_nei_y, b_lat, s, s_x1, s_y1, s_z1, s_x2, s_y2, s_z2);
		}
		
	}
	if (c >= 0 and c <= 11)//соседи справа
	{
		if (s(r, c + 1 ) > 0) 
		{	
			define_lat_energy(r, c + 1, right_nei_x, right_nei_y, left_nei_x, left_nei_y, b_lat, s, s_x1, s_y1, s_z1, s_x2, s_y2, s_z2);
		}

	}
	if (c >=1 and c <= 12)//соседи слева
	{
		if (s(r, c - 1) > 0) 
		{	
			define_lat_energy(r, c - 1, right_nei_x, right_nei_y, left_nei_x, left_nei_y, b_lat, s, s_x1, s_y1, s_z1, s_x2, s_y2, s_z2);
		}
	}
	for (int i = r + 1; i < leng(c); i++)
	{
		s_x1(i, c) = s_x2(i - 1, c) + (s_d(i, c) + 4)*sin(total_tetta + s_tetta(i, c))*cos(total_fi + s_fi(i, c)); //вместо матрицы можно передавать число
		s_y1(i, c) = s_y2(i - 1, c) + (s_d(i, c) + 4)*sin(total_tetta + s_tetta(i, c))*sin(total_fi + s_fi(i, c));	
		s_z1(i, c) = s_z2(i - 1, c) + (s_d(i, c) + 4)*cos(total_tetta + s_tetta(i, c));
		s_x2(i, c) = s_x2(i - 1, c) + (s_d(i, c) + 8)*sin(total_tetta + s_tetta(i, c))*cos(total_fi + s_fi(i, c));
		s_y2(i, c) = s_y2(i - 1, c) + (s_d(i, c) + 8)*sin(total_tetta + s_tetta(i, c))*sin(total_fi + s_fi(i, c));	
		s_z2(i, c) = s_z2(i - 1, c) + (s_d(i, c) + 8)*cos(total_tetta + s_tetta(i, c));
		define_all_energies(i, c, right_nei_x, right_nei_y, left_nei_x, left_nei_y, b_lat, b_curl, b_long, s_tetta, s_fi, s_d, s_x1, s_y1, s_z1, s_x2, s_y2, s_z2, s);
		if (c == 12)
		{			
			if (s(i + 1, 0) > 0) // верхний мономер ПФ № 0  ряд r+1
			{	
				define_lat_energy(i + 1, 0, right_nei_x, right_nei_y, left_nei_x, left_nei_y, b_lat, s, s_x1, s_y1, s_z1, s_x2, s_y2, s_z2);
			}
			if (s(i + 2, 0) > 0) // нижний мономер ПФ № 0  ряд r+2
			{	
				define_lat_energy(i + 2, 0, right_nei_x, right_nei_y, left_nei_x, left_nei_y, b_lat, s, s_x1, s_y1, s_z1, s_x2, s_y2, s_z2);
			}
			
		}
		if (c == 0) // смотрим соседей слева
		{			
			if (s(i - 2, 12) > 0) // верхний мономер ПФ № 12  ряд r-2 (нижний монососед димера ПФ номер 0)
			{	
				define_lat_energy(i - 2, 12, right_nei_x, right_nei_y, left_nei_x, left_nei_y, b_lat, s, s_x1, s_y1, s_z1, s_x2, s_y2, s_z2);
			}
			if (s(i - 1, 12) > 0) // нижний мономер ПФ № 12  ряд r - 1 (верхний сосед)
			{	
				define_lat_energy(i - 1, 12, right_nei_x, right_nei_y, left_nei_x, left_nei_y, b_lat, s, s_x1, s_y1, s_z1, s_x2, s_y2, s_z2);
			}
			
		}
		if (c >= 0 and c <= 11)//соседи справа
		{
			if (s(i, c + 1) > 0) 
			{	
				define_lat_energy(i, c + 1, right_nei_x, right_nei_y, left_nei_x, left_nei_y, b_lat, s, s_x1, s_y1, s_z1, s_x2, s_y2, s_z2);
			}

		}
		if (c >=1 and c <= 12)//соседи слева
		{
			if (s(i, c - 1) > 0) 
			{	
				define_lat_energy(i, c - 1, right_nei_x, right_nei_y, left_nei_x, left_nei_y, b_lat, s, s_x1, s_y1, s_z1, s_x2, s_y2, s_z2);
			}
		}
		
		total_tetta += s_tetta(i, c);
		total_fi += s_fi(i, c);
	}
}


void event(int N, int M, vec& right_nei_x, vec& right_nei_y, vec& left_nei_x, vec& left_nei_y, mat& b_lat, mat& b_curl, mat& b_long, mat& s_tetta, mat& s_fi, mat& s_d, mat& s_x1, mat& s_y1, mat& s_z1, mat& s_x2, mat& s_y2, mat& s_z2, int x, vec& l, vec& t, vec& leng, double& tm, mat& s)
{
  cube final_cube = cube(N, 13, 3);
  final_cube.fill(datum::inf);
  for (int c = 0; c < 13; c++) //on event
  {
	 final_cube(leng(c), c, 0) = -log(1 - double(rand()/double(RAND_MAX)))/(k_plus*conc);
	 
  }
//  cout << "plus " << k_plus*conc << "\n";
  for (int c = 0; c < 13; c++) //off event
  {
	double energy = 0;
	
    for (int r = leng(c) - 1; r >= 0; r--)
    {  
    	double free_energy = energy + b_lat(r, c) + b_curl(r, c);
    	energy = free_energy;
      double k_minus = (k_plus)/exp (-(free_energy + b_long(r, c))/kbt);
//      if (r > 9)
//      {
//    	  cout << "ener "<< r<<" " <<c << " " << -(free_energy + b_long(r, c))/kbt << "\n";
//      }
//      if (r > 9 and (c == 0 or c == 1))
//      {
////    	  cout << "k_minus "<< r<< " " << c << " " << b_long(r, c)<< " " << free_energy<< " " << free_energy/kbt<< " " << exp (-(free_energy + b_long(r, c))/kbt) << " " <<k_minus<<"\n";
////    	  cout << "k_minus " << -log(1 - double(rand()/double(RAND_MAX))) << " " << k_minus << " " << -log(1 - double(rand()/double(RAND_MAX)))/k_minus<<"\n";
//      }
	  final_cube(r, c, 1) = -log(1 - double(rand()/double(RAND_MAX)))/k_minus;   
    }
  }
  for (int c = 0; c < 13; c++) //Hydrolysis
  {
	  for (int r = 0; r < leng(c); r++)
		{
		  if(s(r, c) == 2 && r > 9) 
				  {
				    final_cube(r, c, 2) = -log(1 - double(rand()/double(RAND_MAX)))/(k_hydr);
				  }
		}
  }
      
  double cube_minimum;
  cube_minimum = final_cube.min();
//  cout << cube_minimum<< "\n";
  uword min_row, min_col, min_slice;
  final_cube.min(min_row, min_col, min_slice);
//  cout << "s" << s<< "\n";
//  cout << "sobi " << min_row<< min_col<<min_slice << "\n";
//  cout << "lat " << b_lat << "\n";
//  cout << "long " << b_long << "\n";
//  cout << "curl " << b_curl << "\n";
//  
//  cout << "fi" << final_cube<<"\n";
  
  if (min_slice == 1)
  {
	  
    for (int i = min_row; i < leng(min_col) + 1; i++)
    {
    s(i, min_col) = 0;
    s_x1(i, min_col) = 0;
    s_y1(i, min_col) = 0;
    s_z1(i, min_col) = 0;
    s_x2(i, min_col) = 0;
    s_y2(i, min_col) = 0;
    s_z2(i, min_col) = 0;
    b_lat(i, min_col) = 0;
    b_curl(i, min_col) = 0;
    b_long(i, min_col) = 0;
    s_tetta(i, min_col) = 0;
    s_fi(i, min_col) = 0;
    s_d(i, min_col) = 0;
    if (min_col == 12)
	{			
		if (s(i + 1, 0) > 0) // верхний мономер ПФ № 0  ряд r+1
		{	
			define_lat_energy(i + 1, 0, right_nei_x, right_nei_y, left_nei_x, left_nei_y, b_lat, s, s_x1, s_y1, s_z1, s_x2, s_y2, s_z2);
		}
		if (s(i + 2, 0) > 0) // нижний мономер ПФ № 0  ряд r+2
		{	
			define_lat_energy(i + 2, 0, right_nei_x, right_nei_y, left_nei_x, left_nei_y, b_lat, s, s_x1, s_y1, s_z1, s_x2, s_y2, s_z2);
		}
		
	}
	if (min_col == 0) // смотрим соседей слева
	{			
		
		if (s(i - 2, 12) > 0) // верхний мономер ПФ № 12  ряд r-2 (нижний монососед димера ПФ номер 0)
		{	
			define_lat_energy(i - 2, 12, right_nei_x, right_nei_y, left_nei_x, left_nei_y, b_lat, s, s_x1, s_y1, s_z1, s_x2, s_y2, s_z2);
		}
		if (s(i - 1, 12) > 0) // нижний мономер ПФ № 12  ряд r - 1 (верхний сосед)
		{	
			define_lat_energy(i - 1, 12, right_nei_x, right_nei_y, left_nei_x, left_nei_y, b_lat, s, s_x1, s_y1, s_z1, s_x2, s_y2, s_z2);
		}
		
	}
	if (min_col >= 0 and min_col <= 11)//соседи справа
	{
		if (s(i, min_col + 1 ) > 0) 
		{	
			define_lat_energy(i, min_col + 1, right_nei_x, right_nei_y, left_nei_x, left_nei_y, b_lat, s, s_x1, s_y1, s_z1, s_x2, s_y2, s_z2);
		}

	}
	if (min_col >=1 and min_col <= 12)//соседи слева
	{
		if (s(i, min_col - 1) > 0) 
		{	
			define_lat_energy(i, min_col - 1, right_nei_x, right_nei_y, left_nei_x, left_nei_y, b_lat, s, s_x1, s_y1, s_z1, s_x2, s_y2, s_z2);
		}
	}
    }
    leng(min_col) = min_row;
    int non_zeroo = 0;
     for (int c = 0; c < 13; c++)
     {
       for (int r = 0; r < N; r++)
       {
         if (s(r, c) > 0)
         {
        	 non_zeroo += 1;
         }
       }
     }
    vec non_z = vec(non_zeroo, fill::zeros);
    double k = 0;
    for (int c = 0; c < 13; c++)
	 {
	   for (int r = 0; r < N; r++)
	   {
		 if (s(r, c) > 0)
		 {
			 non_z(k) = r * 13 + c;
			 k += 1;
		 }
	   }
	 }
    int numb, row_min, col_min;
    for (int u = 0; u < 3*non_zeroo; u++)
    {
    	numb = rand() % non_zeroo;
    	row_min = int(double(non_z(numb))/double(13));
    	col_min = non_z(numb) - 13 * row_min;
    	if (row_min > 9)
    	{
    		minimize_energy(row_min, col_min, right_nei_x, right_nei_y, left_nei_x, left_nei_y, b_lat, b_curl, b_long, s_tetta, s_fi, s_d, s_x1, s_y1, s_z1, s_x2, s_y2, s_z2, s, leng);
    	}
    }
    
  }
  else if (min_slice == 0)
  {
    s(min_row, min_col) = 2;
    minimize_energy(min_row, min_col, right_nei_x, right_nei_y, left_nei_x, left_nei_y, b_lat, b_curl, b_long, s_tetta, s_fi, s_d, s_x1, s_y1, s_z1, s_x2, s_y2, s_z2, s, leng);
    if (min_col == 12)
	{			
		if (s(min_row + 1, 0) > 0) // верхний мономер ПФ № 0  ряд r+1
		{	
			define_lat_energy(min_row + 1, 0, right_nei_x, right_nei_y, left_nei_x, left_nei_y, b_lat, s, s_x1, s_y1, s_z1, s_x2, s_y2, s_z2);
		}
		if (s(min_row + 2, 0) > 0) // нижний мономер ПФ № 0  ряд r+2
		{	
			define_lat_energy(min_row + 2, 0, right_nei_x, right_nei_y, left_nei_x, left_nei_y, b_lat, s, s_x1, s_y1, s_z1, s_x2, s_y2, s_z2);
		}
		
	}
	if (min_col == 0) // смотрим соседей слева
	{			
		if (s(min_row - 2, 12) > 0) // верхний мономер ПФ № 12  ряд r-2 (нижний монососед димера ПФ номер 0)
		{	
			define_lat_energy(min_row - 2, 12, right_nei_x, right_nei_y, left_nei_x, left_nei_y, b_lat, s, s_x1, s_y1, s_z1, s_x2, s_y2, s_z2);
		}
		if (s(min_row - 1, 12) > 0) // нижний мономер ПФ № 12  ряд r - 1 (верхний сосед)
		{	
			define_lat_energy(min_row - 1, 12, right_nei_x, right_nei_y, left_nei_x, left_nei_y, b_lat, s, s_x1, s_y1, s_z1, s_x2, s_y2, s_z2);
		}
		
	}
	if (min_col >= 0 and min_col <= 11)//соседи справа
	{
		if (s(min_row, min_col + 1) > 0) 
		{	
			define_lat_energy(min_row, min_col + 1, right_nei_x, right_nei_y, left_nei_x, left_nei_y, b_lat, s, s_x1, s_y1, s_z1, s_x2, s_y2, s_z2);
		}

	}
	if (min_col >=1 and min_col <= 12)//соседи слева
	{
		if (s(min_row, min_col - 1) > 0) 
		{	
			define_lat_energy(min_row, min_col - 1, right_nei_x, right_nei_y, left_nei_x, left_nei_y, b_lat, s, s_x1, s_y1, s_z1, s_x2, s_y2, s_z2);
		}
	}
    leng(min_col) += 1;
    int non_zeroo = 0;
	 for (int c = 0; c < 13; c++)
	 {
	   for (int r = 0; r < N; r++)
	   {
		 if (s(r, c) > 0)
		 {
			 non_zeroo += 1;
		 }
	   }
	 }
	vec non_z = vec(non_zeroo, fill::zeros);
	double k = 0;
	for (int c = 0; c < 13; c++)
	 {
	   for (int r = 0; r < N; r++)
	   {
		 if (s(r, c) > 0)
		 {
			 non_z(k) = r * 13 + c;
			 k += 1;
		 }
	   }
	 }
	int numb, row_min, col_min;
	for (int u = 0; u < 3*non_zeroo; u++)
	{
		numb = rand() % non_zeroo;
		row_min = int(double(non_z(numb))/double(13));
		col_min = non_z(numb) - 13 * row_min;
		if (row_min > 9)
		{
			minimize_energy(row_min, col_min, right_nei_x, right_nei_y, left_nei_x, left_nei_y, b_lat, b_curl, b_long, s_tetta, s_fi, s_d, s_x1, s_y1, s_z1, s_x2, s_y2, s_z2, s, leng);
		}
	}
  }
  else if (min_slice == 2)
    {
	  s(min_row, min_col) = 1;
    }  
//  cout << "s " << s<<"\n";
//    cout << "c" << b_curl<< "\n";
//    cout << "la " << b_lat<<"\n";
//  if (x == 4)
//  {
//	  ofstream outfile;
//	  std::string fxname = "/home/marge/vis2/1x" + to_string(x) + ".txt";
//	  std::string fyname = "/home/marge/vis2/1y" + to_string(x) + ".txt";
//	  std::string fzname = "/home/marge/vis2/1z" + to_string(x) + ".txt";
//	  std::string sxname = "/home/marge/vis2/2x" + to_string(x) + ".txt";
//	  std::string syname = "/home/marge/vis2/2y" + to_string(x) + ".txt";
//	  std::string szname = "/home/marge/vis2/2z" + to_string(x) + ".txt";
//    std::string sname = "/home/marge/vis2/s" + to_string(x) + ".txt";
//	  outfile.open(fxname,  std::ios_base::app);
//	       for (int r = 0; r < N; r++)
//	       {
//	    	   for (int c = 0; c < 12; c++)
//	    	   {
//	    		   outfile << s_x1(r, c) << ",";
//	    	   }
//	    	   outfile << s_x1(r, 12) << "\n";
//	       }
//	       outfile.close();
//	       outfile.open(fyname,  std::ios_base::app);
//	           for (int r = 0; r < N; r++)
//	           {
//	        	   for (int c = 0; c < 12; c++)
//	        	   {
//	        		   outfile << s_y1(r, c) << ",";
//	        	   }
//	        	   outfile << s_y1(r, 12) << "\n";
//	           }
//	           outfile.close();
//	           outfile.open(fzname,  std::ios_base::app);
//	               for (int r = 0; r < N; r++)
//	               {
//	            	   for (int c = 0; c < 12; c++)
//	            	   {
//	            		   outfile << s_z1(r, c) << ",";
//	            	   }
//	            	   outfile << s_z1(r, 12) << "\n";
//	               }
//	               outfile.close();
//	               outfile.open(sxname,  std::ios_base::app);
//	                   for (int r = 0; r < N; r++)
//	                   {
//	                	   for (int c = 0; c < 12; c++)
//	                	   {
//	                		   outfile << s_x2(r, c) << ",";
//	                	   }
//	                	   outfile << s_x2(r, 12) << "\n";
//	                   }
//	                   outfile.close();
//	                   
//	                   outfile.open(syname,  std::ios_base::app);
//	                       for (int r = 0; r < N; r++)
//	                       {
//	                    	   for (int c = 0; c < 12; c++)
//	                    	   {
//	                    		   outfile << s_y2(r, c) << ",";
//	                    	   }
//	                    	   outfile << s_y2(r, 12) << "\n";
//	                       }
//	                       outfile.close();
//	                       outfile.open(szname,  std::ios_base::app);
//	                           for (int r = 0; r < N; r++)
//	                           {
//	                        	   for (int c = 0; c < 12; c++)
//	                        	   {
//	                        		   outfile << s_z2(r, c) << ",";
//	                        	   }
//	                        	   outfile << s_z2(r, 12) << "\n";
//	                           }
//	                           outfile.close();
//  outfile.open(sname,  std::ios_base::app);
//  	       for (int r = 0; r < N; r++)
//  	       {
//  	    	   for (int c = 0; c < 12; c++)
//  	    	   {
//  	    		   outfile << s(r, c) << ",";
//  	    	   }
//  	    	   outfile << s(r, 12) << "\n";
//  	       }
//  	       outfile.close();
//  }
 t(x) = cube_minimum;
 tm += cube_minimum;
 double non_zer = 0;
 for (int c = 0; c < 13; c++)
 {
   for (int r = 0; r < N; r++)
   {
     if (s(r, c) > 0)
     {
	non_zer += 1;
     }
   }
 }
 l(x) = (double)(non_zer / 13.0f - 10.0f)*0.008;
 cout << tm << "," << l(x) << "\n";
 
 
 
		  
}


int main(int argc, char** argv)
{
	if( !mclInitializeApplication(NULL,0) )
	  {
		  std::cout << "Could not initialize the application.\n";
	  	  exit(1);
	  }
	  if (!libminimInitialize())
	  {
	  	  std::cout << "Could not initialize the library.\n";
	  	  exit(1);
	  }
  cout << "time,length" << "\n";
  struct timeb tp;
  ftime(&tp);
  srand(static_cast<unsigned int>(getpid()) ^
  static_cast<unsigned int>(pthread_self()) ^
  static_cast<unsigned int >(tp.millitm));
  //srand(time(NULL));
  
  int N = atoi(argv[1]);
  int M = atoi(argv[2]);
  vec l = vec(M, fill::zeros);
  vec t = vec(M, fill::zeros);
  vec leng = vec(13, fill::zeros);
  for (int a = 0; a < 13; a++)
  {
	  leng(a) = 15; // 10
  }
  double tm = 0.0f;
  mat s = mat(N, 13, fill::zeros);
  mat b_lat = mat(N, 13, fill::zeros);
  mat b_curl = mat(N, 13, fill::zeros);
  mat b_long = mat(N, 13, fill::zeros);
  mat s_tetta = mat(N, 13, fill::zeros);
  mat s_fi = mat(N, 13, fill::zeros);
  mat s_d = mat(N, 13, fill::zeros);
  mat s_x1 = mat(N, 13, fill::zeros);
  mat s_y1 = mat(N, 13, fill::zeros);
  mat s_z1 = mat(N, 13, fill::zeros);
  mat s_x2 = mat(N, 13, fill::zeros);
  mat s_y2 = mat(N, 13, fill::zeros);
  mat s_z2 = mat(N, 13, fill::zeros);
  for (int i = 0; i < 13; i++)
  {
    for (int n = 0; n < 10; n++)
    {
      s(n, i) = 2;    
    }
  }
  for (int i = 0; i < 13; i++) //этот блок убрать 
    {
      for (int n = 10; n < 15; n++)
      {
        s(n, i) = 1;    
      }
    }
//  for (int i = 0; i < 13; i++)
//    {
//      for (int n = 0; n < 10; n++)
//      {
//        b_lat(n, i) = 2*g_lat;
//      }
//    }
  for (int i = 0; i < 13; i++)
      {
        for (int n = 0; n < 15; n++) //10
        {
          b_long(n, i) = g_long*kbt;
        }
      }
  for (int n = 0; n < 15; n++)   //n < 10
    {
      for (int i = 0; i < 13; i++)
      {
        if (i == 0)
        {
        	s_x1(n, i) = sqrt(16 - double(144)/169)/double(2)/sin(pi/13);
        	s_y1(n, i) = 0;
        	s_x2(n, i) = sqrt(16 - double(144)/169)/double(2)/sin(pi/13);
        	s_y2(n, i) = 0;
        	
        }
        else
        {
        	s_x1(n, i) = s_x1(n, i - 1) * cos(2*pi/13) + s_y1(n, i - 1) * sin(2*pi/13);
        	s_y1(n, i) = -s_x1(n, i - 1) * sin(2*pi/13) + s_y1(n, i - 1) * cos(2*pi/13);
        	s_x2(n, i) = s_x2(n, i - 1) * cos(2*pi/13) + s_y2(n, i - 1) * sin(2*pi/13);
        	s_y2(n, i) = -s_x2(n, i - 1) * sin(2*pi/13) + s_y2(n, i - 1) * cos(2*pi/13);
        	
        }
        s_z1(n, i) = n * 8 + i * double(12)/13;
        s_z2(n, i) = n * 8 + 4 + i * double(12)/13;
      }
    }
  vec right_nei_x = vec(13, fill::zeros);
  vec left_nei_x = vec(13, fill::zeros);
  vec right_nei_y = vec(13, fill::zeros);
  vec left_nei_y = vec(13, fill::zeros);
  right_nei_x(12) = (s_x2(1, 0) - s_x1(0, 12))/2;
  right_nei_y(12) = (s_y2(1, 0) - s_y1(0, 12))/2;
  left_nei_x(0) = (s_x2(0, 12) - s_x1(2, 0))/2;
  left_nei_y(0) = (s_y2(0, 12) - s_y1(2, 0))/2;
  for (int i = 0; i < 12; i++)
  {
	  right_nei_x(i) = (s_x1(0, i + 1) - s_x1(0, i))/2; 
	  right_nei_y(i) = (s_y1(0, i + 1) - s_y1(0, i))/2;
  }
  for (int i = 1; i < 13; i++)
  {
	  left_nei_x(i) = (s_x1(0, i - 1) - s_x1(0, i))/2;
	  left_nei_y(i) = (s_y1(0, i - 1) - s_y1(0, i))/2;
  }
  for (int i = 0; i < 13; i++)
  {
  for (int n = 0; n < 10; n++)
	  {
		b_lat(n, i) = -6.4*kbt;
	  }
  }
  for (int i = 0; i < 13; i++) //этот блок убрать 
    {
    for (int n = 10; n < 15; n++)
  	  {
  		b_lat(n, i) = -6.4*kbt;
  	  }
    }
  for (int i = 0; i < 13; i++) //этот блок убрать 
  {
  for (int n = 10; n < 15; n++)
	  {
		b_curl(n, i) = -g_curl*kbt;
	  }
  }
  b_lat(0, 0) += 3.2*kbt;
  b_lat(13, 12) += 1.6*kbt; //8
  b_lat(14, 12) += 3.2*kbt; //9
  int x;
  for (x = 0; x < M; x++)
  {
	  event(N, M, right_nei_x, right_nei_y, left_nei_x, left_nei_y, b_lat, b_curl, b_long, s_tetta, s_fi, s_d, s_x1, s_y1, s_z1, s_x2, s_y2, s_z2, x, l, t, leng, tm, s);
  }
  libminimTerminate();
  mclTerminateApplication();
  return 0;
}





