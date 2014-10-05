#include <iostream>
#include <cstdlib>
#include <armadillo>
#include <time.h>

using namespace std;
using namespace arma;

const double k_minus_GTP = 0.02;
const double k_minus_GDP = 20;
const double k_plus = 12.5;
const double k_break1 = 800;
const double k_break2 = 180;
const double k_break3 = 140;
const double k_break4 = 400;
const double k_break5 = 90;
const double k_break6 = 70;
const double k_bond = 100;
const double k_hydr = 0.70;
const int minus_event = 0;
const int plus_event = 1;
const int lateral_break = 2;
const int lateral_bond = 3;
const int hydrolysis = 4;


void event(int N, umat& s, umat& b, int x, vec& l, vec& t, double& tm)
{
  //cout << s << '\n';
  //cout << b << '\n';
  cube final_cube = cube(N, 13, 5);
  final_cube.fill(datum::inf);
  for (int c = 0; c < 13; c++)
  {
    for (int r = 0; r < N; r++)
    {
      bool k = false; //off events
      if (s(r, c) == 0)
      {
	  break;
      }
      if (s(r, c) > 0 && b(r, c) == 0)
      {
        if (c > 0)
        {
          if (b(r, c-1) == 0)
          {
              k = true;
          }
        }
        else
        {
          if (r > 2)
          {
            if (b(r - 2, 12) < 2 && b(r - 1, 12) == 0)
            {
              k = true;
            }
          }
        }
      }
      if (k)
      {
        if (s(r - 1, c) == 2)
	{
          final_cube(r, c, 0) = -log(double(rand()) / double((RAND_MAX)))/k_minus_GTP;
	 }
	 else
	 {
           final_cube(r, c, 0) = -log(double(rand()) / double((RAND_MAX)))/k_minus_GDP;
         }
        }
      if (r < N - 1) //on event
      {
        if (s(r, c) > 0 && s(r + 1, c) == 0) {
	  final_cube(r + 1, c, 1) = -log(double(rand()) / double((RAND_MAX)))/k_plus;
        }
      }
      if (b(r, c) > 0 && b(r + 1, c) == 0) // break lateral event
      {
	double p = 1;
	double k_break;
        if (c == 12)
        {
          if (b(r, c) == 1) // one monomer bond
          {
            if (b(r + 1, 0) == 1 && b(r, 11) == 1) // check if neighbor lateral bonds exist
            {
              p = 1000;
            }
            if (s(r, c) == 1 && s(r + 1, 0) == 1) // both GDP
            {
              k_break = k_break1;
            }
            else if (s(r, c) == 2 && s(r + 1, 0) == 2) //both GTP
            {
              k_break = k_break3;
            }
            else if ((s(r, c) == 2 && s(r + 1, 0) == 1) or (s(r, c) == 1 && s(r + 1, 0) == 2)) //different
            {
              k_break = k_break2;
            }
          }
          else // both monomer bonds
          {
            if (b(r + 2, 0) == 1 && b(r, 11) == 1)
            {
              p = 1000;
            }
            if (s(r, c) == 1 && s(r + 2, 0) == 1) // both GDP
            {
              k_break = k_break1;
            }
            else if (s(r, c) == 2 && s(r + 2, 0) == 2) //both GTP
            {
              k_break = k_break3;
            }
            else if ((s(r, c) == 2 && s(r + 2, 0) == 1) or (s(r, c) == 1 && s(r + 2, 0) == 2)) //different
            {
              k_break = k_break2;
            }
          }
        }
        else if (c == 0)
        {
          if (r > 2)
          {
            if (b(r - 2, 12) == 2 && b(r, 1) == 1 && b(r - 1, 12) >= 1) //if neighbor lat bonds exist
            {
              p = 1000;
            }
            if (s(r, c) == 1 && s(r, 1) == 1) //both GDP
            {
              k_break = k_break4;
            }
            else if ((s(r, c) == 1 && s(r, 1) == 2) || (s(r, c) == 2 && s(r, 1) == 1))
            {
              k_break = k_break5;
            }
            else if (s(r, c) == 2 && s(r, 1) == 2) // both GTP
            {
              k_break = k_break6;
            }
          }
        }
        else if (c == 11)//PF 12
        {
          if (b(r, 12) == 2 && b(r, 10) == 1)
          {
            p = 1000;
          }
          if (s(r, c) == 1 && s(r, 12) == 1) //both GDP
          {
            k_break = k_break4;
          }
          else if ((s(r, c) == 2 && s(r, 12) == 1) || (s(r, c) == 1 && s(r, 12) == 2))
          {
            k_break = k_break5;
          }
          else if (s(r, c) == 2 && s(r, 12) == 2)
          {
            k_break = k_break6;
          }
        }
        else //PF 2-11
        {
          if (b(r, c - 1) == 1 && b(r, c + 1) == 1)
          {
            p = 1000;
          }
          if (s(r, c) == 1 && s(r, c + 1) == 1) //both GDP
          {
            k_break = k_break4;
          }
          else if ((s(r, c) == 1 && s(r, c + 1) == 2) || (s(r, c) == 2 && s(r, c + 1) == 1))
          {
            k_break = k_break5;
          }
          else if (s(r, c) == 2 && s(r, c + 1) == 2)
          {
            k_break = k_break6;
          }
        }
//        cout << r << c << k_break/p << '\n';
        if (r > 2)
        {
           final_cube(r, c, 2) = -log(double(rand()) / double((RAND_MAX)))/(k_break/p);
        }
      }
      if (r > 1)
      {
	if (c < 12 && b(r - 1, c) > 0) // bond lateral event
	{
	  if (b(r, c) == 0 && s(r, c + 1) > 0) //no lat. bond; right neighbor exists; there's lat. bond below
	  {
	     final_cube(r, c, 3) = -log(double(rand()) / double((RAND_MAX)))/(k_bond);
	  }
	}
	else if (c == 12 && b(r - 1, c) == 2)//13th PF
	{
	  if (r > 2 && ((b(r, c) == 0 && s(r + 1, 0) > 0) || (b(r, c) == 1 && s(r + 2, 0) > 0)))
	  {
	    final_cube(r, c, 3) = -log(double(rand()) / double((RAND_MAX)))/(k_bond);
	  }
	}
      }
	if (r > 9) //Hydrolysis event
	{
	  if(s(r, c) == 2 && s(r + 1, c) != 0)
	  {
	    final_cube(r, c, 4) = -log(double(rand()) / double((RAND_MAX)))/(k_hydr);
	  }
	}
      }
    }
  double cube_minimum;
  cube_minimum = final_cube.min();
  uword min_row, min_col, min_slice;
  final_cube.min(min_row, min_col, min_slice);
    //cout << min_slice << ' ' << min_row << ' ' << min_col << '\n' ;
  if (min_slice == 0)
  {
    for (int i = min_row; i < N; i++)
    {
      s(i, min_col) = 0;
    }
  }
  else if (min_slice == 1)
  {
    s(min_row, min_col) = 2;
  }
  else if (min_slice == 2)
  {
    b(min_row, min_col) = b(min_row, min_col) - 1;
  }
  else if (min_slice == 3)
  {
    b(min_row, min_col) = b(min_row, min_col) + 1;
  }
  else
  {
    s(min_row, min_col) = 1;
  }
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
  l(x) = (double)(non_zer / 13.0f - 10.0f);
  cout << tm << "," << l(x) << "\n";
//  cout << final_cube << '\n';
}


int main(int argc, char** argv)
{
  srand(time(NULL));
  int N = atoi(argv[1]);
  int M = atoi(argv[2]);
  vec l = vec(M, fill::zeros);
  vec t = vec(M, fill::zeros);
  double tm = 0.0f;
  umat s_in = umat(N, 13, fill::zeros);
  umat b_in = umat(N, 13, fill::zeros);
  for (int i = 0; i < 13; i++)
  {
    for (int n = 0; n < 10; n++)
    {
      s_in(n, i) = 2;
      b_in(n, i) = (i == 12) ? 2 : 1;
    }
  }
  int x;
  for (x = 0; x < M; x++)
  {
    event(N, s_in, b_in, x, l, t, tm);
  }
//  s_in.print();
  return 0;
}

