#include <iostream>
#include <cstdlib>
#include <armadillo>
#include <time.h>
#include <random>
#include <sys/timeb.h>
#include <fstream>
#include <math.h> 
#include <stdexcept> 

using namespace std;
using namespace arma;


const double k_plus = 40; 
const double g_lat = -5.7;
const double g_long = -6.8;
const double g_lat_gdp = -3.2;
const double k_hydr = 0.95;
//const int minus_event = 0;
//const int plus_event = 1;
//const int hydrolysis = 2;

double neighbors(umat& s, umat& s_tagged, int R, int C) //count GTP-neighbors in 4 positions (to define if T is tagged or not)
{
	double neibor_tag = 0;
	if (C == 0) //1st PF
	              {
	                if (s(R, 1) == 2 and s_tagged(R, 1) == 0)
	                {
	                	neibor_tag += 1;
	                }
	                if (s(R + 1, 1) == 2 and s_tagged(R + 1, 1) == 0)
	                {
	                	neibor_tag += 1;
	                }
	                if (R >= 2)
	                {
	                if (s(R - 1, 12) == 2  and s_tagged(R - 1, 12) == 0) 
	                {
	                	neibor_tag += 1;
	                }
	                if (s(R, 12) == 2 and s_tagged(R, 12) == 0 )
	                {
	                	neibor_tag += 1;
	                }
	                }
	                
	            	  
	              }

	              if (C == 12) //13th PF
	              {
	            	  if (s(R, 11)  == 2 and s_tagged(R, 11) == 0)
	            	  {
	            		  neibor_tag += 1;
	            	  }
	            	  if (s(R + 1, 11) == 2 and s_tagged(R + 1, 11) == 0)
	            	  {
	            		  neibor_tag += 1;
	            	  }
	            	  if (s(R + 2, 0) == 2 and s_tagged(R + 2, 0) == 0) 
	            	  {
	            		  neibor_tag += 1;
	            	  }
	            	  if (s(R + 3, 0) == 2 and s_tagged(R + 3, 0) == 0)
	                  {
	                  	neibor_tag += 1;
	                  }
	              }
	              if (C < 12 && C > 0)
	              {
	            	  if (s(R, C - 1) == 2 and s_tagged(R, C - 1) == 0)
	            	  {
	            		  neibor_tag += 1;
	            	  }
	            	  if (s(R + 1, C - 1) == 2 and s_tagged(R + 1, C - 1) == 0)
	            	  {
						  neibor_tag += 1;
					  }
	            	  if (s(R, C + 1) == 2 and s_tagged(R, C + 1) == 0)
	            	  {
	            		  neibor_tag += 1;
	            	  }
	            	  if (s(R + 1, C + 1) == 2 and s_tagged(R + 1, C + 1) == 0)
					  {
						  neibor_tag += 1;
					  }
	              }
	              return neibor_tag;
}



double neibors_lat(umat& s, int R, int C) //count lateral neighbors
{
  double neibor = 0;
  if (C == 0) //1st PF
  {
	if (s(R, 1) > 0)
	{
		neibor += 1;
	}
	
	if (R >= 2)
	{
	if (s(R - 1, 12) > 0) 
	{
		neibor += 1;
	}
	if (s(R - 1, 12) == 0 && s(R - 2, 12) > 0)
	{
		neibor += 0.5;
	}
	}
	
	  
  }

  if (C == 12) //13th PF
  {
	  if (s(R, 11)  > 0)
	  {
		  neibor += 1;
	  }
	  
	  if (s(R + 2, 0) > 0) 
	  {
		  neibor += 1;
	  }
	  if (s(R + 2, 0) == 0 && s(R + 1, 0) > 0)
	  {
		neibor += 0.5;
	  }
  }
  if (C < 12 && C > 0)
  {
	  if (s(R, C - 1) > 0)
	  {
		  neibor += 1;
	  }
	  
	  if (s(R, C + 1) > 0)
	  {
		  neibor += 1;
	  }
	 
  }
  return neibor;
}

void change_neighbors_state(umat& s, mat& b, umat& s_tagged, int min_row, int min_col, vec& leng) //change matrix of bonds, dimers, tagged dimers in case of tubulin assotiation/dissotiation
{
	double lat_bond;
	if (min_col == 0 && min_row > 1)
	{
	   if (neighbors(s, s_tagged, min_row, 1) >= 2)
	   {
		   s_tagged(min_row, 1) = 0;
	   }
	   else
	   {
		   if (s(min_row, 1) == 2 and s(min_row - 1, 1) == 1)
		   {
			   s_tagged(min_row, 1) = 1;
		   }
	   }
	   if (neighbors(s, s_tagged, min_row - 1, 1) >= 2)
		  {
		   s_tagged(min_row - 1, 1) = 0;
		  }
	   else
	   	   {
	   		   if (s(min_row - 1, 1) == 2 and s(min_row - 2, 1) == 1)
	   		   {
	   			   s_tagged(min_row - 1, 1) = 1;
	   		   }
	   	   }
	   for (int ch = min_row - 1; ch < leng(1); ch++)
	   {
		   if ((s(ch, 1) == 2 && s_tagged(ch, 1) == 1) or s(ch, 1) == 1)
		   {
			   lat_bond = g_lat_gdp;
		   }
		   if (s(ch, 1) == 2 && s_tagged(ch, 1) == 0)
		   {
			   lat_bond = g_lat;
		   }
		   if (s(ch, 1) > 0)
		   {
			   b(ch, 1) = neibors_lat(s, ch, 1) * lat_bond;
		   }
	   }
	   if (neighbors(s, s_tagged, min_row - 2, 12) >= 2)
			  {
			   s_tagged(min_row - 2, 12) = 0;
			  }
	   else
	   	   	   {
	   	   		   if (s(min_row - 2, 12) == 2 and s(min_row - 3, 12) == 1)
	   	   		   {
	   	   			   s_tagged(min_row - 2, 12) = 1;
	   	   		   }
	   	   	   }
	   if (neighbors(s, s_tagged, min_row - 1, 12) >= 2)
			  {
			   s_tagged(min_row - 1, 12) = 0;
			  }
	   else
	   	   	   {
	   	   		   if (s(min_row - 1, 12) == 2 and s(min_row - 2, 12) == 1)
	   	   		   {
	   	   			   s_tagged(min_row - 1, 12) = 1;
	   	   		   }
	   	   	   }
	   if (neighbors(s, s_tagged, min_row - 3, 12) >= 2)
	   			  {
	   			   s_tagged(min_row - 3, 12) = 0;
	   			  }
	   else
	   	   	   {
	   	   		   if (s(min_row - 3, 12) == 2 and s(min_row - 4, 12) == 1)
	   	   		   {
	   	   			   s_tagged(min_row - 3, 12) = 1;
	   	   		   }
	   	   	   }
	   for (int ch = min_row - 3; ch < leng(12); ch++)
			  {
			   if ((s(ch, 12) == 2 && s_tagged(ch, 12) == 1) or s(ch, 12) == 1)
			   {
				   lat_bond = g_lat_gdp;
			   }
			   if (s(ch, 12) == 2 && s_tagged(ch, 12) == 0)
			   {
				   lat_bond = g_lat;
			   }
			   if (s(ch, 12) > 0)
			   {
				   b(ch, 12) = neibors_lat(s, ch, 12) * lat_bond;
			   }
			  }
	}
	if (min_col == 12 && min_row > 0)
	{
		if (neighbors(s, s_tagged, min_row, 11) >= 2)
			   {
				   s_tagged(min_row, 11) = 0;
			   }
		else
		   {
			   if (s(min_row, 11) == 2 and s(min_row - 1, 11) == 1)
			   {
				   s_tagged(min_row, 11) = 1;
			   }
		   }
		if (neighbors(s, s_tagged, min_row - 1, 11) >= 2)
			   {
				   s_tagged(min_row - 1, 11) = 0;
			   }
		else
		   {
			   if (s(min_row - 1, 11) == 2 and s(min_row - 2, 11) == 1)
			   {
				   s_tagged(min_row - 1, 11) = 1;
			   }
		   }
		for (int ch = min_row - 1; ch < leng(11); ch++)
			   {
				   if ((s(ch, 11) == 2 && s_tagged(ch, 11) == 1) or s(ch, 11) == 1)
				   {
					   lat_bond = g_lat_gdp;
				   }
				   if (s(ch, 11) == 2 && s_tagged(ch, 11) == 0)
				   {
					   lat_bond = g_lat;
				   } 
				   if (s(ch, 11) > 0)
				   {
					   b(ch, 11) = neibors_lat(s, ch, 11) * lat_bond;
				   }
			   }
		if (neighbors(s, s_tagged, min_row + 1, 0) >= 2)
			   {
				   s_tagged(min_row + 1, 0) = 0;
			   }
		else
		   {
			   if (s(min_row + 1, 0) == 2 and s(min_row, 0) == 1)
			   {
				   s_tagged(min_row + 1, 0) = 1;
			   }
		   }
		if (neighbors(s, s_tagged, min_row + 2, 0) >= 2)
			   {
				   s_tagged(min_row + 2, 0) = 0;
			   }
		else
		   {
			   if (s(min_row + 2, 0) == 2 and s(min_row + 1, 0) == 1)
			   {
				   s_tagged(min_row + 2, 0) = 1;
			   }
		   }
		if (neighbors(s, s_tagged, min_row, 0) >= 2)
					  {
					   s_tagged(min_row, 0) = 0;
					  }
		else
		   {
			   if (s(min_row, 0) == 2 and s(min_row - 1, 0) == 1)
			   {
				   s_tagged(min_row, 0) = 1;
			   }
		   }
		for (int ch = min_row; ch < leng(0); ch++)
			   {
				   if ((s(ch, 0) == 2 && s_tagged(ch, 0) == 1) or s(ch, 0) == 1)
				   {
					   lat_bond = g_lat_gdp;
				   }
				   if (s(ch, 0) == 2 && s_tagged(ch, 0) == 0)
				   {
					   lat_bond = g_lat;
				   }
				   if (s(ch, 0) > 0)
				   {
					   b(ch, 0) = neibors_lat(s, ch, 0) * lat_bond;
				   }
			   }
	}
	if (min_col < 12 && min_col > 0)
	{
		if (min_row > 0)
		{
			if (neighbors(s, s_tagged, min_row, min_col - 1) >= 2)
				   {
					   s_tagged(min_row, min_col - 1) = 0;
				   }
			else
			   {
				   if (s(min_row, min_col - 1) == 2 and s(min_row - 1, min_col - 1) == 1)
				   {
					   s_tagged(min_row, min_col - 1) = 1;
				   }
			   }
			if (neighbors(s, s_tagged, min_row - 1, min_col - 1) >= 2)
				   {
					   s_tagged(min_row - 1, min_col - 1) = 0;
				   }
			else
				   {
					   if (s(min_row - 1, min_col - 1) == 2 and s(min_row - 2, min_col - 1) == 1)
					   {
						   s_tagged(min_row - 1, min_col - 1) = 1;
					   }
				   }
			for (int ch = min_row - 1; ch < leng(min_col - 1); ch++) 
				   {
					   if ((s(ch, min_col - 1) == 2 && s_tagged(ch, min_col - 1) == 1) or s(ch, min_col - 1) == 1)
					   {
						   lat_bond = g_lat_gdp;
					   }
					   if (s(ch, min_col - 1) == 2 && s_tagged(ch, min_col - 1) == 0)
					   {
						   lat_bond = g_lat;
					   }
					   if (s(ch, min_col - 1) > 0)
					   {
						   b(ch, min_col - 1) = neibors_lat(s, ch, min_col - 1) * lat_bond;
					   }////????
				   }
			if (neighbors(s, s_tagged, min_row, min_col + 1) >= 2)
				   {
					   s_tagged(min_row, min_col + 1) = 0;
				   }
			else
				   {
					   if (s(min_row, min_col + 1) == 2 and s(min_row - 1, min_col + 1) == 1)
					   {
						   s_tagged(min_row, min_col + 1) = 1;
					   }
				   }
			if (neighbors(s, s_tagged, min_row - 1, min_col + 1) >= 2)
				   {
					   s_tagged(min_row - 1, min_col + 1) = 0;
				   }
			else
				   {
					   if (s(min_row - 1, min_col + 1) == 2 and s(min_row - 2, min_col + 1) == 1)
					   {
						   s_tagged(min_row - 1, min_col + 1) = 1;
					   }
				   }
			for (int ch = min_row - 1; ch < leng(min_col + 1); ch++)
			   {
				   if ((s(ch, min_col + 1) == 2 && s_tagged(ch, min_col + 1) == 1) or s(ch, min_col + 1) == 1)
				   {
					   lat_bond = g_lat_gdp;
				   }
				   if (s(ch, min_col + 1) == 2 && s_tagged(ch, min_col + 1) == 0)
				   {
					   lat_bond = g_lat;
				   }    	
				   if (s(ch, min_col + 1) > 0)
				   {
					   b(ch, min_col + 1) = neibors_lat(s, ch, min_col + 1) * lat_bond;
				   }
			   }
		}
	}
}

void event(int N, int M, umat& s, mat& b, umat& s_tagged, int x, vec& l, vec& t, vec& leng, double& tm)
{
  //cout << s << '\n';
  //cout << b << '\n';
  cube final_cube = cube(N, 13, 3);
  final_cube.fill(datum::inf); // cube for all probabilities to get minimum. 0 slice - assotiation, 1 - dissotiation, 2 - hydrolysis
  for (int c = 0; c < 13; c++) //on event probabilities
  {
	 final_cube(leng(c), c, 0) = -log(double(rand()/double(RAND_MAX)))/k_plus;
  }
  for (int c = 0; c < 13; c++) //off event probabilities
  {
	double energy = 0;
    for (int r = leng(c) - 1; r >= 0; r--)
    {  
    	double free_energy = energy + b(r, c); //go from top to bottom of a protofilament, summing lateral energies
    	energy = free_energy;
//      for (dimers_above = r; dimers_above < leng(c); dimers_above ++)
//      {
//    	  free_energy_change += b(dimers_above, c);
//      }
      double k_minus = (k_plus*100000)/exp (-(free_energy + g_long));    
	  final_cube(r, c, 1) = -log(double(rand()/double(RAND_MAX)))/k_minus;   
    }
  }
  for (int c = 0; c < 13; c++) //Hydrolysis probabilities
  {
	  for (int r = 0; r < leng(c) - 1; r++)
		{
		  if(s(r, c) == 2 && r > 9) 
				  {
				    final_cube(r, c, 2) = -log(double(rand()/double(RAND_MAX)))/(k_hydr);
				  }
		}
  }
      
  double cube_minimum;
  cube_minimum = final_cube.min();
  uword min_row, min_col, min_slice;
  final_cube.min(min_row, min_col, min_slice);
  if (min_slice == 1) // if dimer dissotiates
  {
	  
    for (int i = min_row; i < leng(min_col) + 1; i++)
    {
      s(i, min_col) = 0; //its position and all above turn to 0
      b(i, min_col) = 0; // bonds turn to 0
      s_tagged(i, min_col) = 0;//tags turn to 0
      
    }
    for (int nei = min_row; nei < leng(min_col); nei ++)
    {
    	change_neighbors_state(s, b, s_tagged, nei, min_col, leng); // change lateral bonds and tags of dimers of neighbor protofilaments
    }
    leng(min_col) = min_row; // change length of PF 
  }
  else if (min_slice == 0) // if dimer assotiates
  {
    s(min_row, min_col) = 2; 
    leng(min_col) += 1;
    if (neighbors(s, s_tagged, min_row, min_col) < 2 && s(min_row - 1, min_col) == 1 && min_row > 9)
    {
    	s_tagged(min_row, min_col) = 1; //check if it is tagged
    }
    if (s_tagged(min_row, min_col) == 1)
    {
    	b(min_row, min_col) = g_lat_gdp * neibors_lat(s, min_row, min_col); // keep it lateral bonds
    }
    if (s_tagged(min_row, min_col) == 0)
    {
    	b(min_row, min_col) = g_lat * neibors_lat(s, min_row, min_col);
    }
    change_neighbors_state(s, b, s_tagged, min_row, min_col, leng); // change lateral bonds and tags of dimers of neighbor PF's
    
  }
  else if (min_slice == 2) //if hydrolysis
    {
  	 s(min_row, min_col) = 1;
  	 s_tagged(min_row, min_col) = 0;
  	 if (neighbors(s, s_tagged, min_row + 1, min_col) < 2 and s(min_row + 1, min_col) == 2) // if above dimer have less than 2 GTP-neighbors we mark it as tagged
  	 {
  		 s_tagged(min_row + 1, min_col) = 1;
  		 for (int ch = min_row + 1; ch < leng(min_col); ch++) //check lateral bonds of all above dimers
  		 {
  			double lat_bond = 0;
  			if ((s(ch, min_col) == 2 && s_tagged(ch, min_col) == 1) or s(ch, min_col) == 1)
			   {
				   lat_bond = g_lat_gdp;
			   }
			   if (s(ch, min_col) == 2 && s_tagged(ch, min_col) == 0)
			   {
				   lat_bond = g_lat;
			   }    	
			   if (s(ch, min_col) > 0)
			   {
				   b(ch, min_col) = neibors_lat(s, ch, min_col) * lat_bond;
			   }
  		 }
  	 }
  	 
//  	 for (int ch = min_row; ch < min_row + 1; ch++)
//  	 {
//  		 double lat_bond;
//  		 if ((s(ch, min_col) == 2 && s_tagged(ch, min_col) == 1) or s(ch, min_col) == 1)
//  		 {
//  			 lat_bond = g_lat_gdp;
//  		 }
//  		 if (s(ch, min_col) == 2 && s_tagged(ch, min_col) == 0)
//  		 {
//  			 lat_bond = g_lat;
//  		 }
//  		 b(ch, min_col) = g_long + neibors_lat(s, ch, min_col) * g_lat;
//  	 }
  	change_neighbors_state(s, b, s_tagged, min_row, min_col, leng);
  	b(min_row, min_col) = neibors_lat(s, min_row, min_col) * g_lat_gdp;
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
 l(x) = (double)(non_zer / 13.0f - 10.0f)*0.008;
 cout << tm << "," << l(x) << "\n";
// ofstream outfile;  // file to write cap size
// outfile.open("cap.txt",  std::ios_base::app); 
// double cap = 0;
// for (int c = 0; c < 13; c++) // cap
//   {
// 	  for (int r = 0; r < leng(c); r++)
// 		{
// 		  if(s(r, c) == 2 && r > 9) 
// 				  {
// 			  	  	  cap += 1;
// 				  }
// 		}
//   }
// outfile << tm << "," << cap << "\n";
// outfile.close();
// if (x > 396074 and x < 397518)
// {
//	 std::string numb = std::to_string(x + 1);
//   ofstream outfile;
////   std::string matrixname = "/home/marge/mtmatlab2/matrix" + numb + ".txt";
//   std::string bondsname = "/home/marge/mtmatlab3/bonds" + numb + ".txt";
//   std::string matrsname = "/home/marge/mtmatlab3/tagged" + numb + ".txt";
////   outfile.open(matrixname,  std::ios_base::app);
////   for (int r = 1040; r < 1075; r++)
////   {
////	   for (int c = 0; c < 12; c++)
////	   {
////		   outfile << s(r, c) << ",";
////	   }
////	   outfile << s(r, 12) << "\n";
////   }
////   for (int c = 0; c < 12; c++)
////   	   {
////   		   outfile << s(1075, c) << ",";
////   	   }
////   outfile << s(1075, 12);
////   outfile.close();
//   outfile.open(bondsname,  std::ios_base::app);
//   for (int r = 1605; r > 1565; r--)
//      {
//   	   for (int c = 0; c < 12; c++)
//   	   {
//   		   outfile << b(r, c) << ",";
//   	   }
//   	   outfile << b(r, 12) << "\n";
//      }
//      for (int c = 0; c < 12; c++)
//      	   {
//      		   outfile << b(1565, c) << ",";
//      	   }
//      outfile << b(1565, 12);
//   outfile.close();
//   outfile.open(matrsname,  std::ios_base::app);
//   for (int r = 1605; r > 1565; r--)
//      {
//   	   for (int c = 0; c < 12; c++)
//   	   {
// 	 	 	   if (s_tagged(r, c) == 1)
// 	 	 	   {
// 	 	 		   outfile << "2,";
// 	 	 	   }
// 	 	 	   if (s_tagged(r, c) == 0 and s(r, c) == 2)
// 	 	 	   {
// 	 	 		   outfile << "3,";
// 	 	 	   }
// 	 	 	   if (s_tagged(r, c) == 0 and s(r, c) == 1)
// 	 	 	   {
// 	 	 		   outfile << "1,";
// 	 	 	   }
// 	 	 	   if (s(r, c) == 0)
// 	 	 	   {
// 	 	 		   outfile << "0,";
// 	 	 	   }
//////   		   outfile << s_tagged(r, c) << ",";
//   	   }
// 	 	 	if (s_tagged(r, 12) == 1)
//			   {
//				   outfile << "2" << "\n";
//			   }
//			   if (s_tagged(r, 12) == 0 and s(r, 12) == 2)
//			   {
//				   outfile << "3" << "\n";
//			   }
//			   if (s_tagged(r, 12) == 0 and s(r, 12) == 1)
//			   {
//				   outfile << "1" << "\n";
//			   }
//			   if (s(r, 12) == 0)
//			   {
//				   outfile << "0" << "\n";
//			   }
//////   	   outfile << s_tagged(r, 12) << "\n";
//      }
//      for (int c = 0; c < 12; c++)
//      	   {
//			   if (s_tagged(1565, c) == 1)
//				   {
//					   outfile << "2,";
//				   }
//				   if (s_tagged(1565, c) == 0 and s(1565, c) == 2)
//				   {
//					   outfile << "3,";
//				   }
//				   if (s_tagged(1565, c) == 0 and s(1565, c) == 1)
//				   {
//					   outfile << "1,";
//				   }
//				   if (s(1565, c) == 0)
//				   {
//					   outfile << "0,";
//				   }
////      	//	   outfile << s_tagged(1055, c) << ",";
//      	   }
//      if (s_tagged(1565, 12) == 1)
//		   {
//			   outfile << "2";
//		   }
//		   if (s_tagged(1565, 12) == 0 and s(1565, 12) == 2)
//		   {
//			   outfile << "3";
//		   }
//		   if (s_tagged(1565, 12) == 0 and s(1565, 12) == 1)
//		   {
//			   outfile << "1";
//		   }
//		   if (s(1565, 12) == 0)
//		   {
//			   outfile << "0";
//		   }
//////      outfile << s_tagged(1055, 12);
//   outfile.close();
// }

// cout << min_row<< " " <а<  min_col << " " << min_slice <<"\n";
//  std::string numb = std::to_string(x + 1);
//
  
//  outfile.open(bondsname,  std::ios_base::app);
//  outfile << b << '\n';
//  outfile.close();
 
//  cout << b << '\n';
//  cout << s_tagged << '\n';
		  
}


int main(int argc, char** argv)
{
  cout << "time,length" << "\n";
  struct timeb tp;
  ftime(&tp);
  srand(static_cast<unsigned int>(getpid()) ^
  static_cast<unsigned int>(pthread_self()) ^
  static_cast<unsigned int >(tp.millitm));
//  srand(883634);
  int N = atoi(argv[1]);
  int M = atoi(argv[2]);
  vec l = vec(M, fill::zeros);
  vec t = vec(M, fill::zeros);
  vec leng = vec(13, fill::zeros);
  for (int a = 0; a < 13; a++)
  {
	  leng(a) = 10;
  }
  double tm = 0.0f;
  umat s_in = umat(N, 13, fill::zeros);
  mat b_in = mat(N, 13, fill::zeros);
  umat s_tagged = umat(N, 13, fill::zeros);
  for (int i = 0; i < 13; i++)
  {
    for (int n = 0; n < 10; n++)
    {
      s_in(n, i) = 2;
    }
  }
//  for (int i = 0; i < 13; i++)
//  {
//    for (int n = 10; n < 100; n++)
//    {
//      s_in(n, i) = 1;
//    }
//  }
  for (int i = 0; i < 13; i++)
      {
	  for (int n = 0; n < 10; n++)
	      {
		    b_in(n, i) = - 11.4;
	      }
      }
//  for (int i = 0; i < 13; i++)
//        {
//  	  for (int n = 10; n < 100; n++)
//  	      {
//  		    b_in(n, i) = - 6.4;
//  	      }
//        }
  b_in(0, 0) += 5.7;
  b_in(8, 12) += 2.85;
  b_in(9, 12) += 5.7;
  s_tagged(9, 12) = 1;
  int x;
  for (x = 0; x < M; x++)
  {
	  event(N, M, s_in, b_in, s_tagged, x, l, t, leng, tm);
  }
  return 0;
}





