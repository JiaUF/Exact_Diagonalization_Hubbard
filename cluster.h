#ifndef _H_cluster
#define _H_cluster

#include <vector>

using namespace std;

class Cluster{
  private:
   const static int Num=20; 
   int nghb[Num][4], xn[Num], yn[Num];
   double xph[Num], yph[Num];
   double x[Num], y[Num]; //coordinates
   int xx[Num], yy[Num];
   
  public:
   //constructor
   Cluster(int n, int dim){
     int N=n;
     if(dim==2){
      switch(N){
         case 4:
            /*       
                   0 --- 1
                   |     |
                   |     |
                   2 --- 3
            */
            // [][0] up, [][1] dn, [][2] right, [][3] left
            nghb[0][0]=2; nghb[0][1]=2; nghb[0][2]=1; nghb[0][3]=1;
            nghb[1][0]=3; nghb[1][1]=3; nghb[1][2]=0; nghb[1][3]=0;
            nghb[2][0]=0; nghb[2][1]=0; nghb[2][2]=3; nghb[2][3]=3; 
            nghb[3][0]=1; nghb[3][1]=1; nghb[3][2]=2; nghb[3][3]=2; 
          
            //neighbor in the x direction : xn
            //neighbor in the x direction : yn
            xn[0]=1; xn[1]=0; xn[2]=3; xn[3]=2;           
            yn[0]=2; yn[1]=3; yn[2]=0; yn[3]=1;           

            break;

         case 8:
            /*           
                      0 --- 1       
                      |     |
                      |     |
                2 --- 3 --- 4 --- 5 --- 2 --- 3 --- 4     
                      |     |     |     |     |     |
                      |     |     |     |     |     |
                      6 --- 7 --- 0 --- 1 --- 6 --- 7
                            |     |     |     |
                            |     |     |     |
                            2 --- 3 --- 4 --- 5
            */
            // [][0] up, [][1] dn, [][2] right, [][3] left
            nghb[0][0]=5; nghb[0][1]=3; nghb[0][2]=1; nghb[0][3]=7;
            nghb[1][0]=2; nghb[1][1]=4; nghb[1][2]=6; nghb[1][3]=0;
            nghb[2][0]=7; nghb[2][1]=1; nghb[2][2]=3; nghb[2][3]=5; 
            nghb[3][0]=0; nghb[3][1]=6; nghb[3][2]=4; nghb[3][3]=2; 
            nghb[4][0]=1; nghb[4][1]=7; nghb[4][2]=5; nghb[4][3]=3;
            nghb[5][0]=6; nghb[5][1]=0; nghb[5][2]=2; nghb[5][3]=4;
            nghb[6][0]=3; nghb[6][1]=5; nghb[6][2]=7; nghb[6][3]=1; 
            nghb[7][0]=4; nghb[7][1]=2; nghb[7][2]=0; nghb[7][3]=6; 
            
            xn[0]=1; xn[1]=6; xn[2]=3; xn[3]=4; 
            xn[4]=5; xn[5]=2; xn[6]=7; xn[7]=0;
            yn[0]=5; yn[1]=2; yn[2]=7; yn[3]=0;
            yn[4]=1; yn[5]=6; yn[6]=3; yn[7]=4;

            xph[0]=1.0;  xph[1]=-1.0; xph[2]=1.0;  xph[3]=1.0;
            xph[4]=1.0;  xph[5]=-1.0; xph[6]=1.0;  xph[7]=-1.0;
            yph[0]=-1.0; yph[1]=-1.0; yph[2]=-1.0; yph[3]=1.0;
            yph[4]=1.0;  yph[5]=-1.0; yph[6]=1.0;  yph[7]=1.0;
            
            x[0] = 0.0;  y[0] = 0.0;
            x[1] = 0.25; y[1] = 0.25;
            x[2] = 0.0;  y[2] = 0.5;
            x[3] = 0.25; y[3] = 0.75;
            x[4] = 0.5;  y[4] = 0.0;
            x[5] = 0.75; y[5] = 0.25;
            x[6] = 0.5;  y[6] = 0.5;
            x[7] = 0.75; y[7] = 0.75;

            xx[0] = 0;  yy[0] = 0;
            xx[1] = 1;  yy[1] = 0;
            xx[2] = 1;  yy[2] = 1;
            xx[3] = 2;  yy[3] = 1;
            xx[4] = 3;  yy[4] = 1;
            xx[5] = 2;  yy[5] = -1;
            xx[6] = 2;  yy[6] = 0;
            xx[7] = 3;  yy[7] = 0;

            break;

         case 10:
           /*           0 --- 1 --- 2 --- 3 --- 4
                        |     |     |     |
                        |     |     |     |
                        7 --- 8 --- 9 --- 0 --- 1                     7    8    9
                        |     |     |     |
                        |     |     |     |
                        4 --- 5 --- 6  -- 7 --- 8                     4    5    6
                        |     |     |     |
                        |     |     |     |
                  0 --- 1 --- 2 --- 3 --- 4 --- 5               0     1    2    3
                        |     |     |     |
                        |     |     |     |
                        8 --- 9 --- 0 --- 1 --- 2

           */ 
            // [][0] up, [][1] dn, [][2] right, [][3] left
            nghb[0][0]=3;  nghb[0][1]=7; nghb[0][2]=1;  nghb[0][3]=9;  
            nghb[1][0]=4;  nghb[1][1]=8; nghb[1][2]=2;  nghb[1][3]=0;  
            nghb[2][0]=5;  nghb[2][1]=9; nghb[2][2]=3;  nghb[2][3]=1;  
            nghb[3][0]=6;  nghb[3][1]=0; nghb[3][2]=4;  nghb[3][3]=2;  
            nghb[4][0]=7;  nghb[4][1]=1; nghb[4][2]=5;  nghb[4][3]=3;  
            nghb[5][0]=8;  nghb[5][1]=2; nghb[5][2]=6;  nghb[5][3]=4;  
            nghb[6][0]=9;  nghb[6][1]=3; nghb[6][2]=7;  nghb[6][3]=5;  
            nghb[7][0]=0;  nghb[7][1]=4; nghb[7][2]=8;  nghb[7][3]=6;  
            nghb[8][0]=1;  nghb[8][1]=5; nghb[8][2]=9;  nghb[8][3]=7; 
            nghb[9][0]=2;  nghb[9][1]=6; nghb[9][2]=0;  nghb[9][3]=8;  

            xn[0]=1; xn[1]=2; xn[2]=3; xn[3]=4; 
            xn[4]=5; xn[5]=6; xn[6]=7; xn[7]=8; 
            xn[8]=9; xn[9]=0; 
            yn[0]=3; yn[1]=4; yn[2]=5; yn[3]=6; 
            yn[4]=7; yn[5]=8; yn[6]=9; yn[7]=0; 
            yn[8]=1; yn[9]=2;  
            
            xph[0]=1.0;  xph[1]=1.0;  xph[2]=1.0;  xph[3]=-1.0;
            xph[4]=1.0;  xph[5]=1.0;  xph[6]=-1.0; xph[7]=1.0;
            xph[8]=1.0;  xph[9]=-1.0; 
            yph[0]=-1.0; yph[1]=1.0;  yph[2]=1.0;  yph[3]=1.0;
            yph[4]=1.0;  yph[5]=1.0;  yph[6]=1.0;  yph[7]=-1.0;
            yph[8]=-1.0; yph[9]=-1.0; 
            
            x[0] = 0.0;  y[0] = 0.0;
            x[1] = 0.3;  y[1] = 0.1;
            x[2] = 0.6;  y[2] = 0.2;
            x[3] = 0.9;  y[3] = 0.3;
            x[4] = 0.2;  y[4] = 0.4;
            x[5] = 0.5;  y[5] = 0.5;
            x[6] = 0.8;  y[6] = 0.6;
            x[7] = 0.1;  y[7] = 0.7;
            x[8] = 0.4;  y[8] = 0.8;
            x[9] = 0.7;  y[9] = 0.9;

            xx[0] = 0;  yy[0] = 0;
            xx[1] = 1;  yy[1] = 0;
            xx[2] = 2;  yy[2] = 0;
            xx[3] = 3;  yy[3] = 0;
            xx[4] = 1;  yy[4] = 1;
            xx[5] = 2;  yy[5] = 1;
            xx[6] = 3;  yy[6] = 1;
            xx[7] = 1;  yy[7] = 2;
            xx[8] = 2;  yy[8] = 2;
            xx[9] = 3;  yy[9] = 2;

            break;                       

         case 12:
           /*           0 --- 1 --- 2 --- 3 --- 0
                        |     |     |     |
                        |     |     |     |
                        8 --- 9 --- 10--- 11--- 8              8    9   10   11
                        |     |     |     |
                        |     |     |     |
                        4 --- 5 --- 6  -- 7 --- 4              4    5    6    7
                        |     |     |     |
                        |     |     |     |
                  0 --- 1 --- 2 --- 3 --- 0 --- 1         0    1    2    3
                        |     |     |     |
                        |     |     |     |
                        9 --- 10--- 11--- 8 --- 9

           */ 
            // [][0] up, [][1] dn, [][2] right, [][3] left
            nghb[0][0]=7;  nghb[0][1]=8; nghb[0][2]=1;  nghb[0][3]=3;  
            nghb[1][0]=4;  nghb[1][1]=9; nghb[1][2]=2;  nghb[1][3]=0;  
            nghb[2][0]=5;  nghb[2][1]=10;nghb[2][2]=3;  nghb[2][3]=1;  
            nghb[3][0]=6;  nghb[3][1]=11;nghb[3][2]=0;  nghb[3][3]=2;  
            nghb[4][0]=8;  nghb[4][1]=1; nghb[4][2]=5;  nghb[4][3]=7;  
            nghb[5][0]=9;  nghb[5][1]=2; nghb[5][2]=6;  nghb[5][3]=4;  
            nghb[6][0]=10; nghb[6][1]=3; nghb[6][2]=7;  nghb[6][3]=5;  
            nghb[7][0]=11; nghb[7][1]=0; nghb[7][2]=4;  nghb[7][3]=6;  
            nghb[8][0]=0;  nghb[8][1]=4; nghb[8][2]=9;  nghb[8][3]=11; 
            nghb[9][0]=1;  nghb[9][1]=5; nghb[9][2]=10; nghb[9][3]=8;  
            nghb[10][0]=2; nghb[10][1]=6;nghb[10][2]=11;nghb[10][3]=9; 
            nghb[11][0]=3; nghb[11][1]=7;nghb[11][2]=8; nghb[11][3]=10;

            xn[0]=1; xn[1]=2; xn[2]=3; xn[3]=0; 
            xn[4]=5; xn[5]=6; xn[6]=7; xn[7]=4; 
            xn[8]=9; xn[9]=10; xn[10]=11; xn[11]=8; 
            yn[0]=7; yn[1]=4; yn[2]=5; yn[3]=6; 
            yn[4]=8; yn[5]=9; yn[6]=10; yn[7]=11; 
            yn[8]=0; yn[9]=1; yn[10]=2; yn[11]=3; 

            xph[0]=1.0;  xph[1]=1.0;  xph[2]=1.0;  xph[3]=-1.0;
            xph[4]=1.0;  xph[5]=1.0;  xph[6]=1.0;  xph[7]=-1.0;
            xph[8]=1.0;  xph[9]=1.0;  xph[10]=1.0; xph[11]=-1.0; 
            yph[0]=-1.0; yph[1]=1.0;  yph[2]=1.0;  yph[3]=1.0;
            yph[4]=1.0;  yph[5]=1.0;  yph[6]=1.0;  yph[7]=1.0;
            yph[8]=-1.0; yph[9]=-1.0; yph[10]=-1.0;yph[11]=-1.0;
            
            
            x[0] = 0.0;  y[0] = 0.0;
            x[1] = 0.25; y[1] = 0.0;
            x[2] = 0.5;  y[2] = 0.0;
            x[3] = 0.75; y[3] = 0.0;
            x[4] = 0.167;  y[4] = 0.333;
            x[5] = 0.417;  y[5] = 0.333;
            x[6] = 0.567;  y[6] = 0.333;
            x[7] = 0.917;  y[7] = 0.333;
            x[8] = 0.083;  y[8] = 0.667;
            x[9] = 0.333;  y[9] = 0.667;
            x[10] = 0.583; y[10] = 0.667;
            x[11] = 0.833; y[11] = 0.667;

            xx[0] = 0;  yy[0] = 0;
            xx[1] = 1;  yy[1] = 0;
            xx[2] = 2;  yy[2] = 0;
            xx[3] = 3;  yy[3] = 0;
            xx[4] = 1;  yy[4] = 1;
            xx[5] = 2;  yy[5] = 1;
            xx[6] = 3;  yy[6] = 1;
            xx[7] = 4;  yy[7] = 1;
            xx[8] = 1;  yy[8] = 2;
            xx[9] = 2;  yy[9] = 2;
            xx[10] = 3;  yy[10] = 2;
            xx[11] = 4;  yy[11] = 2;

            break;                       


 	 case 16:
                    
           /*           0 --- 1 --- 2 --- 3 --- 0
                        |     |     |     |
                        |     |     |     |
                        4 --- 5 --- 6 --- 7 --- 4   
                        |     |     |     |
                        |     |     |     |
                        8 --- 9 --- 10 -- 11--- 8
                        |     |     |     |
                        |     |     |     |
                        12--- 13--- 14--- 15--- 12
                        |     |     |     |
                        |     |     |     |
                        0 --- 1 --- 2 --- 3 --- 0

           */ 
            // [][0] up, [][1] dn, [][2] right, [][3] left
            nghb[0][0]=12; nghb[0][1]=4;   nghb[0][2]=1;   nghb[0][3]=3;
            nghb[1][0]=13; nghb[1][1]=5;   nghb[1][2]=2;   nghb[1][3]=0;
            nghb[2][0]=14; nghb[2][1]=6;   nghb[2][2]=3;   nghb[2][3]=1;
            nghb[3][0]=15; nghb[3][1]=7;   nghb[3][2]=0;   nghb[3][3]=2;
            nghb[4][0]=0;  nghb[4][1]=8;   nghb[4][2]=5;   nghb[4][3]=7;
            nghb[5][0]=1;  nghb[5][1]=9;   nghb[5][2]=6;   nghb[5][3]=4;
            nghb[6][0]=2;  nghb[6][1]=10;  nghb[6][2]=7;   nghb[6][3]=5;
            nghb[7][0]=3;  nghb[7][1]=11;  nghb[7][2]=4;   nghb[7][3]=6;
            nghb[8][0]=4;  nghb[8][1]=12;  nghb[8][2]=9;   nghb[8][3]=11;
            nghb[9][0]=5;  nghb[9][1]=13;  nghb[9][2]=10;  nghb[9][3]=8;
            nghb[10][0]=6; nghb[10][1]=14; nghb[10][2]=11; nghb[10][3]=9;
            nghb[11][0]=7; nghb[11][1]=15; nghb[11][2]=8;  nghb[11][3]=10;
            nghb[12][0]=8; nghb[12][1]=0;  nghb[12][2]=13; nghb[12][3]=15;
            nghb[13][0]=9; nghb[13][1]=1;  nghb[13][2]=14; nghb[13][3]=12;
            nghb[14][0]=10; nghb[14][1]=2; nghb[14][2]=15; nghb[14][3]=13;
            nghb[15][0]=11; nghb[15][1]=3; nghb[15][2]=12; nghb[15][3]=14;

            xn[0]=1; xn[1]=2; xn[2]=3; xn[3]=0; 
            xn[4]=5; xn[5]=6; xn[6]=7; xn[7]=4; 
            xn[8]=9; xn[9]=10; xn[10]=11; xn[11]=8; 
            xn[12]=13; xn[13]=14; xn[14]=15; xn[15]=12; 
            yn[0]=12; yn[1]=13; yn[2]=14; yn[3]=15; 
            yn[4]=0; yn[5]=1; yn[6]=2; yn[7]=3; 
            yn[8]=4; yn[9]=5; yn[10]=6; yn[11]=7; 
            yn[12]=8; yn[13]=9; yn[14]=10; yn[15]=11; 

            break;
         default:
            break;
      }
     } else if(dim==1){
      switch(N){
         case 4:
            /*
                   0 --- 1 --- 2 --- 3
            */
            // [][0] up, [][1] dn, [][2] right, [][3] left
            nghb[0][0]=-1; nghb[0][1]=-1; nghb[0][2]=1; nghb[0][3]=3;
            nghb[1][0]=-1; nghb[1][1]=-1; nghb[1][2]=2; nghb[1][3]=0;
            nghb[2][0]=-1; nghb[2][1]=-1; nghb[2][2]=3; nghb[2][3]=1;
            nghb[3][0]=-1; nghb[3][1]=-1; nghb[3][2]=0; nghb[3][3]=2;

            //neighbor in the x direction : xn
            //neighbor in the y direction : yn
            xn[0]=1; xn[1]=2; xn[2]=3; xn[3]=0;
            yn[0]=-1; yn[1]=-1; yn[2]=-1; yn[3]=-1;

            xph[0]=1.0;  xph[1]=1.0; xph[2]=1.0;  xph[3]=1.0;
            yph[0]=1.0;  yph[1]=1.0; yph[2]=1.0; yph[3]=1.0;

            break;

         case 6:
            /*
                   0 --- 1 --- 2 --- 3 --- 4 --- 5
            */
            // [][0] up, [][1] dn, [][2] right, [][3] left
            nghb[0][0]=-1; nghb[0][1]=-1; nghb[0][2]=1; nghb[0][3]=5;
            nghb[1][0]=-1; nghb[1][1]=-1; nghb[1][2]=2; nghb[1][3]=0;
            nghb[2][0]=-1; nghb[2][1]=-1; nghb[2][2]=3; nghb[2][3]=1;
            nghb[3][0]=-1; nghb[3][1]=-1; nghb[3][2]=4; nghb[3][3]=2;
            nghb[4][0]=-1; nghb[4][1]=-1; nghb[4][2]=5;  nghb[4][3]=3;  
            nghb[5][0]=-1; nghb[5][1]=-1; nghb[5][2]=0;  nghb[5][3]=4;  

            //neighbor in the x direction : xn
            //neighbor in the y direction : yn
            xn[0]=1; xn[1]=2; xn[2]=3; xn[3]=4; xn[4]=5; xn[5]=0;
            yn[0]=-1; yn[1]=-1; yn[2]=-1; yn[3]=-1; yn[4]=-1; yn[5]=-1;

            xph[0]=1.0;  xph[1]=1.0; xph[2]=1.0;  xph[3]=1.0;
            xph[4]=1.0;  xph[5]=1.0; 
            yph[0]=1.0;  yph[1]=1.0; yph[2]=1.0; yph[3]=1.0;
            yph[4]=1.0;  yph[5]=1.0; 

            break;

         case 8:
            /*
                   0 --- 1 --- 2 --- 3 --- 4 --- 5 --- 6 --- 7
            */
            // [][0] up, [][1] dn, [][2] right, [][3] left
            nghb[0][0]=-1; nghb[0][1]=-1; nghb[0][2]=1; nghb[0][3]=7;
            nghb[1][0]=-1; nghb[1][1]=-1; nghb[1][2]=2; nghb[1][3]=0;
            nghb[2][0]=-1; nghb[2][1]=-1; nghb[2][2]=3; nghb[2][3]=1;
            nghb[3][0]=-1; nghb[3][1]=-1; nghb[3][2]=4; nghb[3][3]=2;
            nghb[4][0]=-1; nghb[4][1]=-1; nghb[4][2]=5;  nghb[4][3]=3;  
            nghb[5][0]=-1; nghb[5][1]=-1; nghb[5][2]=6;  nghb[5][3]=4;  
            nghb[6][0]=-1;  nghb[6][1]=-1;  nghb[6][2]=7;   nghb[6][3]=5;
            nghb[7][0]=-1;  nghb[7][1]=-1;  nghb[7][2]=0;   nghb[7][3]=6;

            //neighbor in the x direction : xn
            //neighbor in the y direction : yn
            xn[0]=1; xn[1]=2; xn[2]=3; xn[3]=4; xn[4]=5; xn[5]=6; xn[6]=7; xn[7]=0;
            yn[0]=-1; yn[1]=-1; yn[2]=-1; yn[3]=-1; yn[4]=-1; yn[5]=-1;yn[6]=-1;yn[7]=-1;

            xph[0]=1.0;  xph[1]=1.0; xph[2]=1.0;  xph[3]=1.0;
            xph[4]=1.0;  xph[5]=1.0; xph[6]=1.0;  xph[7]=1.0;
            yph[0]=1.0;  yph[1]=1.0; yph[2]=1.0; yph[3]=1.0;
            yph[4]=1.0;  yph[5]=1.0; yph[6]=1.0;  yph[7]=1.0;

            break;

         case 10:
            /*
                   0 --- 1 --- 2 --- 3 --- 4 --- 5 --- 6 --- 7 --- 8 --- 9
            */
            // [][0] up, [][1] dn, [][2] right, [][3] left
            nghb[0][0]=-1; nghb[0][1]=-1; nghb[0][2]=1; nghb[0][3]=9;
            nghb[1][0]=-1; nghb[1][1]=-1; nghb[1][2]=2; nghb[1][3]=0;
            nghb[2][0]=-1; nghb[2][1]=-1; nghb[2][2]=3; nghb[2][3]=1;
            nghb[3][0]=-1; nghb[3][1]=-1; nghb[3][2]=4; nghb[3][3]=2;
            nghb[4][0]=-1; nghb[4][1]=-1; nghb[4][2]=5;  nghb[4][3]=3;  
            nghb[5][0]=-1; nghb[5][1]=-1; nghb[5][2]=6;  nghb[5][3]=4;  
            nghb[6][0]=-1;  nghb[6][1]=-1;  nghb[6][2]=7;   nghb[6][3]=5;
            nghb[7][0]=-1;  nghb[7][1]=-1;  nghb[7][2]=8;   nghb[7][3]=6;
            nghb[8][0]=-1;  nghb[8][1]=-1;  nghb[8][2]=9;   nghb[8][3]=7;
            nghb[9][0]=-1;  nghb[9][1]=-1;  nghb[9][2]=0;  nghb[9][3]=8;

            //neighbor in the x direction : xn
            //neighbor in the y direction : yn
            xn[0]=1; xn[1]=2; xn[2]=3; xn[3]=4; xn[4]=5; xn[5]=6; xn[6]=7; xn[7]=8; xn[8]=9;xn[9]=0;
            yn[0]=-1; yn[1]=-1; yn[2]=-1; yn[3]=-1; yn[4]=-1; yn[5]=-1;yn[6]=-1;yn[7]=-1;yn[8]=-1;yn[9]=-1;

            xph[0]=1.0;  xph[1]=1.0; xph[2]=1.0;  xph[3]=1.0;
            xph[4]=1.0;  xph[5]=1.0; xph[6]=1.0;  xph[7]=1.0;
	    xph[8]=1.0;  xph[9]=1.0;
            yph[0]=1.0;  yph[1]=1.0; yph[2]=1.0; yph[3]=1.0;
            yph[4]=1.0;  yph[5]=1.0; yph[6]=1.0;  yph[7]=1.0;
	    yph[8]=1.0;  yph[9]=1.0;

            break;
         default:
            break;

      }
     }
   }

   vector<int> neighbor(int n){
      vector<int> il(4);
      vector<int>::iterator itr = il.begin();
      while(itr != il.end()){
         *itr = nghb[n][itr-il.begin()];
         itr++;
      }
      return il;
   }

   vector<int> nnneighbor(int n){
      vector<int> il(4);
      vector<int>::iterator itr = il.begin();
      *itr = nghb[nghb[n][2]][0];
      itr++;
      *itr = nghb[nghb[n][2]][1];
      itr++;
      *itr = nghb[nghb[n][3]][0];
      itr++;
      *itr = nghb[nghb[n][3]][1];
      return il;
   }

   vector<int> nnnneighbor(int n){
      vector<int> il(4);
      vector<int>::iterator itr = il.begin();
      *itr = nghb[nghb[n][2]][2];
      itr++;
      *itr = nghb[nghb[n][1]][1];
      itr++;
      *itr = nghb[nghb[n][0]][0];
      itr++;
      *itr = nghb[nghb[n][3]][3];
      return il;
   }

   int xneighbor(int n){ return xn[n];}
   int yneighbor(int n){ return yn[n];}
   double xphase(int n){ return xph[n];}
   double yphase(int n){ return yph[n];}
   int xcoord(int n) {return xx[n];}
   int ycoord(int n) {return yy[n];}
};

#endif
