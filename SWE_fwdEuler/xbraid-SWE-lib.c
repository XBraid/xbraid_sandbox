#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

#ifdef M_PI
#define PI M_PI
#else
#define PI 3.14159265358979
#endif


/* Print accumulated info on space-time grids visited during the simulation */
void
print_sc_info(double *sc_info, 
              int     max_levels)
{
   int i;
   double dx, dt; 

   printf("\n-----------------------------------------------------------------\n"); 
   printf("-----------------------------------------------------------------\n\n"); 
   printf( " Per level diagnostic information \n\n");
   printf("level       dx          dt        dt/dx\n"); 
   printf("-----------------------------------------------------------------\n"); 
   for( i = 0; i < max_levels; i++)
   {
      dx = sc_info[i*2];
      dt = sc_info[i*2+1];
      if (dx == -1){
         break;
      }
      printf(" %2d   |   %1.2e    %1.2e    %1.2e\n", i, dx, dt, dt/(dx) ); 
   }
   printf( "\n" );
}


