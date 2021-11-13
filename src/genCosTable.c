#include <stdio.h>
#include <stdlib.h>
#include <math.h>

main()
{
    int ip, ic;

    for (ip = 0; ip < 8; ip++)
    {
        for (ic = 0; ic < 8; ic ++)
        {
            printf("%30.27lf\n",
                   cos( ((double) (((2*ip)+1)*ic)) * ((double) M_PI) / ((double) 16.0) ));
        }
    }
}

