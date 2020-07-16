
/*
* -----------------------------------------------------------------
*  Seed Generator --- seed-gen.c
*  Version: 1.0
*  Date: Nov 16, 2010
* ----------------------------------------------------------------- 
*  Programmer: Americo Barbosa da Cunha Junior
*              americo.cunhajr@gmail.com
* -----------------------------------------------------------------
*  Copyright (c) 2010 by Americo Barbosa da Cunha Junior
*
*  This program is free software: you can redistribute it and/or
*  modify it under the terms of the GNU General Public License as
*  published by the Free Software Foundation, either version 3 of
*  the License, or (at your option) any later version.
*
*  This program is distributed in the hope that it will be useful,
*  but WITHOUT ANY WARRANTY; without even the implied warranty of
*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*  GNU General Public License for more details.
*
*  A copy of the GNU General Public License is available in
*  LICENSE.txt or http://www.gnu.org/licenses/.
* -----------------------------------------------------------------
*  This program generates seeds for a randon number generator.
* -----------------------------------------------------------------
*/




#include <time.h>
#include <stdio.h>





int main(void)
{
    unsigned long int seed;
    
    seed = time(NULL);
    printf("\n seed: %lu\n", seed);
    
    return 0;
}
