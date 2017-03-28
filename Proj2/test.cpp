#include <pthread.h>
#include <stdio.h>
#include <iostream>
#include <stdlib.h>

using namespace std;

#define NUM_THREADS     4

int main (int argc, char *argv[])
{
	for(int resX = 0 ; resX  < 100 ; resX++)
		cout << min(NUM_THREADS, resX - 2) << endl;
	return	1;
}
