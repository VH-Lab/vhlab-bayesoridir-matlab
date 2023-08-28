#include "mystats.h"
#include <stdio.h>


int main() {

	double y1 = normpdf(0,0,1);
	double y2 = normpdf(1,0,1);

	printf("%f %f\n", y1, y2);
}

