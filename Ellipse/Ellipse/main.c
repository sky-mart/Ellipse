#include <stdio.h>
#include "tests.h"

int main()
{
	//blasTest(1);
	//gesvTest(1);
	//cbrtTest(1);
	//equTest(1);
	//distToEllipseTest(1);
	//ellipseFittingTest(1);
	ellipseFittingCircleTest(1);
	ellipseFittingNoiseTest(1);
	printf("\nPress ENTER key to continue...\n");
	getchar();

	return 0;
}