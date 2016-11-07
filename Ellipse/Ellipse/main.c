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
	ellipseFittingTestWrapper(30, -1.0, -1.0, 1.0, 6.0, 0.8, 0, 1e-3, 1);
	//ellipseFittingNoiseTest(1);
	//ellipseFittingLoopTest(1);
	printf("\nPress ENTER key to continue...\n");
	getchar();

	return 0;
}