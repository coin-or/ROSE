#include <iostream>
#include <string>
#include <strstream>

int main(int argc, char** argv) {
	char r[100];
	ostrstream outr(r, sizeof(r));
	outr << 2.3 << ends;
	cout <<	r << endl;
	return 0;
}
