#include <iostream>
#include <strstream>
#include <cstdio>
#include <cstdlib>
int main(int argc, char** argv) {
	if (argc == 1) {
		cerr << "give me a string, any string [on the command line]\n";
		exit(1);
	}
	istrstream* input = new istrstream(argv[1]);
	char ch;
	while (input->get(ch)) {
		cout << ch;
	}
	cout << endl;
	return 0;
}

