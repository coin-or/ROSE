#include <iostream>
#include <cstring>

using namespace std;

int main(int argc, char** argv) {

	char buffer[1024];
	cin >> buffer;
	char* nt;
	char* t = buffer;
	nt = strtok(t, " ");
	//cout << t << "; " << nt - t << endl;
	while(nt) {
		cout << nt << endl;
		nt = strtok(0, " ");
	}
	return 0;	
}
