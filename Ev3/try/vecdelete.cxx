#include <iostream>
#include <vector>
int main(int argc, char** argv) {
	vector<int> vi;
	vi.push_back(1);
	vi.push_back(2);
	vi.push_back(3);
	vi.erase(vi.begin() + 0);
	vi.erase(vi.begin() + 1);
	cout << vi.size() << endl;
	cout << vi[0] << endl;
/*	cout << vi[1] << endl;
	cout << vi[2] << endl; */
	return 0;
}
