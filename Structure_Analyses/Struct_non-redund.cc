#include <iostream>
#include <fstream>
#include <vector>
#include <unordered_set>
#include <algorithm>
#include <cctype>
#include <sstream>

using namespace std;


unordered_set<string> read_matrix (string file)
{
	ifstream f (file);
	unordered_set<string> sp;
	string s;
	
	while (getline (f,s) && (s[0]=='#')) { }
	
	while (f >> s){
		s.pop_back();
		sp.insert (s);
	}
	return sp;
}

string idupper (string l)
{
	stringstream ss (l);
	
	while (getline (ss, l, '/')) { }
	
	l = l.substr (0, l.size()-4);
	
	for (int i =0; i < l.size(); ++i) l[i] = toupper (l[i]);

	return l;
}


int main (int argc, char *argv[])
{
	unordered_set<string> usnr = read_matrix ("non_redundant_VAST.txt");

	string p, d;
	
	for (int i = 1; argv[i] != NULL; ++i){
		p = argv[i];
		d = idupper (p);
//cout <<p <<"\n" <<d <<endl;		
		if (usnr.find (d) != usnr.end())
			system (("cp " + p + " ./non-redundant_PDB/" + d + ".pdb").c_str());
	}
	
	cout <<"OK!" <<endl;

	return 0;
}