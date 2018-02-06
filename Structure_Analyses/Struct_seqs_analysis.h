#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <set>
#include <regex>
#include <sstream>
#include <unordered_map>
#include <tuple>
#include <boost/math/distributions/chi_squared.hpp>

using namespace std;

pair<string,string> psplit (string l)
{
	stringstream ss (l);
	string m;
	
	ss >> l >> m;
	
	return make_pair(l,m);
}

vector<string> vsplit (string l)
{
	stringstream ss (l);
	vector<string> v;
	
	while (ss >> l)
		v.push_back(l);
	
	return v;
}

string get_id (string l)
{
	stringstream ss (l);
	
	ss >> l;
	l.erase (0,1);
	
	return l;
}

unordered_map<string,string> read_seq (string file)
{
	ifstream f (file);
	string s;
	string seq="", id="";
	unordered_map<string,string> umpp;
	
	while (getline (f,s))
		if (s[0]=='>'){
			if (seq!=""){
				umpp[id] = seq;
				seq="";
			}
			id = get_id (s);
		}
		else
			seq += s;
	
	umpp[id] = seq;
	
	return umpp;
}
		