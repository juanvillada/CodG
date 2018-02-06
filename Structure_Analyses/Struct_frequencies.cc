#include "seqs_analysis.h"

string get_PDB_ID (string file)
{
	ifstream f (file);
	string s;
	regex e ("(Y\\D{2}\\d{3}\\D{1})");
	smatch sm;
	
	while (getline (f,s))
		if ((s.find ("SOURCE") != string::npos)&&(s.find ("GENE:") !=string::npos))
			if (regex_search (s,sm,e)){
				return sm[1];
			}			
	f.close();
	
	return "";
}

unordered_map<string, string> get_codons (string file)
{
	ifstream f (file);
	string s, d;
	pair<string,string> pp;
	unordered_map<string,string> ums;
	regex e ("\\D-");
	
	getline (f, s); //get rid head of file
	
	while (getline (f, s)){
		pp = psplit (s);
		d = regex_replace (pp.first, e, "");
		ums[d] = pp.second;
	}
	f.close();
	
	return ums;
}

inline string get_first_tag (string l)
{
	stringstream ss (l);
	ss >> l;
	
	return l;
}

void get_PDB_data (string pdbfile, string cds, vector<tuple<string,int,int>> &vHelix, vector<tuple<string,int,int>> &vSheet)
{

	ifstream f (pdbfile);
	if (!f.is_open()) { cerr <<"Error in open the PDB file: " <<pdbfile <<endl; exit(1); }
	string s;
	int pi, pe;
	vector<string> vseq, vt;
	map<int,int> mpHelix;
		
  try{
	
	while (getline (f, s))
		if (get_first_tag(s) == "HELIX"){
			vt = vsplit (s);
			pi = (atoi (vt[5].c_str())-1)*3;
			pe = (atoi (vt[8].c_str())-1)*3;
			if ((pe+3 < cds.size()) && (pi > 0) && (pe > 0) && (pi < (pe+3)))
				vHelix.push_back (make_tuple(cds.substr (pi, pe-pi+3), pi, pe));
		}
		else if (get_first_tag(s) == "SHEET"){
			vt = vsplit (s);
			if ((vt.size() == 11)||(vt.size()>18)){
				pi = (atoi (vt[6].c_str())-1)*3;
				pe = (atoi (vt[9].c_str())-1)*3;
				if ((pe+3 < cds.size()) && (pi > 0) && (pe > 0) && (pi < (pe+3)))
					vSheet.push_back (make_tuple(cds.substr (pi, pe-pi+3), pi, pe));
			}
		}

}catch (const exception &e)
{
	cerr <<e.what() <<"\t" <<pi <<"\t" <<pe <<endl;
}
	f.close();
}	

//cout <<"chegou!" <<endl; ofstream fo ("test_seq.txt"); for (auto q: vseq) fo <<q <<endl; fo.close(); cin.get(); 
//cout <<sq.substr (pi, pe-pi+3) <<"\t" <<pi <<"\t"<<pe  <<endl; 
//cout <<sq.substr (pi-1, 3) <<"\n" <<sq.substr (pe, 3) <<"\n" <<pi <<"\t" <<pe <<endl;
//cin.get();

vector<string> codons (string l)
{
	vector<string> vcd;
	
	for (int i = 0; i < l.size(); i+=3)
		vcd.push_back (l.substr (i, 3));
	
	return vcd;
}

map<string, int> freq_codons (vector<tuple<string,int,int>> vhs, unordered_map<string,string> &umcd)
{
	vector<string> vcodons, vc;

	for (auto g: vhs){
		vc = codons (get<0> (g));
		vcodons.insert (vcodons.cend(), vc.begin(), vc.end());
	}
	
	map <string, int> mpc;
	string cd;
	
	for (auto r: vcodons){
		if (umcd.find (r) != umcd.end()){
			cd = umcd[r];
			mpc[cd] += 1;
		}
	}

	return mpc;
}


double chisq_test (map<string,vector<int>> mpi, ofstream &of)
{
	vector<float> vtot;
	float sumtotal = 0.0;
	float sum = 0.0;
	float n = float (mpi.size());

	for (auto a: mpi){
		for (auto b: a.second)
			sum += float(b);
		vtot.push_back (sum);
		sumtotal += sum;
		sum = 0.0;
	}
	of <<"Total" <<flush;
	for (auto t: vtot)
		of <<"\t" <<t <<flush;
	of<<endl;
	
	of <<"Degree of freedom: " <<n-1.0 <<endl;
	
	float e = sumtotal/n;
	of <<"Expected Proportion: " <<e/n <<endl;
	
	float chi = 0.0;

	of <<"H0: All codons have the same proportion: 1:1:1:1" <<endl;
	
	of <<"Observed Proportion:" <<flush;
	
	for (auto o: vtot){
		chi += ((o - e)*(o - e))/e;
		of <<"\t" <<o/n <<flush;
	}
	of <<endl;
	
	boost::math::chi_squared mydist(n-1.0);
	
	return boost::math::cdf(complement(mydist, chi));	
}

int main (int argc, char *argv[])
{
	if (argc==1) {cerr <<"Missing PDB files! Please, choose the directory with the parameter *.pdb" <<endl; return 0; }
	unordered_map<string, string> umCods = get_codons ("codons_Sc.txt");
	unordered_map<string,string> umcds = read_seq ("orf_coding.fasta");
	unordered_map<string,string>::iterator it;
	
cout <<"Codon information about CDS sequences stored!" <<endl;
	
	vector<string> vfiles;
	
	for (int i = 1; argv[i] != NULL; ++i)
		vfiles.push_back (argv[i]);
		
	string sid;
	map<string,string> mppdb;
	
	for (auto a: vfiles)
		if ((sid = get_PDB_ID(a)) != "")
			mppdb[sid] = a;

	cout <<"Number of PDB proteins: " <<mppdb.size() <<endl;
	
	map<string, vector<tuple<string,int,int>>> mpHelix;
	map<string, vector<tuple<string,int,int>>> mpSheet;
	vector<tuple<string,int,int>> vHelix;
	vector<tuple<string,int,int>> vSheet;
	
	
	for (auto p: mppdb){
cout <<"Analysing PDB file: " <<p.second <<"\tID: " <<p.first <<endl;
		if ((it=umcds.find (p.first)) != umcds.end()){
			get_PDB_data (p.second, it->second, vHelix, vSheet);
			if (!vHelix.empty())
				mpHelix[p.first] = vHelix;
			if (!vSheet.empty())
				mpSheet[p.first] = vSheet;
		}
		else
			cout <<p.first <<" not found in CDS sequences. Ignoring the PDB file: " << p.second <<endl;
		
		vHelix.clear();
		vSheet.clear();
	}

cout <<"PDB sequences processed! \t" <<mpHelix.size() <<" helix\t" <<mpSheet.size() <<" Sheets processed!" <<endl;

	map<string,vector<int>> mpvi;

	ofstream fo ("Helix_table_Scerevisae.tsv");
	fo <<"SGD_ID\tFreNO\tFreO\tRareNO\tRareO"<<endl;
	
	map<string,int> mpf;
	map<string,vector<int>> mpi;

	for (auto h: mpHelix){
		fo <<h.first <<flush;
		mpf =  freq_codons (h.second, umCods);
		for (auto m: mpf){
			fo <<"\t" <<m.second <<flush;
			mpi[m.first].push_back (m.second);
		}
		fo <<endl;
	}
	
	double pval = chisq_test (mpi, fo);	
	fo <<"Chi-square P-value = " <<pval <<endl;

	fo.close();
	
	fo.open ("Sheet_table_Scerevisae.tsv");
	fo <<"SGD_ID\tFreNO\tFreO\tRareNO\tRareO"<<endl;

	for (auto b: mpSheet){
		fo <<b.first <<flush;
		mpf =  freq_codons (b.second, umCods);
		for (auto m: mpf){
			fo <<"\t" <<m.second <<flush;
			mpi[m.first].push_back (m.second);
		}
		fo <<endl;
	}

	pval = chisq_test (mpi, fo);	
	fo <<"Chi-square P-value = " <<pval <<endl;;

	fo.close();
	cout <<"OK!" <<endl;
	
	return 0;
}







































