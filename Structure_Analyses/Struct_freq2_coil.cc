#include "../seqs_analysis.h"

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


void chisq_freqs (ofstream &of, map<int,tuple<int,int,int,int>> mpe)
{
	of <<"\n\n\n* * * * * Comparison of frequencies among the first, second and thrid codons * * * * * * *\n\n" <<endl;
	double chi;
	double pval;
	
	map<int,tuple<int,int,int,int>>::iterator it, jt;
	boost::math::chi_squared mydist(3.0);
	

	for (it = mpe.begin(); it != mpe.end(); ++it)
		for (jt = it; jt != mpe.end(); ++jt)
		if (it!=jt){
	
		of <<"______________________________________________\n Chi Squared test\n" <<"Codons " <<it->first <<" vs " <<jt->first <<"\n______________________________________________\n\n";
   	
   		chi = pow((get<0>(it->second)-get<0>(jt->second)),2)/get<0>(jt->second) + pow((get<1>(it->second)-get<1>(jt->second)),2)/get<1>(jt->second) + 
   		  	  pow((get<2>(it->second)-get<2>(jt->second)),2)/get<2>(jt->second) + pow((get<3>(it->second)-get<3>(jt->second)),2)/get<3>(jt->second);

		pval = boost::math::cdf(complement(mydist, chi));
		of << "Degree of Freedom = 3.0" <<endl;
		of << "Chi-square = " <<chi <<endl;
		of <<"p-value = " <<pval <<" => " <<((pval < 0.05)? "Reject H0 - proportions are different!\n" : "Accepted H0 - the proportions are the same!\n") <<"\n" <<endl;
	}
}



void freq_codons (string outfile, vector<vector<string>> vvs)
{
	ofstream fo (outfile);
	ofstream fo2 ("Raw_table_"+outfile);
	
	int FreNO, FreO, RareNO, RareO;
	int total;
	
	double E;
	double chi;
	double pval;
	
	map<int,tuple<int,int,int,int>> mpentre;

	for (int i = 1; i < 4; ++i){
		fo2 <<"Table related to Codon " <<i <<endl;
		fo <<"Table related to Codon " <<i <<endl;
		fo2 <<"SGD_ID]\tType of Codons: \tFreNO\tFreO\tRareNO\tRareO"<<endl;
		
		FreNO = FreO = RareNO = RareO =0;
		
		for (auto x: vvs){
			fo2 <<x[0] <<"\t" <<x[i] <<endl;
		
			if (x[i]=="FreNO") ++FreNO;
			else if (x[i]=="FreO") ++FreO; 
			else if (x[i]=="RareNO") ++RareNO;
			else if (x[i]=="RareO") ++RareO;
		}		
		fo2 <<"Total\t" <<FreNO <<"\t" <<FreO <<"\t" <<RareNO <<"\t" <<RareO <<endl <<endl;
		fo <<"Total\t" <<FreNO <<"\t" <<FreO <<"\t" <<RareNO <<"\t" <<RareO <<endl;
		
		mpentre[i] = make_tuple (FreNO, FreO, RareNO, RareO);
		
		total = FreNO+FreO+RareNO+RareO;
		E = double (total)/4.0;

		fo <<"______________________________________________\n"
   			 "Chi Squared test for sample standard deviation\n"
   			 "______________________________________________\n\n";

		chi = pow((FreNO-E),2)/E + pow((FreO-E),2)/E + pow(RareNO-E, 2)/E + pow(RareO-E,2)/E;

		fo <<"Degree of freedom: 3.0" <<endl;
		fo <<"Expected value: " <<E <<endl;
		fo <<"H0: All types of codons have the same proportion: 1:1:1:1" <<endl;
		fo <<"Significance of p-value below 0.05" <<endl;

		boost::math::chi_squared mydist(3.0);
		pval = boost::math::cdf(complement(mydist, chi));
		
		fo <<"p-value = " <<pval <<" => " <<((pval < 0.05)? "Reject H0 - proportions are different!\n" : "Accepted H0 - the proportions are the same!\n") <<endl;
	}

cout <<"Analysing frequencies among the types of codons.\n" <<endl;
	
	chisq_freqs (fo, mpentre);

	fo.close();
	fo2.close();

}


vector<vector<string>> find_cod (string id, vector<tuple<string,int,int>> vtup, unordered_map<string, string> &umCod)
{
	vector<vector<string>> vvm;
	
	for (auto t: vtup)
		if (get<0>(t).size() > 9){
			string c1 = get<0>(t).substr(0,3);
			string c2 = get<0>(t).substr(3,3);
			string c3 = get<0>(t).substr(6,3);
	
			if (!((umCod.find(c1) == umCod.end())||(umCod.find(c2) == umCod.end())||(umCod.find(c3) == umCod.end())))
				vvm.push_back (vector<string>({id, umCod[c1],umCod[c2],umCod[c3]}));
		}

	
	return vvm;
}

inline void get_regions (const vector<tuple<string,int,int>> &vt, map<int,int> &mpp)
{
	for (auto t: vt)
		mpp[get<1>(t)] = get<2>(t);	
}

map<string, vector<tuple<string,int,int>>> get_Coil_region (const map<string, vector<tuple<string,int,int>>> &mph, const map<string, vector<tuple<string,int,int>>> &mpsh, const unordered_map<string,string> &ucds)
{
	map<int,int> mppos;
	map<string, vector<tuple<string,int,int>>> mpc;
	vector<tuple<string,int,int>> vtup;
	
	map<string, vector<tuple<string,int,int>>>::const_iterator it;
	map<int,int>::iterator xt, yt;
	int p1, p2;
	string cds;

	unordered_map<string,string>::const_iterator ut;

	for (auto h: mph){
		get_regions (h.second, mppos);
		if ((it=mpsh.find (h.first)) != mpsh.end())
			get_regions (it->second, mppos);

		if ((ut=ucds.find(h.first))==ucds.end()) { cout <<"Critical error, cds is not found!" <<endl; exit(1); }
		cds = ut->second;

//cout <<h.first <<endl;

		yt = mppos.begin();
		for (xt = mppos.begin(); yt != mppos.end(); ++xt){
			yt = xt; ++yt;
		
			p1 = xt->second + 1;
			p2 = yt->first -1;
			
			//cout <<xt->first <<"\t" <<xt->second <<"\n" <<yt->first <<"\t" <<yt->second << "\t" <<p2-p1 <<endl; cin.get();
		
			if ((p1+3 < cds.size()) && (p1 > 0) && (p2 > 0) && (p1 < (p2+3)))
				vtup.push_back (make_tuple(cds.substr (p1, p2-p1+3), p1, p2));
	
		}
		
		mpc[h.first] = vtup;
		vtup.clear();
		mppos.clear();
		
	}
	

	return mpc;
}


int main (int argc, char *argv[])
{
	if (argc==1) {cerr <<"Missing PDB files! Please, choose the directory with the parameter *.pdb" <<endl; return 0; }
	unordered_map<string, string> umCods = get_codons ("../codons_Sc.txt");
	unordered_map<string,string> umcds = read_seq ("../orf_coding.fasta");
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


cout <<"Obtaining sequences between the helix and sheet structures" <<endl;
	map<string, vector<tuple<string,int,int>>> mpCoil = get_Coil_region (mpHelix, mpSheet, umcds);


cout <<"Analizing Coil Structure . . ." <<flush;	
	vector<vector <string>> vmc;
	vector<vector<string>> vmat;

	for (auto h: mpCoil){		
		
		vmc = find_cod (h.first, h.second, umCods);
		
		if (!vmc.empty()){
			vmat.insert (vmat.cend(), vmc.begin(), vmc.end());
		}
	}
cout <<"Size of structures analysed: " <<vmat.size() <<endl;
		
	freq_codons ("Coil_by_first_second_third_codons_chisq.tsv", vmat);

cout <<"OK" <<endl;	

	return 0;
}







































