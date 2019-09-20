/*
 * write.cpp
 *
 *  Created on: 6 nov. 2018
 *      Author: scholz
 */


#include<iostream>
#include<fstream>
#include<vector>
#include <string>
#include<cmath>
#include "write.h"

using namespace std;
//builds output files.
//empties obsolete informations.

void resInit(std::vector<Arc*> *res)
{
	//empties arc* vector
	int s=(*res).size();
	for(int i=0; i<s; i++)
	{
		(*res).pop_back();
	}
}

void windInit(std::vector<std::vector<Arc*>> *wind)
{
	//empties windows vector
	int s=(*wind).size();
	for(int i=0; i<s; i++)
	{
		(*wind).pop_back();
	}
}

void SciPlot(int n, std::string_view query, std::vector<std::string> ref, std::vector<Arc> branches, std::vector<std::vector<Arc*>> res, int shift)
{
	//SciLab executable file for printing curves.
	//Not optimized; temporary (until plot possible from csv)
	int w=res.size();
	if(w>0)
	{
		int m=res[0].size();
		int b=branches.size();
		int a=0;
		std::string ad="/home/guillaume/Documents/LIRMM/Samples/plot_" + to_string(n) + ".sce";
		ofstream plot(ad);
		if(plot)
		{
			//plot << "// query: " << query << endl;
			plot << "clear;" << endl;
			plot << "W=0:" << w-1 << endl;
			plot << "M=%nan*ones(" << b-1 <<"," << w << ")" << endl;
			plot << "check=zeros(1," << b-1 << ")" << endl;
			for(int i=0; i<w; i++)
			{
				for(int j=0; j<m; j++)
				{
					a=(*res[i][j]).getPlace();
					plot << "M(" << a << "," << i+shift <<")=" << (*res[i][j]).getScore() << endl;
					plot << "check(" << a << ")=1" << endl;
				}
			}
			plot << "for i=1:" << b-1 << endl;
			plot << "if check(i)>0" << endl;
			plot << "plot2d(W,M(i,1:" << w << "),style=modulo(i," << ref.size() << ")+1)" << endl;
			plot << "end" << endl << "end" << endl;
			plot << "xtitle(\"" << query << "\",\"Window\",\"Score\")" << endl;
			/*plot << "legend(";
			for(int i=1; i<ref.size(); i++)
			{
				plot << "\"" << ref[i];
				if(i<ref.size()-1)
				{
					plot << "\", ";
				}
			}
			plot << "\")" << endl;*/
		}
	}
}

void Csv(int n, std::string query, std::vector<std::string> ref, std::vector<std::vector<Arc*>> res, int shift)
{
	//Writes Csv file; Window, Group, Score.
	int w=res.size();
	int a=0;
	if(w>0)
	{
		int m=res[0].size();
		std::string ad="/home/guillaume/Documents/LIRMM/Samples/HIV/holidays-data/res_" + query + ".csv";
		ofstream plot(ad);
		if(plot)
		{
			//plot << "// query: " << query << endl;
			for(int i=0; i<w; i++)
			{
				for(int j=0; j<m; j++)
				{
					a=(*res[i][j]).getPlace();
					if(a<ref.size())
					{
						plot << i+shift << "," << ref[a] << "," << (*res[i][j]).getScore() << endl;
					}
					else
					{
						plot << i+shift << "," << a << "," << (*res[i][j]).getScore() << endl;
					}
				}
			}

		}
	}
}

int stars(std::string gr)
{
	int res=0;
	for(int i=0; i<gr.length(); i++)
	{
		if(gr[i]=='*')
		{
			res++;
		}
	}
	return res;
}

std::string color(std::string gr)
{
	std::string res="";
	if(gr.length()<3)
	{
		res=gr;
	}
	else
	{
		int s=stars(gr);
		if(s>3)
		{
			res="black";
		}
		else
		{
			std::string ex=to_string(100/(s+1));
			for(int i=0; i<gr.length(); i++)
			{
				if(gr[i]=='*')
				{
					res=res+"!" + ex + "!";
				}
				else
				{
					res.push_back(gr[i]);
				}
			}
		}
	}
    return res;
}

std::string div(std::string bp)
{
	std::string res="0.000";
	int l=bp.length();
	if(l==4)
	{
		res[0]=bp[0];
	}
	for(int i=0; i<4; i++)
	{
		if (i<=l)
		{
		res[5-i]=bp[l-i];
		}
	}
	return res;
}

void tikzLine(std::string rec, int y)
//takes line of infos or output; returns command line to plot the info in tikz (needs post-treatment)
{
	std::vector<std::string> infos(0);
	std::string tmp="";
	for(int i=0; i<rec.length(); i++)
	{
		if(rec[i] !=',')
		{
			tmp+=rec[i];
		}
		else
		{
			infos.push_back(tmp);
			tmp="";
		}
	}
	infos.push_back(tmp);
	int lines=(infos.size()-1)/2;
	for(int j=0; j<lines; j++)
	{
		cout << "\\fill[color=" << color(infos[2*j+1]) << "](" << div(infos[2*j]) << "," << y << ") rectangle(" << div(infos[2*(j+1)]) << "," << y+1 << ");" << endl;
	}
	cout << endl;
}

void tikzDoc(std::string ref_R, std::string ref_A, std::string ref_B)
{
	std::string line="";
	std::vector<std::string> nm (0);
	std::vector<std::string> real (0);
	std::vector<std::string> resA (0);
	std::vector<std::string> resB (0);
	ifstream dataR(ref_R);
	if(dataR)
	{
		std::string prefix="0,";
		while(getline(dataR,line)) //-> does not work on pc ??
		{
			if(line[0] !='>')
			{
				line=prefix+line;
				//cout << line << endl;
				real.push_back(line);
				cout << endl;
			}
			else
			{
				if(line[6]=='_')
				{
					line[6]=' ';
				}
				nm.push_back(line);
			}
		}
		dataR.close();
	}
	ifstream dataA(ref_A);
	if(dataA)
	{
		while(getline(dataA,line)) //-> does not work on pc ??
		{
			if(line[0] !='>')
			{
				//cout << line << endl;
				resA.push_back(line);
				cout << endl;
			}
		}
		dataA.close();
	}
	ifstream dataB(ref_B);
	if(dataB)
	{
		while(getline(dataB,line)) //-> does not work on pc ??
		{
			if(line[0] !='>')
			{
				//cout << line << endl;
				resB.push_back(line);
				cout << endl;
			}
		}
		dataB.close();
	}
	if(real.size()>0 && resA.size()>0 && resB.size() > 0)
	{
		for(int i=0; i<resA.size(); i++)
		{
			cout << "\\begin{tikzpicture}" << endl;
			cout << "\\draw (0,2.5) node{};" << endl;
			cout << "\\draw(0,1.5) node{" << nm[i] <<"};" << endl;
			cout << "\\draw(-0.5,0.5) node{R};" << endl;
			cout << "\\draw(-0.5,-1.5) node{1.};" << endl;
			cout << "\\draw(-0.5,-3.5) node{2.};" << endl;
			tikzLine(real[i],0);
			tikzLine(resA[i],-2);
			tikzLine(resB[i],-4);
			cout << "\\end{tikzpicture}" << endl << endl;
		}
	}
	else
	{
		cout << "Something went wrong when reading files" << endl;
		cout << real.size() << "; " << resA.size() << "; " << resB.size() << endl; 
	}
}
