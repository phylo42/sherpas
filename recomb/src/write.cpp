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

void tikzLine(std::string rec, double y, std::vector<std::string> *gr)
//takes line of infos or output; returns command line to plot the info in tikz (needs post-treatment)
{
	std::vector<std::string> infos(0);
	std::string tmp="";
	std:string col="";
	int check=0;
	int p=0;
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
		check=0;
		p=0;
		col=color(infos[2*j+1]);
		cout << "\\fill[color=" << col << "](" << div(infos[2*j]) << "," << y << ") rectangle(" << div(infos[2*(j+1)]) << "," << y+0.5 << ");" << endl;
		while(check==0 && p<(*gr).size())
		{
			check=(col==(*gr)[p]);
			p++;
		}
		if(check==0)
		{
			(*gr).push_back(col);
		}
		
	}
	cout << endl;
}

void tikzLegend(std::vector<std::string> gr)
//takes list of groups, plots the legend (name & colored squares)
{
	double x=0;
	for(int i=0; i<gr.size(); i++)
	{
		cout << "\\fill[color=" << gr[i] << "](" << x << "," << -3 << ") rectangle(" << x+1 << "," << -2.5 << ");" << endl;
		cout << "\\draw (" << x-0.5 << ",-3) node[above]{" << gr[i] << "};" << endl;
		x+=2;
		
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
	std::vector<std::string> gr (0);
	ifstream dataR(ref_R);
	if(dataR)
	{
		std::string prefix="0,";
		while(getline(dataR,line)) //-> does not work on pc ??
		{
			if(line[0] !='>')
			{
				//line=prefix+line;
				//cout << line << endl;
				real.push_back(line);
			}
			else
			{
				/*if(line[6]=='_')
				{
					line[6]=' ';
				}*/
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
			}
		}
		dataB.close();
	}
	if(real.size()>0 && resA.size()>0 && resB.size() > 0)
	{
		for(int i=0; i<resA.size(); i++)
		{
			//if(i+1==75 || i+1==83 || i+1==152 || i+1==191 || i+1==319 || i+1==383 || i+1==465 || i+1==480 || i+1==487 || i+1==562 || i+1==706 || i+1==782 || i+1==940 || i+1==963 || i+1==994 || i+1==1111 || i+1==1150 || i+1==1243 || i+1==1326 || i+1==1333 || i+1==1397 || i+1==1403 || i+1==1474 || i+1==1477 || i+1==1493 || i+1==1564 || i+1==1566 || i+1==1587 || i+1==1650 || i+1==1776 || i+1==1877 || i+1==1893 || i+1==1926 || i+1==1931 || i+1==1979 || i+1==2098 || i+1==2105 || i+1==2112 || i+1==2133 || i+1==2444 || i+1==2541 || i+1==2546 || i+1==2677 || i+1==2705 || i+1==2739 || i+1==2777 || i+1==2787 || i+1==2984 || i+1==3031 || i+1==3082 || i+1==3189 || i+1==3220 || i+1==3326 || i+1==3358 || i+1==3375 || i+1==3459 || i+1==3537 || i+1==3644 || i+1==3686 || i+1==3940 || i+1==4074 || i+1==4141 || i+1==4215 || i+1==4249 || i+1==4309 || i+1==4372 || i+1==4475 || i+1==4528 || i+1==4731 || i+1==4758 || i+1==4764 || i+1==4941 || i+1==4968 || i+1==4991 || i+1==5064 || i+1==5097 || i+1==5128 || i+1==5184 || i+1==5341 || i+1==5348 || i+1==5367 || i+1==5535 || i+1==5538 || i+1==5778 || i+1==5784 || i+1==5874 || i+1==6017 || i+1==6030 || i+1==6250 || i+1==6314 || i+1==6325 || i+1==6413 || i+1==6509 || i+1==6563 || i+1==6684 || i+1==6977 || i+1==7043 || i+1==7094 || i+1==7107 || i+1==7153 || i+1==7221 || i+1==7273 || i+1==7321 || i+1==7339 || i+1==7360 || i+1==7460 || i+1==7528 || i+1==7571 || i+1==7592 || i+1==7599 || i+1==7883 || i+1==7931 || i+1==8136 || i+1==8221 || i+1==8251 || i+1==8298 || i+1==8301 || i+1==8513 || i+1==8608 || i+1==8620 || i+1==8715 || i+1==8793 || i+1==8890 || i+1==8932 || i+1==8980 || i+1==9015 || i+1==9068 || i+1==9168 || i+1==9169 || i+1==9202 || i+1==9222 || i+1==9348 || i+1==9396 || i+1==9566 || i+1==9630 || i+1==9675 || i+1==9698 || i+1==9755 || i+1==9861 || i+1==9905 || i+1==9930 || i+1==9942 || i+1==9955)
			if(i>1990 && i<2021)
			{
				cout << "\\begin{tikzpicture}" << endl;
				cout << "\\draw (0,2) node{};" << endl;
				cout << "\\draw(0,1) node{>query " << i+1 <<"};" << endl;
				cout << "\\draw(-0.5,0) node[above]{R};" << endl;
				cout << "\\draw(-0.5,-1) node[above]{1.};" << endl;
				cout << "\\draw(-0.5,-2) node[above]{2.};" << endl;
				tikzLine(real[i],0, &gr);
				tikzLine(resA[i],-1, &gr);
				tikzLine(resB[i],-2, &gr);
				tikzLegend(gr);
				gr.clear();
				cout << "\\end{tikzpicture}" << endl;
				cout << "\\quad" << endl;
			}
		}
	}
	else
	{
		cout << "Something went wrong when reading files" << endl;
		cout << real.size() << "; " << resA.size() << "; " << resB.size() << endl; 
	}
}
