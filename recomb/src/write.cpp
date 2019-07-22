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

int findArc(int a, std::string t)
{
	int res=-1;
	int s=t.length();
	int i=0;
	int j=0;
	int test=-1;
	char c;
	std::string sub="";
	while(i<s && res==-1)
	{
		c=t[i];
		if(c==':')
		{
			j=i;
		}
		if(c=='{')
		{
			sub=t.substr(i+1);
			test=std::stoi(sub,nullptr);
			if(test==a)
			{
				res=j+1;
			}
		}

		i++;
	}
	return res;
}

double lengthArc(int pos, std::string t)
{
	double res=-1;
	if (pos>-1)
	{
		std::string sub=t.substr(pos);
		res=std::stod(sub,nullptr);
	}
	return res;
}

std::string subtree(int pos, std::string t)
{
	std::string sub="(";
	if(pos>-1)
	{
		int p=pos-2;
		if(t[p]==')')
		{
			int brackets=1;
			p+=-1;
			while(brackets !=0)
			{
				brackets+=(t[p]==')')-(t[p]=='(');
				p+=-1;
			}
			sub+=t.substr(p+2,pos-p-4);
		}
		else
		{
			while(t[p]!=',' && t[p]!='(')
			{
				p+=-1;
			}
			sub+=t.substr(p+1,pos-p-2);
		}
	}
	sub+=");";
	return sub;
}

double distNode(int pos, std::string t)
{
	double res=0;
	if(pos>-1)
	{
		res+=lengthArc(pos,t);
		int s=t.length();
		for(int p=pos; p<s; p++)
		{
			if(t[p]==')' && t[p+1]!=';')
			{
				res+=lengthArc(p+2,t);
				p++;
			}
		}
	}
	return res;
}

int nbLeaves(std::string t)
{
	double res=0;
	int s=t.length();
	for(int p=0; p<s; p++)
	{
		if((t[p]=='(' || t[p]== ',') && t[p+1]!='(')
		{
			res++;
		}
	}
	return res;
}

double sumPaths(std::string t)
{
	double res=0;
	int s=t.length();
	int q=0;
	for(int p=0; p<s; p++)
	{
		if((t[p]=='(' || t[p]== ',') && t[p+1]!='(')
		{
			q=p;
			while(t[q]!=':' && q<s)
			{
				q++;
			}
			if(q<s)
			{
				res+=distNode(q+1,t);
			}
			p=q;
		}
	}
	return res;
}

double lengthPendant(int pos, std::string t)
{
	double res=0;
	std::string sub=subtree(pos,t);
	res=sumPaths(sub)/nbLeaves(sub)+lengthArc(pos,t)/2;
	return res;
}

void resInit(std::vector<Arc*> *res)
{
	int s=(*res).size();
	for(int i=0; i<s; i++)
	{
		(*res).pop_back();
	}
}

void windInit(std::vector<std::vector<Arc*>> *wind)
{
	int s=(*wind).size();
	for(int i=0; i<s; i++)
	{
		(*wind).pop_back();
	}
}

/*void jplace(std::vector<std::string> q, std::string t, std::vector<std::vector<Kmer*>> list,std::vector<Arc>* branches, std::vector<std::string> infos, std::vector<Arc*> *res)
{
	ofstream doc("placement.jplace");
	if(doc)
	{
		doc << "\"tree\":\"" << t << "\"" << endl << "[";
		int s=q.size();
		Arc a(0,0);
		std::vector<Arc*> read(0);
		for(int i=0; i<s; i++)
		{
			resInit(res);
			doc << "{" << endl;
			doc << "\"p\":" << endl;
			readQuery(q[i], list, branches, infos, res);
			a=*(*res)[bestArc(*res)];
			doc << "[[" << a.getPlace() << ", "; //Name of winning Arc
			doc << a.getScore() << ", "; // Score of that Arc
			int pos=findArc(a.getPlace(), t);
			doc << lengthArc(pos, t)/2 << ", "; // Place of the ghost node along that Arc (half the length)
			doc << lengthPendant(pos,t) << "]]" << endl; // Length of the arc from the ghost node to the ghost leaf (mean of the lengths from the ghost nodes to all leaves below)
			doc << "\"nm\":" << endl;
			doc << "[[\"query" << i+1 << "\",1]]" << endl;
			doc << "}" << endl;
		}
		doc << "]" << endl;
	}
}*/

void SciPlot(int n, std::string_view query, std::vector<std::string> ref, std::vector<Arc> branches, std::vector<std::vector<Arc*>> res, int shift)
{
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
	int w=res.size();
	int a=0;
	if(w>0)
	{
		int m=res[0].size();
		std::string ad="/home/guillaume/Documents/LIRMM/Samples/res_" + query + ".csv";
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
