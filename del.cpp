#include <algorithm>
#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>

#define REAL double
 
extern "C" double doublerand(); 
extern "C" REAL orient2dexact(REAL* pa, REAL* pb, REAL* pc);
extern "C" REAL incircleexact(REAL* pa, REAL* pb, REAL* pc, REAL* pd); 
using namespace std;

struct point
{
	double coor[2];
	bool operator < (const point& p2) const
	{
		if(coor[0] < p2.coor[0])
			{
				return true;
			}
		else if (coor[0] > p2.coor[0])
		{
			return false;
		}
		else
		{
			return (coor[1] < p2.coor[1]);
		}
	} 
};

class QuadEdge;

class edge
{
	friend QuadEdge; 
	friend void Splice(edge*, edge*);
private:
	int num; 
	edge *next;
	point *data;
public:
	edge() {data = 0;}
	edge* Rot();
	edge* invRot();
	edge* Sym();
	
	edge* Onext();
	edge* Oprev();
	edge* Dnext();
	edge* Dprev();
	edge* Lnext();
	edge* Lprev();
	edge* Rnext();
	edge* Rprev();

	point* Org();
	point* Dest();
	void EndPoints(point*, point*);
	QuadEdge* Qedge() {return (QuadEdge*)(this-num);}
	void changeorg(point*);
	void changedest(point*);

};

class QuadEdge 
{
	friend edge *MakeEdge();
private:
	edge e[4];
public:
	QuadEdge();


};

QuadEdge::QuadEdge()
{
	e[0].num = 0, e[1].num = 1, e[2].num = 2, e[3].num = 3;
	e[0].next = &(e[0]); e[1].next = &(e[3]);
	e[2].next = &(e[2]); e[3].next = &(e[1]);

}
edge* edge::Rot()
{
	return (num<3)? this +1 : this -3;
}
edge* edge::invRot()
{
	return (num >0)? this -1: this +3;
}
edge* edge::Sym()
{
	return (num <2) ? this +2 : this - 2;
}
edge* edge::Onext()
{
	return next; 
}
edge* edge::Oprev()
{
	return Rot() -> Onext() -> Rot();
}

edge* edge::Dnext()
{
	return Sym() -> Onext() -> Sym();
}
edge* edge::Dprev()
{
	return invRot() -> Onext() -> invRot();
}
edge* edge::Lnext()
{
	return invRot() -> Onext() -> Rot();
}
edge* edge::Lprev()
{
	return Onext() -> Sym();
}
edge* edge::Rnext()
{
	return Rot() -> Onext() -> invRot();
}
edge* edge::Rprev()
{
	return Sym() -> Onext();
}

void edge::EndPoints(point* origin, point* de)
{
	data = origin;
	Sym() -> data = de;
}

point* edge::Org()
{
	return data;
}
point* edge::Dest()
{
	return Sym()-> data;
}

void edge::changeorg(point* p)
{
	data = p;
}

void edge::changedest(point* p)
{
	Sym() -> data = p;
}

edge* MakeEdge()
{
	QuadEdge *q1 = new QuadEdge; 
	return q1->e; 
}

void Splice(edge* a, edge* b)
{
	edge* alpha = a->Onext() -> Rot();
	edge* beta = b->Onext() -> Rot();

	edge* t1 = b->Onext();
	edge* t2 = a-> Onext();
	edge* t3 = beta -> Onext();
	edge* t4 = alpha -> Onext();

	a->next = t1; 
	b->next = t2; 
	alpha->next = t3;
	beta-> next = t4;
}
void DeleteEdge(edge* e)
{
	Splice(e, e->Oprev());
	Splice(e->Sym(), e-> Sym()->Oprev());
	delete e-> Qedge();
}



edge* Connect(edge* a, edge* b)
{
	edge* e = MakeEdge();
	e -> changeorg(a -> Dest());
	e -> changedest(b -> Org());
	Splice(e, a-> Lnext());
	Splice(e->Sym(), b);
	//e-> EndPoints(a->Dest(), b-> Org());
	return e;
}

void Swap(edge* e)
{
	edge* a = e->Oprev();
	edge* b = e-> Sym() -> Oprev();
	Splice(e,a);
	Splice(e->Sym(), b);
	Splice(e, a->Lnext());
	Splice(e->Sym(), b->Lnext());
	//e->EndPoints(a->Dest(), b->Dest());
	e -> changeorg(a -> Dest());
	e -> changedest(b ->Dest());
}

bool RightOf(point* p, edge* e)
{
	return (orient2dexact(p -> coor, e-> Dest()->coor, e->Org()->coor) > 0);
}
bool LeftOf(point* p, edge* e)
{
	return (orient2dexact(p->coor, e->Org()->coor, e->Dest()-> coor) > 0);
}

bool Valid(edge* e, edge* base1){
	return RightOf(e->Dest(), base1);
} 


struct edgepair
{
edge* le;
edge* re;
};

void printep(edgepair ep)
{
	double x1 = (ep.le -> Dest() -> coor)[0];
	double y1 = (ep.le -> Dest() -> coor)[1];
	double x2 = (ep.le -> Org() -> coor)[0];
 	double y2 = (ep.le -> Org() -> coor)[1];

 	double x3 = (ep.re -> Dest() -> coor)[0];
 	double y3 = (ep.re -> Dest() -> coor)[1];
 	double x4 = (ep.re -> Org() -> coor)[0];
 	double y4 = (ep.re -> Org() -> coor)[1];

 	cout << x1 << " " << y1 << " to " << x2 << " " << y2 << endl;
 	cout << x3 << " " << y3 << " to " << x4 << " " << y4 << endl;
 	
}

void printedge(edge* e)
{
	double x1 = (e -> Dest() -> coor)[0];
	double y1 = (e -> Dest() -> coor)[1];
	double x2 = (e -> Org() -> coor)[0];
 	double y2 = (e -> Org() -> coor)[1];
 	cout << x1 << " " << y1 << " to " << x2 << " " << y2 << endl;
}

void printpoints(vector<point> &s)
{
	cout << "Points In Order" << endl;
	for(int i=0; i< s.size(); i++)
	{	
		double x = ((s[i].coor)[0]);
		double y = ((s[i].coor)[1]);
		cout << x << " " << y << endl;
	}


}


//--------------------------------- entire delaunay program
edgepair delaunay(vector<point> &s, int begin, int end) //begin = 0, end = length -1
{

if ( (end - begin +1) == 2 )
{
edge* a = MakeEdge();
a->changeorg(&s[begin]);
a->changedest(&s[begin+1]); 
edgepair ep;
ep.le = a;
ep.re = a->Sym();
//------------------print
printedge(a);
//------------------print
return ep;
}//-----------end of if size = 2
else if((end - begin +1) == 3)
{
	edge* a = MakeEdge(); 
	edge* b = MakeEdge();
	Splice(a->Sym(), b);
	a->changeorg(&s[begin]); //a,Org <- s1
	a->changedest(&s[begin+1]); //a.Dest <--b.Org <--- s2;s
	b->changeorg(&s[begin+1]);
	b->changedest(&s[begin+2]);

	//------------------------------print
	printedge(a);
	printedge(b);
	//------------------------------print

//--------[Now Close the Triangle]--
	if(orient2dexact(s[begin].coor, s[begin+1].coor, s[begin+2].coor) >0)
	{
		edge* c = Connect(b,a);
		edgepair ep; 
		ep.le = a; 
		ep.re = b -> Sym();

		//---------print
		printedge(c);
		//---------print
		
		return ep;
	}
	else if(orient2dexact(s[begin].coor, s[begin+2].coor, s[begin+1].coor) >0)
	{
		edge* c = Connect(b,a);
		edgepair ep;
		ep.le = c->Sym();
		ep.re = c;
		//-----------------------------
		printedge(c);

		//-----------------------------
		return ep;
	}
	else
	{	edgepair ep;
		ep.le = a;
		ep.re = b->Sym();
		//-----------------
		
		return ep;

	} //---- three points are collinear
} //---end of |S| = 3
else //|S| >= 4--------------------
{	//[ldo, ldi] = [epl.le, epl.re]
	//[rdi, rdo] = [epr.le, epr.re]
	edgepair epl = delaunay(s, 0, end/2);
	edgepair epr = delaunay(s, end/2+1, end);
	if(LeftOf(epr.le ->Org(), epl.re)) //LeftOf[rdi.Org, ldi]
	{
		epl.re = epl.re->Lnext();
	}
	else if (RightOf(epl.re -> Org(), epr.le))
	{
		epl.le = epl.le -> Rprev();
	}
	else
	{
		// exit IF---
	}
	edge* base1 = Connect(epr.le -> Sym(), epl.re);

	//--------------------------------
	printedge(base1);
	//--------------------------------
	if(epl.re -> Org() == epl.le -> Org())
	{
		epl.le = base1 ->Sym();
	}
	if(epr.le -> Org() == epr.re->Org())
	{
		epr.re = base1;
	}


//DO this is the merge loop
	edge* lcand = base1 -> Sym() -> Onext();
	if(Valid(lcand, base1))
	{
		while( incircleexact(base1->Dest() -> coor, base1 ->Org() -> coor, lcand ->Dest() -> coor, lcand ->Onext() -> Dest() -> coor  )  > 0 )
		{
			edge* t = lcand -> Onext();
			DeleteEdge(lcand);
			lcand = t;
		}
	} 
	edge* rcand = base1 -> Oprev();
	if(Valid(rcand, base1))
	{
		while(incircleexact(base1->Dest() -> coor, base1 ->Org() -> coor, rcand ->Dest() -> coor, rcand ->Oprev() -> Dest() -> coor  )  > 0)
		{
			edge* t = rcand -> Oprev();
			DeleteEdge(rcand);
			rcand = t;
		}
	}
	if (!Valid(lcand, base1) and !Valid(rcand, base1))
	{
		//exit FI
	}
	if ( !Valid(lcand, base1) or (!Valid(rcand, base1) and incircleexact(lcand-> Dest()-> coor, lcand->Org() -> coor, rcand ->Org() -> coor, rcand ->Dest() -> coor  )  > 0))
	{
		base1 = Connect(rcand, base1 -> Sym());

		//----------------------------------
		printedge(base1);
		//----------------------------------
	}
	else
	{
		base1 = Connect(base1 -> Sym(), lcand ->Sym());
		//-----------------------------------
		printedge(base1);
		//-------------------------------------
	}
	//OD 
	edgepair output;
	output.le = epl.le;
	output.re = epr.re;
	//------------------------------
	printep(output);
	return output;
} //---end of |S| >=4



} //---------end of delaunay



int main()
{
string line;
vector<point> s; 
ifstream myfile("4.node");
if(myfile.is_open())
{
	getline(myfile, line);
	int N = line[0] - '0';
	
	while(getline(myfile, line))
	{
		istringstream is(line);
		double temp;
		double x;
		double y;
		is >> temp;
		is >> x;
		is >> y;
		point p; 
		p.coor[0] = x;
		p.coor[1] = y;
		s.push_back(p);
	}
myfile.close();
}
else
	{
		cout <<"Unable to open file";
	}

sort(s.begin(), s.end());
printpoints(s);


edgepair eppp = delaunay(s, 0, s.size()-1);


}