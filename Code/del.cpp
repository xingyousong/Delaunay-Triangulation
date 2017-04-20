#include <algorithm>
#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <chrono>
#include <map>
#include <tuple>
#include <set>

#define REAL double
#define MILLI_PER_NANO 0.000001

extern "C" REAL orient2dexact(REAL* pa, REAL* pb, REAL* pc);
extern "C" REAL incircleexact(REAL* pa, REAL* pb, REAL* pc, REAL* pd);
/*
REAL orient2dexact(REAL* pa, REAL* pb, REAL* pc)
{
	REAL acx, bcx, acy, bcy;

	acx = pa[0] - pc[0];
	bcx = pb[0] - pc[0];
	acy = pa[1] - pc[1];
	bcy = pb[1] - pc[1];
	return acx * bcy - acy * bcx;

}

REAL incircleexact(REAL* pa, REAL* pb, REAL* pc, REAL* pd)
{
	REAL adx, ady, bdx, bdy, cdx, cdy;
	REAL abdet, bcdet, cadet;
	REAL alift, blift, clift;

	adx = pa[0] - pd[0];
	ady = pa[1] - pd[1];
	bdx = pb[0] - pd[0];
	bdy = pb[1] - pd[1];
	cdx = pc[0] - pd[0];
	cdy = pc[1] - pd[1];

	abdet = adx * bdy - bdx * ady;
	bcdet = bdx * cdy - cdx * bdy;
	cadet = cdx * ady - adx * cdy;
	alift = adx * adx + ady * ady;
	blift = bdx * bdx + bdy * bdy;
	clift = cdx * cdx + cdy * cdy;

	return alift * bcdet + blift * cadet + clift * abdet;
}
code so that visual studio works with this 
*/

using namespace std;

struct point
{	int id; 
	double coor[2];
	
};
bool x_first(point p1, point p2)
{
	return p1.coor[0] < p2.coor[0] || (p1.coor[0] == p2.coor[0] && p1.coor[1] < p2.coor[1]);
}
bool y_first(point p1, point p2)
{
	return p1.coor[1] > p2.coor[1] || (p1.coor[1] == p2.coor[1] && p1.coor[0] > p2.coor[0]);
}





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
	return (Rot() -> Onext()) -> Rot();
}
edge* edge::Dnext()
{
	return (Sym() -> Onext()) -> Sym();
}
edge* edge::Dprev()
{
	return (invRot() -> Onext()) -> invRot();
}
edge* edge::Lnext()
{
	return (invRot() -> Onext()) -> Rot();
}
edge* edge::Lprev()
{
	return Onext() -> Sym();
}
edge* edge::Rnext()
{
	return (Rot() -> Onext()) -> invRot();
}
edge* edge::Rprev()
{
	return Sym() -> Onext();
}
void edge::EndPoints(point* origin, point* de)
{
	data = origin;
	(Sym() -> data) = de;
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
	(Sym() -> data) = p;
}

edge* MakeEdge()
{
	QuadEdge *q1 = new QuadEdge;
	return q1->e;
}
//-------------------------------------------------- end of edge class-------------------------------
//---------------------------------------------------topological functions-------------------------------
void Splice(edge* a, edge* b)
{
	edge* alpha = (a->Onext()) -> Rot();
	edge* beta = (b->Onext()) -> Rot();

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
	Splice(e->Sym(), (e-> Sym())->Oprev());
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

void Swap(edge* e) // I never used this function
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

bool Valid(edge* e, edge* base1)
{
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

//-------------------------------------





class subdivision //this is a basic container that will output the triangles
{
public:
	void addEdge(int, int);
	vector<point> s;
	point* getpoint (int);
	void addpoint(point);
	void addEdgepoint(point*, point*);
	void killdupedge();
	vector< pair< int, int  > > edgelist; //supposed to be private
	vector< tuple<int, int, int> > trianglelist; // supposed to be private
	void neighborsort(); 
	void maketriangle(); // main meat
	void printalledge();
private:
	map< int, point*> pointid; 
	map< int, vector<int> > adjlist;  

};
//------------------------------------
// algorithm to make a triangluation output will be sort each adjacency list for each vertex, 
// as we add more vertices, if there is a collision, check consecutive
void subdivision::addEdge(int a, int b)
{	if (a==b)
	{
		return;
	}
	if (a > b) {swap(b,a);}
	pair < int, int> temp (a, b ); //NEED UNIQUENESS 

	edgelist.push_back(temp);

	adjlist[a].push_back(b);
	adjlist[b].push_back(a);
}
point* subdivision::getpoint(int i)
{
	return pointid[i];
}

void subdivision::addpoint(point p)
{	
	
	s.push_back(p);
	pointid[p.id] = &p;
	vector<int> temp; 
	adjlist[p.id] = temp;
	
}
bool tuplecomp(tuple<int, int> e1, tuple<int, int> e2 )
{
	int p1 = get<0>(e1); 
	int p2 = get<1>(e1);
	int p3 = get<0>(e2); 
	int p4 = get<1>(e2);

	if( (p1 == p3 && p2 == p4) || (p1 == p4 && p2 == p3) )
	{
		return true; 
	}
	else
	{
		return false;
	}
}
/*
bool custom(tuple<int, int> e1, tuple<int, int> e2)
{
	return p1.coor[0] < p2.coor[0] || (p1.coor[0] == p2.coor[0] && p1.coor[1] < p2.coor[1]);
}
*/

void subdivision::killdupedge()
{	/*
	vector<pair<int, int>>::iterator it;
	it = unique (edgelist.begin(), edgelist.end());
	edgelist.resize(distance(edgelist.begin(), it));
	*/
	
	set< pair<int, int> > temp(edgelist.begin(), edgelist.end());
	edgelist.assign(temp.begin(), temp.end());
	
}
void subdivision::addEdgepoint(point* a, point* b)
{
	addEdge(a -> id, b -> id );
}
/*
bool anglecomp ( int p1, int p2, int center )
{

}

void subdivision::neighborsort() // for each node, sorts his neighbors in cc order 
{

}
*/
void subdivision::printalledge()
{
	for (int i = 0; i<edgelist.size(); i++)
	{	int temp1 = get<0>(edgelist[i]);
		int temp2 = get<1>(edgelist[i]);
		cout << temp1 << " to " << temp2 << endl;
	}
}



//---- end of subdivision container



//--------------------------------- entire delaunay program
edgepair delaunay(subdivision &sub, int begin, int end, bool vertical, bool alternate) //begin = 0, end = length 
{

if ( (end - begin ) == 2 )
{	
	if(vertical)
	{
		sort((sub.s).begin() + begin, (sub.s).begin() + end , x_first);
	} 
	else
	{
		sort((sub.s).begin() + begin, (sub.s).begin() + end , y_first);
	}



edge* a = MakeEdge();
a->changeorg(&(sub.s[begin]));
a->changedest(&(sub.s[begin+1]));
edgepair ep;
ep.le = a;
ep.re = a->Sym();
//------------------print
printedge(a);
cout << "a" << endl;
sub.addEdgepoint( a-> Org(), a->Dest() );
//------------------print
return ep;
}//-----------end of if size = 2
else if((end - begin ) == 3)
{	
	if(vertical)
	{
		sort(sub.s.begin() + begin, sub.s.begin() + end, x_first);
	} 
	else
	{
		sort(sub.s.begin() + begin, sub.s.begin() + end, y_first);
	}
	
	edge* a = MakeEdge();
	edge* b = MakeEdge();
	
	a->changeorg(&(sub.s[begin])); //a,Org <- s1
	b->changeorg(&(sub.s[begin+1]));
	a->changedest(b -> Org()); //a.Dest <--b.Org <--- s2;s
	b->changedest(&(sub.s[begin+2]));
	Splice(a->Sym(), b);
	//------------------------------print
	printedge(a);
	printedge(b);
	cout << "a, b" << endl;
	sub.addEdgepoint( a-> Org(), a->Dest() );
	sub.addEdgepoint( b-> Org(), b->Dest() );
	//------------------------------print

//--------[Now Close the Triangle]--
	if(orient2dexact(sub.s[begin].coor, sub.s[begin+1].coor, sub.s[begin+2].coor) >0)
	{
		edge* c = Connect(b,a);
		edgepair ep;
		ep.le = a;
		ep.re = b -> Sym();

		//---------print
		printedge(c);
		cout << "c" << endl;
		sub.addEdgepoint( c-> Org(), c->Dest() );
		//---------print

		return ep;
	}
	else if(orient2dexact(sub.s[begin].coor, sub.s[begin+2].coor, sub.s[begin+1].coor) >0)
	{
		edge* c = Connect(b,a);
		edgepair ep;
		ep.le = c->Sym();
		ep.re = c;
		//-----------------------------
		printedge(c);
		cout << "c" << endl;
		sub.addEdgepoint( c-> Org(), c->Dest() );
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
{	int middle = (begin + end) / 2;
	
	if (vertical) 
	{
        std::nth_element(
            sub.s.begin() + begin,
            sub.s.begin() + middle,
            sub.s.begin() + end,
            x_first
        );
    } else 
    {
        std::nth_element(
            sub.s.begin() + begin,
            sub.s.begin() + middle,
            sub.s.begin() + end,
            y_first
        );

    }

    edgepair epl;
    edgepair epr; 
	if (alternate) 
	{
	epl = delaunay(sub, begin, middle, !vertical, alternate);
	epr = delaunay(sub, middle, end, !vertical, alternate);
	}
	else
	{
	epl = delaunay(sub, begin, middle, vertical, alternate);
	epr = delaunay(sub, middle, end, vertical, alternate);
	}
	
	//[ldo, ldi] = [epl.le, epl.re]
	//[rdi, rdo] = [epr.le, epr.re]
	edge* ldo = epl.le;
	edge* ldi = epl.re;
	edge* rdi = epr.le;
	edge* rdo = epr.re;

// alternate switches 

	if (alternate) 
	{
        // If vertical is True, then we need leftmost + rightmost 
        // If vertical is False, then we need topmost + bottommost
       
        bool (*comparator)(point, point ) = vertical ? x_first : y_first; // selection which sorting

        // ldo is ccw edge from topmost/leftmost vertex
        
        if (vertical) { // rest of loops are to move to the correct positions 
            while (comparator(*(ldo->Rprev()->Org()), *(ldo->Org()))) {
                ldo = ldo->Rprev();
            }
            
            while (comparator(*(ldi->Org()), *(ldi->Lprev()->Org()))) {
                ldi = ldi->Lprev();
            }
            // and the same logic for rdi, rdo
            while (comparator(*(rdi->Rprev()->Org()), *(rdi->Org()))) {
                rdi = rdi->Rprev();
            }
            while (comparator(*(rdo->Org()), *(rdo->Lprev()->Org()))) {
                rdo = rdo->Lprev();
            }
        } else {
            
            while (comparator(*(ldo->Rnext()->Org()), *(ldo->Org()))) {
                ldo = ldo->Rnext();
            }
           
            while (comparator(*(ldi->Org()), *(ldi->Lnext()->Org()))) {
                ldi = ldi->Lnext();
            }
          
            while (comparator(*(rdi->Rnext()->Org()), *(rdi->Org()))) {
                rdi = rdi->Rnext();
            }
            while (comparator(*(rdo->Org()), *(rdo->Lnext()->Org()))) {
                rdo = rdo->Lnext();
            }

        }
    }

//alternate switches end
	while(true)
	{

		if(LeftOf(rdi ->Org(), ldi)) //LeftOf[rdi.Org, ldi]
		{
			ldi = ldi->Lnext();
		}
		else if (RightOf(ldi -> Org(), rdi))
		{
			rdi = rdi -> Rprev();
		}
		else
		{
			break;
		}
	}
	edge* base1 = Connect(rdi -> Sym(), ldi);

	//--------------------------------
	printedge(base1);
	cout << "base1-------" << endl; 
	sub.addEdgepoint( base1-> Org(), base1->Dest() );
	//--------------------------------
	if ((ldi-> Org() ) == (ldo -> Org())) {ldo = base1 ->Sym();}
	if( (rdi -> Org()) == (rdo->Org())) {rdo = base1;}

while(true)
{
//DO this is the merge loop
	edge* lcand = (base1 -> Sym()) -> Onext();
	if(Valid(lcand, base1))
	{
		while( incircleexact(base1->Dest() -> coor, base1 ->Org() -> coor, lcand ->Dest() -> coor, (lcand ->Onext()) -> Dest() -> coor  )  > 0 )
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
	if (!Valid(lcand, base1) && !Valid(rcand, base1))
	{
		break;
	}
	if ( !Valid(lcand, base1) || (Valid(rcand, base1) && incircleexact(lcand-> Dest()-> coor, lcand->Org() -> coor, rcand ->Org() -> coor, rcand ->Dest() -> coor  )  > 0))
	{
		base1 = Connect(rcand, base1 -> Sym());

		//----------------------------------
		printedge(base1);
		cout << "base1x" << endl;
		sub.addEdgepoint( base1-> Org(), base1->Dest() );
		//----------------------------------
	}
	else
	{
		base1 = Connect(base1 -> Sym(), lcand ->Sym());
		//-----------------------------------
		printedge(base1);
		cout << "base1y" << endl;
		sub.addEdgepoint( base1-> Org(), base1->Dest() );
		//-------------------------------------
	}
	//OD
}
	edgepair output;
	output.le = ldo;
	output.re = rdo;
	//------------------------------
	printep(output);
	cout << "output" << endl;
	sub.addEdgepoint( ldo-> Org(), ldo->Dest() );
	sub.addEdgepoint( rdo-> Org(), rdo->Dest() );
	return output;
} //---end of |S| >=4



} //---------end of delaunay

//--------------------------

int main()
{
string line;
string filename;




cout <<"Name of file? Ex: 4.node" << endl;
//cin >>  filename;
//char* S = filename.c_str();

subdivision sub; 


ifstream myfile;
myfile.open("box.node");



if(myfile.is_open())
{
	getline(myfile, line);
	int N = line[0] - '0';

	while(getline(myfile, line))
	{
		istringstream is(line);
		double id;
		double x;
		double y;
		is >> id;
		is >> x;
		is >> y;
		point p;
		p.id = id;
		p.coor[0] = x;
		p.coor[1] = y;
		sub.addpoint(p);
	}
myfile.close();

printpoints(sub.s);

bool vertical = true;
bool alternate = false; 
edgepair eppp = delaunay(sub, 0, (sub.s).size(), vertical, alternate);
sub.killdupedge(); //not working????
cout << sub.edgelist.size() << endl; 
sub.printalledge();


}
else
	{
		cout <<"Unable to open file" << endl;
	}



}
