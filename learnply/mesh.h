#ifndef __MESH_H__
#define __MESH_H__

#include "icVector.H"
#include <stdio.h>
#include <stddef.h>

#include <vector>


using std::vector;

class Face;
class Mesh;
class Edge;
class Mat4x4;

class Mat4x4 {
public:
	double val[4][4];// = { 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 };

	Mat4x4() {
		val[0][0] = 0.;
		val[0][1] = 0.;
		val[0][2] = 0.;
		val[0][3] = 0.;

		val[1][0] = 0.;
		val[1][1] = 0.;
		val[1][2] = 0.;
		val[1][3] = 0.;

		val[2][0] = 0.;
		val[2][1] = 0.;
		val[2][2] = 0.;
		val[2][3] = 0.;

		val[3][0] = 0.;
		val[3][1] = 0.;
		val[3][2] = 0.;
		val[3][3] = 0.;
	}

	Mat4x4(double aa) {
		val[0][0] = aa;
		val[0][1] = aa;
		val[0][2] = aa;
		val[0][3] = aa;

		val[1][0] = aa;
		val[1][1] = aa;
		val[1][2] = aa;
		val[1][3] = aa;

		val[2][0] = aa;
		val[2][1] = aa;
		val[2][2] = aa;
		val[2][3] = aa;

		val[3][0] = aa;
		val[3][1] = aa;
		val[3][2] = aa;
		val[3][3] = aa;
	}

	Mat4x4(double aa, double ab, double ac, double ad,
		double ba, double bb, double bc, double bd,
		double ca, double cb, double cc, double cd,
		double da, double db, double dc, double dd) {
		val[0][0] = aa;
		val[0][1] = ab;
		val[0][2] = ac;
		val[0][3] = ad;

		val[1][0] = ba;
		val[1][1] = bb;
		val[1][2] = bc;
		val[1][3] = bd;

		val[2][0] = ca;
		val[2][1] = cb;
		val[2][2] = cc;
		val[2][3] = cd;

		val[3][0] = da;
		val[3][1] = db;
		val[3][2] = dc;
		val[3][3] = dd;

	}

	double at(int a, int b) {
		return val[a][b];
	}

	double a(int a, int b) {
		return val[a - 1][b - 1];
	}

	inline Mat4x4 operator+ (Mat4x4& B) {

		Mat4x4 tmp(val[0][0], val[0][1], val[0][2], val[0][3],
			val[1][0], val[1][1], val[1][2], val[1][3],
			val[2][0], val[2][1], val[2][2], val[2][3],
			val[3][0], val[3][1], val[3][2], val[3][3]);

		Mat4x4 res = Mat4x4(at(0, 0) + B.at(0, 0), at(0, 1) + B.at(0, 1), at(0, 2) + B.at(0, 2), at(0, 3) + B.at(0, 3),
			at(1, 0) + B.at(1, 0), at(1, 1) + B.at(1, 1), at(1, 2) + B.at(1, 2), at(1, 3) + B.at(1, 3),
			at(2, 0) + B.at(2, 0), at(2, 1) + B.at(2, 1), at(2, 2) + B.at(2, 2), at(2, 3) + B.at(2, 3),
			at(3, 0) + B.at(3, 0), at(3, 1) + B.at(3, 1), at(3, 2) + B.at(3, 2), at(3, 3) + B.at(3, 3) );
		return res;
	}

	inline Mat4x4 &Mat4x4::operator += (Mat4x4 & B) {
		val[0][0] += B.at(0, 0); 
		val[0][1] += B.at(0, 1); 
		val[0][2] += B.at(0, 2);
		val[0][3] += B.at(0, 3);
		val[1][0] += B.at(1, 0);
		val[1][1] += B.at(1, 1);
		val[1][2] += B.at(1, 2);
		val[1][3] += B.at(1, 3);
		val[2][0] += B.at(2, 0);
		val[2][1] += B.at(2, 1);
		val[2][2] += B.at(2, 2);
		val[2][3] += B.at(2, 3);
		val[3][0] += B.at(3, 0);
		val[3][1] += B.at(3, 1);
		val[3][2] += B.at(3, 2);
		val[3][3] += B.at(3, 3);
		
		return (*this);
	}

	double determinant() {
		return a(1, 4) * a(2, 3) * a(3, 2) * a(4, 1) - a(1, 3) * a(2, 4) * a(3, 2) * a(4, 1) 
			- a(1, 4) * a(2, 2) * a(3, 3) * a(4, 1) + a(1, 2) * a(2, 4) * a(3, 3) * a(4, 1) 
			+ a(1, 3) * a(2, 2) * a(3, 4) * a(4, 1) - a(1, 2) * a(2, 3) * a(3, 4) * a(4, 1) 
			- a(1, 4) * a(2, 3) * a(3, 1) * a(4, 2) + a(1, 3) * a(2, 4) * a(3, 1) * a(4, 2) 
			+ a(1, 4) * a(2, 1) * a(3, 3) * a(4, 2) - a(1, 1) * a(2, 4) * a(3, 3) * a(4, 2) 
			- a(1, 3) * a(2, 1) * a(3, 4) * a(4, 2) + a(1, 1) * a(2, 3) * a(3, 4) * a(4, 2) 
			+ a(1, 4) * a(2, 2) * a(3, 1) * a(4, 3) - a(1, 2) * a(2, 4) * a(3, 1) * a(4, 3) 
			- a(1, 4) * a(2, 1) * a(3, 2) * a(4, 3) + a(1, 1) * a(2, 4) * a(3, 2) * a(4, 3) 
			+ a(1, 2) * a(2, 1) * a(3, 4) * a(4, 3) - a(1, 1) * a(2, 2) * a(3, 4) * a(4, 3) 
			- a(1, 3) * a(2, 2) * a(3, 1) * a(4, 4) + a(1, 2) * a(2, 3) * a(3, 1) * a(4, 4) 
			+ a(1, 3) * a(2, 1) * a(3, 2) * a(4, 4) - a(1, 1) * a(2, 3) * a(3, 2) * a(4, 4) 
			- a(1, 2) * a(2, 1) * a(3, 3) * a(4, 4) + a(1, 1) * a(2, 2) * a(3, 3) * a(4, 4);
	}

	Mat4x4 inverse() {
		double det = determinant();
		Mat4x4 res = Mat4x4(at(0, 0) / det, at(0, 1) / det, at(0, 2) / det, at(0, 3) / det,
			at(1, 0) / det, at(1, 1) / det, at(1, 2) / det, at(1, 3) / det,
			at(2, 0) / det, at(2, 1) / det, at(2, 2) / det, at(2, 3) / det,
			at(3, 0) / det, at(3, 1) / det, at(3, 2) / det, at(3, 3) / det);
		return res;
	}

};

class Vertex {
public:
	double x, y, z;
	int index;

	Mat4x4 Q;
	Mesh* mesh;
	vector <Face*> faces;
	vector <Edge*> edges;

public:
	Vertex(double xx, double yy, double zz) { x = xx; y = yy; z = zz; Q = Mat4x4(0.); }
	double distance(Vertex * v1);
	bool makesEdgeWith(Vertex* v1);
};

class Edge {
public:
	Vertex* verts[2];
	vector <Face*> faces;
	Mesh* mesh;

	int index;

public:
	Edge(Vertex* v0, Vertex* v1, Mesh* meshIn) { verts[0] = v0; verts[1] = v1; mesh = meshIn; };
};

class Face {
public:
	Vertex* verts[3];
	Edge* edges[3];
	Mesh* mesh;

	Mat4x4 Quadric;
	icVector3 normal;
	int index;

public:
	Face(Vertex* v0, Vertex* v1, Vertex* v2, icVector3 norm) { verts[0] = v0; verts[1] = v1; verts[2] = v2; normal = norm; }
	void calcQuadric();
};

class Pair {
public:
	Vertex* v0;
	Vertex* v1;
	Edge* edge = NULL;
	int index;
	Mat4x4 Quadric;
	double Error;

public:
	Pair(Vertex* v0In, Vertex* v1In) { v0 = v0In; v1 = v1In; }
	Pair(Vertex* v0In, Vertex* v1In, Edge* edgeIn) { v0 = v0In; v1 = v1In; edge = edgeIn; }

	void updateQuadric();
	void updateError();
	icVector3 Vector();
	icVector3 QuadricVector();
	double QuadricError(icVector3 v);

};

class PairComparison
{
	bool reverse;
public:
	PairComparison(const bool& revparam = false)
	{
		reverse = revparam;
	}
	bool operator() ( Pair* lhs, Pair* rhs)
	{
		if (reverse) return (lhs->Error > rhs->Error);
		else return (lhs->Error < rhs->Error);
	}
};

class Mesh {
public:

	vector <Face*> flist;
	vector <Vertex*> vlist;
	vector <Edge*> elist;
	unsigned char orientation; //0=ccw, 1=cw
	icVector3 center;
	double radius;
	int maxFaces;

	vector <Pair*> validPairs;

public:

	/*constructors*/
	Mesh();
	Mesh(char*);

	/*utilties*/
	void seedValidPairs(double maxDistForPairs);
	void seedInitialQuadrics();

	/*initialization functions*/
	void initialize(double maxDistForPairs);
	void finalize();

	void simplify(int target);

};



#endif