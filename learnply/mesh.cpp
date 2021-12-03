#define _CRT_SECURE_NO_WARNINGS

#include <math.h>
#include <stdio.h>
#include <fstream>
#include "icVector.H"
#include "icMatrix.H"
#include "mesh.h"

#include <vector>
#include <algorithm>
#include <queue>


// delimiters for parsing the obj file:

#define OBJDELIMS		" \t"

char* ReadRestOfLine(FILE*);
void ReadObjVTN(char*, int*, int*, int*);
void mjbCross(double[3], double[3], double[3]);
double mjbUnit(double[3]);
double mjbUnit(double[3], double[3]);

struct face {
	int v, n, t;
};


Mesh::Mesh(char* file) {

	char* cmd;		// the command string
	char* str;		// argument string

	double xmin = std::numeric_limits<double>::max();
	double ymin = std::numeric_limits<double>::max();
	double zmin = std::numeric_limits<double>::max();
	double xmax = std::numeric_limits<double>::min();
	double ymax = std::numeric_limits<double>::min();
	double zmax = std::numeric_limits<double>::min();

	std::vector <Vertex*> Vertices;
	std::vector <Face*> Faces;
	std::vector <Edge*> Edges;

	Vertices.clear();
	Faces.clear();

	// open the input file:

	FILE* fp = fopen(file, "r");
	if (fp == NULL) {
		fprintf(stderr, "Cannot open .obj file '%s'\n", file);
	}
	else {
		for (; ; ) {
			char* line = ReadRestOfLine(fp);
			if (line == NULL)
				break;

			// skip this line if it is a comment:
			if (line[0] == '#')
				continue;

			// skip this line if it is something we don't feel like handling today:
			if (line[0] == 'g' || line[0] == 'm' || line[0] == 's' || line[0] == 'u')
				continue;

			// get the command string:
			cmd = strtok(line, OBJDELIMS);


			// skip this line if it is empty:
			if (cmd == NULL) continue;

			// skip vertex normals
			if (strcmp(cmd, "vn") == 0) continue;

			// skip texture coords
			if (strcmp(cmd, "vt") == 0) continue;

			if (strcmp(cmd, "v") == 0) {
				str = strtok(NULL, OBJDELIMS);
				double x = atof(str);

				str = strtok(NULL, OBJDELIMS);
				double y = atof(str);

				str = strtok(NULL, OBJDELIMS);
				double z = atof(str);

				Vertex* newVert = new Vertex(x, y, z);

				Vertices.push_back(newVert);

				if (x < xmin) xmin = x;
				if (x > xmax) xmax = x;
				if (y < ymin) ymin = y;
				if (y > ymax) ymax = y;
				if (z < zmin) zmin = z;
				if (z > zmax) zmax = z;

				continue;
			}

			if (strcmp(cmd, "f") == 0)
			{
				struct face vertices[10];
				for (int i = 0; i < 10; i++)
				{
					vertices[i].v = 0;
					vertices[i].n = 0;
					vertices[i].t = 0;
				}

				int sizev = (int)Vertices.size();

				int numVertices = 0;
				bool valid = true;
				int vtx = 0;
				char* str;
				while ((str = strtok(NULL, OBJDELIMS)) != NULL)
				{
					int v, n, t;
					ReadObjVTN(str, &v, &t, &n);

					// if v, n, or t are negative, they are wrt the end of their respective list:

					if (v < 0)
						v += (sizev + 1);

					// be sure we are not out-of-bounds (<vector> will abort):

					if (v > sizev)
					{
						if (v != 0)
							fprintf(stderr, "Read vertex coord %d, but only have %d so far\n", v, sizev);
						v = 0;
						valid = false;
					}

					vertices[vtx].v = v;
					//vertices[vtx].n = n;
					//vertices[vtx].t = t;
					vtx++;

					if (vtx >= 10)
						break;

					numVertices++;
				}

				// if vertices are invalid, don't draw anything this time:
				if (!valid || (numVertices < 3)) continue;

				// list the vertices:
				int numTriangles = numVertices - 2;

				for (int it = 0; it < numTriangles; it++) {
					int vv[3];
					vv[0] = 0;
					vv[1] = it + 1;
					vv[2] = it + 2;

					// get the planar normal, in case vertex normals are not defined:

					Vertex* v0 = Vertices[vertices[vv[0]].v - 1];
					Vertex* v1 = Vertices[vertices[vv[1]].v - 1];
					Vertex* v2 = Vertices[vertices[vv[2]].v - 1];

					double v01[3], v02[3], norm[3];
					v01[0] = v1->x - v0->x;
					v01[1] = v1->y - v0->y;
					v01[2] = v1->z - v0->z;
					v02[0] = v2->x - v0->x;
					v02[1] = v2->y - v0->y;
					v02[2] = v2->z - v0->z;
					mjbCross(v01, v02, norm);
					mjbUnit(norm, norm);

					Face* newFace = new Face(v0, v1, v2, icVector3(norm));
					Edge* e01 = NULL;
					Edge* e12 = NULL;
					Edge* e20 = NULL;

					for (auto e : Edges) {
						if ((e->verts[0] == v0 && e->verts[1] == v1) || (e->verts[1] == v0 && e->verts[0] == v1))
							e01 = e;
						if ((e->verts[0] == v1 && e->verts[1] == v2) || (e->verts[1] == v1 && e->verts[0] == v2))
							e12 = e;
						if ((e->verts[0] == v2 && e->verts[1] == v0) || (e->verts[1] == v2 && e->verts[0] == v0))
							e20 = e;
					}

					if (e01 == NULL) {
						e01 = new Edge(v0, v1, this);
						Edges.push_back(e01);
					}
					if (e12 == NULL) { 
						e12 = new Edge(v1, v2, this);
						Edges.push_back(e12);
					}
					if (e20 == NULL) {
						e20 = new Edge(v2, v0, this);
						Edges.push_back(e20);
					}

					v0->mesh = this;
					v1->mesh = this;
					v2->mesh = this;

					v0->edges.push_back(e01);
					v0->edges.push_back(e20);
					v1->edges.push_back(e01);
					v1->edges.push_back(e12);
					v2->edges.push_back(e12);
					v2->edges.push_back(e20);

					v0->faces.push_back(newFace);
					v1->faces.push_back(newFace);
					v2->faces.push_back(newFace);

					e01->faces.push_back(newFace);
					e12->faces.push_back(newFace);
					e20->faces.push_back(newFace);

					newFace->edges[0] = e01;
					newFace->edges[1] = e12;
					newFace->edges[2] = e20;

					newFace->mesh = this;

					Faces.push_back(newFace);
				}
				continue;
			}


			if (strcmp(cmd, "s") == 0)
			{
				continue;
			}

		}

		orientation = 0;
		flist = Faces;
		vlist = Vertices;
		elist = Edges;
		center = icVector3((xmin + xmax)/2., (ymin + ymax)/2., (zmin + zmax)/2.);
		radius = std::max({ abs(center.x - xmax), abs(center.y - ymax), abs(center.z - zmax) });

		fclose(fp);
	}
}

void Mesh::initialize(double maxDistForPairs) {
	//set orientation to clockwise (CW). Used for normals for lighting
	orientation = 0;
	maxFaces = flist.size();
	//set indexes for vlist and flist, not strictly required, but may be convienient later
	for (int i = 0; i < flist.size(); i++) flist.at(i)->index = i;
	for (int i = 0; i < vlist.size(); i++) vlist.at(i)->index = i;
	for (int i = 0; i < elist.size(); i++) elist.at(i)->index = i;
	//some initializing for QEM mesh simplification
	seedInitialQuadrics();
	seedValidPairs(maxDistForPairs);
	
}

void Mesh::finalize() {
	for (auto f : flist) free(f->verts);
	flist.clear();
	for (auto v : vlist) free(v);
	vlist.clear();
	for (auto e : elist) free(e);
	elist.clear();
	for (auto p : validPairs) free(p);
	validPairs.clear();

}

void Mesh::seedValidPairs(double t) {
	//Pair is valid if the verts make an edge
	for (auto e : elist) {
		validPairs.push_back(new Pair( e->verts[0], e->verts[1], e));
	}
	//Pair is valid if the verts their distance is less than t
	for (auto v0 : vlist) {
		for (auto v1 : vlist) {
			if (v0 != v1 && !v0->makesEdgeWith(v1) && v0->distance(v1) < t) {
				validPairs.push_back(new Pair(v0, v1));
			}
		}
	}
	for (auto p : validPairs) {
		p->updateError();
	}
}

void Mesh::seedInitialQuadrics() {
	for (auto f : flist) {
		f->calcQuadric();
		//for (int i = 0; i < 3; i++) {
		//	f->verts[i]->Q += f->Quadric;
		//}
	}
	for (auto v : vlist) {
		for (auto f : v->faces) {
			v->Q += f->Quadric;
		}
	}
}

void Mesh::simplify(int target) {
	if (target > flist.size() || target < 4) return;

	for (int i = 0; i < validPairs.size(); i++) {
		validPairs.at(i)->index = i;
	}

	std::sort(validPairs.begin(), validPairs.end(), [](Pair* a, Pair* b) {return a->Error < b->Error; });

	while (flist.size() > target) {

		Pair* p = validPairs.front();

		for (int i = 0; i < flist.size(); i++) {
			flist.at(i)->index = i;
		}

		if (p->edge != NULL) {

			icVector3 newVert = p->Vector();
			Vertex* v0 = p->v0;
			Vertex* v1 = p->v1;

			Vertex* vBar = new Vertex(newVert.x, newVert.y, newVert.z);
			vBar->mesh = v0->mesh;
			vBar->index = v0->index;
			vBar->Q = v0->Q;
			for (int i = 0; i < flist.size(); i++) {
				flist.at(i)->index = i;
			}
			Face* f1 = NULL;
			Face* f2 = NULL;
			for (auto f : v0->faces) {
				bool correctFace = false;
				for (int i = 0; i < 3; i++) {
					if (f->edges[i] == p->edge) {
						correctFace = true;
						if (f->verts[i] == v0 || f->verts[i] == v1) f->verts[i] = vBar;
					}
				}
				if(correctFace) {
					if (f1 == NULL) f1 = f;
					else if (f2 == NULL) f2 = f;	
				}
			}
			std::vector<Face*>::iterator pos1 = std::find(flist.begin(), flist.end(), f1);
			if (pos1 != flist.end()) {
				//Other things to fix when removing face?
				this->flist.erase(pos1);
			} else {
				int g = 2;
			}
			std::vector<Face*>::iterator pos2 = std::find(flist.begin(), flist.end(), f2);
			if (pos2 != flist.end()) {
				//Other things to fix when removing face?
				this->flist.erase(pos2);
			} else {
				int g = 2;
			}



			//perform edge contraction

			/*
			Mesh:
				Edge from elist
				Remove redundant Vertex from vlist
				Remove 2x faces from flist
				Remove pair from validPairs

			Face:
				Delete 2x faces
				All adjacent faces Need to be updated:
					One of the verts for each Face is moved (recompute normal? and Quadric?)
					The other edges from deleted faces are now shared, remove redundancies?
			
			Edge:
				Some edges now overlap

			Vertex:

			Pair:
			
			*/


		} else {
			//perform non-edge contraction
			icVector3 newVert = p->Vector();
			Vertex* v0 = p->v0;
			Vertex* v1 = p->v1;
			Vertex* vBar = new Vertex(newVert.x, newVert.y, newVert.z);
			vBar->mesh = v0->mesh;
			vBar->index = v0->index;
			vBar->Q = v0->Q;

			for (auto f : v0->faces) {
				for (int i = 0; i < 3; i++) {
					if (f->verts[i] == v0)
						f->verts[i] = vBar;
				}
			}

			for (auto f : v1->faces) {
				for (int i = 0; i < 3; i++) {
					if (f->verts[i] == v1)
						f->verts[i] = vBar;
				}
			}

		}
		validPairs.erase(validPairs.begin());

	}
}

double Vertex::distance(Vertex * v1) {
	Vertex* v0 = this;
	double dist = sqrt(pow(v1->x-v0->x, 2.)+ pow(v1->y - v0->y, 2.)+ pow(v1->z - v0->z, 2.));
	return dist;
}

bool Vertex::makesEdgeWith(Vertex* v1) {
	Vertex* v0 = this;
	//If they are the same, then they do not make an edge
	if (v0 == v1) return false;
	//for all edges, if v1 is not one of the verts, return false
	for (auto e : v0->edges) {
		if (e->verts[0] == v1 || e->verts[1] == v0) return true;
	}
	return false;
}

void Face::calcQuadric() {
	Vertex* v1 = verts[0];
	Vertex* v2 = verts[1];
	Vertex* v3 = verts[2];

	double x = v1->x;
	double y = v1->y;
	double z = v1->z;

	double a = normal.x;
	double b = normal.y;
	double c = normal.z;
	double d = -a * x - b * y - c * z;

	Quadric = glm::mat4x4(a*a, a*b, a*c, a*d, a*b, b*b, b*c, b*d, a*c, b*c, c*c, c*d, a*d, b*d, c*d, d*d);
}

void Pair::updateQuadric() {
	Quadric = v0->Q + v1->Q;
}

icVector3 Pair::QuadricVector() {
	updateQuadric();
	glm::mat4x4 q = Quadric;
	q[3][0] = 0.;
	q[3][1] = 0.;
	q[3][2] = 0.; 
	q[3][3] = 1.;
	q = glm::inverse(q);
	return icVector3(q[0][3], q[1][3], q[2][3]);
}

icVector3 Pair::Vector() {
	updateQuadric();
	glm::mat4x4 q = Quadric;
	if (abs(glm::determinant(q) > .001)) {
		return QuadricVector();
	}
	//if cannot be computed from matrix, look along edge
	//For now this is just halfway along the edge 
	//TODO: Fix
	return icVector3((v0->x + v1->x) / 2., (v0->y + v1->y) / 2, (v0->z + v1->z) / 2.);
}

double Pair::QuadricError(icVector3 v) {
	glm::mat4x4 q = Quadric;
	return 
		(v.x * q[0][0] * v.x + v.y * q[1][0] * v.x + v.z * q[2][0] * v.x + q[3][0] * v.x +
		v.x * q[0][1] * v.y + v.y * q[1][1] * v.y + v.z * q[2][1] * v.y + q[3][1] * v.y +
		v.x * q[0][2] * v.z + v.y * q[1][2] * v.z + v.z * q[2][2] * v.z + q[3][2] * v.z +
		v.x * q[0][3] + v.y * q[1][3] + v.z * q[2][3] + q[3][3]);
}

void Pair::updateError() {
	//could be cached and updated as needed, but we'll just calc it in full each time for now
	Error = QuadricError(QuadricVector());
}

char* ReadRestOfLine(FILE* fp){
	static char* line;
	std::vector<char> tmp(1000);
	tmp.clear();

	for (; ; )
	{
		int c = getc(fp);

		if (c == EOF && tmp.size() == 0)
		{
			return NULL;
		}

		if (c == EOF || c == '\n')
		{
			delete[] line;
			line = new char[tmp.size() + 1];
			for (int i = 0; i < (int)tmp.size(); i++)
			{
				line[i] = tmp[i];
			}
			line[tmp.size()] = '\0';	// terminating null
			return line;
		}
		else
		{
			tmp.push_back(c);
		}
	}

	return "";
}

void ReadObjVTN(char* str, int* v, int* t, int* n){
	// can be one of v, v//n, v/t, v/t/n:

	if (strstr(str, "//"))				// v//n
	{
		*t = 0;
		sscanf(str, "%d//%d", v, n);
		return;
	}
	else if (sscanf(str, "%d/%d/%d", v, t, n) == 3)	// v/t/n
	{
		return;
	}
	else
	{
		*n = 0;
		if (sscanf(str, "%d/%d", v, t) == 2)		// v/t
		{
			return;
		}
		else						// v
		{
			*n = *t = 0;
			sscanf(str, "%d", v);
		}
	}
}

void mjbCross(double v1[3], double v2[3], double vout[3]) {

	double tmp[3];

	tmp[0] = v1[1] * v2[2] - v2[1] * v1[2];
	tmp[1] = v2[0] * v1[2] - v1[0] * v2[2];
	tmp[2] = v1[0] * v2[1] - v2[0] * v1[1];

	vout[0] = tmp[0];
	vout[1] = tmp[1];
	vout[2] = tmp[2];
}

double mjbUnit(double v[3]) {
	double dist = v[0] * v[0] + v[1] * v[1] + v[2] * v[2];

	if (dist > 0.0) {
		dist = sqrt(dist);
		v[0] /= dist;
		v[1] /= dist;
		v[2] /= dist;
	}

	return dist;
}

double mjbUnit(double vin[3], double vout[3]) {
	double dist = vin[0] * vin[0] + vin[1] * vin[1] + vin[2] * vin[2];

	if (dist > 0.0) {
		dist = sqrt(dist);
		vout[0] = vin[0] / dist;
		vout[1] = vin[1] / dist;
		vout[2] = vin[2] / dist;
	} else {
		vout[0] = vin[0];
		vout[1] = vin[1];
		vout[2] = vin[2];
	}

	return dist;
}