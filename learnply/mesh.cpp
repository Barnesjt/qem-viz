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
#include <unordered_set>


// delimiters for parsing the obj file:

#define OBJDELIMS		" \t"

char* ReadRestOfLine(FILE*);
void ReadObjVTN(char*, int*, int*, int*);
void mjbCross(double[3], double[3], double[3]);
double mjbUnit(double[3]);
double mjbUnit(double[3], double[3]);
double mjbDot(double v1[3], double v2[3]);

struct face {
	int v, n, t;
};

//constructs a Mesh from obj file, adapted from Bailey's 450/550 loadobjfile.cpp
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
	Edges.clear();

	// open the input file:
	FILE* fp = fopen(file, "r");
	if (fp == NULL) fprintf(stderr, "Cannot open .obj file '%s'\n", file);
	else {
		for (; ; ) {
			char* line = ReadRestOfLine(fp);
			if (line == NULL) break;

			// skip this line if it is a comment:
			if (line[0] == '#') continue;

			// skip this line if it is something we don't feel like handling today:
			if (line[0] == 'g' || line[0] == 'm' || line[0] == 's' || line[0] == 'u') continue;

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

				bool uniqueVert = true;

				Vertex* newVert = NULL;

				if (newVert == NULL) {
					newVert = new Vertex(x, y, z);
				}
				
				Vertices.push_back(newVert);

				if (x < xmin) xmin = x;
				if (x > xmax) xmax = x;
				if (y < ymin) ymin = y;
				if (y > ymax) ymax = y;
				if (z < zmin) zmin = z;
				if (z > zmax) zmax = z;

				continue;
			}

			if (strcmp(cmd, "f") == 0) {
				struct face vertices[10];
				for (int i = 0; i < 10; i++) {
					vertices[i].v = 0;
					vertices[i].n = 0;
					vertices[i].t = 0;
				}

				int sizev = (int)Vertices.size();

				int numVertices = 0;
				bool valid = true;
				int vtx = 0;
				char* str;
				while ((str = strtok(NULL, OBJDELIMS)) != NULL) {
					int v, n, t;
					ReadObjVTN(str, &v, &t, &n);

					// if v, n, or t are negative, they are wrt the end of their respective list:
					if (v < 0) v += (sizev + 1);

					// be sure we are not out-of-bounds (<vector> will abort):
					if (v > sizev) {
						if (v != 0) fprintf(stderr, "Read vertex coord %d, but only have %d so far\n", v, sizev);
						v = 0;
						valid = false;
					}

					vertices[vtx].v = v;
					vtx++;

					if (vtx >= 10) break;

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


			if (strcmp(cmd, "s") == 0) continue;

		}

		//Remove duplicates from vectors (not strictly needed, but seems like a good idea?)
		std::unordered_set<Face*> sf;
		for (auto f : Faces)
			sf.insert(f);
		Faces.assign(sf.begin(), sf.end());
		std::sort(Faces.begin(), Faces.end());

		std::unordered_set<Vertex*> sv;
		for (auto v : Vertices)
			sv.insert(v);
		Vertices.assign(sv.begin(), sv.end());
		std::sort(Vertices.begin(), Vertices.end());

		std::unordered_set<Edge*> se;
		for (auto e : Edges)
			se.insert(e);
		Edges.assign(se.begin(), se.end());
		std::sort(Edges.begin(), Edges.end());

		orientation = 0; //Clockwise
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
	for (unsigned int i = 0; i < flist.size(); i++) flist.at(i)->index = i;
	for (unsigned int i = 0; i < vlist.size(); i++) vlist.at(i)->index = i;
	for (unsigned int i = 0; i < elist.size(); i++) elist.at(i)->index = i;
	//some initializing for QEM mesh simplification
	seedInitialQuadrics();
	seedValidPairs(maxDistForPairs);
}

void Mesh::finalize() {}

void Mesh::seedValidPairs(double t) {

	//remove duplicates from edge list first
	std::unordered_set<Edge*> se;
	for (auto f : elist)
		se.insert(f);
	elist.assign(se.begin(), se.end());
	std::sort(elist.begin(), elist.end()); //dunno if this is needed?

	//also remove duplicates from vertex list
	std::unordered_set<Vertex*> sv;
	for (auto f : vlist)
		sv.insert(f);
	vlist.assign(sv.begin(), sv.end());
	std::sort(vlist.begin(), vlist.end()); //dunno if this is needed?

	//Pair is valid if the verts make an edge
	for (auto e : elist) {
		Pair* p = new Pair(e->verts[0], e->verts[1], e);
		e->verts[0]->pairs.push_back(p);
		e->verts[1]->pairs.push_back(p);
		validPairs.push_back(p);
	}
	//Pair is valid if the verts their distance is less than t
	for (auto v0 : vlist) {
		for (auto v1 : vlist) {
			if (v0 != v1 && ((v0->x != v1->x) && (v0->y != v1->y) && (v0->z != v1->z)) &&   !v0->makesEdgeWith(v1) && v0->distance(v1) < t) {
				Pair* p = new Pair(v0, v1);
				v0->pairs.push_back(p);
				v1->pairs.push_back(p);
				validPairs.push_back(p);
			}
		}
	}
	for (auto p : validPairs) p->updateError();
}

void Mesh::seedInitialQuadrics() {
	for (auto f : flist) f->calcQuadric();
	for (auto v : vlist) {
		for (auto f : v->faces) {
			v->Q += f->Quadric;
		}
	}
}

// Seems this may be broken in some way
// Symptom: Face with 3 identical verts observed in Skeletal Foot object
// Unsure if this is an issue with the input file or our implementation,
//   but the symptom is only observed in this specific mesh
void Mesh::nonEdgeContract(Pair* target) {
	icVector3 newVert = target->Vector();
	Vertex* vBar = new Vertex(newVert.x, newVert.y, newVert.z);
	Vertex* remV0 = target->v0;
	Vertex* remV1 = target->v1;

	//now combine the 2 verticies face lists
	vBar->faces.insert(vBar->faces.end(), remV0->faces.begin(), remV0->faces.end());
	vBar->faces.insert(vBar->faces.end(), remV1->faces.begin(), remV1->faces.end());

	//remove dups
	std::unordered_set<Face*> sf;
	for (auto f : vBar->faces)
		sf.insert(f);
	vBar->faces.assign(sf.begin(), sf.end());
	std::sort(vBar->faces.begin(), vBar->faces.end());

	//now combine the 2 verticies edge lists
	vBar->edges.insert(vBar->edges.end(), remV0->edges.begin(), remV0->edges.end());
	vBar->edges.insert(vBar->edges.end(), remV1->edges.begin(), remV1->edges.end());

	//remove dups
	std::unordered_set<Edge*> se;
	for (auto f : vBar->edges)
		se.insert(f);
	vBar->edges.assign(se.begin(), se.end());
	std::sort(vBar->edges.begin(), vBar->edges.end());

	//next combine the 2 verticies pair lists
	vBar->pairs.insert(vBar->pairs.end(), remV0->pairs.begin(), remV0->pairs.end());
	vBar->pairs.insert(vBar->pairs.end(), remV1->pairs.begin(), remV1->pairs.end());

	//remove the dups
	std::unordered_set<Pair*> sp;
	for (auto f : vBar->pairs)
		sp.insert(f);
	vBar->pairs.assign(sp.begin(), sp.end());
	std::sort(vBar->pairs.begin(), vBar->pairs.end());

	//fix any pairs to point to the new vert
	for (auto p : vBar->pairs) {
		if (p->v0 == remV0 || p->v0 == remV1) p->v0 = vBar;
		if (p->v1 == remV0 || p->v1 == remV1) p->v1 = vBar;
	}

	//and remove pair that will be cut
	vBar->pairs.erase(std::remove(vBar->pairs.begin(), vBar->pairs.end(), target), vBar->pairs.end());

	//for all the faces, replace v0 and v1 with vBar
	//and if edge is remEF0b or remEF1b, replace with remEF0a or remEF0a
	for (auto f : vBar->faces) {
		for (int i = 0; i < 3; i++) {
			if (f->verts[i] == remV0 || f->verts[i] == remV1) {
				f->verts[i] = vBar;
			}
		}
		//now calc new normal, and new Quadric for each face
		f->calcNormal();
		f->calcQuadric();
	}
	//calc new Q for each vertex on those faces
	for (auto f : vBar->faces) {
		for (int i = 0; i < 3; i++) {
			auto v = f->verts[i];
			v->Q = glm::mat4x4(0.0);
			for (auto vF : v->faces) {
				v->Q += vF->Quadric;
			}
			//for each vert, loop through it's pairs and update
			for (auto p : v->pairs) {
				p->updateError();
			}
		}
	}

	//add new vert
	vlist.push_back(vBar);

	//from vlist: remV0 remV1
	vlist.erase(std::remove(vlist.begin(), vlist.end(), remV0), vlist.end());
	vlist.erase(std::remove(vlist.begin(), vlist.end(), remV1), vlist.end());

	//from validPairs: target
	validPairs.erase(std::remove(validPairs.begin(), validPairs.end(), target), validPairs.end());
}
/*
	icVector3 newVert = target->Vector();
	Vertex* vBar = new Vertex(newVert.x, newVert.y, newVert.z);
	Vertex* remV0 = target->v0;
	Vertex* remV1 = target->v1;

	//now combine the 2 verticies face lists
	vBar->faces.insert(vBar->faces.end(), remV0->faces.begin(), remV0->faces.end());
	vBar->faces.insert(vBar->faces.end(), remV1->faces.begin(), remV1->faces.end());

	//remove dups
	std::unordered_set<Face*> sf;
	for (auto f : vBar->faces)
		sf.insert(f);
	vBar->faces.assign(sf.begin(), sf.end());
	std::sort(vBar->faces.begin(), vBar->faces.end());

	//now combine the 2 verticies edge lists
	vBar->edges.insert(vBar->edges.end(), remV0->edges.begin(), remV0->edges.end());
	vBar->edges.insert(vBar->edges.end(), remV1->edges.begin(), remV1->edges.end());

	//remove dups
	std::unordered_set<Edge*> se;
	for (auto f : vBar->edges)
		se.insert(f);
	vBar->edges.assign(se.begin(), se.end());
	std::sort(vBar->edges.begin(), vBar->edges.end());

	//next combine the 2 verticies pair lists
	vBar->pairs.insert(vBar->pairs.end(), remV0->pairs.begin(), remV0->pairs.end());
	vBar->pairs.insert(vBar->pairs.end(), remV1->pairs.begin(), remV1->pairs.end());

	//remove the dups
	std::unordered_set<Pair*> sp;
	for (auto f : vBar->pairs)
		sp.insert(f);
	vBar->pairs.assign(sp.begin(), sp.end());
	std::sort(vBar->pairs.begin(), vBar->pairs.end());

	//fix any pairs to point to the new vert
	for (auto p : vBar->pairs) {
		if (p->v0 == remV0 || p->v0 == remV1) p->v0 = vBar;
		if (p->v1 == remV0 || p->v1 == remV1) p->v1 = vBar;
	}

	//and remove pair that will be cut
	vBar->pairs.erase(std::remove(vBar->pairs.begin(), vBar->pairs.end(), target), vBar->pairs.end());

	//for all the faces, replace v0 and v1 with vBar
	for (auto f : vBar->faces) {
		for (int i = 0; i < 3; i++) {
			if (f->verts[i] == remV0 || f->verts[i] == remV1) {
				f->verts[i] = vBar;
			}
		}
		//now calc new normal, and new Quadric for each face
		f->calcNormal();
		f->calcQuadric();
	}

	//add new vert
	vlist.push_back(vBar);

	//from vlist: remV0 remV1
	vlist.erase(std::remove(vlist.begin(), vlist.end(), remV0), vlist.end());
	vlist.erase(std::remove(vlist.begin(), vlist.end(), remV1), vlist.end());

	//from validPairs: target
	validPairs.erase(std::remove(validPairs.begin(), validPairs.end(), target), validPairs.end());

}*/

void Mesh::edge2fContract(Pair* target) {

	icVector3 newVert = target->Vector();
	Vertex* vBar = new Vertex(newVert.x, newVert.y, newVert.z);
	Vertex* remV0 = target->v0;
	Vertex* remV1 = target->v1;
	Edge* remE = target->edge;
	Face* remF0 = remE->faces.at(0);
	Face* remF1 = remE->faces.at(1);
	Edge* remEF0a = NULL;
	Edge* remEF0b = NULL;
	Edge* remEF1a = NULL;
	Edge* remEF1b = NULL;
	Vertex* oppVertF0 = NULL;
	Vertex* oppVertF1 = NULL;

	//find the opposite vertex for F0 and F1
	for (int i = 0; i < 3; i++) {
		if (remF0->verts[i] != remV0 || remF0->verts[i] != remV1) oppVertF0 = remF0->verts[i];
		if (remF1->verts[i] != remV0 || remF1->verts[i] != remV1) oppVertF1 = remF1->verts[i];
	}

	//find the 4 edges for the faces that will be removed
	for (int i = 0; i < 3; i++) {
		if (remF0->edges[i] != remE) {
			if (remEF0a == NULL) remEF0a = remF0->edges[i];
			else remEF0b = remF0->edges[i];
		}
		if (remF1->edges[i] != remE) {
			if (remEF1a == NULL) remEF1a = remF1->edges[i];
			else remEF1b = remF1->edges[i];
		}
	}

	//make new edges
	Edge* repEdgeF0 = new Edge(oppVertF0, vBar, this);
	Edge* repEdgeF1 = new Edge(oppVertF1, vBar, this);

	//in oppVertF0 remove F0, remEF0a, remEF0b. add repEdgeF0
	oppVertF0->faces.erase(std::remove(oppVertF0->faces.begin(), oppVertF0->faces.end(), remF0), oppVertF0->faces.end());
	oppVertF0->edges.erase(std::remove(oppVertF0->edges.begin(), oppVertF0->edges.end(), remEF0a), oppVertF0->edges.end());
	oppVertF0->edges.erase(std::remove(oppVertF0->edges.begin(), oppVertF0->edges.end(), remEF0b), oppVertF0->edges.end());
	oppVertF0->edges.push_back(repEdgeF0);

	//in oppVertF1 remove remF1, remEF1a, remEF1b. Add repEdgeF1
	oppVertF1->faces.erase(std::remove(oppVertF1->faces.begin(), oppVertF1->faces.end(), remF1), oppVertF1->faces.end());
	oppVertF1->edges.erase(std::remove(oppVertF1->edges.begin(), oppVertF1->edges.end(), remEF1a), oppVertF1->edges.end());
	oppVertF1->edges.erase(std::remove(oppVertF1->edges.begin(), oppVertF1->edges.end(), remEF1b), oppVertF1->edges.end());
	oppVertF0->edges.push_back(repEdgeF1);

	//find the adjacent faces to f0 and f1, add those faces to the replacement edges, add the replacement edge to those faces
	for(auto f : remEF0a->faces){
		if (f != remF0) {
			repEdgeF0->faces.push_back(f);
			for (int i = 0; i < 3; i++) {
				if (f->edges[i] == remEF0a) f->edges[i] = repEdgeF0;
			}
		}
	}
	for (auto f : remEF0b->faces) {
		if (f != remF0) {
			repEdgeF0->faces.push_back(f);
			for (int i = 0; i < 3; i++) {
				if (f->edges[i] == remEF0b) f->edges[i] = repEdgeF0;
			}
		}
	}
	for (auto f : remEF1a->faces) {
		if (f != remF1) {
			repEdgeF1->faces.push_back(f);
			for (int i = 0; i < 3; i++) {
				if (f->edges[i] == remEF1a) f->edges[i] = repEdgeF1;
			}
		}
	}
	for (auto f : remEF1b->faces) {
		if (f != remF1) {
			repEdgeF1->faces.push_back(f);
			for (int i = 0; i < 3; i++) {
				if (f->edges[i] == remEF1b) f->edges[i] = repEdgeF1;
			}
		}
	}

	//now combine the 2 verticies face lists
	vBar->faces.insert(vBar->faces.end(), remV0->faces.begin(), remV0->faces.end());
	vBar->faces.insert(vBar->faces.end(), remV1->faces.begin(), remV1->faces.end());

	//remove dups
	std::unordered_set<Face*> sf;
	for (auto f : vBar->faces)
		sf.insert(f);
	vBar->faces.assign(sf.begin(), sf.end());
	std::sort(vBar->faces.begin(), vBar->faces.end());

	//now combine the 2 verticies edge lists
	vBar->edges.insert(vBar->edges.end(), remV0->edges.begin(), remV0->edges.end());
	vBar->edges.insert(vBar->edges.end(), remV1->edges.begin(), remV1->edges.end());

	//remove dups
	std::unordered_set<Edge*> se;
	for (auto f : vBar->edges)
		se.insert(f);
	vBar->edges.assign(se.begin(), se.end());
	std::sort(vBar->edges.begin(), vBar->edges.end());

	//also remove the faces that will be cut
	vBar->faces.erase(std::remove(vBar->faces.begin(), vBar->faces.end(), remF0), vBar->faces.end());
	vBar->faces.erase(std::remove(vBar->faces.begin(), vBar->faces.end(), remF1), vBar->faces.end());

	//next combine the 2 verticies pair lists
	vBar->pairs.insert(vBar->pairs.end(), remV0->pairs.begin(), remV0->pairs.end());
	vBar->pairs.insert(vBar->pairs.end(), remV1->pairs.begin(), remV1->pairs.end());

	//remove the dups
	std::unordered_set<Pair*> sp;
	for (auto f : vBar->pairs)
		sp.insert(f);
	vBar->pairs.assign(sp.begin(), sp.end());
	std::sort(vBar->pairs.begin(), vBar->pairs.end());

	//fix any pairs to point to the new vert
	for (auto p : vBar->pairs) {
		if (p->v0 == remV0 || p->v0 == remV1) p->v0 = vBar;
		if (p->v1 == remV0 || p->v1 == remV1) p->v1 = vBar;
	}

	//and remove pair that will be cut
	vBar->pairs.erase(std::remove(vBar->pairs.begin(), vBar->pairs.end(), target), vBar->pairs.end());


	//for all the faces, replace v0 and v1 with vBar
	//and if edge is remEF0b or remEF1b, replace with remEF0a or remEF0a
	for (auto f : vBar->faces) {
		for (int i = 0; i < 3; i++) {
			if (f->verts[i] == remV0 || f->verts[i] == remV1){
				f->verts[i] = vBar;
			}
			if (f->edges[i] == remEF0b) {
				f->edges[i] = remEF0a;
			} else if (f->edges[i] == remEF1b) {
				f->edges[i] = remEF1a;
			}
		}
		//now calc new normal, and new Quadric for each face
		f->calcNormal();
		f->calcQuadric();
	}
	//calc new Q for each vertex on those faces
	for (auto f : vBar->faces) {
		for (int i = 0; i < 3; i++) {
			auto v = f->verts[i];
			v->Q = glm::mat4x4(0.0);
			for (auto vF : v->faces) {
				v->Q += vF->Quadric;
			}
			//for each vert, loop through it's pairs and update
			for (auto p : v->pairs) {
				p->updateError();
			}
		}
	}

	//add new vert
	vlist.push_back(vBar);

	//add new edge
	elist.push_back(repEdgeF0);
	elist.push_back(repEdgeF1);

	//remove degenerate geometry from Mesh lists: flist, elist, vlist, validPairs
	//from flist: remF0 remF1
	flist.erase(std::remove(flist.begin(), flist.end(), remF0), flist.end());
	flist.erase(std::remove(flist.begin(), flist.end(), remF1), flist.end());

	//from vlist: remV0 remV1
	vlist.erase(std::remove(vlist.begin(), vlist.end(), remV0), vlist.end());
	vlist.erase(std::remove(vlist.begin(), vlist.end(), remV1), vlist.end());

	//from elist: remE remEF0b remEF1b
	elist.erase(std::remove(elist.begin(), elist.end(), remE), elist.end());
	elist.erase(std::remove(elist.begin(), elist.end(), remEF0a), elist.end());
	elist.erase(std::remove(elist.begin(), elist.end(), remEF1a), elist.end());
	elist.erase(std::remove(elist.begin(), elist.end(), remEF0b), elist.end());
	elist.erase(std::remove(elist.begin(), elist.end(), remEF1b), elist.end());

	//from validPairs: target
	validPairs.erase(std::remove(validPairs.begin(), validPairs.end(), target), validPairs.end());

}

void Mesh::simplify(unsigned int target) {
	//don't perform simplify if target is too small
	//lower limit set to 64 as this is the minimum face count shown in the paper's figures
	if (target < 64) return;

	while (flist.size() > target) {
		//sort validPairs based on error cost
		std::sort(validPairs.begin(), validPairs.end(), PairComp());
		
		//take the best candidate
		Pair* p = validPairs.front();

		//if it's an edge do edge contraction, otherwise non-edge contraction
		if (p->edge != NULL) {
			if (p->edge->faces.size() == 2) edge2fContract(p);
			else validPairs.erase(std::remove(validPairs.begin(), validPairs.end(), p), validPairs.end()); //just remove the pair for now
				
		} else nonEdgeContract(p);
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

void Face::calcNormal() {
	double v01[3], v02[3], norm[3];
	v01[0] = verts[1]->x - verts[0]->x;
	v01[1] = verts[1]->y - verts[0]->y;
	v01[2] = verts[1]->z - verts[0]->z;
	v02[0] = verts[2]->x - verts[0]->x;
	v02[1] = verts[2]->y - verts[0]->y;
	v02[2] = verts[2]->z - verts[0]->z;

	mjbCross(v01, v02, norm);
	mjbUnit(norm, norm);

	normal.x = norm[0];
	normal.y = norm[1];
	normal.z = norm[2];


	/* Sanity check on the new normal, may be bad for rougher geometry (large normal changes)
	double oldNorm[3] = { normal.x, normal.y, normal.z };
	double dot = mjbDot(oldNorm, norm);
	if (dot > 0.) {
		normal.x = norm[0];
		normal.y = norm[1];
		normal.z = norm[2];
	} else {
		normal.x = norm[0] * -1.;
		normal.y = norm[1] * -1.;
		normal.z = norm[2] * -1.;
	}*/
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
	if (abs(glm::determinant(q) > .001)) { //may want a smaller value to check against?
		return QuadricVector();
	}
	//if cannot be computed from matrix, look along edge
	// our implementation just picks the midpoint, this could be optimized
	return icVector3((v0->x + v1->x) / 2., (v0->y + v1->y) / 2, (v0->z + v1->z) / 2.);
}

double Pair::QuadricError(icVector3 v) {
	glm::mat4x4 q = Quadric;
	//return 
	//	(v.x * q[0][0] * v.x + v.y * q[1][0] * v.x + v.z * q[2][0] * v.x + q[3][0] * v.x +
	//	v.x * q[0][1] * v.y + v.y * q[1][1] * v.y + v.z * q[2][1] * v.y + q[3][1] * v.y +
	//	v.x * q[0][2] * v.z + v.y * q[1][2] * v.z + v.z * q[2][2] * v.z + q[3][2] * v.z +
	//	v.x * q[0][3] + v.y * q[1][3] + v.z * q[2][3] + q[3][3]);
	return
		(v.x * q[0][0] * v.x + v.y * q[0][1] * v.x + v.z * q[0][2] * v.x + q[0][3] * v.x +
			v.x * q[1][0] * v.y + v.y * q[1][1] * v.y + v.z * q[1][2] * v.y + q[1][3] * v.y +
			v.x * q[2][0] * v.z + v.y * q[2][1] * v.z + v.z * q[2][2] * v.z + q[2][3] * v.z +
			v.x * q[3][0] + v.y * q[3][1] + v.z * q[3][2] + q[3][3]);
}

void Pair::updateError() {
	Error = QuadricError(QuadricVector());
}

char* ReadRestOfLine(FILE* fp){
	static char* line;
	std::vector<char> tmp(1000);
	tmp.clear();
	for (; ; ) {
		int c = getc(fp);
		if (c == EOF && tmp.size() == 0) return NULL;

		if (c == EOF || c == '\n') {
			delete[] line;
			line = new char[tmp.size() + 1];
			for (int i = 0; i < (int)tmp.size(); i++) line[i] = tmp[i];
			line[tmp.size()] = '\0';	// terminating null
			return line;
		}
		else tmp.push_back(c);
	}
	return "";
}

void ReadObjVTN(char* str, int* v, int* t, int* n){
	// can be one of v, v//n, v/t, v/t/n:

	if (strstr(str, "//")) {			// v//n
		*t = 0;
		sscanf(str, "%d//%d", v, n);
		return;
	} else if (sscanf(str, "%d/%d/%d", v, t, n) == 3){	// v/t/n
		return;
	} else{
		*n = 0;
		if (sscanf(str, "%d/%d", v, t) == 2) return;		// v/t
		else {						// v
			*n = *t = 0;
			sscanf(str, "%d", v);
		}
	}
}

double mjbDot(double v1[3], double v2[3]) {
	return v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2];
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