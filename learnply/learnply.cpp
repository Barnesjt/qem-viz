#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <fstream>
#include <vector>
#include <string>
#include <iostream>

#include <Eigen/Core>
#include <Eigen/Eigenvalues>
#include "glm/mat4x4.hpp"
#include "glm/gtc/type_ptr.hpp"
#include "glm/gtc/matrix_transform.hpp"

#include "glError.h"
#include "gl/glew.h"
#include "gl/freeglut.h"
#include "icVector.H"
#include "icMatrix.H"
#include "trackball.h"
#include "tmatrix.h"

#include "drawUtil.h"

// Data structure to hold Mesh data, includes code adapted from: http://web.engr.oregonstate.edu/~mjb/cs550/loadobjfile.cpp
#include "mesh.h"

//Added include, specifically for reading directory contents
#include <Windows.h>

//Added natural sort from: https://github.com/scopeInfinity/NaturalSort
//For sorting file names into best order (e.g. r2.ply before r10.ply)
#include "natural_sort.hpp"


//Vectors and Strings
using std::vector;
using std::string;

//Eigen
using namespace Eigen;

//Vector of Meshes instead (Because why not just load everything at once)
//Individual vectors for different data, this could be rolled into a single class, maybe later
vector<string> all_obj_files;
vector<Mesh*> all_meshes;

//The marker for which ply we're looking at, used to index into the vectors and cycle through them
unsigned int curr_mesh = 0;
unsigned int face_target = 0;
double t_ratio = 15.;
bool view_error = false;

struct Ellipsoid {
	double x,y,z;
	double xS,yS,zS,wS;
	glm::dmat4x4 totalT;
	glm::dmat4x4 rotT;
	glm::dmat4x4 scaleT;
	glm::dmat4x4 transT;

};

vector<Ellipsoid*> errorEllip;

//To load all geometry use 
string OBJ_PATH = "../data/geometry/";

//string OBJ_PATH = "../data/geometry/cow/";

Mesh* mesh;

/*scene related variables*/
const float zoomspeed = 0.9f;
int win_width = 1024;
int win_height = 1024;
float aspectRatio = (float)win_width / (float)win_height;
const int view_mode = 0;		// 0 = othogonal, 1=perspective
const double radius_factor = 0.9;

int display_mode = 1;

/*User Interaction related variabes*/
float s_old, t_old;
float rotmat[4][4];
double zoom = 1.0;
double translation[2] = { 0, 0 };
int mouse_mode = -2;	// -1 = no action, 1 = translate y, 2 = rotate

/******************************************************************************
Forward declaration of functions
******************************************************************************/

void init(void);

/*glut attaching functions*/
void keyboard(unsigned char key, int x, int y);
void motion(int x, int y);
void display(void);
void mouse(int button, int state, int x, int y);
void mousewheel(int wheel, int direction, int x, int y);
void reshape(int width, int height);

/*display vis results*/
void display_mesh(Mesh* mesh);

/*added functions*/
string dispModeString(int mode);
void DoRasterString(float x, float y, float z, char* s);
std::vector<std::string> get_all_obj_in_folder(std::string folder);

void generate_ellipsoids(Mesh* mesh);
void ellipsoid_transformation(Vertex* vert, Ellipsoid*);
void display_ellipsoid(Ellipsoid* ellip);


/******************************************************************************
Main program.
******************************************************************************/

int main(int argc, char* argv[]) {
	//Get all obj files in directory
	all_obj_files = get_all_obj_in_folder(OBJ_PATH);

	//Natural sorting of obj files
	std::sort(all_obj_files.begin(), all_obj_files.end(), SI::natural::compare<string>);

	for (std::string i : all_obj_files) {
		all_meshes.push_back(new Mesh(const_cast<char*>((OBJ_PATH + i).c_str())));
		all_meshes.back()->initialize(all_meshes.back()->radius/t_ratio);
		//all_meshes.back()->initialize(0.0);
	}

	//if no ply files, then exit with error. otherwise set the first ply
	if (all_meshes.empty()) {
		return 0; //no meshes, no program
	} else {
		mesh = all_meshes.at(0);
		curr_mesh = 0;
		face_target = mesh->flist.size();
	}

	/*init glut and create window*/
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB);
	glutInitWindowPosition(20, 20);
	glutInitWindowSize(win_width, win_height);
	glutCreateWindow("CS553 - Group 1 - QEM w/ Ellipsoids");

	/*initialize openGL*/
	init();
	
	/*the render function and callback registration*/
	glutKeyboardFunc(keyboard);
	glutReshapeFunc(reshape);
	glutDisplayFunc(display);
	glutMotionFunc(motion);
	glutMouseFunc(mouse);
	glutMouseWheelFunc(mousewheel);
	
	/*event processing loop*/
	glutMainLoop();
	
	/*clear memory before exit*/
	mesh->finalize();	// finalize everything
	return 0;
}

/******************************************************************************
Get all .obj files

Function adapted from: https://stackoverflow.com/a/20847429
******************************************************************************/
// Function adapted from: https://stackoverflow.com/a/20847429
std::vector<std::string> get_all_obj_in_folder(std::string folder) {
	std::vector<std::string> names;
	std::string search_path = folder + "/*.obj";
	WIN32_FIND_DATA fd;
	HANDLE hFind = ::FindFirstFile(search_path.c_str(), &fd);
	if (hFind != INVALID_HANDLE_VALUE) {
		do {
			// read all (real) files in current folder
			// , delete '!' read other 2 default folder . and ..
			if (!(fd.dwFileAttributes & FILE_ATTRIBUTE_DIRECTORY)) {
				names.push_back(fd.cFileName);
			}
		} while (::FindNextFile(hFind, &fd));
		::FindClose(hFind);
	}
	return names;
}

/******************************************************************************
Use glut to display a string of characters using a raster font

Borrowed from CS450 Sample Code
******************************************************************************/
void
DoRasterString(float x, float y, float z, char* s) {
	glRasterPos3f((GLfloat)x, (GLfloat)y, (GLfloat)z);
	char c;			// one character to print
	for (; (c = *s) != '\0'; s++) {
		glutBitmapCharacter(GLUT_BITMAP_TIMES_ROMAN_24, c);
	}
}

/******************************************************************************
Convert a display mode number to a string
******************************************************************************/
string dispModeString(int mode) {
	switch (mode) {
	case 1:
		return "Solid";
	case 2:
		return "Solid w/ Error Ellipsoids";
	default:
		return "Unrecognized";
	}
}

/******************************************************************************
Set projection mode
******************************************************************************/

void set_view(GLenum mode) {
	GLfloat light_ambient0[] = { 0.3f, 0.3f, 0.3f, 1.0f };
	GLfloat light_diffuse0[] = { 0.7f, 0.7f, 0.7f, 1.0f };
	GLfloat light_specular0[] = { 0.0f, 0.0f, 0.0f, 1.0f };

	GLfloat light_ambient1[] = { 0.0f, 0.0f, 0.0f, 1.0f };
	GLfloat light_diffuse1[] = { 0.5f, 0.5f, 0.5f, 1.0f };
	GLfloat light_specular1[] = { 0.0f, 0.0f, 0.0f, 1.0f };

	glLightfv(GL_LIGHT0, GL_AMBIENT, light_ambient0);
	glLightfv(GL_LIGHT0, GL_DIFFUSE, light_diffuse0);
	glLightfv(GL_LIGHT0, GL_SPECULAR, light_specular0);

	glLightfv(GL_LIGHT1, GL_AMBIENT, light_ambient1);
	glLightfv(GL_LIGHT1, GL_DIFFUSE, light_diffuse1);
	glLightfv(GL_LIGHT1, GL_SPECULAR, light_specular1);


	glMatrixMode(GL_PROJECTION);
	if (mode == GL_RENDER)
		glLoadIdentity();

	if (aspectRatio >= 1.0) {
		if (view_mode == 0)
			glOrtho(-radius_factor * zoom * aspectRatio, radius_factor * zoom * aspectRatio, -radius_factor * zoom, radius_factor * zoom, -1000, 1000);
		else
			glFrustum(-radius_factor * zoom * aspectRatio, radius_factor * zoom * aspectRatio, -radius_factor * zoom, radius_factor * zoom, 0.1, 1000);
	}
	else {
		if (view_mode == 0)
			glOrtho(-radius_factor * zoom, radius_factor * zoom, -radius_factor * zoom / aspectRatio, radius_factor * zoom / aspectRatio, -1000, 1000);
		else
			glFrustum(-radius_factor * zoom, radius_factor * zoom, -radius_factor * zoom / aspectRatio, radius_factor * zoom / aspectRatio, 0.1, 1000);
	}

	GLfloat light_position[3];
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	light_position[0] = 5.5f;
	light_position[1] = 0.0f;
	light_position[2] = 0.0f;
	glLightfv(GL_LIGHT0, GL_POSITION, light_position);
	light_position[0] = -0.1f;
	light_position[1] = 0.0f;
	light_position[2] = 0.0f;
	glLightfv(GL_LIGHT2, GL_POSITION, light_position);
}

/******************************************************************************
Update the scene
******************************************************************************/

void set_scene(GLenum mode, Mesh* mesh){
	glTranslatef((float)translation[0], (float)translation[1], -3.0f);

	/*multiply rotmat to current mat*/
	{
		int i, j, index = 0;

		GLfloat mat[16];

		for (i = 0; i < 4; i++)
			for (j = 0; j < 4; j++)
				mat[index++] = rotmat[i][j];

		glMultMatrixf(mat);
	}

	//glScalef(0.9, 0.9, 0.9);
	//glTranslatef(0.0, 0.0, 0.0);

	glScalef(0.9f / (float)mesh->radius, 0.9f / (float)mesh->radius, 0.9f / (float)mesh->radius);
	glTranslatef((float)-mesh->center.entry[0], (float)-mesh->center.entry[1], (float)-mesh->center.entry[2]);
}

/******************************************************************************
Init scene
******************************************************************************/

void init(void) {

	mat_ident(rotmat);

	/* select clearing color */
	glClearColor(0.0f, 0.0f, 0.0f, 0.0f);  // background
	glShadeModel(GL_FLAT);
	glPolygonMode(GL_FRONT, GL_FILL);

	glDisable(GL_DITHER);
	glEnable(GL_DEPTH_TEST);
	glDepthFunc(GL_LESS);
	
	//set pixel storage modes
	glPixelStorei(GL_PACK_ALIGNMENT, 1);
	
	glEnable(GL_NORMALIZE);
	if (mesh->orientation == 0)
		glFrontFace(GL_CW);
	else
		glFrontFace(GL_CCW);
}

/******************************************************************************
Process a keyboard action.  In particular, exit the program when an
"escape" is pressed in the window.
******************************************************************************/

void keyboard(unsigned char key, int x, int y) {

	switch (key) {
	case 27:	// set excape key to exit program
		for (auto m : all_meshes) m->finalize(); // finalize_everything
		exit(0);
		break;

	case 'q':
		unsigned int newFaceTarget;
		std::cout << "Input New Face Target: ";
		std::cin >> newFaceTarget;
		if (newFaceTarget < 4) std::cout << "Too few faces, press 'q' to try again";
		else if (newFaceTarget > mesh->maxFaces) std::cout << "Too many faces, press 'q' to try again";
		else face_target = newFaceTarget;

		if (face_target > mesh->flist.size()) {
			free(all_meshes.at(curr_mesh));
			all_meshes.at(curr_mesh) = new Mesh(const_cast<char*>((OBJ_PATH + all_obj_files.at(curr_mesh)).c_str()));
			mesh = all_meshes.at(curr_mesh);
			mesh->initialize(mesh->radius / t_ratio);
			//mesh->initialize(0.0); //edge contractions only: t=0.0
			//generate_ellipsoids(mesh);
		}
		if (face_target < mesh->flist.size()){
			mesh->simplify(face_target);
			generate_ellipsoids(mesh);
		}
		glutPostRedisplay();
		break;

	case 'e':
		view_error = !view_error;
		glutPostRedisplay();
		break;

	case 'r':	// reset rotation and transformation
		mat_ident(rotmat);
		translation[0] = 0;
		translation[1] = 0;
		zoom = 1.0;
		glutPostRedisplay();
		break;

	case 't':
		//toggle through all ply files, increments and modulo on the vector size
		curr_mesh = ++curr_mesh % all_meshes.size();

		//reset the poly we are displaying, w/ min and max
		mesh = all_meshes.at(curr_mesh);

		face_target = mesh->flist.size();

		//recall keyboard function to make sure the correct display mode is set
		keyboard('0' + display_mode, x, y);
		glutPostRedisplay();
		break;
	}
}

/******************************************************************************
Callback function for dragging mouse
******************************************************************************/

void motion(int x, int y) {
	float r[4];
	float s, t;

	s = (2.0f * (float)x - (float)win_width) / (float)win_width;
	t = (2.0f * ((float)win_height - (float)y) - (float)win_height) / (float)win_height;

	if ((s == s_old) && (t == t_old))
		return;

	switch (mouse_mode) {
	case 2:

		float rvec[4];

		mat_to_quat(rotmat, rvec);
		trackball(r, s_old, t_old, s, t);
		add_quats(r, rvec, rvec);
		quat_to_mat(rvec, rotmat);

		s_old = s;
		t_old = t;

		display();
		break;

	case 1:

		translation[0] += ((double)s - (double)s_old);
		translation[1] += ((double)t - (double)t_old);

		s_old = s;
		t_old = t;

		display();
		break;
	}
}

/******************************************************************************
Callback function for mouse clicks
******************************************************************************/

void mouse(int button, int state, int x, int y) {

	int key = glutGetModifiers();

	if (button == GLUT_LEFT_BUTTON || button == GLUT_RIGHT_BUTTON) {
		
		if (state == GLUT_DOWN) {
			float xsize = (float)win_width;
			float ysize = (float)win_height;

			float s = (2.0f * (float) x - (float)win_width) / (float)win_width;
			float t = (2.0f * ((float)win_height - (float)y) - (float)win_height) / (float)win_height;

			s_old = s;
			t_old = t;

			/*translate*/
			if (button == GLUT_LEFT_BUTTON)
			{
				mouse_mode = 1;
			}

			/*rotate*/
			if (button == GLUT_RIGHT_BUTTON)
			{
				mouse_mode = 2;
			}
		}
		else if (state == GLUT_UP) {

			mouse_mode = -1;
		}
	}
}

/******************************************************************************
Callback function for mouse wheel scroll
******************************************************************************/

void mousewheel(int wheel, int direction, int x, int y) {
	if (direction == 1) {
		zoom *= zoomspeed;
		glutPostRedisplay();
	}
	else if (direction == -1) {
		zoom /= zoomspeed;
		glutPostRedisplay();
	}
}

/******************************************************************************
Callback function for window reshaping
******************************************************************************/

void reshape(int width, int height) {
	win_width = width;
	win_height = height;

	aspectRatio = (float)width / (float)height;

	glViewport(0, 0, width, height);

	set_view(GL_RENDER);
}


/******************************************************************************
Callback function for scene display
******************************************************************************/

void display(void) {
	glClearColor(1.0, 1.0, 1.0, 1.0);  // background for rendering color coding and lighting

	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	set_view(GL_RENDER);
	set_scene(GL_RENDER, mesh);

	/*display the mesh*/
	display_mesh(mesh);

	glFlush();
	glutSwapBuffers();
	glFinish();

	CHECK_GL_ERROR();
}

/******************************************************************************
Diaplay the polygon with visualization results
******************************************************************************/

void display_mesh(Mesh* mesh) {

	glEnable(GL_POLYGON_OFFSET_FILL);
	glPolygonOffset(1., 1.);

	glEnable(GL_DEPTH_TEST);
	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
	glShadeModel(GL_SMOOTH);

	glEnable(GL_LIGHTING);
	glEnable(GL_LIGHT0);
	glEnable(GL_LIGHT1);

	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
	GLfloat mat_diffuse[4] = { 0.24f, 0.4f, 0.47f, 0.0f };
	GLfloat mat_specular[] = { 1.0f, 1.0f, 1.0f, 1.0f };
	glMaterialfv(GL_FRONT, GL_DIFFUSE, mat_diffuse);
	glMaterialfv(GL_FRONT, GL_SPECULAR, mat_specular);
	glMaterialf(GL_FRONT, GL_SHININESS, 50.0f);

	glBegin(GL_TRIANGLES);
	glColor3f(0.5f, 0.5f, 0.5f);
	for (auto f : mesh->flist) {
		glNormal3d(f->normal.x, f->normal.y, f->normal.z);
		for (int i = 0; i < 3; i++) {
			glVertex3d(f->verts[i]->x, f->verts[i]->y, f->verts[i]->z);
		}

	}

	glEnd();

	GLfloat mat_diffuse_ellip[4] = { 0.f, 0.4f, 0.f, 0.0f };
	glMaterialfv(GL_FRONT, GL_DIFFUSE, mat_diffuse_ellip);

	if (view_error) {
		for (auto e : errorEllip) {
			display_ellipsoid(e);
		}
	}

}

void generate_ellipsoids(Mesh* mesh) {

	errorEllip.clear();

	for (auto v : mesh->vlist) {
		if (!v->faces.empty()) {
			Ellipsoid* e = new Ellipsoid();
			e->x = v->x;
			e->y = v->y;
			e->z = v->z;
			ellipsoid_transformation(v,e);
			errorEllip.push_back(e);
		}
		
	}

}


void ellipsoid_transformation(Vertex* vert, Ellipsoid* ellip) {
		
	glm::dmat4x4 dmat = vert->Q;
	Matrix4d eigMat;

	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			eigMat(i,j) = dmat[i][j];
		}
	}

	EigenSolver<Matrix4d> es(eigMat);

	Matrix4cd D = es.eigenvalues().asDiagonal();
	Matrix4cd R = es.eigenvectors();
	Matrix4cd Rinv = es.eigenvectors().inverse();

	Matrix4cd RinvD = Rinv*D;
	glm::dmat4x4 totalT = glm::dmat4x4(1.0);
	glm::dmat4x4 rotT = glm::dmat4x4(1.0);
	glm::dmat4x4 scaleT = glm::dmat4x4(1.0);
	glm::dmat4x4 transT = glm::dmat4x4(1.0);
	transT = glm::translate(glm::mat4x4(1.), glm::vec3(vert->x, vert->y, vert->z));
	scaleT = glm::scale(glm::mat4x4(1.), glm::vec3(D(0,0).real(), D(1, 1).real(), D(2, 2).real()));

	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			//totalT[i][j] = RinvD(i, j).real();
			rotT[i][j] = Rinv(i, j).real();
			//scaleT[i][j] = D(j, i).real();
		}
	}
	ellip->xS = D(0, 0).real();
	ellip->yS = D(1, 1).real();
	ellip->zS = D(2, 2).real();
	//ellip->wS = D(3, 3).real();
	
	ellip->transT = transT;
	ellip->rotT = glm::transpose(rotT);
	ellip->scaleT = scaleT*.001;//glm::transpose(scaleT*.001);
	ellip->totalT = scaleT * transT;//* glm::transpose(rotT);
}

void display_ellipsoid(Ellipsoid* ellip) {
	glPushMatrix();

//		glMultMatrixd(glm::value_ptr(ellip->rotT));
		glTranslated(ellip->x, ellip->y, ellip->z);
		glScaled(ellip->xS/500., ellip->yS / 100., ellip->zS / 100.);
//		glMultMatrixd(glm::value_ptr(ellip->totalT));
//		glTranslated(ellip->x, ellip->y, ellip->z);
//		glScaled(ellip->xS / 1000., ellip->yS / 100., ellip->zS / 100.);
//		glMultMatrixd(glm::value_ptr(ellip->transT * ellip->scaleT* ellip->rotT));
//		glMultMatrixd(glm::value_ptr(ellip->transformation));
//		glMultMatrixd(glm::value_ptr(ellip->scaleT * ellip->rotT));
//		glMultMatrixd(glm::value_ptr(ellip->scaleT));
//		glMultMatrixd(glm::value_ptr(ellip->totalT));
			
		glutSolidSphere(1.0, 20, 20);
	glPopMatrix();
	
}
