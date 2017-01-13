#include <GL/glut.h>
#include "assignment3.h"
#include "init.h"
#include <iostream>
#include <math.h>

using namespace std;

vector<vector<Pt>> drawingVector;
vector<Matrix> transformMatrices;
vector<Pt> condensationPoly;

float pi = atan(1.0) * 4.0;

Matrix translate ( Vec v )
{
	Matrix rvalue;

	rvalue.data[0][0] = 1;
	rvalue.data[1][1] = 1;
	rvalue.data[2][2] = 1;
	rvalue.data[0][2] = v.x;
	rvalue.data[1][2] = v.y;

	return rvalue;
}

Matrix rotate ( Pt o, float theta )
{
	Matrix rvalue;

	rvalue.data[0][0] = cos(theta);
	rvalue.data[0][1] = -1.0*sin(theta);
	rvalue.data[0][2] = o.x + o.y*sin(theta) - o.x*cos(theta);
	rvalue.data[1][0] = sin(theta);
	rvalue.data[1][1] = cos(theta);
	rvalue.data[1][2] = o.y - o.y*cos(theta) - o.x*sin(theta);
	rvalue.data[2][0] = 0;
	rvalue.data[2][1] = 0;
	rvalue.data[2][2] = 1;

	return rvalue;
}

Matrix scale ( Pt o, float alpha )
{
	Matrix rvalue;

	rvalue.data[0][0] = alpha;
	rvalue.data[1][1] = alpha;
	rvalue.data[2][2] = 1;

	rvalue.data[0][2] = (1-alpha)*o.x;
	rvalue.data[1][2] = (1-alpha)*o.y;

	return rvalue;
}

Matrix nscale ( Pt o, Vec v, float alpha )
{
	Matrix rvalue;

	rvalue.data[0][0] = 1 + (alpha - 1)*v.x*v.x;
	rvalue.data[0][1] = (alpha - 1)*v.x*v.y;
	rvalue.data[0][2] = v.x*(o.x*v.x + o.y*v.y)*(1-alpha);
	rvalue.data[1][0] = (alpha - 1)*v.x*v.y;
	rvalue.data[1][1] = 1 + (alpha - 1)*v.y*v.y;
	rvalue.data[1][2] = v.y*(o.x*v.x + o.y*v.y)*(1 - alpha);
	rvalue.data[2][0] = 0;
	rvalue.data[2][1] = 0;
	rvalue.data[2][2] = 1;

	return rvalue;
}

Matrix image ( Pt p1, Pt p2, Pt p3, Pt q1, Pt q2, Pt q3 )
{
	Matrix rvalue;

	Matrix m1, m2;

	m1.data[0][0] = q1.x;
	m1.data[0][1] = q2.x;
	m1.data[0][2] = q3.x;
	m1.data[1][0] = q1.y;
	m1.data[1][1] = q2.y;
	m1.data[1][2] = q3.y;
	m1.data[2][0] = 1;
	m1.data[2][1] = 1;
	m1.data[2][2] = 1;

	m2.data[0][0] = p1.x;
	m2.data[0][1] = p2.x;
	m2.data[0][2] = p3.x;
	m2.data[1][0] = p1.y;
	m2.data[1][1] = p2.y;
	m2.data[1][2] = p3.y;
	m2.data[2][0] = 1;
	m2.data[2][1] = 1;
	m2.data[2][2] = 1;

	m2 = inverse(m2);

	rvalue = compose(m1,m2);

	return rvalue;
}

Matrix inverse(Matrix m1) {

	Matrix out;

	float det = m1.data[0][0] * m1.data[1][1] * m1.data[2][2]
		+ m1.data[1][0] * m1.data[2][1] * m1.data[0][2]
		+ m1.data[2][0] * m1.data[0][1] * m1.data[1][2]
		- m1.data[0][0] * m1.data[2][1] * m1.data[1][2]
		- m1.data[2][0] * m1.data[1][1] * m1.data[0][2]
		- m1.data[1][0] * m1.data[0][1] * m1.data[2][2];
	
	det = 1.0 / det;

	out.data[0][0] = det*(m1.data[1][1] * m1.data[2][2] - m1.data[1][2] * m1.data[2][1]);
	out.data[0][1] = det*(m1.data[0][2] * m1.data[2][1] - m1.data[0][1] * m1.data[2][2]);
	out.data[0][2] = det*(m1.data[0][1] * m1.data[1][2] - m1.data[0][2] * m1.data[1][1]);
	out.data[1][0] = det*(m1.data[1][2] * m1.data[2][0] - m1.data[1][0] * m1.data[2][2]);
	out.data[1][1] = det*(m1.data[0][0] * m1.data[2][2] - m1.data[0][2] * m1.data[2][0]);
	out.data[1][2] = det*(m1.data[0][2] * m1.data[1][0] - m1.data[0][0] * m1.data[1][2]);
	out.data[2][0] = det*(m1.data[1][0] * m1.data[2][1] - m1.data[1][1] * m1.data[2][0]);
	out.data[2][1] = det*(m1.data[0][1] * m1.data[2][0] - m1.data[0][0] * m1.data[2][1]);
	out.data[2][2] = det*(m1.data[0][0] * m1.data[1][1] - m1.data[0][1] * m1.data[1][0]);
	
	return out;
}

Matrix compose ( Matrix m1, Matrix m2 )
{
	Matrix rvalue;

	for (int x = 0; x < 3; x++) {
		for (int y = 0; y < 3; y++) {
			rvalue.data[x][y] = m1.data[x][0] * m2.data[0][y] 
							  + m1.data[x][1] * m2.data[1][y] 
							  + m1.data[x][2] * m2.data[2][y];
		}
	}

	return rvalue;
}

void setCondensationSet ( vector<Pt> pts )
{
	condensationPoly = pts;
}

void setIATTransformations ( vector<Matrix> transformations )
{
	transformMatrices = transformations;
}

// Draws the current IAT
void display(void){

	//SOME DEBUGGING CODE

	cout << "OUR MAGIC FRACTAL ELVES ARE PREPARING YOUR FRACTAL, PLEASE STAND BY......" << endl;

	/*	cout << "Drawing New shape:        " << endl;
	cout << "----------------------" << endl;

	for (Matrix x : transformMatrices) {

		cout << x.data[0][0] << " " << x.data[0][1] << " " << x.data[0][2] << endl;
		cout << x.data[1][0] << " " << x.data[1][1] << " " << x.data[1][2] << endl;
		cout << x.data[2][0] << " " << x.data[2][1] << " " << x.data[2][2] << endl;
		cout << "----------------------" << endl;

	}
*/
	//DETERMINE AMOUNT OF ITERATIONS

	int k = 3;

	switch (transformMatrices.size()) {
		case 1: k = 50; break;
		case 2: k = 10; break;
		case 3: k = 8; break;
		case 4: k = 7; break;
		case 5: k = 6; break;
		case 6: k = 5; break;
		case 7: k = 4; break;
	}

	//POPULATE CURRENT ITERATIONS

	drawingVector.clear();

	Pt origin;
	vector<Pt> addPolyLine;
	addPolyLine.push_back(origin);
	drawingVector.push_back(addPolyLine);

	//START THE ACTUAL LOOP

	vector<vector<Pt>> currentIteration;
	currentIteration = drawingVector;

	for (int i = 0; i < k; i++) {
		if (condensationPoly.size() > 0) {
			currentIteration.push_back(condensationPoly);
		}
		for (int x = 0; x < currentIteration.size(); x++) {
			for (int m = 0; m < transformMatrices.size(); m++) {
				addPolyLine.clear();
				for (int y = 0; y < currentIteration[x].size(); y++) {
					Pt add = currentIteration[x][y];
					float xtemp;
					float ytemp;

					xtemp = transformMatrices[m].data[0][0] * add.x + transformMatrices[m].data[0][1] * add.y + transformMatrices[m].data[0][2];
					ytemp = transformMatrices[m].data[1][0] * add.x + transformMatrices[m].data[1][1] * add.y + transformMatrices[m].data[1][2];

					add.x = xtemp;
					add.y = ytemp;

					addPolyLine.push_back(add);
				}
				drawingVector.push_back(addPolyLine);
			}
		}
		currentIteration = drawingVector;
		drawingVector.clear();
	}
	drawingVector = currentIteration;
	drawingVector.push_back(condensationPoly);
	//DRAW THE FRACTAL

	for (vector<Pt> vec : drawingVector) {
		if (vec.size() == 1) {
			glBegin(GL_POINTS);
			glVertex2f(vec[0].x, vec[0].y);
			glEnd();
		}
		if (vec.size() == 2) {
			glBegin(GL_LINES);
			glVertex2f(vec[0].x, vec[0].y);
			glVertex2f(vec[1].x, vec[1].y);
			glEnd();
		}
		if (vec.size() > 2) {
			glBegin(GL_LINE_LOOP);
			for (Pt i : vec) {
				glVertex2f(i.x, i.y);
			}
			glEnd();
		}
	}

	glFlush ( );
}

/* do not modify the reshape function */
void reshape ( int width, int height )
{
	glViewport ( 0, 0, width, height );
	glMatrixMode ( GL_PROJECTION );
	glLoadIdentity ( );    
	gluOrtho2D (-1, 1, -1, 1);
	glMatrixMode ( GL_MODELVIEW );
    glLoadIdentity ( );
}

void keyboard(unsigned char key, int x, int y)
{
	glClear(GL_COLOR_BUFFER_BIT);

	if (key == '1') {
		vector<Matrix> iat;

		glColor4f(1, 0, 1, 1);

		iat.push_back(image(
			Pt(-1, 1),
			Pt(-1,-1),
			Pt( 1,-1),
			Pt( 0, 0),
			Pt( 0,-1),
			Pt( 1,-1)
		));
		iat.push_back(image(
			Pt(-1, 1),
			Pt(-1,-1),
			Pt( 1,-1),
			Pt( 0, 1),
			Pt(-1, 1),
			Pt(-1, 0)
		));
		iat.push_back(image(
			Pt(-1, 1),
			Pt(-1,-1),
			Pt( 1,-1),
			Pt(-1, 0),
			Pt(-1,-1),
			Pt( 0,-1)
		));


		setIATTransformations(iat);

		vector<Pt> pts;
		// no condensation set
		//pts.push_back(Pt(0.25,0.25));
		//pts.push_back(Pt(0.25,0.75));
		//pts.push_back(Pt(0.75,0.75));
		//pts.push_back(Pt(0.75,0.25));

		setCondensationSet(pts);
	}
	if (key == '2') {

		glColor4f(0, 1, 0, 1);

		vector<Matrix> iat;

		iat.push_back(scale(Pt(1, 1), 0.5));
		iat.push_back(scale(Pt(-1, -1), 0.5));

		setIATTransformations(iat);

		vector<Pt> pts;
		// no condensation set
		pts.push_back(Pt(1,-1));
		pts.push_back(Pt(0,-1));
		pts.push_back(Pt(0, 0));
		pts.push_back(Pt(1, 0));

		setCondensationSet(pts);
	}
	if (key == '3') {

		vector<Matrix> iat;

		glColor4f(1, 0, 0, 1);

		Matrix rot = rotate(Pt(0,0),3.14);
		Matrix sca = scale(Pt(0, 0), 0.38);

		iat.push_back(compose(sca,rot));
		iat.push_back(scale(Pt(       0,     .8), 0.38));
		iat.push_back(scale(Pt(	.760845, .24721), 0.38));
		iat.push_back(scale(Pt(-.760845, .24721), 0.38));
		iat.push_back(scale(Pt( .470230,-.64721), 0.38));
		iat.push_back(scale(Pt(-.470230,-.64721), 0.38));
		setIATTransformations(iat);

		vector<Pt> pts;
		// no condensation set
		//pts.push_back(Pt(1, -1));
		//pts.push_back(Pt(0, -1));
		//pts.push_back(Pt(0, 0));
		//pts.push_back(Pt(1, 0));

		setCondensationSet(pts);
	}
	if (key == '4') {

		vector<Matrix> iat;

		glColor4f(0, 0, 1, 1);

		iat.push_back(scale(Pt(-.5, .866), 0.3333));
		iat.push_back(scale(Pt(.5, .866), 0.3333));
		iat.push_back(scale(Pt(-.5, -.866), 0.3333));
		iat.push_back(scale(Pt(.5, -.866), 0.3333));
		iat.push_back(scale(Pt(1, 0), 0.3333));
		iat.push_back(scale(Pt(-1, 0), 0.3333));
		setIATTransformations(iat);

		vector<Pt> pts;
		// no condensation set
		setCondensationSet(pts);
	}
	if (key == '5') {

		vector<Matrix> iat;

		glColor4f(0.3, 1, 0.3, 1);

		iat.push_back(compose(
			scale(Pt(1, 0), 0.58),
			rotate(Pt(0, 0), pi *7.0 / 6.0)
		));
		iat.push_back(compose(
			scale(Pt(-1, -0), 0.58),
			rotate(Pt(0, 0), -pi*7.0 / 6.0)
		));

		setIATTransformations(iat);

		vector<Pt> pts;
		// no condensation set
		pts.push_back(Pt(0, .1));
		pts.push_back(Pt(.0866, .05));
		pts.push_back(Pt(.0866, -.05));
		pts.push_back(Pt(0, -.1));
		pts.push_back(Pt(-.0866, -.05));
		pts.push_back(Pt(-.0866, .05));

		setCondensationSet(pts);
	}
	if (key == '6') {

		vector<Matrix> iat;

		glColor4f(0, 1, 1, 1);

		iat.push_back(image(
			Pt(-1, 1),
			Pt(-1, -1),
			Pt(1, -1),
			Pt(1, 0),
			Pt(1, -1),
			Pt(0, -1)
		));
		iat.push_back(image(
			Pt(-1, 1),
			Pt(-1, -1),
			Pt(1, -1),
			Pt(-1, 0),
			Pt(-1, -1),
			Pt(0, -1)
		));
		iat.push_back(image(
			Pt(-1, 1),
			Pt(-1, -1),
			Pt(1, -1),
			Pt(0, 1),
			Pt(-1, 1),
			Pt(-1, 0)
		));

		setIATTransformations(iat);

		vector<Pt> pts;
		// no condensation set
		pts.push_back(Pt(0, .25));
		pts.push_back(Pt(0.25, 0));
		pts.push_back(Pt(0, -0.25));
		pts.push_back(Pt(-0.25, 0));

		setCondensationSet(pts);
	}
	if (key == '7') {

		vector<Matrix> iat;

		glColor4f(1, 1, 0, 1);

		iat.push_back(compose(scale(Pt(-0.9, -0.9), 0.44), rotate(Pt(0, 0), pi / 6.0)));
		iat.push_back(compose(scale(Pt( 0.9, -0.9), 0.44), rotate(Pt(0, 0), -pi / 6.0)));
		iat.push_back(compose(scale(Pt(-0.9,  0.9), 0.44), rotate(Pt(0, 0), -pi*7 / 6.0 )));
		iat.push_back(compose(scale(Pt( 0.9,  0.9), 0.44), rotate(Pt(0, 0), pi*7 / 6.0)));
		//iat.push_back(scale(Pt(0, 0), 0.45));


		setIATTransformations(iat);

		vector<Pt> pts;
		// no condensation set
		//pts.push_back(Pt(0, -0.2));
		//pts.push_back(Pt(-0.2,  0));
		//pts.push_back(Pt(0,  0.2));
		//pts.push_back(Pt(0.2, 0));

		setCondensationSet(pts);
	}

	glutPostRedisplay();
}


int main ( int argc, char** argv )
{
	glutInit ( &argc, argv );

	glutInitDisplayMode ( GLUT_SINGLE | GLUT_RGB );
	glutInitWindowSize ( 500, 500 );
	glutInitWindowPosition ( 100, 100 );
	glutCreateWindow ( "Tunc Gocay - Homework 3" );
	glClear(GL_COLOR_BUFFER_BIT); 
	glColor4f(1, 1, 1, 1);
	init ( );	

	glutDisplayFunc ( display );
	glutReshapeFunc ( reshape );
	glutKeyboardFunc( keyboard);
	glutMainLoop ( );


	return 0;
}
