#include "gb.h"





vectormath operator+(vectormath a, vectormath b) {
	vectormath temp;
	temp.x = a.x + b.x;
	temp.y = a.y + b.y;
	temp.z = a.z + b.z;
	return temp;
}

vectormath operator-(vectormath a, vectormath b) {
	vectormath temp;
	temp.x = a.x - b.x;
	temp.y = a.y - b.y;
	temp.z = a.z - b.z;
	return temp;
}

vectormath operator*(vectormath a, double scalar) {
	vectormath temp;
	temp.x = a.x * scalar;
	temp.y = a.y * scalar;
	temp.z = a.z * scalar;
	return temp;
}

double operator*(vectormath a, vectormath b) {
	return a.x*b.x + a.y*b.y + a.z*b.z;
}


vectormath operator%(vectormath a, vectormath b) {
	vectormath d;
	d.x=(a.y*b.z - a.z*b.y);
	d.y=(a.z*b.x - a.x*b.z);
	d.z=(a.x*b.y - a.y*b.x);
	return d;
}

double abs(vectormath a) {
	return sqrt(a.x*a.x + a.y*a.y + a.z*a.z);
}

vectormath operator/(vectormath a, double V) {
	vectormath d;
	d.x=a.x/V;d.y=a.y/V;d.z=a.z/V;
	return d;
}

void vectormath::print() {
	printf("(%f,%f,%f) \n",x,y,z);
}

void vectormath::set(db r[3]) {
	x = r[0]; y = r[1]; z = r[2];
}