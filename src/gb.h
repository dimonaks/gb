//Alglib is turned off to save time
#include <iostream>//
#include <cstdlib>// div
#include <algorithm>
#include <fstream>
#include <ctime>
#define _USE_MATH_DEFINES
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

//#include <cmath>
#include <new>
#include <complex>
#include <cstring>
#include <vector>
//#include <solvers.h> // Alglib library
//#include <ap.h> // Alglib library
//#include <alglibinternal.h> // Alglib library
#define type double
#define db double
#define max_at 8000
#define n_sort 2
//#define to_angstrom 0.52917721092
//#define len_units " Bohr"
#define to_angstrom 1 //switch to internal angstroms units
#define len_units " Angstrom"
using namespace std;

typedef complex<double> compln;

struct limits{
    db xmax; db xmin; db ymax,ymin,zmax,zmin;
};
struct vector_dec{
    double x;double y;double z; double length();
};

class vectormath {
    public:
    double x,y,z;
    public:
    vectormath(double vX=0, double vY=0, double vZ=0) { x=vX; y=vY; z=vZ; } // конструктор
    // ~vector();
    friend vectormath operator+(vectormath a, vectormath b); //сложение векторов
    friend vectormath operator-(vectormath a, vectormath b); //вычитание векторов
    friend vectormath operator*(vectormath a, double scalar); //умножение вектора на скаляр
    friend double operator*(vectormath a, vectormath b); //скалярное произведение векторов
    friend vectormath operator%(vectormath a, vectormath b); //векторное произведение векторов
    friend vectormath operator/(vectormath a, double V);//деление компонент вектора на скаляр
    friend double abs(vectormath a); //модуль вектора
    void print();//вывод вектора в терминале
    void set(db r[3]);
};


class CrystalCell{
    db rprim[3][3];// Первый индекс - номер вектора. Второй индекс - компоненты. Наоборот по отношению к Abinit, так как в C сначала меняется второй индекс
    db acell[3];
    db znucl[100]; //Should be the same in fortran
    db xred[max_at][3],xcart[max_at][3];

    db firotatez;
    int typat[max_at];
    // sort[][0] means the grain; sort[][0] = 0 - grain A, 1 - grain B;
    // sort[][1] means the closed packed layer; sort[][1] = 0 - layer A, 1 - layer B, 2 - layer C 
    int sort[max_at][n_sort]; 
    db mul[3];
    int periodcell[6];

    db spacing; // spacing between hkl1 planes
    vectormath recip[3], normal_to_boundary; 
    db deldist, shiftcell[2], shiftgrain[2];
    db slicepms[3]; //parameters of slicing
    db shiftatoms[4]; //shifting atoms
    db hex_a,hex_c; //only for case of hex lattices, needed to remember initial lat. constants
    db gbpos; //position of first grain boundary plane along x
    string calctype;
public:
    string name;
    int natom;
    int version;
    db  a_c_conv[4], sixvalues[6]; //used in create_scaled_versions()
    vectormath rprimd[3];
    int hkl1[3], uvw1[3], uvw2[3], uvw3[3];//,uvtw1[3], uvtw2[3], uvtw3[3];
public:
    CrystalCell();
    void choose_new_rprimd(vectormath &r1, vectormath &r2, vectormath &r3);
    void readin(string inname);//read from file all characteristics of cell
    void periodicfill_and_cut(vectormath r1, vectormath r2, vectormath r3,db bob,db upb);//Заполняет пространство в соответствии с параметром periodcell
    void periodicfill(int nx, int ny, int nz);
    void rotate(vectormath &r1, vectormath &r2, vectormath &r3); //Поворачивает оси координат на угол firotatez относительно оси z
    void cut(limits lim);//Вырезает параллелипед в соответствии с заданными пределами lim
    void mirrory();//Отражает ячейку в плоскости перпендикулярной оси y.
    db make_boundary(db delta_v); //delta_v - specific volume added to each grain boundary
    int delnear();//Удаляет близколежащие атомы
    void delnear_by_xred();//Удаляет близколежащие атомы
    void write_xyz(string additional);
    void write_abinit(string additional,  const CrystalCell & init = CrystalCell() );   //const - быстрее работает, но нельзя изменять
    void outcell();
    void findpores(db ri);
    void cpxyz(int ,int ,int); //copy nx,ny,nz times by x, y,z
    void return_atoms_to_cell();
    void shift_atoms(); //control input parameter shiftatoms
    void shift_cell();
    void shift_grain();
    void scale_hcp_acell(db a,db c);
    void center();
    void cut_rectangular(db xl, db xr, db yl, db yr , db zl, db zr, string cuttyp);
};

extern string read_plus(string, string);



//Inline functions
inline db red_prec(db x, db prec = 1000.) {x = x * prec; x = round(x) / prec; return x; }
std::vector<double> linspace(double a, double b, int n);
//Constants:
extern int fnlen;     //! maximum length of file name variables
extern int maxstrlen; //! maximum length of input string
static db rprimd_hcp;
extern int test;

