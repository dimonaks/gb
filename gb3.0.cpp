//Not adding impurities at all


#include <iostream>//
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
#define type double
#define db double
using namespace std;

typedef complex<double> compln;

struct limits{
db xmax; db xmin; db ymax,ymin,zmax,zmin;
};
struct vector_dec{
double x;double y;double z; double length();
};
double vector_dec::length(){return(sqrt(x*x+y*y+z*z));}
vector_dec calculate_xred(double L[3][3],double x,double y, double z);

void distance(bool *at_to_delete, int n_at, db *x_cell, db *y_cell, db *z_cell,db *acell);
void write_abinit(char *name, bool *at_to_delete, int n_at, \
		db *x_cell, db *y_cell, db *z_cell,db *acell, int n_d_at,int ntypat);
void cut_cell(int &n_k_nodes,db*,db*,db*,db*,db*,db*,limits lim,int s);
void write_xyz(char *name,int n_at,int *i_znucl_at,db *x,db *y, db *z);



int h[100000],k[100000],l[100000];
//This source file was used to create 54 ti atom samples with grain boundary, жестко забито, что всего 1 атом примеси в конце
int main(){
double a[3][3];
int nbasis,i,j,ki,s;;
double basis_1[100],basis_2[100],basis_3[100],basis_dec[3];
compln im_unit(0,1);

int node1,node2,node0;

db acell_i[3];
ifstream in; // input
in.open("real.in", ios::in);
in >> acell_i[0]>>acell_i[1]>>acell_i[2];
for(int i=0;i<=2; i++){
in >> a[i][0] >> a[i][1] >> a[i][2] ;
cout << a[i][0]<<" " << a[i][1]<<" " << a[i][2] <<"-rprim \n";}
for(int j=0;j<=2; j++){
a[0][j]=a[0][j]*acell_i[0];
a[1][j]=a[1][j]*acell_i[1];
a[2][j]=a[2][j]*acell_i[2];
}
for(int i=0;i<=2; i++){
cout << a[i][0]<<" " << a[i][1]<<" " << a[i][2] <<"-rpimd \n";}

in >> nbasis;
cout << nbasis<<"\n";
int typat[nbasis];

for(i=0;i<nbasis; i++)
{
in >> basis_1[i]>>basis_2[i]>>basis_3[i];
cout << basis_1[i]<<" " << basis_2[i]<<" " << basis_3[i] <<" \n";
}
for(i=0;i<nbasis; i++)
{
in >> typat[i];
cout << "typat ="<<typat[i]<< " ";
}
cout << "\n";
in.close();

int ni=14,nj=14,nk=4,N;//Dimensions of the system, можно выбирать с запасом, чтобы потом вырезать нужный кусок после поворота
N=ni*nj*nk*8*nbasis;
double x[N],y[N],z[N],d[N];
int i_znucl_at[N];
int fordopes[100];//массив для номеров атомов, возле которых будут вставляться примеси


//cout <<"N =" <<N;
s=0;
int i_at=0;
//Используя заданные максимальные размеры системы (сделать: Вынести размеры системы во входной файл конфигурации) пройти
//по всем атомам,определить номера узловых точек ячейки и рассчитать декартовы координаты всех атомов.
for(i=-ni;i<ni;i++)
	for(j=-nj;j<nj;j++)
		for(ki=-nk;ki<nk;ki++)
{
//cout << "s= "<<s << "\n";
h[s]=i;
k[s]=j;
l[s]=ki;
int U[3];
U[0]=i;U[1]=j;U[2]=ki;
double sum_bi[3]={0,0,0},r_length=0;
for(int i=0;i<3;i++)
for(int k=0;k<3;k++)
sum_bi[i]+=U[k]*a[k][i];// координаты векторов
for(int i=0;i<3;i++)
r_length+=sum_bi[i]*sum_bi[i];
r_length=sqrt(r_length); //длина векторов для теста
cout <<"r_length= "<<r_length<<"\n";
printf( "%i %i %i \n",h[s],k[s],l[s] );
//Выбор угловых узлов решетки, которые будут определять размеры ячейки
//в соответствии рисунком из работы [] и моим выбором осей, угловыми узлами будут [320] и [-140]
//Запомнить номера узлов.
if(U[0]==0&&U[1]==0&&U[2]==0) node0=s;
if(U[0]==3&&U[1]==2&&U[2]==0) node1=s;
if(U[0]==-1&&U[1]==4&&U[2]==0) node2=s;


//Расчитать декартовы координаты всех атомов
for(int i=0;i<nbasis;i++){
for(int j=0;j<3;j++)
{
//cout<<"basis_1[i] "<<basis_1[i]<<" basis_2[i]"<<basis_2[i]<<"basis_3[i]"<<basis_3[i]<<"\n";
basis_dec[j]=basis_1[i]*a[0][j]+basis_2[i]*a[1][j]+basis_3[i]*a[2][j];//Вычисление декартовых координат атомов базиса
//cout<<"basis_dec[j]"<<basis_dec[j]<<"\n";
}
x[s+i]=h[s]*a[0][0]+k[s]*a[1][0]+l[s]*a[2][0]+basis_dec[0];
y[s+i]=h[s]*a[0][1]+k[s]*a[1][1]+l[s]*a[2][1]+basis_dec[1];
z[s+i]=h[s]*a[0][2]+k[s]*a[1][2]+l[s]*a[2][2]+basis_dec[2];
d[s+i]=sqrt(x[s+i]*x[s+i]+y[s+i]*y[s+i]+z[s+i]*z[s+i]);//расстояние из точки 0 до каждого узла (по сути атома)
i_at++;
}
s+=nbasis;
}


int n_at;
int ntypat=1;
//Поворот ячейки, чтобы выбрать нужную ориентацию ячейки относительно осей координат
double fi,tan_fi,fi2,x_new[N],y_new[N],z_new[N];
tan_fi=y[node1]/x[node1];
fi=atan(tan_fi);
fi2=atan(x[node2]/y[node2]);//угол  между вектором, направленным к узлу node2 и осью y
//повернуть против часовой стрелки всю ячейку так, чтобы node2 начал лежать на оси y
for(int i=0;i<s;i++){
x_new[i]=x[i]*cos(fi2)-y[i]*sin(fi2);
y_new[i]=x[i]*sin(fi2)+y[i]*cos(fi2);
z_new[i]=z[i];
}
cout<<"fi= "<<fi/M_PI*180.*2<<" degrees"<<"\n";
cout<<"fi2= "<<fi2/M_PI*180.*2<<" degrees"<<"\n";





//Вырезаем нужную ячйку
double x_cell[N*2],y_cell[N*2],z_cell[N*2];
int n_k_nodes;
//задаем границы ячейки
limits lim;
lim.ymin=0;
lim.ymax=d[node2]*2;
lim.xmin=0;
lim.xmax=d[node1];
lim.zmin=0;
lim.zmax=a[2][2];

cut_cell(n_k_nodes,x_cell,y_cell,z_cell,x_new,y_new,z_new,lim,s);

//Задаем acell для абинита
double acell[3];
acell[0]=d[node1];
acell[1]=lim.ymax*2.;
acell[2]=lim.zmax;

bool at_to_delete[N*2]; memset(at_to_delete, 0, N*2);
char outname5[80];


//Строим ячейку с границей путем зеркального отражения
for(int i=0;i<n_k_nodes;i++){
y_cell[i+n_k_nodes]=acell[1]-y_cell[i];
x_cell[i+n_k_nodes]=x_cell[i];
z_cell[i+n_k_nodes]=z_cell[i];}
n_at=n_k_nodes*2;

//Определим расстояния между атомами и запишем в массив at_to_delete номера слишком близко лежащих атомов.
distance(at_to_delete, n_at,x_cell, y_cell, z_cell,acell);

int n_at_end = 0, n_d_at = 0,iend = 0;

strcat(outname5,"0.yesgb.geo");

//Пока вручную указываем ядерный заряд атомов
//Также удаляем (пропускаем при записи в конечный массив) слишком близко расположенные атомы
db x_cell_end[N],y_cell_end[N],z_cell_end[N];
for (int i=0; i<n_at; i++){
	if(at_to_delete[i]) continue;
	i_znucl_at[iend] = 22;
	x_cell_end[iend] = x_cell[i];
	y_cell_end[iend] = y_cell[i];
	z_cell_end[iend] = z_cell[i];
	iend++;
}
n_at_end = iend;
iend = 0;
//вывод координат в формате .xyz
write_xyz(outname5,n_at_end,i_znucl_at,x_cell_end,y_cell_end,z_cell_end);
strcat(outname5,".in");

write_abinit(outname5, at_to_delete, n_at, \
		x_cell, y_cell, z_cell,acell, n_d_at,ntypat);
cout <<"\n\n\n\n";












//строим ячейку без границы обычной трансляцией!
cut_cell(n_k_nodes,x_cell,y_cell,z_cell,x_new,y_new,z_new,lim,s);
ntypat=1; n_d_at=0;
for(int i=0;i<n_k_nodes;i++){
	x_cell[i+n_k_nodes]=x_cell[i];
	y_cell[i+n_k_nodes]=acell[1]*0.5+y_cell[i];
	z_cell[i+n_k_nodes]=z_cell[i];
//	cout<<"i="<<i<<" coord nogb:"<<x_cell[i]<<" "<<y_cell[i]<<" "<<z_cell[i]<<endl;
	}
n_at=n_k_nodes*2;



distance(at_to_delete, n_at,x_cell, y_cell, z_cell,acell);
char outname1[80];
strcpy(outname1,"1nogb.geo.in");
write_abinit(outname1, at_to_delete, n_at, \
		x_cell, y_cell, z_cell,acell, n_d_at,ntypat);

//write_xyz в отличии от write_abinit не может удалить атомы, содержащиеся в массиве at_to_delete,
// поэтому делаем это здесь
n_at_end = 0; iend = 0;

for (int i=0; i<n_at; i++){
	if(at_to_delete[i]) continue;
	i_znucl_at[iend] = 22;
	x_cell_end[iend] = x_cell[i];
	y_cell_end[iend] = y_cell[i];
	z_cell_end[iend] = z_cell[i];
	iend++;
}
n_at_end = iend;
//вывод координат в формате .xyz
write_xyz(outname1,n_at_end,i_znucl_at,x_cell_end,y_cell_end,z_cell_end);




//Ячейка без границы вырезанная сразу нужных размеров для теста
lim.ymax=d[node2]*2;
cout << lim.ymax<<endl;
cut_cell(n_k_nodes,x_cell,y_cell,z_cell,x_new,y_new,z_new,lim,s);
distance(at_to_delete, n_k_nodes,x_cell, y_cell, z_cell,acell);
n_at=n_k_nodes;
char outname2[80];
strcpy(outname2,"test1nogb.geo.in");
write_abinit(outname2, at_to_delete, n_k_nodes, \
		x_cell, y_cell, z_cell,acell, n_d_at,ntypat);












}//end main

void write_abinit(char *name, bool *at_to_delete, int n_at, \
		db *x_cell, db *y_cell, db *z_cell,db *acell, int n_d_at,int ntypat){
	//Выводим файл геометрии для abinit
	ofstream out;
	out.open(name, ios::out);
	out.precision(12);
	int natom=0;
	out<<"acell ";
	for(int i=0;i<3;i++)
	out<<acell[i]<<" ";
	out<<"\n";
	out<<"rprim 1 0 0 0 1 0 0 0 1\n";
    db L[3][3]={{acell[0],0,0},{0,acell[1],0},{0,0,acell[2]}};



    vector_dec xred;
//	out<<"xcart ";
	for(int i=0;i<n_at;i++){
	if(!at_to_delete[i]){//Исключаем вывод ранее отмеченных к удалению атомов
	natom++;

//	out << x_cell[i]<<" "<<y_cell[i]<<" "<<z_cell[i]<<" "<<"\n";
	}
	}

//Пишем xred
    out << "xred ";
	for(int i=0;i<n_at;i++){
	if(!at_to_delete[i]){//Исключаем вывод ранее отмеченных к удалению атомов
    xred=calculate_xred(L,x_cell[i],y_cell[i],z_cell[i]);
	out << xred.x<<" "<<xred.y<<" "<<xred.z<<" "<<"\n";
	}
	}






	out<<"ntypat "<<ntypat<<endl;
	out<<"typat "<<natom-n_d_at<<"*1 ";	if(n_d_at>0&&ntypat>1) out<<n_d_at<<"*2"<<endl;
	out<<endl;
	out<<"natom "<<natom<<endl;
	out.close();

}





void distance(bool *at_to_delete, int n_at, db *x_cell, db *y_cell, db *z_cell,db *acell){
bool checked[n_at];
memset(checked, 0, n_at);
memset(at_to_delete, 0, n_at);
db sizex=acell[0];
db sizey=acell[1];
db sizez=acell[2];
db sizex05=sizex*0.5;
db sizey05=sizey*0.5;
db sizez05=sizez*0.5;
for (int i=0;i<n_at;i++){
	 checked[i]=1;
	for (int j=0;j<n_at;j++){
if(!checked[j]){

 db xx=x_cell[i]-x_cell[j];
 db yy=y_cell[i]-y_cell[j];
 db zz=z_cell[i]-z_cell[j];

    if(xx > sizex05) xx=xx-sizex;
	if(xx < -sizex05) xx=xx+sizex;
    if(yy > sizey05) yy=yy-sizey;
	if(yy < -sizey05) yy=yy+sizey;
    if(zz > sizez05) zz=zz-sizez;
	if(zz < -sizez05) zz=zz+sizez;
db	rr=xx*xx+yy*yy+zz*zz;
db	d=sqrt(rr);
db r=3;
if(d<r){ //если расстояние между атомами меньше 2 Бор

	//if(abs(x_cell[i]-x_cell[j])>r)x_cell[i]=(x_cell[i]+x_cell[j]-sizex)*0.5;
//if(abs(y_cell[i]-y_cell[j])>r)y_cell[i]=(y_cell[i]+y_cell[j]-sizey)*0.5;
//if(abs(z_cell[i]-z_cell[j])>r)z_cell[i]=(z_cell[i]+z_cell[j]-sizez)*0.5;
//x_cell[i]=(x_cell[i]+x_cell[j])/2.;
//y_cell[i]=(y_cell[i]+y_cell[j])/2.;
//z_cell[i]=(z_cell[i]+z_cell[j])/2.;//т.е. атом i теперь занимает центральное положение между бывшими положениями
	cout<<"i="<<i<<" coord was:"<<x_cell[i]<<" "<<y_cell[i]<<" "<<z_cell[i]<<endl;
	x_cell[i]=xx*0.5+x_cell[j];//Позволяет найти среднее положение, то что выше закоментено работает неправильно
	y_cell[i]=yy*0.5+y_cell[j];
	z_cell[i]=zz*0.5+z_cell[j];
	at_to_delete[j]=1;
	cout<<"i="<<i<<" coord become:"<<x_cell[i]<<" "<<y_cell[i]<<" "<<z_cell[i]<<endl;
cout<<"atom j="<<j<<" will be deleted, because it is to close to atom i="<<i<<" ,d="<<d<<endl;
cout<<"j  coord:"<<x_cell[j]<<" "<<y_cell[j]<<" "<<z_cell[j]<<endl;
cout<<"i coord:"<<x_cell[i]<<" "<<y_cell[i]<<" "<<z_cell[i]<<endl;
 	 	 }
				}
								}}

}








void cut_cell(int &n_k_nodes, db *x_cell,db *y_cell,db *z_cell,db *x_new,db *y_new,db *z_new,limits lim,int s){
int k_nodes=0;
for(int i=0;i<s;i++){
//Выбираем только те атомы, которые попадают в обозначенную ячейку. По оси z выбраются от 0 до вектора с-т.е всего два ГПУ слоя.
//Здесь по x и z не включаются узлы через период, однако добавлен узел по периоду y, так как он потом используется
//Верхние коментарии теперь мало имеют значения, так теперь с учетом переиодических условий
//совпадающие атомы удаляются функцией distance
if(y_new[i]>=lim.ymin&&y_new[i]<=lim.ymax\
		&&x_new[i]>=lim.xmin-1e-11&&x_new[i]<=lim.xmax+1e-11\
		&&z_new[i]>=lim.zmin&&z_new[i]<lim.zmax){
//cout<<"[uvw], s="<<i<<" k_nodes=" <<k_nodes<<"\n";
//cout<<h[i]<<" "<<k[i]<<" "<<l[i]<<"\n";//индексы hkl есть только для узлов (четные i!)
//cout<<x_new[i]<<" "<<y_new[i]<<"\n";//В данном случае для всех значений i есть координаты
x_cell[k_nodes]=x_new[i];
y_cell[k_nodes]=y_new[i];
z_cell[k_nodes]=z_new[i];
k_nodes++;
}
}
n_k_nodes=k_nodes;
//cout<<"x[node1]"<<x[node1]<<"\n";
//cout<<"y[node1]"<<y[node1]<<"\n";
//cout<<"x[node2]"<<x[node2]<<"\n";
//cout<<"y[node2]"<<y[node2]<<"\n";
//cout<<"x_new[node2]"<<x_new[node2]<<"\n";
//cout<<"y_new[node2]"<<y_new[node2]<<"\n";
//cout<<"y_new[0]= "<<y_new[0]<<"\n";
//cout<<"d[node2]"<<d[node2]<<"\n";
//cout<<"d[node1]"<<d[node1]<<"\n";
//cout<<"n_k_nodes= "<<n_k_nodes<<"\n";
}


void write_xyz(char *name,int n_at,int *i_znucl_at,db *x,db *y, db *z){
char name2[80];
strcpy(name2,name);
ofstream out_xyz; // открываем файл для записи
out_xyz.open(strcat(name2,".xyz"), ios::out);
out_xyz<<n_at<<endl;
out_xyz<<name<<endl;
for(int i=0;i<n_at;i++){
if (i_znucl_at[i] == 22) out_xyz<<"Ti ";
if (i_znucl_at[i] == 6) out_xyz<<"C ";
out_xyz<<x[i]<<" "<<y[i]<<" "<<z[i]<<endl;}
out_xyz.close();
}



vector_dec calculate_xred(double L[3][3],double x,double y, double z){

db a1=-L[0][1]/1./L[0][0];
if(L[0][1]==0)a1=0;
db c=-L[0][2]/1./L[0][1];
if(L[0][2]==0)c=0;
db a11=a1*L[1][0]+L[1][1];
db a12=a1*L[2][0]+L[2][1];
db c11=c*L[1][1]+L[1][2];
db c12=c*L[2][1]+L[2][2];
db A=-c11/a11;
if(c11==0)A=0;
db n12=a1*x+y;
db n23=c*y+z;
vector_dec xred;
 xred.z=(n12*A+n23)/(A*a12+c12);
 xred.y=(n12-xred.z*a12)/a11;
 xred.x=(x-xred.y*L[1][0]-xred.z*L[2][0])/L[0][0];
return xred;
}

