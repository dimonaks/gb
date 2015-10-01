#include "gb.h"
vectormath calculate_xred(vectormath rprimd[3], double x,double y, double z);

void four_index_to_three(int uvtw[3], int *uvw);
void plane_to_direction(int *uvw, int hkl[3], vectormath rprimd[3]);
void combine_z_with(vectormath vector_to_combine_with, vectormath *arrv, int size_of_arrv);
void show_possible_r1(vectormath rprimd[3], vectormath r2, vectormath r3);
//void readplus(string);
extern "C" {void readinnatom_(char* filnamin,\
int *natom);};

extern "C" {void readin_(char* filnamin,\
db *acell,db rprim[3][3],db xred[max_at][3],int *typat, db *znucl,
int *periodcell, db *firotatez,\
int *hkl1, int *uvtw1, int *uvtw2, int *uvtw3, int *uvw1, int *uvw2, int *uvw3, db *mul, int *test, \
db *slicepms, db *shiftatoms, db *deldist, db *shiftcell, db *shiftgrain,  db *a_c_conv, db *sixvalues);};

vectormath vector_dec_for_UVW(int d[3], vectormath rprimd[3], db mul);



//Constants:
int fnlen = 264;     //! maximum length of file name variables
int maxstrlen = 20000; //! maximum length of input string
int test;
//Functions for classes:
double vector_dec::length(){return(sqrt(x*x+y*y+z*z));};




//Functions
vectormath vector_dec_for_UVW(int d[3], vectormath rprimd[3], db mul) {
    //Calculate decart coordinates of vector for direction d
    vectormath r;
    db U, V, W;

    U = d[0] * mul; V = d[1] * mul; W=d[2] * mul; 

    r.x = U*rprimd[0].x + V*rprimd[1].x + W*rprimd[2].x;
    r.y = U*rprimd[0].y + V*rprimd[1].y + W*rprimd[2].y;
    r.z = U*rprimd[0].z + V*rprimd[1].z + W*rprimd[2].z;

    return r;
}

void four_index_to_three(int uvtw[3], int *uvw) {
 
    if( (uvtw[0] + uvtw[1]) != -uvtw[2] ) {cerr <<"\nError: four_index_to_three, direction is not correct\n"; throw 20;}

    uvw[0] = 2 * uvtw[0] + uvtw[1];
    uvw[1] = 2 * uvtw[1] + uvtw[0];
    uvw[2] = uvtw[3];
}
void three2four(int *uvtw, int uvw[3]) {
 
    //db u, v, t, w;
    uvtw[0] = (2*uvw[0] - uvw[1]);
    uvtw[1] = (2*uvw[1] - uvw[0]);
    uvtw[2] = -(uvtw[0]+uvtw[1]);
    uvtw[3] =  3 * uvw[2];

    if( (uvtw[0] + uvtw[1]) != -uvtw[2] ) {cerr <<"\nthree2four Error, direction is not correct\n"; throw 20;}
    
    if(uvtw[0]%3 == 0 && uvtw[1]%3 == 0 && uvtw[3]%3 == 0 )
        for (int i=0;i<4;i++)
            uvtw[i] = uvtw[i]/3;
    

}

void plane_to_direction(int *uvw, int hkl[3], vectormath rprimd[3]) {
    //Empty function for future - to find perpendicular to plane through reciprocal space
    //and solving of linear system using matrix solvel in library alglib.

 //find description in alglib/manual.cpp.html#sub_rmatrixsolve
    //first find reciprocal rprimd[3]
//    alglib::rmatrixsolve(
//    real_2d_array a,
//    ae_int_t n,
//    real_1d_array b,
//    ae_int_t& info,
//    densesolverreport& rep,
//    real_1d_array& x);
    return;
}
void show_possible_r1(vectormath rprimd[3], vectormath r2, vectormath r3) {
    cout << "\nHere is the list of possible r1:\n";   
    int u, v, w, max = 13;
    vectormath r1;
    int mul = 1;
    int d[3];
    int uvtw[4], uvw[3];
    db alfa, fi, maxlength = 60,temp ;//40 len_units
    //cout.precision(20);
    for(u = -max;u < max;u++)
        for(v = -max;v < max;v++)
            for(w = -max;w < max;w++) {
                d[0] = u; 
                d[1] = v; 
                d[2] = w;
                r1 = vector_dec_for_UVW(d, rprimd, mul);
                if(abs(r1) > maxlength || abs(r1) == 0) continue;

                temp = r1 * r3 / abs(r1) / abs(r3);
                //reduce the precision of temp 
                temp = roundf(temp * 10000.0f) / 10000.0f; //because temp can be 1.000000002 which cause nan for acos

                alfa = acos(temp) * 180 / M_PI;
                if(alfa < 89 || alfa > 91) continue;
 
                temp = r1 * r2 / abs(r1) / abs(r2);
                temp = roundf(temp * 10000.0f) / 10000.0f;  
 
                fi = acos(temp) * 180 / M_PI;
                if(fi < 70 || fi > 110) continue;
                //temp = r1 * r2 / abs(r1) / abs(r2);
                //cout << "\nfi "<<fi<<" " << temp<< "\n";
                uvw[0] = u; uvw[1] = v; uvw[2] = w;
                three2four(uvtw, uvw);
                cout << "\nVector [" << u << v << w << "] has length = " << abs(r1) \
                <<"; in four index system is ["<<uvtw[0]<<" "<<uvtw[1]<<" "<<uvtw[2]<<" "<<uvtw[3]<<"]\n";
                cout << "\nAngle with r3 is " << alfa << " , angle with r2 is " << fi<<endl;



            }


}



void show_possible_r2(vectormath rprimd[3], vectormath norm, vectormath r3) {
    //the function is identical with show_possible_r1; they can be merged into one

    //norm is normal_to_boundary
    cout << "\nHere is the list of possible r2:\n";   
    int u, v, w, max = 12;
    vectormath r2;
    int mul = 1;
    int d[3];
    int uvtw[4], uvw[3];
    db alfa, fi, maxlength = 40,temp ;//40 len_units
    //cout.precision(20);
    for(u = -max;u < max;u++)
        for(v = -max;v < max;v++)
            for(w = -max;w < max;w++) {
                d[0] = u; 
                d[1] = v; 
                d[2] = w;
                r2 = vector_dec_for_UVW(d, rprimd, mul);
                if(abs(r2) > maxlength || abs(r2) == 0) continue;

                temp = r2 * r3 / abs(r2) / abs(r3);
                //reduce the precision of temp 
                temp = roundf(temp * 10000.0f) / 10000.0f; //because temp can be 1.000000002 which cause nan for acos

                alfa = acos(temp) * 180 / M_PI;
                if(alfa < 87 || alfa > 93) continue; //between r2 and r3
 
                temp = norm * r2 / abs(norm) / abs(r2);
                temp = roundf(temp * 10000.0f) / 10000.0f;  
 
                fi = acos(temp) * 180 / M_PI;
                if(fi < 87 || fi > 93) continue; //between r2 and normal
                //temp = r1 * r2 / abs(r1) / abs(r2);
                //cout << "\nfi "<<fi<<" " << temp<< "\n";
                uvw[0] = u; uvw[1] = v; uvw[2] = w;
                three2four(uvtw, uvw);
                cout << "\nVector [" << u << v << w << "] has length = " << abs(r2) \
                <<"; in for index system is ["<<uvtw[0]<<" "<<uvtw[1]<<" "<<uvtw[2]<<" "<<uvtw[3]<<"]\n";
                cout << "\nAngle with r3 is " << alfa << " , angle with normal_to_boundary is " << fi<<endl;



            }


}


void combine_z_with(vectormath  r3, vectormath *arrv, int size_of_arrv) {

    //1. Here a rotation in space around arbitary axis perpendicular to z is used (rotate_axis)
    //Rotate vectors arrv in such a way to combine r3 with z.
    db fi, cosfi, sinfi, omc, rh, rk, rl;
    vector_dec rotate_axis;
    db x[3], y[3], z[3];
    //Determine angle for rotation
    fi = acos( r3.z / abs(r3) );
    if (test == 1) cout << "\nAngle between uvw3 and z is " << fi*180/M_PI;

    //Additional parameters
    cosfi = cos(fi); sinfi = sin(fi); omc = 1-cosfi;
    rotate_axis.x = r3.y; rotate_axis.y = -r3.x; rotate_axis.z=0;
    rh = rotate_axis.x / rotate_axis.length();
    rk = rotate_axis.y / rotate_axis.length();
    rl = 0;
    if(rotate_axis.length() == 0) {rh = 0;rk = 0;}
    //cout << "\nParameters for rotation: rh, rk = " << rh << ", " << rk;
    if (test == 1) cout << "\nLength of rotation axis = " << rotate_axis.length();

    //Rotation
    if (test == 1) cout << "\nrprimd before rotation to z:\n";
    for(int i=0;i<size_of_arrv;i++) {
        x[i] = arrv[i].x;
        y[i] = arrv[i].y;
        z[i] = arrv[i].z;
         
        if (test == 1) cout << arrv[i].x << " " << arrv[i].y << " " << arrv[i].z << endl;
    }

    if (test == 1) cout << "\nrprimd after rotation to z:\n";
    for(int i=0;i<size_of_arrv;i++) {
        arrv[i].x = x[i] * (cosfi + omc * rh * rh) + y[i] * omc * rh * rk + z[i] * sinfi * rk;
        arrv[i].y = x[i] * omc * rh * rk + y[i] * (cosfi + omc * rk * rk) - z[i] * sinfi * rh;
        arrv[i].z = -x[i] * (sinfi * rk) + y[i] * sinfi * rh + z[i] * cosfi;
        if (test == 1) cout << arrv[i].x << " " << arrv[i].y << " " << arrv[i].z << endl;
    }


}



//Constructor
CrystalCell::CrystalCell() {
    periodcell[0] = 0;
    periodcell[1] = 1;
    periodcell[2] = 0;
    periodcell[3] = 1;
    periodcell[4] = 0;
    periodcell[5] = 1;
    firotatez = 0;
    name = "";
  
}  


//void test_directions_and_volume



void CrystalCell::choose_new_rprimd(vectormath &r1, vectormath &r2, vectormath &r3) {
    //Method choose free vectors r1, r2 and r3 for cell with grain (GC) according to three input directions uvtw
    //Also calculate distance between planes parralel to boundary plane



    //Determination of recip and normal_to_boundary
    db Vi = rprimd[0] * (rprimd[1] % rprimd[2]); //initial volume
    recip[0] = rprimd[1] % rprimd[2];
    recip[1] = rprimd[2] % rprimd[0];
    recip[2] = rprimd[0] % rprimd[1];
    recip[0] = recip[0] / Vi;
    recip[1] = recip[1] / Vi;
    recip[2] = recip[2] / Vi;
    normal_to_boundary = vector_dec_for_UVW(hkl1, recip, 1);
    vectormath n = normal_to_boundary;

    //db Vi = rprimd[0] * (rprimd[1] % rprimd[2]); //initial volume
    cout <<"\nThree vectors in 3-index system before mul:";  
    cout << "\nuvw1 is "<<uvw1[0]<<uvw1[1]<<uvw1[2]; 
    cout << "\nuvw2 is "<<uvw2[0]<<uvw2[1]<<uvw2[2]; 
    cout << "\nuvw3 is "<<uvw3[0]<<uvw3[1]<<uvw3[2]<<endl; 

    cout << "\nVolume of initial cell is " << Vi<<endl;

    r1 = vector_dec_for_UVW(uvw1, rprimd, mul[0]);
    r2 = vector_dec_for_UVW(uvw2, rprimd, mul[1]);
    r3 = vector_dec_for_UVW(uvw3, rprimd, mul[2]);

    //Check r3 and normal_to_boundary
    cout <<"\nAngle between r3 and normal is "<< acos( r3 * n / abs(r3) / abs(n) ) * 180 / M_PI;
    //Check r2 and r3
    if(abs(r2 * r3) > 0.00001) {
        cout << "\nWarning! r2 and r3 not perpendicular; Scalar product is "<< r2 * r3 << endl;
    }

    //Show possible r2
    if(abs(r2) == 0) {
        cout << "\nWarning! r2 is zero.";
        show_possible_r2(rprimd, normal_to_boundary, r3);
        cout << "\n\nPlease choose one of the direction for r2 from above\n";
        exit(1);
    }    

    //Show possible r1
    if(abs(r1) == 0) {
        cout << "\nWarning! r1 is zero.";
        show_possible_r1(rprimd,r2,r3);
        cout << "\n\nPlease choose one of the direction for r1 from above\n";
        exit(1);
    }

    cout <<"\nThe following vectors have been chosen:"<<endl;
    cout << r1.x << " " << r1.y << " " << r1.z << endl;
    cout << r2.x << " " << r2.y << " " << r2.z << endl;
    cout << r3.x << " " << r3.y << " " << r3.z << endl;
    //r2.y = r2.y+0.1;
    //r3.z = r3.z+0.1;


    db V = r1 * (r2 % r3); //Volume of new cell
    cout << "\nVolume of constructed cell is " << V<<endl;
    if(V == 0 ) cerr << "\nError!, volume is zero. Something wrong with new vectors, try to change directions";
    if(V <= 0 ) cerr << "\nWarning!, Volume is negative. May be the vectors form left system";

    int natom_new = round( V / (Vi / natom) ); // round up
    if(test == 1 ) {
 
        cout << "\nHere is possible sizes of the one grain. Use mul parameter to change sizes\n";
        cout << "\nLength of r1 " << abs(r1) << len_units;
        cout << "\nLength of r2 " << abs(r2) << len_units;
        cout << "\nLength of r3 " << abs(r3) << len_units;
        db fi = acos( r1 * r2 / abs(r1) / abs(r2) );
        cout << "\nAngle between r1 and r2 is " << fi * 180 / M_PI;
        cout << "\nDistance between grain boundary planes due to PBC is "<< abs(r1) * cos( fi- (M_PI / 2) );
        cout << "\nApproximate number of atoms: " << natom_new;

        cout << "\nTest 1 finished. Change test parameter to 0 for GB construction\n";
    }





}


void CrystalCell::readin(string inname) {
    char filnamin[fnlen];
    strcpy(filnamin,inname.c_str());

    cout <<"\nName of input file "<<filnamin<<"\n";

    int len_filnamin = strlen(filnamin);
//    db acell[4];
//    int natom;
    readinnatom_(filnamin,\
    &natom);
//    cout << "natom" <<natom;
    if(natom>max_at) cout <<"Error, natom > max_at";
//    xred=new db*[4]; for(int i=0; i < 4; i++) xred[i]=new db[natom+1]; //allocate memory for xred

    int uvtw1[3], uvtw2[3] , uvtw3[3];

    readin_(filnamin,\
    acell,rprim,xred,typat, znucl,
    periodcell,&firotatez,\
    hkl1, uvtw1, uvtw2, uvtw3, uvw1, uvw2, uvw3, mul, &test,\
    slicepms, shiftatoms, &deldist, shiftcell, shiftgrain, a_c_conv, sixvalues);



    version = (int)a_c_conv[0]; 
    hex_a = acell[0]; //works only for special case of hcp
    hex_c = acell[2];

 
    cout <<"\nuvw1 is "<<uvw1[0]<<uvw1[1]<<uvw1[2];
    cout <<"\nuvw2 is "<<uvw2[0]<<uvw2[1]<<uvw2[2];
    cout <<"\nuvw3 is "<<uvw3[0]<<uvw3[1]<<uvw3[2];
    if(test == 1) {
        cout<<"\nshiftatoms ";
        for(int i=0;i<4;i++)
            cout<<shiftatoms[i]<<" ";
        cout<<"\nshiftgrain "<<shiftgrain[0]<<" "<<shiftgrain[1];
    }
    if (uvw1[0] == 0 && uvw1[1] == 0 && uvw1[2] == 0 ) four_index_to_three(uvtw1, uvw1);
    if (uvw2[0] == 0 && uvw2[1] == 0 && uvw2[2] == 0 ) four_index_to_three(uvtw2, uvw2);
    if (uvw3[0] == 0 && uvw3[1] == 0 && uvw3[2] == 0 ) four_index_to_three(uvtw3, uvw3);
    //Calculate rprimd
    for(int i=0; i<3; i++) {
        rprimd[i].x = rprim[i][0] * acell[i];
        rprimd[i].y = rprim[i][1] * acell[i];
        rprimd[i].z = rprim[i][2] * acell[i];
    }



    //Calculate direction from plane if direction is zero
//    if(uvw1[0] == 0 && uvw1[1] == 0 && uvw1[2] == 0)
        //        plane_to_direction(uvw1, hkl1, rprimd); //this may not to work. 

    //Fill in sort array
    //Should be carefully checked for other lattices 
    cout <<"\nWarning! Closed packed layers in sort[][1] filled in correctly only in the particular case of \
    input hcp lattice with two atoms!\n";
    for(int i=0;i<natom;i++) {
        if( i/2 * 2 == i ) sort[i][1] = 0;
        if( i/2 * 2 != i ) sort[i][1] = 1;
    }

    //Cope with name
    name = string(filnamin);
    name.resize(100); //to remove some strange symbols at the end
    //cout <<"\nName of output structure before trimming "<<name<<";\n";
    name.erase( remove(name.begin(), name.end(), ' '), name.end() );
    name = read_plus( inname, "name" );
    cout <<"\nName of output structure "<<name<<";\n";

    calctype = read_plus( inname, "calctype" ); //read type of calculation
    return;
}

void CrystalCell::periodicfill_and_cut(vectormath r1, vectormath r2, vectormath r3,db bob,db upb) {
//Вычисляет новые xcart и natom после копирования ячейки periodcell[] раз.
    //Makes grain A
    int nh = periodcell[1] - periodcell[0];
    int nk = periodcell[3] - periodcell[2];
    int nl = periodcell[5] - periodcell[4];
    if(nh < 0 || nk < 0 || nl < 0) cerr <<"Error periodcell[] are wrong!!!";
    double basis_dec[3];
    int i_typat[max_at], i_sort[max_at][n_sort];
    vectormath t_xcart, tx;
    vectormath r[3]; r[0] = r1; r[1] = r2; r[2] = r3;
    int i, u, v, w;//Номера узлов решетки
    int i_at = 0, tat; // счетчик по новым атомам
    //db bob = 0.001, upb = 1.02;
    db i_xred[max_at][3];
    for(i=0;i<natom;i++) {
        i_xred[i][0] = xred[i][0];
        i_xred[i][1] = xred[i][1];
        i_xred[i][2] = xred[i][2];
        i_typat[i] = typat[i];
        //cout << i_typat[i] <<" ";
        i_sort[i][0] = 0; //Grain A
        i_sort[i][1] = sort[i][1];
    }
    //Используя заданные максимальные размеры системы пройти
    //по всем атомам,определить номера узловых точек ячейки и рассчитать декартовы координаты всех атомов.
    for(u = periodcell[0];u < periodcell[1];u++)
        for(v = periodcell[2];v < periodcell[3];v++)
            for(w = periodcell[4];w < periodcell[5];w++)

    //for(u = 0;u <= uvw1[0];u++)
    //    for(v = 0;v <= uvw[1];v++)
    //        for(w = periodcell[4];w < periodcell[5];w++)
                for(i=0;i<natom;i++) {
                    //Вычисление декартовых координат атомов базиса
                    basis_dec[0] = i_xred[i][0] * rprimd[0].x + i_xred[i][1] * rprimd[1].x + i_xred[i][2] * rprimd[2].x;
                    basis_dec[1] = i_xred[i][0] * rprimd[0].y + i_xred[i][1] * rprimd[1].y + i_xred[i][2] * rprimd[2].y;
                    basis_dec[2] = i_xred[i][0] * rprimd[0].z + i_xred[i][1] * rprimd[1].z + i_xred[i][2] * rprimd[2].z;
                    //Вычисление новых координат
                    t_xcart.x = u * rprimd[0].x + v * rprimd[1].x + w * rprimd[2].x + basis_dec[0];
                    t_xcart.y = u * rprimd[0].y + v * rprimd[1].y + w * rprimd[2].y + basis_dec[1];
                    t_xcart.z = u * rprimd[0].z + v * rprimd[1].z + w * rprimd[2].z + basis_dec[2];
                    
                    //Включаем, только те атомы, условные координаты которых попадают в новую ячейку
                    tx = calculate_xred(r, t_xcart.x, t_xcart.y, t_xcart.z);
                    tat = i_typat[i];


                    //cout << tat <<" ";
                    if( tx.x >= bob && tx.x < upb && \
                        tx.y >= bob && tx.y < upb && \
                        tx.z >= bob && tx.z < upb ) {


                        xcart[i_at][0] = t_xcart.x;
                        xcart[i_at][1] = t_xcart.y;
                        xcart[i_at][2] = t_xcart.z;
                        xred[i_at][0] = tx.x;
                        xred[i_at][1] = tx.y;
                        xred[i_at][2] = tx.z;
                        typat[i_at] = tat;
                        //cout << tat <<" "<<endl;
                        sort[i_at][0] = i_sort[i][0];
                        sort[i_at][1] = i_sort[i][1];

                        i_at++;
                        if(i_at > max_at) cerr <<"Error periodcell[] are too big! Encrease max_at at cpp and in fortran files\n";

                    }
                    //cout << endl;
                    //d[s+i]=sqrt(x[s+i]*x[s+i]+y[s+i]*y[s+i]+z[s+i]*z[s+i]);//расстояние из точки 0 до каждого узла (по сути атома)
                }
    natom = i_at;

    if(test == 1) cout << "\n'periodicfill_and_cut': New cell contains " << natom << " atoms.\n";
    //Calculate rprimd
    rprimd[0] = r1;
    rprimd[1] = r2;
    rprimd[2] = r3;
}





void CrystalCell::rotate(vectormath &r1, vectormath &r2, vectormath &r3) {
    //Two step combination of directions with axes
    //1. Combine uvw3 with z
    //2. Try to combine uvw2 with y


    //1. Here is a rotation in space around arbitary axis perpendicular to z is used (rotate_axis)

    combine_z_with(r3, rprimd, 3);

    r1 = vector_dec_for_UVW(uvw1, rprimd, mul[0]);
    r2 = vector_dec_for_UVW(uvw2, rprimd, mul[1]);
    r3 = vector_dec_for_UVW(uvw3, rprimd, mul[2]); 

    if(test == 1) {
    cout << "\nTest of rotation 1, Lengths of vectors after first rotation obtained from uvw and new rprimd, \
     should be the same \n";
    cout << "\nLength of r1 " << abs(r1) << len_units;
    cout << "\nLength of r2 " << abs(r2) << len_units;
    cout << "\nLength of r3 " << abs(r3) << len_units;
    }


    //Determination of recip and normal to boundary after rotation
    db Vi = rprimd[0] * (rprimd[1] % rprimd[2]); //initial volume
    if(test == 1)    cout << "\nVolume of initial cell after first rotation is " << Vi;
    recip[0] = rprimd[1] % rprimd[2];
    recip[1] = rprimd[2] % rprimd[0];
    recip[2] = rprimd[0] % rprimd[1];
    recip[0] = recip[0] / Vi;
    recip[1] = recip[1] / Vi;
    recip[2] = recip[2] / Vi;
    normal_to_boundary = vector_dec_for_UVW(hkl1, recip, 1);

    //Determination of spacing. This is 1/length, where length is length of hkl1
    spacing = 1 / abs(normal_to_boundary);
    if(test == 1) cout << "\nSpacing between (" << hkl1[0]<<hkl1[1]<<hkl1[2] <<") planes is " << spacing <<len_units;



    //2. Here rotation around z is used.
    vectormath x, y, z;  
    db fi;
    
    
    //Check if it is possible to  match r2 and y by using triple product. If triple product is
    //not equal to zero then r2 lies outside xy plane and it is impossible to match it with
    //y by rotation around z.
    x.x = 1; x.y = 0; x.z = 0;
    y.x = 0; y.y = 1; y.z = 0;
    z.x = 0; z.y = 0; z.z = 1;
    if( abs(x * (r2 % y)) > 1e-10 ) {
        cout << "\nr1, r2, r3 after first rotation to z:\n";
        cout << r1.x << " " << r1.y << " " << r1.z << endl;
        cout << r2.x << " " << r2.y << " " << r2.z << endl;
        cout << r3.x << " " << r3.y << " " << r3.z << endl;
        cout << "\nThe triple product of x, r2 and y is "<< x * (r2 % y);
        cout << "\nWarning! Imposible to match r2 and y, as r2 does not lie in xy plane. Nothing done\n";

        //return;
    }//may be no more needed, check!!!

    //Additional check if normal_to_boundary is parrallel to r2 % r3
    db par = abs( normal_to_boundary % (r2 % r3) );
    vectormath n = normal_to_boundary;
    if ( abs(par) > 0.0001) 
        cout <<"\nWarning! normal_to_boundary is not parrallel to double product of r2 and r3; par is " << par <<"\n";
    
    db temp = r1 * n / abs(r1) / abs(n) ;
    temp = roundf(temp * 10000.0f) / 10000.0f; //because temp can be 1.000000002 which cause nan for acos
    cout <<"\nAngle between r1 and normal before last rotation is "<< acos( temp ) * 180 / M_PI;
    cout <<"\nAngle between r2 and normal before last rotation is "<< acos( r2 * n / abs(r2) / abs(n) ) * 180 / M_PI;
    cout <<"\nAngle between r3 and normal before last rotation is "<< acos( r3 * n / abs(r3) / abs(n) ) * 180 / M_PI;
    

    //Here we match normal_to_boundary with x
    //Determine angle between x and normal_to_boundary:
    //fi = -acos( normal_to_boundary.x / abs(normal_to_boundary) );
    if (normal_to_boundary.y < 0)
        //Знак угла зависит от того в какой четверти находится вектор; для углов более 180 градусов поворот работает не верно. Знак углов похоже не совпадает с общепринятым (отсчитывается по часовой стрелке, а должен против)
        fi =  acos( normal_to_boundary.x / abs(normal_to_boundary) );
    else
        fi =  -acos( normal_to_boundary.x / abs(normal_to_boundary) );        

    if (test == 1) cout << "\nnormal_to_boundary before rotation is "; n.print();
    if (test == 1) cout << "\nAngle between normal_to_boundary and x before rotation is " << fi*180/M_PI;

    //rotate
    db xt,yt,zt, cosfi = cos(fi), sinfi = sin(fi);
    for(int i=0;i<3;i++) {
        xt = rprimd[i].x;  yt = rprimd[i].y;
        rprimd[i].x = xt * cosfi - yt * sinfi;
        rprimd[i].y = xt * sinfi + yt * cosfi;
    }


    r1 = vector_dec_for_UVW(uvw1, rprimd, mul[0]);
    r2 = vector_dec_for_UVW(uvw2, rprimd, mul[1]);
    r3 = vector_dec_for_UVW(uvw3, rprimd, mul[2]); 

    //Determination of recip and normal to boundary after second rotation
    Vi = rprimd[0] * (rprimd[1] % rprimd[2]); //initial volume
    recip[0] = rprimd[1] % rprimd[2];
    recip[1] = rprimd[2] % rprimd[0];
    recip[2] = rprimd[0] % rprimd[1];
    recip[0] = recip[0] / Vi;
    recip[1] = recip[1] / Vi;
    recip[2] = recip[2] / Vi;
    normal_to_boundary = vector_dec_for_UVW(hkl1, recip, 1);
    n = normal_to_boundary;


    temp = r1 * n / abs(r1) / abs(n) ;
    temp = roundf(temp * 10000.0f) / 10000.0f; //because temp can be 1.000000002 which cause nan for acos
    cout <<"\nAngle between r1 and normal after rotation is "<< acos( temp ) * 180 / M_PI;
    cout <<"\nAngle between r2 and normal after rotation is "<< acos( r2 * n / abs(r2) / abs(n) ) * 180 / M_PI;
    cout <<"\nAngle between r3 and normal after rotation is "<< acos( r3 * n / abs(r3) / abs(n) ) * 180 / M_PI;
    
    if (test == 1) cout << "\nnormal_to_boundary after rotation is "; n.print();
    //Determine angle between x and normal_to_boundary after rotation:
    fi = -acos( normal_to_boundary.x / abs(normal_to_boundary) );
    if (test == 1) cout << "\nAngle between normal_to_boundary and x after rotation is " << fi*180/M_PI;











    if(test == 1) {
        cout << "\nTest of rotation 2, Lengths of vectors after second rotation obtained from uvw and new rprimd, \
        should be the same \n";
        cout << "\nLength of r1 " << abs(r1) << len_units;
        cout << "\nLength of r2 " << abs(r2) << len_units;
        cout << "\nLength of r3 " << abs(r3) << len_units;

        cout << "\nrprimd after second rotation around z:\n";
        for(int i=0;i<3;i++) 
            cout << rprimd[i].x << " " << rprimd[i].y << " " << rprimd[i].z << endl;
        Vi = rprimd[0] * (rprimd[1] % rprimd[2]); //initial volume
        cout << "\nVolume of initial cell after second rotation is " << Vi;
        cout.precision(5);
        cout << std::fixed;
        cout << "\nr1, r2, r3 after second rotation around z:\n";
        cout << r1.x << " " << r1.y << " " << r1.z << endl;
        cout << r2.x << " " << r2.y << " " << r2.z << endl;
        cout << r3.x << " " << r3.y << " " << r3.z << endl;
    }


}







void CrystalCell::cut(limits lim)
{
//Вырезаем нужную ячйку

int i2=0,j;
db xcart2[max_at][3],typat2[max_at];

for(int i=0;i<natom;i++)
    {
//Выбираем только те атомы, которые попадают в обозначенную ячейку. По оси z выбраются от 0 до вектора с-т.е всего два ГПУ слоя.
//Здесь по x и z не включаются узлы через период, однако добавлен узел по периоду y, так как он потом используется
//Верхние коментарии теперь мало имеют значения, так теперь с учетом переиодических условий

    if(xcart[i][0] >= lim.xmin-1e-11 \
	&& xcart[i][0] <= lim.xmax+1e-11 \
    && xcart[i][1] >= lim.ymin \
    && xcart[i][1] <= lim.ymax \
	&& xcart[i][2] >= lim.zmin \
	&& xcart[i][2] <  lim.zmax)
	    {
//cout<<"[uvw], s="<<i<<" k_nodes=" <<k_nodes<<"\n";
//cout<<h[i]<<" "<<k[i]<<" "<<l[i]<<"\n";//индексы hkl есть только для узлов (четные i!)
//cout<<x_new[i]<<" "<<y_new[i]<<"\n";//В данном случае для всех значений i есть координаты
        for(j=0; j<3; j++)  
            xcart2[i2][j] = xcart[i][j];
//        cout <<xcart2[i2][0]<<" "<<xcart2[i2][1]<<" "<<xcart2[i2][2]<<endl;             
        typat2[i2] = typat[i];
        i2++;
        }
    }

natom=i2;
for(int i=0; i<natom; i++)
    {
    typat[i] = typat2[i];
    for(int j=0; j<3; j++)  
        xcart[i][j] = xcart2[i][j];
    }

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



//Задаем acell and rprim
acell[0]=lim.xmax - lim.xmin;
acell[1]=lim.ymax - lim.ymin;
acell[2]=lim.zmax - lim.zmin;

for(int i=0; i<3; i++)
    for(int j=0; j<3; j++)
        {
        if(i==j) rprim[i][j] = 1;
        if(i!=j) rprim[i][j] = 0;
        } 

    for(int i=0; i<3; i++) {
        rprimd[i].x = rprim[i][0] * acell[i];
        rprimd[i].y = rprim[i][1] * acell[i];
        rprimd[i].z = rprim[i][2] * acell[i];
    }    

}












db CrystalCell::make_boundary(db delta_v = 0) {
    //Строим ячейку с границей путем зеркального отражения
    //Input: delta_v - this specific volume added to each grain boundary
    //n_slice is number of planes which should be sliced from each parts of grain A
    //n_spacing is number of minimum plane spacings inside spacing
    
    //n_relative is integer number for choosing relative position of two boundaries; 
    //in fact makes slice only from rigth side; influence on size of grain
    
    //Makes reflection relative to yz plane and use rprimd[0].x size



    //First step is to make the same slices from both sides of 
    //initial grain
    db n_slice = slicepms[0];
    db n_relative = slicepms[1]; 
    db n_spacing = slicepms[2];
    int i_at = 0;
    db new_xcart[max_at][3];
    int new_typat[max_at], new_sort[max_at][n_sort];
    db min_spacing = spacing / n_spacing; //Minimum spacing between planes. Should be choosen by hands.
    db slice = n_slice* min_spacing; //
    db width_of_slice_left = slice + 0 * min_spacing - 0.001; // was 0.000001
    //at right we have additional slice parameter to make slicing from both sides symmetrical
    db width_of_slice_right = slice + n_relative*min_spacing - 0.001;//0.0001;  // little additional reduce to the layer of atoms was 0.01 

    //Calculate volume per one atom in the initial ideal bulk grain
    db V_1at = ( rprimd[0] * (rprimd[1] % rprimd[2]) ) / (natom*1.)  ; 
    //cout<<"\ndelta_v inside make_boundary "<<delta_v;


    //Make sliceing from left and right
    //write new coordinates in new_xcart
    for(int i=0;i<natom;i++) {
        if( xcart[i][0] < width_of_slice_left || \
            xcart[i][0] > (rprimd[0].x - width_of_slice_right) ) continue; 
        
        new_xcart[i_at][0] = xcart[i][0];
        new_xcart[i_at][1] = xcart[i][1];
        new_xcart[i_at][2] = xcart[i][2];
        new_typat[i_at] = typat[i];
        new_sort[i_at][1] = sort[i][1];      
        i_at++;
    }
    natom = i_at;
    rprimd[0].x = rprimd[0].x - width_of_slice_right;
    if(test == 1) cout << "\nmake_boundary: After slicing cell contains " << natom << " atoms.\n";

    //Additional volume for first boundary; the most right boundary
    gbpos = rprimd[0].x + delta_v / 4.; // geometrical position of first gb plane

    rprimd[0].x = rprimd[0].x + delta_v / 2.; //because the cell will be doubled, we add only half here     





    //Make mirroring of initial grain along x.
    rprimd[0].x = rprimd[0].x * 2  ;
    //cout << "\nCheck "; rprimd[0].print();
    if(natom * 2 > max_at) cerr << "\nError! Two much atoms. Check Mirrorx";
    for(int i=0;i<natom;i++) {
        //For grain A
        xcart[i][0] = new_xcart[i][0];
        xcart[i][1] = new_xcart[i][1];
        xcart[i][2] = new_xcart[i][2];
        typat[i] = new_typat[i];
        sort[i][0] = 0; //Grain A.
        sort[i][1] = new_sort[i][1]; //atoms of closed packed layer 
        //For grain B - mirror of A
        xcart[i + natom][0] = rprimd[0].x - new_xcart[i][0];
        xcart[i + natom][1] = new_xcart[i][1];
        xcart[i + natom][2] = new_xcart[i][2];
        typat[i + natom] = new_typat[i];
        sort[i + natom][0] = 1; //Grain B.
        sort[i + natom][1] = new_sort[i][1];
    }
    natom = natom * 2;

    //Here second adjusting of cell size to make the second boundary the same;
    rprimd[0].x = rprimd[0].x + 0 * min_spacing - 2 * slice;
    rprimd[0].y = 0.0; // because mirror  is perpendicular to x. What about rprimd[0].z ???

    //Additional volume for second boundary
    rprimd[0].x = rprimd[0].x + delta_v; //here we add the whole delta_v  


    //Calculate xred after the selection 
    vectormath new_xred;
    for(int i=0;i<natom;i++) {    
        new_xred = calculate_xred(rprimd, xcart[i][0], xcart[i][1], xcart[i][2]);
        //cout <<  new_xred.x<<" "<<new_xred.y<<" "<<new_xred.z<<" "<<"\n";
        xred[i][0] = new_xred.x;
        xred[i][1] = new_xred.y;
        xred[i][2] = new_xred.z;
    }
    if(test == 1) cout << "\nAfter mirroring cell contains " << natom << " atoms.\n";

    shift_grain();
    
    delnear();//delete close lying atoms
    delnear();
    
    if (shiftatoms[2] != 0 && shiftatoms[3] != 0)
        //Only shift in the case then the shift vector is provided and abs value of shift is > 0
        shift_atoms();// input parameter shift_atoms
    


    shift_cell(); //
    return_atoms_to_cell();

    //Calculate total volume of cell with grain boundary
    db V_gb = ( rprimd[0] * (rprimd[1] % rprimd[2]) )  ; 

    //Calculate specific volume of grain boundary =  difference between total volumes of cell with boundary and ideal cell divided by the surface area of boundary

    db v_gb = (V_gb - V_1at * natom) / abs(rprimd[1] % rprimd[2]) / 2.;
    cout.precision(6);
    //cout <<"\nV_gb inside make_boundary "<<V_gb;
    //cout <<"\nV_1at "<<V_1at;
    cout <<"\nnatom "<<natom;
    //cout <<"\nSurface area " <<abs(rprimd[1] % rprimd[2]);
    //cout <<"\nv_gb inside make_boundary "<<v_gb<<endl;
    return v_gb;
}



void CrystalCell::shift_grain() {
    //Формат во входном файле shiftgrain y z
    //Grain A is shifted by:
    //y - along r2 direction (parts of r2)
    //z - along r3 direction (parts of r3)
    // int y = shiftgrain[0];
    // int z = shiftgrain[1];
    cout <<"\nshift_grain(): The grain A will be shifted by " << shiftgrain[0]*100 << " and "<< shiftgrain[1]*100 <<" % along r2 and r3\n";


    for(int i=0;i<natom;i++) {
        //cout << sort[i][1]<<endl;
        if(sort[i][0] == 0 ) {
            xred[i][1] = shiftgrain[0] + xred[i][1];
            xred[i][2] = shiftgrain[1] + xred[i][2];
            //cout << "\natom shifted";
            //calculate xcart from xred
            xcart[i][0] = xred[i][0] * rprimd[0].x + xred[i][1] * rprimd[1].x + xred[i][2] * rprimd[2].x;
            xcart[i][1] = xred[i][0] * rprimd[0].y + xred[i][1] * rprimd[1].y + xred[i][2] * rprimd[2].y;
            xcart[i][2] = xred[i][0] * rprimd[0].z + xred[i][1] * rprimd[1].z + xred[i][2] * rprimd[2].z;
        }
    }

    return;
}


void CrystalCell::mirrory() {
//Строим ячейку с границей путем зеркального отражения
    //Only for cubic cell
    cout << "\nWarning! mirrory() works only for cubic cells";
rprimd[1].y = rprimd[1].y * 2;
for(int i=0;i<natom;i++)
    {
    xcart[i + natom][0] = xcart[i][0];
    xcart[i + natom][1] = rprimd[1].y-xcart[i][1];
    xcart[i + natom][2] = xcart[i][2];
    typat[i + natom] = typat[i];
    }
natom = natom * 2;
}



int CrystalCell::delnear() {
    //works only for cubic cells
    //"minimum_distance" is minimum distance between atoms 
bool checked[max_at];
bool at_to_delete[max_at];
memset(checked, 0, natom);
memset(at_to_delete, 0, natom);
db sizex = rprimd[0].x, sizey = rprimd[1].y, sizez = rprimd[2].z;
db sizex05=sizex*0.5, sizey05=sizey*0.5, sizez05=sizez*0.5;
db rr, d, r = deldist; //Class attribute
db xx,yy,zz;
int natom2 = natom;
int j=0;
for (int i=0;i<natom;i++)
    {
    checked[i]=1;
	for (j=0;j<natom;j++)
	    {
        if(!checked[j])
            {
            xx = xcart[i][0] - xcart[j][0]; yy=xcart[i][1] - xcart[j][1]; zz=xcart[i][2] - xcart[j][2];
            if(xx > sizex05) xx=xx-sizex;
    	    if(xx < -sizex05) xx=xx+sizex;
            if(yy > sizey05) yy=yy-sizey;
	        if(yy < -sizey05) yy=yy+sizey;
            if(zz > sizez05) zz=zz-sizez;
	        if(zz < -sizez05) zz=zz+sizez;
            rr=xx*xx+yy*yy+zz*zz;
            d=sqrt(rr);
            if(d < r)
                { //если расстояние между атомами меньше 3 Бор
//	            cout<<"i="<<i<<" coord was:"<<x_cell[i]<<" "<<y_cell[i]<<" "<<z_cell[i]<<endl;
	            xcart[i][0] = xx*0.5 + xcart[j][0];//Позволяет найти среднее положение, т
	            xcart[i][1] = yy*0.5 + xcart[j][1];
	            xcart[i][2] = zz*0.5 + xcart[j][2];
	            at_to_delete[j] = 1;
	            if(test == 1) cout <<"atom "<<j<<" was deleted"<<endl;
//	            natom2--;
//	cout<<"i="<<i<<" coord become:"<<x_cell[i]<<" "<<y_cell[i]<<" "<<z_cell[i]<<endl;
//cout<<"atom j="<<j<<" will be deleted, because it is to close to atom i="<<i<<" ,d="<<d<<endl;
//cout<<"j  coord:"<<x_cell[j]<<" "<<y_cell[j]<<" "<<z_cell[j]<<endl;
//cout<<"i coord:"<<x_cell[i]<<" "<<y_cell[i]<<" "<<z_cell[i]<<endl;
 	 	        }
			}
		}
	}

int i2=0;
for (int i=0;i<natom;i++) {
    if(at_to_delete[i]) continue;
    xcart[i2][0] = xcart[i][0]; 
    xcart[i2][1] = xcart[i][1]; 
    xcart[i2][2] = xcart[i][2];
    typat[i2] = typat[i];
    sort[i2][0] =  sort[i][0];
    sort[i2][1] =  sort[i][1];
    i2++;
}
//if(i2!=natom2) cout << "Error in distance\n";
natom = i2;
    //calculate xred;
    vectormath new_xred;
    for(int i=0;i<natom;i++) {    
        new_xred = calculate_xred(rprimd, xcart[i][0], xcart[i][1], xcart[i][2]);
        //cout <<  new_xred.x<<" "<<new_xred.y<<" "<<new_xred.z<<" "<<"\n";
        xred[i][0] = new_xred.x;
        xred[i][1] = new_xred.y;
        xred[i][2] = new_xred.z;
    }
    if(test == 1) cout << "\nAfter deleting near atoms cell contains " << natom << " atoms.\n";


    //cout <<"\nTable for latex, sizes in Angstroms and number of atoms\n";
    cout.precision(3);
    if(test == 1) cout <<name<<" & " <<abs(rprimd[0]) * to_angstrom <<" & "<< abs(rprimd[1]) * to_angstrom << \
                    " & " << abs(rprimd[2]) * to_angstrom << " & " << natom << "\n\n";
return natom;
}




















void CrystalCell::write_xyz(string add) {
    //Determine number of atoms of different sorts
    //Only for grains, turned off
    int choosensort = 0;
    int nstates = 1;
    int value[2] = {0,1};
    int num_of_atoms[2] = {0,0};
    int z;
    for(int j=0;j<nstates;j++)
        for(int i=0;i<natom;i++)
            if(sort[i][choosensort] == value[j] || 1) num_of_atoms[j]++;
    //cout << "\nNumber of atoms in Grain A " << num_of_atoms[0]<<endl;

    ofstream out_xyz; // открываем файл для записи
    //cout << "out/"+name+".xyz\n";

    //string folder = "out/"+name+"/"+add+"/xyz/";
    string folder = "xyz/";
    //system( ("mkdir -p "+folder).c_str() );
    //out_xyz.open(( folder+name+add+"."+to_string(version)+".xyz").c_str(), ios::out);

    out_xyz.open(( folder+name+add+"."+to_string(version)+".xyz").c_str(), ios::out);

    for(int j=0;j<nstates;j++) {
        out_xyz << num_of_atoms[j] << endl;
        out_xyz << "Model "<< j+1 << endl;
        for(int i=0;i<natom;i++) {

            //C переводом в Ангстремы
            if (sort[i][choosensort] == value[j] || 1) {
                z = (int)znucl[typat[i]-1];
                if      (z == 22) out_xyz    <<"Ti ";
                else if (z == 74) out_xyz    <<"W ";
                else if (z == 26) out_xyz    <<"Fe ";
                else if (z == 27) out_xyz    <<"Co ";
                else if (z == 8)  out_xyz    <<"O ";
                else if (z == 6)  out_xyz    <<"C ";
                else              out_xyz    <<"Pu ";
                
                out_xyz << xcart[i][0]*to_angstrom <<" "<< \
                           xcart[i][1]*to_angstrom <<" "<< \
                           xcart[i][2]*to_angstrom << endl;
            }
        }
    }
    out_xyz.close();

    //Write for jmol additional information about grains
    for(int i=0;i<natom;i++) {
        //if(sort[i][0] == 0) cout << "Ti" << i+1 << ",";
        ;
    }
    cout<<endl;

}



void CrystalCell::outcell() {

    cout << "\nperiodcell = "<<periodcell[0]<<" "<<periodcell[1]<<" "<<periodcell[2]<<" "<<periodcell[3]<<" "<<\
    periodcell[4]<<" "<<periodcell[5]<<" \n";
    cout << "\nacell = "<<acell[0]<<" "<<acell[1]<<" "<<acell[2]<<endl;

    cout << "\nrprim  = \n";
    cout <<rprim[0][0]<<" "<<rprim[0][1]<<" "<<rprim[0][2]<<endl;
    cout <<rprim[1][0]<<" "<<rprim[1][1]<<" "<<rprim[1][2]<<endl;
    cout <<rprim[2][0]<<" "<<rprim[2][1]<<" "<<rprim[2][2]<<endl;

    cout << "natom = " <<natom<<endl;

//    cout <<"\nxred[i,j]\n";
//    for(int i=0; i<natom;i++)
//        cout <<xred[i][0]<<" "<<xred[i][1]<<" "<<xred[i][2]<<endl;

    cout <<"\nxcart[i,j]\n";
    for(int i=0; i<natom;i++)
        cout <<xcart[i][0]<<" "<<xcart[i][1]<<" "<<xcart[i][2]<<endl;

    cout <<"\ntypat[i]\n";
        for(int i=0; i<natom;i++)
        cout <<typat[i]<<" ";    

    //cout xcoordinate of type 3 objects (intersitial sites or atoms)
    cout <<"\nType 3 xcart:\n";
    for(int i=0; i<natom;i++)
        if(typat[i] == 3) {
        //if(xcart[i][1]> (acell[1]*0.25 - 1) && xcart[i][1]< (acell[1]*0.25 + 1))//Кординаты между границами
        //          cout <<xcart[i][0]<<" "<<xcart[i][1]<<" "<<xcart[i][2]<<endl;

        //        if(xcart[i][1]> (11) && xcart[i][1]< (12)) //поиск 601 положения
        //          cout <<xcart[i][0]<<" "<<xcart[i][1]<<" "<<xcart[i][2]<<endl;
            if(xcart[i][1]>  (acell[1]*0.5 - 3) && xcart[i][1]< (acell[1]*0.5 + 3)) //поиск междоузлий на границе 601 положения
                cout <<xcart[i][0]<<" "<<xcart[i][1]<<" "<<xcart[i][2]<<endl;
        }


    cout <<"\ntest is: "<<test;
    cout << "\nuvw1 is "<<uvw1[0]<<uvw1[1]<<uvw1[2]; 
    cout << "\nuvw2 is "<<uvw2[0]<<uvw2[1]<<uvw2[2]; 
    cout << "\nuvw3 is "<<uvw3[0]<<uvw3[1]<<uvw3[2]; 
    cout <<"\n__________________________________\n";   
}

void CrystalCell::findpores(db ri)
{
db maxri;//не используется, была попытка вывести все поры меньше определенного размера,
//однако это не так оказалось просто, так как из-за сдвига центра, большая пора
//распознается как маленькая и её нельзя отсечь без дополнительной проверки, которая оказывается
//достаточно сложной. Поэтому в текущем виде программа позволяет найти поры, радиус которых
//не меньше заданного.
//В титане Примерный радиус атома титана ~ 2.8 Бор. углерода 1.1 Бор
db rm = 2.7;//,ri = 0.75;//1.15//Радиус атомов матрицы, радиус атомов примеси.
//набор rm = 2.7 ri = 1.15 и scans = 0.1 позволяет достаточно хорошо находить октапоры в чистом ГПУ.
limits bounds; 
bounds.xmin = 0;
bounds.xmax = acell[0];
bounds.ymin = 0;//acell[1]*0.5 - 6*rm;//примерно По два слоя атомов с каждой стороны границы
bounds.ymax = acell[1];//*0.5 + 6*rm;
bounds.zmin = 0;
bounds.zmax = acell[2];

int i,j,k,i2=0;
//Пока закоментировано, так как работает не всегда правильно. Отрезает атомы, у которых были координаты
//меньше нуля. Нужно сделать так, чтобы область выбиралась всегда правильно!
/*for(i=0;i<natom;i++)
    {
    if(xcart[i][0] >= bounds.xmin \
    && xcart[i][0] <= bounds.xmax \
    && xcart[i][1] >= bounds.ymin \
    && xcart[i][1] <= bounds.ymax \
    && xcart[i][2] >= bounds.zmin \
    && xcart[i][2] <= bounds.zmax)
	    {
        for(j=0; j<3; j++)  
            xcart[i2][j] = xcart[i][j];
        typat[i2] = typat[i];
        i2++;
        }
    }
natom = i2;
*/

db sizex=acell[0], sizey=acell[1], sizez=acell[2];
db sizex05=sizex*0.5, sizey05=sizey*0.5, sizez05=sizez*0.5;
//Теперь пройтись по выбранной области с малым шагом
db scans = 0.1;//Шаг сканирования в Бор ; 1 //0.2

//Размеры сетки
db lx =  bounds.xmax - bounds.xmin; //нужны период условия, пока жесткие стенки
db ly =  bounds.ymax - bounds.ymin;// - 2*rm;//-2rm - Так как с каждой стороны стенки могли 
//разрезать атомы и в итоге добавить свободного места
db lz =  bounds.zmax - bounds.zmin; //нужны период условия, пока жесткие стенки
int nx = lx/scans;
int ny = ly/scans;
int nz = lz/scans;
cout <<"\nfindpore: nx, ny, nz = "<<nx<<" "<<ny<<" "<<nz<<endl;
//cout <<"\nfindpore: natom before= "<<natom<<endl;
//Пока без периодических условий
int n,ii=0,natom2=natom,npore=0;
db x,y,z,xx,yy,zz,dsqr,dsqrmin=100,refsqr,maxrefsqr;
refsqr = (rm+ri)*(rm+ri);//Минимальное расстояние между атомом матрицы и атомом примеси
maxrefsqr = (maxri+rm) * (maxri+rm);
//Попытка сделать один массив ни к чему не привела, ускорение не заметно, но требует больших массивов
//Динамическое создание больших массивов занимает очень много времени 
cout <<"\nfindpore: maxrefsqr (Bohr^2) = "<<maxrefsqr<<endl;
cout <<"\nfindpore: refsqr (Bohr^2) = "<<refsqr<<endl;
for(i=0; i<nx; i++) 
    for(j=0; j<ny; j++) 
        for(k=0; k<nz; k++)
            {
            dsqrmin=100;
            //для ускорения счета операцию суммирования можно убрать, сдвинув все атомы до циклов
            x = i*scans + bounds.xmin;
            y = j*scans + bounds.ymin;// + rm;
            z = k*scans + bounds.zmin;
//cout << x << " "<< y << " "<< z << " \n";
            for(n=0; n<natom; n++)
                {
                xx = x - xcart[n][0];
                yy = y - xcart[n][1];
                zz = z - xcart[n][2];
                if(xx > sizex05) xx=xx-sizex;
    	        if(xx < -sizex05) xx=xx+sizex;
                if(yy > sizey05) yy=yy-sizey;
	            if(yy < -sizey05) yy=yy+sizey;
                if(zz > sizez05) zz=zz-sizez;
	            if(zz < -sizez05) zz=zz+sizez;
                dsqr = xx*xx + yy*yy + zz*zz;
//                if(dsqr < refsqr)cout << sqrt(dsqr)<<endl;
                if(dsqr < dsqrmin) dsqrmin = dsqr;
                }
//cout << sqrt(dsqrmin)<<endl;

//цикл для отсечки пор - доделать, но получается по ресурсам не очень, лучше переработать алгоритм поиска пор -
//путем объеденения точек и нахождения их центра и радиуса

            if(dsqrmin > refsqr)
                {
cout <<"\nfindpore: dsqrmin (Bohr^2) = "<<dsqrmin<<endl;
                
            //ii++;
           //Здесь можно оставить natom, тогда в для следующего набора hkl будет учитываться атом углерода и атомы
           // не будут друг на друга садиться. А можно  поставить natom2 - что позволит запихивать атомы углерода друг на друга
                xcart[natom][0]=x;
                xcart[natom][1]=y;
                xcart[natom][2]=z;
                typat[natom]=2;
                natom++;
                npore++;
               cout << x << " "<< y << " "<< z << " \n";
                } 
        
            } 
//natom=natom2;            
cout <<"\nfindpore: For pore radius (Bohr) = "<<ri<<endl;
cout <<"\nfindpore: Pores found = "<<npore<<endl;
} 



void CrystalCell::cpxyz(int nx,int ny,int nz) {
//Копирует ячейку по x,y,z nx,ny,nz раз. nx,ny,nz должны быть >=0!
int nn[3],nnc[3];nn[0] = nx;nn[1]=ny;nn[2]=nz;
int natom2,i,j,k;
natom2 = nx * ny * nz * natom;
if(natom2 > max_at) cout <<" Error, too many atoms after copy";

for(j=0; j<3; j++) //Цикл для копирования по различным осям
    {
    nnc[0] = 0;nnc[1] = 0;nnc[2] = 0;//если равно 1, т
//    о по соответствующей оси будет прибавлен acell к кординатам скопированных атомов
    nnc[j] = 1;
    for(k=1; k<nn[j]; k++) //Цикл копирует ячейку по одной оси nn[j] раз
        for(i=0; i<natom; i++)
            { 
            xcart[natom * k + i][0] = xcart[i][0] + rprimd[0].x * k * nnc[0];
            xcart[natom * k + i][1] = xcart[i][1] + rprimd[1].y * k * nnc[1];
            xcart[natom * k + i][2] = xcart[i][2] + rprimd[2].z * k * nnc[2];
            typat[natom * k + i]    = typat[i];
            }
    natom = natom*nn[j];
    if(j == 0) rprimd[j].x = rprimd[j].x * nn[j];
    if(j == 1) rprimd[j].y = rprimd[j].y * nn[j];
    if(j == 2) rprimd[j].z = rprimd[j].z * nn[j];
    }
//natom = natom2;

}





void CrystalCell::delnear_by_xred() {
    //Need to be written. should work for any cell,
bool checked[max_at];
bool at_to_delete[max_at];
memset(checked, 0, natom);
memset(at_to_delete, 0, natom);
db sizex=acell[0], sizey=acell[1], sizez=acell[2];
db sizex05=sizex*0.5, sizey05=sizey*0.5, sizez05=sizez*0.5;
db rr,d,r = 3;
db xx,yy,zz;
int natom2 = natom;
int j=0;
for (int i=0;i<natom;i++)
    {
    checked[i]=1;
    for (j=0;j<natom;j++)
        {
        if(!checked[j])
            {
            xx = xcart[i][0] - xcart[j][0]; yy=xcart[i][1] - xcart[j][1]; zz=xcart[i][2] - xcart[j][2];
            if(xx > sizex05) xx=xx-sizex;
            if(xx < -sizex05) xx=xx+sizex;
            if(yy > sizey05) yy=yy-sizey;
            if(yy < -sizey05) yy=yy+sizey;
            if(zz > sizez05) zz=zz-sizez;
            if(zz < -sizez05) zz=zz+sizez;
            rr=xx*xx+yy*yy+zz*zz;
            d=sqrt(rr);
            if(d < r)
                { //если расстояние между атомами меньше 3 Бор
//              cout<<"i="<<i<<" coord was:"<<x_cell[i]<<" "<<y_cell[i]<<" "<<z_cell[i]<<endl;
                xcart[i][0] = xx*0.5 + xcart[j][0];//Позволяет найти среднее положение, т
                xcart[i][1] = yy*0.5 + xcart[j][1];
                xcart[i][2] = zz*0.5 + xcart[j][2];
                at_to_delete[j] = 1;
                cout <<"atom "<<j<<" was deleted"<<endl;
//              natom2--;
//  cout<<"i="<<i<<" coord become:"<<x_cell[i]<<" "<<y_cell[i]<<" "<<z_cell[i]<<endl;
//cout<<"atom j="<<j<<" will be deleted, because it is to close to atom i="<<i<<" ,d="<<d<<endl;
//cout<<"j  coord:"<<x_cell[j]<<" "<<y_cell[j]<<" "<<z_cell[j]<<endl;
//cout<<"i coord:"<<x_cell[i]<<" "<<y_cell[i]<<" "<<z_cell[i]<<endl;
                }
            }
        }
    }

int i2=0;
for (int i=0;i<natom;i++) {
    if(at_to_delete[i]) continue;
    xcart[i2][0] = xcart[i][0]; 
    xcart[i2][1] = xcart[i][1]; 
    xcart[i2][2] = xcart[i][2];
    typat[i2] = typat[i]; 
    i2++;
}
//if(i2!=natom2) cout << "Error in distance\n";
natom = i2;
}







void CrystalCell::periodicfill(int nx, int ny, int nz) {
//Вычисляет новые xcart и natom после копирования ячейки periodcell[] раз.
    int nh = nx;
    int nk = ny;
    int nl = nz;
    if(nh < 0 || nk < 0 || nl < 0) cerr <<"Error periodcell[] are wrong!!!";
    double basis_dec[3];
    vectormath t_xcart, tx;
    int i, u, v, w;//Номера узлов решетки
    int i_at = 0; // счетчик по новым атомам
    int i_typat[max_at], i_sort[max_at][n_sort];
    db i_xred[max_at][3];
    for(i=0;i<natom;i++) {
        i_xred[i][0] = xred[i][0];
        i_xred[i][1] = xred[i][1];
        i_xred[i][2] = xred[i][2];
        i_typat[i] = typat[i];
        i_sort[i][0] = sort[i][0];
        i_sort[i][1] = sort[i][1];
    }
    //Calculate rprimd
    vectormath new_rprimd[3];
    new_rprimd[0] = rprimd[0] * nx;
    new_rprimd[1] = rprimd[1] * ny;
    new_rprimd[2] = rprimd[2] * nz;

    //Используя заданные максимальные размеры системы пройти
    //по всем атомам,определить номера узловых точек ячейки и рассчитать декартовы координаты всех атомов.
    for(u = 0;u < nx;u++)
        for(v = 0;v < ny;v++)
            for(w = 0;w < nz;w++)
                for(i=0;i<natom;i++) {
                    //Вычисление декартовых координат атомов базиса
                    basis_dec[0] = i_xred[i][0] * rprimd[0].x + i_xred[i][1] * rprimd[1].x + i_xred[i][2] * rprimd[2].x;
                    basis_dec[1] = i_xred[i][0] * rprimd[0].y + i_xred[i][1] * rprimd[1].y + i_xred[i][2] * rprimd[2].y;
                    basis_dec[2] = i_xred[i][0] * rprimd[0].z + i_xred[i][1] * rprimd[1].z + i_xred[i][2] * rprimd[2].z;
                    //Вычисление новых координат
                    t_xcart.x = u * rprimd[0].x + v * rprimd[1].x + w * rprimd[2].x + basis_dec[0];
                    t_xcart.y = u * rprimd[0].y + v * rprimd[1].y + w * rprimd[2].y + basis_dec[1];
                    t_xcart.z = u * rprimd[0].z + v * rprimd[1].z + w * rprimd[2].z + basis_dec[2];
                    
                    tx = calculate_xred(new_rprimd, t_xcart.x, t_xcart.y, t_xcart.z);

                    xcart[i_at][0] = t_xcart.x;
                    xcart[i_at][1] = t_xcart.y;
                    xcart[i_at][2] = t_xcart.z;
                    xred[i_at][0] = tx.x;
                    xred[i_at][1] = tx.y;
                    xred[i_at][2] = tx.z;
                    typat[i_at] = i_typat[i];
                    sort[i_at][0] = i_sort[i][0];
                    sort[i_at][1] = i_sort[i][1];


                    i_at++;
                    if(i_at > max_at) cerr <<"Error periodcell[] are too big! Encrease max_at at cpp and fortran";

                    //d[s+i]=sqrt(x[s+i]*x[s+i]+y[s+i]*y[s+i]+z[s+i]*z[s+i]);//расстояние из точки 0 до каждого узла (по сути атома)
                }
    natom = i_at;

    cout << "\n'periodicfill': New cell contains " << natom << " atoms.\n";

    rprimd[0] = new_rprimd[0];
    rprimd[1] = new_rprimd[1];
    rprimd[2] = new_rprimd[2];

}





void CrystalCell::return_atoms_to_cell() {
    //To make: take into account that xred can be more than 2
    db upb = 0.999, bob = -0.001;
    int returned =0;
    for(int i=0;i<natom;i++) {
        //cout <<xred[i][0]<<" "<<xred[i][1]<<" "<<xred[i][2]<<endl;
        for(int j=0;j<3;j++) {
            if(xred[i][j] <= bob) xred[i][j] = xred[i][j] - truncf(xred[i][j]) + 1;
            if(xred[i][j] > upb) xred[i][j] = xred[i][j] - truncf(xred[i][j]);
        }
        xcart[i][0] = xred[i][0] * rprimd[0].x + xred[i][1] * rprimd[1].x + xred[i][2] * rprimd[2].x;
        xcart[i][1] = xred[i][0] * rprimd[0].y + xred[i][1] * rprimd[1].y + xred[i][2] * rprimd[2].y;
        xcart[i][2] = xred[i][0] * rprimd[0].z + xred[i][1] * rprimd[1].z + xred[i][2] * rprimd[2].z;
        //cout <<xred[i][0]<<" "<<xred[i][1]<<" "<<xred[i][2]<<endl;
        //returned++;
    }
    //if(test == 1) cout << "\n " << returned << " atoms were returned to cell.\n";
}

void CrystalCell::shift_atoms(){
//void CrystalCell::shift_atoms(int grain, int sort_of_layer, vectormath r, db mul){
    // Формат во входном файле shiftatoms grain sort_of_layer r mul
    //vectormath r - вектор вдоль которого делается сдвиг, пока возможны варианты выбора только оси rprimd; в файле должно быть 1, 2 или 3
    //mul - определяет длину на которую нужно сдвинуть
    //grain - зерно, в котором нужно сдвинуть атомы; 0 - grain A; 1 - grain B
    //sort_of_layer - слой, который нужно сдвинуть;  0 - layer A, 1 - layer B, 2 - layer C 
    int grain         = shiftatoms[0];
    int sort_of_layer = shiftatoms[1];
    int k             = shiftatoms[2] - 1;
    vectormath r      = rprimd[k];

    db mul            = shiftatoms[3];

    r = r * mul;
    cout <<"\nShift_atoms(): Length of shift vector is " << abs(r);

    for(int i=0;i<natom;i++) {
        //cout << sort[i][1]<<endl;
        if(sort[i][0] == grain && sort[i][1] == sort_of_layer) {
            xcart[i][0] = r.x + xcart[i][0];
            xcart[i][1] = r.y + xcart[i][1];
            xcart[i][2] = r.z + xcart[i][2];
            //cout << "\natom shifted";

        }
    }
    //calculate xred;
    vectormath new_xred;
    for(int i=0;i<natom;i++) {    
        new_xred = calculate_xred(rprimd, xcart[i][0], xcart[i][1], xcart[i][2]);
        //cout <<  new_xred.x<<" "<<new_xred.y<<" "<<new_xred.z<<" "<<"\n";
        xred[i][0] = new_xred.x;
        xred[i][1] = new_xred.y;
        xred[i][2] = new_xred.z;
    }

}

 

void CrystalCell::shift_cell(){
    //direction - in which 0,1,2 direction to shift 
    //mul - what part of rpmid vector
    int direction = (int)shiftcell[0];
    db mul = shiftcell[1];
    vectormath r = rprimd[direction] * mul;
    if(test == 1) cout <<"\nShift_cell(): Length of shift vector is " << abs(r);

    for(int i=0;i<natom;i++) {
            xcart[i][0] = r.x + xcart[i][0];
            xcart[i][1] = r.y + xcart[i][1];
            xcart[i][2] = r.z + xcart[i][2];

    }
    gbpos = gbpos + r.x;


    //calculate xred;
    vectormath new_xred;
    for(int i=0;i<natom;i++) {    
        new_xred = calculate_xred(rprimd, xcart[i][0], xcart[i][1], xcart[i][2]);
        //cout <<  new_xred.x<<" "<<new_xred.y<<" "<<new_xred.z<<" "<<"\n";
        xred[i][0] = new_xred.x;
        xred[i][1] = new_xred.y;
        xred[i][2] = new_xred.z;
    }

}

db maxval(const int array[], int length) {
    int max, t;
    max = array[0];
    for(int i=0;i<length;i++) {
        t = array[i];
        if(t>max) max = t;
    }
    return max;
}
void CrystalCell::write_abinit(string add,  const CrystalCell &init  ){
    //Выводим файл геометрии для abinit
    ofstream out;
    string folder = "out/"+name+"/"+add+"/";
    system( ("mkdir -p "+folder).c_str() );


    out.open(( folder+name+add+"."+to_string(version)+".in.geo").c_str(), ios::out);
    out.precision(12);
    out << fixed;
    if (len_units == " Angstrom") out <<"len_units Angstrom\n";
    out << "\nversion "<<version;
    out << "\nhex_a "<< hex_a * to_angstrom <<" angstrom";
    out << "\nhex_c "<< hex_c * to_angstrom <<" angstrom";
    out << "\ngbpos " << gbpos;   
    out << "\nacell  1 1 1\n";
    out <<"rprim  ";
    for(int i=0;i<3;i++)
        out<<rprimd[i].x<<" "<<rprimd[i].y<<" "<<rprimd[i].z<<endl;
    out<<endl;
//Пишем xred
    out << "xred  ";
    for(int i=0;i<natom;i++)
        out << xred[i][0]<<" "<<xred[i][1]<<" "<<xred[i][2]<<" "<<"\n";
    out<<endl;

    out << "xcart  ";
    for(int i=0;i<natom;i++)
        out << xcart[i][0]<<" "<<xcart[i][1]<<" "<<xcart[i][2]<<" "<<"\n";
    out<<endl;
    int ntypat;
    ntypat = maxval(typat, natom);
    out<<"ntypat  "<< ntypat<<endl;
    out<<endl;
    out<<"typat  ";
    int step = 20;
    for(int i=0;i<natom;i++) {
        out<<typat[i]<<" ";
        if(i >= step-1) {
            out <<endl;
            step = step + 20;
        }
    }

    out<<endl;
    out<<"natom  "<<natom<<endl;
    out<<"\nznucl  ";
    for(int i=0;i<ntypat;i++)
        out<<znucl[i]<<" ";


    if(init.name != "" && version == 1) {
        //Write information about how to build this structure using this programm
        //Only if init exist and version == 1
        out<<"\n\n\n#BEGIN BUILD INFORMATION!!!\n";
        //out <<"\nbuild_name  "    <<init.name;
        out <<"\ncalctype "<<init.calctype;
        out <<"\na_c_conv "<<init.a_c_conv[0]<<" "<<init.a_c_conv[1]<<" "<<init.a_c_conv[2]<<" "<<init.a_c_conv[3];
        out <<"\nbuild_natom  "  <<init.natom<<endl;
        out <<"\nbuild_acell "<<init.acell[0]<<" "<<init.acell[1]<<" "<<init.acell[2];    
        out <<"\nbuild_rprim\n";
        for(int i=0;i<3;i++)
            out<<init.rprim[i][0]<<" "<<init.rprim[i][1]<<" "<<init.rprim[i][2]<<endl;

        out << "build_xred  \n";
        for(int i=0;i<init.natom;i++)
            out << init.xred[i][0]<<" "<<init.xred[i][1]<<" "<<init.xred[i][2]<<" "<<"\n";

        ntypat = maxval(init.typat, init.natom);
        out<<"\nbuild_ntypat  "<< ntypat;
        out<<"\nbuild_typat  ";
        step = 20;
        for(int i=0;i<init.natom;i++) {
            out<<init.typat[i]<<" ";
            if(i >= step-1) {
                out <<endl;
                step = step + 20;
            }
        }
        out<<"\nbuild_znucl  ";
        for(int i=0;i<ntypat;i++)
            out<<init.znucl[i]<<" ";
        
        out<<"\nperiodcell  ";
        for(int i=0;i<6;i++)
            out<<init.periodcell[i]<<" ";

        out<<"\nhkl1  "<<init.hkl1[0]<<" "<<init.hkl1[1]<<" "<<init.hkl1[2];    
        out<<"\nuvw1  "<<init.uvw1[0]<<" "<<init.uvw1[1]<<" "<<init.uvw1[2];   
        out<<"\nuvw2  "<<init.uvw2[0]<<" "<<init.uvw2[1]<<" "<<init.uvw2[2]; 
        out<<"\nuvw3  "<<init.uvw3[0]<<" "<<init.uvw3[1]<<" "<<init.uvw3[2];
        out<<"\nmul   "<<init.mul [0]<<" "<<init.mul [1]<<" "<<init.mul [2];  
        out<<"\n#END BUILD INFORMATION!!!\n";
    }

    out.close();

}


void CrystalCell::scale_hcp_acell(db a, db c) {
    //cout << a;
    acell[0] += acell[0] * a / 100;
    acell[1] += acell[1] * a / 100;
    acell[2] += acell[2] * c / 100; 
    cout <<"acell "<< acell[0]<<" " <<acell[1]<<" " <<acell[2]<<endl;
    hex_a = acell[0];
    hex_c = acell[2];
    for(int i=0; i<3; i++) {
        rprimd[i].x = rprim[i][0] * acell[i];
        rprimd[i].y = rprim[i][1] * acell[i];
        rprimd[i].z = rprim[i][2] * acell[i];
    }
}



void CrystalCell::center() {
    for(int i=0; i<3; i++) {
        ;
        //rprimd[i].x = rprim[i][0] * acell[i];
        //rprimd[i].y = rprim[i][1] * acell[i];
        //rprimd[i].z = rprim[i][2] * acell[i];
    }
    if(test == 1) cout <<"\nStart to center... "<<name<<endl;
    db xc=0, yc=0, zc=0;
    for(int i=0;i<natom;i++) {
        xc+=xcart[i][0];
        yc+=xcart[i][1];
        zc+=xcart[i][2];
    }
    xc = xc / natom;
    yc = yc / natom;
    zc = zc / natom;

    for(int i=0;i<natom;i++) {
        xcart[i][0]-=xc;
        xcart[i][1]-=yc;
        xcart[i][2]-=zc;
    }

}


void CrystalCell::cut_rectangular(db xl, db xr, db yl, db yr, db zl, db zr, string cuttyp) {
    //Cut cell by planes perpendicular to x, y and z.
    int i_new = 0;
    if(cuttyp == "xred") {
    if(test == 1) cout <<"\nStart to cut... "<<cuttyp<<" "<<xl<<" "<<xr<<endl;

        for(int i=0;i<natom;i++) {
            

            if (xred[i][0] < xl or xred[i][0] > xr or
                xred[i][1] < yl or xred[i][1] > yr or
                xred[i][2] < zl or xred[i][2] > zr) 
                continue;
                

            xred[i_new][0] = xred[i][0];
            xred[i_new][1] = xred[i][1];
            xred[i_new][2] = xred[i][2];
            typat[i_new]   = typat[i];
            
            i_new++;
           
        }
        natom = i_new;
        if(test == 1) cout <<"\nNew number of atoms is  "<<natom<<" "<<endl;

        for(int i=0;i<natom;i++) {
            //Вычисление декартовых координат атомов базиса
            xcart[i][0] = xred[i][0] * rprimd[0].x + xred[i][1] * rprimd[1].x + xred[i][2] * rprimd[2].x;
            xcart[i][1] = xred[i][0] * rprimd[0].y + xred[i][1] * rprimd[1].y + xred[i][2] * rprimd[2].y;
            xcart[i][2] = xred[i][0] * rprimd[0].z + xred[i][1] * rprimd[1].z + xred[i][2] * rprimd[2].z;
        }

    }
    else if (cuttyp == "xcart")

        if(test == 1) cout <<"\nStart to cut... "<<cuttyp<<" "<<xl<<" "<<xr<<" "<<yl<<" "<<yr<<endl;
        
        for(int i=0;i<natom;i++) {
            

            if (xcart[i][0] < xl or xcart[i][0] > xr or
                xcart[i][1] < yl or xcart[i][1] > yr or
                xcart[i][2] < zl or xcart[i][2] > zr) 
                continue;
                

            xcart[i_new][0] = xcart[i][0];
            xcart[i_new][1] = xcart[i][1];
            xcart[i_new][2] = xcart[i][2];
            typat[i_new]   = typat[i];
            
            i_new++;
           
        }
        natom = i_new;
        if(test == 1) cout <<"\nNew number of atoms is  "<<natom<<" "<<endl;

}



//Functions:

void create_scaled_versions(string inname, int uniform_scale = 0, db nva = 0, db nvc = 0, int nsteps = 0 ) {
//Function creates sets of structures with different lattice constants
//for manual search of equillibrium values.
//Input:
//    inname - name of base geometry
//    nva - value in percents; lattice parameter a will be changed in [-nva, +nva] range
//    nvc - the same for c parameter
//    nsteps - number of values between ranges
    // uniform_scale - if 1 then a and c increased uniformly - needed for volume fit


    cout <<"Start of function create_scaled_versions()"<<endl;
    ofstream out;
    CrystalCell hcp, grainA, bulk_cell;     vectormath r1, r2, r3;
    int version = 1; //number of start version
    std::vector<double> a, c;

    nsteps = 0; //default argument does not work for some reason
    //cout << nsteps;

    out.open("a_and_c", ios::out);


    hcp.readin( inname );
    cout <<"Inintial version is "<<hcp.version<<endl;

    if (nsteps == 0) { //Determine parameters from input file
        version = (int)hcp.a_c_conv[0]; 
        nva     = hcp.a_c_conv[1];
        nvc     = hcp.a_c_conv[2];
        nsteps  = (int)hcp.a_c_conv[3];
        //cout << nsteps;
    }
    
    a = linspace(-nva, nva, nsteps);
    c = linspace(-nvc, nvc, nsteps);
    //for (int i=0; i<nsteps;i++)
    //    cout << nsteps << " ";



    db i = 0, j=0;
    for( auto &a_i : a ) {
        i++; j = 0;
        for( auto &c_i : c ) {
            j++;
            if (uniform_scale && j != i) continue;
            out <<"a, c = "<<a_i<<" "<<c_i<<" percents"<<endl;
            if (j == i) out <<"version " <<version<<" corresponds to uniform extension of cell "<<endl;
            //cout <<"c_i "<<c_i<<endl;
            grainA = hcp;
            grainA.version = version;
            grainA.scale_hcp_acell(a_i, c_i);
            grainA.choose_new_rprimd(r1, r2, r3);
            grainA.rotate(r1, r2, r3);
            grainA.periodicfill_and_cut(r1, r2, r3,-0.0001,0.9999);
            grainA.write_abinit("grainA_s", hcp);
            grainA.write_xyz("grainA_s");            
            bulk_cell = grainA;
            bulk_cell.periodicfill(2,1,1);
            bulk_cell.write_abinit("bulk_s");
            version++;
        }
    }
    
    out.close();
    if(test == 1) cout << "\na_c_conv: "<<version<<" "<<nva<<" "<<nvc<<" "<<nsteps<<endl; 
    cout << "\nScaled versions of the bulk cell were succesfully created ";    
    




    
}


std::vector<double> linspace(double a, double b, int n) {
    std::vector<db> array;
    db step = (b-a) / (n-1);
    if (n == 1) step = 0;
    for(int i=0; i<n; i++) 
        array.push_back(a + step * i); 
    return array;
}




























void matr3inv(db aa[3][3], db ait[3][3]) {
    //Arguments ------------------------------------
    //arrays


    //Local variables-------------------------------
    //scalars
    db   dd, t0, t1, t2;

    // *************************************************************************

     t0 = aa[1][1] * aa[2][2] - aa[2][1] * aa[1][2];
     t1 = aa[2][1] * aa[0][2] - aa[0][1] * aa[2][2];
     t2 = aa[0][1] * aa[1][2] - aa[1][1] * aa[0][2];
     dd = 1. / (aa[0][0] * t0 + aa[1][0] * t1 + aa[2][0] * t2);
     ait[0][0] = t0 * dd;
     ait[1][0] = t1 * dd;
     ait[2][0] = t2 * dd;
     ait[0][1] = ( aa[2][0] * aa[1][2] - aa[1][0] * aa[2][2] ) * dd; 
     ait[1][1] = ( aa[0][0] * aa[2][2] - aa[2][0] * aa[0][2] ) * dd; 
     ait[2][1] = ( aa[1][0] * aa[0][2] - aa[0][0] * aa[1][2] ) * dd; 
     ait[0][2] = ( aa[1][0] * aa[2][1] - aa[2][0] * aa[1][1] ) * dd; 
     ait[1][2] = ( aa[2][0] * aa[0][1] - aa[0][0] * aa[2][1] ) * dd; 
     ait[2][2] = ( aa[0][0] * aa[1][1] - aa[1][0] * aa[0][1] ) * dd; 

    // for(int i=0;i<3;i++)
    //     for(int j=0;j<3;j++)
    //         cout<< ait[i][j] <<" ";

}


vectormath calculate_xred(vectormath rprimd[3],double x,double y, double z) {
    //rprimd[3] are vectors in which reduced coordinates are calculated. x,y,z are decart coordinates of atom
    //The function still have problems then several components are zero
    db L[3][3]; vectormath tx_xcart;
    db prec = 10000.0f;
    for(int i=0;i<3;i++) {
       L[i][0] =  rprimd[i].x;
       L[i][1] =  rprimd[i].y;
       L[i][2] =  rprimd[i].z;
    }
    //Start determination of xred
    // db a1 = -L[0][1] / 1. / L[0][0];
    // if(L[0][1]==0) a1 = 0;
    // db c = -L[0][2] / 1. / L[0][1];
    // if(L[0][2]==0) c = 0;
    // db a11 = a1 * L[1][0] + L[1][1];
    // db a12 = a1 * L[2][0] + L[2][1];
    // db c11 = c * L[1][1] + L[1][2];
    // db c12 = c * L[2][1] + L[2][2];
    // db A = -c11 / a11;
    // if(c11==0) A = 0;
    // db n12 = a1 * x + y;
    // db n23 = c * y + z;
    
    // vectormath xred_o;
    
    // xred_o.z = (n12 * A + n23)/(A * a12 + c12);
    // xred_o.y = (n12 - xred_o.z * a12) / a11;
    // xred_o.x = (x - xred_o.y * L[1][0] - xred_o.z * L[2][0]) / L[0][0];


    //Calculate xred using inversion of matrix
    db gprimd[3][3];
    vectormath xred;

    matr3inv(L, gprimd);
    
    // for(int i=0;i<3;i++)
    //     for(int j=0;j<3;j++)
    //         cout<< L[i][j] <<" ";
    // cout<<"\n";

    // for(int i=0;i<3;i++)
    //     for(int j=0;j<3;j++)
    //         cout<< gprimd[i][j] <<" ";
    
    xred.x = gprimd[0][0] * x + gprimd[0][1] * y + gprimd[0][2] * z;
    xred.y = gprimd[1][0] * x + gprimd[1][1] * y + gprimd[1][2] * z;
    xred.z = gprimd[2][0] * x + gprimd[2][1] * y + gprimd[2][2] * z;

    //if(test == 1) { cout << "\nReduced coord using old alg "; xred_o.print(); }
    //if(test == 1) { cout << "Reduced coord using matr3inv"; xred.print(); }





    //Для теста, проверка tx
    //if(test == 1) cout << "\nDecart coord is "; t_xcart.print();
    //if(test == 1) cout << "\nReduced coord is "; tx.print();
    tx_xcart.x = xred.x * rprimd[0].x + xred.y * rprimd[1].x + xred.z * rprimd[2].x;
    tx_xcart.y = xred.x * rprimd[0].y + xred.y * rprimd[1].y + xred.z * rprimd[2].y;
    tx_xcart.z = xred.x * rprimd[0].z + xred.y * rprimd[1].z + xred.z * rprimd[2].z;
    //if(test == 1) cout << "\nDecart coord after check is "; tx_xcart.print();
    //Округляется с точностью до prec, так после повторного вычисления xcart есть небольшая ошибка
    if(roundf(x * prec) / prec != roundf(tx_xcart.x * prec) / prec || \
       roundf(y * prec) / prec != roundf(tx_xcart.y * prec) / prec || \
       roundf(z * prec) / prec != roundf(tx_xcart.z * prec) / prec ) {
        cerr << "\nError! calculate_xred works bad.\n";
        exit(1);
    }







    return xred;
}


string read_plus(string inname, string token) {
    //nelements - number of elements to read
    std::ifstream infile( inname );
    std::stringstream buffer; buffer << infile.rdbuf(); string source = buffer.str(); //read file to string
    
    while(source.find("#") != std::string::npos) { // remove comments
        size_t Beg = source.find("#");
        source.erase(Beg, source.find("\n", Beg) - Beg);
    }
    buffer.str("");
    buffer << source;

    //cout <<source;
    
    //cout << instr;
    //cout << inname<<endl;
    //std::size_t found = instr.find("calctype");
    string word, token_value;
    //bool read = 1;
    while(buffer >> word)
        if (word == token) {
            buffer >> word; //get the next word after calc type
            //cout << word <<endl;
            token_value = word;
            break;
        }
    //exit(1);
    return token_value;
}