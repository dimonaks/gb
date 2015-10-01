//Not adding impurities at all
//Во входном можно пока давать только xred! координаты
//Таким должен быть periodcell 0 1 0 1 0 1, если ячейку не нужно транслировать
//координаты пор пишуться в текущий объект с типом 3.
//Единицы измерения?
//Сделать что то с max_at, потому что сейчас в фортране и си это разные константы

/*
Для увеличения производительности:
Сейчас calculate_xred вызывается для каждого атомам - нужно вызывать один раз для всех атомов передавая их сразу

Для увеличения числа атомов:
Перевести все большие массивы на тип vector, тем самым используя динамическое выделение памяти только там где это нужно.

Для читаемости и ошибкоустойчивости кода:

все таки лучше иметь xcart[0][0] чем xcart[0].x, так как для первого легче организовывать циклы, тем самым уменьшая
количество строк и ошибок.

*/

#include "gb.h"

extern void write_abinit(char *name, bool *at_to_delete, int n_at, \
		db *x_cell, db *y_cell, db *z_cell,db *acell, int n_d_at,int ntypat);
extern void write_xyz(char *name,int n_at,int *i_znucl_at,db *x,db *y, db *z);
void create_scaled_versions(string inname, int uniform_scale = 0, db nva=0, db nvc=0, int nsteps=0 );

extern vector_dec vector_dec_for_UVW(int d[3], db rprimd[3][3], int mul);




int main(int argc, char *argv[]){
    
    string inname = string(argv[1]);
    if(argc!=2) {
        cout<<"Not enouth input data! The format must be: ./program <name of input file>\n";return 0;
    }
    cout << "Name of input file = "<<inname<<endl;  

    CrystalCell incell; //unit cell, big cell    
    vectormath r1, r2, r3; 
    string calctype = "empty";
    //std::vector<db> shift_atoms;
    //shift_atoms = read_plus_db( inname, "shift_atoms", 4 ); //read type of calculation

    calctype = read_plus( inname, "calctype" ); //read type of calculation
    incell.readin( inname ); //read main parameters with parser
    //return 1;
    cout << "You have choosen "<< calctype <<endl;
    //Working version for creating one structure
    if( calctype == "simple" ) {
        CrystalCell hcp_ucell, hcp_bcell, grainA;
        hcp_ucell.readin( inname ); 
        //system( ("mkdir out/"+hcp_ucell.name).c_str() ); 
        hcp_ucell.outcell();
        hcp_ucell.choose_new_rprimd(r1, r2, r3);
        hcp_bcell = hcp_ucell;
        hcp_bcell.rotate(r1, r2, r3);
        hcp_bcell.periodicfill_and_cut(r1, r2, r3,-0.0001,0.9999);
        grainA = hcp_bcell;
        hcp_bcell.make_boundary(0);
        //hcp_bcell.delnear();
        //if(inname.find("T2") != string::npos) hcp_bcell.shift_atoms(1, 1, r3, -2./6.);//For T2!
        //hcp_bcell.shift_cell(); moved to make boundary
        //hcp_bcell.return_atoms_to_cell();

        hcp_bcell.periodicfill(1,2,1);
        //hcp_bcell.cpxyz(2,1,1);
        //hcp_bcell.outcell();
        hcp_bcell.write_xyz("gb");
        grainA.write_xyz("bulk");
        grainA.periodicfill(2,1,1);
        grainA.write_xyz("bulk");
        grainA.write_abinit("bulk");
        //cout <<truncf(-2.4);
    }
    else if (calctype == "bulk_scale")
        create_scaled_versions(inname);

    else if (calctype == "bulk_uniform_scale") //
        create_scaled_versions(inname, 1);

    else if (calctype == "gb_scale") {
        CrystalCell gbcell, grainA, grainA_sc, bulkcell; //unit cell, big cell   
        grainA = incell;
        grainA.choose_new_rprimd(r1, r2, r3);
        grainA.rotate(r1, r2, r3);
        grainA.periodicfill_and_cut(r1, r2, r3,-0.0001,0.9999);
        grainA.write_abinit("grainA_r");
        bulkcell = grainA;
        bulkcell.periodicfill(2,1,1);
        bulkcell.write_abinit("bulk_r");
        //create boundary with its initial volume
        gbcell = grainA;
        db v_gb_i = gbcell.make_boundary(0); int n_at_i = gbcell.natom; //initial volume of gb and number of atoms
        gbcell.write_xyz("initial");
        gbcell.write_abinit("initial");
        cout<<"\nInitial specific volume is "<<v_gb_i * to_angstrom <<len_units<<endl;

        //Read target volumes
        db v_min = incell.a_c_conv[1];
        db v_max = incell.a_c_conv[2];
        int nsteps = (int)incell.a_c_conv[3];
        //cout << nsteps;
        std::vector<double> volumes = linspace(v_min, v_max, nsteps);
        //db v_target = 0; //Angstrom
        
        //Make cycle by volumes
        db delta_v, v_gb_t; int n_at_t;
        for( auto &v_target : volumes ) {
            v_target = v_target / to_angstrom; // return to Bohr 
            delta_v = v_target - v_gb_i;
            grainA_sc = grainA;
            grainA_sc.rprimd[0].x = grainA_sc.rprimd[0].x + 2*v_target;
            grainA_sc.write_abinit("grainA_sc");            
            gbcell = grainA;
            v_gb_t = gbcell.make_boundary(delta_v); n_at_t = gbcell.natom; //cell with target volume of gb and resulted numаber of atoms
            cout<<"\nObtained target specific volume is "<< v_gb_t * to_angstrom<<len_units<<endl;;
            if (red_prec(v_gb_t) != red_prec(v_target) ) cout << "\nI cant create cell with target "<<red_prec(v_target) <<" gb specific volume; I get "<<red_prec(v_gb_t)<<"\n";
            if (n_at_t != n_at_i) cout << "\nWarning! You have asked for specific volume, which caused additional removing of atoms and conseqently increase of specific volume\n";
            
            //if(inname.find("T2") != string::npos) gbcell.shift_atoms(1, 1, r3, -2./6.);//For T2!
            gbcell.write_abinit("target");
            gbcell.periodicfill(1,4,1);
            gbcell.write_xyz("target");
            
            grainA.version++;
        cout <<"\nI created "<< grainA.version-1<<" versions\n";
        }
    }
    else if (calctype == "planes") {
            CrystalCell planes, scell; //unit cell, big cell   
            scell = incell;
            scell.choose_new_rprimd(r1, r2, r3);
            scell.rotate(r1, r2, r3);
            scell.periodicfill_and_cut(r1, r2, r3,-0.0001,0.9999);

            scell.center();            
            planes = scell;





            db xl = incell.sixvalues[0];
            db xr = incell.sixvalues[1];
            db yl = incell.sixvalues[2];
            db yr = incell.sixvalues[3];
            db zl = incell.sixvalues[4];
            db zr = incell.sixvalues[5];


            planes.cut_rectangular(xl, xr, yl, yr, zl, zr, "xcart");       
            //planes.center();
            string name = "("+ to_string(planes.hkl1[0]) + to_string(planes.hkl1[1]) + to_string(planes.hkl1[2])+")";
            name = name +"_["+to_string(planes.uvw3[0])+to_string(planes.uvw3[1])+to_string(planes.uvw3[2])+"]";
            
            //scell.write_abinit("supercell");
            //scell.write_xyz(name+"_sc");
            planes.write_xyz(name+"_cut");
    }
    else cerr << "\nError!, you did not set calctype variable\n ";
    




    cout << endl;
    return 0;
}//end main



