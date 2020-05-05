#include <iostream>
#include <string>
#include <vector>
#include <thread>
#include "defs.h"
#include "graph.h"
#include "coloring_rainbow.h"
#include "matrix.h"
#include "domination.h"


using namespace std;

int N_THREADS;

/**
 * Computes Independent r-rainbow domination numbers for Petersen graphs P(n,k)
 * @param k parameter of P(n,k)
 * @param r number of colors
 */
void main_rainbow(u K, u r) {

    // MATRIX FILENAMES
    auto mono_filename = filename_data("rainbow",r,"petersen",K,"mono");
    auto mod_filenames = filename_datas("rainbow", r, "petersen", K, "mod", 2*K );

    // GENERATE MONOGRAPHS
    Graph *monograph = PetersenMonograph(K, K); // one-side monograph
    Graph *duograph = PetersenMonograph(K, 2*K); // double monograph
    Graph *quatrograph = PetersenMonograph(K, 4*K); // four monographs
    Graph *pentagraph = PetersenMonograph(K, 5*K); // five monographs
    Graph *start_pentagraph[K]; // starting monographs (mod)

    for (int m = 0; m < K; m++)
        start_pentagraph[m] = PetersenMonograph(K, 5*K+m);

    vvvu ex_colorings_in, ex_colorings_out; // store coloring extensions
    vvu colorings; // coloring of monograph

    // compute colorings

    vvu ex_colorings = Partial_Independent_Rainbow_Domination_Colorings(monograph, vector<u>(), vector<u>(), r, CUBIC);
    colorings = Middle_Colorings(quatrograph, ex_colorings, ex_colorings, r, CUBIC, Partial_Independent_Rainbow_Domination_Colorings, ex_colorings_in, ex_colorings_out);
    cout << colorings.size() << endl;

    // compute matrices
    matrix * mat =  MatrixEx(duograph, pentagraph,colorings, ex_colorings_in, ex_colorings_out, r, CUBIC, GoodPIRDQ);
    matrix * start_mat[K];
    start_mat[0] = mat;

    cout << "Monomatrix: " << *mat << endl;
    for (int m=1; m < K; m++) {
        start_mat[m] = MatrixEx(duograph, start_pentagraph[m], colorings, ex_colorings_in, ex_colorings_out, r, CUBIC, GoodPIRDQ);
        cout << "Matrix mod " << m << ": " << *start_mat[m] << endl;
    }

    // compute matrix powers
    for (u m = 0; m < K; m++) {
        u n0 = monograph->n/2 + m, nm = monograph->n/2;
        compute_matrix_powers(*start_mat[m], *mat,
                              [n0, nm, m, K](int i) { return "P("+to_string(n0 + i*nm)+","+to_string(K)+")";},
                              [n0, nm, m, K, r](int i) { return I2RD_conjecture(n0 + i*nm, K, r);}
        );
    }
}




int main() {


#if _WIN32 || _WIN64
#if _WIN64
    cout << "System is 64 bit." << endl;
#else
    cout << "System is 32 bit" << endl;
#endif
#endif

    N_THREADS = std::thread::hardware_concurrency()/2;

    cout << "Number of supported threads: " << std::thread::hardware_concurrency() << endl;
    cout << "Number of used threads:      " << N_THREADS << endl;

    timer_start();
    main_rainbow(2, 2);

    cout << "Time: " <<timer_value() << endl;
    return 0;
}
