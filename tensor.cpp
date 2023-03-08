#include <iostream>
#include <string>
#include <random>
#include <math.h>
#include <fstream>

#include "dais_exc.h"
#include "tensor.h"

#define PI 3.141592654
#define FLT_MAX 3.402823466e+38F /* max value */
#define FLT_MIN 1.175494351e-38F /* min positive value */

using namespace std;

/**
 * Random Initialization
 *
 * Perform a random initialization of the tensor
 *
 * @param mean The mean
 * @param std  Standard deviation
 */
void Tensor::init_random(float mean, float std) {
    if (data) {

        std::default_random_engine generator;
        std::normal_distribution<float> distribution(mean, std);

        for (int i = 0;i < r;i++) {
            for (int j = 0;j < c;j++) {
                for (int k = 0;k < d;k++) {
                    this->operator()(i, j, k) = distribution(generator);
                }
            }
        }

    }
    else {
        throw(tensor_not_initialized());
    }
}

Tensor::Tensor() {
    data = nullptr;
}

Tensor::Tensor(int r, int c, int d, float v) {   //costruttore che inizializza tutta la matrice con il valore v
    init(r, c, d, v);
}

float& Tensor::operator()(int i, int j, int k) {        //i riga, j colonna , k dimensione
    if (i >= 0 and i < r and j >= 0 and j < c and k >= 0 and k < d) //verifico se è nel range
        return data[i][j][k];   //se è nel range ritorno l'indirizzo di memoria posizione 
    else
        throw (index_out_of_bound());     //se non è nel range lancio un'eccezione
}

float Tensor::operator()(int i, int j, int k) const {    //i riga, j colonna , k dimensione
    if (i >= 0 and i < r and j >= 0 and j < c and k >= 0 and k < d) //verifico se è nel range
        return data[i][j][k];   //se è nel range ritorno l'elemento corrispondente alla posizione
    else
        throw (index_out_of_bound());     //se non è nel range lancio un'eccezione
}

void Tensor::init(int r, int c, int d, float v) {

    this->d = d;  //assegno le dimensioni
    this->r = r;
    this->c = c;

    this->data = new float** [r];  //alloco la memoria di data dove verranno contenute le matrici
    //per scorrere la matrice uso tre for(uno per la dimensione ,uno per la riga ed uno, per la colonna)
    for (int i = 0; i < r; i++) {       //i per la dimensione
        this->data[i] = new float* [c];
        for (int j = 0; j < c; j++) {   //j per le righe
            this->data[i][j] = new float[d];
            for (int k = 0; k < d; k++) {   //k per le colonne
                data[i][j][k] = v;
            }
        }
    }
}

Tensor::~Tensor() { //distruttore

    for (int i = 0; i < r; i++) {       //i per la dimensione
        for (int j = 0; j < c; j++) {   //j per le righe
            delete[] data[i][j];       //dealloco riga per riga
        }
        delete[] data[i];  //dealloco la dimensione
    }
    delete[] data; //dealloco le dimensioni
}

Tensor::Tensor(const Tensor& that) {    //costruttore per copia

    this->r = that.r;
    this->c = that.c;
    this->d = that.d;
    data = new float** [r];
    for (int i = 0;i < r;i++) {       //i per la dimensione
        data[i] = new float* [c];
        for (int j = 0;j < c;j++) {   //j per le righe
            data[i][j] = new float[d];
            for (int k = 0;k < d;k++) {   //k per le colonne
                this->operator()(i, j, k) = that.operator()(i, j, k);
            }
        }
    }
}

int Tensor::depth() const {    //ritorna il numero delle dimensioni
    return this->d;
}
int Tensor::rows() const {     //ritorna il numero delle righe
    return this->r;
}

int Tensor::cols() const {     //ritorna il numero delle colonne
    return this->c;
}


/**
 * Operator overloading ==
 *
 * It performs the point-wise equality check between two Tensors.
 *
 * The equality check between floating points cannot be simply performed using the
 * operator == but it should take care on their approximation.
 *
 * This approximation is known as rounding (do you remember "Architettura degli Elaboratori"?)
 *
 * For example, given a=0.1232 and b=0.1233 they are
 * - the same, if we consider a rounding with 1, 2 and 3 decimals
 * - different when considering 4 decimal points. In this case b>a
 *
 * So, given two floating point numbers "a" and "b", how can we check their equivalence?
 * through this formula:
 *
 * a ?= b if and only if |a-b|<EPSILON
 *
 * where EPSILON is fixed constant (defined at the beginning of this header file)
 *
 * Two tensors A and B are the same if:
 * A[i][j][k] == B[i][j][k] for all i,j,k
 * where == is the above formula.
 *
 * The two tensors must have the same size otherwise throw a dimension_mismatch()
 *
 * @return returns true if all their entries are "floating" equal
 */

bool Tensor::operator==(const Tensor& rhs) const {
    if (r != rhs.r or c != rhs.c or d != rhs.d) { //se i due Tensor sono uguali fa la sottrazione altrimenti da un eccezione
        throw (dimension_mismatch());
    }

    bool out_bool = true;
    int i = 0;

    while (i < this->rows() and out_bool) {
        int j = 0;
        while (j < this->cols() and out_bool) {
            int k = 0;
            while (k < this->depth() and out_bool) {
                if ((abs(this->operator()(i, j, k) - rhs.operator()(i, j, k)) >= EPSILON)) //come sopra: a==b if and only if |a-b|<EPSILON
                    out_bool = false;
                ++k;
            }
            ++j;
        }
        ++i;
    }
    return out_bool;
}

/**
* Operator overloading +
*
* It performs the point-wise sum between two Tensors.
*
* result(i,j,k)=this(i,j,k)+rhs(i,j,k)
*
* The two tensors must have the same size otherwise throw a dimension_mismatch()
*
* @return returns a new Tensor containing the result of the operation
*/
Tensor Tensor::operator +(const Tensor& rhs)const {
    if (this->r != rhs.r or this->c != rhs.c or this->d != rhs.d) {
        throw (dimension_mismatch());
    }

    Tensor out(*this);
    for (int i = 0; i < r; i++) {
        for (int j = 0; j < c; j++) {
            for (int k = 0; k < d; k++) {
                out.operator()(i, j, k) += rhs.operator()(i, j, k);
            }
        }
    }
    return out;
}

Tensor Tensor::operator-(const Tensor& rhs) const {
    if (r != rhs.r or c != rhs.c or d != rhs.d) { //se i due Tensor sono uguali fa la sottrazione altrimenti da un eccezione
        throw (dimension_mismatch());
    }

    Tensor out(*this);  //Tensor da dare in output (utilizzando il copy constructor)

    for (int i = 0; i < r; i++) {
        for (int j = 0; j < c; j++) {
            for (int k = 0; k < d; k++) {
                out.operator()(i, j, k) -= rhs.operator()(i, j, k);
            }
        }
    }
    return out;
}

Tensor Tensor::operator*(const Tensor& rhs) const {
    if (r != rhs.r or c != rhs.c or d != rhs.d) { //se i due Tensor sono uguali fa la sottrazione altrimenti da un eccezione
        throw (dimension_mismatch());
    }

    Tensor out(*this);  //Tensor da dare in output

    for (int i = 0; i < r; i++) {
        for (int j = 0; j < c; j++) {
            for (int k = 0; k < d; k++) {
                out.operator()(i, j, k) *= rhs.operator()(i, j, k);
            }
        }
    }
    return out;
}

Tensor Tensor::operator/(const Tensor& rhs) const {
    if (r != rhs.r or c != rhs.c or d != rhs.d) { //se i due Tensor sono uguali fa la sottrazione altrimenti da un eccezione
        throw (dimension_mismatch());
    }

    Tensor out(*this);  //Tensor da dare in output

    for (int i = 0; i < r; i++) {
        for (int j = 0; j < c; j++) {
            for (int k = 0; k < d; k++) {
                if (rhs.operator()(i, j, k) != 0) {
                    out.operator()(i, j, k) /= rhs.operator()(i, j, k);
                }
                else {
                    throw (unknown_exception());
                }
            }
        }
    }
    return out;
}

/**
* Operator overloading +
*
* It performs the point-wise sum between a Tensor and a constant
*
* result(i,j,k)=this(i,j,k)+rhs
*
* @return returns a new Tensor containing the result of the operation
*/
Tensor Tensor::operator+(const float& rhs)const {

    Tensor out(*this);

    for (int i = 0; i < r; i++) {
        for (int j = 0; j < c; j++) {
            for (int k = 0; k < d; k++) {
                out.operator()(i, j, k) += rhs;
            }
        }
    }
    return out;
}

Tensor Tensor::operator-(const float& rhs) const {
    Tensor out(*this);  //Tensor da dare in output

    for (int i = 0; i < r; i++) {
        for (int j = 0; j < c; j++) {
            for (int k = 0; k < d; k++) {
                out.operator()(i, j, k) -= rhs;
            }
        }
    }

    return out;
}

Tensor Tensor::operator*(const float& rhs) const {
    Tensor out(*this);  //Tensor da dare in output

    for (int i = 0; i < r; i++) {
        for (int j = 0; j < c; j++) {
            for (int k = 0; k < d; k++) {
                out.operator()(i, j, k) *= rhs;
            }
        }
    }
    return out;
}

Tensor Tensor::operator/(const float& rhs) const {
    Tensor out(*this);  //Tensor da dare in output

    for (int i = 0; i < r; i++) {
        for (int j = 0; j < c; j++) {
            for (int k = 0; k < d; k++) {
                if (rhs != 0) {
                    out.operator()(i, j, k) /= rhs;
                }
                else {
                    throw (unknown_exception());
                }
            }
        }
    }
    return out;
}

float Tensor::getMin(int k) const {
    float min = this->operator()(0, 0, k);
    for (int j = 0;j < this->rows();j++) {    //uso la j per scorrere le righe
        for (int z = 0;z < this->cols();z++) {    //uso la z per scorrere le colonne
            if (min > this->operator()(j, z, k))
                min = this->operator()(j, z, k);
        }
    }
    return min;
}

float Tensor::getMax(int k) const {    //ritorna il numero più alto nella dimensione k
    float max = this->operator()(0, 0, k);
    for (int j = 0;j < this->rows();j++) {    //uso la j per scorrere le righe
        for (int z = 0;z < this->cols();z++) {    //uso la z per scorrere le colonne
            if (max < this->operator()(j, z, k))
                max = this->operator()(j, z, k);
        }
    }
    return max;
}

void Tensor::showSize() const {    //output le dimensioni: righe x colonne x dimensioni
    cout << this->rows() << "x" << this->cols() << "x" << this->depth() << endl;
}

void Tensor::clamp(float low, float high) {

    for (int i = 0; i < r; i++) {
        for (int j = 0; j < c; j++) {
            for (int k = 0; k < d; k++) {
                if (this->operator()(i, j, k) < low)
                    this->operator()(i, j, k) = low;
                if (this->operator()(i, j, k) > high)
                    this->operator()(i, j, k) = high;
            }
        }
    }
}

/**
* Tensor Rescaling
*
* Rescale the value of the tensor following this rule:
*
* newvalue(i,j,k) = ((data(i,j,k)-min(i))/(max(i)-min(i)))*new_max
*
* where max(i) and min(i) are the maximum and minimum value in the i-th channel.
*
* new_max is the new value for the maximum
*
* @param new_max New maximum vale
*/
void Tensor::rescale(float new_max) {

    for (int k = 0; k < this->d; k++) { //controllo ogni canale
        float max = this->getMax(k);
        float min = this->getMin(k);

        for (int i = 0; i < this->r; i++) {
            for (int j = 0; j < this->c; j++) {
                if (max == min)
                    this->operator()(i, j, k) = new_max;
                else
                    this->operator()(i, j, k) = (this->operator()(i, j, k) - min) / (max - min) * new_max;
            }
        }
    }
}

/**
 * Subset a tensor
 *
 * retuns a part of the tensor having the following indices:
 * row_start <= i < row_end
 * col_start <= j < col_end
 * depth_start <= k < depth_end
 *
 * The right extrema is NOT included
 *
 * @param row_start
 * @param row_end
 * @param col_start
 * @param col_end
 * @param depth_start
 * @param depth_end
 * @return the subset of the original tensor
 */
Tensor Tensor::subset(unsigned int row_start, unsigned int row_end, unsigned int col_start, unsigned int col_end, unsigned int depth_start, unsigned int depth_end)const {

    if ((this->rows() < row_start or row_start < 0) or (this->rows() < row_end or row_end < 0) or (this->cols() < col_start or col_start < 0) or (this->cols() < col_end or col_end < 0) or (this->depth() < depth_start or depth_start < 0) or (this->depth() < depth_end or depth_end < 0)) throw (index_out_of_bound()); //controllo gli estremi (se uno non va bene allora lancio un'eccezione)

    if (row_end < row_start or col_end < col_start or depth_end < depth_start) throw (index_out_of_bound()); //controllo che gli "end" siano maggiori degli "start"

    Tensor out(row_end - row_start, col_end - col_start, depth_end - depth_start, 0);

    for (int i = 0; i < row_end - row_start; i++)
        for (int j = 0; j < col_end - col_start; j++)
            for (int k = 0; k < depth_end - depth_start; k++)
                out(i, j, k) = this->operator()(i + row_start, j + col_start, k + depth_start);

    return out;
}

/**
    * Tensor padding
    *
    * Zero pad a tensor in height and width, the new tensor will have the following dimensions:
    *
    * (rows+2*pad_h) x (cols+2*pad_w) x (depth)
    *
    * @param pad_h the height padding
    * @param pad_w the width padding
    * @return the padded tensor
*/
Tensor Tensor::padding(int pad_h, int pad_w) const {
    Tensor out(this->rows() + (2 * pad_h), this->cols() + (2 * pad_w), this->depth(), 0);    //inizializzo il Tensor da dare in output con dimensioni : (rows+2*pad_h) x (cols+2*pad_w) x (depth)

    for (int i = pad_h; i < out.rows() - pad_h; ++i) {
        for (int j = pad_w; j < out.cols() - pad_w; ++j) {
            for (int z = 0; z < out.depth(); ++z) {
                out.operator()(i, j, z) = this->operator()(i - pad_h, j - pad_w, z);
            }
        }
    }
    return out;
}


/**
 * Operator overloading = (assignment)
 *
 * Perform the assignment between this object and another
 *
 * @return a reference to the receiver object
 */

Tensor& Tensor::operator=(const Tensor& other) {
    if (this == &other) return *this; // se stesso oggetto

    if (this->data) { // se tensore già inizializzato
        if (this->r != other.r or this->c != other.c or this->d != other.d) { // se dimensioni dei tensori diverse
            this->~Tensor();
            this->data = nullptr;
            init(other.r, other.c, other.d); //inizializza
        }

        // fa la copia dei tensori
        for (int z = 0; z < this->d; ++z)
            for (int i = 0; i < this->r; ++i)
                for (int j = 0; j < this->c; ++j)
                    this->operator()(i, j, z) = other.operator()(i, j, z);
    }
    else { // se tensore vuoto
        init(other.r, other.c, other.d); // inizializza
        for (int z = 0; z < this->d; ++z) // popola il tensore
            for (int i = 0; i < this->r; ++i)
                for (int j = 0; j < this->c; ++j)
                    this->operator()(i, j, z) = other.operator()(i, j, z);

    }
    return *this;
}

/**
 * Concatenate
 *
 * The function concatenates two tensors along a give axis
 *
 * Example: this is of size 10x5x6 and rhs is of 25x5x6
 *
 * if concat on axis 0 (row) the result will be a new Tensor of size 35x5x6
 *
 * if concat on axis 1 (columns) the operation will fail because the number
 * of rows are different (10 and 25).
 *
 * In order to perform the concatenation is mandatory that all the dimensions
 * different from the axis should be equal, other wise throw concat_wrong_dimension().
 *
 * @param rhs The tensor to concatenate with
 * @param axis The axis along which perform the concatenation
 * @return a new Tensor containing the result of the concatenation
 */

Tensor Tensor::concat(const Tensor& rhs, int axis) const {
    if (axis == 0) { //concatena righe
        if (this->c != rhs.c or this->d != rhs.d) throw concat_wrong_dimension();

        Tensor out(this->r + rhs.r, this->c, this->d);

        for (int i = 0; i < this->c; ++i) {
            for (int j = 0; j < this->d; ++j) {
                for (int k = 0; k < out.r; ++k) {
                    if (k < this->r)
                        out(k, i, j) = (*this)(k, i, j);
                    else
                        out(k, i, j) = rhs(k - this->r, i, j);
                }
            }
        }
        return out;

    }
    else if (axis == 1) { //concatena colonne
        if (this->r != rhs.r or this->d != rhs.d) throw concat_wrong_dimension();

        Tensor out(this->r, this->c + rhs.c, this->d);

        for (int i = 0; i < this->r; ++i) {
            for (int j = 0; j < this->d; ++j) {
                for (int k = 0; k < out.c; ++k) {
                    if (k < this->c)
                        out(i, k, j) = (*this)(i, k, j);
                    else
                        out(i, k, j) = rhs(i, k - this->c, j);
                }
            }
        }
        return out;

    }
    else if (axis == 2) { //concatena canali
        if (this->r != rhs.r or this->c != rhs.c) throw concat_wrong_dimension();

        Tensor out(this->r, this->c, this->d + rhs.d);

        for (int i = 0; i < this->r; ++i) {
            for (int j = 0; j < this->c; ++j) {
                for (int k = 0; k < out.d; ++k) {
                    if (k < this->d)
                        out(i, j, k) = (*this)(i, j, k);
                    else
                        out(i, j, k) = rhs(i, j, k - this->d);
                }
            }
        }
        return out;
    }
    else throw (unknown_exception());

    return *this;
}

/**
 * Convolution
 *
 * This function performs the convolution of the Tensor with a filter.
 *
 * The filter f must have odd dimensions and same depth.
 *
 * Remeber to apply the padding before running the convolution
 *
 * @param f The filter
 * @return a new Tensor containing the result of the convolution
 */


Tensor Tensor::convolve(const Tensor& f) const {
    if (this->depth() != f.depth()) throw (dimension_mismatch());

    if ((f.rows() % 2) == 0 or (f.cols() % 2) == 0) throw (filter_odd_dimensions());

    int pad_r = f.rows() / 2;
    int pad_c = f.cols() / 2;

    Tensor conv(this->padding(pad_r, pad_c));
    Tensor out(this->rows(), this->cols(), this->depth());

    // ciclo tensore allargato
    for (int i = 0; i < out.depth(); ++i) {										// ciclo profondità tensore
        for (int j = 0; j < out.rows(); ++j) {									// ciclo righe tensore allargato
            for (int k = 0; k < out.cols(); ++k) {								// ciclo colonne tensore allargato
                // ciclo il filtro
                for (int m = 0; m < f.rows(); ++m) {							// ciclo righe tensore filtro
                    for (int n = 0; n < f.cols(); ++n) {						// ciclo colonne tensore filtro
                        out(j, k, i) += conv(j + m, k + n, i) * f(m, n, i);	// aggiorno tensore destinazione
                    }
                }

            }
        }
    }
    return out;
}


/**
 * Operator overloading <<
 *
 * Use the overaloading of << to show the content of the tensor.
 *
 * You are free to chose the output format, btw we suggest you to show the tensor by layer.
 *
 * [..., ..., 0]
 * [..., ..., 1]
 * ...
 * [..., ..., k]
 */
ostream& operator<<(ostream& stream, const Tensor& obj) {
    for (int k = 0; k < obj.d; ++k) {
        stream << "[";
        for (int j = 0; j < obj.rows(); ++j) {
            if (j != 0) stream << " ";
            for (int i = 0; i < obj.cols(); ++i) {
                stream << obj.operator()(j, i, k);
                if (i != (obj.cols() - 1)) stream << " ";
            }
            stream << ",";
        }
        stream << " ch " << k << "]" << endl;
    }

    return stream;
}


void Tensor::read_file(string filename) {
    ifstream ifs{ filename }; // stream input file

    if (ifs.eof() or ifs.bad() or ifs.fail()) throw unable_to_read_file(); // se non ci sono righe fallisce

    int val;

    // preleva da file prima dimensione tensore
    ifs >> val; // try to read rows
    if (val <= 0 or ifs.fail() or ifs.eof()) throw unable_to_read_file(); // se non trova intero lancia eccezione
    int dim_r = val;

    // preleva da file seconda dimensione tensore
    ifs >> val; // try to read cols
    //se tensore di dimensione 1,1,1
    if (val <= 0 or ifs.fail() or ifs.eof()) {
        init(1, 1, 1, dim_r);
        return;
    }
    int dim_c = val;

    // preleva da file terza dimensione tensore
    ifs >> val; // try to read depth
    if (val <= 0 or ifs.fail() or ifs.eof()) throw unable_to_read_file(); // se non trova intero lancia eccezione
    int dim_d = val;

    init(dim_r, dim_c, dim_d);


    float field;

    // ciclo per leggere linee colore: data(x,x..,x)
    for (int i = 0; i < (this->r * this->c * this->d); ++i) {
        ifs >> field;
        if (ifs.fail() or ifs.eof()) throw unable_to_read_file(); // se non trova intero lancia eccezione

        int ir = (i / this->c) % this->r;
        int ic = i % this->c;
        int id = i / (this->c * this->r);

        this->operator()(ir, ic, id) = field; // mette il valore nel tensore
    }
    ifs >> field;
    if (not ifs.eof()) throw unable_to_read_file(); // se non trova intero lancia eccezione

}

/**
 * Write the tensor to a file
 *
 * Write the content of a tensor to a textual file.
 *
 * The file should have this structure: the first three lines provide the dimensions while
 * the following lines contains the actual data by channel.
 *
 * For example, a tensor of size 4x3x2 will have the following structure:
 * 4
 * 3
 * 2
 * data(0,0,0)
 * data(0,1,0)
 * data(0,2,0)
 * data(1,0,0)
 * data(1,1,0)
 * .
 * .
 * .
 * data(3,1,1)
 * data(3,2,1)
 *
 * if the file is not reachable throw unable_to_read_file()
 *
 */
void Tensor::write_file(string filename) {
    ofstream ofs{ filename };
    ofs << this->rows() << endl;
    ofs << this->cols() << endl;
    ofs << this->depth() << endl;

    for (int k = 0; k < this->depth(); ++k) {
        for (int i = 0; i < this->rows(); ++i) {
            for (int j = 0; j < this->cols(); ++j) {
                ofs << this->operator()(i, j, k) << endl;
            }
        }

    }
}