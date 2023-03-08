#include <iostream>
#include <string>

#include "dais_exc.h"
#include "tensor.h"
#include "libbmp.h"
#include "DAISGram.h"

using namespace std;

/**
 * Load a bitmap from file
 *
 * @param filename String containing the path of the file
 */
void DAISGram::load_image(string filename) {
	BmpImg img = BmpImg();

	img.read(filename.c_str());

	const int h = img.get_height();
	const int w = img.get_width();

	data = Tensor(h, w, 3, 0.0);

	for (int i = 0;i < img.get_height();i++) {
		for (int j = 0;j < img.get_width();j++) {
			data(i, j, 0) = (float)img.red_at(j, i);
			data(i, j, 1) = (float)img.green_at(j, i);
			data(i, j, 2) = (float)img.blue_at(j, i);
		}
	}
}

/**
 * Save a DAISGram object to a bitmap file.
 *
 * Data is clamped to 0,255 before saving it.
 *
 * @param filename String containing the path where to store the image.
 */
void DAISGram::save_image(string filename) {

	data.clamp(0.0, 255.0);

	BmpImg img = BmpImg(getCols(), getRows());

	img.init(getCols(), getRows());

	for (int i = 0;i < getRows();i++) {
		for (int j = 0;j < getCols();j++) {
			img.set_pixel(j, i, (unsigned char)data(i, j, 0), (unsigned char)data(i, j, 1), (unsigned char)data(i, j, 2));
		}
	}

	img.write(filename);

}

DAISGram::DAISGram() {}

DAISGram::~DAISGram() {}

/**
* Get rows
*
* @return returns the number of rows in the image
*/
int DAISGram::getRows() {
	return data.rows();
}

/**
* Get columns
*
* @return returns the number of columns in the image
*/
int DAISGram::getCols() {
	return data.cols();
}

/**
* Get depth
*
* @return returns the number of channels in the image
*/
int DAISGram::getDepth() {
	return data.depth();
}

/**
* Brighten the image
*
* It sums the bright variable to all the values in the image.
*
* Before returning the image, the corresponding tensor should be clamped in [0,255]
*
* @param bright the amount of bright to add (if negative the image gets darker)
* @return returns a new DAISGram containing the modified object
*/
DAISGram DAISGram::brighten(float bright) {
	DAISGram out;
	out.data = this->data + bright;
	out.data.clamp(0.0, 255.0);
	return out;
}

/**
* Create a grayscale version of the object
*
* A grayscale image is produced by substituting each pixel with its average on all the channel
*
* @return returns a new DAISGram containing the modified object
*/
DAISGram DAISGram::grayscale() {
	DAISGram grayscale{};
	grayscale.data = this->data;
	float value{};
	for (int i = 0; i < (*this).getRows(); i++)
		for (int j = 0; j < (*this).getCols(); j++) {
			value = 0.0;
			for (int k = 0; k < (*this).getDepth(); k++) {
				value += this->data(i, j, k);
			}
			for (int k = 0; k < (*this).getDepth(); k++) {
				grayscale.data(i, j, k) = value / (*this).getDepth();
			}
		}
	return grayscale;
}

/**
 * Create a Warhol effect on the image
 *
 * This function returns a composition of 4 different images in which the:
 * - top left is the original image
 * - top right is the original image in which the Red and Green channel are swapped
 * - bottom left is the original image in which the Blue and Green channel are swapped
 * - bottom right is the original image in which the Red and Blue channel are swapped
 *
 * The output image is twice the dimensions of the original one.
 *
 * @return returns a new DAISGram containing the modified object
 */
DAISGram DAISGram::warhol() {
	DAISGram out{};

	//creo 3 copie dell'originale
	Tensor top_right(this->data); //copy constructor
	Tensor bottom_left(this->data);
	Tensor bottom_right(this->data);

	// swapping channel color
	for (int i = 0; i < this->data.rows(); ++i) {
		for (int j = 0; j < this->data.cols(); ++j) {
			for (int k = 0; k < this->data.depth(); ++k) {
				if (k == 0) { //RED
					top_right(i, j, k) = this->data(i, j, k + 1); //metto verde
					bottom_right(i, j, k) = this->data(i, j, k + 2); //metto blu
					bottom_left(i, j, k) = this->data(i, j, k); //lascio rosso
				}
				else if (k == 1) { //GREEN
					top_right(i, j, k) = this->data(i, j, k - 1); //metto rosso
					bottom_right(i, j, k) = this->data(i, j, k); //lascio verde
					bottom_left(i, j, k) = this->data(i, j, k + 1); //metto blu
				}
				else if (k == 2) { //BLUE
					top_right(i, j, k) = this->data(i, j, k); //lascio blu
					bottom_right(i, j, k) = this->data(i, j, k - 2); //metto rosso
					bottom_left(i, j, k) = this->data(i, j, k - 1); //metto verde
				}
			}
		}
	}
	//creo Tensore per la parte alta
	Tensor top = this->data.concat(top_right, 1);

	//creo Tensore per la parte bassa
	Tensor bottom = bottom_left.concat(bottom_right, 1);

	//creo Tensore con parte alta e bassa concatenata
	Tensor outcome = top.concat(bottom, 0);

	out.data = outcome;

	return out;
}

/**
 * Sharpen the image
 *
 * This function makes the image sharper by convolving it with a sharp filter
 *
 * filter[3][3]
 *    0  -1  0
 *    -1  5 -1
 *    0  -1  0
 *
 * Before returning the image, the corresponding tensor should be clamped in [0,255]
 *
 * @return returns a new DAISGram containing the modified object
 */
DAISGram DAISGram::sharpen() {
	DAISGram out{};
	Tensor filter(3, 3, this->getDepth());

	for (int i = 0; i < 3; ++i) {
		for (int j = 0; j < 3; ++j) {
			for (int k = 0; k < this->getDepth(); ++k) {
				if (((i + j) % 2) != 0)
					filter(i, j, k) = -1;
				else if (i == 1 and j == 1)
					filter(i, j, k) = 5.0;
			}
		}
	}
	out.data = this->data.convolve(filter);
	out.data.clamp(0.0, 255.0);
	return out;
}

/**
 * Emboss the image
 *
 * This function makes the image embossed (a light 3D effect) by convolving it with an
 * embossing filter
 *
 * filter[3][3]
 *    -2 -1  0
 *    -1  1  1
 *     0  1  2
 *
 * Before returning the image, the corresponding tensor should be clamped in [0,255]
 *
 * @return returns a new DAISGram containing the modified object
 */
DAISGram DAISGram::emboss() {
	DAISGram out{};

	Tensor filter(3, 3, this->getDepth());

	for (int i = 0; i < 3; ++i) {
		for (int j = 0; j < 3; ++j) {
			for (int k = 0; k < this->getDepth(); ++k) {
				if (i == 1 and j == 1)
					filter(i, j, k) = 1;
				else
					filter(i, j, k) = -2 + i + j;
			}
		}
	}

	out.data = this->data.convolve(filter);
	out.data.clamp(0.0, 255.0);
	return out;
}

/**
 * Smooth the image
 *
 * This function remove the noise in an image using convolution and an average filter
 * of size h*h:
 *
 * c = 1/(h*h)
 *
 * filter[3][3]
 *    c c c
 *    c c c
 *    c c c
 *
 * @param h the size of the filter
 * @return returns a new DAISGram containing the modified object
 */
DAISGram DAISGram::smooth(int h) {
	DAISGram out{};

	Tensor filter(h, h, this->getDepth());
	float c = 1.0 / (h * h);

	for (int i = 0; i < h; ++i) {
		for (int j = 0; j < h; ++j) {
			for (int k = 0; k < this->getDepth(); ++k) {
				filter(i, j, k) = c;
			}
		}
	}

	out.data = this->data.convolve(filter);
	return out;
}

/**
* Edges of an image
*
* This function extract the edges of an image by using the convolution
* operator and the following filter
*
*
* filter[3][3]
* -1  -1  -1
* -1   8  -1
* -1  -1  -1
*
* Remeber to convert the image to grayscale before running the convolution.
*
* Before returning the image, the corresponding tensor should be clamped in [0,255]
*
* @return returns a new DAISGram containing the modified object
*/
DAISGram DAISGram::edge() { //MODIFICA all'infice (messo k al posto di d_), ho messo 8.0 al posto di 8

	DAISGram edge{};
	DAISGram gray = (*this).grayscale();

	Tensor filter(3, 3, (*this).getDepth(), -1);
	for (int k = 0; k < filter.depth(); k++) {
		filter(1, 1, k) = 8.0;
	}
	edge.data = gray.data.convolve(filter);
	edge.data.clamp(0.0, 255.0);
	return edge;
}

/**
 * Blend with anoter image
 *
 * This function generate a new DAISGram which is the composition
 * of the object and another DAISGram object
 *
 * The composition follows this convex combination:
 * results = alpha*this + (1-alpha)*rhs
 *
 * rhs and this obejct MUST have the same dimensions.
 *
 * @param rhs The second image involved in the blending
 * @param alpha The parameter of the convex combination
 * @return returns a new DAISGram containing the blending of the two images.
 */
DAISGram DAISGram::blend(const DAISGram& rhs, float alpha) {
	DAISGram out{};
	if (alpha > 1.0 or alpha < 0.0) throw unknown_exception();
	if (this->data.rows() != rhs.data.rows() or this->data.cols() != rhs.data.cols() or this->data.depth() != rhs.data.depth()) throw dimension_mismatch();

	Tensor mix(this->getRows(), this->getCols(), this->getDepth());

	for (int i = 0; i < this->getRows(); ++i) {
		for (int j = 0; j < this->getCols(); ++j) {
			for (int k = 0; k < this->getDepth(); ++k) {
				mix(i, j, k) = (alpha * this->data(i, j, k)) + ((1 - alpha) * rhs.data(i, j, k));
			}
		}
	}
	out.data = mix;
	return out;
}

/**
 * Green Screen
 *
 * This function substitutes a pixel with the corresponding one in a background image
 * if its colors are in the surrounding (+- threshold) of a given color (rgb).
 *
 * (rgb - threshold) <= pixel <= (rgb + threshold)
 *
 *
 * @param bkg The second image used as background
 * @param rgb[] The color to substitute (rgb[0] = RED, rgb[1]=GREEN, rgb[2]=BLUE)
 * @param threshold[] The threshold to add/remove for each color (threshold[0] = RED, threshold[1]=GREEN, threshold[2]=BLUE)
 * @return returns a new DAISGram containing the result.
 */
 /*
 *	[Parametri]
 *   - rgb[100, 180, 50]
 *   - threshold[60, 60, 50]
 */
DAISGram DAISGram::greenscreen(DAISGram& bkg, int rgb[], float threshold[]) {
	if (this->data.rows() != bkg.data.rows() or this->data.cols() != bkg.data.cols() or this->data.depth() != bkg.data.depth()) throw dimension_mismatch();

	DAISGram out{};

	out.data = this->data;

	for (int i = 0; i < this->data.rows(); ++i) {
		for (int j = 0; j < this->data.cols(); ++j) {
			if ((this->data(i, j, 0) >= (rgb[0] - threshold[0])) and
				(this->data(i, j, 0) <= (rgb[0] + threshold[0])) and
				(this->data(i, j, 1) >= (rgb[1] - threshold[1])) and
				(this->data(i, j, 1) <= (rgb[1] + threshold[1])) and
				(this->data(i, j, 2) >= (rgb[2] - threshold[2])) and
				(this->data(i, j, 2) <= (rgb[2] + threshold[2]))) {
				out.data(i, j, 0) = bkg.data(i, j, 0);
				out.data(i, j, 1) = bkg.data(i, j, 1);
				out.data(i, j, 2) = bkg.data(i, j, 2);
			}
		}
	}

	return out;
}

/**
 * Equalize
 *
 * Stretch the distribution of colors of the image in order to use the full range of intesities.
 *
 * See https://it.wikipedia.org/wiki/Equalizzazione_dell%27istogramma
 *
 * @return returns a new DAISGram containing the equalized image.
 */
DAISGram DAISGram::equalize() {
	Tensor T_distro_freq(256, 1, 3);		// tensore per funzione distribuzione di frequenza
	Tensor T_distro_freq_cumula(256, 1, 3);		// tensore per funzione distribuzione di frequenza cumulata
	Tensor T_distro_freq_cumula_scaled(256, 1, 3);	// tensore per funzione distribuzione val colore scalato

	// distribuzione di frequenza per canale
	for (int i = 0; i < data.rows(); ++i)
		for (int j = 0; j < data.cols(); ++j)
			for (int k = 0; k < data.depth(); ++k)
				// aggiunge uno al contatore ogni volta che nel canale trova un certo valore
				T_distro_freq(round(this->data(i, j, k)), 0, k) = T_distro_freq(round(this->data(i, j, k)), 0, k) + 1;

	// distribuzione di frequenza cumulativa per canale
	Tensor T_cumula_min(1, 1, 3); // tensore per i minimi cumulati di canale
	for (int k = 0; k < this->getDepth(); ++k) {
		T_distro_freq_cumula(0, 0, k) = T_distro_freq(0, 0, k); // inizializza il primo valore
		T_cumula_min(0, 0, k) = T_distro_freq_cumula(0, 0, k);
		for (int i = 1; i < 256; ++i) {
			// somma il valore di frequenza al cumulato precedente
			T_distro_freq_cumula(i, 0, k) = T_distro_freq(i, 0, k) + T_distro_freq_cumula(i - 1, 0, k);
			// trovo il minimo
			if (T_cumula_min(0, 0, k) > T_distro_freq_cumula(i, 0, k) and T_distro_freq_cumula(i, 0, k) > 0) T_cumula_min(0, 0, k) = T_distro_freq_cumula(i, 0, k);
		}
	}

	// distribuzione di frequenza cumulativa scalata per canale
	for (int k = 0; k < this->getDepth(); ++k)
		for (int i = 0; i < 256; ++i)
			// applico la funzione che rimodula il range dei valori colore
			T_distro_freq_cumula_scaled(i, 0, k) = round((T_distro_freq_cumula(i, 0, k) - T_cumula_min(0, 0, k)) / ((data.rows() * data.cols()) - T_cumula_min(0, 0, k)) * (256 - 1)); //256 livelli di grigio

	// aggiorno valori colore sul tensore di destinazione
	DAISGram out;
	out.data = Tensor(data.rows(), data.cols(), data.depth());

	// aggiorno tensore di output
	for (int i = 0; i < data.rows(); ++i)
		for (int j = 0; j < data.cols(); ++j)
			for (int k = 0; k < data.depth(); ++k)
				out.data(i, j, k) = T_distro_freq_cumula_scaled(round(this->data(i, j, k)), 0, k);


	return out;
}

/**
 * Generate Random Image
 *
 * Generate a random image from nois
 *
 * @param h height of the image
 * @param w width of the image
 * @param d number of channels
 * @return returns a new DAISGram containing the generated image.
 */
void DAISGram::generate_random(int h, int w, int d) {
	data = Tensor(h, w, d, 0.0);
	data.init_random(128.0, 50.0);
	data.rescale(255.0);
}